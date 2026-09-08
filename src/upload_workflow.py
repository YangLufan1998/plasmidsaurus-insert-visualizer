from __future__ import annotations

from dataclasses import dataclass
import io
import json
from pathlib import Path, PurePosixPath
import tempfile
import zipfile

from .analysis import AnalysisResult, Thresholds, run_analysis, write_outputs


MAX_ZIP_UNCOMPRESSED_BYTES = 250 * 1024 * 1024


@dataclass
class BatchArtifacts:
    results: list[AnalysisResult]
    metadata: dict
    summary_csv: bytes
    inserts_fasta: bytes
    results_zip: bytes


def run_uploaded_batch(
    candidate_file: tuple[str, bytes] | None,
    sample_map_file: tuple[str, bytes],
    backbone_file: tuple[str, bytes],
    sequencing_files: list[tuple[str, bytes]],
    thresholds: Thresholds,
    replacement_marker: str,
    flank_length: int,
) -> BatchArtifacts:
    with tempfile.TemporaryDirectory(prefix="insert-visualizer-") as temp_name:
        temp_root = Path(temp_name)
        candidate_root = temp_root / "candidates"
        sample_root = temp_root / "sequencing"
        results_dir = temp_root / "results"
        candidate_root.mkdir()
        sample_root.mkdir()

        candidate_path = write_upload(candidate_root, *candidate_file) if candidate_file else None
        sample_map_path = write_upload(temp_root, *sample_map_file)
        backbone_path = write_upload(temp_root, *backbone_file)
        for name, data in sequencing_files:
            if name.lower().endswith(".zip"):
                extract_zip(data, sample_root)
            else:
                write_upload(sample_root, name, data)

        results, metadata = run_analysis(
            candidate_root,
            sample_root,
            sample_map_path,
            results_dir,
            thresholds,
            backbone_path=backbone_path,
            replacement_marker=replacement_marker,
            flank_length=flank_length,
        )
        sanitize_output_paths(
            results,
            metadata,
            candidate_path.name if candidate_path else "not supplied",
            backbone_path.name,
        )
        write_outputs(results, results_dir, metadata.get("missing_notes", []), metadata)
        (results_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        return BatchArtifacts(
            results=results,
            metadata=metadata,
            summary_csv=(results_dir / "summary.csv").read_bytes(),
            inserts_fasta=(results_dir / "extracted_inserts.fasta").read_bytes(),
            results_zip=zip_directory(results_dir),
        )


def write_upload(directory: Path, name: str, data: bytes) -> Path:
    safe_name = Path(name.replace("\\", "/")).name
    if not safe_name or safe_name in {".", ".."}:
        raise ValueError("uploaded file has an invalid name")
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / safe_name
    counter = 2
    while target.exists():
        target = directory / f"{Path(safe_name).stem}_{counter}{Path(safe_name).suffix}"
        counter += 1
    target.write_bytes(data)
    return target


def extract_zip(data: bytes, destination: Path) -> None:
    try:
        archive = zipfile.ZipFile(io.BytesIO(data))
    except zipfile.BadZipFile as exc:
        raise ValueError("a sequencing-results ZIP file is invalid") from exc
    with archive:
        total_size = sum(item.file_size for item in archive.infolist() if not item.is_dir())
        if total_size > MAX_ZIP_UNCOMPRESSED_BYTES:
            raise ValueError("sequencing-results ZIP expands beyond the 250 MB limit")
        for item in archive.infolist():
            if item.is_dir():
                continue
            relative = PurePosixPath(item.filename)
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError("sequencing-results ZIP contains an unsafe path")
            relative_parts = [part for part in relative.parts if part not in {"", "."}]
            if not relative_parts:
                continue
            target = destination.joinpath(*relative_parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(archive.read(item))


def sanitize_output_paths(
    results: list[AnalysisResult], metadata: dict, candidate_name: str, backbone_name: str
) -> None:
    for result in results:
        result.source_file = Path(result.source_file).name if result.source_file else ""
        for row in result.ranking:
            if row.get("source_file"):
                row["source_file"] = Path(str(row["source_file"])).name
    metadata.update(
        {
            "candidate_root": candidate_name,
            "sample_root": "uploaded sequencing files",
            "results_dir": "browser session download",
            "backbone_path": backbone_name,
        }
    )


def zip_directory(directory: Path) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(directory.rglob("*")):
            if path.is_file():
                archive.write(path, path.relative_to(directory).as_posix())
    return buffer.getvalue()
