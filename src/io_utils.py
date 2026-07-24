from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable


SUPPORTED_EXTENSIONS = {".fa", ".fasta", ".fna", ".gb", ".gbk", ".fastq", ".fq", ".txt", ".ape"}
DNA_RE = re.compile(r"^[ACGTRYSWKMBDHVNacgtryswkmbdhvn\-\s]+$")
HIT_RE = re.compile(r"H[Ii][Tt][_-]?0*(\d+)")


@dataclass(frozen=True)
class SequenceRecord:
    id: str
    name: str
    sequence: str
    source_file: str
    description: str = ""

    @property
    def hit_id(self) -> str | None:
        return normalize_hit_id(self.name) or normalize_hit_id(self.id) or normalize_hit_id(self.description)


def normalize_hit_id(text: str | None) -> str | None:
    if not text:
        return None
    match = HIT_RE.search(text)
    if not match:
        return None
    return f"HIT_{int(match.group(1)):06d}"


def normalize_sample_id(text: str | None) -> str | None:
    if not text:
        return None
    match = re.search(r"sample[_-]?(\d+)", text, re.IGNORECASE)
    if match:
        return f"sample_{int(match.group(1))}"
    return None


def normalize_nt(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq).upper().replace("U", "T")


def parse_sample_map(path: str | Path) -> dict[str, dict[str, str]]:
    import pandas as pd

    table = pd.read_csv(path, sep="\t")
    required = {"sample_id", "colony_name", "expected_orf"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"sample map is missing columns: {', '.join(sorted(missing))}")
    out: dict[str, dict[str, str]] = {}
    for row in table.to_dict("records"):
        sample_id = normalize_sample_id(row["sample_id"]) or str(row["sample_id"])
        if not sample_id or sample_id.lower() == "nan":
            raise ValueError("sample map contains an empty sample_id")
        if sample_id in out:
            raise ValueError(f"sample map contains duplicate sample_id: {sample_id}")
        colony_name = str(row["colony_name"])
        expected_orf = normalize_hit_id(str(row["expected_orf"])) or str(row["expected_orf"])
        if expected_orf.lower() == "nan":
            raise ValueError(f"sample map contains an empty expected_orf for {sample_id}")
        out[sample_id] = {
            "sample_id": sample_id,
            "colony_name": colony_name,
            "expected_orf": expected_orf,
        }
    return out


def discover_sequence_files(root: str | Path) -> list[Path]:
    root = Path(root).expanduser()
    if root.is_file():
        return [root] if root.suffix.lower() in SUPPORTED_EXTENSIONS else []
    if not root.exists():
        return []
    files = [p for p in root.rglob("*") if p.is_file() and p.suffix.lower() in SUPPORTED_EXTENSIONS]
    return sorted(files)


def candidate_files(root: str | Path) -> list[Path]:
    files = discover_sequence_files(root)
    preferred = [p for p in files if "predicted_orf" in p.name.lower()]
    if preferred:
        return preferred
    return [p for p in files if "pza11" not in p.name.lower() and "backbone" not in p.name.lower()]


def sample_files(root: str | Path) -> list[Path]:
    files = discover_sequence_files(root)
    preferred = [p for p in files if "fasta-file" in str(p).lower() or p.suffix.lower() in {".fa", ".fasta", ".fna", ".fastq", ".fq", ".gb", ".gbk"}]
    named = [p for p in preferred if normalize_sample_id(p.name) or normalize_sample_id(p.stem)]
    return sorted(named or preferred)


def parse_sequence_file(path: str | Path) -> list[SequenceRecord]:
    path = Path(path).expanduser()
    ext = path.suffix.lower()
    if ext in {".fa", ".fasta", ".fna"}:
        return _parse_fasta(path)
    if ext in {".fastq", ".fq"}:
        return _parse_fastq(path)
    if ext in {".gb", ".gbk", ".ape"}:
        records = _parse_with_biopython(path, "genbank")
        if records:
            return records
        return _parse_fasta(path)
    if ext == ".txt":
        return _parse_text_sequence(path)
    return []


def load_backbone_sequence(path: str | Path | None) -> SequenceRecord | None:
    if not path:
        return None
    path = Path(path).expanduser()
    if not path.exists() or not path.is_file():
        return None
    records = parse_sequence_file(path)
    if records:
        return records[0]
    text = path.read_text(errors="ignore")
    match = re.search(r"ORIGIN\s*(.*?)//", text, re.IGNORECASE | re.DOTALL)
    if not match:
        return None
    seq = normalize_nt(match.group(1))
    if not seq:
        return None
    return SequenceRecord(id=path.stem, name=path.stem, sequence=seq, source_file=str(path), description=path.name)


def load_candidates(root: str | Path) -> dict[str, SequenceRecord]:
    candidates: dict[str, SequenceRecord] = {}
    for path in candidate_files(root):
        for record in parse_sequence_file(path):
            hit_id = record.hit_id
            if hit_id and hit_id not in candidates:
                candidates[hit_id] = SequenceRecord(
                    id=hit_id,
                    name=record.name,
                    sequence=record.sequence,
                    source_file=record.source_file,
                    description=record.description,
                )
    return dict(sorted(candidates.items()))


def load_samples(root: str | Path) -> dict[str, list[SequenceRecord]]:
    samples: dict[str, list[SequenceRecord]] = {}
    for path in sample_files(root):
        records = parse_sequence_file(path)
        sample_id = normalize_sample_id(path.name)
        for record in records:
            sid = normalize_sample_id(record.id) or normalize_sample_id(record.name) or sample_id
            if sid:
                samples.setdefault(sid, []).append(record)
    return dict(sorted(samples.items()))


def _parse_fasta(path: Path) -> list[SequenceRecord]:
    records: list[SequenceRecord] = []
    name = ""
    desc = ""
    chunks: list[str] = []
    for line in path.read_text(errors="ignore").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name and chunks:
                records.append(_record(path, name, desc, chunks))
            desc = line[1:].strip()
            name = desc.split()[0] if desc else path.stem
            chunks = []
        else:
            chunks.append(line)
    if name and chunks:
        records.append(_record(path, name, desc, chunks))
    return records


def _parse_fastq(path: Path) -> list[SequenceRecord]:
    records: list[SequenceRecord] = []
    lines = path.read_text(errors="ignore").splitlines()
    i = 0
    while i + 3 < len(lines):
        header = lines[i].strip()
        seq = lines[i + 1].strip()
        plus = lines[i + 2].strip()
        if header.startswith("@") and plus.startswith("+"):
            desc = header[1:].strip()
            name = desc.split()[0] if desc else path.stem
            records.append(_record(path, name, desc, [seq]))
            i += 4
        else:
            i += 1
    return records


def _parse_text_sequence(path: Path) -> list[SequenceRecord]:
    text = path.read_text(errors="ignore")
    if ">" in text:
        return _parse_fasta(path)
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    dna_lines = [line for line in lines if DNA_RE.match(line)]
    seq = normalize_nt("".join(dna_lines))
    if len(seq) < 50:
        return []
    return [SequenceRecord(id=path.stem, name=path.stem, sequence=seq, source_file=str(path), description=path.name)]


def _parse_with_biopython(path: Path, fmt: str) -> list[SequenceRecord]:
    try:
        from Bio import SeqIO
    except Exception:
        return []
    records: list[SequenceRecord] = []
    try:
        with path.open() as handle:
            for item in SeqIO.parse(handle, fmt):
                seq = normalize_nt(str(item.seq))
                if seq:
                    records.append(
                        SequenceRecord(
                            id=item.id or path.stem,
                            name=item.name or item.id or path.stem,
                            sequence=seq,
                            source_file=str(path),
                            description=item.description or "",
                        )
                    )
    except Exception:
        return []
    return records


def _record(path: Path, name: str, desc: str, chunks: Iterable[str]) -> SequenceRecord:
    return SequenceRecord(id=name, name=name, sequence=normalize_nt("".join(chunks)), source_file=str(path), description=desc)
