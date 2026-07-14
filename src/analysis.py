from __future__ import annotations

from dataclasses import asdict, dataclass, field
from html import escape
import argparse
import json
from pathlib import Path

import pandas as pd

from .alignment import LocalAlignment, align_best_global_orientation, align_best_orientation, global_align, reverse_complement
from .io_utils import SequenceRecord, load_backbone_sequence, load_candidates, load_samples, parse_sample_map


DNA_BASES = set("ACGT")
STOP_CODONS = {"TAA", "TAG", "TGA"}
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


@dataclass
class Thresholds:
    pass_identity: float = 99.5
    pass_ref_coverage: float = 98.0
    warning_identity: float = 95.0
    warning_ref_coverage: float = 85.0
    backbone_pass_identity: float = 99.9
    backbone_pass_coverage: float = 99.5
    backbone_warning_identity: float = 98.0
    backbone_warning_coverage: float = 95.0


@dataclass
class AnalysisResult:
    sample_id: str
    colony_name: str
    expected_orf: str
    best_matching_orf: str
    status: str
    observed_insert_length: int
    expected_orf_length: int
    length_difference: int
    percent_identity: float
    query_coverage: float
    reference_coverage: float
    orientation: str
    matches: int
    mismatches: int
    insertions: int
    deletions: int
    ambiguous_bases: int
    starts_with_atg: bool
    length_divisible_by_3: bool
    frameshift: bool
    internal_stop_codon_count: int
    expected_stop_codon_position: str
    observed_stop_codon_positions: str
    amino_acid_translation_summary: str
    notes: str
    verdict: str
    source_file: str
    extraction_method: str = "local_alignment"
    backbone_status: str = "NOT CHECKED"
    backbone_identity: float = 0.0
    backbone_coverage: float = 0.0
    observed_backbone_length: int = 0
    expected_backbone_length: int = 0
    backbone_length_difference: int = 0
    backbone_matches: int = 0
    backbone_mismatches: int = 0
    backbone_insertions: int = 0
    backbone_deletions: int = 0
    backbone_ambiguous_bases: int = 0
    backbone_notes: str = ""
    backbone_events: list[dict] = field(default_factory=list)
    backbone_runs: list[dict] = field(default_factory=list)
    backbone_alignment: dict = field(default_factory=dict)
    alignment: dict = field(default_factory=dict)
    events: list[dict] = field(default_factory=list)
    stops: list[dict] = field(default_factory=list)
    ends: dict = field(default_factory=dict)
    ranking: list[dict] = field(default_factory=list)
    runs: list[dict] = field(default_factory=list)


def run_analysis(
    candidate_root: str | Path,
    sample_root: str | Path,
    sample_map_path: str | Path,
    results_dir: str | Path,
    thresholds: Thresholds | None = None,
    backbone_path: str | Path | None = None,
    replacement_marker: str = "GGGCCCCCCCT",
    flank_length: int = 40,
) -> tuple[list[AnalysisResult], dict]:
    thresholds = thresholds or Thresholds()
    validate_thresholds(thresholds)
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    candidates = load_candidates(candidate_root)
    samples = load_samples(sample_root)
    sample_map = parse_sample_map(sample_map_path)
    backbone = load_backbone_sequence(backbone_path)
    missing_notes: list[str] = []
    if not candidates:
        raise RuntimeError(f"No candidate ORF sequences found under {candidate_root}")
    if not samples:
        raise RuntimeError(f"No Plasmidsaurus sample sequences found under {sample_root}")
    missing_orfs = sorted({meta["expected_orf"] for meta in sample_map.values()} - set(candidates))
    if missing_orfs:
        raise RuntimeError(f"sample map references ORFs absent from the candidate FASTA: {', '.join(missing_orfs)}")

    results: list[AnalysisResult] = []
    for sample_id, meta in sample_map.items():
        records = samples.get(sample_id, [])
        if not records:
            missing_notes.append(f"{sample_id}: no parseable sample sequence found")
            results.append(missing_result(sample_id, meta, "No parseable sample sequence found."))
            continue
        result = analyze_sample(records, candidates, meta, thresholds, backbone, replacement_marker, flank_length)
        results.append(result)

    metadata = {
        "candidate_count": len(candidates),
        "sample_count": len(samples),
        "missing_notes": missing_notes,
        "candidate_root": str(candidate_root),
        "sample_root": str(sample_root),
        "results_dir": str(results_dir),
        "backbone_path": str(backbone_path) if backbone_path else "",
        "backbone_loaded": backbone is not None,
        "replacement_marker": replacement_marker,
        "flank_length": flank_length,
        "thresholds": asdict(thresholds),
    }
    write_outputs(results, results_dir, missing_notes, metadata)
    (results_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return results, metadata


def analyze_sample(
    records: list[SequenceRecord],
    candidates: dict[str, SequenceRecord],
    meta: dict[str, str],
    thresholds: Thresholds,
    backbone: SequenceRecord | None = None,
    replacement_marker: str = "GGGCCCCCCCT",
    flank_length: int = 40,
) -> AnalysisResult:
    expected_orf = meta["expected_orf"]
    colony_name = meta["colony_name"]
    sample_id = meta["sample_id"]
    ranked: list[tuple[dict, LocalAlignment, SequenceRecord, dict | None]] = []
    for record in records:
        analysis_record = record
        extraction = None
        if backbone:
            extraction = extract_insert_from_backbone(record.sequence, backbone.sequence, replacement_marker, flank_length)
            if extraction:
                analysis_record = SequenceRecord(
                    id=f"{record.id}_extracted_insert",
                    name=f"{record.name} extracted insert",
                    sequence=extraction["insert"],
                    source_file=record.source_file,
                    description=f"insert extracted with {flank_length} bp backbone flanks from {backbone.name}",
                )
        for candidate_id, candidate in candidates.items():
            aln = (
                align_best_global_orientation(candidate_id, candidate.sequence, analysis_record.sequence)
                if extraction
                else align_best_orientation(candidate_id, candidate.sequence, analysis_record.sequence)
            )
            counts = count_alignment(candidate.sequence, aln)
            ref_cov = 100.0 if extraction else safe_pct(aln.ref_end - aln.ref_start + 1, len(candidate.sequence))
            aligned_obs = sum(1 for c in aln.obs_aln if c != "-")
            query_cov = safe_pct(aligned_obs, max(1, len(aln.observed_insert)))
            comparable = counts["matches"] + counts["mismatches"] + counts["ambiguous_bases"]
            identity = safe_pct(counts["matches"], comparable)
            ranked.append(
                (
                    {
                        "orf": candidate_id,
                        "record_id": record.id,
                        "orientation": aln.orientation,
                        "score": round(aln.score, 2),
                        "percent_identity": round(identity, 3),
                        "reference_coverage": round(ref_cov, 3),
                        "query_coverage": round(query_cov, 3),
                        "observed_insert_length": len(aln.observed_insert),
                        "source_file": record.source_file,
                        "extraction_method": "backbone_flank" if extraction else "local_alignment",
                    },
                    aln,
                    analysis_record,
                    extraction,
                )
            )
    ranked.sort(key=lambda x: (x[0]["score"], x[0]["percent_identity"], x[0]["reference_coverage"]), reverse=True)
    best_meta, best_aln, best_record, _ = ranked[0]
    selected = next((item for item in ranked if item[0]["orf"] == expected_orf), ranked[0])
    selected_meta, selected_aln, selected_record, selected_extraction = selected
    selected_candidate = candidates.get(expected_orf) or candidates[selected_meta["orf"]]
    counts = count_alignment(selected_candidate.sequence, selected_aln)
    events, runs = build_events(selected_candidate.sequence, selected_aln, sample_id, expected_orf)
    expected_stops = stop_codons(selected_candidate.sequence)
    observed_stops = stop_codons(selected_aln.observed_insert)
    internal_stops = [s for s in observed_stops if not s["terminal"]]
    net_indel = counts["insertions"] - counts["deletions"]
    frameshift = (len(selected_aln.observed_insert) % 3 != 0) or (net_indel % 3 != 0)
    comparable = counts["matches"] + counts["mismatches"] + counts["ambiguous_bases"]
    identity = safe_pct(counts["matches"], comparable)
    ref_cov = 100.0 if selected_meta.get("extraction_method") == "backbone_flank" else safe_pct(selected_aln.ref_end - selected_aln.ref_start + 1, len(selected_candidate.sequence))
    aligned_obs = sum(1 for c in selected_aln.obs_aln if c != "-")
    query_cov = safe_pct(aligned_obs, max(1, len(selected_aln.observed_insert)))
    status, notes = call_status(
        expected_orf=expected_orf,
        best_orf=best_meta["orf"],
        selected_meta=selected_meta,
        top_meta=best_meta,
        identity=identity,
        ref_cov=ref_cov,
        orientation=selected_aln.orientation,
        counts=counts,
        frameshift=frameshift,
        internal_stop_count=len(internal_stops),
        thresholds=thresholds,
        candidate_count=len(candidates),
    )
    backbone_qc = analyze_backbone(selected_extraction, thresholds, backbone is not None)
    if backbone_qc["status"] == "FAIL":
        status = "FAIL"
        notes = "; ".join(filter(None, [notes, f"backbone QC failed: {backbone_qc['notes']}"]))
    elif backbone_qc["status"] == "WARNING" and status == "PASS":
        status = "WARNING"
        notes = "; ".join(filter(None, [notes, f"backbone QC warning: {backbone_qc['notes']}"]))
    verdict = verdict_text(status, sample_id, colony_name, expected_orf, best_meta["orf"], identity, ref_cov, selected_aln.orientation, notes)
    verdict += f" Backbone QC {backbone_qc['status']}: {backbone_qc['identity']:.3f}% identity, {backbone_qc['coverage']:.3f}% coverage; {backbone_qc['notes']}."
    exp_aa = translate(selected_candidate.sequence)
    obs_aa = translate(selected_aln.observed_insert)
    aa_summary = f"expected {len(exp_aa)} aa, observed {len(obs_aa)} aa; common prefix {common_prefix(exp_aa, obs_aa)} aa"
    return AnalysisResult(
        sample_id=sample_id,
        colony_name=colony_name,
        expected_orf=expected_orf,
        best_matching_orf=best_meta["orf"],
        status=status,
        observed_insert_length=len(selected_aln.observed_insert),
        expected_orf_length=len(selected_candidate.sequence),
        length_difference=len(selected_aln.observed_insert) - len(selected_candidate.sequence),
        percent_identity=round(identity, 3),
        query_coverage=round(query_cov, 3),
        reference_coverage=round(ref_cov, 3),
        orientation=selected_aln.orientation,
        matches=counts["matches"],
        mismatches=counts["mismatches"],
        insertions=counts["insertions"],
        deletions=counts["deletions"],
        ambiguous_bases=counts["ambiguous_bases"],
        starts_with_atg=selected_aln.observed_insert[:3] == "ATG",
        length_divisible_by_3=(len(selected_aln.observed_insert) % 3 == 0),
        frameshift=frameshift,
        internal_stop_codon_count=len(internal_stops),
        expected_stop_codon_position=format_stop_positions(expected_stops),
        observed_stop_codon_positions=format_stop_positions(observed_stops),
        amino_acid_translation_summary=aa_summary,
        notes=notes,
        verdict=verdict,
        source_file=selected_record.source_file,
        extraction_method=str(selected_meta.get("extraction_method", "local_alignment")),
        backbone_status=backbone_qc["status"],
        backbone_identity=backbone_qc["identity"],
        backbone_coverage=backbone_qc["coverage"],
        observed_backbone_length=backbone_qc["observed_length"],
        expected_backbone_length=backbone_qc["expected_length"],
        backbone_length_difference=backbone_qc["length_difference"],
        backbone_matches=backbone_qc["matches"],
        backbone_mismatches=backbone_qc["mismatches"],
        backbone_insertions=backbone_qc["insertions"],
        backbone_deletions=backbone_qc["deletions"],
        backbone_ambiguous_bases=backbone_qc["ambiguous_bases"],
        backbone_notes=backbone_qc["notes"],
        backbone_events=backbone_qc["events"],
        backbone_runs=backbone_qc["runs"],
        backbone_alignment=backbone_qc["alignment"],
        alignment={
            "ref_start": selected_aln.ref_start,
            "ref_end": selected_aln.ref_end,
            "obs_start": selected_aln.obs_start,
            "obs_end": selected_aln.obs_end,
            "ref_aln": selected_aln.ref_aln,
            "obs_aln": selected_aln.obs_aln,
            "expected_sequence": selected_candidate.sequence,
            "observed_insert": selected_aln.observed_insert,
            "observed_record_id": selected_record.id,
        },
        events=events,
        stops=[{"sequence": "observed", **s} for s in observed_stops] + [{"sequence": "expected", **s} for s in expected_stops],
        ends={
            "observed_first60": selected_aln.observed_insert[:60],
            "expected_first60": selected_candidate.sequence[:60],
            "observed_last60": selected_aln.observed_insert[-60:],
            "expected_last60": selected_candidate.sequence[-60:],
            "observed_contig_start": selected_aln.obs_start,
            "observed_contig_end": selected_aln.obs_end,
        },
        ranking=[item[0] | {"rank": i + 1} for i, item in enumerate(ranked[:50])],
        runs=runs,
    )


def count_alignment(reference: str, aln: LocalAlignment) -> dict[str, int]:
    counts = {"matches": 0, "mismatches": 0, "insertions": 0, "deletions": 0, "ambiguous_bases": 0}
    for ref_base, obs_base in zip(aln.ref_aln, aln.obs_aln):
        if ref_base == "-":
            counts["insertions"] += 1
        elif obs_base == "-":
            counts["deletions"] += 1
        elif ref_base not in DNA_BASES or obs_base not in DNA_BASES:
            counts["ambiguous_bases"] += 1
        elif ref_base == obs_base:
            counts["matches"] += 1
        else:
            counts["mismatches"] += 1
    return counts


def extract_insert_from_backbone(observed: str, backbone: str, replacement_marker: str, flank_length: int = 40) -> dict | None:
    marker = replacement_marker.upper()
    marker_start = backbone.upper().find(marker)
    if marker_start < 0:
        raise RuntimeError(f"replacement marker {replacement_marker!r} was not found in the backbone sequence")
    left = backbone[:marker_start]
    right = backbone[marker_start + len(marker) :]
    if len(left) < flank_length or len(right) < flank_length:
        raise RuntimeError("backbone flanks are shorter than the requested flank length")
    forward_left = left[-flank_length:].upper()
    forward_right = right[:flank_length].upper()
    backbone_upper = backbone.upper()
    if backbone_upper.count(forward_left) != 1 or backbone_upper.count(forward_right) != 1:
        raise RuntimeError("backbone flanks are not unique; increase the flank length or choose a unique replacement marker")
    expected_backbone = (right + left).upper()
    orientations = [
        ("same-as-backbone", observed.upper()),
        ("reverse-complement-of-backbone", reverse_complement(observed)),
    ]
    candidates: list[dict] = []
    for orientation, oriented_observed in orientations:
        circular = oriented_observed + oriented_observed
        for found in extract_between_flanks(circular, len(oriented_observed), forward_left, forward_right):
            observed_backbone = circular[
                found["right_index"] : found["left_index"] + len(oriented_observed) + len(forward_left)
            ]
            backbone_aln = global_align("backbone-anchor-candidate", expected_backbone, observed_backbone, "forward")
            counts = count_alignment(expected_backbone, backbone_aln)
            variant_burden = counts["mismatches"] + counts["insertions"] + counts["deletions"] + counts["ambiguous_bases"]
            candidates.append({
                "insert": found["insert"],
                "sample_orientation": orientation,
                "left_flank_start": found["left_flank_start"],
                "right_flank_start": found["right_flank_start"],
                "observed_backbone": observed_backbone,
                "expected_backbone": expected_backbone,
                "selection_key": (
                    variant_burden,
                    abs(len(observed_backbone) - len(expected_backbone)),
                    -counts["matches"],
                    -backbone_aln.score,
                ),
            })
    if not candidates:
        return None
    candidates.sort(key=lambda item: item["selection_key"])
    best = candidates[0]
    tied = [item for item in candidates if item["selection_key"] == best["selection_key"]]
    distinct_ties = {(item["insert"], item["observed_backbone"]) for item in tied}
    if len(distinct_ties) > 1:
        raise RuntimeError(
            f"multiple flank pairs ({len(tied)}) have indistinguishable backbone scores; increase the flank length"
        )
    best.pop("selection_key", None)
    best["anchor_candidate_count"] = len(candidates)
    best["selection_notes"] = (
        f"selected the best full-backbone match from {len(candidates)} flank pair candidates"
        if len(candidates) > 1
        else "one flank pair candidate found"
    )
    return best


def extract_between_flanks(circular: str, observed_length: int, left_flank: str, right_flank: str) -> list[dict]:
    candidates: list[dict] = []
    left_start = 0
    while left_start < observed_length:
        left_index = circular.find(left_flank, left_start, observed_length + len(left_flank) - 1)
        if left_index < 0:
            break
        insert_start = left_index + len(left_flank)
        right_start = insert_start
        right_limit = left_index + observed_length
        while right_start < right_limit:
            right_index = circular.find(right_flank, right_start, right_limit + len(right_flank))
            if right_index < 0 or right_index >= right_limit:
                break
            insert = circular[insert_start:right_index]
            candidates.append({
                "insert": insert,
                "left_flank_start": left_index + 1,
                "right_flank_start": (right_index % observed_length) + 1,
                "left_index": left_index,
                "right_index": right_index,
            })
            right_start = right_index + 1
        left_start = left_index + 1
    return candidates


def analyze_backbone(extraction: dict | None, thresholds: Thresholds, reference_loaded: bool) -> dict:
    empty = {
        "status": "NOT CHECKED" if not reference_loaded else "FAIL",
        "identity": 0.0,
        "coverage": 0.0,
        "observed_length": 0,
        "expected_length": 0,
        "length_difference": 0,
        "matches": 0,
        "mismatches": 0,
        "insertions": 0,
        "deletions": 0,
        "ambiguous_bases": 0,
        "notes": "no backbone reference was supplied" if not reference_loaded else "insert-flanking anchors were not both found; full backbone could not be verified",
        "events": [],
        "runs": [],
        "alignment": {},
    }
    if not extraction:
        return empty

    expected = extraction["expected_backbone"]
    observed = extraction["observed_backbone"]
    aln = align_best_global_orientation("backbone", expected, observed)
    counts = count_alignment(expected, aln)
    comparable = counts["matches"] + counts["mismatches"] + counts["ambiguous_bases"]
    identity = safe_pct(counts["matches"], comparable)
    coverage = safe_pct(comparable, len(expected))
    events, runs = build_backbone_events(aln)
    notes: list[str] = []
    if identity < thresholds.backbone_warning_identity:
        notes.append(f"low backbone identity: {identity:.3f}%")
    elif identity < thresholds.backbone_pass_identity:
        notes.append(f"backbone identity below PASS threshold: {identity:.3f}%")
    if coverage < thresholds.backbone_warning_coverage:
        notes.append(f"incomplete backbone coverage: {coverage:.3f}%")
    elif coverage < thresholds.backbone_pass_coverage:
        notes.append(f"backbone coverage below PASS threshold: {coverage:.3f}%")
    if counts["mismatches"]:
        notes.append(f"{counts['mismatches']} backbone mismatch(es)")
    if counts["insertions"] or counts["deletions"]:
        notes.append(f"{counts['insertions']} insertion(s), {counts['deletions']} deletion(s) in backbone")
    if counts["ambiguous_bases"]:
        notes.append(f"{counts['ambiguous_bases']} ambiguous backbone base(s)")

    if identity < thresholds.backbone_warning_identity or coverage < thresholds.backbone_warning_coverage:
        status = "FAIL"
    elif (
        identity < thresholds.backbone_pass_identity
        or coverage < thresholds.backbone_pass_coverage
        or counts["mismatches"]
        or counts["insertions"]
        or counts["deletions"]
        or counts["ambiguous_bases"]
    ):
        status = "WARNING"
    else:
        status = "PASS"
    if not notes:
        notes.append("full backbone matches the reference outside the replacement marker")
    if extraction.get("anchor_candidate_count", 1) > 1:
        notes.append(extraction["selection_notes"])
    return {
        "status": status,
        "identity": round(identity, 3),
        "coverage": round(coverage, 3),
        "observed_length": len(observed),
        "expected_length": len(expected),
        "length_difference": len(observed) - len(expected),
        **counts,
        "notes": "; ".join(notes),
        "events": events,
        "runs": runs,
        "alignment": {
            "ref_aln": aln.ref_aln,
            "obs_aln": aln.obs_aln,
            "expected_sequence": expected,
            "observed_sequence": observed,
            "sample_orientation": extraction["sample_orientation"],
        },
    }


def build_backbone_events(aln: LocalAlignment) -> tuple[list[dict], list[dict]]:
    events: list[dict] = []
    statuses: list[tuple[int, str]] = []
    ref_pos = 0
    obs_pos = 0
    for ref_base, obs_base in zip(aln.ref_aln, aln.obs_aln):
        if ref_base != "-":
            ref_pos += 1
        if obs_base != "-":
            obs_pos += 1
        if ref_base == "-":
            status = "insertion"
            event_type = "insertion"
        elif obs_base == "-":
            status = "deletion"
            event_type = "deletion"
        elif ref_base not in DNA_BASES or obs_base not in DNA_BASES:
            status = "ambiguous"
            event_type = "ambiguous base"
        elif ref_base == obs_base:
            status = "match"
            event_type = "match"
        else:
            status = "mismatch"
            event_type = "mismatch"
        if ref_base != "-":
            statuses.append((ref_pos, status))
        if event_type != "match":
            events.append(
                {
                    "event_type": event_type,
                    "backbone_coordinate": max(1, ref_pos),
                    "observed_coordinate": obs_pos if obs_base != "-" else None,
                    "expected_base": ref_base,
                    "observed_base": obs_base,
                }
            )
    return events, compress_runs(statuses)


def build_events(reference: str, aln: LocalAlignment, sample_id: str, expected_orf: str) -> tuple[list[dict], list[dict]]:
    events: list[dict] = []
    statuses: list[tuple[int, str]] = []
    ref_pos = aln.ref_start - 1
    obs_pos = aln.obs_start - 1 if aln.orientation == "forward" else 0
    oriented_obs_pos = 0
    codon_obs_bases: dict[int, list[str]] = {}
    for ref_base, obs_base in zip(aln.ref_aln, aln.obs_aln):
        if ref_base != "-":
            ref_pos += 1
        if obs_base != "-":
            oriented_obs_pos += 1
            if aln.orientation == "forward":
                obs_pos += 1
            else:
                obs_pos = aln.obs_end - oriented_obs_pos + 1
        if ref_base != "-":
            codon_index = (ref_pos - 1) // 3 + 1
            if obs_base != "-":
                codon_obs_bases.setdefault(codon_index, []).append(obs_base)
        if ref_base == "-":
            anchor = max(1, ref_pos)
            events.append(event_dict(sample_id, expected_orf, "insertion", anchor, obs_pos, "-", obs_base, reference, aln.observed_insert))
        elif obs_base == "-":
            statuses.append((ref_pos, "deletion"))
            events.append(event_dict(sample_id, expected_orf, "deletion", ref_pos, None, ref_base, "-", reference, aln.observed_insert))
        elif ref_base not in DNA_BASES or obs_base not in DNA_BASES:
            statuses.append((ref_pos, "ambiguous"))
            events.append(event_dict(sample_id, expected_orf, "ambiguous base", ref_pos, obs_pos, ref_base, obs_base, reference, aln.observed_insert))
        elif ref_base == obs_base:
            statuses.append((ref_pos, "match"))
        else:
            statuses.append((ref_pos, "mismatch"))
            ev = event_dict(sample_id, expected_orf, "mismatch", ref_pos, obs_pos, ref_base, obs_base, reference, aln.observed_insert)
            ev.update(codon_context(reference, codon_obs_bases, ref_pos))
            events.append(ev)
    return events, compress_runs(statuses)


def event_dict(sample_id: str, expected_orf: str, event_type: str, exp_pos: int, obs_pos: int | None, exp: str, obs: str, reference: str, observed: str) -> dict:
    return {
        "sample_id": sample_id,
        "expected_orf": expected_orf,
        "event_type": event_type,
        "expected_coordinate": exp_pos,
        "observed_coordinate": obs_pos,
        "expected_base": exp,
        "observed_base": obs,
        "nearby_expected_context": context(reference, exp_pos),
        "nearby_observed_context": context(observed, obs_pos) if obs_pos else "",
        "codon_position": ((exp_pos - 1) % 3) + 1 if exp_pos else "",
    }


def codon_context(reference: str, codon_obs_bases: dict[int, list[str]], ref_pos: int) -> dict:
    codon_index = (ref_pos - 1) // 3 + 1
    start = (codon_index - 1) * 3
    expected_codon = reference[start : start + 3]
    observed_codon = "".join(codon_obs_bases.get(codon_index, []))
    out = {"codon_index": codon_index, "expected_codon": expected_codon, "observed_codon": observed_codon}
    if len(expected_codon) == 3 and len(observed_codon) == 3:
        exp_aa = translate_codon(expected_codon)
        obs_aa = translate_codon(observed_codon)
        out["amino_acid_consequence"] = f"{exp_aa}{codon_index}{obs_aa}" if exp_aa != obs_aa else "synonymous"
    return out


def compress_runs(statuses: list[tuple[int, str]]) -> list[dict]:
    if not statuses:
        return []
    runs: list[dict] = []
    start, status = statuses[0]
    end = start
    count = 1
    for pos, st in statuses[1:]:
        if st == status and pos == end + 1:
            end = pos
            count += 1
        else:
            runs.append({"start": start, "end": end, "status": status, "count": count})
            start, end, status, count = pos, pos, st, 1
    runs.append({"start": start, "end": end, "status": status, "count": count})
    return runs


def call_status(expected_orf: str, best_orf: str, selected_meta: dict, top_meta: dict, identity: float, ref_cov: float, orientation: str, counts: dict, frameshift: bool, internal_stop_count: int, thresholds: Thresholds, candidate_count: int) -> tuple[str, str]:
    notes: list[str] = []
    if expected_orf != best_orf:
        notes.append(f"best match is {best_orf}, not expected {expected_orf}")
    if orientation != "forward":
        notes.append("insert aligns in reverse-complement orientation")
    if ref_cov < thresholds.warning_ref_coverage:
        notes.append(f"strong truncation or partial match: {ref_cov:.2f}% reference coverage")
    if identity < thresholds.warning_identity:
        notes.append(f"low identity: {identity:.2f}%")
    if frameshift:
        notes.append("frameshift risk from length or net indel phase")
    if internal_stop_count:
        notes.append(f"{internal_stop_count} premature/internal stop codon(s)")
    if counts["insertions"] or counts["deletions"]:
        notes.append(f"{counts['insertions']} insertion(s), {counts['deletions']} deletion(s)")
    if counts["ambiguous_bases"]:
        notes.append(f"{counts['ambiguous_bases']} ambiguous aligned base(s)")
    top_margin = top_meta["score"] - selected_meta["score"]
    if expected_orf == best_orf and candidate_count > 1 and top_margin < 3 and len(notes) == 0:
        notes.append("best ORF margin is small; candidate distinction may be ambiguous")

    fail_conditions = [
        expected_orf != best_orf and top_margin >= 10,
        orientation != "forward",
        ref_cov < thresholds.warning_ref_coverage,
        identity < thresholds.warning_identity,
        frameshift,
        internal_stop_count > 0,
    ]
    if any(fail_conditions):
        return "FAIL", "; ".join(notes)
    warning_conditions = [
        expected_orf != best_orf,
        ref_cov < thresholds.pass_ref_coverage,
        identity < thresholds.pass_identity,
        counts["insertions"] > 0,
        counts["deletions"] > 0,
        counts["mismatches"] > 0,
        counts["ambiguous_bases"] > 0,
        top_margin < 3,
    ]
    if any(warning_conditions):
        if not notes:
            notes.append("mostly correct insert but not perfect under configured PASS thresholds")
        return "WARNING", "; ".join(notes)
    return "PASS", "expected ORF is the best forward full-length match with no serious variants"


def verdict_text(status: str, sample_id: str, colony: str, expected: str, best: str, identity: float, coverage: float, orientation: str, notes: str) -> str:
    action = "Keep as a candidate colony." if status == "PASS" else "Review manually before keeping." if status == "WARNING" else "Reject for this expected insert."
    return f"{status}: {sample_id} ({colony}) expected {expected}; best match {best}. Identity {identity:.2f}%, reference coverage {coverage:.2f}%, orientation {orientation}. {notes}. {action}"


def missing_result(sample_id: str, meta: dict[str, str], note: str) -> AnalysisResult:
    return AnalysisResult(
        sample_id=sample_id,
        colony_name=meta["colony_name"],
        expected_orf=meta["expected_orf"],
        best_matching_orf="",
        status="WARNING",
        observed_insert_length=0,
        expected_orf_length=0,
        length_difference=0,
        percent_identity=0,
        query_coverage=0,
        reference_coverage=0,
        orientation="unknown",
        matches=0,
        mismatches=0,
        insertions=0,
        deletions=0,
        ambiguous_bases=0,
        starts_with_atg=False,
        length_divisible_by_3=False,
        frameshift=False,
        internal_stop_codon_count=0,
        expected_stop_codon_position="",
        observed_stop_codon_positions="",
        amino_acid_translation_summary="",
        notes=note,
        verdict=f"WARNING: {note}",
        source_file="",
    )


def stop_codons(seq: str) -> list[dict]:
    stops: list[dict] = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        if codon in STOP_CODONS:
            codon_index = i // 3 + 1
            terminal = i + 3 >= len(seq) - (len(seq) % 3)
            stops.append({"codon_index": codon_index, "nt_start": i + 1, "nt_end": i + 3, "codon": codon, "terminal": terminal})
    return stops


def format_stop_positions(stops: list[dict]) -> str:
    return "; ".join(f"aa{item['codon_index']} nt{item['nt_start']}-{item['nt_end']} {item['codon']}" for item in stops)


def translate(seq: str) -> str:
    return "".join(translate_codon(seq[i : i + 3]) for i in range(0, len(seq) - 2, 3))


def translate_codon(codon: str) -> str:
    return CODON_TABLE.get(codon.upper(), "X")


def common_prefix(a: str, b: str) -> int:
    count = 0
    for x, y in zip(a, b):
        if x != y:
            break
        count += 1
    return count


def context(seq: str, pos: int | None, radius: int = 15) -> str:
    if not pos:
        return ""
    start = max(0, pos - radius - 1)
    end = min(len(seq), pos + radius)
    return seq[start:end]


def safe_pct(num: float, den: float) -> float:
    return 0.0 if den == 0 else 100.0 * num / den


def validate_thresholds(thresholds: Thresholds) -> None:
    if thresholds.warning_identity > thresholds.pass_identity or thresholds.warning_ref_coverage > thresholds.pass_ref_coverage:
        raise ValueError("insert WARNING floors cannot exceed PASS thresholds")
    if thresholds.backbone_warning_identity > thresholds.backbone_pass_identity or thresholds.backbone_warning_coverage > thresholds.backbone_pass_coverage:
        raise ValueError("backbone WARNING floors cannot exceed PASS thresholds")


def overview_dataframe(results: list[AnalysisResult]) -> pd.DataFrame:
    rows = []
    for result in results:
        row = asdict(result)
        for key in ["alignment", "events", "stops", "ends", "ranking", "runs", "backbone_events", "backbone_runs", "backbone_alignment"]:
            row.pop(key, None)
        rows.append(row)
    return pd.DataFrame(rows)


def write_outputs(results: list[AnalysisResult], results_dir: Path, missing_notes: list[str], metadata: dict | None = None) -> None:
    overview = overview_dataframe(results)
    overview.to_csv(results_dir / "summary.csv", index=False)
    if missing_notes:
        (results_dir / "missing_or_ambiguous_files.txt").write_text("\n".join(missing_notes) + "\n", encoding="utf-8")
    else:
        (results_dir / "missing_or_ambiguous_files.txt").write_text("No missing or unparseable expected samples detected.\n", encoding="utf-8")
    for result in results:
        prefix = results_dir / result.sample_id
        pd.DataFrame(result.events).to_csv(prefix.with_suffix(".events.tsv"), sep="\t", index=False)
        pd.DataFrame(result.backbone_events).to_csv(prefix.with_suffix(".backbone_events.tsv"), sep="\t", index=False)
        pd.DataFrame(result.stops).to_csv(prefix.with_suffix(".stops.tsv"), sep="\t", index=False)
        pd.DataFrame(result.ranking).to_csv(prefix.with_suffix(".ranking.tsv"), sep="\t", index=False)
        prefix.with_suffix(".report.html").write_text(sample_report_html(result), encoding="utf-8")
    try:
        from .visualization import dashboard_html

        (results_dir / "index.html").write_text(dashboard_html(results, metadata or {}), encoding="utf-8")
    except Exception as exc:
        (results_dir / "dashboard_error.txt").write_text(str(exc) + "\n", encoding="utf-8")


def sample_report_html(result: AnalysisResult) -> str:
    event_rows = "".join(
        f"<tr><td>{escape(str(e.get('event_type', '')))}</td><td>{e.get('expected_coordinate','')}</td><td>{e.get('observed_coordinate','')}</td><td>{escape(str(e.get('expected_base','')))}</td><td>{escape(str(e.get('observed_base','')))}</td><td>{escape(str(e.get('amino_acid_consequence','')))}</td></tr>"
        for e in result.events[:500]
    )
    backbone_event_rows = "".join(
        f"<tr><td>{escape(str(e.get('event_type', '')))}</td><td>{e.get('backbone_coordinate','')}</td><td>{e.get('observed_coordinate','')}</td><td>{escape(str(e.get('expected_base','')))}</td><td>{escape(str(e.get('observed_base','')))}</td></tr>"
        for e in result.backbone_events[:500]
    )
    return f"""<!doctype html>
<html><head><meta charset="utf-8"><title>{escape(result.sample_id)} insert report</title>
<style>body{{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;margin:24px;color:#20242c}}table{{border-collapse:collapse;width:100%}}td,th{{border-bottom:1px solid #d9dee8;padding:6px;text-align:left}}.badge{{font-weight:700}}</style></head>
<body><h1>{escape(result.sample_id)} · {escape(result.expected_orf)}</h1>
<p class="badge">{escape(result.status)}</p><p>{escape(result.verdict)}</p>
<table><tbody>
<tr><th>Colony</th><td>{escape(result.colony_name)}</td></tr>
<tr><th>Best matching ORF</th><td>{escape(result.best_matching_orf)}</td></tr>
<tr><th>Identity</th><td>{result.percent_identity}%</td></tr>
<tr><th>Reference coverage</th><td>{result.reference_coverage}%</td></tr>
	<tr><th>Observed insert length</th><td>{result.observed_insert_length}</td></tr>
	<tr><th>Backbone status</th><td>{escape(result.backbone_status)}</td></tr>
	<tr><th>Backbone identity</th><td>{result.backbone_identity}%</td></tr>
	<tr><th>Backbone coverage</th><td>{result.backbone_coverage}%</td></tr>
	<tr><th>Backbone length</th><td>{result.observed_backbone_length} / {result.expected_backbone_length} bp</td></tr>
	<tr><th>Backbone notes</th><td>{escape(result.backbone_notes)}</td></tr>
	<tr><th>Notes</th><td>{escape(result.notes)}</td></tr>
	</tbody></table>
	<h2>Variant Events</h2><table><thead><tr><th>Event</th><th>Expected nt</th><th>Observed nt</th><th>Expected</th><th>Observed</th><th>AA consequence</th></tr></thead><tbody>{event_rows or '<tr><td colspan="6">No events.</td></tr>'}</tbody></table>
	<h2>Backbone Events</h2><table><thead><tr><th>Event</th><th>Backbone nt</th><th>Observed nt</th><th>Expected</th><th>Observed</th></tr></thead><tbody>{backbone_event_rows or '<tr><td colspan="5">No backbone events.</td></tr>'}</tbody></table>
</body></html>"""


def main() -> None:
    parser = argparse.ArgumentParser(description="Run insert-level ORF analysis and write reports.")
    project_root = Path(__file__).resolve().parents[1]
    input_root = project_root / "input"
    parser.add_argument("--candidate-root", default=str(input_root))
    parser.add_argument("--sample-root", default=str(input_root / "sequencing_results"))
    parser.add_argument("--sample-map", default=str(input_root / "sample_map.tsv"))
    parser.add_argument("--backbone", default=str(input_root / "backbone.ape"))
    parser.add_argument("--replacement-marker", default="GGGCCCCCCCT")
    parser.add_argument("--flank-length", type=int, default=40)
    parser.add_argument("--results-dir", default=str(Path(__file__).resolve().parents[1] / "results"))
    args = parser.parse_args()
    results, metadata = run_analysis(
        args.candidate_root,
        args.sample_root,
        args.sample_map,
        args.results_dir,
        backbone_path=args.backbone,
        replacement_marker=args.replacement_marker,
        flank_length=args.flank_length,
    )
    print(f"Loaded {metadata['candidate_count']} candidate ORFs and {metadata['sample_count']} samples.")
    print(f"Wrote {len(results)} sample reports to {metadata['results_dir']}")
    print(overview_dataframe(results)[["sample_id", "expected_orf", "best_matching_orf", "status", "percent_identity", "backbone_status", "backbone_identity", "backbone_coverage"]].to_string(index=False))


if __name__ == "__main__":
    main()
