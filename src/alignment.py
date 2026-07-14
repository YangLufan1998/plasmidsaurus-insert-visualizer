from __future__ import annotations

from dataclasses import dataclass


COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


@dataclass
class LocalAlignment:
    ref_id: str
    orientation: str
    score: float
    ref_aln: str
    obs_aln: str
    ref_start: int
    ref_end: int
    obs_start: int
    obs_end: int
    observed_insert: str


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1].upper()


def align_best_orientation(ref_id: str, reference: str, observed: str) -> LocalAlignment:
    forward = local_align(ref_id, reference, observed, "forward")
    reverse = local_align(ref_id, reference, reverse_complement(observed), "reverse-complement")
    if reverse.score > forward.score:
        n = len(observed)
        reverse.obs_start, reverse.obs_end = n - reverse.obs_end + 1, n - reverse.obs_start + 1
        return reverse
    return forward


def align_best_global_orientation(ref_id: str, reference: str, observed: str) -> LocalAlignment:
    forward = global_align(ref_id, reference, observed, "forward")
    reverse = global_align(ref_id, reference, reverse_complement(observed), "reverse-complement")
    if reverse.score > forward.score:
        n = len(observed)
        reverse.obs_start, reverse.obs_end = n - reverse.obs_end + 1, n - reverse.obs_start + 1
        return reverse
    return forward


def global_align(ref_id: str, reference: str, observed_oriented: str, orientation: str) -> LocalAlignment:
    try:
        from Bio.Align import PairwiseAligner
    except Exception as exc:
        raise RuntimeError("Biopython is required for global alignment. Install with `pip install -r requirements.txt`.") from exc

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -4.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(reference, observed_oriented)[0]
    ref_aln, obs_aln = gapped_strings(reference, observed_oriented, alignment.coordinates)
    return LocalAlignment(
        ref_id=ref_id,
        orientation=orientation,
        score=float(alignment.score),
        ref_aln=ref_aln,
        obs_aln=obs_aln,
        ref_start=1,
        ref_end=len(reference),
        obs_start=1,
        obs_end=len(observed_oriented),
        observed_insert=observed_oriented,
    )


def local_align(ref_id: str, reference: str, observed_oriented: str, orientation: str) -> LocalAlignment:
    try:
        from Bio.Align import PairwiseAligner
    except Exception as exc:
        raise RuntimeError("Biopython is required for local alignment. Install with `pip install -r requirements.txt`.") from exc

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -4.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(reference, observed_oriented)[0]
    ref_aln, obs_aln = gapped_strings(reference, observed_oriented, alignment.coordinates)
    ref_start, ref_end, obs_start, obs_end = aligned_bounds(alignment.coordinates)
    observed_insert = observed_oriented[max(0, obs_start - 1) : obs_end]
    return LocalAlignment(
        ref_id=ref_id,
        orientation=orientation,
        score=float(alignment.score),
        ref_aln=ref_aln,
        obs_aln=obs_aln,
        ref_start=ref_start,
        ref_end=ref_end,
        obs_start=obs_start,
        obs_end=obs_end,
        observed_insert=observed_insert,
    )


def gapped_strings(reference: str, observed: str, coordinates) -> tuple[str, str]:
    ref_parts: list[str] = []
    obs_parts: list[str] = []
    coords = coordinates.tolist() if hasattr(coordinates, "tolist") else coordinates
    for i in range(len(coords[0]) - 1):
        r0, r1 = int(coords[0][i]), int(coords[0][i + 1])
        q0, q1 = int(coords[1][i]), int(coords[1][i + 1])
        r_len = r1 - r0
        q_len = q1 - q0
        if r_len and q_len:
            ref_parts.append(reference[r0:r1])
            obs_parts.append(observed[q0:q1])
        elif r_len:
            ref_parts.append(reference[r0:r1])
            obs_parts.append("-" * r_len)
        elif q_len:
            ref_parts.append("-" * q_len)
            obs_parts.append(observed[q0:q1])
    return "".join(ref_parts), "".join(obs_parts)


def aligned_bounds(coordinates) -> tuple[int, int, int, int]:
    coords = coordinates.tolist() if hasattr(coordinates, "tolist") else coordinates
    ref_points = [int(x) for x in coords[0]]
    obs_points = [int(x) for x in coords[1]]
    return min(ref_points) + 1, max(ref_points), min(obs_points) + 1, max(obs_points)
