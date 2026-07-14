from __future__ import annotations

import unittest
import random

from src.alignment import reverse_complement
from src.analysis import Thresholds, analyze_backbone, extract_insert_from_backbone


MARKER = "GGGCCCCCCCT"


def dna(seed: int, length: int) -> str:
    return "".join(random.Random(seed).choices("ACGT", k=length))


LEFT = dna(11, 80)
RIGHT = dna(29, 80)
INSERT = "ATG" + "GAA" * 20 + "TAA"
BACKBONE = LEFT + MARKER + RIGHT


def rotate(sequence: str, offset: int) -> str:
    return sequence[offset:] + sequence[:offset]


class BackboneQcTests(unittest.TestCase):
    def qc(self, plasmid: str) -> dict:
        extraction = extract_insert_from_backbone(plasmid, BACKBONE, MARKER, flank_length=12)
        self.assertIsNotNone(extraction)
        return analyze_backbone(extraction, Thresholds(), reference_loaded=True)

    def test_exact_rotated_plasmid_passes(self) -> None:
        plasmid = rotate(LEFT + INSERT + RIGHT, 37)
        qc = self.qc(plasmid)
        self.assertEqual(qc["status"], "PASS")
        self.assertEqual(qc["identity"], 100.0)
        self.assertEqual(qc["coverage"], 100.0)

    def test_reverse_complement_plasmid_passes(self) -> None:
        plasmid = reverse_complement(rotate(LEFT + INSERT + RIGHT, 51))
        qc = self.qc(plasmid)
        self.assertEqual(qc["status"], "PASS")

    def test_single_backbone_substitution_warns(self) -> None:
        right = RIGHT[:40] + ("A" if RIGHT[40] != "A" else "C") + RIGHT[41:]
        qc = self.qc(LEFT + INSERT + right)
        self.assertEqual(qc["status"], "WARNING")
        self.assertEqual(qc["mismatches"], 1)

    def test_large_backbone_deletion_fails(self) -> None:
        right = RIGHT[:35] + RIGHT[55:]
        qc = self.qc(LEFT + INSERT + right)
        self.assertEqual(qc["status"], "FAIL")
        self.assertEqual(qc["deletions"], 20)

    def test_repeated_right_flank_inside_insert_selects_exact_backbone(self) -> None:
        repeated_insert = INSERT[:18] + RIGHT[:12] + INSERT[18:]
        plasmid = reverse_complement(rotate(LEFT + repeated_insert + RIGHT, 43))
        extraction = extract_insert_from_backbone(plasmid, BACKBONE, MARKER, flank_length=12)
        self.assertIsNotNone(extraction)
        self.assertGreater(extraction["anchor_candidate_count"], 1)
        self.assertEqual(extraction["insert"], repeated_insert)
        qc = analyze_backbone(extraction, Thresholds(), reference_loaded=True)
        self.assertEqual(qc["status"], "PASS")
        self.assertEqual(qc["observed_length"], len(BACKBONE) - len(MARKER))


if __name__ == "__main__":
    unittest.main()
