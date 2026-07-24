from __future__ import annotations

import io
from pathlib import Path
import random
import tempfile
import unittest
import zipfile

from src.analysis import Thresholds
from src.upload_workflow import extract_zip, run_uploaded_batch


def dna(seed: int, length: int) -> str:
    return "".join(random.Random(seed).choices("ACGT", k=length))


class UploadWorkflowTests(unittest.TestCase):
    def test_uploaded_batch_returns_download_bundle(self) -> None:
        marker = "GGGCCCCCCCT"
        left = dna(101, 80)
        right = dna(202, 80)
        insert = "ATG" + "GAA" * 20 + "TAA"
        sequencing_zip = io.BytesIO()
        with zipfile.ZipFile(sequencing_zip, "w") as archive:
            archive.writestr("fasta-files/run_sample_1.fasta", f">sample_1\n{left}{insert}{right}\n")

        artifacts = run_uploaded_batch(
            ("candidate_orfs.fna", f">HIT_000001\n{insert}\n".encode()),
            ("sample_map.tsv", b"sample_id\tcolony_name\texpected_orf\nsample_1\tcolony_1\tHIT_000001\n"),
            ("backbone.ape", f">backbone\n{left}{marker}{right}\n".encode()),
            [("plasmidsaurus_results.zip", sequencing_zip.getvalue())],
            Thresholds(),
            marker,
            12,
        )

        self.assertEqual(len(artifacts.results), 1)
        self.assertEqual(artifacts.results[0].backbone_status, "PASS")
        self.assertEqual(artifacts.results[0].source_file, "run_sample_1.fasta")
        with zipfile.ZipFile(io.BytesIO(artifacts.results_zip)) as archive:
            names = set(archive.namelist())
        self.assertIn("summary.csv", names)
        self.assertIn("index.html", names)
        self.assertIn("sample_1.report.html", names)
        self.assertNotIn("candidate_orfs.fna", names)
        self.assertNotIn("run_sample_1.fasta", names)

    def test_zip_path_traversal_is_rejected(self) -> None:
        payload = io.BytesIO()
        with zipfile.ZipFile(payload, "w") as archive:
            archive.writestr("../outside.fasta", ">sample_1\nACGT\n")
        with self.assertRaisesRegex(ValueError, "unsafe path"):
            with tempfile.TemporaryDirectory() as temp_name:
                extract_zip(payload.getvalue(), Path(temp_name))

    def test_missing_candidate_fasta_extracts_inserts_only(self) -> None:
        marker = "GGGCCCCCCCT"
        left = dna(303, 80)
        right = dna(404, 80)
        insert = "ATG" + "CAA" * 18 + "TGA"

        artifacts = run_uploaded_batch(
            None,
            ("sample_map.tsv", b"sample_id\tcolony_name\nsample_1\tcolony_1\n"),
            ("backbone.ape", f">backbone\n{left}{marker}{right}\n".encode()),
            [("run_sample_1.fasta", f">sample_1\n{left}{insert}{right}\n".encode())],
            Thresholds(),
            marker,
            12,
        )

        result = artifacts.results[0]
        self.assertEqual(artifacts.metadata["analysis_mode"], "extraction_only")
        self.assertEqual(result.analysis_mode, "extraction_only")
        self.assertEqual(result.alignment["observed_insert"], insert)
        self.assertEqual(result.observed_insert_length, len(insert))
        self.assertEqual(result.backbone_status, "PASS")
        self.assertIn(f">sample_1 colony=colony_1 length={len(insert)}", artifacts.inserts_fasta.decode())
        self.assertIn(insert, artifacts.inserts_fasta.decode())
        summary = artifacts.summary_csv.decode()
        self.assertIn("extracted_insert_sequence", summary)
        self.assertNotIn("expected_orf_length", summary)
        with zipfile.ZipFile(io.BytesIO(artifacts.results_zip)) as archive:
            names = archive.namelist()
            self.assertIn("extracted_inserts.fasta", names)
            self.assertNotIn("sample_1.events.tsv", names)
            self.assertNotIn("sample_1.ranking.tsv", names)
            self.assertNotIn("dashboard_error.txt", names)


if __name__ == "__main__":
    unittest.main()
