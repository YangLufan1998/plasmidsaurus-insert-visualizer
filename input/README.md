# Local Input Directory

Add or replace these files for each Plasmidsaurus batch:

- `sample_map.tsv`: columns must be `sample_id`, `colony_name`, `expected_orf`.
- `candidate_orfs.fna`: predicted ORF FASTA with `HIT_000000`-style IDs in headers.
- `backbone.ape`: backbone ApE/GenBank file containing the replacement marker exactly once. The selected left and right flanks must be unique in the reference. Repeated flank matches in an insert are resolved using the complete-backbone alignment.
- `sequencing_results/`: sample FASTA/FASTQ/GenBank files. Filenames or headers should include `sample_1`, `sample_2`, etc.

Default replacement marker: `GGGCCCCCCCT`.

Run from the project root:

```bash
python -m src.analysis
```

Open `results/index.html` for the unified interactive report.

The output reports insert status and full-backbone status separately. A backbone PASS means all sequence outside the replacement marker matches the supplied backbone under the configured thresholds and contains no called mismatch or indel.

Files in this directory are ignored by Git, except for this README.
