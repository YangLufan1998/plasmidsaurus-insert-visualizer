# Plasmidsaurus Insert Visualizer

A local Streamlit app for checking colony sequencing results against an expected insert and plasmid backbone.

The app reports:

- the best matching candidate ORF,
- insert mismatches, indels, coverage, frame, and stop codons,
- full-backbone identity and coverage,
- separate insert, backbone, and overall PASS/WARNING/FAIL results,
- interactive sequence tracks and downloadable tables.

## Install

Python 3.10 or newer is recommended.

```bash
git clone git@github.com:YangLufan1998/plasmidsaurus-insert-visualizer.git
cd plasmidsaurus-insert-visualizer
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Conda users can instead run:

```bash
conda env create -f environment.yml
conda activate insert_visualizer
```

## Prepare Input Files

Create the following files inside `input/`:

```text
input/
  sample_map.tsv
  candidate_orfs.fna
  backbone.ape
  sequencing_results/
```

### Sample map

`input/sample_map.tsv` must be tab-separated:

```text
sample_id	colony_name	expected_orf
sample_1	colony_1	HIT_000001
sample_2	colony_2	HIT_000002
```

### Candidate ORFs

`input/candidate_orfs.fna` must be a FASTA file. Each header must contain an ID such as `HIT_000001`.

### Backbone

`input/backbone.ape` must be an ApE or GenBank file containing the replacement marker. The default marker is:

```text
GGGCCCCCCCT
```

### Sequencing results

Place Plasmidsaurus FASTA, FASTQ, or GenBank files inside `input/sequencing_results/`. Filenames or sequence headers must contain sample IDs such as `sample_1`.

## Run the App

```bash
streamlit run app.py
```

Open [http://localhost:8501](http://localhost:8501). Replace input files for a new batch, then click **Refresh analysis**.

## Run Without the App

```bash
python -m src.analysis
```

Results are written to `results/`:

- `index.html`: interactive batch report,
- `summary.csv`: sample-level summary,
- `sample_*.report.html`: individual reports,
- `sample_*.events.tsv`: insert variants,
- `sample_*.backbone_events.tsv`: backbone variants.

## Backbone Check

The program finds the insert using sequence on both sides of the replacement marker. If a flank also appears inside the insert, all valid flank pairs are tested and the pair with the best complete-backbone match is selected.

The backbone check supports circular sequence rotation and reverse-complement sequencing output.

## Tests

```bash
python -m unittest discover -s tests -v
```

## Data Privacy

Input sequences and generated results are ignored by Git. They stay on the local computer and should not be committed to the repository.
