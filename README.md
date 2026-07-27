# Plasmidsaurus Insert Visualizer

A Streamlit web app for extracting inserts from Plasmidsaurus colony sequencing results and checking the plasmid backbone.

When candidate ORFs are supplied, the app also reports ORF identity, insert variants, frame, and stop codons.

## Use the Web App

Upload three required input groups:

1. Sample mapping TSV
2. Backbone ApE or GenBank file
3. Plasmidsaurus sequence files, or a ZIP of the result folder

Candidate ORF FASTA is optional. When it is omitted, the app runs in extraction-only mode and returns the insert sequences without ORF comparison.

Click **Analyze uploaded files**. Results can be downloaded as extracted-insert FASTA, summary CSV, or a complete report ZIP.

The sample mapping file must be tab-separated:

```text
sample_id	colony_name	expected_orf
sample_1	colony_1	HIT_000001
sample_2	colony_2	HIT_000002
```

Sample filenames or FASTA headers must contain matching IDs such as `sample_1`. Candidate FASTA headers must contain IDs such as `HIT_000001`.

In extraction-only mode, `expected_orf` may be omitted:

```text
sample_id	colony_name
sample_1	colony_1
sample_2	colony_2
```

## Run Locally

Python 3.10 or newer is recommended.

```bash
git clone git@github.com:YangLufan1998/plasmidsaurus-insert-visualizer.git
cd plasmidsaurus-insert-visualizer
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

Open [http://localhost:8501](http://localhost:8501).

For a local Conda environment instead, run
`conda env create -f environment-local.yml`. The nonstandard filename keeps
Streamlit Community Cloud on the faster `requirements.txt` installation path.

The command-line workflow is also available for files stored under `input/`:

```bash
python -m src.analysis
```

## Deploy on Streamlit Community Cloud

1. Sign in at [share.streamlit.io](https://share.streamlit.io) with GitHub.
2. Create an app from `YangLufan1998/plasmidsaurus-insert-visualizer`.
3. Select branch `main` and entrypoint `app.py`.
4. Deploy, then copy the generated `streamlit.app` URL to the lab resource page.

The GitHub repository may remain private. To let anyone use the app, set the deployed app to public in its Streamlit sharing settings.

## Data Privacy

Input data and generated results are not stored in GitHub. Web uploads are processed in a temporary server directory and removed after analysis. Results remain in the active browser session until cleared or the session ends.

For sensitive or unpublished data, use a private deployment or run the app locally.

## Tests

```bash
python -m unittest discover -s tests -v
```
