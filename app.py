from __future__ import annotations

from html import escape

import pandas as pd
import streamlit as st

from src.analysis import Thresholds, overview_dataframe
from src.upload_workflow import BatchArtifacts, run_uploaded_batch
from src.visualization import alignment_blocks, backbone_track, coordinate_track, status_badge, summary_cards, translation_view


st.set_page_config(page_title="Insert Sequencing Visualizer", layout="wide")

st.html(
    """
<style>
body { color: #20242c; }
.block-container { padding-top: 1.2rem; }
.badge { border-radius: 999px; padding: 2px 8px; font-size: 12px; font-weight: 800; white-space: nowrap; }
.summary-grid { display:grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 10px; margin: 0 0 12px; }
.metric-card { border: 1px solid #d9dee8; border-radius: 8px; padding: 10px; background: #fbfcfe; min-height: 74px; }
.metric-value { font-size: 21px; font-weight: 760; line-height: 1.2; overflow-wrap: anywhere; }
.metric-name { color: #667085; font-size: 12px; margin-top: 4px; }
.verdict { border-left: 4px solid #d9dee8; padding: 12px 14px; background: #fbfcfe; border-radius: 6px; line-height: 1.45; }
.legend { display:grid; grid-template-columns: 1fr 1fr; gap: 6px 10px; font-size: 12px; color: #667085; margin: 8px 0 12px; }
.swatch { display:inline-block; width: 13px; height: 9px; border-radius: 2px; margin-right: 6px; }
.alignment-block { white-space: pre; overflow-x: auto; border: 1px solid #d9dee8; border-radius: 8px; padding: 10px; background: #fff; font-size: 12px; line-height: 1.45; }
@media (max-width: 900px) { .summary-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); } }
</style>
"""
)


def sidebar_controls() -> tuple[Thresholds, float, bool]:
    st.sidebar.title("Review settings")
    with st.sidebar.expander("Insert QC thresholds"):
        pass_identity = st.slider("PASS identity threshold (%)", 90.0, 100.0, 99.5, 0.1)
        pass_coverage = st.slider("PASS reference coverage (%)", 80.0, 100.0, 98.0, 0.5)
        warn_identity = st.slider("WARNING identity floor (%)", 50.0, 100.0, 95.0, 0.5)
        warn_coverage = st.slider("WARNING coverage floor (%)", 50.0, 100.0, 85.0, 0.5)
    with st.sidebar.expander("Backbone QC thresholds"):
        backbone_pass_identity = st.slider("Backbone PASS identity (%)", 95.0, 100.0, 99.9, 0.1)
        backbone_pass_coverage = st.slider("Backbone PASS coverage (%)", 90.0, 100.0, 99.5, 0.1)
        backbone_warn_identity = st.slider("Backbone WARNING identity floor (%)", 80.0, 100.0, 98.0, 0.1)
        backbone_warn_coverage = st.slider("Backbone WARNING coverage floor (%)", 80.0, 100.0, 95.0, 0.5)
    with st.sidebar.expander("Track display"):
        zoom = st.slider("Track zoom", 0.7, 3.0, 1.15, 0.05)
        show_matches = st.checkbox("Show matching segments", True)
    st.sidebar.markdown(
        """
<div class="legend">
<div><span class="swatch" style="background:#2e9d63"></span>Match</div>
<div><span class="swatch" style="background:#d63b3b"></span>Mismatch</div>
<div><span class="swatch" style="background:#f29f05"></span>Deletion</div>
<div><span class="swatch" style="background:#7b61ff"></span>Insertion</div>
<div><span class="swatch" style="background:#0ea5e9"></span>Ambiguous</div>
<div><span class="swatch" style="background:#111827"></span>Stop codon</div>
</div>
""",
        unsafe_allow_html=True,
    )
    thresholds = Thresholds(
        pass_identity,
        pass_coverage,
        warn_identity,
        warn_coverage,
        backbone_pass_identity,
        backbone_pass_coverage,
        backbone_warn_identity,
        backbone_warn_coverage,
    )
    return thresholds, zoom, show_matches


def upload_panel(thresholds: Thresholds) -> BatchArtifacts | None:
    st.subheader("Upload a sequencing batch")
    left, right = st.columns(2)
    with left:
        candidate = st.file_uploader(
            "Candidate ORF FASTA (optional)",
            type=["fasta", "fa", "fna"],
            key="candidate_upload",
            help="Leave empty to extract inserts without ORF comparison.",
        )
        sample_map = st.file_uploader(
            "Sample mapping TSV",
            type=["tsv", "txt"],
            key="sample_map_upload",
        )
    with right:
        backbone = st.file_uploader(
            "Backbone ApE or GenBank file",
            type=["ape", "gb", "gbk", "fasta", "fa"],
            key="backbone_upload",
        )
        sequencing = st.file_uploader(
            "Plasmidsaurus results",
            type=["zip", "fasta", "fa", "fna", "fastq", "fq", "gb", "gbk", "txt"],
            accept_multiple_files=True,
            key="sequencing_upload",
        )
    settings_left, settings_right = st.columns(2)
    with settings_left:
        replacement_marker = st.text_input("Replacement marker", "GGGCCCCCCCT")
    with settings_right:
        flank_length = st.number_input("Backbone flank length", min_value=12, max_value=120, value=40, step=1)

    st.download_button(
        "Download sample-map template",
        data=b"sample_id\tcolony_name\texpected_orf\nsample_1\tcolony_1\tHIT_000001\n",
        file_name="sample_map_template.tsv",
        mime="text/tab-separated-values",
    )
    ready = sample_map is not None and backbone is not None and bool(sequencing)
    if st.button("Analyze uploaded files", type="primary", disabled=not ready, use_container_width=True):
        st.session_state.pop("batch_artifacts", None)
        try:
            with st.spinner("Analyzing inserts and full backbones..."):
                artifacts = run_uploaded_batch(
                    (candidate.name, candidate.getvalue()) if candidate else None,
                    (sample_map.name, sample_map.getvalue()),
                    (backbone.name, backbone.getvalue()),
                    [(item.name, item.getvalue()) for item in sequencing],
                    thresholds,
                    replacement_marker,
                    int(flank_length),
                )
            st.session_state["batch_artifacts"] = artifacts
        except Exception as exc:
            st.error(f"Analysis could not run: {exc}")
    st.caption("Uploads are processed in a temporary workspace and removed after analysis.")
    return st.session_state.get("batch_artifacts")


def result_downloads(artifacts: BatchArtifacts) -> None:
    left, middle, right, clear = st.columns(4)
    with left:
        st.download_button(
            "Download summary CSV",
            data=artifacts.summary_csv,
            file_name="summary.csv",
            mime="text/csv",
            use_container_width=True,
        )
    with middle:
        st.download_button(
            "Download extracted inserts",
            data=artifacts.inserts_fasta,
            file_name="extracted_inserts.fasta",
            mime="text/plain",
            use_container_width=True,
        )
    with right:
        st.download_button(
            "Download all reports",
            data=artifacts.results_zip,
            file_name="insert_visualizer_results.zip",
            mime="application/zip",
            use_container_width=True,
        )
    with clear:
        st.button("Clear session results", use_container_width=True, on_click=clear_session)


def clear_session() -> None:
    for key in ["batch_artifacts", "candidate_upload", "sample_map_upload", "backbone_upload", "sequencing_upload"]:
        st.session_state.pop(key, None)


def render_results(artifacts: BatchArtifacts, zoom: float, show_matches: bool) -> None:
    results = artifacts.results
    metadata = artifacts.metadata
    if not results:
        st.warning("No mapped samples were available for review.")
        return
    result_downloads(artifacts)
    if metadata.get("analysis_mode") == "extraction_only":
        render_extraction_results(artifacts, zoom, show_matches)
        return
    result_by_label = {
        f"{r.sample_id} · {r.expected_orf} · {r.status} · {r.percent_identity:.2f}%": r for r in results
    }
    selected_label = st.sidebar.radio("Sample list", list(result_by_label))
    result = result_by_label[selected_label]
    counts = pd.Series([item.status for item in results]).value_counts()
    backbone_counts = pd.Series([item.backbone_status for item in results]).value_counts()
    st.caption(
        f"Batch: {len(results)} mapped samples · Insert PASS {counts.get('PASS', 0)}, WARNING {counts.get('WARNING', 0)}, FAIL {counts.get('FAIL', 0)} · "
        f"Backbone PASS {backbone_counts.get('PASS', 0)}, WARNING {backbone_counts.get('WARNING', 0)}, FAIL {backbone_counts.get('FAIL', 0)}"
    )
    st.html(
        f"<div>{status_badge(result.status)} <strong>{escape(result.sample_id)}</strong> · {escape(result.colony_name)} · expected <code>{escape(result.expected_orf)}</code> · best <code>{escape(result.best_matching_orf)}</code></div>"
    )
    st.html(summary_cards(result))
    st.html(f"<div class='verdict'>{escape(result.verdict)}</div>")

    st.subheader("Expected ORF Coordinate Track")
    st.plotly_chart(coordinate_track(result, show_matches=show_matches, zoom=zoom), use_container_width=False)

    tab_backbone, tab_events, tab_stops, tab_ends, tab_alignment, tab_translation, tab_ranking, tab_overview = st.tabs(
        ["Backbone QC", "Insert Variants", "Stop Codons", "Sequence Ends", "Alignment Viewer", "ORF Translation", "Best ORF Ranking", "Overview"]
    )
    with tab_backbone:
        st.html(
            "<div class='summary-grid'>"
            f"<div class='metric-card'><div class='metric-value'>{escape(result.backbone_status)}</div><div class='metric-name'>Backbone status</div></div>"
            f"<div class='metric-card'><div class='metric-value'>{result.backbone_identity:.3f}%</div><div class='metric-name'>Identity</div></div>"
            f"<div class='metric-card'><div class='metric-value'>{result.backbone_coverage:.3f}%</div><div class='metric-name'>Coverage</div></div>"
            f"<div class='metric-card'><div class='metric-value'>{result.observed_backbone_length} / {result.expected_backbone_length} bp</div><div class='metric-name'>Observed / expected backbone</div></div>"
            "</div>"
        )
        st.caption(result.backbone_notes)
        st.plotly_chart(backbone_track(result, show_matches=show_matches, zoom=zoom), use_container_width=False)
        backbone_events = pd.DataFrame(result.backbone_events)
        st.dataframe(
            backbone_events if not backbone_events.empty else pd.DataFrame(columns=["event_type", "backbone_coordinate", "observed_coordinate", "expected_base", "observed_base"]),
            width="stretch",
            hide_index=True,
        )
    with tab_events:
        events = pd.DataFrame(result.events)
        st.dataframe(events if not events.empty else pd.DataFrame(columns=["event_type", "expected_coordinate", "observed_coordinate"]), width="stretch", hide_index=True)
    with tab_stops:
        st.dataframe(pd.DataFrame(result.stops), width="stretch", hide_index=True)
    with tab_ends:
        st.dataframe(pd.DataFrame([result.ends]), width="stretch", hide_index=True)
    with tab_alignment:
        event_options = ["All events"] + [f"{e['event_type']} at expected nt {e['expected_coordinate']}" for e in result.events[:500]]
        choice = st.selectbox("Focus", event_options)
        focus = None if choice == "All events" else int(choice.rsplit(" ", 1)[-1])
        st.html(alignment_blocks(result, focus_coordinate=focus))
    with tab_translation:
        st.write(
            {
                "observed_in_frame": not result.frameshift,
                "length_divisible_by_3": result.length_divisible_by_3,
                "internal_stop_codon_count": result.internal_stop_codon_count,
                "summary": result.amino_acid_translation_summary,
            }
        )
        st.dataframe(translation_view(result), width="stretch", hide_index=True)
    with tab_ranking:
        st.dataframe(pd.DataFrame(result.ranking), width="stretch", hide_index=True)
    with tab_overview:
        st.dataframe(overview_dataframe(results), width="stretch", hide_index=True)
    st.caption(f"Loaded {metadata['candidate_count']} candidate ORFs and {metadata['sample_count']} sequenced samples.")


def render_extraction_results(artifacts: BatchArtifacts, zoom: float, show_matches: bool) -> None:
    results = artifacts.results
    result_by_label = {
        f"{result.sample_id} · {result.observed_insert_length} bp · {result.status}": result for result in results
    }
    selected_label = st.sidebar.radio("Sample list", list(result_by_label))
    result = result_by_label[selected_label]
    counts = pd.Series([item.status for item in results]).value_counts()
    st.caption(
        f"Extraction-only mode · {len(results)} mapped samples · "
        f"PASS {counts.get('PASS', 0)}, WARNING {counts.get('WARNING', 0)}, FAIL {counts.get('FAIL', 0)}"
    )
    st.html(
        f"<div>{status_badge(result.status)} <strong>{escape(result.sample_id)}</strong> · "
        f"{escape(result.colony_name)} · extracted insert</div>"
    )
    backbone_variants = result.backbone_mismatches + result.backbone_insertions + result.backbone_deletions
    st.html(
        "<div class='summary-grid'>"
        f"<div class='metric-card'><div class='metric-value'>{result.observed_insert_length} bp</div><div class='metric-name'>Extracted insert</div></div>"
        f"<div class='metric-card'><div class='metric-value'>{escape(result.orientation)}</div><div class='metric-name'>Sample orientation</div></div>"
        f"<div class='metric-card'><div class='metric-value'>{escape(result.backbone_status)}</div><div class='metric-name'>Backbone status</div></div>"
        f"<div class='metric-card'><div class='metric-value'>{result.backbone_identity:.3f}%</div><div class='metric-name'>Backbone identity</div></div>"
        f"<div class='metric-card'><div class='metric-value'>{result.backbone_coverage:.3f}%</div><div class='metric-name'>Backbone coverage</div></div>"
        f"<div class='metric-card'><div class='metric-value'>{backbone_variants}</div><div class='metric-name'>Backbone variants</div></div>"
        "</div>"
    )
    st.html(f"<div class='verdict'>{escape(result.verdict)}</div>")

    insert_sequence = result.alignment.get("observed_insert", "")
    if insert_sequence:
        fasta = f">{result.sample_id} length={len(insert_sequence)}\n{insert_sequence}\n".encode()
        st.download_button(
            "Download selected insert",
            data=fasta,
            file_name=f"{result.sample_id}_insert.fasta",
            mime="text/plain",
        )

    tab_sequence, tab_backbone, tab_overview = st.tabs(["Extracted Sequence", "Backbone QC", "Overview"])
    with tab_sequence:
        st.code(insert_sequence or "No insert was extracted.", language="text", wrap_lines=True)
    with tab_backbone:
        st.caption(result.backbone_notes)
        st.plotly_chart(backbone_track(result, show_matches=show_matches, zoom=zoom), use_container_width=False)
        backbone_events = pd.DataFrame(result.backbone_events)
        st.dataframe(
            backbone_events if not backbone_events.empty else pd.DataFrame(
                columns=["event_type", "backbone_coordinate", "observed_coordinate", "expected_base", "observed_base"]
            ),
            width="stretch",
            hide_index=True,
        )
    with tab_overview:
        rows = [
            {
                "sample_id": item.sample_id,
                "colony_name": item.colony_name,
                "status": item.status,
                "insert_length": item.observed_insert_length,
                "sample_orientation": item.orientation,
                "backbone_status": item.backbone_status,
                "backbone_identity": item.backbone_identity,
                "backbone_coverage": item.backbone_coverage,
                "source_file": item.source_file,
            }
            for item in results
        ]
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)


def main() -> None:
    thresholds, zoom, show_matches = sidebar_controls()
    st.title("Insert Sequencing Visualizer")
    st.caption("Extract inserts from a Plasmidsaurus batch, review the complete backbone, and optionally compare candidate ORFs.")
    artifacts = upload_panel(thresholds)
    if artifacts is None:
        st.info("Upload sample mapping, backbone, and sequencing files. Candidate ORF FASTA is optional.")
        return
    st.divider()
    render_results(artifacts, zoom, show_matches)


if __name__ == "__main__":
    main()
