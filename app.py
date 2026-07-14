from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from src.analysis import Thresholds, overview_dataframe, run_analysis
from src.visualization import alignment_blocks, backbone_track, coordinate_track, status_badge, summary_cards, translation_view


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT_DIR = ROOT / "input"
DEFAULT_CANDIDATE_ROOT = DEFAULT_INPUT_DIR
DEFAULT_SAMPLE_ROOT = DEFAULT_INPUT_DIR / "sequencing_results"
DEFAULT_SAMPLE_MAP = DEFAULT_INPUT_DIR / "sample_map.tsv"
DEFAULT_BACKBONE = DEFAULT_INPUT_DIR / "backbone.ape"
DEFAULT_RESULTS_DIR = ROOT / "results"


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
""",
)


def sidebar_controls() -> tuple[str, str, Path, Path, Path, str, int, Thresholds, float, bool, bool]:
    st.sidebar.title("Samples")
    with st.sidebar.expander("Input files", expanded=True):
        candidate_root = st.text_input("Candidate ORF directory", str(DEFAULT_CANDIDATE_ROOT))
        sample_root = st.text_input("Plasmidsaurus results directory", str(DEFAULT_SAMPLE_ROOT))
        sample_map = Path(st.text_input("Sample map", str(DEFAULT_SAMPLE_MAP)))
        backbone = Path(st.text_input("Backbone ApE/GenBank file", str(DEFAULT_BACKBONE)))
        replacement_marker = st.text_input("Replacement marker", "GGGCCCCCCCT")
        flank_length = st.number_input("Backbone flank length", min_value=12, max_value=120, value=40, step=1)
        results_dir = Path(st.text_input("Results directory", str(DEFAULT_RESULTS_DIR)))
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
    refresh = st.sidebar.button("Refresh analysis", type="primary")
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
    return candidate_root, sample_root, sample_map, results_dir, backbone, replacement_marker, int(flank_length), thresholds, zoom, show_matches, refresh


@st.cache_data(show_spinner=False)
def cached_analysis(
    candidate_root: str,
    sample_root: str,
    sample_map: str,
    results_dir: str,
    backbone: str,
    replacement_marker: str,
    flank_length: int,
    thresholds: Thresholds,
):
    return run_analysis(
        candidate_root,
        sample_root,
        sample_map,
        results_dir,
        thresholds,
        backbone_path=backbone,
        replacement_marker=replacement_marker,
        flank_length=flank_length,
    )


def main() -> None:
    candidate_root, sample_root, sample_map, results_dir, backbone, replacement_marker, flank_length, thresholds, zoom, show_matches, refresh = sidebar_controls()
    if refresh:
        cached_analysis.clear()
    try:
        with st.spinner("Loading sequences and aligning inserts and backbones..."):
            results, metadata = cached_analysis(candidate_root, sample_root, str(sample_map), str(results_dir), str(backbone), replacement_marker, flank_length, thresholds)
    except Exception as exc:
        st.error(f"Analysis could not run: {exc}")
        st.stop()
    if not results:
        st.warning("No mapped samples were available for review.")
        st.stop()
    result_by_label = {
        f"{r.sample_id} · {r.expected_orf} · {r.status} · {r.percent_identity:.2f}%": r for r in results
    }
    selected_label = st.sidebar.radio("Sample list", list(result_by_label), label_visibility="collapsed")
    result = result_by_label[selected_label]
    st.title("Insert Sequencing Visualizer")
    st.caption("Insert/ORF review plus full-length backbone QC outside the replacement marker. The default input directory is editable under ./input.")
    counts = pd.Series([item.status for item in results]).value_counts()
    backbone_counts = pd.Series([item.backbone_status for item in results]).value_counts()
    st.caption(
        f"Batch: {len(results)} mapped samples · Insert PASS {counts.get('PASS', 0)}, WARNING {counts.get('WARNING', 0)}, FAIL {counts.get('FAIL', 0)} · "
        f"Backbone PASS {backbone_counts.get('PASS', 0)}, WARNING {backbone_counts.get('WARNING', 0)}, FAIL {backbone_counts.get('FAIL', 0)}"
    )
    st.html(
        f"<div>{status_badge(result.status)} <strong>{result.sample_id}</strong> · {result.colony_name} · expected <code>{result.expected_orf}</code> · best <code>{result.best_matching_orf}</code></div>"
    )
    st.html(summary_cards(result))
    st.html(f"<div class='verdict'>{result.verdict}</div>")

    st.subheader("Expected ORF Coordinate Track")
    st.plotly_chart(coordinate_track(result, show_matches=show_matches, zoom=zoom), use_container_width=False)

    tab_backbone, tab_events, tab_stops, tab_ends, tab_alignment, tab_translation, tab_ranking, tab_overview = st.tabs(
        ["Backbone QC", "Insert Variants", "Stop Codons", "Sequence Ends", "Alignment Viewer", "ORF Translation", "Best ORF Ranking", "Overview"]
    )
    with tab_backbone:
        st.html(
            "<div class='summary-grid'>"
            f"<div class='metric-card'><div class='metric-value'>{result.backbone_status}</div><div class='metric-name'>Backbone status</div></div>"
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
        focus = None
        if choice != "All events":
            focus = int(choice.rsplit(" ", 1)[-1])
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

    st.caption(
        f"Loaded {metadata['candidate_count']} candidate ORFs and {metadata['sample_count']} sequenced samples. Unified report: {Path(metadata['results_dir']) / 'index.html'}."
    )


if __name__ == "__main__":
    main()
