from __future__ import annotations

from dataclasses import asdict
from html import escape
import json
from pathlib import Path

import pandas as pd

from .analysis import AnalysisResult, translate


COLORS = {
    "match": "#2e9d63",
    "mismatch": "#d63b3b",
    "deletion": "#f29f05",
    "insertion": "#7b61ff",
    "ambiguous": "#0ea5e9",
    "stop": "#111827",
    "frameshift": "#b54708",
    "end": "#475467",
}


STATUS_COLORS = {"PASS": "#087443", "WARNING": "#b54708", "FAIL": "#b42318"}


def status_badge(status: str) -> str:
    color = STATUS_COLORS.get(status, "#344054")
    bg = {"PASS": "#dcfae6", "WARNING": "#fef0c7", "FAIL": "#fee4e2"}.get(status, "#f2f4f7")
    return f"<span class='badge' style='color:{color};background:{bg}'>{escape(status)}</span>"


def metric_card(name: str, value: str, accent: str = "#344054") -> str:
    return (
        "<div class='metric-card'>"
        f"<div class='metric-value' style='color:{accent}'>{escape(str(value))}</div>"
        f"<div class='metric-name'>{escape(name)}</div>"
        "</div>"
    )


def summary_cards(result: AnalysisResult) -> str:
    accent = STATUS_COLORS.get(result.status, "#344054")
    return "<div class='summary-grid'>" + "".join(
        [
            metric_card("Status", result.status, accent),
            metric_card("Identity", f"{result.percent_identity:.2f}%"),
            metric_card("Observed / Expected", f"{result.observed_insert_length} / {result.expected_orf_length} bp"),
            metric_card("Mismatches", str(result.mismatches)),
            metric_card("Indels", str(result.insertions + result.deletions)),
            metric_card("Internal Stops", str(result.internal_stop_codon_count)),
            metric_card("Frameshift", "YES" if result.frameshift else "NO", "#b54708" if result.frameshift else "#087443"),
            metric_card("Orientation", result.orientation),
            metric_card("Backbone", result.backbone_status, STATUS_COLORS.get(result.backbone_status, "#344054")),
            metric_card("Backbone Identity", f"{result.backbone_identity:.3f}%"),
            metric_card("Backbone Coverage", f"{result.backbone_coverage:.3f}%"),
            metric_card("Backbone Variants", str(result.backbone_mismatches + result.backbone_insertions + result.backbone_deletions)),
        ]
    ) + "</div>"


def coordinate_track(result: AnalysisResult, show_matches: bool = True, zoom: float = 1.0) -> go.Figure:
    import plotly.graph_objects as go

    width_hint = max(900, min(2400, int(result.expected_orf_length * zoom + 250)))
    fig = go.Figure()
    y = 0
    for run in result.runs:
        status = run["status"]
        if status == "match" and not show_matches:
            continue
        fig.add_trace(
            go.Scatter(
                x=[run["start"], run["end"]],
                y=[y, y],
                mode="lines",
                line={"color": COLORS.get(status, "#667085"), "width": 12},
                name=status,
                legendgroup=status,
                showlegend=not any(t.name == status for t in fig.data),
                hovertemplate=f"{status}<br>expected nt {run['start']}-{run['end']}<br>{run['count']} bp<extra></extra>",
            )
        )
    event_df = pd.DataFrame(result.events)
    for event_type, symbol, color in [
        ("mismatch", "circle", COLORS["mismatch"]),
        ("deletion", "diamond", COLORS["deletion"]),
        ("ambiguous base", "x", COLORS["ambiguous"]),
    ]:
        if not event_df.empty:
            sub = event_df[event_df["event_type"] == event_type]
            if not sub.empty:
                fig.add_trace(
                    go.Scatter(
                        x=sub["expected_coordinate"],
                        y=[0.16] * len(sub),
                        mode="markers",
                        marker={"symbol": symbol, "size": 10, "color": color, "line": {"width": 1, "color": "#ffffff"}},
                        name=event_type,
                        customdata=sub[["sample_id", "expected_orf", "observed_coordinate", "expected_base", "observed_base", "nearby_expected_context", "expected_codon", "observed_codon", "amino_acid_consequence"]].fillna("").to_numpy()
                        if {"expected_codon", "observed_codon", "amino_acid_consequence"}.issubset(sub.columns)
                        else sub.assign(expected_codon="", observed_codon="", amino_acid_consequence="")[["sample_id", "expected_orf", "observed_coordinate", "expected_base", "observed_base", "nearby_expected_context", "expected_codon", "observed_codon", "amino_acid_consequence"]].fillna("").to_numpy(),
                        hovertemplate=(
                            "event: " + event_type
                            + "<br>sample: %{customdata[0]}"
                            + "<br>expected ORF: %{customdata[1]}"
                            + "<br>expected nt: %{x}"
                            + "<br>observed nt: %{customdata[2]}"
                            + "<br>expected base: %{customdata[3]}"
                            + "<br>observed base: %{customdata[4]}"
                            + "<br>context: %{customdata[5]}"
                            + "<br>expected codon: %{customdata[6]}"
                            + "<br>observed codon: %{customdata[7]}"
                            + "<br>AA: %{customdata[8]}<extra></extra>"
                        ),
                    )
                )
    if not event_df.empty:
        insertions = event_df[event_df["event_type"] == "insertion"]
        for _, row in insertions.iterrows():
            x = row["expected_coordinate"]
            fig.add_trace(
                go.Scatter(
                    x=[x, x],
                    y=[-0.28, 0.28],
                    mode="lines",
                    line={"color": COLORS["insertion"], "width": 3},
                    name="insertion",
                    legendgroup="insertion",
                    showlegend=not any(t.name == "insertion" for t in fig.data),
                    hovertemplate=(
                        f"insertion<br>expected nt anchor: {x}<br>observed nt: {row.get('observed_coordinate','')}"
                        f"<br>observed base: {escape(str(row.get('observed_base','')))}<br>context: {escape(str(row.get('nearby_observed_context','')))}<extra></extra>"
                    ),
                )
            )
    stop_df = pd.DataFrame(result.stops)
    if not stop_df.empty:
        obs_stops = stop_df[stop_df["sequence"] == "observed"]
        if not obs_stops.empty:
            fig.add_trace(
                go.Scatter(
                    x=obs_stops["nt_start"].clip(upper=max(1, result.expected_orf_length)),
                    y=[0.42] * len(obs_stops),
                    mode="markers",
                    marker={"symbol": "triangle-down", "size": 13, "color": COLORS["stop"]},
                    name="stop codon",
                    customdata=obs_stops[["codon_index", "nt_start", "nt_end", "codon", "terminal"]].to_numpy(),
                    hovertemplate="observed stop codon<br>aa %{customdata[0]}<br>nt %{customdata[1]}-%{customdata[2]}<br>codon %{customdata[3]}<br>terminal: %{customdata[4]}<extra></extra>",
                )
            )
    fig.add_trace(
        go.Scatter(
            x=[1, result.expected_orf_length],
            y=[-0.48, -0.48],
            mode="markers+text",
            marker={"symbol": "line-ns-open", "size": 12, "color": COLORS["end"]},
            text=["expected start", "expected end"],
            textposition="bottom center",
            name="expected ends",
            hovertemplate="expected ORF %{text}: nt %{x}<extra></extra>",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[max(1, result.alignment.get("ref_start", 1)), min(result.expected_orf_length, result.alignment.get("ref_end", result.expected_orf_length))],
            y=[-0.66, -0.66],
            mode="markers+text",
            marker={"symbol": "line-ns-open", "size": 12, "color": "#101828"},
            text=["observed aligned start", "observed aligned end"],
            textposition="bottom center",
            name="observed ends",
            hovertemplate="observed aligned %{text}<br>expected coordinate %{x}<extra></extra>",
        )
    )
    fig.update_layout(
        height=360,
        width=width_hint,
        margin={"l": 30, "r": 30, "t": 20, "b": 60},
        plot_bgcolor="#ffffff",
        paper_bgcolor="#ffffff",
        xaxis={"title": "Expected ORF nucleotide coordinate", "range": [1, max(1, result.expected_orf_length)], "showgrid": True, "gridcolor": "#eef2f6"},
        yaxis={"visible": False, "range": [-0.9, 0.65]},
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.02, "xanchor": "left", "x": 0},
    )
    return fig


def backbone_track(result: AnalysisResult, show_matches: bool = True, zoom: float = 1.0) -> go.Figure:
    import plotly.graph_objects as go

    length = result.expected_backbone_length
    width_hint = max(900, min(2400, int(max(1, length) * zoom * 0.55 + 250)))
    fig = go.Figure()
    if not length or not result.backbone_runs:
        fig.add_annotation(
            text="Backbone could not be extracted and compared.",
            x=0.5,
            y=0.5,
            xref="paper",
            yref="paper",
            showarrow=False,
            font={"color": "#667085", "size": 14},
        )
    for run in result.backbone_runs:
        status = run["status"]
        if status == "match" and not show_matches:
            continue
        fig.add_trace(
            go.Scatter(
                x=[run["start"], run["end"]],
                y=[0, 0],
                mode="lines",
                line={"color": COLORS.get(status, "#667085"), "width": 13},
                name=status,
                legendgroup=f"backbone-{status}",
                showlegend=not any(t.name == status for t in fig.data),
                hovertemplate=f"{status}<br>backbone nt {run['start']}-{run['end']}<br>{run['count']} bp<extra></extra>",
            )
        )
    events = pd.DataFrame(result.backbone_events)
    for event_type, symbol, color in [
        ("mismatch", "circle", COLORS["mismatch"]),
        ("deletion", "diamond", COLORS["deletion"]),
        ("ambiguous base", "x", COLORS["ambiguous"]),
    ]:
        if events.empty:
            continue
        subset = events[events["event_type"] == event_type]
        if subset.empty:
            continue
        fig.add_trace(
            go.Scatter(
                x=subset["backbone_coordinate"],
                y=[0.22] * len(subset),
                mode="markers",
                marker={"symbol": symbol, "size": 10, "color": color, "line": {"width": 1, "color": "#ffffff"}},
                name=event_type,
                customdata=subset[["observed_coordinate", "expected_base", "observed_base"]].fillna("").to_numpy(),
                hovertemplate=(
                    f"{event_type}<br>backbone nt %{{x}}<br>observed nt %{{customdata[0]}}"
                    "<br>expected %{customdata[1]}<br>observed %{customdata[2]}<extra></extra>"
                ),
            )
        )
    if not events.empty:
        for _, event in events[events["event_type"] == "insertion"].iterrows():
            x = event["backbone_coordinate"]
            fig.add_trace(
                go.Scatter(
                    x=[x, x],
                    y=[-0.25, 0.25],
                    mode="lines",
                    line={"color": COLORS["insertion"], "width": 3},
                    name="insertion",
                    legendgroup="backbone-insertion",
                    showlegend=not any(t.name == "insertion" for t in fig.data),
                    hovertemplate=f"insertion after backbone nt {x}<br>observed base {escape(str(event.get('observed_base', '')))}<extra></extra>",
                )
            )
    fig.update_layout(
        height=280,
        width=width_hint,
        margin={"l": 30, "r": 30, "t": 20, "b": 55},
        plot_bgcolor="#ffffff",
        paper_bgcolor="#ffffff",
        xaxis={"title": "Expected backbone coordinate (replacement marker removed)", "range": [1, max(1, length)], "showgrid": True, "gridcolor": "#eef2f6"},
        yaxis={"visible": False, "range": [-0.45, 0.5]},
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.02, "xanchor": "left", "x": 0},
    )
    return fig


def alignment_blocks(result: AnalysisResult, focus_coordinate: int | None = None, flank: int = 80) -> str:
    ref = result.alignment.get("ref_aln", "")
    obs = result.alignment.get("obs_aln", "")
    if not ref or not obs:
        return "<p>No alignment available.</p>"
    positions: list[int | None] = []
    ref_pos = result.alignment.get("ref_start", 1) - 1
    for base in ref:
        if base != "-":
            ref_pos += 1
            positions.append(ref_pos)
        else:
            positions.append(None)
    start = 0
    end = len(ref)
    if focus_coordinate:
        indexes = [i for i, pos in enumerate(positions) if pos == focus_coordinate]
        if indexes:
            start = max(0, indexes[0] - flank)
            end = min(len(ref), indexes[0] + flank)
    chunks = []
    for i in range(start, end, 100):
        r = ref[i : i + 100]
        o = obs[i : i + 100]
        marker = "".join("|" if a == b and a != "-" else " " for a, b in zip(r, o))
        chunks.append(f"<pre class='alignment-block'>REF {escape(r)}\n    {marker}\nOBS {escape(o)}</pre>")
    return "".join(chunks)


def translation_view(result: AnalysisResult) -> pd.DataFrame:
    expected = translate(result.alignment.get("expected_sequence", ""))
    observed = translate(result.alignment.get("observed_insert", ""))
    rows = []
    for i, (exp, obs) in enumerate(zip(expected, observed), start=1):
        if exp != obs or obs == "*":
            rows.append({"aa_position": i, "expected_aa": exp, "observed_aa": obs, "note": "stop" if obs == "*" else "changed"})
    return pd.DataFrame(rows)


def dashboard_html(results: list[AnalysisResult], metadata: dict) -> str:
    records = []
    for result in results:
        row = asdict(result)
        row["category"] = {"PASS": "near-perfect", "WARNING": "frameshift", "FAIL": "failed"}.get(result.status, "review")
        row["observed_stops"] = [item for item in result.stops if item.get("sequence") == "observed"]
        row["expected_stops"] = [item for item in result.stops if item.get("sequence") == "expected"]
        row["identity_expected_space"] = round(100.0 * result.matches / max(1, result.expected_orf_length), 2)
        records.append(row)
    data_json = json.dumps(records)
    metadata_json = json.dumps(metadata)
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>Insert Sequencing Review</title>
<style>
:root {{ --bg:#f7f8fb; --panel:#fff; --ink:#20242c; --muted:#667085; --line:#d9dee8; --match:#2e9d63; --mismatch:#d63b3b; --deletion:#f29f05; --insertion:#7b61ff; --ambiguous:#0ea5e9; --stop:#111827; --failed:#b42318; --warn:#b54708; --pass:#087443; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; color:var(--ink); background:var(--bg); }}
header {{ padding:24px 28px 14px; border-bottom:1px solid var(--line); background:#fff; }}
h1 {{ margin:0 0 8px; font-size:24px; letter-spacing:0; }}
.subtitle {{ color:var(--muted); font-size:14px; line-height:1.45; max-width:1120px; }}
main {{ padding:18px 28px 36px; display:grid; grid-template-columns:360px minmax(0,1fr); gap:18px; }}
.panel {{ background:var(--panel); border:1px solid var(--line); border-radius:8px; }}
.sidebar {{ align-self:start; position:sticky; top:14px; max-height:calc(100vh - 28px); overflow:auto; }}
.controls {{ padding:14px; border-bottom:1px solid var(--line); display:grid; gap:12px; }}
label {{ font-size:12px; color:var(--muted); display:grid; gap:6px; }}
input[type="range"] {{ width:100%; }}
.legend {{ display:grid; grid-template-columns:1fr 1fr; gap:7px 10px; font-size:12px; color:var(--muted); }}
.swatch {{ display:inline-block; width:14px; height:9px; margin-right:6px; border-radius:2px; vertical-align:middle; }}
.sample-list {{ display:grid; }}
.sample-row {{ all:unset; cursor:pointer; display:grid; gap:6px; padding:12px 14px; border-bottom:1px solid var(--line); }}
.sample-row:hover,.sample-row.active {{ background:#f1f5f9; }}
.sample-title {{ font-size:13px; font-weight:650; display:flex; justify-content:space-between; gap:8px; }}
.badge {{ font-size:11px; border-radius:999px; padding:2px 8px; font-weight:800; white-space:nowrap; }}
.badge.failed {{ color:var(--failed); background:#fee4e2; }} .badge.frameshift {{ color:var(--warn); background:#fef0c7; }} .badge.near-perfect {{ color:var(--pass); background:#dcfae6; }}
.metric-line {{ font-size:12px; color:var(--muted); }}
.content {{ display:grid; gap:18px; min-width:0; }}
.summary-grid {{ display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:10px; padding:14px; }}
.metric {{ border:1px solid var(--line); border-radius:8px; padding:10px; background:#fbfcfe; }}
.metric .value {{ font-size:20px; font-weight:750; overflow-wrap:anywhere; }}
.metric .name {{ color:var(--muted); font-size:12px; margin-top:2px; }}
.section {{ padding:14px; }}
.section h2 {{ margin:0 0 10px; font-size:15px; }}
.verdict {{ border-left:4px solid var(--line); padding:10px 12px; background:#fbfcfe; font-size:14px; line-height:1.45; }}
.verdict.failed {{ border-left-color:var(--failed); }} .verdict.frameshift {{ border-left-color:var(--warn); }} .verdict.near-perfect {{ border-left-color:var(--pass); }}
.track-wrap {{ overflow-x:auto; border:1px solid var(--line); border-radius:8px; background:#fff; }}
svg {{ display:block; }}
.event-table {{ width:100%; border-collapse:collapse; font-size:12px; }}
th,td {{ text-align:left; padding:7px 8px; border-bottom:1px solid var(--line); vertical-align:top; }}
th {{ color:var(--muted); font-weight:650; background:#fbfcfe; position:sticky; top:0; }}
.code {{ font-family:ui-monospace,SFMono-Regular,Menlo,monospace; font-size:12px; word-break:break-all; }}
.tabs {{ display:flex; gap:6px; padding:10px 14px 0; }}
.tab {{ all:unset; cursor:pointer; font-size:12px; padding:6px 10px; border-radius:6px; color:var(--muted); }}
.tab.active {{ background:#e8eef7; color:var(--ink); font-weight:650; }}
.tooltip {{ position:fixed; pointer-events:none; background:#111827; color:#fff; padding:8px 10px; border-radius:6px; font-size:12px; line-height:1.35; max-width:320px; z-index:10; display:none; }}
.footer-note {{ color:var(--muted); font-size:12px; line-height:1.45; }}
@media (max-width:900px) {{ main {{ grid-template-columns:1fr; padding:14px; }} .sidebar {{ position:static; max-height:none; }} .summary-grid {{ grid-template-columns:repeat(2,minmax(0,1fr)); }} }}
</style>
</head>
<body>
<header>
  <h1>Insert Sequencing Review</h1>
  <div class="subtitle">Interactive comparison of each sequenced insert against its expected ORF plus full-length backbone QC outside the replacement marker. Green = identical, red = mismatch, orange = deletion, purple = insertion, blue = ambiguous base, black = stop codon.</div>
</header>
<main>
  <aside class="panel sidebar">
    <div class="controls">
      <label>Zoom <input id="zoom" type="range" min="0.7" max="2.8" value="1.15" step="0.05"></label>
      <label><input id="showMatches" type="checkbox" checked> Show matching segments</label>
      <div class="legend">
        <div><span class="swatch" style="background:var(--match)"></span>Match</div>
        <div><span class="swatch" style="background:var(--mismatch)"></span>Mismatch</div>
        <div><span class="swatch" style="background:var(--deletion)"></span>Deletion</div>
        <div><span class="swatch" style="background:var(--insertion)"></span>Insertion</div>
        <div><span class="swatch" style="background:var(--ambiguous)"></span>Ambiguous</div>
        <div><span class="swatch" style="background:var(--stop)"></span>Stop codon</div>
      </div>
    </div>
    <div id="sampleList" class="sample-list"></div>
  </aside>
  <section class="content">
    <div class="panel"><div id="summary" class="summary-grid"></div></div>
    <div class="panel section"><h2 id="sampleHeading"></h2><div id="verdict" class="verdict"></div></div>
    <div class="panel section"><h2>Predicted ORF Coordinate Track</h2><div id="track" class="track-wrap"></div></div>
    <div class="panel section"><h2>Full Backbone QC</h2><div id="backboneVerdict" class="verdict"></div><div id="backboneTrack" class="track-wrap" style="margin-top:10px"></div></div>
    <div class="panel">
      <div class="tabs">
        <button class="tab active" data-tab="events">Variant Events</button>
        <button class="tab" data-tab="stops">Stop Codons</button>
        <button class="tab" data-tab="ends">Sequence Ends</button>
      </div>
      <div id="tabBody" class="section"></div>
    </div>
    <div class="panel section footer-note" id="sourceInfo"></div>
  </section>
</main>
<div id="tooltip" class="tooltip"></div>
<script>
const DATA = {data_json};
const SOURCE = {metadata_json};
let selected = 0;
let activeTab = 'events';
const colors = {{ match:'#2e9d63', mismatch:'#d63b3b', deletion:'#f29f05', insertion:'#7b61ff', ambiguous:'#0ea5e9', stop:'#111827' }};
const list = document.getElementById('sampleList');
const tooltip = document.getElementById('tooltip');
function category(d) {{ if (d.status === 'PASS') return 'near-perfect'; if (d.status === 'WARNING') return 'frameshift'; return 'failed'; }}
function badgeText(d) {{ return d.status === 'PASS' ? 'PASS' : d.status === 'WARNING' ? 'WARNING' : 'FAIL'; }}
function esc(s) {{ return String(s ?? '').replace(/[&<>"']/g, m => ({{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}}[m])); }}
function renderList() {{
  list.innerHTML = DATA.map((d,i)=>`<button class="sample-row ${{i===selected?'active':''}}" data-i="${{i}}">
    <div class="sample-title"><span>${{esc(d.sample_id)}} · ${{esc(d.expected_orf)}}</span><span class="badge ${{category(d)}}">${{badgeText(d)}}</span></div>
    <div class="metric-line">${{esc(d.colony_name)}}</div>
    <div class="metric-line">Insert ${{d.observed_insert_length}} bp vs expected ${{d.expected_orf_length}} bp · ${{d.identity_expected_space}}% expected-space match</div>
  </button>`).join('');
  list.querySelectorAll('.sample-row').forEach(btn=>btn.onclick=()=>{{ selected=Number(btn.dataset.i); renderAll(); }});
}}
function metric(name,value) {{ return `<div class="metric"><div class="value">${{esc(value)}}</div><div class="name">${{esc(name)}}</div></div>`; }}
function renderSummary(d) {{
  document.getElementById('summary').innerHTML = [
    metric('Status', d.status), metric('Observed insert length', `${{d.observed_insert_length}} bp`), metric('Predicted ORF length', `${{d.expected_orf_length}} bp`), metric('Length delta', `${{d.length_difference>0?'+':''}}${{d.length_difference}} bp`),
    metric('Expected-space match', `${{d.identity_expected_space}}%`), metric('Mismatches', d.mismatches), metric('Indels', d.insertions + d.deletions), metric('Observed internal stops', d.internal_stop_codon_count),
    metric('Backbone status', d.backbone_status), metric('Backbone identity', `${{d.backbone_identity}}%`), metric('Backbone coverage', `${{d.backbone_coverage}}%`), metric('Backbone variants', d.backbone_mismatches + d.backbone_insertions + d.backbone_deletions)
  ].join('');
}}
function eventTitle(e) {{
  if (e.event_type === 'mismatch') return `Mismatch at predicted nt ${{e.expected_coordinate}}: expected ${{e.expected_base}}, observed ${{e.observed_base}}`;
  if (e.event_type === 'deletion') return `Deletion at predicted nt ${{e.expected_coordinate}}: expected ${{e.expected_base}}, observed gap`;
  if (e.event_type === 'insertion') return `Insertion after predicted nt ${{e.expected_coordinate}}: observed ${{e.observed_base}} at observed nt ${{e.observed_coordinate}}`;
  if (e.event_type === 'ambiguous base') return `Ambiguous base at predicted nt ${{e.expected_coordinate}}`;
  return '';
}}
function showTip(evt, text) {{ tooltip.innerHTML = esc(text); tooltip.style.display='block'; tooltip.style.left=(evt.clientX+12)+'px'; tooltip.style.top=(evt.clientY+12)+'px'; }}
function hideTip() {{ tooltip.style.display='none'; }}
function renderTrack(d) {{
  const zoom = Number(document.getElementById('zoom').value);
  const showMatches = document.getElementById('showMatches').checked;
  const padL=48, padR=22, y=60, h=22;
  const width = Math.max(900, Math.ceil(d.expected_orf_length * zoom + padL + padR));
  let svg = `<svg width="${{width}}" height="150" viewBox="0 0 ${{width}} 150" xmlns="http://www.w3.org/2000/svg">`;
  svg += `<text x="${{padL}}" y="24" font-size="12" fill="#667085">Predicted ORF coordinates: 1-${{d.expected_orf_length}} nt</text>`;
  svg += `<line x1="${{padL}}" y1="${{y+h+18}}" x2="${{width-padR}}" y2="${{y+h+18}}" stroke="#98a2b3" stroke-width="1"/>`;
  const tickStep = d.expected_orf_length > 800 ? 100 : 50;
  for (let p=0; p<=d.expected_orf_length; p+=tickStep) {{
    const x=padL+p*zoom;
    svg += `<line x1="${{x}}" x2="${{x}}" y1="${{y+h+14}}" y2="${{y+h+22}}" stroke="#98a2b3"/><text x="${{x}}" y="${{y+h+36}}" font-size="10" fill="#667085" text-anchor="middle">${{p}}</text>`;
  }}
  for (const r of d.runs) {{
    if (r.status==='match' && !showMatches) continue;
    const x=padL+(r.start-1)*zoom;
    const w=Math.max(1,(r.end-r.start+1)*zoom);
    svg += `<rect class="track-event" data-tip="${{esc(r.status)}}: predicted nt ${{r.start}}-${{r.end}} (${{r.count}} bp)" x="${{x.toFixed(2)}}" y="${{y}}" width="${{w.toFixed(2)}}" height="${{h}}" fill="${{colors[r.status] || '#667085'}}" opacity="${{r.status==='match'?0.78:0.95}}"/>`;
  }}
  for (const e of d.events.filter(x=>x.event_type==='insertion')) {{
    const x=padL+(e.expected_coordinate)*zoom;
    svg += `<line class="track-event" data-tip="${{esc(eventTitle(e))}}" x1="${{x}}" x2="${{x}}" y1="${{y-12}}" y2="${{y+h+8}}" stroke="${{colors.insertion}}" stroke-width="3"/>`;
  }}
  for (const s of d.observed_stops) {{
    const x=padL+(Math.min(s.nt_start,d.expected_orf_length)-1)*zoom;
    svg += `<path class="track-event" data-tip="${{s.terminal?'Terminal':'Internal'}} stop codon in observed insert: ${{s.codon}} at observed nt ${{s.nt_start}}-${{s.nt_end}} (aa ${{s.codon_index}})" d="M ${{x}} ${{y-17}} l 7 13 h -14 z" fill="${{colors.stop}}" opacity="${{s.terminal?0.65:1}}"/>`;
  }}
  svg += `</svg>`;
  document.getElementById('track').innerHTML=svg;
  document.querySelectorAll('.track-event').forEach(el=>{{ el.addEventListener('mousemove', e=>showTip(e, el.dataset.tip)); el.addEventListener('mouseleave', hideTip); }});
}}
function renderBackboneTrack(d) {{
  const zoom = Number(document.getElementById('zoom').value);
  const showMatches = document.getElementById('showMatches').checked;
  const length=d.expected_backbone_length || 0, padL=48, padR=22, y=48, h=20;
  const width=Math.max(900, Math.ceil(length*zoom*0.55+padL+padR));
  let svg=`<svg width="${{width}}" height="125" viewBox="0 0 ${{width}} 125" xmlns="http://www.w3.org/2000/svg">`;
  svg += `<text x="${{padL}}" y="22" font-size="12" fill="#667085">Expected backbone (replacement marker removed): 1-${{length}} nt</text>`;
  if (!length || !d.backbone_runs.length) svg += `<text x="${{padL}}" y="62" font-size="13" fill="#667085">Backbone could not be extracted and compared.</text>`;
  const scale=zoom*0.55;
  for (const r of d.backbone_runs) {{
    if (r.status==='match' && !showMatches) continue;
    const x=padL+(r.start-1)*scale, w=Math.max(1,(r.end-r.start+1)*scale);
    svg += `<rect class="track-event" data-tip="${{esc(r.status)}}: backbone nt ${{r.start}}-${{r.end}} (${{r.count}} bp)" x="${{x.toFixed(2)}}" y="${{y}}" width="${{w.toFixed(2)}}" height="${{h}}" fill="${{colors[r.status] || '#667085'}}"/>`;
  }}
  for (const e of d.backbone_events.filter(x=>x.event_type==='insertion')) {{
    const x=padL+e.backbone_coordinate*scale;
    svg += `<line class="track-event" data-tip="Insertion after backbone nt ${{e.backbone_coordinate}}: ${{esc(e.observed_base)}}" x1="${{x}}" x2="${{x}}" y1="${{y-10}}" y2="${{y+h+8}}" stroke="${{colors.insertion}}" stroke-width="3"/>`;
  }}
  svg += `</svg>`;
  document.getElementById('backboneTrack').innerHTML=svg;
  const v=document.getElementById('backboneVerdict'); v.className='verdict '+category({{status:d.backbone_status}}); v.textContent=`${{d.backbone_status}}: identity ${{d.backbone_identity}}%, coverage ${{d.backbone_coverage}}%, observed/expected ${{d.observed_backbone_length}}/${{d.expected_backbone_length}} bp. ${{d.backbone_notes}}`;
  document.querySelectorAll('.track-event').forEach(el=>{{ el.addEventListener('mousemove', e=>showTip(e, el.dataset.tip)); el.addEventListener('mouseleave', hideTip); }});
}}
function renderEvents(d) {{
  const rows = d.events.filter(e=>e.event_type!=='match').slice(0,500).map(e=>`<tr><td>${{esc(e.event_type)}}</td><td>${{e.expected_coordinate}}</td><td>${{e.observed_coordinate ?? ''}}</td><td class="code">${{esc(e.expected_base)}}</td><td class="code">${{esc(e.observed_base)}}</td><td>${{esc(e.amino_acid_consequence || eventTitle(e))}}</td></tr>`).join('');
  return `<table class="event-table"><thead><tr><th>Event</th><th>Predicted nt</th><th>Observed nt</th><th>Expected</th><th>Observed</th><th>Note</th></tr></thead><tbody>${{rows || '<tr><td colspan="6">No variant events.</td></tr>'}}</tbody></table>`;
}}
function renderStops(d) {{
  const rows = d.stops.map(s=>`<tr><td>${{esc(s.sequence)}}</td><td>${{s.codon_index}}</td><td>${{s.nt_start}}-${{s.nt_end}}</td><td class="code">${{esc(s.codon)}}</td><td>${{s.terminal?'Terminal stop':'Internal stop'}}</td></tr>`).join('');
  return `<table class="event-table"><thead><tr><th>Sequence</th><th>AA position</th><th>NT position</th><th>Codon</th><th>Type</th></tr></thead><tbody>${{rows || '<tr><td colspan="5">No in-frame stop codon found.</td></tr>'}}</tbody></table>`;
}}
function renderEnds(d) {{
  return `<table class="event-table"><tbody>
    <tr><th>Observed first 60 nt</th><td class="code">${{esc(d.ends.observed_first60)}}</td></tr>
    <tr><th>Predicted first 60 nt</th><td class="code">${{esc(d.ends.expected_first60)}}</td></tr>
    <tr><th>Observed last 60 nt</th><td class="code">${{esc(d.ends.observed_last60)}}</td></tr>
    <tr><th>Predicted last 60 nt</th><td class="code">${{esc(d.ends.expected_last60)}}</td></tr>
  </tbody></table>`;
}}
function renderTab(d) {{ document.getElementById('tabBody').innerHTML = activeTab==='events' ? renderEvents(d) : activeTab==='stops' ? renderStops(d) : renderEnds(d); }}
function renderAll() {{
  const d=DATA[selected];
  renderList(); renderSummary(d);
  document.getElementById('sampleHeading').textContent = `${{d.sample_id}} · ${{d.colony_name}} · expected ${{d.expected_orf}}`;
  const v=document.getElementById('verdict'); v.className='verdict '+category(d); v.textContent=d.verdict;
  renderTrack(d); renderBackboneTrack(d); renderTab(d);
  document.getElementById('sourceInfo').innerHTML = `Candidate root: <span class="code">${{esc(SOURCE.candidate_root)}}</span><br>Sample root: <span class="code">${{esc(SOURCE.sample_root)}}</span><br>Backbone: <span class="code">${{esc(SOURCE.backbone_path || 'not used')}}</span><br>Replacement marker: <span class="code">${{esc(SOURCE.replacement_marker || '')}}</span> · Flank length: ${{SOURCE.flank_length || ''}} bp · Extraction loaded: ${{SOURCE.backbone_loaded}}`;
}}
document.getElementById('zoom').addEventListener('input', ()=>{{ renderTrack(DATA[selected]); renderBackboneTrack(DATA[selected]); }});
document.getElementById('showMatches').addEventListener('change', ()=>{{ renderTrack(DATA[selected]); renderBackboneTrack(DATA[selected]); }});
document.querySelectorAll('.tab').forEach(t=>t.onclick=()=>{{ activeTab=t.dataset.tab; document.querySelectorAll('.tab').forEach(x=>x.classList.toggle('active', x===t)); renderTab(DATA[selected]); }});
renderAll();
</script>
</body>
</html>"""
