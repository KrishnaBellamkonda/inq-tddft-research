#!/usr/bin/env python3
"""Workflow-diagram template — .drawio source + matched matplotlib PNG.

The project has NO headless draw.io renderer, so a workflow diagram is two
artifacts sharing ONE node/edge definition (see scientific-figures/SKILL.md §7):

  1. a .drawio page (editable source, grounded module names, MathJax equations)
  2. a matched 1920x1080 matplotlib PNG in the report idiom (the deliverable)

This is a trimmed, copy-me starting point distilled from
docs/diagrams/{build_contribution_page.py, render_contribution_png.py}. Define
your NODES/EDGES once, then both emit_drawio() and render_png() consume them.

Run with venv python:
  /local/data/public/skcb2/tddft/venv/bin/python3 workflow_render_template.py
"""
from pathlib import Path
import html
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ---- report idiom: Calibri-metric sans, muted restrained palette ----
mpl.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": ["Calibri", "Carlito", "Liberation Sans", "DejaVu Sans"],
})
W, H = 1920, 1080
FILL, STROKE, FONT = "#E9E9E9", "#8A8A8A", "#2A2A2A"      # neutral box
ACC_FILL, ACC_STROKE, ACC_FONT = "#CDE5C4", "#4E8A3A", "#21501A"  # accent box
CONT_FILL, CONT_STROKE = "#EFF6EB", "#6FA85B"             # swimlane container
EDGE = "#555555"

# ---- ONE definition: nodes (id -> dict) and edges (src -> dst) ----
# x,y,w,h in a 0..1920 / 0..1080 canvas (y-down, draw.io convention).
NODES = {
    "in":  dict(x=120, y=460, w=360, h=150, title="Region average",
                sub="mean density over the slab region", accent=False),
    "n":   dict(x=560, y=460, w=360, h=150, title="Choose N (even)",
                sub="round to nearest even electron count", accent=True),
    "pois":dict(x=1000, y=460, w=380, h=150, title="Poisson solve",
                sub=r"$v_{bg}=-\,\mathrm{poisson}(n_+)$", accent=True),
    "add": dict(x=1460, y=460, w=340, h=150, title="Add to KS potential",
                sub="every Hamiltonian rebuild", accent=True),
}
EDGES = [("in", "n", ""), ("n", "pois", ""), ("pois", "add", "")]


# ---------- matplotlib render ----------
def _box(ax, nd):
    fill, stroke, font = (ACC_FILL, ACC_STROKE, ACC_FONT) if nd["accent"] else (FILL, STROKE, FONT)
    ax.add_patch(FancyBboxPatch((nd["x"], nd["y"]), nd["w"], nd["h"],
                 boxstyle="round,pad=0,rounding_size=10", linewidth=1.6,
                 edgecolor=stroke, facecolor=fill, zorder=4))
    cx = nd["x"] + nd["w"] / 2
    ax.text(cx, nd["y"] + nd["h"] / 2 - 16, nd["title"], fontsize=24,
            fontweight="bold", color=font, ha="center", va="center", zorder=5)
    if nd.get("sub"):
        ax.text(cx, nd["y"] + nd["h"] / 2 + 22, nd["sub"], fontsize=16,
                color=font, ha="center", va="center", zorder=5)


def _edge(ax, a, b, label=""):
    p0 = (a["x"] + a["w"], a["y"] + a["h"] / 2)
    p1 = (b["x"], b["y"] + b["h"] / 2)
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=18,
                 linewidth=2.2, color=EDGE, connectionstyle="arc3,rad=0",
                 shrinkA=3, shrinkB=4, zorder=6))
    if label:
        mx, my = (p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2
        ax.text(mx, my - 14, label, fontsize=15, color=EDGE, style="italic",
                ha="center", va="center", zorder=7,
                bbox=dict(fc="white", ec="none", alpha=0.9, pad=1.5))


def render_png(out_png, title="Workflow"):
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0)
    ax.axis("off"); fig.patch.set_facecolor("white")
    ax.text(W / 2, 60, title, fontsize=28, fontweight="bold", ha="center",
            va="center", color="#1a1a1a")  # presentation: title PRESENT
    for nd in NODES.values():
        _box(ax, nd)
    for s, d, lab in EDGES:
        _edge(ax, NODES[s], NODES[d], lab)
    fig.savefig(out_png, dpi=100, facecolor="white"); plt.close(fig)
    print(f"[png] {out_png}")


# ---------- .drawio emit ----------
def _esc(s):
    return html.escape(s, quote=True)


def emit_drawio(out_drawio, page_name="workflow"):
    cells = ['        <mxCell id="0" />', '        <mxCell id="1" parent="0" />']
    for cid, nd in NODES.items():
        fill, stroke = (ACC_FILL, ACC_STROKE) if nd["accent"] else (FILL, STROKE)
        val = nd["title"] + ("<br>" + _esc(nd["sub"]) if nd.get("sub") else "")
        style = f"rounded=1;whiteSpace=wrap;html=1;fillColor={fill};strokeColor={stroke};"
        cells.append(
            f'        <mxCell id="{cid}" value="{val}" style="{style}" vertex="1" parent="1">\n'
            f'          <mxGeometry x="{nd["x"]}" y="{nd["y"]}" width="{nd["w"]}" height="{nd["h"]}" as="geometry"/>\n'
            f'        </mxCell>')
    for s, d, lab in EDGES:
        cells.append(
            f'        <mxCell id="e_{s}_{d}" value="{_esc(lab)}" style="edgeStyle=orthogonalEdgeStyle;rounded=0;html=1;" '
            f'edge="1" parent="1" source="{s}" target="{d}"><mxGeometry relative="1" as="geometry"/></mxCell>')
    body = "\n".join(cells)
    xml = (f'<mxfile host="app">\n  <diagram id="wf" name="{page_name}">\n'
           f'    <mxGraphModel dx="1422" dy="800" grid="1" gridSize="10" guides="1" '
           f'page="1" pageWidth="{W}" pageHeight="{H}" math="1" background="none">\n'
           f'      <root>\n{body}\n      </root>\n    </mxGraphModel>\n  </diagram>\n</mxfile>\n')
    Path(out_drawio).write_text(xml)
    import xml.etree.ElementTree as ET
    ET.fromstring(xml)  # validate well-formed
    print(f"[drawio] {out_drawio}")


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    emit_drawio(here / "_example_workflow.drawio", "example")
    render_png(here / "_example_workflow.png", "Localised-background injection")
