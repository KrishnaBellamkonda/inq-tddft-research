#!/usr/bin/env python3
"""Render the 'contribution' workflow page to a report-standard 1920x1080 PNG.

No draw.io CLI on this host, so this reproduces the 5th page of
INQ-flow-chart.drawio (Misc/build_contribution_page.py) directly with matplotlib,
in the report-1 visual idiom: serif (Computer Modern) typography, muted strokes,
restrained palette, full-canvas balanced layout. Grey = pre-existing INQ base;
green = my inqkit/inqview contribution.

Output: Misc/INQ-flow-chart-contribution.png  (+ copy in draft3_wp).
Run: /local/data/public/skcb2/tddft/venv/bin/python3 Misc/render_contribution_png.py
"""
from pathlib import Path
import shutil
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ---- typography: Calibri (Carlito = metric-compatible open clone of Calibri) ----
mpl.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": ["Calibri", "Carlito", "Liberation Sans", "DejaVu Sans"],
})

W, H = 1920, 1080
GREY_FILL, GREY_STROKE, GREY_FONT = "#E9E9E9", "#8A8A8A", "#2A2A2A"
GRN_FILL, GRN_STROKE, GRN_FONT = "#CDE5C4", "#4E8A3A", "#21501A"
GRNBG_FILL, GRNBG_STROKE = "#EFF6EB", "#6FA85B"
EDGE = "#555555"
MOD = "#6B6B6B"
OUT = Path("/local/data/public/skcb2/tddft/Misc/INQ-flow-chart-contribution.png")
OUT2 = Path("/local/data/public/skcb2/tddft/docs/presentations/assets/draft3_wp/INQ-flow-chart-contribution.png")

fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(H, 0)        # y-down to match draw.io coords
ax.axis("off")
fig.patch.set_facecolor("white")


# Font sizes are deliberately large: the 19.2in figure is shown at ~13.3in on a
# 16:9 slide (x0.69), so on-slide pt = these x0.69 (e.g. box title 27 -> ~19pt).
TITLE_SZ, SUB_SZ, MOD_SZ = 24, 17, 13


def container(x, y, w, h, title):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=12",
                 linewidth=1.6, edgecolor=GRNBG_STROKE, facecolor=GRNBG_FILL,
                 linestyle=(0, (7, 5)), zorder=1))
    ax.text(x + 22, y + 32, title, fontsize=19.5, fontweight="bold", style="italic",
            color=GRN_FONT, ha="left", va="center", zorder=3)


def box(x, y, w, h, title, sub, mod, fill, stroke, font, *, tsize=TITLE_SZ, rounded=True):
    bs = "round,pad=0,rounding_size=10" if rounded else "square,pad=0"
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=bs, linewidth=1.6,
                 edgecolor=stroke, facecolor=fill, zorder=4))
    cx = x + w / 2
    lines = [(title, tsize, "bold", font, "normal")]
    if sub:
        parts = [s.strip() for s in sub.split("•") if s.strip()]
        if "•" in sub:
            lines += [("•  " + s, SUB_SZ, "normal", font, "normal") for s in parts]
        else:
            lines.append((sub, SUB_SZ, "normal", font, "normal"))
    # module-name annotations (mod) intentionally omitted (user, 2026-06-03)
    lh = [sz * 2.05 for (_, sz, *_r) in lines]
    cur = y + h / 2 - sum(lh) / 2
    for (txt, sz, wt, col, st), dh in zip(lines, lh):
        ax.text(cx, cur + dh / 2, txt, fontsize=sz, fontweight=wt, color=col,
                style=st, ha="center", va="center", zorder=5)
        cur += dh


def edge(p0, p1, label="", color=EDGE, conn="arc3,rad=0", lw=2.2, double=False):
    style = "<|-|>" if double else "-|>"
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=18,
                 linewidth=lw, color=color, connectionstyle=conn,
                 shrinkA=3, shrinkB=4, zorder=6))
    if label:
        mx, my = (p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2
        ax.text(mx, my - 14, label, fontsize=15.5, color=color, style="italic",
                ha="center", va="center", zorder=7,
                bbox=dict(fc="white", ec="none", alpha=0.9, pad=1.5))


# ---- title + legend ----
ax.text(W / 2, 50, "rt-TDDFT wave-packet pipeline:  INQ base  +  inqkit / inqview contribution",
        fontsize=28, fontweight="bold", ha="center", va="center", color="#1a1a1a")
# legend sits on its own row below the title, right-aligned, clear of the spine
ax.add_patch(FancyBboxPatch((1300, 104), 32, 28, boxstyle="square,pad=0", linewidth=1.6,
             edgecolor=GREY_STROKE, facecolor=GREY_FILL, zorder=4))
ax.text(1340, 118, "pre-existing INQ", fontsize=16, ha="left", va="center")
ax.add_patch(FancyBboxPatch((1556, 104), 32, 28, boxstyle="square,pad=0", linewidth=1.6,
             edgecolor=GRN_STROKE, facecolor=GRN_FILL, zorder=4))
ax.text(1596, 118, "my contribution", fontsize=16, ha="left", va="center")

# ---- grey spine (INQ base) ----
box(80, 158, 380, 168, "Input", "geometry + pseudopotentials",
    "inq::systems::ions · pseudopod", GREY_FILL, GREY_STROKE, GREY_FONT, tsize=21)
box(540, 158, 380, 168, "Ground-state SCF", "self-consistent Kohn–Sham",
    "inq::ground_state::calculator", GREY_FILL, GREY_STROKE, GREY_FONT, tsize=21)
box(1000, 158, 440, 168, "Real-time TDDFT", "TDKS time-stepping loop",
    "inq::real_time::propagate", GREY_FILL, GREY_STROKE, GREY_FONT, tsize=21)
edge((460, 235), (540, 235))
edge((920, 235), (1000, 235))

# ---- green lane A: inqkit (in-run, C++) ----
container(80, 392, 1080, 596, "My work · inqkit    —    in-run (C++)")
box(330, 452, 560, 150, "Per-step callback", "hooks the propagation loop",
    "inqkit::real_time_session · step_context", GRN_FILL, GRN_STROKE, GRN_FONT)
box(120, 672, 490, 300, "Wave-packet",
    "• injection • orthonormalisation • per-run validation",
    "wavepacket · orbital_overlap · injection_report", GRN_FILL, GRN_STROKE, GRN_FONT)
box(640, 672, 490, 300, "In-run extraction & I/O",
    "• fields → VTI • observables → CSV • manifest / run_summary",
    "io::vti_image_data_writer · observables_writer", GRN_FILL, GRN_STROKE, GRN_FONT)
# g3 (RT propagation) <-> callback ; callback -> wavepacket (inject) & I/O (extract)
edge((1180, 320), (760, 452), "per-step callback", color=GRN_STROKE,
     conn="arc3,rad=0.16")
edge((470, 602), (360, 672), "inject WP", color=GRN_STROKE, conn="arc3,rad=0.12")
edge((760, 602), (880, 672), "extract & write", color=GRN_STROKE, conn="arc3,rad=-0.12")

# ---- green lane B: inqview (post-run, Python) ----
container(1240, 392, 600, 596, "My work · inqview    —    post-run (Python)")
box(1280, 470, 520, 230, "Direct observables", "• densities • wavefunctions",
    "inqview: density · orbitals · fields · vti", GRN_FILL, GRN_STROKE, GRN_FONT)
box(1280, 730, 520, 250, "Derived observables",
    "• loss function • momentum distribution • stopping power",
    "inqview: lindhard · stopping · momentum · spectral_weight",
    GRN_FILL, GRN_STROKE, GRN_FONT)
edge((1130, 770), (1280, 585), "VTI · CSV · manifest", color=EDGE, conn="arc3,rad=-0.1")
edge((1130, 850), (1280, 855), color=EDGE, conn="arc3,rad=0.05")

# ---- save + validate ----
fig.savefig(OUT, dpi=100, facecolor="white")
shutil.copy(OUT, OUT2)
plt.close(fig)
from PIL import Image
im = Image.open(OUT)
print(f"[check] saved {OUT}  size={im.size} (expect (1920, 1080))")
assert im.size == (W, H)
print(f"[check] copy -> {OUT2}")
