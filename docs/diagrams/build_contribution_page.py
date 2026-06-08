#!/usr/bin/env python3
"""Generate and splice a 5th draw.io page showcasing the inqkit/inqview contribution.

Adds a clean, slide-level (16:9, 1920x1080) page to Misc/INQ-flow-chart.drawio:
  * GREY spine = pre-existing INQ base (Input -> Ground-state SCF -> Real-time TDDFT).
  * GREEN = my contribution, split into two lanes:
      - inqkit (in-run, C++):  wave-packet injection/orthonormalisation/validation,
        per-step callback into inq::real_time::propagate, in-run extraction & I/O.
      - inqview (post-run, Python): direct observables (densities, wavefunctions) and
        derived observables (loss function, momentum distribution, stopping power).
The existing 4 pages are left untouched; we insert one <diagram> before </mxfile>.

Colours: grey #EEEEEE/#9E9E9E = INQ base; green #D5E8D4/#82B366 = contribution.
Module names under each box are grounded in inq-stack/include/inqkit and
inq-stack/python/inqview.

Run: /local/data/public/skcb2/tddft/venv/bin/python3 Misc/build_contribution_page.py
"""
from pathlib import Path
import html

DRAWIO = Path("/local/data/public/skcb2/tddft/Misc/INQ-flow-chart.drawio")
PAGE_ID = "contrib_inqkit_inqview"
PAGE_NAME = "contribution"

# palette
GREY_FILL, GREY_STROKE, GREY_FONT = "#EEEEEE", "#9E9E9E", "#333333"
GRN_FILL, GRN_STROKE, GRN_FONT = "#D5E8D4", "#82B366", "#1B5E20"
GRNBG_FILL, GRNBG_STROKE = "#F1F8E9", "#82B366"
EDGE = "#666666"

cells = []


def esc(s):
    return html.escape(s, quote=True)


def xml_esc(s):
    """XML-attribute-escape a (possibly HTML) label so <font ...> tags and quotes
    survive inside value="...". drawio stores HTML labels this way."""
    return (s.replace("&", "&amp;").replace("<", "&lt;")
             .replace(">", "&gt;").replace('"', "&quot;"))


def label(title, sub=None, mod=None, *, tcolor="#000000", tsize=17, bold=True):
    """HTML label: bold title, optional bullet sub-line, optional small mono module line."""
    out = f'<font style="font-size:{tsize}px;" color="{tcolor}">'
    out += (f"<b>{esc(title)}</b>" if bold else esc(title)) + "</font>"
    if sub:
        out += f'<br><font style="font-size:13px;" color="{tcolor}">{esc(sub)}</font>'
    if mod:
        out += (f'<br><font style="font-size:10px;" color="#6B6B6B">'
                f'<i>{esc(mod)}</i></font>')
    return out


def box(cid, x, y, w, h, value, fill, stroke, font, *, rounded=True, dashed=False):
    r = "rounded=1;arcSize=8;" if rounded else "rounded=0;"
    d = "dashed=1;" if dashed else ""
    style = (f"{r}whiteSpace=wrap;html=1;{d}fillColor={fill};strokeColor={stroke};"
             f"fontColor={font};strokeWidth=2;align=center;verticalAlign=middle;"
             f"spacingLeft=6;spacingRight=6;")
    cells.append(
        f'        <mxCell id="{cid}" value="{xml_esc(value)}" style="{style}" vertex="1" parent="1">\n'
        f'          <mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry" />\n'
        f'        </mxCell>')


def container(cid, x, y, w, h, value, fill, stroke):
    style = (f"rounded=1;arcSize=4;whiteSpace=wrap;html=1;fillColor={fill};"
             f"strokeColor={stroke};strokeWidth=2;dashed=1;verticalAlign=top;align=left;"
             f"fontSize=18;fontStyle=1;fontColor={GRN_FONT};spacingLeft=14;spacingTop=10;")
    cells.append(
        f'        <mxCell id="{cid}" value="{xml_esc(value)}" style="{style}" vertex="1" parent="1">\n'
        f'          <mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry" />\n'
        f'        </mxCell>')


def edge(cid, src, dst, value="", *, color=EDGE, dashed=False, width=2):
    d = "dashed=1;" if dashed else ""
    style = (f"edgeStyle=orthogonalEdgeStyle;rounded=1;html=1;{d}strokeColor={color};"
             f"strokeWidth={width};fontColor={color};fontSize=12;endArrow=block;"
             f"endFill=1;")
    val = esc(value) if value else ""
    cells.append(
        f'        <mxCell id="{cid}" value="{val}" style="{style}" edge="1" parent="1" '
        f'source="{src}" target="{dst}">\n'
        f'          <mxGeometry relative="1" as="geometry" />\n'
        f'        </mxCell>')


# ---- title + legend ----
box("c_title", 40, 24, 1840, 56,
    label("rt-TDDFT wave-packet pipeline      —      INQ base  +  my inqkit / inqview contribution",
          tsize=24),
    "none", "none", "#000000", rounded=False)
# legend swatches
box("c_leg_g", 1430, 96, 30, 26, "", GREY_FILL, GREY_STROKE, GREY_FONT, rounded=False)
box("c_leg_gt", 1465, 96, 160, 26, label("pre-existing INQ", tsize=13, bold=False),
    "none", "none", "#000000", rounded=False)
box("c_leg_c", 1640, 96, 30, 26, "", GRN_FILL, GRN_STROKE, GRN_FONT, rounded=False)
box("c_leg_ct", 1675, 96, 200, 26, label("my contribution", tsize=13, bold=False),
    "none", "none", "#000000", rounded=False)

# ---- grey spine (INQ base) ----
box("g1", 60, 170, 300, 150,
    label("Input", "geometry + pseudopotentials",
          "inq::systems::ions · pseudopod"),
    GREY_FILL, GREY_STROKE, GREY_FONT)
box("g2", 430, 170, 300, 150,
    label("Ground-state SCF", "self-consistent Kohn–Sham",
          "inq::ground_state::calculator"),
    GREY_FILL, GREY_STROKE, GREY_FONT)
box("g3", 800, 170, 380, 170,
    label("Real-time TDDFT propagation", "TDKS time-stepping loop",
          "inq::real_time::propagate"),
    GREY_FILL, GREY_STROKE, GREY_FONT)
edge("e_g12", "g1", "g2")
edge("e_g23", "g2", "g3")

# ---- green lane A: inqkit (in-run, C++) ----
container("contA", 40, 380, 1150, 430,
          "My work · inqkit    —    in-run (C++)", GRNBG_FILL, GRNBG_STROKE)
box("a1", 820, 430, 340, 110,
    label("Per-step callback", "hooks into the propagation loop",
          "inqkit::real_time_session · step_context"),
    GRN_FILL, GRN_STROKE, GRN_FONT)
box("a2", 80, 560, 350, 220,
    label("Wave-packet", "• injection • orthonormalisation "
          "• per-run validation",
          "wavepacket · observables::orbital_overlap · injection_report"),
    GRN_FILL, GRN_STROKE, GRN_FONT)
box("a3", 470, 560, 330, 220,
    label("In-run extraction & I/O",
          "• fields → VTI • observables → CSV "
          "• manifest / run_summary",
          "io::vti_image_data_writer · observables_writer · manifest_writer"),
    GRN_FILL, GRN_STROKE, GRN_FONT)
edge("e_a2g3", "a2", "g3", "inject WP orbital", color=GRN_STROKE, width=2)
edge("e_g3a1", "g3", "a1", "per-step callback", color=GRN_STROKE, width=2)
edge("e_a1a3", "a1", "a3", "extract &amp; write", color=GRN_STROKE, width=2)

# ---- green lane B: inqview (post-run, Python) ----
container("contB", 1240, 380, 640, 430,
          "My work · inqview    —    post-run (Python)", GRNBG_FILL, GRNBG_STROKE)
box("b1", 1270, 470, 580, 140,
    label("Direct observables", "• densities • wavefunctions",
          "inqview: density · orbitals · fields · vti"),
    GRN_FILL, GRN_STROKE, GRN_FONT)
box("b2", 1270, 640, 580, 150,
    label("Derived observables",
          "• loss function • momentum distribution • stopping power",
          "inqview: lindhard · stopping · momentum · spectral_weight"),
    GRN_FILL, GRN_STROKE, GRN_FONT)
edge("e_a3b1", "a3", "b1", "VTI · CSV · manifest", color=EDGE, width=2)
edge("e_a3b2", "a3", "b2", "", color=EDGE, width=2)


# ---- assemble page ----
diagram = (
    f'  <diagram id="{PAGE_ID}" name="{PAGE_NAME}">\n'
    f'    <mxGraphModel dx="1422" dy="800" grid="1" gridSize="10" guides="1" '
    f'tooltips="1" connect="1" arrows="1" fold="1" page="1" pageScale="1" '
    f'pageWidth="1920" pageHeight="1080" math="1" shadow="0">\n'
    f'      <root>\n'
    f'        <mxCell id="0" />\n'
    f'        <mxCell id="1" parent="0" />\n'
    + "\n".join(cells) + "\n"
    f'      </root>\n'
    f'    </mxGraphModel>\n'
    f'  </diagram>\n'
)

text = DRAWIO.read_text()
assert text.rstrip().endswith("</mxfile>"), "unexpected drawio tail"
assert PAGE_ID not in text, "page already inserted; restore from .bak first"
new = text.rstrip()[: -len("</mxfile>")] + diagram + "</mxfile>\n"
DRAWIO.write_text(new)

# ---- validation: well-formed XML + page count ----
import xml.etree.ElementTree as ET
root = ET.fromstring(new)
pages = root.findall("diagram")
print(f"[check] XML well-formed; {len(pages)} pages (was 4, expect 5)")
print(f"[check] new page '{pages[-1].get('name')}' id={pages[-1].get('id')}")
ncells = len(cells)
print(f"[check] contribution page cells added: {ncells}")
assert len(pages) == 5 and pages[-1].get("id") == PAGE_ID
print("[ok] spliced contribution page into INQ-flow-chart.drawio")
