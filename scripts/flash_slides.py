"""Flash-pitch slides (1-minute pitch, 2 slides, 16:9) for poster #5 at the
1st Workshop on Neuromuscular Robotics, IROS 2026. Output:
docs/poster/05_AlanRoyceGabriel_Samuel.pptx (filename format from the
organisers' email: XY_name_surname.pptx)."""

from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Emu, Inches, Pt

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "docs" / "poster" / "figures"
OUT = ROOT / "docs" / "poster" / "05_AlanRoyceGabriel_Samuel.pptx"

INK, MUTED, ACCENT = RGBColor(0x34, 0x49, 0x5E), RGBColor(0x6B, 0x7A, 0x89), RGBColor(0xEC, 0x00, 0x8B)
FONT = "Arial"


def text(slide, x, y, w, h, runs, size=18, color=INK, bold=False, align=PP_ALIGN.LEFT, spacing=1.1):
    """runs: list of paragraphs; each paragraph a list of (text, bold, color) tuples or a str."""
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    for i, para in enumerate(runs):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.line_spacing = spacing
        p.space_after = Pt(6)
        if isinstance(para, str):
            para = [(para, bold, color)]
        for t, b, c in para:
            r = p.add_run()
            r.text = t
            r.font.size = Pt(size)
            r.font.bold = b
            r.font.name = FONT
            r.font.color.rgb = c
    return tb


def badge(slide, x, y, w, h, label, size=40):
    shp = slide.shapes.add_shape(1, Inches(x), Inches(y), Inches(w), Inches(h))  # 1 = rectangle
    shp.fill.solid(); shp.fill.fore_color.rgb = ACCENT
    shp.line.fill.background()
    tf = shp.text_frame
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
    r = p.add_run(); r.text = label; r.font.size = Pt(size); r.font.bold = True
    r.font.name = FONT; r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)


def picture(slide, name, x, y, w=None, h=None):
    kw = {}
    if w: kw["width"] = Inches(w)
    if h: kw["height"] = Inches(h)
    return slide.shapes.add_picture(str(FIG / name), Inches(x), Inches(y), **kw)


prs = Presentation()
prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
blank = prs.slide_layouts[6]

# ---------------- slide 1: the claim ----------------
s = prs.slides.add_slide(blank)
badge(s, 11.9, 0.35, 1.1, 1.1, "5", size=44)
text(s, 0.5, 0.35, 11.2, 1.3, [[("Non-Smooth Equilibria in Hill-Type Muscle Models: ", True, INK),
                                ("a Quantified Consequence for Linearization-Based Control", True, ACCENT)]], size=28)
text(s, 0.5, 1.45, 11.2, 0.5, [[("Alan Royce Gabriel Samuel", True, INK),
                                ("   Indian Institute of Technology Madras   ·   poster 5", False, MUTED)]], size=16)
picture(s, "fig1_kinks.png", 0.4, 2.15, w=7.2)
text(s, 0.5, 5.05, 7.0, 0.5, ["Both switches sit exactly on every equilibrium: u* = a* and zero velocity."],
     size=15, color=MUTED)
text(s, 7.85, 2.15, 5.1, 4.9, [
    [("Equilibrium linearization assumes ", False, INK), ("one", True, INK), (" local linear model.", False, INK)],
    [("A Hill-type muscle has ", False, INK), ("four", True, ACCENT), (": two literature-cited switches (activation time constant, force-velocity) sit on every equilibrium. Verified exactly as a direct sum of two bimodal piecewise-linear switches.", False, INK)],
    [("Does it matter? ", True, INK), ("Wrong branch = ", False, INK), ("14 to 75%", True, ACCENT), (" of the predicted signal; both wrong, at least 23% at every equilibrium tested.", False, INK)],
    [("In closed loop: ", True, INK), ("8 to 9x the raising-step overshoot", True, ACCENT), (" for a blended linearization; smoothing only the activation switch still gives 5x the sign test's.", False, INK)],
], size=17, spacing=1.15)
text(s, 0.5, 6.85, 12.3, 0.4, ["Thelen 2003 · Holzbaur 2005 · Murray 1995 · Camlibel, Heemels & Schumacher 2008"],
     size=11, color=MUTED)

# ---------------- slide 2: the three results ----------------
s = prs.slides.add_slide(blank)
badge(s, 11.9, 0.35, 1.1, 1.1, "5", size=44)
text(s, 0.5, 0.4, 11.2, 0.8, [[("One message: ", True, INK), ("track which side of the switch you are on", True, ACCENT)]], size=26)
picture(s, "fig2_prediction.png", 0.3, 1.5, w=4.25)
picture(s, "fig4_closedloop.png", 4.55, 1.5, w=4.25)
picture(s, "fig5_pair_2panel.png", 9.05, 1.5, h=2.1)
text(s, 0.35, 3.75, 4.15, 1.6, [[("Open loop. ", True, ACCENT), ("A 0.01 excitation step: the matched branch tracks the truth; the other three miss by 14 to 75% of the swing within 50 to 500 ms.", False, INK)]], size=14)
text(s, 4.6, 3.75, 4.15, 1.6, [[("Closed loop. ", True, ACCENT), ("Only the linearization changes. Raising step: central differences overshoot 0.65 vs 0.07 deg; smoothing only activation, 0.42. Lowering: sign test 0.37, relinearizing at the measured state 0.02 deg.", False, INK)]], size=14)
text(s, 8.85, 3.75, 4.2, 1.6, [[("Antagonist pair. ", True, ACCENT), ("Structure survives (3 switches, 8 modes). Co-contraction shrinks the velocity kink from 2 to within 7% of 1 for both triceps moment arms (purple: x0.81).", False, INK)]], size=14)
text(s, 0.5, 5.35, 12.3, 1.0, [[("Take-home: ", True, INK), ("never central-difference through the kink; relinearize at the measured state and read the branch from it; compensate delay with a predictor.", False, INK)]], size=18)
text(s, 0.5, 6.3, 12.3, 0.5, [[("Come find me at poster 5.  ", True, ACCENT), ("Alan Royce Gabriel Samuel, IIT Madras", False, MUTED)]], size=16)

prs.save(OUT)
print("saved", OUT)
