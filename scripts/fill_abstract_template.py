"""Fill the IROS 2026 Neuromuscular Robotics workshop poster-abstract
template with this project's Phase 5 content (CLAUDE.md section 4.2,
docs/abstract.md). Run once, locally (not part of the package/tests),
produces docs/abstract_IROS2026_neuromuscular_robotics.docx.

Author/affiliation/email filled per the user's explicit instruction
(single author, single affiliation, so the template's two-author/
two-affiliation fields are collapsed to one rather than left with a
dangling blank second entry).
"""

import copy
from pathlib import Path

import docx
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = ROOT / "abstract_template.docx"
FIGURE = ROOT / "results" / "phase4_activating_step_small.png"
# Submission naming requirement: POSTER_NAME_SURNAME.pdf. PDF itself is
# produced separately (LibreOffice conversion, see CLAUDE.md section 9
# Phase 5 note) -- this controls the .docx base name so both stay
# consistent.
OUT = ROOT / "docs" / "POSTER_AlanRoyceGabriel_Samuel.docx"

TITLE = ("Non-Smooth Equilibria in Hill-Type Muscle Models: A Quantified "
         "Consequence for Linearization-Based Control")

BACKGROUND = (
    "Controllers for muscle-actuated systems (prosthetics, exosuits, "
    "neuromuscular robots) routinely rely on equilibrium-linearization-"
    "based design (LQR, successive-linearization MPC), which implicitly "
    "assumes a single well-defined local linear model at each operating "
    "point. Whether Hill-type muscle models actually admit one has not "
    "been examined, distinct from Yeo et al.'s (2023) negative-stiffness "
    "eigenvalue-instability finding."
)

OBJECTIVES = (
    "Characterize the local dynamics of a Hill-type muscle-actuated joint "
    "at its operating equilibria, and determine whether any resulting "
    "model ambiguity is large enough to matter for a realistic control "
    "design."
)

METHODS = (
    "Single-joint elbow model (brachialis, Thelen (2003) Hill dynamics, "
    "Holzbaur et al. (2005), Murray et al. (1995), de Leva (1996) "
    "parameters), validated against published torque-angle data. "
    "Linearized at an isometric equilibrium and ZOH-discretized at a "
    "representative 10 ms control rate; the resulting equilibrium "
    "structure was formalized via bimodal piecewise-linear system theory "
    "and quantified against one-step MPC-style predictions."
)

RESULTS = (
    "Every equilibrium of this model sits exactly on two independent "
    "non-smooth switches: an activation-dynamics asymmetry and a "
    "force-velocity asymmetry, verified exactly (not approximately) "
    "against Camlibel et al.'s (2008) bimodal piecewise-linear continuity "
    "condition, with disjoint effects on the local Jacobian. In the "
    "small-signal regime, using the wrong local linear model produces "
    "joint-angle prediction error equal to 13-88% of the true signal "
    "across 50-500 ms horizons."
)

CONCLUSIONS = (
    "Hill-type muscle-actuated joints have an overlooked non-smoothness "
    "at every equilibrium that is practically significant, not a "
    "curiosity, for linearization-based control, directly relevant to "
    "neuromuscular and muscle-driven robotic systems. Future work: "
    "antagonist pair, closed-loop MPC test of branch-tracking against a "
    "naive fixed linearization."
)

FIGURE_CAPTION = (
    "True joint-angle trajectory vs. the four branch-linearizations' "
    "one-step predictions. The matched branch tracks truth almost "
    "exactly; the others diverge within 50-100 ms."
)

AUTHOR_NAME = "Alan Royce Gabriel Samuel"
AFFILIATION = "Indian Institute of Technology Madras, Chennai, India"
EMAIL = "alanroyce2010@gmail.com"


def set_paragraph_text(paragraph, text):
    """Replace a single-run paragraph's text, keeping its run's formatting."""
    if len(paragraph.runs) != 1:
        raise ValueError(f"expected exactly 1 run, got {len(paragraph.runs)}: {paragraph.text!r}")
    paragraph.runs[0].text = text


def insert_paragraph_after(paragraph, text="", style=None):
    """Standard python-docx recipe: insert a new paragraph immediately
    after the given one (python-docx has no public API for this)."""
    new_p = copy.deepcopy(paragraph._p)
    for child in list(new_p):
        new_p.remove(child)
    paragraph._p.addnext(new_p)
    new_paragraph = docx.text.paragraph.Paragraph(new_p, paragraph._parent)
    if style is not None:
        new_paragraph.style = style
    if text:
        new_paragraph.add_run(text)
    return new_paragraph


def fill_author_table(d):
    """Author(s)/Affiliation(s)/Presenting Author Email live in the
    template's single table cell, not the body paragraphs -- single
    author/affiliation here, so the template's second author and second
    affiliation lines are removed rather than left blank."""
    cell = d.tables[0].rows[0].cells[0]
    p_author, p_affiliation_label, p_aff1, p_aff2, p_email = cell.paragraphs

    # "Author(s): First Name Last Name¹, First Name Last Name²"
    # -> "Author(s): Alan Royce Gabriel Samuel" (no superscript needed,
    # only one affiliation).
    p_author.runs[1].text = AUTHOR_NAME

    # "¹ Department, Institution, City, Country" -> just the affiliation,
    # dropping the now-unneeded superscript-number prefix run.
    p_aff1.runs[0].text = ""
    p_aff1.runs[1].text = AFFILIATION

    # Second affiliation line is unused with one author -- remove outright
    # (same reasoning as the blank-spacer removal below: an empty line
    # still costs vertical space even with zeroed paragraph spacing).
    p_aff2._p.getparent().remove(p_aff2._p)

    p_email.add_run(EMAIL)


def main():
    d = docx.Document(str(TEMPLATE))
    fill_author_table(d)

    paras = d.paragraphs
    set_paragraph_text(paras[0], TITLE)
    set_paragraph_text(paras[5], BACKGROUND)
    set_paragraph_text(paras[8], OBJECTIVES)
    set_paragraph_text(paras[11], METHODS)
    set_paragraph_text(paras[14], RESULTS)
    set_paragraph_text(paras[17], CONCLUSIONS)
    # paras[18] was "Please state relevance to the workshop theme." --
    # folded into CONCLUSIONS above (explicitly addresses workshop
    # relevance), so clear this now-redundant instruction line.
    set_paragraph_text(paras[18], "")
    set_paragraph_text(paras[21], FIGURE_CAPTION)

    fig_paragraph = insert_paragraph_after(paras[21])
    run = fig_paragraph.add_run()
    if FIGURE.exists():
        run.add_picture(str(FIGURE), width=Inches(2.6))
    else:
        raise FileNotFoundError(f"figure not found: {FIGURE}")

    # Tighten spacing to reclaim room for the figure on page 1. Zeroing
    # space_before/after alone (tried first) barely moved the render --
    # each blank spacer paragraph still costs a full line at the base
    # font size regardless of its before/after spacing. Deleting them
    # outright is what actually reclaims the room (verified by
    # re-rendering after this change, not assumed).
    for p in d.paragraphs:
        p.paragraph_format.space_before = Pt(0)
        p.paragraph_format.space_after = Pt(4)
    for p in list(d.paragraphs):
        if not p.text.strip() and p._p is not fig_paragraph._p:
            p._p.getparent().remove(p._p)

    d.save(str(OUT))
    print(f"saved: {OUT}")


if __name__ == "__main__":
    main()
