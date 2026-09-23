"""Flash-pitch slides built from the project video (1-minute pitch, poster 5,
1st Workshop on Neuromuscular Robotics, IROS 2026).

Source video: results/submission_video/Hill_ICRA27.mp4 (144 s, 1080p, silent).
Its seven chapters are cut into four clips whose lengths match the pitch script
in docs/poster/flash_pitch.md (34 + 9 + 12 + 6 s), each embedded in one slide
and set to start automatically when the slide appears. Slides advance on click
(three clicks in 60 s; cues in the speaker notes). A timed auto-advance was
tried and dropped: PowerPoint honours the times on slides without auto-started
media, but with an auto-playing clip it advanced slides early and
non-deterministically (measured with PowerPoint's own slideshow position).

Output: docs/poster/05_AlanRoyceGabriel_Samuel.pptx (organisers' XY_name_surname
rule). Clips and poster frames: results/submission_video/flash_clips/.

Run with a python that has python-pptx and lxml; FFmpeg must be on PATH.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from lxml import etree
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "results" / "submission_video" / "Hill_ICRA27.mp4"
CLIPS = ROOT / "results" / "submission_video" / "flash_clips"
OUT = Path(os.environ.get("FLASH_OUT", ROOT / "docs" / "poster" / "05_AlanRoyceGabriel_Samuel.pptx"))
# FLASH_TIMING=simple keeps python-pptx's own media timing with a zero delay
# (diagnostic variant; PowerPoint reads it as play-on-entry=false, i.e. no
# autoplay). Default writes PowerPoint's own "start automatically" timing tree,
# which PowerPoint reads as play-on-entry=true.
TIMING_MODE = os.environ.get("FLASH_TIMING", "mediacall")

# Palette of the video itself (cream canvas, navy ink, red accent), so the
# embedded frames do not sit as a cream rectangle on a white slide.
BG = RGBColor(0xFD, 0xF0, 0xD5)
INK = RGBColor(0x00, 0x30, 0x49)
MUTED = RGBColor(0x40, 0x59, 0x67)
ACCENT = RGBColor(0xC1, 0x12, 0x1F)
FONT = "Arial"

# Chapter boundaries of Hill_ICRA27.mp4 (seconds), verified frame by frame.
CH = {"intro": (0, 12), "switches": (12, 34), "mismatch": (34, 56),
      "raise": (56, 80), "lower": (80, 100), "track": (100, 124),
      "pair": (124, 142)}

REBUILD_CLIPS = "--recut" in sys.argv

NS = {"p": "http://schemas.openxmlformats.org/presentationml/2006/main",
      "a": "http://schemas.openxmlformats.org/drawingml/2006/main"}


def run(cmd):
    subprocess.run(cmd, check=True)


def cut(name, start, end, speed=1.0, pad_to=None):
    """Cut [start, end) from the source, play it `speed` times faster, and
    hold the last frame until `pad_to` seconds if given. Returns the path."""
    dest = CLIPS / f"{name}.mp4"
    if dest.exists() and not REBUILD_CLIPS:
        return dest
    dur = (end - start) / speed
    vf = f"setpts=PTS/{speed},fps=30"
    if pad_to and pad_to > dur:
        vf += f",tpad=stop_mode=clone:stop_duration={pad_to - dur:.3f}"
    run(["ffmpeg", "-v", "error", "-y", "-ss", str(start), "-t", str(end - start),
         "-i", str(SRC), "-vf", vf, "-an", "-c:v", "libx264", "-crf", "22",
         "-preset", "slow", "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(dest)])
    return dest


def cut_pair(name, left, right, speed, pad_to):
    """Two chapters side by side (each 1920x1080, output 3840x1080)."""
    dest = CLIPS / f"{name}.mp4"
    if dest.exists() and not REBUILD_CLIPS:
        return dest
    (l0, l1), (r0, r1) = left, right
    ld, rd = (l1 - l0) / speed, (r1 - r0) / speed
    fl = f"[0:v]setpts=PTS/{speed},fps=30,tpad=stop_mode=clone:stop_duration={max(0, pad_to - ld):.3f}[l]"
    fr = f"[1:v]setpts=PTS/{speed},fps=30,tpad=stop_mode=clone:stop_duration={max(0, pad_to - rd):.3f}[r]"
    run(["ffmpeg", "-v", "error", "-y",
         "-ss", str(l0), "-t", str(l1 - l0), "-i", str(SRC),
         "-ss", str(r0), "-t", str(r1 - r0), "-i", str(SRC),
         "-filter_complex", f"{fl};{fr};[l][r]hstack=inputs=2,trim=duration={pad_to}",
         "-an", "-c:v", "libx264", "-crf", "22", "-preset", "slow", "-pix_fmt", "yuv420p",
         "-movflags", "+faststart", str(dest)])
    return dest


def cut_seq(name, parts):
    """Chapters played one after another, full frame, each at its own speed."""
    dest = CLIPS / f"{name}.mp4"
    if dest.exists() and not REBUILD_CLIPS:
        return dest
    pieces = []
    for k, ((t0, t1), speed) in enumerate(parts):
        piece = CLIPS / f"_{name}_{k}.mp4"
        run(["ffmpeg", "-v", "error", "-y", "-ss", str(t0), "-t", str(t1 - t0), "-i", str(SRC),
             "-vf", f"setpts=PTS/{speed},fps=30", "-an", "-c:v", "libx264", "-crf", "22",
             "-preset", "slow", "-pix_fmt", "yuv420p", str(piece)])
        pieces.append(piece)
    lst = CLIPS / f"_{name}.txt"
    lst.write_text("".join(f"file '{p}'\n" for p in pieces))
    run(["ffmpeg", "-v", "error", "-y", "-f", "concat", "-safe", "0", "-i", str(lst), "-c:v", "libx264",
         "-crf", "22", "-preset", "slow", "-pix_fmt", "yuv420p", "-movflags", "+faststart", "-an", str(dest)])
    for f in pieces + [lst]:
        f.unlink()
    return dest


def poster_frame(clip, t=0.2):
    png = clip.with_suffix(".png")
    run(["ffmpeg", "-v", "error", "-y", "-ss", str(t), "-i", str(clip), "-frames:v", "1", str(png)])
    return png


def duration(clip):
    out = subprocess.run(["ffprobe", "-v", "error", "-show_entries", "format=duration",
                          "-of", "csv=p=0", str(clip)], capture_output=True, text=True, check=True)
    return float(out.stdout.strip())


def text(slide, x, y, w, h, runs, size=18, color=INK, bold=False, align=PP_ALIGN.LEFT, spacing=1.1):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = 0
    for i, para in enumerate(runs):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.line_spacing = spacing
        p.space_after = Pt(4)
        if isinstance(para, str):
            para = [para]
        para = [(t, bold, color) if isinstance(t, str) else t for t in para]
        for t, b, c in para:
            r = p.add_run()
            r.text = t
            r.font.size = Pt(size)
            r.font.bold = b
            r.font.name = FONT
            r.font.color.rgb = c
    return tb


def badge(slide, x, y, w, h, label, size=40):
    shp = slide.shapes.add_shape(1, Inches(x), Inches(y), Inches(w), Inches(h))
    shp.fill.solid()
    shp.fill.fore_color.rgb = ACCENT
    shp.line.fill.background()
    shp.shadow.inherit = False
    p = shp.text_frame.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    r = p.add_run()
    r.text = label
    r.font.size = Pt(size)
    r.font.bold = True
    r.font.name = FONT
    r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)


TIMING = """<p:timing xmlns:p="{p}" xmlns:a="{a}">
  <p:tnLst>
    <p:par>
      <p:cTn id="1" dur="indefinite" restart="never" nodeType="tmRoot">
        <p:childTnLst>
          <p:seq concurrent="1" nextAc="seek">
            <p:cTn id="2" dur="indefinite" nodeType="mainSeq">
              <p:childTnLst>
                <p:par>
                  <p:cTn id="3" fill="hold">
                    <p:stCondLst>
                      <p:cond delay="indefinite"/>
                      <p:cond evt="onBegin" delay="0"><p:tn val="2"/></p:cond>
                    </p:stCondLst>
                    <p:childTnLst>
                      <p:par>
                        <p:cTn id="4" fill="hold">
                          <p:stCondLst><p:cond delay="0"/></p:stCondLst>
                          <p:childTnLst>
                            <p:par>
                              <p:cTn id="5" presetID="1" presetClass="mediacall" presetSubtype="0" fill="hold" nodeType="withEffect">
                                <p:stCondLst><p:cond delay="0"/></p:stCondLst>
                                <p:childTnLst>
                                  <p:cmd type="call" cmd="playFrom(0.0)">
                                    <p:cBhvr>
                                      <p:cTn id="6" dur="{ms}" fill="hold"/>
                                      <p:tgtEl><p:spTgt spid="{spid}"/></p:tgtEl>
                                    </p:cBhvr>
                                  </p:cmd>
                                </p:childTnLst>
                              </p:cTn>
                            </p:par>
                          </p:childTnLst>
                        </p:cTn>
                      </p:par>
                    </p:childTnLst>
                  </p:cTn>
                </p:par>
              </p:childTnLst>
            </p:cTn>
            <p:prevCondLst>
              <p:cond evt="onPrev" delay="0"><p:tgtEl><p:sldTgt/></p:tgtEl></p:cond>
            </p:prevCondLst>
            <p:nextCondLst>
              <p:cond evt="onNext" delay="0"><p:tgtEl><p:sldTgt/></p:tgtEl></p:cond>
            </p:nextCondLst>
          </p:seq>
          <p:video>
            <p:cMediaNode vol="80000">
              <p:cTn id="7" fill="hold" display="0">
                <p:stCondLst><p:cond delay="indefinite"/></p:stCondLst>
                <p:endCondLst>
                  <p:cond evt="onNext" delay="0"><p:tgtEl><p:sldTgt/></p:tgtEl></p:cond>
                </p:endCondLst>
              </p:cTn>
              <p:tgtEl><p:spTgt spid="{spid}"/></p:tgtEl>
            </p:cMediaNode>
          </p:video>
        </p:childTnLst>
      </p:cTn>
    </p:par>
  </p:tnLst>
</p:timing>"""


def movie(slide, clip, x, y, w, h):
    """Embed the clip and make it start automatically on slide entry."""
    png = poster_frame(clip)
    shp = slide.shapes.add_movie(str(clip), Inches(x), Inches(y), Inches(w), Inches(h),
                                 poster_frame_image=str(png), mime_type="video/mp4")
    spid = shp.shape_id
    ms = int(round(duration(clip) * 1000))
    sld = slide._element
    if TIMING_MODE == "simple":
        for cond in sld.xpath(".//p:video/p:cMediaNode/p:cTn/p:stCondLst/p:cond"):
            cond.set("delay", "0")
        timing = sld.find("p:timing", NS)
        sld.remove(timing)
    else:
        for old in sld.findall("p:timing", NS):
            sld.remove(old)
        timing = etree.fromstring(TIMING.format(p=NS["p"], a=NS["a"], ms=ms, spid=spid))
    # p:transition precedes p:timing; both come after p:clrMapOvr.
    trans = etree.SubElement(sld, "{%s}transition" % NS["p"])
    trans.set("spd", "fast")
    trans.set("advClick", "1")
    sld.append(timing)
    return shp


def main():
    CLIPS.mkdir(parents=True, exist_ok=True)
    c1 = cut("s1_question", CH["intro"][0], CH["switches"][1])                   # 34 s
    c2 = cut("s2_open_loop", *CH["mismatch"], speed=2.5, pad_to=9.0)            # 8.8 -> 9 s
    c3 = cut_seq("s3_closed_loop_seq", [(CH["raise"], 3.0), (CH["lower"], 4.0)])       # 8 + 5 s
    c4 = cut("s4_pair", *CH["pair"], speed=3.0, pad_to=6.0)                      # 6 s

    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
    blank = prs.slide_layouts[6]

    def new_slide(notes):
        s = prs.slides.add_slide(blank)
        s.background.fill.solid()
        s.background.fill.fore_color.rgb = BG
        badge(s, 11.9, 0.3, 1.1, 1.1, "5", size=44)
        s.notes_slide.notes_text_frame.text = notes
        return s

    # Slide 1 (0-35 s): the question and the two switches.
    s = new_slide(
        "0-35 s. Hi, I'm Alan from IIT Madras, poster five. If you control anything driven "
        "by muscle, a prosthesis, an exosuit, an FES stimulator, your controller starts from "
        "a model of how the limb answers a small push. Real muscle doesn't answer the same "
        "way in both directions: it switches on three times faster than it switches off, "
        "and it fights being stretched twice as hard as shortening. A limb at rest sits "
        "exactly where both of these flip. So there isn't one local model. There are four. "
        "[Clip: the elbow, then the two switches at the resting posture. Click at 35 s.]")
    text(s, 0.5, 0.3, 11.2, 0.55, [
        [("Non-Smooth Equilibria in Hill-Type Muscle Models", True, INK)]], size=26)
    text(s, 0.5, 0.82, 11.2, 0.35, [
        [("Alan Royce Gabriel Samuel", True, INK),
         ("   Indian Institute of Technology Madras   ·   poster 5", False, MUTED)]], size=13)
    text(s, 0.5, 1.2, 11.2, 0.65, [
        [("Prosthesis, exosuit, FES or muscle-driven robot: ", True, ACCENT),
         ("your controller assumes one model of how the limb answers a small push. Real muscle "
          "switches on 3x faster than off and resists stretch 2x harder than shortening. "
          "At rest it sits exactly where both flip.", False, INK)]], size=14)
    movie(s, c1, 1.82, 1.9, 9.7, 5.46)

    # Slide 2 (35-44 s): open-loop consequence.
    s = new_slide(
        "35-44 s. Does it matter? Pick the wrong one and a small movement is mispredicted "
        "by up to 75 percent. [Clip: mismatch study at 2.5x speed, 9 s. Click at 44 s.]")
    text(s, 0.5, 0.3, 11.2, 1.1, [
        [("Does it matter?  ", True, INK),
         ("Pick the wrong model and a small movement is mispredicted by 14 to 75%", True, ACCENT)]], size=24)
    text(s, 0.5, 1.28, 11.2, 0.4, [
        ["Same tiny command, four textbook-valid models, one true answer: only the matched branch tracks it"]], size=14, color=MUTED)
    movie(s, c2, 1.77, 1.8, 9.8, 5.51)

    # Slide 3 (44-56 s): closed loop, raising and lowering side by side.
    s = new_slide(
        "44-56 s. Feedback keeps you stable. But average the models, which is what smoothing "
        "does, and a small lift overshoots eight times more. The fix is free: relinearize at "
        "the measured state. [Clip: raising then lowering (right) steps, four rules, "
        "2x speed, 12 s. Click at 56 s.]")
    text(s, 0.5, 0.3, 11.2, 1.1, [
        [("Feedback keeps you stable.  ", True, INK),
         ("Averaging the models, which is what smoothing does, gives 8 to 9x the overshoot", True, ACCENT)]], size=24)
    text(s, 0.5, 1.28, 11.2, 0.4, [
        ["The fix is free: relinearize at the measured state. Same MPC, same weights, no new sensors, one QP per step"]], size=14, color=MUTED)
    movie(s, c3, 2.22, 1.8, 8.9, 5.0)
    text(s, 0.5, 6.85, 6.0, 0.55, [
        [("Small lift (first). ", True, ACCENT),
         ("The averaged (blended) model overshoots 0.65 deg; the others 0.07 to 0.08. Same controller, "
          "only the model changed.", False, INK)]], size=12)
    text(s, 6.83, 6.85, 6.0, 0.55, [
        [("Small lowering (second). ", True, ACCENT),
         ("Reading the switch at the reference posture: 0.37 deg overshoot. Relinearizing at the "
          "measured state: 0.02 deg.", False, INK)]], size=12)

    # Slide 4 (56-60 s): antagonist pair and the invitation.
    s = new_slide(
        "56-60 s. With two muscles, co-contraction shrinks the problem. Poster five: come and "
        "tell me where this bites in your system. [Clip: co-contraction sweep at 3x, 6 s.]")
    text(s, 0.5, 0.3, 11.2, 1.1, [
        [("Two muscles: ", True, INK),
         ("co-contraction, the way people stiffen a joint, shrinks the velocity switch", True, ACCENT)]], size=24)
    text(s, 0.5, 1.28, 11.2, 0.4, [
        ["Brachialis + triceps: three switches, eight models, and a knob the user already has; the activation switch remains"]], size=14, color=MUTED)
    movie(s, c4, 2.31, 1.8, 8.71, 4.9)
    text(s, 0.5, 6.85, 12.3, 0.5, [
        [("Poster 5: come and tell me where this bites in your system.  ", True, ACCENT),
         ("Alan Royce Gabriel Samuel, IIT Madras", False, MUTED)]], size=17)

    prs.save(OUT)
    print("saved", OUT, f"{OUT.stat().st_size / 1e6:.1f} MB")
    for c in (c1, c2, c3, c4):
        print(f"  {c.name}: {duration(c):.1f} s")


if __name__ == "__main__":
    main()
