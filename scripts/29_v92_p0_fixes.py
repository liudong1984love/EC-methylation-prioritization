"""v9.2: respond to P0 review with data-backed edits."""
import re
from pathlib import Path

from docx import Document

BASE = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19")
SRC = BASE / "Manuscript_EC_methylation_panel_v9.1.docx"
DST = BASE / "Manuscript_EC_methylation_panel_v9.2.docx"

doc = Document(SRC)
log = []

def set_text(p, text):
    for r in list(p.runs):
        r._element.getparent().remove(r._element)
    p.add_run(text)

def sub_para(p, old, new, tag):
    if old in p.text:
        set_text(p, p.text.replace(old, new))
        log.append(f"OK   {tag}")
        return True
    return False

def try_all(old, new, tag, first_only=True):
    n = 0
    for p in doc.paragraphs:
        if sub_para(p, old, new, tag):
            n += 1
            if first_only:
                break
    if n == 0:
        log.append(f"MISS {tag}: {old[:60]!r}")

# ---------- P0-1：队列依赖数字分解 ----------
try_all(
    "Background estimates were two- to three-fold higher in symptomatic controls than in healthy donors.",
    "Background estimates were higher in symptomatic controls than in healthy donors; although the raw cross-platform contrast suggested a two- to three-fold difference, a decomposition controlling for platform and normalisation (uniform EPIC/GMQN pipeline) confirmed a consistent elevation of smaller magnitude (e.g., cancer-free vs symptomatic endometrial controls: P95 beta 0.059 vs 0.100 for cg24589459; 0.266 vs 0.427 for cg07495363; 0.252 vs 0.372 for CDO1).",
    "abstract-2to3x")

try_all(
    "Background P95 β values in benign endometrium from symptomatic women (n=347) were two- to three-fold higher than in healthy premenopausal donors (n=17), with a continuous distribution from 0.05 to 0.70;",
    "Background P95 β values in benign endometrium from symptomatic women (n=347) were higher than in healthy premenopausal donors (n=17), with a continuous distribution from 0.05 to 0.70. Because this contrast spans platform (450K versus EPIC) and normalisation (noob versus GMQN), we decomposed it: on a uniform EPIC/GMQN pipeline, symptomatic controls remained consistently elevated relative to cancer-free endometrial controls (n=13; e.g., P95 β 0.100 vs 0.059 for cg24589459, 0.427 vs 0.266 for cg07495363, 0.372 vs 0.252 for CDO1, 0.425 vs 0.319 for GHSR), i.e., approximately 1.2- to 1.6-fold; normalisation alone accounted for only +0.01-0.04 β on the same healthy samples. The qualitative conclusion is therefore platform-independent, although the magnitude of the raw two- to three-fold contrast was inflated by processing differences;",
    "results-decompose")

try_all(
    "The cohort dependence we quantified (two- to three-fold background inflation in symptomatic versus healthy endometrium)",
    "The cohort dependence we quantified (background inflation in symptomatic versus healthy endometrium; 1.2- to 1.6-fold on a uniform processing pipeline, up to two- to three-fold in the raw cross-platform contrast)",
    "discussion-2to3x")

# ---------- P0-2：ZSCAN12 小队列验证 + 增生佐证 ----------
try_all(
    "cg27577527 (ZSCAN12) remains unevaluable on EPIC and requires EPIC-compatible probe designs.",
    "cg27577527 (ZSCAN12) remains unevaluable on EPIC and requires EPIC-compatible probe designs; as a partial substitute, we assessed it in an additional small 450K cohort not used in discovery (GSE93589, n=9 EC), where it was positive in 8/9 tumours (mean β=0.62), and the two BOLL-region probes were positive in 9/9 each. Independently, the 8 hyperplasia samples of GSE67116 (450K) also showed elevated background at these loci (P95 β 0.44–0.69), corroborating the neoplasia-gradient pattern observed in GSE136791.",
    "gse93589-zscan12")

# ---------- P0-3：OR 规则特异度余量声明 ----------
try_all(
    "formal panel-rule definition and threshold locking are deferred to wet-lab validation.",
    "formal panel-rule definition and threshold locking are deferred to wet-lab validation; we note that the false-positive estimate rests on 13 controls (0/13; upper 95% confidence bound ≈23%) and that the partner probe's control P95 (0.270) sits close to the illustrative threshold, so the specificity margin of this rule requires confirmation in larger benign series.",
    "or-rule-caveat")

# ---------- P0-5/其他：局限性增补 ----------
try_all(
    "Several limitations warrant emphasis.",
    "Several limitations warrant emphasis. The healthy-donor background estimates rest on 17 samples, so P95 statistics are unstable (approximating the third-largest value); menopause status in TCGA was taken from the UCSC Xena clinical matrix as clinically annotated and was not histologically confirmed.",
    "limitations-n17")

# ---------- 小问题 ----------
try_all("an NMPA-approved three-gene adjunctive kit and international tests (WID-qEC; CISENDO) now exist",
        "an NMPA-approved three-gene adjunctive kit and international programmes (the WID-qEC test; the CISENDO cohort) now exist",
        "cisendo-wording")
try_all("and and cg27577527", "and cg27577527", "typo-andand")

# BOLL 生物学合理性（讨论段，cohort-dependence 段之后插入由段内追加实现）
try_all(
    "a gap only purpose-built cohorts can fill.",
    "a gap only purpose-built cohorts can fill. Regarding biological plausibility, BOLL (boule homolog) encodes a germ-cell meiotic regulator that is normally silenced in somatic tissues; cg24589459 lies within a 12-CpG, ~1.8-kb differentially methylated region on chromosome 2 that is hypermethylated across three independent EC cohorts. Whether this hypermethylation represses or merely marks the locus is unknown—expression-correlation analysis and functional follow-up are warranted before mechanistic claims are made.",
    "boll-plausibility")

# 23 探针明细引用
try_all(
    "A 23-probe panel (15 fully compliant plus 8 backup candidates), all free of common SNPs and 22/23 embedded in multi-CpG DMRs of 263 bp–3.4 kb, is proposed for targeted bisulfite sequencing validation",
    "A 23-probe panel (15 fully compliant plus 8 backup candidates; full probe list, genomic coordinates and primer-design annotations in Supplementary Data), all free of common SNPs and 22/23 embedded in multi-CpG DMRs of 263 bp–3.4 kb, is proposed for targeted bisulfite sequencing validation",
    "23probe-supp-ref")

# in silico / in-silico 统一为 in silico
cnt = 0
for p in doc.paragraphs:
    if "in-silico" in p.text:
        set_text(p, p.text.replace("in-silico", "in silico"))
        cnt += 1
log.append(f"OK   in-silico unify ({cnt} paras)")

doc.save(DST)
print("saved:", DST)
print("\n".join(log))
