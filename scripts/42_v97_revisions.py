# -*- coding: utf-8 -*-
"""v9.6 -> v9.7: implement third review.
New analyses inserted: menopause-stratified + interaction (A), matched-pair/endometrioid/
subtype sensitivity (B), threshold retention scan with Figure 5 (C), bootstrap/Wilson CIs (D),
bumphunter cross-check (E). Text fixes: GSE223817 renaming, threshold-contradiction correction,
pan-compartment wording, third-largest-value fix, FFPE alternative explanations, intended-use
statement for hyperplasia, 99.8% weakening, Stouffer correlation disclosure, deduplication.
"""
from pathlib import Path
from docx import Document
from docx.shared import Inches

BASE = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19")
SRC = BASE / "Manuscript_EC_methylation_panel_v9.6.docx"
DST = BASE / "Manuscript_EC_methylation_panel_v9.7.docx"

doc = Document(SRC)
log = []

def set_text(p, text):
    for r in list(p.runs):
        r._element.getparent().remove(r._element)
    p.add_run(text)

def try_all(old, new, tag, first_only=True, count_expected=None):
    n = 0
    for p in doc.paragraphs:
        if old in p.text:
            set_text(p, p.text.replace(old, new))
            n += 1
            log.append(f"OK   {tag}")
            if first_only:
                break
    if n == 0:
        log.append(f"MISS {tag}: {old[:70]!r}")
    elif count_expected and n != count_expected:
        log.append(f"WARN {tag}: expected {count_expected}, got {n}")

def replace_everywhere(old, new, tag):
    n = 0
    for p in doc.paragraphs:
        if old in p.text:
            set_text(p, p.text.replace(old, new))
            n += 1
    log.append(f"{'OK  ' if n else 'MISS'} {tag} x{n}")

# ============ 1. Abstract Methods: rename GSE223817 + mention new analyses ============
try_all(
    "benign endometrium from symptomatic women (GSE223817, 347 controls, 637 endometriosis), paired menstrual-cycle phases",
    "benign endometrium from an endometriosis-enriched cohort approximating the intended-use population (GSE223817, 347 controls, 637 endometriosis; per-sample menopausal status and AUB presentation not annotated), paired menstrual-cycle phases",
    "abs-methods-gse223817")
try_all(
    "candidates were scored for sampling-compartment suitability (intrauterine brush, cervical scrape/self-sampling, urine/ctDNA).",
    "candidates were scored for sampling-compartment suitability (intrauterine brush, cervical scrape/self-sampling, urine/ctDNA). Menopause-stratified, patient-matched paired, histology- and molecular-subtype-restricted sensitivity analyses were performed in TCGA; candidate-retention was scanned across background-stringency thresholds; lead regions were cross-checked with the permutation-based bumphunter; and background P95 and positivity estimates are reported with bootstrap and Wilson confidence intervals.",
    "abs-methods-newanalyses")

# ============ 2. Abstract Results: rename + add key new findings (kept compact) ============
try_all(
    "Background estimates were higher in symptomatic controls than in healthy donors (approximately",
    "Background estimates were higher in the benign endometrial cohort than in healthy donors (approximately",
    "abs-results-rename")
try_all(
    "Sampling-compartment analysis identified three pan-compartment candidates on assessable dimensions (ZSCAN12, BOLL, ARL5C);",
    "Three candidates (ZSCAN12, BOLL, ARL5C) showed a favourable in silico background profile across all assessed compartments;",
    "abs-results-pancomp")
try_all(
    "A quantitative hyperplasia-to-carcinoma methylation gradient was observed at panel loci in two cohorts.",
    "A quantitative hyperplasia-to-carcinoma methylation gradient was observed at panel loci in two cohorts. Tumour methylation at panel loci did not differ significantly between premenopausal (n=28) and postmenopausal (n=345) TCGA tumours, was confirmed in 33 patient-matched tumour\u2013normal pairs (all probes, FDR<0.01), and held across endometrioid-only and all four molecular-subtype strata. Strict all-dimension compliance was, however, threshold-dependent: only ZSCAN12 survived a background gate below 0.15, whereas BOLL required a gate of \u22650.15. Both lead regions were independently recovered by bumphunter (FWER=0).",
    "abs-results-newfindings")

# ============ 3. Methods: DMR paragraph — correlation disclosure + bumphunter ============
try_all(
    "Probe identities were additionally cross-checked against the Zhou et al. hg38 re-annotation",
    "Unlike comb-p, the implementation does not model inter-probe correlation when combining P values, which can be anti-conservative at region level; the two lead regions were therefore re-tested with bumphunter (B=500 permutations; maxGap 500 bp; chromosome-wide null on chr2 and chr6). Probe identities were additionally cross-checked against the Zhou et al. hg38 re-annotation",
    "methods-dmr-bumphunter")

# ============ 4. Methods: triage paragraph — rename + threshold correction ============
try_all(
    "(ii) benign endometrium from symptomatic women and endometriosis patients (GSE223817) [19]\u2014per-sample menopausal status,",
    "(ii) benign endometrium from an endometriosis-enriched surgical cohort approximating the intended-use population (GSE223817) [19]\u2014per-sample menopausal status,",
    "methods-triage-rename")
try_all(
    "candidate rankings were robust across gates of 0.10\u20130.20 (Results, threshold sensitivity)",
    "candidate rankings were stable across gates of 0.10\u20130.20, although strict all-dimension pass/fail status was threshold-dependent below 0.15 (Results, threshold sensitivity)",
    "methods-triage-threshold")

# ============ 5. Methods: statistics paragraph — new methods ============
try_all(
    "because of this processing mismatch, all thresholds (entry \u0394\u03b2\u22650.30 or positivity \u22650.50, background P95 \u03b2>0.15, weighted composite score) were set empirically, and final cut-off design and threshold sensitivity analyses are planned on uniformly processed data (see Limitations).",
    "because of this processing mismatch, all thresholds (entry \u0394\u03b2\u22650.30 or positivity \u22650.50, background P95 \u03b2>0.15, weighted composite score) were set empirically, and final cut-off design and threshold sensitivity analyses are planned on uniformly processed data (see Limitations). Background P95 confidence intervals were estimated by 10,000-draw non-parametric bootstrap over per-sample \u03b2 values; positivity rates are reported with Wilson score intervals. Menopause-stratified analyses compared 28 premenopausal with 345 postmenopausal TCGA tumours (Welch t on M-values; Fisher exact on positivity; tissue\u00d7menopause interaction by ordinary least squares). Thirty-three patient-matched tumour\u2013adjacent-normal pairs (matched on the 12-character patient barcode) were analysed by Wilcoxon signed-rank tests; molecular subtypes were taken from the PanCanAtlas annotations via cBioPortal. Candidate retention was scanned across background P95 gates of 0.05\u20130.30.",
    "methods-stats-new")

# ============ 6. Results: cohort-dependence paragraph — rename symptomatic ============
try_all(
    "Background P95 \u03b2 values in benign endometrium from symptomatic women (n=347) were higher than in healthy premenopausal donors (n=17),",
    "Background P95 \u03b2 values in the benign endometrial cohort approximating the intended-use population (n=347) were higher than in healthy premenopausal donors (n=17),",
    "results-cohortdep-rename")
try_all(
    "symptomatic controls remained consistently elevated relative to cancer-free endometrial controls",
    "the benign cohort remained consistently elevated relative to cancer-free endometrial controls",
    "results-cohortdep-rename2")

# ============ 7. Results: final candidates — 99.8% weakening + bootstrap CI ============
try_all(
    "with overall coverage reaching 99.8% after two additions\u2014a discovery-cohort estimate that is optimistic by construction.",
    "with overall coverage reaching 99.8% after two additions\u2014a within-training-set estimate that is optimistic by construction and must not be read as expected validation performance.",
    "results-998")
try_all(
    "cg24589459 (BOLL, assessable in all five cohorts; background P95 \u03b2 0.087\u20130.141, \u0394\u03b2=0.533, premenopausal positivity 89.3%)",
    "cg24589459 (BOLL, assessable in all five cohorts; background P95 \u03b2 0.087\u20130.141 [bootstrap 95% CI at the limiting cervical dimension 0.100\u20130.158], \u0394\u03b2=0.533, premenopausal positivity 89.3% [Wilson 95% CI 72.8\u201396.3%])",
    "results-boll-ci")

# ============ 8. Results: compartment — pan-compartment wording ============
try_all(
    "Only three candidates were compatible with all assessed compartments: ZSCAN12 (cervical 0.044, blood 0.058), BOLL (0.141, 0.112), and ARL5C (0.041, 0.052) (Table 7).",
    "Three candidates showed a favourable in silico background profile across all assessed compartments: ZSCAN12 (cervical 0.044, blood 0.058), BOLL (0.141, 0.112), and ARL5C (0.041, 0.052) (Table 7); this profile reflects purified-reference backgrounds and does not establish suitability for real self-collected samples, in which tumour-DNA fraction, cellular mixing and degradation dominate.",
    "results-compartment")
try_all(
    "a panel anchored on pan-compartment markers offers a potential route",
    "a panel anchored on markers with a favourable multi-compartment background profile offers a potential route",
    "results-compartment2")
try_all(
    "ZSCAN12's pan-compartment label rests on the four measurable dimensions;",
    "ZSCAN12's multi-compartment compatibility rests on the four measurable dimensions;",
    "results-compartment3")

# ============ 9. Results: threshold sensitivity — CORRECT THE CONTRADICTION ============
for p in doc.paragraphs:
    if p.text.startswith("Threshold sensitivity analyses showed that the two fully compliant candidates"):
        set_text(p,
            "Threshold sensitivity analysis distinguished ranking stability from strict gate compliance (Figure 5). "
            "Candidate rankings were stable across background gates of 0.10\u20130.20 and composite-score weight "
            "perturbations (7\u20139 of the top-10 probes unchanged), but strict all-dimension compliance was "
            "threshold-dependent: cg27577527 (ZSCAN12) passed every assessable dimension down to a gate of 0.09, "
            "whereas cg24589459 (BOLL) entered only at gates \u22650.15, limited by its cervical-scrape P95 "
            "(0.141; bootstrap 95% CI 0.100\u20130.158); at a gate of 0.10, only ZSCAN12 would have been retained. "
            "We therefore correct the earlier blanket statement that both candidates tolerate a 0.10\u20130.20 gate: "
            "rank stability and all-gate compliance are distinct properties, gate choice materially affects BOLL's "
            "pass/fail status, and final cut-offs must be set empirically in the validation cohort. The elite pool "
            "was also sensitive to extreme entry stringency (n=1,722 at \u0394\u03b2\u22650.40/positivity \u22650.80 "
            "versus n=67 at \u0394\u03b2\u22650.60/positivity \u22650.95).")
        log.append("OK   results-threshold-rewrite")
        break
else:
    log.append("MISS results-threshold-rewrite")

# ============ 10. Results: insert menopause/paired/subtype subsection after threshold para ============
NEW_RES = (
    "Menopause-stratified, matched-pair and subtype sensitivity analyses in TCGA. "
    "Three analyses tested whether the panel signal depends on menopausal status or on unpaired-design confounding. "
    "First, mean \u03b2 at the 26 lead and panel probes was essentially identical between premenopausal (n=28) and "
    "postmenopausal (n=345) tumours (per-probe difference \u22120.13 to +0.05; no probe significant after "
    "Benjamini\u2013Hochberg correction, minimum FDR 0.17), positivity rates did not differ (all Fisher P\u22650.045, "
    "none significant after correction), and tissue\u00d7menopause interaction terms were non-significant for every "
    "probe (FDR\u22650.96 for the lead probes), although the interaction analysis is underpowered because only two "
    "adjacent normals were premenopausal. Second, in 33 patient-matched tumour\u2013normal pairs all 26 probes were "
    "hypermethylated in tumours (mean paired \u0394\u03b2 0.46\u20130.70; \u226591% of pairs concordant per probe; "
    "Wilcoxon FDR<0.01 throughout), excluding confounding by patient-level factors. Third, the signal persisted in "
    "the endometrioid-only subset (n=333; \u0394\u03b2 0.44\u20130.67, FDR<0.01) and across all four TCGA molecular "
    "subtypes (POLE, MSI, copy-number-low and copy-number-high; \u0394\u03b2 0.34\u20130.71), with copy-number-high "
    "tumours showing the smallest, yet still strongly positive, effect sizes. Both lead regions were independently "
    "recovered by bumphunter from the chromosome-wide permutation null: BOLL chr2:197,785,878\u2013197,786,773 "
    "(18-probe cluster, area 22.4, FWER=0) and ZSCAN12 chr6:28,399,501\u201328,400,120 (10-probe cluster, area 33.4, "
    "FWER=0), closely matching the comb-p-like boundaries. These case-only and tissue-only comparisons show that "
    "the tumour signal is not driven by menopausal status, histology or molecular subtype; they cannot establish "
    "premenopausal assay specificity, which requires the purpose-built cohort described below."
)
for p in doc.paragraphs:
    if p.text.startswith("Threshold sensitivity analysis distinguished ranking stability"):
        newp = p.insert_paragraph_before("")
        p._p.addnext(newp._p)  # move the new empty paragraph AFTER p
        set_text(newp, NEW_RES)
        log.append("OK   insert-menopause-subsection")
        break
else:
    log.append("MISS insert-menopause-subsection")

# ============ 11. Results: external validation — 0/13 Wilson bound ============
try_all(
    "the 13 endometrial controls showed a markedly low background for cg24589459 (P95 \u03b2=0.025)",
    "the 13 endometrial controls showed a markedly low background for cg24589459 (P95 \u03b2=0.025; no control exceeded \u03b2 0.3 for either BOLL probe, 0/13, Wilson 95% CI 0\u201322.8%)",
    "results-013-wilson")

# ============ 12. Results: hyperplasia gradient — intended-use statement ============
try_all(
    "They nevertheless motivate quantitative, region-level assay readouts that could stratify risk rather than merely flag carcinoma.",
    "Because the intended use is triage of abnormal uterine bleeding\u2014in which detection of atypical hyperplasia is clinically valuable rather than a false positive\u2014we interpret hyperplasia positivity as detection of the neoplastic continuum and explicitly avoid carcinoma-specificity claims; if the intended use were instead carcinoma-versus-benign discrimination, the observed hyperplasia positivity would limit specificity. These data motivate quantitative, region-level assay readouts that could stratify risk rather than merely flag carcinoma.",
    "results-hyperplasia-intendeduse")

# ============ 13. Discussion: limitations fixes ============
try_all(
    "The healthy-donor background estimates rest on 17 samples, so P95 statistics are unstable (approximating the third-largest value);",
    "The healthy-donor background estimates rest on 17 samples, so P95 statistics are unstable and driven by the top one to two observations (bootstrap 95% CI for cg07495363 spans 0.099\u20130.203);",
    "disc-thirdlargest")
try_all(
    "consistent with formaldehyde-fixation-related methylation artefacts; fresh-frozen estimates were used for cross-cohort comparison,",
    "compatible with formaldehyde-fixation-related methylation artefacts, although cohort, tissue-composition and processing-batch differences offer alternative explanations; fresh-frozen estimates were used for cross-cohort comparison,",
    "disc-ffpe")
try_all(
    "Urine and ctDNA suitability was inferred from blood-cell background without direct validation.",
    "Urine and ctDNA suitability was inferred from blood-cell background without direct validation. The comb-p-like region merger pools same-direction significant probes without modelling inter-probe correlation and may yield optimistic region-level P values; both lead regions were therefore confirmed with the permutation-based bumphunter (FWER=0). True self-sampling performance additionally depends on tumour-DNA shedding, cellular mixing, DNA degradation and bisulfite conversion efficiency, none of which is captured by purified-reference backgrounds.",
    "disc-stouffer")

# ============ 14. Figure 3 legend wording + Figure 5 legend ============
try_all(
    "Green: P95 \u03b2 \u22640.15 in both cervix and blood (pan-compartment); red: compartment-restricted.",
    "Green: P95 \u03b2 \u22640.15 in both cervix and blood (favourable in silico background profile across assessed compartments); red: compartment-restricted.",
    "fig3-legend")
# add Figure 5 legend after Figure 4 legend
for p in doc.paragraphs:
    if p.text.startswith("Figure 4. In-silico greedy incremental coverage"):
        newp = p.insert_paragraph_before("")
        p._p.addnext(newp._p)
        set_text(newp,
            "Figure 5. Candidate retention across background-stringency thresholds. Line: number of elite probes "
            "passing every assessable background P95 dimension at each gate; step indicators: pass/fail status of "
            "the three lead probes. At the primary gate (0.15), cg24589459 (BOLL) and cg27577527 (ZSCAN12) comply; "
            "below 0.15 only ZSCAN12 remains. Rank stability across gates is distinct from strict all-gate compliance.")
        log.append("OK   fig5-legend")
        break
else:
    log.append("MISS fig5-legend")

# ============ 15. Insert Figure 5 image after its legend ============
img_path = BASE / "results" / "v97_threshold_retention_curve.png"
for p in doc.paragraphs:
    if p.text.startswith("Figure 5. Candidate retention"):
        newp = p.insert_paragraph_before("")
        run = newp.add_run()
        run.add_picture(str(img_path), width=Inches(5.8))
        log.append("OK   fig5-image")
        break
else:
    log.append("MISS fig5-image")

# ============ 16. residual 'symptomatic' cleanup (selective) ============
replace_everywhere("symptomatic-benign background could not be assessed",
                   "background in the approximating cohort could not be assessed", "cleanup-symp1")
replace_everywhere("its symptomatic-benign compliance",
                   "its compliance in the approximating cohort", "cleanup-symp2")

doc.save(DST)
print("\n".join(log))
miss = [l for l in log if l.startswith("MISS")]
print(f"\n{len(log)} ops, {len(miss)} MISS")
