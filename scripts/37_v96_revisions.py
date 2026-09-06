"""v9.5 -> v9.6: implement second (statistical) review.
Hard fixes: weights statement; stale expression sentence; cg07495363 28/28 clarification;
fully-compliant definition conflict; cross-cohort AUC removal.
Tone-down: title, abstract, conclusions, population-matched -> aware, pan-compartment
qualifiers, mandatory/striking wording, screening -> detection/triage, three-concept
distinction, code availability at submission, TCGA unpaired design note, GSE223817
annotation caveat, ZSCAN12 EPIC-compatible substitutes (new Results paragraph).
"""
from pathlib import Path

from docx import Document

BASE = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19")
SRC = BASE / "Manuscript_EC_methylation_panel_v9.5.docx"
DST = BASE / "Manuscript_EC_methylation_panel_v9.6.docx"

doc = Document(SRC)
log = []

def set_text(p, text):
    for r in list(p.runs):
        r._element.getparent().remove(r._element)
    p.add_run(text)

def try_all(old, new, tag, first_only=True):
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

# ---------- 1. Title ----------
try_all(
    "In silico prioritization of BOLL and ZSCAN12 as candidate DNA methylation markers for premenopausal endometrial cancer detection: a multi-cohort, population-matched background triage analysis",
    "Population- and compartment-aware in silico prioritization of DNA methylation markers for endometrial cancer: BOLL and ZSCAN12 as candidates for premenopausal validation",
    "title")

# ---------- 2. Abstract Results (compressed, EPIC substitutes in, AUC note) ----------
for p in doc.paragraphs:
    if p.text.startswith("Results: TCGA and GSE67116 yielded"):
        set_text(p,
            "Results: TCGA and GSE67116 yielded 23,099 and 292 DMRs; 50 DMRs and 476 elite probes passed "
            "cross-cohort consistency and entry criteria. Background estimates were higher in symptomatic "
            "controls than in healthy donors (approximately 1.2- to 1.6-fold on a uniform EPIC/GMQN pipeline, "
            "a descriptive cross-cohort contrast; the raw two- to three-fold difference was inflated by "
            "processing). Two probes passed every assessable background dimension: cg24589459 (BOLL; "
            "background P95 β 0.087–0.141, Δβ=0.533, premenopausal positivity 89.3%) and cg27577527 "
            "(ZSCAN12; P95 β 0.044–0.060 where measurable, Δβ=0.614, positivity 92.9%). The ZSCAN12 probe "
            "is absent from the EPIC platform, so its symptomatic-benign background could not be assessed "
            "directly, but three EPIC-compatible probes within the same DMR passed all background gates. "
            "These in silico estimates derive from the discovery cohort and should be regarded as upper "
            "bounds pending independent validation. Literature convergence was substantial: ZSCAN12 is the "
            "core of WID-qEC, CDO1 anchors the approved Chinese kit and showed premenopausal sensitivity in "
            "CISENDO, and GHSR is a component of the Dutch nine-marker panel, whereas BOLL has no prior "
            "report for EC detection. Sampling-compartment analysis identified three pan-compartment "
            "candidates on assessable dimensions (ZSCAN12, BOLL, ARL5C); the approved-kit genes AJAP1 and "
            "GALR1 showed elevated blood-cell background unsuitable for cervical self-sampling. Probe-level "
            "external validation in independent cohorts partially supported the panel and revealed a "
            "dissociation within the BOLL region—moderate tumour positivity of the anchor probe with an "
            "exceptionally clean control background—supporting two-probe, region-level interpretation; no "
            "cross-cohort ROC metric is reported, because case and control series come from different "
            "studies. A quantitative hyperplasia-to-carcinoma methylation gradient was observed at panel "
            "loci in two cohorts.")
        log.append("OK   abstract-results")
        break

# ---------- 3. Abstract Conclusions ----------
for p in doc.paragraphs:
    if p.text.startswith("Conclusions: Marker background is cohort- and compartment-dependent"):
        set_text(p,
            "Conclusions: Marker background is cohort- and compartment-dependent; benign controls should "
            "come from the target clinical population. BOLL and ZSCAN12 are prioritised as candidates for "
            "premenopausal EC detection, with a sampling-compartment profile compatible in principle with "
            "cervical or self-collected formats; they do not yet constitute a validated panel. Confirmation "
            "requires wet-lab validation in a purpose-built premenopausal AUB cohort with benign "
            "uterine-condition controls.")
        log.append("OK   abstract-conclusions")
        break

# ---------- 4. Keywords ----------
try_all("premenopausal; screening; sampling compartment",
        "premenopausal; diagnostic triage; sampling compartment", "keywords")

# ---------- 5. TCGA unpaired design ----------
try_all(
    "Contrasts: TCGA tumour versus adjacent normal; GSE67116 primary tumour versus hyperplasia.",
    "Contrasts: TCGA tumour versus adjacent normal (unpaired design, 431 tumours versus 46 adjacent "
    "normals, without adjustment for tumour purity or molecular subtype, because effect-size filtering "
    "rather than statistical significance was the primary criterion); GSE67116 primary tumour versus "
    "hyperplasia.",
    "tcga-design")

# ---------- 6. Weights ----------
try_all(
    "A weighted composite score combined benign background (25%), premenopausal tumour positivity (25%), "
    "overall tumour positivity (20%), cross-cohort consistency (15%) and amplicon suitability (5%).",
    "A weighted composite score combined benign background (weight 0.25), premenopausal tumour positivity "
    "(0.25), overall tumour positivity (0.20), cross-cohort consistency (0.15) and amplicon suitability "
    "(0.05). The consistency and amplicon terms were constant across candidates at this stage (all "
    "candidates had already passed cross-cohort consistency and amplicon design), so ranking was driven by "
    "the first three terms; the score is an unnormalised prioritisation index (maximum 0.90), not an "
    "absolute performance metric.",
    "weights")

# ---------- 7. GSE223817 annotation caveat ----------
try_all(
    "(ii) benign endometrium from symptomatic women and endometriosis patients (GSE223817) [19];",
    "(ii) benign endometrium from symptomatic women and endometriosis patients (GSE223817) [19]—per-sample "
    "menopausal status, AUB presentation and benign uterine diagnoses (polyps, leiomyomas, adenomyosis) are "
    "not annotated in this cohort, which therefore approximates, rather than strictly matches, the "
    "intended-use population;",
    "gse223817-caveat")

# ---------- 8. Methods: three independent cohorts wording ----------
try_all(
    "External validation was performed at the probe level in three independent EPIC-array cohorts not used in discovery:",
    "External validation was performed at the probe level in three independent cohorts with EPIC-array "
    "data, none used in discovery:",
    "methods-epic-wording")

# ---------- 9. P34 descriptive qualifier ----------
try_all(
    "i.e., approximately 1.2- to 1.6-fold;",
    "i.e., approximately 1.2- to 1.6-fold; we treat this as a descriptive cross-cohort comparison rather "
    "than a confirmed population-level difference, because cohort, centre and residual processing effects "
    "cannot be fully separated and the cancer-free reference comprises only 13 samples;",
    "descriptive-qualifier")

# ---------- 10. Union-rule naming ----------
try_all(
    "A union rule requiring P95 β≤0.15 in every cohort left two fully compliant probes.",
    "A union rule requiring P95 β≤0.15 in every assessable cohort (hereafter 'five-dimension "
    "union-compliant') left two probes.",
    "union-naming")

# ---------- 11. cg07495363 28/28 + 99.8% qualifier ----------
try_all(
    "Greedy coverage among compliant probes showed that a BOLL-region probe alone detected all 28 "
    "premenopausal tumours and 94.0% of all TCGA tumours in silico, with overall coverage reaching 99.8% "
    "after two additions.",
    "Greedy coverage among compliant probes showed that the high-sensitivity BOLL-region probe cg07495363 "
    "alone detected all 28 premenopausal tumours and 94.0% of all TCGA tumours in silico (the anchor probe "
    "cg24589459 detected 25 of 28, 89.3%), with overall coverage reaching 99.8% after two additions—a "
    "discovery-cohort estimate that is optimistic by construction.",
    "greedy-clarify")

# ---------- 12. 23-probe pool rewording ----------
try_all(
    "A 23-probe panel (15 fully compliant plus 8 backups; full probe list, genomic coordinates and "
    "primer-design annotations in Supplementary Data), all SNP-free and 22/23 embedded in multi-CpG DMRs "
    "(263 bp–3.4 kb), is proposed for validation by targeted bisulfite sequencing.",
    "A 23-probe experimental candidate pool—7 tier-A/B lead candidates plus 16 region-coverage probes, of "
    "which 15 additionally passed the narrower cervix–blood–cycle triage gate (a distinct criterion from "
    "the five-dimension union rule above); full probe list, genomic coordinates and primer-design "
    "annotations in Supplementary Data—all SNP-free and 22/23 embedded in multi-CpG DMRs (263 bp–3.4 kb), "
    "is proposed for validation by targeted bisulfite sequencing.",
    "panel23-reword")

# ---------- 13. P41 pan-compartment qualifier + soften ----------
try_all(
    "a panel anchored on pan-compartment markers offers a route that existing kits cannot take.",
    "a panel anchored on pan-compartment markers offers a potential route that the current "
    "intrauterine-brush kit does not take. ZSCAN12's pan-compartment label rests on the four measurable "
    "dimensions; its symptomatic-benign endometrial background was verified only through the "
    "EPIC-compatible region probes described below.",
    "pancomp-qualifier")

# ---------- 14. P45 full-series wording ----------
try_all(
    "In the full 103-sample GSE136791 series, 13 of 21 evaluable probes reached at least 80% tumour "
    "positivity (β>0.3);",
    "In the full 103-sample GSE136791 series (carcinoma and hyperplasia combined), 13 of 21 evaluable "
    "probes reached at least 80% positivity (β>0.3);",
    "gse136791-full-series")

# ---------- 15. Results: three independent cohorts wording ----------
try_all(
    "Probe-level external validation was performed in three independent EPIC cohorts not used in discovery:",
    "Probe-level external validation was performed in three independent cohorts with EPIC-array data, none "
    "used in discovery:",
    "results-epic-wording")

# ---------- 16. Remove cross-cohort AUC ----------
try_all(
    "as an illustrative discrimination estimate, the two-probe BOLL-region score (maximum of the two "
    "probes) separated carcinomas from the 13 endometrial controls (GSE155760—the only cancer-free "
    "endometrial controls among the validation cohorts) with an AUC of 0.96 (GSE136791) and 0.98 "
    "(GSE155760); formal panel-rule definition and threshold locking are deferred to wet-lab validation;",
    "we do not report ROC/AUC estimates, because the only cancer-free endometrial controls among the "
    "validation cohorts (n=13, GSE155760) come from a different study than the GSE136791 carcinomas—so a "
    "case–control ROC would conflate disease status with cohort and batch effects—and because 13 controls "
    "are too few for a stable within-cohort estimate; formal panel-rule definition and threshold locking "
    "are deferred to wet-lab validation;",
    "auc-removal")

# ---------- 17. mandatory x2 ----------
try_all(
    "region-level (two-probe) interpretation should therefore be regarded as mandatory rather than "
    "optional in the wet-lab panel.",
    "we therefore pre-specify region-level (two-probe) interpretation in the wet-lab panel.",
    "mandatory-1")
try_all(
    "further supporting mandatory region-level interpretation;",
    "further supporting region-level interpretation;",
    "mandatory-2")

# ---------- 18. New Results paragraph: ZSCAN12 EPIC-compatible substitutes ----------
NEW_PARA = (
    "To address the EPIC-platform absence of cg27577527, we screened all nine EPIC-manifest probes "
    "spanning the same ZSCAN12 DMR (hg38 chr6:28,399,501–28,400,362, which contains cg27577527 at "
    "28,399,766) against the decisive symptomatic-benign background (GSE223817) and the external tumour "
    "cohorts. Three EPIC-compatible probes—cg25666433, cg23164203 and cg20275132—passed the P95 β≤0.15 "
    "gate in symptomatic controls (P95 0.088, 0.120 and 0.099; positivity ≤0.3%) and showed low background "
    "across all remaining dimensions (healthy premenopausal endometrium P95 0.143–0.151; cervical scrapes "
    "0.060–0.095; blood cells 0.080–0.092; menstrual-cycle drift 0.023–0.029), while retaining high tumour "
    "positivity (TCGA 86.8–90.3% overall and 89.3–96.4% premenopausal; GSE136791 carcinoma 91.3–95.7%; "
    "GSE155760 endometrioid EC 65.2–73.9%; full per-cohort metrics in Supplementary Data). The ZSCAN12 "
    "region therefore remains a viable pan-platform candidate: its symptomatic-benign compliance, not "
    "directly assessable for cg27577527 itself, is supported at the region level by EPIC-compatible "
    "readouts that can be incorporated into the targeted assay.")
inserted = False
for p in doc.paragraphs:
    if p.text.strip().startswith("Hyperplasia-to-carcinoma gradient and BOLL expression analysis"):
        p.insert_paragraph_before(NEW_PARA)
        inserted = True
        log.append("OK   zscan12-substitutes-para")
        break
if not inserted:
    log.append("MISS zscan12-substitutes-para")

# ---------- 19. P47 cross-cohort reference caveat ----------
try_all(
    "same platform, used as a cross-cohort benign reference)",
    "same platform, used as a cross-cohort benign reference; cohort effects cannot be excluded, and the "
    "within-cohort orderings below are the primary evidence)",
    "gradient-caveat")

# ---------- 20. Discussion tone ----------
try_all("a candidate set designed for premenopausal women,",
        "a candidate set prioritised for premenopausal physiology,", "disc-tone")

# ---------- 21. Delete stale expression sentence ----------
try_all(
    "Whether this hypermethylation represses or merely marks the locus is unknown—expression-correlation "
    "analysis and functional follow-up are warranted before mechanistic claims are made.",
    "Functional follow-up is warranted before any mechanistic claim is made.",
    "stale-expression")

# ---------- 22. Limitations additions ----------
try_all(
    "Finally, no wet-lab validation has been performed;",
    "The benign endometrial controls in the external validation cohorts number only 13, so background and "
    "false-positive estimates from that group are unstable and cross-cohort contrasts involving them are "
    "descriptive. GSE223817 lacks per-sample annotation of menopausal status, AUB presentation and benign "
    "uterine diagnoses, so the symptomatic-benign dimension approximates rather than strictly matches the "
    "intended-use population. The 28 premenopausal TCGA tumours are the only explicitly premenopausal "
    "cancers analysed; external cohorts were not menopause-stratified, so any premenopausal-specific "
    "advantage remains to be demonstrated prospectively. Finally, no wet-lab validation has been "
    "performed;",
    "limitations-add")

# ---------- 23. Conclusion with three-concept distinction ----------
for p in doc.paragraphs:
    if p.text.startswith("In conclusion, multi-cohort discovery"):
        set_text(p,
            "In conclusion, multi-cohort discovery with population- and compartment-aware background "
            "triage prioritised BOLL and ZSCAN12 as the most robust in silico candidates for premenopausal "
            "EC detection and reconciled the composition of existing tests with their sampling designs. We "
            "distinguish three outputs of this study: two lead candidate regions (BOLL and ZSCAN12), a "
            "23-probe experimental pool for wet-lab adjudication, and a final locked diagnostic panel, "
            "which does not yet exist and must be defined in a purpose-built premenopausal AUB cohort. "
            "Within these bounds, the framework and tiered candidates provide a reproducible foundation "
            "for wet-lab validation toward a premenopausal-suitable, minimally invasive EC detection test.")
        log.append("OK   conclusion")
        break

# ---------- 24. Code availability ----------
try_all(
    "Analysis scripts, the DMR pipeline and the reference-correction workflow are available from the "
    "corresponding author upon reasonable request and will be deposited in a public repository upon "
    "acceptance.",
    "Analysis scripts, the DMR pipeline (including the comb-p-like region-merging implementation) and the "
    "key result tables are provided as Supplementary Material with this submission and will be deposited "
    "in a public repository upon acceptance.",
    "code-availability")

# ---------- 25. background screening -> triage ----------
try_all(
    "the audit demonstrates the attrition of tumour-signal-driven candidates under clinically relevant "
    "background screening.",
    "the audit demonstrates the attrition of tumour-signal-driven candidates under clinically relevant "
    "background triage.",
    "screening-word")

doc.save(DST)
print("\n".join(log))
print("saved:", DST.name)
