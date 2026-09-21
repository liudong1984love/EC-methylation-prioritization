"""v9.16 -> v9.17: 11th review round.
Fix real numeric errors (Table 8 ctrl13 column, 0.025), unify pipeline language,
title/positioning, WID-qEC hedging, BOLL OR-rule full-cohort disclosure,
pairing caveat, surrogate independence, figure legends harmonised.
"""
from docx import Document

SRC = r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Manuscript_EC_methylation_panel_v9.16.docx"
DST = r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Manuscript_EC_methylation_panel_v9.17.docx"

doc = Document(SRC)
log = []

def set_text(p, new):
    if not p.runs:
        p.add_run(new)
        return
    p.runs[0].text = new
    for r in p.runs[1:]:
        r.text = ""

def sub(p, old, new, tag):
    if old in p.text:
        set_text(p, p.text.replace(old, new))
        log.append("OK   " + tag)
        return True
    log.append("MISS " + tag)
    return False

# ---------------- Title ----------------
for p in doc.paragraphs:
    if p.text.startswith("Population- and compartment-informed in silico prioritization"):
        set_text(p, "Population- and compartment-informed prioritization of DNA methylation markers "
                    "for endometrial cancer: an in silico evaluation of BOLL and ZSCAN12")
        log.append("OK   title")
        break

# ---------------- Abstract ----------------
for p in doc.paragraphs:
    t = p.text
    if t.startswith("Background:"):
        sub(p, "We aimed to build and triage methylation candidates for premenopausal EC.",
            "We aimed to prioritise methylation candidates for future validation in premenopausal EC triage.", "abstract aim")
    if t.startswith("Methods:"):
        sub(p, "benign endometrium from a surgical cohort used as an intended-use proxy (menopausal status and AUB presentation unavailable)",
            "benign endometrium from a surgical cohort (individual menopausal status and AUB presentation unavailable)", "abstract methods proxy")
    if t.startswith("Results:"):
        sub(p, "Background estimates were descriptively higher in the intended-use-proxy cohort than in a cancer-free endometrial reference cohort.",
            "Background estimates were descriptively higher in the benign surgical cohort than in cancer-free endometrial reference controls.", "abstract results proxy")

# ---------------- Methods ----------------
for p in doc.paragraphs:
    t = p.text
    if t.startswith("Candidates were triaged against five dimensions"):
        sub(p, "(ii) benign endometrium from an endometriosis-enriched surgical cohort approximating the intended-use population (GSE223817) [19]—per-sample menopausal status, AUB presentation and benign uterine diagnoses (polyps, leiomyomas, adenomyosis) are not annotated, so this cohort approximates rather than strictly matches the intended-use population; (iii) paired early/mid-secretory endometrium (GSE90060; 17 donor pairs reconstructed by methylation fingerprinting, within-pair r=0.991);",
            "(ii) benign endometrium from an endometriosis-enriched surgical cohort (GSE223817) [19], used as an approximate intended-use proxy—per-sample menopausal status, AUB presentation and benign uterine diagnoses (polyps, leiomyomas, adenomyosis) are not annotated, so this cohort approximates rather than strictly matches the intended-use population; (iii) early/mid-secretory endometrium (GSE90060; the source publication reports paired donor sampling [13], but donor identifiers are not annotated in the accessed records, so 17 pairs were reconstructed by genome-wide methylation similarity—within-pair r median 0.992, best cross-donor r 0.991, 15/17 pairs age-concordant—and drift estimates were confirmed on the age-concordant subset);", "methods dimensions")
    if t.startswith("Redundant probes (tumour β correlation"):
        sub(p, "Literature convergence was evaluated by querying Europe PMC for each candidate and recent EC methylation panels [8,23,24].",
            "Literature convergence was evaluated by querying Europe PMC for each candidate and recent EC methylation panels [8,23,24] (last searched 6 September 2026).", "methods lit date")
    if t.startswith("Analyses used R 4.6.0"):
        sub(p, "Discovery-layer β values combined Xena-precomputed TCGA data with noob-normalised GSE67116 data, whereas background triage used GMQN-normalised values; because of this mismatch all thresholds were set empirically, and final thresholds will be established in a prospectively collected cohort processed with one predefined pipeline (see Limitations).",
            "The primary pipeline was prespecified per data source: noob-normalised IDAT processing for GSE73949 and GSE67116; GMQN-normalised β values as distributed by the NGDC EWAS Data Hub for the other GEO cohorts; and Xena-precomputed β values for TCGA-UCEC. All background estimates in the main text use these primary pipelines; raw/GMQN reprocessing of the GSE73949 and GSE67116 IDATs was used only in pipeline-sensitivity analyses. Because pipelines differ across sources, all thresholds were set empirically, and final thresholds will be established in a prospectively collected cohort processed with one predefined pipeline (see Limitations).", "methods primary pipeline")
        sub(p, "Thirty-three patient-matched tumour–normal pairs were analysed by Wilcoxon signed-rank tests",
            "Thirty-three tumour–normal pairs matched by TCGA patient barcode were analysed by Wilcoxon signed-rank tests", "methods tcga barcode")

# ---------------- Results ----------------
for p in doc.paragraphs:
    t = p.text
    if t.startswith("Background P95 β values in the benign endometrial cohort"):
        sub(p, "in the benign endometrial cohort approximating the intended-use population (n=347)",
            "in the benign surgical endometrial cohort (n=347)", "results cohort name")
        sub(p, "for cg27577527 (ZSCAN12), the intended-use-proxy dimension was not directly assessable",
            "for cg27577527 (ZSCAN12), the benign-surgical dimension was not directly assessable", "results zscan12 dim")
    if t.startswith("The two lead probes differ in evidentiary status"):
        sub(p, "(P95 β 0.044–0.060; Δβ=0.614, same contrast; positivity 92.9%)",
            "(P95 β 0.044–0.060 under the primary pipeline of each dimension; Δβ=0.614, same contrast; positivity 92.9%)", "zscan12 range pipeline")
        sub(p, "the decisive intended-use-proxy dimension is not directly assessable",
            "the decisive benign-surgical-cohort dimension is not directly assessable", "zscan12 dim2")
        sub(p, "Menstrual-cycle drift of the anchor probe was minimal (mean |Δβ| 0.026 between fingerprint-paired early- and mid-secretory samples; experimental-pool range 0.014–0.042); proliferative-phase, menstrual and anovulatory states were not assayed.",
            "Menstrual-cycle drift of the anchor probe was minimal (mean |Δβ| 0.026 between fingerprint-inferred early/mid-secretory pairs, 0.028 on the 15 age-concordant pairs; experimental-pool range 0.014–0.042, concordant-subset 0.013–0.044; pairing is inference—see Methods); proliferative-phase, menstrual and anovulatory states were not assayed.", "cycle drift concordant")
    if t.startswith("Cross-pipeline consistency."):
        sub(p, "changed pass/fail status at the 0.15 gate for 11 of 26 probes between noob and GMQN—all near-gate probes (P95 0.125–0.178 across pipelines).",
            "changed pass/fail status at the 0.15 gate for 11 of 26 probes between noob and GMQN—all near-gate probes (P95 0.125–0.178 across pipelines). We therefore report gate compliance in three grades: stable (point-estimate pass under every tested pipeline; both lead probes), pipeline-dependent (the 11 near-gate probes, including the three ZSCAN12 substitutes) and not passing; only the stable grade supports pipeline-independent candidate claims.", "grading")
        sub(p, "(P95 0.160–0.167)", "(P95 0.158–0.169)", "subs gmqn range p46")
    if t.startswith("Of the 23 experimental-pool probes, 21 are present"):
        new49 = ("Of the 23 experimental-pool probes, 21 are present on the EPIC manifest (exceptions cg27577527 [ZSCAN12] and cg20998319); cg27577527 was therefore not evaluable externally. "
        "In GSE136791, 13 of 21 evaluable probes reached ≥80% positivity (β>0.3) across all 103 samples and 16 of 21 in the 69 carcinomas, including six at ≥94% and the BOLL anchor cg24589459 at 87.0% (Table 8). "
        "In the 23 endometrioid EC samples of GSE155760, positivity was broadly reproducible (91.3–95.7%) but lower for cg24589459 (65.2%), whereas the 13 endometrial controls showed a markedly lower background for cg24589459 (P95 β=0.059) than for the other tiered candidates (P95 0.177–0.488), including cg07495363 (P95 β=0.266); no control exceeded β 0.3 for either BOLL probe (0/13) (Table 8). "
        "This dissociation indicates complementary behaviour within the region—cg24589459 low-background, cg07495363 high-sensitivity—and constitutes a post hoc, hypothesis-generating observation made after inspecting these external data. "
        "As a post hoc rule derived from these same data (reported for transparency, not as a performance claim), positivity of either BOLL-region probe yielded 65/69 (94.2%) in GSE136791 carcinomas and 22/23 (95.7%) in GSE155760 with 0/13 positives among the endometrial controls; we do not report ROC/AUC, because the only cancer-free controls (n=13) come from a different study than the GSE136791 carcinomas and are too few for a specificity estimate. "
        "Critically, the same OR rule is frequently positive in benign endometrium—159/347 (45.8%) in the benign surgical cohort and 214/637 (33.6%) in its endometriosis subset, driven by cg07495363 (anchor alone 2/347 and 0/637)—although it was negative in every healthy-endometrium (0/17), HPV-negative cervical (0/20), purified blood-cell (0/42) and endometrial-control (0/13) sample. "
        "A simple β>0.3 OR rule therefore inherits the partner probe's benign-tissue background, and the 0/13 result cannot be read as evidence of high specificity; any two-probe region-level rule must be developed and locked in an independent threshold-setting cohort containing symptomatic benign samples before any validation. "
        "In GSE178610, positivity was preservation-dependent (higher in FFPE than fresh-frozen; Supplementary Data), consistent with fixation- or batch-related inflation; fresh-frozen estimates were used for cross-cohort comparison. "
        "cg27577527 was positive in 8/9 tumours (mean β=0.62) in a small 450K cohort not used in discovery (GSE93589). "
        "In the 34 GSE136791 hyperplasia samples, positivity remained high (58.8–100%), indicating a neoplasia-gradient rather than carcinoma-specific profile (atypia unannotated; see below).")
        set_text(p, new49)
        log.append("OK   p49 rewrite")
    if t.startswith("To address the EPIC-platform absence of cg27577527"):
        sub(p, "against the decisive intended-use-proxy background (GSE223817)",
            "against the decisive benign-surgical-cohort background (GSE223817)", "subs background name")
        sub(p, "passed the P95 β≤0.15 gate in the approximating cohort",
            "passed the P95 β≤0.15 gate in the benign surgical cohort", "subs cohort name")
        sub(p, "(P95 0.160–0.167)", "(P95 0.158–0.169)", "subs gmqn range p50")
        sub(p, "we do not interpret sub-0.01 point differences (e.g., 0.149 versus 0.151) as meaningful.",
            "we do not interpret sub-0.01 point differences (e.g., 0.149 versus 0.151) as meaningful. Because the external tumour cohorts were used to select the three substitutes, they no longer constitute fully independent validation for them.", "subs independence")

# ---------------- Discussion ----------------
for p in doc.paragraphs:
    t = p.text
    if t.startswith("Methylation-based EC detection has entered clinical reality"):
        sub(p, "a candidate set prioritised for premenopausal physiology",
            "a candidate set prioritised for future validation in premenopausal AUB triage", "disc opening 1")
        sub(p, "a reproducible audit framework that reconciles our candidates with—and explains the design choices of—existing tests",
            "a reproducible audit framework that generates testable hypotheses about marker choice and sampling design", "disc opening 2")
    if t.startswith("Re-identification of established markers supports"):
        sub(p, "These probe-level findings do not evaluate the exact CpGs, amplicons or classification algorithms used in the commercial or published assays: WID-qEC absorbs GYPC's probe-level background through region-level ΣPMR quantification, and the approved kit deploys AJAP1/GALR1 via an intrauterine brush with cohort-specific cut-offs.",
            "These probe-level findings do not evaluate the exact CpGs, amplicons or classification algorithms used in the commercial or published assays. One plausible—but untested—interpretation is that region-level quantification (such as ΣPMR) tolerates probe-level background that eliminates single probes, and that intrauterine-brush deployment bypasses blood-cell-dominated compartments; we did not assay the actual amplicons or classifiers, and these probe-level observations cannot be used to evaluate the performance of existing test products.", "wid-qec hedge")
    if t.startswith("No EC-detection report for BOLL"):
        sub(p, "was identified in our literature search.",
            "was identified in our literature search (Europe PMC; last searched 6 September 2026).", "boll search date")
    if t.startswith("The cohort dependence we quantified"):
        sub(p, "background inflation in the intended-use-proxy cohort versus healthy endometrium",
            "background inflation in the benign surgical cohort versus healthy endometrium", "disc cohort name")
    if t.startswith("Several limitations warrant emphasis"):
        set_text(p, p.text.rstrip() +
            " Menstrual-cycle pairing in GSE90060 was inferred from genome-wide methylation similarity because donor identifiers are not annotated; within-pair and best cross-donor correlations are close (median 0.992 versus 0.991), so pairing is probabilistic, although drift estimates changed little on the 15 age-concordant pairs. "
            "The post hoc two-probe BOLL rule was derived from the same external data and inherits the partner probe's benign-endometrium background (45.8% rule positivity in the benign surgical cohort); it is a hypothesis to be locked and tested in a threshold-setting cohort, not a validated rule.")
        log.append("OK   limitations append")
    if t.startswith("These limitations define a concrete validation agenda"):
        sub(p, "and prespecified handling of atypical hyperplasia.",
            "and prespecified enrolment of both atypical and non-atypical hyperplasia in adequate numbers.", "agenda hyperplasia")
    if t.startswith("In conclusion, multi-cohort discovery"):
        sub(p, "and reconciled the composition of existing tests with their sampling designs",
            "and generated hypotheses linking the composition of existing tests to their sampling designs", "conclusion hedge")

# ---------------- Figure & table legends ----------------
for p in doc.paragraphs:
    t = p.text
    if t.startswith("Figure 3."):
        set_text(p, "Figure 3. Sampling-compartment suitability of candidate markers (data from Table 7; bar colours use the Table 7 categories). Green: P95 β ≤0.15 (✓); yellow: >0.15–0.30 (△); red: >0.30 (✗). *Borderline: point-estimate pass whose bootstrap 95% confidence interval crosses the 0.15 gate (BOLL anchor).")
        log.append("OK   fig3 legend")
    if t.startswith("Figure 4."):
        set_text(p, "Figure 4. Cross-cohort performance of the BOLL and ZSCAN12 regions. (A) Background P95 β of the three region probes across background cohorts (primary pipeline per cohort; dashed line: 0.15 gate; n/a: probe absent from the platform). (B) Tumour-tissue positivity (β>0.3) across discovery and corroboration cohorts (fresh-frozen stratum shown for GSE178610). The benign-surgical-cohort background of cg07495363 and the lower external positivity of cg24589459 illustrate the complementary behaviour motivating region-level assay design.")
        log.append("OK   fig4 legend")
    if t.startswith("Table 7. Sampling-compartment suitability"):
        set_text(p, t.rstrip() + " Figure 3 uses the same colour categories.")
        log.append("OK   table7 footnote")

# ---------------- Table 8 ctrl13 column ----------------
FIX = {"cg18507379": "0.352", "cg10109500": "0.319", "cg15790037": "0.177",
       "cg23180938": "0.252", "cg07495363": "0.266", "cg24589459": "0.059", "cg16439198": "0.488"}
t8 = doc.tables[7]
fixed = 0
for row in t8.rows:
    probe = row.cells[0].text.strip()
    if probe in FIX:
        set_text(row.cells[5].paragraphs[0], FIX[probe])
        fixed += 1
log.append(f"OK   table8 ctrl13 cells fixed: {fixed}/7")

doc.save(DST)
print("saved:", DST)
for l in log:
    print(l)
print("MISSES:", sum(1 for l in log if l.startswith("MISS")))
