"""v9.11 revisions: repositioning (7th review round).
- population-aware -> population-informed (title, conclusion)
- external validation -> independent tissue-cohort corroboration
- BOLL vs ZSCAN12 evidentiary asymmetry explicit
- post hoc labelling of two-probe rule; intended-use definition
- new: adjusted-model + cross-pipeline results; heuristic threshold framing
- GSE73949 noob recompute (bootstrap CI fix); P90/max reporting
- light compression
"""
import shutil
from docx import Document

SRC = 'Manuscript_EC_methylation_panel_v9.10.docx'
DST = 'Manuscript_EC_methylation_panel_v9.11.docx'
shutil.copy(SRC, DST)
doc = Document(DST)

LOG = []

def set_text(par, text):
    for r in list(par.runs):
        r._element.getparent().remove(r._element)
    par.add_run(text)

def replace_segment(anchor, old, new, must=True):
    for p in doc.paragraphs:
        if anchor in p.text:
            t = p.text
            if old not in t:
                if must: LOG.append(f'MISS-SEG: {old[:60]}')
                return False
            set_text(p, t.replace(old, new))
            LOG.append(f'OK seg: {old[:50]}...')
            return True
    if must: LOG.append(f'MISS-ANCHOR: {anchor[:60]}')
    return False

def replace_para(anchor, new_text):
    for p in doc.paragraphs:
        if anchor in p.text:
            set_text(p, new_text)
            LOG.append(f'OK para: {anchor[:50]}')
            return True
    LOG.append(f'MISS-PARA: {anchor[:60]}')
    return False

def insert_before(anchor, new_text):
    for p in doc.paragraphs:
        if anchor in p.text:
            np = p.insert_paragraph_before(new_text)
            np.style = p.style
            LOG.append(f'OK insert-before: {anchor[:50]}')
            return True
    LOG.append(f'MISS-INSERT: {anchor[:60]}')
    return False

def append_to(anchor, addition):
    for p in doc.paragraphs:
        if anchor in p.text:
            set_text(p, p.text + addition)
            LOG.append(f'OK append: {anchor[:50]}')
            return True
    LOG.append(f'MISS-APPEND: {anchor[:60]}')
    return False

# ---------------- 1. Title ----------------
replace_segment('in silico prioritization of DNA methylation markers for endometrial cancer',
    'Population- and compartment-aware', 'Population- and compartment-informed')

# ---------------- 2. Abstract Methods ----------------
replace_segment('Sensitivity analyses comprised',
    'Sensitivity analyses comprised menopause-stratified, patient-matched and subtype-restricted analyses in TCGA, threshold scanning, bumphunter cross-checking and bootstrap/Wilson confidence intervals.',
    'Sensitivity analyses comprised menopause-stratified, patient-matched and subtype-restricted analyses in TCGA, confounder-adjusted models (age and histological type), cross-pipeline consistency checks (raw, noob and GMQN preprocessing), threshold scanning, bumphunter cross-checking and bootstrap/Wilson confidence intervals.')

# ---------------- 3. Abstract Results ----------------
replace_segment('Two probes passed every assessable background dimension',
    'Two probes passed every assessable background dimension: cg24589459 (BOLL; Δβ=0.533, premenopausal positivity 89.3%) and cg27577527 (ZSCAN12; Δβ=0.614, positivity 92.9%). The ZSCAN12 probe is absent from the EPIC platform, but three EPIC-compatible probes within the same DMR passed (two) or marginally exceeded (one, borderline) the background gates at the region level. ZSCAN12 showed the most consistently low background across measurable dimensions, whereas BOLL met the prespecified 0.15 background gate by point estimate but its bootstrap confidence interval crossed the threshold (borderline). Probe-level external validation in three independent cohorts—performed on tumour tissue, not the intended sampling compartments—partially supported the candidates and, for BOLL, motivated a two-probe region-level interpretation to be locked prospectively; no cross-cohort ROC metric is reported because case and control series come from different studies.',
    'One probe was assessable in all five background dimensions: cg24589459 (BOLL; Δβ=0.533, premenopausal positivity 89.3%), which met the prespecified 0.15 background gate by point estimate, although its bootstrap confidence interval crossed the threshold (borderline). cg27577527 (ZSCAN12; Δβ=0.614, positivity 92.9%) passed the four dimensions measurable at probe level and showed the most consistently low background; the intended-use-proxy dimension is not directly assessable because the probe is absent from the EPIC platform, and support on that dimension is indirect and region-level only, through EPIC-compatible substitute probes whose own background compliance was pipeline-dependent. Probe-level corroboration in three independent tumour-tissue cohorts—not the intended sampling compartments—partially supported the candidates and, for BOLL, generated a post hoc two-probe region-level hypothesis to be locked prospectively; no cross-cohort ROC metric is reported because case and control series come from different studies.')
replace_segment('although the analysis was underpowered and equivalence was not established. Carcinoma',
    'although the analysis was underpowered and equivalence was not established. Carcinoma',
    'although the analysis was underpowered and equivalence was not established. Tumour–normal differences were essentially unchanged after adjustment for age and histological type (maximum change in Δβ 0.03). Background point estimates and near-gate pass/fail status were pipeline-dependent (pass/fail changed for 11 of 26 probes between noob and GMQN; both lead probes robust), so the fixed β thresholds are heuristic screening rules, not cross-platform constants. Carcinoma')

# ---------------- 4. Abstract Conclusions ----------------
replace_segment('they do not yet constitute a validated panel. Confirmation requires',
    'they do not yet constitute a validated panel. Confirmation requires',
    'they do not yet constitute a validated panel. Population matching was approximated rather than achieved, and external evidence derives from tumour tissue, not the intended sampling compartments. Confirmation requires')

# ---------------- 5. Introduction: intended use ----------------
append_to('benchmarking the final candidates against the contemporary literature',
    ' The intended use guiding every design choice is triage of premenopausal women with abnormal uterine bleeding—identifying patients who warrant further evaluation—with endometrial carcinoma and atypical hyperplasia as the target condition and histopathology as the reference standard. Because no public cohort matches this population, we treat population matching as a validation-cohort design requirement rather than a property claimed for the present data.')

# ---------------- 6. Five dimensions: heuristic thresholds + tail stats ----------------
replace_segment('The Δβ≥0.30 entry criterion matches effect sizes',
    'The Δβ≥0.30 entry criterion matches effect sizes',
    'The β>0.3 positivity cut-off and the P95 β≤0.15 background gate are heuristic screening rules for candidate prioritisation, not cross-platform biological constants; final thresholds will be locked in the validation cohort. Because P95 estimates in cohorts of 17–20 samples are governed by one or two extreme observations, we additionally report the maximum, P90 and the proportion of samples above 0.15 for every cohort–probe pair (Supplementary Data), and we treat near-gate pass/fail status as unstable. The Δβ≥0.30 entry criterion matches effect sizes')

# ---------------- 7. Statistical analysis: new methods ----------------
append_to('Candidate retention was scanned across background P95 gates of 0.05–0.30.',
    ' Confounder sensitivity was assessed with limma models on M-values adjusting for age and histological type (three-level tissue factor: adjacent normal, endometrioid tumour, non-endometrioid tumour), an endometrioid-only age-adjusted model, and a tumour-only model testing molecular-subtype association. Pipeline sensitivity was assessed by reprocessing the GSE73949 and GSE67116 IDATs under raw and noob pipelines and comparing background quantiles and tumour-signal rankings with the GMQN-normalised values used in triage.')

# ---------------- 8-10. Section renames ----------------
replace_para('External validation', 'Independent tissue-cohort corroboration')  # Methods heading (first exact match paragraph)
replace_segment('GSE178610 (80 endometrioid EC; 49 fresh-frozen, 31 FFPE).',
    'External validation was performed at the probe level in three independent EPIC cohorts, none used in discovery:',
    'Probe-level corroboration was performed in three independent EPIC cohorts, none used in discovery:')
append_to('GSE178610 (80 endometrioid EC; 49 fresh-frozen, 31 FFPE).',
    ' These cohorts provide tumour-tissue evidence only: they corroborate tumour association but do not validate the intended sampling scenario (cervical scraping or self-sampling), which requires purpose-collected specimens.')
replace_para('External validation in independent EPIC cohorts', 'Independent tissue-cohort corroboration in EPIC cohorts')

# ---------------- 11. Lead-probe evidentiary asymmetry ----------------
replace_segment('Two probes passed all measurable background dimensions',
    'Two probes passed all measurable background dimensions: cg24589459 (BOLL, assessable in all five cohorts; background P95 β 0.087–0.141 [bootstrap 95% CI at the limiting cervical dimension 0.100–0.158; a point-estimate pass whose confidence interval crosses the 0.15 gate—borderline], Δβ=0.533, premenopausal positivity 89.3% [Wilson 95% CI 72.8–96.3%]) and cg27577527 (ZSCAN12; P95 β 0.044–0.060 in the four measurable dimensions, Δβ=0.614, positivity 92.9%; absent from EPIC, so not assessable in the approximating cohort).',
    'The two lead probes differ in evidentiary status. cg24589459 (BOLL) is the only probe assessable in all five dimensions: background P95 β 0.100–0.141 across dimensions (bootstrap 95% CI at the limiting cervical dimension 0.100–0.158; a point-estimate pass whose confidence interval crosses the 0.15 gate—borderline), Δβ=0.533, premenopausal positivity 89.3% (Wilson 95% CI 72.8–96.3%). cg27577527 (ZSCAN12) passed the four dimensions measurable at probe level (P95 β 0.044–0.060; Δβ=0.614; positivity 92.9%) and showed the most consistently low background; the decisive intended-use-proxy dimension is not directly assessable because the probe is absent from EPIC arrays, and its compliance on that dimension is supported only indirectly, at region level, by the EPIC-compatible substitutes described below. We therefore do not claim five-dimension compliance for ZSCAN12.')

# ---------------- 12. New Results paragraph: cross-pipeline consistency ----------------
insert_before('Menopause-stratified, matched-pair and subtype sensitivity analyses in TCGA.',
    'Cross-pipeline consistency. Because discovery and triage combined Xena-precomputed TCGA β values, noob-processed IDATs and GMQN-normalised matrices, we tested how preprocessing affects the fixed-threshold gates. Reprocessing the GSE73949 IDATs (n=17) under raw and noob pipelines and comparing with GMQN values from the EWAS Data Hub shifted background P95 point estimates by a median of 0.022 (maximum 0.045) and changed pass/fail status at the 0.15 gate for 11 of 26 probes between noob and GMQN—all of them near-gate probes (P95 0.125–0.178 across pipelines). Both lead probes were robust: cg24589459 and cg27577527 passed under every pipeline evaluated (noob P95 0.110 and 0.060; GMQN 0.131 and 0.085). The three ZSCAN12 substitute probes passed under noob but exceeded the gate under GMQN (P95 0.160–0.167), so their healthy-endometrium compliance is pipeline-dependent and unresolved. In GSE67116, tumour–hyperplasia Δβ rankings were highly consistent between raw and noob (Spearman ρ=0.98). Quantile normalisation materially distorted probe-level β values in these small, compositionally extreme designs and was not used further. Candidate ranking is thus stable across mainstream pipelines, but absolute pass/fail status near the gate is pipeline-bound—reinforcing that the fixed thresholds are heuristic screening rules and that final cut-offs must be locked in the validation cohort under one predefined pipeline.')

# ---------------- 13. Adjusted-model results ----------------
replace_segment('The signal also persisted in the endometrioid-only subset',
    'The signal also persisted in the endometrioid-only subset (n=333; Δβ 0.44–0.67) and across all four molecular subtypes (Δβ 0.34–0.71; Supplementary Data).',
    'The signal also persisted in the endometrioid-only subset (n=333; Δβ 0.44–0.67) and across all four molecular subtypes (Δβ 0.34–0.71; Supplementary Data). Adjustment for measured confounders changed these estimates minimally: all 26 probes remained significant (FDR<0.05) in age-adjusted and age-plus-histology-adjusted limma models (adjusted versus unadjusted logFC Pearson r=0.97–0.99; maximum change in Δβ 0.03), and an endometrioid-only age-adjusted model was concordant (Δβ 0.46–0.67). Within tumours, methylation magnitude varied by molecular subtype (subtype F-test FDR<0.05 for all 26 probes), consistent with known subtype-associated global methylation differences; tumour purity, cell composition and batch remain unmeasured residual confounders.')

# ---------------- 14. Post hoc labelling ----------------
replace_segment('motivates our proposal to prospectively specify and lock',
    'This dissociation indicates complementary behaviour within the region—cg24589459 low-false-positive, cg07495363 high-sensitivity—and, being exploratory, motivates our proposal to prospectively specify and lock a two-probe region-level rule before wet-lab validation.',
    'This dissociation indicates complementary behaviour within the region—cg24589459 low-false-positive, cg07495363 high-sensitivity—and constitutes a post hoc, hypothesis-generating observation made after inspecting these external data; it motivates our proposal to pre-specify and lock a two-probe region-level rule in an independent threshold-setting cohort before any validation.')
replace_segment('As an illustrative, non-locked rule',
    'As an illustrative, non-locked rule, positivity of either BOLL-region probe yielded',
    'As a post hoc, hypothesis-generating rule derived from these same data (reported for transparency, not as a performance claim), positivity of either BOLL-region probe yielded')

# ---------------- 15. Substitute pipeline dependence ----------------
append_to('whereas cg25666433 was borderline (0.151, 95% CI 0.124–0.155, crossing the gate).',
    ' These healthy-endometrium estimates are pipeline-dependent: under GMQN-normalised values from the same 17 samples, all three substitutes exceed the gate (P95 0.160–0.167), so substitute compliance at this dimension is unresolved and will depend on the locked validation pipeline; we do not interpret sub-0.01 point differences (e.g., 0.149 versus 0.151) as meaningful.')

# ---------------- 16. Intended use in gradient paragraph ----------------
replace_segment('because the intended use is triage of abnormal uterine bleeding',
    'because the intended use is triage of abnormal uterine bleeding—where detecting atypical hyperplasia is clinically valuable—we interpret hyperplasia positivity as detection of the neoplastic continuum and avoid carcinoma-specificity claims.',
    'the intended use is triage of premenopausal abnormal uterine bleeding (target condition: EC and atypical hyperplasia; reference standard: histopathology), in which detecting atypical hyperplasia is clinically valuable; because atypia is unannotated, hyperplasia positivity is ambiguous—premalignant detection, non-atypical false positivity and cohort or technical bias cannot be separated—so we interpret it as a neoplastic-continuum signal and avoid carcinoma-specificity claims.')

# ---------------- 17. Limitations updates ----------------
replace_segment('bootstrap 95% CI for cg07495363 spans',
    '(bootstrap 95% CI for cg07495363 spans 0.099–0.203)',
    '(bootstrap 95% CI for cg07495363 spans 0.069–0.182; per-probe maximum, P90 and proportion above 0.15 in Supplementary Data—for several probes a minority of healthy donors exceed 0.15 even where the P95 passes, so the gate is a screening heuristic rather than a specificity estimate)')
replace_segment('Discovery and triage layers used different normalisation',
    'Discovery and triage layers used different normalisation (noob versus GMQN); final thresholds will be set in a prospectively collected cohort processed with one predefined pipeline.',
    'Discovery and triage layers used different normalisation (noob versus GMQN); background point estimates shifted by up to 0.045 between pipelines and near-gate pass/fail status changed for 11 of 26 probes (both lead probes robust), so thresholds are heuristic and pipeline-bound; final thresholds will be set in a prospectively collected cohort processed with one predefined pipeline.')
replace_segment('External validation used tumour tissue',
    'External validation used tumour tissue, not the intended sampling compartments',
    'Tissue-cohort corroboration used tumour tissue, not the intended sampling compartments')

# ---------------- 18. Conclusion ----------------
replace_segment('multi-cohort discovery with population- and compartment-aware',
    'population- and compartment-aware background triage',
    'population- and compartment-informed background triage')
replace_segment('remained borderline at the cervical-background dimension.',
    'remained borderline at the cervical-background dimension.',
    'remained borderline at the cervical-background dimension, and support for ZSCAN12 at the intended-use-proxy dimension is region-level and pipeline-dependent.')

# ---------------- 19. Compress cohort-dependence paragraph ----------------
replace_para('Background P95 β values in the benign endometrial cohort approximating the intended-use population (n=347) were higher than in healthy premenopausal donors',
    'Background P95 β values in the benign endometrial cohort approximating the intended-use population (n=347) were higher than in healthy premenopausal donors (n=17). Because this contrast spans platform (450K versus EPIC) and normalisation (noob versus GMQN), we decomposed it in an EPIC/GMQN-aligned comparison (same platform and reported normalisation framework, different cohorts and centres): the benign cohort remained consistently elevated relative to cancer-free endometrial controls (n=13; e.g., P95 β 0.100 versus 0.059 for cg24589459, 0.372 versus 0.252 for CDO1), approximately 1.2- to 1.6-fold, so the raw two- to three-fold contrast is inflated by processing while the qualitative difference persists (cohort, centre and residual processing effects cannot be fully separated; reference n=13). Endometriosis samples (n=637) showed per-probe P95 values comparable to—and consistently slightly below—those of the 347 controls (per-group values in Supplementary Data), so the control-based gate is the conservative choice for this dimension. A union rule requiring P95 β≤0.15 in every assessable cohort left two probes (below); for cg27577527 (ZSCAN12), the intended-use-proxy dimension was not directly assessable because the probe is absent from EPIC arrays.')

# ---------------- 20. Compress kit-audit paragraph ----------------
replace_para('Auditing the NMPA-kit and WID-qEC components through the same funnel',
    'Auditing the NMPA-kit and WID-qEC components through the same funnel (Table 6) showed that CDO1 combined the lowest cervical (0.055) and blood (0.080) backgrounds among kit genes, consistent with its dual use in cervical (CISENDO) and intrauterine formats. AJAP1 and GALR1 showed elevated blood-cell backgrounds (0.27–0.44), and GYPC showed the heaviest background across all dimensions (cervical 0.37–0.68, blood up to 0.84) and was eliminated by our first gate (per-probe values in Supplementary Data). CELF4, a principal marker of the Taiwanese series [8,25], has no 450K probe and was invisible to our pipeline. Several literature-supported genes were near-misses (GRIA4, SOX11, HS3ST2: Δβ≈0.5, acceptable background, sub-90% positivity) and are recommended as wet-lab controls.')

# ---------------- 21. Compress discussion re-identification paragraph ----------------
replace_para('Re-identification of established markers supports the face validity',
    'Re-identification of established markers supports the face validity of the prioritisation framework: without any prior gene list, the funnel re-identified the core of WID-qEC (ZSCAN12 [5,6]), the anchor of the only approved national kit (CDO1 [8,10]) and a Dutch-panel component (GHSR [24]). Other components did not meet our prespecified probe-level background criteria (GYPC, AJAP1, GALR1; Results), and CELF4 is structurally invisible to 450K-based discovery. These probe-level findings do not evaluate the exact CpGs, amplicons or classification algorithms used in the commercial or published assays: WID-qEC absorbs GYPC\'s probe-level background through region-level ΣPMR quantification, and the approved kit deploys AJAP1/GALR1 via an intrauterine brush with cohort-specific cut-offs. Marker choice and sampling design are inseparable—a point we formalised as compartment-suitability scoring.')

# ---------------- 22. Compress BOLL novelty paragraph ----------------
replace_segment('plausibly reflects the candidate-gene path-dependence',
    'plausibly reflects the candidate-gene path-dependence of the field, the dilution of probe-level signals in gene-level rankings, and the fact that our ranking rewards low background—a dimension most studies do not weight. Absence from the literature may also reflect unpublished negative results; whether this novelty',
    'plausibly reflects the candidate-gene path-dependence of the field, dilution of probe-level signals in gene-level rankings, unpublished negative results, and the fact that our ranking rewards low background—a dimension most studies do not weight. Whether this novelty')

# ---------------- 23. Table 1 role cells + Table 8 caption + Figure 4 legend ----------------
n_cells = 0
for t in doc.tables:
    for row in t.rows:
        for c in row.cells:
            if 'External validation' in c.text:
                for p in c.paragraphs:
                    if 'External validation' in p.text:
                        set_text(p, p.text.replace('External validation', 'Tissue-cohort corroboration'))
                        n_cells += 1
LOG.append(f'table cells renamed: {n_cells}')
replace_segment('Table 8. Probe-level external validation',
    'Table 8. Probe-level external validation in independent EPIC cohorts.',
    'Table 8. Probe-level corroboration in independent EPIC tumour-tissue cohorts.')
replace_segment('independent probe-level validation is reported in the External validation section',
    'independent probe-level validation is reported in the External validation section',
    'independent probe-level corroboration is reported in the Independent tissue-cohort corroboration section')

doc.save(DST)
print('\n'.join(LOG))
misses = [l for l in LOG if l.startswith('MISS')]
print(f'\n=== {len(misses)} MISSES ===')
for m in misses: print(m)
