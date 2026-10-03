# EC-methylation-prioritization

Analysis code accompanying the manuscript:

> **Population- and compartment-informed prioritization of DNA methylation markers for endometrial cancer: an in silico evaluation of BOLL and ZSCAN12** (manuscript v9.29)

The study re-analyses public 450K/EPIC methylation data to prioritise DNA methylation
regions for prospective validation in premenopausal abnormal uterine bleeding (AUB)
cohorts, triaging candidates against five background dimensions (healthy premenopausal
endometrium, a symptomatic/benign endometrial proxy, menstrual cycle, cervical scrapes
and blood-reference fractions). It is a computational prioritisation study: the outputs are
candidate regions (BOLL, ZSCAN12) and a 23-probe experimental pool (plus a 26-probe
extended sensitivity set), **not** a validated diagnostic panel.

## Repository contents

```
README.md                       this file
LICENSE                         MIT
CITATION.cff                    citation metadata
requirements.txt                Python dependencies
sessionInfo.txt                 R session info (versions)
config/analysis_parameters.yml  all thresholds, gates and seeds
data/accession_manifest.tsv     datasets, composition, role, access route
data/README.md                  how to obtain the raw data (nothing is redistributed)
scripts/                        analysis pipeline, numbered in execution order
results/candidate_probes.tsv    23-probe pool + 3 EPIC substitutes (extended set of 26)
results/candidate_regions.tsv   BOLL and ZSCAN12 lead regions (hg38 + hg19 coordinates)
results/v917_*.csv              v9.17 verification results (see below)
results/v917_figure3_compartment.png, v917_figure4_crosscohort.png
results/v921_figure3_compartment.png  Figure 3 with corrected Table 6 legend (v9.21+)
results/v922_selection_funnel_bootstrap.csv  per-probe bootstrap, v9.22 design (superseded)
results/v923_selection_funnel_bootstrap.csv  per-probe conditional selection-stability (v9.23, S2f)
results/v928_gate_proportion_ci_*.csv        background-gate proportions + Wilson CIs (v9.27, S2g/S2h inputs)
results/v930_expanded_background_476.csv      five additional benign/background cohorts (v9.27, S2h)
results/v931_pipeline_sensitivity_gse73949.csv  noob vs GMQN pipeline sensitivity (v9.28, S2i)
results/v932_boll_rule_*.csv                  BOLL five-rule comparison, B=2000 + benign LOO (v9.28, S2j/S2k)
results/v933_threshold_grid_summary.csv       discovery-threshold grid, 27 cells (v9.29, S2l/S2m; superseded by v2)
results/v936_threshold_grid_summary_v2.csv    corrected grid v2: cache-complete gates + official score (v9.30, S2l)
results/v936_new_entrants_background_v2.csv   new-entrant gates incl. endometriosis background (v9.30, S2q)
results/v936_reconciliation.csv               cross-analysis reconciliation of counts/ranks (v9.30, S2r/Table S14)
results/v934_sampling_paths_*.csv             application-specific gating paths (v9.29, S2n/S2o)
results/v935_retention_ledger.csv             per-probe retention ledger, 23-pool (v9.29, S2p/Table S10)
results/v929_figure1_funnel.png               Figure 1 funnel, corrected branch counts (v9.29)
results/stability_cache.npz     cached analysis matrices for scripts 81-84
```

## Environment

- R 4.6.0 (minfi 1.58.0, limma 3.68.4, bumphunter 1.54.0; see `sessionInfo.txt`)
- Python 3.13 (pandas, numpy, matplotlib, python-docx; see `requirements.txt`)

## Data access

All datasets are public; no raw data are redistributed here. Accession numbers,
composition, study role and download route for each dataset are in
`data/accession_manifest.tsv`. TCGA-UCEC beta values were taken precomputed from
UCSC Xena; GSE67116 and GSE73949 IDATs were obtained via the EBI BioStudies mirror
and processed with minfi (noob); remaining cohorts were obtained as normalised beta
matrices from the EWAS Data Hub or GEO. Sample inclusion/exclusion rules are
described in the manuscript (Table 1 and Methods).

## Running order

Scripts are numbered in execution order (`scripts/01_...` to `scripts/60_...`).
They assume the working directory is the project root with `data/` and `results/`
subfolders. Key stages:

| Stage | Scripts |
|---|---|
| Preprocessing (minfi/noob from IDAT) | 01-03 |
| DMP/DMR discovery and candidate ranking | 04-12 |
| Background triage (five dimensions) | 05, 09-11 |
| Panel selection and literature/kit audit | 06-08, 12 |
| Sensitivity analyses (pairs, subtypes, menopause, threshold scan) | 18-32, 40 |
| bumphunter cross-check | 41 |
| Tissue-cohort corroboration | 20-22 |
| Confounder-adjusted models (v9.11) | 52-53 |
| Cross-pipeline consistency (v9.11) | 54b |
| Background tail statistics + noob recompute (v9.11) | 55-56 |
| OR-rule joint positivity in benign cohorts (v9.17) | 70, 70c |
| Pairing-fingerprint diagnostics, age-concordant drift (v9.17) | 70b, 70c |
| Table 8 control-column recomputation (v9.17) | 70c |
| Figure generation (v9.17) | 70d, 70e |
| Conditional selection-stability analysis (v9.22/v9.23) | 81-84 |
| Background-gate proportions + Wilson CI, cervix denominator fix (v9.27) | 86, 86b |
| Expanded benign-background cohorts via EWAS Data Hub (v9.27) | 87 |
| Pipeline sensitivity noob vs GMQN, GSE73949 (v9.28) | 85 |
| BOLL five-rule comparison, bootstrap + benign LOO (v9.28) | 88 |
| Discovery-threshold sensitivity grid (v9.29) | 89 |
| Corrected grid v2: cache-complete gates, official score, reconciliation (v9.30) | 89b |
| Application-specific gating paths (v9.29) | 90 |
| Per-probe retention ledger (v9.29) | 91 |
| Figure 1 funnel redraw (v9.29) | 92 |

## v9.29 rule-closure revisions (pre-submission)

Four pre-submission analyses close the rule loop requested during internal audit:

- **Discovery-threshold sensitivity** (`89_c4_threshold_sensitivity.py`, v9.29, S2l/S2m):
  a 27-cell grid over the three discovery cut-offs (delta-beta 0.30/0.40/0.50, overall
  positivity 0.80/0.85/0.90, premenopausal positivity 0.75/0.80/0.85), re-derived
  genome-wide (485,512 probes). ZSCAN12 ranks first in 27/27 cells; BOLL anchor passes
  triage in all cells (rank 3-11). No new low-background candidates enter under the
  original >=0.90 positivity requirement. CDO1 passes all discovery criteria (it is
  eliminated only at the background layer); AJAP1 is excluded by the overall-positivity
  cut-off (84.2%); GALR1/GYPC enter in no cell.
- **Application-specific gating paths** (`90_c5_sampling_paths.py`, v9.29, S2n/S2o):
  strict four-gate (2/476), endometrial-brush path (11/476) and cervix/self-sampling
  path (43/476) compared on both the 476-probe universe and the 23-probe pool. CDO1
  passes the cervix path (0.055/0.080), consistent with its approved cervical kit role.
- **Per-probe retention ledger** (`91_retention_ledger.py`, v9.29, Table S10/S2p):
  all 23 pool probes annotated with entry path (strict / literature-reference /
  region-coverage), per-gate results, evidence completeness and retention rationale;
  8/23 are region-coverage probes retained despite failing triage, now disclosed
  consistently in Figure 1, Table S10 and the main text.
- **Figure 1 corrections** (`92_figure1_v6.py`): triage elimination count corrected to
  461/476 (was 457, inconsistent with Table 2); region-coverage branch corrected to
  15 probes (was 16, which summed to 24 across the four branches vs the 23-probe pool).

## v9.30 reconciliation corrections (pre-submission)

A pre-submission numeric audit found that the v1 threshold grid (89) read gate values from
the elite audit table, whose unevaluated cells were pipeline-stage completeness artefacts:
cg15285250 (NXPH1) lacked cervix/blood values in that table but had measured values
(0.186/0.276) in the stability cache, so v1 wrongly let it pass at the original cell and
ranked BOLL 3/3 with a simplified score. The corrected grid (`89b_c4_threshold_sensitivity_v2.py`)
uses cache-complete four-gate values for the 476, de novo four-cohort evaluation for new
entrants (including the endometriosis background used by the score), and the official
score_v6 formula. The original cell now reproduces the main analysis exactly: 479 candidates
(476 + 3 annotation-filtered probes), 2 triage survivors, ZSCAN12 rank 1/2, BOLL rank 2/2;
BOLL ranks 2-6 across the 27-cell grid; ZSCAN12 ranks first in 27/27 cells. The
cross-analysis reconciliation table (main funnel 476->65->15->2; C5 paths; C4 v2; pipeline
sensitivity 15->10) is `results/v936_reconciliation.csv` (manuscript Table S14, Note S2.17).

Scripts 85-95b are analysis, revision and verification utilities for v9.24-v9.29.

## v9.23 conditional selection-stability analysis

Scripts 81-84 recompute the complete per-probe single-probe funnel (discovery criteria,
four background gates, effective three-component score, re-ranking) inside 2,000
stratified bootstrap replicates, conditional on the fixed 476-probe audit universe.
Because the 28 premenopausal tumours are a subset of the 431 tumours, the two disjoint
strata (28 premenopausal / 403 other tumours) are resampled separately and pooled per
replicate (nested design; v9.23). Genome-wide significance is not re-derived: this is a
conditional selection-stability analysis, not a re-run of the genome-wide screen.

- `81_stability_extract.py`: builds the analysis matrices (TCGA beta for the 476 probes;
  the four background cohorts) and **validates exact reproduction** of every audit point
  estimate (max |diff| <= 0.0005) before any resampling. Note the data sources per gate:
  gate 1 = minfi-noob `GSE73949_betas.rds`; gate 2 = GSE223817 controls (n=347);
  gate 3 = GSE46306 HPV-negative cervix (n=20); gate 4 = GSE35069 all 60 samples,
  per-fraction P95 then max across the 10 blood-reference fractions.
- `82_stability_bootstrap.py`: the main nested-resampling loop with per-condition failure
  decomposition -> `results/v923_selection_funnel_bootstrap.csv`.
- `83_stability_pre_influence.py`: premenopausal-case influence (leave-one-out jackknife
  over the 28 premenopausal tumours).
- `84_stability_pre_removal_v2.py`: nested-design variant without the premenopausal
  criterion/score term.

Headline results (Supplementary Note S2.9, v9.23): ZSCAN12 anchor Tier A in 53.1% of
replicates (rank 1 in all retained replicates; attrition = overall positivity at the
>=0.90 cut-off, 45.1% failures); BOLL anchor Tier A in 48.3% (rank <=5 in 99.8% when
retained; decomposition: overall positivity 23.1% failures, premenopausal positivity
16.1%, cervical gate 28.4% failures, other gates 100% stable); leading reserve candidate
cg06014401 (EOMES, 38.3%). Removing the premenopausal criterion: BOLL 57.1%,
ZSCAN12 54.9%.

## v9.17 verification results

- `results/v917_boll_or_rule_background.csv`: positivity of the post hoc two-probe
  BOLL OR rule (beta > 0.3) across benign background cohorts — 159/347 (45.8%) in the
  benign surgical cohort and 214/637 (33.6%) in its endometriosis subset (driven by
  cg07495363; anchor alone 2/347 and 0/637); 0 in healthy endometrium (17), HPV-negative
  cervix (20), purified blood cells (42) and the 13 endometrial controls.
- `results/v917_table8_ctrl13_recompute.csv` / `v917_table8_ctrl13_all21.csv`: recomputed
  P95 beta for the 13 GSE155760 endometrial controls (manuscript Table 8 control column).
- `results/v917_cycle_drift_concordant.csv`: cycle drift on all 17 inferred pairs vs the
  15 age-concordant pairs (GSE90060).

Scripts 57-60, 58b-58g, 66-68 and 71 are manuscript-revision and verification utilities
(docx editing and checklists), included for provenance; they are not part of the
analysis pipeline.

## Key parameters

All gates and seeds are in `config/analysis_parameters.yml`: positivity cut-off
beta > 0.3, background gate P95 beta <= 0.15 (union rule), entry effect size
delta-beta >= 0.30, bumphunter B = 500 permutations (chromosome-wise, seed 42),
bootstrap 10,000 draws (seeds 42 / 20260906), quantile type 7. These thresholds
are heuristic screening rules, not cross-platform constants.

## Licence

Code: MIT (see `LICENSE`). No patient-level data are included; all datasets remain
subject to their original terms of use.
