# EC-methylation-prioritization

Analysis code accompanying the manuscript:

> **Population- and compartment-informed prioritization of DNA methylation markers for endometrial cancer: an in silico evaluation of BOLL and ZSCAN12** (manuscript v9.17)

The study re-analyses public 450K/EPIC methylation data to prioritise DNA methylation
regions for prospective validation in premenopausal abnormal uterine bleeding (AUB)
cohorts, triaging candidates against five background dimensions (healthy premenopausal
endometrium, a symptomatic/benign endometrial proxy, menstrual cycle, cervical scrapes
and purified blood cells). It is a computational prioritisation study: the outputs are
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
