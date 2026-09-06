# Data directory

No raw or processed methylation data are redistributed in this repository.
All datasets are public; see `accession_manifest.tsv` for accession numbers,
composition, study role and the exact access route used for each dataset.

## Access routes (summary)

| Source | Datasets | Route |
|---|---|---|
| UCSC Xena | TCGA-UCEC | Precomputed 450K beta values (`TCGA.UCEC.sampleMap_HumanMethylation450`); clinical matrix from Xena; molecular subtypes from PanCanAtlas via cBioPortal |
| EBI BioStudies mirror | GSE67116, GSE73949 | IDAT files (`E-GEOD-67116`, `E-GEOD-73949`), processed with minfi/noob — used because the local network blocks NCBI/GEO directly |
| EWAS Data Hub | GSE223817, GSE90060, GSE46306, GSE35069 | Normalised beta matrices (GMQN for GSE223817) |
| GEO | GSE136791, GSE178610, GSE155760, GSE93589 | Published beta matrices |

## Sample inclusion / exclusion

Cohort composition and the role of each dataset are listed in
`accession_manifest.tsv` and in Table 1 of the manuscript. Probe-level
filtering for IDAT-processed cohorts: samples retained when <10% of probes had
detection P>0.01; probes excluded for detection failure, sex chromosomes or
common SNPs. For GSE223817 the union-rule gate was computed on the 347 control
samples; the 637 endometriosis samples were evaluated descriptively.

## Note on processing heterogeneity

Discovery-layer beta values combined Xena-precomputed TCGA data with
noob-normalised GSE67116 data, whereas background triage used GMQN-normalised
values. Cross-pipeline sensitivity of every gate is quantified in the
manuscript (scripts 54b-56; Results/v911_pipeline_*.csv in the supplementary
package). All fixed beta thresholds are heuristic screening rules.
