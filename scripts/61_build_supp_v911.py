"""Build EC_methylation_supplementary_v911.zip = v910 contents + v9.11 review-round additions."""
import zipfile, shutil, os

SRC = 'EC_methylation_supplementary_v910.zip'
DST = 'EC_methylation_supplementary_v911.zip'

NEW_RESULTS = [
    ('results/v911_tcga_adjusted_models.csv',       'Results/v911_tcga_adjusted_models.csv'),
    ('results/v911_tcga_eec_only_adjusted.csv',     'Results/v911_tcga_eec_only_adjusted.csv'),
    ('results/v911_tcga_tumour_only_subtype.csv',   'Results/v911_tcga_tumour_only_subtype.csv'),
    ('results/v911_background_quantiles_final.csv', 'Results/v911_background_quantiles_final.csv'),
    ('results/v911_gse73949_noob_stats.csv',        'Results/v911_gse73949_noob_stats.csv'),
    ('results/v911_pipeline_p95_gse73949.csv',      'Results/v911_pipeline_p95_gse73949.csv'),
    ('results/v911_pipeline_discordant_gse73949.csv','Results/v911_pipeline_discordant_gse73949.csv'),
    ('results/v911_pipeline_tumour_gse67116.csv',   'Results/v911_pipeline_tumour_gse67116.csv'),
]

NEW_CODE = [
    ('scripts/52_v911_tcga_prep.py',               'Code/52_v911_tcga_prep.py'),
    ('scripts/53_v911_tcga_adjusted.R',            'Code/53_v911_tcga_adjusted.R'),
    ('scripts/54b_v911_gse67116_pipelines.R',      'Code/54b_v911_gse67116_pipelines.R'),
    ('scripts/55_v911_background_quantiles.py',    'Code/55_v911_background_quantiles.py'),
    ('scripts/56_v911_gse73949_noob_recompute.R',  'Code/56_v911_gse73949_noob_recompute.R'),
    ('scripts/57_v911_revisions.py',               'Code/57_v911_revisions.py'),
    ('scripts/58_v911_checklist.py',               'Code/58_v911_checklist.py'),
    ('scripts/59_v911_compression.py',             'Code/59_v911_compression.py'),
    ('scripts/60_v911_compression2.py',            'Code/60_v911_compression2.py'),
]

INDEX_ADD = """## Added in v9.11 (review round: repositioning + confounder/pipeline sensitivity)
- Code/52-53: TCGA confounder-adjusted sensitivity (limma models adjusting for age and
  histological type; endometrioid-only age-adjusted model; tumour-only molecular-subtype model)
- Code/54b: cross-pipeline consistency (GSE73949 raw/noob vs GMQN background P95;
  GSE67116 raw vs noob tumour-hyperplasia delta-beta ranking)
- Code/55-56: background tail statistics (max, P90, P95, proportion above 0.15 for every
  cohort-probe pair) and GSE73949 noob recomputation from IDAT with 10,000-draw bootstrap CIs
- Code/57-60: manuscript v9.10 -> v9.11 edits, verification checklist, compression passes
- Results/v911_tcga_adjusted_models.csv: adjusted vs unadjusted effects for all 26 probes
- Results/v911_tcga_eec_only_adjusted.csv, v911_tcga_tumour_only_subtype.csv
- Results/v911_background_quantiles_final.csv: max/P90/P95/proportion>0.15 per cohort-probe
  (GSE73949 values from noob-processed IDAT, not GMQN)
- Results/v911_gse73949_noob_stats.csv: noob P95 with bootstrap 95% CIs (n=17)
- Results/v911_pipeline_*.csv: per-pipeline P95, discordant probes (11/26), GSE67116 rankings

"""

zin = zipfile.ZipFile(SRC)
missing = [s for s, _ in NEW_RESULTS + NEW_CODE if not os.path.exists(s)]
if missing:
    print('MISSING SOURCE FILES:'); [print(' ', m) for m in missing]; raise SystemExit(1)

with zipfile.ZipFile(DST, 'w', zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == 'INDEX.md':
            text = data.decode('utf-8')
            text = text.replace('# Supplementary package v9.6', '# Supplementary package v9.11', 1)
            # insert the v9.11 section right after the header block (before first "## Added")
            idx = text.find('## Added in v9.6')
            text = text[:idx] + INDEX_ADD + text[idx:]
            data = text.encode('utf-8')
        zout.writestr(item, data)
    for src, arc in NEW_RESULTS + NEW_CODE:
        zout.write(src, arc)

zin.close()
z = zipfile.ZipFile(DST)
print(DST, '->', len(z.namelist()), 'files (v910 had 70)')
bad = z.testzip()
print('zip integrity:', 'OK' if bad is None else f'BAD: {bad}')
