"""Build EC_methylation_supplementary_v917.zip = v911 contents + v9.12-v9.17 additions."""
import os
import zipfile

SRC = 'EC_methylation_supplementary_v911.zip'
DST = 'EC_methylation_supplementary_v917.zip'

NEW_RESULTS = [
    ('results/v917_boll_or_rule_background.csv',   'Results/v917_boll_or_rule_background.csv'),
    ('results/v917_table8_ctrl13_recompute.csv',   'Results/v917_table8_ctrl13_recompute.csv'),
    ('results/v917_table8_ctrl13_all21.csv',       'Results/v917_table8_ctrl13_all21.csv'),
    ('results/v917_cycle_drift_concordant.csv',    'Results/v917_cycle_drift_concordant.csv'),
    ('results/v917_figure3_compartment.png',       'Figures/v917_figure3_compartment.png'),
    ('results/v917_figure4_crosscohort.png',       'Figures/v917_figure4_crosscohort.png'),
]

NEW_CODE = [
    ('scripts/58b_v912_checklist.py',  'Code/58b_v912_checklist.py'),
    ('scripts/58c_v913_checklist.py',  'Code/58c_v913_checklist.py'),
    ('scripts/66_v914_revisions.py',   'Code/66_v914_revisions.py'),
    ('scripts/58d_v914_checklist.py',  'Code/58d_v914_checklist.py'),
    ('scripts/67_v915_revisions.py',   'Code/67_v915_revisions.py'),
    ('scripts/58e_v915_checklist.py',  'Code/58e_v915_checklist.py'),
    ('scripts/68_v916_revisions.py',   'Code/68_v916_revisions.py'),
    ('scripts/58f_v916_checklist.py',  'Code/58f_v916_checklist.py'),
    ('scripts/70_or_rule_background.py',      'Code/70_or_rule_background.py'),
    ('scripts/70b_cycle_pairing_check.py',    'Code/70b_cycle_pairing_check.py'),
    ('scripts/70c_v917_verification.py',      'Code/70c_v917_verification.py'),
    ('scripts/70d_figure4_crosscohort.py',    'Code/70d_figure4_crosscohort.py'),
    ('scripts/70e_figure3_harmonised.py',     'Code/70e_figure3_harmonised.py'),
    ('scripts/71_v917_revisions.py',          'Code/71_v917_revisions.py'),
    ('scripts/58g_v917_checklist.py',         'Code/58g_v917_checklist.py'),
]

INDEX_ADD = """## Added in v9.17 (review round: numeric verification + repositioning)
- Code/70, 70b, 70c: verification analyses — BOLL two-probe OR-rule joint positivity across all
  benign background cohorts; GSE90060 pairing-fingerprint diagnostics (within-pair vs best
  cross-donor correlation, age concordance); Table 8 control-column recomputation from per-sample files
- Code/70d, 70e: regenerated Figure 4 (cross-cohort performance of the two lead regions) and
  Figure 3 (colour categories harmonised with Table 7)
- Code/71 + 58g: manuscript v9.16 -> v9.17 edits and verification checklist
- Results/v917_boll_or_rule_background.csv: OR-rule positivity per background cohort
  (159/347 = 45.8% in the benign surgical cohort; 0 in all others)
- Results/v917_table8_ctrl13_recompute.csv / _all21.csv: recomputed P95 for the 13 GSE155760
  endometrial controls (corrects the v9.16 Table 8 control column)
- Results/v917_cycle_drift_concordant.csv: cycle drift on all 17 pairs vs 15 age-concordant pairs
- Figures/v917_figure3_compartment.png, v917_figure4_crosscohort.png
- Code/66-68 + 58b-58f: earlier text-revision scripts and checklists for v9.12-v9.16
  (abstract compression; target-condition wording; GSE93589/RNA-seq methods; panel terminology)

"""

zin = zipfile.ZipFile(SRC)
missing = [s for s, _ in NEW_RESULTS + NEW_CODE if not os.path.exists(s)]
if missing:
    print('MISSING SOURCE FILES:')
    for m in missing:
        print(' ', m)
    raise SystemExit(1)

with zipfile.ZipFile(DST, 'w', zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == 'INDEX.md':
            text = data.decode('utf-8')
            text = text.replace('# Supplementary package v9.11', '# Supplementary package v9.17', 1)
            idx = text.find('## Added in v9.11')
            text = text[:idx] + INDEX_ADD + text[idx:]
            data = text.encode('utf-8')
        zout.writestr(item, data)
    for src, arc in NEW_RESULTS + NEW_CODE:
        zout.write(src, arc)

zin.close()
z = zipfile.ZipFile(DST)
print(DST, '->', len(z.namelist()), 'files')
bad = z.testzip()
print('zip integrity:', 'OK' if bad is None else f'BAD: {bad}')
