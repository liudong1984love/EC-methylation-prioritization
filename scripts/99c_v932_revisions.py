# fix_ms_v932.py — Manuscript v9.31 -> v9.32
# 1) Methods: donor-level blood resampling (v3 design)
# 2) Methods: correct cross-cohort consistency statement (code-verified: GSE67116 direction filter applied)
# 3) Results: stability numbers updated to donor-level v3 rerun
# 4) "purified blood cells" -> "blood-reference fractions"
# 5) Abstract: OR rule = two BOLL-region probes
# 6) Intro: soften "no public cohort matches"
# 7) BOLL expression wording softened
# 8) CDO1 dual-use claims tightened (x2)
# 9) Limitations 15->10 simplified
# 10) Figure 4 caption + image replacement
import zipfile, os, shutil

SRC = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Manuscript_EC_methylation_panel_v9.31.docx'
DST = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Manuscript_EC_methylation_panel_v9.32.docx'
FIG4 = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\results\v932_figure4_crosscohort.png'

zin = zipfile.ZipFile(SRC)
xml = zin.read('word/document.xml').decode('utf-8')

GE = '\u2265'   # ≥
DB = '\u0394\u03b2'  # Δβ
MD = '\u2014'   # em dash

REPL = [
 # --- 1) Methods: donor-level blood resampling ---
 ('(GSE35069 within cell-type fractions)',
  '(GSE35069 at donor level: the six verified donors resampled with replacement, each drawn donor '
  'contributing all ten blood-reference fractions' + MD + 'per-fraction P95 then maximum' + MD +
  'preserving cross-fraction donor correlation; donor identities verified against ArrayExpress sample '
  'metadata, Supplementary Note S2.18)'),
 # --- 2) Methods: cross-cohort consistency correction ---
 ('imposed no cross-cohort consistency requirement (Supplementary Data S2b)',
  'additionally required cross-cohort direction consistency (concordant '
  'tumour' + MD + 'hyperplasia ' + DB + '>0 in GSE67116) and presence in the GSE73949 noob matrix '
  '(Supplementary Data S2b)'),
 # --- 3) Results: v3 stability numbers ---
 ('re-entered Tier A in 53.1% and 48.3% of replicates',
  're-entered Tier A in 52.8% and 48.7% of replicates'),
 ('criterion (45.1% of replicates)',
  'criterion (45.9% of replicates)'),
 ('overall positivity below ' + GE + '0.90 in 23.1% of replicates (observed 0.909), premenopausal '
  'positivity below ' + GE + '0.85 in 16.1%, and the borderline cervical gate failing in 28.4% '
  '(71.6% compliance)',
  'overall positivity below ' + GE + '0.90 in 21.7% of replicates (observed 0.909), premenopausal '
  'positivity below ' + GE + '0.85 in 16.8%, and the borderline cervical gate failing in 26.6% '
  '(73.5% compliance)'),
 # --- Results stability preamble: note donor-level design ---
 ('(2,000 nested bootstrap replicates; premenopausal and other tumours resampled as separate strata '
  'and pooled, so that both positivity metrics derive from the same draw; Supplementary Note S2.9)',
  '(2,000 nested bootstrap replicates; premenopausal and other tumours resampled as separate strata '
  'and pooled, so that both positivity metrics derive from the same draw, and blood-reference '
  'fractions resampled at the donor level; Supplementary Note S2.9)'),
 # --- 4) purified blood cells ---
 ('cervical scrapes and purified blood cells)',
  'cervical scrapes and blood-reference fractions)'),
 # --- 5) Abstract OR rule ---
 ('and a simple two-probe OR combination was frequently positive',
  'and an OR rule combining the two BOLL-region probes was frequently positive'),
 # --- 6) Intro softening ---
 ('Because no public cohort matches this population',
  'Because no eligible cohort identified in our search matched this population'),
 # --- 7) BOLL expression wording ---
 ('BOLL expression was essentially absent from endometrial tissue irrespective of methylation state',
  'BOLL expression was low in most endometrial tissue samples irrespective of methylation state'),
 # --- 8) CDO1 dual-use tightening (Discussion) ---
 ('consistent with its dual use in cervical (CISENDO) and intrauterine formats',
  'supporting its cervical-format use (CISENDO); its intrauterine application rests on prior '
  'clinical studies, as its benign endometrial background in our cohorts did not meet the gate'),
 # --- 8b) CDO1 dual-use tightening (Results) ---
 ('CDO1 showed a dual-compartment profile (cervical 0.055, blood 0.080) consistent with its use in '
  'both cervical and intrauterine formats.',
  'CDO1 showed low cervical (0.055) and blood (0.080) backgrounds supporting cervical-pathway use; '
  'its intrauterine application derives from prior clinical studies rather than from the '
  'benign-endometrial background assessed here.'),
 # --- 9) Limitations 15->10 simplified ---
 ('the stage-2 triage-passing set' + MD + 'an auxiliary stage upstream of the strict all-four-gate '
  'Tier-A designation, not the Tier-A survivors' + MD + 'from 15 to 10',
  'the auxiliary stage-2 triage-passing set (upstream of the strict Tier-A designation, not the '
  'Tier-A survivors) from 15 to 10'),
 # --- 10) Figure 4 caption ---
 ('Figure 4. Cross-cohort performance of the BOLL and ZSCAN12 regions.',
  'Figure 4. Cross-cohort methylation backgrounds and tumour-tissue positivity of the BOLL and '
  'ZSCAN12 regions.'),
]

for old, new in REPL:
    n = xml.count(old)
    print(n, '|', old[:70])
    assert n == 1, f'anchor count {n} != 1: {old[:60]}'
    xml = xml.replace(old, new, 1)

# residual scan
for bad in ['within cell-type fractions', 'purified blood', 'no cross-cohort consistency requirement',
            '53.1%', '48.3%', '45.1%', '23.1% of replicates', 'two-probe OR combination',
            'essentially absent', 'dual-compartment profile', 'no public cohort matches']:
    c = xml.count(bad)
    if c: print('RESIDUAL', bad, c)

fig4_bytes = open(FIG4, 'rb').read()
with zipfile.ZipFile(DST + '.tmp', 'w', zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == 'word/document.xml':
            data = xml.encode('utf-8')
        elif item.filename == 'word/media/image4.png':
            data = fig4_bytes
        zout.writestr(item, data)
zin.close()
os.replace(DST + '.tmp', DST)

import xml.etree.ElementTree as ET
z = zipfile.ZipFile(DST)
assert z.testzip() is None
ET.fromstring(z.read('word/document.xml'))
print('image4 size:', z.getinfo('word/media/image4.png').file_size, '(expect', len(fig4_bytes), ')')
print('OK ->', DST, os.path.getsize(DST))
