# fix_supp_v932.py — Supplementary Notes v9.31 -> v9.32
# S2.9 rewritten to donor-level blood resampling (v3); all frequencies updated;
# ablation (nopre) updated; EOMES reserve updated; within-fraction vs donor-level robustness note added.
import zipfile, os

SRC = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Supplementary_Notes_EC_methylation_panel_v9.31.docx'
DST = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Supplementary_Notes_EC_methylation_panel_v9.32.docx'

zin = zipfile.ZipFile(SRC)
xml = zin.read('word/document.xml').decode('utf-8')

GE = '\u2265'
DB = '\u0394\u03b2'
MD = '\u2014'

REPL = [
 # --- Design: donor-level blood resampling ---
 ('GSE35069 within each of the 10 blood-reference fractions, n=6 per fraction)',
  'GSE35069 at donor level: the six verified donors' + MD + 'identities taken from ArrayExpress '
  'sample metadata (Note S2.18)' + MD + 'resampled with replacement, each drawn donor contributing '
  'all ten blood-reference fractions, with per-fraction P95 then maximum across fractions, so that '
  'cross-fraction donor correlation is preserved)'),
 # --- Tier-A membership numbers (v3) ---
 ('cg27577527 (ZSCAN12 anchor) re-entered Tier A in 53.1% of replicates',
  'cg27577527 (ZSCAN12 anchor) re-entered Tier A in 52.8% of replicates'),
 ('fell below the ' + GE + '0.90 criterion in 45.1% of replicates and premenopausal positivity below '
  + GE + '0.85 in 4.9%',
  'fell below the ' + GE + '0.90 criterion in 45.9% of replicates and premenopausal positivity below '
  + GE + '0.85 in 4.6%'),
 ('it ranked first in all 1,061 retained replicates',
  'it ranked first in all 1,055 retained replicates'),
 ('cg24589459 (BOLL anchor) re-entered Tier A in 48.3% of replicates, with attrition',
  'cg24589459 (BOLL anchor) re-entered Tier A in 48.7% of replicates, with attrition'),
 ('overall tumour positivity fell below ' + GE + '0.90 in 23.1% of replicates (observed 0.909), '
  'premenopausal positivity below ' + GE + '0.85 in 16.1%, ' + DB + ' below ' + GE + '0.50 in 0.4%, '
  'and the borderline cervical gate failed in 28.4% (71.6% compliance',
  'overall tumour positivity fell below ' + GE + '0.90 in 21.7% of replicates (observed 0.909), '
  'premenopausal positivity below ' + GE + '0.85 in 16.8%, ' + DB + ' below ' + GE + '0.50 in 0.2%, '
  'and the borderline cervical gate failed in 26.6% (73.5% compliance'),
 ('The premenopausal criterion and the cervical gate jointly held in 60.0% of replicates' + MD +
  'consistent with the Fr\u00e9chet lower bound (55.6%) for their marginal rates' + MD + 'and the '
  'remaining gap to the 48.3% Tier-A frequency',
  'The premenopausal criterion and the cervical gate jointly held in 61.2% of replicates' + MD +
  'consistent with the Fr\u00e9chet lower bound (56.7%) for their marginal rates' + MD + 'and the '
  'remaining gap to the 48.7% Tier-A frequency'),
 ('it ranked ' + '\u2264' + '2 in 70.0% and ' + '\u2264' + '5 in 99.8% of its 966 retained replicates',
  'it ranked ' + '\u2264' + '2 in 69.5% and ' + '\u2264' + '5 in 99.7% of its 974 retained replicates'),
 ('Both leads were jointly Tier A and occupied the top two ranks in 15.6% of replicates',
  'Both leads were jointly Tier A and occupied the top two ranks in 14.8% of replicates'),
 # --- Ablation (nopre, donor-level) ---
 ('raised the Tier-A frequency of cg24589459 from 48.3% to 57.1% and of cg27577527 from 53.1% to '
  '54.9% (joint top-2 frequency 15.6%' + '\u2192' + '23.9%)',
  'raised the Tier-A frequency of cg24589459 from 48.7% to 55.4% and of cg27577527 from 52.8% to '
  '52.9% (joint top-2 frequency 14.8%' + '\u2192' + '21.0%)'),
 # --- EOMES reserve ---
 ('Tier A in 38.3% of replicates, limited mainly by the healthy-endometrium gate at 38.5% compliance',
  'Tier A in 37.5% of replicates, limited mainly by the healthy-endometrium gate at 37.7% compliance'),
 # --- Robustness note: within-fraction variant ---
 ('its background evidence is therefore not deterministic either. These frequencies quantify '
  'conditional selection stability',
  'its background evidence is therefore not deterministic either. A within-fraction blood-resampling '
  'variant (six samples drawn independently inside each fraction, the design used before donor '
  'identities were verified) gave nearly identical lead retention (53.1% versus 52.8% for '
  'cg27577527; 48.3% versus 48.7% for cg24589459), so the stability results are insensitive to how '
  'residual blood-donor correlation is handled. These frequencies quantify conditional selection '
  'stability'),
]

for old, new in REPL:
    n = xml.count(old)
    print(n, '|', old[:70])
    assert n == 1, f'anchor count {n} != 1: {old[:60]}'
    xml = xml.replace(old, new, 1)

for bad in ['53.1% of replicates', '48.3% of replicates', '45.1% of replicates', '23.1% of replicates',
            '16.1%', '28.4%', '71.6% compliance', '15.6%', '57.1%', '54.9%', '23.9%', '38.3%',
            '38.5% compliance', '1,061', '966 retained', 'within each of the 10 blood-reference']:
    c = xml.count(bad)
    if c: print('RESIDUAL', bad, c)

with zipfile.ZipFile(DST + '.tmp', 'w', zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == 'word/document.xml':
            data = xml.encode('utf-8')
        zout.writestr(item, data)
zin.close()
os.replace(DST + '.tmp', DST)

import xml.etree.ElementTree as ET
z = zipfile.ZipFile(DST)
assert z.testzip() is None
ET.fromstring(z.read('word/document.xml'))
print('OK ->', DST, os.path.getsize(DST))
