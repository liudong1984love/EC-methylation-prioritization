# -*- coding: utf-8 -*-
"""Manuscript v9.27 -> v9.28: three insertions reporting items 1 & 3.
1. Discussion para (BOLL novelty): append region-mean pre-specification after the OR-rule sentence.
2. Limitations para: quantify pipeline sensitivity after the "pipeline-bound" sentence.
3. Limitations final sentence: update the post-hoc BOLL rule verdict (region-mean now pre-specified;
   45.8% / 159/347 confirmed by recomputation on 3-decimal source betas).
4. Figure 4 legend: point to S2.13.
"""
import zipfile, re

SRC = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Manuscript_EC_methylation_panel_v9.27.docx'
DST = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Manuscript_EC_methylation_panel_v9.28.docx'

def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

def run(text):
    return f'<w:r><w:t xml:space="preserve">{esc(text)}</w:t></w:r>'

zin = zipfile.ZipFile(SRC)
xml = zin.read('word/document.xml').decode('utf-8')

EDITS = []

# --- 1. BOLL discussion paragraph (para 56): after OR-rule sentence ---
a1 = 'a simple OR rule was frequently positive in benign surgical endometrium. The expression analysis'
ins1 = ('a simple OR rule was frequently positive in benign surgical endometrium. To replace this '
        'post-hoc concept with a pre-specifiable rule, five candidate region-level rules (anchor only, '
        'partner only, OR, AND and region mean; positivity \u03b2>0.3) were compared by stratified '
        'bootstrap (B=2000) and benign leave-one-out selection (Supplementary Note S2.13): the '
        'region-mean rule was selected in 65.1% of bootstrap replicates and 100% of benign leave-one-out '
        'folds (anchor-only 34.9%; OR, AND and partner-only never selected), retaining 92.9% premenopausal '
        'tumour positivity with 2.0% benign surgical positivity, versus 45.8% benign positivity under the '
        'OR rule (159/347, partner-driven). The region-mean rule is therefore pre-specified as the BOLL '
        'region-level readout for the prospective phase, and the OR rule is abandoned. The expression analysis')
assert xml.count(a1) == 1, 'anchor 1 not unique: %d' % xml.count(a1)
# anchor spans two runs? check; if single text node, simple replace works
EDITS.append((a1, ins1))

# --- 2. Limitations: after pipeline-bound sentence ---
a2 = 'and near-gate pass/fail status was pipeline-bound (Results), so the fixed thresholds are heuristic.'
ins2 = ('and near-gate pass/fail status was pipeline-bound (Results), so the fixed thresholds are '
        'heuristic. Recomputing the healthy-endometrium gate under the GMQN normalisation used for the '
        'other background dimensions (same 17 donors) shifted the P95 upward for a minority of '
        'probes\u201424/476 pass\u2013fail flips, all in the conservative (noob-pass-to-GMQN-fail) '
        'direction, Cohen\u2019s \u03ba=0.747\u2014reducing the healthy-pass set from 65 to 41 and the '
        'finalists from 15 to 10; all five finalist losses were marginal (GMQN P95 0.156\u20130.168) and '
        'included CDO1 (cg23180938) and the BOLL partner (cg07495363), whereas both leads passed under '
        'either pipeline (Supplementary Note S2.12).')
assert xml.count(a2) == 1, 'anchor 2 not unique: %d' % xml.count(a2)
EDITS.append((a2, ins2))

# --- 3. Limitations final BOLL-rule sentence ---
a3 = ("The post hoc two-probe BOLL rule was derived from the same external data and inherits the partner "
      "probe's benign-endometrium background (45.8% rule positivity in the benign surgical cohort); "
      "it is a hypothesis to be locked and tested in a threshold-setting cohort, not a validated rule.")
ins3 = ('The two-probe BOLL rule was derived from the same external data and inherits the partner '
        'probe\u2019s benign-endometrium background (45.8% OR-rule positivity in the benign surgical '
        'cohort, 159/347); resampling-based selection over five pre-declared candidate rules favoured the '
        'region-mean rule (bootstrap 65.1%, benign leave-one-out 100%; Supplementary Note S2.13), which '
        'is now pre-specified for the threshold-setting cohort and remains to be validated there.')
assert xml.count(a3) == 1, 'anchor 3 not unique: %d' % xml.count(a3)
EDITS.append((a3, ins3))

# --- 4. Figure 4 legend tail ---
a4 = ('a simple OR rule was frequently positive in benign surgical endometrium, so region-level assay '
      'design remains a hypothesis to be tested prospectively.')
ins4 = ('a simple OR rule was frequently positive in benign surgical endometrium (45.8% benign '
        'positivity); the region-mean rule\u2014stably selected over OR, AND and single-probe '
        'alternatives by bootstrap and leave-one-out resampling (Supplementary Note S2.13)\u2014is '
        'pre-specified for prospective testing.')
assert xml.count(a4) == 1, 'anchor 4 not unique: %d' % xml.count(a4)
EDITS.append((a4, ins4))

for old, new in EDITS:
    # anchors must live inside a single w:t node for safe replacement
    m = re.search(re.escape(old), xml)
    s = xml.rfind('<w:t', 0, m.start())
    e = xml.find('</w:t>', m.end())
    seg = xml[s:e]
    assert '</w:t>' not in seg[len('<w:t'):], 'anchor crosses text nodes: %s' % old[:60]
    xml = xml.replace(old, new, 1)

with zipfile.ZipFile(DST, 'w', zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == 'word/document.xml':
            data = xml.encode('utf-8')
        zout.writestr(item, data)
zin.close()

z = zipfile.ZipFile(DST)
assert z.testzip() is None
import xml.etree.ElementTree as ET
ET.fromstring(z.read('word/document.xml'))
x = z.read('word/document.xml').decode('utf-8')
for k in ['S2.12', 'S2.13', 'region-mean rule', '65.1%', '45.8%', '\u03ba=0.747', '159/347']:
    print(k, '->', k in x)
print('saved:', DST)
