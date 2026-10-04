# -*- coding: utf-8 -*-
"""Supplementary Notes v9.30 -> v9.31.
S1  S2.15: uniform-filter description + version-3 reproduction statement.
S2  Table S11 caption: version 3 note; replace 27 grid rows with v938 values.
S3  S2.17 text: 'adds three annotation-filtered probes' -> '479 before / 476 after'.
S4  Table S14 row 3 cells updated (v938 reconciliation wording).
S5  New note S2.18: verified GSE35069 donor identities (Data S2s).
"""
import zipfile, re
import xml.etree.ElementTree as ET
import pandas as pd

SRC = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Supplementary_Notes_EC_methylation_panel_v9.30.docx'
DST = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Supplementary_Notes_EC_methylation_panel_v9.31.docx'
RES = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\results'

MD = '\u2014'; LE = '\u2264'; GE = '\u2265'; DL = '\u0394'; B = '\u03b2'

def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

zin = zipfile.ZipFile(SRC)
xml = zin.read('word/document.xml').decode('utf-8')

REPL = []

# ---- S1a: uniform filters in S2.15 methods ----
REPL.append(('S1a',
  'with FDR and the adjacent-normal criterion held fixed.',
  'with FDR, the adjacent-normal criterion and the audited cross-cohort-consistency and '
  'annotation filters held fixed; the filters (GSE67116 directional consistency ' + DL + B +
  '>0, and presence in the noob-processed GSE73949 background matrix' + MD + 'through which '
  'the minfi sex-chromosome, SNP and detection-P exclusions enter) were applied identically '
  'to every grid cell, so loosened cells cannot admit probes the main pipeline would have '
  'excluded.'))

# ---- S1b: version-3 reproduction statement ----
REPL.append(('S1b',
  'The corrected analysis reproduces the main analysis exactly at the original cell: 479 '
  'candidates (the 476-probe audit universe plus three probes excluded by the audited '
  'annotation filters, none passing triage), two triage survivors (ZSCAN12 rank 1/2, BOLL '
  'anchor rank 2/2).',
  'The corrected analysis (version 3, with the uniform filters) reproduces the main '
  'analysis exactly at the original cell: 479 candidates before the audited '
  'cross-cohort-consistency and annotation filters and 476 after' + MD + 'the exact audit '
  'universe' + MD + 'with two triage survivors (ZSCAN12 rank 1/2, BOLL anchor rank 2/2); '
  'the three filter-excluded probes are cg17837162 (chrY) and cg09687738 (chrX), absent '
  'from the minfi-processed GSE67116/GSE73949 matrices, and cg18661205 (GSE67116 ' + DL +
  B + ' 0.108>0 but absent from the GSE73949 matrix).'))

# ---- S2a: Table S11 caption ----
REPL.append(('S2a',
  'gate values are cache-complete for the 476 and de novo four-cohort evaluations for new '
  'entrants (corrected version 2).',
  'gate values are cache-complete for the 476 and de novo four-cohort evaluations for new '
  'entrants (corrected version 3; uniform cross-cohort-consistency and annotation filters '
  'applied to every cell).'))

# ---- S3: S2.17 text ----
REPL.append(('S3',
  'the starting set (476 audited probes versus genome-wide re-derivation, which adds three '
  'annotation-filtered probes)',
  'the starting set (476 audited probes versus genome-wide re-derivation, which yields 479 '
  'probes before the audited cross-cohort-consistency and annotation filters and the same '
  '476 after them)'))

for tag, old, new in REPL:
    n = xml.count(old)
    assert n == 1, '%s anchor count = %d' % (tag, n)
    xml = xml.replace(old, new, 1)
    print(tag, 'replaced')

with zipfile.ZipFile(DST, 'w', zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == 'word/document.xml':
            data = xml.encode('utf-8')
        zout.writestr(item, data)
zin.close()

# ---------- python-docx stage: Table S11 rows, Table S14 row, new S2.18 ----------
import docx
from docx.oxml.ns import qn
doc = docx.Document(DST)

g = pd.read_csv(RES + r'\v938_threshold_grid_summary_v3.csv')
recon = pd.read_csv(RES + r'\v938_reconciliation.csv')

# Table S11 = tables[11] (28 rows x 9 cols). Rebuild rows.
t11 = doc.tables[11]
assert len(t11.rows) == 28 and len(t11.columns) == 9, 'S11 shape %d x %d' % (len(t11.rows), len(t11.columns))
hdr = [c.text.strip() for c in t11.rows[0].cells]
print('S11 header:', hdr)
# map grid columns to table columns by position (assumes same order as v930 fill)
for ri in range(1, 28):
    rec = g.iloc[ri - 1]
    vals = [str(rec['delta_beta']), str(rec['case_pos']), str(rec['pre_pos']),
            str(int(rec['n_candidates'])), str(int(rec['n_new_vs_476'])),
            str(int(rec['n_triage_pass'])), str(int(rec['n_new_lowbg'])),
            str(rec['BOLL_anchor']), str(rec['ZSCAN12_anchor'])]
    for ci, v in enumerate(vals):
        cell = t11.rows[ri].cells[ci]
        # keep first paragraph, rewrite runs
        p = cell.paragraphs[0]
        for r in list(p.runs):
            r._r.getparent().remove(r._r)
        p.add_run(v)
print('S11 rows replaced')

# Table S14 = tables[14] (5 rows x 6 cols); row 3 = threshold grid original cell
t14 = doc.tables[14]
print('S14 row3 first cell:', t14.rows[3].cells[0].text.strip()[:60])
r3 = recon.iloc[2]
vals = [str(r3[c]) for c in recon.columns]
for ci, v in enumerate(vals):
    cell = t14.rows[3].cells[ci]
    p = cell.paragraphs[0]
    for r in list(p.runs):
        r._r.getparent().remove(r._r)
    p.add_run(v)
# also update row 4 (pipeline sensitivity: 'finalists 15 -> 10' wording)
r4 = recon.iloc[3]
vals = [str(r4[c]) for c in recon.columns]
for ci, v in enumerate(vals):
    cell = t14.rows[4].cells[ci]
    p = cell.paragraphs[0]
    for r in list(p.runs):
        r._r.getparent().remove(r._r)
    p.add_run(v)
print('S14 rows 3-4 updated')

# S14 row 1 (main funnel) also got updated wording in v938 (start_set/discovery) - update rows 1-2 too for consistency
for rix in [1, 2]:
    rec = recon.iloc[rix - 1]
    vals = [str(rec[c]) for c in recon.columns]
    for ci, v in enumerate(vals):
        cell = t14.rows[rix].cells[ci]
        p = cell.paragraphs[0]
        for r in list(p.runs):
            r._r.getparent().remove(r._r)
        p.add_run(v)
print('S14 rows 1-2 refreshed')

# ---- S5: new note S2.18 at end of document ----
sec_h = doc.add_paragraph()
run = sec_h.add_run('S2.18 Verified donor identities in GSE35069')
run.bold = True
body = doc.add_paragraph(
  'GSE35069 comprises ten blood fractions from the same six donors. Donor correspondence was '
  'recovered from the ArrayExpress record E-GEOD-35069 (sample-data-relationship file): sample '
  'titles encode the fraction and donor (e.g. CD14+_1 ' + MD + ' CD14+_6), yielding a complete, '
  'duplicate-free 6 ' + '\u00d7' + ' 10 donor-by-fraction mapping. The verified mapping is provided '
  'in Supplementary Data S2s. Re-running the 10,000-draw donor-clustered bootstrap '
  '(Methods; Supplementary Data S2c) with the verified mapping reproduced the earlier '
  'order-inferred intervals exactly (point estimates and all interval bounds identical to three '
  'decimals) and changed no gate decision. The earlier within-fraction order inference is '
  'therefore superseded by verified metadata; the clustered intervals remain sensitivity-level '
  'because six donors cannot support precise cluster-level uncertainty estimates.')
print('S2.18 added')

doc.save(DST)

# verify
import docx as dx
d2 = dx.Document(DST)
t = d2.tables[11]
print('S11 row1:', [c.text for c in t.rows[1].cells][:7])
print('S11 row27:', [c.text for c in t.rows[27].cells][:7])
t14b = d2.tables[14]
print('S14 r3:', t14b.rows[3].cells[0].text[:60], '|', t14b.rows[3].cells[4].text[:60])
z = zipfile.ZipFile(DST)
assert z.testzip() is None
ET.fromstring(z.read('word/document.xml'))
x = z.read('word/document.xml').decode('utf-8')
for k in ['applied identically to every grid cell', '479 candidates before the audited',
          'corrected version 3', '479 probes before the audited', 'S2.18 Verified donor identities',
          'Supplementary Data S2s']:
    print(k[:55], '->', k in x)
print('saved:', DST)
