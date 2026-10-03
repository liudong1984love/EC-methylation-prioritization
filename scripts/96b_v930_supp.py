# -*- coding: utf-8 -*-
"""Supplementary Notes v9.29 -> v9.30.
S1 S2.11: soften systematic-null wording (eligibility-based).
S2 S2.13: 'defined a priori' -> pre-declared before resampling (post hoc exploratory).
S3 S2.15: methods text updated to v2 (cache-complete gates, official score_v6, NXPH1 artefact).
S4 Table S11: replace 27 data rows with C4 v2 values; header note updated.
S5 S2.16: remove AUB-overestimate/conservative inference.
S6 New S2.17 + Table S14 (cross-analysis reconciliation), inserted after Table S13.
"""
import re
import pandas as pd
import docx
from docx.oxml.ns import qn

SRC = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Supplementary_Notes_EC_methylation_panel_v9.29.docx'
DST = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Supplementary_Notes_EC_methylation_panel_v9.30.docx'
RES = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\results'

doc = docx.Document(SRC)

def set_para_text(p, text):
    """Replace paragraph text, keeping paragraph style and first run's format."""
    runs = p.runs
    if not runs:
        p.add_run(text)
        return
    runs[0].text = text
    for r in runs[1:]:
        r.text = ''

def find_para(substr, start=0):
    for i, p in enumerate(doc.paragraphs):
        if i >= start and substr in p.text:
            return i, p
    raise AssertionError('paragraph not found: ' + substr[:60])

# ---------- S1: S2.11 soften ----------
i, p = find_para('No publicly available 450K/EPIC cohorts of endometrial polyps')
t = p.text
old = ('No publicly available 450K/EPIC cohorts of endometrial polyps, leiomyomas, adenomyosis '
       'or PCOS exist within the reachable repositories; this systematic null result '
       'substantiates the statement in the main text that intended-use benign backgrounds must '
       'be collected prospectively.')
new = ('Within the repositories searched, we identified no eligible publicly available 450K/EPIC '
       'cohort providing premenopausal benign-lesion (polyp, leiomyoma, adenomyosis or PCOS) '
       'endometrial methylation profiles with the required benign-diagnosis annotations '
       '(databases, search date and eligibility criteria as stated above); lesion tissue itself '
       'and endometrium from patients with these conditions are distinguished in this '
       'eligibility definition. This systematic null result substantiates the statement in the '
       'main text that intended-use benign backgrounds must be collected prospectively.')
assert old in t, 'S1 anchor'
set_para_text(p, t.replace(old, new))
print('S1 done')

# ---------- S2: S2.13 'a priori' fix ----------
i, p = find_para('five candidate region-level rules were defined a priori')
t = p.text
old = ('To replace it with a pre-specifiable rule, five candidate region-level rules were '
       'defined a priori and compared:')
new = ('To replace it with a pre-specifiable rule, five candidate region-level rules '
       '(post hoc exploratory rules, formulated after inspection of the external cohorts and '
       'pre-declared before this resampling comparison) were compared:')
assert old in t, 'S2 anchor'
set_para_text(p, t.replace(old, new))
print('S2 done')

# ---------- S3: S2.15 methods text -> v2 ----------
i, p = find_para('Rationale and methods. Single-probe discovery required FDR<0.05')
t = p.text
old = ('At each grid cell, the candidate set was re-derived genome-wide from the TCGA discovery '
       'statistics (485,512 probes), new entrants versus the strict 476 were identified, and all '
       'new entrants were evaluated against the same four background gates (union rule, n/a '
       'excluded; per-probe gate values in Supplementary Data S2q). Reproduction check: at the '
       'original criteria the genome-wide re-derivation returns 479 candidates, 476 of which are '
       'the audited elite set (the three additional probes map to chrX/chrY or a PCDHB-region '
       'shore and were excluded by the audited annotation filters).')
new = ('At each grid cell, the candidate set was re-derived genome-wide from the TCGA discovery '
       'statistics (485,512 probes), new entrants versus the strict 476 were identified, and '
       'candidates were evaluated against the same four background gates (union rule; n/a = not '
       'assessable on the platform, not failed). Gate values for the 476 audit-universe probes '
       'are the complete four-gate evaluations from the stability cache; new entrants were '
       'evaluated de novo in all four cohorts plus the endometriosis background used by the '
       'score (per-probe values in Supplementary Data S2q). Ranks use the official effective '
       'score (0.25' + '\u00d7' + 's_bg + 0.25' + '\u00d7' + 'premenopausal positivity + '
       '0.20' + '\u00d7' + 'overall positivity + 0.20; s_bg from the worst of the three '
       'endometrial backgrounds). An earlier version of this grid read gate values from the '
       'elite audit table, in which unevaluated cells (pipeline-stage completeness artefacts, '
       'e.g. cg15285250/NXPH1 lacking cervix and blood values that were in fact measured at '
       '0.186/0.276) were wrongly treated as missing, and used a simplified score; that version '
       'let NXPH1 pass at the original cell and ranked BOLL 3/3. The corrected analysis '
       'reproduces the main analysis exactly at the original cell: 479 candidates (the '
       '476-probe audit universe plus three probes excluded by the audited annotation filters, '
       'none passing triage), two triage survivors (ZSCAN12 rank 1/2, BOLL anchor rank 2/2). '
       'Cross-analysis reconciliation is provided in Supplementary Note S2.17 (Table S14).')
assert old in t, 'S3 anchor'
set_para_text(p, t.replace(old, new))
print('S3 done')

# ---------- S4: Table S11 values -> v2 ----------
hdr_i, hdr_p = find_para('Table S11. Discovery-threshold grid (27 combinations).')
old_hdr = 'Ranks are by the three-term effective score among triage survivors.'
new_hdr = 'Ranks are by the official effective score (see text) among triage survivors; gate values are cache-complete for the 476 and de novo four-cohort evaluations for new entrants (corrected version 2).'
assert old_hdr in hdr_p.text
set_para_text(hdr_p, hdr_p.text.replace(old_hdr, new_hdr))

g = pd.read_csv(RES + r'\v936_threshold_grid_summary_v2.csv')
# find Table S11: header row contains 'Candidates' and 'ZSCAN12 rank'
t11 = None
for t in doc.tables:
    h = [c.text.strip() for c in t.rows[0].cells]
    if 'Candidates' in h and any('ZSCAN12' in x for x in h):
        t11 = t
        break
assert t11 is not None and len(t11.rows) == 28, 'Table S11 not found or wrong size: %s' % (None if t11 is None else len(t11.rows))
def rank_of(s, lead):
    v = s[lead]
    m = re.search(r'rank (\d+/\d+)', v)
    return m.group(1) if m else 'n/a'
for ri, (_, r) in enumerate(g.iterrows(), start=1):
    vals = [f"{r['delta_beta']:.2f}", f"{r['case_pos']:.2f}", f"{r['pre_pos']:.2f}",
            str(int(r['n_candidates'])), str(int(r['n_new_vs_476'])),
            str(int(r['n_triage_pass'])), str(int(r['n_new_lowbg'])),
            rank_of(r, 'ZSCAN12_anchor'), rank_of(r, 'BOLL_anchor')]
    for ci, v in enumerate(vals):
        cell = t11.rows[ri].cells[ci]
        cell.paragraphs[0].runs[0].text = v if cell.paragraphs[0].runs else None
        if not cell.paragraphs[0].runs:
            cell.paragraphs[0].add_run(v)
        for extra in cell.paragraphs[0].runs[1:]:
            extra.text = ''
print('S4 Table S11 updated (27 rows)')

# ---------- S5: S2.16 conservative inference removed ----------
i, p = find_para('the endometriosis-enriched surgical cohort likely overestimates')
t = p.text
old = ('while its endometrial-brush suitability remains unresolved, because the '
       'endometriosis-enriched surgical cohort likely overestimates an AUB-screening '
       'background. Endometrial-path retentions are thus conservative, and cervical-path '
       'passes do not establish endometrial suitability.')
new = ('while its endometrial-brush suitability remains unresolved: because the surgical '
       'cohort lacks AUB annotation and differs from the intended-use population, the '
       'direction and magnitude of background bias relative to a premenopausal AUB cohort '
       'remain unknown, and cervical-path passes do not establish endometrial suitability.')
assert old in t, 'S5 anchor'
set_para_text(p, t.replace(old, new))
print('S5 done')

# ---------- S6: new S2.17 + Table S14 after Table S13 ----------
recon = pd.read_csv(RES + r'\v936_reconciliation.csv')

# build content at document end, then move after the Table S13 table
h = doc.add_paragraph()
r = h.add_run('S2.17 Cross-analysis reconciliation of candidate counts and ranks')
r.bold = True
p1 = doc.add_paragraph(
    'Different analyses in this study report different candidate and survivor counts because '
    'they ask different questions of the same data. Table S14 reconciles them: the main funnel '
    '(476-probe audit universe), the application-specific gating paths (S2.16), the '
    'discovery-threshold grid (S2.15) and the pipeline sensitivity analysis (S2.12). The key '
    'distinctions are the starting set (476 audited probes versus genome-wide re-derivation, '
    'which adds three annotation-filtered probes), the background rule applied at each stage '
    '(healthy-endometrium screen; triage on healthy/cervix/blood with benign-surgical '
    'descriptive; strict four-gate union with missing = not failed) and the score used for '
    'ranking (official effective score throughout the corrected analyses). With these '
    'definitions applied consistently, the original threshold-grid cell reproduces the main '
    'analysis exactly, and the "finalists" of the pipeline analysis are the stage-2 triage '
    'survivors of the main funnel (15), not the strict-union Tier A set (2).')
cap = doc.add_paragraph()
r = cap.add_run('Table S14. Cross-analysis reconciliation of candidate counts, background rules, survivor counts and lead ranks.')
r.bold = True
t14 = doc.add_table(rows=1 + len(recon), cols=len(recon.columns))
hdr = ['Analysis', 'Starting set', 'Discovery criteria', 'Background rule', 'Survivors', 'Lead ranks']
for j, txt in enumerate(hdr):
    c = t14.rows[0].cells[j]
    c.paragraphs[0].add_run(txt).bold = True
for i2, (_, row) in enumerate(recon.iterrows(), start=1):
    for j, col in enumerate(recon.columns):
        t14.rows[i2].cells[j].paragraphs[0].add_run(str(row[col]))

# move the four new elements after the table that follows the 'Table S13.' caption
cap_i, cap_p = find_para('Table S13. The 23-probe experimental pool under the three gating rules')
body = cap_p._p.getparent()
node = cap_p._p
t13_tbl = None
while node is not None:
    node = node.getnext()
    if node is not None and node.tag == qn('w:tbl'):
        t13_tbl = node
        break
assert t13_tbl is not None, 'Table S13 table element not found'
# reuse Table S13's table properties for the new table
import copy
t13_tblPr = t13_tbl.find(qn('w:tblPr'))
t14_tblPr = t14._tbl.find(qn('w:tblPr'))
if t13_tblPr is not None and t14_tblPr is not None:
    t14._tbl.replace(t14_tblPr, copy.deepcopy(t13_tblPr))
for el in [t14._tbl, cap._p, p1._p, h._p]:
    t13_tbl.addnext(el)
print('S6 S2.17 + Table S14 inserted')

doc.save(DST)

# ---------- verification ----------
d2 = docx.Document(DST)
full = '\n'.join(p.text for p in d2.paragraphs)
for k in ['we identified no eligible publicly available 450K/EPIC cohort',
          'pre-declared before this resampling comparison',
          'reproduces the main analysis exactly at the original cell',
          'NXPH1', 'S2.17 Cross-analysis reconciliation',
          'Table S14', 'corrected version 2']:
    print(k[:55], '->', k in full)
print('residual likely overestimates:', full.count('likely overestimates'))
print('residual defined a priori:', full.count('defined a priori'))
print('tables now:', len(d2.tables))
# verify Table S11 first/last data rows
for t in d2.tables:
    hcells = [c.text.strip() for c in t.rows[0].cells]
    if 'Candidates' in hcells and any('ZSCAN12' in x for x in hcells):
        print('S11 row1:', [c.text for c in t.rows[1].cells])
        print('S11 row27:', [c.text for c in t.rows[27].cells])
print('saved:', DST)
