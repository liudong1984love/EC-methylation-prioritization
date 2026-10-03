# -*- coding: utf-8 -*-
"""Manuscript v9.28 -> v9.29: rule-closure revision.
E1 Methods: resampling description -> nested/conditional (was 'full funnel ... independently').
E2 Results: blood OR-rule denominator 0/42 -> 0/60 (verified on all ten fractions).
E3 Results: pool accounting (7+1+15) + cg05639246 rationale (overall, not premenopausal, coverage).
E4 Discussion: BOLL readout positioning unified (anchor primary exploratory; region mean candidate,
   internal stability only; 2.0% benign proportion != P95 gate pass).
E5 Limitations: region-mean pre-specification scope (internal stability, not validation).
E6 Results: discovery-threshold sensitivity paragraph (new, exploratory).
E7 Results: sampling-site parallel paths paragraph (new, exploratory).
E8 Figure 1 legend: point to retention ledger (S2.14/Table S10/S2p).
E9 Replace embedded Figure 1 image (corrected counts: 461/476; 15 region-coverage, 14/15).
"""
import zipfile, re

SRC = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Manuscript_EC_methylation_panel_v9.28.docx'
DST = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Manuscript_EC_methylation_panel_v9.29.docx'
FIG = r'C:\Users\ld\WorkBuddy\2026-09-25-20-09-18\fig1_v6.png'

B = '\u03b2'; LE = '\u2264'; GE = '\u2265'; MD = '\u2014'; DQ = '\u2019'

def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

zin = zipfile.ZipFile(SRC)
xml = zin.read('word/document.xml').decode('utf-8')

REPL = []

# ---- E1: Methods resampling wording ----
old1 = ('Stability of the multi-step selection itself was assessed by resampling the full funnel: '
        '2,000 stratified bootstrap replicates in which TCGA tumours, adjacent normals and premenopausal '
        'tumours were resampled independently, each background cohort was resampled at sample level '
        '(GSE35069 within cell-type fractions), and discovery criteria, background gates, the effective '
        'three-component score and candidate ranks were recomputed per replicate, conditional on the '
        '476-probe audit universe (Supplementary Note S2.9; Supplementary Data S2f).')
new1 = ('Stability of the multi-step selection itself was assessed by conditional resampling within the '
        'fixed 476-probe audit universe (not a re-run of genome-wide discovery): 2,000 bootstrap '
        'replicates in which premenopausal (n=28) and other (n=403) tumours were resampled as separate '
        'strata and pooled' + MD + 'so that overall and premenopausal positivity derive from the same '
        'nested draw' + MD + 'adjacent normals and each background cohort were resampled at sample level '
        '(GSE35069 within cell-type fractions), and discovery criteria, background gates, the effective '
        'three-component score and candidate ranks were recomputed per replicate (Supplementary Note '
        'S2.9; Supplementary Data S2f).')
REPL.append(('E1', old1, new1))

# ---- E2: blood denominator ----
old2 = 'purified blood-cell fractions (0/42; six donors, none positive in any fraction)'
new2 = ('blood-reference fractions (0/60; ten fractions ' + '\u00d7' + ' six donors, none positive in '
        'any fraction)')
REPL.append(('E2', old2, new2))

# ---- E3: pool accounting + cg05639246 rationale ----
old3 = ('A/B leads plus 16 region-coverage probes '
        'retained for regional characterisation and incremental training-set coverage, 15 of which '
        'additionally showed cervical and blood P95 ' + B + ' ' + LE + '0.30 and small between-phase '
        'differences' + MD + 'a descriptive post-assembly audit, not an applied gate; the exception, '
        'cg05639246, cervical P95 0.441, was retained because its region contributed incremental '
        'premenopausal-tumour coverage in the greedy assembly; Supplementary Note S2.6;')
new3 = ('A/B leads, the BOLL partner probe and 15 '
        'region-coverage probes retained for regional characterisation and incremental training-set '
        'coverage; of these 16 regional and coverage probes, 15 additionally showed cervical and blood '
        'P95 ' + B + ' ' + LE + '0.30 and small between-phase differences' + MD + 'a descriptive '
        'post-assembly audit, not an applied gate; the exception, cg05639246, cervical P95 0.441, was '
        'retained because its region contributed incremental overall-tumour coverage in the greedy '
        'assembly, premenopausal coverage being already complete before its addition; per-probe '
        'retention ledger in Supplementary Note S2.14 and Table S10; Supplementary Note S2.6;')
REPL.append(('E3', old3, new3))

# ---- E4: BOLL readout positioning (Discussion) ----
old4 = ('The region-mean rule is therefore pre-specified as the BOLL region-level readout for the '
        'prospective phase, and the OR rule is abandoned.')
new4 = ('Accordingly, the BOLL anchor remains the primary exploratory readout, and the region mean is '
        'carried forward as a post-hoc-derived candidate region-level readout' + MD + 'selected on '
        'internal resampling stability only (the bootstrap and leave-one-out frequencies measure '
        'rule-selection stability, not independent validation)' + MD + 'to be evaluated together with '
        'its threshold in the independent development cohort and locked before any validation; the OR '
        'rule is abandoned. The region-mean rule' + DQ + 's 2.0% benign positivity is a rule-level '
        + B + '>0.3 proportion and does not constitute a pass of the probe-level P95 ' + LE +
        '0.15 background gate; it belongs to the region-level exploratory path, not to the strict '
        'probe funnel.')
REPL.append(('E4', old4, new4))

# ---- E5: Limitations scope of pre-specification ----
old5 = ('which is now pre-specified for the threshold-setting cohort and remains to be validated there.')
new5 = ('which is pre-specified only as the candidate readout to be evaluated and locked, with its '
        'threshold, in the threshold-setting cohort; the resampling frequencies measure internal '
        'rule-selection stability, not independent validation.')
REPL.append(('E5', old5, new5))

# ---- E8: Figure 1 legend -> ledger pointer ----
old8 = 'n/a marks gates not assessable for a probe (missing, not failed).'
new8 = ('n/a marks gates not assessable for a probe (missing, not failed). The per-probe retention '
        'ledger (entry path, gate outcomes, evidence completeness, retention rationale, intended use) '
        'is provided in Supplementary Note S2.14 (Table S10; Supplementary Data S2p).')
REPL.append(('E8', old8, new8))

for tag, old, new in REPL:
    n = xml.count(old)
    assert n == 1, '%s anchor count = %d' % (tag, n)
    m = re.search(re.escape(old), xml)
    s = xml.rfind('<w:t', 0, m.start())
    e = xml.find('</w:t>', m.end())
    seg = xml[s:e]
    assert '</w:t>' not in seg[len('<w:t'):], '%s crosses text nodes' % tag
    xml = xml.replace(old, new, 1)
    print(tag, 'replaced')

# ---- paragraph insertions (E6, E7) ----
def para_xml(lead, rest):
    lead_run = (f'<w:r><w:rPr><w:b w:val="1"/></w:rPr><w:t xml:space="preserve">{esc(lead)}</w:t></w:r>')
    return f'<w:p>{lead_run}<w:r><w:t xml:space="preserve">{esc(rest)}</w:t></w:r></w:p>'

E6 = para_xml('Discovery-threshold sensitivity (exploratory). ',
    'Because the leads sit close to the discovery cut-offs, the cut-offs themselves were scanned over a '
    'pre-specified, one-sided 27-combination grid loosening ' + '\u0394' + B + ' (' + GE + '0.50 to '
    + GE + '0.30), overall tumour positivity (' + GE + '90% to ' + GE + '80%) and premenopausal '
    'positivity (' + GE + '85% to ' + GE + '75%), with FDR and the adjacent-normal P95 criterion held '
    'fixed (Supplementary Note S2.15). Relaxation enlarged the candidate set from 479 to 2,983 probes, '
    'but the four-gate background triage eliminated nearly all entrants: at the original ' + GE +
    '90% overall-positivity requirement no new low-background candidate entered in any grid cell, and '
    'at the loosest corner only 17 of 2,507 new entrants passed triage. cg27577527 (ZSCAN12) ranked '
    'first among triage survivors in all 27 combinations; cg24589459 (BOLL anchor) survived triage in '
    'all combinations (rank 3' + MD + '11, reflecting its threshold-adjacent cervical background rather '
    'than signal loss). Approved-kit and reference markers showed distinct exclusion modes: CDO1 passes '
    'all discovery criteria and is eliminated only at background triage; AJAP1 is excluded by the '
    'overall-positivity criterion (84.2%); GALR1 and GYPC enter at no grid cell. Single-probe positivity '
    'below 90% therefore does not imply absence of complementary value' + MD + 'but relaxing the '
    'positivity requirement does not rescue the kit comparators either (Supplementary Data S2l, S2m).')

E7 = para_xml('Application-specific gating paths (exploratory). ',
    'The strict four-gate union selects cross-compartment low-background markers and is more stringent '
    'than any single intended application requires. In a parallel analysis that does not alter the main '
    'rule (Supplementary Note S2.16), two specimen-matched gating paths were evaluated: an '
    'endometrial-brush path (gates: healthy and benign-surgical endometrium; blood reported as '
    'bleeding-contamination risk) retained 11 of 476 candidates, and a cervical/self-sampling path '
    '(gates: HPV-negative cervix and blood) retained 43, versus 2 under the strict union. The paths '
    'partition the experimental pool informatively: CDO1' + MD + 'eliminated by the benign-surgical '
    'dimension (P95 0.372)' + MD + 'passes the cervical path (cervix 0.055, blood 0.080), consistent '
    'with its established cervical application, whereas six pool probes show the mirror pattern '
    '(endometrium-passing, cervix/blood-elevated; Table S13; Supplementary Data S2n). Because the '
    'benign-surgical cohort is endometriosis-enriched and lacks AUB annotation, its background likely '
    'overestimates an AUB-screening background; endometrial-path retentions are therefore conservative, '
    'and cervical-path passes do not establish endometrial suitability. Path assignments guide '
    'specimen-matched assay development; the strict four-gate rule remains the main analysis.')

def insert_after_para(anchor_tail, new_para, tag):
    global xml
    n = xml.count(anchor_tail)
    assert n == 1, '%s anchor count = %d' % (tag, n)
    i = xml.find(anchor_tail)
    e = xml.find('</w:p>', i)
    assert e > 0
    xml = xml[:e + len('</w:p>')] + new_para + xml[e + len('</w:p>'):]
    print(tag, 'inserted')

insert_after_para('Final cut-offs must be established in a separate development cohort and locked '
                  'before evaluation in an independent validation cohort.', E6, 'E6')
insert_after_para('this probe-level audit does not evaluate the exact CpGs, amplicons or '
                  'classification algorithms used in commercial or published assays.', E7, 'E7')

# ---- E9: swap Figure 1 image ----
fig_bytes = open(FIG, 'rb').read()
with zipfile.ZipFile(DST, 'w', zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == 'word/document.xml':
            data = xml.encode('utf-8')
        elif item.filename == 'word/media/image1.png':
            data = fig_bytes
        zout.writestr(item, data)
zin.close()

z = zipfile.ZipFile(DST)
assert z.testzip() is None
import xml.etree.ElementTree as ET
ET.fromstring(z.read('word/document.xml'))
x = z.read('word/document.xml').decode('utf-8')
for k in ['conditional resampling within the fixed 476-probe audit universe', '0/60',
          'incremental overall-tumour coverage', 'post-hoc-derived candidate region-level readout',
          'Discovery-threshold sensitivity (exploratory)', 'Application-specific gating paths',
          'per-probe retention ledger', 'S2.14', 'S2.15', 'S2.16', 'Table S10', 'Table S13']:
    print(k, '->', k in x)
print('saved:', DST)
