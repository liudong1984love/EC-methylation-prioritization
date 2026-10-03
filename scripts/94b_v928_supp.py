# -*- coding: utf-8 -*-
"""Supplementary v9.27 -> v9.28: append S2.12 (cross-pipeline gate sensitivity)
and S2.13 (pre-specified BOLL region combination rule)."""
import zipfile, shutil, re, os
import pandas as pd

SRC = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Supplementary_Notes_EC_methylation_panel_v9.27.docx'
DST = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Supplementary_Notes_EC_methylation_panel_v9.28.docx'
RES = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\results'
OUTDIR = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19'

def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

def heading(text):
    return ('<w:p><w:pPr><w:keepNext w:val="1"/><w:keepLines w:val="1"/>'
            '<w:spacing w:before="100" w:beforeAutospacing="1" w:after="100" w:afterAutospacing="1"/>'
            '<w:outlineLvl w:val="1"/><w:rPr><w:sz w:val="32"/><w:szCs w:val="32"/></w:rPr></w:pPr>'
            f'<w:r><w:rPr><w:b w:val="1"/><w:sz w:val="32"/></w:rPr>'
            f'<w:t xml:space="preserve">{esc(text)}</w:t></w:r></w:p>')

def para(lead, rest):
    lead_run = f'<w:r><w:rPr><w:b w:val="1"/></w:rPr><w:t xml:space="preserve">{esc(lead)}</w:t></w:r>' if lead else ''
    return f'<w:p>{lead_run}<w:r><w:t xml:space="preserve">{esc(rest)}</w:t></w:r></w:p>'

def caption(text):
    return ('<w:p><w:pPr><w:keepNext w:val="1"/></w:pPr>'
            f'<w:r><w:rPr><w:b w:val="1"/></w:rPr><w:t xml:space="preserve">{esc(text)}</w:t></w:r></w:p>')

TBLPR = ('<w:tblPr><w:tblW w:w="0" w:type="auto"/><w:tblBorders>'
         '<w:top w:val="single" w:sz="12" w:space="0" w:color="000000"/>'
         '<w:left w:val="none"/><w:bottom w:val="single" w:sz="12" w:space="0" w:color="000000"/>'
         '<w:right w:val="none"/><w:insideH w:val="none"/><w:insideV w:val="none"/></w:tblBorders>'
         '<w:tblLayout w:type="autofit"/><w:tblCellMar>'
         '<w:top w:w="15" w:type="dxa"/><w:left w:w="15" w:type="dxa"/>'
         '<w:bottom w:w="15" w:type="dxa"/><w:right w:w="15" w:type="dxa"/></w:tblCellMar></w:tblPr>')

def table(rows):
    ncol = len(rows[0])
    grid = ''.join('<w:gridCol w:w="%d"/>' % (8306 // ncol) for _ in range(ncol))
    out = ['<w:tbl>', TBLPR, '<w:tblGrid>', grid, '</w:tblGrid>']
    for i, row in enumerate(rows):
        out.append('<w:tr>')
        for cell in row:
            if i == 0:
                tcpr = ('<w:tcPr><w:tcBorders><w:top w:val="none"/><w:left w:val="none"/>'
                        '<w:bottom w:val="single" w:sz="6" w:space="0" w:color="000000"/>'
                        '<w:right w:val="none"/></w:tcBorders></w:tcPr>')
                run = f'<w:r><w:rPr><w:b w:val="1"/></w:rPr><w:t xml:space="preserve">{esc(str(cell))}</w:t></w:r>'
            else:
                tcpr = '<w:tcPr/>'
                run = f'<w:r><w:t xml:space="preserve">{esc(str(cell))}</w:t></w:r>'
            out.append(f'<w:tc>{tcpr}<w:p>{run}</w:p></w:tc>')
        out.append('</w:tr>')
    out.append('</w:tbl>')
    out.append('<w:p/>')
    return ''.join(out)

B = '\u03b2'; LE = '\u2264'; KAPPA = '\u03ba'; GT = '\u2265'

# ---------- data for tables ----------
v931 = pd.read_csv(RES + r'\v931_pipeline_sensitivity_gse73949.csv')
audit = pd.read_csv(OUTDIR + r'\Supplementary_Data_S2a_single_probe_audit.csv')
tri = audit.loc[audit['funnel_pass_triage'] == True, 'probe'].tolist()
elite = pd.read_csv(RES + r'\elite_ranked_v6.csv').set_index('probe')
def gshort(p):
    g = elite.loc[p, 'gene']
    return g.split(';')[0] if isinstance(g, str) else '(unannotated)'

s7 = [['Probe', 'Gene', 'P95 noob', 'P95 GMQN', 'k/n noob', 'k/n GMQN', 'Gate-1 under single pipeline']]
for p in tri:
    r = v931[v931.probe == p].iloc[0]
    s7.append([p, gshort(p), f'{r.p95_noob:.4f}', f'{r.p95_gmqn:.4f}',
               f'{int(r.k_noob)}/17', f'{int(r.k_gmqn)}/17',
               'pass' if r.pass_gmqn else 'FAIL (marginal)'])

met = pd.read_csv(RES + r'\v932_boll_rule_metrics.csv')
s8 = [['Region-level rule', 'Tumour, all (n=431)', 'Tumour, premenopausal (n=28)',
       'Benign surgical (n=347)', 'GSE87621 benign (n=5)', 'GSE73948 endometriosis (n=26)', 'GSE73949 healthy (n=17)']]
label = {'R1_anchor_only': 'R1 anchor only', 'R2_partner_only': 'R2 partner only',
         'R3_OR': 'R3 OR (either probe)', 'R4_AND': 'R4 AND (both probes)',
         'R5_region_mean': 'R5 region mean'}
for _, r in met.iterrows():
    s8.append([label[r['rule']], f"{r['tumour_all']:.3f}", f"{r['tumour_pre']:.3f}",
               f"{r['benign347_GSE223817']:.3f}", f"{r['GSE87621_ctrl5_450K']:.3f}",
               f"{r['GSE73948_endo26_450K']:.3f}", f"{r['GSE73949_healthy17_GMQN']:.3f}"])

sel = pd.read_csv(RES + r'\v932_boll_rule_selection.csv')
s9 = [['Region-level rule', 'Bootstrap, criterion A', 'Bootstrap, criterion B', 'Benign LOOCV, criterion B']]
for _, r in sel.iterrows():
    s9.append([label[r['rule']], f"{r['boot_critA_freq']:.3f}", f"{r['boot_critB_freq']:.3f}",
               f"{r['loocv_critB_freq']:.3f}"])

parts = []

# ---------------- S2.12 ----------------
parts.append(heading('S2.12 Cross-pipeline sensitivity of the healthy-endometrium gate (noob versus GMQN)'))
parts.append(para('Rationale and methods. ',
    'The background layer of the funnel was assembled from two normalisation pipelines: gate 1 (healthy '
    'premenopausal endometrium, GSE73949, n=17) used minfi-noob betas reprocessed from IDATs, whereas '
    'gates 2-4 used EWAS Data Hub GMQN betas. To quantify how much of the gate-1 verdict is pipeline-bound, '
    'gate 1 was recomputed for all 476 elite probes from the EWAS GMQN per-sample betas of the same 17 '
    'donors (sample-level identity verified by GSM accession; noob recomputation reproduced the audit '
    'values with max |difference| 0.0005). Pass/fail uses the same rule as the main text (P95 ' + B + ' ' + LE + ' 0.15).'))
parts.append(para('Results. ',
    'GMQN shifted the healthy-endometrium P95 upward for a minority of probes: 24 of 476 probes (5.0%) '
    'changed pass/fail status, and every flip was in the conservative direction (pass under noob, fail '
    'under GMQN; zero flips in the opposite direction; pass/pass 41, fail/fail 411, Cohen\'s ' + KAPPA + ' = 0.747). '
    'The healthy-pass set therefore shrinks from 65 (mixed pipeline) to 41 (single GMQN pipeline), with no '
    'newly passing probes. Among the 15 finalists, five are lost under the single pipeline, all by marginal '
    'exceedances (GMQN P95 0.156-0.168 versus the 0.15 gate): cg23180938 (CDO1), cg07495363 (BOLL partner), '
    'cg25340966 (TBX15), cg08189989 and cg10812708 (Table S7). Both leads are pipeline-robust: cg27577527 '
    '(ZSCAN12; noob 0.060, GMQN 0.085) and cg24589459 (BOLL anchor; noob 0.110, GMQN 0.131) pass under '
    'either pipeline. Within the 23-probe experimental pool, 9 of 23 probes flip at gate 1 (all '
    'conservatively). Per-probe values for all 476 candidates: Supplementary Data S2i '
    '(Supplementary_Data_S2i_pipeline_sensitivity_gse73949.csv).'))
parts.append(caption('Table S7. The 15 finalists under the two gate-1 pipelines (k/n = samples exceeding '
    + B + ' 0.15). The single-pipeline column applies the GMQN gate used for all other background dimensions.'))
parts.append(table(s7))
parts.append(para('Interpretation. ',
    'The mixed-pipeline funnel is the more liberal of the two; the single-pipeline GMQN rerun is strictly '
    'more conservative and leaves the study\'s central claims intact: both leads and 10 of 15 finalists '
    'survive, and every loss is a marginal, threshold-adjacent exceedance rather than a qualitative '
    'reversal. CDO1\'s gate-1 pass, in addition to its benign-surgical failure, is pipeline-dependent, '
    'which further weakens its tier-B standing; the BOLL partner\'s marginal gate-1 failure is consistent '
    'with its contextual (non-anchor) role. For the prospective phase this analysis is moot by design: '
    'thresholds will be locked on a single quantitative assay (targeted bisulfite sequencing or qMSP), '
    'eliminating array-normalisation choice altogether.'))

# ---------------- S2.13 ----------------
parts.append(heading('S2.13 Pre-specification of the BOLL region-level combination rule'))
parts.append(para('Rationale and methods. ',
    'The two-probe OR concept for the BOLL region (anchor cg24589459 + partner cg07495363, 270 bp apart) '
    'was formulated after inspecting external data and is therefore post hoc. To replace it with a '
    'pre-specifiable rule, five candidate region-level rules were defined a priori and compared: R1 '
    'anchor only, R2 partner only, R3 OR (either probe positive), R4 AND (both positive), R5 region mean '
    '(mean of the two probe betas positive). Positivity uses the study-wide tumour-positivity threshold '
    '(' + B + ' > 0.3). Rules were scored on TCGA tumours (all, n=431; premenopausal, n=28) and on benign '
    'backgrounds (GSE223817 surgical controls, n=347; GSE87621 benign, n=5; GSE73948 eutopic endometriosis, '
    'n=26; GSE73949 healthy, n=17). Selection stability was assessed by stratified bootstrap (B=2000; '
    'resampling tumours, premenopausal tumours and benign controls independently) under two pre-declared '
    'criteria-(A) maximise premenopausal positivity subject to benign positivity ' + LE + ' 0.05, and '
    '(B) maximise Youden J = premenopausal positivity - benign positivity-and by leave-one-out over the '
    '347 benign controls under criterion B.'))
parts.append(caption('Table S8. Positivity (' + B + ' > 0.3) of the five candidate region-level rules in '
    'tumour and background cohorts.'))
parts.append(table(s8))
parts.append(caption('Table S9. Selection frequencies of the five candidate rules (B=2000 bootstrap; '
    '347-fold benign leave-one-out).'))
parts.append(table(s9))
parts.append(para('Results and decision. ',
    'The OR rule reached the highest tumour positivity (94.9% overall, 100% premenopausal) but was '
    'positive in 45.8% of benign surgical endometrium (159/347), driven entirely by the partner probe '
    '(161/347 partner-only-positive versus 0/347 anchor-only-positive); AND inherited the anchor\'s '
    'sensitivity without adding anything. The region-mean rule retained nearly all of the sensitivity '
    '(93.0% overall, 92.9% premenopausal) at 2.0% benign positivity and was stably selected: 65.1% of '
    'bootstrap replicates under both criteria (anchor-only 34.9%; OR, AND and partner-only never selected) '
    'and 100% of benign leave-one-out folds. The region-mean rule is therefore PRE-SPECIFIED as the '
    'region-level readout for the BOLL locus in the prospective phase; the OR rule is abandoned. The '
    'partner\'s benign elevation was not observed in the three smaller 450K background cohorts (0/5, '
    '0/26, 0/17; Note S2.11), so the elevation may be specific to the endometriosis-enriched surgical '
    'population; this does not change the decision, because that cohort is the largest and most '
    'stringent benign dimension. Per-rule values: Supplementary Data S2j '
    '(Supplementary_Data_S2j_boll_rule_metrics.csv); selection frequencies: Supplementary Data S2k '
    '(Supplementary_Data_S2k_boll_rule_selection.csv).'))

new_xml = ''.join(parts)

zin = zipfile.ZipFile(SRC)
xml = zin.read('word/document.xml').decode('utf-8')
i = xml.rfind('<w:sectPr')
assert i > 0, 'sectPr not found'
xml_new = xml[:i] + new_xml + xml[i:]

with zipfile.ZipFile(DST, 'w', zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == 'word/document.xml':
            data = xml_new.encode('utf-8')
        zout.writestr(item, data)
zin.close()

z = zipfile.ZipFile(DST)
assert z.testzip() is None
import xml.etree.ElementTree as ET
ET.fromstring(z.read('word/document.xml'))
x = z.read('word/document.xml').decode('utf-8')
print('S2.12 present:', 'S2.12 Cross-pipeline' in x)
print('S2.13 present:', 'S2.13 Pre-specification' in x)
print('tables now:', x.count('<w:tbl>'), '(was 7)')
print('Table S7/S8/S9:', all(f'Table S{n}' in x for n in [7, 8, 9]))
print('saved:', DST)

shutil.copy(RES + r'\v931_pipeline_sensitivity_gse73949.csv',
            OUTDIR + r'\Supplementary_Data_S2i_pipeline_sensitivity_gse73949.csv')
shutil.copy(RES + r'\v932_boll_rule_metrics.csv',
            OUTDIR + r'\Supplementary_Data_S2j_boll_rule_metrics.csv')
shutil.copy(RES + r'\v932_boll_rule_selection.csv',
            OUTDIR + r'\Supplementary_Data_S2k_boll_rule_selection.csv')
print('S2i/S2j/S2k CSVs written')
