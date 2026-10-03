# -*- coding: utf-8 -*-
"""Supplementary v9.24 -> v9.27: append S2.10 (proportion+CI) and S2.11 (extended benign search)."""
import zipfile, shutil, re, os

SRC = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Supplementary_Notes_EC_methylation_panel_v9.24.docx'
DST = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Supplementary_Notes_EC_methylation_panel_v9.27.docx'

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
    # spacer paragraph after table (Word requires para between consecutive tables/end)
    out.append('<w:p/>')
    return ''.join(out)

B = '\u03b2'      # beta
LE = '\u2264'     # <=
GE = '\u2265'     # >=

parts = []

# ---------------- S2.10 ----------------
parts.append(heading('S2.10 Proportion-based background assessment with Wilson confidence intervals'))
parts.append(para('Rationale. ',
    'On small reference cohorts a P95 point estimate is dominated by one or two extreme samples, and a '
    'point-estimate pass does not certify that the underlying exceedance rate sits below the gate. As a '
    'complement to the P95 gate audit, every cohort-by-probe pair among the 476 elite candidates was '
    'summarised by the proportion of samples exceeding ' + B + ' 0.15 (k/n) with a two-sided Wilson 95% '
    'confidence interval. Pairs that pass the gate by point estimate but whose confidence interval crosses '
    '0.15 are classified as uncertain and are not treated as certified passes until larger reference samples '
    'narrow the interval. The cervical analysis uses the 20 HPV-negative control scrapes of GSE46306, '
    'identical to the gate definition in the main text; the blood analysis uses the worst per-fraction value '
    'on the six independent donors. Per-probe values for all 476 candidates: Supplementary Data S2g '
    '(Supplementary_Data_S2g_background_proportion_wilson_ci.csv).'))
parts.append(caption('Table S3. Exceedance proportions k/n at ' + B + ' 0.15 (Wilson 95% upper bound in '
    'parentheses) for the lead and core-set probes across the four background dimensions.'))
parts.append(table([
    ['Probe', 'Healthy endometrium (n=17)', 'Benign surgical (n=347)', 'HPV-negative cervix (n=20)', 'Blood, worst fraction (n=6)'],
    ['cg24589459 (BOLL anchor)', '0/17 (0.184)', '5/347 (0.033)', '1/20 (0.236)', '0/6 (0.390)'],
    ['cg27577527 (ZSCAN12)', '0/17 (0.184)', 'n/a (absent from EPIC)', '0/20 (0.161)', '0/6 (0.390)'],
    ['cg07495363 (BOLL partner)', '1/17 (0.270)', '327/347 (fail)', '0/20 (0.161)', '0/6 (0.390)'],
    ['cg23180938 (CDO1, tier-B)', '1/17 (0.270)', '270/347 (fail)', '0/20 (0.161)', '0/6 (0.390)'],
]))
parts.append(para('Lead-probe profile. ',
    'cg24589459 (BOLL anchor) is the only probe whose benign-surgical pass is certified at interval level '
    '(5/347 exceedances, upper 95% bound 0.033), because n=347 is the only background cohort large enough to '
    'narrow the interval below the gate; its cervical gate remains borderline under both methods (1/20, upper '
    'bound 0.236; P95 0.141), consistent with the bootstrap P95 interval reported in the main text. '
    'cg27577527 (ZSCAN12) shows zero exceedances in every probe-assessable dimension (0/17 healthy, 0/20 '
    'cervical, 0/6 worst blood fraction), the cleanest point-estimate profile among all candidates; at these '
    'sample sizes the intervals are wide, so this is point-estimate-level evidence rather than a certified '
    'specificity bound. The benign-surgical dimension remains not assessable for this probe on EPIC '
    '(first 450K benign assessments in Note S2.11).'))
parts.append(caption('Table S4. Minimum reference-sample size n for the Wilson 95% upper bound to fall at or '
    'below 5% or 10%, by the observed number of exceedances k.'))
parts.append(table([
    ['Observed exceedances k', 'n for upper bound ' + LE + ' 5%', 'n for upper bound ' + LE + ' 10%'],
    ['0', '59', '29'], ['1', '93', '46'], ['2', '124', '61'], ['3', '153', '75'], ['5', '209', '102'],
]))
parts.append(para('Sample-size implications. ',
    'Even with zero observed exceedances, 59 reference samples are required before the interval excludes a '
    '5% exceedance rate. The healthy (n=17), cervical (n=20) and per-fraction blood (n=6) gates therefore '
    'certify nothing beyond their point estimates at current sample sizes; pilot cohorts of about ten samples '
    'per arm function as rule-out screens, and definitive benign-specificity estimation will require on the '
    'order of 60-100 annotated benign samples. These requirements are carried into the design of the planned '
    'premenopausal AUB cohort.'))

# ---------------- S2.11 ----------------
parts.append(heading('S2.11 Extended benign background search (EWAS Data Hub)'))
parts.append(para('Search. ',
    'To determine whether the benign-endometrium background could be strengthened with public data, the EWAS '
    'Data Hub was searched exhaustively for endometrial methylation datasets (tissue field "uterus - '
    'endometrium"; 18 projects, 1,746 samples, the platform\'s entire endometrium collection as of '
    '29 September 2026). After removing cohorts already used in this study and non-informative entries '
    '(FFPE sarcoma series, molar pregnancies, metastasis panels), five candidate cohorts (93 samples) were '
    'downloaded and evaluated (Table S5). No publicly available 450K/EPIC cohorts of endometrial polyps, '
    'leiomyomas, adenomyosis or PCOS exist within the reachable repositories; this systematic null result '
    'substantiates the statement in the main text that intended-use benign backgrounds must be collected '
    'prospectively.'))
parts.append(caption('Table S5. Disposition of the five newly identified endometrial methylation cohorts.'))
parts.append(table([
    ['Cohort', 'Platform', 'Composition', 'Characterisation', 'Disposition'],
    ['GSE87621', '450K', '5 controls + 4 endometriosis', 'benign tissue controls', 'retained (benign dimension)'],
    ['GSE73948', '450K', '26 endometriosis eutopic endometrium (secretory phase)',
     'same study series as the healthy reference GSE73949', 'retained (symptomatic background, not a benign control)'],
    ['GSE122148', 'EPIC', '8 tissue samples (5 normal-weight, 3 obese)',
     'probe-specific background elevations (Table S6)', 'sensitivity cohort only, not a gate'],
    ['GSE130028', '450K', '48 (4 donors x 4 hormone treatments)',
     'in-vitro hormone-treated cultured cells', 'excluded'],
    ['GSE108143', '450K', '2 (decidualized / non-decidualized)',
     'in-vitro decidualization of cultured stromal cells', 'excluded'],
]))
parts.append(caption('Table S6. Background P95 ' + B + ' (exceedances k/n at ' + B + ' 0.15) for the lead and '
    'core-set probes in the new cohorts, with the GSE223817 surgical-benign cohort for contrast.'))
parts.append(table([
    ['Probe', 'GSE87621 benign (n=5)', 'GSE73948 endometriosis (n=26)', 'GSE122148 obesity-enriched (n=8)', 'GSE223817 surgical benign (n=347)'],
    ['cg27577527 (ZSCAN12)', '0.056 (0/5)', '0.108 (1/26)', 'n/a (absent from EPIC)', 'n/a (absent from EPIC)'],
    ['cg24589459 (BOLL anchor)', '0.136 (0/5)', '0.139 (1/26)', '0.044 (0/8)', '0.100 (5/347)'],
    ['cg07495363 (BOLL partner)', '0.093 (0/5)', '0.127 (0/26)', '0.253 (8/8)', '0.427 (327/347)'],
    ['cg23180938 (CDO1, tier-B)', '0.033 (0/5)', '0.149 (2/26)', '0.217 (3/8)', '0.372 (270/347)'],
]))
parts.append(para('Findings. ',
    'First, cg27577527 (ZSCAN12), not assessable in the EPIC-based surgical cohort, showed low background in '
    'two independent 450K cohorts: 0/5 exceedances (P95 0.056) in benign endometrium and 1/26 (P95 0.108) in '
    'eutopic endometrium from endometriosis patients. These are the first public benign-tissue assessments of '
    'this probe; at these sample sizes they provide point-estimate support only, and we continue to refrain '
    'from claiming four-gate compliance. Second, the BOLL partner probe cg07495363 and the tier-B comparator '
    'cg23180938 (CDO1) were clean in both new 450K cohorts (P95 ' + LE + ' 0.149), in contrast to their high '
    'background in GSE223817 (P95 0.427 and 0.372); this supports the interpretation that the surgical-cohort '
    'elevation is a property of its endometriosis-enriched, symptomatic population rather than of benign '
    'endometrium in general, and is consistent with the grading of GSE223817 as the most stringent background '
    'dimension. Third, in the obesity-enriched EPIC cohort, cg07495363 (8/8) and cg23180938 (3/8) exceeded '
    'the gate while cg24589459 remained low (0/8); the pattern is probe-specific rather than a global shift, '
    'its origin (adiposity-related field effects, platform or cohort composition) cannot be resolved at n=8, '
    'and body-mass index should be recorded as a covariate in future background cohorts.'))
parts.append(para('Data. ',
    'Per-probe values for all 476 elite candidates across the new cohorts: Supplementary Data S2h '
    '(Supplementary_Data_S2h_extended_benign_background.csv). Accession-to-publication mapping for the new '
    'cohorts is pending because NCBI services were unreachable from the analysis environment; cohort '
    'characterisations above are taken from the curated EWAS Data Hub sample metadata.'))

new_xml = ''.join(parts)

# ---- splice into document ----
zin = zipfile.ZipFile(SRC)
xml = zin.read('word/document.xml').decode('utf-8')
i = xml.rfind('<w:sectPr')
assert i > 0, 'sectPr not found'
xml_new = xml[:i] + new_xml + xml[i:]

os.makedirs(os.path.dirname(DST), exist_ok=True)
with zipfile.ZipFile(DST, 'w', zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == 'word/document.xml':
            data = xml_new.encode('utf-8')
        zout.writestr(item, data)
zin.close()

# ---- verify ----
z = zipfile.ZipFile(DST)
assert z.testzip() is None
import xml.etree.ElementTree as ET
ET.fromstring(z.read('word/document.xml'))
x = z.read('word/document.xml').decode('utf-8')
print('S2.10 present:', 'S2.10 Proportion-based' in x)
print('S2.11 present:', 'S2.11 Extended benign' in x)
print('tables now:', x.count('<w:tbl>'), '(was 3)')
print('Table S3/S4/S5/S6:', all(f'Table S{n}' in x for n in [3,4,5,6]))
print('saved:', DST)

# ---- copy Supplementary Data CSVs ----
res = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\results'
outdir = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19'
shutil.copy(res + r'\v928_gate_proportion_ci_476.csv',
            outdir + r'\Supplementary_Data_S2g_background_proportion_wilson_ci.csv')
shutil.copy(res + r'\v930_expanded_background_476.csv',
            outdir + r'\Supplementary_Data_S2h_extended_benign_background.csv')
print('S2g/S2h CSVs written')
