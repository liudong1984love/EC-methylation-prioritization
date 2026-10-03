# -*- coding: utf-8 -*-
"""Supplementary v9.28 -> v9.29:
1. Fix S2.13 resampling wording (disjoint groups; overall-tumour draw not used).
2. Append S2.14 retention ledger & terminology (Table S10).
3. Append S2.15 discovery-threshold sensitivity (Tables S11, S12).
4. Append S2.16 sampling-site parallel gating (Table S13).
5. Copy Supplementary Data S2l-S2q CSVs.
"""
import zipfile, shutil
import pandas as pd

SRC = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Supplementary_Notes_EC_methylation_panel_v9.28.docx'
DST = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Supplementary_Notes_EC_methylation_panel_v9.29.docx'
RES = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\results'
OUTDIR = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19'

B = '\u03b2'; LE = '\u2264'; GE = '\u2265'; ARR = '\u2192'; MD = '\u2014'

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

def table(rows, fsz=16):
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
                run = (f'<w:r><w:rPr><w:b w:val="1"/><w:sz w:val="{fsz}"/></w:rPr>'
                       f'<w:t xml:space="preserve">{esc(str(cell))}</w:t></w:r>')
            else:
                tcpr = '<w:tcPr/>'
                run = (f'<w:r><w:rPr><w:sz w:val="{fsz}"/></w:rPr>'
                       f'<w:t xml:space="preserve">{esc(str(cell))}</w:t></w:r>')
            out.append(f'<w:tc>{tcpr}<w:p>{run}</w:p></w:tc>')
        out.append('</w:tr>')
    out.append('</w:tbl>')
    out.append('<w:p/>')
    return ''.join(out)

# ---------- load data ----------
ledger = pd.read_csv(RES + r'\v935_retention_ledger.csv')
grid = pd.read_csv(RES + r'\v933_threshold_grid_summary.csv')
known = pd.read_csv(RES + r'\v933_known_markers_discovery.csv')
paths = pd.read_csv(RES + r'\v934_sampling_paths_pool23.csv')

def na(v):
    return 'n/a' if v != v or v == '' else v

# Table S10: ledger
s10 = [['Probe', 'Gene', 'Entry path', 'G1 healthy', 'G2 benign-surg', 'G3 cervix', 'G4 blood',
        'Compl.', 'Intended use']]
for _, r in ledger.iterrows():
    s10.append([r['probe'], r['gene'], r['entry_path'],
                na(r['g1_healthy17_p95']), na(r['g2_benign347_p95']),
                na(r['g3_cervix20_p95']), na(r['g4_blood60_p95']),
                r['evidence_completeness'].split(' ')[0], r['intended_use']])

# Table S11: threshold grid
s11 = [['\u0394' + B + ' ' + GE, 'Overall pos ' + GE, 'Premeno pos ' + GE, 'Candidates', 'New vs 476',
        'Triage pass', 'New low-bg', 'ZSCAN12 rank', 'BOLL anchor rank']]
for _, r in grid.iterrows():
    def rk(s):
        # 'in-set, triage-PASS, rank 11/20' -> '11/20'
        if 'rank' in s:
            return s.split('rank ')[1]
        return 'not in set'
    s11.append([f"{r['delta_beta']:.2f}", f"{r['case_pos']:.2f}", f"{r['pre_pos']:.2f}",
                int(r['n_candidates']), int(r['n_new_vs_476']), int(r['n_triage_pass']),
                int(r['n_new_lowbg']), rk(r['ZSCAN12_anchor']), rk(r['BOLL_anchor'])])

# Table S12: known markers
s12 = [['Marker', '\u0394' + B, 'Overall pos', 'Premeno pos', 'Ctrl P95', 'In 476',
        'Strict-criterion failure', 'Enters loosest grid']]
for _, r in known.iterrows():
    s12.append([r['marker'], f"{r['delta_beta']:.3f}", f"{r['case_pos03']:.3f}", f"{r['pre_pos03']:.3f}",
                f"{r['ctrl_p95']:.3f}", 'yes' if r['in_strict_476'] else 'no',
                r['strict_fail'], 'yes' if r['enters_at_loosest_grid'] else 'no'])

# Table S13: pool23 path matrix
s13 = [['Probe', 'Gene', 'G1 healthy', 'G2 benign-surg', 'G3 cervix', 'G4 blood', 'Compl.',
        'Strict 4-gate', 'Endo-brush', 'Cervix/self']]
def tf(v):
    return 'pass' if v else 'FAIL'
for _, r in paths.iterrows():
    s13.append([r['probe'], r['gene'], na(round(r['g1_healthy17'], 3)) if r['g1_healthy17'] == r['g1_healthy17'] else 'n/a',
                na(round(r['g2_benign347'], 3)) if r['g2_benign347'] == r['g2_benign347'] else 'n/a',
                na(round(r['g3_cervix20'], 3)) if r['g3_cervix20'] == r['g3_cervix20'] else 'n/a',
                na(round(r['g4_blood60'], 3)) if r['g4_blood60'] == r['g4_blood60'] else 'n/a',
                f"{int(r['n_gates_assessable'])}/4", tf(r['STRICT4']), tf(r['ENDO_brush']), tf(r['CERV_selfsample'])])

parts = []

# ---------------- S2.14 ----------------
parts.append(heading('S2.14 Per-probe retention ledger and terminology for the 23-probe experimental pool'))
parts.append(para('Purpose. ',
    'The 23-probe experimental pool is NOT the product of a single strict funnel: it is a heterogeneous, '
    'deliberately assembled library combining three retention paths. This note gives the per-probe ledger '
    '(Table S10; machine-readable: Supplementary Data S2p, Supplementary_Data_S2p_retention_ledger.csv) and '
    'fixes the terminology used across the main text and these Notes.'))
parts.append(para('Accounting. ',
    '476 elite single probes ' + ARR + ' 65 pass the healthy-endometrium screen ' + ARR + ' 15 pass the '
    'four-gate triage ("finalists") ' + ARR + ' 2 Tier A leads with favourable background across all '
    'assessable gates. The 23-probe pool = 7 tier-A/B leads (2 Tier A + 5 Tier B comparators) + 1 BOLL '
    'partner probe (regional characterisation) + 15 region-coverage probes (7 triage-passing; 8 retained '
    'despite triage failure for region representation or incremental training-set coverage). Under the '
    'single-pipeline GMQN rerun of gate 1 (Note S2.12), the healthy-pass set shrinks 65 ' + ARR + ' 41 and '
    'the finalists 15 ' + ARR + ' 10; this is a pipeline variant of the same set, not a different set, and '
    'both leads pass under either pipeline. The strict "two leads" statement and the "15 finalists" '
    'statement therefore refer to different, explicitly defined stages of one funnel.'))
parts.append(para('Evidence completeness. ',
    'The two leads are NOT evidence-equivalent. cg24589459 (BOLL anchor) is fully assessed (4/4 gates) with '
    'a point-estimate pass whose cervical bootstrap CI crosses the 0.15 gate (borderline). cg27577527 '
    '(ZSCAN12 anchor) is assessed at 3/3 probe-assessable gates with the decisive benign-surgical '
    'dimension pending (probe absent from EPIC; 450K support in Note S2.11 is point-estimate-level and is '
    'not treated as a gate pass). On the three commonly assessable dimensions (healthy endometrium, '
    'cervix, blood), ZSCAN12 shows P95 ' + B + ' 0.044' + MD + '0.060 versus BOLL 0.110' + MD + '0.141; '
    'both pass, and ZSCAN12 ranks first on the effective score in all 27 discovery-threshold grid cells '
    '(Note S2.15). Total scores are therefore compared only on commonly assessable dimensions.'))
parts.append(caption('Table S10. Retention ledger of the 23 experimental-pool probes: entry path, gate P95 '
    + B + ' values (G1 healthy premenopausal endometrium n=17; G2 benign surgical endometrium n=347; '
    'G3 HPV-negative cervix n=20; G4 blood worst-fraction n=60), evidence completeness and intended use. '
    'n/a = dimension not assessable for the probe (missing, not failed).'))
parts.append(table(s10))

# ---------------- S2.15 ----------------
parts.append(heading('S2.15 Discovery-threshold sensitivity analysis (exploratory)'))
parts.append(para('Rationale and methods. ',
    'Single-probe discovery required FDR<0.05, ' + '\u0394' + B + ' ' + GE + ' 0.50, overall tumour '
    'positivity ' + GE + ' 90% and premenopausal positivity ' + GE + ' 85% (adjacent-normal P95 '
    + LE + ' 0.40 throughout). Because both leads sit close to these cut-offs, the cut-offs themselves '
    'were scanned over a pre-specified, one-sided (looser) 27-combination grid: ' + '\u0394' + B + ' '
    + GE + ' {0.30, 0.40, 0.50}, overall positivity ' + GE + ' {80%, 85%, 90%}, premenopausal positivity '
    + GE + ' {75%, 80%, 85%}, with FDR and the adjacent-normal criterion held fixed. The stricter '
    'direction is uninformative by construction (the leads sit at the boundary). Thresholds were not '
    'chosen to protect the current leads. At each grid cell, the candidate set was re-derived genome-wide '
    'from the TCGA discovery statistics (485,512 probes), new entrants versus the strict 476 were '
    'identified, and all new entrants were evaluated against the same four background gates (union rule, '
    'n/a excluded; per-probe gate values in Supplementary Data S2q). Reproduction check: at the original '
    'criteria the genome-wide re-derivation returns 479 candidates, 476 of which are the audited elite '
    'set (the three additional probes map to chrX/chrY or a PCDHB-region shore and were excluded by the '
    'audited annotation filters).'))
parts.append(caption('Table S11. Discovery-threshold grid (27 combinations). Candidates = genome-wide '
    'probes meeting the criteria; New vs 476 = candidates outside the strict elite set; Triage pass = '
    'candidates passing all assessable background gates; New low-bg = triage-passing new entrants. Ranks '
    'are by the three-term effective score among triage survivors.'))
parts.append(table(s11))
parts.append(para('Results. ',
    'Relaxation enlarged the candidate set from 479 to 2,983 probes, but background triage eliminated '
    'nearly all entrants: at the original ' + GE + '90% overall-positivity requirement no new '
    'low-background candidate entered in any grid cell, and at the loosest corner only 17 of 2,507 new '
    'entrants passed triage. cg27577527 (ZSCAN12) ranked first among triage survivors in all 27 '
    'combinations; cg24589459 (BOLL anchor) survived triage in all combinations (rank 3' + MD + '11), its '
    'rank sliding at looser thresholds because lower-background entrants out-score its borderline '
    'cervical dimension, not because its signal is lost. Approved-kit and reference markers show distinct '
    'exclusion modes (Table S12): CDO1 passes every discovery criterion and is eliminated only at '
    'background triage; AJAP1 is excluded by the overall-positivity criterion (84.2% versus ' + GE +
    '90%); GALR1 (overall positivity 63' + MD + '70%) and GYPC (' + '\u0394' + B + ' 0.45, '
    'adjacent-normal P95 0.52) enter at no grid cell. Single-probe positivity below 90% therefore does '
    'not imply absence of complementary value ' + MD + ' but relaxing the positivity requirement does '
    'not rescue the kit comparators either. This analysis is exploratory and does not change the main '
    'funnel.'))
parts.append(caption('Table S12. Discovery-metric decomposition for approved-kit and reference markers '
    '(positivity at ' + B + '>0.3; ctrl P95 = adjacent-normal P95).'))
parts.append(table(s12))
parts.append(para('Data. ',
    'Grid summary: Supplementary Data S2l (Supplementary_Data_S2l_threshold_grid_summary.csv); marker '
    'decomposition: Supplementary Data S2m (Supplementary_Data_S2m_known_markers_discovery.csv); '
    'per-probe background gate values for all 2,507 grid entrants: Supplementary Data S2q '
    '(Supplementary_Data_S2q_new_entrants_background.csv).'))

# ---------------- S2.16 ----------------
parts.append(heading('S2.16 Sampling-site parallel gating analysis (exploratory)'))
parts.append(para('Rationale and methods. ',
    'The strict four-gate union selects cross-compartment low-background markers, which is more '
    'stringent than any single intended application requires. Without altering the main rule, two '
    'application-specific gating paths were evaluated on the 476 elite probes and on the 23-probe pool '
    '(same gate values and union logic as the main analysis; blood reported separately in the '
    'endometrial path as bleeding-contamination risk): (i) ENDO ' + MD + ' physician-collected '
    'endometrial brush: gates = healthy premenopausal endometrium (GSE73949) + benign surgical '
    'endometrium (GSE223817); (ii) CERV ' + MD + ' cervical scraping / self-sampling: gates = '
    'HPV-negative cervix (GSE46306) + blood fractions (GSE35069).'))
parts.append(para('Results. ',
    'Among the 476 elite probes, the strict union retains 2, the ENDO path 11 (3 if blood is also '
    'required), and the CERV path 43; 52 probes pass at least one path. Within the 23-probe pool '
    '(Table S13), the two leads pass every path; 15 pool probes ' + MD + ' including CDO1, the BOLL '
    'partner, GHSR and ARL5C ' + MD + ' fail the benign-surgical dimension (P95 0.25' + MD + '0.56) but '
    'pass the cervical path, whereas six pool probes show the mirror pattern (endometrium-passing, '
    'cervix/blood-elevated). CDO1 therefore illustrates the value of the parallel analysis: eliminated '
    'by the benign-surgical dimension, it passes the cervical path (0.055/0.080), consistent with its '
    'established cervical application ' + MD + ' while its endometrial-brush suitability remains '
    'unresolved, because the endometriosis-enriched surgical cohort likely overestimates an '
    'AUB-screening background. Endometrial-path retentions are thus conservative, and cervical-path '
    'passes do not establish endometrial suitability. Path assignments guide specimen-matched assay '
    'development only; the strict four-gate rule remains the main analysis.'))
parts.append(caption('Table S13. The 23-probe experimental pool under the three gating rules (P95 ' + B +
    ' per dimension; pass = P95 ' + LE + ' 0.15 on all assessable dimensions of the rule).'))
parts.append(table(s13))
parts.append(para('Data. ',
    'Pool matrix: Supplementary Data S2n (Supplementary_Data_S2n_sampling_paths_pool23.csv); per-probe '
    'path results for all 476 elite probes: Supplementary Data S2o '
    '(Supplementary_Data_S2o_sampling_paths_476.csv).'))

new_xml = ''.join(parts)

zin = zipfile.ZipFile(SRC)
xml = zin.read('word/document.xml').decode('utf-8')

# ---- S2.13 wording fix (disjoint-group resampling; overall-tumour draw not a selection input) ----
old13 = ('Selection stability was assessed by stratified bootstrap (B=2000; '
         'resampling tumours, premenopausal tumours and benign controls independently) under two '
         'pre-declared criteria')
new13 = ('Selection stability was assessed by bootstrap (B=2000; premenopausal tumours and benign '
         'controls resampled independently of each other, the two groups being disjoint cohorts; '
         'overall-tumour positivity was not an input to rule selection) under two pre-declared criteria')
assert xml.count(old13) == 1, 'S2.13 anchor: %d' % xml.count(old13)
xml = xml.replace(old13, new13, 1)

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
for k in ['S2.14 Per-probe retention ledger', 'S2.15 Discovery-threshold', 'S2.16 Sampling-site',
          'Table S10', 'Table S11', 'Table S12', 'Table S13', 'S2l', 'S2m', 'S2n', 'S2o', 'S2p', 'S2q']:
    print(k, '->', k in x)
print('saved:', DST)

copies = [
    ('v933_threshold_grid_summary.csv', 'Supplementary_Data_S2l_threshold_grid_summary.csv'),
    ('v933_known_markers_discovery.csv', 'Supplementary_Data_S2m_known_markers_discovery.csv'),
    ('v934_sampling_paths_pool23.csv', 'Supplementary_Data_S2n_sampling_paths_pool23.csv'),
    ('v934_sampling_paths_476.csv', 'Supplementary_Data_S2o_sampling_paths_476.csv'),
    ('v935_retention_ledger.csv', 'Supplementary_Data_S2p_retention_ledger.csv'),
    ('v933_new_entrants_background.csv', 'Supplementary_Data_S2q_new_entrants_background.csv'),
]
for src_f, dst_f in copies:
    shutil.copy(RES + '\\' + src_f, OUTDIR + '\\' + dst_f)
print('S2l-S2q CSVs written')
