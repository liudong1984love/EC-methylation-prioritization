# -*- coding: utf-8 -*-
"""Manuscript v9.30 -> v9.31: consistency closure (reviewer round 3).
V1  Abstract Background: drop regulatory/consensus opening.
V2  Abstract Methods: drop DMR technical clause.
V3  Methods: compress historical score-constant detail (pointer to S2.3).
V4  Methods: donor correspondence VERIFIED (ArrayExpress titles; Data S2s).
V5  Results [45]: top-10 stability caveat + pointer to 48-53% retention.
V6  Results [46]: C4 v3 rewrite (479 before filters / 476 after; 2,942; 2,466; 17/19).
V7  Results [51]: region-rule selection criteria + cross-cohort caveat.
V8  Discussion [55]: BOLL expression numbers + softened mechanism statement.
V9  Discussion [60]: delete absolute 'devoid' sentence; keep qualified search.
V10 Limitations [62]: 'finalists 15 to 10' -> stage-named auxiliary set.
V11 [42]: 'purified-reference backgrounds' -> 'reference-fraction backgrounds'.
V12 [19]: blood fractions donor identities verified note.
"""
import zipfile, re
import xml.etree.ElementTree as ET

SRC = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Manuscript_EC_methylation_panel_v9.30.docx'
DST = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Manuscript_EC_methylation_panel_v9.31.docx'

B = '\u03b2'; LE = '\u2264'; GE = '\u2265'; MD = '\u2014'; DL = '\u0394'; RHO = '\u03c1'; MINUS = '\u2212'

def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

zin = zipfile.ZipFile(SRC)
xml = zin.read('word/document.xml').decode('utf-8')

def para_rebuild(xml, marker, new_runs, tag):
    for m in re.finditer(r'<w:p\b[^>]*>.*?</w:p>', xml, re.S):
        block = m.group(0)
        if marker not in block:
            continue
        ppr = re.search(r'<w:pPr>.*?</w:pPr>', block, re.S)
        ppr_s = ppr.group(0) if ppr else ''
        runs = ''.join(
            ('<w:r><w:rPr><w:b w:val="1"/></w:rPr>' if bold else '<w:r>') +
            '<w:t xml:space="preserve">' + esc(text) + '</w:t></w:r>'
            for text, bold in new_runs)
        return xml[:m.start()] + '<w:p>' + ppr_s + runs + '</w:p>' + xml[m.end():]
    raise AssertionError(tag + ': paragraph not found')

REPL = []

# ---- V1: abstract Background ----
REPL.append(('V1',
  'is maturing rapidly, with an NMPA-approved three-gene kit, international programmes and a '
  '2026 Chinese expert consensus. However, existing markers were developed',
  'is maturing rapidly, but existing markers were developed'))

# ---- V2: abstract Methods ----
REPL.append(('V2',
  'with single-probe differential-methylation filtering as the candidate-generating track; a '
  'parallel DMR analysis of TCGA and GSE67116 provided regional context but contributed no '
  'candidate to the final pool. Candidates were triaged',
  'with single-probe differential-methylation filtering as the candidate-generating track. '
  'Candidates were triaged'))

# ---- V3: Methods score constants compressed ----
REPL.append(('V3',
  'The historical implementation appended two further terms (cross-cohort consistency, 0.15; '
  'amplicon suitability, 0.05) that were constant across single-probe-track candidates at this '
  'stage and therefore cannot affect ranking; they are documented in Note S2.3 and are not '
  'counted as evidence here.',
  'Two further historical score terms were constant across single-probe-track candidates at '
  'this stage and therefore cannot affect ranking (Note S2.3).'))

# ---- V4: donor correspondence verified (first text node only) ----
REPL.append(('V4',
  'with donor correspondence inferred from within-fraction sample order and not independently '
  'verified' + MD + 'yielded similar intervals and changed no ',
  'with donor correspondence verified against ArrayExpress sample metadata (E-GEOD-35069 '
  'sample titles encode donor identifiers; verified mapping in Supplementary Data S2s)' + MD +
  'yielded intervals identical to the earlier order-inferred analysis and changed no '))

# ---- V5: top-10 caveat ----
REPL.append(('V5',
  'cannot support a top-10 comparison; Supplementary Data)',
  'cannot support a top-10 comparison; because the score formula is fixed in this analysis, '
  'unchanged rankings partly follow from the analysis definition rather than evidencing '
  'stability of the gated selection itself' + MD + 'the direct stability measure for the '
  'gated pipeline is the conditional-resampling retention of the two leads (48' + MD + '53%; '
  'Limitations); Supplementary Data)'))

# ---- V8: BOLL expression (anchor mid-paragraph; leading B in prior run) ----
REPL.append(('V8',
  'OLL expression is essentially absent from endometrial tissue irrespective of methylation '
  'state (172 TCGA tumours with paired RNA-seq; full analysis in Supplementary Data), so BOLL '
  'hypermethylation in EC is best interpreted as ectopic hypermethylation of a repressed '
  'germline-specific locus, not silencing of an expressed gene.',
  'OLL expression was essentially absent from endometrial tissue irrespective of methylation '
  'state (median log2(x+1) expression 0 in 172 TCGA tumours with paired RNA-seq and in 24 '
  'adjacent-normal samples; expression above 1 in 16% of tumours versus 0% of adjacent '
  'normals; methylation' + MD + 'expression correlation weakly negative, Spearman ' + RHO +
  ' ' + MINUS + '0.19 for cg24589459; full analysis in Supplementary Data). Low BOLL '
  'expression in both tumour and adjacent endometrial tissue therefore does not support '
  'transcriptional silencing of an actively expressed gene as the explanation for the '
  'observed methylation signal.'))

# ---- V9: delete absolute 'devoid' claim ----
REPL.append(('V9',
  'Public data remain devoid of premenopausal benign AUB methylation profiles (polyps, '
  'leiomyomas, adenomyosis), a gap only purpose-built cohorts can fill. A systematic search '
  'of the EWAS Data Hub (18 endometrium projects, 1,746 samples) confirmed this absence '
  'within the resources searched: we identified no eligible publicly available 450K/EPIC '
  'cohort providing premenopausal benign-lesion (polyp, leiomyoma, adenomyosis or PCOS) '
  'endometrial methylation profiles with the required benign-diagnosis annotations '
  '(databases, search date and eligibility criteria in Supplementary Note S2.11).',
  'A systematic search of the EWAS Data Hub (18 endometrium projects, 1,746 samples) '
  'identified no eligible publicly available 450K/EPIC cohort providing premenopausal '
  'benign-lesion (polyp, leiomyoma, adenomyosis or PCOS) endometrial methylation profiles '
  'with the required benign-diagnosis annotations in the resources searched (databases, '
  'search date and eligibility criteria in Supplementary Note S2.11)' + MD + 'a gap only '
  'purpose-built cohorts can fill.'))

# ---- V10: finalists -> stage-named auxiliary set ----
REPL.append(('V10',
  'reducing the healthy-pass set from 65 to 41 and the finalists from 15 to 10; all five '
  'finalist losses were marginal',
  'reducing the healthy-endometrium-pass set from 65 to 41 and the stage-2 triage-passing '
  'set' + MD + 'an auxiliary stage upstream of the strict all-four-gate Tier-A designation, '
  'not the Tier-A survivors' + MD + 'from 15 to 10; all five losses were marginal'))

# ---- V11: purified-reference -> reference-fraction ----
REPL.append(('V11',
  'these profiles reflect purified-reference backgrounds',
  'these profiles reflect reference-fraction backgrounds'))

# ---- V12: blood fractions donor identities verified ----
REPL.append(('V12',
  'whole blood, PBMC and granulocytes' + MD + 'each from the same six male donors) [21]',
  'whole blood, PBMC and granulocytes' + MD + 'each from the same six male donors, donor '
  'identities verified from ArrayExpress sample titles) [21]'))

for tag, old, new in REPL:
    n = xml.count(old)
    assert n == 1, '%s anchor count = %d' % (tag, n)
    xml = xml.replace(old, new, 1)
    print(tag, 'replaced')

# ---- V6: E6 paragraph rewritten with C4 v3 numbers ----
xml = para_rebuild(xml, 'Discovery-threshold sensitivity (exploratory).',
  [('Discovery-threshold sensitivity (exploratory). ', True),
   ('Because the leads sit close to the discovery cut-offs, the cut-offs themselves were '
    'scanned over a pre-specified, one-sided 27-combination grid loosening ' + DL + B +
    ' (' + GE + '0.50 to ' + GE + '0.30), overall tumour positivity (' + GE + '90% to ' +
    GE + '80%) and premenopausal positivity (' + GE + '85% to ' + GE + '75%), with FDR, the '
    'adjacent-normal P95 criterion and the audited cross-cohort-consistency and annotation '
    'filters held fixed and applied identically to every grid cell (Supplementary Note '
    'S2.15). At the original cut-offs the genome-wide re-derivation yields 479 candidates '
    'before these filters and 476 after' + MD + 'exactly the main-analysis audit '
    'universe' + MD + 'with the same two Tier-A survivors: ZSCAN12 rank 1/2 and BOLL anchor '
    'rank 2/2 by the effective score. Relaxation enlarged the filtered candidate set to up '
    'to 2,942 probes, but the four-gate background triage eliminated nearly all entrants: '
    'under the original ' + GE + '90% overall-positivity requirement no new low-background '
    'candidate entered in any grid cell, and at the loosest corner only 17 of 2,466 new '
    'entrants passed triage (19 survivors in total). cg27577527 (ZSCAN12) ranked first '
    'among triage survivors in all 27 combinations; cg24589459 (BOLL anchor) survived '
    'triage in all combinations (rank 2' + MD + '6 by the effective score, the variation '
    'reflecting higher-scoring low-background entrants at looser thresholds rather than '
    'signal loss). Approved-kit and reference markers showed distinct exclusion modes: '
    'CDO1 passes all discovery criteria and is eliminated only at background triage; '
    'AJAP1 is excluded by the overall-positivity criterion (84.2%); GALR1 and GYPC enter '
    'at no grid cell. Single-probe positivity below 90% therefore does not imply absence '
    'of complementary value' + MD + 'but relaxing the positivity requirement does not '
    'rescue the kit comparators either (Supplementary Data S2l, S2m; cross-analysis '
    'reconciliation in Supplementary Note S2.17, Table S14).', False)], 'V6')
print('V6 E6 rewritten')

# ---- V7: region-rule paragraph with selection criteria ----
xml = para_rebuild(xml, 'Region-level rules for BOLL (exploratory).',
  [('Region-level rules for BOLL (exploratory). ', True),
   ('Comparing the five pre-declared candidate rules (Methods; Supplementary Note S2.13) '
    'under two pre-declared selection criteria' + MD + 'maximising premenopausal tumour '
    'positivity subject to benign positivity ' + LE + '5% (fallback: minimum benign '
    'positivity), and maximising the Youden index (premenopausal minus benign '
    'positivity)' + MD + 'the region-mean rule was the most frequently selected under '
    'both criteria (65.1% of bootstrap replicates and 100% of benign leave-one-out folds; '
    'anchor-only 34.9%; OR, AND and partner-only never selected), retaining 92.9% '
    'premenopausal tumour positivity with 2.0% benign surgical positivity, versus 45.8% '
    'under the OR rule (159/347, partner-driven). Tumour and benign readouts derive from '
    'different cohorts (TCGA versus GSE223817), so resampling cannot remove '
    'cohort' + MD + 'disease confounding: these figures describe rule behaviour within '
    'these data, measure internal rule-selection stability rather than independent '
    'validation, and do not by themselves establish cross-cohort transportability; the '
    'rule and its threshold remain to be evaluated and locked in the independent '
    'development cohort.', False)], 'V7')
print('V7 region rules rewritten')

with zipfile.ZipFile(DST, 'w', zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == 'word/document.xml':
            data = xml.encode('utf-8')
        zout.writestr(item, data)
zin.close()

z = zipfile.ZipFile(DST)
assert z.testzip() is None
ET.fromstring(z.read('word/document.xml'))
x = z.read('word/document.xml').decode('utf-8')
checks = ['is maturing rapidly, but existing markers',
          'candidate-generating track. Candidates were triaged',
          'Two further historical score terms',
          'verified against ArrayExpress sample metadata',
          'partly follow from the analysis definition',
          '479 candidates before these filters and 476 after',
          'only 17 of 2,466 new entrants',
          'up to 2,942 probes',
          'maximising the Youden index',
          'cross-cohort transportability',
          'in 24 adjacent-normal samples',
          'does not support transcriptional silencing of an actively expressed gene',
          'A systematic search of the EWAS Data Hub (18 endometrium projects, 1,746 samples) identified',
          'stage-2 triage-passing set',
          'reference-fraction backgrounds',
          'donor identities verified from ArrayExpress sample titles']
for k in checks:
    print(k[:60], '->', k in x)
print('residual NMPA-approved three-gene kit, international:', x.count('international programmes'))
print('residual devoid:', x.count('devoid'))
print('residual finalists from 15 to 10:', x.count('finalists from 15 to 10'))
print('residual ectopic:', x.count('ectopic'))
print('residual inferred from within-fraction sample order:', x.count('inferred from within-fraction sample order'))
print('residual 2,983:', x.count('2,983'), '| residual 2,507:', x.count('2,507'))
print('residual purified-reference:', x.count('purified-reference'))
print('saved:', DST)
