# -*- coding: utf-8 -*-
"""Build the per-probe retention ledger for the 23-probe experimental pool
(v935): entry path, gate outcomes, evidence completeness, retention rationale,
intended use. Sources: audit S2a (dispositions), pool23 gate table (v928),
wetlab panel annotation. Output: results/v935_retention_ledger.csv"""
import numpy as np
import pandas as pd

BASE = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19'
RES = BASE + r'\results'

audit = pd.read_csv(BASE + r'\Supplementary_Data_S2a_single_probe_audit.csv').set_index('probe')
pool = pd.read_csv(RES + r'\wetlab_panel_23_primer_annotation.csv')
ci = pd.read_csv(RES + r'\v928_gate_proportion_ci_pool23.csv')

CO = {'g1': 'healthy17_GSE73949_noob', 'g2': 'benign347_GSE223817',
      'g3': 'cervix20_HPVneg_GSE46306', 'g4': 'blood60_GSE35069'}
piv = {}
for (pr, co), g in ci.groupby(['probe', 'cohort']):
    piv.setdefault(pr, {})[co] = float(g['p95'].iloc[0])

CORE = {'cg27577527', 'cg23180938', 'cg24589459', 'cg07495363'}
ROLE = {
    'cg27577527': 'Core set (ZSCAN12 anchor)',
    'cg23180938': 'Core set (CDO1; clinical reference)',
    'cg24589459': 'Core set (BOLL anchor)',
    'cg07495363': 'Core set (BOLL partner; regional readout candidate)',
}

def entry_path(p, disp, triage):
    if disp.startswith('Tier A lead'):
        return 'Strict funnel: Tier A lead'
    if disp.startswith('Tier B comparator'):
        if 'CDO1' in disp:
            return 'Strict funnel: Tier B comparator (clinical reference: approved-kit/CISENDO evidence)'
        return 'Strict funnel: Tier B comparator'
    if disp.startswith('BOLL partner'):
        return 'Strict funnel: BOLL regional characterisation (partner)'
    if disp.startswith('Region-coverage'):
        if triage is True:
            return 'Strict funnel: region-coverage (triage-passing)'
        return 'Region-coverage retention despite triage failure (region representation / incremental training-set coverage)'
    return disp

RATIONALE = {
    'cg27577527': 'Lowest background across assessable gates; discovery metrics at threshold but stable (rank 1 in 27/27 threshold-grid cells)',
    'cg24589459': 'Only probe passing all four gates by point estimate (cervical gate borderline, CI crosses 0.15)',
    'cg23180938': 'Kit anchor; lowest cervix/blood background among kit genes; fails benign-surgical dimension (endometriosis-enriched cohort)',
    'cg07495363': '100% premenopausal positivity; benign-surgical elevation; regional readout via pre-specified region-mean rule',
    'cg15790037': 'Favourable background in all assessed compartments; overall positivity near threshold',
    'cg16439198': 'Literature-supported; specimen-relevant background profile',
    'cg18507379': 'PCDHGA cluster representative; literature convergence',
    'cg10109500': 'Strong discovery metrics; benign-surgical background elevated (symptomatic cohort)',
}
rows = []
for _, r in pool.iterrows():
    p = r['probe']
    a = audit.loc[p]
    disp = a['final_disposition']
    triage = a['funnel_pass_triage']
    v = {k: piv.get(p, {}).get(c, np.nan) for k, c in CO.items()}
    n_as = sum(x == x for x in v.values())
    missing = [k for k, x in v.items() if x != x]
    compl = f'{n_as}/4' + (' (n/a: ' + ','.join(missing) + ')' if missing else '')
    rows.append(dict(
        probe=p, gene=str(r['gene']).split(';')[0] if isinstance(r['gene'], str) else '(unannotated)',
        chr=r['chr'], pos=r['pos'],
        entry_path=entry_path(p, disp, triage),
        audit_disposition=disp,
        g1_healthy17_p95=round(v['g1'], 4) if v['g1'] == v['g1'] else '',
        g2_benign347_p95=round(v['g2'], 4) if v['g2'] == v['g2'] else 'n/a',
        g3_cervix20_p95=round(v['g3'], 4) if v['g3'] == v['g3'] else '',
        g4_blood60_p95=round(v['g4'], 4) if v['g4'] == v['g4'] else '',
        evidence_completeness=compl,
        retention_rationale=RATIONALE.get(p, 'Region-coverage: ' + ('triage-passing' if triage is True else 'retained for region representation/incremental coverage')),
        intended_use=ROLE.get(p, 'Reserve (expansion conditional on core-set performance)'),
    ))
df = pd.DataFrame(rows)
df.to_csv(RES + r'\v935_retention_ledger.csv', index=False)
pd.set_option('display.width', 250)
print(df[['probe', 'gene', 'entry_path', 'evidence_completeness']].to_string(index=False))
print('\nentry path counts:')
print(df['entry_path'].value_counts())
print('\nsum check:', len(df))
