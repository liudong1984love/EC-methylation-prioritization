# -*- coding: utf-8 -*-
"""Item 2: proportion + Wilson CI as primary background-gate criterion.
Recomputes per-cohort exceedance stats for the 476 elite probes using the
same beta sources as the audited gates, cross-checks P95 vs elite_ranked_v6.csv.
"""
import os, math, json
import pandas as pd
import numpy as np

BASE = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19'
RES = BASE + r'\results'
DATA = BASE + r'\data'
THR = 0.15
Z = 1.959964

elite = pd.read_csv(RES + r'\elite_ranked_v6.csv')
probes = elite['probe'].tolist()
pset = set(probes)
print('elite probes:', len(probes))

def wilson(k, n):
    if n == 0:
        return (float('nan'),) * 3
    p = k / n
    d = 1 + Z * Z / n
    c = p + Z * Z / (2 * n)
    h = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n))
    return p, max(0.0, (c - h) / d), min(1.0, (c + h) / d)

def stats_row(vals):
    v = np.asarray([x for x in vals if x == x], dtype=float)  # drop NaN
    n = len(v)
    k = int((v > THR).sum())
    p, lo, hi = wilson(k, n)
    return dict(n=n, k=k, prop=p, ci_lo=lo, ci_hi=hi,
                median=float(np.median(v)) if n else np.nan,
                p90=float(np.percentile(v, 90)) if n else np.nan,
                p95=float(np.percentile(v, 95)) if n else np.nan,
                max=float(v.max()) if n else np.nan)

def read_gsm_files(gse, gsm_keep=None):
    """Return dict probe -> {gsm: beta} restricted to pset."""
    d = os.path.join(DATA, gse, 'ewas_betas')
    out = {p: {} for p in probes}
    files = [f for f in os.listdir(d) if f.endswith('.txt')]
    n_done = 0
    for f in files:
        gsm = f.split('.')[0]
        if gsm_keep is not None and gsm not in gsm_keep:
            continue
        with open(os.path.join(d, f), encoding='utf-8', errors='ignore') as fh:
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) != 2:
                    continue
                pr, val = parts
                if pr in pset:
                    try:
                        out[pr][gsm] = float(val)
                    except ValueError:
                        pass
        n_done += 1
        if n_done % 100 == 0:
            print(f'  {gse}: {n_done} files...')
    return out

rows = []

# ---- Gate 1: healthy premenopausal endometrium, minfi-noob RDS (audit-consistent)
import pyreadr
r = pyreadr.read_r(RES + r'\GSE73949_betas.rds')
key = list(r.keys())[0]
b1 = r[key]
print('gate1 rds:', key, b1.shape)
b1 = b1.loc[b1.index.intersection(probes)]
for pr in probes:
    vals = b1.loc[pr].values if pr in b1.index else []
    s = stats_row(vals)
    rows.append(dict(probe=pr, cohort='healthy17_GSE73949_noob', group='all', **s))
print('gate1 done')

# ---- Gate 2: benign surgical endometrium, GSE223817 controls
meta = pd.read_csv(RES + r'\GSE223817_ewas_meta.csv')
ctrl = set(meta.loc[meta['sample type'] == 'control', 'sample id'])
print('gate2 controls:', len(ctrl))
b2 = read_gsm_files('GSE223817', ctrl)
for pr in probes:
    s = stats_row(list(b2[pr].values()))
    rows.append(dict(probe=pr, cohort='benign347_GSE223817', group='all', **s))
print('gate2 done')

# ---- Gate 3: cervix HPV-neg, GSE46306
b3 = read_gsm_files('GSE46306', None)
for pr in probes:
    s = stats_row(list(b3[pr].values()))
    rows.append(dict(probe=pr, cohort='cervix44_GSE46306', group='all', **s))
print('gate3 done')

# ---- Gate 4: blood fractions, GSE35069 (per tissue fraction + worst)
m4 = pd.read_csv(RES + r'\GSE35069_ewas_meta.csv')
frac = dict(zip(m4['sample id'], m4['tissue']))
b4 = read_gsm_files('GSE35069', None)
fractions = sorted(set(frac.values()))
print('gate4 fractions:', len(fractions))
for pr in probes:
    per = {}
    for gsm, val in b4[pr].items():
        per.setdefault(frac.get(gsm, '?'), []).append(val)
    worst_hi = -1.0
    for fr in fractions:
        s = stats_row(per.get(fr, []))
        rows.append(dict(probe=pr, cohort='blood60_GSE35069', group=fr, **s))
        if s['ci_hi'] == s['ci_hi']:
            worst_hi = max(worst_hi, s['ci_hi'])
    # fraction-worst summary row (mirrors current 'max over fractions' logic)
    p95s = [stats_row(per.get(fr, []))['p95'] for fr in fractions]
    ks = [stats_row(per.get(fr, []))['k'] for fr in fractions]
    rows.append(dict(probe=pr, cohort='blood60_GSE35069', group='WORST_FRACTION',
                     n=6, k=max(ks), prop=np.nan, ci_lo=np.nan, ci_hi=worst_hi,
                     median=np.nan, p90=np.nan, p95=max(p95s), max=np.nan))
print('gate4 done')

df = pd.DataFrame(rows)
df.to_csv(RES + r'\v928_gate_proportion_ci_476.csv', index=False)
print('saved v928_gate_proportion_ci_476.csv', df.shape)

# ---- cross-check P95 vs elite_v6 current gate columns
chk = elite.set_index('probe')
def cmp(col, cohort, group='all'):
    mine = df[(df.cohort == cohort) & (df.group == group)].set_index('probe')['p95']
    ref = chk[col]
    both = pd.concat([mine, ref], axis=1, keys=['mine', 'ref']).dropna()
    diff = (both['mine'] - both['ref']).abs()
    print(f'{cohort} vs {col}: n={len(both)} maxdiff={diff.max():.4f} meandiff={diff.mean():.4f}')
cmp('ne_bg_p95', 'healthy17_GSE73949_noob')
cmp('ne347_p95', 'benign347_GSE223817')
cmp('cervix_p95', 'cervix44_GSE46306')
cmp('blood_max_p95', 'blood60_GSE35069', 'WORST_FRACTION')

# ---- new-criterion summary: pass iff Wilson CI upper <= 0.05
CUT = 0.05
summ = []
for cohort, group in [('healthy17_GSE73949_noob', 'all'), ('benign347_GSE223817', 'all'),
                      ('cervix44_GSE46306', 'all'), ('blood60_GSE35069', 'WORST_FRACTION')]:
    sub = df[(df.cohort == cohort) & (df.group == group)].set_index('probe')
    pt_pass = sub['p95'] <= THR
    ci_pass = sub['ci_hi'] <= CUT
    summ.append(dict(cohort=cohort, n_samples=int(sub['n'].max()),
                     point_pass=int(pt_pass.sum()), ci_pass=int(ci_pass.sum()),
                     point_pass_ci_fail=int((pt_pass & ~ci_pass).sum())))
sdf = pd.DataFrame(summ)
print(sdf.to_string(index=False))
sdf.to_csv(RES + r'\v928_gate_proportion_ci_summary.csv', index=False)

# ---- flips among the 23-probe pool and the 15/65 sets
pool = pd.read_csv(RES + r'\wetlab_panel_23_primer_annotation.csv')['probe'].tolist()
pool_rows = df[df.probe.isin(pool) & (df.group.isin(['all', 'WORST_FRACTION']))]
pool_rows.to_csv(RES + r'\v928_gate_proportion_ci_pool23.csv', index=False)
print('pool23 rows saved')
