# -*- coding: utf-8 -*-
"""Fix cervix gate: use 20 HPV- controls only (audit-consistent, script 07 logic)."""
import os, math
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

def wilson(k, n):
    if n == 0:
        return (float('nan'),) * 3
    p = k / n
    d = 1 + Z * Z / n
    c = p + Z * Z / (2 * n)
    h = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n))
    return p, max(0.0, (c - h) / d), min(1.0, (c + h) / d)

def stats_row(vals):
    v = np.asarray([x for x in vals if x == x], dtype=float)
    n = len(v)
    k = int((v > THR).sum())
    p, lo, hi = wilson(k, n)
    return dict(n=n, k=k, prop=p, ci_lo=lo, ci_hi=hi,
                median=float(np.median(v)) if n else np.nan,
                p90=float(np.percentile(v, 90)) if n else np.nan,
                p95=float(np.percentile(v, 95)) if n else np.nan,
                max=float(v.max()) if n else np.nan)

meta = pd.read_csv(RES + r'\GSE46306_ewas_meta.csv')
keep = set(meta.loc[(meta['infection'] == 'HPV-') & (meta['sample type'] == 'control'), 'sample id'])
print('HPV- controls:', len(keep))

d = os.path.join(DATA, 'GSE46306', 'ewas_betas')
out = {p: {} for p in probes}
for f in os.listdir(d):
    if not f.endswith('.txt'):
        continue
    gsm = f.split('.')[0]
    if gsm not in keep:
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

rows = []
for pr in probes:
    s = stats_row(list(out[pr].values()))
    rows.append(dict(probe=pr, cohort='cervix20_HPVneg_GSE46306', group='all', **s))
new = pd.DataFrame(rows)

# cross-check vs elite_v6
chk = elite.set_index('probe')['cervix_p95']
both = pd.concat([new.set_index('probe')['p95'], chk], axis=1, keys=['mine', 'ref']).dropna()
diff = (both['mine'] - both['ref']).abs()
print(f'cervix20 vs cervix_p95: n={len(both)} maxdiff={diff.max():.4f} meandiff={diff.mean():.4f}')

# patch the main CSV
df = pd.read_csv(RES + r'\v928_gate_proportion_ci_476.csv')
df = df[df.cohort != 'cervix44_GSE46306']
df = pd.concat([df, new], ignore_index=True)
df.to_csv(RES + r'\v928_gate_proportion_ci_476.csv', index=False)
print('patched main CSV:', df.shape)

# redo summary
CUT = 0.05
summ = []
for cohort, group in [('healthy17_GSE73949_noob', 'all'), ('benign347_GSE223817', 'all'),
                      ('cervix20_HPVneg_GSE46306', 'all'), ('blood60_GSE35069', 'WORST_FRACTION')]:
    sub = df[(df.cohort == cohort) & (df.group == group)].set_index('probe')
    pt = sub['p95'] <= THR
    ci = sub['ci_hi'] <= CUT
    summ.append(dict(cohort=cohort, n_samples=int(sub['n'].max()),
                     point_pass=int(pt.sum()), ci_pass=int(ci.sum()),
                     point_pass_ci_fail=int((pt & ~ci).sum())))
sdf = pd.DataFrame(summ)
print(sdf.to_string(index=False))
sdf.to_csv(RES + r'\v928_gate_proportion_ci_summary.csv', index=False)

pool = pd.read_csv(RES + r'\wetlab_panel_23_primer_annotation.csv')['probe'].tolist()
pool_rows = df[df.probe.isin(pool) & (df.group.isin(['all', 'WORST_FRACTION']))]
pool_rows.to_csv(RES + r'\v928_gate_proportion_ci_pool23.csv', index=False)
print('pool23 updated')
