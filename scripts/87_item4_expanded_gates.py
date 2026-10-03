# -*- coding: utf-8 -*-
"""Item 4 follow-up: expanded background gates with newly censused cohorts.
New benign surgical controls (450K/GMQN): GSE130028(16) + GSE87621(5) + GSE108143(2) = +23
New healthy epithelial (EPIC): GSE122148(8, obese-annotated)
Disease-background expansion: GSE73948 (26 eutopic endometriosis, 450K)
"""
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
    return dict(n=n, k=k, prop=round(p, 4), ci_lo=round(lo, 4), ci_hi=round(hi, 4),
                p95=round(float(np.percentile(v, 95)), 4) if n else np.nan)

def load(gse, gsm_keep=None):
    d = os.path.join(DATA, gse, 'ewas_betas')
    out = {p: {} for p in probes}
    if not os.path.isdir(d):
        print('MISSING DIR', d)
        return out
    for f in os.listdir(d):
        if not f.endswith('.txt'):
            continue
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
    return out

def ctrl_set(gse):
    m = pd.read_csv(RES + f'\\{gse}_ewas_meta.csv')
    return set(m.loc[m['sample type'] == 'control', 'sample id'])

rows = []
SETS = []

# new benign surgical controls (450K)
b13 = load('GSE130028', ctrl_set('GSE130028'))
b87 = load('GSE87621', ctrl_set('GSE87621'))
b10 = load('GSE108143', ctrl_set('GSE108143'))
SETS.append(('benign_new16_GSE130028', b13))
SETS.append(('benign_new5_GSE87621', b87))
SETS.append(('benign_new2_GSE108143', b10))
# combined new 450K benign
bnew = {p: dict(list(b13[p].items()) + list(b87[p].items()) + list(b10[p].items())) for p in probes}
SETS.append(('benign_new23_pooled', bnew))

# GSE223817 existing controls for pooled benign
m22 = pd.read_csv(RES + r'\GSE223817_ewas_meta.csv')
c22 = set(m22.loc[m22['sample type'] == 'control', 'sample id'])
b22 = load('GSE223817', c22)
ball = {p: dict(list(b22[p].items()) + list(bnew[p].items())) for p in probes}
SETS.append(('benign_all370_pooled', ball))

# healthy epithelial EPIC (obese-annotated)
m12 = pd.read_csv(RES + r'\GSE122148_ewas_meta.csv')
print('GSE122148 meta:', m12.to_string()[:600])
b12 = load('GSE122148', None)
SETS.append(('healthy_epithelial8_GSE122148', b12))

# eutopic endometriosis disease background
b48 = load('GSE73948', None)
SETS.append(('eutopic_endo26_GSE73948', b48))

for name, mat in SETS:
    for pr in probes:
        s = stats_row(list(mat[pr].values()))
        rows.append(dict(probe=pr, cohort=name, **s))
    n0 = len(next(iter(mat.values()))) if mat else 0
    print(name, 'samples:', n0)

df = pd.DataFrame(rows)
df.to_csv(RES + r'\v930_expanded_background_476.csv', index=False)
print('saved v930_expanded_background_476.csv', df.shape)

# headline: anchors + partner across new cohorts
key_probes = ['cg27577527', 'cg24589459', 'cg07495363', 'cg23180938', 'cg23164203']
sub = df[df.probe.isin(key_probes)]
pd.set_option('display.width', 250)
print(sub.pivot_table(index='probe', columns='cohort', values='p95').to_string())
print()
print(sub.pivot_table(index='probe', columns='cohort', values='k').to_string())
print()
print(sub.pivot_table(index='probe', columns='cohort', values='ci_hi').to_string())

# CI-pass counts per cohort
summ = []
for name, _ in SETS:
    s = df[df.cohort == name]
    pt = s['p95'] <= THR
    ci = s['ci_hi'] <= 0.05
    summ.append(dict(cohort=name, n=int(s['n'].max()), point_pass=int(pt.sum()),
                     ci_pass=int(ci.sum()), uncertain=int((pt & ~ci).sum())))
print()
print(pd.DataFrame(summ).to_string(index=False))
