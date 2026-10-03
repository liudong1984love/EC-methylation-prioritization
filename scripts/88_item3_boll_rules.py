# -*- coding: utf-8 -*-
"""Item 3: replace the post-hoc BOLL OR rule with resampling-based rule selection.

Candidate region-level rules for the BOLL locus (anchor cg24589459 + partner
cg07495363), positivity threshold beta > 0.3 (the manuscript's tumour-positivity
definition):
  R1 anchor-only | R2 partner-only | R3 OR | R4 AND | R5 region-mean
Evaluated on:
  - TCGA tumours: all (431) and premenopausal (28) positivity
  - benign backgrounds: GSE223817 controls (347, EPIC), GSE87621 (5, 450K),
    GSE73948 eutopic endometriosis (26, 450K), GSE73949 healthy (17, GMQN)
Selection stability:
  - B=2000 bootstrap: premenopausal tumours (28) and benign controls (347) are
    resampled independently of each other (disjoint cohorts). The overall-tumour
    draw (431) is kept in the RNG stream for seed compatibility but is NOT used:
    overall-tumour positivity is not an input to rule selection.
    two pre-specified criteria:
      (A) max premenopausal positivity s.t. benign positivity <= 0.05
          (fallback: min benign positivity)
      (B) max Youden J = pre_positivity - benign_positivity
  - LOOCV over the 347 benign samples under criterion (B)
Outputs: results/v932_boll_rule_metrics.csv
         results/v932_boll_rule_selection.csv
"""
import os
import numpy as np
import pandas as pd

BASE = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19'
RES = BASE + r'\results'
WS = r'C:\Users\ld\WorkBuddy\2026-09-25-20-09-18'
THR = 0.3
B = 2000
rng = np.random.default_rng(20260930)

ANCHOR, PARTNER = 'cg24589459', 'cg07495363'

# ---------- load tumour + benign-347 from stability cache ----------
z = np.load(WS + r'\stability_cache.npz', allow_pickle=True)
probes = z['probes'].tolist()
ia, ip = probes.index(ANCHOR), probes.index(PARTNER)
X = z['tcga'].astype(float)           # 476 x samples
tumor, pre = z['tumor'], z['pre_tumor']
g2 = z['g2'].astype(float)            # 476 x 347 GSE223817 controls
# float32 storage in the npz turns exact 3-decimal EWAS betas (e.g. 0.300) into
# 0.30000001, flipping strict >0.3 positivity counts at the boundary (+4 false
# positives on the partner probe). EWAS Data Hub betas are 3-decimal: round back
# to source precision (verified identical to per-GSM source files; 159/347 OR).
g2 = np.round(g2, 3)
print('tumours:', int(tumor.sum()), '| pre:', int(pre.sum()), '| benign347:', g2.shape[1])

tum_a, tum_p = X[ia, tumor], X[ip, tumor]
pre_a, pre_p = X[ia, pre], X[ip, pre]
ben_a, ben_p = g2[ia], g2[ip]

def rules(a, p):
    a = np.asarray(a, float); p = np.asarray(p, float)
    return {
        'R1_anchor_only': a > THR,
        'R2_partner_only': p > THR,
        'R3_OR': (a > THR) | (p > THR),
        'R4_AND': (a > THR) & (p > THR),
        'R5_region_mean': (a + p) / 2 > THR,
    }

def pos(rule_dict):
    return {k: float(np.nanmean(v)) for k, v in rule_dict.items()}

# ---------- point estimates ----------
rows = []
t_all = pos(rules(tum_a, tum_p)); t_pre = pos(rules(pre_a, pre_p)); b_347 = pos(rules(ben_a, ben_p))
print('\n=== point estimates (positivity) ===')
for r in t_all:
    print(f'{r:16s} tumour-all {t_all[r]:.3f} | tumour-pre {t_pre[r]:.3f} | benign347 {b_347[r]:.3f}')

# ---------- other background cohorts (EWAS GMQN per-sample) ----------
def read2(gse, keep=None):
    d = os.path.join(BASE, 'data', gse, 'ewas_betas')
    out = {}
    for f in sorted(os.listdir(d)):
        if not f.endswith('.txt'):
            continue
        gsm = f.split('.')[0]
        if keep is not None and gsm not in keep:
            continue
        vals = {}
        with open(os.path.join(d, f), encoding='utf-8', errors='ignore') as fh:
            for line in fh:
                p10 = line[:10]
                if p10 == ANCHOR or p10 == PARTNER:
                    pr, v = line.rstrip('\n').split('\t')[:2]
                    try:
                        vals[pr] = float(v)
                    except ValueError:
                        pass
        out[gsm] = vals
    return out

m87 = pd.read_csv(RES + r'\GSE87621_ewas_meta.csv')
keep87 = set(m87.loc[m87['sample type'] == 'control', 'sample id'])
cohorts = {
    'GSE223817_ctrl347_EPIC': (ben_a, ben_p),
    'GSE87621_ctrl5_450K': None,
    'GSE73948_endo26_450K': None,
    'GSE73949_healthy17_GMQN': None,
}
ext = {'GSE87621_ctrl5_450K': read2('GSE87621', keep87),
       'GSE73948_endo26_450K': read2('GSE73948'),
       'GSE73949_healthy17_GMQN': read2('GSE73949')}
bg_extra = {}
for name, dat in ext.items():
    a = np.array([v.get(ANCHOR, np.nan) for v in dat.values()])
    p = np.array([v.get(PARTNER, np.nan) for v in dat.values()])
    bg_extra[name] = (a, p)

rule_names = list(t_all.keys())
for r in rule_names:
    row = {'rule': r, 'tumour_all': round(t_all[r], 4), 'tumour_pre': round(t_pre[r], 4),
           'benign347_GSE223817': round(b_347[r], 4)}
    for name, (a, p) in bg_extra.items():
        row[name] = round(pos(rules(a, p))[r], 4)
    rows.append(row)
met = pd.DataFrame(rows)
pd.set_option('display.width', 250)
print('\n', met.to_string(index=False))
met.to_csv(RES + r'\v932_boll_rule_metrics.csv', index=False)

# ---------- bootstrap selection ----------
def boot_pos(a, p, idx):
    return rules(a[idx], p[idx])

nT, nP, nB = tumor.sum(), pre.sum(), g2.shape[1]
selA = {r: 0 for r in rule_names}; selB = {r: 0 for r in rule_names}
fallbackA = 0
for b in range(B):
    ti = rng.integers(0, nT, nT)
    pi = rng.integers(0, nP, nP)
    bi = rng.integers(0, nB, nB)
    pre_r = pos(rules(pre_a[pi], pre_p[pi]))
    ben_r = pos(rules(ben_a[bi], ben_p[bi]))
    ok = [r for r in rule_names if ben_r[r] <= 0.05]
    if ok:
        selA[max(ok, key=lambda r: pre_r[r])] += 1
    else:
        fallbackA += 1
        selA[min(rule_names, key=lambda r: ben_r[r])] += 1
    selB[max(rule_names, key=lambda r: pre_r[r] - ben_r[r])] += 1
print('\n=== bootstrap selection frequencies (B=%d) ===' % B)
print('criterion A (max pre pos | benign<=0.05; fallback min benign):')
for r in rule_names:
    print(f'  {r:16s} {selA[r]/B:.3f}')
print('  (fallback used in', fallbackA, 'replicates)')
print('criterion B (max Youden J = pre - benign):')
for r in rule_names:
    print(f'  {r:16s} {selB[r]/B:.3f}')

# ---------- LOOCV over benign under criterion B ----------
selL = {r: 0 for r in rule_names}
pre_full = pos(rules(pre_a, pre_p))
for i in range(nB):
    mask = np.ones(nB, bool); mask[i] = False
    ben_r = pos(rules(ben_a[mask], ben_p[mask]))
    selL[max(rule_names, key=lambda r: pre_full[r] - ben_r[r])] += 1
print('\n=== benign LOOCV selection (criterion B, %d folds) ===' % nB)
for r in rule_names:
    print(f'  {r:16s} {selL[r]/nB:.3f}')

sel = pd.DataFrame([
    {'rule': r, 'boot_critA_freq': round(selA[r]/B, 4),
     'boot_critB_freq': round(selB[r]/B, 4),
     'loocv_critB_freq': round(selL[r]/nB, 4)} for r in rule_names])
sel.to_csv(RES + r'\v932_boll_rule_selection.csv', index=False)
print('\nsaved v932_boll_rule_metrics.csv + v932_boll_rule_selection.csv')

# ---------- anchor vs partner discordance in benign (mechanism check) ----------
both_hi = ((ben_a > THR) & (ben_p > THR)).sum()
only_p = ((ben_a <= THR) & (ben_p > THR)).sum()
only_a = ((ben_a > THR) & (ben_p <= THR)).sum()
print(f'\nbenign347 discordance: both>0.3: {both_hi} | partner-only: {only_p} | anchor-only: {only_a} | neither: {nB-both_hi-only_p-only_a}')
print('benign partner median/p95:', round(float(np.nanmedian(ben_p)),3), round(float(np.nanpercentile(ben_p,95)),3),
      '| anchor median/p95:', round(float(np.nanmedian(ben_a)),3), round(float(np.nanpercentile(ben_a,95)),3))
