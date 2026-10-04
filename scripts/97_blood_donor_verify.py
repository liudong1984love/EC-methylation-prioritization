# -*- coding: utf-8 -*-
"""Verified-donor re-run of the GSE35069 blood bootstrap (S2c) + BOLL expression n.

1. Blood: reproduce Supplementary_Data_S2c with donor identities VERIFIED from
   ArrayExpress E-GEOD-35069 SDRF sample titles (v937_gse35069_verified_donors.csv).
   - pooled_p95 over all 60 samples; naive 10,000-draw bootstrap CI
   - donor-clustered bootstrap (resample 6 donors, keep all their fractions)
   - gate statistic = max per-fraction P95; donor-clustered CI
   Compare with existing S2c (order-inferred pairing).

2. Expression: from UCEC_HiSeqV2.gz count tumour/adjacent-normal samples with
   BOLL expression and medians (sample-type from TCGA barcode positions 14-15),
   restricted to barcodes present in the 450K methylation TSV header.
Outputs: results/v937_blood_donorcl_verified.csv, results/v937_boll_expr_n.csv
"""
import gzip
import numpy as np
import pandas as pd

BASE = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19'
RES = BASE + r'\results'
DATA = BASE + r'\data'

don = pd.read_csv(RES + r'\v937_gse35069_verified_donors.csv')
probes = pd.read_csv(BASE + r'\Supplementary_Data_S2c_blood_donor_clustered_bootstrap.csv')['probe'].tolist()
print('probes:', probes)

# read betas for the 9 probes across 60 GSM
import os
vals = {p: {} for p in probes}
dd = os.path.join(DATA, 'GSE35069', 'ewas_betas')
PSET = set(probes)
for f in os.listdir(dd):
    gsm = f.split('.')[0]
    with open(os.path.join(dd, f), encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            if line[:10] in PSET:
                parts = line.rstrip('\n').split('\t')
                if parts[0] in PSET:
                    vals[parts[0]][gsm] = float(parts[1])

def p95(v):
    return float(np.percentile(np.asarray(v, float), 95))

rng = np.random.default_rng(20261004)
B = 10000
rows = []
donors = sorted(don['donor'].unique())
fracs = sorted(don['fraction'].unique())
# fast lookup: gsm -> (donor index, fraction index)
didx = {k: i for i, k in enumerate(donors)}
fidx = {k: i for i, k in enumerate(fracs)}
gsm_di = dict(zip(don['gsm'], don['donor'].map(didx)))
gsm_fi = dict(zip(don['gsm'], don['fraction'].map(fidx)))
gsms = don['gsm'].tolist()
for p in probes:
    # grid[d, f] = beta
    grid = np.full((6, 10), np.nan)
    for g in gsms:
        if g in vals[p]:
            grid[gsm_di[g], gsm_fi[g]] = vals[p][g]
    v60 = grid.flatten()
    pooled = p95(v60)
    boot = np.percentile(rng.choice(v60, size=(B, 60), replace=True), 95, axis=1)
    nlo, nhi = np.percentile(boot, [2.5, 97.5])
    gate = float(np.nanpercentile(grid, 95, axis=0).max())
    # donor-clustered bootstrap: resample donor rows, keep all fractions
    cp = np.empty(B); cg = np.empty(B)
    for b in range(B):
        sd = rng.integers(0, 6, 6)
        sub = grid[sd, :]                      # 6 x 10
        cg[b] = np.percentile(sub, 95, axis=0).max()
        cp[b] = np.percentile(sub, 95)
    rows.append(dict(probe=p, pooled_p95=round(pooled, 4),
                     naive_lo=round(nlo, 4), naive_hi=round(nhi, 4),
                     donorcl_pooled_lo=round(np.percentile(cp, 2.5), 4),
                     donorcl_pooled_hi=round(np.percentile(cp, 97.5), 4),
                     gate_maxtype_p95=round(gate, 4),
                     gate_donorcl_lo=round(np.percentile(cg, 2.5), 4),
                     gate_donorcl_hi=round(np.percentile(cg, 97.5), 4)))
new = pd.DataFrame(rows)
new.to_csv(RES + r'\v937_blood_donorcl_verified.csv', index=False)
old = pd.read_csv(BASE + r'\Supplementary_Data_S2c_blood_donor_clustered_bootstrap.csv')
cmp = new.merge(old, on='probe', suffixes=('_verified', '_s2c'))
for c in ['pooled_p95', 'gate_maxtype_p95', 'gate_donorcl_hi']:
    cmp[c + '_diff'] = (cmp[c + '_verified'] - cmp[c + '_s2c']).round(4)
print(cmp[['probe', 'pooled_p95_diff', 'gate_maxtype_p95_diff', 'gate_donorcl_hi_verified', 'gate_donorcl_hi_s2c']].to_string(index=False))
print('max |gate P95 diff|:', cmp['gate_maxtype_p95_diff'].abs().max())
print('any gate CI crossing 0.15 (verified):', bool((new['gate_donorcl_hi'] > 0.15).any() & (new['gate_maxtype_p95'] <= 0.15).any()))

# ---------- BOLL expression n ----------
hdr = None
with gzip.open(DATA + r'\UCEC_HiSeqV2.gz', 'rt') as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
print('expr samples:', len(hdr) - 1)
codes = [h[13:15] for h in hdr[1:]]
from collections import Counter
print('sample-type codes:', Counter(codes))
# methylation TSV header
with open(DATA + r'\TCGA_UCEC_methylation450_beta.tsv', encoding='utf-8', errors='ignore') as fh:
    mhdr = fh.readline().rstrip('\n').split('\t')
mbc = set(h[:15] for h in mhdr[1:])
ebc = [h[:15] for h in hdr[1:]]
paired = [i for i, b in enumerate(ebc) if b in mbc]
pt = [i for i in paired if ebc[i][13:15] == '01']
pn = [i for i in paired if ebc[i][13:15] == '11']
print('paired tumour:', len(pt), '| paired adjacent normal:', len(pn))
# BOLL expression row
boll = None
with gzip.open(DATA + r'\UCEC_HiSeqV2.gz', 'rt') as fh:
    fh.readline()
    for line in fh:
        if line.startswith('BOLL'):
            parts = line.rstrip('\n').split('\t')
            if parts[0] == 'BOLL':
                boll = np.array([float(x) if x not in ('', 'NA') else np.nan for x in parts[1:]])
                break
if boll is not None:
    import numpy as np
    tv = boll[pt]; nv = boll[pn]
    print('BOLL expr tumour: median %.3f | >1 frac %.3f' % (np.nanmedian(tv), float((tv > 1).mean())))
    print('BOLL expr normal: median %.3f | >1 frac %.3f' % (np.nanmedian(nv), float((nv > 1).mean())))
    pd.DataFrame([
        dict(item='paired_tumours', value=len(pt)),
        dict(item='paired_adjacent_normals', value=len(pn)),
        dict(item='BOLL_expr_tumour_median_log2x1', value=round(float(np.nanmedian(tv)), 3)),
        dict(item='BOLL_expr_normal_median_log2x1', value=round(float(np.nanmedian(nv)), 3)),
        dict(item='BOLL_expr_gt1_tumour_frac', value=round(float((tv > 1).mean()), 3)),
        dict(item='BOLL_expr_gt1_normal_frac', value=round(float((nv > 1).mean()), 3)),
    ]).to_csv(RES + r'\v937_boll_expr_n.csv', index=False)
print('DONE')
