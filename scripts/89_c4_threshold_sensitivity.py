# -*- coding: utf-8 -*-
"""Item C4: discovery-threshold sensitivity analysis (exploratory).

Pre-specified grid around the original strict discovery criteria
(FDR<0.05 and adjacent-normal P95<=0.40 held fixed):
  delta_beta  in {0.30, 0.40, 0.50}   (original 0.50)
  case_pos03  in {0.80, 0.85, 0.90}   (original 0.90)
  pre_pos03   in {0.75, 0.80, 0.85}   (original 0.85)
For each of 27 combos: candidate count, leads retained + score rank among
background-triage survivors, new entrants vs the strict 476, and how many new
entrants pass the four-gate background triage (union rule, n/a excluded).
Also decomposes which discovery criterion excludes known markers
(CDO1/AJAP1/GALR1/GYPC).

NOTE: thresholds were NOT chosen to protect the current leads; the grid is
symmetric looser-than-original (stricter direction is uninformative: the leads
sit at the boundary, so any stricter step removes them by construction).

Outputs: results/v933_threshold_grid_summary.csv
         results/v933_new_entrants_background.csv
         results/v933_known_markers_discovery.csv
"""
import os, math
import numpy as np
import pandas as pd

BASE = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19'
RES = BASE + r'\results'
DATA = BASE + r'\data'
THR = 0.15

d = pd.read_csv(RES + r'\TCGA_dmp.csv')
elite = pd.read_csv(RES + r'\elite_ranked_v6.csv')
pset = set(elite['probe'])
d = d.set_index('probe')

GRID_DB = [0.30, 0.40, 0.50]
GRID_CS = [0.80, 0.85, 0.90]
GRID_PR = [0.75, 0.80, 0.85]
LEADS = {'cg24589459': 'BOLL anchor', 'cg27577527': 'ZSCAN12 anchor'}

base = d[(d['adj.P.Val'] < 0.05) & (d['ctrl_p95'] <= 0.40)]

def cand_set(db, cs, pr):
    s = base[(base['delta_beta'] >= db) & (base['case_pos03'] >= cs) & (base['pre_pos03'] >= pr)]
    return set(s.index)

# ---- union of all grid candidates needing background evaluation ----
union_new = set()
for db in GRID_DB:
    for cs in GRID_CS:
        for pr in GRID_PR:
            union_new |= (cand_set(db, cs, pr) - pset)
print('union new entrants across grid:', len(union_new), flush=True)
EVAL = sorted(union_new)
ESET = set(EVAL)

# ---- background gates for the union set ----
def p95(v):
    v = np.asarray([x for x in v if x == x], float)
    return float(np.percentile(v, 95)) if len(v) else np.nan

bg = {p: [np.nan] * 4 for p in EVAL}

# gate 1: GSE73949 healthy premenopausal endometrium (minfi-noob RDS, full genome)
import pyreadr
g1full = list(pyreadr.read_r(RES + r'\GSE73949_betas.rds').values())[0]
sub = g1full.reindex(g1full.index.intersection(EVAL))
for pr in sub.index:
    bg[pr][0] = p95(sub.loc[pr].values)
print('gate1 done:', int(sub.notna().any(axis=1).sum()), 'probes matched', flush=True)

def read_cohort(gse, gsm_keep, gi, frac_map=None):
    dd = os.path.join(DATA, gse, 'ewas_betas')
    per = {p: {} for p in EVAL}
    files = [f for f in os.listdir(dd) if f.endswith('.txt')]
    n = 0
    for f in files:
        gsm = f.split('.')[0]
        if gsm_keep is not None and gsm not in gsm_keep:
            continue
        with open(os.path.join(dd, f), encoding='utf-8', errors='ignore') as fh:
            for line in fh:
                p10 = line[:10]
                if p10 in ESET:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) >= 2 and parts[0] in ESET:
                        try:
                            per[parts[0]][gsm] = float(parts[1])
                        except ValueError:
                            pass
        n += 1
        if n % 100 == 0:
            print(f'  {gse}: {n} files...', flush=True)
    if frac_map is None:
        for pr in EVAL:
            bg[pr][gi] = p95(list(per[pr].values()))
    else:
        for pr in EVAL:
            byf = {}
            for gsm, v in per[pr].items():
                byf.setdefault(frac_map.get(gsm, '?'), []).append(v)
            bg[pr][gi] = max((p95(v) for v in byf.values()), default=np.nan)
    print(gse, 'done', flush=True)

m223 = pd.read_csv(RES + r'\GSE223817_ewas_meta.csv')
ctrl223 = set(m223.loc[m223['sample type'] == 'control', 'sample id'])
read_cohort('GSE223817', ctrl223, 1)
m463 = pd.read_csv(RES + r'\GSE46306_ewas_meta.csv')
hpvneg = set(m463.loc[(m463['infection'] == 'HPV-') & (m463['sample type'] == 'control'), 'sample id'])
read_cohort('GSE46306', hpvneg, 2)
m350 = pd.read_csv(RES + r'\GSE35069_ewas_meta.csv')
frac_map = dict(zip(m350['sample id'], m350['tissue']))
read_cohort('GSE35069', set(m350['sample id']), 3, frac_map)

bgdf = pd.DataFrame([(p, *bg[p]) for p in EVAL],
                    columns=['probe', 'g1_healthy17_p95', 'g2_benign347_p95',
                             'g3_cervix20_p95', 'g4_blood60_worst_p95'])
bgdf.to_csv(RES + r'\v933_new_entrants_background.csv', index=False)
print('background table saved', flush=True)

bgx = bgdf.set_index('probe')
def triage_pass(probe):
    if probe in pset:
        row = elite.set_index('probe').loc[probe]
        vals = [row['ne_bg_p95'], row['ne347_p95'], row['cervix_p95'], row['blood_max_p95']]
    elif probe in bgx.index:
        r = bgx.loc[probe]
        vals = [r['g1_healthy17_p95'], r['g2_benign347_p95'], r['g3_cervix20_p95'], r['g4_blood60_worst_p95']]
    else:
        return False
    ok = [(v <= THR) for v in vals if v == v]
    return all(ok) if ok else False

def worst_bg(probe):
    if probe in pset:
        row = elite.set_index('probe').loc[probe]
        vals = [row['ne_bg_p95'], row['ne347_p95'], row['cervix_p95'], row['blood_max_p95']]
    else:
        r = bgx.loc[probe]
        vals = [r['g1_healthy17_p95'], r['g2_benign347_p95'], r['g3_cervix20_p95'], r['g4_blood60_worst_p95']]
    vv = [v for v in vals if v == v]
    return max(vv) if vv else np.nan

# ---- grid loop ----
rows = []
for db in GRID_DB:
    for cs in GRID_CS:
        for pr in GRID_PR:
            S = cand_set(db, cs, pr)
            new = S - pset
            tri = {p for p in S if triage_pass(p)}
            new_tri = tri - pset
            # score rank of leads among triage survivors
            rec = {}
            if tri:
                sc = {}
                for p in tri:
                    wb = worst_bg(p)
                    sbg = max(0.0, 1 - wb / THR) if wb == wb else 0.0
                    sc[p] = 0.25 * sbg + 0.25 * d.loc[p, 'pre_pos03'] + 0.20 * d.loc[p, 'case_pos03']
                rk = pd.Series(sc).rank(ascending=False, method='min')
                for lp, nm in LEADS.items():
                    if lp in S:
                        rec[nm] = ('in-set, triage-%s, rank %d/%d'
                                   % ('PASS' if lp in tri else 'fail',
                                      int(rk[lp]) if lp in tri else -1, len(tri)))
                    else:
                        rec[nm] = 'NOT in discovery set'
            rows.append(dict(delta_beta=db, case_pos=cs, pre_pos=pr,
                             n_candidates=len(S), n_new_vs_476=len(new),
                             n_triage_pass=len(tri), n_new_lowbg=len(new_tri),
                             **{nm.replace(' ', '_'): rec.get(nm, 'n/a') for nm in LEADS.values()}))
            print(rows[-1], flush=True)

g = pd.DataFrame(rows)
g.to_csv(RES + r'\v933_threshold_grid_summary.csv', index=False)

# ---- known-marker discovery decomposition ----
KNOWN = {'cg23180938': 'CDO1 (EcPredict kit)', 'cg13495205': 'AJAP1 (kit, top-dBeta)',
         'cg11835068': 'AJAP1 (kit, low-bg)', 'cg15484532': 'AJAP1 (kit)',
         'cg03502002': 'GALR1 (kit, top-dBeta)', 'cg04534765': 'GALR1 (kit)',
         'cg10390058': 'GALR1 (kit)', 'cg25390440': 'GYPC (WID-qEC)',
         'cg27577527': 'ZSCAN12 (lead)', 'cg24589459': 'BOLL anchor (lead)',
         'cg07495363': 'BOLL partner', 'cg10109500': 'GHSR'}
kr = []
for p, nm in KNOWN.items():
    if p not in d.index:
        continue
    r = d.loc[p]
    fails = []
    if not r['adj.P.Val'] < 0.05: fails.append('FDR')
    if not r['delta_beta'] >= 0.50: fails.append('dBeta<0.50')
    if not r['case_pos03'] >= 0.90: fails.append('case_pos<0.90')
    if not r['pre_pos03'] >= 0.85: fails.append('pre_pos<0.85')
    if not r['ctrl_p95'] <= 0.40: fails.append('ctrl_p95>0.40')
    # loosest grid point membership
    loose = (r['adj.P.Val'] < 0.05 and r['ctrl_p95'] <= 0.40 and r['delta_beta'] >= 0.30
             and r['case_pos03'] >= 0.80 and r['pre_pos03'] >= 0.75)
    kr.append(dict(probe=p, marker=nm, adj_P=float(r['adj.P.Val']),
                   delta_beta=round(float(r['delta_beta']), 3),
                   case_pos03=round(float(r['case_pos03']), 3),
                   pre_pos03=round(float(r['pre_pos03']), 3),
                   ctrl_p95=round(float(r['ctrl_p95']), 3),
                   in_strict_476=p in pset,
                   strict_fail=';'.join(fails) if fails else '(passes discovery)',
                   enters_at_loosest_grid=bool(loose)))
kdf = pd.DataFrame(kr)
kdf.to_csv(RES + r'\v933_known_markers_discovery.csv', index=False)
print(kdf.to_string(index=False), flush=True)
print('DONE C4', flush=True)
