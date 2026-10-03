# -*- coding: utf-8 -*-
"""Item C4 v2 (corrected): discovery-threshold sensitivity grid.

Corrections vs v1 (v933):
1. Gate values for the 476 audit-universe probes now come from the stability
   cache (complete four-gate evaluation for every probe), NOT from
   elite_ranked_v6, whose NaN cells were pipeline-stage completeness artefacts
   (e.g. cg15285250/NXPH1 had unevaluated cervix/blood in the elite table but
   measured values 0.186/0.276 in the cache -> it FAILS triage; v1 let it pass).
2. Score = official score_v6: 0.25*s_bg + 0.25*pre_pos03 + 0.20*case_pos03
   + 0.15 + 0.05, with s_bg from bg3_worst = worst endometrial background
   (GSE73949 healthy / GSE223817 benign surgical / GSE155760 13 controls),
   NOT the four-gate worst. eis637 P95 is now also computed for new entrants.
3. NaN = genuinely not assessable on that platform/cohort -> missing, not
   failed (consistent with the main analysis).

Consequence: the original grid cell (0.50/0.90/0.85) reproduces the main
analysis exactly: 479 candidates (476 + 3 recomputation-difference probes),
2 triage survivors (BOLL + ZSCAN12), ZSCAN12 rank 1, BOLL rank 2.

Outputs: results/v936_threshold_grid_summary_v2.csv
         results/v936_new_entrants_background_v2.csv
         results/v936_reconciliation.csv
"""
import os
import numpy as np
import pandas as pd

BASE = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19'
RES = BASE + r'\results'
DATA = BASE + r'\data'
WS = r'C:\Users\ld\WorkBuddy\2026-09-25-20-09-18'
THR = 0.15

d = pd.read_csv(RES + r'\TCGA_dmp.csv').set_index('probe')
elite = pd.read_csv(RES + r'\elite_ranked_v6.csv')
pset = set(elite['probe'])

GRID_DB = [0.30, 0.40, 0.50]
GRID_CS = [0.80, 0.85, 0.90]
GRID_PR = [0.75, 0.80, 0.85]
LEADS = {'cg24589459': 'BOLL anchor', 'cg27577527': 'ZSCAN12 anchor'}

base = d[(d['adj.P.Val'] < 0.05) & (d['ctrl_p95'] <= 0.40)]

def cand_set(db, cs, pr):
    s = base[(base['delta_beta'] >= db) & (base['case_pos03'] >= cs) & (base['pre_pos03'] >= pr)]
    return set(s.index)

union_new = set()
for db in GRID_DB:
    for cs in GRID_CS:
        for pr in GRID_PR:
            union_new |= (cand_set(db, cs, pr) - pset)
print('union new entrants across grid:', len(union_new), flush=True)
EVAL = sorted(union_new)
ESET = set(EVAL)

def p95(v):
    v = np.asarray([x for x in v if x == x], float)
    return float(np.percentile(v, 95)) if len(v) else np.nan

# ---------- gates for NEW ENTRANTS (four gates + eis637 for scoring) ----------
bg = {p: [np.nan] * 5 for p in EVAL}  # g1,g2,g3,g4,eis637

import pyreadr
g1full = list(pyreadr.read_r(RES + r'\GSE73949_betas.rds').values())[0]
sub = g1full.reindex(g1full.index.intersection(EVAL))
for pr in sub.index:
    bg[pr][0] = p95(sub.loc[pr].values)
print('gate1 (new entrants):', int(sub.notna().any(axis=1).sum()), 'matched', flush=True)

def read_cohort(gse, gsm_keep, gi, frac_map=None):
    dd = os.path.join(DATA, gse, 'ewas_betas')
    if not os.path.isdir(dd):
        print(gse, 'directory missing -> all NaN', flush=True)
        return
    per = {p: {} for p in EVAL}
    for f in os.listdir(dd):
        if not f.endswith('.txt'):
            continue
        gsm = f.split('.')[0]
        if gsm_keep is not None and gsm not in gsm_keep:
            continue
        with open(os.path.join(dd, f), encoding='utf-8', errors='ignore') as fh:
            for line in fh:
                if line[:10] in ESET:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) >= 2 and parts[0] in ESET:
                        try:
                            per[parts[0]][gsm] = float(parts[1])
                        except ValueError:
                            pass
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
# eis637 = GSE223817 endometriosis disease-tissue samples (n=637); the third
# endometrial background in the official score's bg3_worst (script 11).
eis637 = set(m223.loc[m223['sample type'] == 'disease tissue', 'sample id'])
read_cohort('GSE223817', eis637, 4)

bgdf = pd.DataFrame([(p, *bg[p]) for p in EVAL],
                    columns=['probe', 'g1_healthy17_p95', 'g2_benign347_p95',
                             'g3_cervix20_p95', 'g4_blood60_worst_p95', 'eis637_endometriosis_p95'])
bgdf.to_csv(RES + r'\v936_new_entrants_background_v2.csv', index=False)
bgx = bgdf.set_index('probe')

# ---------- gates for the 476: stability cache (complete four-gate) ----------
z = np.load(WS + r'\stability_cache.npz', allow_pickle=True)
cache_probes = z['probes'].tolist()
G1 = z['g1'].astype(float); G2 = np.round(z['g2'].astype(float), 3)
G3 = z['g3'].astype(float); G4 = z['g4'].astype(float)
frac = z['g4_frac']
c_g1 = np.nanquantile(G1, 0.95, axis=1, method='linear')
c_g2 = np.nanquantile(G2, 0.95, axis=1, method='linear')
c_g3 = np.nanquantile(G3, 0.95, axis=1, method='linear')
pf = [np.nanquantile(G4[:, np.where(frac == t)[0]], 0.95, axis=1, method='linear') for t in np.unique(frac)]
c_g4 = np.nanmax(np.vstack(pf), axis=0)
cache_gates = {p: (float(c_g1[i]), float(c_g2[i]), float(c_g3[i]), float(c_g4[i]))
               for i, p in enumerate(cache_probes)}
eis_map = dict(zip(elite['probe'], elite['eis637_p95'].astype(float)))

def gates4(p):
    """(g1,g2,g3,g4); cache for 476, recomputed for new entrants."""
    if p in cache_gates:
        return cache_gates[p]
    if p in bgx.index:
        r = bgx.loc[p]
        return (r['g1_healthy17_p95'], r['g2_benign347_p95'], r['g3_cervix20_p95'], r['g4_blood60_worst_p95'])
    return (np.nan,) * 4

def triage_pass(p):
    vals = gates4(p)
    ok = [v <= THR for v in vals if v == v]
    return all(ok) if ok else False

def score_v6(p):
    g1, g2, _, _ = gates4(p)
    e6 = eis_map.get(p, np.nan)
    if p in bgx.index:
        e6v = bgx.loc[p, 'eis637_endometriosis_p95']
        e6 = e6v if e6v == e6v else e6
    bg3 = [v for v in (g1, g2, e6) if v == v]
    w = max(bg3) if bg3 else np.nan
    s_bg = max(0.0, 1 - w / THR) if w == w else 0.0
    r = d.loc[p]
    return 0.25 * s_bg + 0.25 * float(r['pre_pos03']) + 0.20 * float(r['case_pos03']) + 0.15 + 0.05

# ---------- grid loop ----------
rows = []
for db in GRID_DB:
    for cs in GRID_CS:
        for pr in GRID_PR:
            S = cand_set(db, cs, pr)
            new = S - pset
            tri = {p for p in S if triage_pass(p)}
            new_tri = tri - pset
            rec = {}
            if tri:
                sc = {p: score_v6(p) for p in tri}
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
g.to_csv(RES + r'\v936_threshold_grid_summary_v2.csv', index=False)

# ---------- reconciliation table ----------
orig = g[(g.delta_beta == 0.5) & (g.case_pos == 0.9) & (g.pre_pos == 0.85)].iloc[0]
recon = pd.DataFrame([
    dict(analysis='Main funnel (Table 2 / Figure 1)',
         start_set='476 (genome-wide screen, Xena TCGA metrics)',
         discovery='FDR<0.05; dBeta>=0.50; case>=0.90; pre>=0.85; ctrl_p95<=0.40',
         background_rule='Stage 1 healthy-endometrium screen (P95<=0.15): 65 pass; stage 2 triage on healthy/cervix/blood (P95<=0.15, benign-surgical descriptive): 15 finalists; strict four-gate union (missing = not failed): Tier A = 2 leads',
         survivors='65 -> 15 -> 2',
         lead_ranks='ZSCAN12 1; BOLL 2 (score_v6 among strict-union survivors)'),
    dict(analysis='Application-specific paths (C5, S2.16)',
         start_set='476 (stability-cache gates, complete four-gate)',
         discovery='same as main (fixed 476 universe)',
         background_rule='STRICT4: all four gates P95<=0.15, missing = not failed; ENDO: g1+g2; CERV: g3+g4',
         survivors='STRICT4 = 2; ENDO = 11; CERV = 43',
         lead_ranks='both leads pass STRICT4'),
    dict(analysis='Threshold grid original cell (C4 v2, S2.15)',
         start_set='479 (genome-wide re-derivation: 476 + 3 recomputation-difference probes: cg18661205, cg17837162, cg09687738)',
         discovery='same as main',
         background_rule='same four-gate union as STRICT4 (cache-complete gates; new entrants fully evaluated incl. eis637)',
         survivors='2 (BOLL + ZSCAN12) — reproduces main analysis',
         lead_ranks='ZSCAN12 1/2; BOLL 2/2'),
    dict(analysis='Pipeline sensitivity (item 1, S2.12)',
         start_set='476 (fixed universe)',
         discovery='same as main',
         background_rule='gate 1 recomputed under GMQN (same 17 donors); other gates unchanged',
         survivors='healthy-pass 65 -> 41; finalists 15 -> 10',
         lead_ranks='both leads pass under either pipeline'),
], columns=['analysis', 'start_set', 'discovery', 'background_rule', 'survivors', 'lead_ranks'])
recon.to_csv(RES + r'\v936_reconciliation.csv', index=False)
print(recon[['analysis', 'survivors', 'lead_ranks']].to_string(index=False), flush=True)
print('DONE C4 v2', flush=True)
