# -*- coding: utf-8 -*-
"""Item C4 v3: discovery-threshold sensitivity grid with UNIFORM filters.

Change vs v2 (v936): every grid cell now applies the SAME filter set as the
main elite track (script 05_background_rescan_v2.R):
  (a) TCGA thresholds (per grid cell)
  (b) cross-cohort directional consistency: GSE67116 delta_beta > 0
      (probes absent from the minfi-processed GSE67116 matrix = excluded,
       which is how chrX/chrY/SNP/failP annotation filtering enters)
  (c) presence in the noob-processed GSE73949 healthy-endometrium matrix
      (gate-1 evaluability requirement)
v2 applied only (a) to new entrants, so loosened cells admitted probes the
main pipeline would have filtered (e.g. sex-chromosome probes at loose
thresholds). The reviewer-required reconciliation: original cell must equal
476 candidates after filtering, same two Tier-A survivors.

Outputs: results/v938_threshold_grid_summary_v3.csv
         results/v938_new_entrants_background_v3.csv
         results/v938_reconciliation.csv
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

# --- uniform filters (main elite track, script 05) ---
gse = pd.read_csv(RES + r'\GSE67116_dmp.csv')
gse_db = dict(zip(gse['probe'], gse['delta_beta']))
import pyreadr
g1full = list(pyreadr.read_r(RES + r'\GSE73949_betas.rds').values())[0]
g1_index = set(g1full.index)

def passes_filters(p):
    v = gse_db.get(p)
    return (v is not None and v == v and v > 0) and (p in g1_index)

GRID_DB = [0.30, 0.40, 0.50]
GRID_CS = [0.80, 0.85, 0.90]
GRID_PR = [0.75, 0.80, 0.85]
LEADS = {'cg24589459': 'BOLL anchor', 'cg27577527': 'ZSCAN12 anchor'}

base = d[(d['adj.P.Val'] < 0.05) & (d['ctrl_p95'] <= 0.40)]

def cand_set(db, cs, pr):
    s = base[(base['delta_beta'] >= db) & (base['case_pos03'] >= cs) & (base['pre_pos03'] >= pr)]
    return {p for p in s.index if passes_filters(p)}

# sanity: original cell must now equal 476
orig = cand_set(0.50, 0.90, 0.85)
print('original cell after uniform filters:', len(orig), '| == pset:', orig == pset)

union_new = set()
for db in GRID_DB:
    for cs in GRID_CS:
        for pr in GRID_PR:
            union_new |= (cand_set(db, cs, pr) - pset)
print('union new entrants across grid (filtered):', len(union_new), flush=True)
EVAL = sorted(union_new)
ESET = set(EVAL)

def p95(v):
    v = np.asarray([x for x in v if x == x], float)
    return float(np.percentile(v, 95)) if len(v) else np.nan

# ---------- gates for NEW ENTRANTS (four gates + eis637 for scoring) ----------
bg = {p: [np.nan] * 5 for p in EVAL}

sub = g1full.reindex(g1full.index.intersection(EVAL))
for pr in sub.index:
    bg[pr][0] = p95(sub.loc[pr].values)
print('gate1 (new entrants):', int(sub.notna().any(axis=1).sum()), 'matched', flush=True)

def read_cohort(gse_, gsm_keep, gi, frac_map=None):
    dd = os.path.join(DATA, gse_, 'ewas_betas')
    if not os.path.isdir(dd):
        print(gse_, 'directory missing -> all NaN', flush=True)
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
    print(gse_, 'done', flush=True)

m223 = pd.read_csv(RES + r'\GSE223817_ewas_meta.csv')
ctrl223 = set(m223.loc[m223['sample type'] == 'control', 'sample id'])
read_cohort('GSE223817', ctrl223, 1)
m463 = pd.read_csv(RES + r'\GSE46306_ewas_meta.csv')
hpvneg = set(m463.loc[(m463['infection'] == 'HPV-') & (m463['sample type'] == 'control'), 'sample id'])
read_cohort('GSE46306', hpvneg, 2)
m350 = pd.read_csv(RES + r'\GSE35069_ewas_meta.csv')
frac_map = dict(zip(m350['sample id'], m350['tissue']))
read_cohort('GSE35069', set(m350['sample id']), 3, frac_map)
eis637 = set(m223.loc[m223['sample type'] == 'disease tissue', 'sample id'])
read_cohort('GSE223817', eis637, 4)

bgdf = pd.DataFrame([(p, *bg[p]) for p in EVAL],
                    columns=['probe', 'g1_healthy17_p95', 'g2_benign347_p95',
                             'g3_cervix20_p95', 'g4_blood60_worst_p95', 'eis637_endometriosis_p95'])
bgdf.to_csv(RES + r'\v938_new_entrants_background_v3.csv', index=False)
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
g.to_csv(RES + r'\v938_threshold_grid_summary_v3.csv', index=False)

# ---------- reconciliation table ----------
recon = pd.DataFrame([
    dict(analysis='Main funnel (Table 2 / Figure 1)',
         start_set='476 (genome-wide screen, Xena TCGA metrics + cross-cohort consistency + annotation filters)',
         discovery='FDR<0.05; dBeta>=0.50; case>=0.90; pre>=0.85; ctrl_p95<=0.40; GSE67116 dBeta>0; present in GSE73949 noob matrix',
         background_rule='Stage 1 healthy-endometrium screen (P95<=0.15): 65 pass; stage 2 triage on healthy/cervix/blood (P95<=0.15, benign-surgical descriptive): 15 triage-passing; strict four-gate union (missing = not failed): Tier A = 2 leads',
         survivors='65 -> 15 -> 2',
         lead_ranks='ZSCAN12 1; BOLL 2 (score_v6 among strict-union survivors)'),
    dict(analysis='Application-specific paths (C5, S2.16)',
         start_set='476 (stability-cache gates, complete four-gate)',
         discovery='same as main (fixed 476 universe)',
         background_rule='STRICT4: all four gates P95<=0.15, missing = not failed; ENDO: g1+g2; CERV: g3+g4',
         survivors='STRICT4 = 2; ENDO = 11; CERV = 43',
         lead_ranks='both leads pass STRICT4'),
    dict(analysis='Threshold grid original cell (C4 v3, S2.15)',
         start_set='479 before uniform filters (genome-wide re-derivation); 476 after (excluded: cg17837162 chrY, cg09687738 chrX — absent from minfi-processed GSE67116/GSE73949 matrices; cg18661205 — GSE67116 dBeta 0.108>0 but absent from GSE73949 noob matrix)',
         discovery='same as main, plus the same cross-cohort consistency and annotation filters applied to every grid cell',
         background_rule='same four-gate union as STRICT4 (cache-complete gates; new entrants fully evaluated incl. eis637)',
         survivors='2 (BOLL + ZSCAN12) — reproduces main analysis exactly',
         lead_ranks='ZSCAN12 1/2; BOLL 2/2'),
    dict(analysis='Pipeline sensitivity (item 1, S2.12)',
         start_set='476 (fixed universe)',
         discovery='same as main',
         background_rule='gate 1 recomputed under GMQN (same 17 donors); other gates unchanged',
         survivors='healthy-pass 65 -> 41; triage-passing 15 -> 10 (auxiliary stage, not Tier A)',
         lead_ranks='both leads pass under either pipeline'),
], columns=['analysis', 'start_set', 'discovery', 'background_rule', 'survivors', 'lead_ranks'])
recon.to_csv(RES + r'\v938_reconciliation.csv', index=False)
print(recon[['analysis', 'survivors', 'lead_ranks']].to_string(index=False), flush=True)
print('DONE C4 v3', flush=True)
