# -*- coding: utf-8 -*-
"""Item C5 (lightweight): sampling-site parallel-gating analysis (exploratory).

Main analysis stays the strict four-gate rule. In parallel, two
application-specific gating paths are evaluated WITHOUT re-running the funnel:
  - ENDO  (physician-collected endometrial brush): gates = healthy premenopausal
          endometrium (GSE73949) + benign surgical endometrium (GSE223817);
          blood reported separately as bleeding-contamination risk
          (variant ENDO+BLOOD also reported).
  - CERV  (cervical / self-sampling): gates = HPV-negative cervix (GSE46306)
          + blood fractions (GSE35069).
Applied to (i) the 476 elite probes and (ii) the 23-probe experimental pool.
Outputs: results/v934_sampling_paths_476.csv
         results/v934_sampling_paths_pool23.csv
"""
import numpy as np
import pandas as pd

BASE = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19'
RES = BASE + r'\results'
WS = r'C:\Users\ld\WorkBuddy\2026-09-25-20-09-18'
THR = 0.15

# ---------- 476 elite probes from stability cache ----------
z = np.load(WS + r'\stability_cache.npz', allow_pickle=True)
probes = z['probes'].tolist()
audit = pd.DataFrame(z['audit'])
G = [z['g1'].astype(float), z['g2'].astype(float), z['g3'].astype(float), z['g4'].astype(float)]
g2m = np.round(G[1], 3)  # float32 storage artefact fix (3-decimal source precision)
frac = z['g4_frac']

def p95_axis(M):
    return np.nanquantile(M, 0.95, axis=1, method='linear')

g1p = p95_axis(G[0]); g2p = p95_axis(g2m); g3p = p95_axis(G[2])
pf = []
for t in np.unique(frac):
    cols = np.where(frac == t)[0]
    pf.append(np.nanquantile(G[3][:, cols], 0.95, axis=1, method='linear'))
g4p = np.nanmax(np.vstack(pf), axis=0)

def passrule(v):
    return True if v != v else (v <= THR)

rows = []
for i, p in enumerate(probes):
    vals = [g1p[i], g2p[i], g3p[i], g4p[i]]
    strict = all(passrule(v) for v in vals)
    endo = passrule(vals[0]) and passrule(vals[1])
    endo_blood = endo and passrule(vals[3])
    cerv = passrule(vals[2]) and passrule(vals[3])
    n_assess = sum(v == v for v in vals)
    rows.append(dict(probe=p, gene=audit.loc[i, 'gene'] if 'gene' in audit else '',
                     g1_healthy17=round(vals[0], 3) if vals[0] == vals[0] else np.nan,
                     g2_benign347=round(vals[1], 3) if vals[1] == vals[1] else np.nan,
                     g3_cervix20=round(vals[2], 3) if vals[2] == vals[2] else np.nan,
                     g4_blood60_worst=round(vals[3], 3) if vals[3] == vals[3] else np.nan,
                     n_gates_assessable=n_assess,
                     STRICT4=strict, ENDO_brush=endo, ENDO_plus_blood=endo_blood, CERV_selfsample=cerv))
df = pd.DataFrame(rows)
df.to_csv(RES + r'\v934_sampling_paths_476.csv', index=False)

n = len(df)
print('=== 476 elite probes: path retention ===')
for c in ['STRICT4', 'ENDO_brush', 'ENDO_plus_blood', 'CERV_selfsample']:
    print(f'  {c:16s} {int(df[c].sum()):3d}/{n}')
print('ENDO-only gain vs STRICT4 (cervix/blood-limited):', int((df.ENDO_brush & ~df.STRICT4).sum()))
print('CERV-only gain vs STRICT4 (endometrium-limited):', int((df.CERV_selfsample & ~df.STRICT4).sum()))
print('both paths pass:', int((df.ENDO_brush & df.CERV_selfsample).sum()))
print('either path pass:', int((df.ENDO_brush | df.CERV_selfsample).sum()))

# ---------- 23-probe experimental pool ----------
pool = pd.read_csv(RES + r'\wetlab_panel_23_primer_annotation.csv')
ci = pd.read_csv(RES + r'\v928_gate_proportion_ci_pool23.csv')
piv = {}
for (pr, co), g in ci.groupby(['probe', 'cohort']):
    piv.setdefault(pr, {})[co] = float(g['p95'].iloc[0])
CO = {'g1': 'healthy17_GSE73949_noob', 'g2': 'benign347_GSE223817',
      'g3': 'cervix20_HPVneg_GSE46306', 'g4': 'blood60_GSE35069'}
elite_probes = set(probes)
prow = []
for _, r in pool.iterrows():
    pr = r['probe']
    v = {k: piv.get(pr, {}).get(c, np.nan) for k, c in CO.items()}
    n_assess = sum(x == x for x in v.values())
    strict = all(passrule(x) for x in v.values())
    endo = passrule(v['g1']) and passrule(v['g2'])
    cerv = passrule(v['g3']) and passrule(v['g4'])
    prow.append(dict(probe=pr, gene=str(r['gene']).split(';')[0], in_476=pr in elite_probes,
                     g1_healthy17=v['g1'], g2_benign347=v['g2'], g3_cervix20=v['g3'], g4_blood60=v['g4'],
                     n_gates_assessable=n_assess,
                     STRICT4=strict, ENDO_brush=endo, CERV_selfsample=cerv))
pd.set_option('display.width', 250)
pool_df = pd.DataFrame(prow)
pool_df.to_csv(RES + r'\v934_sampling_paths_pool23.csv', index=False)
print('\n=== 23-probe pool: path matrix ===')
print(pool_df.to_string(index=False))
print('\nDONE C5', flush=True)
