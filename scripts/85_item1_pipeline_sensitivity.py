# -*- coding: utf-8 -*-
"""Item 1: single-pipeline (GMQN) rerun of gate 1 + noob/GMQN sensitivity on GSE73949.

The manuscript gate 1 (healthy premenopausal endometrium, n=17) used minfi-noob
betas (RDS) while gates 2-4 used EWAS Data Hub GMQN betas. Here we:
  A) recompute gate 1 under GMQN (EWAS per-sample betas, same 17 donors);
  B) quantify pass/fail concordance (P95 <= 0.15 rule) across all 476 elite probes;
  C) assess impact on the 65-probe healthy-pass set and the 15 finalists;
  D) assemble a single-pipeline (GMQN everywhere) four-gate summary.
Outputs: results/v931_pipeline_sensitivity_gse73949.csv
         results/v931_single_pipeline_allfour.csv
"""
import os, math, re
import numpy as np
import pandas as pd
import pyreadr

BASE = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19'
RES = BASE + r'\results'
THR = 0.15

audit = pd.read_csv(BASE + r'\Supplementary_Data_S2a_single_probe_audit.csv')
elite = pd.read_csv(RES + r'\elite_ranked_v6.csv').set_index('probe')
probes = audit['probe'].tolist()
pset = set(probes)
print('probes:', len(probes))

# ---------- A) GMQN gate 1 from EWAS per-sample files ----------
d = os.path.join(BASE, 'data', 'GSE73949', 'ewas_betas')
gmqn = {}  # probe -> {gsm: beta}
files = [f for f in os.listdir(d) if f.endswith('.txt')]
print('EWAS GSE73949 files:', len(files))
for f in sorted(files):
    gsm = f.split('.')[0]
    with open(os.path.join(d, f), encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            p = line[:10]
            if p in pset:
                pr, val = line.rstrip('\n').split('\t')[:2]
                if pr in pset:
                    try:
                        gmqn.setdefault(pr, {})[gsm] = float(val)
                    except ValueError:
                        pass
gsms = sorted({g for v in gmqn.values() for g in v})
G = pd.DataFrame({g: pd.Series({p: gmqn.get(p, {}).get(g, np.nan) for p in probes}) for g in gsms})
G = G.reindex(probes)
print('GMQN matrix:', G.shape, '| mean probes/sample:', int(G.notna().sum(0).mean()))

# ---------- B) noob gate 1 from RDS (audit source) ----------
rds = list(pyreadr.read_r(RES + r'\GSE73949_betas.rds').values())[0]
rds.columns = [c.split('_')[0] for c in rds.columns]
N = rds.reindex(probes)
common = [g for g in G.columns if g in N.columns]
print('noob matrix:', N.shape, '| common GSMs:', len(common), '| rds-only:', sorted(set(N.columns)-set(G.columns)), '| ewas-only:', sorted(set(G.columns)-set(N.columns)))
assert len(common) == 17, 'sample mismatch'

G = G[common]; N = N[common]

def p95(df):
    return pd.Series(np.nanpercentile(df.values.astype(float), 95, axis=1), index=df.index)

p95_g = p95(G); p95_n = p95(N)
k_g = (G > THR).sum(axis=1); k_n = (N > THR).sum(axis=1)

df = pd.DataFrame({
    'probe': probes,
    'gene': audit['gene'],
    'p95_noob': p95_n.values, 'p95_gmqn': p95_g.values,
    'k_noob': k_n.values, 'k_gmqn': k_g.values, 'n': len(common),
})
df['pass_noob'] = df['p95_noob'] <= THR
df['pass_gmqn'] = df['p95_gmqn'] <= THR
df['flip'] = df['pass_noob'] != df['pass_gmqn']
df['dp95'] = (df['p95_gmqn'] - df['p95_noob']).round(4)

# validate noob recompute vs audit gate1
chk = df.merge(audit[['probe', 'gate1_healthy_endometrium_P95']], on='probe')
dv = (chk['p95_noob'] - chk['gate1_healthy_endometrium_P95']).abs()
print(f'VALIDATE noob vs audit gate1: max|diff|={dv.max():.4f} n>0.01={int((dv>0.01).sum())}')

# ---------- concordance ----------
both = df.dropna(subset=['p95_noob', 'p95_gmqn'])
a = int(((both.pass_noob) & (both.pass_gmqn)).sum())
b = int(((both.pass_noob) & (~both.pass_gmqn)).sum())
c = int(((~both.pass_noob) & (both.pass_gmqn)).sum())
dd = int(((~both.pass_noob) & (~both.pass_gmqn)).sum())
n_tot = a + b + c + dd
po = (a + dd) / n_tot
pe = ((a + b) * (a + c) + (c + dd) * (b + dd)) / n_tot ** 2
kappa = (po - pe) / (1 - pe) if pe < 1 else float('nan')
print(f'concordance: pass/pass={a} noob-only={b} gmqn-only={c} fail/fail={dd} | kappa={kappa:.3f}')
print(f'flips total: {b+c}/{n_tot} ({100*(b+c)/n_tot:.1f}%)')

# ---------- impact on funnel sets ----------
hp = set(audit.loc[audit['funnel_pass_healthy_background'] == True, 'probe'])
tri = set(audit.loc[audit['funnel_pass_triage'] == True, 'probe'])
print('audit healthy-pass set:', len(hp), '| finalists:', len(tri))
g_mq_pass = set(df.loc[df['pass_gmqn'], 'probe'])
n_noob_pass = set(df.loc[df['pass_noob'], 'probe'])
print('healthy-pass under noob:', len(n_noob_pass), '| under GMQN:', len(g_mq_pass))
print('65-set members lost under GMQN:', sorted(hp - g_mq_pass))
print('probes newly passing under GMQN (not in 65):', len(g_mq_pass - hp))
fin = df[df.probe.isin(tri)][['probe', 'gene', 'p95_noob', 'p95_gmqn', 'k_noob', 'k_gmqn', 'pass_noob', 'pass_gmqn']]
print('\n15 finalists under both pipelines:')
pd.set_option('display.width', 250)
print(fin.to_string(index=False))
print('finalists failing gate1 under GMQN:', sorted(tri - g_mq_pass))

# ---------- pool-23 subset ----------
pool = pd.read_csv(RES + r'\wetlab_panel_23_primer_annotation.csv')
pool23 = pool['probe'].tolist()
sub = df[df.probe.isin(pool23)]
print(f'\npool-23 subset: flips {int(sub.flip.sum())}/{len(sub)}')

df.to_csv(RES + r'\v931_pipeline_sensitivity_gse73949.csv', index=False)
print('saved v931_pipeline_sensitivity_gse73949.csv', df.shape)

# ---------- D) single-pipeline (GMQN) four-gate summary ----------
v928 = pd.read_csv(RES + r'\v928_gate_proportion_ci_476.csv')
def get(co, gr):
    s = v928[(v928.cohort == co) & (v928.group == gr)].set_index('probe')['p95']
    return s
gates = pd.DataFrame({
    'probe': probes,
    'g1_healthy17_gmqn': p95_g.values,
    'g2_benign347': get('benign347_GSE223817', 'all').reindex(probes).values,
    'g3_cervix20': get('cervix20_HPVneg_GSE46306', 'all').reindex(probes).values,
    'g4_blood60_worstfrac': get('blood60_GSE35069', 'WORST_FRACTION').reindex(probes).values,
})
for c in gates.columns[1:]:
    gates[c + '_pass'] = gates[c] <= THR
gates['all4_pass_gmqn'] = gates[[c + '_pass' for c in gates.columns[1:5]]].all(axis=1)
gates['is_65set'] = gates.probe.isin(hp)
gates['is_finalist'] = gates.probe.isin(tri)
gates.to_csv(RES + r'\v931_single_pipeline_allfour.csv', index=False)
print('\nsingle-pipeline GMQN all-four pass:', int(gates['all4_pass_gmqn'].sum()))
print('finalists passing all four under single-pipeline GMQN:',
      int(gates.loc[gates.is_finalist, 'all4_pass_gmqn'].sum()), '/', len(tri))
print('saved v931_single_pipeline_allfour.csv', gates.shape)
