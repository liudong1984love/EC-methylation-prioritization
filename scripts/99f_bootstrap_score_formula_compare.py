"""Q1 (NoteGPT pack verification): does the bootstrap rank summary change if the
score uses the MAIN-ANALYSIS background term instead of the four-gate worst?

Main analysis (elite_ranked_v6 / S2.3): s_bg = max(0, 1 - worst3/0.15),
worst3 = max over ASSESSABLE {g1 healthy17, g2 benign347, eis637 endometriosis}.
Bootstrap v3 (99a): sbg = worst over the four RESAMPLED gates (g1,g2,g3,g4);
eis637 never resampled (conditional stability design).

This script reruns the identical 2000-replicate donor-level bootstrap (same seed,
same draws) and, within each replicate, computes BOTH scores and BOTH rankings of
Tier-A probes. eis637 enters the main-formula score as its OBSERVED (fixed) P95,
consistent with the not-resampled design.

Outputs: results/v935_bootstrap_score_formula_compare.csv + console summary.
"""
import time

import numpy as np
import pandas as pd

OUT = r"C:/Users/ld/WorkBuddy/2026-09-25-20-09-18"
RES = r"C:/Users/ld/WorkBuddy/2026-08-19-15-50-19/results"
B = 2000
SEED = 20260926
rng = np.random.default_rng(SEED)

d = np.load(OUT + "/stability_cache.npz", allow_pickle=True)
probes = d["probes"]
X = d["tcga"].astype(float)
tumor_idx = np.where(d["tumor"])[0]
normal_idx = np.where(d["normal"])[0]
pre_idx = np.where(d["pre_tumor"])[0]
nonpre_idx = np.setdiff1d(tumor_idx, pre_idx)
assert len(pre_idx) == 28 and len(nonpre_idx) == 403 and len(tumor_idx) == 431
G = [d["g1"].astype(float), d["g2"].astype(float), d["g3"].astype(float), d["g4"].astype(float)]
frac = d["g4_frac"].astype(str)
nP = len(probes)
LEADS = {"cg24589459": "BOLL anchor", "cg27577527": "ZSCAN12 anchor"}
lead_i = {p: int(np.where(probes == p)[0][0]) for p in LEADS}

meta = pd.read_csv(RES + r"/GSE35069_ewas_meta.csv")
don = pd.read_csv(RES + r"/v937_gse35069_verified_donors.csv")
gsm2donor = dict(zip(don["gsm"], don["donor"]))
assert G[3].shape[1] == len(meta) == 60
col_donor = np.array([gsm2donor[g] for g in meta["sample id"]])
col_frac = frac
fractions = sorted(set(col_frac))
pair_col = {}
for j in range(60):
    pair_col[(int(col_donor[j]), col_frac[j])] = j
assert len(pair_col) == 60
donors = sorted(set(int(x) for x in col_donor))
assert donors == [1, 2, 3, 4, 5, 6]

# ---- observed eis637 P95 per probe (fixed, not resampled) ----
v6 = pd.read_csv(RES + r"/elite_ranked_v6.csv")
eis_map = dict(zip(v6["probe"], v6["eis637_p95"].astype(float)))
eis_obs = np.array([eis_map.get(p, np.nan) for p in probes], dtype=float)
print("eis637 available for", int(np.isfinite(eis_obs).sum()), "of", nP, "probes")
for p, i in lead_i.items():
    print(f"  {p}: eis637 P95 = {eis_obs[i]}")

def tcga_stats(t_idx, n_idx, p_idx):
    db = np.nanmean(X[:, t_idx], 1) - np.nanmean(X[:, n_idx], 1)
    case = np.nanmean(X[:, t_idx] > 0.3, 1)
    pre = np.nanmean(X[:, p_idx] > 0.3, 1)
    ctrl = np.nanquantile(np.where(np.isin(np.arange(X.shape[1]), n_idx)[None, :], X, np.nan),
                          0.95, axis=1, method="linear")
    return db, case, pre, ctrl

def gates123_p95(sel_list):
    out = []
    for k, sel in enumerate(sel_list):
        out.append(np.nanquantile(G[k][:, sel], 0.95, axis=1, method="linear"))
    return np.vstack(out).T

def blood_gate_donor_boot(drawn):
    pf = []
    for f in fractions:
        cols = [pair_col[(int(dn), f)] for dn in drawn]
        pf.append(np.nanquantile(G[3][:, cols], 0.95, axis=1, method="linear"))
    return np.nanmax(np.vstack(pf), axis=0)

COV = np.vstack([~np.isnan(G[k]).all(axis=1) for k in range(4)]).T
COV3 = np.column_stack([COV[:, 0], COV[:, 1], np.isfinite(eis_obs)])

rank_rec4 = {p: [] for p in LEADS}   # current bootstrap formula (4-gate worst)
rank_rec3 = {p: [] for p in LEADS}   # main-analysis formula (3 endometrial)
top2_both_4 = 0
top2_both_3 = 0
t0 = time.time()

for b in range(B):
    p_i = rng.choice(pre_idx, size=len(pre_idx), replace=True)
    q_i = rng.choice(nonpre_idx, size=len(nonpre_idx), replace=True)
    t_i = np.concatenate([p_i, q_i])
    n_i = rng.choice(normal_idx, size=len(normal_idx), replace=True)
    db, case, pre, ctrl = tcga_stats(t_i, n_i, p_i)
    disc = (db >= 0.50) & (case >= 0.90) & (pre >= 0.85) & (ctrl <= 0.40)
    g123 = gates123_p95([rng.integers(0, G[0].shape[1], G[0].shape[1]),
                         rng.integers(0, G[1].shape[1], G[1].shape[1]),
                         rng.integers(0, G[2].shape[1], G[2].shape[1])])
    g4b = blood_gate_donor_boot(rng.choice(donors, size=6, replace=True))
    gs = np.column_stack([g123, g4b])
    gp = np.where(COV, gs <= 0.15, True)
    triage = gp.all(1)
    tA = disc & triage
    if tA.sum() >= 1:
        # variant 4: worst of four resampled gates (current 99a)
        worst4 = np.where(COV, gs, -np.inf).max(1)
        sbg4 = np.maximum(0, 1 - worst4 / 0.15)
        score4 = 0.25 * sbg4 + 0.25 * pre + 0.20 * case
        # variant 3: main-analysis formula: worst of {g1_boot, g2_boot, eis637_obs}
        g3endo = np.column_stack([g123[:, 0], g123[:, 1], eis_obs])
        worst3 = np.where(COV3, g3endo, -np.inf).max(1)
        sbg3 = np.maximum(0, 1 - worst3 / 0.15)
        score3 = 0.25 * sbg3 + 0.25 * pre + 0.20 * case
        idxs = np.where(tA)[0]
        r4 = dict(zip(idxs, pd.Series(score4[tA]).rank(ascending=False, method="min").values))
        r3 = dict(zip(idxs, pd.Series(score3[tA]).rank(ascending=False, method="min").values))
        for p, i in lead_i.items():
            if tA[i]:
                rank_rec4[p].append(int(r4[i]))
                rank_rec3[p].append(int(r3[i]))
        li4 = [i for p, i in lead_i.items() if tA[i] and r4[i] <= 2]
        li3 = [i for p, i in lead_i.items() if tA[i] and r3[i] <= 2]
        if len(li4) == 2:
            top2_both_4 += 1
        if len(li3) == 2:
            top2_both_3 += 1
    if (b + 1) % 500 == 0:
        print(f"  {b+1}/{B}  ({time.time()-t0:.0f}s)", flush=True)

print("\n=== Q1 RESULT: rank summaries under two score formulas ===")
rows = []
for p, name in LEADS.items():
    a4 = np.array(rank_rec4[p]); a3 = np.array(rank_rec3[p])
    n4, n3 = len(a4), len(a3)
    assert n4 == n3
    row = dict(probe=p, retained_reps=n4,
               v4_rank1=round(float((a4 == 1).mean()), 4),
               v4_top2=round(float((a4 <= 2).mean()), 4),
               v4_top5=round(float((a4 <= 5).mean()), 4),
               v3_rank1=round(float((a3 == 1).mean()), 4),
               v3_top2=round(float((a3 <= 2).mean()), 4),
               v3_top5=round(float((a3 <= 5).mean()), 4),
               rank_changed=round(float((a4 != a3).mean()), 4))
    rows.append(row)
    print(row)
print(f"joint top-2 frequency: four-gate score {top2_both_4/B:.4f} | three-bg score {top2_both_3/B:.4f}")
pd.DataFrame(rows).to_csv(RES + r"/v935_bootstrap_score_formula_compare.csv", index=False)
print("saved results/v935_bootstrap_score_formula_compare.csv")
