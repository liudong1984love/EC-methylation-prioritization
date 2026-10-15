"""Q1 companion: premenopausal-criterion-removal sensitivity (99b design, seed 20260927)
with BOTH score background-term variants recorded per replicate:
  v4 = worst of four resampled gates (as run in 99b)
  v3 = main-analysis formula: worst of {g1_boot, g2_boot, eis637_observed-fixed}
Discovery here excludes the premenopausal criterion and the premenopausal score term.
Outputs: results/v935_bootstrap_nopre_score_formula_compare.csv + console summary.
"""
import time

import numpy as np
import pandas as pd

OUT = r"C:/Users/ld/WorkBuddy/2026-09-25-20-09-18"
RES = r"C:/Users/ld/WorkBuddy/2026-08-19-15-50-19/results"
B = 2000
SEED = 20260927
rng = np.random.default_rng(SEED)

d = np.load(OUT + "/stability_cache.npz", allow_pickle=True)
probes = d["probes"]
X = d["tcga"].astype(float)
tumor_idx = np.where(d["tumor"])[0]
normal_idx = np.where(d["normal"])[0]
pre_idx = np.where(d["pre_tumor"])[0]
nonpre_idx = np.setdiff1d(tumor_idx, pre_idx)
G = [d["g1"].astype(float), d["g2"].astype(float), d["g3"].astype(float), d["g4"].astype(float)]
frac = d["g4_frac"].astype(str)
nP = len(probes)
LEADS = {"cg24589459": "BOLL anchor", "cg27577527": "ZSCAN12 anchor"}
lead_i = {p: int(np.where(probes == p)[0][0]) for p in LEADS}

meta = pd.read_csv(RES + r"/GSE35069_ewas_meta.csv")
don = pd.read_csv(RES + r"/v937_gse35069_verified_donors.csv")
gsm2donor = dict(zip(don["gsm"], don["donor"]))
col_donor = np.array([gsm2donor[g] for g in meta["sample id"]])
col_frac = frac
fractions = sorted(set(col_frac))
pair_col = {(int(col_donor[j]), col_frac[j]): j for j in range(60)}
donors = sorted(set(int(x) for x in col_donor))

v6 = pd.read_csv(RES + r"/elite_ranked_v6.csv")
eis_map = dict(zip(v6["probe"], v6["eis637_p95"].astype(float)))
eis_obs = np.array([eis_map.get(p, np.nan) for p in probes], dtype=float)

def tcga_stats(t_idx, n_idx, p_idx):
    db = np.nanmean(X[:, t_idx], 1) - np.nanmean(X[:, n_idx], 1)
    case = np.nanmean(X[:, t_idx] > 0.3, 1)
    pre = np.nanmean(X[:, p_idx] > 0.3, 1)
    ctrl = np.nanquantile(np.where(np.isin(np.arange(X.shape[1]), n_idx)[None, :], X, np.nan),
                          0.95, axis=1, method="linear")
    return db, case, pre, ctrl

def gates123_p95(sel_list):
    return np.vstack([np.nanquantile(G[k][:, sel], 0.95, axis=1, method="linear")
                      for k, sel in enumerate(sel_list)]).T

def blood_gate_donor_boot(drawn):
    pf = []
    for f in fractions:
        cols = [pair_col[(int(dn), f)] for dn in drawn]
        pf.append(np.nanquantile(G[3][:, cols], 0.95, axis=1, method="linear"))
    return np.nanmax(np.vstack(pf), axis=0)

COV = np.vstack([~np.isnan(G[k]).all(axis=1) for k in range(4)]).T
COV3 = np.column_stack([COV[:, 0], COV[:, 1], np.isfinite(eis_obs)])

# ---------- observed ranks under BOTH formulas (no resampling; Q1 item 3) ----------
obs_db, obs_case, obs_pre, obs_ctrl = tcga_stats(tumor_idx, normal_idx, pre_idx)
obs_g = np.column_stack([
    gates123_p95([np.arange(G[0].shape[1]), np.arange(G[1].shape[1]), np.arange(G[2].shape[1])]),
    np.nanmax(np.vstack([np.nanquantile(G[3][:, np.where(col_frac == f)[0]], 0.95, axis=1, method="linear")
                         for f in fractions]), axis=0)])
obs_tA = (obs_db >= 0.50) & (obs_case >= 0.90) & (obs_pre >= 0.85) & (obs_ctrl <= 0.40) \
         & np.where(COV, obs_g <= 0.15, True).all(1)
print("observed Tier-A probes:", [str(probes[i]) for i in np.where(obs_tA)[0]])
w4 = np.where(COV, obs_g, -np.inf).max(1)
s4 = np.maximum(0, 1 - w4 / 0.15)
sc4 = 0.25 * s4 + 0.25 * obs_pre + 0.20 * obs_case
g3e = np.column_stack([obs_g[:, 0], obs_g[:, 1], eis_obs])
w3 = np.where(COV3, g3e, -np.inf).max(1)
s3 = np.maximum(0, 1 - w3 / 0.15)
sc3 = 0.25 * s3 + 0.25 * obs_pre + 0.20 * obs_case
for p, i in lead_i.items():
    print(f"observed {p}: score_4gate={sc4[i]:.4f} score_3bg={sc3[i]:.4f}")
r4_obs = pd.Series(sc4[obs_tA]).rank(ascending=False, method="min").values
r3_obs = pd.Series(sc3[obs_tA]).rank(ascending=False, method="min").values
print("observed ranks 4-gate:", r4_obs.tolist(), "| 3-bg:", r3_obs.tolist())

# ---------- no-premenopausal-criterion bootstrap, dual formula ----------
rank_rec4 = {p: [] for p in LEADS}
rank_rec3 = {p: [] for p in LEADS}
tierA_cnt = np.zeros(nP, dtype=int)
top2_4 = 0
top2_3 = 0
t0 = time.time()
for b in range(B):
    p_i = rng.choice(pre_idx, size=len(pre_idx), replace=True)
    q_i = rng.choice(nonpre_idx, size=len(nonpre_idx), replace=True)
    t_i = np.concatenate([p_i, q_i])
    n_i = rng.choice(normal_idx, size=len(normal_idx), replace=True)
    db, case, pre, ctrl = tcga_stats(t_i, n_i, p_i)
    disc = (db >= 0.50) & (case >= 0.90) & (ctrl <= 0.40)   # no premenopausal criterion
    g123 = gates123_p95([rng.integers(0, G[0].shape[1], G[0].shape[1]),
                         rng.integers(0, G[1].shape[1], G[1].shape[1]),
                         rng.integers(0, G[2].shape[1], G[2].shape[1])])
    g4b = blood_gate_donor_boot(rng.choice(donors, size=6, replace=True))
    gs = np.column_stack([g123, g4b])
    tA = disc & np.where(COV, gs <= 0.15, True).all(1)
    tierA_cnt += tA
    if tA.sum() >= 1:
        w4b = np.where(COV, gs, -np.inf).max(1)
        sc4b = 0.25 * np.maximum(0, 1 - w4b / 0.15) + 0.20 * case   # no pre term
        g3eb = np.column_stack([g123[:, 0], g123[:, 1], eis_obs])
        w3b = np.where(COV3, g3eb, -np.inf).max(1)
        sc3b = 0.25 * np.maximum(0, 1 - w3b / 0.15) + 0.20 * case
        idxs = np.where(tA)[0]
        r4 = dict(zip(idxs, pd.Series(sc4b[tA]).rank(ascending=False, method="min").values))
        r3 = dict(zip(idxs, pd.Series(sc3b[tA]).rank(ascending=False, method="min").values))
        for p, i in lead_i.items():
            if tA[i]:
                rank_rec4[p].append(int(r4[i]))
                rank_rec3[p].append(int(r3[i]))
        if len([i for p, i in lead_i.items() if tA[i] and r4[i] <= 2]) == 2:
            top2_4 += 1
        if len([i for p, i in lead_i.items() if tA[i] and r3[i] <= 2]) == 2:
            top2_3 += 1
    if (b + 1) % 500 == 0:
        print(f"  {b+1}/{B}  ({time.time()-t0:.0f}s)", flush=True)

print("\n=== no-pre variant: Tier-A frequencies (score-independent) ===")
rows = []
for p, name in LEADS.items():
    i = lead_i[p]
    a4 = np.array(rank_rec4[p]); a3 = np.array(rank_rec3[p])
    row = dict(probe=p, P_TierA=round(tierA_cnt[i] / B, 4), retained=len(a4),
               v4_top2=round(float((a4 <= 2).mean()), 4) if len(a4) else None,
               v3_top2=round(float((a3 <= 2).mean()), 4) if len(a3) else None,
               v3_rank1=round(float((a3 == 1).mean()), 4) if len(a3) else None,
               v4_rank1=round(float((a4 == 1).mean()), 4) if len(a4) else None)
    rows.append(row)
    print(row)
print(f"joint top-2 (no-pre): four-gate {top2_4/B:.4f} | three-bg {top2_3/B:.4f}")
pd.DataFrame(rows).to_csv(RES + r"/v935_bootstrap_nopre_score_formula_compare.csv", index=False)
print("DONE")
