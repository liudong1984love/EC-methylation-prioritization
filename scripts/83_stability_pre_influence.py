"""Stability analysis step 3: premenopausal-case influence.
(a) jackknife over the 28 premenopausal tumours (deterministic)
(b) bootstrap variant with the premenopausal criterion and score term removed
"""
import numpy as np
import pandas as pd

OUT = r"C:/Users/ld/WorkBuddy/2026-08-19-15-50-19/results"  # results dir of the analysis base
d = np.load(OUT + "/stability_cache.npz", allow_pickle=True)
probes = d["probes"]
X = d["tcga"].astype(float)
tumor_idx = np.where(d["tumor"])[0]
normal_idx = np.where(d["normal"])[0]
pre_idx = np.where(d["pre_tumor"])[0]
G = [d["g1"].astype(float), d["g2"].astype(float), d["g3"].astype(float), d["g4"].astype(float)]
frac = d["g4_frac"]
frac_groups = [np.where(frac == t)[0] for t in np.unique(frac)]
LEADS = {"cg24589459": "BOLL anchor", "cg27577527": "ZSCAN12 anchor"}
lead_i = {p: int(np.where(probes == p)[0][0]) for p in LEADS}
COV = np.vstack([~np.isnan(G[k]).all(axis=1) for k in range(4)]).T

def gates_p95_obs():
    out = [np.nanquantile(G[k], 0.95, axis=1, method="linear") for k in range(3)]
    pf = [np.nanquantile(G[3][:, g], 0.95, axis=1, method="linear") for g in frac_groups]
    out.append(np.nanmax(np.vstack(pf), axis=0))
    return np.vstack(out).T

obs_gates = gates_p95_obs()
worst = np.where(COV, obs_gates, -np.inf).max(1)
sbg = np.maximum(0, 1 - worst / 0.15)
triage = np.where(COV, obs_gates <= 0.15, True).all(1)
db = np.nanmean(X[:, tumor_idx], 1) - np.nanmean(X[:, normal_idx], 1)
case = np.nanmean(X[:, tumor_idx] > 0.3, 1)
ctrl = np.nanquantile(X[:, normal_idx], 0.95, axis=1, method="linear")

print("=== (a) jackknife over 28 premenopausal tumours ===")
for p, i in lead_i.items():
    pv = X[i, pre_idx] > 0.3
    npos = int(pv.sum())
    print(f"{LEADS[p]} ({p}): pre_pos03 = {npos}/28 = {npos/28:.4f}")
    lo = []
    for j in range(28):
        rest = np.delete(pre_idx, j)
        lo.append(np.nanmean(X[i, rest] > 0.3))
    lo = np.array(lo)
    print(f"  leave-one-out pre_pos03: min {lo.min():.4f} max {lo.max():.4f} | any < 0.85: {(lo < 0.85).any()}")
    # does any single omission change the score ranking of the two leads?
    sc = 0.25 * sbg[i] + 0.25 * lo + 0.20 * case[i]
    print(f"  leave-one-out score range: {sc.min():.4f}-{sc.max():.4f}")

print()
print("=== (b) bootstrap variant: no premenopausal criterion / score term ===")
B = 2000
rng = np.random.default_rng(20260926)
tierA_cnt = {p: 0 for p in LEADS}
disc_cnt = {p: 0 for p in LEADS}
rank_rec = {p: [] for p in LEADS}
both_top2 = 0
for b in range(B):
    t_i = rng.choice(tumor_idx, size=len(tumor_idx), replace=True)
    n_i = rng.choice(normal_idx, size=len(normal_idx), replace=True)
    db_b = np.nanmean(X[:, t_i], 1) - np.nanmean(X[:, n_i], 1)
    case_b = np.nanmean(X[:, t_i] > 0.3, 1)
    ctrl_b = np.nanquantile(X[:, n_i], 0.95, axis=1, method="linear")
    disc = (db_b >= 0.50) & (case_b >= 0.90) & (ctrl_b <= 0.40)  # no pre criterion
    gs = np.vstack([np.nanquantile(G[0][:, rng.integers(0, 17, 17)], 0.95, axis=1, method="linear"),
                    np.nanquantile(G[1][:, rng.integers(0, 347, 347)], 0.95, axis=1, method="linear"),
                    np.nanquantile(G[2][:, rng.integers(0, 20, 20)], 0.95, axis=1, method="linear"),
                    np.nanmax(np.vstack([np.nanquantile(G[3][:, g[rng.integers(0, len(g), len(g))]],
                                                        0.95, axis=1, method="linear")
                                         for g in frac_groups]), axis=0)]).T
    tri = np.where(COV, gs <= 0.15, True).all(1)
    tA = disc & tri
    if tA.sum():
        worst_b = np.where(COV, gs, -np.inf).max(1)
        score = 0.25 * np.maximum(0, 1 - worst_b / 0.15) + 0.20 * case_b  # no pre term
        r = pd.Series(score[tA]).rank(ascending=False, method="min").values
        rmap = dict(zip(np.where(tA)[0], r))
        for p, i in lead_i.items():
            if tA[i]:
                rank_rec[p].append(int(rmap[i]))
        if sum(1 for p, i in lead_i.items() if tA[i] and rmap[i] <= 2) == 2:
            both_top2 += 1
    for p, i in lead_i.items():
        disc_cnt[p] += disc[i]
        tierA_cnt[p] += tA[i]

for p, i in lead_i.items():
    rr = np.array(rank_rec[p])
    print(f"{LEADS[p]}: P(discovery, no-pre) = {disc_cnt[p]/B:.4f} | P(Tier A, no-pre) = {tierA_cnt[p]/B:.4f}",
          f"| rank median {np.median(rr):.0f} P(rank<=2) {(rr<=2).mean():.3f}" if len(rr) else "")
print(f"P(both leads Tier A & top-2, no-pre) = {both_top2/B:.4f}")
print("DONE")
