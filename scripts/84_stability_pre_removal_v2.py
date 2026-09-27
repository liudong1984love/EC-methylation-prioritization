"""v2 companion: premenopausal-criterion-removal sensitivity under the same
nested-strata design as run_bootstrap_v2.py.

Discovery criteria WITHOUT the premenopausal-positivity requirement and WITHOUT
the premenopausal score term: disc = dBeta>=0.50 & overall pos>=0.90 & ctrl<=0.40;
score = 0.25*s_bg + 0.20*overall positivity (renormalisation unchanged, as in v1).
"""
import numpy as np
import pandas as pd
import time

OUT = r"C:/Users/ld/WorkBuddy/2026-08-19-15-50-19/results"  # results dir of the analysis base
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
frac = d["g4_frac"]
frac_groups = [np.where(frac == t)[0] for t in np.unique(frac)]
nP = len(probes)
LEADS = {"cg24589459": "BOLL anchor", "cg27577527": "ZSCAN12 anchor"}
lead_i = {p: int(np.where(probes == p)[0][0]) for p in LEADS}
COV = np.vstack([~np.isnan(G[k]).all(axis=1) for k in range(4)]).T

def gates_p95(gidx_list):
    out = []
    for k, sel in enumerate(gidx_list):
        M = G[k]
        if k == 3:
            pf = [np.nanquantile(M[:, g], 0.95, axis=1, method="linear") for g in sel if len(g) > 0]
            out.append(np.nanmax(np.vstack(pf), axis=0))
        else:
            out.append(np.nanquantile(M[:, sel], 0.95, axis=1, method="linear"))
    return np.vstack(out).T

tierA_cnt = np.zeros(nP, dtype=int)
top2_both = 0
t0 = time.time()
for b in range(B):
    p_i = rng.choice(pre_idx, size=len(pre_idx), replace=True)
    q_i = rng.choice(nonpre_idx, size=len(nonpre_idx), replace=True)
    t_i = np.concatenate([p_i, q_i])
    n_i = rng.choice(normal_idx, size=len(normal_idx), replace=True)
    db = np.nanmean(X[:, t_i], 1) - np.nanmean(X[:, n_i], 1)
    case = np.nanmean(X[:, t_i] > 0.3, 1)
    ctrl = np.nanquantile(np.where(np.isin(np.arange(X.shape[1]), n_i)[None, :], X, np.nan),
                          0.95, axis=1, method="linear")
    disc = (db >= 0.50) & (case >= 0.90) & (ctrl <= 0.40)
    gs = gates_p95([rng.integers(0, G[0].shape[1], G[0].shape[1]),
                    rng.integers(0, G[1].shape[1], G[1].shape[1]),
                    rng.integers(0, G[2].shape[1], G[2].shape[1]),
                    [g[rng.integers(0, len(g), len(g))] for g in frac_groups]])
    triage = np.where(COV, gs <= 0.15, True).all(1)
    tA = disc & triage
    tierA_cnt += tA
    if tA.sum() >= 1:
        worst_b = np.where(COV, gs, -np.inf).max(1)
        sbg = np.maximum(0, 1 - worst_b / 0.15)
        score = 0.25 * sbg + 0.20 * case
        r = pd.Series(score[tA]).rank(ascending=False, method="min").values
        rmap = dict(zip(np.where(tA)[0], r))
        li = [i for p, i in lead_i.items() if tA[i] and rmap[i] <= 2]
        if len(li) == 2:
            top2_both += 1
    if (b + 1) % 500 == 0:
        print(f"  {b+1}/{B} ({time.time()-t0:.0f}s)", flush=True)

for p, name in LEADS.items():
    i = lead_i[p]
    print(f"{name} ({p}): P(Tier A | premenopausal criterion removed) = {tierA_cnt[i]/B:.4f}")
print("P(both leads Tier A and top-2 | removed):", round(top2_both / B, 4))
print("DONE")
