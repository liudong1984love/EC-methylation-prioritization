"""Stability analysis step 2: bootstrap the full selection funnel (B replicates).

Per replicate:
  - resample TCGA strata (tumors / adjacent normals / premenopausal tumors)
  - resample each background cohort (blood: within cell-type fractions)
  - recompute discovery criteria, four background gates (P95, type-7),
    strict triage (all assessable gates <= 0.15), effective 3-term score
  - re-rank probes passing strict triage

Outputs: results_stability_bootstrap.csv (per-probe), stability_summary.txt
"""
import numpy as np
import pandas as pd
import time

OUT = r"C:/Users/ld/WorkBuddy/2026-08-19-15-50-19/results"  # results dir of the analysis base
B = 2000
SEED = 20260926
rng = np.random.default_rng(SEED)

d = np.load(OUT + "/stability_cache.npz", allow_pickle=True)
probes = d["probes"]
audit = pd.DataFrame(d["audit"])
X = d["tcga"].astype(float)
tumor_idx = np.where(d["tumor"])[0]
normal_idx = np.where(d["normal"])[0]
pre_idx = np.where(d["pre_tumor"])[0]
G = [d["g1"].astype(float), d["g2"].astype(float), d["g3"].astype(float), d["g4"].astype(float)]
frac = d["g4_frac"]
frac_groups = [np.where(frac == t)[0] for t in np.unique(frac)]
nP = len(probes)
LEADS = {"cg24589459": "BOLL anchor", "cg27577527": "ZSCAN12 anchor"}
lead_i = {p: int(np.where(probes == p)[0][0]) for p in LEADS}

# ---------- observed reference ----------
def tcga_stats(t_idx, n_idx, p_idx):
    db = np.nanmean(X[:, t_idx], 1) - np.nanmean(X[:, n_idx], 1)
    case = np.nanmean(X[:, t_idx] > 0.3, 1)
    pre = np.nanmean(X[:, p_idx] > 0.3, 1)
    ctrl = np.nanquantile(np.where(np.isin(np.arange(X.shape[1]), n_idx)[None, :], X, np.nan),
                          0.95, axis=1, method="linear")
    return db, case, pre, ctrl

def gates_p95(gidx_list):
    """gidx_list: per cohort, either column indices (pooled) or list-of-groups (blood)."""
    out = []
    for k, sel in enumerate(gidx_list):
        M = G[k]
        if k == 3:  # blood: per-fraction P95 then max
            pf = [np.nanquantile(M[:, g], 0.95, axis=1, method="linear") for g in sel if len(g) > 0]
            out.append(np.nanmax(np.vstack(pf), axis=0))
        else:
            out.append(np.nanquantile(M[:, sel], 0.95, axis=1, method="linear"))
    return np.vstack(out).T  # probes x 4

def assessable(gate_mat):
    """probes x 4 bool: gate assessable = cohort covers probe (not all-NaN)."""
    cov = []
    for k in range(4):
        cov.append(~np.isnan(G[k]).all(axis=1))
    return np.vstack(cov).T

COV = assessable(None)

obs_db, obs_case, obs_pre, obs_ctrl = tcga_stats(tumor_idx, normal_idx, pre_idx)
obs_gates = gates_p95([np.arange(G[0].shape[1]), np.arange(G[1].shape[1]),
                       np.arange(G[2].shape[1]), frac_groups])
obs_disc = (obs_db >= 0.50) & (obs_case >= 0.90) & (obs_pre >= 0.85) & (obs_ctrl <= 0.40)
obs_gate_pass = np.where(COV, obs_gates <= 0.15, True)
obs_triage = obs_gate_pass.all(1)
obs_tierA = obs_disc & obs_triage
worst = np.where(COV, obs_gates, -np.inf).max(1)
obs_sbg = np.maximum(0, 1 - worst / 0.15)
obs_score = 0.25 * obs_sbg + 0.25 * obs_pre + 0.20 * obs_case
obs_rank = pd.Series(obs_score[obs_tierA]).rank(ascending=False, method="min")
obs_rank_map = {int(idx): int(r) for idx, r in zip(np.where(obs_tierA)[0], obs_rank)}
print("observed: discovery-retained:", int(obs_disc.sum()), "/", nP, flush=True)
print("observed strict triage pass:", int(obs_triage.sum()), "| Tier A:", int(obs_tierA.sum()), flush=True)
for p, i in lead_i.items():
    print(f"  {LEADS[p]}: disc={obs_disc[i]} triage={obs_triage[i]} "
          f"rank={obs_rank_map.get(i, 'n/a')} "
          f"gates={np.round(obs_gates[i], 3)} cov={COV[i]}", flush=True)

# ---------- bootstrap ----------
tierA_cnt = np.zeros(nP, dtype=int)
disc_cnt = np.zeros(nP, dtype=int)
gatepass_cnt = np.zeros((nP, 4), dtype=int)
gate_p95_sum = np.zeros((nP, 4))
rank_rec = {p: [] for p in LEADS}
tierA_n = []
pre_below = {p: 0 for p in LEADS}
pre_vals = {p: [] for p in LEADS}
top2_both = 0
t0 = time.time()

for b in range(B):
    t_i = rng.choice(tumor_idx, size=len(tumor_idx), replace=True)
    n_i = rng.choice(normal_idx, size=len(normal_idx), replace=True)
    p_i = rng.choice(pre_idx, size=len(pre_idx), replace=True)
    db, case, pre, ctrl = tcga_stats(t_i, n_i, p_i)
    disc = (db >= 0.50) & (case >= 0.90) & (pre >= 0.85) & (ctrl <= 0.40)
    gs = gates_p95([rng.integers(0, G[0].shape[1], G[0].shape[1]),
                    rng.integers(0, G[1].shape[1], G[1].shape[1]),
                    rng.integers(0, G[2].shape[1], G[2].shape[1]),
                    [g[rng.integers(0, len(g), len(g))] for g in frac_groups]])
    gp = np.where(COV, gs <= 0.15, True)
    triage = gp.all(1)
    tA = disc & triage
    disc_cnt += disc
    gatepass_cnt += np.where(COV, gs <= 0.15, -1)  # -1 marks n/a
    gate_p95_sum += np.where(COV, gs, 0)
    tierA_cnt += tA
    tierA_n.append(int(tA.sum()))
    if tA.sum() >= 1:
        worst_b = np.where(COV, gs, -np.inf).max(1)
        sbg = np.maximum(0, 1 - worst_b / 0.15)
        score = 0.25 * sbg + 0.25 * pre + 0.20 * case
        r = pd.Series(score[tA]).rank(ascending=False, method="min").values
        idxs = np.where(tA)[0]
        rmap = dict(zip(idxs, r))
        for p, i in lead_i.items():
            if tA[i]:
                rank_rec[p].append(int(rmap[i]))
        li = [i for p, i in lead_i.items() if tA[i] and rmap[i] <= 2]
        if len(li) == 2:
            top2_both += 1
    for p, i in lead_i.items():
        pre_vals[p].append(pre[i])
        if pre[i] < 0.85:
            pre_below[p] += 1
    if (b + 1) % 200 == 0:
        print(f"  {b+1}/{B}  ({time.time()-t0:.0f}s)", flush=True)

# ---------- per-probe summary ----------
rows = []
for i, p in enumerate(probes):
    rows.append(dict(
        probe=p, gene=audit.loc[i, "gene"], disposition=audit.loc[i, "final_disposition"],
        P_discovery_retained=round(disc_cnt[i] / B, 4),
        P_gate1_pass=(round(gatepass_cnt[i, 0] / B, 4) if COV[i, 0] else None),
        P_gate2_pass=(round(gatepass_cnt[i, 1] / B, 4) if COV[i, 1] else None),
        P_gate3_pass=(round(gatepass_cnt[i, 2] / B, 4) if COV[i, 2] else None),
        P_gate4_pass=(round(gatepass_cnt[i, 3] / B, 4) if COV[i, 3] else None),
        P_TierA=round(tierA_cnt[i] / B, 4)))
res = pd.DataFrame(rows).sort_values("P_TierA", ascending=False)
res.to_csv(OUT + "/results_stability_bootstrap.csv", index=False)

with open(OUT + "/stability_summary.txt", "w", encoding="utf-8") as fh:
    w = lambda *a: (print(*a, flush=True), fh.write(" ".join(str(x) for x in a) + "\n"))
    w("Bootstrap stability of the 476-probe selection funnel | B =", B, "| seed", SEED)
    w("=" * 70)
    w("Tier A count per replicate: median", int(np.median(tierA_n)),
      "IQR", int(np.percentile(tierA_n, 25)), "-", int(np.percentile(tierA_n, 75)),
      "range", min(tierA_n), "-", max(tierA_n))
    w("P(both leads Tier A and ranked top-2):", round(top2_both / B, 4))
    w("")
    for p, name in LEADS.items():
        i = lead_i[p]
        rr = np.array(rank_rec[p])
        w(f"--- {name} ({p}) ---")
        w("  P(discovery retained) =", round(disc_cnt[i] / B, 4))
        for k, gname in enumerate(["healthy endometrium", "benign surgical", "cervix", "blood"]):
            if COV[i, k]:
                w(f"  P(gate {k+1} {gname} <=0.15) =", round(gatepass_cnt[i, k] / B, 4),
                  "| mean bootstrap P95 =", round(gate_p95_sum[i, k] / B, 4))
            else:
                w(f"  gate {k+1} {gname}: n/a (platform)")
        w("  P(Tier A) =", round(tierA_cnt[i] / B, 4))
        if len(rr):
            w("  rank among Tier A: median", int(np.median(rr)), "IQR",
              int(np.percentile(rr, 25)), "-", int(np.percentile(rr, 75)),
              "| P(rank=1)", round((rr == 1).mean(), 3),
              "| P(rank<=2)", round((rr <= 2).mean(), 3),
              "| P(rank<=5)", round((rr <= 5).mean(), 3), f"(n={len(rr)} replicates)")
        pv = np.array(pre_vals[p])
        w("  pre_pos03: mean", round(pv.mean(), 4), "P2.5", round(np.percentile(pv, 2.5), 4),
          "P97.5", round(np.percentile(pv, 97.5), 4), "| P(pre_pos03<0.85) =", round(pre_below[p] / B, 4))
    w("")
    w("Top-15 probes by P(Tier A):")
    w(res.head(15).to_string(index=False))
print("DONE", flush=True)
