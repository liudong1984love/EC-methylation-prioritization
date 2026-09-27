"""Conditional selection-stability analysis within the 476-probe candidate
universe (v2: nested/stratified tumour resampling).

Change vs v1: premenopausal tumours (n=28) are a SUBSET of all tumours (n=431),
so the two positivity metrics share data. v1 resampled the 431 and the 28
independently; v2 resamples the two disjoint strata (28 premenopausal / 403
other tumours) separately and pools them, so each replicate's overall-positivity
and premenopausal-positivity are computed from the SAME nested draw.

Also new: per-condition failure decomposition for the two lead probes
(each discovery criterion and each background gate counted separately),
to explain final Tier-A retention without independence assumptions.

Outputs: results_stability_bootstrap_v2.csv (per-probe), stability_summary_v2.txt
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
nonpre_idx = np.setdiff1d(tumor_idx, pre_idx)
assert len(pre_idx) == 28 and len(nonpre_idx) == 403 and len(tumor_idx) == 431
G = [d["g1"].astype(float), d["g2"].astype(float), d["g3"].astype(float), d["g4"].astype(float)]
frac = d["g4_frac"]
frac_groups = [np.where(frac == t)[0] for t in np.unique(frac)]
nP = len(probes)
LEADS = {"cg24589459": "BOLL anchor", "cg27577527": "ZSCAN12 anchor"}
lead_i = {p: int(np.where(probes == p)[0][0]) for p in LEADS}

def tcga_stats(t_idx, n_idx, p_idx):
    db = np.nanmean(X[:, t_idx], 1) - np.nanmean(X[:, n_idx], 1)
    case = np.nanmean(X[:, t_idx] > 0.3, 1)
    pre = np.nanmean(X[:, p_idx] > 0.3, 1)
    ctrl = np.nanquantile(np.where(np.isin(np.arange(X.shape[1]), n_idx)[None, :], X, np.nan),
                          0.95, axis=1, method="linear")
    return db, case, pre, ctrl

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

COV = np.vstack([~np.isnan(G[k]).all(axis=1) for k in range(4)]).T

obs_db, obs_case, obs_pre, obs_ctrl = tcga_stats(tumor_idx, normal_idx, pre_idx)
obs_gates = gates_p95([np.arange(G[0].shape[1]), np.arange(G[1].shape[1]),
                       np.arange(G[2].shape[1]), frac_groups])
obs_disc = (obs_db >= 0.50) & (obs_case >= 0.90) & (obs_pre >= 0.85) & (obs_ctrl <= 0.40)
obs_triage = np.where(COV, obs_gates <= 0.15, True).all(1)
obs_tierA = obs_disc & obs_triage
print("observed: discovery-retained:", int(obs_disc.sum()), "/", nP,
      "| strict triage pass:", int(obs_triage.sum()), "| Tier A:", int(obs_tierA.sum()), flush=True)

# ---------- bootstrap (nested strata) ----------
tierA_cnt = np.zeros(nP, dtype=int)
disc_cnt = np.zeros(nP, dtype=int)
gatepass_cnt = np.zeros((nP, 4), dtype=int)
gate_p95_sum = np.zeros((nP, 4))
rank_rec = {p: [] for p in LEADS}
tierA_n = []
pre_below = {p: 0 for p in LEADS}
pre_vals = {p: [] for p in LEADS}
top2_both = 0
# per-condition failure decomposition for leads
cond_fail = {p: dict(db=0, case=0, pre=0, ctrl=0, g1=0, g2=0, g3=0, g4=0,
                     disc_ok=0, triage_ok=0, tierA=0,
                     pre_and_cervix_ok=0, pre_ok=0, cervix_ok=0) for p in LEADS}
t0 = time.time()

for b in range(B):
    p_i = rng.choice(pre_idx, size=len(pre_idx), replace=True)          # 28 from pre stratum
    q_i = rng.choice(nonpre_idx, size=len(nonpre_idx), replace=True)    # 403 from other tumours
    t_i = np.concatenate([p_i, q_i])                                    # nested pooled 431
    n_i = rng.choice(normal_idx, size=len(normal_idx), replace=True)
    db, case, pre, ctrl = tcga_stats(t_i, n_i, p_i)
    c_db = db >= 0.50; c_case = case >= 0.90; c_pre = pre >= 0.85; c_ctrl = ctrl <= 0.40
    disc = c_db & c_case & c_pre & c_ctrl
    gs = gates_p95([rng.integers(0, G[0].shape[1], G[0].shape[1]),
                    rng.integers(0, G[1].shape[1], G[1].shape[1]),
                    rng.integers(0, G[2].shape[1], G[2].shape[1]),
                    [g[rng.integers(0, len(g), len(g))] for g in frac_groups]])
    gp = np.where(COV, gs <= 0.15, True)
    triage = gp.all(1)
    tA = disc & triage
    disc_cnt += disc
    gatepass_cnt += np.where(COV, gs <= 0.15, -1)
    gate_p95_sum += np.where(COV, gs, 0)
    tierA_cnt += tA
    tierA_n.append(int(tA.sum()))
    for p, i in lead_i.items():
        cf = cond_fail[p]
        cf["db"] += int(not c_db[i]); cf["case"] += int(not c_case[i])
        cf["pre"] += int(not c_pre[i]); cf["ctrl"] += int(not c_ctrl[i])
        for k, name in enumerate(["g1", "g2", "g3", "g4"]):
            if COV[i, k]:
                cf[name] += int(gs[i, k] > 0.15)
        cf["disc_ok"] += int(disc[i]); cf["triage_ok"] += int(triage[i]); cf["tierA"] += int(tA[i])
        cervix_ok = (gs[i, 2] <= 0.15) if COV[i, 2] else True
        cf["pre_ok"] += int(c_pre[i]); cf["cervix_ok"] += int(cervix_ok)
        cf["pre_and_cervix_ok"] += int(c_pre[i] and cervix_ok)
        pre_vals[p].append(pre[i])
        if pre[i] < 0.85:
            pre_below[p] += 1
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
res.to_csv(OUT + "/results_stability_bootstrap_v2.csv", index=False)

with open(OUT + "/stability_summary_v2.txt", "w", encoding="utf-8") as fh:
    w = lambda *a: (print(*a, flush=True), fh.write(" ".join(str(x) for x in a) + "\n"))
    w("Conditional selection-stability analysis within the 476-probe candidate universe")
    w("B =", B, "| seed", SEED, "| nested strata: 28 pre + 403 other tumours pooled per replicate")
    w("=" * 74)
    w("Tier A count per replicate: median", int(np.median(tierA_n)),
      "IQR", int(np.percentile(tierA_n, 25)), "-", int(np.percentile(tierA_n, 75)),
      "range", min(tierA_n), "-", max(tierA_n))
    w("P(both leads Tier A and ranked top-2):", round(top2_both / B, 4))
    w("")
    for p, name in LEADS.items():
        i = lead_i[p]
        cf = cond_fail[p]
        rr = np.array(rank_rec[p])
        w(f"--- {name} ({p}) ---")
        w("  P(Tier A) =", round(tierA_cnt[i] / B, 4),
          "| P(discovery ok) =", round(cf["disc_ok"] / B, 4),
          "| P(triage ok) =", round(cf["triage_ok"] / B, 4))
        w("  discovery-condition failures:  P(dBeta<0.50) =", round(cf["db"] / B, 4),
          " P(overall pos<0.90) =", round(cf["case"] / B, 4),
          " P(pre pos<0.85) =", round(cf["pre"] / B, 4),
          " P(adjacent-normal P95>0.40) =", round(cf["ctrl"] / B, 4))
        for k, gname in enumerate(["healthy endometrium", "benign surgical", "cervix", "blood"]):
            if COV[i, k]:
                w(f"  gate {k+1} {gname}: P(fail) =", round(cf[f'g{k+1}'] / B, 4),
                  "| P(pass) =", round(gatepass_cnt[i, k] / B, 4),
                  "| mean bootstrap P95 =", round(gate_p95_sum[i, k] / B, 4))
            else:
                w(f"  gate {k+1} {gname}: n/a (platform)")
        w("  joint: P(pre>=0.85) =", round(cf["pre_ok"] / B, 4),
          "| P(cervix pass) =", round(cf["cervix_ok"] / B, 4),
          "| P(pre AND cervix ok) =", round(cf["pre_and_cervix_ok"] / B, 4))
        if len(rr):
            w("  rank among Tier A: median", int(np.median(rr)),
              "| P(rank=1)", round((rr == 1).mean(), 3),
              "| P(rank<=2)", round((rr <= 2).mean(), 3),
              "| P(rank<=5)", round((rr <= 5).mean(), 3), f"(n={len(rr)})")
        pv = np.array(pre_vals[p])
        w("  pre_pos03: mean", round(pv.mean(), 4), "P2.5", round(np.percentile(pv, 2.5), 4),
          "P97.5", round(np.percentile(pv, 97.5), 4))
    w("")
    w("Top-15 probes by P(Tier A):")
    w(res.head(15).to_string(index=False))
print("DONE", flush=True)
