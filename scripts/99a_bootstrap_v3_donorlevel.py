"""Conditional selection-stability analysis v3 (v9.32 review item 2).

Change vs v2: GSE35069 blood gate is resampled at DONOR level instead of
within-fraction. v2 drew 6 samples with replacement inside each of the 10
fractions independently, breaking the cross-fraction correlation structure of
the six donors. v3 draws 6 donors with replacement (verified donor identities
from ArrayExpress E-GEOD-35069 sample titles, Supplementary Data S2s), keeps
ALL fractions of each drawn donor, then computes per-fraction P95 and takes
the maximum across fractions - identical to the observed statistic.

Donor mapping: g4 columns follow the row order of GSE35069_ewas_meta.csv
(see extract_stability_data.py); donors assigned via
v937_gse35069_verified_donors.csv (gsm -> donor 1..6).

Everything else unchanged: nested tumour strata (28 pre + 403 other),
sample-level resampling for adjacent normals and gates 1-3, B=2000, same seed.

Outputs: results_stability_bootstrap_v3.csv, stability_summary_v3.txt
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
audit = pd.DataFrame(d["audit"])
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

# ---- verified donor mapping for the 60 g4 columns ----
meta = pd.read_csv(RES + r"/GSE35069_ewas_meta.csv")
don = pd.read_csv(RES + r"/v937_gse35069_verified_donors.csv")
gsm2donor = dict(zip(don["gsm"], don["donor"]))
assert G[3].shape[1] == len(meta) == 60
col_donor = np.array([gsm2donor[g] for g in meta["sample id"]])
col_frac = frac  # length-60 fraction labels from the cache build
fractions = sorted(set(col_frac))
# col index for each (donor, fraction) pair
pair_col = {}
for j in range(60):
    pair_col[(int(col_donor[j]), col_frac[j])] = j
assert len(pair_col) == 60, "expected a complete 6x10 donor x fraction grid"
donors = sorted(set(int(x) for x in col_donor))
assert donors == [1, 2, 3, 4, 5, 6]


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
    return np.vstack(out).T  # nP x 3


def blood_gate_observed():
    pf = [np.nanquantile(G[3][:, np.where(col_frac == f)[0]], 0.95, axis=1, method="linear")
          for f in fractions]
    return np.nanmax(np.vstack(pf), axis=0)


def blood_gate_donor_boot(drawn):
    """drawn: array of 6 donor ids (with replacement)."""
    pf = []
    for f in fractions:
        cols = [pair_col[(int(dn), f)] for dn in drawn]
        pf.append(np.nanquantile(G[3][:, cols], 0.95, axis=1, method="linear"))
    return np.nanmax(np.vstack(pf), axis=0)


COV = np.vstack([~np.isnan(G[k]).all(axis=1) for k in range(4)]).T

obs_db, obs_case, obs_pre, obs_ctrl = tcga_stats(tumor_idx, normal_idx, pre_idx)
obs_gates = np.column_stack([
    gates123_p95([np.arange(G[0].shape[1]), np.arange(G[1].shape[1]), np.arange(G[2].shape[1])]),
    blood_gate_observed()])
obs_disc = (obs_db >= 0.50) & (obs_case >= 0.90) & (obs_pre >= 0.85) & (obs_ctrl <= 0.40)
obs_triage = np.where(COV, obs_gates <= 0.15, True).all(1)
obs_tierA = obs_disc & obs_triage
print("observed: discovery-retained:", int(obs_disc.sum()), "/", nP,
      "| strict triage pass:", int(obs_triage.sum()), "| Tier A:", int(obs_tierA.sum()), flush=True)

# ---------- bootstrap (nested strata; donor-level blood) ----------
tierA_cnt = np.zeros(nP, dtype=int)
disc_cnt = np.zeros(nP, dtype=int)
gatepass_cnt = np.zeros((nP, 4), dtype=int)
gate_p95_sum = np.zeros((nP, 4))
rank_rec = {p: [] for p in LEADS}
tierA_n = []
pre_below = {p: 0 for p in LEADS}
pre_vals = {p: [] for p in LEADS}
top2_both = 0
cond_fail = {p: dict(db=0, case=0, pre=0, ctrl=0, g1=0, g2=0, g3=0, g4=0,
                     disc_ok=0, triage_ok=0, tierA=0,
                     pre_and_cervix_ok=0, pre_ok=0, cervix_ok=0) for p in LEADS}
t0 = time.time()

for b in range(B):
    p_i = rng.choice(pre_idx, size=len(pre_idx), replace=True)
    q_i = rng.choice(nonpre_idx, size=len(nonpre_idx), replace=True)
    t_i = np.concatenate([p_i, q_i])
    n_i = rng.choice(normal_idx, size=len(normal_idx), replace=True)
    db, case, pre, ctrl = tcga_stats(t_i, n_i, p_i)
    c_db = db >= 0.50; c_case = case >= 0.90; c_pre = pre >= 0.85; c_ctrl = ctrl <= 0.40
    disc = c_db & c_case & c_pre & c_ctrl
    g123 = gates123_p95([rng.integers(0, G[0].shape[1], G[0].shape[1]),
                         rng.integers(0, G[1].shape[1], G[1].shape[1]),
                         rng.integers(0, G[2].shape[1], G[2].shape[1])])
    g4b = blood_gate_donor_boot(rng.choice(donors, size=6, replace=True))
    gs = np.column_stack([g123, g4b])
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
res.to_csv(OUT + "/results_stability_bootstrap_v3.csv", index=False)

# ---- comparison vs v2 (within-fraction resampling) ----
try:
    v2 = pd.read_csv(OUT + "/results_stability_bootstrap_v2.csv").set_index("probe")
    cmp_rows = []
    for p, name in LEADS.items():
        cmp_rows.append(dict(lead=name, probe=p,
                             P_TierA_v2_within_fraction=v2.loc[p, "P_TierA"],
                             P_TierA_v3_donor_level=round(tierA_cnt[lead_i[p]] / B, 4),
                             P_gate4_pass_v2=v2.loc[p, "P_gate4_pass"],
                             P_gate4_pass_v3=round(gatepass_cnt[lead_i[p], 3] / B, 4)))
    cmp_df = pd.DataFrame(cmp_rows)
    cmp_df.to_csv(RES + "/v939_blood_resampling_comparison.csv", index=False)
    print(cmp_df.to_string(index=False), flush=True)
except FileNotFoundError:
    print("v2 results not found for comparison", flush=True)

with open(OUT + "/stability_summary_v3.txt", "w", encoding="utf-8") as fh:
    w = lambda *a: (print(*a, flush=True), fh.write(" ".join(str(x) for x in a) + "\n"))
    w("Conditional selection-stability analysis within the 476-probe candidate universe")
    w("v3: GSE35069 blood resampled at DONOR level (verified identities, S2s);")
    w("6 donors drawn with replacement, all 10 fractions per drawn donor retained,")
    w("per-fraction P95 then max - preserving cross-fraction donor correlation.")
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
