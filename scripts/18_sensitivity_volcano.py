"""Todo B (threshold sensitivity) + Todo F (volcano figure)."""
import csv
import numpy as np
import pandas as pd
from pathlib import Path

RES = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\results")

# ---------- B1: P95 背景阈值敏感性（69 个排雷探针） ----------
tri = pd.read_csv(RES / "triage_closedloop_v4.csv")
bg = pd.read_csv(RES / "bg347_recheck_v5.csv")[["probe", "ne347_p95"]]
v6 = pd.read_csv(RES / "elite_ranked_v6.csv")[["probe", "ne_bg_p95", "eis637_p95"]]
d = tri.merge(bg, on="probe", how="left").merge(v6, on="probe", how="left")
for c in ["ne347_p95", "eis637_p95", "cervix_p95", "blood_max_p95"]:
    d[c] = pd.to_numeric(d[c], errors="coerce")
d["bg_worst"] = d[["ne_bg_p95", "ne347_p95", "eis637_p95"]].max(axis=1)
d["cycle_ok"] = (d["cycle_mean_abs_d"] <= 0.05) & (d["cycle_frac_gt01"] <= 0.25)

print("=== 背景阈值敏感性 ===")
print(f"{'P95 gate':>10} {'union_pass_n':>12} {'BOLL_pass':>10} {'ZSCAN12_pass':>13} {'CDO1_pass':>10}")
for gate in [0.10, 0.15, 0.20]:
    ok = d[(d["bg_worst"] <= gate) & (d["cervix_p95"] <= gate) & (d["blood_max_p95"] <= gate) & d["cycle_ok"]]
    b = "cg24589459" in set(ok["probe"])
    z = "cg27577527" in set(ok["probe"])
    cdo1 = "cg23180938" in set(ok["probe"])
    print(f"{gate:>10.2f} {len(ok):>12} {str(b):>10} {str(z):>13} {str(cdo1):>10}")

# ---------- B2: 入围门槛敏感性（TCGA 精英池规模） ----------
t = pd.read_csv(RES / "TCGA_dmp.csv", usecols=["probe", "adj.P.Val", "delta_beta", "case_pos03", "pre_pos03", "ctrl_p95"])
gse = pd.read_csv(RES / "GSE67116_dmp.csv", usecols=["probe", "delta_beta"]).rename(columns={"delta_beta": "gse_db"})
t = t.merge(gse, on="probe", how="left")
key = ["cg24589459", "cg27577527", "cg23180938", "cg15790037", "cg10109500"]
print("\n=== 入围门槛敏感性 ===")
print(f"{'entry dBeta>=':>13} {'positivity>=':>13} {'elite_n':>8}", " ".join(f"{p[-6:]:>9}" for p in key))
for db in [0.4, 0.5, 0.6]:
    for pos in [0.80, 0.90, 0.95]:
        e = t[(t["adj.P.Val"] < 0.05) & (t["delta_beta"] >= db) & (t["case_pos03"] >= pos)
              & (t["pre_pos03"] >= max(pos - 0.05, 0.80)) & (t["ctrl_p95"] <= 0.40) & (t["gse_db"] > 0)]
        flags = ["  in" if p in set(e["probe"]) else " out" for p in key]
        print(f"{db:>13.2f} {pos:>13.2f} {len(e):>8}", " ".join(f"{f:>9}" for f in flags))

# ---------- B3: 综合评分权重扰动 ----------
v6f = pd.read_csv(RES / "elite_ranked_v6.csv")
for c in ["bg3_worst", "pre_pos03", "case_pos03", "score_v6"]:
    v6f[c] = pd.to_numeric(v6f[c], errors="coerce")
base_top10 = set(v6f.nlargest(10, "score_v6")["probe"])
print("\n=== 权重扰动 ===")
print(f"{'weights(bg,pre,case)':>22} {'top10_overlap_with_base':>24}")
for wb, wp, wc in [(0.20, 0.30, 0.20), (0.30, 0.20, 0.20), (0.25, 0.20, 0.25)]:
    s = (wb * (1 - v6f["bg3_worst"] / 0.15).clip(lower=0) + wp * v6f["pre_pos03"] + wc * v6f["case_pos03"] + 0.20)
    top10 = set(v6f.loc[s.nlargest(10).index, "probe"])
    print(f"{str((wb, wp, wc)):>22} {len(base_top10 & top10):>24}")

# ---------- F: 火山图 ----------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

t2 = pd.read_csv(RES / "TCGA_dmp.csv", usecols=["probe", "adj.P.Val", "delta_beta"])
t2 = t2.dropna(subset=["adj.P.Val", "delta_beta"])
t2["nlp"] = -np.log10(t2["adj.P.Val"].clip(lower=1e-300))
panel = set(pd.read_csv(RES / "wetlab_panel_23_primer_annotation.csv")["probe"])
five = {"cg05016408", "cg01268824", "cg18675097", "cg23180938", "cg01580681"}
sig = (t2["adj.P.Val"] < 0.05) & (t2["delta_beta"].abs() >= 0.3)

plt.figure(figsize=(9, 6.5), dpi=200)
plt.scatter(t2.loc[~sig, "delta_beta"], t2.loc[~sig, "nlp"], s=1, c="#C9D3E0", alpha=0.4, rasterized=True)
plt.scatter(t2.loc[sig, "delta_beta"], t2.loc[sig, "nlp"], s=2, c="#3D7BD9", alpha=0.5, rasterized=True, label="FDR<0.05 & |Δβ|≥0.3")
sel = t2["probe"].isin(panel | five)
colors = ["#D9534F" if p in five else "#E8A33D" for p in t2.loc[sel, "probe"]]
plt.scatter(t2.loc[sel, "delta_beta"], t2.loc[sel, "nlp"], s=26, c=colors, edgecolors="k", linewidths=0.4, zorder=5)
lab = {"cg24589459": "BOLL", "cg27577527": "ZSCAN12", "cg23180938": "CDO1",
       "cg01268824": "ZNF154", "cg18675097": "NKAPL", "cg01580681": "HAND2", "cg05016408": "LOC134466"}
for p, g in lab.items():
    row = t2[t2["probe"] == p]
    if len(row):
        plt.annotate(g, (row["delta_beta"].iloc[0], row["nlp"].iloc[0]),
                     textcoords="offset points", xytext=(6, 6), fontsize=9, fontweight="bold")
from matplotlib.lines import Line2D
leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor="#3D7BD9", markersize=6, label="Significant probes"),
       Line2D([0], [0], marker="o", color="w", markerfacecolor="#E8A33D", markeredgecolor="k", markersize=7, label="23-probe panel"),
       Line2D([0], [0], marker="o", color="w", markerfacecolor="#D9534F", markeredgecolor="k", markersize=7, label="Previous 5 CpGs")]
plt.legend(handles=leg, loc="upper left", frameon=False, fontsize=9)
plt.xlabel("Δβ (tumour − adjacent normal)")
plt.ylabel("−log10(adj. P)")
plt.title("TCGA-UCEC differential methylation (431 tumours vs 46 adjacent normals)")
plt.axvline(0.3, ls="--", lw=0.8, c="grey"); plt.axvline(-0.3, ls="--", lw=0.8, c="grey")
plt.tight_layout()
out = RES / "fig_volcano_tcga.png"
plt.savefig(out)
print("\nvolcano saved:", out)
