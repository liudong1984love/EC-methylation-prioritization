"""Figure 4 (v9.17): cross-cohort performance of BOLL and ZSCAN12 regions.
Panel A: background P95 beta across background cohorts (gate 0.15 dashed).
Panel B: tumour positivity across cohorts.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

probes = ["cg24589459\n(BOLL anchor)", "cg07495363\n(BOLL partner)", "cg27577527\n(ZSCAN12 anchor)"]

# Panel A: background P95 (primary pipeline per cohort; ctrl13 = GSE155760 endometrial mucosa)
bg_cohorts = ["Healthy\nendometrium\n(n=17)", "Benign surgical\nendometrium\n(n=347)", "Cervix HPV-neg\n(n=20)", "Blood purified\n(n=42)", "Endometrial controls\nGSE155760 (n=13)"]
bg = {
    "cg24589459\n(BOLL anchor)": [0.110, 0.100, 0.141, 0.112, 0.059],
    "cg07495363\n(BOLL partner)": [0.133, 0.427, 0.080, 0.124, 0.266],
    "cg27577527\n(ZSCAN12 anchor)": [0.060, np.nan, 0.044, 0.058, np.nan],
}
# Panel B: tumour positivity (%) per cohort
tu_cohorts = ["TCGA\noverall\n(n=431)", "TCGA\npremenop.\n(n=28)", "GSE67116\ncarcinoma\n(n=33)", "GSE136791\ncarcinoma\n(n=69)", "GSE155760\nEEC\n(n=23)", "GSE178610\nFF\n(n=49)", "GSE93589\n(n=9)"]
tu = {
    "cg24589459\n(BOLL anchor)": [91.0, 89.3, 75.4, 87.0, 65.2, 55.1, np.nan],
    "cg07495363\n(BOLL partner)": [94.0, 100.0, 73.8, 94.2, 95.7, 73.5, np.nan],
    "cg27577527\n(ZSCAN12 anchor)": [90.0, 92.9, 87.7, np.nan, np.nan, np.nan, 88.9],
}
colors = ["#1f77b4", "#2ca02c", "#d62728"]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), dpi=300)

ax = axes[0]
x = np.arange(len(bg_cohorts)); w = 0.26
for k, p in enumerate(probes):
    vals = bg[p]
    ax.bar(x + (k - 1) * w, [np.nan_to_num(v) for v in vals], w, label=p.split("\n")[0], color=colors[k])
    for xi, v in zip(x + (k - 1) * w, vals):
        if np.isnan(v):
            ax.text(xi, 0.02, "n/a", ha="center", fontsize=7, color=colors[k], rotation=90)
ax.axhline(0.15, ls="--", lw=1, color="grey")
ax.text(len(bg_cohorts) - 0.5, 0.16, "P95 gate 0.15", fontsize=8, color="grey", ha="right")
ax.set_xticks(x); ax.set_xticklabels(bg_cohorts, fontsize=7.5)
ax.set_ylabel("Background P95 β"); ax.set_ylim(0, 0.5)
ax.set_title("A  Background by cohort (primary pipeline)", fontsize=10)
ax.legend(fontsize=7, loc="upper left")

ax = axes[1]
x = np.arange(len(tu_cohorts))
for k, p in enumerate(probes):
    vals = tu[p]
    ax.bar(x + (k - 1) * w, [np.nan_to_num(v) for v in vals], w, color=colors[k])
    for xi, v in zip(x + (k - 1) * w, vals):
        if np.isnan(v):
            ax.text(xi, 3, "n/a", ha="center", fontsize=7, color=colors[k], rotation=90)
ax.set_xticks(x); ax.set_xticklabels(tu_cohorts, fontsize=7.5)
ax.set_ylabel("Tumour positivity (β>0.3, %)"); ax.set_ylim(0, 112)
ax.set_title("B  Tumour-tissue positivity across cohorts", fontsize=10)

fig.tight_layout()
fig.savefig(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\results\v917_figure4_crosscohort.png", bbox_inches="tight")
print("saved results/v917_figure4_crosscohort.png")
