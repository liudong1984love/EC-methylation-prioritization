# 48_fig3_regen.py — regenerate Figure 3 with three-tier colouring
# green = clear pass (P95<=0.15 both, CI not crossing), yellow = borderline
# (point-estimate pass, bootstrap CI crosses gate; BOLL), red = compartment-restricted
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "C:/Users/ld/WorkBuddy/2026-08-19-15-50-19"
e = pd.read_csv(f"{BASE}/results/elite_ranked_v6.csv")

pts = [
    ("cg27577527", "ZSCAN12", "green"),
    ("cg15790037", "ARL5C", "green"),
    ("cg23180938", "CDO1", "green"),
    ("cg16439198", "CYP1B1", "green"),
    ("cg18507379", "PCDHGA", "green"),
    ("cg10109500", "GHSR", "green"),
    ("cg07495363", "BOLL (partner)", "green"),
    ("cg24589459", "BOLL (anchor)", "yellow"),
    ("cg01268824", "ZNF154", "red"),
    ("cg05016408", "LOC134466", "red"),
    ("cg01580681", "HAND2", "red"),
]
rows = []
for probe, label, tier in pts:
    r = e[e["probe"] == probe].iloc[0]
    rows.append(dict(label=label, cervix=r["cervix_p95"], blood=r["blood_max_p95"], tier=tier))
df = pd.DataFrame(rows)

COL = {"green": "#2e9e5b", "yellow": "#e6a817", "red": "#c0392b"}
fig, ax = plt.subplots(figsize=(7.2, 5.2), dpi=200)
for _, r in df.iterrows():
    ax.scatter(r["cervix"], r["blood"], s=90, c=COL[r["tier"]],
               edgecolors="white", linewidths=1.0, zorder=3)
    dx, dy = 0.008, 0.008
    ha = "left"
    if r["label"] == "LOC134466":
        dx, dy = 0.008, -0.018
    if r["label"] == "ZNF154":
        dx, dy = 0.008, -0.018
    if r["label"] == "BOLL (partner)":
        dx, dy = 0.010, -0.018
    if r["label"] == "ZSCAN12":
        dx, dy = -0.008, 0.012; ha = "right"
    if r["label"] == "ARL5C":
        dx, dy = -0.008, -0.018; ha = "right"
    if r["label"] == "PCDHGA":
        dx, dy = 0.010, 0.012
    if r["label"] == "CYP1B1":
        dx, dy = 0.006, 0.016
    if r["label"] == "CDO1":
        dx, dy = -0.006, 0.014; ha = "right"
    ax.annotate(r["label"], (r["cervix"] + dx, r["blood"] + dy),
                fontsize=8.5, ha=ha, color="#222222")
ax.axvline(0.15, color="#888888", lw=1, ls="--", zorder=1)
ax.axhline(0.15, color="#888888", lw=1, ls="--", zorder=1)
ax.text(0.152, 0.635, "gate 0.15", fontsize=8, color="#666666")
ax.set_xlabel("Cervical-scrape background P95 β (GSE46306)", fontsize=10)
ax.set_ylabel("Blood-cell background max P95 β (GSE35069)", fontsize=10)
ax.set_xlim(-0.02, 0.30)
ax.set_ylim(-0.02, 0.66)
from matplotlib.lines import Line2D
handles = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor=COL["green"], markersize=9,
           label="clear pass (P95 ≤ 0.15 both compartments)"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor=COL["yellow"], markersize=9,
           label="borderline (point-estimate pass, bootstrap CI crosses gate)"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor=COL["red"], markersize=9,
           label="compartment-restricted (P95 > 0.15)"),
]
ax.legend(handles=handles, fontsize=8, loc="upper left", frameon=False)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
out = f"{BASE}/results/v99_figure3_compartment.png"
fig.savefig(out)
print("saved", out)
