"""Figure 3 (v9.17): compartment suitability, colour categories identical to Table 7.
Green: max(cervix,blood) P95 <= 0.15; yellow: 0.15-0.30; red: >0.30. Asterisk: CI-crossing borderline.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# gene: (cervix P95, blood P95)  [manuscript Table 7 values]
data = [
    ("ZSCAN12", 0.044, 0.058),
    ("ARL5C", 0.041, 0.052),
    ("CDO1", 0.055, 0.080),
    ("GHSR", 0.072, 0.149),
    ("PCDHGA", 0.080, 0.116),
    ("CYP1B1", 0.108, 0.122),
    ("BOLL*", 0.141, 0.112),
    ("AJAP1", 0.088, 0.295),
    ("GALR1", 0.226, 0.435),
    ("GYPC", 0.680, 0.838),
]

def cat(v):
    return "#2ca02c" if v <= 0.15 else ("#e6b800" if v <= 0.30 else "#d62728")

genes = [d[0] for d in data]
cervix = [d[1] for d in data]
blood = [d[2] for d in data]
limit = [max(c, b) for c, b in zip(cervix, blood)]
order = np.argsort(limit)
genes = [genes[i] for i in order]; cervix = [cervix[i] for i in order]; blood = [blood[i] for i in order]

x = np.arange(len(genes)); w = 0.38
fig, ax = plt.subplots(figsize=(8.5, 4), dpi=300)
b1 = ax.bar(x - w / 2, cervix, w, label="Cervix (HPV-neg, n=20)")
b2 = ax.bar(x + w / 2, blood, w, label="Blood (purified, n=42)")
for rects in (b1, b2):
    for r, g in zip(rects, genes):
        v = r.get_height()
        r.set_color(cat(v))
        r.set_alpha(0.92)
for xi, g in zip(x, genes):
    if g.endswith("*"):
        ax.text(xi, max(cervix[genes.index(g)], blood[genes.index(g)]) + 0.02, "*", ha="center", fontsize=14, fontweight="bold")
ax.axhline(0.15, ls="--", lw=1, color="grey"); ax.text(len(genes) - 0.4, 0.16, "0.15 gate", fontsize=8, color="grey", ha="right")
ax.axhline(0.30, ls=":", lw=1, color="grey"); ax.text(len(genes) - 0.4, 0.31, "0.30", fontsize=8, color="grey", ha="right")
ax.set_xticks(x); ax.set_xticklabels([g.rstrip("*") for g in genes], rotation=30, ha="right", fontsize=9)
ax.set_ylabel("Background P95 β"); ax.set_ylim(0, 0.95)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(fc="#2ca02c", label="≤0.15 (✓)"), Patch(fc="#e6b800", label="0.15–0.30 (△)"), Patch(fc="#d62728", label=">0.30 (✗)")],
          fontsize=8, loc="upper left", title="Colour = Table 7 category (per bar)", title_fontsize=8)
fig.tight_layout()
fig.savefig(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\results\v917_figure3_compartment.png", bbox_inches="tight")
print("saved results/v917_figure3_compartment.png")
