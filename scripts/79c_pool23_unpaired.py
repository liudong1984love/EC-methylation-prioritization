"""Unpaired early-vs-mid diffs for the exact 23-probe experimental pool (from boll_cycle_drift_v95.csv)."""
import csv
from pathlib import Path

import pandas as pd

BASE = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19")
pool = [r["probe"] for r in csv.DictReader(open(BASE / "results/boll_cycle_drift_v95.csv", encoding="utf-8"))]
print("pool n =", len(pool))
meta = {r["sample id"]: r for r in csv.DictReader(open(BASE / "results/GSE90060_ewas_meta.csv", encoding="utf-8"))}
files = sorted((BASE / "data/GSE90060/ewas_betas").glob("GSM*.txt"))
vals = {p: {"early secretory": [], "mid secretory": []} for p in pool}
for f in files:
    phase = meta[f.stem]["menstrual cycle phase"]
    d = pd.read_csv(f, sep="\t", names=["probe", "beta"], index_col=0)
    d["beta"] = pd.to_numeric(d["beta"], errors="coerce")
    for p in pool:
        if p in d.index:
            vals[p][phase].append(float(d.loc[p, "beta"]))
rows = []
for p in pool:
    e, m = vals[p]["early secretory"], vals[p]["mid secretory"]
    rows.append((p, abs(sum(e) / len(e) - sum(m) / len(m))))
res = pd.DataFrame(rows, columns=["probe", "unpaired_absdiff"]).sort_values("unpaired_absdiff")
res.to_csv(BASE / "results/v918_cycle_unpaired_pool23.csv", index=False, float_format="%.4f")
print(res.to_string(index=False))
print(f"\npool23 range: {res.unpaired_absdiff.min():.4f}-{res.unpaired_absdiff.max():.4f}")
print(f"anchor cg24589459: {res[res.probe=='cg24589459'].unpaired_absdiff.iloc[0]:.4f}")
