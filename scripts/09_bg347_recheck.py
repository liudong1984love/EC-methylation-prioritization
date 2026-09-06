"""Large-scale benign background recheck with GSE223817 controls (n=347, cycle-phased)."""
import csv
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19")
DATA, RES = BASE / "data", BASE / "results"

tri = list(csv.DictReader(open(RES / "triage_closedloop_v4.csv", encoding="utf-8")))
POI = [r["probe"] for r in tri]
meta = {r["sample id"]: r for r in csv.DictReader(open(RES / "GSE223817_ewas_meta.csv", encoding="utf-8"))}
ctrl = [g for g, r in meta.items() if r["sample type"] == "control"]
ctrl = [g for g in ctrl if (DATA / "GSE223817/ewas_betas" / f"{g}.txt").exists()]
print(f"controls loaded: {len(ctrl)}", flush=True)

cols = {}
for g in ctrl:
    s = pd.read_csv(DATA / "GSE223817/ewas_betas" / f"{g}.txt", sep="\t",
                    names=["probe", "beta"], index_col=0, dtype={"beta": str})
    cols[g] = pd.to_numeric(s["beta"], errors="coerce")
df = pd.DataFrame(cols)
print("matrix:", df.shape, flush=True)

phases = {g: meta[g].get("menstrual cycle phases", "") for g in ctrl}
rows = []
for p in POI:
    if p not in df.index:
        rows.append({"probe": p, "ne347_n": 0})
        continue
    v = df.loc[p]
    rec = {"probe": p, "ne347_n": int(v.notna().sum()),
           "ne347_mean": float(v.mean()), "ne347_p95": float(v.quantile(0.95)),
           "ne347_max": float(v.max())}
    ph = {}
    for g in ctrl:
        ph.setdefault(phases[g], []).append(v[g])
    for phase, vals in ph.items():
        vals = pd.Series(vals).dropna()
        if len(vals) >= 5:
            rec[f"mean_{phase.replace(' ', '_')}"] = round(float(vals.mean()), 3)
    pv = [rec[k] for k in rec if k.startswith("mean_")]
    rec["phase_max_mean"] = max(pv) if pv else np.nan
    rec["phase_min_mean"] = min(pv) if pv else np.nan
    rec["phase_swing"] = rec["phase_max_mean"] - rec["phase_min_mean"] if pv else np.nan
    rec["ne347_pass"] = rec["ne347_p95"] <= 0.15
    rows.append(rec)

odf = pd.DataFrame(rows)
tri_df = pd.DataFrame(tri).set_index("probe")
odf = odf.merge(tri_df[["triage_pass"]], left_on="probe", right_index=True, how="left")
odf = odf.sort_values("ne347_p95")
odf.to_csv(RES / "bg347_recheck_v5.csv", index=False)
sel = ["probe", "ne347_n", "ne347_mean", "ne347_p95", "phase_swing", "ne347_pass", "triage_pass"]
pd.set_option("display.width", 200)
print(odf[sel].to_string(index=False), flush=True)
