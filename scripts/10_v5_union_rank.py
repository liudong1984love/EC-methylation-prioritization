"""v5: ne347 background for all elite probes + union-gate re-ranking."""
import csv
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19")
DATA, RES = BASE / "data", BASE / "results"

elite = list(csv.DictReader(open(RES / "elite_probes_v2.csv", encoding="utf-8")))
allp = [r["probe"] for r in elite]
meta = {r["sample id"]: r for r in csv.DictReader(open(RES / "GSE223817_ewas_meta.csv", encoding="utf-8"))}
ctrl = [g for g, r in meta.items()
        if r["sample type"] == "control" and (DATA / "GSE223817/ewas_betas" / f"{g}.txt").exists()]

need = set(allp)
cols = {}
for g in ctrl:
    vals = {}
    with open(DATA / "GSE223817/ewas_betas" / f"{g}.txt") as f:
        for line in f:
            pr, b = line.rstrip("\n").split("\t")
            if pr in need:
                vals[pr] = float(b) if b != "NA" else np.nan
    cols[g] = vals
df = pd.DataFrame(cols)
print("matrix:", df.shape, flush=True)

tri = {r["probe"]: r for r in csv.DictReader(open(RES / "triage_closedloop_v4.csv", encoding="utf-8"))}
out = []
for r in elite:
    p = r["probe"]
    rec = dict(r)
    if p in df.index:
        v = df.loc[p].dropna()
        rec["ne347_n"] = len(v)
        rec["ne347_mean"] = round(float(v.mean()), 3) if len(v) else ""
        rec["ne347_p95"] = round(float(v.quantile(0.95)), 3) if len(v) else ""
    else:
        rec["ne347_n"] = 0; rec["ne347_mean"] = ""; rec["ne347_p95"] = ""
    t = tri.get(p, {})
    for k in ["cycle_mean_abs_d", "cervix_p95", "blood_max_p95", "triage_pass"]:
        rec[k] = t.get(k, "")
    ne17 = float(r["ne_bg_p95"]) if r["ne_bg_p95"] not in ("", None) else np.nan
    ne347 = float(rec["ne347_p95"]) if rec["ne347_p95"] != "" else np.nan
    rec["bg_worst"] = round(np.nanmax([ne17, ne347]), 3)
    tri_ok = str(t.get("triage_pass", "")).upper() == "TRUE"
    rec["union_pass"] = bool(tri_ok and not np.isnan(ne347) and ne347 <= 0.15)
    # v5 评分：低背景项用最坏背景
    s_bg = max(0, 1 - rec["bg_worst"] / 0.15) if not np.isnan(rec["bg_worst"]) else 0
    rec["score_v5"] = round(0.25 * s_bg + 0.25 * float(r["pre_pos03"]) +
                            0.20 * float(r["case_pos03"]) + 0.15 * 1 + 0.05 * 1, 3)
    out.append(rec)

odf = pd.DataFrame(out).sort_values("score_v5", ascending=False)
odf.to_csv(RES / "elite_ranked_v5.csv", index=False)
up = odf[odf["union_pass"] == True]
print(f"elite total={len(odf)}, union_pass={len(up)}", flush=True)
show = ["probe", "gene", "delta_beta", "ne_bg_p95", "ne347_p95", "bg_worst",
        "cervix_p95", "blood_max_p95", "pre_pos03", "union_pass", "score_v5"]
pd.set_option("display.width", 250)
print(odf.head(25)[show].to_string(index=False), flush=True)
