"""v6: endometriosis (n=637) benign-disease background for all elite probes; merge into final table."""
import csv
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19")
DATA, RES = BASE / "data", BASE / "results"

elite = pd.read_csv(RES / "elite_ranked_v5.csv")
need = set(elite["probe"])
meta = {r["sample id"]: r for r in csv.DictReader(open(RES / "GSE223817_ewas_meta.csv", encoding="utf-8"))}
dis = [g for g, r in meta.items()
       if r["sample type"] == "disease tissue" and (DATA / "GSE223817/ewas_betas" / f"{g}.txt").exists()]
print(f"endometriosis samples: {len(dis)}", flush=True)

cols = {}
for g in dis:
    vals = {}
    with open(DATA / "GSE223817/ewas_betas" / f"{g}.txt") as f:
        for line in f:
            pr, b = line.rstrip("\n").split("\t")
            if pr in need:
                vals[pr] = float(b) if b != "NA" else np.nan
    cols[g] = vals
df = pd.DataFrame(cols)
print("matrix:", df.shape, flush=True)

eis_mean = df.mean(axis=1).round(3)
eis_p95 = df.quantile(0.95, axis=1).round(3)
elite["eis637_n"] = elite["probe"].map(df.notna().sum(axis=1)).fillna(0).astype(int)
elite["eis637_mean"] = elite["probe"].map(eis_mean)
elite["eis637_p95"] = elite["probe"].map(eis_p95)

# 三维良性背景最坏值（健康育龄 ne17 / 有症状对照 ne347 / 内异症 eis637）
def worst(row):
    vals = []
    for k in ["ne_bg_p95", "ne347_p95", "eis637_p95"]:
        try:
            v = float(row[k])
            if not np.isnan(v):
                vals.append(v)
        except (TypeError, ValueError):
            pass
    return round(max(vals), 3) if vals else np.nan
elite["bg3_worst"] = elite.apply(worst, axis=1)
elite["tri_ok"] = elite["triage_pass"].astype(str).str.upper() == "TRUE"
elite["union3_pass"] = elite.apply(
    lambda r: bool(r["tri_ok"]) and not np.isnan(r["bg3_worst"]) and r["bg3_worst"] <= 0.15, axis=1)
s_bg = elite["bg3_worst"].apply(lambda w: max(0, 1 - w / 0.15) if not np.isnan(w) else 0)
elite["score_v6"] = (0.25 * s_bg + 0.25 * elite["pre_pos03"].astype(float) +
                     0.20 * elite["case_pos03"].astype(float) + 0.15 * 1 + 0.05 * 1).round(3)
elite = elite.sort_values("score_v6", ascending=False)
elite.to_csv(RES / "elite_ranked_v6.csv", index=False)

up = elite[elite["union3_pass"]]
print(f"elite={len(elite)}, union3_pass={len(up)}", flush=True)
panel7 = ["cg24589459", "cg27577527", "cg15790037", "cg16439198", "cg23180938", "cg18507379", "cg10109500"]
show = ["probe", "gene", "delta_beta", "ne_bg_p95", "ne347_p95", "eis637_p95", "bg3_worst",
        "cervix_p95", "blood_max_p95", "union3_pass", "score_v6"]
pd.set_option("display.width", 260)
print(elite[elite["probe"].isin(panel7)][show].to_string(index=False), flush=True)
print("\n--- union3_pass probes ---")
print(up[["probe", "gene", "delta_beta", "bg3_worst", "pre_pos03", "score_v6"]].to_string(index=False), flush=True)
