"""Corrected external validation:
1) GSE136791: unsupervised clustering to find the hyperplasia subgroup, then carcinoma-only positivity.
2) GSE155760: metadata-verified uterine cancers (23 EEC + 10 USC) vs 36 endometrial controls."""
import csv
import json
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19")
DATA, RES = BASE / "data", BASE / "results"

panel = [r["probe"] for r in csv.DictReader(open(RES / "wetlab_panel_23_primer_annotation.csv", encoding="utf-8"))]
kit = ["cg23180938", "cg13495205", "cg03502002", "cg25390440"]
KEY = ["cg24589459", "cg07495363", "cg23180938", "cg15790037", "cg16439198",
       "cg18507379", "cg10109500", "cg13495205", "cg03502002", "cg25390440"]
GENES = {"cg24589459": "BOLL", "cg07495363": "BOLL(twin)", "cg23180938": "CDO1",
         "cg15790037": "ARL5C", "cg16439198": "CYP1B1", "cg18507379": "PCDHGA",
         "cg10109500": "GHSR", "cg13495205": "AJAP1", "cg03502002": "GALR1", "cg25390440": "GYPC"}
need = set(panel + kit)

def load(gse, probes=None):
    cols = {}
    for f in sorted((DATA / gse / "ewas_betas").glob("GSM*.txt")):
        vals = {}
        with open(f) as fh:
            for line in fh:
                pr, b = line.rstrip("\n").split("\t")
                if probes is None or pr in probes:
                    vals[pr] = float(b) if b != "NA" else np.nan
        cols[f.stem] = vals
    return pd.DataFrame(cols)

# ---------- GSE136791: clustering ----------
df136 = load("GSE136791")  # full probe set for clustering
print("GSE136791 matrix:", df136.shape, flush=True)
sub = df136.sample(n=80000, random_state=7).dropna()
X = (sub - sub.mean(axis=1).values[:, None]).values  # row-center
X = X / (sub.std(axis=1).values[:, None] + 1e-9)
# 简单 2-means（不用 sklearn，直接 numpy）
C = np.corrcoef(X.T)  # 样本间相关
# 用相关矩阵的第一主成分做二分
from numpy.linalg import eigh
w, V = eigh(C)
pc1 = V[:, -1]
med = np.median(pc1)
grp = np.where(pc1 >= med, 1, 0)
sizes = pd.Series(grp).value_counts().to_dict()
print("cluster sizes:", sizes, flush=True)
# 用我方高甲基化标志物判断哪个簇是增生（低信号簇）
cancer_markers = ["cg18507379", "cg10109500", "cg07495363", "cg23180938"]
mean_by_grp = {}
for g in [0, 1]:
    cols = df136.columns[grp == g]
    mean_by_grp[g] = df136.loc[cancer_markers, cols].mean().mean()
print("mean beta of cancer markers by cluster:", mean_by_grp, flush=True)
hyper_grp = min(mean_by_grp, key=mean_by_grp.get)  # 信号低的簇=增生
hyp_cols = list(df136.columns[grp == hyper_grp])
ca_cols = list(df136.columns[grp != hyper_grp])
print(f"hyperplasia-like cluster n={len(hyp_cols)}, carcinoma n={len(ca_cols)}", flush=True)

# ---------- GSE155760: metadata split ----------
def get_meta(gse, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(
                f"https://ngdc.cncb.ac.cn/ewas/datahub/repository/basic?field=project%20id&val={gse}&relationship=&limit=500",
                headers={"User-Agent": "Mozilla/5.0"})
            return json.loads(urllib.request.urlopen(req, timeout=90).read())["content"]
        except Exception:
            time.sleep(3)
    return []
meta155 = {r["sample id"]: r for r in get_meta("GSE155760")}
df155 = load("GSE155760", need)
eec = [c for c in df155.columns if meta155.get(c, {}).get("disease") == "endometrial endometrioid carcinoma"]
usc = [c for c in df155.columns if meta155.get(c, {}).get("disease") == "uterine serous carcinoma"]
ctrl155 = [c for c in df155.columns if meta155.get(c, {}).get("sample type") == "control"]
print(f"GSE155760: EEC={len(eec)}, USC={len(usc)}, controls={len(ctrl155)}", flush=True)

df136k = df136.loc[df136.index.intersection(need)]
rows = []
for p in KEY:
    def stat(dframe, cols):
        if p not in dframe.index or not cols:
            return np.nan
        v = dframe.loc[p, cols].dropna()
        return round(float((v > 0.3).mean()), 3) if len(v) else np.nan
    def ctrl_p95(dframe, cols):
        if p not in dframe.index or not cols:
            return np.nan
        v = dframe.loc[p, cols].dropna()
        return round(float(v.quantile(0.95)), 3) if len(v) else np.nan
    rows.append({
        "probe": p, "gene": GENES[p],
        "GSE136791_all103": stat(df136k, list(df136k.columns)),
        "GSE136791_ca_only": stat(df136k, ca_cols),
        "GSE136791_hyp_like": stat(df136k, hyp_cols),
        "GSE155760_EEC23": stat(df155, eec),
        "GSE155760_USC10": stat(df155, usc),
        "GSE155760_ctrl36_p95": ctrl_p95(df155, ctrl155),
    })
out = pd.DataFrame(rows)
out.to_csv(RES / "external_validation_corrected.csv", index=False)
pd.set_option("display.width", 220)
print(out.to_string(index=False), flush=True)
