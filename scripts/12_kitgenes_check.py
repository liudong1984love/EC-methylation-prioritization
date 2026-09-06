"""Five-dimension triage metrics for the approved-kit genes (CDO1/AJAP1/GALR1) + GYPC probes."""
import csv
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19")
DATA, RES = BASE / "data", BASE / "results"

PROBES = {
    "cg23180938": "CDO1(获批kit)",
    "cg13495205": "AJAP1(获批kit,topΔβ)",
    "cg11835068": "AJAP1(获批kit,低背景)",
    "cg15484532": "AJAP1(获批kit)",
    "cg03502002": "GALR1(获批kit,topΔβ)",
    "cg04534765": "GALR1(获批kit)",
    "cg10390058": "GALR1(获批kit)",
    "cg25390440": "GYPC(WID-qEC)",
    "cg04453971": "GYPC(WID-qEC)",
}
need = set(PROBES)

def scan_dir(path):
    cols = {}
    for f in sorted(Path(path).glob("GSM*.txt")):
        vals = {}
        with open(f) as fh:
            for line in fh:
                pr, b = line.rstrip("\n").split("\t")
                if pr in need:
                    vals[pr] = float(b) if b != "NA" else np.nan
        cols[f.stem] = vals
    return pd.DataFrame(cols)

def p95(v):
    v = pd.Series(v).dropna()
    return float(v.quantile(0.95)) if len(v) >= 3 else np.nan

def mean(v):
    v = pd.Series(v).dropna()
    return float(v.mean()) if len(v) else np.nan

# GSE73949 (育龄正常内膜, noob 450K rds 需 R 转换——这里用 EWAS? GSE73949 不在 EWAS 下载集,用 R 导出的 csv)
# 预检: 若无 csv 则用 R 先导出
ne17_csv = RES / "GSE73949_betas_export.csv"
if not ne17_csv.exists():
    import subprocess
    subprocess.run(["C:/Program Files/R/R-4.6.0/bin/Rscript.exe", "-e",
        "b <- readRDS('C:/Users/ld/WorkBuddy/2026-08-19-15-50-19/results/GSE73949_betas.rds'); "
        "write.csv(b, 'C:/Users/ld/WorkBuddy/2026-08-19-15-50-19/results/GSE73949_betas_export.csv')"],
        check=True, env={"PATH": "C:/Program Files/R/R-4.6.0/bin", "R_LIBS_USER": "C:/Users/ld/Rlibs"})
ne17 = pd.read_csv(ne17_csv, index_col=0)
ne17.index.name = "probe"

meta22 = {r["sample id"]: r for r in csv.DictReader(open(RES / "GSE223817_ewas_meta.csv", encoding="utf-8"))}
df22 = scan_dir(DATA / "GSE223817/ewas_betas")
ctrl22 = [c for c in df22.columns if meta22[c]["sample type"] == "control"]
eis22 = [c for c in df22.columns if meta22[c]["sample type"] == "disease tissue"]

meta46 = {r["sample id"]: r for r in csv.DictReader(open(RES / "GSE46306_ewas_meta.csv", encoding="utf-8"))}
df46 = scan_dir(DATA / "GSE46306/ewas_betas")
norm46 = [c for c in df46.columns if meta46[c].get("infection") == "HPV-" and meta46[c].get("sample type") == "control"]

meta35 = {r["sample id"]: r for r in csv.DictReader(open(RES / "GSE35069_ewas_meta.csv", encoding="utf-8"))}
df35 = scan_dir(DATA / "GSE35069/ewas_betas")
blood_groups = {}
for c in df35.columns:
    blood_groups.setdefault(meta35[c]["tissue"], []).append(c)

meta90 = {r["sample id"]: r for r in csv.DictReader(open(RES / "GSE90060_ewas_meta.csv", encoding="utf-8"))}
# 指纹配对需要全探针矩阵（不做目标过滤）
files90 = sorted((DATA / "GSE90060/ewas_betas").glob("GSM*.txt"))
df90full = pd.DataFrame({f.stem: pd.to_numeric(pd.read_csv(f, sep="\t", names=["probe", "beta"], index_col=0, dtype={"beta": str})["beta"], errors="coerce") for f in files90})
df90 = df90full
early = [c for c in df90.columns if meta90[c]["menstrual cycle phase"] == "early secretory"]
mid = [c for c in df90.columns if meta90[c]["menstrual cycle phase"] == "mid secretory"]
sub90 = df90.sample(n=min(120000, len(df90)), random_state=1)
E = sub90[early].to_numpy(); M = sub90[mid].to_numpy()
def colcorr(A, B):
    An = A - np.nanmean(A, axis=0); Bn = B - np.nanmean(B, axis=0)
    num = np.nansum(An[:, :, None] * Bn[:, None, :], axis=0)
    den = np.sqrt(np.nansum(An**2, axis=0))[:, None] * np.sqrt(np.nansum(Bn**2, axis=0))[None, :]
    return num / den
C = colcorr(E, M)
pairs, used = {}, set()
for i, j in np.dstack(np.unravel_index(np.argsort(-C, axis=None), C.shape))[0]:
    e, m = early[i], mid[j]
    if e not in pairs and m not in used:
        pairs[e] = m; used.add(m)

rows = []
for p, label in PROBES.items():
    r = {"probe": p, "gene": label}
    r["ne17_p95"] = round(p95(ne17.loc[p]) if p in ne17.index else np.nan, 3)
    r["ne347_p95"] = round(p95(df22.loc[p, ctrl22]) if p in df22.index else np.nan, 3)
    r["eis637_p95"] = round(p95(df22.loc[p, eis22]) if p in df22.index else np.nan, 3)
    r["cervix_p95"] = round(p95(df46.loc[p, norm46]) if p in df46.index else np.nan, 3)
    if p in df35.index:
        bt = {t: p95(df35.loc[p, cols]) for t, cols in blood_groups.items()}
        bt = {k: v for k, v in bt.items() if not np.isnan(v)}
        r["blood_max_p95"] = round(max(bt.values()), 3) if bt else np.nan
        r["blood_worst"] = max(bt, key=bt.get) if bt else ""
    if p in df90.index:
        deltas = [df90.loc[p, m] - df90.loc[p, e] for e, m in pairs.items()
                  if pd.notna(df90.loc[p, e]) and pd.notna(df90.loc[p, m])]
        r["cycle_mean_abs"] = round(float(np.mean(np.abs(deltas))), 3) if deltas else np.nan
    rows.append(r)

odf = pd.DataFrame(rows)
odf.to_csv(RES / "kitgenes_fivedim_check.csv", index=False)
pd.set_option("display.width", 250)
print(odf.to_string(index=False))
