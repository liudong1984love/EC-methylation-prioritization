"""Triage closed loop v4: cycle stability (GSE90060, fingerprint-paired),
cervix background (GSE46306 HPV- normals), blood cell background (GSE35069)."""
import csv
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19")
DATA, RES = BASE / "data", BASE / "results"

# ---- probes of interest ----
elite = list(csv.DictReader(open(RES / "elite_probes_v2.csv", encoding="utf-8")))
passed = [r for r in elite if r["ne_pass"].upper() == "TRUE"]
top_elite = [r["probe"] for r in passed]  # 全部 65 个背景通过精英探针
panel = ["cg20935165", "cg01589129", "cg14768785", "cg15790037", "cg16439198", "cg23180938"]
five = ["cg05016408", "cg01268824", "cg18675097", "cg23180938", "cg01580681"]
POI = list(dict.fromkeys(panel + five + top_elite))
print(f"probes of interest: {len(POI)}", flush=True)

def load_betas(gse):
    files = sorted((DATA / gse / "ewas_betas").glob("GSM*.txt"))
    cols = {}
    for f in files:
        s = pd.read_csv(f, sep="\t", names=["probe", "beta"], index_col=0, dtype={"beta": str})
        cols[f.stem] = pd.to_numeric(s["beta"], errors="coerce")
    return pd.DataFrame(cols)

# ---- 1. GSE90060: fingerprint pairing + paired cycle deltas ----
meta = {r["sample id"]: r for r in csv.DictReader(open(RES / "GSE90060_ewas_meta.csv", encoding="utf-8"))}
df = load_betas("GSE90060")
early = [c for c in df.columns if meta[c]["menstrual cycle phase"] == "early secretory"]
mid = [c for c in df.columns if meta[c]["menstrual cycle phase"] == "mid secretory"]
sub = df.sample(n=120000, random_state=1)
corr = sub[early].corrwith(sub[mid].set_axis(early, axis=1), axis=0)  # placeholder
# proper: cor between each early col and each mid col
E = sub[early].values; Mm = sub[mid].values
Ec = E - np.nanmean(E, axis=0); Mc = Mm - np.nanmean(Mm, axis=0)
def colcorr(A, B):
    An = A - np.nanmean(A, axis=0); Bn = B - np.nanmean(B, axis=0)
    num = np.nansum(An[:, :, None] * Bn[:, None, :], axis=0)
    den = np.sqrt(np.nansum(An**2, axis=0))[:, None] * np.sqrt(np.nansum(Bn**2, axis=0))[None, :]
    return num / den
C = colcorr(E, Mm)
pairs = {}
used_mid = set()
order = np.dstack(np.unravel_index(np.argsort(-C, axis=None), C.shape))[0]
for i, j in order:
    e, m = early[i], mid[j]
    if e not in pairs and m not in used_mid:
        pairs[e] = m; used_mid.add(m)
print(f"fingerprint pairs: {len(pairs)}/17; min pair corr={min(C[early.index(e), mid.index(m)] for e,m in pairs.items()):.4f}", flush=True)
within = [C[early.index(e), mid.index(m)] for e, m in pairs.items()]
C_off = C.copy()
for k, (e, m) in enumerate(pairs.items()):
    C_off[early.index(e), mid.index(m)] = np.nan
print(f"within-pair corr mean={np.mean(within):.4f} vs off-pair max={np.nanmax(C_off):.4f}", flush=True)

poi_idx = df.index.intersection(POI)
cyc = {}
for p in poi_idx:
    deltas = []
    for e, m in pairs.items():
        a, b = df.loc[p, e], df.loc[p, m]
        if pd.notna(a) and pd.notna(b):
            deltas.append(b - a)  # mid - early
    deltas = np.array(deltas)
    cyc[p] = dict(cycle_mean_abs_d=float(np.mean(np.abs(deltas))) if len(deltas) else np.nan,
                  cycle_max_abs_d=float(np.max(np.abs(deltas))) if len(deltas) else np.nan,
                  cycle_frac_gt01=float(np.mean(np.abs(deltas) > 0.1)) if len(deltas) else np.nan)

# ---- 2. GSE46306: HPV- normal cervix background ----
meta4 = {r["sample id"]: r for r in csv.DictReader(open(RES / "GSE46306_ewas_meta.csv", encoding="utf-8"))}
df4 = load_betas("GSE46306")
norm4 = [c for c in df4.columns if meta4[c].get("infection") == "HPV-" and meta4[c].get("sample type") == "control"]
print(f"cervix normals: {len(norm4)}", flush=True)
cerv = {}
for p in df4.index.intersection(POI):
    v = df4.loc[p, norm4].dropna()
    cerv[p] = dict(cervix_mean=float(v.mean()), cervix_p95=float(v.quantile(0.95)) if len(v) else np.nan)

# ---- 3. GSE35069: per-cell-type blood background ----
meta3 = {r["sample id"]: r for r in csv.DictReader(open(RES / "GSE35069_ewas_meta.csv", encoding="utf-8"))}
df3 = load_betas("GSE35069")
ctypes = {}
for c in df3.columns:
    ctypes.setdefault(meta3[c]["tissue"], []).append(c)
print("blood cell types:", {k: len(v) for k, v in ctypes.items()}, flush=True)
blood = {}
for p in df3.index.intersection(POI):
    per_type = {}
    for t, cols in ctypes.items():
        v = df3.loc[p, cols].dropna()
        if len(v):
            per_type[t] = float(v.quantile(0.95))
    blood[p] = dict(blood_max_p95=max(per_type.values()) if per_type else np.nan,
                    blood_worst_type=max(per_type, key=per_type.get) if per_type else "")

# ---- merge ----
out = []
for p in POI:
    row = {"probe": p}
    row.update(cyc.get(p, {})); row.update(cerv.get(p, {})); row.update(blood.get(p, {}))
    row["triage_pass"] = (row.get("cervix_p95", 1) <= 0.15 and row.get("blood_max_p95", 1) <= 0.15
                          and row.get("cycle_mean_abs_d", 1) <= 0.05 and row.get("cycle_frac_gt01", 1) <= 0.25)
    out.append(row)
odf = pd.DataFrame(out)
odf.to_csv(RES / "triage_closedloop_v4.csv", index=False)
pd.set_option("display.width", 250)
print(odf[["probe", "cycle_mean_abs_d", "cycle_frac_gt01", "cervix_p95", "blood_max_p95", "blood_worst_type", "triage_pass"]].to_string(index=False), flush=True)
