"""Review-round-2 (Li Yunhui comments) analyses:
A) Menstrual-cycle drift for all final candidates (GSE90060, fingerprint-paired early vs mid secretory)  -> boll_cycle_drift_v95.csv
B) Hyperplasia->carcinoma quantitative gradient (GSE136791 EPIC: 34 hyperplasia vs 69 carcinoma vs 13 endometrial controls from GSE155760;
   GSE67116 450K: 8 hyperplasia vs 33 primary EC) -> hyperplasia_gradient_v95.csv
"""
import csv
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

BASE = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19")
DATA, RES = BASE / "data", BASE / "results"
DL = Path(r"C:/Users/ld/Downloads")

panel_rows = list(csv.DictReader(open(RES / "wetlab_panel_23_primer_annotation.csv", encoding="utf-8")))
PANEL = [r["probe"] for r in panel_rows]
GENE = {r["probe"]: r["gene"].split(";")[0] for r in panel_rows}
EXTRA = {"cg23180938": "CDO1(kit)", "cg07495363": "BOLL(twin)"}
PROBES = list(dict.fromkeys(PANEL + list(EXTRA)))
for p, g in EXTRA.items():
    GENE.setdefault(p, g)

def p95(x):
    x = pd.Series(x).dropna()
    return float(np.percentile(x, 95)) if len(x) else np.nan

def load_probes(gse, probes):
    cols = {}
    for f in sorted((DATA / gse / "ewas_betas").glob("GSM*.txt")):
        vals = {}
        with open(f) as fh:
            for line in fh:
                pr, b = line.rstrip("\n").split("\t")
                if pr in probes:
                    vals[pr] = float(b) if b != "NA" else np.nan
        cols[f.stem] = vals
    return pd.DataFrame(cols)

def read_list(fn):
    return [m.group(0) for line in open(DL / fn, encoding="utf-8") for m in [re.search(r"GSM\d+", line)] if m]

# ============ A) Cycle drift (GSE90060) ============
print("== A) cycle drift ==", flush=True)
meta90 = {r["sample id"]: r for r in csv.DictReader(open(RES / "GSE90060_ewas_meta.csv", encoding="utf-8"))}
files90 = sorted((DATA / "GSE90060/ewas_betas").glob("GSM*.txt"))
df90 = pd.DataFrame({f.stem: pd.to_numeric(pd.read_csv(f, sep="\t", names=["probe", "beta"], index_col=0, dtype={"beta": str})["beta"], errors="coerce") for f in files90})
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
print(f"paired {len(pairs)} early/mid samples", flush=True)

rows = []
for p in PROBES:
    if p not in df90.index:
        rows.append({"probe": p, "gene": GENE.get(p, ""), "cycle_drift": np.nan, "early_mean": np.nan, "mid_mean": np.nan, "note": "probe absent on 450K"})
        continue
    deltas = [df90.loc[p, m] - df90.loc[p, e] for e, m in pairs.items()
              if pd.notna(df90.loc[p, e]) and pd.notna(df90.loc[p, m])]
    rows.append({"probe": p, "gene": GENE.get(p, ""),
                 "cycle_drift": round(float(np.mean(np.abs(deltas))), 3),
                 "early_mean": round(float(df90.loc[p, early].mean()), 3),
                 "mid_mean": round(float(df90.loc[p, mid].mean()), 3), "note": ""})
cd = pd.DataFrame(rows).sort_values("cycle_drift", na_position="last")
cd.to_csv(RES / "boll_cycle_drift_v95.csv", index=False)
print(cd.to_string(index=False), flush=True)

# ============ B) Hyperplasia gradient ============
print("== B) hyperplasia gradient ==", flush=True)
ca69 = read_list("NoteGPT_GSE136791_carcinoma_GSM_list.txt")
hyp34 = read_list("NoteGPT_GSE136791_hyperplasia_GSM_list.txt")
ctrl13 = read_list("NoteGPT_GSE155760_endometrial_controls_GSM_list.txt")
print(f"GSE136791 ca={len(ca69)} hyp={len(hyp34)}; GSE155760 ctrl={len(ctrl13)}", flush=True)

EPIC_ABSENT = {"cg27577527", "cg20998319"}
EPIC_PROBES = [p for p in PROBES if p not in EPIC_ABSENT]
df136 = load_probes("GSE136791", set(EPIC_PROBES))
df155 = load_probes("GSE155760", set(EPIC_PROBES))

grad_rows = []
for p in EPIC_PROBES:
    if p not in df136.index:
        continue
    v_ca = df136.loc[p, [c for c in ca69 if c in df136.columns]].dropna()
    v_hyp = df136.loc[p, [c for c in hyp34 if c in df136.columns]].dropna()
    v_ct = df155.loc[p, [c for c in ctrl13 if c in df155.columns]].dropna() if p in df155.index else pd.Series(dtype=float)
    w_hyp_ca = stats.mannwhitneyu(v_hyp, v_ca, alternative="less").pvalue if len(v_hyp) and len(v_ca) else np.nan
    grad_rows.append({
        "probe": p, "gene": GENE.get(p, ""),
        "ctrl_mean": round(v_ct.mean(), 3) if len(v_ct) else np.nan,
        "hyp_mean": round(v_hyp.mean(), 3), "ca_mean": round(v_ca.mean(), 3),
        "hyp_pos03": round(float((v_hyp > 0.3).mean()), 3), "ca_pos03": round(float((v_ca > 0.3).mean()), 3),
        "gradient_ctrl_hyp_ca": bool(len(v_ct) and v_ct.mean() < v_hyp.mean() < v_ca.mean()),
        "wilcoxon_hyp_lt_ca_p": round(w_hyp_ca, 4) if pd.notna(w_hyp_ca) else np.nan,
    })
g1 = pd.DataFrame(grad_rows)
print("GSE136791 gradient:", flush=True)
print(g1.to_string(index=False), flush=True)

# GSE67116 (450K): betas from rds via pyreadr
g2 = pd.DataFrame()
try:
    import pyreadr
    b67116 = pyreadr.read_r(str(RES / "GSE67116_betas.rds"))
    df67 = list(b67116.values())[0]
    df67.index = df67.index.astype(str)
    pheno = pd.read_csv(RES / "GSE67116_pheno.csv")
    hyp67 = [g for g in pheno.loc[pheno.sample_type == "hyperplasia", "GSM"] if g in df67.columns]
    ec67 = [g for g in pheno.loc[pheno.sample_type.str.contains("primary|cancer|tumou?r", case=False, na=False), "GSM"] if g in df67.columns]
    if not ec67:
        print("GSE67116 sample_type values:", pheno.sample_type.unique(), flush=True)
        ec67 = [g for g in pheno.loc[pheno.sample_type != "hyperplasia", "GSM"] if g in df67.columns]
    rows67 = []
    for p in PROBES:
        if p not in df67.index:
            continue
        v_h = pd.to_numeric(df67.loc[p, hyp67], errors="coerce").dropna()
        v_c = pd.to_numeric(df67.loc[p, ec67], errors="coerce").dropna()
        if not len(v_h) or not len(v_c):
            continue
        rows67.append({"probe": p, "gene": GENE.get(p, ""),
                       "hyp_mean": round(v_h.mean(), 3), "ec_mean": round(v_c.mean(), 3),
                       "hyp_pos03": round(float((v_h > 0.3).mean()), 3), "ec_pos03": round(float((v_c > 0.3).mean()), 3),
                       "gradient_hyp_lt_ec": bool(v_h.mean() < v_c.mean())})
    g2 = pd.DataFrame(rows67)
    print(f"GSE67116 gradient (hyp n={len(hyp67)}, EC n={len(ec67)}):", flush=True)
    print(g2.to_string(index=False), flush=True)
except Exception as ex:
    print("GSE67116 rds read failed:", ex, flush=True)

with pd.ExcelWriter(RES / "hyperplasia_gradient_v95.xlsx") as xw:
    g1.to_excel(xw, sheet_name="GSE136791_EPIC", index=False)
    if len(g2):
        g2.to_excel(xw, sheet_name="GSE67116_450K", index=False)
g1.to_csv(RES / "hyperplasia_gradient_GSE136791_v95.csv", index=False)
if len(g2):
    g2.to_csv(RES / "hyperplasia_gradient_GSE67116_v95.csv", index=False)
print("DONE", flush=True)
