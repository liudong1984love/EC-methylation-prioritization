"""v9.17 verification bundle:
1) OR-rule joint positivity on manuscript-defined background subsets
2) Table 8 ctrl13 P95 recomputation for all listed probes
3) Cycle drift on age-concordant pairs (sensitivity)
"""
import csv
import re
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19")
DL = Path(r"C:/Users/ld/Downloads")
RES = BASE / "results"
BOLL2 = {"cg24589459", "cg07495363"}

def read_sample(path, probes):
    vals = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            pr, b = line.rstrip("\n").split("\t")[:2]
            if pr in probes and b not in ("NA", "", "NaN"):
                vals[pr] = float(b)
    return vals

def cohort_df(gse, gsm_list, probes):
    rows = {}
    for g in gsm_list:
        f = BASE / "data" / gse / "ewas_betas" / f"{g}.txt"
        if f.exists():
            rows[g] = read_sample(f, probes)
    return pd.DataFrame.from_dict(rows, orient="index")

def or_report(df, label, out):
    df = df.dropna(subset=["cg24589459", "cg07495363"])
    n = len(df)
    a = int((df["cg24589459"] > 0.3).sum())
    p = int((df["cg07495363"] > 0.3).sum())
    o = int(((df["cg24589459"] > 0.3) | (df["cg07495363"] > 0.3)).sum())
    print(f"{label}: n={n} anchor>0.3:{a} partner>0.3:{p} OR>0.3:{o} ({o/n*100:.1f}%)", flush=True)
    out.append(dict(cohort=label, n=n, anchor_pos=a, partner_pos=p, or_pos=o))

out = []
ctrl13 = [re.search(r"GSM\d+", l).group(0) for l in open(DL / "NoteGPT_GSE155760_endometrial_controls_GSM_list.txt", encoding="utf-8") if re.search(r"GSM\d+", l)]
or_report(cohort_df("GSE155760", ctrl13, BOLL2), "GSE155760_ctrl13", out)

g49 = [f.stem for f in (BASE / "data/GSE73949/ewas_betas").glob("GSM*.txt")]
or_report(cohort_df("GSE73949", g49, BOLL2), "GSE73949_healthy17", out)

meta817 = pd.read_csv(RES / "GSE223817_ewas_meta.csv")
ctrl347 = meta817.loc[meta817["sample type"] == "control", "sample id"].tolist()
endo = meta817.loc[meta817["sample type"] != "control", "sample id"].tolist()
print(f"GSE223817 ctrl={len(ctrl347)} endo={len(endo)}", flush=True)
or_report(cohort_df("GSE223817", ctrl347, BOLL2), "GSE223817_ctrl347", out)
or_report(cohort_df("GSE223817", endo, BOLL2), "GSE223817_endometriosis", out)

meta306 = pd.read_csv(RES / "GSE46306_ewas_meta.csv")
hpvneg = meta306.loc[(meta306["infection"] == "HPV-") & (meta306["sample type"] == "control"), "sample id"].tolist()
or_report(cohort_df("GSE46306", hpvneg, BOLL2), "GSE46306_cervix_HPVneg", out)

meta069 = pd.read_csv(RES / "GSE35069_ewas_meta.csv")
pur = meta069.loc[~meta069["tissue"].isin(["whole blood", "peripheral blood mononuclear cell", "granulocyte"]), "sample id"].tolist()
or_report(cohort_df("GSE35069", pur, BOLL2), "GSE35069_blood_purified", out)

pd.DataFrame(out).to_csv(RES / "v917_boll_or_rule_background.csv", index=False)
print("saved v917_boll_or_rule_background.csv", flush=True)

# ---- 2) Table 8 ctrl13 P95 for all probes listed ----
t8_probes = ["cg18507379", "cg10109500", "cg15790037", "cg23180938", "cg07495363", "cg24589459", "cg16439198"]
d13 = cohort_df("GSE155760", ctrl13, set(t8_probes))
rows = []
for p in t8_probes:
    x = d13[p].dropna().values
    rows.append(dict(probe=p, n=len(x), p95=round(float(np.quantile(x, 0.95, method="linear")), 4),
                     max=round(float(np.max(x)), 4), prop_gt_0_15=round(float(np.mean(x > 0.15)), 4),
                     n_gt_0_3=int((x > 0.3).sum())))
r = pd.DataFrame(rows)
print(r.to_string(index=False), flush=True)
r.to_csv(RES / "v917_table8_ctrl13_recompute.csv", index=False)

# ---- 3) cycle drift on age-concordant pairs ----
meta90 = {row["sample id"]: row for row in csv.DictReader(open(RES / "GSE90060_ewas_meta.csv", encoding="utf-8"))}
panel = pd.read_csv(RES / "wetlab_panel_23_primer_annotation.csv")["probe"].tolist()
probe_set = set(panel)
files90 = sorted((BASE / "data/GSE90060/ewas_betas").glob("GSM*.txt"))
df90 = {}
for f in files90:
    df90[f.stem] = read_sample(f, probe_set)
df90 = pd.DataFrame.from_dict(df90, orient="index").T

# reuse pairing from 70b logic but on 120k genome-wide sample requires full matrix; instead reuse saved pairs via age rule:
# reconstruct pairs with the same greedy algorithm requires full matrix; here we reuse the pairs printed by 70b via re-run on probe subset is invalid.
# So: hardcode nothing — rerun full pairing quickly using all probe rows from per-sample files is heavy; instead load pairing from a dump if exists.
# Simpler: recompute pairing using the panel probes is NOT comparable. -> Do full re-run with sampled probes (same as 70b).
def full_matrix():
    cols = {}
    for f in files90:
        cols[f.stem] = pd.to_numeric(pd.read_csv(f, sep="\t", names=["probe", "beta"], index_col=0, dtype={"beta": str})["beta"], errors="coerce")
    return pd.DataFrame(cols)

M = full_matrix()
early = [c for c in M.columns if meta90[c]["menstrual cycle phase"] == "early secretory"]
mid = [c for c in M.columns if meta90[c]["menstrual cycle phase"] == "mid secretory"]
sub = M.sample(n=min(120000, len(M)), random_state=1)
E = sub[early].to_numpy(); MM = sub[mid].to_numpy()
def colcorr(A, B):
    An = A - np.nanmean(A, axis=0); Bn = B - np.nanmean(B, axis=0)
    num = np.nansum(An[:, :, None] * Bn[:, None, :], axis=0)
    den = np.sqrt(np.nansum(An**2, axis=0))[:, None] * np.sqrt(np.nansum(Bn**2, axis=0))[None, :]
    return num / den
C = colcorr(E, MM)
pairs, used = {}, set()
for i, j in np.dstack(np.unravel_index(np.argsort(-C, axis=None), C.shape))[0]:
    e, m = early[i], mid[j]
    if e not in pairs and m not in used:
        pairs[e] = m; used.add(m)
concord = {e: m for e, m in pairs.items()
           if abs(float(meta90[e]["age (year)"]) - float(meta90[m]["age (year)"])) <= 1.0}
print(f"pairs={len(pairs)}, age-concordant={len(concord)}", flush=True)

drift = []
for p in panel:
    if p not in df90.index:
        continue
    d_all = [abs(df90.loc[p, m] - df90.loc[p, e]) for e, m in pairs.items()
             if pd.notna(df90.loc[p, e]) and pd.notna(df90.loc[p, m])]
    d_con = [abs(df90.loc[p, m] - df90.loc[p, e]) for e, m in concord.items()
             if pd.notna(df90.loc[p, e]) and pd.notna(df90.loc[p, m])]
    drift.append(dict(probe=p,
                      drift_all17=round(float(np.mean(d_all)), 4) if d_all else np.nan,
                      drift_concordant=round(float(np.mean(d_con)), 4) if d_con else np.nan,
                      n_concordant=len(d_con)))
dd = pd.DataFrame(drift)
dd.to_csv(RES / "v917_cycle_drift_concordant.csv", index=False)
key = dd[dd.probe.isin(["cg24589459", "cg27577527", "cg07495363"])]
print(key.to_string(index=False), flush=True)
print("pool range (concordant):", round(dd.drift_concordant.min(), 4), "-", round(dd.drift_concordant.max(), 4), flush=True)
