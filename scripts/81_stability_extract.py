"""Stability analysis step 1: extract TCGA + background betas for the 476 audited
single-probe candidates and VALIDATE against the audit CSV (S2a).
Caches matrices to results/stability_cache.npz for the bootstrap step.
"""
import os, glob, sys
import numpy as np
import pandas as pd

BASE = r"C:/Users/ld/WorkBuddy/2026-08-19-15-50-19"
os.chdir(BASE)
OUT = r"C:/Users/ld/WorkBuddy/2026-08-19-15-50-19/results"  # results dir of the analysis base

audit = pd.read_csv("Supplementary_Data_S2a_single_probe_audit.csv")
PROBES = audit["probe"].tolist()
PSET = set(PROBES)
print("audit probes:", len(PROBES), flush=True)

# ---------------- TCGA ----------------
print("scanning TCGA beta tsv ...", flush=True)
with open("data/TCGA_UCEC_methylation450_beta.tsv") as fh:
    header = fh.readline().rstrip("\n").split("\t")
samples = header[1:]
rows = {}
with open("data/TCGA_UCEC_methylation450_beta.tsv") as fh:
    fh.readline()
    for line in fh:
        pid = line[:10]
        if pid in PSET:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in PSET:
                rows[parts[0]] = [float(x) if x not in ("NA", "", "NaN") else np.nan for x in parts[1:]]
print("TCGA probes found:", len(rows), "samples:", len(samples), flush=True)
tcga = pd.DataFrame(rows, index=samples).T.reindex(PROBES)

tumor = np.array([s.endswith("-01") for s in samples])
normal = np.array([s.endswith("-11") for s in samples])
clin = pd.read_csv("data/UCEC_clinicalMatrix.tsv", sep="\t", low_memory=False)
meno = dict(zip(clin["sampleID"].str[:15], clin["menopause_status"]))
is_pre = np.array([str(meno.get(s[:15], "")).startswith("Pre") for s in samples])
pre_tumor = tumor & is_pre
print(f"tumors={tumor.sum()} normals={normal.sum()} pre_tumors={pre_tumor.sum()}", flush=True)

X = tcga.values  # probes x samples
delta_beta = np.nanmean(X[:, tumor], 1) - np.nanmean(X[:, normal], 1)
case_pos03 = np.nanmean(X[:, tumor] > 0.3, 1)
pre_pos03 = np.nanmean(X[:, pre_tumor] > 0.3, 1)
ctrl_p95 = np.nanquantile(np.where(normal[None, :], X, np.nan), 0.95, axis=1, method="linear")

val = pd.DataFrame({"probe": PROBES, "delta_beta": delta_beta, "case_pos03": case_pos03,
                    "pre_pos03": pre_pos03, "ctrl_p95": ctrl_p95}).merge(
    audit[["probe", "delta_beta", "case_pos03", "pre_pos03"]], on="probe", suffixes=("_new", "_audit"))
for c in ["delta_beta", "case_pos03", "pre_pos03"]:
    d = (val[f"{c}_new"] - val[f"{c}_audit"]).abs()
    print(f"VALIDATE {c}: max|diff|={d.max():.4f} mean|diff|={d.mean():.5f} n>{0.01}={int((d>0.01).sum())}", flush=True)
print("ctrl_p95>0.40 count:", int((ctrl_p95 > 0.40).sum()), flush=True)

# ---------------- background cohorts ----------------
def read_probes(path, pset):
    out = {}
    with open(path) as fh:
        for line in fh:
            p = line[:10]
            if p in pset:
                pr, b = line.rstrip("\n").split("\t")[:2]
                if b not in ("NA", "", "NaN"):
                    out[pr] = float(b)
    return out

def cohort_matrix(gsms, cdir, tag):
    recs = {}
    for i, g in enumerate(gsms):
        fp = os.path.join(cdir, g + ".txt")
        if os.path.exists(fp):
            recs[g] = read_probes(fp, PSET)
        if (i + 1) % 100 == 0:
            print(f"  {tag}: {i+1}/{len(gsms)}", flush=True)
    m = pd.DataFrame(recs, index=PROBES)  # probes x samples (NaN absent)
    print(tag, m.shape, "mean probes/sample:", int(m.notna().sum(0).mean()), flush=True)
    return m

# gate1 GSE73949: audit used the minfi-noob reprocessed betas (full-genome RDS;
# the CSV export only covers 9 panel probes)
import pyreadr
_g1full = list(pyreadr.read_r("results/GSE73949_betas.rds").values())[0]
g1 = _g1full.reindex(PROBES)
print("GSE73949 noob RDS", g1.shape, "matched:", int(g1.notna().any(axis=1).sum()), flush=True)
# gate2 GSE223817 controls
m223 = pd.read_csv("results/GSE223817_ewas_meta.csv")
ctrl223 = m223.loc[m223["sample type"] == "control", "sample id"].tolist()
print("GSE223817 controls:", len(ctrl223), flush=True)
g2 = cohort_matrix(ctrl223, "data/GSE223817/ewas_betas", "GSE223817")
# gate3 GSE46306 HPV-neg controls
m463 = pd.read_csv("results/GSE46306_ewas_meta.csv")
hpvneg = m463.loc[(m463["infection"] == "HPV-") & (m463["sample type"] == "control"), "sample id"].tolist()
print("GSE46306 HPV-neg:", len(hpvneg), flush=True)
g3 = cohort_matrix(hpvneg, "data/GSE46306/ewas_betas", "GSE46306")
# gate4 GSE35069: audit used ALL 60 samples (10 fractions incl. granulocyte/PBMC/whole blood),
# per-fraction P95 (n=6) then max across fractions
m350 = pd.read_csv("results/GSE35069_ewas_meta.csv")
pur_meta = m350
print("GSE35069 all:", len(pur_meta))
print(pur_meta["tissue"].value_counts(), flush=True)
g4 = cohort_matrix(pur_meta["sample id"].tolist(), "data/GSE35069/ewas_betas", "GSE35069")
frac = pur_meta.set_index("sample id")["tissue"]

# ---------------- validate gates vs audit ----------------
def p95(mat):
    return np.nanquantile(mat.values.astype(float), 0.95, axis=1, method="linear")

gate_vals = {"gate1_healthy_endometrium_P95": p95(g1),
             "gate2_benign_surgical_P95": p95(g2),
             "gate3_cervix_P95": p95(g3)}
# blood: pooled vs per-fraction max
blood_pooled = p95(g4)
pf = []
for t in frac.unique():
    cols = [c for c in g4.columns if frac.get(c) == t]
    if cols:
        pf.append(np.nanquantile(g4[cols].values.astype(float), 0.95, axis=1, method="linear"))
blood_fracmax = np.nanmax(np.vstack(pf), axis=0)
gate_vals["gate4_blood_P95_pooled"] = blood_pooled
gate_vals["gate4_blood_P95_fracmax"] = blood_fracmax

chk = pd.DataFrame({"probe": PROBES})
for k, v in gate_vals.items():
    chk[k] = v
chk = chk.merge(audit[["probe", "gate1_healthy_endometrium_P95", "gate2_benign_surgical_P95",
                       "gate3_cervix_P95", "gate4_blood_P95"]].rename(columns={
    "gate1_healthy_endometrium_P95": "audit_gate1", "gate2_benign_surgical_P95": "audit_gate2",
    "gate3_cervix_P95": "audit_gate3", "gate4_blood_P95": "audit_gate4"}), on="probe")
for new, old in [("gate1_healthy_endometrium_P95", "audit_gate1"),
                 ("gate2_benign_surgical_P95", "audit_gate2"),
                 ("gate3_cervix_P95", "audit_gate3"),
                 ("gate4_blood_P95_pooled", "audit_gate4"),
                 ("gate4_blood_P95_fracmax", "audit_gate4")]:
    d = (chk[new] - chk[old]).abs().dropna()
    print(f"VALIDATE {new} vs audit: n={len(d)} max|diff|={d.max():.4f} mean|diff|={d.mean():.5f} n>0.01={int((d>0.01).sum())}", flush=True)

np.savez_compressed(os.path.join(OUT, "stability_cache.npz"),
                    probes=np.array(PROBES), tcga=X.astype(np.float32),
                    tumor=tumor, normal=normal, pre_tumor=pre_tumor,
                    g1=g1.values.astype(np.float32), g2=g2.values.astype(np.float32),
                    g3=g3.values.astype(np.float32), g4=g4.values.astype(np.float32),
                    g4_frac=frac.reindex(g4.columns).astype(str).values,
                    audit=audit.to_records(index=False))
chk.to_csv(os.path.join(OUT, "stability_validation.csv"), index=False)
print("cache saved", flush=True)
