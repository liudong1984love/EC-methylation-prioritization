# 36_substitute_fulldim.py
# Five-dimension background + TCGA profile for EPIC-compatible ZSCAN12-DMR
# substitute probes (cg25666433, cg20275132, cg23164203; reference cg27577527).
# All four are on HM450, so 450K-based cohorts (TCGA, GSE73949, GSE90060,
# GSE46306, GSE35069) are queryable; GSE223817/GSE136791/GSE155760 already in
# results/zscan12_epic_substitutes_v96.csv.
import gzip
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path("C:/Users/ld/WorkBuddy/2026-08-19-15-50-19")
RES = BASE / "results"
PROBES = ["cg25666433", "cg20275132", "cg23164203", "cg27577527"]

rows = {p: {"probe": p} for p in PROBES}

# ---- TCGA-UCEC: tumour positivity (all + premenopausal), adjacent-normal P95
clin = pd.read_csv(BASE / "data/UCEC_clinicalMatrix.tsv", sep="\t", low_memory=False)
clin = clin.set_index(clin.columns[0])
# original pipeline definition: startswith("Pre") on full 15-char sample ID
pre_ids = set(clin.index[clin["menopause_status"].astype(str).str.startswith("Pre")])
hdr = None
betas = {}
with open(BASE / "data/TCGA_UCEC_methylation450_beta.tsv") as f:
    hdr = f.readline().rstrip("\n").split("\t")[1:]
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if parts[0] in PROBES:
            betas[parts[0]] = parts[1:]
    # file is probe-per-line; read fully
tum_cols = [i for i, s in enumerate(hdr) if s.endswith("-01")]
nor_cols = [i for i, s in enumerate(hdr) if s.endswith("-11")]
pre_cols = [i for i in tum_cols if hdr[i] in pre_ids or hdr[i][:-4] in pre_ids]
print(f"TCGA: tumour={len(tum_cols)}, normal={len(nor_cols)}, premenopausal tumour={len(pre_cols)}")
for p, vals in betas.items():
    v = pd.to_numeric(pd.Series(vals), errors="coerce")
    vt = v.iloc[tum_cols].dropna(); vn = v.iloc[nor_cols].dropna(); vp = v.iloc[pre_cols].dropna()
    rows[p].update({
        "tcga_tumour_pos03": round(float((vt > 0.3).mean()), 3),
        "tcga_tumour_mean": round(vt.mean(), 3),
        "tcga_pre_pos03": round(float((vp > 0.3).mean()), 3) if len(vp) else np.nan,
        "tcga_pre_n": len(vp),
        "tcga_adjnormal_p95": round(np.percentile(vn, 95), 3) if len(vn) else np.nan,
        "tcga_adjnormal_mean": round(vn.mean(), 3) if len(vn) else np.nan,
    })

# ---- GSE73949 healthy premenopausal endometrium (full 450K betas rds)
import pyreadr
g = list(pyreadr.read_r(str(RES / "GSE73949_betas.rds")).values())[0]
g.index = g.index.astype(str)
for p in PROBES:
    if p in g.index:
        v = pd.to_numeric(g.loc[p], errors="coerce").dropna()
        rows[p].update({"gse73949_mean": round(v.mean(), 3),
                        "gse73949_p95": round(np.percentile(v, 95), 3),
                        "gse73949_n": len(v)})

# ---- ewas_betas cohorts: GSE90060 (cycle), GSE46306 (cervix), GSE35069 (blood)
def ewas_load(folder, gsms):
    out = {}
    for gsm in gsms:
        fp = BASE / "data" / folder / "ewas_betas" / f"{gsm}.txt"
        if not fp.exists():
            continue
        d = {}
        with open(fp) as fh:
            for line in fh:
                p, val = line.rstrip("\n").split("\t")
                if p in PROBES:
                    d[p] = float(val) if val not in ("", "NA", "nan") else np.nan
        out[gsm] = d
    return pd.DataFrame(out)

# GSE90060: use metadata to get all samples (drift uses pairs; here background level)
m90 = pd.read_csv(RES / "GSE90060_ewas_meta.csv")
col = [c for c in m90.columns if "sample" in c.lower() and "id" in c.lower()][0]
b90 = ewas_load("GSE90060", m90[col].tolist())
for p in PROBES:
    if p in b90.index:
        v = pd.to_numeric(b90.loc[p], errors="coerce").dropna()
        rows[p].update({"gse90060_mean": round(v.mean(), 3), "gse90060_p95": round(np.percentile(v, 95), 3), "gse90060_n": len(v)})

m46 = pd.read_csv(RES / "GSE46306_ewas_meta.csv")
col46 = [c for c in m46.columns if "sample" in c.lower() and "id" in c.lower()][0]
# manuscript definition: HPV-negative normal cervical scrapes (n=20)
norm46 = m46.loc[(m46["sample type"] == "control") & (m46["infection"] == "HPV-"), col46].tolist()
print(f"GSE46306 HPV- controls: {len(norm46)}")
b46 = ewas_load("GSE46306", norm46)
for p in PROBES:
    if p in b46.index:
        v = pd.to_numeric(b46.loc[p], errors="coerce").dropna()
        rows[p].update({"gse46306_cervix_mean": round(v.mean(), 3), "gse46306_cervix_p95": round(np.percentile(v, 95), 3), "gse46306_n": len(v)})

m35 = pd.read_csv(RES / "GSE35069_ewas_meta.csv")
col35 = [c for c in m35.columns if "sample" in c.lower() and "id" in c.lower()][0]
b35 = ewas_load("GSE35069", m35[col35].tolist())
# per-cell-type P95, then max across cell types (matches triage definition)
for p in PROBES:
    if p in b35.index:
        v = pd.to_numeric(b35.loc[p], errors="coerce")
        per_type = []
        for tis, grp in m35.groupby("tissue"):
            vv = v[[g for g in grp[col35] if g in b35.columns]].dropna()
            if len(vv):
                per_type.append(np.percentile(vv, 95))
        allv = v.dropna()
        rows[p].update({"gse35069_blood_max_p95": round(max(per_type), 3) if per_type else np.nan,
                        "gse35069_blood_mean": round(allv.mean(), 3), "gse35069_n": len(allv)})

out = pd.DataFrame(rows.values())
out.to_csv(RES / "zscan12_substitutes_fulldim_v96.csv", index=False)
print(out.to_string(index=False))
