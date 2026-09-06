# 35_zscan12_epic_substitutes.py
# Find EPIC-compatible substitute probes in the ZSCAN12 DMR and benchmark them
# on the decisive background dimension (GSE223817 symptomatic controls) and on
# external tumour cohorts (GSE136791 carcinomas, GSE155760 EEC + endometrial controls).
import csv, gzip
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path("C:/Users/ld/WorkBuddy/2026-08-19-15-50-19")
RES = BASE / "results"
DL = Path("C:/Users/ld/Downloads")

# hg38 DMR window: cg27577527 hg38 = chr6:28,399,766 (HM450 manifest);
# hg19 DMR 28,367,279-28,367,898 -> hg38 approx 28,399,501-28,400,120 (+32,222 shift)
DMR_LO, DMR_HI = 28399300, 28400400  # slightly widened

# --- EPIC probes in region (already extracted: 13); re-derive with coords ---
region = []
with gzip.open(BASE / "data/EPIC.hg38.manifest.tsv.gz", "rt") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    for line in f:
        c = line.rstrip("\n").split("\t")
        if c[0] == "chr6" and DMR_LO <= int(c[1]) <= DMR_HI:
            region.append({"probe": c[8], "pos": int(c[1]), "mapq": c[12]})
region = pd.DataFrame(region).sort_values("pos")
print(f"EPIC probes in window: {len(region)}")

# which are also on HM450 (checkable in TCGA/450K cohorts)
hm450 = set()
with gzip.open(BASE / "data/HM450.hg38.manifest.tsv.gz", "rt") as f:
    f.readline()
    for line in f:
        hm450.add(line.split("\t")[8])
region["on_450K"] = region.probe.isin(hm450)
print(region.to_string(index=False))

probes = region.probe.tolist()

def load_betas(folder, gsms, probes):
    """Return DataFrame probe x gsm of betas for given sample files."""
    out = {}
    for g in gsms:
        fp = BASE / "data" / folder / "ewas_betas" / f"{g}.txt"
        if not fp.exists():
            continue
        d = {}
        with open(fp) as fh:
            for line in fh:
                p, v = line.rstrip("\n").split("\t")
                if p in probes:
                    d[p] = float(v) if v not in ("", "NA", "nan") else np.nan
        out[g] = d
    return pd.DataFrame(out)

# --- GSE223817 symptomatic controls ---
meta = pd.read_csv(RES / "GSE223817_ewas_meta.csv")
ctrl = meta.loc[meta["sample type"] == "control", "sample id"].tolist()
eis = meta.loc[meta["sample type"] == "disease tissue", "sample id"].tolist()
print(f"GSE223817 controls={len(ctrl)}, endometriosis={len(eis)}")
b_ctrl = load_betas("GSE223817", ctrl, probes)
b_eis = load_betas("GSE223817", eis, probes)
print(f"loaded: ctrl {b_ctrl.shape}, eis {b_eis.shape}")

# --- GSE136791 carcinomas (GEO-annotated list) ---
ca136 = [l.split("\t")[0].strip() for l in open(DL / "NoteGPT_GSE136791_carcinoma_GSM_list.txt", encoding="utf-8", errors="ignore") if l.strip().startswith("GSM")]
hyp136 = [l.split("\t")[0].strip() for l in open(DL / "NoteGPT_GSE136791_hyperplasia_GSM_list.txt", encoding="utf-8", errors="ignore") if l.strip().startswith("GSM")]
b_ca = load_betas("GSE136791", ca136, probes)
b_hyp = load_betas("GSE136791", hyp136, probes)
print(f"GSE136791: ca {b_ca.shape}, hyp {b_hyp.shape}")

# --- GSE155760: GEO-annotated EEC and endometrial-control lists ---
def read_gsm_list(fp):
    gsms = []
    for l in open(fp, encoding="utf-8", errors="ignore"):
        l = l.strip()
        if l.startswith("GSM"):
            gsms.append(l.split("\t")[0].split()[0])
    return gsms

eec155 = read_gsm_list(DL / "NoteGPT_GSE155760_EEC_GSM_list.txt")
ctrl155 = read_gsm_list(DL / "NoteGPT_GSE155760_endometrial_controls_GSM_list.txt")
b_eec = load_betas("GSE155760", eec155, probes)
b_c155 = load_betas("GSE155760", ctrl155, probes)
print(f"GSE155760: eec {b_eec.shape}, ctrl {b_c155.shape}")

rows = []
for p in probes:
    r = {"probe": p, "pos": int(region.loc[region.probe == p, "pos"].iloc[0]),
         "on_450K": bool(region.loc[region.probe == p, "on_450K"].iloc[0])}
    for name, df in [("gse223817_ctrl", b_ctrl), ("gse223817_eis", b_eis),
                     ("gse136791_ca", b_ca), ("gse136791_hyp", b_hyp),
                     ("gse155760_eec", b_eec), ("gse155760_ctrl", b_c155)]:
        if p in df.index and df.shape[1] > 0:
            v = pd.to_numeric(df.loc[p], errors="coerce").dropna()
            r[f"{name}_n"] = len(v)
            r[f"{name}_mean"] = round(v.mean(), 3)
            r[f"{name}_p95"] = round(np.percentile(v, 95), 3) if len(v) else np.nan
            r[f"{name}_pos03"] = round(float((v > 0.3).mean()), 3)
        else:
            r[f"{name}_n"] = 0
    rows.append(r)
out = pd.DataFrame(rows)
out.to_csv(RES / "zscan12_epic_substitutes_v96.csv", index=False)
print(out.to_string(index=False))
