"""External validation of panel probes in GSE136791 (103 EC tumors, EPIC) and GSE155760 (EPIC, tumors+controls)."""
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
kit = ["cg23180938", "cg13495205", "cg03502002", "cg25390440"]  # CDO1, AJAP1, GALR1, GYPC
need = set(panel + kit)

def get_meta(gse, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(
                f"https://ngdc.cncb.ac.cn/ewas/datahub/repository/basic?field=project%20id&val={gse}&relationship=&limit=500",
                headers={"User-Agent": "Mozilla/5.0"})
            return json.loads(urllib.request.urlopen(req, timeout=90).read())["content"]
        except Exception as e:
            print("meta retry", i + 1, e); time.sleep(3)
    return []

def load(gse):
    cols = {}
    for f in sorted((DATA / gse / "ewas_betas").glob("GSM*.txt")):
        vals = {}
        with open(f) as fh:
            for line in fh:
                pr, b = line.rstrip("\n").split("\t")
                if pr in need:
                    vals[pr] = float(b) if b != "NA" else np.nan
        cols[f.stem] = vals
    return pd.DataFrame(cols)

out = []
for gse in ["GSE136791", "GSE155760"]:
    meta = {r["sample id"]: r for r in get_meta(gse)}
    df = load(gse)
    tumors = [c for c in df.columns if meta.get(c, {}).get("sample type") == "disease tissue"]
    ctrls = [c for c in df.columns if meta.get(c, {}).get("sample type") == "control"]
    if not tumors:  # metadata missing -> treat all as tumor for GSE136791
        if gse == "GSE136791":
            tumors = list(df.columns)
    print(f"{gse}: files={df.shape[1]}, tumors={len(tumors)}, controls={len(ctrls)}", flush=True)
    for p in need:
        if p not in df.index:
            out.append({"gse": gse, "probe": p, "n_tumor": 0, "pos03": np.nan,
                        "mean_beta": np.nan, "ctrl_p95": np.nan})
            continue
        vt = df.loc[p, tumors].dropna() if tumors else pd.Series(dtype=float)
        vc = df.loc[p, ctrls].dropna() if ctrls else pd.Series(dtype=float)
        out.append({"gse": gse, "probe": p, "n_tumor": len(vt),
                    "pos03": round(float((vt > 0.3).mean()), 3) if len(vt) else np.nan,
                    "mean_beta": round(float(vt.mean()), 3) if len(vt) else np.nan,
                    "ctrl_p95": round(float(vc.quantile(0.95)), 3) if len(vc) else np.nan})

odf = pd.DataFrame(out)
odf.to_csv(RES / "external_validation_v7.csv", index=False)
key = ["cg24589459", "cg27577527", "cg23180938", "cg15790037", "cg16439198", "cg18507379",
       "cg10109500", "cg07495363", "cg13495205", "cg03502002", "cg25390440"]
pd.set_option("display.width", 200)
print(odf[odf.probe.isin(key)].to_string(index=False), flush=True)
n_pass = odf[(odf.gse == "GSE136791") & (odf.pos03 >= 0.8)]
print(f"\nGSE136791 中 pos03>=0.8 的 panel 探针数: {len(n_pass)}/23")
