"""GSE90060 pair-fingerprint diagnostics: within-pair r, best cross-donor r, age concordance."""
import csv
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19")
meta90 = {r["sample id"]: r for r in csv.DictReader(open(BASE / "results/GSE90060_ewas_meta.csv", encoding="utf-8"))}
files90 = sorted((BASE / "data/GSE90060/ewas_betas").glob("GSM*.txt"))
df90 = pd.DataFrame({f.stem: pd.to_numeric(pd.read_csv(f, sep="\t", names=["probe", "beta"], index_col=0, dtype={"beta": str})["beta"], errors="coerce") for f in files90})
early = [c for c in df90.columns if meta90[c]["menstrual cycle phase"] == "early secretory"]
mid = [c for c in df90.columns if meta90[c]["menstrual cycle phase"] == "mid secretory"]
print(f"early={len(early)}, mid={len(mid)}", flush=True)
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

within = [C[early.index(e), mid.index(m)] for e, m in pairs.items()]
# best remaining cross-donor correlation: max off-pair entry
mask = np.ones_like(C, dtype=bool)
for e, m in pairs.items():
    mask[early.index(e), mid.index(m)] = False
cross = C[mask]
print(f"pairs={len(pairs)}; within-pair r: min={min(within):.4f}, median={float(np.median(within)):.4f}, max={max(within):.4f}", flush=True)
print(f"best cross-donor r (off-pair max)={cross.max():.4f}; off-pair median={float(np.median(cross)):.4f}", flush=True)
# age concordance
bad = [(e, m, meta90[e]["age (year)"], meta90[m]["age (year)"]) for e, m in pairs.items()
       if abs(float(meta90[e]["age (year)"]) - float(meta90[m]["age (year)"])) > 1.0]
print(f"age-mismatched pairs (>1 yr): {len(bad)}", flush=True)
for b in bad:
    print("  ", b, flush=True)
ages = sorted(set(float(meta90[g]["age (year)"]) for g in early + mid))
print("distinct ages:", ages, flush=True)
