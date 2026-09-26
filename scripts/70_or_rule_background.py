"""OR-rule joint positivity of BOLL two-probe rule across all benign background cohorts.
Early-exit file scan (both probes found -> stop). GMQN per-sample files from EWAS hub.
"""
import re
from pathlib import Path

BASE = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19")
DL = Path(r"C:/Users/ld/Downloads")
PROBES = {"cg07495363", "cg24589459"}

def load(gse, ids=None):
    d = BASE / "data" / gse / "ewas_betas"
    if not d.exists():
        return None
    files = sorted(d.glob("GSM*.txt"))
    if ids is not None:
        keep = set(ids)
        files = [f for f in files if f.stem in keep]
    out = {}
    for f in files:
        found = {}
        with open(f) as fh:
            for line in fh:
                pr, b = line.rstrip("\n").split("\t")
                if pr in PROBES:
                    if b != "NA":
                        found[pr] = float(b)
                    if len(found) == 2 or (pr == "cg24589459" and len(found) >= 1):
                        if len(found) == 2:
                            break
        out[f.stem] = found
    return out

def or_stats(v, label):
    if v is None:
        print(f"{label}: NO DATA", flush=True)
        return None
    ids = [g for g, d in v.items() if "cg24589459" in d and "cg07495363" in d]
    n = len(ids)
    b1 = sum(1 for g in ids if v[g]["cg24589459"] > 0.3)
    b2 = sum(1 for g in ids if v[g]["cg07495363"] > 0.3)
    either = sum(1 for g in ids if v[g]["cg24589459"] > 0.3 or v[g]["cg07495363"] > 0.3)
    print(f"{label}: n={n}, anchor>0.3: {b1}, partner>0.3: {b2}, OR>0.3: {either} ({either/n*100:.2f}%)", flush=True)
    return label, n, b1, b2, either

ctrl13 = [re.search(r"GSM\d+", l).group(0) for l in open(DL / "NoteGPT_GSE155760_endometrial_controls_GSM_list.txt", encoding="utf-8") if re.search(r"GSM\d+", l)]

results = []
results.append(or_stats(load("GSE155760", ctrl13), "GSE155760_ctrl13_EPIC_GMQN"))
results.append(or_stats(load("GSE73949"), "GSE73949_healthy17_450K_GMQN"))
results.append(or_stats(load("GSE46306"), "GSE46306_cervix_450K_GMQN"))
results.append(or_stats(load("GSE35069"), "GSE35069_blood_450K_GMQN"))
results.append(or_stats(load("GSE223817"), "GSE223817_benign_EPIC_GMQN"))

import csv
with open(BASE / "results" / "v917_boll_or_rule_background.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["cohort", "n", "anchor_pos", "partner_pos", "or_rule_pos"])
    for r in results:
        if r:
            w.writerow(r)
print("saved results/v917_boll_or_rule_background.csv", flush=True)
