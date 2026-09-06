"""P0-1 confound decomposition + GSE93589 ZSCAN12 validation."""
import re
import subprocess
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19")
KEY = ["cg24589459", "cg07495363", "cg27577527", "cg23180938", "cg15790037",
       "cg16439198", "cg18507379", "cg10109500", "cg01268824", "cg18675097"]
need = set(KEY)

def dl(gse):
    url = f"https://download.cncb.ac.cn/ewas/datahub/EWAS_db/{gse}/"
    html = urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"}), timeout=60).read().decode("utf-8", "replace")
    files = sorted(set(re.findall(r'href="(GSM\d+)\.txt"', html)))
    outdir = BASE / "data" / gse / "ewas_betas"
    outdir.mkdir(parents=True, exist_ok=True)
    def fetch(gsm):
        dest = outdir / f"{gsm}.txt"
        if dest.exists() and dest.stat().st_size > 1_000_000:
            return True
        tmp = dest.with_suffix(".part")
        try:
            subprocess.run(["curl", "-sL", "--max-time", "600", "-o", str(tmp), f"{url}{gsm}.txt"], check=True, timeout=660)
            if tmp.exists() and tmp.stat().st_size > 1_000_000:
                tmp.replace(dest)
                return True
        except Exception:
            pass
        return False
    ok = 0
    with ThreadPoolExecutor(max_workers=10) as ex:
        for r in as_completed([ex.submit(fetch, g) for g in files]):
            ok += r.result()
    print(f"{gse}: {ok}/{len(files)}", flush=True)
    return files

f49 = dl("GSE73949")
f67 = dl("GSE67116")

def probe_vals(gse, ids):
    out = {p: [] for p in KEY}
    for g in ids:
        f = BASE / "data" / gse / "ewas_betas" / f"{g}.txt"
        if not f.exists():
            continue
        with open(f) as fh:
            for line in fh:
                pr, b = line.rstrip("\n").split("\t")
                if pr in need and b != "NA":
                    out[pr].append(float(b))
    return out

# GSE67116 分组：前 8 个为增生（按 GEO 系列顺序：hyperplasia 在前）——用元数据确认
import json, time
def meta(gse, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(
                f"https://ngdc.cncb.ac.cn/ewas/datahub/repository/basic?field=project%20id&val={gse}&relationship=&limit=500",
                headers={"User-Agent": "Mozilla/5.0"})
            return json.loads(urllib.request.urlopen(req, timeout=90).read())["content"]
        except Exception:
            time.sleep(3)
    return []
m67 = {r["sample id"]: r for r in meta("GSE67116")}
from collections import Counter
print("GSE67116 disease:", Counter(str(r.get("disease", ""))[:40] for r in m67.values()).most_common(6), flush=True)
hyp67 = [g for g in f67 if "hyperplasia" in str(m67.get(g, {}).get("disease", "")).lower()]
ca67 = [g for g in f67 if "adenocarcinoma" in str(m67.get(g, {}).get("disease", "")).lower() or "carcinoma" in str(m67.get(g, {}).get("disease", "")).lower()]
print(f"GSE67116: hyperplasia={len(hyp67)}, carcinoma={len(ca67)}", flush=True)

v49 = probe_vals("GSE73949", f49)          # 健康 450K GMQN
v67h = probe_vals("GSE67116", hyp67)       # 增生 450K GMQN
v67c = probe_vals("GSE67116", ca67)        # 癌 450K GMQN

# 已有数据（前次分析）：GSE223817 ctrl347、GSE155760 ctrl13、GSE73949-noob——从 v6 精英表取
v6 = pd.read_csv(BASE / "results/elite_ranked_v6.csv").set_index("probe")
ctrl13 = [l.strip() for l in open(r"C:/Users/ld/Downloads/NoteGPT_GSE155760_endometrial_controls_GSM_list.txt", encoding="utf-8")]
ctrl13 = [re.search(r"GSM\d+", l).group(0) for l in ctrl13 if re.search(r"GSM\d+", l)]
v155 = probe_vals("GSE155760", ctrl13)     # 无癌内膜 EPIC GMQN

def p95(v, p):
    s = v.get(p, [])
    return round(float(pd.Series(s).quantile(0.95)), 3) if s else np.nan
def pos(v, p):
    s = v.get(p, [])
    return round(float(np.mean(np.array(s) > 0.3)), 3) if s else np.nan

rows = []
for p in KEY:
    rows.append({
        "probe": p,
        "healthy_450K_GMQN_p95": p95(v49, p),
        "hyp_450K_GMQN_p95": p95(v67h, p),
        "ca_450K_GMQN_pos03": pos(v67c, p),
        "cancerfree_EPIC_GMQN_p95": p95(v155, p),
        "symptomatic_EPIC_GMQN_p95(v6)": (round(float(v6.loc[p, "ne347_p95"]), 3) if p in v6.index and pd.notna(v6.loc[p].get("ne347_p95", np.nan)) else np.nan),
        "healthy_450K_noob_p95(v6)": (round(float(v6.loc[p, "ne_bg_p95"]), 3) if p in v6.index and pd.notna(v6.loc[p].get("ne_bg_p95", np.nan)) else np.nan),
    })
df = pd.DataFrame(rows)
df.to_csv(BASE / "results/confound_decomposition.csv", index=False)
pd.set_option("display.width", 240)
print(df.to_string(index=False), flush=True)

# GSE93589 ZSCAN12 等阳性率（450K，9 例 EC）
f93 = sorted((BASE / "data/GSE93589/ewas_betas").glob("GSM*.txt"))
v93 = probe_vals("GSE93589", [f.stem for f in f93])
print("\nGSE93589 (9 EC, 450K) positivity:", flush=True)
for p in KEY:
    print(f"  {p}: pos03={pos(v93, p)}, mean={round(float(np.mean(v93[p])),3) if v93[p] else 'NA'}", flush=True)
