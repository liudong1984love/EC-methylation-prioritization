"""Fetch sample metadata from EBI BioSamples (mirrors NCBI) for EWAS-downloaded datasets."""
import json
import time
import urllib.request
import urllib.parse
import csv
from pathlib import Path

DATA = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\data")
OUT = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\results")
DATASETS = ["GSE90060", "GSE46306", "GSE35069"]

def query_gsm(gsm, retries=3):
    q = urllib.parse.quote(gsm)
    url = f"https://www.ebi.ac.uk/biosamples/samples?text={q}&size=25"
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                js = json.loads(r.read().decode("utf-8"))
            for s in js.get("_embedded", {}).get("samples", []):
                if s.get("name") == gsm:
                    chars = {}
                    for k, v in s.get("characteristics", {}).items():
                        if v and isinstance(v, list):
                            chars[k] = v[0].get("text", "")
                    return chars
            return None
        except Exception as e:
            print(f"  retry {i+1} {gsm}: {e}", flush=True)
            time.sleep(3)
    return None

def main():
    for gse in DATASETS:
        d = DATA / gse / "ewas_betas"
        gsms = sorted(f.stem for f in d.glob("GSM*.txt"))
        print(f"=== {gse}: {len(gsms)} samples ===", flush=True)
        rows = []
        for i, gsm in enumerate(gsms):
            chars = query_gsm(gsm)
            if chars is None:
                print(f"  MISS {gsm}", flush=True)
                chars = {}
            rows.append({"GSM": gsm, **chars})
            if (i + 1) % 20 == 0:
                print(f"  {i+1}/{len(gsms)}", flush=True)
            time.sleep(0.2)
        keys = ["GSM"]
        for r in rows:
            for k in r:
                if k not in keys:
                    keys.append(k)
        out = OUT / f"{gse}_metadata.csv"
        with open(out, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
        print(f"  saved {out}", flush=True)
    print("ALL DONE", flush=True)

if __name__ == "__main__":
    main()
