"""Parallel GEO downloader via EBI BioStudies mirror using curl subprocesses (urllib proved unstable)."""
import json
import subprocess
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

BASE = "https://www.ebi.ac.uk/biostudies"
DATA = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\data")

DATASETS = {
    "GSE67116": "E-GEOD-67116",
    "GSE90060": "E-GEOD-90060",
    "GSE73949": "E-GEOD-73949",
    "GSE46306": "E-GEOD-46306",
    "GSE35069": "E-GEOD-35069",
}
WORKERS = 10

def list_files(acc):
    files, offset = [], 0
    while True:
        url = f"{BASE}/api/v1/studies/{acc}/files?offset={offset}&limit=25"
        js = None
        for attempt in range(5):
            try:
                with urllib.request.urlopen(url, timeout=60) as r:
                    js = json.loads(r.read().decode("utf-8"))
                break
            except Exception as e:
                print(f"  list retry {attempt+1}: {e}", flush=True)
                time.sleep(5)
        if js is None:
            raise RuntimeError(f"cannot list {acc}")
        items = js.get("items", [])
        files.extend(items)
        total = js.get("pagination", {}).get("total", len(files))
        offset += len(items)
        if not items or offset >= total:
            return files

def fetch(url, dest, expected_size):
    if dest.exists() and expected_size and dest.stat().st_size == expected_size:
        return True
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(4):
        try:
            subprocess.run(
                ["curl", "-sL", "--max-time", "600", "--retry", "2", "-o", str(tmp), url],
                check=True, timeout=660)
            if tmp.exists() and (not expected_size or tmp.stat().st_size == expected_size):
                tmp.replace(dest)
                return True
        except Exception as e:
            print(f"  retry {attempt+1} {dest.name}: {e}", flush=True)
        time.sleep(3)
    if tmp.exists():
        tmp.unlink()
    return False

def main():
    only = sys.argv[1:] or list(DATASETS)
    for gse in only:
        acc = DATASETS[gse]
        outdir = DATA / gse
        outdir.mkdir(parents=True, exist_ok=True)
        print(f"=== {gse} ({acc}) ===", flush=True)
        try:
            files = [f for f in list_files(acc)
                     if f.get("isDirectory") == "false" and f["path"].lower().endswith(".idat")]
        except RuntimeError as e:
            print(f"  SKIP: {e}", flush=True)
            continue
        print(f"  {len(files)} files", flush=True)
        ok = fail = 0
        with ThreadPoolExecutor(max_workers=WORKERS) as ex:
            futs = {}
            for f in files:
                name = Path(f["path"]).name
                url = f"{BASE}/files/{acc}/{f['path']}"
                size = int(f.get("Size") or 0)
                futs[ex.submit(fetch, url, outdir / name, size)] = name
            for fut in as_completed(futs):
                if fut.result():
                    ok += 1
                else:
                    fail += 1
                if (ok + fail) % 20 == 0:
                    print(f"  progress {ok+fail}/{len(files)}", flush=True)
        print(f"  done: {ok} ok, {fail} failed", flush=True)
    print("ALL DONE", flush=True)

if __name__ == "__main__":
    main()
