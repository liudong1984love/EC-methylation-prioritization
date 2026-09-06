"""Download per-sample normalized beta txt files from NGDC EWAS Data Hub (NCBI-free route)."""
import re
import subprocess
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

BASE = "https://download.cncb.ac.cn/ewas/datahub/EWAS_db"
DATA = Path(r"C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\data")
WORKERS = 10

def list_txt(gse):
    url = f"{BASE}/{gse}/"
    html = ""
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                html = r.read().decode("utf-8", errors="replace")
            break
        except Exception as e:
            print(f"  list retry {attempt+1}: {e}", flush=True)
            time.sleep(5)
    return sorted(set(re.findall(r'href="(GSM\d+\.txt)"', html)))

def fetch(url, dest):
    if dest.exists() and dest.stat().st_size > 1_000_000:
        return True
    tmp = dest.with_suffix(".part")
    for attempt in range(4):
        try:
            subprocess.run(["curl", "-sL", "--max-time", "600", "--retry", "2",
                            "-o", str(tmp), url], check=True, timeout=660)
            if tmp.exists() and tmp.stat().st_size > 1_000_000:
                tmp.replace(dest)
                return True
        except Exception as e:
            print(f"  retry {attempt+1} {dest.name}: {e}", flush=True)
        time.sleep(3)
    if tmp.exists():
        tmp.unlink()
    return False

def main():
    for gse in sys.argv[1:]:
        outdir = DATA / gse / "ewas_betas"
        outdir.mkdir(parents=True, exist_ok=True)
        files = list_txt(gse)
        print(f"=== {gse}: {len(files)} samples ===", flush=True)
        if not files:
            continue
        ok = fail = 0
        with ThreadPoolExecutor(max_workers=WORKERS) as ex:
            futs = {ex.submit(fetch, f"{BASE}/{gse}/{fn}", outdir / fn): fn for fn in files}
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
