# -*- coding: utf-8 -*-
"""v9.7 reviewer-3 analyses:
A. TCGA pre- vs postmenopausal comparison + tissue x menopause interaction
B. TCGA paired tumour-normal, endometrioid subset, molecular-subtype sensitivity
C. Background-threshold retention scan (0.05-0.30) + curve figure
D. Bootstrap CIs for background P95 (lead probes, 4 cohorts) + Wilson CIs for positivity
E. Export chr2/chr6 M-value matrix for R bumphunter cross-check
"""
import os, subprocess, sys, gzip
import numpy as np
import pandas as pd
from scipy import stats

BASE = r"C:/Users/ld/WorkBuddy/2026-08-19-15-50-19"
D = os.path.join(BASE, "data"); R = os.path.join(BASE, "results")
os.makedirs(R, exist_ok=True)

# ---------- probe sets ----------
panel = pd.read_csv(os.path.join(R, "wetlab_panel_23_primer_annotation.csv"))
PANEL23 = panel["probe"].tolist()
SUBS = ["cg25666433", "cg23164203", "cg20275132"]           # ZSCAN12 EPIC substitutes (on 450K)
LEADS = ["cg24589459", "cg07495363", "cg27577527",          # BOLL x2, ZSCAN12
         "cg23180938", "cg10109500", "cg15790037"] + SUBS   # CDO1, GHSR, ARL5C + subs
FOCUS = sorted(set(PANEL23 + SUBS))

# ---------- manifest: chr2 / chr6 probe list for bumphunter ----------
man = pd.read_csv(os.path.join(D, "HM450.hg38.manifest.tsv.gz"), sep="\t", comment="#")
man = man[["Probe_ID", "CpG_chrm", "CpG_beg"]].dropna()
man.columns = ["probe", "chr", "pos"]
man["pos"] = pd.to_numeric(man["pos"], errors="coerce")
man = man.dropna()
chr26 = man[man["chr"].isin(["chr2", "chr6"])].copy()
KEEP = set(FOCUS) | set(chr26["probe"])
print(f"manifest: {len(man)} probes; chr2+chr6 {len(chr26)}; keep {len(KEEP)}")

# ---------- stream TCGA beta matrix ----------
rows = {}
with open(os.path.join(D, "TCGA_UCEC_methylation450_beta.tsv")) as f:
    header = f.readline().rstrip("\n").split("\t")
    samples = header[1:]
    for line in f:
        p = line.split("\t", 1)
        if p[0] in KEEP:
            rows[p[0]] = p[1].rstrip("\n").split("\t")
beta = pd.DataFrame.from_dict(rows, orient="index", columns=samples)
beta = beta.apply(pd.to_numeric, errors="coerce")
print("beta extracted:", beta.shape)
beta.to_csv(os.path.join(R, "v97_tcga_beta_focus.csv"))

# ---------- clinical ----------
clin = pd.read_csv(os.path.join(D, "UCEC_clinicalMatrix.tsv"), sep="\t", low_memory=False)
clin["bar15"] = clin["sampleID"].astype(str).str[:15]
clin = clin.drop_duplicates("bar15").set_index("bar15")
common = [s for s in beta.columns if s in clin.index]
clin = clin.loc[common]
beta = beta[common]
meno = clin["menopause_status"].astype(str)
clin["meno_grp"] = np.where(meno.str.startswith("Pre"), "Pre",
                    np.where(meno.str.startswith("Post"), "Post", "Other"))
clin["is_tumor"] = clin["sample_type"].eq("Primary Tumor")
clin["is_normal"] = clin["sample_type"].eq("Solid Tissue Normal")
clin["patient"] = clin.index.str[:12]
print(clin.groupby(["meno_grp", "sample_type"]).size())

# molecular subtype (cBioPortal pan-can-atlas patient file)
sub = pd.read_csv(os.path.join(D, "ucec_cbioportal_clinical_patient.txt"),
                  sep="\t", skiprows=4, low_memory=False)
sub = sub[["PATIENT_ID", "SUBTYPE"]].dropna()
submap = dict(zip(sub["PATIENT_ID"], sub["SUBTYPE"]))
clin["subtype"] = clin["patient"].map(submap)
print("subtype coverage:", clin["subtype"].notna().sum(), "/", len(clin))
print(clin.loc[clin["is_tumor"], "subtype"].value_counts(dropna=False))

def mval(b):
    b = np.clip(b, 1e-4, 1 - 1e-4)
    return np.log2(b / (1 - b))

M = mval(beta)

# ============ PART A: pre vs post + interaction ============
tu = clin[clin["is_tumor"]]
pre_s = tu[tu["meno_grp"] == "Pre"].index
post_s = tu[tu["meno_grp"] == "Post"].index
norm_s = clin[clin["is_normal"]].index
print(f"tumours Pre {len(pre_s)}  Post {len(post_s)}  normals {len(norm_s)}")

recs = []
for pr in FOCUS:
    if pr not in beta.index:
        continue
    b_pre = beta.loc[pr, pre_s].dropna(); b_post = beta.loc[pr, post_s].dropna()
    b_norm = beta.loc[pr, norm_s].dropna()
    t, p = stats.ttest_ind(M.loc[pr, pre_s].dropna(), M.loc[pr, post_s].dropna(), equal_var=False)
    pos_pre = (b_pre > 0.3).mean(); pos_post = (b_post > 0.3).mean()
    # Fisher on positivity
    a, b_ = int((b_pre > 0.3).sum()), int((b_pre <= 0.3).sum())
    c, d_ = int((b_post > 0.3).sum()), int((b_post <= 0.3).sum())
    _, pf = stats.fisher_exact([[a, b_], [c, d_]])
    recs.append(dict(probe=pr, n_pre=len(b_pre), n_post=len(b_post),
                     mean_pre=b_pre.mean(), mean_post=b_post.mean(),
                     d_pre_vs_post=b_pre.mean() - b_post.mean(),
                     welch_p=p, pos_pre=pos_pre, pos_post=pos_post, fisher_p=pf,
                     dbeta_pre=b_pre.mean() - b_norm.mean(),
                     dbeta_post=b_post.mean() - b_norm.mean()))
A = pd.DataFrame(recs)
A["welch_fdr"] = stats.false_discovery_control(A["welch_p"])
A.to_csv(os.path.join(R, "v97_pre_vs_post_tcga.csv"), index=False)
print("\n[A] pre vs post tumours (focus probes):")
print(A[["probe", "mean_pre", "mean_post", "d_pre_vs_post", "welch_fdr",
         "pos_pre", "pos_post", "fisher_p"]].round(3).to_string())

# interaction: M ~ tissue * menopause on samples with Pre/Post known
import statsmodels.formula.api as smf
irec = []
use = clin[clin["meno_grp"].isin(["Pre", "Post"]) & (clin["is_tumor"] | clin["is_normal"])]
for pr in FOCUS:
    if pr not in M.index:
        continue
    df = pd.DataFrame({"M": M.loc[pr, use.index],
                       "tissue": np.where(use["is_tumor"], "T", "N"),
                       "meno": use["meno_grp"]})
    df = df.dropna()
    try:
        fit = smf.ols("M ~ tissue * meno", data=df).fit()
        key = "tissue[T.T]:meno[T.Pre]"
        irec.append(dict(probe=pr, int_coef=fit.params.get(key, np.nan),
                         int_p=fit.pvalues.get(key, np.nan),
                         n_pre_norm=int(((df.meno == "Pre") & (df.tissue == "N")).sum())))
    except Exception as e:
        irec.append(dict(probe=pr, int_coef=np.nan, int_p=np.nan, n_pre_norm=np.nan))
I = pd.DataFrame(irec)
I["int_fdr"] = stats.false_discovery_control(I["int_p"].fillna(1))
I.to_csv(os.path.join(R, "v97_interaction_tcga.csv"), index=False)
print("\n[A] tissue x menopause interaction:")
print(I.round(4).to_string())

# ============ PART B: paired / endometrioid / subtype ============
tum_pat = set(clin[clin["is_tumor"]]["patient"])
nor_pat = set(clin[clin["is_normal"]]["patient"])
pair_pats = sorted(tum_pat & nor_pat)
pair_t = clin[(clin["is_tumor"]) & clin["patient"].isin(pair_pats)].groupby("patient").apply(lambda g: g.index[0], include_groups=False)
pair_n = clin[(clin["is_normal"]) & clin["patient"].isin(pair_pats)].groupby("patient").apply(lambda g: g.index[0], include_groups=False)
print(f"\n[B] matched pairs: {len(pair_pats)}")

prec = []
for pr in FOCUS:
    if pr not in beta.index:
        continue
    btv = beta.loc[pr, pair_t.values].to_numpy(dtype=float)
    bnv = beta.loc[pr, pair_n.values].to_numpy(dtype=float)
    ok = ~np.isnan(btv) & ~np.isnan(bnv)
    d = btv[ok] - bnv[ok]
    if ok.sum() >= 8 and (d != 0).any():
        w, pw = stats.wilcoxon(d)
    else:
        pw = np.nan
    prec.append(dict(probe=pr, n_pairs=int(ok.sum()),
                     mean_paired_dbeta=float(np.mean(d)) if len(d) else np.nan,
                     median_paired_dbeta=float(np.median(d)) if len(d) else np.nan,
                     wilcoxon_p=pw,
                     frac_pairs_tumour_higher=float((d > 0).mean()) if len(d) else np.nan))
B1 = pd.DataFrame(prec)
fdr_ok = B1["wilcoxon_p"].notna()
B1.loc[fdr_ok, "wilcoxon_fdr"] = stats.false_discovery_control(B1.loc[fdr_ok, "wilcoxon_p"])
B1.to_csv(os.path.join(R, "v97_paired_tcga.csv"), index=False)
print(B1.round(4).to_string())

# endometrioid subset
eec = tu[tu["histological_type"].astype(str).str.contains("Endometrioid", case=False)].index
erec = []
for pr in FOCUS:
    if pr not in beta.index:
        continue
    b_e = beta.loc[pr, eec].dropna(); b_n = beta.loc[pr, norm_s].dropna()
    t, p = stats.ttest_ind(M.loc[pr, eec].dropna(), M.loc[pr, norm_s].dropna(), equal_var=False)
    erec.append(dict(probe=pr, n_eec=len(b_e), dbeta_eec=b_e.mean() - b_n.mean(),
                     welch_p=p, pos03_eec=(b_e > 0.3).mean()))
B2 = pd.DataFrame(erec)
B2["welch_fdr"] = stats.false_discovery_control(B2["welch_p"])
B2.to_csv(os.path.join(R, "v97_endometrioid_tcga.csv"), index=False)
print(f"\n[B] endometrioid-only (n={len(eec)}):")
print(B2[["probe", "n_eec", "dbeta_eec", "welch_fdr", "pos03_eec"]].round(3).to_string())

# molecular subtype sensitivity (tumours with subtype, vs all normals)
srec = []
for st in ["UCEC_POLE", "UCEC_MSI", "UCEC_CN_LOW", "UCEC_CN_HIGH"]:
    ss = tu[tu["subtype"] == st].index
    if len(ss) < 5:
        continue
    for pr in FOCUS:
        if pr not in beta.index:
            continue
        b_s = beta.loc[pr, ss].dropna(); b_n = beta.loc[pr, norm_s].dropna()
        srec.append(dict(subtype=st, n=len(b_s), probe=pr,
                         dbeta=b_s.mean() - b_n.mean(), pos03=(b_s > 0.3).mean()))
B3 = pd.DataFrame(srec)
B3.to_csv(os.path.join(R, "v97_subtype_tcga.csv"), index=False)
piv = B3.pivot_table(index="probe", columns="subtype", values="dbeta")
print("\n[B] subtype-stratified dBeta:")
print(piv.round(3).to_string())

# ============ PART C: threshold retention scan ============
# recompute GSE73949 (n=17) P95 from raw for elite 476 probes; reuse elite table for other cohorts
elite = pd.read_csv(os.path.join(R, "elite_ranked_v6.csv"))
g73949_dir = os.path.join(D, "GSE73949", "ewas_betas")
elite_set = set(elite["probe"])
vals = {pr: [] for pr in elite_set}
for fn in os.listdir(g73949_dir):
    with open(os.path.join(g73949_dir, fn)) as f:
        for line in f:
            pr, v = line.rstrip("\n").split("\t")
            if pr in elite_set and v not in ("NA", "", "NaN"):
                vals[pr].append(float(v))
g7 = {pr: (np.percentile(v, 95) if v else np.nan) for pr, v in vals.items()}
elite["g73949_p95_raw"] = elite["probe"].map(g7)
print("\n[C] GSE73949 recomputed P95 for leads:",
      {p: round(g7.get(p, np.nan), 4) for p in LEADS if p in g7})

# per-dimension P95 table for scan
dims = pd.DataFrame({
    "probe": elite["probe"],
    "g73949": elite["g73949_p95_raw"],
    "ne347": elite["ne347_p95"],
    "cervix": elite["cervix_p95"],
    "blood": elite["blood_max_p95"],
})
dims["worst"] = dims[["g73949", "ne347", "cervix", "blood"]].max(axis=1)
dims["ne_pass"] = elite["ne_pass"].values
scan = []
thr_grid = np.round(np.arange(0.05, 0.301, 0.01), 2)
for thr in thr_grid:
    passed = dims[(dims["ne_pass"]) & (dims["worst"] <= thr)]
    scan.append(dict(threshold=thr, n_pass=len(passed),
                     boll_cg24589459=bool((passed["probe"] == "cg24589459").any()),
                     boll_cg07495363=bool((passed["probe"] == "cg07495363").any()),
                     zscan12_cg27577527=bool((passed["probe"] == "cg27577527").any()),
                     arl5c=bool((passed["probe"] == "cg15790037").any()),
                     ghsr=bool((passed["probe"] == "cg10109500").any())))
S = pd.DataFrame(scan)
S.to_csv(os.path.join(R, "v97_threshold_scan.csv"), index=False)
print(S.to_string())

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(S["threshold"], S["n_pass"], marker="o", color="#2c7fb8", lw=2, label="All elite probes passing")
for pr, col, lab in [("boll_cg24589459", "#d7301f", "BOLL cg24589459"),
                     ("boll_cg07495363", "#fc9272", "BOLL cg07495363"),
                     ("zscan12_cg27577527", "#238b45", "ZSCAN12 cg27577527")]:
    ax.step(S["threshold"], S[pr].astype(int), where="post", color=col, lw=1.8, label=lab)
ax.axvline(0.15, color="grey", ls="--", lw=1); ax.text(0.152, ax.get_ylim()[1]*0.9, "primary gate 0.15", fontsize=8, color="grey")
ax.set_xlabel("Background P95 gate (beta)"); ax.set_ylabel("Candidates retained (count / indicator)")
ax.set_title("Candidate retention across background-stringency thresholds")
ax.legend(fontsize=8); fig.tight_layout()
fig.savefig(os.path.join(R, "v97_threshold_retention_curve.png"), dpi=200)
print("figure saved")

# ============ PART D: bootstrap P95 CIs + Wilson CIs ============
# extract per-sample values for LEADS from raw cohort files (awk one-pass for GSE223817)
def extract_leads(cohort_dir, out_tsv):
    probe_re = "|".join(LEADS)
    cmd = f'awk -F"\\t" -v OFS="\\t" \'BEGIN{{split("{probe_re}",a,"|");for(i in a)keep[a[i]]=1}} ($1 in keep){{print FILENAME,$1,$2}}\' "{cohort_dir}"/*.txt'
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    with open(out_tsv, "w") as f:
        f.write(r.stdout)
    return out_tsv

cohorts = {
    "GSE73949_healthy17": dict(dir=os.path.join(D, "GSE73949", "ewas_betas"),
                               meta=None, filt=None),
    "GSE223817_benign347": dict(dir=os.path.join(D, "GSE223817", "ewas_betas"),
                                meta=os.path.join(R, "GSE223817_ewas_meta.csv"),
                                filt=("sample type", "control")),
    "GSE46306_cervix_HPVneg": dict(dir=os.path.join(D, "GSE46306", "ewas_betas"),
                                   meta=os.path.join(R, "GSE46306_ewas_meta.csv"),
                                   filt=("infection", "HPV-")),
    "GSE35069_blood60": dict(dir=os.path.join(D, "GSE35069", "ewas_betas"),
                             meta=None, filt=None),
}
rng = np.random.default_rng(42)
ci_rows = []
for cname, cfg in cohorts.items():
    tsv = extract_leads(cfg["dir"], os.path.join(R, f"_v97_{cname}_leads.tsv"))
    df = pd.read_csv(tsv, sep="\t", names=["file", "probe", "beta"])
    df["gsm"] = df["file"].str.replace(".txt", "", regex=False).str.replace(".*[/\\\\]", "", regex=True)
    df = df[df["beta"] != "NA"]
    df["beta"] = df["beta"].astype(float)
    if cfg["meta"]:
        meta = pd.read_csv(cfg["meta"])
        keep = set(meta.loc[meta[cfg["filt"][0]] == cfg["filt"][1], "sample id"])
        df = df[df["gsm"].isin(keep)]
    if cname == "GSE46306_cervix_HPVneg":
        pass
    for pr in LEADS:
        v = df.loc[df["probe"] == pr, "beta"].dropna().values
        if len(v) < 5:
            continue
        p95 = np.percentile(v, 95)
        boot = np.percentile(rng.choice(v, size=(10000, len(v)), replace=True), 95, axis=1)
        lo, hi = np.percentile(boot, [2.5, 97.5])
        ci_rows.append(dict(cohort=cname, probe=pr, n=len(v), p95=round(p95, 4),
                            boot_lo=round(lo, 4), boot_hi=round(hi, 4)))
        if cname == "GSE35069_blood60":
            pass
CI = pd.DataFrame(ci_rows)
CI.to_csv(os.path.join(R, "v97_background_bootstrap_ci.csv"), index=False)
print("\n[D] bootstrap P95 CIs:")
print(CI.to_string())

# per-cell-type worst blood P95 for leads (GSE35069)
meta35 = pd.read_csv(os.path.join(R, "GSE35069_ewas_meta.csv"))
df35 = pd.read_csv(os.path.join(R, "_v97_GSE35069_blood60_leads.tsv"), sep="\t",
                   names=["file", "probe", "beta"])
df35["gsm"] = df35["file"].str.replace(".txt", "", regex=False).str.replace(".*[/\\\\]", "", regex=True)
df35 = df35[df35["beta"] != "NA"]; df35["beta"] = df35["beta"].astype(float)
df35 = df35.merge(meta35[["sample id", "tissue"]], left_on="gsm", right_on="sample id")
ct = df35.groupby(["probe", "tissue"])["beta"].agg(n="count", p95=lambda x: np.percentile(x, 95)).reset_index()
ct.to_csv(os.path.join(R, "v97_blood_celltype_p95.csv"), index=False)
print(ct[ct["probe"].isin(LEADS[:3])].to_string())

# Wilson CIs for key positivity counts
def wilson(k, n):
    z = 1.96; p = k / n
    den = 1 + z**2 / n
    ctr = (p + z**2 / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / den
    return round(ctr - half, 4), round(ctr + half, 4)

wrows = []
# TCGA pre positivity from data
for pr in ["cg07495363", "cg24589459", "cg27577527", "cg23180938", "cg10109500", "cg15790037"]:
    if pr in beta.index:
        b = beta.loc[pr, pre_s].dropna(); k = int((b > 0.3).sum()); n = len(b)
        lo, hi = wilson(k, n)
        wrows.append(dict(group="TCGA premenopausal tumours", probe=pr, k=k, n=n,
                          rate=round(k / n, 4), wilson_lo=lo, wilson_hi=hi))
# paired tumours
for pr in ["cg07495363", "cg24589459"]:
    bt = beta.loc[pr, pair_t.values].dropna(); k = int((bt > 0.3).sum()); n = len(bt)
    lo, hi = wilson(k, n)
    wrows.append(dict(group="TCGA paired tumours", probe=pr, k=k, n=n,
                      rate=round(k / n, 4), wilson_lo=lo, wilson_hi=hi))
# external key counts (from v9.6 manuscript)
for grp, pr, k, n in [("GSE155760 benign cervix controls", "cg07495363", 0, 13),
                      ("GSE155760 benign cervix controls", "cg24589459", 0, 13),
                      ("GSE136791 carcinoma", "cg07495363", None, 69),
                      ("GSE136791 carcinoma", "cg24589459", None, 69)]:
    if k is None:
        ext = pd.read_csv(os.path.join(R, "external_validation_corrected.csv"))
        rate = float(ext.loc[ext["probe"] == pr, "GSE136791_ca_only"].iloc[0])
        k = int(round(rate * n))
    lo, hi = wilson(k, n)
    wrows.append(dict(group=grp, probe=pr, k=k, n=n, rate=round(k / n, 4),
                      wilson_lo=lo, wilson_hi=hi))
W = pd.DataFrame(wrows)
W.to_csv(os.path.join(R, "v97_positivity_wilson_ci.csv"), index=False)
print("\n[D] Wilson CIs:")
print(W.to_string())

# ============ PART E: export chr2/chr6 M matrix for bumphunter ============
sub_beta = beta.loc[beta.index.isin(set(chr26["probe"]))]
annot = chr26.set_index("probe").loc[sub_beta.index.intersection(chr26["probe"])]
outM = mval(sub_beta)
outM.insert(0, "chr", annot["chr"]); outM.insert(1, "pos", annot["pos"].astype(int))
outM.to_csv(os.path.join(R, "v97_bumphunter_chr26_M.tsv"), sep="\t")
pheno_bh = pd.DataFrame({"sample": outM.columns[2:],
                         "group": np.where(clin.loc[outM.columns[2:], "is_tumor"], "T", "N")})
pheno_bh.to_csv(os.path.join(R, "v97_bumphunter_pheno.tsv"), sep="\t", index=False)
print("\n[E] bumphunter input:", outM.shape, "pheno groups:", pheno_bh.group.value_counts().to_dict())
print("ALL DONE")
