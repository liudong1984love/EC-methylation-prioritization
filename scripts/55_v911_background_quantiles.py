"""v9.11 Analysis B: background distribution summaries per cohort x probe.
Reports n, median, P90, P95, max, and proportion > 0.15 (R quantile type 7).
Output: results/v911_background_quantiles.csv
"""
import glob, os
import numpy as np
import pandas as pd

PANEL23 = pd.read_csv('results/wetlab_panel_23_primer_annotation.csv')['probe'].tolist()
SUBS = ['cg23164203', 'cg20275132', 'cg25666433']
PROBES = set(PANEL23 + SUBS)

def read_sample_probes(path):
    vals = {}
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            pid, b = parts[0], parts[1]
            if pid in PROBES and b not in ('NA', '', 'NaN'):
                vals[pid] = float(b)
    return vals

def cohort_matrix(gsm_list, cohort_dir):
    rows = {}
    for gsm in gsm_list:
        path = os.path.join(cohort_dir, f'{gsm}.txt')
        if os.path.exists(path):
            rows[gsm] = read_sample_probes(path)
    return pd.DataFrame.from_dict(rows, orient='index')

def summarise(mat, cohort_name):
    recs = []
    for p in PANEL23 + SUBS:
        if p not in mat.columns:
            continue
        x = mat[p].dropna().values
        if len(x) == 0:
            continue
        recs.append(dict(cohort=cohort_name, probe=p, n=int(len(x)),
                         median=round(float(np.median(x)), 4),
                         p90=round(float(np.quantile(x, 0.9, method='linear')), 4),
                         p95=round(float(np.quantile(x, 0.95, method='linear')), 4),
                         max=round(float(np.max(x)), 4),
                         prop_gt_0_15=round(float(np.mean(x > 0.15)), 4)))
    return recs

out = []

# GSE73949: all 17 healthy
gsms = [os.path.basename(f)[:-4] for f in glob.glob('data/GSE73949/ewas_betas/*.txt')]
m = cohort_matrix(gsms, 'data/GSE73949/ewas_betas')
out += summarise(m, 'GSE73949_healthy17')
print('GSE73949', m.shape)

# GSE223817: 347 controls
meta = pd.read_csv('results/GSE223817_ewas_meta.csv')
ctrl = meta.loc[meta['sample type'] == 'control', 'sample id'].tolist()
print('GSE223817 controls:', len(ctrl))
m = cohort_matrix(ctrl, 'data/GSE223817/ewas_betas')
out += summarise(m, 'GSE223817_benign347')
print('GSE223817', m.shape)

# GSE46306: HPV-negative normal cervix
meta = pd.read_csv('results/GSE46306_ewas_meta.csv')
hpvneg = meta.loc[(meta['infection'] == 'HPV-') & (meta['sample type'] == 'control'), 'sample id'].tolist()
print('GSE46306 HPV-neg controls:', len(hpvneg))
m = cohort_matrix(hpvneg, 'data/GSE46306/ewas_betas')
out += summarise(m, 'GSE46306_cervix_HPVneg20')
print('GSE46306', m.shape)

# GSE35069: purified blood cell types (exclude whole blood/PBMC/granulocyte composites)
meta = pd.read_csv('results/GSE35069_ewas_meta.csv')
pur = meta.loc[~meta['tissue'].isin(['whole blood', 'peripheral blood mononuclear cell', 'granulocyte']), 'sample id'].tolist()
print('GSE35069 purified:', len(pur))
m = cohort_matrix(pur, 'data/GSE35069/ewas_betas')
out += summarise(m, 'GSE35069_blood_purified42')
print('GSE35069', m.shape)

df = pd.DataFrame(out)
df.to_csv('results/v911_background_quantiles.csv', index=False)
print(df.groupby('cohort')['probe'].count())
# key check: which probes exceed 0.15 at P95 / max per cohort
print(df[df.p95 > 0.15][['cohort','probe','p95','max','prop_gt_0.15']].to_string())
