"""v9.11 Analysis A prep: extract TCGA beta for 26 probes + phenotype table.
Output: results/v911_tcga_adj_input.tsv (samples x [probe betas + group/age/histology/subtype/menopause])
"""
import pandas as pd

PANEL23 = pd.read_csv('results/wetlab_panel_23_primer_annotation.csv')['probe'].tolist()
SUBS = ['cg23164203', 'cg20275132', 'cg25666433']
PROBES = PANEL23 + SUBS
print('probes:', len(PROBES))

# --- beta ---
beta = pd.read_csv('data/TCGA_UCEC_methylation450_beta.tsv', sep='\t', index_col=0)
beta = beta.loc[beta.index.isin(PROBES)]
print('beta probes found:', beta.shape[0], 'samples:', beta.shape[1])
missing = [p for p in PROBES if p not in beta.index]
print('missing probes:', missing)
beta = beta.T  # samples x probes
beta.index.name = 'sampleID'
beta = beta.reset_index()

# sample type from barcode position 14-15 (TCGA-XX-YYYY-01)
beta['code'] = beta['sampleID'].str.slice(13, 15)
beta['group'] = beta['code'].map(lambda c: 'tumour' if c == '01' else ('normal' if c == '11' else 'other'))
print(beta['group'].value_counts().to_dict())

# --- clinical ---
clin = pd.read_csv('data/UCEC_clinicalMatrix.tsv', sep='\t', low_memory=False)
clin = clin[['sampleID', 'age_at_initial_pathologic_diagnosis', 'histological_type', 'menopause_status']]
clin['age'] = pd.to_numeric(clin['age_at_initial_pathologic_diagnosis'], errors='coerce')

# --- subtype (patient level) ---
sub = pd.read_csv('data/ucec_cbioportal_clinical_patient.txt', sep='\t', skiprows=4)
sub = sub[['PATIENT_ID', 'SUBTYPE']]
beta['patient'] = beta['sampleID'].str.slice(0, 12)

df = beta.merge(clin, on='sampleID', how='left').merge(sub, left_on='patient', right_on='PATIENT_ID', how='left')
df['histology_simple'] = df['histological_type'].map(
    lambda h: 'endometrioid' if isinstance(h, str) and 'endometrioid' in h.lower() and 'mixed' not in h.lower()
    else ('other' if isinstance(h, str) else None))
df['subtype_simple'] = df['SUBTYPE'].map(lambda s: s.replace('UCEC_', '') if isinstance(s, str) else None)
df['menopause_simple'] = df['menopause_status'].map(
    lambda m: 'pre' if isinstance(m, str) and m.lower().startswith('pre')
    else ('peri' if isinstance(m, str) and m.lower().startswith('peri')
    else ('post' if isinstance(m, str) and m.lower().startswith('post') else None)))

keep = ['sampleID', 'group', 'age', 'histology_simple', 'subtype_simple', 'menopause_simple'] + [p for p in PROBES if p in df.columns]
out = df[keep]
out = out[out['group'].isin(['tumour', 'normal'])]
out.to_csv('results/v911_tcga_adj_input.tsv', sep='\t', index=False)
print(out.shape)
print(out[['group']].value_counts().to_dict())
print('tumours with age:', out.loc[out.group=='tumour','age'].notna().sum(),
      '| histology:', out.loc[out.group=='tumour','histology_simple'].value_counts().to_dict(),
      '| subtype:', out.loc[out.group=='tumour','subtype_simple'].value_counts(dropna=False).to_dict(),
      '| menopause:', out.loc[out.group=='tumour','menopause_simple'].value_counts(dropna=False).to_dict())
print('normals with age:', out.loc[out.group=='normal','age'].notna().sum())
