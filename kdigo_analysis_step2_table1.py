#!/usr/bin/env python3
"""
KDIGO AKI 分期 vs 预后 — Step 2: Table 1 + 合并症 + 生命体征
==============================================================
产出:
  - table1_baseline.csv: 按 KDIGO 分期的基线特征表
  - table1_formatted.txt: 可直接用于论文的格式化 Table 1
"""
import pandas as pd
import numpy as np
import os
import gc
from datetime import datetime

print("=" * 70)
print("KDIGO-AKI Step 2: Table 1 — Baseline Characteristics")
print(f"Start: {datetime.now()}")
print("=" * 70)

BASE = "E:/mimic-iv"
HOSP = os.path.join(BASE, "v3.1/physionet.org/files/mimiciv/3.1/hosp")
ICU = os.path.join(BASE, "v3.1/physionet.org/files/mimiciv/3.1/icu")

def read_gz(path, **kwargs):
    return pd.read_csv(path, compression='gzip', low_memory=False, **kwargs)

# ============================================================
# 1. Load KDIGO cohort
# ============================================================
print("[1] Loading KDIGO cohort...")

IN = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
cohort = pd.read_csv(os.path.join(IN, "kdigo_cohort.csv"), low_memory=False)
print(f"  Cohort: {len(cohort):,} rows")

# Stage as categorical
cohort['stage_label'] = cohort['kdigo_stage'].map({0: 'No AKI', 1: 'Stage 1', 2: 'Stage 2', 3: 'Stage 3'})

# ============================================================
# 2. Comorbidities (Elixhauser / Charlson via ICD)
# ============================================================
print("[2] Building comorbidities...")

diagnoses = read_gz(os.path.join(HOSP, "diagnoses_icd.csv.gz"))
d_icd = read_gz(os.path.join(HOSP, "d_icd_diagnoses.csv.gz"))

# Merge comorbidity status into cohort
cohort_sids = set(cohort['subject_id'].unique())
diag_cohort = diagnoses[diagnoses['subject_id'].isin(cohort_sids)]

# Elixhauser-like categories (simplified for MIMIC ICD-9/10 mix)
comorb_map = {
    'Congestive Heart Failure': r'^(I50|I11\.0|I13\.[02]|I25\.5|I42|428|398\.91|402\.[01xi]|404\.[01x]3)',
    'Cardiac Arrhythmias': r'^(I4[4-9]|427|426)',
    'Valvular Disease': r'^(I0[5-8]|I3[4-9]|39[3-4]|424)',
    'Hypertension': r'^(I1[0-5]|40[1-5])',
    'Diabetes (uncomplicated)': r'^(E1[0-4][0-9]|250[0-3])',
    'Diabetes (complicated)': r'^(E1[0-4][2-9]|250[4-9])',
    'Chronic Pulmonary': r'^(J4[0-7]|J6[0-7]|49[0-6]|500|501|502|503|504|505|506\.4)',
    'Liver Disease': r'^(K7[0-7]|B18|571|572|456\.[0-2])',
    'Coagulopathy': r'^(D6[5-9]|286|287\.[3-5])',
    'Fluid/Electrolyte': r'^(E86|E87|276)',
    'Metastatic Cancer': r'^(C7[7-9]|C80|19[6-9])',
    'Solid Tumor (non-met)': r'^(C[0-7][0-6]|1[4-9]\d|20[0-8])',
    'Obesity': r'^(E66|278\.0)',
    'Peripheral Vascular': r'^(I7[0-9]|44[0-8]|09[3-4])',
    'Cerebrovascular': r'^(I6[0-9]|43[0-8]|G45|433|434|435|436|437|438)',
    'Sepsis': r'^(A4[01]|R65\.[2]|038|995\.9[12]|785\.52)',
    'Myocardial Infarction': r'^(I2[1-2]|I25\.2|410|412)',
    'Peptic Ulcer': r'^(K2[5-8]|53[1-4])',
}

comorb_status = {}
for comorb, pattern in comorb_map.items():
    matched = diag_cohort[diag_cohort['icd_code'].str.match(pattern, na=False)]
    has_comorb = set(matched['subject_id'].unique())
    comorb_status[comorb] = has_comorb
    print(f"  {comorb}: {len(has_comorb):,} patients")

# Apply to cohort
for comorb, patients_set in comorb_status.items():
    cohort[f'com_{comorb.replace(" ", "_").replace("(", "").replace(")", "")}'] = \
        cohort['subject_id'].isin(patients_set).astype(int)

print(f"  Comorbidities added: {len(comorb_status)} conditions")

# ============================================================
# 3. ICU Admission Vital Signs
# ============================================================
print("[3] Extracting ICU admission vital signs...")

chartevents_path = os.path.join(ICU, "chartevents.csv.gz")

# Vitals itemids
vital_itemids = {
    'heart_rate': [220045],
    'sbp': [220050],          # non-invasive SBP
    'dbp': [220051],          # non-invasive DBP
    'mbp': [220052],          # non-invasive MBP
    'resp_rate': [220210],
    'temperature': [223761],  # Temperature Fahrenheit
    'spo2': [220277],
}

# For each ICU stay, get first 24h vitals
# Strategy: read chartevents and filter by our cohort's stay_ids
cohort_stay_ids = set(cohort['stay_id'].unique()) if 'stay_id' in cohort.columns else set()
cohort_hadm_ids = set(cohort['hadm_id'].unique())

vital_data = []
chunk_size = 2000000
total_processed = 0

for i, chunk in enumerate(pd.read_csv(chartevents_path, compression='gzip', chunksize=chunk_size, low_memory=False)):
    # Filter: only relevant itemids
    all_vital_ids = [item for sublist in vital_itemids.values() for item in sublist]
    chunk = chunk[chunk['itemid'].isin(all_vital_ids)].copy()
    
    # Filter: only cohort patients
    chunk = chunk[chunk['hadm_id'].isin(cohort_hadm_ids)].copy()
    
    if len(chunk) > 0:
        chunk['valuenum'] = pd.to_numeric(chunk['valuenum'], errors='coerce')
        chunk = chunk.dropna(subset=['valuenum'])
        vital_data.append(chunk)
    
    total_processed += chunk_size if hasattr(chunk, 'shape') else 0
    if (i + 1) % 5 == 0:
        print(f"  ... processed ~{total_processed/1e6:.0f}M rows, kept {sum(len(v) for v in vital_data):,} vital records")

vitals_all = pd.concat(vital_data, ignore_index=True)
del vital_data; gc.collect()
print(f"  Total vital records: {len(vitals_all):,}")

# Merge with ICU admission time
vitals_all['charttime'] = pd.to_datetime(vitals_all['charttime'])

# Get first 24h vitals per ICU stay
# Merge with intime for time filtering
vitals_with_time = vitals_all.merge(cohort[['hadm_id', 'stay_id', 'intime']].rename(columns={'stay_id': 'cohort_stay_id'}), on='hadm_id', how='inner')
vitals_with_time['intime'] = pd.to_datetime(vitals_with_time['intime'])
vitals_with_time['hours_from_icu'] = (vitals_with_time['charttime'] - vitals_with_time['intime']).dt.total_seconds() / 3600

# First 24h only
vitals_24h = vitals_with_time[(vitals_with_time['hours_from_icu'] >= 0) & (vitals_with_time['hours_from_icu'] <= 24)]

# Summary per stay_id
vital_summary = {}
for name, itemids in vital_itemids.items():
    sub = vitals_24h[vitals_24h['itemid'].isin(itemids)]
    if name == 'temperature':
        # Convert F to C
        sub['valuenum'] = (sub['valuenum'] - 32) * 5 / 9
    agg = sub.groupby('stay_id')['valuenum'].agg(['mean', 'min', 'max', 'first']).reset_index()
    agg.columns = ['stay_id', f'{name}_mean', f'{name}_min', f'{name}_max', f'{name}_first']
    vital_summary[name] = agg
    print(f"  {name}: {len(sub):,} records, {len(agg):,} stays")

# Merge vitals into cohort
for name, df in vital_summary.items():
    cohort = cohort.merge(df, on='stay_id', how='left')

print(f"  Vital signs merged. Missing rates:")
for name in vital_itemids:
    col = f'{name}_mean'
    if col in cohort.columns:
        miss = cohort[col].isna().mean() * 100
        print(f"    {name}: {miss:.1f}% missing")

# ============================================================
# 4. APACHE severity scores (if available)
# ============================================================
print("[4] Extracting severity scores...")

# Try loading Apache scores directly
try:
    apsiii = read_gz(os.path.join(ICU, "apsiii.csv.gz"))
    print(f"  apsiii: {len(apsiii):,} rows")
    # Keep first row per stay
    apsiii = apsiii.sort_values('charttime').drop_duplicates('stay_id', keep='first')
    cohort = cohort.merge(apsiii[['stay_id', 'apsiii', 'apsiii_prob', 'apsiii_hr']],
                          on='stay_id', how='left')
    print(f"  APS III merged: {cohort['apsiii'].notna().sum():,} stays")
except Exception as e:
    print(f"  APS III not available: {e}")

# ============================================================
# 5. Lab values on ICU admission (beyond Cr)
# ============================================================
print("[5] Extracting key admission labs...")

labevents_path = os.path.join(HOSP, "labevents.csv.gz")

key_labs = {
    'bun': 51006,           # Urea Nitrogen (BUN)
    'egfr': 50920,          # eGFR (MDRD)
    'sodium': 50983,        # Sodium
    'potassium': 50971,     # Potassium
    'chloride': 50902,      # Chloride
    'bicarbonate': 50882,   # Bicarbonate
    'lactate': 50813,       # Lactate
    'wbc': 51301,           # WBC
    'hemoglobin': 51222,    # Hemoglobin
    'platelet': 51265,      # Platelet
    'bilirubin_total': 50885, # Total Bilirubin
    'albumin': 50862,       # Albumin
    'glucose': 50931,       # Glucose
    'ph': 50820,           # pH (arterial)
}

all_lab_ids = list(key_labs.values())
lab_data = []
chunk_size = 2000000

for i, chunk in enumerate(pd.read_csv(labevents_path, compression='gzip', chunksize=chunk_size, low_memory=False)):
    chunk = chunk[chunk['itemid'].isin(all_lab_ids)].copy()
    chunk = chunk[chunk['hadm_id'].isin(cohort_hadm_ids)].copy()
    if len(chunk) > 0:
        chunk['valuenum'] = pd.to_numeric(chunk['valuenum'], errors='coerce')
        chunk = chunk.dropna(subset=['valuenum'])
        lab_data.append(chunk[['subject_id', 'hadm_id', 'itemid', 'valuenum', 'charttime']])
    if (i + 1) % 10 == 0:
        print(f"  ... processed ~{chunk_size*(i+1)/1e6:.0f}M rows")

labs_all = pd.concat(lab_data, ignore_index=True)
del lab_data; gc.collect()
labs_all['charttime'] = pd.to_datetime(labs_all['charttime'])

# Get first measurement within +/- 24h of ICU admission per stay
labs_with_time = labs_all.merge(cohort[['hadm_id', 'stay_id', 'intime']], on='hadm_id', how='inner')
labs_with_time['intime'] = pd.to_datetime(labs_with_time['intime'])
labs_with_time['hours'] = (labs_with_time['charttime'] - labs_with_time['intime']).dt.total_seconds() / 3600
labs_admission = labs_with_time[(labs_with_time['hours'] >= -24) & (labs_with_time['hours'] <= 24)]
labs_admission = labs_admission.sort_values('hours')

for name, itemid in key_labs.items():
    sub = labs_admission[labs_admission['itemid'] == itemid]
    first_val = sub.groupby('stay_id').first().reset_index()
    first_val = first_val[['stay_id', 'valuenum']].rename(columns={'valuenum': f'lab_{name}'})
    cohort = cohort.merge(first_val, on='stay_id', how='left')
    n_avail = cohort[f'lab_{name}'].notna().sum()
    print(f"  lab_{name}: {n_avail:,}/{len(cohort):,} available ({n_avail/len(cohort)*100:.1f}%)")

# ============================================================
# 6. Generate Table 1
# ============================================================
print("\n[5.5] Saving enriched cohort (before Table 1)...")

# Add derived columns for Table 1
cohort['male'] = cohort['gender'].map({'M': 1, 'F': 0}).fillna(0).astype(int) if 'gender' in cohort.columns else 0
cohort['race_white'] = cohort['race'].str.contains('White', na=False).astype(int) if 'race' in cohort.columns else 0

# Save enriched cohort now (in case Table 1 fails)
cohort.to_csv(os.path.join(IN, "kdigo_cohort_enriched.csv"), index=False)
print(f"  Enriched cohort saved: {len(cohort):,} rows, {len(cohort.columns)} columns")

# ============================================================
# 6. Generate Table 1
# ============================================================
print("\n[6] Generating Table 1...")

def describe_var(df, col, is_categorical=False, decimals=1):
    """Describe variable by KDIGO stage. col can be a string column name."""
    results = {}
    for stage in [0, 1, 2, 3]:
        sub = df[df['kdigo_stage'] == stage]
        if is_categorical:
            n = sub[col].sum()
            pct = n / len(sub) * 100 if len(sub) > 0 else 0
            results[f'Stage_{stage}'] = f"{n} ({pct:.{decimals}f}%)"
        else:
            vals = sub[col].dropna()
            if len(vals) == 0:
                results[f'Stage_{stage}'] = "N/A"
            else:
                results[f'Stage_{stage}'] = f"{vals.mean():.{decimals}f} ± {vals.std():.{decimals}f}"
    # Overall
    if is_categorical:
        n = df[col].sum()
        pct = n / len(df) * 100
        results['Overall'] = f"{n} ({pct:.{decimals}f}%)"
    else:
        vals = df[col].dropna()
        results['Overall'] = f"{vals.mean():.{decimals}f} ± {vals.std():.{decimals}f}"
    return results

def row_cat(name, col):
    s = describe_var(cohort, col, True)
    return (name, s['Stage_0'], s['Stage_1'], s['Stage_2'], s['Stage_3'], s['Overall'])

def row_cont(name, col):
    s = describe_var(cohort, col, False, decimals=1)
    return (name, s['Stage_0'], s['Stage_1'], s['Stage_2'], s['Stage_3'], s['Overall'])

table1_rows = []
table1_rows.append(('Demographics', '', '', '', '', ''))
table1_rows.append(('  N', len(cohort[cohort['kdigo_stage']==0]), len(cohort[cohort['kdigo_stage']==1]),
                    len(cohort[cohort['kdigo_stage']==2]), len(cohort[cohort['kdigo_stage']==3]), len(cohort)))
table1_rows.append(row_cont('  Age (years)', 'anchor_age'))
table1_rows.append(row_cat('  Male, n (%)', 'male'))
table1_rows.append(row_cat('  Race — White, n (%)', 'race_white'))

# Cr metrics
table1_rows.append(('Renal Function', '', '', '', '', ''))
table1_rows.append(row_cont('  Baseline Cr (mg/dL)', 'baseline_cr'))
table1_rows.append(row_cont('  Peak Cr (mg/dL)', 'peak_cr'))
table1_rows.append(row_cont('  Cr ratio (peak/baseline)', 'cr_ratio'))
table1_rows.append(row_cat('  CKD (ICD), n (%)', 'ckd_icd'))

# Comorbidities
table1_rows.append(('Comorbidities', '', '', '', '', ''))
comorb_cols = [c for c in cohort.columns if c.startswith('com_')]
for col in sorted(comorb_cols):
    name = col.replace('com_', '').replace('_', ' ')
    table1_rows.append(row_cat(f'  {name}', col))

# Vitals
table1_rows.append(('Vital Signs (24h mean)', '', '', '', '', ''))
for vname in ['heart_rate', 'sbp', 'dbp', 'mbp', 'resp_rate', 'temperature', 'spo2']:
    col = f'{vname}_mean'
    if col in cohort.columns:
        table1_rows.append(row_cont(f'  {vname}', col))

# Labs
table1_rows.append(('Admission Labs', '', '', '', '', ''))
for lname in ['bun', 'sodium', 'potassium', 'bicarbonate', 'lactate', 'wbc', 'hemoglobin', 'platelet', 'bilirubin_total', 'albumin', 'glucose']:
    col = f'lab_{lname}'
    if col in cohort.columns:
        table1_rows.append(row_cont(f'  {lname}', col))

# Outcomes
table1_rows.append(('Outcomes', '', '', '', '', ''))
table1_rows.append(row_cat('  Hospital mortality, n (%)', 'hospital_expire_flag'))
table1_rows.append(row_cont('  ICU LOS (hours)', 'icu_los_hours'))
table1_rows.append(row_cont('  Hospital LOS (days)', 'hosp_los_days'))

# Save Table 1
df_table1 = pd.DataFrame(table1_rows, columns=['Variable', 'No AKI', 'Stage 1', 'Stage 2', 'Stage 3', 'Overall'])
df_table1.to_csv(os.path.join(IN, 'table1_baseline.csv'), index=False)
print(f"  Table 1 saved: {len(df_table1)} rows")

print("\n" + "=" * 70)
print("Step 2 Complete!")
print("=" * 70)
