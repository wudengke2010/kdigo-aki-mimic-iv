#!/usr/bin/env python3
"""
KDIGO AKI 分期 vs 预后 — Step 1: 构建队列
=============================================
基于 KDIGO 肌酐标准对 ICU 患者进行 AKI 分期:
  Stage 1: Cr 1.5-1.9x baseline OR ≥0.3 mg/dL increase
  Stage 2: Cr 2.0-2.9x baseline
  Stage 3: Cr ≥3.0x baseline OR Cr ≥4.0 OR RRT initiation

产出: kdigo_cohort.csv — 每个 ICU stay 一行，含分期+结局
"""
import pandas as pd
import numpy as np
import os, sys, gc, json
from datetime import datetime

print("=" * 70)
print("KDIGO-AKI Cohort Builder")
print(f"Start: {datetime.now()}")
print("=" * 70)

BASE = "E:/mimic-iv"
HOSP = os.path.join(BASE, "v3.1/physionet.org/files/mimiciv/3.1/hosp")
ICU = os.path.join(BASE, "v3.1/physionet.org/files/mimiciv/3.1/icu")

OUT = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"

def read_gz(path, **kwargs):
    return pd.read_csv(path, compression='gzip', low_memory=False, **kwargs)

# ============================================================
# 1. Load base tables
# ============================================================
print("\n[1] Loading base tables...")

patients = read_gz(os.path.join(HOSP, "patients.csv.gz"))
admissions = read_gz(os.path.join(HOSP, "admissions.csv.gz"))
icustays = read_gz(os.path.join(ICU, "icustays.csv.gz"))
diagnoses = read_gz(os.path.join(HOSP, "diagnoses_icd.csv.gz"))
d_icd = read_gz(os.path.join(HOSP, "d_icd_diagnoses.csv.gz"))
d_lab = read_gz(os.path.join(HOSP, "d_labitems.csv.gz"))

print(f"  patients: {len(patients):,}")
print(f"  admissions: {len(admissions):,}")
print(f"  icustays: {len(icustays):,}")

# ============================================================
# 2. Find renal patients via ICD codes
# ============================================================
print("\n[2] Identifying renal patients via ICD...")

aki10 = diagnoses[(diagnoses['icd_version'] == 10) & diagnoses['icd_code'].str.match(r'^N17', na=False)]
aki9 = diagnoses[(diagnoses['icd_version'] == 9) & diagnoses['icd_code'].str.match(r'^584', na=False)]
ckd10 = diagnoses[(diagnoses['icd_version'] == 10) & diagnoses['icd_code'].str.match(r'^N18', na=False)]
ckd9 = diagnoses[(diagnoses['icd_version'] == 9) & diagnoses['icd_code'].str.match(r'^585', na=False)]

aki_ids = set(aki10['subject_id'].unique()) | set(aki9['subject_id'].unique())
ckd_ids = set(ckd10['subject_id'].unique()) | set(ckd9['subject_id'].unique())

print(f"  AKI patients (ICD): {len(aki_ids):,}")
print(f"  CKD patients (ICD): {len(ckd_ids):,}")

# ============================================================
# 3. ICU stays — merge with admission info
# ============================================================
print("\n[3] Preparing ICU cohort...")

# Only keep first ICU stay per hospital admission
icu = icustays.sort_values('intime').drop_duplicates(subset=['hadm_id'], keep='first').copy()
print(f"  Unique ICU stays (first per hadm): {len(icu):,}")

# Merge admission info
adm_cols = ['hadm_id', 'subject_id', 'admittime', 'dischtime', 'deathtime',
            'admission_type', 'admission_location', 'discharge_location',
            'insurance', 'language', 'marital_status', 'race',
            'hospital_expire_flag', 'edregtime', 'edouttime']
icu = icu.merge(admissions[adm_cols], on=['subject_id', 'hadm_id'], how='left')

# Merge patient info
pat_cols = ['subject_id', 'gender', 'anchor_age', 'anchor_year', 'dod']
icu = icu.merge(patients[pat_cols], on='subject_id', how='left')

# Time columns
for col in ['intime', 'outtime', 'admittime', 'dischtime', 'deathtime', 'dod']:
    if col in icu.columns:
        icu[col] = pd.to_datetime(icu[col])

# ICU LOS
icu['icu_los_hours'] = (icu['outtime'] - icu['intime']).dt.total_seconds() / 3600
icu['hosp_los_days'] = (icu['dischtime'] - icu['admittime']).dt.total_seconds() / 86400

# Group labels
def get_renal_group(sid):
    if sid in aki_ids and sid in ckd_ids:
        return 'AKI+CKD'
    elif sid in aki_ids:
        return 'AKI'
    elif sid in ckd_ids:
        return 'CKD'
    else:
        return 'Non-renal'

icu['renal_group_icd'] = icu['subject_id'].apply(get_renal_group)

# AKI diagnosis flag (ICD-based, for reference only)
icu['aki_icd'] = icu['subject_id'].isin(aki_ids).astype(int)
icu['ckd_icd'] = icu['subject_id'].isin(ckd_ids).astype(int)

print(f"  ICU cohort size: {len(icu):,}")
print(icu['renal_group_icd'].value_counts())

# ============================================================
# 4. Extract Creatinine from labevents
# ============================================================
print("\n[4] Extracting creatinine from labevents...")

creatinine_itemid = 50912  # 血清肌酐

# Strategy: chunked read, filter only Cr, only patients with ICU stays
icu_subject_ids = set(icu['subject_id'].unique())
print(f"  ICU patients to filter: {len(icu_subject_ids):,}")

labevents_path = os.path.join(HOSP, "labevents.csv.gz")

cr_chunks = []
chunk_size = 2000000
total_rows = 0
matched_rows = 0

# First pass: get all Cr measurements for ICU patients
for i, chunk in enumerate(pd.read_csv(labevents_path, compression='gzip', chunksize=chunk_size, low_memory=False)):
    total_rows += len(chunk)
    
    # Filter: only creatinine itemid
    chunk = chunk[chunk['itemid'] == creatinine_itemid].copy()
    
    # Filter: only ICU patients
    chunk = chunk[chunk['subject_id'].isin(icu_subject_ids)].copy()
    
    if len(chunk) > 0:
        # Keep essential columns
        chunk = chunk[['subject_id', 'hadm_id', 'charttime', 'valuenum']].copy()
        chunk['valuenum'] = pd.to_numeric(chunk['valuenum'], errors='coerce')
        chunk = chunk.dropna(subset=['valuenum'])
        chunk = chunk[chunk['valuenum'] > 0]
        chunk = chunk[chunk['valuenum'] < 100]  # remove absurd values
        cr_chunks.append(chunk)
        matched_rows += len(chunk)
    
    if (i + 1) % 10 == 0:
        print(f"  ... processed {total_rows/1e6:.1f}M lab rows, kept {matched_rows:,} Cr measurements")

print(f"  Total labevents rows scanned: {total_rows:,}")
print(f"  Total Cr measurements kept: {matched_rows:,}")

cr_all = pd.concat(cr_chunks, ignore_index=True)
del cr_chunks
gc.collect()

# Parse charttime
cr_all['charttime'] = pd.to_datetime(cr_all['charttime'])
cr_all = cr_all.sort_values(['subject_id', 'hadm_id', 'charttime'])

print(f"  Creatinine measurements: {len(cr_all):,}")
print(f"  Unique patients with Cr: {cr_all['subject_id'].nunique():,}")
print(f"  Unique hospital stays with Cr: {cr_all['hadm_id'].nunique():,}")

# ============================================================
# 5. Calculate baseline and peak Cr per ICU stay
# ============================================================
print("\n[5] Calculating baseline & peak Cr per ICU stay...")

# Merge Cr with ICU stay times
cr_icu = cr_all.merge(icu[['hadm_id', 'intime', 'outtime']], on='hadm_id', how='inner')

# Baseline Cr: minimum Cr within the first 7 days before ICU admission OR first measurement
# For MIMIC, we use the first hospital Cr measurement before/at ICU admission as baseline

# Filter Cr measurements before or at ICU intime (for baseline)
# And all Cr during ICU stay (for peak)
cr_icu['hours_from_icu_in'] = (cr_icu['charttime'] - cr_icu['intime']).dt.total_seconds() / 3600

# Baseline Cr: first measurement within the hospital stay (before/during ICU)
# Use hosp admission Cr (the first measurement per hospital stay)
baseline_cr = cr_all.groupby('hadm_id').first().reset_index()
baseline_cr = baseline_cr[['hadm_id', 'valuenum']].rename(columns={'valuenum': 'baseline_cr'})

# ICU peak Cr: highest during ICU stay (including first 24h)
peak_cr = cr_icu.groupby('hadm_id')['valuenum'].max().reset_index()
peak_cr.columns = ['hadm_id', 'peak_cr']

# First Cr of ICU stay
first_icu_cr = cr_icu[cr_icu['hours_from_icu_in'] >= -24].groupby('hadm_id').first().reset_index()
first_icu_cr = first_icu_cr[['hadm_id', 'valuenum']].rename(columns={'valuenum': 'icu_first_cr'})

# Last Cr of ICU stay
last_icu_cr = cr_icu.groupby('hadm_id').last().reset_index()
last_icu_cr = last_icu_cr[['hadm_id', 'valuenum']].rename(columns={'valuenum': 'icu_last_cr'})

# Count of Cr measurements per stay
cr_count = cr_icu.groupby('hadm_id').size().reset_index(name='n_cr_measurements')

# Merge all Cr metrics
icu = icu.merge(baseline_cr, on='hadm_id', how='left')
icu = icu.merge(peak_cr, on='hadm_id', how='left')
icu = icu.merge(first_icu_cr, on='hadm_id', how='left')
icu = icu.merge(last_icu_cr, on='hadm_id', how='left')
icu = icu.merge(cr_count, on='hadm_id', how='left')

print(f"  Patients with baseline Cr: {icu['baseline_cr'].notna().sum():,}")
print(f"  Patients with peak Cr: {icu['peak_cr'].notna().sum():,}")

# ============================================================
# 6. KDIGO Staging (Creatinine criteria)
# ============================================================
print("\n[6] Applying KDIGO creatinine staging...")

# Need at least baseline Cr to stage
valid = icu[icu['baseline_cr'].notna() & icu['peak_cr'].notna()].copy()
print(f"  Stays with valid Cr data: {len(valid):,}")

valid['cr_ratio'] = valid['peak_cr'] / valid['baseline_cr']
valid['cr_increase'] = valid['peak_cr'] - valid['baseline_cr']

# KDIGO stage
conditions = [
    # Stage 3: Cr >= 3x baseline OR Cr >= 4.0 OR RRT
    (valid['cr_ratio'] >= 3.0) | (valid['peak_cr'] >= 4.0),
    # Stage 2: 2.0-2.9x baseline
    (valid['cr_ratio'] >= 2.0) & (valid['cr_ratio'] < 3.0),
    # Stage 1: 1.5-1.9x baseline OR >=0.3 increase
    ((valid['cr_ratio'] >= 1.5) & (valid['cr_ratio'] < 2.0)) | (valid['cr_increase'] >= 0.3),
]
choices = [3, 2, 1]
valid['kdigo_stage'] = np.select(conditions, choices, default=0)

print(f"\n  KDIGO Stage distribution:")
for s in [0, 1, 2, 3]:
    n = (valid['kdigo_stage'] == s).sum()
    pct = n / len(valid) * 100
    death_n = (valid[valid['kdigo_stage'] == s]['hospital_expire_flag'] == 1).sum()
    death_pct = death_n / n * 100 if n > 0 else 0
    print(f"    Stage {s}: {n:,} ({pct:.1f}%), mortality={death_n} ({death_pct:.1f}%)")

# ============================================================
# 7. CRRT/Dialysis flag for Stage 3 upgrade
# ============================================================
print("\n[7] Adding RRT status for Stage 3...")

proc_events = read_gz(os.path.join(ICU, "procedureevents.csv.gz"))
d_items = read_gz(os.path.join(ICU, "d_items.csv.gz"))

# Dialysis itemids
dial_items = d_items[d_items['label'].str.lower().str.contains(
    'dialysis - crrt|dialysis - cvvhd|dialysis - cvvhdf|hemodialysis|dialysis - scuf',
    na=False)]
dial_itemids = set(dial_items['itemid'].tolist())
print(f"  RRT itemids: {dial_itemids}")

rrt_events = proc_events[proc_events['itemid'].isin(dial_itemids)]
rrt_hadms = set(rrt_events['hadm_id'].unique())
print(f"  RRT events: {len(rrt_events):,}")
print(f"  RRT patients (ICU): {rrt_events['subject_id'].nunique():,}")
print(f"  RRT hospital stays: {len(rrt_hadms):,}")

# Upgrade to Stage 3 if RRT was used
valid['rrt'] = valid['hadm_id'].isin(rrt_hadms).astype(int)
valid.loc[valid['rrt'] == 1, 'kdigo_stage'] = valid.loc[valid['rrt'] == 1, 'kdigo_stage'].clip(lower=3)

print(f"\n  After RRT upgrade:")
for s in [0, 1, 2, 3]:
    n = (valid['kdigo_stage'] == s).sum()
    n_rrt = ((valid['kdigo_stage'] == s) & (valid['rrt'] == 1)).sum()
    pct = n / len(valid) * 100
    death_n = (valid[valid['kdigo_stage'] == s]['hospital_expire_flag'] == 1).sum()
    death_pct = death_n / n * 100 if n > 0 else 0
    print(f"    Stage {s}: {n:,} ({pct:.1f}%) [RRT={n_rrt}], mortality={death_n} ({death_pct:.1f}%)")

# ============================================================
# 8. Save cohort
# ============================================================
print("\n[8] Saving cohort...")

output_path = os.path.join(OUT, "kdigo_cohort.csv")
valid.to_csv(output_path, index=False)
print(f"  Saved to: {output_path}")
print(f"  Rows: {len(valid):,}")
print(f"  Columns: {list(valid.columns)}")

# Summary stats
summary = {
    'total_icu_stays': len(valid),
    'stage_0': int((valid['kdigo_stage'] == 0).sum()),
    'stage_1': int((valid['kdigo_stage'] == 1).sum()),
    'stage_2': int((valid['kdigo_stage'] == 2).sum()),
    'stage_3': int((valid['kdigo_stage'] == 3).sum()),
    'overall_mortality': float(valid['hospital_expire_flag'].mean()),
    'stage_mortality': {
        '0': float(valid[valid['kdigo_stage']==0]['hospital_expire_flag'].mean()),
        '1': float(valid[valid['kdigo_stage']==1]['hospital_expire_flag'].mean()),
        '2': float(valid[valid['kdigo_stage']==2]['hospital_expire_flag'].mean()),
        '3': float(valid[valid['kdigo_stage']==3]['hospital_expire_flag'].mean()),
    },
    'mean_baseline_cr': float(valid['baseline_cr'].mean()),
    'mean_peak_cr': float(valid['peak_cr'].mean()),
    'n_rrt': int(valid['rrt'].sum()),
}

with open(os.path.join(OUT, "kdigo_summary.json"), 'w') as f:
    json.dump(summary, f, indent=2, default=str)

print("\n" + "=" * 70)
print("Step 1 Complete!")
print("=" * 70)
print(json.dumps(summary, indent=2))
print(f"\nEnd: {datetime.now()}")
