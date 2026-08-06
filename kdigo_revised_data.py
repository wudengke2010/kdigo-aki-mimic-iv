#!/usr/bin/env python3
"""
KDIGO AKI 分期 — 修订版全量数据提取
=====================================
审稿修订内容:
  Fix #1: 生命体征生理范围过滤 (HR 20-250, SBP 40-250, DBP 20-200, etc.)
  Fix #2: 种族 White 解析修复
  Fix #3: 不再同时放入 KDIGO 分期 + Cr ratio (此脚本只提取数据)
  Fix #4: RRT-stratified Stage 3 分层标记
  Fix #5: MDRD 估算基线肌酐 (eGFR=75 反推)
  Fix #6: 实验室 CKD 定义 (CKD-EPI eGFR < 60)
  Fix #7: GCS 评分提取
  Fix #8: 机械通气状态
  Fix #9: 血管活性药物使用
  Fix #10: 患者筛选流程图数据
  Fix #11: 缺失数据处理报告

产出: kdigo_cohort_revised.csv + kdigo_cohort_revised_enriched.csv + table1_revised.csv
"""
import pandas as pd
import numpy as np
import os, json, gc
from datetime import datetime

print("=" * 80)
print("KDIGO-AKI REVISED Data Extraction (All Reviewer Fixes)")
print(f"Start: {datetime.now()}")
print("=" * 80)

BASE = "E:/mimic-iv"
HOSP = os.path.join(BASE, "v3.1/physionet.org/files/mimiciv/3.1/hosp")
ICU = os.path.join(BASE, "v3.1/physionet.org/files/mimiciv/3.1/icu")
OUT = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"

def read_gz(path, **kwargs):
    return pd.read_csv(path, compression='gzip', low_memory=False, **kwargs)

# ============================================================
# PHASE 0: 患者筛选流程图
# ============================================================
print("\n" + "=" * 60)
print("PHASE 0: Patient Selection Flowchart")
print("=" * 60)

icustays_all = read_gz(os.path.join(ICU, "icustays.csv.gz"))
patients_all = read_gz(os.path.join(HOSP, "patients.csv.gz"))
admissions_all = read_gz(os.path.join(HOSP, "admissions.csv.gz"))

n_total_icu = len(icustays_all)
n_total_patients = icustays_all['subject_id'].nunique()
n_total_hadm = icustays_all['hadm_id'].nunique()

print(f"  Total ICU stays: {n_total_icu:,}")
print(f"  Total ICU patients: {n_total_patients:,}")
print(f"  Total hospital admissions: {n_total_hadm:,}")

# 排除: 年龄 < 18（不修改 patients_all，避免后续 join 膨胀）
print(f"  (ICU patients before age filter: {n_total_patients:,})")
admittime_ref = admissions_all.groupby('subject_id')['admittime'].first().reset_index()
admittime_ref.columns = ['subject_id', 'first_admittime']
pat_age = patients_all[['subject_id', 'anchor_age', 'anchor_year']].drop_duplicates('subject_id')
pat_age = pat_age.merge(admittime_ref, on='subject_id', how='left')
pat_age['first_admittime'] = pd.to_datetime(pat_age['first_admittime'])
pat_age['age_at_first'] = pat_age['anchor_age'] + (pat_age['first_admittime'].dt.year - pat_age['anchor_year'])
icu_patient_ids = set(icustays_all['subject_id'].unique())
adult_non_icu = set(pat_age[pat_age['age_at_first'] >= 18]['subject_id'].unique())
adult_ids = adult_non_icu & icu_patient_ids  # only ICU patients
n_adult_icu = len(adult_ids)
n_excluded_age = n_total_patients - n_adult_icu
print(f"  Excluded age < 18: {n_excluded_age:,} ICU patients")

# First ICU per hadm
icustays_all['intime'] = pd.to_datetime(icustays_all['intime'])
icu_first = icustays_all.sort_values('intime').drop_duplicates('hadm_id', keep='first')
n_first_icu = len(icu_first)
n_excluded_duplicate = icustays_all['hadm_id'].nunique() - n_first_icu  # stays beyond first per hadm
print(f"  After first ICU per hadm: {n_first_icu:,} stays (excluded {n_excluded_duplicate:,} duplicate stays)")

# Filter adults
icu_adult = icu_first[icu_first['subject_id'].isin(adult_ids)]
n_adult_stays = len(icu_adult)
n_excluded_age_stays = n_first_icu - n_adult_stays
print(f"  After age >= 18: {n_adult_stays:,} stays (excluded {n_excluded_age_stays:,})")

# ============================================================
# PHASE 1: 基础数据加载 + ICD 编码
# ============================================================
print("\n" + "=" * 60)
print("PHASE 1: Base tables + ICD coding")
print("=" * 60)

diagnoses = read_gz(os.path.join(HOSP, "diagnoses_icd.csv.gz"))
d_icd = read_gz(os.path.join(HOSP, "d_icd_diagnoses.csv.gz"))
d_lab = read_gz(os.path.join(HOSP, "d_labitems.csv.gz"))

# AKI/CKD ICD
aki10 = diagnoses[(diagnoses['icd_version'] == 10) & diagnoses['icd_code'].str.match(r'^N17', na=False)]
aki9 = diagnoses[(diagnoses['icd_version'] == 9) & diagnoses['icd_code'].str.match(r'^584', na=False)]
ckd10 = diagnoses[(diagnoses['icd_version'] == 10) & diagnoses['icd_code'].str.match(r'^N18', na=False)]
ckd9 = diagnoses[(diagnoses['icd_version'] == 9) & diagnoses['icd_code'].str.match(r'^585', na=False)]

aki_ids = set(aki10['subject_id'].unique()) | set(aki9['subject_id'].unique())
ckd_icd_ids = set(ckd10['subject_id'].unique()) | set(ckd9['subject_id'].unique())

print(f"  AKI (ICD): {len(aki_ids):,} patients")
print(f"  CKD (ICD): {len(ckd_icd_ids):,} patients")

# ============================================================
# PHASE 2: 合并 ICU + 入院 + 患者信息
# ============================================================
print("\n" + "=" * 60)
print("PHASE 2: Build ICU cohort")
print("=" * 60)

icu = icu_adult.copy()

adm_cols = ['hadm_id', 'subject_id', 'admittime', 'dischtime', 'deathtime',
            'admission_type', 'admission_location', 'discharge_location',
            'insurance', 'language', 'marital_status', 'race',
            'hospital_expire_flag', 'edregtime', 'edouttime']
icu = icu.merge(admissions_all[adm_cols], on=['subject_id', 'hadm_id'], how='left')

pat_cols = ['subject_id', 'gender', 'anchor_age', 'anchor_year', 'dod']
# Use clean unduplicated patients table (patients_all was NOT mutated above)
icu = icu.merge(patients_all[pat_cols].drop_duplicates('subject_id'), on='subject_id', how='left')

# Time columns
for col in ['intime', 'outtime', 'admittime', 'dischtime', 'deathtime', 'dod']:
    if col in icu.columns:
        icu[col] = pd.to_datetime(icu[col])

# LOS
icu['icu_los_hours'] = (icu['outtime'] - icu['intime']).dt.total_seconds() / 3600
icu['hosp_los_days'] = (icu['dischtime'] - icu['admittime']).dt.total_seconds() / 86400

# Renal group labels
icu['aki_icd'] = icu['subject_id'].isin(aki_ids).astype(int)
icu['ckd_icd'] = icu['subject_id'].isin(ckd_icd_ids).astype(int)

n_before_cr = len(icu)
print(f"  Adult ICU cohort (before Cr filtering): {n_before_cr:,}")

# ============================================================
# PHASE 3: 肌酐提取 + KDIGO 分期
# ============================================================
print("\n" + "=" * 60)
print("PHASE 3: Creatinine extraction + KDIGO staging")
print("=" * 60)

creatinine_itemid = 50912
icu_subject_ids = set(icu['subject_id'].unique())
labevents_path = os.path.join(HOSP, "labevents.csv.gz")

cr_chunks = []
chunk_size = 2000000
total_rows = 0
matched_rows = 0

for i, chunk in enumerate(pd.read_csv(labevents_path, compression='gzip', chunksize=chunk_size, low_memory=False)):
    total_rows += len(chunk)
    chunk = chunk[(chunk['itemid'] == creatinine_itemid) & chunk['subject_id'].isin(icu_subject_ids)].copy()
    if len(chunk) > 0:
        chunk = chunk[['subject_id', 'hadm_id', 'charttime', 'valuenum']].copy()
        chunk['valuenum'] = pd.to_numeric(chunk['valuenum'], errors='coerce')
        chunk = chunk.dropna(subset=['valuenum'])
        chunk = chunk[(chunk['valuenum'] > 0.1) & (chunk['valuenum'] < 100)]
        cr_chunks.append(chunk)
        matched_rows += len(chunk)
    if (i + 1) % 10 == 0:
        print(f"  ... scanned {total_rows/1e6:.1f}M rows, kept {matched_rows:,} Cr")

cr_all = pd.concat(cr_chunks, ignore_index=True)
del cr_chunks; gc.collect()
cr_all['charttime'] = pd.to_datetime(cr_all['charttime'])
cr_all = cr_all.sort_values(['subject_id', 'hadm_id', 'charttime'])
print(f"  Total Cr measurements: {len(cr_all):,}")

# Merge with ICU times
cr_icu = cr_all.merge(icu[['hadm_id', 'intime', 'outtime']], on='hadm_id', how='inner')

# Baseline Cr: first hospital measurement
baseline_cr = cr_all.groupby('hadm_id').first().reset_index()[['hadm_id', 'valuenum']]
baseline_cr.columns = ['hadm_id', 'baseline_cr']

# Peak Cr during ICU
peak_cr = cr_icu.groupby('hadm_id')['valuenum'].max().reset_index()
peak_cr.columns = ['hadm_id', 'peak_cr']

# First/last ICU Cr
cr_icu['hours_from_icu_in'] = (cr_icu['charttime'] - cr_icu['intime']).dt.total_seconds() / 3600
first_icu_cr = cr_icu[cr_icu['hours_from_icu_in'] >= -24].groupby('hadm_id').first().reset_index()[['hadm_id', 'valuenum']]
first_icu_cr.columns = ['hadm_id', 'icu_first_cr']
last_icu_cr = cr_icu.groupby('hadm_id').last().reset_index()[['hadm_id', 'valuenum']]
last_icu_cr.columns = ['hadm_id', 'icu_last_cr']
cr_count = cr_icu.groupby('hadm_id').size().reset_index(name='n_cr_measurements')

# Merge
icu = icu.merge(baseline_cr, on='hadm_id', how='left')
icu = icu.merge(peak_cr, on='hadm_id', how='left')
icu = icu.merge(first_icu_cr, on='hadm_id', how='left')
icu = icu.merge(last_icu_cr, on='hadm_id', how='left')
icu = icu.merge(cr_count, on='hadm_id', how='left')

# KDIGO staging
valid = icu[icu['baseline_cr'].notna() & icu['peak_cr'].notna()].copy()
valid['cr_ratio'] = valid['peak_cr'] / valid['baseline_cr']
valid['cr_increase'] = valid['peak_cr'] - valid['baseline_cr']

conditions = [
    (valid['cr_ratio'] >= 3.0) | (valid['peak_cr'] >= 4.0),
    (valid['cr_ratio'] >= 2.0) & (valid['cr_ratio'] < 3.0),
    ((valid['cr_ratio'] >= 1.5) & (valid['cr_ratio'] < 2.0)) | (valid['cr_increase'] >= 0.3),
]
choices = [3, 2, 1]
valid['kdigo_stage'] = np.select(conditions, choices, default=0)

n_after_cr = len(valid)
n_excluded_no_cr = n_before_cr - n_after_cr
print(f"  Cohort with valid Cr: {n_after_cr:,} (excluded {n_excluded_no_cr:,} missing Cr)")

# Stage 3 criteria breakdown
valid['stage3_cr_criteria'] = ((valid['cr_ratio'] >= 3.0) | (valid['peak_cr'] >= 4.0)).astype(int)

# ============================================================
# PHASE 4: RRT + Stage 3 分层
# ============================================================
print("\n" + "=" * 60)
print("PHASE 4: RRT + Stage 3 stratification")
print("=" * 60)

proc_events = read_gz(os.path.join(ICU, "procedureevents.csv.gz"))
d_items = read_gz(os.path.join(ICU, "d_items.csv.gz"))

dial_items = d_items[d_items['label'].str.lower().str.contains(
    'dialysis - crrt|dialysis - cvvhd|dialysis - cvvhdf|hemodialysis|dialysis - scuf',
    na=False)]
dial_itemids = set(dial_items['itemid'].tolist())
print(f"  RRT itemids: {len(dial_itemids)}")

rrt_events = proc_events[proc_events['itemid'].isin(dial_itemids)]
rrt_hadms = set(rrt_events['hadm_id'].unique())
rrt_subjects = set(rrt_events['subject_id'].unique())
print(f"  RRT hospital stays: {len(rrt_hadms):,}")
print(f"  RRT patients: {len(rrt_subjects):,}")

valid['rrt'] = valid['hadm_id'].isin(rrt_hadms).astype(int)

# Stage 3 entry mechanism
valid['stage3_by_rrt_only'] = 0  # RRT but no Cr criteria
valid['stage3_by_cr_only'] = 0   # Cr criteria but no RRT
valid['stage3_by_both'] = 0      # Both

# Identify Stage 3 by mechanism (before RRT upgrade)
stage3_pre_rrt = valid['kdigo_stage'].copy()  # Cr-based before RRT

# Apply RRT upgrade
valid.loc[valid['rrt'] == 1, 'kdigo_stage'] = valid.loc[valid['rrt'] == 1, 'kdigo_stage'].clip(lower=3)

# Now classify Stage 3 patients
stage3_mask = valid['kdigo_stage'] == 3
valid.loc[stage3_mask & (valid['rrt'] == 1) & (valid['stage3_cr_criteria'] == 0), 'stage3_by_rrt_only'] = 1
valid.loc[stage3_mask & (valid['rrt'] == 0) & (valid['stage3_cr_criteria'] == 1), 'stage3_by_cr_only'] = 1
valid.loc[stage3_mask & (valid['rrt'] == 1) & (valid['stage3_cr_criteria'] == 1), 'stage3_by_both'] = 1

print(f"\n  Stage 3 stratification:")
print(f"    RRT only (no Cr criteria): {(valid['stage3_by_rrt_only']==1).sum():,}")
print(f"    Cr only (no RRT): {(valid['stage3_by_cr_only']==1).sum():,}")
print(f"    Both RRT + Cr: {(valid['stage3_by_both']==1).sum():,}")

# ============================================================
# PHASE 5: MDRD baseline Cr + Lab-based CKD
# ============================================================
print("\n" + "=" * 60)
print("PHASE 5: MDRD baseline + Lab-based CKD")
print("=" * 60)

# MDRD: eGFR = 175 × SCr^(-1.154) × Age^(-0.203) × [0.742 if female] × [1.212 if Black]
# Back-calculate: baseline_SCr = (75 / (175 × Age^(-0.203) × [0.742 if F] × [1.212 if B]))^(-1/1.154)
# This is equivalent to: SCr = (75 / (175 * A * R * S))^(-1/1.154)
# where A = age^(-0.203), R = 1.212 if Black else 1, S = 0.742 if female else 1

def mdrd_baseline_cr(age, sex, race_str):
    age_factor = max(age, 18) ** (-0.203)
    sex_factor = 0.742 if sex == 'F' else 1.0
    race_factor = 1.212 if (isinstance(race_str, str) and 'black' in race_str.lower()) else 1.0
    denominator = 175 * age_factor * race_factor * sex_factor
    if denominator > 0:
        return (75.0 / denominator) ** (-1.0 / 1.154)
    return np.nan

valid['mdrd_baseline_cr'] = valid.apply(
    lambda r: mdrd_baseline_cr(r['anchor_age'], r.get('gender', 'M'), r.get('race', '')),
    axis=1
)
mdrd_available = valid['mdrd_baseline_cr'].notna().sum()
print(f"  MDRD baseline Cr available: {mdrd_available:,}/{len(valid):,}")

# MDRD-based KDIGO staging
valid['cr_ratio_mdrd'] = valid['peak_cr'] / valid['mdrd_baseline_cr']
valid['kdigo_stage_mdrd'] = np.select([
    (valid['cr_ratio_mdrd'] >= 3.0) | (valid['peak_cr'] >= 4.0),
    (valid['cr_ratio_mdrd'] >= 2.0) & (valid['cr_ratio_mdrd'] < 3.0),
    ((valid['cr_ratio_mdrd'] >= 1.5) & (valid['cr_ratio_mdrd'] < 2.0)) | ((valid['peak_cr'] - valid['mdrd_baseline_cr']) >= 0.3),
], [3, 2, 1], default=0)

# CKD-EPI eGFR from first SCr
def ckd_epi_egfr(cr, age, sex, race_str):
    """CKD-EPI 2021 (race-free)"""
    if pd.isna(cr) or pd.isna(age) or cr <= 0:
        return np.nan
    kappa = 0.7 if sex == 'F' else 0.9
    alpha = -0.241 if sex == 'F' else -0.302
    sex_factor = 1.012 if sex == 'F' else 1.0
    egfr = 142 * (min(cr / kappa, 1.0) ** alpha) * (max(cr / kappa, 1.0) ** (-1.200)) * (0.9938 ** age) * sex_factor
    return egfr

valid['egfr_ckdepi'] = valid.apply(
    lambda r: ckd_epi_egfr(r['baseline_cr'], r['anchor_age'], r.get('gender', 'M'), r.get('race', '')),
    axis=1
)
valid['ckd_lab'] = (valid['egfr_ckdepi'] < 60).astype(int)
ckd_lab_n = valid['ckd_lab'].sum()
ckd_icd_n = valid['ckd_icd'].sum()
ckd_both = ((valid['ckd_icd'] == 1) & (valid['ckd_lab'] == 1)).sum()
ckd_icd_only = ((valid['ckd_icd'] == 1) & (valid['ckd_lab'] == 0)).sum()
ckd_lab_only = ((valid['ckd_icd'] == 0) & (valid['ckd_lab'] == 1)).sum()

print(f"  CKD by ICD: {ckd_icd_n:,} ({ckd_icd_n/len(valid)*100:.1f}%)")
print(f"  CKD by lab (eGFR<60): {ckd_lab_n:,} ({ckd_lab_n/len(valid)*100:.1f}%)")
print(f"  Both: {ckd_both:,} | ICD only: {ckd_icd_only:,} | Lab only: {ckd_lab_only:,}")
print(f"  ICD sensitivity (vs lab): {ckd_both/ckd_lab_n*100:.1f}%" if ckd_lab_n > 0 else "")

# ============================================================
# PHASE 6: 生命体征 (physiological range filtered)
# ============================================================
print("\n" + "=" * 60)
print("PHASE 6: Vital signs (physiological range filtered)")
print("=" * 60)

chartevents_path = os.path.join(ICU, "chartevents.csv.gz")
cohort_hadm_ids = set(valid['hadm_id'].unique())

# Itemids
vital_itemids = {
    'heart_rate': [220045],
    'sbp': [220050],
    'dbp': [220051],
    'mbp': [220052],
    'resp_rate': [220210],
    'temperature': [223761],
    'spo2': [220277],
}

# Physiological range filters
vital_ranges = {
    'heart_rate': (20, 250),
    'sbp': (40, 250),
    'dbp': (20, 200),
    'mbp': (30, 250),
    'resp_rate': (5, 60),
    'temperature': (30, 42),  # Celsius (after F→C conversion)
    'spo2': (50, 100),
}

vital_outlier_counts = {name: {'raw': 0, 'kept': 0, 'excluded': 0} for name in vital_itemids}

all_vital_ids = [item for sublist in vital_itemids.values() for item in sublist]
vital_data = []
chunk_size = 2000000

for i, chunk in enumerate(pd.read_csv(chartevents_path, compression='gzip', chunksize=chunk_size, low_memory=False)):
    chunk = chunk[(chunk['itemid'].isin(all_vital_ids)) & (chunk['hadm_id'].isin(cohort_hadm_ids))].copy()
    if len(chunk) > 0:
        chunk['valuenum'] = pd.to_numeric(chunk['valuenum'], errors='coerce')
        chunk = chunk.dropna(subset=['valuenum'])
        vital_data.append(chunk[['hadm_id', 'stay_id', 'itemid', 'valuenum', 'charttime']])
    if (i + 1) % 5 == 0:
        print(f"  ... processed ~{chunk_size*(i+1)/1e6:.0f}M chartevents rows, kept {sum(len(v) for v in vital_data):,} vital records")

vitals_all = pd.concat(vital_data, ignore_index=True)
del vital_data; gc.collect()
print(f"  Raw vital records: {len(vitals_all):,}")

# Merge with ICU times (chartevents already has stay_id; merge on both to avoid column suffix)
vitals_all['charttime'] = pd.to_datetime(vitals_all['charttime'])
vitals_wtime = vitals_all.merge(valid[['hadm_id', 'stay_id', 'intime']], on=['hadm_id', 'stay_id'], how='inner')
vitals_wtime['intime'] = pd.to_datetime(vitals_wtime['intime'])
vitals_wtime['hours_from_icu'] = (vitals_wtime['charttime'] - vitals_wtime['intime']).dt.total_seconds() / 3600

# First 24h only
vitals_24h = vitals_wtime[(vitals_wtime['hours_from_icu'] >= 0) & (vitals_wtime['hours_from_icu'] <= 24)].copy()

# Apply physiological range filters + aggregate
vital_summary = {}
for name, itemids in vital_itemids.items():
    sub = vitals_24h[vitals_24h['itemid'].isin(itemids)].copy()
    raw_count = len(sub)
    
    # Temperature conversion F->C
    if name == 'temperature':
        sub['valuenum'] = (sub['valuenum'] - 32) * 5 / 9
    
    # Apply range filter
    lo, hi = vital_ranges[name]
    sub_filtered = sub[(sub['valuenum'] >= lo) & (sub['valuenum'] <= hi)]
    kept_count = len(sub_filtered)
    
    vital_outlier_counts[name]['raw'] = raw_count
    vital_outlier_counts[name]['kept'] = kept_count
    vital_outlier_counts[name]['excluded'] = raw_count - kept_count
    
    # Aggregate per stay
    agg = sub_filtered.groupby('stay_id')['valuenum'].agg(['mean', 'min', 'max', 'first', 'count']).reset_index()
    agg.columns = ['stay_id', f'{name}_mean', f'{name}_min', f'{name}_max', f'{name}_first', f'{name}_count']
    vital_summary[name] = agg
    
    excl_pct = (raw_count - kept_count) / max(raw_count, 1) * 100
    print(f"  {name}: {kept_count:,}/{raw_count:,} kept ({excl_pct:.1f}% excluded)")

# Merge vitals into cohort
for name, df in vital_summary.items():
    valid = valid.merge(df, on='stay_id', how='left')

# ============================================================
# PHASE 7: GCS + 机械通气 + 血管活性药物
# ============================================================
print("\n" + "=" * 60)
print("PHASE 7: GCS + Ventilation + Vasopressors")
print("=" * 60)

# GCS from chartevents (itemids: 220739 Eye, 223900 Verbal, 223901 Motor)
gcs_itemids = [220739, 223900, 223901]
gcs_components = {220739: 'gcs_eye', 223900: 'gcs_verbal', 223901: 'gcs_motor'}

gcs_data = []
for i, chunk in enumerate(pd.read_csv(chartevents_path, compression='gzip', chunksize=chunk_size, low_memory=False)):
    chunk = chunk[(chunk['itemid'].isin(gcs_itemids)) & (chunk['hadm_id'].isin(cohort_hadm_ids))].copy()
    if len(chunk) > 0:
        chunk['valuenum'] = pd.to_numeric(chunk['valuenum'], errors='coerce')
        chunk = chunk.dropna(subset=['valuenum'])
        chunk = chunk[chunk['valuenum'].between(1, 5)]
        gcs_data.append(chunk[['hadm_id', 'stay_id', 'itemid', 'valuenum', 'charttime']])
    if (i + 1) % 5 == 0:
        print(f"  GCS scan: processed ~{chunk_size*(i+1)/1e6:.0f}M rows, kept {sum(len(g) for g in gcs_data):,} GCS records")

gcs_all = pd.concat(gcs_data, ignore_index=True)
del gcs_data; gc.collect()
gcs_all['charttime'] = pd.to_datetime(gcs_all['charttime'])

# Merge with ICU times, first 24h
gcs_wtime = gcs_all.merge(valid[['hadm_id', 'stay_id', 'intime']], on=['hadm_id', 'stay_id'], how='inner')
gcs_wtime['intime'] = pd.to_datetime(gcs_wtime['intime'])
gcs_wtime['hours'] = (gcs_wtime['charttime'] - gcs_wtime['intime']).dt.total_seconds() / 3600
gcs_24h = gcs_wtime[(gcs_wtime['hours'] >= 0) & (gcs_wtime['hours'] <= 24)]

# First GCS per component
for itemid, comp_name in gcs_components.items():
    comp = gcs_24h[gcs_24h['itemid'] == itemid]
    first_comp = comp.sort_values('hours').groupby('stay_id').first().reset_index()
    valid = valid.merge(first_comp[['stay_id', 'valuenum']].rename(columns={'valuenum': f'{comp_name}_first'}), on='stay_id', how='left')

# GCS total (sum of eye + verbal + motor, first measurement)
gcs_cols = [c for c in valid.columns if c.startswith('gcs_') and c.endswith('_first')]
if len(gcs_cols) >= 3:
    valid['gcs_total'] = valid[gcs_cols].sum(axis=1, min_count=3)
    gcs_avail = valid['gcs_total'].notna().sum()
    print(f"  GCS total available: {gcs_avail:,}/{len(valid):,} ({gcs_avail/len(valid)*100:.1f}%)")
    print(f"  GCS mean: {valid['gcs_total'].mean():.1f} ± {valid['gcs_total'].std():.1f}")

# Ventilation
vent_items = d_items[d_items['label'].str.lower().str.contains(
    'invasive ventilation|mechanically ventilated|ventilated at any time',
    na=False)]
vent_itemids = set(vent_items['itemid'].tolist())
print(f"  Ventilation itemids: {len(vent_itemids)}")

# Combine procedureevents (vent from procedures) + chartevents (vent from chart)
vent_patients = set()

# From procedureevents
vent_proc = proc_events[proc_events['itemid'].isin(vent_itemids)]
vent_patients |= set(vent_proc['subject_id'].unique())

# From chartevents
vent_chart = []
for i, chunk in enumerate(pd.read_csv(chartevents_path, compression='gzip', chunksize=chunk_size, low_memory=False)):
    chunk = chunk[(chunk['itemid'].isin(vent_itemids)) & (chunk['hadm_id'].isin(cohort_hadm_ids))]
    if len(chunk) > 0:
        vent_patients |= set(chunk['subject_id'].unique())
    if (i + 1) % 10 == 0:
        print(f"  Vent scan: ~{chunk_size*(i+1)/1e6:.0f}M rows")

valid['mech_vent'] = valid['subject_id'].isin(vent_patients).astype(int)
vent_n = valid['mech_vent'].sum()
print(f"  Mechanical ventilation: {vent_n:,}/{len(valid):,} ({vent_n/len(valid)*100:.1f}%)")

# Vasopressors from inputevents
vaso_itemids = [221906, 221289, 221662, 222315, 221749, 221653]  # NE, Epi, Dopamine, Vasopressin, Phenylephrine, Dobutamine
inputevents_path = os.path.join(ICU, "inputevents.csv.gz")

vaso_patients = set()
for i, chunk in enumerate(pd.read_csv(inputevents_path, compression='gzip', chunksize=chunk_size, low_memory=False)):
    chunk = chunk[(chunk['itemid'].isin(vaso_itemids)) & (chunk['hadm_id'].isin(cohort_hadm_ids))]
    if len(chunk) > 0:
        vaso_patients |= set(chunk['subject_id'].unique())
    if (i + 1) % 5 == 0:
        print(f"  Vaso scan: ~{chunk_size*(i+1)/1e6:.0f}M rows")

valid['vasopressor'] = valid['subject_id'].isin(vaso_patients).astype(int)
vaso_n = valid['vasopressor'].sum()
print(f"  Vasopressor use: {vaso_n:,}/{len(valid):,} ({vaso_n/len(valid)*100:.1f}%)")

# ============================================================
# PHASE 8: 实验室指标
# ============================================================
print("\n" + "=" * 60)
print("PHASE 8: Admission labs")
print("=" * 60)

key_labs = {
    'bun': 51006, 'egfr': 50920, 'sodium': 50983, 'potassium': 50971,
    'bicarbonate': 50882, 'lactate': 50813, 'wbc': 51301, 'hemoglobin': 51222,
    'platelet': 51265, 'bilirubin_total': 50885, 'albumin': 50862, 'glucose': 50931,
}

all_lab_ids = list(key_labs.values())
lab_data = []

for i, chunk in enumerate(pd.read_csv(labevents_path, compression='gzip', chunksize=chunk_size, low_memory=False)):
    chunk = chunk[(chunk['itemid'].isin(all_lab_ids)) & (chunk['hadm_id'].isin(cohort_hadm_ids))].copy()
    if len(chunk) > 0:
        chunk['valuenum'] = pd.to_numeric(chunk['valuenum'], errors='coerce')
        chunk = chunk.dropna(subset=['valuenum'])
        lab_data.append(chunk[['subject_id', 'hadm_id', 'itemid', 'valuenum', 'charttime']])
    if (i + 1) % 10 == 0:
        print(f"  Lab scan: ~{chunk_size*(i+1)/1e6:.0f}M rows")

labs_all = pd.concat(lab_data, ignore_index=True)
del lab_data; gc.collect()
labs_all['charttime'] = pd.to_datetime(labs_all['charttime'])

# First measurement within +/- 24h
# labevents has no stay_id column; merge on hadm_id (1:1 since first ICU per hadm)
labs_wtime = labs_all.merge(valid[['hadm_id', 'stay_id', 'intime']], on='hadm_id', how='inner')
labs_wtime['intime'] = pd.to_datetime(labs_wtime['intime'])
labs_wtime['hours'] = (labs_wtime['charttime'] - labs_wtime['intime']).dt.total_seconds() / 3600
labs_48h = labs_wtime[(labs_wtime['hours'] >= -24) & (labs_wtime['hours'] <= 24)]
labs_48h = labs_48h.sort_values('hours')

for name, itemid in key_labs.items():
    sub = labs_48h[labs_48h['itemid'] == itemid]
    first_val = sub.groupby('stay_id').first().reset_index()
    valid = valid.merge(first_val[['stay_id', 'valuenum']].rename(columns={'valuenum': f'lab_{name}'}), on='stay_id', how='left')
    n_avail = valid[f'lab_{name}'].notna().sum()
    print(f"  lab_{name}: {n_avail:,}/{len(valid):,} ({n_avail/len(valid)*100:.1f}%)")

# ============================================================
# PHASE 9: 合并症 (Elixhauser)
# ============================================================
print("\n" + "=" * 60)
print("PHASE 9: Comorbidities (Elixhauser)")
print("=" * 60)

valid_sids = set(valid['subject_id'].unique())
diag_cohort = diagnoses[diagnoses['subject_id'].isin(valid_sids)]

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

for comorb, patients_set in comorb_status.items():
    col_name = f'com_{comorb.replace(" ", "_").replace("(", "").replace(")", "")}'
    valid[col_name] = valid['subject_id'].isin(patients_set).astype(int)

print(f"  Comorbidities added: {len(comorb_status)} conditions")

# ============================================================
# PHASE 10: Derived columns + Table 1 generation
# ============================================================
print("\n" + "=" * 60)
print("PHASE 10: Derived columns + Table 1")
print("=" * 60)

# Demographics
valid['male'] = valid['gender'].map({'M': 1, 'F': 0}).fillna(0).astype(int)

# Race: robust parsing
def parse_race(race_str):
    if pd.isna(race_str):
        return 'Unknown'
    race_str = str(race_str).strip().upper()
    if 'WHITE' in race_str:
        return 'White'
    elif 'BLACK' in race_str or 'AFRICAN' in race_str:
        return 'Black'
    elif 'ASIAN' in race_str:
        return 'Asian'
    elif 'HISPANIC' in race_str or 'LATINO' in race_str:
        return 'Hispanic'
    else:
        return 'Other/Unknown'

valid['race_parsed'] = valid['race'].apply(parse_race)
valid['race_white'] = (valid['race_parsed'] == 'White').astype(int)
valid['race_black'] = (valid['race_parsed'] == 'Black').astype(int)
valid['race_asian'] = (valid['race_parsed'] == 'Asian').astype(int)
valid['race_hispanic'] = (valid['race_parsed'] == 'Hispanic').astype(int)

print(f"  Race distribution:")
for r in ['White', 'Black', 'Asian', 'Hispanic', 'Other/Unknown']:
    n = (valid['race_parsed'] == r).sum()
    print(f"    {r}: {n:,} ({n/len(valid)*100:.1f}%)")

# SOFA components (simplified)
# Respiratory: mechanical ventilation = 3-4 (approximate)
# Coagulation: platelet < 20k = 4, <50k = 3, <100k = 2, <150k = 1
# Liver: bilirubin > 12 = 4, 6-12 = 3, 2-6 = 2, 1.2-2 = 1
# Cardiovascular: vasopressor = 3-4 (approximate), MAP < 70 = 1
# Neurological: GCS 13-14 = 1, 10-12 = 2, 6-9 = 3, <6 = 4
# Renal: Cr > 5 = 4, 3.5-5 = 3, 2-3.5 = 2, 1.2-2 = 1 (or urine output)

# We'll derive a simplified SOFA for model adjustment
if 'gcs_total' in valid.columns:
    valid['sofa_neuro'] = np.select([
        valid['gcs_total'] < 6, valid['gcs_total'].between(6, 9),
        valid['gcs_total'].between(10, 12), valid['gcs_total'].between(13, 14),
    ], [4, 3, 2, 1], default=0)

if 'lab_platelet' in valid.columns:
    valid['sofa_coag'] = np.select([
        valid['lab_platelet'] < 20, valid['lab_platelet'].between(20, 49),
        valid['lab_platelet'].between(50, 99), valid['lab_platelet'].between(100, 149),
    ], [4, 3, 2, 1], default=0)

if 'lab_bilirubin_total' in valid.columns:
    valid['sofa_liver'] = np.select([
        valid['lab_bilirubin_total'] > 12, valid['lab_bilirubin_total'].between(6, 12),
        valid['lab_bilirubin_total'].between(2, 6), valid['lab_bilirubin_total'].between(1.2, 2),
    ], [4, 3, 2, 1], default=0)

# Cardiovascular: vasopressor = min 3, MAP < 70 = 1
valid['sofa_cv'] = 0
valid.loc[(valid['vasopressor'] == 1), 'sofa_cv'] = 3
valid.loc[(valid['vasopressor'] == 0) & (valid['mbp_mean'].notna()) & (valid['mbp_mean'] < 70), 'sofa_cv'] = 1

# Renal: based on peak Cr
valid['sofa_renal'] = np.select([
    valid['peak_cr'] > 5.0, valid['peak_cr'].between(3.5, 5.0),
    valid['peak_cr'].between(2.0, 3.5), valid['peak_cr'].between(1.2, 2.0),
], [4, 3, 2, 1], default=0)
valid.loc[valid['rrt'] == 1, 'sofa_renal'] = 4

# Resp: mechanical ventilation >= 3
valid['sofa_resp'] = 0
valid.loc[valid['mech_vent'] == 1, 'sofa_resp'] = 3

# Total SOFA
sofa_cols = ['sofa_neuro', 'sofa_coag', 'sofa_liver', 'sofa_cv', 'sofa_renal', 'sofa_resp']
valid['sofa_total'] = valid[sofa_cols].sum(axis=1, min_count=3)
sofa_avail = valid['sofa_total'].notna().sum()
print(f"  SOFA total available: {sofa_avail:,}/{len(valid):,} ({sofa_avail/len(valid)*100:.1f}%)")
print(f"  SOFA mean: {valid['sofa_total'].mean():.1f} ± {valid['sofa_total'].std():.1f}")

# Save revised cohort
valid.to_csv(os.path.join(OUT, "kdigo_cohort_revised.csv"), index=False)
print(f"\n  Revised cohort saved: {len(valid):,} rows, {len(valid.columns)} columns")

# ============================================================
# PHASE 11: Table 1 generation (with fixes)
# ============================================================
print("\n" + "=" * 60)
print("PHASE 11: Table 1")
print("=" * 60)

def describe_var_cont(df, col, decimals=1, drop_outliers=False):
    results = {}
    for stage in [0, 1, 2, 3]:
        sub = df[df['kdigo_stage'] == stage]
        vals = sub[col].dropna()
        if drop_outliers:
            vals = vals[(vals >= vals.quantile(0.01)) & (vals <= vals.quantile(0.99))]
        if len(vals) < 5:
            results[f'Stage_{stage}'] = "N/A"
        else:
            results[f'Stage_{stage}'] = f"{vals.mean():.{decimals}f} ± {vals.std():.{decimals}f}"
    vals = df[col].dropna()
    if drop_outliers:
        vals = vals[(vals >= vals.quantile(0.01)) & (vals <= vals.quantile(0.99))]
    results['Overall'] = f"{vals.mean():.{decimals}f} ± {vals.std():.{decimals}f}"
    return results

def describe_var_cat(df, col, decimals=1):
    results = {}
    for stage in [0, 1, 2, 3]:
        sub = df[df['kdigo_stage'] == stage]
        n = sub[col].sum()
        pct = n / len(sub) * 100 if len(sub) > 0 else 0
        results[f'Stage_{stage}'] = f"{int(n)} ({pct:.{decimals}f}%)"
    n = df[col].sum()
    pct = n / len(df) * 100
    results['Overall'] = f"{int(n)} ({pct:.{decimals}f}%)"
    return results

def row(label, s0, s1, s2, s3, ovr):
    return [label, s0, s1, s2, s3, ovr]

table1 = []

# Demographics
table1.append(row('Demographics', '', '', '', '', ''))
table1.append(row('N', str((valid['kdigo_stage']==0).sum()), str((valid['kdigo_stage']==1).sum()),
                  str((valid['kdigo_stage']==2).sum()), str((valid['kdigo_stage']==3).sum()), str(len(valid))))
d = describe_var_cont(valid, 'anchor_age')
table1.append(row('Age (years)', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))
d = describe_var_cat(valid, 'male')
table1.append(row('Male, n (%)', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))
d = describe_var_cat(valid, 'race_white')
table1.append(row('Race — White, n (%)', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))
d = describe_var_cat(valid, 'race_black')
table1.append(row('Race — Black, n (%)', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))

# Severity scores
table1.append(row('Severity Scores', '', '', '', '', ''))
d = describe_var_cont(valid, 'sofa_total', decimals=1)
table1.append(row('SOFA score', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))
d = describe_var_cont(valid, 'gcs_total', decimals=1)
table1.append(row('GCS', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))
d = describe_var_cat(valid, 'mech_vent')
table1.append(row('Mechanical ventilation, n (%)', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))
d = describe_var_cat(valid, 'vasopressor')
table1.append(row('Vasopressor use, n (%)', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))

# Renal
table1.append(row('Renal Function', '', '', '', '', ''))
d = describe_var_cont(valid, 'baseline_cr', decimals=2)
table1.append(row('Baseline Cr (mg/dL)', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))
d = describe_var_cont(valid, 'peak_cr', decimals=2)
table1.append(row('Peak Cr (mg/dL)', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))
d = describe_var_cont(valid, 'cr_ratio', decimals=2)
table1.append(row('Cr ratio (peak/baseline)', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))
d = describe_var_cat(valid, 'ckd_icd')
table1.append(row('CKD (ICD), n (%)', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))
d = describe_var_cat(valid, 'ckd_lab')
table1.append(row('CKD (lab, eGFR<60), n (%)', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))
d = describe_var_cont(valid, 'egfr_ckdepi', decimals=1)
table1.append(row('eGFR (CKD-EPI)', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))

# Comorbidities
table1.append(row('Comorbidities', '', '', '', '', ''))
comorb_cols = sorted([c for c in valid.columns if c.startswith('com_')])
for col in comorb_cols:
    name = col.replace('com_', '').replace('_', ' ')
    d = describe_var_cat(valid, col)
    table1.append(row(f'  {name}, n (%)', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))

# Vital Signs (filtered)
table1.append(row('Vital Signs (first 24h, filtered)', '', '', '', '', ''))
for vname in ['heart_rate', 'sbp', 'dbp', 'mbp', 'resp_rate', 'temperature', 'spo2']:
    col = f'{vname}_mean'
    if col in valid.columns:
        d = describe_var_cont(valid, col, decimals=1)
        table1.append(row(f'  {vname}', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))

# Labs
table1.append(row('Admission Labs', '', '', '', '', ''))
for lname in ['bun', 'sodium', 'potassium', 'bicarbonate', 'lactate', 'wbc', 'hemoglobin', 'platelet', 'bilirubin_total', 'albumin', 'glucose']:
    col = f'lab_{lname}'
    if col in valid.columns:
        d = describe_var_cont(valid, col, decimals=1)
        table1.append(row(f'  {lname}', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))

# Outcomes
table1.append(row('Outcomes', '', '', '', '', ''))
d = describe_var_cat(valid, 'hospital_expire_flag')
table1.append(row('Hospital mortality, n (%)', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))
d = describe_var_cont(valid, 'icu_los_hours', decimals=1)
table1.append(row('ICU LOS (hours)', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))
d = describe_var_cont(valid, 'hosp_los_days', decimals=1)
table1.append(row('Hospital LOS (days)', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))
d = describe_var_cat(valid, 'rrt')
table1.append(row('RRT, n (%)', d['Stage_0'], d['Stage_1'], d['Stage_2'], d['Stage_3'], d['Overall']))

# Save Table 1
df_t1 = pd.DataFrame(table1, columns=['Variable', 'No AKI', 'Stage 1', 'Stage 2', 'Stage 3', 'Overall'])
df_t1.to_csv(os.path.join(OUT, 'table1_revised.csv'), index=False)
print(f"  Table 1 saved: {len(df_t1)} rows")

# ============================================================
# PHASE 12: 缺失数据报告
# ============================================================
print("\n" + "=" * 60)
print("PHASE 12: Missing data report")
print("=" * 60)

key_vars = ['sofa_total', 'gcs_total', 'mech_vent', 'vasopressor',
            'heart_rate_mean', 'sbp_mean', 'dbp_mean', 'mbp_mean', 'resp_rate_mean', 'spo2_mean',
            'lab_sodium', 'lab_potassium', 'lab_bicarbonate', 'lab_lactate',
            'lab_wbc', 'lab_hemoglobin', 'lab_platelet', 'lab_bilirubin_total', 'lab_albumin', 'lab_glucose',
            'egfr_ckdepi', 'ckd_lab', 'mdrd_baseline_cr']

missing_report = []
for var in key_vars:
    if var in valid.columns:
        miss_n = valid[var].isna().sum()
        miss_pct = miss_n / len(valid) * 100
        missing_report.append({'variable': var, 'missing_n': int(miss_n), 'missing_pct': round(miss_pct, 1)})
        print(f"  {var}: {miss_n:,} missing ({miss_pct:.1f}%)")

missing_df = pd.DataFrame(missing_report)
missing_df.to_csv(os.path.join(OUT, 'missing_data_report.csv'), index=False)

# ============================================================
# PHASE 13: 流程图数据 + Final save
# ============================================================
print("\n" + "=" * 60)
print("PHASE 13: Flowchart data + Final save")
print("=" * 60)

flowchart = {
    'total_icu_stays': int(n_total_icu),
    'total_hadm': int(n_total_hadm),
    'total_patients': int(n_total_patients),
    'excluded_age_lt_18': int(n_excluded_age),
    'excluded_duplicate_stays': int(n_excluded_duplicate),
    'after_first_icu_adult': int(n_adult_stays),
    'excluded_missing_cr': int(n_excluded_no_cr),
    'final_cohort': int(len(valid)),
    'aki_n': int((valid['kdigo_stage'] >= 1).sum()),
    'aki_pct': round((valid['kdigo_stage'] >= 1).sum() / len(valid) * 100, 1),
    'stage_counts': {
        '0': int((valid['kdigo_stage'] == 0).sum()),
        '1': int((valid['kdigo_stage'] == 1).sum()),
        '2': int((valid['kdigo_stage'] == 2).sum()),
        '3': int((valid['kdigo_stage'] == 3).sum()),
    },
    'stage3_stratification': {
        'rrt_only': int((valid['stage3_by_rrt_only'] == 1).sum()),
        'cr_only': int((valid['stage3_by_cr_only'] == 1).sum()),
        'both': int((valid['stage3_by_both'] == 1).sum()),
    },
    'vital_outlier_counts': vital_outlier_counts,
}

with open(os.path.join(OUT, 'flowchart_data.json'), 'w') as f:
    json.dump(flowchart, f, indent=2, default=str)

# Save enriched cohort for stats
valid.to_csv(os.path.join(OUT, "kdigo_cohort_revised_enriched.csv"), index=False)

# Final summary
print("\n" + "=" * 80)
print("DATA EXTRACTION COMPLETE")
print("=" * 80)
print(f"\nFinal cohort: {len(valid):,} ICU stays")
print(f"AKI prevalence: {flowchart['aki_n']:,} ({flowchart['aki_pct']}%)")
print(f"\nKDIGO distribution:")
for s in [0, 1, 2, 3]:
    n = flowchart['stage_counts'][str(s)]
    mort = valid[valid['kdigo_stage'] == s]['hospital_expire_flag'].mean() * 100
    print(f"  Stage {s}: {n:,} ({n/len(valid)*100:.1f}%), mortality={mort:.2f}%")

print(f"\nFiles saved:")
print(f"  - kdigo_cohort_revised.csv")
print(f"  - kdigo_cohort_revised_enriched.csv")
print(f"  - table1_revised.csv")
print(f"  - flowchart_data.json")
print(f"  - missing_data_report.csv")
print(f"\nEnd: {datetime.now()}")
