#!/usr/bin/env python3
"""
eICU-CRD External Validation of KDIGO AKI Staging → ICU Mortality
================================================================
Replicates the KDIGO AKI staging analysis from MIMIC-IV v3.1 (N=84,167)
in the eICU-CRD database (N=~200,000) as an independent validation cohort.

KDIGO staging: creatinine-based only (consistent with MIMIC-IV paper)
- Stage 1: 1.5-1.9x baseline OR >=0.3 mg/dL increase in 48h
- Stage 2: 2.0-2.9x baseline
- Stage 3: >=3.0x baseline OR >=4.0 mg/dL OR RRT

Author: Dengke Wu
Date: 2026-08-08
"""

import duckdb
import pandas as pd
import numpy as np
import json
import os
import warnings
warnings.filterwarnings('ignore')

# Paths
EICU_PATH = 'E:/mimic-iv/eicu/physionet.org/files/eicu-crd/2.0'
OUTPUT_DIR = 'C:/Users/admin/WorkBuddy/2026-07-06-13-22-19'

print("=" * 70)
print("eICU-CRD External Validation of KDIGO AKI Staging")
print("=" * 70)

# ============================================================
# Phase 1: Cohort Selection
# ============================================================
print("\n[Phase 1] Cohort Selection...")

# Get patient table with key fields
patient_query = f"""
SELECT
    patientunitstayid,
    gender,
    age,
    ethnicity,
    hospitalid,
    unittype,
    unitdischargeoffset,
    hospitaldischargeoffset,
    unitdischargestatus,
    hospitaldischargestatus,
    admissionweight,
    unitstaytype
FROM read_csv_auto('{EICU_PATH}/patient.csv.gz', sample_size=100000)
WHERE age IS NOT NULL
  AND age != ''
  AND age != '> 89'
"""

# Handle age > 89 (encoded as '> 89' in eICU)
patient_df = duckdb.query(patient_query).df()
# Convert age to numeric, set '> 89' to 90
patient_df['age_num'] = pd.to_numeric(patient_df['age'], errors='coerce')
patient_df.loc[patient_df['age_num'].isna() & patient_df['age'].astype(str).str.contains('> 89', na=False), 'age_num'] = 90

# Filter: adult (>=18), has valid age
cohort = patient_df[patient_df['age_num'] >= 18].copy()
print(f"  Total eICU patients: {len(patient_df):,}")
print(f"  Adult (>=18): {len(cohort):,}")

# Mortality
cohort['mortality'] = (cohort['hospitaldischargestatus'] == 'Expired').astype(int)
print(f"  Hospital mortality: {cohort['mortality'].sum():,} ({cohort['mortality'].mean()*100:.1f}%)")

# Unit discharge offset is in minutes; convert to days for ICU LOS
cohort['icu_los_days'] = pd.to_numeric(cohort['unitdischargeoffset'], errors='coerce') / 1440.0
cohort = cohort[cohort['icu_los_days'] >= 0.5].copy()  # At least 12h stay (consistent with MIMIC)
print(f"  After >=12h ICU stay filter: {len(cohort):,}")

# ============================================================
# Phase 2: Creatinine Data Extraction & KDIGO Staging
# ============================================================
print("\n[Phase 2] Creatinine Extraction & KDIGO Staging...")

# Extract all creatinine values with offsets (in minutes from admission)
cr_query = f"""
SELECT
    patientunitstayid,
    labresultoffset,
    labresult
FROM read_csv_auto('{EICU_PATH}/lab.csv.gz', sample_size=1000000)
WHERE labname = 'creatinine'
  AND labresult IS NOT NULL
  AND labresult > 0
  AND labresult < 50
"""
cr_df = duckdb.query(cr_query).df()
print(f"  Total creatinine measurements: {len(cr_df):,}")
print(f"  Patients with creatinine data: {cr_df['patientunitstayid'].nunique():,}")

# Filter to cohort patients
cr_df = cr_df[cr_df['patientunitstayid'].isin(cohort['patientunitstayid'])]
print(f"  In cohort: {len(cr_df):,} measurements, {cr_df['patientunitstayid'].nunique():,} patients")

# Sort by patient and offset
cr_df = cr_df.sort_values(['patientunitstayid', 'labresultoffset'])

# Calculate baseline creatinine per patient
# Strategy: minimum SCr within first 24h (1440 min) of ICU stay
# If no SCr in first 24h, use first available SCr
# If no SCr at all, use MDRD back-calculation (assume eGFR=75)

def calculate_baseline_cr(group):
    """Calculate baseline creatinine for a patient."""
    vals_24h = group[group['labresultoffset'] <= 1440]
    if len(vals_24h) > 0:
        return vals_24h['labresult'].min()
    return group['labresult'].iloc[0]  # First available

baseline_cr = cr_df.groupby('patientunitstayid').apply(calculate_baseline_cr).reset_index()
baseline_cr.columns = ['patientunitstayid', 'baseline_cr']

# MDRD back-calculation for patients without any SCr
# eGFR = 175 * SCr^(-1.154) * age^(-0.203) * (0.742 if female)
# SCr = (175 * age^(-0.203) * gender_factor / 75)^(1/1.154)
no_cr_patients = cohort[~cohort['patientunitstayid'].isin(cr_df['patientunitstayid'])].copy()
print(f"  Patients without SCr (will use MDRD): {len(no_cr_patients):,}")

if len(no_cr_patients) > 0:
    # MDRD back-calculation
    no_cr_patients['gender_factor'] = np.where(no_cr_patients['gender'] == 'Female', 0.742, 1.0)
    no_cr_patients['baseline_cr_mdrd'] = (
        175 * (no_cr_patients['age_num'] ** (-0.203)) * no_cr_patients['gender_factor'] / 75
    ) ** (1 / 1.154)
    # Round to typical SCr values
    no_cr_patients['baseline_cr'] = no_cr_patients['baseline_cr_mdrd'].clip(0.4, 3.0)

    # For these patients, KDIGO stage = 0 (no SCr measurements to stage)
    no_cr_patients['kdigo_stage'] = 0
    no_cr_patients['cr_data_source'] = 'mdrd'

# For patients with SCr data, apply KDIGO staging
print("\n  Applying KDIGO creatinine criteria...")

def apply_kdigo_staging(group, baseline):
    """Apply KDIGO creatinine-based staging."""
    if baseline is None or pd.isna(baseline) or baseline <= 0:
        return 0  # Cannot stage

    max_cr = group['labresult'].max()

    # Check for >=0.3 mg/dL increase in 48h window
    # Need to check consecutive measurements within 48h
    group_sorted = group.sort_values('labresultoffset')
    cr_values = group_sorted['labresult'].values
    offsets = group_sorted['labresultoffset'].values

    stage_1_48h = False
    for i in range(len(cr_values)):
        for j in range(i + 1, len(cr_values)):
            if offsets[j] - offsets[i] <= 2880:  # 48h = 2880 min
                if cr_values[j] - cr_values[i] >= 0.3:
                    stage_1_48h = True
                    break
        if stage_1_48h:
            break

    ratio = max_cr / baseline

    # Stage 3: >=3x baseline OR >=4.0 mg/dL OR RRT
    if ratio >= 3.0 or max_cr >= 4.0:
        return 3
    # Stage 2: 2-2.9x baseline
    elif ratio >= 2.0:
        return 2
    # Stage 1: 1.5-1.9x baseline OR >=0.3 increase in 48h
    elif ratio >= 1.5 or stage_1_48h:
        return 1
    else:
        return 0

# Apply KDIGO staging (this may take a minute for large datasets)
print("  Processing KDIGO staging (may take 1-2 minutes)...")
kdigo_results = []
for pid, group in cr_df.groupby('patientunitstayid'):
    base_val = baseline_cr.loc[baseline_cr['patientunitstayid'] == pid, 'baseline_cr'].values
    if len(base_val) > 0:
        stage = apply_kdigo_staging(group, base_val[0])
    else:
        stage = 0
    kdigo_results.append({'patientunitstayid': pid, 'kdigo_stage': stage, 'cr_data_source': 'measured'})

kdigo_df = pd.DataFrame(kdigo_results)

# Merge with MDRD-based patients
if len(no_cr_patients) > 0:
    mdrd_stages = no_cr_patients[['patientunitstayid', 'kdigo_stage', 'cr_data_source']]
    kdigo_df = pd.concat([kdigo_df, mdrd_stages], ignore_index=True)

print(f"\n  KDIGO staging complete:")
print(f"  Stage 0: {(kdigo_df['kdigo_stage'] == 0).sum():,}")
print(f"  Stage 1: {(kdigo_df['kdigo_stage'] == 1).sum():,}")
print(f"  Stage 2: {(kdigo_df['kdigo_stage'] == 2).sum():,}")
print(f"  Stage 3: {(kdigo_df['kdigo_stage'] == 3).sum():,}")

# ============================================================
# Phase 3: RRT Data
# ============================================================
print("\n[Phase 3] RRT Extraction...")

rrt_query = f"""
SELECT DISTINCT patientunitstayid
FROM read_csv_auto('{EICU_PATH}/treatment.csv.gz', sample_size=100000)
WHERE treatmentstring ILIKE '%dialysis%'
  OR treatmentstring ILIKE '%CRRT%'
  OR treatmentstring ILIKE '%hemodialysis%'
  OR treatmentstring ILIKE '%ultrafiltration%'
"""
rrt_patients = duckdb.query(rrt_query).df()
rrt_set = set(rrt_patients['patientunitstayid'].values)
print(f"  Patients with RRT: {len(rrt_set):,}")

# Override KDIGO stage to 3 for RRT patients
kdigo_df.loc[kdigo_df['patientunitstayid'].isin(rrt_set), 'kdigo_stage'] = 3
print(f"  KDIGO Stage 3 (after RRT override): {(kdigo_df['kdigo_stage'] == 3).sum():,}")

# ============================================================
# Phase 4: Severity Adjustment Variables
# ============================================================
print("\n[Phase 4] Severity Adjustment Variables...")

# APACHE score (equivalent to SOFA in MIMIC)
apache_query = f"""
SELECT
    patientunitstayid,
    apachescore,
    apacheversion
FROM read_csv_auto('{EICU_PATH}/apachePatientResult.csv.gz', sample_size=100000)
WHERE apacheversion = 'IVa'
"""
apache_df = duckdb.query(apache_query).df()
# Keep one row per patient (take max if duplicates)
apache_df = apache_df.groupby('patientunitstayid').first().reset_index()
print(f"  Patients with APACHE IVa: {len(apache_df):,}")

# Mechanical ventilation (ventday1 from apachePredVar)
vent_query = f"""
SELECT
    patientunitstayid,
    ventday1
FROM read_csv_auto('{EICU_PATH}/apachePredVar.csv.gz', sample_size=100000)
"""
vent_df = duckdb.query(vent_query).df()
vent_df = vent_df.groupby('patientunitstayid').first().reset_index()
print(f"  Patients with vent data: {len(vent_df):,}")

# CKD from diagnosis
ckd_query = f"""
SELECT DISTINCT patientunitstayid
FROM read_csv_auto('{EICU_PATH}/diagnosis.csv.gz', sample_size=100000)
WHERE diagnosisstring ILIKE '%chronic kidney disease%'
"""
ckd_patients = duckdb.query(ckd_query).df()
ckd_set = set(ckd_patients['patientunitstayid'].values)
print(f"  Patients with CKD diagnosis: {len(ckd_set):,}")

# ============================================================
# Phase 5: Merge All Data
# ============================================================
print("\n[Phase 5] Merging All Data...")

# Start with cohort
analysis_df = cohort.merge(kdigo_df, on='patientunitstayid', how='inner')
analysis_df = analysis_df.merge(apache_df[['patientunitstayid', 'apachescore']], on='patientunitstayid', how='left')
analysis_df = analysis_df.merge(vent_df[['patientunitstayid', 'ventday1']], on='patientunitstayid', how='left')
analysis_df['ckd'] = analysis_df['patientunitstayid'].isin(ckd_set).astype(int)
analysis_df['rrt'] = analysis_df['patientunitstayid'].isin(rrt_set).astype(int)

# Create mechanical ventilation binary variable
analysis_df['mech_vent'] = (analysis_df['ventday1'] > 0).astype(int)
analysis_df.loc[analysis_df['ventday1'].isna(), 'mech_vent'] = np.nan

# Impute missing APACHE with median
apache_median = analysis_df['apachescore'].median()
analysis_df['apachescore'] = analysis_df['apachescore'].fillna(apache_median)
analysis_df['apache_imputed'] = (analysis_df['apachescore'] == apache_median).astype(int)

# Impute missing mech_vent with 0 (most likely not ventilated)
analysis_df['mech_vent'] = analysis_df['mech_vent'].fillna(0)

print(f"\n  Final analysis cohort: {len(analysis_df):,}")
print(f"  Mortality: {analysis_df['mortality'].sum():,} ({analysis_df['mortality'].mean()*100:.1f}%)")

# Mortality by KDIGO stage
print("\n  Mortality by KDIGO Stage:")
for stage in [0, 1, 2, 3]:
    subset = analysis_df[analysis_df['kdigo_stage'] == stage]
    mort = subset['mortality'].mean() * 100 if len(subset) > 0 else 0
    print(f"    Stage {stage}: n={len(subset):,}, mortality={mort:.1f}%")

# ============================================================
# Phase 6: Statistical Analysis
# ============================================================
print("\n[Phase 6] Statistical Analysis...")
import statsmodels.api as sm
from lifelines import CoxPHFitter

# Prepare variables
analysis_df['kdigo_1'] = (analysis_df['kdigo_stage'] == 1).astype(int)
analysis_df['kdigo_2'] = (analysis_df['kdigo_stage'] == 2).astype(int)
analysis_df['kdigo_3'] = (analysis_df['kdigo_stage'] == 3).astype(int)

analysis_df['gender_male'] = (analysis_df['gender'] == 'Male').astype(int)

# --- Model A: Unadjusted ---
print("\n  [Model A] Unadjusted Logistic Regression:")
X_a = analysis_df[['kdigo_1', 'kdigo_2', 'kdigo_3']].astype(float)
X_a = sm.add_constant(X_a)
y = analysis_df['mortality'].astype(float)
model_a = sm.Logit(y, X_a).fit(disp=0)

results_a = {}
for var in ['kdigo_1', 'kdigo_2', 'kdigo_3']:
    or_val = np.exp(model_a.params[var])
    ci_low = np.exp(model_a.conf_int().loc[var, 0])
    ci_high = np.exp(model_a.conf_int().loc[var, 1])
    p = model_a.pvalues[var]
    results_a[var] = {'OR': round(or_val, 3), 'CI_low': round(ci_low, 3), 'CI_high': round(ci_high, 3), 'p': p}
    print(f"    {var}: OR={or_val:.3f} ({ci_low:.3f}-{ci_high:.3f}), p={p:.4f}")

# --- Model B: Adjusted (APACHE + mech_vent + age + gender) ---
print("\n  [Model B] Adjusted Logistic Regression:")
X_b = analysis_df[['kdigo_1', 'kdigo_2', 'kdigo_3', 'apachescore', 'mech_vent', 'age_num', 'gender_male']].astype(float)
X_b = sm.add_constant(X_b)
model_b = sm.Logit(y, X_b).fit(disp=0)

results_b = {}
for var in ['kdigo_1', 'kdigo_2', 'kdigo_3']:
    or_val = np.exp(model_b.params[var])
    ci_low = np.exp(model_b.conf_int().loc[var, 0])
    ci_high = np.exp(model_b.conf_int().loc[var, 1])
    p = model_b.pvalues[var]
    results_b[var] = {'OR': round(or_val, 3), 'CI_low': round(ci_low, 3), 'CI_high': round(ci_high, 3), 'p': p}
    print(f"    {var}: OR={or_val:.3f} ({ci_low:.3f}-{ci_high:.3f}), p={p:.4f}")

# AUC comparison
from sklearn.metrics import roc_auc_score
auc_a = roc_auc_score(y, model_a.predict(X_a))
auc_b = roc_auc_score(y, model_b.predict(X_b))
print(f"\n  AUC Model A: {auc_a:.3f}")
print(f"  AUC Model B: {auc_b:.3f}")

# --- Cox Regression ---
print("\n  [Cox Regression] Model B (adjusted):")
# Use unitdischargeoffset as time, mortality as event
# Remove patients with invalid time
cox_df = analysis_df[analysis_df['icu_los_days'] > 0].copy()
cox_df['duration'] = cox_df['icu_los_days']
cox_df['event'] = cox_df['mortality']

# Cap duration at 30 days (like MIMIC)
cox_df['duration'] = cox_df['duration'].clip(upper=30)

cox_vars = ['kdigo_1', 'kdigo_2', 'kdigo_3', 'apachescore', 'mech_vent', 'age_num', 'gender_male']
cox_data = cox_df[['duration', 'event'] + cox_vars].dropna().astype(float)

cph = CoxPHFitter()
cph.fit(cox_data, duration_col='duration', event_col='event')

results_cox = {}
for var in ['kdigo_1', 'kdigo_2', 'kdigo_3']:
    hr_val = np.exp(cph.params_[var])
    ci_low = np.exp(cph.confidence_intervals_.loc[var, '95% lower-bound'])
    ci_high = np.exp(cph.confidence_intervals_.loc[var, '95% upper-bound'])
    p = cph.summary.loc[var, 'p']
    results_cox[var] = {'HR': round(hr_val, 3), 'CI_low': round(ci_low, 3), 'CI_high': round(ci_high, 3), 'p': p}
    print(f"    {var}: HR={hr_val:.3f} ({ci_low:.3f}-{ci_high:.3f}), p={p:.4f}")

# ============================================================
# Phase 7: CKD × AKI Interaction (AKI Paradox)
# ============================================================
print("\n[Phase 7] CKD × AKI Interaction Test:")

analysis_df['ckd_kdigo1'] = analysis_df['ckd'] * analysis_df['kdigo_1']
analysis_df['ckd_kdigo2'] = analysis_df['ckd'] * analysis_df['kdigo_2']
analysis_df['ckd_kdigo3'] = analysis_df['ckd'] * analysis_df['kdigo_3']

X_int = analysis_df[['kdigo_1', 'kdigo_2', 'kdigo_3', 'ckd',
                      'ckd_kdigo1', 'ckd_kdigo2', 'ckd_kdigo3',
                      'apachescore', 'mech_vent', 'age_num', 'gender_male']].astype(float)
X_int = sm.add_constant(X_int)
model_int = sm.Logit(y, X_int).fit(disp=0)

print("  Interaction ORs:")
interaction_results = {}
for var in ['ckd_kdigo1', 'ckd_kdigo2', 'ckd_kdigo3']:
    or_val = np.exp(model_int.params[var])
    ci_low = np.exp(model_int.conf_int().loc[var, 0])
    ci_high = np.exp(model_int.conf_int().loc[var, 1])
    p = model_int.pvalues[var]
    interaction_results[var] = {'OR': round(or_val, 3), 'CI_low': round(ci_low, 3), 'CI_high': round(ci_high, 3), 'p': p}
    print(f"    {var}: OR={or_val:.3f} ({ci_low:.3f}-{ci_high:.3f}), p={p:.4f}")

# ============================================================
# Phase 8: Save Results
# ============================================================
print("\n[Phase 8] Saving Results...")

# Summary
summary = {
    'database': 'eICU-CRD v2.0',
    'n_total': int(len(analysis_df)),
    'kdigo_distribution': {
        '0': int((analysis_df['kdigo_stage'] == 0).sum()),
        '1': int((analysis_df['kdigo_stage'] == 1).sum()),
        '2': int((analysis_df['kdigo_stage'] == 2).sum()),
        '3': int((analysis_df['kdigo_stage'] == 3).sum()),
    },
    'aki_prevalence_pct': round(analysis_df['kdigo_stage'].gt(0).mean() * 100, 1),
    'mortality_by_stage': {
        '0': round(analysis_df[analysis_df['kdigo_stage'] == 0]['mortality'].mean() * 100, 1),
        '1': round(analysis_df[analysis_df['kdigo_stage'] == 1]['mortality'].mean() * 100, 1),
        '2': round(analysis_df[analysis_df['kdigo_stage'] == 2]['mortality'].mean() * 100, 1),
        '3': round(analysis_df[analysis_df['kdigo_stage'] == 3]['mortality'].mean() * 100, 1),
    },
    'overall_mortality_pct': round(analysis_df['mortality'].mean() * 100, 1),
    'mean_age': round(analysis_df['age_num'].mean(), 1),
    'male_pct': round((analysis_df['gender'] == 'Male').mean() * 100, 1),
    'mean_apache': round(analysis_df['apachescore'].mean(), 1),
    'mech_vent_pct': round(analysis_df['mech_vent'].mean() * 100, 1),
    'ckd_pct': round(analysis_df['ckd'].mean() * 100, 1),
    'rrt_pct': round(analysis_df['rrt'].mean() * 100, 1),
    'kdigo_or_model_a': {k: v for k, v in results_a.items()},
    'kdigo_or_model_b': {k: v for k, v in results_b.items()},
    'kdigo_hr_model_b': {k: v for k, v in results_cox.items()},
    'ckd_interaction': interaction_results,
    'auc_model_a': round(auc_a, 3),
    'auc_model_b': round(auc_b, 3),
}

# Save to JSON
output_path = os.path.join(OUTPUT_DIR, 'eicu_validation_results.json')
with open(output_path, 'w', encoding='utf-8') as f:
    json.dump(summary, f, indent=2, ensure_ascii=False, default=str)
print(f"\n  Results saved to: {output_path}")

# Save analysis data for further use
analysis_df.to_csv(os.path.join(OUTPUT_DIR, 'eicu_validation_data.csv'), index=False)
print(f"  Analysis data saved to: eicu_validation_data.csv")

# ============================================================
# Phase 9: Comparison with MIMIC-IV
# ============================================================
print("\n" + "=" * 70)
print("COMPARISON: MIMIC-IV (Discovery) vs eICU-CRD (Validation)")
print("=" * 70)

# Load MIMIC results
mimic_path = os.path.join(OUTPUT_DIR, 'summary_revised.json')
if os.path.exists(mimic_path):
    with open(mimic_path, 'r', encoding='utf-8') as f:
        mimic = json.load(f)

    print(f"\n{'Metric':<25} {'MIMIC-IV':>15} {'eICU-CRD':>15}")
    print("-" * 55)
    print(f"{'N':<25} {mimic['n_total']:>15,} {len(analysis_df):>15,}")
    print(f"{'AKI prevalence (%)':<25} {mimic['aki_prevalence_pct']:>15.1f} {summary['aki_prevalence_pct']:>15.1f}")

    for stage in ['0', '1', '2', '3']:
        m_mort = mimic['mortality_by_stage'][stage]
        e_mort = summary['mortality_by_stage'][stage]
        print(f"{'Stage ' + stage + ' mortality (%)':<25} {m_mort:>15.1f} {e_mort:>15.1f}")

    print(f"\n{'OR Model B':<25} {'MIMIC-IV':>15} {'eICU-CRD':>15}")
    print("-" * 55)
    for stage in ['1', '2', '3']:
        key = f'kdigo_{stage}'
        m_or = mimic['kdigo_or_model_b'][stage]['OR']
        e_or = results_b.get(key, {}).get('OR', 'N/A')
        m_str = f"{m_or:.3f}"
        e_str = f"{e_or:.3f}" if isinstance(e_or, float) else str(e_or)
        print(f"{'  Stage ' + stage + ' OR':<25} {m_str:>15} {e_str:>15}")

    print(f"\n{'HR Model B':<25} {'MIMIC-IV':>15} {'eICU-CRD':>15}")
    print("-" * 55)
    for stage in ['1', '2', '3']:
        key = f'kdigo_{stage}'
        m_hr = mimic['kdigo_hr_model_b'][stage]['HR']
        e_hr = results_cox.get(key, {}).get('HR', 'N/A')
        m_str = f"{m_hr:.3f}"
        e_str = f"{e_hr:.3f}" if isinstance(e_hr, float) else str(e_hr)
        print(f"{'  Stage ' + stage + ' HR':<25} {m_str:>15} {e_str:>15}")

    print(f"\n{'AUC Model B':<25} {'0.828':>15} {auc_b:>15.3f}")

    print(f"\n{'CKD×Stage3 Interaction':<25} {'MIMIC-IV':>15} {'eICU-CRD':>15}")
    print("-" * 55)
    m_int = mimic.get('ckd_interaction', {}).get('stage3_interaction', {}).get('OR', 'N/A')
    e_int = interaction_results.get('ckd_kdigo3', {}).get('OR', 'N/A')
    print(f"{'  OR':<25} {m_int:>15} {e_int:>15}")

print("\n" + "=" * 70)
print("eICU-CRD External Validation COMPLETE")
print("=" * 70)
