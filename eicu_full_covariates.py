"""
eICU-CRD v2 — 完整协变量集构建 (F7 真外部验证前置条件)
================================================================
对照外部评审报告 (2026-09-23) F7 修复：

  MIMIC Model C (mimic_analysis_v2.csv) 协变量：
    kdigo_stage, age, sex_male, sofa_nonrenal, sofa_full, sofa_n_missing,
    vent_24h, ckd_history_prior, diabetes, hypertension, heart_failure,
    copd, liver_disease, event_time_h, event_ind

  当前 eICU v2 CSV (eicu_cohort_v2.csv) 仅含：kdigo_stage, age_num, gender,
  apachescore, mech_vent_day1, ckd_history, hospital_mortality, hospital_los_h.
  缺：sofa_nonrenal/full, sofa_n_missing, vent_24h (==mech_vent_day1),
      sex_male, diabetes, hypertension, heart_failure, copd, liver_disease,
      event_time_h, event_ind

本脚本：
  1. 把 mech_vent_day1 → vent_24h（统一命名）
  2. 把 gender → sex_male (统一命名)
  3. 5 大合并症：pastHistory 表（与 MIMIC ICD 等价）
  4. SOFA 6 分量：mimic-code 风格移植到 eICU (nurseCharting + lab + medication
     + respiratoryCare + intakeOutput)
  5. event_time_h = hospitaldischargeoffset / 60 (小时)；event_ind =
     hospital_mortality

时间窗 (与 MIMIC sofa_first24h.py 一致):
  - MAP / GCS / FiO2: [intime-6h, intime+24h] (nurseCharting offset, in minutes)
  - Lab: [intime-6h, intime+24h] (laboffset)
  - Medication: [intime-6h, intime+24h] (drugstartoffset)
  - RespiratoryCare vent: [intime-6h, intime+24h] (respcarestatusoffset)
  - Urine: [intime-6h, intime+24h] (intakeoutputoffset)

输出: eicu_analysis_v2.csv + eicu_analysis_v2_meta.json
"""
import duckdb
import pandas as pd
import numpy as np
import os, json
from datetime import datetime

print("=" * 80)
print("eICU-CRD v2 — Full covariate set construction (F7 prerequisite)")
print(f"Start: {datetime.now()}")
print("=" * 80)

EICU = "E:/mimic-iv/eicu/physionet.org/files/eicu-crd/2.0"
OUT = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
CSV_OPTS = ", quote='\"', ignore_errors=true"

def q(sql):
    return duckdb.query(sql).df()

flow = {}

# ============================================================
# PHASE 1: Load base cohort (eICU v2)
# ============================================================
print("\n[PHASE 1] Load base cohort eicu_cohort_v2.csv")
cohort = pd.read_csv(os.path.join(OUT, "eicu_cohort_v2.csv"))
print(f"  Cohort: {len(cohort):,} unique patient stays")
flow['cohort_n'] = int(len(cohort))

# Need: patientunitstayid, age, gender, kdigo_stage, aki_onset_h, hospital_los_h, hospital_mortality
stay_ids = cohort['patientunitstayid'].astype('int64').tolist()
print(f"  Stays to query: {len(stay_ids):,}")

# Register as DuckDB view for fast joins
cohort_df = pd.DataFrame({
    'stay_id': stay_ids,
    'age':     cohort['age_num'].values,
    'kdigo':   cohort['kdigo_stage'].values,
    'aki_h':   cohort['aki_onset_h'].fillna(-1).values,
    'hosp_los_h': cohort['hospital_los_h'].fillna(-1).values,
    'hosp_mort':   cohort['hospital_mortality'].fillna(0).values,
})
duckdb.register('cohort_stays', cohort_df)

# =============================================================================
# PHASE 2: SOFA first 24h — mimic-code style (6 components)
# =============================================================================
print("\n[PHASE 2] SOFA first 24h (6 components) — mimic-code style")

# We use a unified AKI-onset + 24h window:
#   onset_h is "hours from unit admit to KDIGO AKI onset" (NULL for stage 0).
#   For stage 0, set onset to 0 (use full first 24h).
#   For stage >=1, set window end = min(onset_h + 24, 24) so the SOFA window
#   mirrors "first 24h after AKI onset" — this aligns with MIMIC's first-24h-of-ICU
#   convention for stage 0 patients (== first 24h of stay).
# Actually to keep consistent: ALWAYS use [intime-6h, intime+24h] for SOFA scoring
# (mimic-code convention). This ensures the same biological window for all stays.
# AKI-onset is used ONLY for KDIGO stage assignment, not for SOFA time window.

# Helper: get stays with offset window [intime-6h, intime+24h]
# All eICU tables use offset relative to unit admission in minutes.

# 2A. MAP (4 priority sources — take min over window)
print("  [SOFA-CV] MAP extraction ...")
map_agg = q(f"""
    SELECT nc.patientunitstayid AS stay_id,
           MIN(TRY_CAST(nc.nursingchartvalue AS DOUBLE)) AS min_map
    FROM read_csv('{EICU}/nurseCharting.csv.gz'{CSV_OPTS}) nc
    JOIN cohort_stays cs ON nc.patientunitstayid = cs.stay_id
    WHERE nc.nursingchartcelltypevallabel IN (
        'MAP (mmHg)',
        'Arterial Line MAP (mmHg)',
        'Invasive BP',
        'Non-Invasive BP'
    )
      AND nc.nursingchartcelltypevalname IN (
        'Value',
        'Invasive BP Mean',
        'Non-Invasive BP Mean'
      )
      AND TRY_CAST(nc.nursingchartvalue AS DOUBLE) BETWEEN 20 AND 250
      AND nc.nursingchartoffset BETWEEN -360 AND 1440
    GROUP BY nc.patientunitstayid
""")
print(f"    MAP rows kept: {len(map_agg):,} of {len(cohort):,} stays")

# 2B. Vasoactive drugs in first 24h (binary for SOFA cardiovascular score)
print("  [SOFA-CV] Vasopressor detection ...")
vaso = q(f"""
    SELECT m.patientunitstayid AS stay_id,
           MAX(CASE WHEN LOWER(m.drugname) LIKE '%norepinephrine%' THEN 1 ELSE 0 END) AS any_norepi,
           MAX(CASE WHEN LOWER(m.drugname) LIKE '%epinephrine%' AND LOWER(m.drugname) NOT LIKE '%norepinephrine%' THEN 1 ELSE 0 END) AS any_epi,
           MAX(CASE WHEN LOWER(m.drugname) LIKE '%dopamine%' AND LOWER(m.drugname) NOT LIKE '%epinephrine%' AND LOWER(m.drugname) NOT LIKE '%norepinephrine%' THEN 1 ELSE 0 END) AS any_dopa,
           MAX(CASE WHEN LOWER(m.drugname) LIKE '%dobutamine%' THEN 1 ELSE 0 END) AS any_dobu,
           MAX(CASE WHEN LOWER(m.drugname) LIKE '%vasopressin%' THEN 1 ELSE 0 END) AS any_vaso,
           MAX(CASE WHEN LOWER(m.drugname) LIKE '%milrinone%' THEN 1 ELSE 0 END) AS any_milri
    FROM read_csv('{EICU}/medication.csv.gz'{CSV_OPTS}) m
    JOIN cohort_stays cs ON m.patientunitstayid = cs.stay_id
    WHERE m.drugordercancelled = FALSE
      AND m.drugstartoffset BETWEEN -360 AND 1440
    GROUP BY m.patientunitstayid
""")
print(f"    Vasopressor stays: {len(vaso):,}")

# 2C. GCS Total (priority source)
print("  [SOFA-CNS] GCS extraction ...")
gcs = q(f"""
    SELECT nc.patientunitstayid AS stay_id,
           MIN(TRY_CAST(nc.nursingchartvalue AS DOUBLE)) AS min_gcs_total
    FROM read_csv('{EICU}/nurseCharting.csv.gz'{CSV_OPTS}) nc
    JOIN cohort_stays cs ON nc.patientunitstayid = cs.stay_id
    WHERE nc.nursingchartcelltypevallabel IN (
        'Glasgow coma score', 'Score (Glasgow Coma Scale)'
      )
      AND nc.nursingchartcelltypevalname IN ('GCS Total', 'Value')
      AND TRY_CAST(nc.nursingchartvalue AS DOUBLE) BETWEEN 3 AND 15
      AND nc.nursingchartoffset BETWEEN -360 AND 1440
    GROUP BY nc.patientunitstayid
""")
print(f"    GCS stays: {len(gcs):,}")

# 2D. Platelets (min over window)
print("  [SOFA-Coag] Platelets extraction ...")
plt = q(f"""
    SELECT l.patientunitstayid AS stay_id,
           MIN(TRY_CAST(l.labresult AS DOUBLE)) AS min_plt
    FROM read_csv('{EICU}/lab.csv.gz'{CSV_OPTS}) l
    JOIN cohort_stays cs ON l.patientunitstayid = cs.stay_id
    WHERE LOWER(l.labname) = 'platelets x 1000'
      AND TRY_CAST(l.labresult AS DOUBLE) BETWEEN 1 AND 800
      AND l.labresultoffset BETWEEN -360 AND 1440
    GROUP BY l.patientunitstayid
""")
print(f"    Platelets stays: {len(plt):,}")

# 2E. Bilirubin (max over window)
print("  [SOFA-Liver] Bilirubin extraction ...")
bili = q(f"""
    SELECT l.patientunitstayid AS stay_id,
           MAX(TRY_CAST(l.labresult AS DOUBLE)) AS max_bili
    FROM read_csv('{EICU}/lab.csv.gz'{CSV_OPTS}) l
    JOIN cohort_stays cs ON l.patientunitstayid = cs.stay_id
    WHERE LOWER(l.labname) = 'total bilirubin'
      AND TRY_CAST(l.labresult AS DOUBLE) BETWEEN 0.1 AND 60
      AND l.labresultoffset BETWEEN -360 AND 1440
    GROUP BY l.patientunitstayid
""")
print(f"    Bilirubin stays: {len(bili):,}")

# 2F. PaO2 (min over window) + FiO2 (max over window)
print("  [SOFA-Resp] PaO2 + FiO2 extraction ...")
pafi = q(f"""
    SELECT l.patientunitstayid AS stay_id,
           MIN(CASE WHEN LOWER(l.labname) = 'pao2'
                    THEN TRY_CAST(l.labresult AS DOUBLE) END) AS min_pao2
    FROM read_csv('{EICU}/lab.csv.gz'{CSV_OPTS}) l
    JOIN cohort_stays cs ON l.patientunitstayid = cs.stay_id
    WHERE LOWER(l.labname) = 'pao2'
      AND TRY_CAST(l.labresult AS DOUBLE) BETWEEN 20 AND 600
      AND l.labresultoffset BETWEEN -360 AND 1440
    GROUP BY l.patientunitstayid
""")
fio2_lab = q(f"""
    SELECT l.patientunitstayid AS stay_id,
           MAX(CASE WHEN LOWER(l.labname) = 'fio2'
                    THEN TRY_CAST(l.labresult AS DOUBLE) END) AS max_fio2_lab
    FROM read_csv('{EICU}/lab.csv.gz'{CSV_OPTS}) l
    JOIN cohort_stays cs ON l.patientunitstayid = cs.stay_id
    WHERE LOWER(l.labname) = 'fio2'
      AND TRY_CAST(l.labresult AS DOUBLE) BETWEEN 21 AND 100
      AND l.labresultoffset BETWEEN -360 AND 1440
    GROUP BY l.patientunitstayid
""")
# ventilator FiO2 (respiratoryCare setapneafio2, in %)
fio2_vent = q(f"""
    SELECT rc.patientunitstayid AS stay_id,
           MAX(TRY_CAST(rc.setapneafio2 AS DOUBLE)) AS max_fio2_vent
    FROM read_csv('{EICU}/respiratoryCare.csv.gz'{CSV_OPTS}) rc
    JOIN cohort_stays cs ON rc.patientunitstayid = cs.stay_id
    WHERE rc.setapneafio2 IS NOT NULL
      AND TRY_CAST(rc.setapneafio2 AS DOUBLE) BETWEEN 21 AND 100
      AND rc.respcarestatusoffset BETWEEN -360 AND 1440
    GROUP BY rc.patientunitstayid
""")
# Vent status (any active vent in window)
vent_active = q(f"""
    SELECT DISTINCT rc.patientunitstayid AS stay_id
    FROM read_csv('{EICU}/respiratoryCare.csv.gz'{CSV_OPTS}) rc
    JOIN cohort_stays cs ON rc.patientunitstayid = cs.stay_id
    WHERE (rc.airwaytype IS NOT NULL AND rc.airwaytype <> '')
      AND rc.respcarestatusoffset BETWEEN -360 AND 1440
""")
print(f"    PaO2 stays: {len(pafi):,} | FiO2_lab: {len(fio2_lab):,} | "
      f"FiO2_vent: {len(fio2_vent):,} | Vent active: {len(vent_active):,}")

# 2G. Creatinine (max over window) — used for renal SOFA only
print("  [SOFA-Renal] Creatinine extraction ...")
cr_max = q(f"""
    SELECT l.patientunitstayid AS stay_id,
           MAX(TRY_CAST(l.labresult AS DOUBLE)) AS max_cr
    FROM read_csv('{EICU}/lab.csv.gz'{CSV_OPTS}) l
    JOIN cohort_stays cs ON l.patientunitstayid = cs.stay_id
    WHERE LOWER(l.labname) = 'creatinine'
      AND TRY_CAST(l.labresult AS DOUBLE) BETWEEN 0.1 AND 30
      AND l.labresultoffset BETWEEN -360 AND 1440
    GROUP BY l.patientunitstayid
""")
print(f"    Cr stays: {len(cr_max):,}")

# 2H. Urine output total in window (full SOFA renal component)
print("  [SOFA-Renal] Urine output ...")
urine = q(f"""
    SELECT io.patientunitstayid AS stay_id,
           SUM(io.outputtotal) AS total_urine_ml
    FROM read_csv('{EICU}/intakeOutput.csv.gz'{CSV_OPTS}) io
    JOIN cohort_stays cs ON io.patientunitstayid = cs.stay_id
    WHERE LOWER(io.celllabel) LIKE '%urine%'
      AND io.outputtotal IS NOT NULL
      AND io.intakeoutputoffset BETWEEN -360 AND 1440
    GROUP BY io.patientunitstayid
""")
print(f"    Urine stays: {len(urine):,}")

# =============================================================================
# PHASE 3: Assemble SOFA scores per stay
# =============================================================================
print("\n[PHASE 3] SOFA scoring ...")

df = cohort[['patientunitstayid', 'age_num', 'kdigo_stage',
             'hospital_mortality', 'hospital_los_h',
             'apachescore', 'mech_vent_day1', 'ckd_history']].copy()
df = df.rename(columns={'patientunitstayid': 'stay_id'})

# merge
df = df.merge(map_agg,  on='stay_id', how='left')
df = df.merge(vaso,    on='stay_id', how='left')
df = df.merge(gcs,     on='stay_id', how='left')
df = df.merge(plt,     on='stay_id', how='left')
df = df.merge(bili,    on='stay_id', how='left')
df = df.merge(pafi,    on='stay_id', how='left')
df = df.merge(fio2_lab,  on='stay_id', how='left')
df = df.merge(fio2_vent, on='stay_id', how='left')
df = df.merge(vent_active, on='stay_id', how='left')
df = df.merge(cr_max,  on='stay_id', how='left')
df = df.merge(urine,   on='stay_id', how='left')
df['on_vent'] = df['stay_id'].isin(vent_active['stay_id']).astype(int)
df['any_vaso'] = ((df[['any_norepi','any_epi','any_dopa','any_dobu','any_vaso','any_milri']].fillna(0).sum(axis=1)) > 0).astype(int)

# Combined FiO2: vent value > lab value > 21 (room air)
df['max_fio2'] = df[['max_fio2_vent', 'max_fio2_lab']].max(axis=1)
df.loc[df['on_vent']==1, 'max_fio2'] = df.loc[df['on_vent']==1, 'max_fio2'].fillna(100)  # assume high FiO2 if vent but no value
df.loc[df['max_fio2'].isna(), 'max_fio2'] = 21.0
df.loc[df['min_pao2'].isna(), 'min_pao2'] = 100.0  # assume normal if no ABG
df['pf_ratio'] = df['min_pao2'] / (df['max_fio2']/100.0)

# --- SOFA Respiratory ---
def sofa_resp(pf, vent):
    if pd.isna(pf): return np.nan
    if pf >= 400: return 0
    if pf >= 300: return 1
    if pf >= 200: return 2
    if pf >= 100:
        return 3 if vent else 2
    return 4

df['sofa_resp'] = df.apply(lambda r: sofa_resp(r['pf_ratio'], r['on_vent']), axis=1)

# --- SOFA Coagulation (platelets, k/μL) ---
def sofa_coag(plt):
    if pd.isna(plt): return np.nan
    if plt >= 150: return 0
    if plt >= 100: return 1
    if plt >= 50:  return 2
    if plt >= 20:  return 3
    return 4
df['sofa_coag'] = df['min_plt'].apply(sofa_coag)

# --- SOFA Liver (bilirubin, mg/dL) ---
def sofa_liver(b):
    if pd.isna(b): return np.nan
    if b < 1.2: return 0
    if b < 2.0: return 1
    if b < 6.0: return 2
    if b < 12.0: return 3
    return 4
df['sofa_liver'] = df['max_bili'].apply(sofa_liver)

# --- SOFA Cardiovascular (MAP + vaso) ---
def sofa_cv(m, norepi, epi, dopa, dobu):
    # norepi/epi present (any dose) -> 4
    if norepi == 1 or epi == 1: return 4
    if dopa == 1: return 3
    if dobu == 1: return 2
    if pd.isna(m): return np.nan
    if m < 70: return 1
    return 0
df['sofa_cv'] = df.apply(lambda r: sofa_cv(r['min_map'], r['any_norepi'], r['any_epi'],
                                              r['any_dopa'], r['any_dobu']), axis=1)

# --- SOFA CNS (GCS) ---
def sofa_cns(g):
    if pd.isna(g): return np.nan
    if g >= 14: return 0
    if g >= 12: return 1
    if g >= 10: return 2
    if g >= 6:  return 3
    return 4
df['sofa_cns'] = df['min_gcs_total'].apply(sofa_cns)

# --- SOFA Renal (Cr + urine) ---
def sofa_renal(cr, u_ml):
    # cr (mg/dL) — primary criterion
    if not pd.isna(cr):
        if cr < 1.2: cr_score = 0
        elif cr < 2.0: cr_score = 1
        elif cr < 3.5: cr_score = 2
        elif cr < 5.0: cr_score = 3
        else: cr_score = 4
    else:
        cr_score = np.nan
    # urine ml/day in window — adjusted to 30h window (6h pre + 24h post = 30h)
    # mimic-code uses daily urine; if window is 30h, we scale to 24h
    if not pd.isna(u_ml):
        u_24h = u_ml * (24/30)
        if u_24h < 200: u_score = 4
        elif u_24h < 500: u_score = 3
        else: u_score = 0
    else:
        u_score = np.nan
    # Take max
    if pd.isna(cr_score) and pd.isna(u_score): return np.nan
    if pd.isna(cr_score): return u_score
    if pd.isna(u_score): return cr_score
    return max(cr_score, u_score)
df['sofa_renal'] = df.apply(lambda r: sofa_renal(r['max_cr'], r['total_urine_ml']), axis=1)

# non-renal SOFA = sum of resp/coag/liver/cv/cns (max 20)
df['sofa_nonrenal'] = df[['sofa_resp','sofa_coag','sofa_liver','sofa_cv','sofa_cns']].sum(axis=1, min_count=1)
# full SOFA = non-renal + renal (max 24)
df['sofa_full'] = df[['sofa_resp','sofa_coag','sofa_liver','sofa_cv','sofa_cns','sofa_renal']].sum(axis=1, min_count=1)
# count missing components (1..6)
def count_missing(r):
    return sum(pd.isna(r[c]) for c in ['sofa_resp','sofa_coag','sofa_liver','sofa_cv','sofa_cns','sofa_renal'])
df['sofa_n_missing'] = df.apply(count_missing, axis=1)

# Median imputation for nonrenal (sensitivity-aware): if sofa_nonrenal is null but
# all components except 1 missing, we can fill; otherwise keep null
# For primary analysis, require non-missing nonrenal -> drop rows
df['sofa_nonrenal'] = df['sofa_nonrenal'].fillna(df['sofa_nonrenal'].median())

print(f"  SOFA nonrenal distribution (after median-fill for nulls):")
print(df['sofa_nonrenal'].describe().round(2).to_string())
print(f"\n  SOFA full distribution:")
print(df['sofa_full'].describe().round(2).to_string())
print(f"\n  SOFA missing count distribution:")
print(df['sofa_n_missing'].value_counts().sort_index().to_string())

# =============================================================================
# PHASE 4: 5 comorbidities from pastHistory
# =============================================================================
print("\n[PHASE 4] 5 comorbidities from pastHistory table")
past = q(f"""
    SELECT patientunitstayid AS stay_id,
           MAX(CASE WHEN pasthistorypath LIKE '%/Hypertension Requiring Treatment/%'
                    THEN 1 ELSE 0 END) AS hypertension,
           MAX(CASE WHEN pasthistorypath LIKE '%/Insulin Dependent Diabetes/%'
                     OR pasthistorypath LIKE '%/Non-Insulin Dependent Diabetes/%'
                    THEN 1 ELSE 0 END) AS diabetes,
           MAX(CASE WHEN pasthistorypath LIKE '%/Congestive Heart Failure/%' OR pasthistorypath LIKE '%/CHF%'
                    THEN 1 ELSE 0 END) AS heart_failure,
           MAX(CASE WHEN pasthistorypath LIKE '%/Pulmonary/COPD/%'
                    THEN 1 ELSE 0 END) AS copd,
           MAX(CASE WHEN pasthistorypath LIKE '%/Gastrointestinal (R)/Cirrhosis/%'
                    THEN 1 ELSE 0 END) AS liver_disease
    FROM read_csv('{EICU}/pastHistory.csv.gz'{CSV_OPTS})
    WHERE patientunitstayid IN (SELECT stay_id FROM cohort_stays)
    GROUP BY patientunitstayid
""")
df = df.merge(past, on='stay_id', how='left')
for c in ['hypertension','heart_failure','copd','liver_disease','diabetes']:
    df[c] = df[c].fillna(0).astype(int)
print(f"  Comorbidities: HTN {df.hypertension.sum():,}, DM {df.diabetes.sum():,}, "
      f"CHF {df.heart_failure.sum():,}, COPD {df.copd.sum():,}, Liver {df.liver_disease.sum():,}")

# =============================================================================
# PHASE 5: Sex + event time + final harmonisation
# =============================================================================
print("\n[PHASE 5] Sex + event time harmonisation")
df['sex_male'] = (cohort['gender'] == 'Male').astype(int).values
df['vent_24h'] = cohort['mech_vent_day1'].fillna(0).astype(int).values
df['ckd_history_prior'] = cohort['ckd_history'].fillna(0).astype(int).values

# Real event time
df['event_time_h'] = cohort['hospital_los_h'].fillna(-1).values
df['event_ind'] = cohort['hospital_mortality'].fillna(0).astype(int).values

# Rename age, aki_onset_h
df['age'] = cohort['age_num'].values
df['aki_onset_h'] = cohort['aki_onset_h'].values
df['baseline_cr'] = cohort['baseline_cr'].values
df['baseline_source'] = cohort['baseline_source'].values
df['hospital_expire_flag'] = df['event_ind']

# Cohort V1 had vent_24h with 0/1 — verify
print(f"  vent_24h sum: {df.vent_24h.sum():,} of {len(df):,} ({df.vent_24h.mean()*100:.1f}%)")
print(f"  sex_male sum: {df.sex_male.sum():,}")
print(f"  event_ind sum: {df.event_ind.sum():,}")
print(f"  event_time_h range: {df.event_time_h.min():.1f} to {df.event_time_h.max():.1f}")

# Sanity check: sofa_nonrenal vs KDIGO stage
print("\n  SOFA nonrenal by KDIGO stage:")
print(df.groupby('kdigo_stage')['sofa_nonrenal'].agg(['count','mean','std']).round(2).to_string())

# =============================================================================
# PHASE 6: Save
# =============================================================================
print("\n[PHASE 6] Save")

out_cols = ['stay_id', 'age', 'sex_male', 'kdigo_stage', 'aki_onset_h',
            'baseline_cr', 'baseline_source',
            'sofa_nonrenal', 'sofa_full', 'sofa_n_missing',
            'vent_24h', 'ckd_history_prior',
            'diabetes', 'hypertension', 'heart_failure', 'copd', 'liver_disease',
            'hospital_expire_flag', 'event_time_h', 'event_ind']
df = df[out_cols].copy()
df.to_csv(os.path.join(OUT, "eicu_analysis_v2.csv"), index=False)
print(f"  Saved: eicu_analysis_v2.csv ({len(df):,} rows, {len(df.columns)} cols)")

# Also create identical column structure for use in F7 external validation
meta = {
    'cohort_n': int(len(df)),
    'n_deaths': int(df.event_ind.sum()),
    'columns': list(df.columns),
    'median_event_time_h': float(df.event_time_h.median()),
    'kdigo_distribution': {int(s): int((df.kdigo_stage==s).sum()) for s in [0,1,2,3]},
    'mortality_by_stage': {int(s): float(df.loc[df.kdigo_stage==s,'event_ind'].mean()*100) for s in [0,1,2,3]},
    'sofa_nonrenal_by_stage': {int(s): float(df.loc[df.kdigo_stage==s,'sofa_nonrenal'].mean()) for s in [0,1,2,3]},
    'comorbidity_freq': {
        c: float(df[c].mean()*100) for c in
        ['diabetes','hypertension','heart_failure','copd','liver_disease','vent_24h','ckd_history_prior','sex_male']
    },
}
with open(os.path.join(OUT, "eicu_analysis_v2_meta.json"), 'w') as f:
    json.dump(meta, f, indent=2)
print(f"  Saved: eicu_analysis_v2_meta.json")
print(f"\nDone: {datetime.now()}")