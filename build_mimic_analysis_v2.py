"""
构建 MIMIC v2 分析就绪数据集
合并: kdigo_cohort_v2 + SOFA(已就位) + 协变量(性别/通气24h/合并症) + 真实事件时间
输出: mimic_analysis_v2.csv + mimic_analysis_v2_meta.json
"""
import json
import os
import duckdb
import pandas as pd
import numpy as np

HOSP = "E:/mimic-iv/v3.1/physionet.org/files/mimiciv/3.1/hosp"
ICU = "E:/mimic-iv/v3.1/physionet.org/files/mimiciv/3.1/icu"
WORK = r"C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"

def q(sql):
    return duckdb.query(sql).df()

print("=" * 70)
print("MIMIC v2 分析数据集构建")
print("=" * 70)

# --- Base cohort ---
cohort = pd.read_csv(os.path.join(WORK, "kdigo_cohort_v2.csv"))
print(f"\n[1] Base cohort: {len(cohort):,}")

# --- Sex (gender M/F -> 1/0) ---
print("\n[2] Encode sex")
cohort['sex_male'] = (cohort['gender'] == 'M').astype(int)
print(f"  M: {(cohort.sex_male==1).sum():,}, F: {(cohort.sex_male==0).sum():,}")

# --- Ventilation in first 24h of ICU (procedureevents, time-filtered) ---
# F4 fix: restrict to procedures STARTED within [intime, intime+24h] of the index stay.
# (Previous version had no time filter -> whole-stay ventilation.)
print("\n[3] Mechanical ventilation first 24h of ICU stay")
duckdb.register('cohort_stays',
                cohort[['stay_id', 'intime']].copy())
vent_raw = q(f"""
    SELECT p.stay_id, p.starttime, c.intime
    FROM read_csv_auto('{ICU}/procedureevents.csv.gz') p
    INNER JOIN cohort_stays c USING (stay_id)
    WHERE p.itemid IN (
        225792, -- Invasive Ventilation
        224422, -- Intubation
        227011, -- Ventilation (CareVue)
        225468, -- Ventilator Care
        224385  -- Intubation (CareVue)
    )
""")
vent_raw['starttime'] = pd.to_datetime(vent_raw['starttime'])
vent_raw['intime']    = pd.to_datetime(vent_raw['intime'])
dt_h = (vent_raw['starttime'] - vent_raw['intime']).dt.total_seconds() / 3600.0
vent_24h_ids = set(vent_raw.loc[(dt_h >= 0) & (dt_h <= 24), 'stay_id'])
cohort['vent_24h'] = cohort['stay_id'].isin(vent_24h_ids).astype(int)
print(f"  ventilated within first 24h: {cohort['vent_24h'].sum():,} "
      f"({cohort['vent_24h'].mean()*100:.1f}%)  "
      f"[any-time during stay would be {cohort['stay_id'].isin(set(vent_raw['stay_id'])).sum():,}]")

# --- Comorbidities from index-admission diagnoses ---
# Elixhauser-style simplified set used in literature:
#   Diabetes: ICD-9 250* / ICD-10 E10*, E11*
#   Hypertension: ICD-9 401*-405* / ICD-10 I10*-I15*
#   Heart Failure: ICD-9 428* / ICD-10 I50*
#   COPD: ICD-9 490*-496* / ICD-10 J44*
#   Liver disease (moderate-severe): ICD-9 570*, 571*, 572* / ICD-10 K70*-K74*
print("\n[4] Comorbidities from index-admission ICD codes")
adm_ids_list = cohort['hadm_id'].astype('int64').unique().tolist()
print(f"  index hadm_ids: {len(adm_ids_list):,}")
# Use a temp view to avoid giant IN list
duckdb.register('cohort_adm', pd.DataFrame({'hadm_id': adm_ids_list}))

diag = q(f"""
    SELECT d.hadm_id, d.icd_version, d.icd_code
    FROM read_csv_auto('{HOSP}/diagnoses_icd.csv.gz') d
    INNER JOIN cohort_adm c USING (hadm_id)
""")
print(f"  diagnosis rows: {len(diag):,}")

def diag_match(pat9, pat10):
    return ((diag.icd_version==9)  & diag.icd_code.str.match(pat9,  na=False)) | \
           ((diag.icd_version==10) & diag.icd_code.str.match(pat10, na=False))

flags = {}
flags['diabetes']      = diag_match(r'^250',  r'^E1[01]')
flags['heart_failure']  = diag_match(r'^428',  r'^I50')
flags['copd']          = diag_match(r'^49[0-6]', r'^J44')
flags['hypertension']  = diag_match(r'^40[1-5]', r'^I1[0-5]')
flags['liver_disease'] = diag_match(r'^5(70|71|72)', r'^K7[0-4]')

for c, mask in flags.items():
    hadm_set = set(diag.loc[mask, 'hadm_id'])
    cohort[c] = cohort['hadm_id'].isin(hadm_set).astype(int)
    print(f"  {c}: {cohort[c].sum():,} ({cohort[c].mean()*100:.1f}%)")

# --- Real event times (deathtime-admittime for deaths; dischtime for survivors) ---
print("\n[5] Real event times (deathtime for deaths, dischtime for survivors)")
cohort['intime']    = pd.to_datetime(cohort['intime'])
cohort['admittime'] = pd.to_datetime(cohort['admittime'])
cohort['deathtime'] = pd.to_datetime(cohort['deathtime'])
cohort['dischtime'] = pd.to_datetime(cohort['dischtime'])

cohort['event_time_h'] = np.where(
    cohort['hospital_expire_flag']==1,
    (cohort['deathtime'] - cohort['intime']).dt.total_seconds()/3600,
    (cohort['dischtime'] - cohort['intime']).dt.total_seconds()/3600
)
cohort['event_ind'] = (cohort['hospital_expire_flag']==1).astype(int)

# Validate deathtime within hospitalization (rare anomalies)
bad = ((cohort['hospital_expire_flag']==1) &
       (cohort['deathtime'] > cohort['dischtime'] + pd.Timedelta(hours=2))).sum()
print(f"  deaths with deathtime AFTER discharge (anomaly): {bad}")

# Time must be positive
print(f"  event_time_h range: {cohort['event_time_h'].min():.1f} to {cohort['event_time_h'].max():.1f}")
print(f"  deaths: {cohort['event_ind'].sum():,}")
print(f"  median time_h: {cohort['event_time_h'].median():.1f}")

# --- Save ---
out_cols = ['subject_id','hadm_id','stay_id','age','sex_male',
            'kdigo_stage','aki_onset_h','baseline_cr','baseline_source',
            'sofa_nonrenal','sofa_full','sofa_n_missing',
            'vent_24h','ckd_history_prior','diabetes','hypertension',
            'heart_failure','copd','liver_disease',
            'hospital_expire_flag','event_time_h','event_ind']
out = cohort[out_cols].copy()
out.to_csv(os.path.join(WORK, "mimic_analysis_v2.csv"), index=False)
print(f"\n  mimic_analysis_v2.csv saved ({len(out):,} rows, {len(out.columns)} cols)")

# Column reference JSON for the survival script
meta = {
    'cohort_n': int(len(out)),
    'n_deaths': int(out.event_ind.sum()),
    'columns': list(out.columns),
    'median_event_time_h': float(out.event_time_h.median()),
    'kdigo_distribution': {
            int(s): int((out.kdigo_stage==s).sum()) for s in [0,1,2,3]
        },
    'mortality_by_stage': {
            int(s): float(out.loc[out.kdigo_stage==s,'event_ind'].mean()*100) for s in [0,1,2,3]
        },
}
with open(os.path.join(WORK, "mimic_analysis_v2_meta.json"), 'w') as f:
    json.dump(meta, f, indent=2)
print(f"  meta saved")