#!/usr/bin/env python3
"""
eICU-CRD KDIGO AKI 分期 v2 — 队列构建 (与 MIMIC kdigo_cohort_v2.py 规则统一)
================================================================================
对照外部评审报告 (2026-09-23) F2/F5/F7 修复，消除两库流程不一致：

  F2-a: 队列单位 = 患者 — 每位 uniquepid 仅保留首次合格 unit stay
        (v1 为 200,859 全部住院级 stay，含同一患者多次入组)
  F2-b: 无 SCr 患者明确排除并计数 (不再 MDRD 反推 + 强制 Stage 0)
  F2-c: ESRD/慢性透析排除 — ICD 585.6/N18.6/V45.1/Z99.2 (本次住院诊断)
        + pastHistory 慢性透析 (hemodialysis/peritoneal dialysis)
        + pastHistory "renal failure - not currently dialyzed" (ESRD 等价)
  F5:   真实随访时间 — unitdischargeoffset / hospitaldischargeoffset (分钟级
        真实偏移)，不再用 ICU LOS 代理死亡时间；同时记录 unit/hospital 两级死亡率
  F7:   age '> 89' → 90 (v1 直接排除 7,081 例)
        AKI 分期采用与 MIMIC v2 完全相同的严格算法 (48h 滚动窗 + 7d 比值窗
        + SCr>=4.0 需急性变化 + RRT 仅计 index stay 前 7 天)

基线肌酐层级 (eICU 无住院前门诊化验，无法复制 MIMIC 第 1 层)：
  1. ICU 入科前 7 天内 SCr 中位数 (labresultoffset ∈ [max(hospitaladmitoffset,
     -10080), 0)，即本院住院期至入科前)
  2. ICU 入科后 24h 内首值 (offset ∈ [0, 1440])
  → baseline_source 记录来源；两层均无 → 排除

产出: eicu_cohort_v2.csv + eicu_cohort_v2_flow.json
"""
import duckdb
import pandas as pd
import numpy as np
import os, json, gc
from datetime import datetime

print("=" * 80)
print("eICU-CRD KDIGO-AKI v2: Cohort Construction (harmonised with MIMIC v2)")
print(f"Start: {datetime.now()}")
print("=" * 80)

EICU = "E:/mimic-iv/eicu/physionet.org/files/eicu-crd/2.0"
OUT = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"

CR_MIN, CR_MAX = 0.1, 30.0        # physiologic plausibility bounds (mg/dL)
WINDOW_RATIO_MIN = 7 * 1440       # ratio criteria window: 7 days from ICU admission (min)
WINDOW_RISE_MIN = 48 * 60         # rolling window for +0.3 mg/dL criterion (min)
MIN_ICU_LOS_MIN = 12 * 60         # minimum ICU stay: 12h (min)
STAGE3_ABS_CUT = 4.0              # mg/dL absolute Stage 3 threshold

# eICU CSV 引号/脏行处理 (pastHistory 等表存在不规范行)
CSV_OPTS = ", quote='\"', ignore_errors=true"

def q(sql):
    return duckdb.query(sql).df()

flow = {}

# ============================================================
# PHASE 1: patient 表 — 成人 + LOS>=12h + 患者级去重 (F2-a, F7)
# ============================================================
print("\n" + "=" * 60)
print("PHASE 1: patient table — adults, LOS>=12h, first eligible stay per PATIENT")
print("=" * 60)

patient = q(f"""
    SELECT patientunitstayid, uniquepid, patienthealthsystemstayid,
           gender, age, ethnicity, hospitalid, unittype, unitstaytype,
           hospitaladmitoffset, unitdischargeoffset, hospitaldischargeoffset,
           unitdischargestatus, hospitaldischargestatus
    FROM read_csv('{EICU}/patient.csv.gz'{CSV_OPTS})
""")

flow['total_unit_stays'] = len(patient)
flow['total_patients'] = patient['uniquepid'].nunique()
print(f"  Total unit stays: {flow['total_unit_stays']:,} "
      f"({flow['total_patients']:,} unique patients)")

# Age: '> 89' -> 90 (F7); blank/null -> excluded (cannot verify adult)
age_raw = patient['age'].astype(str).str.strip()
age_num = pd.to_numeric(age_raw, errors='coerce')
n_gt89 = (age_raw == '> 89').sum()
age_num[age_raw == '> 89'] = 90
patient['age_num'] = age_num
flow['age_gt89_recoded'] = int(n_gt89)
flow['excluded_age_missing'] = int(patient['age_num'].isna().sum())

adult = patient[patient['age_num'] >= 18].copy()
flow['after_adult'] = len(adult)
flow['excluded_age_lt18'] = int((patient['age_num'] < 18).sum())
print(f"  Adults (>=18, '> 89' recoded to 90): {flow['after_adult']:,} "
      f"(excluded: {flow['excluded_age_lt18']:,} minors, "
      f"{flow['excluded_age_missing']:,} missing age)")

# ICU LOS >= 12h (unitdischargeoffset is minutes from unit admission)
adult['icu_los_h'] = pd.to_numeric(adult['unitdischargeoffset'], errors='coerce') / 60.0
eligible = adult[adult['icu_los_h'] >= MIN_ICU_LOS_MIN / 60].copy()
flow['after_los_12h'] = len(eligible)
flow['excluded_los_lt12h'] = flow['after_adult'] - flow['after_los_12h']
print(f"  After ICU LOS >= 12h: {flow['after_los_12h']:,} stays "
      f"(excluded {flow['excluded_los_lt12h']:,})")

# --- F2-a: first ELIGIBLE unit stay per uniquepid ---
# eICU has no absolute timestamps; patientunitstayid is assigned in
# chronological order of unit admission — used as the ordering key.
eligible = eligible.sort_values(['uniquepid', 'patientunitstayid'])
first_stay = eligible.drop_duplicates('uniquepid', keep='first').copy()
flow['after_patient_dedup'] = len(first_stay)
flow['excluded_repeat_stays'] = len(eligible) - len(first_stay)
print(f"  First eligible unit stay per patient: {flow['after_patient_dedup']:,} "
      f"(dropped {flow['excluded_repeat_stays']:,} later stays)")

cohort = first_stay.copy()
cohort_stays = set(cohort['patientunitstayid'])

# ============================================================
# PHASE 2: ESRD/慢性透析排除 (F2-c) + CKD 病史标记
# ============================================================
print("\n" + "=" * 60)
print("PHASE 2: ESRD/chronic dialysis exclusion + CKD history flag")
print("=" * 60)

# (a) Diagnosis ICD codes (ICD-9,ICD-10 pairs stored together, e.g. '585.6, N18.6')
esrd_icd = q(f"""
    SELECT DISTINCT patientunitstayid
    FROM read_csv('{EICU}/diagnosis.csv.gz'{CSV_OPTS})
    WHERE regexp_matches(icd9code, '585[.]6|N18[.]6|V45[.]1|Z99[.]2')
""")
esrd_icd_set = set(esrd_icd['patientunitstayid']) & cohort_stays

# (b) pastHistory: chronic dialysis / ESRD-equivalent (strictly PRE-EXISTING conditions)
esrd_ph = q(f"""
    SELECT DISTINCT patientunitstayid
    FROM read_csv('{EICU}/pastHistory.csv.gz'{CSV_OPTS})
    WHERE pasthistorypath ILIKE '%renal failure - hemodialysis%'
       OR pasthistorypath ILIKE '%renal failure - peritoneal dialysis%'
       OR pasthistorypath ILIKE '%renal failure- not currently dialyzed%'
""")
esrd_ph_set = set(esrd_ph['patientunitstayid']) & cohort_stays

esrd_all = esrd_icd_set | esrd_ph_set
cohort = cohort[~cohort['patientunitstayid'].isin(esrd_all)].copy()
flow['after_esrd_exclusion'] = len(cohort)
flow['excluded_esrd_dialysis'] = flow['after_patient_dedup'] - len(cohort)
flow['esrd_by_icd'] = len(esrd_icd_set)
flow['esrd_by_pasthistory'] = len(esrd_ph_set)
print(f"  After ESRD/chronic-dialysis exclusion: {flow['after_esrd_exclusion']:,} "
      f"(excluded {flow['excluded_esrd_dialysis']:,}: "
      f"{flow['esrd_by_icd']:,} by ICD, {flow['esrd_by_pasthistory']:,} by past history)")

# CKD stage 3-5 history: ICD 585.3-585.5/N18.3-N18.5 + pastHistory baseline SCr 3-5
ckd_icd = q(f"""
    SELECT DISTINCT patientunitstayid
    FROM read_csv('{EICU}/diagnosis.csv.gz'{CSV_OPTS})
    WHERE regexp_matches(icd9code, '585[.][345]|N18[.][345]')
""")
ckd_ph = q(f"""
    SELECT DISTINCT patientunitstayid
    FROM read_csv('{EICU}/pastHistory.csv.gz'{CSV_OPTS})
    WHERE pasthistorypath ILIKE '%creatinine 3-4%'
       OR pasthistorypath ILIKE '%creatinine 4-5%'
       OR pasthistorypath ILIKE '%creatinine > 5%'
""")
ckd_set = (set(ckd_icd['patientunitstayid']) | set(ckd_ph['patientunitstayid'])) & set(cohort['patientunitstayid'])
cohort['ckd_history'] = cohort['patientunitstayid'].isin(ckd_set).astype(int)
print(f"  CKD stage 3-5 history: {cohort['ckd_history'].sum():,} "
      f"({cohort['ckd_history'].mean()*100:.1f}%)")

cohort_stays = set(cohort['patientunitstayid'])
del esrd_icd, esrd_ph, ckd_icd, ckd_ph; gc.collect()

# ============================================================
# PHASE 3: 肌酐化验提取 (F1 数据基础)
# ============================================================
print("\n" + "=" * 60)
print("PHASE 3: lab — all serum creatinine for cohort stays")
print("=" * 60)

cr_df = q(f"""
    SELECT patientunitstayid, labresultoffset, labresult
    FROM read_csv('{EICU}/lab.csv.gz'{CSV_OPTS})
    WHERE labname = 'creatinine'
      AND labresult >= {CR_MIN}
      AND labresult <= {CR_MAX}
      AND labresultoffset IS NOT NULL
""")
cr_df = cr_df[cr_df['patientunitstayid'].isin(cohort_stays)]
cr_df['labresultoffset'] = cr_df['labresultoffset'].astype(float)
cr_df = cr_df.sort_values(['patientunitstayid', 'labresultoffset'])
print(f"  Creatinine rows for cohort stays: {len(cr_df):,} "
      f"across {cr_df['patientunitstayid'].nunique():,} stays")

# ============================================================
# PHASE 4: 基线肌酐层级 (F1)
# ============================================================
print("\n" + "=" * 60)
print("PHASE 4: Baseline creatinine hierarchy")
print("=" * 60)

# Attach hospital admission offset (negative = before ICU admission)
cr_idx = cr_df.merge(
    cohort[['patientunitstayid', 'hospitaladmitoffset']].rename(
        columns={'hospitaladmitoffset': 'hosp_adm_off'}),
    on='patientunitstayid', how='inner')

off = cr_idx['labresultoffset']
hosp = cr_idx['hosp_adm_off'].astype(float)
# pre-ICU window: [max(hospital admit offset, -7d), 0)
lo = np.maximum(hosp.fillna(-WINDOW_RATIO_MIN), -WINDOW_RATIO_MIN)
pre_icu = cr_idx[(off >= lo) & (off < 0)]
# ICU first 24h
icu_24h = cr_idx[(off >= 0) & (off <= 1440)]

bl_pre = pre_icu.groupby('patientunitstayid')['labresult'].median()
bl_24 = icu_24h.sort_values('labresultoffset').groupby('patientunitstayid')['labresult'].first()

cohort['baseline_cr'] = np.nan
cohort['baseline_source'] = None
for src, series in [('preICU_0_7d', bl_pre), ('icu_first24h', bl_24)]:
    need = cohort['baseline_cr'].isna()
    cohort.loc[need, 'baseline_cr'] = cohort.loc[need, 'patientunitstayid'].map(series)
    cohort.loc[need & cohort['baseline_cr'].notna(), 'baseline_source'] = src

n_bl = int(cohort['baseline_cr'].notna().sum())
print(f"  Baseline available: {n_bl:,}/{len(cohort):,}")
print("  Baseline source distribution:")
for src, cnt in cohort['baseline_source'].value_counts(dropna=False).items():
    print(f"    {src}: {cnt:,} ({cnt/len(cohort)*100:.1f}%)")
flow['after_baseline_available'] = n_bl
flow['excluded_no_baseline'] = int(len(cohort) - n_bl)

# ============================================================
# PHASE 5: ICU 入科后 7 天窗口 + 严格 KDIGO 分期 (F1)
# ============================================================
print("\n" + "=" * 60)
print("PHASE 5: Strict KDIGO staging in [0, +7d] from unit admission")
print("=" * 60)

cr_win = cr_idx[(cr_idx['labresultoffset'] >= 0)
                & (cr_idx['labresultoffset'] <= WINDOW_RATIO_MIN)][
    ['patientunitstayid', 'labresultoffset', 'labresult']].copy()
n_win = cr_win['patientunitstayid'].nunique()
print(f"  Stays with >=1 SCr in first 7 ICU days: {n_win:,}")
flow['after_scr_in_window'] = int(n_win)
flow['excluded_no_window_scr'] = int(n_bl - n_win)

stageable = cohort[cohort['baseline_cr'].notna()
                   & cohort['patientunitstayid'].isin(set(cr_win['patientunitstayid']))].copy()
flow['final_cohort'] = len(stageable)
print(f"  Final stageable cohort: {flow['final_cohort']:,} patients")

cr_win = cr_win[cr_win['patientunitstayid'].isin(set(stageable['patientunitstayid']))]
cr_win = cr_win.sort_values(['patientunitstayid', 'labresultoffset'])

bl_lookup = stageable.set_index('patientunitstayid')['baseline_cr']

def stage_stay(gr, base):
    """Strict KDIGO staging (identical algorithm to MIMIC v2, minutes)."""
    t = gr['labresultoffset'].values
    v = gr['labresult'].values
    n = len(v)

    max_ratio = 0.0
    max_rise = 0.0
    stage_traj = np.zeros(n, dtype=int)
    stage3_abs = False
    stage1_ratio = False

    for j in range(n):
        ratio = v[j] / base if base > 0 else 0.0
        if ratio > max_ratio:
            max_ratio = ratio
        if 1.5 <= ratio < 2.0:
            stage1_ratio = True
        s_ratio = 3 if ratio >= 3.0 else (2 if ratio >= 2.0 else (1 if ratio >= 1.5 else 0))

        # +0.3 mg/dL within rolling 48h window ending at t_j
        rise = 0.0
        i = j
        while i >= 0 and (t[j] - t[i]) <= WINDOW_RISE_MIN:
            d = v[j] - v[i]
            if d > rise:
                rise = d
            i -= 1
        if rise > max_rise:
            max_rise = rise
        s_rise = 1 if rise >= 0.3 else 0

        # absolute Stage 3 requires acute change (F1-c)
        s_abs = 0
        if v[j] >= STAGE3_ABS_CUT and (s_ratio >= 1 or s_rise >= 1):
            s_abs = 3
            stage3_abs = True

        stage_traj[j] = max(s_ratio, s_rise, s_abs)

    stage = int(stage_traj.max()) if n > 0 else 0
    onset_h = np.nan
    if stage > 0:
        first_idx = int(np.argmax(stage_traj >= stage))
        onset_h = t[first_idx] / 60.0

    return pd.Series({'kdigo_stage': stage,
                      'aki_onset_h': onset_h,
                      'max_cr_7d': v.max(),
                      'max_ratio_7d': max_ratio,
                      'max_rise48h': max_rise,
                      'stage3_abs_aki': int(stage3_abs),
                      'stage1_ratio_met': int(stage1_ratio),
                      'n_cr_window': n})

results = []
staged_ids = []
for sid, gr in cr_win.groupby('patientunitstayid'):
    staged_ids.append(sid)
    results.append(stage_stay(gr, bl_lookup.loc[sid]))
    if len(results) % 20000 == 0:
        print(f"    ... staged {len(results):,} stays")

staged = pd.DataFrame(results, index=pd.Index(staged_ids, name='patientunitstayid'))
print(f"  Staged: {len(staged):,} stays")

# ============================================================
# PHASE 6: RRT — 仅 index stay 入科后 7 天内 (F1/F4)
# ============================================================
print("\n" + "=" * 60)
print("PHASE 6: RRT within first 7 ICU days (index stay only)")
print("=" * 60)

# Dialysis treatments; exclude catheter-insertion/radiology entries (not RRT itself)
rrt_df = q(f"""
    SELECT patientunitstayid, treatmentoffset
    FROM read_csv('{EICU}/treatment.csv.gz'{CSV_OPTS})
    WHERE (treatmentstring ILIKE '%dialysis%' OR treatmentstring ILIKE '%crrt%'
           OR treatmentstring ILIKE '%ultrafiltration%')
      AND treatmentstring NOT ILIKE '%catheter%'
      AND treatmentstring NOT ILIKE '%radiology%'
""")
rrt_df = rrt_df[rrt_df['patientunitstayid'].isin(set(stageable['patientunitstayid']))]
rrt_df['treatmentoffset'] = pd.to_numeric(rrt_df['treatmentoffset'], errors='coerce')
rrt_7d = rrt_df[(rrt_df['treatmentoffset'] >= 0)
                & (rrt_df['treatmentoffset'] <= WINDOW_RATIO_MIN)]
rrt_ids = set(rrt_7d['patientunitstayid'])
print(f"  RRT within first 7 ICU days: {len(rrt_ids):,} stays "
      f"(any-time RRT would have been {rrt_df['patientunitstayid'].nunique():,})")

cohort_v2 = stageable.merge(staged, left_on='patientunitstayid',
                            right_index=True, how='left')
cohort_v2['rrt_7d'] = cohort_v2['patientunitstayid'].isin(rrt_ids).astype(int)
rrt_mask = (cohort_v2['rrt_7d'] == 1) & (cohort_v2['kdigo_stage'] < 3)
cohort_v2.loc[rrt_mask, 'kdigo_stage'] = 3
cohort_v2['stage3_by_rrt'] = ((cohort_v2['rrt_7d'] == 1) &
                              (cohort_v2['stage3_abs_aki'] == 0) &
                              (cohort_v2['max_ratio_7d'] < 3.0)).astype(int)

# ============================================================
# PHASE 7: 协变量 + 结局 (F5: 真实随访时间)
# ============================================================
print("\n" + "=" * 60)
print("PHASE 7: Covariates + outcomes (real follow-up times)")
print("=" * 60)

# APACHE IVa score (no imputation here — downstream stats decide)
apache = q(f"""
    SELECT patientunitstayid, apachescore
    FROM read_csv('{EICU}/apachePatientResult.csv.gz'{CSV_OPTS})
    WHERE apacheversion = 'IVa' AND apachescore IS NOT NULL
""")
apache = apache.groupby('patientunitstayid')['apachescore'].first().reset_index()
cohort_v2 = cohort_v2.merge(apache, on='patientunitstayid', how='left')

# Ventilation on ICU day 1 (index stay)
vent = q(f"""
    SELECT patientunitstayid, ventday1
    FROM read_csv('{EICU}/apachePredVar.csv.gz'{CSV_OPTS})
""")
vent = vent.groupby('patientunitstayid')['ventday1'].first().reset_index()
cohort_v2 = cohort_v2.merge(vent, on='patientunitstayid', how='left')
cohort_v2['mech_vent_day1'] = np.where(cohort_v2['ventday1'].isna(), np.nan,
                                       (cohort_v2['ventday1'] > 0).astype(float))

# Outcomes (F5): real follow-up times in hours + unit/hospital mortality
cohort_v2['unit_mortality'] = (cohort_v2['unitdischargestatus'] == 'Expired').astype(int)
cohort_v2['hospital_mortality'] = (cohort_v2['hospitaldischargestatus'] == 'Expired').astype(int)
cohort_v2['hospital_los_h'] = pd.to_numeric(
    cohort_v2['hospitaldischargeoffset'], errors='coerce') / 60.0

# ============================================================
# PHASE 8: 汇总 + 保存
# ============================================================
print("\n" + "=" * 60)
print("PHASE 8: Summary + save")
print("=" * 60)

final = cohort_v2[cohort_v2['kdigo_stage'].notna()].copy()
final['kdigo_stage'] = final['kdigo_stage'].astype(int)

print("\n  FINAL eICU COHORT (patient-level):")
print(f"  N = {len(final):,}")
print("\n  KDIGO stage distribution:")
summary = final.groupby('kdigo_stage').agg(
    n=('kdigo_stage', 'size'),
    unit_mort=('unit_mortality', 'mean'),
    hosp_mort=('hospital_mortality', 'mean'))
for stage, row in summary.iterrows():
    print(f"    Stage {stage}: {row['n']:>7,} ({row['n']/len(final)*100:5.1f}%)  "
          f"unit mort {row['unit_mort']*100:5.1f}%  hosp mort {row['hosp_mort']*100:5.1f}%")
print(f"  Overall unit mortality: {final['unit_mortality'].mean()*100:.1f}%")
print(f"  Overall hospital mortality: {final['hospital_mortality'].mean()*100:.1f}%")
print(f"  AKI (stage >= 1): {(final['kdigo_stage']>=1).mean()*100:.1f}%")
print(f"  APACHE IVa available: {final['apachescore'].notna().mean()*100:.1f}%")
print(f"  Vent day-1 available: {final['mech_vent_day1'].notna().mean()*100:.1f}%")

out_cols = ['patientunitstayid', 'uniquepid', 'patienthealthsystemstayid',
            'gender', 'age_num', 'ethnicity', 'hospitalid', 'unittype', 'unitstaytype',
            'unitdischargeoffset', 'hospitaldischargeoffset',
            'unitdischargestatus', 'hospitaldischargestatus',
            'unit_mortality', 'hospital_mortality', 'icu_los_h', 'hospital_los_h',
            'baseline_cr', 'baseline_source', 'ckd_history',
            'n_cr_window', 'max_cr_7d', 'max_ratio_7d', 'max_rise48h',
            'kdigo_stage', 'aki_onset_h', 'stage3_abs_aki', 'stage1_ratio_met',
            'rrt_7d', 'stage3_by_rrt',
            'apachescore', 'mech_vent_day1']
final[out_cols].to_csv(os.path.join(OUT, "eicu_cohort_v2.csv"), index=False)

flow['final_cohort'] = len(final)
flow['kdigo_distribution'] = {int(k): int(v) for k, v in
                              final['kdigo_stage'].value_counts().sort_index().items()}
flow['unit_mortality'] = float(final['unit_mortality'].mean())
flow['hospital_mortality'] = float(final['hospital_mortality'].mean())
flow['aki_prevalence'] = float((final['kdigo_stage'] >= 1).mean())
flow['mortality_by_stage_hospital'] = {
    int(k): float(v) for k, v in
    final.groupby('kdigo_stage')['hospital_mortality'].mean().items()}
flow['mortality_by_stage_unit'] = {
    int(k): float(v) for k, v in
    final.groupby('kdigo_stage')['unit_mortality'].mean().items()}
flow['baseline_source_dist'] = {k: int(v) for k, v in
                                final['baseline_source'].value_counts().items()}

with open(os.path.join(OUT, "eicu_cohort_v2_flow.json"), 'w') as f:
    json.dump(flow, f, indent=2)

print(f"\n  Saved: eicu_cohort_v2.csv ({len(final):,} rows)")
print(f"  Saved: eicu_cohort_v2_flow.json")
print(f"\nDone: {datetime.now()}")
