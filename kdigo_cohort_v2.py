#!/usr/bin/env python3
"""
KDIGO AKI 分期 v2 — 队列构建 + KDIGO 严格重实现 (F1+F2 修复)
================================================================
对照外部评审报告 (2026-09-23) 的阻断性问题修复:

  F1-a: KDIGO 严格时间窗 — Stage 1 的 +0.3 mg/dL 标准限制在任意 48h
        滚动窗内；比值标准 (1.5x/2x/3x) 限制在 ICU 入科后 7 天内
  F1-b: peak SCr 不再取住院全程最大值 — 只用 [intime, intime+7d] 内测量
  F1-c: SCr >= 4.0 mg/dL 的 Stage 3 绝对标准必须伴随急性变化
        (>=1.5x baseline 或 48h 内升高 >=0.3)，慢性高 SCr 不再直接进 Stage 3
  F2-a: 队列单位 = 患者 — 每位患者仅保留首次合格 ICU stay (subject_id 去重)
  F2-b: 无 SCr 患者明确排除并计数 (不再 MDRD 反推 + 强制 Stage 0)
  F2-c: 排除 ESRD/慢性透析 (N18.6/585.6/Z99.2/V45.1, 限定 index 入院及之前)

基线肌酐层级 (审稿建议优先序):
  1. 入院前 7-365 天 SCr 中位数 (优先门诊值, labevents hadm_id 为空)
  2. 入院前 0-7 天 SCr 中位数
  3. 入院后 24h 内首值
  → baseline_source 记录来源, 供敏感性分析

RRT: 仅 index ICU stay 入科后 7 天内的透析事件 (procedureevents, 按
     starttime 限定时间窗)，不再按整次住院或患者历史提升分期

CKD: 仅用 index 入院【之前】的住院诊断 (N18.3-5 / 585.3-5)，
     消除未来住院信息泄漏

产出: kdigo_cohort_v2.csv + cohort_v2_flow.json
"""
import pandas as pd
import numpy as np
import os, json, gc
from datetime import datetime

print("=" * 80)
print("KDIGO-AKI v2: Cohort Construction + Strict KDIGO Staging (F1+F2)")
print(f"Start: {datetime.now()}")
print("=" * 80)

BASE = "E:/mimic-iv"
HOSP = os.path.join(BASE, "v3.1/physionet.org/files/mimiciv/3.1/hosp")
ICU = os.path.join(BASE, "v3.1/physionet.org/files/mimiciv/3.1/icu")
OUT = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"

CREATININE_ITEMID = 50912          # serum creatinine (verified against d_labitems)
CR_MIN, CR_MAX = 0.1, 30.0         # physiologic plausibility bounds (mg/dL)
WINDOW_RATIO_DAYS = 7              # ratio criteria window from ICU admission
WINDOW_RISE_HOURS = 48             # rolling window for +0.3 mg/dL criterion
MIN_ICU_LOS_HOURS = 12             # minimum ICU stay (harmonised with eICU pipeline)
STAGE3_ABS_CUT = 4.0               # mg/dL absolute Stage 3 threshold

def read_gz(path, **kwargs):
    return pd.read_csv(path, compression='gzip', low_memory=False, **kwargs)

flow = {}

# ============================================================
# PHASE 1: 核心表 + 成人筛选 + 患者级去重 (F2-a)
# ============================================================
print("\n" + "=" * 60)
print("PHASE 1: Core tables, adults, first eligible ICU stay per PATIENT")
print("=" * 60)

icustays = read_gz(os.path.join(ICU, "icustays.csv.gz"))
patients = read_gz(os.path.join(HOSP, "patients.csv.gz"))
admissions = read_gz(os.path.join(HOSP, "admissions.csv.gz"),
                     usecols=['subject_id','hadm_id','admittime','dischtime',
                              'deathtime','hospital_expire_flag','admission_type',
                              'race','insurance','marital_status'])

for df, cols in [(icustays, ['intime','outtime']),
                 (admissions, ['admittime','dischtime','deathtime'])]:
    for c in cols:
        df[c] = pd.to_datetime(df[c])

flow['total_icu_stays'] = len(icustays)
flow['total_icu_patients'] = icustays['subject_id'].nunique()
print(f"  Total ICU stays: {flow['total_icu_stays']:,} "
      f"({flow['total_icu_patients']:,} patients)")

# Age at ICU admission (anchor_age; >89 encoded as 91 in MIMIC-IV)
icustays = icustays.merge(
    patients[['subject_id','anchor_age','anchor_year','gender']],
    on='subject_id', how='left')
icustays['age'] = icustays['anchor_age'] + (
    icustays['intime'].dt.year - icustays['anchor_year'])

adult = icustays[icustays['age'] >= 18].copy()
flow['after_adult'] = len(adult)
flow['excluded_age_lt18'] = flow['total_icu_stays'] - flow['after_adult']
print(f"  After age >= 18: {flow['after_adult']:,} stays "
      f"(excluded {flow['excluded_age_lt18']:,})")

# Valid admission link
adult = adult[adult['hadm_id'].notna()].copy()
# Minimum ICU LOS (harmonised with eICU >= 12h rule)
adult['icu_los_h'] = (adult['outtime'] - adult['intime']).dt.total_seconds() / 3600
eligible = adult[adult['icu_los_h'] >= MIN_ICU_LOS_HOURS].copy()
flow['after_los_12h'] = len(eligible)
flow['excluded_los_lt12h'] = len(adult) - len(eligible)
print(f"  After ICU LOS >= 12h: {flow['after_los_12h']:,} stays "
      f"(excluded {flow['excluded_los_lt12h']:,})")

# --- F2-a: first ELIGIBLE ICU stay per PATIENT (subject_id, not hadm_id) ---
eligible = eligible.sort_values(['subject_id', 'intime', 'stay_id'])
first_stay = eligible.drop_duplicates('subject_id', keep='first').copy()
flow['after_patient_dedup'] = len(first_stay)
flow['excluded_repeat_stays'] = len(eligible) - len(first_stay)
print(f"  First eligible ICU stay per patient: {flow['after_patient_dedup']:,} "
      f"(dropped {flow['excluded_repeat_stays']:,} later stays)")

# Link admission-level info
first_stay = first_stay.merge(
    admissions[['subject_id','hadm_id','admittime','dischtime','deathtime',
                'hospital_expire_flag','race','admission_type']],
    on=['subject_id','hadm_id'], how='left', validate='one_to_one')

# ============================================================
# PHASE 2: ESRD / 慢性透析排除 (F2-c) + CKD 病史 (仅既往)
# ============================================================
print("\n" + "=" * 60)
print("PHASE 2: ESRD/chronic dialysis exclusion + prior-admission CKD flag")
print("=" * 60)

diagnoses = read_gz(os.path.join(HOSP, "diagnoses_icd.csv.gz"),
                    usecols=['subject_id','hadm_id','icd_version','icd_code'])

# Admission times for every hadm of cohort subjects (to date-restrict diagnoses)
subj_hadms = admissions[admissions['subject_id'].isin(
    set(first_stay['subject_id']))][['subject_id','hadm_id','admittime']]
diag_dated = diagnoses.merge(subj_hadms, on=['subject_id','hadm_id'], how='inner')

def icd_match(df, pattern10, pattern9):
    # MIMIC-IV stores ICD codes WITHOUT the decimal point (e.g. 'N186' = N18.6)
    return df[
        ((df['icd_version'] == 10) & df['icd_code'].str.match(pattern10, na=False)) |
        ((df['icd_version'] == 9) & df['icd_code'].str.match(pattern9, na=False))
    ]

# ESRD / chronic dialysis: N18.6, Z99.2 / 585.6, V45.1  (index admission or earlier)
esrd_pat = icd_match(diag_dated, r'^(N186|Z992)', r'^(5856|V451)')

# restrict to diagnoses dated on/before the index admission
idx_adm = first_stay[['subject_id','admittime']].rename(columns={'admittime':'idx_admittime'})
esrd_pat = esrd_pat.merge(idx_adm, on='subject_id', how='left')
esrd_pre = set(esrd_pat[esrd_pat['admittime'] <= esrd_pat['idx_admittime']]['subject_id'])

cohort = first_stay[~first_stay['subject_id'].isin(esrd_pre)].copy()
flow['after_esrd_exclusion'] = len(cohort)
flow['excluded_esrd_dialysis'] = len(first_stay) - len(cohort)
print(f"  After ESRD/chronic-dialysis exclusion: {flow['after_esrd_exclusion']:,} "
      f"(excluded {flow['excluded_esrd_dialysis']:,})")

# CKD history: N18.3-N18.5 / 585.3-585.5 — STRICTLY BEFORE index admission
ckd_pat = icd_match(diag_dated, r'^N18[345]', r'^585[345]')
ckd_pat = ckd_pat.merge(idx_adm, on='subject_id', how='left')
ckd_prior = set(ckd_pat[ckd_pat['admittime'] < ckd_pat['idx_admittime']]['subject_id'])
ckd_incl_index = set(ckd_pat[ckd_pat['admittime'] <= ckd_pat['idx_admittime']]['subject_id'])

cohort['ckd_history_prior'] = cohort['subject_id'].isin(ckd_prior).astype(int)
cohort['ckd_incl_index'] = cohort['subject_id'].isin(ckd_incl_index).astype(int)
print(f"  CKD stage 3-5 history (prior admissions only): "
      f"{cohort['ckd_history_prior'].sum():,} ({cohort['ckd_history_prior'].mean()*100:.1f}%)")

del diag_dated, subj_hadms, esrd_pat, ckd_pat; gc.collect()

# ============================================================
# PHASE 3: labevents 全量扫描 — 队列患者的所有血清肌酐 (F1 数据基础)
# ============================================================
print("\n" + "=" * 60)
print("PHASE 3: labevents scan — all serum creatinine for cohort patients")
print("=" * 60)

# Verify creatinine itemid label
d_lab = read_gz(os.path.join(HOSP, "d_labitems.csv.gz"))
lab_label = d_lab[d_lab['itemid'] == CREATININE_ITEMID]['label'].iloc[0]
print(f"  itemid {CREATININE_ITEMID} = '{lab_label}'")
del d_lab

cohort_subjects = set(cohort['subject_id'])
labevents_path = os.path.join(HOSP, "labevents.csv.gz")

cr_parts = []
chunk_size = 5_000_000
n_chunks = 0
for chunk in pd.read_csv(labevents_path, compression='gzip', chunksize=chunk_size,
                         usecols=['subject_id','itemid','charttime','valuenum','hadm_id'],
                         dtype={'subject_id':'int64','itemid':'int64',
                                'valuenum':'float64'},
                         low_memory=False):
    n_chunks += 1
    m = chunk[(chunk['itemid'] == CREATININE_ITEMID)
              & (chunk['subject_id'].isin(cohort_subjects))]
    if len(m) > 0:
        cr_parts.append(m)
    if n_chunks % 5 == 0:
        print(f"    ... scanned {n_chunks*chunk_size/1e6:.0f}M rows")

cr_all = pd.concat(cr_parts, ignore_index=True)
del cr_parts, chunk; gc.collect()
cr_all['charttime'] = pd.to_datetime(cr_all['charttime'])
cr_all = cr_all[(cr_all['valuenum'] >= CR_MIN) & (cr_all['valuenum'] <= CR_MAX)]
cr_all = cr_all.sort_values(['subject_id','charttime'])
print(f"  Creatinine rows for cohort subjects: {len(cr_all):,} "
      f"across {cr_all['subject_id'].nunique():,} patients")

# ============================================================
# PHASE 4: 基线肌酐层级 (F1 — 参考入院时间, 非住院首值)
# ============================================================
print("\n" + "=" * 60)
print("PHASE 4: Baseline creatinine hierarchy")
print("=" * 60)

# Attach index admission time to each creatinine row
cr_idx = cr_all.merge(
    cohort[['subject_id','admittime','intime']].rename(columns={'admittime':'idx_adm','intime':'idx_icu'}),
    on='subject_id', how='inner')

days_before_adm = (cr_idx['idx_adm'] - cr_idx['charttime']).dt.total_seconds() / 86400

prior_outpatient = cr_idx[(days_before_adm >= 7) & (days_before_adm <= 365)
                          & cr_idx['hadm_id'].isna()]
prior_any        = cr_idx[(days_before_adm >= 7) & (days_before_adm <= 365)]
prior_0_7        = cr_idx[(days_before_adm >= 0) & (days_before_adm < 7)]
hosp_first24h    = cr_idx[(cr_idx['charttime'] >= cr_idx['idx_adm'])
                          & (cr_idx['charttime'] <= cr_idx['idx_adm'] + pd.Timedelta(hours=24))]

bl_out = prior_outpatient.groupby('subject_id')['valuenum'].median()
bl_any = prior_any.groupby('subject_id')['valuenum'].median()
bl_07  = prior_0_7.groupby('subject_id')['valuenum'].median()
bl_24  = hosp_first24h.sort_values('charttime').groupby('subject_id')['valuenum'].first()

cohort['baseline_cr'] = np.nan
cohort['baseline_source'] = None
for src, series in [('prior_outpt_7_365d', bl_out), ('prior_any_7_365d', bl_any),
                    ('prior_0_7d', bl_07), ('hosp_first24h', bl_24)]:
    need = cohort['baseline_cr'].isna()
    cohort.loc[need, 'baseline_cr'] = cohort.loc[need, 'subject_id'].map(series)
    cohort.loc[need & cohort['baseline_cr'].notna(), 'baseline_source'] = src

n_bl = cohort['baseline_cr'].notna().sum()
print(f"  Baseline available: {n_bl:,}/{len(cohort):,}")
print("  Baseline source distribution:")
for src, cnt in cohort['baseline_source'].value_counts(dropna=False).items():
    print(f"    {src}: {cnt:,} ({cnt/len(cohort)*100:.1f}%)")

flow['after_baseline_available'] = int(n_bl)
flow['excluded_no_baseline'] = int(len(cohort) - n_bl)

# ============================================================
# PHASE 5: ICU 入科后 7 天窗口内的肌酐 + 严格 KDIGO 分期 (F1)
# ============================================================
print("\n" + "=" * 60)
print("PHASE 5: Strict KDIGO staging in [intime, intime+7d]")
print("=" * 60)

win_end = cohort['intime'] + pd.Timedelta(days=WINDOW_RATIO_DAYS)
cr_win_mask = (cr_idx['charttime'] >= cr_idx['idx_icu']) & \
              (cr_idx['charttime'] <= cr_idx['idx_icu'] + pd.Timedelta(days=WINDOW_RATIO_DAYS))
cr_win = cr_idx[cr_win_mask][['subject_id','charttime','valuenum']].copy()
n_win_subj = cr_win['subject_id'].nunique()
print(f"  Patients with >=1 SCr in [intime, intime+7d]: {n_win_subj:,}")

flow['after_scr_in_window'] = int(n_win_subj)
flow['excluded_no_window_scr'] = int(n_bl - n_win_subj)

# Stageable cohort: baseline + in-window SCr
stageable_ids = cohort[cohort['baseline_cr'].notna()
                       & cohort['subject_id'].isin(set(cr_win['subject_id']))]
flow['final_cohort'] = len(stageable_ids)
print(f"  Final stageable cohort: {flow['final_cohort']:,} patients")

cr_win = cr_win[cr_win['subject_id'].isin(set(stageable_ids['subject_id']))]
cr_win = cr_win.sort_values(['subject_id','charttime'])

# Baseline lookup
def stage_patient(gr, intime, base):
    """Strict KDIGO staging for one patient. Returns Series of stage, onset_h,
    max_ratio, max_rise48, stage3_abs flag, stage1_ratio flag, n measurements."""
    t = gr['charttime'].values
    v = gr['valuenum'].values
    n = len(v)
    t0 = np.datetime64(intime)

    max_ratio = 0.0
    max_rise = 0.0
    stage_traj = np.zeros(n, dtype=int)   # stage achieved at each time point
    stage3_abs = False
    stage1_ratio = False

    for j in range(n):
        # ratio criteria (vs baseline)
        ratio = v[j] / base if base > 0 else 0.0
        if ratio > max_ratio:
            max_ratio = ratio
        if 1.5 <= ratio < 2.0:
            stage1_ratio = True
        s_ratio = 3 if ratio >= 3.0 else (2 if ratio >= 2.0 else (1 if ratio >= 1.5 else 0))

        # +0.3 mg/dL within rolling 48h window ending at t_j
        rise = 0.0
        i = j
        while i >= 0 and (t[j] - t[i]) <= np.timedelta64(WINDOW_RISE_HOURS, 'h'):
            d = v[j] - v[i]
            if d > rise:
                rise = d
            i -= 1
        if rise > max_rise:
            max_rise = rise
        s_rise = 1 if rise >= 0.3 else 0

        # absolute Stage 3 threshold requires acute change (F1-c)
        s_abs = 0
        if v[j] >= STAGE3_ABS_CUT and (s_ratio >= 1 or s_rise >= 1):
            s_abs = 3
            stage3_abs = True

        stage_traj[j] = max(s_ratio, s_rise, s_abs)

    stage = int(stage_traj.max()) if n > 0 else 0
    onset_h = np.nan
    if stage > 0:
        first_idx = int(np.argmax(stage_traj >= stage))
        onset_h = (t[first_idx] - t0) / np.timedelta64(1, 'h')

    return pd.Series({'kdigo_stage': stage,
                      'aki_onset_h': onset_h,
                      'max_cr_7d': v.max(),
                      'max_ratio_7d': max_ratio,
                      'max_rise48h': max_rise,
                      'stage3_abs_aki': int(stage3_abs),
                      'stage1_ratio_met': int(stage1_ratio),
                      'n_cr_window': n})

# Group-level staging: intime/baseline passed explicitly (no attrs dependency)
cr_win_grouped = cr_win.groupby('subject_id')
intime_map = cohort.set_index('subject_id')['intime']
bl_lookup = cohort.set_index('subject_id')['baseline_cr']

results = []
staged_ids = []
for sid, gr in cr_win_grouped:
    staged_ids.append(sid)
    results.append(stage_patient(gr, intime_map.loc[sid], bl_lookup.loc[sid]))
    if len(results) % 10000 == 0:
        print(f"    ... staged {len(results):,} patients")

staged = pd.DataFrame(results, index=pd.Index(staged_ids, name='subject_id'))
print(f"  Staged: {len(staged):,} patients")

# ============================================================
# PHASE 6: RRT — 仅 index stay 入科后 7 天内 (F1/F4)
# ============================================================
print("\n" + "=" * 60)
print("PHASE 6: RRT within first 7 ICU days (index stay only)")
print("=" * 60)

d_items = read_gz(os.path.join(ICU, "d_items.csv.gz"))
dial_items = d_items[d_items['label'].str.lower().str.contains(
    'dialysis - crrt|dialysis - cvvhd|dialysis - cvvhdf|hemodialysis|dialysis - scuf',
    na=False)]
dial_itemids = set(dial_items['itemid'].tolist())
print(f"  RRT itemids: {len(dial_itemids)}")
del d_items

proc = read_gz(os.path.join(ICU, "procedureevents.csv.gz"),
               usecols=['subject_id','hadm_id','stay_id','starttime','itemid'])
proc['starttime'] = pd.to_datetime(proc['starttime'])
proc = proc[proc['itemid'].isin(dial_itemids)]

proc = proc[proc['stay_id'].isin(set(cohort['stay_id']))]
proc = proc.merge(cohort[['stay_id','intime']], on='stay_id', how='left')
rrt_7d = proc[(proc['starttime'] >= proc['intime']) &
              (proc['starttime'] <= proc['intime'] + pd.Timedelta(days=WINDOW_RATIO_DAYS))]
rrt_ids = set(rrt_7d['stay_id'])
cohort['rrt_7d'] = cohort['stay_id'].isin(rrt_ids).astype(int)
print(f"  RRT within first 7 ICU days: {cohort['rrt_7d'].sum():,} patients")

# RRT upgrade to Stage 3 (time-restricted, index stay only)
cohort_v2 = cohort.merge(staged, left_on='subject_id', right_index=True, how='left')
rrt_mask = (cohort_v2['rrt_7d'] == 1) & (cohort_v2['kdigo_stage'] < 3)
cohort_v2.loc[rrt_mask, 'kdigo_stage'] = 3
cohort_v2['stage3_by_rrt'] = ((cohort_v2['rrt_7d'] == 1) &
                              (cohort_v2['stage3_abs_aki'] == 0) &
                              (cohort_v2['max_ratio_7d'] < 3.0)).astype(int)

# ============================================================
# PHASE 7: 汇总 + 保存
# ============================================================
print("\n" + "=" * 60)
print("PHASE 7: Summary + save")
print("=" * 60)

# Restrict to final stageable cohort
final = cohort_v2[cohort_v2['kdigo_stage'].notna()].copy()
final['kdigo_stage'] = final['kdigo_stage'].astype(int)

print("\n  FINAL COHORT (patient-level):")
print(f"  N = {len(final):,}")
print("\n  KDIGO stage distribution:")
mortality_by_stage = final.groupby('kdigo_stage').agg(
    n=('kdigo_stage','size'),
    hospital_death=('hospital_expire_flag','mean'))
for stage, row in mortality_by_stage.iterrows():
    print(f"    Stage {stage}: {row['n']:>7,} ({row['n']/len(final)*100:5.1f}%)  "
          f"hospital mortality {row['hospital_death']*100:5.1f}%")
print(f"  Overall hospital mortality: "
      f"{final['hospital_expire_flag'].mean()*100:.1f}%")
print(f"  AKI (stage >= 1): {(final['kdigo_stage']>=1).mean()*100:.1f}%")

# Save
out_cols = ['subject_id','hadm_id','stay_id','gender','age','race','admission_type',
            'intime','outtime','admittime','dischtime','deathtime',
            'hospital_expire_flag','icu_los_h',
            'baseline_cr','baseline_source','ckd_history_prior','ckd_incl_index',
            'n_cr_window','max_cr_7d','max_ratio_7d','max_rise48h',
            'kdigo_stage','aki_onset_h','stage3_abs_aki','stage1_ratio_met',
            'rrt_7d','stage3_by_rrt']
final[out_cols].to_csv(os.path.join(OUT, "kdigo_cohort_v2.csv"), index=False)

flow['final_cohort'] = len(final)
flow['kdigo_distribution'] = {int(k): int(v) for k, v in
                              final['kdigo_stage'].value_counts().sort_index().items()}
flow['overall_mortality'] = float(final['hospital_expire_flag'].mean())
flow['aki_prevalence'] = float((final['kdigo_stage'] >= 1).mean())
flow['baseline_source_dist'] = {k: int(v) for k, v in
                                final['baseline_source'].value_counts().items()}

with open(os.path.join(OUT, "cohort_v2_flow.json"), 'w') as f:
    json.dump(flow, f, indent=2)

print(f"\n  Saved: kdigo_cohort_v2.csv ({len(final):,} rows)")
print(f"  Saved: cohort_v2_flow.json")
print(f"\nDone: {datetime.now()}")
