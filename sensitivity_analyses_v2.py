#!/usr/bin/env python3
"""
增值敏感性分析 v2 (基于 v2 重分析流水线)
================================================
(1) 院前基线敏感性 — 仅纳入有入院前基线SCr的患者重跑 Model C
(2) eGFR 反推基线 — CKD-EPI 2021 race-free, eGFR=75 反推基线SCr,
    对无院前基线者插补后重新 KDIGO 分期 + Model C
(3) 亚组 + 交互检验 — age>=65 / 性别 / CKD / 糖尿病 / 通气 / SOFA三分位
    的 Stage 3 效应森林图 + 交互 p
(4) 一过性 vs 持续性 AKI 表型 — AKI 发生后 72h 内 SCr 是否回落
    (<1.5x 基线 且 < 基线+0.3 mg/dL) 分类, Model C + 30d 死亡率

前置: extract_scr_series_v2.py 已生成
  mimic_scr_series.csv.gz (subject_id, charttime, valuenum)
  eicu_scr_series.csv.gz  (patientunitstayid, labresultoffset, labresult)

预处理与 survival_mimic_v2.py / survival_eicu_v2.py 完全一致:
  MIMIC: event_time 下限 0.5h; 死亡事件封顶 dischtime+24h; landmark>=24h;
         time_from_lm 截断 696h; Model C penalizer=0.0
  eICU : event_time_h = hospitaldischargeoffset/60; landmark>=24h; 截断 696h

输出: v2_outputs/sensitivity_results.json
      v2_outputs/sensitivity_report.md
      v2_outputs/fig_subgroup_forest.png/pdf
      v2_outputs/fig_phenotype_forest.png/pdf
"""
import json
import os
import warnings
import pandas as pd
import numpy as np
from scipy.optimize import brentq
from lifelines import CoxPHFitter
warnings.filterwarnings('ignore')

WORK = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
OUTD = os.path.join(WORK, "v2_outputs")
os.makedirs(OUTD, exist_ok=True)

LANDMARK = 24
TRUNC_H = 24 * 30
WINDOW_RISE_H = 48.0
WINDOW_RATIO_H = 168.0      # 7 days
PHENO_RESOLVE_H = 72.0      # transient/persistent cutoff
STAGE3_ABS_CUT = 4.0

MODEL_C = ['kdigo_1', 'kdigo_2', 'kdigo_3', 'age', 'sex_male',
           'sofa_nonrenal', 'vent_24h', 'ckd_history_prior', 'diabetes',
           'hypertension', 'heart_failure', 'copd', 'liver_disease']

results = {'generated': str(pd.Timestamp.now())}

print("=" * 70)
print("增值敏感性分析 (v2 pipeline)")
print("=" * 70)

# ============================================================
# 0. 数据加载 + 与生存脚本一致的预处理
# ============================================================
print("\n[0] Load + landmark preprocessing (identical to survival scripts)")

# ---- MIMIC ----
mim = pd.read_csv(os.path.join(WORK, "mimic_analysis_v2.csv"))
mim['event_time_h'] = mim['event_time_h'].clip(lower=0.5)
coh_m = pd.read_csv(os.path.join(WORK, "kdigo_cohort_v2.csv"),
                    usecols=['stay_id', 'subject_id', 'intime', 'dischtime',
                             'deathtime', 'hospital_expire_flag', 'rrt_7d',
                             'baseline_cr', 'baseline_source'])
coh_m['intime'] = pd.to_datetime(coh_m['intime'])
coh_m['dischtime'] = pd.to_datetime(coh_m['dischtime'])
mim = mim.merge(coh_m, on=['stay_id', 'hospital_expire_flag'], how='left',
                suffixes=('', '_coh'))
mim['event_time_h'] = np.where(
    mim['hospital_expire_flag'] == 1,
    np.minimum(mim['event_time_h'],
               (mim['dischtime'] - mim['intime']).dt.total_seconds() / 3600 + 24),
    mim['event_time_h'])
mim.loc[mim['event_time_h'] < 0, 'event_time_h'] = 0.5
for s in [1, 2, 3]:
    mim[f'kdigo_{s}'] = (mim['kdigo_stage'] == s).astype(int)
msur = mim[mim['event_time_h'] >= LANDMARK].copy()
msur['time_from_lm'] = (msur['event_time_h'] - LANDMARK).clip(upper=TRUNC_H - LANDMARK)
print(f"  MIMIC landmark cohort: {len(msur):,} (events {msur['event_ind'].sum():,})")

# ---- eICU ----
eic = pd.read_csv(os.path.join(WORK, "eicu_analysis_v2.csv"))
for s in [1, 2, 3]:
    eic[f'kdigo_{s}'] = (eic['kdigo_stage'] == s).astype(int)
esur = eic[eic['event_time_h'] >= LANDMARK].copy()
esur['time_from_lm'] = (esur['event_time_h'] - LANDMARK).clip(upper=TRUNC_H - LANDMARK)
print(f"  eICU landmark cohort: {len(esur):,} (events {esur['event_ind'].sum():,})")


def fit_model_c(df, extra_covars=None, label=''):
    """Model C fit; returns dict of n, C, variables(kdigo only + extras)."""
    covars = MODEL_C + (extra_covars or [])
    keep = [v for v in covars if v in df.columns]
    d = df[['time_from_lm', 'event_ind'] + keep].dropna()
    cph = CoxPHFitter(penalizer=0.0)
    cph.fit(d, duration_col='time_from_lm', event_col='event_ind')
    out = {'n': int(len(d)), 'C': float(cph.concordance_index_), 'variables': {}}
    for v in keep:
        r = cph.summary.loc[v]
        out['variables'][v] = {
            'HR': float(np.exp(r['coef'])),
            'CI_lo': float(r['exp(coef) lower 95%']),
            'CI_hi': float(r['exp(coef) upper 95%']),
            'p': float(r['p'])}
    if label:
        print(f"  [{label}] n={out['n']:,} C={out['C']:.3f} | "
              + " ".join(f"S{s}={out['variables'][f'kdigo_{s}']['HR']:.2f}"
                         for s in [1, 2, 3]))
    return out

# ============================================================
# (1) 院前基线敏感性
# ============================================================
print("\n[1] Sensitivity: pre-admission baseline only")
MIMIC_PRE = ['prior_outpt_7_365d', 'prior_any_7_365d', 'prior_0_7d']
EICU_PRE = ['preICU_0_7d']

m_pre = msur[msur['baseline_source'].isin(MIMIC_PRE)]
e_pre = esur[esur['baseline_source'].isin(EICU_PRE)]
results['s1_pre_admission_baseline'] = {
    'mimic_n': int(len(m_pre)),
    'mimic_pct': float(len(m_pre) / len(msur)),
    'eicu_n': int(len(e_pre)),
    'eicu_pct': float(len(e_pre) / len(esur)),
    'mimic_model_C': fit_model_c(m_pre, label='MIMIC pre-adm'),
    'eicu_model_C': fit_model_c(e_pre, label='eICU pre-adm')}

# ============================================================
# (3) 亚组 + 交互检验 (Stage 3 效应)
# ============================================================
print("\n[3] Subgroup analyses + interaction tests (Stage 3 HR)")


def subgroup_analysis(df, cohort_name):
    d0 = df.copy()
    d0['sofa_tert'] = pd.qcut(d0['sofa_nonrenal'], 3,
                              labels=['T1', 'T2', 'T3']).astype(str)
    d0['age_ge65'] = (d0['age'] >= 65).astype(int)
    d0['sofa_t3'] = (d0['sofa_tert'] == 'T3').astype(int)
    subs = [
        ('Age < 65 vs >=65', 'age_ge65'),
        ('Female vs Male', 'sex_male'),
        ('No prior CKD vs CKD', 'ckd_history_prior'),
        ('No diabetes vs diabetes', 'diabetes'),
        ('Not ventilated vs ventilated', 'vent_24h'),
        ('SOFA non-renal T1-2 vs T3', 'sofa_t3'),
    ]
    rows = []
    for label, col in subs:
        d = d0.copy()
        d['sub'] = d[col]
        d['k3_sub'] = d['kdigo_3'] * d['sub']
        # if the subgroup variable is itself a Model C covariate, do NOT add
        # a duplicate 'sub' column (perfect collinearity)
        extra = ['k3_sub'] if col in MODEL_C else ['sub', 'k3_sub']
        covars = MODEL_C + extra
        dd = d[['time_from_lm', 'event_ind'] + covars].dropna()
        try:
            cph = CoxPHFitter(penalizer=0.0)
            cph.fit(dd, duration_col='time_from_lm', event_col='event_ind')
            b3, b3se = cph.params_['kdigo_3'], cph.standard_errors_['kdigo_3']
            bi, bise = cph.params_['k3_sub'], cph.standard_errors_['k3_sub']
            cov = cph.variance_matrix_.loc['kdigo_3', 'k3_sub']
            def hr_ci(b, se):
                return (float(np.exp(b)),
                        float(np.exp(b - 1.96 * se)),
                        float(np.exp(b + 1.96 * se)))
            hr0 = hr_ci(b3, b3se)
            hr1 = hr_ci(b3 + bi, np.sqrt(b3se**2 + bise**2 + 2 * cov))
            z = bi / bise
            from scipy.stats import norm
            p_int = float(2 * norm.sf(abs(z)))
            n0 = int((d.loc[dd.index, col] == 0).sum())
            n1 = int((d.loc[dd.index, col] == 1).sum())
            rows.append({'subgroup': label, 'n_ref': n0, 'n_sub': n1,
                         'hr_ref': hr0[0], 'hr_ref_lo': hr0[1], 'hr_ref_hi': hr0[2],
                         'hr_sub': hr1[0], 'hr_sub_lo': hr1[1], 'hr_sub_hi': hr1[2],
                         'p_interaction': p_int})
        except Exception as ex:
            rows.append({'subgroup': label, 'error': str(ex)})
    print(f"  [{cohort_name}]")
    for r in rows:
        if 'error' not in r:
            print(f"    {r['subgroup']:>34s}: ref HR={r['hr_ref']:.2f} "
                  f"sub HR={r['hr_sub']:.2f} p_int={r['p_interaction']:.3g}")
    return rows


results['s3_subgroups'] = {
    'mimic': subgroup_analysis(msur, 'MIMIC'),
    'eicu': subgroup_analysis(esur, 'eICU')}

# ============================================================
# 共享: SCr 序列加载 + 分期函数
# ============================================================
print("\n[load] SCr series caches")


def stage_series(t, v, base):
    """Strict KDIGO staging (identical to v2 cohort scripts), t in hours."""
    n = len(v)
    stage_traj = np.zeros(n, dtype=int)
    for j in range(n):
        ratio = v[j] / base if base > 0 else 0.0
        s_ratio = 3 if ratio >= 3.0 else (2 if ratio >= 2.0 else (1 if ratio >= 1.5 else 0))
        rise = 0.0
        i = j
        while i >= 0 and (t[j] - t[i]) <= WINDOW_RISE_H:
            d = v[j] - v[i]
            if d > rise:
                rise = d
            i -= 1
        s_rise = 1 if rise >= 0.3 else 0
        s_abs = 3 if (v[j] >= STAGE3_ABS_CUT and (s_ratio >= 1 or s_rise >= 1)) else 0
        stage_traj[j] = max(s_ratio, s_rise, s_abs)
    stage = int(stage_traj.max()) if n > 0 else 0
    if (stage_traj >= 1).any():
        oi = int(np.argmax(stage_traj >= 1))
        onset_h = float(t[oi])
    else:
        onset_h = np.nan
    return stage, stage_traj, onset_h


def classify_phenotype(t, v, base, stage_traj, rrt, event_ind, event_time_h):
    """Transient vs persistent AKI. Resolution = SCr <1.5x baseline AND
    < baseline+0.3 within 72h of AKI onset. RRT or death before resolution
    = persistent. No post-onset SCr and alive = indeterminate."""
    if rrt == 1:
        return 'persistent', np.nan
    if not (stage_traj >= 1).any():
        return 'unclassified', np.nan
    oi = int(np.argmax(stage_traj >= 1))
    onset_h = float(t[oi])
    resolved = False
    n_after = 0
    for j in range(oi + 1, len(v)):
        if t[j] - onset_h > PHENO_RESOLVE_H:
            break
        n_after += 1
        if v[j] / base < 1.5 and v[j] < base + 0.3:
            resolved = True
            break
    if resolved:
        return 'transient', onset_h
    if n_after == 0:
        # no SCr after onset within 72h: death before resolution -> persistent
        if event_ind == 1 and (event_time_h - onset_h) <= PHENO_RESOLVE_H:
            return 'persistent', onset_h
        return 'indeterminate', onset_h
    return 'persistent', onset_h


# ============================================================
# (2) eGFR 反推基线 (CKD-EPI 2021 race-free, eGFR = 75)
# ============================================================
print("\n[2] Sensitivity: eGFR-imputed baseline (CKD-EPI 2021, eGFR=75)")


def ckdepi2021(scr, age, female):
    k = 0.7 if female else 0.9
    a = -0.241 if female else -0.302
    e = 142.0 * min(scr / k, 1.0) ** a * max(scr / k, 1.0) ** (-1.200) \
        * 0.9938 ** age
    if female:
        e *= 1.012
    return e


def scr_from_egfr75(age, female):
    return brentq(lambda s: ckdepi2021(s, age, female) - 75.0, 0.05, 25.0)


# ---- MIMIC restaging ----
ser_m = pd.read_csv(os.path.join(WORK, "mimic_scr_series.csv.gz"))
ser_m['charttime'] = pd.to_datetime(ser_m['charttime'])
ser_m = ser_m.merge(coh_m[['subject_id', 'intime']], on='subject_id', how='inner')
ser_m['t_h'] = (ser_m['charttime'] - ser_m['intime']).dt.total_seconds() / 3600.0
ser_m = ser_m[(ser_m['t_h'] >= 0) & (ser_m['t_h'] <= WINDOW_RATIO_H)]
rrt_m = coh_m.set_index('subject_id')['rrt_7d']

mim2 = mim.copy()
need_m = mim2['baseline_source'] == 'hosp_first24h'
print(f"  MIMIC restage subset (hosp_first24h baseline): {need_m.sum():,}")
mim2['baseline_cr_egfr'] = np.where(
    need_m,
    [scr_from_egfr75(a, f == 0) for a, f in zip(mim2['age'], mim2['sex_male'])],
    mim2['baseline_cr'])
print(f"  eGFR-imputed baseline SCr: mean="
      f"{mim2.loc[need_m, 'baseline_cr_egfr'].mean():.2f} mg/dL")

new_stage = {}
grp = ser_m[ser_m['subject_id'].isin(set(mim2.loc[need_m, 'subject_id']))] \
    .groupby('subject_id')
bl_map = mim2.set_index('subject_id')['baseline_cr_egfr']
for sid, g in grp:
    t = g['t_h'].values; v = g['valuenum'].values
    st, _, _ = stage_series(t, v, float(bl_map.loc[sid]))
    if rrt_m.loc[sid] == 1 and st < 3:
        st = 3
    new_stage[sid] = st
mim2['kdigo_stage_egfr'] = mim2['kdigo_stage']
mim2.loc[need_m, 'kdigo_stage_egfr'] = mim2.loc[need_m, 'subject_id'].map(new_stage)
mim2.loc[mim2['kdigo_stage_egfr'].isna(), 'kdigo_stage_egfr'] = \
    mim2.loc[mim2['kdigo_stage_egfr'].isna(), 'kdigo_stage']

msur2 = mim2[mim2['event_time_h'] >= LANDMARK].copy()
msur2['time_from_lm'] = (msur2['event_time_h'] - LANDMARK).clip(upper=TRUNC_H - LANDMARK)
for s in [1, 2, 3]:
    msur2[f'kdigo_{s}'] = (msur2['kdigo_stage_egfr'] == s).astype(int)
print(f"  MIMIC eGFR-baseline AKI prevalence: "
      f"{(mim2['kdigo_stage_egfr'] >= 1).mean() * 100:.1f}% "
      f"(main: {(mim2['kdigo_stage'] >= 1).mean() * 100:.1f}%)")
results['s2_egfr_baseline'] = {
    'mimic_n_restage': int(need_m.sum()),
    'mimic_imputed_baseline_mean': float(mim2.loc[need_m, 'baseline_cr_egfr'].mean()),
    'mimic_aki_prevalence_main': float((mim2['kdigo_stage'] >= 1).mean()),
    'mimic_aki_prevalence_egfr': float((mim2['kdigo_stage_egfr'] >= 1).mean()),
    'mimic_model_C': fit_model_c(msur2, label='MIMIC eGFR')}
del ser_m, grp; import gc; gc.collect()

# ---- eICU restaging ----
ser_e = pd.read_csv(os.path.join(WORK, "eicu_scr_series.csv.gz"))
ser_e['t_h'] = ser_e['labresultoffset'].astype(float) / 60.0
ser_e = ser_e[(ser_e['t_h'] >= 0) & (ser_e['t_h'] <= WINDOW_RATIO_H)]
coh_e = pd.read_csv(os.path.join(WORK, "eicu_cohort_v2.csv"),
                    usecols=['patientunitstayid', 'rrt_7d'])
rrt_e = coh_e.set_index('patientunitstayid')['rrt_7d']

eic2 = eic.copy()
need_e = eic2['baseline_source'] == 'icu_first24h'
print(f"  eICU restage subset (icu_first24h baseline): {need_e.sum():,}")
eic2['baseline_cr_egfr'] = np.where(
    need_e,
    [scr_from_egfr75(a, f == 0) for a, f in zip(eic2['age'], eic2['sex_male'])],
    eic2['baseline_cr'])
print(f"  eGFR-imputed baseline SCr: mean="
      f"{eic2.loc[need_e, 'baseline_cr_egfr'].mean():.2f} mg/dL")

new_stage_e = {}
grp_e = ser_e[ser_e['patientunitstayid'].isin(set(eic2.loc[need_e, 'stay_id']))] \
    .groupby('patientunitstayid')
bl_map_e = eic2.set_index('stay_id')['baseline_cr_egfr']
for sid, g in grp_e:
    t = g['t_h'].values; v = g['labresult'].values
    st, _, _ = stage_series(t, v, float(bl_map_e.loc[sid]))
    if rrt_e.loc[sid] == 1 and st < 3:
        st = 3
    new_stage_e[sid] = st
eic2['kdigo_stage_egfr'] = eic2['kdigo_stage']
eic2.loc[need_e, 'kdigo_stage_egfr'] = eic2.loc[need_e, 'stay_id'].map(new_stage_e)
eic2.loc[eic2['kdigo_stage_egfr'].isna(), 'kdigo_stage_egfr'] = \
    eic2.loc[eic2['kdigo_stage_egfr'].isna(), 'kdigo_stage']

esur2 = eic2[eic2['event_time_h'] >= LANDMARK].copy()
esur2['time_from_lm'] = (esur2['event_time_h'] - LANDMARK).clip(upper=TRUNC_H - LANDMARK)
for s in [1, 2, 3]:
    esur2[f'kdigo_{s}'] = (esur2['kdigo_stage_egfr'] == s).astype(int)
print(f"  eICU eGFR-baseline AKI prevalence: "
      f"{(eic2['kdigo_stage_egfr'] >= 1).mean() * 100:.1f}% "
      f"(main: {(eic2['kdigo_stage'] >= 1).mean() * 100:.1f}%)")
results['s2_egfr_baseline'].update({
    'eicu_n_restage': int(need_e.sum()),
    'eicu_imputed_baseline_mean': float(eic2.loc[need_e, 'baseline_cr_egfr'].mean()),
    'eicu_aki_prevalence_main': float((eic2['kdigo_stage'] >= 1).mean()),
    'eicu_aki_prevalence_egfr': float((eic2['kdigo_stage_egfr'] >= 1).mean()),
    'eicu_model_C': fit_model_c(esur2, label='eICU eGFR')})

# ============================================================
# (4) 一过性 vs 持续性 AKI 表型
# ============================================================
print("\n[4] AKI phenotype: transient vs persistent (72h resolution)")

# ---- MIMIC ----
ser_m4 = pd.read_csv(os.path.join(WORK, "mimic_scr_series.csv.gz"))
ser_m4['charttime'] = pd.to_datetime(ser_m4['charttime'])
ser_m4 = ser_m4.merge(coh_m[['subject_id', 'intime', 'baseline_cr']],
                      on='subject_id', how='inner')
ser_m4['t_h'] = (ser_m4['charttime'] - ser_m4['intime']).dt.total_seconds() / 3600.0
ser_m4 = ser_m4[(ser_m4['t_h'] >= 0) & (ser_m4['t_h'] <= WINDOW_RATIO_H)]

aki_m = msur[msur['kdigo_stage'] >= 1].copy()
bl4 = coh_m.set_index('subject_id')['baseline_cr']
pheno_m = {}
grp4 = ser_m4[ser_m4['subject_id'].isin(set(aki_m['subject_id']))].groupby('subject_id')
subj_ev = aki_m.set_index('subject_id')[['event_ind', 'event_time_h']]
for sid, g in grp4:
    t = g['t_h'].values; v = g['valuenum'].values
    base = float(bl4.loc[sid])
    st, traj, onset = stage_series(t, v, base)
    ev = subj_ev.loc[sid]
    ph, onset = classify_phenotype(
        t, v, base, traj, int(rrt_m.loc[sid]),
        int(ev['event_ind']), float(ev['event_time_h']))
    pheno_m[sid] = (ph, onset, st)
aki_m['phenotype'] = aki_m['subject_id'].map(lambda s: pheno_m.get(s, ('unclassified',))[0])
aki_m['onset_h2'] = aki_m['subject_id'].map(lambda s: pheno_m.get(s, (None, np.nan))[1])

print("  MIMIC phenotype distribution:")
for ph, cnt in aki_m['phenotype'].value_counts().items():
    print(f"    {ph}: {cnt:,}")

# 30d binary mortality by phenotype (landmark cohort, no-AKI as ref)
m_noaki = msur[msur['kdigo_stage'] == 0].copy(); m_noaki['phenotype'] = 'none'
m_ph = pd.concat([m_noaki[['event_ind', 'event_time_h', 'phenotype']],
                  aki_m[['event_ind', 'event_time_h', 'phenotype']]])
m_ph['mort30'] = ((m_ph['event_ind'] == 1) &
                  (m_ph['event_time_h'] <= TRUNC_H)).astype(int)
mort_m = m_ph.groupby('phenotype')['mort30'].agg(['mean', 'size'])
print("  MIMIC 30d mortality by phenotype:")
print(mort_m.round(3).to_string())

# Cox: no AKI (ref) vs transient vs persistent, + Model C covariates
mph = msur.copy()
mph['phenotype'] = 'none'
mph.loc[mph['kdigo_stage'] >= 1, 'phenotype'] = \
    mph.loc[mph['kdigo_stage'] >= 1, 'subject_id'].map(
        lambda s: pheno_m.get(s, ('unclassified',))[0])
mph = mph[mph['phenotype'].isin(['none', 'transient', 'persistent'])].copy()
mph['aki_transient'] = (mph['phenotype'] == 'transient').astype(int)
mph['aki_persistent'] = (mph['phenotype'] == 'persistent').astype(int)
cov4 = ['aki_transient', 'aki_persistent'] + \
    [c for c in MODEL_C if not c.startswith('kdigo_')]
d = mph[['time_from_lm', 'event_ind'] + cov4].dropna()
cph = CoxPHFitter(penalizer=0.0)
cph.fit(d, duration_col='time_from_lm', event_col='event_ind')
res_ph_m = {'n': int(len(d)), 'C': float(cph.concordance_index_), 'variables': {}}
for v in cov4:
    r = cph.summary.loc[v]
    res_ph_m['variables'][v] = {
        'HR': float(np.exp(r['coef'])),
        'CI_lo': float(r['exp(coef) lower 95%']),
        'CI_hi': float(r['exp(coef) upper 95%']),
        'p': float(r['p'])}
print(f"  [MIMIC pheno] n={res_ph_m['n']:,} C={res_ph_m['C']:.3f} | "
      f"transient HR={res_ph_m['variables']['aki_transient']['HR']:.2f} "
      f"persistent HR={res_ph_m['variables']['aki_persistent']['HR']:.2f}")

# ---- eICU ----
ser_e4 = pd.read_csv(os.path.join(WORK, "eicu_scr_series.csv.gz"))
ser_e4['t_h'] = ser_e4['labresultoffset'].astype(float) / 60.0
ser_e4 = ser_e4[(ser_e4['t_h'] >= 0) & (ser_e4['t_h'] <= WINDOW_RATIO_H)]
ser_e4 = ser_e4.merge(eic[['stay_id', 'baseline_cr']].rename(
    columns={'stay_id': 'patientunitstayid'}), on='patientunitstayid', how='inner')

aki_e = esur[esur['kdigo_stage'] >= 1].copy()
bl4e = eic.set_index('stay_id')['baseline_cr']
pheno_e = {}
grp4e = ser_e4[ser_e4['patientunitstayid'].isin(set(aki_e['stay_id']))] \
    .groupby('patientunitstayid')
subj_ev_e = aki_e.set_index('stay_id')[['event_ind', 'event_time_h']]
for sid, g in grp4e:
    t = g['t_h'].values; v = g['labresult'].values
    base = float(bl4e.loc[sid])
    st, traj, onset = stage_series(t, v, base)
    ev = subj_ev_e.loc[sid]
    ph, onset = classify_phenotype(
        t, v, base, traj, int(rrt_e.loc[sid]),
        int(ev['event_ind']), float(ev['event_time_h']))
    pheno_e[sid] = (ph, onset, st)
aki_e['phenotype'] = aki_e['stay_id'].map(lambda s: pheno_e.get(s, ('unclassified',))[0])
print("  eICU phenotype distribution:")
for ph, cnt in aki_e['phenotype'].value_counts().items():
    print(f"    {ph}: {cnt:,}")

e_noaki = esur[esur['kdigo_stage'] == 0].copy(); e_noaki['phenotype'] = 'none'
e_ph = pd.concat([e_noaki[['event_ind', 'event_time_h', 'phenotype']],
                  aki_e[['event_ind', 'event_time_h', 'phenotype']]])
e_ph['mort30'] = ((e_ph['event_ind'] == 1) &
                  (e_ph['event_time_h'] <= TRUNC_H)).astype(int)
mort_e = e_ph.groupby('phenotype')['mort30'].agg(['mean', 'size'])
print("  eICU 30d mortality by phenotype:")
print(mort_e.round(3).to_string())

eph = esur.copy()
eph['phenotype'] = 'none'
eph.loc[eph['kdigo_stage'] >= 1, 'phenotype'] = \
    eph.loc[eph['kdigo_stage'] >= 1, 'stay_id'].map(
        lambda s: pheno_e.get(s, ('unclassified',))[0])
eph = eph[eph['phenotype'].isin(['none', 'transient', 'persistent'])].copy()
eph['aki_transient'] = (eph['phenotype'] == 'transient').astype(int)
eph['aki_persistent'] = (eph['phenotype'] == 'persistent').astype(int)
d = eph[['time_from_lm', 'event_ind'] + cov4].dropna()
cph = CoxPHFitter(penalizer=0.0)
cph.fit(d, duration_col='time_from_lm', event_col='event_ind')
res_ph_e = {'n': int(len(d)), 'C': float(cph.concordance_index_), 'variables': {}}
for v in cov4:
    r = cph.summary.loc[v]
    res_ph_e['variables'][v] = {
        'HR': float(np.exp(r['coef'])),
        'CI_lo': float(r['exp(coef) lower 95%']),
        'CI_hi': float(r['exp(coef) upper 95%']),
        'p': float(r['p'])}
print(f"  [eICU pheno] n={res_ph_e['n']:,} C={res_ph_e['C']:.3f} | "
      f"transient HR={res_ph_e['variables']['aki_transient']['HR']:.2f} "
      f"persistent HR={res_ph_e['variables']['aki_persistent']['HR']:.2f}")

results['s4_phenotype'] = {
    'mimic_dist': {k: int(v) for k, v in aki_m['phenotype'].value_counts().items()},
    'eicu_dist': {k: int(v) for k, v in aki_e['phenotype'].value_counts().items()},
    'mimic_mort30': {k: float(v) for k, v in mort_m['mean'].items()},
    'mimic_n': {k: int(v) for k, v in mort_m['size'].items()},
    'eicu_mort30': {k: float(v) for k, v in mort_e['mean'].items()},
    'eicu_n': {k: int(v) for k, v in mort_e['size'].items()},
    'mimic_model_C': res_ph_m,
    'eicu_model_C': res_ph_e,
    'definition': ('transient = SCr falls below 1.5x baseline AND below '
                   'baseline+0.3 mg/dL within 72h of AKI onset; RRT within '
                   '7d or death before resolution = persistent; indeterminate '
                   '(no post-onset SCr, alive) excluded')}

# ============================================================
# 保存 JSON
# ============================================================
with open(os.path.join(OUTD, "sensitivity_results.json"), 'w') as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nSaved: {OUTD}/sensitivity_results.json")

# ============================================================
# 图: 亚组森林图 + 表型森林图
# ============================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def forest(ax, rows, title):
    labels, y = [], []
    for i, r in enumerate(rows):
        if 'error' in r:
            continue
        labels.append(r['subgroup'])
        y.append(len(y))
        for hr, lo, hi, off, color in [
                (r['hr_ref'], r['hr_ref_lo'], r['hr_ref_hi'], -0.18, '#1f77b4'),
                (r['hr_sub'], r['hr_sub_lo'], r['hr_sub_hi'], +0.18, '#d62728')]:
            ax.errorbar(hr, len(y) - 1 + off,
                        xerr=[[hr - lo], [hi - hr]], fmt='o',
                        color=color, capsize=2.5, ms=4)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(1.0, color='gray', lw=0.8, ls='--')
    ax.set_xlabel('Stage 3 HR (95% CI)', fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.invert_yaxis()

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharex=True)
forest(axes[0], results['s3_subgroups']['mimic'], 'MIMIC-IV (derivation)')
forest(axes[1], results['s3_subgroups']['eicu'], 'eICU-CRD (validation)')
fig.suptitle('KDIGO Stage 3 vs Stage 0: subgroup hazard ratios '
             '(blue = reference level, red = subgroup)', fontsize=10)
fig.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(OUTD, f"fig_subgroup_forest.{ext}"), dpi=300)
plt.close(fig)
print(f"Saved: fig_subgroup_forest.png/pdf")

# 表型森林图 (transient / persistent, 两队列)
fig, ax = plt.subplots(figsize=(6.5, 3))
rows_plt = []
for cname, res_c in [('MIMIC-IV', res_ph_m), ('eICU-CRD', res_ph_e)]:
    for v, lab in [('aki_transient', 'Transient AKI'), ('aki_persistent', 'Persistent AKI')]:
        r = res_c['variables'][v]
        rows_plt.append((f"{cname} — {lab}", r['HR'], r['CI_lo'], r['CI_hi']))
for i, (lab, hr, lo, hi) in enumerate(rows_plt):
    ax.errorbar(hr, i, xerr=[[hr - lo], [hi - hr]], fmt='s', capsize=3, ms=5,
                color='#d62728' if 'Persistent' in lab else '#1f77b4')
ax.set_yticks(range(len(rows_plt)))
ax.set_yticklabels([r[0] for r in rows_plt], fontsize=9)
ax.axvline(1.0, color='gray', lw=0.8, ls='--')
ax.set_xlabel('Adjusted HR vs no AKI (95% CI)', fontsize=9)
ax.set_title('30-day mortality by AKI phenotype (72h resolution criterion)', fontsize=10)
ax.invert_yaxis()
fig.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(OUTD, f"fig_phenotype_forest.{ext}"), dpi=300)
plt.close(fig)
print(f"Saved: fig_phenotype_forest.png/pdf")

print("\nAll sensitivity analyses done.")
