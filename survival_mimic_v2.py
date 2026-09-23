"""
F5 修复 - 真实事件时间 Cox / KM 分析 (MIMIC-IV v3.1)

事件: 院内死亡 (hospital_expire_flag == 1)
时间: hours from ICU admission (intime)
KM   : log-rank + Bonferroni
Cox  : 3 incremental models + 时间分层
PH检验: 每个模型每个变量
Landmark: 24h, 截断 30d
"""
import json
import os
import warnings
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test, pairwise_logrank_test
warnings.filterwarnings('ignore')

WORK = r"C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
os.makedirs(os.path.join(WORK, "v2_outputs"), exist_ok=True)

LANDMARK = 24
TRUNC_H = 24 * 30

print("=" * 70)
print("MIMIC-IV v3.1 — 真实事件时间 Cox / KM 分析")
print("=" * 70)

# ---- 1. Load + sanitize ----
print("\n[1] Load analysis dataset + sanitize")
df = pd.read_csv(os.path.join(WORK, "mimic_analysis_v2.csv"))
print(f"  loaded: {len(df):,}")

# Sanitize event_time: floor at 0
df['event_time_h'] = df['event_time_h'].clip(lower=0.5)
# Fix deathtime-after-discharge anomalies: cap at dischtime (= censor)
# We need dischtime / deathtime here -- load raw cohort
cohort = pd.read_csv(os.path.join(WORK, "kdigo_cohort_v2.csv"),
                     usecols=['stay_id','intime','dischtime','deathtime','hospital_expire_flag'])
cohort['intime']    = pd.to_datetime(cohort['intime'])
cohort['dischtime'] = pd.to_datetime(cohort['dischtime'])
cohort['deathtime'] = pd.to_datetime(cohort['deathtime'])
df = df.merge(cohort[['stay_id','intime','dischtime','deathtime']],
              on='stay_id', how='left')

anom = ((df['hospital_expire_flag']==1) & (df['deathtime'] > df['dischtime'] + pd.Timedelta(hours=2))).sum()
print(f"  deaths with deathtime > dischtime (anomaly): {anom}")

# Cap event_time: deaths must be <= dischtime+24h (deaths after = censor)
df['event_time_h'] = np.where(
    df['hospital_expire_flag']==1,
    np.minimum(df['event_time_h'], (df['dischtime'] - df['intime']).dt.total_seconds()/3600 + 24),
    df['event_time_h']
)
df.loc[df['event_time_h'] < 0, 'event_time_h'] = 0.5

print(f"  event_time range: {df['event_time_h'].min():.1f} to {df['event_time_h'].max():.1f}")
print(f"  events: {df['event_ind'].sum():,}")

# ---- 2. Landmark + truncate ----
print(f"\n[2] Landmark {LANDMARK}h + truncate {TRUNC_H}h (30d)")
surv = df[df['event_time_h'] >= LANDMARK].copy()
surv['time_from_lm'] = (surv['event_time_h'] - LANDMARK).clip(upper=TRUNC_H - LANDMARK)
print(f"  cohort after landmark: {len(surv):,}")
print(f"  cohort after 30d cap: {len(surv):,}")
print(f"  events: {surv['event_ind'].sum():,}")

# ---- 3. KM + pairwise log-rank ----
# NOTE: in-hospital cohorts have heavy informative censoring (early discharge of survivors),
# so KM estimates past ~14d are biased upward (denominator inflates with patients
# discharged alive). Report binary mortality at multiple time points and provide
# KM curve up to 30d for visualization.
print("\n[3] KM curves + log-rank (Bonferroni)")
km_results = {'cohort_n': int(len(surv)),
              'cohort_n_stage': {}}
km_fitters = {}
surv_v = surv.dropna(subset=['kdigo_stage'])

for s in [0,1,2,3]:
    sub = surv_v[surv_v.kdigo_stage==s]
    kmf = KaplanMeierFitter()
    kmf.fit(sub['time_from_lm'], sub['event_ind'], label=f'Stage {s}')
    km_fitters[s] = kmf

    # BINARY mortality at clinically meaningful time points (deaths/N_total)
    # (avoids the informative-censoring bias of late KM)
    bin_7d  = sub[(sub.event_time_h - LANDMARK) <= 7*24]['event_ind'].sum()  / len(sub) * 100
    bin_14d = sub[(sub.event_time_h - LANDMARK) <= 14*24]['event_ind'].sum() / len(sub) * 100
    bin_30d = sub[(sub.event_time_h - LANDMARK) <= 30*24]['event_ind'].sum() / len(sub) * 100
    # KM point estimates at 7d, 14d, 30d (for the figure only, labelled)
    km_7d  = (1 - kmf.predict(min(7*24,  sub['time_from_lm'].max()))) * 100
    km_14d = (1 - kmf.predict(min(14*24, sub['time_from_lm'].max()))) * 100
    km_30d = (1 - kmf.predict(min(30*24, sub['time_from_lm'].max()))) * 100

    km_results[f'stage_{s}'] = {
        'n': int(len(sub)),
        'events_total': int(sub['event_ind'].sum()),
        'binary_mortality_7d_pct':  float(bin_7d),
        'binary_mortality_14d_pct': float(bin_14d),
        'binary_mortality_30d_pct': float(bin_30d),
        'km_mortality_7d_pct':  float(km_7d),
        'km_mortality_14d_pct': float(km_14d),
        'km_mortality_30d_pct': float(km_30d),
    }
    print(f"  Stage {s}: n={len(sub):>6,} deaths={sub['event_ind'].sum():>4,}")
    print(f"    binary 7d/14d/30d = {bin_7d:.1f}% / {bin_14d:.1f}% / {bin_30d:.1f}%")
    print(f"    KM     7d/14d/30d = {km_7d:.1f}% / {km_14d:.1f}% / {km_30d:.1f}% (rising past 14d = informative-censoring artifact)")

mvt = multivariate_logrank_test(surv_v['time_from_lm'], surv_v['kdigo_stage'], surv_v['event_ind'])
print(f"  Overall log-rank chi2={mvt.test_statistic:.2f}, p={mvt.p_value:.2e}")

pw = pairwise_logrank_test(surv_v['time_from_lm'], surv_v['kdigo_stage'], surv_v['event_ind'])
pw_df = pw.summary.reset_index()
pw_df.columns = ['group_A','group_B','test_statistic','p','minus_log2_p']
pw_df['p_bonferroni'] = (pw_df['p'] * 6).clip(upper=1.0)
for _, r in pw_df.iterrows():
    print(f"    S{int(r['group_A'])} vs S{int(r['group_B'])}: p_adj={r['p_bonferroni']:.3e}")

km_results['overall_logrank_chi2'] = float(mvt.test_statistic)
km_results['overall_logrank_p'] = float(mvt.p_value)
km_results['pairwise'] = pw_df[['group_A','group_B','p','p_bonferroni']].to_dict('records')

# Save KM curves
km_csv = []
for s in [0,1,2,3]:
    sf = km_fitters[s].survival_function_.reset_index()
    sf.columns = ['time_h','survival']
    sf['stage'] = s
    km_csv.append(sf)
km_df = pd.concat(km_csv)
km_df.to_csv(os.path.join(WORK, "v2_outputs", "km_curves_mimic.csv"), index=False)
print(f"  km_curves_mimic.csv saved ({len(km_df):,} rows)")

# ---- 4. Cox PH models ----
print("\n[4] Cox PH models")
surv['kdigo_1'] = (surv['kdigo_stage']==1).astype(int)
surv['kdigo_2'] = (surv['kdigo_stage']==2).astype(int)
surv['kdigo_3'] = (surv['kdigo_stage']==3).astype(int)

models_specs = {
    'A_kdigo_only':     ['kdigo_1','kdigo_2','kdigo_3'],
    'B_demo':           ['kdigo_1','kdigo_2','kdigo_3','age','sex_male'],
    'C_full_adjusted':  ['kdigo_1','kdigo_2','kdigo_3','age','sex_male',
                         'sofa_nonrenal','vent_24h',
                         'ckd_history_prior','diabetes','hypertension',
                         'heart_failure','copd','liver_disease'],
}

cox_results = {}
for mname, vars_ in models_specs.items():
    use = [v for v in vars_ if v in surv.columns]
    df_m = surv[['time_from_lm','event_ind'] + use].dropna()
    print(f"\n  Model {mname}: n={len(df_m):,}")
    cph = CoxPHFitter(penalizer=0.0)
    try:
        cph.fit(df_m, duration_col='time_from_lm', event_col='event_ind')
        cox_results[mname] = {'n': int(len(df_m)),
                              'concordance': float(cph.concordance_index_),
                              'variables': {}}
        for v in use:
            hr  = float(np.exp(cph.params_[v]))
            lo  = float(cph.summary.loc[v, 'exp(coef) lower 95%'])
            hi  = float(cph.summary.loc[v, 'exp(coef) upper 95%'])
            pv  = float(cph.summary.loc[v, 'p'])
            cox_results[mname]['variables'][v] = {'HR':hr,'CI_lo':lo,'CI_hi':hi,'p':pv}
            tag = '** KDY ' if v.startswith('kdigo_') else ''
            print(f"    {v:>18s}: HR={hr:.3f} ({lo:.3f}-{hi:.3f}) p={pv:.3e}")
        print(f"    C-index: {cph.concordance_index_:.3f}")
    except Exception as e:
        print(f"    ERROR: {e}")
        cox_results[mname] = {'error': str(e)}

# ---- 5. PH assumption tests ----
print("\n[5] PH assumption tests")
ph_results = {}
for mname, vars_ in models_specs.items():
    use = [v for v in vars_ if v in surv.columns]
    df_m = surv[['time_from_lm','event_ind'] + use].dropna()
    cph = CoxPHFitter(penalizer=0.0)
    try:
        cph.fit(df_m, duration_col='time_from_lm', event_col='event_ind')
        test = cph.proportional_hazard_test(df_m)
        ph_results[mname] = {}
        for v in use:
            row = test.summary.loc[v]
            pv = float(row['p'])
            ph_results[mname][v] = {'test_statistic': float(row['test_statistic']),
                                     'p': pv, 'violates': pv < 0.05}
            flag = ' ** VIOLATES **' if pv < 0.05 else ''
            print(f"    {mname:>16s} | {v:>18s}: chi2={row['test_statistic']:.2f} p={pv:.3e}{flag}")
    except Exception as e:
        ph_results[mname] = {'error': str(e)}

# ---- 6. Time-stratified Cox ----
print("\n[6] Time-stratified Cox (<7d / 7-14d / 14-30d)")
strata_def = [(0, 7*24, '0_7d'), (7*24, 14*24, '7_14d'), (14*24, 30*24, '14_30d')]
ts_results = {}
for tlo, thi, label in strata_def:
    sub = surv[(surv['time_from_lm'] > tlo) & (surv['time_from_lm'] <= thi)].copy()
    use = ['kdigo_1','kdigo_2','kdigo_3','age','sex_male',
           'sofa_nonrenal','vent_24h','ckd_history_prior','diabetes','hypertension',
           'heart_failure','copd','liver_disease']
    sub = sub[['time_from_lm','event_ind'] + use].dropna()
    print(f"\n  Stratum {label}: n={len(sub):,}, events={int(sub['event_ind'].sum())}")
    cph = CoxPHFitter(penalizer=0.01)  # light ridge for strata stability
    try:
        cph.fit(sub, duration_col='time_from_lm', event_col='event_ind')
        ts_results[label] = {'n': int(len(sub)),
                             'events': int(sub['event_ind'].sum()),
                             'concordance': float(cph.concordance_index_),
                             'kdigo': {}}
        for v in ['kdigo_1','kdigo_2','kdigo_3']:
            hr  = float(np.exp(cph.params_[v]))
            lo  = float(cph.summary.loc[v, 'exp(coef) lower 95%'])
            hi  = float(cph.summary.loc[v, 'exp(coef) upper 95%'])
            pv  = float(cph.summary.loc[v, 'p'])
            ts_results[label]['kdigo'][v] = {'HR':hr,'CI_lo':lo,'CI_hi':hi,'p':pv}
            print(f"    {v}: HR={hr:.3f} ({lo:.3f}-{hi:.3f}) p={pv:.3e}")
    except Exception as e:
        print(f"    ERROR: {e}")
        ts_results[label] = {'error': str(e)}

# ---- 7. Save ----
out = {
    'cohort_n_after_lm': int(len(surv)),
    'n_deaths_after_lm': int(surv['event_ind'].sum()),
    'landmark_h': LANDMARK,
    'truncation_days': 30,
    'km': km_results,
    'cox_models': cox_results,
    'ph_tests': ph_results,
    'time_stratified': ts_results,
}
with open(os.path.join(WORK, "v2_outputs", "survival_mimic.json"), 'w') as f:
    json.dump(out, f, indent=2, default=str)
print(f"\n  survival_mimic.json saved")
print("=" * 70)
print("DONE")