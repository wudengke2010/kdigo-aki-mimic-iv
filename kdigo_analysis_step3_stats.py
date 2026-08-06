#!/usr/bin/env python3
"""
KDIGO AKI 分期 vs 预后 — Step 3: Statistical Analysis + Figures
=================================================================
产出:
  - Logistic regression: OR for hospital mortality by KDIGO stage
  - Kaplan-Meier survival curves
  - Cox proportional hazards
  - Forest plot data
  - All results saved for paper writing
"""
import pandas as pd
import numpy as np
import os, json
from datetime import datetime

print("=" * 70)
print("KDIGO-AKI Step 3: Statistical Analysis")
print(f"Start: {datetime.now()}")
print("=" * 70)

# ============================================================
# 1. Load cohort
# ============================================================
print("[1] Loading enriched cohort...")
IN = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"

# Try enriched first, fall back to basic
if os.path.exists(os.path.join(IN, "kdigo_cohort_enriched.csv")):
    cohort = pd.read_csv(os.path.join(IN, "kdigo_cohort_enriched.csv"), low_memory=False)
    print("  Using enriched cohort")
else:
    cohort = pd.read_csv(os.path.join(IN, "kdigo_cohort.csv"), low_memory=False)
    print("  Using basic cohort (no enrichment)")

print(f"  N = {len(cohort):,}")
print(f"  Columns: {len(cohort.columns)}")

# ============================================================
# 2. Univariate analysis
# ============================================================
print("\n[2] Univariate analysis — Mortality by KDIGO stage...")

# Crude OR
import warnings
warnings.filterwarnings('ignore')

# Mortality rates
rates = []
for s in [0, 1, 2, 3]:
    sub = cohort[cohort['kdigo_stage'] == s]
    n = len(sub)
    deaths = sub['hospital_expire_flag'].sum()
    rate = deaths / n * 100 if n > 0 else 0
    rates.append({'stage': s, 'n': n, 'deaths': int(deaths), 'mortality_rate': rate})
    print(f"  Stage {s}: {deaths}/{n} ({rate:.2f}%)")

rates_df = pd.DataFrame(rates)

# ============================================================
# 3. Multivariable Logistic Regression
# ============================================================
print("\n[3] Multivariable logistic regression (statsmodels)...")

import statsmodels.api as sm
from scipy import stats

# Define outcome and features
y = cohort['hospital_expire_flag'].values

# Features: stage dummies + demographics + key comorbidities
feature_cols = []

# KDIGO stage (categorical, reference = Stage 0)
for s in [1, 2, 3]:
    col = f'kdigo_stage_{s}'
    cohort[col] = (cohort['kdigo_stage'] == s).astype(int)
    feature_cols.append(col)

# Demographics
if 'anchor_age' in cohort.columns:
    feature_cols.append('anchor_age')
if 'gender' in cohort.columns:
    cohort['male'] = (cohort['gender'] == 'M').astype(int)
    feature_cols.append('male')

# Key comorbidities
comorb_cols = [c for c in cohort.columns if c.startswith('com_')]
# Pick top ones
key_comorb = ['com_Sepsis', 'com_Congestive_Heart_Failure', 'com_Diabetes_complicated',
              'com_Chronic_Pulmonary', 'com_Liver_Disease', 'com_Metastatic_Cancer',
              'com_Cerebrovascular', 'com_Coagulopathy', 'com_Myocardial_Infarction']
for c in key_comorb:
    if c in cohort.columns:
        feature_cols.append(c)

# CKD status
if 'ckd_icd' in cohort.columns:
    feature_cols.append('ckd_icd')

# Creatinine ratio (continuous)
if 'cr_ratio' in cohort.columns:
    # Log transform for normality
    cohort['log_cr_ratio'] = np.log(cohort['cr_ratio'].clip(lower=0.1))
    feature_cols.append('log_cr_ratio')

# Keep only complete cases
model_data = cohort[feature_cols + ['hospital_expire_flag', 'hadm_id', 'stay_id']].dropna()
X = model_data[feature_cols].values
y = model_data['hospital_expire_flag'].values

print(f"  Complete cases for model: {len(model_data):,}")
print(f"  Features: {feature_cols}")

# Fit logistic regression with statsmodels (gives OR, CI, p-value directly)
X_with_const = sm.add_constant(X)
logit_model = sm.Logit(y, X_with_const)
result = logit_model.fit(disp=0, maxiter=200)

print(f"\n  Logistic Regression Results (adjusted OR):")
print(f"  {'Variable':<40} {'OR':>8} {'95% CI':>20} {'p-value':>12}")
print(f"  {'-'*80}")

results = []
for i, col in enumerate(feature_cols):
    coef = result.params[i + 1]  # +1 for constant
    se = result.bse[i + 1]
    or_val = np.exp(coef)
    ci_low = np.exp(coef - 1.96 * se)
    ci_high = np.exp(coef + 1.96 * se)
    p_val = result.pvalues[i + 1]

    sig = '***' if p_val < 0.001 else ('**' if p_val < 0.01 else ('*' if p_val < 0.05 else ''))
    print(f"  {col:<40} {or_val:>8.3f} ({ci_low:.3f}-{ci_high:.3f}) {p_val:>10.4f} {sig}")

    results.append({
        'variable': col,
        'OR': round(float(or_val), 3),
        'CI_low': round(float(ci_low), 3),
        'CI_high': round(float(ci_high), 3),
        'p_value': round(float(p_val), 6)
    })

# Also run univariate (crude) models
print(f"\n  Univariate (crude) OR:")
print(f"  {'Variable':<40} {'OR':>8} {'95% CI':>20} {'p-value':>12}")
print(f"  {'-'*80}")
uni_results = []
for i, col in enumerate(feature_cols):
    X_uni = sm.add_constant(X[:, i:i+1])
    try:
        uni_model = sm.Logit(y, X_uni)
        uni_result = uni_model.fit(disp=0, maxiter=100)
        uni_coef = uni_result.params[1]
        uni_se = uni_result.bse[1]
        uni_or = np.exp(uni_coef)
        uni_ci_low = np.exp(uni_coef - 1.96 * uni_se)
        uni_ci_high = np.exp(uni_coef + 1.96 * uni_se)
        uni_p = uni_result.pvalues[1]
        sig = '***' if uni_p < 0.001 else ('**' if uni_p < 0.01 else ('*' if uni_p < 0.05 else ''))
        print(f"  {col:<40} {uni_or:>8.3f} ({uni_ci_low:.3f}-{uni_ci_high:.3f}) {uni_p:>10.4f} {sig}")
        uni_results.append({
            'variable': col,
            'crude_OR': round(float(uni_or), 3),
            'crude_CI_low': round(float(uni_ci_low), 3),
            'crude_CI_high': round(float(uni_ci_high), 3),
            'crude_p': round(float(uni_p), 6)
        })
    except:
        pass

# Merge univariate into results
if uni_results:
    uni_df = pd.DataFrame(uni_results)
    results_df = pd.DataFrame(results).merge(uni_df, on='variable', how='left')
else:
    results_df = pd.DataFrame(results)

results_df.to_csv(os.path.join(IN, 'logistic_results.csv'), index=False)
print(f"\n  Results saved to logistic_results.csv")

# ============================================================
# 4. Kaplan-Meier Analysis (30-day survival)
# ============================================================
print("\n[4] Kaplan-Meier survival analysis...")

# Prepare survival data
cohort['intime_dt'] = pd.to_datetime(cohort['intime'])
cohort['outtime_dt'] = pd.to_datetime(cohort['outtime'])
cohort['deathtime_dt'] = pd.to_datetime(cohort['deathtime'])
cohort['dod_dt'] = pd.to_datetime(cohort['dod'])

# Time to death (from ICU admission), censored at 30 days
cohort['event'] = 0
cohort['time_to_event'] = 30  # default censored at 30d

for idx, row in cohort.iterrows():
    if pd.notna(row['dod_dt']):
        t = (row['dod_dt'] - row['intime_dt']).total_seconds() / 86400
        if t <= 30:
            cohort.at[idx, 'time_to_event'] = t
            cohort.at[idx, 'event'] = 1
    elif row['hospital_expire_flag'] == 1 and pd.notna(row['deathtime_dt']):
        t = (row['deathtime_dt'] - row['intime_dt']).total_seconds() / 86400
        if t <= 30:
            cohort.at[idx, 'time_to_event'] = t
            cohort.at[idx, 'event'] = 1

print(f"  Events (30d death): {cohort['event'].sum():,}")
print(f"  Censored: {(cohort['event'] == 0).sum():,}")

# Generate KM data for plotting
try:
    from lifelines import KaplanMeierFitter
    LIFELINES_OK = True
except:
    LIFELINES_OK = False
    print("  lifelines not available, computing KM manually")

km_data_by_stage = {}
for stage in [0, 1, 2, 3]:
    sub = cohort[cohort['kdigo_stage'] == stage]
    
    if LIFELINES_OK:
        kmf = KaplanMeierFitter()
        kmf.fit(sub['time_to_event'], sub['event'], label=f'Stage {stage}')
        
        # Export survival table
        surv_table = kmf.survival_function_
        times = surv_table.index.values
        surv = surv_table.values.flatten()
        
        km_data_by_stage[str(stage)] = {
            'times': times.tolist(),
            'survival': surv.tolist(),
            'n': len(sub),
            'events': int(sub['event'].sum()),
        }
    else:
        # Manual KM calculation
        df = sub[['time_to_event', 'event']].sort_values('time_to_event')
        times = sorted(df['time_to_event'].unique())
        surv_probs = []
        n_risk = len(df)
        surv = 1.0
        for t in times:
            events_at_t = df[(df['time_to_event'] == t) & (df['event'] == 1)].shape[0]
            if n_risk > 0:
                surv *= (1 - events_at_t / n_risk)
                n_risk -= events_at_t
            surv_probs.append(surv)
        
        km_data_by_stage[str(stage)] = {
            'times': times,
            'survival': surv_probs,
            'n': len(sub),
            'events': int(sub['event'].sum()),
        }

# Log-rank test
if LIFELINES_OK:
    from lifelines.statistics import multivariate_logrank_test
    lr_result = multivariate_logrank_test(
        cohort['time_to_event'], 
        cohort['kdigo_stage'], 
        cohort['event']
    )
    print(f"\n  Log-rank test: chi2={lr_result.test_statistic:.2f}, p={lr_result.p_value:.6f}")
    
    # Pairwise tests
    from lifelines.statistics import logrank_test
    print("  Pairwise log-rank tests:")
    stages = [0, 1, 2, 3]
    for i, s1 in enumerate(stages):
        for s2 in stages[i+1:]:
            g1 = cohort[cohort['kdigo_stage'] == s1]
            g2 = cohort[cohort['kdigo_stage'] == s2]
            lr = logrank_test(g1['time_to_event'], g2['time_to_event'], g1['event'], g2['event'])
            print(f"    Stage {s1} vs Stage {s2}: p={lr.p_value:.6f}")

# Save KM data
with open(os.path.join(IN, 'km_survival_data.json'), 'w') as f:
    json.dump(km_data_by_stage, f, indent=2)

# ============================================================
# 5. Cox Proportional Hazards
# ============================================================
print("\n[5] Cox proportional hazards...")

if LIFELINES_OK:
    from lifelines import CoxPHFitter
    
    cox_data = cohort[feature_cols + ['time_to_event', 'event']].dropna()
    
    cph = CoxPHFitter()
    cph.fit(cox_data, duration_col='time_to_event', event_col='event')
    
    print("\n  Cox PH Results:")
    cph.print_summary()
    
    # Save
    cox_summary = cph.summary
    cox_summary.to_csv(os.path.join(IN, 'cox_results.csv'))
    
    # Hazard ratios
    print("\n  Hazard Ratios:")
    for col in cox_summary.index:
        hr = np.exp(cox_summary.loc[col, 'coef'])
        ci_low = np.exp(cox_summary.loc[col, 'coef lower 95%'])
        ci_high = np.exp(cox_summary.loc[col, 'coef upper 95%'])
        p = cox_summary.loc[col, 'p']
        print(f"    {col:<40} HR={hr:.3f} ({ci_low:.3f}-{ci_high:.3f}), p={p:.4f}")

# ============================================================
# 6. Dose-response: Mortality by Cr ratio decile
# ============================================================
print("\n[6] Dose-response analysis...")

# Bin Cr ratio
cohort['cr_ratio_bin'] = pd.cut(cohort['cr_ratio'], 
    bins=[0, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 100.0],
    labels=['<1.0', '1.0-1.5', '1.5-2.0', '2.0-3.0', '3.0-5.0', '5.0-10.0', '>=10.0'])

dr_data = cohort.groupby('cr_ratio_bin', observed=False).agg(
    n=('stay_id', 'count'),
    deaths=('hospital_expire_flag', 'sum')
).reset_index()
dr_data['mortality_rate'] = dr_data['deaths'] / dr_data['n'] * 100
dr_data = dr_data.dropna()

print("  Mortality by Cr ratio bin:")
for _, row in dr_data.iterrows():
    print(f"    {row['cr_ratio_bin']}: {row['mortality_rate']:.1f}% ({int(row['deaths'])}/{int(row['n'])})")

dr_data.to_csv(os.path.join(IN, 'dose_response.csv'), index=False)

# ============================================================
# 7. Subgroup Analysis: AKI effect by CKD status
# ============================================================
print("\n[7] Subgroup analysis — AKI effect stratified by CKD...")

for ckd_status in [0, 1]:
    sub = cohort[cohort['ckd_icd'] == ckd_status]
    label = 'CKD+' if ckd_status == 1 else 'CKD-'
    if len(sub) > 100:
        rates = {}
        for s in [0, 1, 2, 3]:
            ss = sub[sub['kdigo_stage'] == s]
            if len(ss) > 0:
                rates[s] = ss['hospital_expire_flag'].mean() * 100
        print(f"  {label}: Stage mortality = {rates}")

# ============================================================
# 8. Output summary for paper
# ============================================================
print("\n[8] Generating paper summary...")

# Overall stats
n_total = len(cohort)
n_aki_cr = (cohort['kdigo_stage'] >= 1).sum()
n_aki_icd = cohort['aki_icd'].sum() if 'aki_icd' in cohort.columns else 0

summary_text = f"""
KDIGO-Cr AKI Staging: Results Summary
=====================================

Cohort: {n_total:,} ICU stays with valid creatinine data

KDIGO Distribution (by Cr criteria):
  No AKI:       {(cohort['kdigo_stage'] == 0).sum():,} ({(cohort['kdigo_stage'] == 0).sum()/n_total*100:.1f}%)
  Stage 1:      {(cohort['kdigo_stage'] == 1).sum():,} ({(cohort['kdigo_stage'] == 1).sum()/n_total*100:.1f}%)
  Stage 2:      {(cohort['kdigo_stage'] == 2).sum():,} ({(cohort['kdigo_stage'] == 2).sum()/n_total*100:.1f}%)
  Stage 3:      {(cohort['kdigo_stage'] == 3).sum():,} ({(cohort['kdigo_stage'] == 3).sum()/n_total*100:.1f}%)
  Total AKI:    {n_aki_cr:,} ({n_aki_cr/n_total*100:.1f}%)

Primary Outcome — Hospital Mortality:
  No AKI:   {(cohort[cohort['kdigo_stage'] == 0]['hospital_expire_flag'].mean()*100):.2f}%
  Stage 1:  {(cohort[cohort['kdigo_stage'] == 1]['hospital_expire_flag'].mean()*100):.2f}%
  Stage 2:  {(cohort[cohort['kdigo_stage'] == 2]['hospital_expire_flag'].mean()*100):.2f}%
  Stage 3:  {(cohort[cohort['kdigo_stage'] == 3]['hospital_expire_flag'].mean()*100):.2f}%
  Overall:  {(cohort['hospital_expire_flag'].mean()*100):.2f}%

Secondary Outcomes:
  RRT rate (Stage 3): {(cohort[(cohort['kdigo_stage'] == 3) & (cohort['rrt'] == 1)].shape[0] / max(cohort[cohort['kdigo_stage'] == 3].shape[0], 1) * 100):.1f}%
  
Mean ICU LOS (hours):
  No AKI: {cohort[cohort['kdigo_stage'] == 0]['icu_los_hours'].mean():.1f}
  Stage 1: {cohort[cohort['kdigo_stage'] == 1]['icu_los_hours'].mean():.1f}
  Stage 2: {cohort[cohort['kdigo_stage'] == 2]['icu_los_hours'].mean():.1f}
  Stage 3: {cohort[cohort['kdigo_stage'] == 3]['icu_los_hours'].mean():.1f}

Mean Hospital LOS (days):
  No AKI: {cohort[cohort['kdigo_stage'] == 0]['hosp_los_days'].mean():.1f}
  Stage 1: {cohort[cohort['kdigo_stage'] == 1]['hosp_los_days'].mean():.1f}
  Stage 2: {cohort[cohort['kdigo_stage'] == 2]['hosp_los_days'].mean():.1f}
  Stage 3: {cohort[cohort['kdigo_stage'] == 3]['hosp_los_days'].mean():.1f}

Creatinine metrics:
  Mean baseline Cr: {cohort['baseline_cr'].mean():.2f} mg/dL
  Mean peak Cr: {cohort['peak_cr'].mean():.2f} mg/dL
  Mean Cr ratio: {cohort['cr_ratio'].mean():.2f}
"""

with open(os.path.join(IN, 'paper_summary.txt'), 'w') as f:
    f.write(summary_text)

print(summary_text)
print("\n" + "=" * 70)
print("Step 3 Complete! All results saved.")
print("=" * 70)
