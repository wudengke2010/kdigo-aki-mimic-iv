"""
Supplementary analyses for manuscript revision:
1. Immortal time bias: 24h landmark analysis
2. Time-stratified Cox regression (0-7d, 7-14d, 14-30d)
3. Missing data stratified by KDIGO stage
"""
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
import statsmodels.api as sm
import json

# Load cohort
df = pd.read_csv('kdigo_cohort_revised_enriched.csv', low_memory=False)
print(f"Loaded cohort: {len(df):,} patients")
print(f"Total deaths: {df['hospital_expire_flag'].sum():,} ({df['hospital_expire_flag'].mean()*100:.1f}%)")

# ========================================
# 1. IMMORTAL TIME BIAS: 24h Landmark Analysis
# ========================================
print("\n" + "="*70)
print("1. IMMORTAL TIME BIAS: 24h Landmark Analysis")
print("="*70)

# Need survival time from ICU admission
df['icu_los_days'] = df['icu_los_hours'] / 24.0

# Patients who died within 24h - their KDIGO stages are unreliable
early_death_24h = df[(df['hospital_expire_flag'] == 1) & (df['icu_los_hours'] <= 24)]
print(f"Deaths within 24h: {len(early_death_24h):,} ({len(early_death_24h)/df['hospital_expire_flag'].sum()*100:.1f}% of deaths)")
for s in [0, 1, 2, 3]:
    n = (early_death_24h['kdigo_stage'] == s).sum()
    pct = n / len(early_death_24h) * 100
    print(f"  Stage {s}: {n:,} ({pct:.1f}%)")

# Landmark: exclude patients who died ≤24h, refit logistic regression
df_landmark = df[~((df['hospital_expire_flag'] == 1) & (df['icu_los_hours'] <= 24))].copy()
print(f"\nLandmark cohort (exclude died ≤24h): {len(df_landmark):,}")
print(f"Excluded: {len(df) - len(df_landmark):,}")

# Prepare variables for logistic regression (Model B)
vars_to_use = ['kdigo_stage', 'anchor_age', 'male', 'com_Cardiac_Arrhythmias',
               'com_Cerebrovascular', 'com_Chronic_Pulmonary', 'com_Coagulopathy',
               'com_Congestive_Heart_Failure', 'com_Diabetes_complicated',
               'com_Diabetes_uncomplicated', 'com_Fluid/Electrolyte', 'com_Hypertension',
               'com_Liver_Disease', 'com_Metastatic_Cancer', 'com_Myocardial_Infarction',
               'com_Obesity', 'com_Peptic_Ulcer', 'com_Peripheral_Vascular', 'com_Sepsis',
               'com_Solid_Tumor_non-met', 'com_Valvular_Disease',
               'sofa_total', 'mech_vent', 'vasopressor', 'ckd_icd']

def fit_logistic(data, label):
    """Fit logistic regression Model B"""
    X_data = data[vars_to_use].copy()
    # Create KDIGO dummies, ensure numeric
    for s in [1, 2, 3]:
        X_data[f'kdigo_{s}'] = (X_data['kdigo_stage'] == s).astype(float)
    X_data = X_data.drop('kdigo_stage', axis=1).astype(float)
    
    # Add constant
    X = sm.add_constant(X_data)
    y = data['hospital_expire_flag'].astype(float)
    
    # Drop any constant columns
    non_constant = X.columns[X.nunique() > 1]
    X = X[non_constant]
    
    model = sm.Logit(y, X).fit(disp=0, maxiter=200)
    results = {}
    for stage, prefix in [(1, 'kdigo_1'), (2, 'kdigo_2'), (3, 'kdigo_3')]:
        if prefix in X.columns:
            coef = model.params.get(prefix, np.nan)
            ci_low = model.conf_int().loc[prefix, 0] if prefix in model.conf_int().index else np.nan
            ci_high = model.conf_int().loc[prefix, 1] if prefix in model.conf_int().index else np.nan
            results[f'stage_{stage}'] = {
                'OR': round(float(np.exp(coef)), 3),
                'CI_low': round(float(np.exp(ci_low)), 3),
                'CI_high': round(float(np.exp(ci_high)), 3),
                'p': round(float(model.pvalues.get(prefix, np.nan)), 4)
            }
    results['N'] = len(data)
    results['AUC'] = round(concordance_index(y.values, model.predict(X)), 3)
    return results

full_results = fit_logistic(df, "Full cohort")
landmark_results = fit_logistic(df_landmark, "Landmark (excl. died ≤24h)")

print("\n=== Model B Comparison: Full vs Landmark ===")
for s in [1, 2, 3]:
    full = full_results[f'stage_{s}']
    lm = landmark_results[f'stage_{s}']
    ratio = lm['OR'] / full['OR']
    print(f"Stage {s}: Full OR={full['OR']} ({full['CI_low']}-{full['CI_high']}) vs "
          f"Landmark OR={lm['OR']} ({lm['CI_low']}-{lm['CI_high']}), ratio={ratio:.3f}, AUC: {full_results['AUC']} vs {landmark_results['AUC']}")

landmark_results['full_OR_comparison'] = {
    f'stage_{s}': {'full_OR': full_results[f'stage_{s}']['OR'], 
                   'landmark_OR': landmark_results[f'stage_{s}']['OR'],
                   'ratio': round(landmark_results[f'stage_{s}']['OR'] / full_results[f'stage_{s}']['OR'], 3)}
    for s in [1, 2, 3]
}

with open('landmark_analysis.json', 'w') as f:
    json.dump(landmark_results, f, indent=2)
print("\nSaved: landmark_analysis.json")

# ========================================
# 2. TIME-STRATIFIED COX REGRESSION
# ========================================
print("\n" + "="*70)
print("2. TIME-STRATIFIED COX (0-7d, 7-14d, 14-30d)")
print("="*70)

# Prepare survival data
# Use survival_time from Kaplan-Meier data
from lifelines import KaplanMeierFitter
df_cox = df.copy()

# We need survival time and event indicator
# 30-day survival: time = min(deathtime - admittime, dischtime - admittime, 30d)
# But we don't have deathtime/dischtime directly, use the enriched data
# Check what time columns we have
time_cols = [c for c in df.columns if 'time' in c.lower() or 'surv' in c.lower() or 'duration' in c.lower()]
print(f"Available time columns: {time_cols}")

# Let's check the cohort columns for survival time
surv_cols = [c for c in df.columns if 'surv' in c.lower() or 'duration' in c.lower() or 'follow' in c.lower()]
print(f"Survival-related columns: {surv_cols}")

# For now, approximate: ICU LOS in days as follow-up time, hospital_expire_flag as event
# This is approximate - ideally we'd use exact survival time
# For the time-stratified Cox, we can use ICU LOS truncated at 30d
df_cox['surv_time'] = np.minimum(df_cox['icu_los_hours'] / 24, 30)

print(f"\nUsing ICU LOS (capped at 30d) as survival time proxy")
print(f"Mean follow-up: {df_cox['surv_time'].mean():.1f} days")
print(f"Events: {df_cox['hospital_expire_flag'].sum():,}")

# Time-stratified Cox
cox_cols = ['kdigo_stage_1', 'kdigo_stage_2', 'kdigo_stage_3',
            'anchor_age', 'male_num'] + \
           [c for c in df_cox.columns if c.startswith('com_')] + \
           ['sofa_total', 'mech_vent', 'vasopressor', 'ckd_icd']

# Create kdigo dummy columns
for s in [1, 2, 3]:
    df_cox[f'kdigo_stage_{s}'] = (df_cox['kdigo_stage'] == s).astype(int)

# Full 30-day Cox
cox_full = CoxPHFitter(penalizer=0.01)
cox_cols_available = [c for c in cox_cols if c in df_cox.columns]
cox_data = df_cox[['surv_time', 'hospital_expire_flag'] + cox_cols_available].dropna()
cox_full.fit(cox_data, duration_col='surv_time', event_col='hospital_expire_flag')
print(f"\nFull 30-day Cox (N={len(cox_data):,}):")
for s in [1, 2, 3]:
    col = f'kdigo_stage_{s}'
    if col in cox_full.summary.index:
        hr = cox_full.summary.loc[col, 'exp(coef)']
        ci_low = cox_full.summary.loc[col, 'exp(coef) lower 95%']
        ci_high = cox_full.summary.loc[col, 'exp(coef) upper 95%']
        print(f"  Stage {s}: HR={hr:.3f} ({ci_low:.3f}-{ci_high:.3f})")

# Time-stratified: 0-7d, 7-14d, 14-30d
stratified_results = {}
for label, (t_min, t_max) in [('0-7d', (0, 7)), ('7-14d', (7, 14)), ('14-30d', (14, 30))]:
    # Restrict to patients still at risk at t_min
    at_risk = cox_data[cox_data['surv_time'] > t_min].copy()
    # Truncate follow-up at t_max
    at_risk['surv_trunc'] = np.minimum(at_risk['surv_time'] - t_min, t_max - t_min)
    # Events: death within stratum
    at_risk['event_trunc'] = at_risk['hospital_expire_flag'] & (at_risk['surv_time'] <= t_max)
    
    cox_strat = CoxPHFitter(penalizer=0.01)
    cox_strat.fit(at_risk[['surv_trunc', 'event_trunc'] + cox_cols_available],
                  duration_col='surv_trunc', event_col='event_trunc')
    
    stratified_results[label] = {'N_at_risk': len(at_risk), 'N_events': at_risk['event_trunc'].sum()}
    print(f"\n{label} (N={len(at_risk):,}, events={at_risk['event_trunc'].sum():,}):")
    for s in [1, 2, 3]:
        col = f'kdigo_stage_{s}'
        if col in cox_strat.summary.index:
            hr = cox_strat.summary.loc[col, 'exp(coef)']
            ci_low = cox_strat.summary.loc[col, 'exp(coef) lower 95%']
            ci_high = cox_strat.summary.loc[col, 'exp(coef) upper 95%']
            stratified_results[label][f'stage_{s}'] = {'HR': round(hr, 3), 'CI_low': round(ci_low, 3), 'CI_high': round(ci_high, 3)}
            print(f"  Stage {s}: HR={hr:.3f} ({ci_low:.3f}-{ci_high:.3f})")

# RMST approximation
print(f"\n=== RMST (Restricted Mean Survival Time at 30d) ===")
kmf = KaplanMeierFitter()
for s in [0, 1, 2, 3]:
    mask = df_cox['kdigo_stage'] == s
    kmf.fit(df_cox.loc[mask, 'surv_time'], 
            df_cox.loc[mask, 'hospital_expire_flag'],
            label=f'Stage {s}')
    rmst = kmf.restricted_mean_survival_time(30)
    print(f"  Stage {s}: RMST(30d)={rmst:.1f} days, N={mask.sum():,}")

with open('time_stratified_cox.json', 'w') as f:
    json.dump(stratified_results, f, indent=2)
print("\nSaved: time_stratified_cox.json")

# ========================================
# 3. MISSING DATA STRATIFIED BY KDIGO STAGE
# ========================================
print("\n" + "="*70)
print("3. MISSING DATA STRATIFIED BY KDIGO STAGE")
print("="*70)

missing_vars = ['gcs_total', 'sbp_mean', 'dbp_mean', 'mbp_mean', 'tot_bili', 
                'lactate', 'albumin', 'hr_mean', 'rr_mean', 'temp_mean', 'spo2_mean']

strat_missing = []
for var in missing_vars:
    if var in df.columns:
        row = {'variable': var, 'overall_missing_pct': round(df[var].isna().mean() * 100, 1)}
        for s in [0, 1, 2, 3]:
            mask = df['kdigo_stage'] == s
            row[f'stage_{s}_missing_pct'] = round(df.loc[mask, var].isna().mean() * 100, 1)
        strat_missing.append(row)

strat_missing_df = pd.DataFrame(strat_missing)
print(strat_missing_df.to_string(index=False))
strat_missing_df.to_csv('missing_by_stage.csv', index=False)
print("\nSaved: missing_by_stage.csv")

print("\n" + "="*70)
print("ALL ANALYSES COMPLETE")
print("="*70)
