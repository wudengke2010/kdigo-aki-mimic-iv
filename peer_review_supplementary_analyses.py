"""
Supplementary analyses for peer review response:
1. Cox PH assumption (Schoenfeld residuals)
2. Complete-case SOFA sensitivity analysis
3. VIF collinearity diagnostics
4. Systematic missing data report
5. Temporal trend analysis
6. DCA quantitative results extraction
"""
import pandas as pd
import numpy as np
import json
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Load data
# ============================================================
print("=" * 70)
print("Loading cohort data...")
df = pd.read_csv('kdigo_cohort_revised_enriched.csv')
print(f"Cohort: {len(df)} patients")

# Prepare variables
df['kdigo_1'] = (df['kdigo_stage'] == 1).astype(int)
df['kdigo_2'] = (df['kdigo_stage'] == 2).astype(int)
df['kdigo_3'] = (df['kdigo_stage'] == 3).astype(int)
df['log_cr_ratio'] = np.log(df['cr_ratio'].clip(lower=0.01))
df['age'] = df['anchor_age']
df['male'] = df['male'].astype(int)

# ============================================================
# 1. SYSTEMATIC MISSING DATA REPORT
# ============================================================
print("\n" + "=" * 70)
print("1. SYSTEMATIC MISSING DATA REPORT")
print("=" * 70)

missing_vars = {
    'age': 'Age',
    'male': 'Sex',
    'sofa_total': 'SOFA score',
    'sofa_neuro': 'SOFA neuro (GCS)',
    'sofa_coag': 'SOFA coag (platelet)',
    'sofa_liver': 'SOFA liver (bilirubin)',
    'sofa_cv': 'SOFA cardiovascular (MAP)',
    'sofa_renal': 'SOFA renal (creatinine)',
    'sofa_resp': 'SOFA resp (vent/PaO2)',
    'mech_vent': 'Mechanical ventilation',
    'vasopressor': 'Vasopressor use',
    'baseline_cr': 'Baseline creatinine',
    'peak_cr': 'Peak creatinine',
    'cr_ratio': 'Creatinine ratio',
    'lab_bun': 'BUN',
    'lab_sodium': 'Sodium',
    'lab_potassium': 'Potassium',
    'lab_bicarbonate': 'Bicarbonate',
    'lab_lactate': 'Lactate',
    'lab_wbc': 'WBC',
    'lab_hemoglobin': 'Hemoglobin',
    'lab_platelet': 'Platelet',
    'lab_bilirubin_total': 'Total bilirubin',
    'lab_albumin': 'Albumin',
    'lab_glucose': 'Glucose',
    'heart_rate_mean': 'Heart rate',
    'sbp_mean': 'Systolic BP',
    'dbp_mean': 'Diastolic BP',
    'mbp_mean': 'Mean arterial pressure',
    'resp_rate_mean': 'Respiratory rate',
    'temperature_mean': 'Temperature',
    'spo2_mean': 'SpO2',
    'gcs_total': 'GCS total',
    'ckd_icd': 'CKD (ICD)',
    'ckd_lab': 'CKD (lab)',
    'com_Sepsis': 'Sepsis',
}

missing_report = []
N = len(df)
for col, label in missing_vars.items():
    if col in df.columns:
        n_missing = df[col].isna().sum()
        # For SOFA components, 0 might mean imputed
        pct_missing = n_missing / N * 100
        missing_report.append({
            'Variable': label,
            'Column': col,
            'N_missing': n_missing,
            'Pct_missing': round(pct_missing, 1),
        })

missing_df = pd.DataFrame(missing_report)
print(missing_df.to_string(index=False))
missing_df.to_csv('missing_data_report.csv', index=False)

# Also check SOFA component imputation (0 = imputed for neuro/liver)
print("\nSOFA component zero-frequencies (imputed as normal):")
for comp in ['sofa_neuro', 'sofa_liver', 'sofa_coag', 'sofa_cv', 'sofa_renal', 'sofa_resp']:
    n_zero = (df[comp] == 0).sum()
    print(f"  {comp}: {n_zero} ({n_zero/N*100:.1f}%) with score=0")

# Complete case for SOFA (non-imputed neuro + liver)
# GCS available means sofa_neuro was not imputed
# We use gcs_total not NaN as proxy for "GCS was measured"
gcs_available = df['gcs_total'].notna()
bilirubin_available = df['lab_bilirubin_total'].notna()
complete_sofa = gcs_available & bilirubin_available
print(f"\nComplete-case SOFA (GCS + bilirubin available): {complete_sofa.sum()} ({complete_sofa.sum()/N*100:.1f}%)")

# ============================================================
# 2. COMPLETE-CASE SOFA SENSITIVITY ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("2. COMPLETE-CASE SOFA SENSITIVITY ANALYSIS")
print("=" * 70)

import statsmodels.api as sm

# Model B covariates
comorb_cols = [c for c in df.columns if c.startswith('com_')]
model_b_vars = ['kdigo_1', 'kdigo_2', 'kdigo_3', 'age', 'male'] + comorb_cols + ['sofa_total', 'mech_vent', 'vasopressor']

# Full cohort Model B
df_full = df.copy()
df_full['mortality'] = df_full['hospital_expire_flag']
X_full = sm.add_constant(df_full[model_b_vars])
y_full = df_full['mortality']
try:
    model_b_full = sm.Logit(y_full, X_full).fit(disp=0, maxiter=100)
    print("\nFull cohort Model B (N=%d):" % len(df_full))
    for stage in ['kdigo_1', 'kdigo_2', 'kdigo_3']:
        or_val = np.exp(model_b_full.params[stage])
        ci_low = np.exp(model_b_full.conf_int().loc[stage, 0])
        ci_high = np.exp(model_b_full.conf_int().loc[stage, 1])
        print(f"  {stage}: aOR={or_val:.3f} ({ci_low:.3f}-{ci_high:.3f})")
except Exception as e:
    print(f"Full model error: {e}")

# Complete-case Model B
df_cc = df[complete_sofa].copy()
df_cc['mortality'] = df_cc['hospital_expire_flag']
X_cc = sm.add_constant(df_cc[model_b_vars])
y_cc = df_cc['mortality']
try:
    model_b_cc = sm.Logit(y_cc, X_cc).fit(disp=0, maxiter=100)
    print(f"\nComplete-case Model B (N={len(df_cc)}):")
    for stage in ['kdigo_1', 'kdigo_2', 'kdigo_3']:
        or_val = np.exp(model_b_cc.params[stage])
        ci_low = np.exp(model_b_cc.conf_int().loc[stage, 0])
        ci_high = np.exp(model_b_cc.conf_int().loc[stage, 1])
        print(f"  {stage}: aOR={or_val:.3f} ({ci_low:.3f}-{ci_high:.3f})")
    
    # Compare AUC
    from sklearn.metrics import roc_auc_score
    auc_full = roc_auc_score(y_full, model_b_full.predict(X_full))
    auc_cc = roc_auc_score(y_cc, model_b_cc.predict(X_cc))
    print(f"\nAUC full: {auc_full:.3f}, AUC complete-case: {auc_cc:.3f}")
    
    # Save results
    cc_results = {}
    for stage in ['kdigo_1', 'kdigo_2', 'kdigo_3']:
        cc_results[stage] = {
            'aOR': round(np.exp(model_b_cc.params[stage]), 3),
            'CI_low': round(np.exp(model_b_cc.conf_int().loc[stage, 0]), 3),
            'CI_high': round(np.exp(model_b_cc.conf_int().loc[stage, 1]), 3),
        }
    cc_results['N_complete_case'] = int(len(df_cc))
    cc_results['N_full'] = int(len(df_full))
    cc_results['AUC_full'] = round(auc_full, 3)
    cc_results['AUC_cc'] = round(auc_cc, 3)
    
    # Attenuation comparison
    for stage_name, stage_var in [('Stage 1', 'kdigo_1'), ('Stage 2', 'kdigo_2'), ('Stage 3', 'kdigo_3')]:
        or_full = np.exp(model_b_full.params[stage_var])
        or_cc = np.exp(model_b_cc.params[stage_var])
        att_full = (1 - or_full / np.exp(model_b_full.params[stage_var])) * 100  # placeholder
        print(f"\n{stage_name}: Full aOR={or_full:.3f}, CC aOR={or_cc:.3f}")
    
    with open('complete_case_sensitivity.json', 'w') as f:
        json.dump(cc_results, f, indent=2)
    print("\nResults saved to complete_case_sensitivity.json")
except Exception as e:
    print(f"Complete-case model error: {e}")

# ============================================================
# 3. VIF COLLINEARITY DIAGNOSTICS
# ============================================================
print("\n" + "=" * 70)
print("3. VIF COLLINEARITY DIAGNOSTICS")
print("=" * 70)

from statsmodels.stats.outliers_influence import variance_inflation_factor

# Check VIF for key variables
vif_vars = ['sofa_total', 'mech_vent', 'vasopressor', 'age', 'kdigo_1', 'kdigo_2', 'kdigo_3']
vif_data = df_full[vif_vars].dropna()
X_vif = sm.add_constant(vif_data)

vif_results = {}
print(f"\nVIF for key variables (N={len(vif_data)}):")
for i, col in enumerate(X_vif.columns):
    if col == 'const':
        continue
    vif = variance_inflation_factor(X_vif.values, i)
    vif_results[col] = round(vif, 2)
    print(f"  {col}: VIF = {vif:.2f}")

with open('vif_results.json', 'w') as f:
    json.dump(vif_results, f, indent=2)
print("Results saved to vif_results.json")

# ============================================================
# 4. COX PH ASSUMPTION TEST (Schoenfeld residuals)
# ============================================================
print("\n" + "=" * 70)
print("4. COX PROPORTIONAL HAZARDS ASSUMPTION TEST")
print("=" * 70)

from lifelines import CoxPHFitter

# Prepare survival data
df_surv = df.copy()
df_surv['mortality'] = df_surv['hospital_expire_flag']

# Survival time: ICU stay to death or 30 days, censored at discharge
df_surv['surv_time'] = df_surv['icu_los_hours'] / 24.0  # convert to days
df_surv.loc[df_surv['surv_time'] > 30, 'surv_time'] = 30
df_surv['event'] = ((df_surv['mortality'] == 1) & (df_surv['surv_time'] <= 30)).astype(int)

# Model B covariates for Cox
cox_vars = ['kdigo_1', 'kdigo_2', 'kdigo_3', 'age', 'male', 'sofa_total', 'mech_vent', 'vasopressor']
# Add comorbidities
cox_vars += comorb_cols

df_cox = df_surv[['surv_time', 'event'] + cox_vars].dropna()
print(f"Cox model N: {len(df_cox)}")

try:
    cph = CoxPHFitter()
    cph.fit(df_cox, duration_col='surv_time', event_col='event')
    
    # Schoenfeld residuals test
    ph_test = cph.check_assumptions(df_cox, p_value_threshold=0.05, show_plots=False)
    
    # Manual Schoenfeld test using proportional_hazard_test
    from lifelines.statistics import proportional_hazard_test
    results = proportional_hazard_test(cph, df_cox, time_transform='rank')
    
    print("\nSchoenfeld residuals test results:")
    ph_results = {}
    for var in cox_vars:
        if var in results.summary.index:
            stat = results.summary.loc[var, 'test_statistic']
            p = results.summary.loc[var, 'p']
            ph_results[var] = {'chi2': round(stat, 3), 'p': round(p, 4)}
            status = "VIOLATED" if p < 0.05 else "OK"
            print(f"  {var}: chi2={stat:.3f}, p={p:.4f} [{status}]")
    
    with open('cox_ph_test.json', 'w') as f:
        json.dump(ph_results, f, indent=2)
    print("\nResults saved to cox_ph_test.json")
except Exception as e:
    print(f"Cox PH test error: {e}")
    import traceback
    traceback.print_exc()

# ============================================================
# 5. TEMPORAL TREND ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("5. TEMPORAL TREND ANALYSIS")
print("=" * 70)

# Extract admission year from admittime
df['admit_year'] = pd.to_datetime(df['admittime']).dt.year

# Group by 3-year periods
def year_to_period(y):
    if y <= 2010:
        return '2008-2010'
    elif y <= 2013:
        return '2011-2013'
    elif y <= 2016:
        return '2014-2016'
    elif y <= 2019:
        return '2017-2019'
    else:
        return '2020-2022'

df['period'] = df['admit_year'].apply(year_to_period)

print("\nAKI incidence and mortality by time period:")
temporal = df.groupby('period').agg(
    N=('kdigo_stage', 'count'),
    aki_pct=('kdigo_stage', lambda x: (x > 0).mean() * 100),
    mortality_pct=('hospital_expire_flag', lambda x: x.mean() * 100),
    stage3_pct=('kdigo_stage', lambda x: (x == 3).mean() * 100),
    mean_sofa=('sofa_total', 'mean'),
    rrt_pct=('rrt', lambda x: x.mean() * 100),
).reset_index()

# Also compute AKI mortality by period
aki_df = df[df['kdigo_stage'] > 0]
aki_temporal = aki_df.groupby('period').agg(
    aki_mortality=('hospital_expire_flag', lambda x: x.mean() * 100),
    N_aki=('kdigo_stage', 'count'),
).reset_index()
aki_temporal.columns = ['period', 'aki_mortality_pct', 'N_aki']

temporal = temporal.merge(aki_temporal, on='period')
print(temporal.to_string(index=False))
temporal.to_csv('temporal_trend.csv', index=False)

# Also by individual year
print("\nBy individual year:")
yearly = df.groupby('admit_year').agg(
    N=('kdigo_stage', 'count'),
    aki_pct=('kdigo_stage', lambda x: (x > 0).mean() * 100),
    mortality_pct=('hospital_expire_flag', lambda x: x.mean() * 100),
).reset_index()
print(yearly.to_string(index=False))

# ============================================================
# 6. DCA QUANTITATIVE RESULTS
# ============================================================
print("\n" + "=" * 70)
print("6. DCA QUANTITATIVE RESULTS")
print("=" * 70)

try:
    dca = pd.read_csv('dca_results.csv')
    print(f"DCA data: {len(dca)} rows")
    print(f"Columns: {list(dca.columns)}")
    print(dca.head(10).to_string())
    
    # Find threshold range where Model B provides net benefit
    model_b_cols = [c for c in dca.columns if 'model_b' in c.lower() or 'Model_B' in c]
    print(f"\nModel B columns: {model_b_cols}")
    print(dca.describe().to_string())
except Exception as e:
    print(f"DCA read error: {e}")

print("\n" + "=" * 70)
print("ALL SUPPLEMENTARY ANALYSES COMPLETE")
print("=" * 70)
