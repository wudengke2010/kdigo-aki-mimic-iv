#!/usr/bin/env python3
"""
MICE Multiple Imputation for AKI Paradox Study
================================================
Strategy:
1. Mark SOFA components as NaN where underlying lab values are missing
   (sofa_liver where bili_max missing, sofa_resp where pao2fio2_min missing, etc.)
2. Run MICE with m=20 imputations using statsmodels
3. For each imputed dataset, recalculate sofa_total & sofa_nonrenal
4. Run 4 logistic models (A, B, C, D) on each imputed dataset
5. Pool results using Rubin's rules
6. Compare with original results
"""
import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.imputation.mice import MICEData
import warnings
import json
import copy

warnings.filterwarnings('ignore')

# ============================================================
# 1. LOAD DATA
# ============================================================
print("=" * 70)
print("STEP 1: Loading data...")
print("=" * 70)

df = pd.read_csv(r'C:\Users\admin\WorkBuddy\2026-07-07-19-40-19\aki_paradox_cohort_with_sofa.csv')
aki_df = df[df['max_aki_stage'] >= 1].copy()
print(f"AKI subset: {len(aki_df)} patients")
print(f"  Stage 1: {(aki_df['max_aki_stage']==1).sum()}, Stage 2: {(aki_df['max_aki_stage']==2).sum()}, Stage 3: {(aki_df['max_aki_stage']==3).sum()}")
print(f"  CKD: {aki_df['ckd'].sum()}, Non-CKD: {(aki_df['ckd']==0).sum()}")

# ============================================================
# 2. PREPARE MISSINGNESS
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: Marking missing SOFA components...")
print("=" * 70)

# Mark SOFA components as NaN where underlying lab values are missing
# sofa_liver: based on bili_max
mask_liver = aki_df['bili_max'].isna()
aki_df.loc[mask_liver, 'sofa_liver'] = np.nan
print(f"sofa_liver set to NaN: {mask_liver.sum()} ({mask_liver.sum()/len(aki_df)*100:.1f}%)")

# sofa_resp: based on pao2fio2_min
mask_resp = aki_df['pao2fio2_min'].isna()
aki_df.loc[mask_resp, 'sofa_resp'] = np.nan
print(f"sofa_resp set to NaN: {mask_resp.sum()} ({mask_resp.sum()/len(aki_df)*100:.1f}%)")

# sofa_coag: based on platelet_min
mask_coag = aki_df['platelet_min'].isna()
aki_df.loc[mask_coag, 'sofa_coag'] = np.nan
print(f"sofa_coag set to NaN: {mask_coag.sum()} ({mask_coag.sum()/len(aki_df)*100:.1f}%)")

# sofa_cns: based on gcs_min
mask_cns = aki_df['gcs_min'].isna()
aki_df.loc[mask_cns, 'sofa_cns'] = np.nan
print(f"sofa_cns set to NaN: {mask_cns.sum()} ({mask_cns.sum()/len(aki_df)*100:.1f}%)")

# sofa_cv: based on map_min AND no vasopressors
# If vasopressors present, sofa_cv is calculated from them regardless of MAP
has_vaso = (aki_df['norepi_max_rate'] > 0) | (aki_df['epi_max_rate'] > 0) | \
           (aki_df['dopa_max_rate'] > 0) | (aki_df['doba_max_rate'] > 0)
mask_cv = aki_df['map_min'].isna() & (~has_vaso)
aki_df.loc[mask_cv, 'sofa_cv'] = np.nan
print(f"sofa_cv set to NaN: {mask_cv.sum()} ({mask_cv.sum()/len(aki_df)*100:.1f}%)")

# ============================================================
# 3. PREPARE IMPUTATION DATAFRAME
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: Preparing imputation dataframe...")
print("=" * 70)

# Convert gender to numeric
aki_df['gender_num'] = (aki_df['gender'] == 'M').astype(int)

# Variables for imputation model (include all analysis vars + auxiliaries)
imp_vars = [
    # Analysis variables
    'age', 'gender_num', 'ckd', 'htn', 'dm', 'hf',
    'max_aki_stage', 'hospital_expire_flag',
    # SOFA components (some have NaN to be imputed)
    'sofa_resp', 'sofa_coag', 'sofa_liver', 'sofa_cv', 'sofa_cns', 'sofa_renal',
    # Auxiliary variables (help predict missingness)
    'norepi_max_rate', 'los', 'weight'
]

imp_df = aki_df[imp_vars].copy()
print(f"Imputation dataframe shape: {imp_df.shape}")
print(f"Variables: {imp_vars}")

# Check missingness in imputation dataframe
print("\nMissingness in imputation dataframe:")
for col in imp_df.columns:
    n_miss = imp_df[col].isna().sum()
    if n_miss > 0:
        print(f"  {col}: {n_miss} ({n_miss/len(imp_df)*100:.1f}%)")

# ============================================================
# 4. RUN MICE (m=20)
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: Running MICE (m=20, burn-in=10)...")
print("=" * 70)

np.random.seed(42)
m = 20
burn_in = 10

# Create MICEData object
imp = MICEData(imp_df, k_pmm=5)  # Predictive mean matching with 5 nearest neighbors

# Set up imputers for each variable with missing data
# Use all other variables as predictors
missing_vars = ['sofa_liver', 'sofa_resp', 'sofa_coag', 'sofa_cns', 'sofa_cv']
for var in missing_vars:
    formula = ' + '.join([v for v in imp_vars if v != var])
    imp.set_imputer(var, formula=formula, model_class=sm.OLS, fit_kwds={})

# Generate m imputed datasets
imputed_datasets = []
print(f"Generating {m} imputed datasets...")

for i in range(m):
    # Update the imputation with new random seed
    imp.update_all()
    
    # Get the imputed data
    imp_data = imp.data.copy()
    
    # Clip SOFA components to valid range [0, 4]
    for comp in missing_vars:
        imp_data[comp] = imp_data[comp].clip(0, 4)
        imp_data[comp] = imp_data[comp].round().astype(int)  # SOFA scores are integers
    
    # Recalculate sofa_total and sofa_nonrenal
    imp_data['sofa_total'] = (imp_data['sofa_resp'] + imp_data['sofa_coag'] + 
                              imp_data['sofa_liver'] + imp_data['sofa_cv'] + 
                              imp_data['sofa_cns'] + imp_data['sofa_renal'])
    imp_data['sofa_nonrenal'] = imp_data['sofa_total'] - imp_data['sofa_renal']
    
    # Add back necessary columns
    imp_data['ckd'] = aki_df['ckd'].values
    imp_data['max_aki_stage'] = aki_df['max_aki_stage'].values
    imp_data['hospital_expire_flag'] = aki_df['hospital_expire_flag'].values
    imp_data['age'] = aki_df['age'].values
    imp_data['gender_num'] = (aki_df['gender'] == 'M').astype(int).values
    imp_data['htn'] = aki_df['htn'].values
    imp_data['dm'] = aki_df['dm'].values
    imp_data['hf'] = aki_df['hf'].values
    
    imputed_datasets.append(imp_data)
    
    if (i + 1) % 5 == 0:
        print(f"  Completed imputation {i+1}/{m}")

print(f"Generated {len(imputed_datasets)} imputed datasets")

# Verify imputation quality
print("\nImputation quality check (imputed dataset #1):")
d1 = imputed_datasets[0]
for comp in missing_vars:
    print(f"  {comp}: mean={d1[comp].mean():.3f}, range=[{d1[comp].min()}, {d1[comp].max()}]")
print(f"  sofa_total: mean={d1['sofa_total'].mean():.3f}")
print(f"  sofa_nonrenal: mean={d1['sofa_nonrenal'].mean():.3f}")

# Compare with original (before marking NaN)
print("\nOriginal SOFA (with missing=0):")
orig_liver = aki_df['sofa_liver'].fillna(0)
orig_resp = aki_df['sofa_resp'].fillna(0)
print(f"  sofa_liver: mean={orig_liver.mean():.3f}")
print(f"  sofa_resp: mean={orig_resp.mean():.3f}")
orig_total = (aki_df['sofa_resp'].fillna(0) + aki_df['sofa_coag'].fillna(0) + 
              aki_df['sofa_liver'].fillna(0) + aki_df['sofa_cv'].fillna(0) + 
              aki_df['sofa_cns'].fillna(0) + aki_df['sofa_renal'].fillna(0))
print(f"  sofa_total: mean={orig_total.mean():.3f}")

# ============================================================
# 5. DEFINE AND RUN 4 LOGISTIC MODELS
# ============================================================
print("\n" + "=" * 70)
print("STEP 5: Running 4 logistic models on each imputed dataset...")
print("=" * 70)

def run_model(data, model_name):
    """Run a logistic regression model and return results."""
    # Create interaction terms
    data = data.copy()
    data['stage2'] = (data['max_aki_stage'] == 2).astype(int)
    data['stage3'] = (data['max_aki_stage'] == 3).astype(int)
    data['ckd_stage2'] = data['ckd'] * data['stage2']
    data['ckd_stage3'] = data['ckd'] * data['stage3']
    
    y = data['hospital_expire_flag']
    
    if model_name == 'A':
        # Base: age, sex, comorbidities, total SOFA, AKI stage (no CKD)
        X = data[['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_total', 'stage2', 'stage3']]
    elif model_name == 'B':
        # + CKD main effect
        X = data[['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_total', 'stage2', 'stage3', 'ckd']]
    elif model_name == 'C':
        # + CKD × AKI stage interaction
        X = data[['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_total', 'stage2', 'stage3', 
                  'ckd', 'ckd_stage2', 'ckd_stage3']]
    elif model_name == 'D':
        # Replace total SOFA with 6 individual SOFA components + CKD + interaction
        X = data[['age', 'gender_num', 'htn', 'dm', 'hf', 
                  'sofa_resp', 'sofa_coag', 'sofa_liver', 'sofa_cv', 'sofa_cns', 'sofa_renal',
                  'stage2', 'stage3', 'ckd', 'ckd_stage2', 'ckd_stage3']]
    
    X = sm.add_constant(X)
    
    try:
        model = sm.Logit(y, X)
        result = model.fit(disp=0, maxiter=200)
        
        # Get coefficients and covariance
        params = result.params
        cov = result.cov_params()
        
        # Calculate OR and CI
        or_vals = np.exp(params)
        se = np.sqrt(np.diag(cov))
        ci_lower = np.exp(params - 1.96 * se)
        ci_upper = np.exp(params + 1.96 * se)
        p_values = result.pvalues
        
        # Log-likelihood
        llf = result.llf
        
        # AUC
        from sklearn.metrics import roc_auc_score
        y_pred = result.predict(X)
        auc = roc_auc_score(y, y_pred)
        
        return {
            'params': params,
            'cov': cov,
            'or': or_vals,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'p_values': p_values,
            'llf': llf,
            'auc': auc,
            'converged': result.mle_retvals['converged']
        }
    except Exception as e:
        print(f"  ERROR in {model_name}: {e}")
        return None

# Also run on original data (with missing=0) for comparison
print("Running models on ORIGINAL data (missing SOFA = 0)...")
orig_data = aki_df.copy()
orig_data['gender_num'] = (orig_data['gender'] == 'M').astype(int)
orig_data['sofa_total_orig'] = (orig_data['sofa_resp'].fillna(0) + orig_data['sofa_coag'].fillna(0) + 
                                 orig_data['sofa_liver'].fillna(0) + orig_data['sofa_cv'].fillna(0) + 
                                 orig_data['sofa_cns'].fillna(0) + orig_data['sofa_renal'].fillna(0))
orig_data['sofa_nonrenal_orig'] = orig_data['sofa_total_orig'] - orig_data['sofa_renal'].fillna(0)
# Use original sofa_total for models
orig_data['sofa_total'] = orig_data['sofa_total_orig']
for comp in missing_vars:
    orig_data[comp] = orig_data[comp].fillna(0)

orig_results = {}
for model_name in ['A', 'B', 'C', 'D']:
    result = run_model(orig_data, model_name)
    if result:
        orig_results[model_name] = result
        print(f"  Model {model_name}: AUC={result['auc']:.4f}, LL={result['llf']:.1f}, converged={result['converged']}")

# Run models on each imputed dataset
print(f"\nRunning models on {m} imputed datasets...")
mice_results = {model_name: [] for model_name in ['A', 'B', 'C', 'D']}

for i, imp_data in enumerate(imputed_datasets):
    for model_name in ['A', 'B', 'C', 'D']:
        result = run_model(imp_data, model_name)
        if result:
            mice_results[model_name].append(result)
    
    if (i + 1) % 5 == 0:
        print(f"  Completed {i+1}/{m} datasets")

# Count successful runs
for model_name in ['A', 'B', 'C', 'D']:
    print(f"  Model {model_name}: {len(mice_results[model_name])}/{m} successful runs")

# ============================================================
# 6. POOL RESULTS USING RUBIN'S RULES
# ============================================================
print("\n" + "=" * 70)
print("STEP 6: Pooling results using Rubin's rules...")
print("=" * 70)

def rubins_rules(results_list, param_name):
    """
    Pool results using Rubin's rules.
    
    For coefficients (on log-odds scale):
    - Pooled estimate: mean of estimates
    - Within-imputation variance: mean of variances
    - Between-imputation variance: variance of estimates
    - Total variance: within + (1 + 1/m) * between
    - 95% CI: pooled ± 1.96 * sqrt(total variance)
    
    For p-values: use the Rubin's rules with t-distribution
    """
    m = len(results_list)
    
    # Extract parameter estimates and variances
    estimates = np.array([r['params'][param_name] for r in results_list])
    variances = np.array([r['cov'].loc[param_name, param_name] if hasattr(r['cov'], 'loc') 
                          else r['cov'][list(r['params'].index).index(param_name), list(r['params'].index).index(param_name)]
                          for r in results_list])
    
    # Pooled estimate (on log-odds scale)
    q_bar = np.mean(estimates)
    
    # Within-imputation variance
    u_bar = np.mean(variances)
    
    # Between-imputation variance
    b = np.var(estimates, ddof=1) if m > 1 else 0
    
    # Total variance
    t_total = u_bar + (1 + 1/m) * b
    
    # Standard error
    se_pooled = np.sqrt(t_total)
    
    # Degrees of freedom (Rubin's formula)
    if b > 0:
        df_old = (m - 1) * (1 + u_bar / ((1 + 1/m) * b))**2
    else:
        df_old = float('inf')
    
    # For large samples, use normal approximation
    # OR and CI
    or_pooled = np.exp(q_bar)
    ci_lower = np.exp(q_bar - 1.96 * se_pooled)
    ci_upper = np.exp(q_bar + 1.96 * se_pooled)
    
    # P-value (two-sided)
    from scipy import stats
    z = q_bar / se_pooled
    p_value = 2 * stats.norm.sf(abs(z))
    
    # Fraction of missing information
    if t_total > 0:
        fmi = (1 + 1/m) * b / t_total
    else:
        fmi = 0
    
    return {
        'estimate': q_bar,
        'se': se_pooled,
        'or': or_pooled,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'p_value': p_value,
        'fmi': fmi,
        'within_var': u_bar,
        'between_var': b
    }

def pool_auc(results_list):
    """Pool AUC values (simple mean, no Rubin's rules needed for AUC)."""
    aucs = [r['auc'] for r in results_list]
    return {
        'mean': np.mean(aucs),
        'std': np.std(aucs, ddof=1),
        'min': np.min(aucs),
        'max': np.max(aucs)
    }

# Pool results for each model
pooled_results = {}
for model_name in ['A', 'B', 'C', 'D']:
    if len(mice_results[model_name]) == 0:
        continue
    
    results_list = mice_results[model_name]
    param_names = list(results_list[0]['params'].index)
    
    pooled = {}
    for param in param_names:
        pooled[param] = rubins_rules(results_list, param)
    
    pooled['auc'] = pool_auc(results_list)
    pooled['llf'] = np.mean([r['llf'] for r in results_list])
    
    pooled_results[model_name] = pooled

# ============================================================
# 7. PRINT COMPARISON TABLE
# ============================================================
print("\n" + "=" * 70)
print("STEP 7: Comparison - Original vs MICE-pooled results")
print("=" * 70)

# Key parameters to compare
key_params = {
    'A': ['const', 'age', 'gender_num', 'sofa_total', 'stage3'],
    'B': ['const', 'age', 'ckd', 'sofa_total', 'stage3'],
    'C': ['const', 'age', 'ckd', 'ckd_stage2', 'ckd_stage3', 'sofa_total', 'stage3'],
    'D': ['const', 'age', 'ckd', 'ckd_stage2', 'ckd_stage3', 
          'sofa_resp', 'sofa_coag', 'sofa_liver', 'sofa_cv', 'sofa_cns', 'sofa_renal', 'stage3']
}

for model_name in ['A', 'B', 'C', 'D']:
    print(f"\n{'='*70}")
    print(f"Model {model_name}")
    print(f"{'='*70}")
    
    if model_name not in pooled_results:
        print("  No results available")
        continue
    
    # AUC comparison
    orig_auc = orig_results[model_name]['auc']
    mice_auc = pooled_results[model_name]['auc']
    print(f"\n  AUC: Original={orig_auc:.4f}, MICE={mice_auc['mean']:.4f} (±{mice_auc['std']:.4f})")
    
    # Key parameters
    print(f"\n  {'Parameter':<20} {'Orig OR':>10} {'Orig CI':>20} {'MICE OR':>10} {'MICE CI':>20} {'FMI':>6} {'ΔOR%':>8}")
    print(f"  {'-'*96}")
    
    params_to_show = key_params[model_name]
    for param in params_to_show:
        if param not in pooled_results[model_name]:
            continue
        
        orig_r = orig_results[model_name]
        orig_idx = list(orig_r['params'].index)
        if param in orig_idx:
            orig_or = orig_r['or'][param]
            orig_ci_l = orig_r['ci_lower'][param]
            orig_ci_u = orig_r['ci_upper'][param]
            orig_str = f"{orig_or:.3f} ({orig_ci_l:.3f}-{orig_ci_u:.3f})"
        else:
            orig_or = None
            orig_str = "N/A"
        
        mice_r = pooled_results[model_name][param]
        mice_or = mice_r['or']
        mice_ci_l = mice_r['ci_lower']
        mice_ci_u = mice_r['ci_upper']
        mice_str = f"{mice_or:.3f} ({mice_ci_l:.3f}-{mice_ci_u:.3f})"
        
        delta = ((mice_or - orig_or) / orig_or * 100) if orig_or and orig_or > 0 else 0
        
        print(f"  {param:<20} {orig_or if orig_or else 0:>10.3f} {orig_str:>20} {mice_or:>10.3f} {mice_str:>20} {mice_r['fmi']:>6.3f} {delta:>7.1f}%")

# ============================================================
# 8. SAVE RESULTS TO JSON
# ============================================================
print("\n" + "=" * 70)
print("STEP 8: Saving results to JSON...")
print("=" * 70)

output = {
    'mice_settings': {
        'm': m,
        'burn_in': burn_in,
        'method': 'PMM (k=5)',
        'seed': 42,
        'missing_vars': missing_vars
    },
    'missing_rates': {
        'sofa_liver': f"{mask_liver.sum()/len(aki_df)*100:.1f}%",
        'sofa_resp': f"{mask_resp.sum()/len(aki_df)*100:.1f}%",
        'sofa_coag': f"{mask_coag.sum()/len(aki_df)*100:.1f}%",
        'sofa_cns': f"{mask_cns.sum()/len(aki_df)*100:.1f}%",
        'sofa_cv': f"{mask_cv.sum()/len(aki_df)*100:.1f}%"
    },
    'original_results': {},
    'mice_pooled_results': {}
}

for model_name in ['A', 'B', 'C', 'D']:
    if model_name in orig_results:
        orig_r = orig_results[model_name]
        output['original_results'][model_name] = {
            'auc': float(orig_r['auc']),
            'llf': float(orig_r['llf']),
            'params': {}
        }
        for param in orig_r['params'].index:
            output['original_results'][model_name]['params'][param] = {
                'or': float(orig_r['or'][param]),
                'ci_lower': float(orig_r['ci_lower'][param]),
                'ci_upper': float(orig_r['ci_upper'][param]),
                'p_value': float(orig_r['p_values'][param])
            }
    
    if model_name in pooled_results:
        mice_r = pooled_results[model_name]
        output['mice_pooled_results'][model_name] = {
            'auc_mean': float(mice_r['auc']['mean']),
            'auc_std': float(mice_r['auc']['std']),
            'llf_mean': float(mice_r['llf']),
            'params': {}
        }
        for param in mice_r:
            if param in ['auc', 'llf']:
                continue
            output['mice_pooled_results'][model_name]['params'][param] = {
                'or': float(mice_r[param]['or']),
                'ci_lower': float(mice_r[param]['ci_lower']),
                'ci_upper': float(mice_r[param]['ci_upper']),
                'p_value': float(mice_r[param]['p_value']),
                'fmi': float(mice_r[param]['fmi']),
                'within_var': float(mice_r[param]['within_var']),
                'between_var': float(mice_r[param]['between_var'])
            }

with open(r'C:\Users\admin\WorkBuddy\2026-07-07-19-40-19\mice_results.json', 'w') as f:
    json.dump(output, f, indent=2, default=str)

print("Results saved to mice_results.json")

# Print summary
print("\n" + "=" * 70)
print("SUMMARY: Key findings")
print("=" * 70)

# Model C - CKD × Stage 3 interaction (the core finding)
if 'C' in pooled_results:
    c_orig = orig_results['C']
    c_mice = pooled_results['C']
    
    # Find the ckd_stage3 parameter
    ckd_s3_key = None
    for key in c_mice:
        if 'ckd_stage3' in key.lower() or 'ckd_s3' in key.lower():
            ckd_s3_key = key
            break
    
    if ckd_s3_key:
        orig_or = c_orig['or'].get(ckd_s3_key, None)
        mice_or = c_mice[ckd_s3_key]['or']
        mice_ci = f"({c_mice[ckd_s3_key]['ci_lower']:.3f}-{c_mice[ckd_s3_key]['ci_upper']:.3f})"
        mice_p = c_mice[ckd_s3_key]['p_value']
        fmi = c_mice[ckd_s3_key]['fmi']
        
        print(f"\n  CKD × Stage 3 interaction (Model C):")
        print(f"    Original: OR={orig_or:.3f}")
        print(f"    MICE:     OR={mice_or:.3f}, 95% CI={mice_ci}, P={mice_p:.4f}")
        print(f"    FMI:      {fmi:.3f}")
        print(f"    Change:   {(mice_or-orig_or)/orig_or*100:+.1f}%")

# Model D AUC
if 'D' in pooled_results:
    d_orig_auc = orig_results['D']['auc']
    d_mice_auc = pooled_results['D']['auc']['mean']
    print(f"\n  Model D AUC:")
    print(f"    Original: {d_orig_auc:.4f}")
    print(f"    MICE:     {d_mice_auc:.4f}")
    print(f"    Change:   {d_mice_auc-d_orig_auc:+.4f}")

# Sofa_total comparison
print(f"\n  SOFA total mean:")
print(f"    Original (missing=0): {orig_total.mean():.3f}")
print(f"    MICE (imputed):       {imputed_datasets[0]['sofa_total'].mean():.3f}")

print("\n" + "=" * 70)
print("MICE analysis complete!")
print("=" * 70)
