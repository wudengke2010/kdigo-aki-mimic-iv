#!/usr/bin/env python3
"""
KDIGO AKI — 5 项新增分析
========================
1. RCS 非线性分析 (Restricted Cubic Splines)
2. 模型比较指标 (AUC, NRI, IDI) + ROC + DCA
3. 脓毒症亚组分析 + 交互项
4. CKD x KDIGO 交互项
5. 校准曲线 (Hosmer-Lemeshow)
"""
import pandas as pd
import numpy as np
import os, json, warnings
from datetime import datetime
from scipy import stats
from scipy.stats import binom

warnings.filterwarnings('ignore')
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score, roc_curve, brier_score_loss

print("=" * 80)
print("KDIGO-AKI: 5 Supplementary Analyses")
print(f"Start: {datetime.now()}")
print("=" * 80)

IN = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"

# ============================================================
# Load and prepare data (same as kdigo_revised_stats.py)
# ============================================================
print("[0] Loading cohort...")
cohort = pd.read_csv(os.path.join(IN, "kdigo_cohort_revised_enriched.csv"), low_memory=False)
print(f"  N = {len(cohort):,}, columns = {len(cohort.columns)}")

# KDIGO stage dummies
for s in [1, 2, 3]:
    cohort[f'kdigo_stage_{s}'] = (cohort['kdigo_stage'] == s).astype(int)

# Demographics
cohort['male_num'] = (cohort['gender'] == 'M').astype(int)

# Severity
cohort['sofa_imputed'] = cohort['sofa_total'].fillna(cohort['sofa_total'].median())
cohort['gcs_imputed'] = cohort['gcs_total'].fillna(15)

base_covars = ['anchor_age', 'male_num']
comorb_vars = sorted([c for c in cohort.columns if c.startswith('com_')])
severity_vars = ['sofa_imputed', 'mech_vent', 'vasopressor', 'ckd_icd']

# log_cr_ratio
cohort['log_cr_ratio'] = np.log(cohort['cr_ratio'].clip(lower=0.1))

# Outcome
y_col = 'hospital_expire_flag'

# ============================================================
# ANALYSIS 1: RCS Nonlinear Analysis
# ============================================================
print("\n" + "=" * 80)
print("[1] RCS Nonlinear Analysis: log(Cr ratio) vs Mortality")
print("=" * 80)

# Model B covariates but replace log_cr_ratio with RCS terms
model_b_covars = ['kdigo_stage_1', 'kdigo_stage_2', 'kdigo_stage_3'] + base_covars + comorb_vars + severity_vars

# For RCS, we use log_cr_ratio (continuous) instead of KDIGO stages
# 4-knot RCS at clinically meaningful Cr ratio thresholds
# Cr ratio: 1.0 (baseline), 1.5 (Stage 1), 2.0 (Stage 2), 3.0 (Stage 3)
knots = np.array([0.0, np.log(1.5), np.log(2.0), np.log(3.0)])
print(f"  RCS knots (clinical thresholds): {knots}")
print(f"  Corresponding Cr ratios: {np.exp(knots)}")

# Create RCS basis functions (4 knots = 3 spline terms)
def make_rcs(x, knots):
    """Create restricted cubic spline basis for k knots (k-1 terms)."""
    k = len(knots)
    n = len(x)
    # First term: x
    basis = [x.values.copy()]
    # Additional terms
    for j in range(2, k):
        term = np.zeros(n)
        t_j = knots[j - 1]
        t_k = knots[-1]
        t_1 = knots[0]
        idx = x.values > t_j
        term[idx] = (x.values[idx] - t_j) ** 3 - (
            (x.values[idx] - t_k) ** 3 * (t_j - t_1) - (x.values[idx] - t_1) ** 3 * (t_j - t_k)
        ) / (t_k - t_1)
        # Actually, use simpler formula
        lam_k = (t_k - t_1)
        term2 = np.where(x.values > t_j,
                         (x.values - t_j) ** 3 -
                         (x.values - t_k) ** 3 * (t_j - t_1) / (t_k - t_1) +
                         (x.values - t_1) ** 3 * (t_j - t_k) / (t_k - t_1),
                         0)
        basis.append(term2)
    return basis

# Actually, let me use a cleaner RCS implementation
def rcs_basis(x, knots):
    """
    Restricted cubic spline basis function.
    Returns matrix with k-1 columns where k = len(knots).
    Based on Harrell's formula.
    """
    x = np.asarray(x, dtype=float)
    k = len(knots)
    n = len(x)
    t = sorted(knots)
    
    # First column: linear term
    basis = np.zeros((n, k - 1))
    basis[:, 0] = x
    
    for j in range(1, k - 1):
        tj = t[j]
        tk_last = t[-1]
        t1 = t[0]
        basis[:, j] = np.where(
            x > tj,
            (x - tj) ** 3 + (t1 - tj) ** 3 - 
            (x - tk_last) ** 3 * (t1 - tj) / (t1 - tk_last) +
            (x - t1) ** 3 * (tj - tk_last) / (t1 - tk_last),
            (t1 - tj) ** 3  # When x <= tj, only the constant part matters
        )
        # Actually the standard formula:
        # When x <= tj, the term is 0 (not the constant)
    
    # Recompute with standard Harrell formula
    basis = np.zeros((n, k - 1))
    basis[:, 0] = x
    
    for j in range(1, k - 1):
        tj = t[j]
        tk_last = t[-1]
        t1 = t[0]
        
        d_j = np.where(x > tj, (x - tj) ** 3, 0)
        d_last = np.where(x > tk_last, (x - tk_last) ** 3, 0)
        d_1 = np.where(x > t1, (x - t1) ** 3, 0)
        
        basis[:, j] = d_j - d_last * (tj - t1) / (tk_last - t1) + d_1 * (tj - tk_last) / (tk_last - t1)
    
    return basis

# Build RCS data
rcs_covars = ['log_cr_ratio'] + base_covars + comorb_vars + severity_vars
rcs_data = cohort[rcs_covars + [y_col]].dropna()
print(f"  N for RCS analysis: {len(rcs_data):,}")

# Create RCS basis
x_rcs = rcs_data['log_cr_ratio'].values
rcs_terms = rcs_basis(x_rcs, knots)
rcs_col_names = ['rcs_1', 'rcs_2', 'rcs_3']  # 4 knots -> 3 terms

# Build design matrix
X_rcs = rcs_data[rcs_covars].copy()
# Replace log_cr_ratio with RCS terms
X_rcs = X_rcs.drop(columns=['log_cr_ratio'])
for i, name in enumerate(rcs_col_names):
    X_rcs[name] = rcs_terms[:, i]

y_rcs = rcs_data[y_col]
X_rcs_c = sm.add_constant(X_rcs)

print("  Fitting RCS logistic model...")
m_rcs = sm.Logit(y_rcs, X_rcs_c).fit(disp=0, maxiter=200)

# Test nonlinearity: compare linear-only vs RCS
# Linear model: log_cr_ratio only
X_lin = rcs_data[rcs_covars].copy()
X_lin_c = sm.add_constant(X_lin)
m_lin = sm.Logit(y_rcs, X_lin_c).fit(disp=0, maxiter=200)

# Nonlinearity test: LRT for rcs_2 and rcs_3 (beyond linear rcs_1)
llr_full = m_rcs.llf
llr_linear = m_lin.llf
chi2_nonlin = 2 * (llr_full - llr_linear)
df_nonlin = 2  # 2 extra parameters (rcs_2, rcs_3)
p_nonlin = 1 - stats.chi2.cdf(chi2_nonlin, df_nonlin)

print(f"\n  Nonlinearity test (LRT):")
print(f"    chi2 = {chi2_nonlin:.2f}, df = {df_nonlin}, p = {p_nonlin:.6f}")
print(f"    -> {'Nonlinear' if p_nonlin < 0.05 else 'Linear'} relationship (p {'<' if p_nonlin < 0.05 else '>'} 0.05)")

# Generate predicted OR curve across Cr ratio range
cr_ratio_range = np.linspace(0.5, 15.0, 200)
log_range = np.log(cr_ratio_range)
rcs_pred = rcs_basis(log_range, knots)

# Get predicted log-odds at reference (Cr ratio = 1.0, i.e., log = 0)
ref_rcs = rcs_basis(np.array([0.0]), knots)

# Calculate OR relative to Cr ratio = 1.0
# For each point, compute the change in log-odds from reference
# due to RCS terms only (holding other covariates at mean)

# Extract RCS coefficients
rcs_coefs = np.array([m_rcs.params[f'rcs_{i+1}'] for i in range(3)])
rcs_se = np.array([m_rcs.bse[f'rcs_{i+1}'] for i in range(3)])

# Compute log-odds difference from reference (Cr ratio=1, log=0)
log_odds_ref = rcs_coefs @ ref_rcs[0]
log_odds_range = rcs_pred @ rcs_coefs
delta_log_odds = log_odds_range - log_odds_ref

# OR = exp(delta_log_odds)
or_curve = np.exp(delta_log_odds)

# Approximate CI using delta method (simplified)
# Variance of linear combination: sum(coef_i * basis_i * coef_j * basis_j * Cov_ij)
# Use covariance matrix of RCS terms
cov_rcs = m_rcs.cov_params().loc[rcs_col_names, rcs_col_names].values

# For each point, variance = basis_vec @ cov_rcs @ basis_vec
# But need to account for reference point subtraction
delta_basis = rcs_pred - ref_rcs[0]  # (200, 3) - (3,) -> broadcasting
var_curve = np.array([db @ cov_rcs @ db for db in delta_basis])
se_curve = np.sqrt(var_curve)

or_low = np.exp(delta_log_odds - 1.96 * se_curve)
or_high = np.exp(delta_log_odds + 1.96 * se_curve)

# Save RCS data for plotting
rcs_output = pd.DataFrame({
    'cr_ratio': cr_ratio_range,
    'OR': or_curve,
    'OR_low': or_low,
    'OR_high': or_high,
})
rcs_output.to_csv(os.path.join(IN, 'rcs_results.csv'), index=False)

print(f"\n  RCS OR curve saved to rcs_results.csv")
print(f"  Sample ORs:")
for cr in [1.0, 1.5, 2.0, 3.0, 5.0, 10.0]:
    idx = np.argmin(np.abs(cr_ratio_range - cr))
    print(f"    Cr ratio {cr:.1f}: OR = {or_curve[idx]:.2f} ({or_low[idx]:.2f}-{or_high[idx]:.2f})")

# ============================================================
# ANALYSIS 2: Model Comparison (AUC, NRI, IDI) + ROC + DCA
# ============================================================
print("\n" + "=" * 80)
print("[2] Model Comparison: AUC, NRI, IDI, ROC, DCA")
print("=" * 80)

# Re-fit Models A, B, C (same as stats script)
# Model A
model_a_covars = ['kdigo_stage_1', 'kdigo_stage_2', 'kdigo_stage_3'] + base_covars + comorb_vars
model_a_data = cohort[model_a_covars + [y_col]].dropna()
X_a = model_a_data[model_a_covars].copy()
y_a = model_a_data[y_col]
X_a_c = sm.add_constant(X_a)
m_a = sm.Logit(y_a, X_a_c).fit(disp=0, maxiter=200)

# Model B
model_b_covars = ['kdigo_stage_1', 'kdigo_stage_2', 'kdigo_stage_3'] + base_covars + comorb_vars + severity_vars
model_b_data = cohort[model_b_covars + [y_col]].dropna()
X_b = model_b_data[model_b_covars].copy()
y_b = model_b_data[y_col]
X_b_c = sm.add_constant(X_b)
m_b = sm.Logit(y_b, X_b_c).fit(disp=0, maxiter=200)

# Model C
model_c_covars = ['log_cr_ratio'] + base_covars + comorb_vars + severity_vars
model_c_data = cohort[model_c_covars + [y_col]].dropna()
X_c = model_c_data[model_c_covars].copy()
y_c = model_c_data[y_col]
X_c_c = sm.add_constant(X_c)
m_c = sm.Logit(y_c, X_c_c).fit(disp=0, maxiter=200)

# Predicted probabilities
pred_a = m_a.predict(X_a_c)
pred_b = m_b.predict(X_b_c)
pred_c = m_c.predict(X_c_c)

# AUC
auc_a = roc_auc_score(y_a, pred_a)
auc_b = roc_auc_score(y_b, pred_b)
auc_c = roc_auc_score(y_c, pred_c)

print(f"  AUC (Model A): {auc_a:.3f} (95% CI computed below)")
print(f"  AUC (Model B): {auc_b:.3f}")
print(f"  AUC (Model C): {auc_c:.3f}")

# AUC 95% CI via bootstrap (DeLong approximation using Hanley-McNeil)
def auc_ci(y_true, y_pred, alpha=0.05):
    """Hanley-McNeil AUC confidence interval."""
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    n1 = np.sum(y_true == 1)
    n0 = np.sum(y_true == 0)
    auc = roc_auc_score(y_true, y_pred)
    
    # Hanley-McNeil
    Q1 = auc / (2 - auc)
    Q2 = 2 * auc**2 / (1 + auc)
    se = np.sqrt((auc * (1 - auc) + (n1 - 1) * (Q1 - auc**2) + (n0 - 1) * (Q2 - auc**2)) / (n1 * n0))
    z = stats.norm.ppf(1 - alpha / 2)
    ci_low = auc - z * se
    ci_high = auc + z * se
    return auc, ci_low, ci_high, se

auc_a, ci_a_low, ci_a_high, se_a = auc_ci(y_a, pred_a)
auc_b, ci_b_low, ci_b_high, se_b = auc_ci(y_b, pred_b)
auc_c, ci_c_low, ci_c_high, se_c = auc_ci(y_c, pred_c)

print(f"\n  AUC with 95% CI:")
print(f"    Model A: {auc_a:.3f} ({ci_a_low:.3f}-{ci_a_high:.3f}), SE={se_a:.4f}")
print(f"    Model B: {auc_b:.3f} ({ci_b_low:.3f}-{ci_b_high:.3f}), SE={se_b:.4f}")
print(f"    Model C: {auc_c:.3f} ({ci_c_low:.3f}-{ci_c_high:.3f}), SE={se_c:.4f}")

# Brier score
brier_a = brier_score_loss(y_a, pred_a)
brier_b = brier_score_loss(y_b, pred_b)
brier_c = brier_score_loss(y_c, pred_c)
print(f"\n  Brier score:")
print(f"    Model A: {brier_a:.4f}")
print(f"    Model B: {brier_b:.4f}")
print(f"    Model C: {brier_c:.4f}")

# NRI and IDI: Model A -> Model B
# Use common sample (intersection of A and B)
common_idx = model_a_data.index.intersection(model_b_data.index)
y_common = cohort.loc[common_idx, y_col]

# Get predictions on common sample
pred_a_common = m_a.predict(sm.add_constant(model_a_data.loc[common_idx, model_a_covars]))
pred_b_common = m_b.predict(sm.add_constant(model_b_data.loc[common_idx, model_b_covars]))

# NRI (Continuous)
def calculate_nri_idi(y_true, pred_old, pred_new):
    """Calculate NRI and IDI."""
    y_true = np.asarray(y_true)
    pred_old = np.asarray(pred_old)
    pred_new = np.asarray(pred_new)
    
    events = y_true == 1
    non_events = y_true == 0
    
    # NRI
    # For events: proportion with increased probability
    up_events = np.sum(pred_new[events] > pred_old[events]) / np.sum(events)
    down_events = np.sum(pred_new[events] < pred_old[events]) / np.sum(events)
    
    # For non-events: proportion with decreased probability
    up_nonevents = np.sum(pred_new[non_events] > pred_old[non_events]) / np.sum(non_events)
    down_nonevents = np.sum(pred_new[non_events] < pred_old[non_events]) / np.sum(non_events)
    
    nri = (up_events - down_events) + (down_nonevents - up_nonevents)
    
    # IDI
    mean_pred_new_events = np.mean(pred_new[events])
    mean_pred_old_events = np.mean(pred_old[events])
    mean_pred_new_nonevents = np.mean(pred_new[non_events])
    mean_pred_old_nonevents = np.mean(pred_old[non_events])
    
    integrated_discrimination_events = mean_pred_new_events - mean_pred_old_events
    integrated_discrimination_nonevents = mean_pred_new_nonevents - mean_pred_old_nonevents
    
    idi = integrated_discrimination_events - integrated_discrimination_nonevents
    
    return {
        'NRI': nri,
        'NRI_events': up_events - down_events,
        'NRI_nonevents': down_nonevents - up_nonevents,
        'IDI': idi,
        'IDI_events': integrated_discrimination_events,
        'IDI_nonevents': integrated_discrimination_nonevents,
    }

nri_idi_ab = calculate_nri_idi(y_common, pred_a_common, pred_b_common)
print(f"\n  NRI/IDI (Model A -> Model B):")
print(f"    NRI = {nri_idi_ab['NRI']:.4f} (events: {nri_idi_ab['NRI_events']:.4f}, non-events: {nri_idi_ab['NRI_nonevents']:.4f})")
print(f"    IDI = {nri_idi_ab['IDI']:.4f} (events: {nri_idi_ab['IDI_events']:.4f}, non-events: {nri_idi_ab['IDI_nonevents']:.4f})")

# NRI/IDI: Model B -> RCS model
# Need to get RCS predictions on same sample
rcs_common_idx = rcs_data.index.intersection(model_b_data.index)
y_rcs_common = cohort.loc[rcs_common_idx, y_col]

pred_b_rcs_common = m_b.predict(sm.add_constant(model_b_data.loc[rcs_common_idx, model_b_covars]))
# RCS predictions
rcs_common_data = rcs_data.loc[rcs_common_idx]
x_rcs_common = rcs_common_data['log_cr_ratio'].values
rcs_terms_common = rcs_basis(x_rcs_common, knots)
X_rcs_common = rcs_common_data[rcs_covars].drop(columns=['log_cr_ratio']).copy()
for i, name in enumerate(rcs_col_names):
    X_rcs_common[name] = rcs_terms_common[:, i]
pred_rcs_common = m_rcs.predict(sm.add_constant(X_rcs_common))

# Compare Model B (KDIGO stages) vs RCS (continuous)
# Actually, let's compare Model C (linear log_cr_ratio) vs RCS
# The key comparison is: does RCS (nonlinear) improve over Model C (linear)?
nri_idi_c_rcs = calculate_nri_idi(y_rcs_common, 
    m_c.predict(sm.add_constant(model_c_data.loc[rcs_common_idx, model_c_covars])),
    pred_rcs_common)
print(f"\n  NRI/IDI (Model C linear -> RCS nonlinear):")
print(f"    NRI = {nri_idi_c_rcs['NRI']:.4f}")
print(f"    IDI = {nri_idi_c_rcs['IDI']:.4f}")

# Save ROC data
fpr_a, tpr_a, _ = roc_curve(y_a, pred_a)
fpr_b, tpr_b, _ = roc_curve(y_b, pred_b)
fpr_c, tpr_c, _ = roc_curve(y_c, pred_c)

roc_data = pd.DataFrame({
    'fpr_a': fpr_a, 'tpr_a': tpr_a,
})
# Pad to same length
max_len = max(len(fpr_a), len(fpr_b), len(fpr_c))
roc_df = pd.DataFrame({
    'fpr_a': np.pad(fpr_a, (0, max_len - len(fpr_a)), constant_values=np.nan),
    'tpr_a': np.pad(tpr_a, (0, max_len - len(tpr_a)), constant_values=np.nan),
    'fpr_b': np.pad(fpr_b, (0, max_len - len(fpr_b)), constant_values=np.nan),
    'tpr_b': np.pad(tpr_b, (0, max_len - len(tpr_b)), constant_values=np.nan),
    'fpr_c': np.pad(fpr_c, (0, max_len - len(fpr_c)), constant_values=np.nan),
    'tpr_c': np.pad(tpr_c, (0, max_len - len(tpr_c)), constant_values=np.nan),
})
roc_df.to_csv(os.path.join(IN, 'roc_comparison.csv'), index=False)

# DCA: Decision Curve Analysis
def decision_curve(y_true, pred_prob, thresholds):
    """Calculate net benefit at each threshold."""
    y_true = np.asarray(y_true)
    pred_prob = np.asarray(pred_prob)
    n = len(y_true)
    nb_values = []
    for t in thresholds:
        # Predicted positive: pred >= threshold
        tp = np.sum((pred_prob >= t) & (y_true == 1))
        fp = np.sum((pred_prob >= t) & (y_true == 0))
        nb = tp / n - fp / n * t / (1 - t)
        nb_values.append(nb)
    return np.array(nb_values)

thresholds = np.arange(0.01, 0.99, 0.01)
nb_a = decision_curve(y_a, pred_a, thresholds)
nb_b = decision_curve(y_b, pred_b, thresholds)
nb_c = decision_curve(y_c, pred_c, thresholds)

# Treat-all: predict everyone positive
nb_treat_all = np.array([np.mean(y_a) - (1 - np.mean(y_a)) * t / (1 - t) for t in thresholds])
# Treat-none: 0
nb_treat_none = np.zeros(len(thresholds))

dca_df = pd.DataFrame({
    'threshold': thresholds,
    'nb_model_a': nb_a,
    'nb_model_b': nb_b,
    'nb_model_c': nb_c,
    'nb_treat_all': nb_treat_all,
    'nb_treat_none': nb_treat_none,
})
dca_df.to_csv(os.path.join(IN, 'dca_results.csv'), index=False)

# Save model comparison summary
model_comparison = {
    'AUC': {
        'Model_A': {'auc': round(auc_a, 3), 'ci_low': round(ci_a_low, 3), 'ci_high': round(ci_a_high, 3)},
        'Model_B': {'auc': round(auc_b, 3), 'ci_low': round(ci_b_low, 3), 'ci_high': round(ci_b_high, 3)},
        'Model_C': {'auc': round(auc_c, 3), 'ci_low': round(ci_c_low, 3), 'ci_high': round(ci_c_high, 3)},
    },
    'Brier': {
        'Model_A': round(brier_a, 4),
        'Model_B': round(brier_b, 4),
        'Model_C': round(brier_c, 4),
    },
    'NRI_IDI_A_to_B': {
        'NRI': round(nri_idi_ab['NRI'], 4),
        'NRI_events': round(nri_idi_ab['NRI_events'], 4),
        'NRI_nonevents': round(nri_idi_ab['NRI_nonevents'], 4),
        'IDI': round(nri_idi_ab['IDI'], 4),
        'IDI_events': round(nri_idi_ab['IDI_events'], 4),
        'IDI_nonevents': round(nri_idi_ab['IDI_nonevents'], 4),
    },
    'NRI_IDI_C_to_RCS': {
        'NRI': round(nri_idi_c_rcs['NRI'], 4),
        'IDI': round(nri_idi_c_rcs['IDI'], 4),
    },
    'RCS_nonlinearity': {
        'chi2': round(chi2_nonlin, 2),
        'df': df_nonlin,
        'p_value': round(p_nonlin, 6),
    },
}

with open(os.path.join(IN, 'model_comparison.json'), 'w') as f:
    json.dump(model_comparison, f, indent=2)

print(f"\n  Model comparison saved to model_comparison.json")

# ============================================================
# ANALYSIS 3: Sepsis Subgroup Analysis + Interaction
# ============================================================
print("\n" + "=" * 80)
print("[3] Sepsis Subgroup Analysis + KDIGO x Sepsis Interaction")
print("=" * 80)

sepsis_col = 'com_Sepsis'
print(f"  Sepsis prevalence: {cohort[sepsis_col].sum():,} / {len(cohort):,} ({cohort[sepsis_col].mean()*100:.1f}%)")

# Subgroup analysis: Model B in sepsis vs non-sepsis
# Remove com_Sepsis from covariates (zero variance within subgroups)
sepsis_subgroup_covars = [c for c in model_b_covars if c != 'com_Sepsis']

sepsis_results = {}
for sepsis_status, label in [(1, 'Sepsis'), (0, 'Non-sepsis')]:
    sub = cohort[cohort[sepsis_col] == sepsis_status]
    sub_data = sub[sepsis_subgroup_covars + [y_col]].dropna()
    X_sub = sub_data[sepsis_subgroup_covars].copy()
    y_sub = sub_data[y_col]
    X_sub_c = sm.add_constant(X_sub)
    
    try:
        m_sub = sm.Logit(y_sub, X_sub_c).fit(disp=0, maxiter=200)
        n_sub = len(sub_data)
        auc_sub = roc_auc_score(y_sub, m_sub.predict(X_sub_c))
        
        print(f"\n  {label} (n={n_sub:,}, mortality={y_sub.mean()*100:.1f}%):")
        print(f"    AUC = {auc_sub:.3f}")
        
        kdigo_results = {}
        for s in [1, 2, 3]:
            col = f'kdigo_stage_{s}'
            idx = sepsis_subgroup_covars.index(col)
            coef = m_sub.params[idx + 1]
            se = m_sub.bse[idx + 1]
            or_val = np.exp(coef)
            ci_low = np.exp(coef - 1.96 * se)
            ci_high = np.exp(coef + 1.96 * se)
            p_val = m_sub.pvalues[idx + 1]
            kdigo_results[s] = {
                'OR': round(float(or_val), 3),
                'CI_low': round(float(ci_low), 3),
                'CI_high': round(float(ci_high), 3),
                'p': round(float(p_val), 6),
            }
            print(f"    Stage {s}: OR={or_val:.3f} ({ci_low:.3f}-{ci_high:.3f}), p={p_val:.4f}")
        
        sepsis_results[label] = {
            'n': n_sub,
            'mortality_pct': round(float(y_sub.mean() * 100), 1),
            'AUC': round(float(auc_sub), 3),
            'kdigo_or': kdigo_results,
        }
    except Exception as e:
        print(f"  {label}: Error - {e}")

# Interaction test: KDIGO x Sepsis
print(f"\n  Interaction test: KDIGO x Sepsis")
int_covars = model_b_covars.copy()
# Add interaction terms
for s in [1, 2, 3]:
    int_term = f'kdigo_{s}_x_sepsis'
    cohort[int_term] = cohort[f'kdigo_stage_{s}'] * cohort[sepsis_col]
    int_covars.append(int_term)

int_data = cohort[int_covars + [y_col]].dropna()
X_int = int_data[int_covars].copy()
y_int = int_data[y_col]
X_int_c = sm.add_constant(X_int)
m_int = sm.Logit(y_int, X_int_c).fit(disp=0, maxiter=200)

# LRT: compare with model without interaction
llr_int = m_int.llf
llr_base = m_b.llf
chi2_int = 2 * (llr_int - llr_base)
df_int = 3  # 3 interaction terms
p_int = 1 - stats.chi2.cdf(chi2_int, df_int)

print(f"    LRT: chi2 = {chi2_int:.2f}, df = {df_int}, p = {p_int:.6f}")

for s in [1, 2, 3]:
    col = f'kdigo_{s}_x_sepsis'
    idx = int_covars.index(col)
    coef = m_int.params[idx + 1]
    se = m_int.bse[idx + 1]
    or_val = np.exp(coef)
    ci_low = np.exp(coef - 1.96 * se)
    ci_high = np.exp(coef + 1.96 * se)
    p_val = m_int.pvalues[idx + 1]
    print(f"    Stage {s} x Sepsis: OR={or_val:.3f} ({ci_low:.3f}-{ci_high:.3f}), p={p_val:.4f}")

sepsis_results['interaction'] = {
    'chi2': round(float(chi2_int), 2),
    'df': df_int,
    'p_value': round(float(p_int), 6),
    'stage1_interaction_OR': round(float(np.exp(m_int.params[int_covars.index('kdigo_1_x_sepsis') + 1])), 3),
    'stage2_interaction_OR': round(float(np.exp(m_int.params[int_covars.index('kdigo_2_x_sepsis') + 1])), 3),
    'stage3_interaction_OR': round(float(np.exp(m_int.params[int_covars.index('kdigo_3_x_sepsis') + 1])), 3),
}

with open(os.path.join(IN, 'sepsis_subgroup.json'), 'w') as f:
    json.dump(sepsis_results, f, indent=2)

print(f"\n  Sepsis subgroup results saved to sepsis_subgroup.json")

# ============================================================
# ANALYSIS 4: CKD x KDIGO Interaction (AKI Paradox Formal Test)
# ============================================================
print("\n" + "=" * 80)
print("[4] CKD x KDIGO Interaction (Formal AKI Paradox Test)")
print("=" * 80)

ckd_int_covars = model_b_covars.copy()
# Add interaction terms
for s in [1, 2, 3]:
    int_term = f'kdigo_{s}_x_ckd'
    cohort[int_term] = cohort[f'kdigo_stage_{s}'] * cohort['ckd_icd']
    ckd_int_covars.append(int_term)

ckd_int_data = cohort[ckd_int_covars + [y_col]].dropna()
X_ckd_int = ckd_int_data[ckd_int_covars].copy()
y_ckd_int = ckd_int_data[y_col]
X_ckd_int_c = sm.add_constant(X_ckd_int)
m_ckd_int = sm.Logit(y_ckd_int, X_ckd_int_c).fit(disp=0, maxiter=200)

# LRT
llr_ckd_int = m_ckd_int.llf
llr_ckd_base = m_b.llf
chi2_ckd = 2 * (llr_ckd_int - llr_ckd_base)
df_ckd = 3
p_ckd = 1 - stats.chi2.cdf(chi2_ckd, df_ckd)

print(f"  Interaction test: KDIGO x CKD")
print(f"    LRT: chi2 = {chi2_ckd:.2f}, df = {df_ckd}, p = {p_ckd:.6f}")

ckd_int_results = {'interaction': {
    'chi2': round(float(chi2_ckd), 2),
    'df': df_ckd,
    'p_value': round(float(p_ckd), 6),
}}

for s in [1, 2, 3]:
    col = f'kdigo_{s}_x_ckd'
    idx = ckd_int_covars.index(col)
    coef = m_ckd_int.params[idx + 1]
    se = m_ckd_int.bse[idx + 1]
    or_val = np.exp(coef)
    ci_low = np.exp(coef - 1.96 * se)
    ci_high = np.exp(coef + 1.96 * se)
    p_val = m_ckd_int.pvalues[idx + 1]
    
    print(f"    Stage {s} x CKD: OR={or_val:.3f} ({ci_low:.3f}-{ci_high:.3f}), p={p_val:.4f}")
    
    ckd_int_results[f'stage{s}_interaction'] = {
        'OR': round(float(or_val), 3),
        'CI_low': round(float(ci_low), 3),
        'CI_high': round(float(ci_high), 3),
        'p': round(float(p_val), 6),
    }

# Also test with lab-CKD
print(f"\n  Interaction test: KDIGO x Lab-CKD (eGFR <60)")
lab_ckd_covars = ['kdigo_stage_1', 'kdigo_stage_2', 'kdigo_stage_3'] + base_covars + comorb_vars + ['sofa_imputed', 'mech_vent', 'vasopressor', 'ckd_lab']
for s in [1, 2, 3]:
    int_term = f'kdigo_{s}_x_labckd'
    cohort[int_term] = cohort[f'kdigo_stage_{s}'] * cohort['ckd_lab']
    lab_ckd_covars.append(int_term)

lab_ckd_int_data = cohort[lab_ckd_covars + [y_col]].dropna()
X_lab_ckd_int = lab_ckd_int_data[lab_ckd_covars].copy()
y_lab_ckd_int = lab_ckd_int_data[y_col]
X_lab_ckd_int_c = sm.add_constant(X_lab_ckd_int)
m_lab_ckd_int = sm.Logit(y_lab_ckd_int, X_lab_ckd_int_c).fit(disp=0, maxiter=200)

# Base model with lab CKD (no interaction)
lab_ckd_base_covars = ['kdigo_stage_1', 'kdigo_stage_2', 'kdigo_stage_3'] + base_covars + comorb_vars + ['sofa_imputed', 'mech_vent', 'vasopressor', 'ckd_lab']
lab_ckd_base_data = cohort[lab_ckd_base_covars + [y_col]].dropna()
X_lab_base = lab_ckd_base_data[lab_ckd_base_covars].copy()
y_lab_base = lab_ckd_base_data[y_col]
X_lab_base_c = sm.add_constant(X_lab_base)
m_lab_base = sm.Logit(y_lab_base, X_lab_base_c).fit(disp=0, maxiter=200)

llr_lab_int = m_lab_ckd_int.llf
llr_lab_base = m_lab_base.llf
chi2_lab = 2 * (llr_lab_int - llr_lab_base)
p_lab = 1 - stats.chi2.cdf(chi2_lab, 3)

print(f"    LRT: chi2 = {chi2_lab:.2f}, df = 3, p = {p_lab:.6f}")

for s in [1, 2, 3]:
    col = f'kdigo_{s}_x_labckd'
    idx = lab_ckd_covars.index(col)
    coef = m_lab_ckd_int.params[idx + 1]
    se = m_lab_ckd_int.bse[idx + 1]
    or_val = np.exp(coef)
    ci_low = np.exp(coef - 1.96 * se)
    ci_high = np.exp(coef + 1.96 * se)
    p_val = m_lab_ckd_int.pvalues[idx + 1]
    print(f"    Stage {s} x Lab-CKD: OR={or_val:.3f} ({ci_low:.3f}-{ci_high:.3f}), p={p_val:.4f}")
    
    ckd_int_results[f'stage{s}_labckd_interaction'] = {
        'OR': round(float(or_val), 3),
        'CI_low': round(float(ci_low), 3),
        'CI_high': round(float(ci_high), 3),
        'p': round(float(p_val), 6),
    }

ckd_int_results['labckd_interaction'] = {
    'chi2': round(float(chi2_lab), 2),
    'df': 3,
    'p_value': round(float(p_lab), 6),
}

with open(os.path.join(IN, 'ckd_interaction.json'), 'w') as f:
    json.dump(ckd_int_results, f, indent=2)

print(f"\n  CKD interaction results saved to ckd_interaction.json")

# ============================================================
# ANALYSIS 5: Calibration Curve (Hosmer-Lemeshow)
# ============================================================
print("\n" + "=" * 80)
print("[5] Calibration Curve (Hosmer-Lemeshow)")
print("=" * 80)

def hosmer_lemeshow(y_true, y_pred, n_groups=10):
    """Hosmer-Lemeshow test and calibration data."""
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    n = len(y_true)
    
    # Sort by predicted probability
    order = np.argsort(y_pred)
    y_true_sorted = y_true[order]
    y_pred_sorted = y_pred[order]
    
    # Divide into groups
    groups = np.array_split(np.arange(n), n_groups)
    
    cal_data = []
    obs_list = []
    exp_list = []
    
    for i, idx in enumerate(groups):
        obs = np.mean(y_true_sorted[idx])
        exp = np.mean(y_pred_sorted[idx])
        n_g = len(idx)
        
        # CI for observed proportion
        ci_low = binom.interval(0.95, n_g, obs)[0] / n_g if n_g > 0 else 0
        ci_high = binom.interval(0.95, n_g, obs)[1] / n_g if n_g > 0 else 0
        
        cal_data.append({
            'group': i + 1,
            'n': n_g,
            'observed': round(float(obs), 4),
            'expected': round(float(exp), 4),
            'obs_ci_low': round(float(ci_low), 4),
            'obs_ci_high': round(float(ci_high), 4),
            'mean_pred': round(float(exp), 4),
        })
        obs_list.append(obs * n_g)
        exp_list.append(exp * n_g)
    
    # Hosmer-Lemeshow chi-square
    obs_arr = np.array(obs_list)
    exp_arr = np.array(exp_list)
    # Use groups of approximately equal size
    hl = np.sum((obs_arr - exp_arr) ** 2 / (exp_arr * (1 - exp_arr / [len(g) for g in groups] + 1e-10)))
    
    # Simpler HL formula
    hl_stat = 0
    for i, idx in enumerate(groups):
        n_g = len(idx)
        o1 = obs_list[i]
        e1 = exp_list[i]
        o0 = n_g - o1
        e0 = n_g - e1
        if e1 > 0 and e0 > 0:
            hl_stat += (o1 - e1) ** 2 / e1 + (o0 - e0) ** 2 / e0
    
    p_hl = 1 - stats.chi2.cdf(hl_stat, n_groups - 2)
    
    return cal_data, hl_stat, p_hl

# Calibration for Model B
cal_b, hl_b, p_hl_b = hosmer_lemeshow(y_b, pred_b, n_groups=10)
print(f"  Model B Hosmer-Lemeshow: chi2 = {hl_b:.2f}, df = 8, p = {p_hl_b:.4f}")
print(f"    -> {'Good fit' if p_hl_b > 0.05 else 'Poor fit'} (p {'>' if p_hl_b > 0.05 else '<='} 0.05)")

# Calibration for Model A
cal_a, hl_a, p_hl_a = hosmer_lemeshow(y_a, pred_a, n_groups=10)
print(f"  Model A Hosmer-Lemeshow: chi2 = {hl_a:.2f}, p = {p_hl_a:.4f}")

# Save calibration data
cal_df = pd.DataFrame(cal_b)
cal_df.to_csv(os.path.join(IN, 'calibration_model_b.csv'), index=False)

calibration_results = {
    'Model_B': {'HL_chi2': round(float(hl_b), 2), 'p_value': round(float(p_hl_b), 4),
                'calibration_data': cal_b},
    'Model_A': {'HL_chi2': round(float(hl_a), 2), 'p_value': round(float(p_hl_a), 4)},
}

with open(os.path.join(IN, 'calibration_results.json'), 'w') as f:
    json.dump(calibration_results, f, indent=2)

print(f"\n  Calibration data saved to calibration_model_b.csv")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 80)
print("ALL 5 ANALYSES COMPLETE")
print("=" * 80)
print(f"\nSummary:")
print(f"  1. RCS: P_nonlinearity = {p_nonlin:.6f} ({'significant' if p_nonlin < 0.05 else 'not significant'})")
print(f"  2. AUC: A={auc_a:.3f}, B={auc_b:.3f}, C={auc_c:.3f}")
print(f"     NRI(A->B)={nri_idi_ab['NRI']:.4f}, IDI(A->B)={nri_idi_ab['IDI']:.4f}")
print(f"  3. Sepsis: KDIGO x Sepsis p = {p_int:.6f}")
print(f"     Sepsis Stage 3 OR = {sepsis_results.get('Sepsis', {}).get('kdigo_or', {}).get(3, {}).get('OR', 'N/A')}")
print(f"     Non-sepsis Stage 3 OR = {sepsis_results.get('Non-sepsis', {}).get('kdigo_or', {}).get(3, {}).get('OR', 'N/A')}")
print(f"  4. CKD x KDIGO: p = {p_ckd:.6f}")
print(f"     Stage 3 x CKD OR = {ckd_int_results['stage3_interaction']['OR']}")
print(f"  5. Calibration: HL p = {p_hl_b:.4f} ({'good fit' if p_hl_b > 0.05 else 'poor fit'})")
print(f"\nEnd: {datetime.now()}")
