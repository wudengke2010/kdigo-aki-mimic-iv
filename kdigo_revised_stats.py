#!/usr/bin/env python3
"""
KDIGO AKI 分期 — 修订版统计分析
================================
审稿修订内容:
  Fix #2: 移除 KDIGO stage 与 Cr ratio 共线性 — 分别建模
  Fix #3: 加入 SOFA/GCS/Vent/Vasopressor 严重度评分
  Fix #4: RRT-stratified Stage 3 AKI paradox 分析
  Fix #5: MDRD baseline 敏感性分析
  Fix #6: Lab-based CKD 敏感性分析
  Fix #8: 报告全部协变量 (Cox + Logistic)
  Fix #9: 患者筛选流程图数据
  Fix #10: 缺失数据报告
  Fix #12: Bonferroni 校正 pairwise log-rank
  Fix #13: 剂量-反应表含 95% CI
"""
import pandas as pd
import numpy as np
import os, json, warnings
from datetime import datetime

warnings.filterwarnings('ignore')

print("=" * 80)
print("KDIGO-AKI REVISED Statistical Analysis")
print(f"Start: {datetime.now()}")
print("=" * 80)

IN = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"

# ============================================================
# 1. Load data
# ============================================================
print("[1] Loading revised cohort...")
cohort_path = os.path.join(IN, "kdigo_cohort_revised_enriched.csv")
if not os.path.exists(cohort_path):
    cohort_path = os.path.join(IN, "kdigo_cohort_revised.csv")

cohort = pd.read_csv(cohort_path, low_memory=False)
print(f"  N = {len(cohort):,}, columns = {len(cohort.columns)}")

# ============================================================
# 2. Descriptive stats
# ============================================================
print("\n[2] Descriptive statistics by KDIGO stage...")
for s in [0, 1, 2, 3]:
    sub = cohort[cohort['kdigo_stage'] == s]
    n = len(sub)
    deaths = sub['hospital_expire_flag'].sum()
    rate = deaths / n * 100 if n > 0 else 0
    print(f"  Stage {s}: {n:,} stays, {deaths:,} deaths ({rate:.2f}%)")

# ============================================================
# 3. Logistic Regression — Multiple Models
# ============================================================
print("\n[3] Logistic Regression — Model comparison...")
import statsmodels.api as sm

# Prepare KDIGO stage dummies
for s in [1, 2, 3]:
    cohort[f'kdigo_stage_{s}'] = (cohort['kdigo_stage'] == s).astype(int)

# Prepare demographics
cohort['male_num'] = (cohort['gender'] == 'M').astype(int) if 'gender' in cohort.columns else cohort.get('male', 0)

# Severity
cohort['has_sofa'] = cohort['sofa_total'].notna().astype(int)
cohort['sofa_imputed'] = cohort['sofa_total'].fillna(cohort['sofa_total'].median())
cohort['gcs_imputed'] = cohort['gcs_total'].fillna(15)

# Demographics + comorbidity covariates (always included)
base_covars = ['anchor_age', 'male_num']
comorb_vars = sorted([c for c in cohort.columns if c.startswith('com_')])
severity_vars = ['sofa_imputed', 'mech_vent', 'vasopressor', 'ckd_icd']

# === MODEL A: KDIGO stages only (NO creatinine ratio) ===
print("\n--- MODEL A: KDIGO stages + demographics + comorbidities ---")
model_a_covars = ['kdigo_stage_1', 'kdigo_stage_2', 'kdigo_stage_3'] + base_covars + comorb_vars

model_a_data = cohort[model_a_covars + ['hospital_expire_flag']].dropna()
X_a = model_a_data[model_a_covars].copy()
y_a = model_a_data['hospital_expire_flag']

print(f"  N = {len(model_a_data):,}")

X_a_c = sm.add_constant(X_a)
m_a = sm.Logit(y_a, X_a_c).fit(disp=0, maxiter=200)

# Extract KDIGO stage ORs
print(f"  {'Variable':<40} {'OR':>8} {'95% CI':>20} {'p-value':>12}")
print(f"  {'-'*80}")
a_results = []
for i, col in enumerate(model_a_covars):
    coef = m_a.params[i + 1]
    se = m_a.bse[i + 1]
    or_val = np.exp(coef)
    ci_low = np.exp(coef - 1.96 * se)
    ci_high = np.exp(coef + 1.96 * se)
    p_val = m_a.pvalues[i + 1]
    sig = '***' if p_val < 0.001 else ('**' if p_val < 0.01 else ('*' if p_val < 0.05 else ''))
    if col.startswith('kdigo') or col in ['sofa_imputed', 'mech_vent', 'vasopressor', 'ckd_icd'] + base_covars or p_val < 0.001:
        print(f"  {col:<40} {or_val:>8.3f} ({ci_low:.3f}-{ci_high:.3f}) {p_val:>10.4f} {sig}")
    a_results.append({'model': 'A', 'variable': col, 'OR': round(float(or_val),3),
                      'CI_low': round(float(ci_low),3), 'CI_high': round(float(ci_high),3),
                      'p_value': round(float(p_val),6)})

# === MODEL B: KDIGO + severity scores ===
print("\n--- MODEL B: KDIGO stages + severity scores ---")
model_b_covars = ['kdigo_stage_1', 'kdigo_stage_2', 'kdigo_stage_3'] + base_covars + comorb_vars + severity_vars

model_b_data = cohort[model_b_covars + ['hospital_expire_flag']].dropna()
X_b = model_b_data[model_b_covars].copy()
y_b = model_b_data['hospital_expire_flag']

print(f"  N = {len(model_b_data):,}")

X_b_c = sm.add_constant(X_b)
m_b = sm.Logit(y_b, X_b_c).fit(disp=0, maxiter=200)

print(f"  {'Variable':<40} {'OR':>8} {'95% CI':>20} {'p-value':>12}")
print(f"  {'-'*80}")
b_results = []
for i, col in enumerate(model_b_covars):
    coef = m_b.params[i + 1]
    se = m_b.bse[i + 1]
    or_val = np.exp(coef)
    ci_low = np.exp(coef - 1.96 * se)
    ci_high = np.exp(coef + 1.96 * se)
    p_val = m_b.pvalues[i + 1]
    sig = '***' if p_val < 0.001 else ('**' if p_val < 0.01 else ('*' if p_val < 0.05 else ''))
    if col.startswith('kdigo') or col in severity_vars + base_covars or p_val < 0.001:
        print(f"  {col:<40} {or_val:>8.3f} ({ci_low:.3f}-{ci_high:.3f}) {p_val:>10.4f} {sig}")
    b_results.append({'model': 'B', 'variable': col, 'OR': round(float(or_val),3),
                      'CI_low': round(float(ci_low),3), 'CI_high': round(float(ci_high),3),
                      'p_value': round(float(p_val),6)})

# Compare AIC
print(f"\n  Model A (no severity) AIC: {m_a.aic:.1f}")
print(f"  Model B (with severity) AIC: {m_b.aic:.1f}")
print(f"  Delta AIC: {m_a.aic - m_b.aic:.1f}")

# === MODEL C: Creatinine ratio only (no KDIGO) ===
print("\n--- MODEL C: Creatinine ratio only (sensitivity) ---")
cohort['log_cr_ratio'] = np.log(cohort['cr_ratio'].clip(lower=0.1))

model_c_covars = ['log_cr_ratio'] + base_covars + comorb_vars + severity_vars
model_c_data = cohort[model_c_covars + ['hospital_expire_flag']].dropna()
X_c = model_c_data[model_c_covars].copy()
y_c = model_c_data['hospital_expire_flag']

X_c_c = sm.add_constant(X_c)
m_c = sm.Logit(y_c, X_c_c).fit(disp=0, maxiter=200)

print(f"  log_cr_ratio OR: {np.exp(m_c.params[1]):.3f} (95%CI {np.exp(m_c.params[1]-1.96*m_c.bse[1]):.3f}-{np.exp(m_c.params[1]+1.96*m_c.bse[1]):.3f}), p={m_c.pvalues[1]:.6f}")
print(f"  Model C AIC: {m_c.aic:.1f}")

# Save all logistic results
all_logistic = a_results + b_results
pd.DataFrame(all_logistic).to_csv(os.path.join(IN, 'logistic_results_revised.csv'), index=False)

# Extract KDIGO OR comparison for paper
print("\n  KDIGO Stage OR comparison:")
print(f"  {'Stage':<12} {'Model A OR':>12} {'Model B OR':>12} {'Delta':>10}")
for stage_num in [1, 2, 3]:
    col = f'kdigo_stage_{stage_num}'
    or_a = [r for r in a_results if r['variable'] == col][0]
    or_b = [r for r in b_results if r['variable'] == col][0]
    delta = (or_a['OR'] - or_b['OR']) / or_a['OR'] * 100
    print(f"  Stage {stage_num:<6} {or_a['OR']:>12.3f} {or_b['OR']:>12.3f} {delta:>9.1f}%")

# ============================================================
# 4. Kaplan-Meier + Log-rank (with Bonferroni)
# ============================================================
print("\n[4] Kaplan-Meier survival analysis...")

# Prepare 30-day survival
cohort['intime_dt'] = pd.to_datetime(cohort['intime'])
cohort['deathtime_dt'] = pd.to_datetime(cohort['deathtime'])
cohort['dod_dt'] = pd.to_datetime(cohort['dod'])

cohort['event_30d'] = 0
cohort['time_30d'] = 30.0

for idx in cohort.index:
    dod = cohort.at[idx, 'dod_dt']
    intime = cohort.at[idx, 'intime_dt']
    if pd.notna(dod) and pd.notna(intime):
        t = (dod - intime).total_seconds() / 86400
        if 0 < t <= 30:
            cohort.at[idx, 'time_30d'] = t
            cohort.at[idx, 'event_30d'] = 1

print(f"  30-day deaths: {cohort['event_30d'].sum():,}")

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test

km_data_by_stage = {}
for stage in [0, 1, 2, 3]:
    sub = cohort[cohort['kdigo_stage'] == stage]
    kmf = KaplanMeierFitter()
    kmf.fit(sub['time_30d'], sub['event_30d'], label=f'Stage {stage}')
    surv_table = kmf.survival_function_
    km_data_by_stage[str(stage)] = {
        'times': surv_table.index.tolist(),
        'survival': surv_table.values.flatten().tolist(),
        'n': len(sub),
        'events': int(sub['event_30d'].sum()),
        'n_at_risk': [len(sub[sub['time_30d'] >= t]) for t in [0, 7, 14, 21, 30]],
        'surv_at': {
            'day7': float(kmf.predict(7)),
            'day14': float(kmf.predict(14)),
            'day30': float(kmf.predict(30)),
        }
    }
    print(f"  Stage {stage}: 30d surv = {km_data_by_stage[str(stage)]['surv_at']['day30']:.4f}, "
          f"n={len(sub)}, events={int(sub['event_30d'].sum())}")

# Multivariate log-rank
lr_all = multivariate_logrank_test(cohort['time_30d'], cohort['kdigo_stage'], cohort['event_30d'])
print(f"\n  Overall log-rank: chi2={lr_all.test_statistic:.2f}, p={lr_all.p_value:.8f}")

# Pairwise with Bonferroni (6 comparisons)
print("\n  Pairwise log-rank (Bonferroni-corrected):")
stages = [0, 1, 2, 3]
pairwise_results = []
for i, s1 in enumerate(stages):
    for j, s2 in enumerate(stages[i+1:], i+1):
        g1 = cohort[cohort['kdigo_stage'] == s1]
        g2 = cohort[cohort['kdigo_stage'] == s2]
        lr = logrank_test(g1['time_30d'], g2['time_30d'], g1['event_30d'], g2['event_30d'])
        p_raw = lr.p_value
        p_bonf = min(p_raw * 6, 1.0)  # 6 pairwise comparisons
        sig = '***' if p_bonf < 0.001 else ('**' if p_bonf < 0.01 else ('*' if p_bonf < 0.05 else 'ns'))
        print(f"    Stage {s1} vs Stage {s2}: p_raw={p_raw:.6f}, p_bonf={p_bonf:.6f} {sig}")
        pairwise_results.append({
            'comparison': f'{s1}_vs_{s2}',
            'p_raw': round(float(p_raw), 8),
            'p_bonferroni': round(float(p_bonf), 8),
        })

with open(os.path.join(IN, 'km_revised.json'), 'w') as f:
    json.dump({'km_data': km_data_by_stage, 'pairwise': pairwise_results,
               'overall_chi2': round(float(lr_all.test_statistic), 2),
               'overall_p': round(float(lr_all.p_value), 8)}, f, indent=2)

# ============================================================
# 5. Cox Proportional Hazards (full covariate reporting)
# ============================================================
print("\n[5] Cox Proportional Hazards...")

from lifelines import CoxPHFitter

# Model A: KDIGO + base + comorbidities
cox_a_covars = ['kdigo_stage_1', 'kdigo_stage_2', 'kdigo_stage_3'] + base_covars + comorb_vars
cox_a_data = cohort[cox_a_covars + ['time_30d', 'event_30d']].dropna()
cph_a = CoxPHFitter()
cph_a.fit(cox_a_data, duration_col='time_30d', event_col='event_30d')

print("\n  --- Cox Model A: KDIGO + demographics + comorbidities ---")
print(f"  N = {cph_a._n_examples:,}, Events = {cox_a_data['event_30d'].sum():,}")
print(f"  C-index = {cph_a.concordance_index_:.3f}")
for col in cox_a_covars:
    hr = cph_a.summary.loc[col, 'exp(coef)']
    ci_low = cph_a.summary.loc[col, 'exp(coef) lower 95%']
    ci_high = cph_a.summary.loc[col, 'exp(coef) upper 95%']
    p = cph_a.summary.loc[col, 'p']
    if col.startswith('kdigo') or p < 0.001:
        print(f"    {col:<40} HR={hr:.3f} ({ci_low:.3f}-{ci_high:.3f}), p={p:.4f}")

# Model B: + severity
cox_b_covars = ['kdigo_stage_1', 'kdigo_stage_2', 'kdigo_stage_3'] + base_covars + comorb_vars + severity_vars
cox_b_data = cohort[cox_b_covars + ['time_30d', 'event_30d']].dropna()
cph_b = CoxPHFitter()
cph_b.fit(cox_b_data, duration_col='time_30d', event_col='event_30d')

print(f"\n  --- Cox Model B: KDIGO + severity scores ---")
print(f"  N = {cph_b._n_examples:,}, Events = {cox_b_data['event_30d'].sum():,}")
print(f"  C-index = {cph_b.concordance_index_:.3f}")
print(f"  ALL covariates:")
for col in cox_b_covars:
    hr = cph_b.summary.loc[col, 'exp(coef)']
    ci_low = cph_b.summary.loc[col, 'exp(coef) lower 95%']
    ci_high = cph_b.summary.loc[col, 'exp(coef) upper 95%']
    p = cph_b.summary.loc[col, 'p']
    sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else ''))
    marker = ' <-- ' if col.startswith('kdigo') or col in severity_vars else ''
    print(f"    {col:<42} HR={hr:.3f} ({ci_low:.3f}-{ci_high:.3f}), p={p:.4f} {sig}{marker}")

# Save full Cox tables
cph_a.summary.to_csv(os.path.join(IN, 'cox_a_full.csv'))
cph_b.summary.to_csv(os.path.join(IN, 'cox_b_full.csv'))

# KDIGO HR comparison
print(f"\n  KDIGO Stage HR comparison:")
for stage_num in [1, 2, 3]:
    col = f'kdigo_stage_{stage_num}'
    hr_a = np.exp(cph_a.params_[col])
    hr_b = np.exp(cph_b.params_[col])
    delta = (hr_a - hr_b) / hr_a * 100
    print(f"    Stage {stage_num}: Model A HR={hr_a:.3f}, Model B HR={hr_b:.3f}, attenuated {delta:.1f}%")

# ============================================================
# 6. Dose-Response with 95% CI
# ============================================================
print("\n[6] Dose-Response analysis...")

# Bin Cr ratio
cohort['cr_ratio_bin'] = pd.cut(cohort['cr_ratio'],
    bins=[0, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 100.0],
    labels=['<1.0', '1.0-1.5', '1.5-2.0', '2.0-3.0', '3.0-5.0', '5.0-10.0', '>=10.0'])

dr_data = []
for label in ['<1.0', '1.0-1.5', '1.5-2.0', '2.0-3.0', '3.0-5.0', '5.0-10.0', '>=10.0']:
    sub = cohort[cohort['cr_ratio_bin'] == label]
    n = len(sub)
    deaths = int(sub['hospital_expire_flag'].sum())
    rate = deaths / n * 100 if n > 0 else 0
    # Binomial 95% CI
    from scipy.stats import binom
    ci_low = binom.interval(0.95, n, rate/100)[0] / n * 100 if n > 0 else 0
    ci_high = binom.interval(0.95, n, rate/100)[1] / n * 100 if n > 0 else 0
    dr_data.append({
        'cr_ratio_bin': label, 'n': n, 'deaths': deaths,
        'mortality_pct': round(rate, 1),
        'CI_low': round(float(ci_low), 1),
        'CI_high': round(float(ci_high), 1),
    })
    print(f"  {label}: n={n:,}, mortality={rate:.1f}% ({ci_low:.1f}-{ci_high:.1f}%)")

pd.DataFrame(dr_data).to_csv(os.path.join(IN, 'dose_response_revised.csv'), index=False)

# ============================================================
# 7. AKI Paradox: RRT-stratified Stage 3
# ============================================================
print("\n[7] AKI Paradox — RRT-stratified Stage 3 analysis...")

for ckd_status, ckd_label in [(0, 'CKD-'), (1, 'CKD+')]:
    sub = cohort[cohort['ckd_icd'] == ckd_status]
    print(f"\n  {ckd_label} (n={len(sub):,}):")
    for stage in [0, 1, 2, 3]:
        ss = sub[sub['kdigo_stage'] == stage]
        mort = ss['hospital_expire_flag'].mean() * 100 if len(ss) > 0 else 0
        print(f"    Stage {stage}: n={len(ss):,}, mortality={mort:.1f}%")

# RRT-stratified Stage 3
print("\n  Stage 3 stratification by entry mechanism:")
for ckd_status, ckd_label in [(0, 'CKD-'), (1, 'CKD+')]:
    print(f"\n  {ckd_label} Stage 3:")
    s3 = cohort[(cohort['kdigo_stage'] == 3) & (cohort['ckd_icd'] == ckd_status)]
    s3_rrt = s3[s3['rrt'] == 1]
    s3_cr_only = s3[s3['rrt'] == 0]
    
    print(f"    Total Stage 3: {len(s3):,}, mortality={s3['hospital_expire_flag'].mean()*100:.1f}%")
    print(f"    RRT received: {len(s3_rrt):,} ({len(s3_rrt)/max(len(s3),1)*100:.1f}%), "
          f"mortality={s3_rrt['hospital_expire_flag'].mean()*100:.1f}%")
    print(f"    Cr criteria only (no RRT): {len(s3_cr_only):,} ({len(s3_cr_only)/max(len(s3),1)*100:.1f}%), "
          f"mortality={s3_cr_only['hospital_expire_flag'].mean()*100:.1f}%")
    
    # By substrata
    for substratum, sublabel in [('stage3_by_rrt_only', 'RRT only'), 
                                  ('stage3_by_cr_only', 'Cr only'), 
                                  ('stage3_by_both', 'Both')]:
        ss = s3[s3[substratum] == 1]
        if len(ss) > 0:
            print(f"    {sublabel}: n={len(ss):,}, mortality={ss['hospital_expire_flag'].mean()*100:.1f}%")

# Key AKI paradox test: Cr-only Stage 3, CKD+ vs CKD-
print("\n  AKI Paradox — Cr-only Stage 3 (most stringent test):")
cr_only_s3 = cohort[(cohort['kdigo_stage'] == 3) & (cohort['rrt'] == 0)]
for ckd_status, ckd_label in [(0, 'CKD-'), (1, 'CKD+')]:
    ss = cr_only_s3[cr_only_s3['ckd_icd'] == ckd_status]
    print(f"    {ckd_label}: n={len(ss):,}, mortality={ss['hospital_expire_flag'].mean()*100:.1f}%")

# ============================================================
# 8. MDRD Sensitivity Analysis
# ============================================================
print("\n[8] MDRD baseline sensitivity analysis...")

mdrd_valid = cohort[cohort['mdrd_baseline_cr'].notna()]
print(f"  Patients with MDRD baseline: {len(mdrd_valid):,}")
print(f"  KDIGO stage distribution (MDRD baseline):")
for s in [0, 1, 2, 3]:
    n = (mdrd_valid['kdigo_stage_mdrd'] == s).sum()
    mort = mdrd_valid[mdrd_valid['kdigo_stage_mdrd'] == s]['hospital_expire_flag'].mean() * 100 if n > 0 else 0
    print(f"    Stage {s}: {n:,} ({n/len(mdrd_valid)*100:.1f}%), mortality={mort:.1f}%")

# Concordance between original and MDRD staging
agreement = (cohort['kdigo_stage'] == cohort['kdigo_stage_mdrd']).sum()
kappa_denom = len(mdrd_valid)
print(f"  Agreement: {agreement}/{kappa_denom} ({agreement/kappa_denom*100:.1f}%)")

# Logistic with MDRD
mdrd_data = mdrd_valid.copy()
for s in [1, 2, 3]:
    mdrd_data[f'kdigo_mdrd_{s}'] = (mdrd_data['kdigo_stage_mdrd'] == s).astype(int)

mdrd_covars = ['kdigo_mdrd_1', 'kdigo_mdrd_2', 'kdigo_mdrd_3'] + base_covars + comorb_vars + severity_vars
mdrd_model_data = mdrd_data[mdrd_covars + ['hospital_expire_flag']].dropna()
X_mdrd = mdrd_model_data[mdrd_covars]
y_mdrd = mdrd_model_data['hospital_expire_flag']
X_mdrd_c = sm.add_constant(X_mdrd)
m_mdrd = sm.Logit(y_mdrd, X_mdrd_c).fit(disp=0, maxiter=200)

print(f"\n  MDRD-based KDIGO OR (Model B equivalent):")
for s in [1, 2, 3]:
    col = f'kdigo_mdrd_{s}'
    idx = mdrd_covars.index(col)
    coef = m_mdrd.params[idx + 1]
    se = m_mdrd.bse[idx + 1]
    or_val = np.exp(coef)
    ci_low = np.exp(coef - 1.96 * se)
    ci_high = np.exp(coef + 1.96 * se)
    p_val = m_mdrd.pvalues[idx + 1]
    print(f"    Stage {s}: OR={or_val:.3f} ({ci_low:.3f}-{ci_high:.3f}), p={p_val:.4f}")

# ============================================================
# 9. Lab-based CKD sensitivity
# ============================================================
print("\n[9] Lab-based CKD sensitivity analysis...")

# Re-run AKI paradox with lab-based CKD
print("  AKI Paradox (Lab-based CKD):")
for ckd_val, ckd_label in [(0, 'CKD- (lab)'), (1, 'CKD+ (lab)')]:
    sub = cohort[cohort['ckd_lab'] == ckd_val]
    print(f"\n  {ckd_label} (n={len(sub):,}):")
    for stage in [0, 1, 2, 3]:
        ss = sub[sub['kdigo_stage'] == stage]
        mort = ss['hospital_expire_flag'].mean() * 100 if len(ss) > 0 else 0
        print(f"    Stage {stage}: n={len(ss):,}, mortality={mort:.1f}%")

# Logistic with lab CKD instead of ICD CKD
lab_ckd_covars = ['kdigo_stage_1', 'kdigo_stage_2', 'kdigo_stage_3'] + base_covars + comorb_vars + ['sofa_imputed', 'mech_vent', 'vasopressor', 'ckd_lab']
lab_ckd_data = cohort[lab_ckd_covars + ['hospital_expire_flag']].dropna()
X_lab = lab_ckd_data[lab_ckd_covars]
y_lab = lab_ckd_data['hospital_expire_flag']
X_lab_c = sm.add_constant(X_lab)
m_lab = sm.Logit(y_lab, X_lab_c).fit(disp=0, maxiter=200)

print(f"\n  Lab-based CKD vs ICD-CKD in full model:")
for ckd_var in ['ckd_icd', 'ckd_lab']:
    if ckd_var in model_b_covars or ckd_var in lab_ckd_covars:
        continue
# Actually, let's just print comparison of CKD OR
# ICD-CKD from Model B
ckd_icd_or = np.exp(m_b.params[list(model_b_covars).index('ckd_icd') + 1])
ckd_icd_p = m_b.pvalues[list(model_b_covars).index('ckd_icd') + 1]
# Lab-CKD
ckd_lab_or = np.exp(m_lab.params[list(lab_ckd_covars).index('ckd_lab') + 1])
ckd_lab_p = m_lab.pvalues[list(lab_ckd_covars).index('ckd_lab') + 1]
print(f"  ICD-CKD: OR={ckd_icd_or:.3f}, p={ckd_icd_p:.4f}")
print(f"  Lab-CKD: OR={ckd_lab_or:.3f}, p={ckd_lab_p:.4f}")

# ============================================================
# 10. Final Summary Output
# ============================================================
print("\n[10] Generating final summary...")

# Extract KDIGO ORs from Model B for paper
kdigo_or = {}
for s in [1, 2, 3]:
    col = f'kdigo_stage_{s}'
    idx = model_b_covars.index(col)
    coef = m_b.params[idx + 1]
    se = m_b.bse[idx + 1]
    kdigo_or[s] = {
        'OR': round(float(np.exp(coef)), 3),
        'CI_low': round(float(np.exp(coef - 1.96 * se)), 3),
        'CI_high': round(float(np.exp(coef + 1.96 * se)), 3),
        'p': round(float(m_b.pvalues[idx + 1]), 6),
    }

# KDIGO HRs from Cox B
kdigo_hr = {}
for s in [1, 2, 3]:
    col = f'kdigo_stage_{s}'
    hr = cph_b.summary.loc[col, 'exp(coef)']
    ci_low = cph_b.summary.loc[col, 'exp(coef) lower 95%']
    ci_high = cph_b.summary.loc[col, 'exp(coef) upper 95%']
    p = cph_b.summary.loc[col, 'p']
    kdigo_hr[s] = {'HR': round(float(hr), 3), 'CI_low': round(float(ci_low), 3),
                   'CI_high': round(float(ci_high), 3), 'p': round(float(p), 6)}

summary = {
    'n_total': int(len(cohort)),
    'kdigo_distribution': {
        '0': int((cohort['kdigo_stage'] == 0).sum()),
        '1': int((cohort['kdigo_stage'] == 1).sum()),
        '2': int((cohort['kdigo_stage'] == 2).sum()),
        '3': int((cohort['kdigo_stage'] == 3).sum()),
    },
    'aki_prevalence_pct': round((cohort['kdigo_stage'] >= 1).sum() / len(cohort) * 100, 1),
    'mortality_by_stage': {
        '0': round(float(cohort[cohort['kdigo_stage'] == 0]['hospital_expire_flag'].mean() * 100), 2),
        '1': round(float(cohort[cohort['kdigo_stage'] == 1]['hospital_expire_flag'].mean() * 100), 2),
        '2': round(float(cohort[cohort['kdigo_stage'] == 2]['hospital_expire_flag'].mean() * 100), 2),
        '3': round(float(cohort[cohort['kdigo_stage'] == 3]['hospital_expire_flag'].mean() * 100), 2),
    },
    'kdigo_or_model_b': kdigo_or,
    'kdigo_hr_model_b': kdigo_hr,
    'model_a_aic': round(float(m_a.aic), 1),
    'model_b_aic': round(float(m_b.aic), 1),
    'cox_b_cindex': round(float(cph_b.concordance_index_), 3),
    'logrank_chi2': round(float(lr_all.test_statistic), 2),
    'logrank_p': round(float(lr_all.p_value), 8),
    'stage3_rrt_stratification': {
        'rrt_only': int((cohort['stage3_by_rrt_only'] == 1).sum()),
        'cr_only': int((cohort['stage3_by_cr_only'] == 1).sum()),
        'both': int((cohort['stage3_by_both'] == 1).sum()),
    },
    'ckd_icd_vs_lab': {
        'ckd_icd_n': int(cohort['ckd_icd'].sum()),
        'ckd_lab_n': int(cohort['ckd_lab'].sum()),
        'both': int(((cohort['ckd_icd'] == 1) & (cohort['ckd_lab'] == 1)).sum()),
        'icd_sensitivity_vs_lab': round(float(ckd_both_n) / max(int(cohort['ckd_lab'].sum()), 1) * 100, 1) if 'ckd_both_n' in dir() else 0,
    },
}

# Fix CKD sensitivity
ckd_both_n = int(((cohort['ckd_icd'] == 1) & (cohort['ckd_lab'] == 1)).sum())
summary['ckd_icd_vs_lab']['both'] = ckd_both_n
summary['ckd_icd_vs_lab']['icd_sensitivity_vs_lab'] = round(ckd_both_n / max(int(cohort['ckd_lab'].sum()), 1) * 100, 1)

with open(os.path.join(IN, 'summary_revised.json'), 'w') as f:
    json.dump(summary, f, indent=2)

print("\n" + "=" * 80)
print("STATISTICAL ANALYSIS COMPLETE")
print("=" * 80)
print(f"\nKey results (Model B = KDIGO + severity):")
for s in [1, 2, 3]:
    or_d = kdigo_or[s]
    hr_d = kdigo_hr[s]
    print(f"  Stage {s}: OR={or_d['OR']} ({or_d['CI_low']}-{or_d['CI_high']}), "
          f"HR={hr_d['HR']} ({hr_d['CI_low']}-{hr_d['CI_high']})")
print(f"\n  Model AIC: {m_a.aic:.0f} -> {m_b.aic:.0f} (delta={m_a.aic-m_b.aic:.0f})")
print(f"  Cox C-index: {cph_b.concordance_index_:.3f}")
print(f"\nEnd: {datetime.now()}")
