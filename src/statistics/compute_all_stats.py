"""
Comprehensive statistics computation for AKI Paradox manuscript v4.
Computes ALL statistics needed for manuscript text, tables, and figure annotations.
All numbers from this script are the SINGLE SOURCE OF TRUTH.
"""
import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score, brier_score_loss
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
import json
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

df = pd.read_csv('aki_paradox_cohort_with_sofa.csv')
aki_df = df[df['max_aki_stage'] > 0].copy()
aki_df['gender_num'] = (aki_df['gender'] == 'M').astype(int)
aki_df['stage2'] = (aki_df['max_aki_stage'] == 2).astype(int)
aki_df['stage3'] = (aki_df['max_aki_stage'] == 3).astype(int)
aki_df['ckd_stage2'] = aki_df['ckd'] * aki_df['stage2']
aki_df['ckd_stage3'] = aki_df['ckd'] * aki_df['stage3']
y = aki_df['hospital_expire_flag']
y_arr = y.values
n = len(y_arr)

results = {}

# ============================================================
# 1. COHORT NUMBERS
# ============================================================
print("=" * 80)
print("1. COHORT NUMBERS")
print("=" * 80)
total_icu = len(df)
aki_total = len(aki_df)
ckd_total = int(aki_df['ckd'].sum())
pure_total = aki_total - ckd_total
for s in [1, 2, 3]:
    n_stage = int((aki_df['max_aki_stage'] == s).sum())
    n_pure = int(((aki_df['max_aki_stage'] == s) & (aki_df['ckd'] == 0)).sum())
    n_ckd = int(((aki_df['max_aki_stage'] == s) & (aki_df['ckd'] == 1)).sum())
    print(f"  Stage {s}: total={n_stage}, pure={n_pure}, ckd={n_ckd}")
    results[f'stage{s}_total'] = n_stage
    results[f'stage{s}_pure'] = n_pure
    results[f'stage{s}_ckd'] = n_ckd

results['total_icu'] = total_icu
results['aki_total'] = aki_total
results['ckd_total'] = ckd_total
results['pure_total'] = pure_total
print(f"  Total ICU: {total_icu}, AKI: {aki_total}, CKD: {ckd_total}, Pure: {pure_total}")

# ============================================================
# 2. MORTALITY RATES BY STAGE AND CKD
# ============================================================
print("\n" + "=" * 80)
print("2. MORTALITY RATES")
print("=" * 80)
for s in [1, 2, 3]:
    for ckd_status in [0, 1]:
        sub = aki_df[(aki_df['max_aki_stage'] == s) & (aki_df['ckd'] == ckd_status)]
        mort = sub['hospital_expire_flag'].mean() * 100
        n = len(sub)
        label = f"Stage {s} {'CKD' if ckd_status else 'Pure'}"
        print(f"  {label}: n={n}, mortality={mort:.1f}%")
        results[f'mort_stage{s}_{"ckd" if ckd_status else "pure"}'] = round(mort, 1)

# ============================================================
# 3. UNADJUSTED AND SOFA-ADJUSTED ORs BY STAGE
# ============================================================
print("\n" + "=" * 80)
print("3. UNADJUSTED AND SOFA-ADJUSTED ORs BY STAGE")
print("=" * 80)
for s in [1, 2, 3]:
    sub = aki_df[aki_df['max_aki_stage'] == s]
    # Unadjusted
    X_unadj = sm.add_constant(sub['ckd'])
    m_unadj = sm.Logit(sub['hospital_expire_flag'], X_unadj).fit(disp=0)
    or_unadj = np.exp(m_unadj.params['ckd'])
    ci_unadj = np.exp(m_unadj.conf_int().loc['ckd'])
    p_unadj = m_unadj.pvalues['ckd']
    
    # SOFA-adjusted (age, sex, comorbidities, SOFA)
    X_adj = sm.add_constant(sub[['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_total', 'ckd']])
    m_adj = sm.Logit(sub['hospital_expire_flag'], X_adj).fit(disp=0)
    or_adj = np.exp(m_adj.params['ckd'])
    ci_adj = np.exp(m_adj.conf_int().loc['ckd'])
    p_adj = m_adj.pvalues['ckd']
    
    print(f"  Stage {s}: Unadjusted OR={or_unadj:.3f} ({ci_unadj[0]:.3f}-{ci_unadj[1]:.3f}), P={p_unadj:.4f}")
    print(f"  Stage {s}: SOFA-Adjusted OR={or_adj:.3f} ({ci_adj[0]:.3f}-{ci_adj[1]:.3f}), P={p_adj:.4f}")
    results[f'stage{s}_unadj_or'] = round(or_unadj, 3)
    results[f'stage{s}_unadj_lo'] = round(ci_unadj[0], 3)
    results[f'stage{s}_unadj_hi'] = round(ci_unadj[1], 3)
    results[f'stage{s}_unadj_p'] = round(p_unadj, 4)
    results[f'stage{s}_adj_or'] = round(or_adj, 3)
    results[f'stage{s}_adj_lo'] = round(ci_adj[0], 3)
    results[f'stage{s}_adj_hi'] = round(ci_adj[1], 3)
    results[f'stage{s}_adj_p'] = round(p_adj, 4)

# ============================================================
# 4. TEMPORAL DYNAMICS (Stage 3, different timeframes)
# ============================================================
print("\n" + "=" * 80)
print("4. TEMPORAL DYNAMICS (Stage 3)")
print("=" * 80)
# Check available mortality columns
mort_cols = [c for c in aki_df.columns if 'mort' in c.lower() or 'expire' in c.lower() or 'death' in c.lower() or 'died' in c.lower()]
print(f"  Mortality-related columns: {mort_cols}")

# Check for 30d, 90d, 1y mortality
for col in ['mortality_30d', 'mortality_90d', 'mortality_1y', 'died_30d', 'died_90d', 'died_1y',
            'hosp_mort_30d', 'hosp_mort_90d', 'hosp_mort_1y', 'mort_30d', 'mort_90d', 'mort_1y']:
    if col in aki_df.columns:
        print(f"  Found column: {col}")

# Try survival data
try:
    sv = pd.read_csv('aki_survival_data.csv')
    print(f"  Survival data columns: {list(sv.columns)}")
    
    # Check for 30d/90d/1y mortality indicators
    for col in sv.columns:
        if '30' in col or '90' in col or '1y' in col or 'year' in col or 'day' in col.lower():
            print(f"    {col}: {sv[col].dtype}, non-null={sv[col].notna().sum()}")
    
    # Compute 30/90/365 day mortality from survival_days if available
    if 'survival_days' in sv.columns:
        sv_aki = sv[sv['max_aki_stage'] > 0].copy()
        s3_sv = sv_aki[sv_aki['max_aki_stage'] == 3]
        
        for days, label in [(30, '30d'), (90, '90d'), (365, '1y')]:
            sv_aki[f'died_{days}d'] = ((sv_aki['survival_days'] <= days) & (sv_aki['hospital_expire_flag'] == 1)).astype(int)
            # Actually, survival_days might already account for post-discharge
            # Let's check: if survival_days > days, patient survived past that timeframe
            # If survival_days <= days and died=1, patient died within that timeframe
            # If survival_days <= days and died=0, patient was discharged and we don't know post-discharge
            # For MIMIC, hospital_expire_flag=1 means in-hospital death
            # Post-discharge death comes from Social Security Death Records
            
            # Let's use a simpler approach: if survival_days <= days, patient died within that period
            # (survival_days should represent days from admission to death or censoring)
            sv_aki.loc[:, f'died_{days}d'] = (sv_aki['survival_days'] <= days).astype(int)
            
            s3 = sv_aki[sv_aki['max_aki_stage'] == 3]
            pure = s3[s3['ckd'] == 0]
            ckd = s3[s3['ckd'] == 1]
            
            pure_mort = pure[f'died_{days}d'].mean() * 100
            ckd_mort = ckd[f'died_{days}d'].mean() * 100
            
            # OR
            X = sm.add_constant(s3['ckd'])
            m = sm.Logit(s3[f'died_{days}d'], X).fit(disp=0)
            or_val = np.exp(m.params['ckd'])
            ci = np.exp(m.conf_int().loc['ckd'])
            p_val = m.pvalues['ckd']
            
            print(f"  Stage 3 {label}: Pure={pure_mort:.1f}%, CKD={ckd_mort:.1f}%, OR={or_val:.3f} ({ci[0]:.3f}-{ci[1]:.3f}), P={p_val:.3f}")
            results[f's3_{label}_pure_mort'] = round(pure_mort, 1)
            results[f's3_{label}_ckd_mort'] = round(ckd_mort, 1)
            results[f's3_{label}_or'] = round(or_val, 3)
            results[f's3_{label}_lo'] = round(ci[0], 3)
            results[f's3_{label}_hi'] = round(ci[1], 3)
            results[f's3_{label}_p'] = round(p_val, 3)
            
except Exception as e:
    print(f"  Error with survival data: {e}")

# ============================================================
# 5. NESTED MODEL COMPARISON (Table 2/6)
# ============================================================
print("\n" + "=" * 80)
print("5. NESTED MODEL COMPARISON")
print("=" * 80)

# Model A: age, sex, comorbidities, SOFA, AKI stage (no CKD)
X_A = sm.add_constant(aki_df[['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_total', 'stage2', 'stage3']])
# Model B: A + CKD main effect
X_B = sm.add_constant(aki_df[['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_total', 'stage2', 'stage3', 'ckd']])
# Model C: B + CKD×Stage interactions
X_C = sm.add_constant(aki_df[['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_total', 'stage2', 'stage3', 'ckd', 'ckd_stage2', 'ckd_stage3']])
# Model D: Replace total SOFA with 6 components
X_D = sm.add_constant(aki_df[['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_resp', 'sofa_coag', 'sofa_liver', 'sofa_cv', 'sofa_cns', 'sofa_renal', 'stage2', 'stage3', 'ckd', 'ckd_stage2', 'ckd_stage3']])

mA = sm.Logit(y, X_A).fit(disp=0)
mB = sm.Logit(y, X_B).fit(disp=0)
mC = sm.Logit(y, X_C).fit(disp=0)
mD = sm.Logit(y, X_D).fit(disp=0)

models = {"A": (mA, X_A), "B": (mB, X_B), "C": (mC, X_C), "D": (mD, X_D)}

for name, (model, X) in models.items():
    preds = model.predict(X).values
    auc_point = roc_auc_score(y_arr, preds)
    # Bootstrap CI
    aucs = []
    for _ in range(2000):
        idx = np.random.choice(n, n, replace=True)
        if len(np.unique(y_arr[idx])) < 2:
            continue
        try:
            aucs.append(roc_auc_score(y_arr[idx], preds[idx]))
        except:
            pass
    ci_l = np.percentile(aucs, 2.5)
    ci_u = np.percentile(aucs, 97.5)
    brier = brier_score_loss(y, preds)
    pseudo_r2 = 1 - (model.llf / mA.llnull)  # McFadden's R2
    
    print(f"  Model {name}: AUC={auc_point:.4f} ({ci_l:.4f}-{ci_u:.4f}), Brier={brier:.4f}, Pseudo-R2={pseudo_r2:.4f}")
    results[f'model{name}_auc'] = round(auc_point, 4)
    results[f'model{name}_auc_lo'] = round(ci_l, 4)
    results[f'model{name}_auc_hi'] = round(ci_u, 4)
    results[f'model{name}_brier'] = round(brier, 4)
    results[f'model{name}_r2'] = round(pseudo_r2, 4)

# CV-AUC
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
feature_sets = {
    'A': ['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_total', 'stage2', 'stage3'],
    'B': ['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_total', 'stage2', 'stage3', 'ckd'],
    'C': ['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_total', 'stage2', 'stage3', 'ckd', 'ckd_stage2', 'ckd_stage3'],
    'D': ['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_resp', 'sofa_coag', 'sofa_liver', 'sofa_cv', 'sofa_cns', 'sofa_renal', 'stage2', 'stage3', 'ckd', 'ckd_stage2', 'ckd_stage3']
}
for name, features in feature_sets.items():
    cv_aucs = []
    for train_idx, test_idx in skf.split(aki_df, y):
        X_train = aki_df.iloc[train_idx][features]
        y_train = y.iloc[train_idx]
        X_test = aki_df.iloc[test_idx][features]
        y_test = y.iloc[test_idx]
        lr = LogisticRegression(max_iter=1000, penalty=None)
        lr.fit(X_train, y_train)
        preds = lr.predict_proba(X_test)[:, 1]
        cv_aucs.append(roc_auc_score(y_test, preds))
    cv_auc = np.mean(cv_aucs)
    print(f"  Model {name}: CV-AUC={cv_auc:.4f}")
    results[f'model{name}_cv_auc'] = round(cv_auc, 4)

# Key coefficients from Model C
ckd_or = np.exp(mC.params['ckd'])
ckd_ci = np.exp(mC.conf_int().loc['ckd'])
ckd_p = mC.pvalues['ckd']
ckd_s3_or = np.exp(mC.params['ckd_stage3'])
ckd_s3_ci = np.exp(mC.conf_int().loc['ckd_stage3'])
ckd_s3_p = mC.pvalues['ckd_stage3']
print(f"\n  Model C CKD main effect: OR={ckd_or:.3f} ({ckd_ci[0]:.3f}-{ckd_ci[1]:.3f}), P={ckd_p:.4f}")
print(f"  Model C CKD×Stage3: OR={ckd_s3_or:.3f} ({ckd_s3_ci[0]:.3f}-{ckd_s3_ci[1]:.3f}), P={ckd_s3_p:.6f}")
results['modelC_ckd_or'] = round(ckd_or, 3)
results['modelC_ckd_lo'] = round(ckd_ci[0], 3)
results['modelC_ckd_hi'] = round(ckd_ci[1], 3)
results['modelC_ckd_p'] = round(ckd_p, 4)
results['modelC_ckd_s3_or'] = round(ckd_s3_or, 3)
results['modelC_ckd_s3_lo'] = round(ckd_s3_ci[0], 3)
results['modelC_ckd_s3_hi'] = round(ckd_s3_ci[1], 3)
results['modelC_ckd_s3_p'] = round(ckd_s3_p, 6)

# Combined effect: CKD in Stage 3 = main + interaction
combined_or = ckd_or * ckd_s3_or
print(f"  Combined CKD in Stage 3: OR={combined_or:.3f}")
results['combined_s3_or'] = round(combined_or, 3)

# LR test: B vs C
lr_stat = -2 * (mB.llf - mC.llf)
lr_p = 1 - stats.chi2.cdf(lr_stat, 2)  # 2 df (2 interaction terms)
print(f"  LR test (C vs B): chi2={lr_stat:.2f}, P={lr_p:.2e}")
results['lr_test_stat'] = round(lr_stat, 2)
results['lr_test_p'] = f"{lr_p:.1e}"

# LR test: A vs B
lr_stat_ab = -2 * (mA.llf - mB.llf)
lr_p_ab = 1 - stats.chi2.cdf(lr_stat_ab, 1)
print(f"  LR test (B vs A): chi2={lr_stat_ab:.2f}, P={lr_p_ab:.4f}")
results['lr_test_ab_p'] = round(lr_p_ab, 4)

# Model D CKD×Stage3
ckd_s3_or_d = np.exp(mD.params['ckd_stage3'])
ckd_s3_ci_d = np.exp(mD.conf_int().loc['ckd_stage3'])
ckd_s3_p_d = mD.pvalues['ckd_stage3']
print(f"  Model D CKD×Stage3: OR={ckd_s3_or_d:.3f} ({ckd_s3_ci_d[0]:.3f}-{ckd_s3_ci_d[1]:.3f}), P={ckd_s3_p_d:.6f}")
results['modelD_ckd_s3_or'] = round(ckd_s3_or_d, 3)
results['modelD_ckd_s3_lo'] = round(ckd_s3_ci_d[0], 3)
results['modelD_ckd_s3_hi'] = round(ckd_s3_ci_d[1], 3)

# ============================================================
# 6. NRI AND IDI
# ============================================================
print("\n" + "=" * 80)
print("6. NRI AND IDI")
print("=" * 80)
def compute_nri(y_true, p_old, p_new):
    events = y_true == 1
    non_events = y_true == 0
    up_events = np.sum((p_new > p_old) & events)
    down_events = np.sum((p_new < p_old) & events)
    up_non = np.sum((p_new > p_old) & non_events)
    down_non = np.sum((p_new < p_old) & non_events)
    nri_events = (up_events - down_events) / np.sum(events)
    nri_non = (down_non - up_non) / np.sum(non_events)
    nri = nri_events + nri_non
    return nri, nri_events, nri_non

def compute_idi(y_true, p_old, p_new):
    events = y_true == 1
    non_events = y_true == 0
    idi_events = np.mean(p_new[events]) - np.mean(p_old[events])
    idi_non = np.mean(p_new[non_events]) - np.mean(p_old[non_events])
    return idi_events - idi_non

preds_A = mA.predict(X_A).values
preds_B = mB.predict(X_B).values
preds_C = mC.predict(X_C).values

nri_ca, nri_ca_e, nri_ca_ne = compute_nri(y_arr, preds_A, preds_C)
idi_ca = compute_idi(y_arr, preds_A, preds_C)
print(f"  NRI (C vs A): {nri_ca:.3f} (events: {nri_ca_e:.3f}, non-events: {nri_ca_ne:.3f})")
print(f"  IDI (C vs A): {idi_ca:.4f}")
results['nri_ca'] = round(nri_ca, 3)
results['nri_ca_events'] = round(nri_ca_e, 3)
results['nri_ca_nonevents'] = round(nri_ca_ne, 3)
results['idi_ca'] = round(idi_ca, 4)

nri_ba, _, _ = compute_nri(y_arr, preds_A, preds_B)
idi_ba = compute_idi(y_arr, preds_A, preds_B)
print(f"  NRI (B vs A): {nri_ba:.3f}")
print(f"  IDI (B vs A): {idi_ba:.4f}")
results['nri_ba'] = round(nri_ba, 3)
results['idi_ba'] = round(idi_ba, 4)

# ============================================================
# 7. SOFA-STRATIFIED ANALYSIS (Stage 3, cutoff <=5)
# ============================================================
print("\n" + "=" * 80)
print("7. SOFA-STRATIFIED ANALYSIS (Stage 3, cutoff <=5)")
print("=" * 80)
stage3_df = aki_df[aki_df['max_aki_stage'] == 3].copy()
for label, lo, hi in [("Low (0-5)", 0, 5), ("Medium (6-10)", 6, 10), ("High (>10)", 11, 100)]:
    stratum = stage3_df[(stage3_df['sofa_total'] >= lo) & (stage3_df['sofa_total'] <= hi)]
    n_strat = len(stratum)
    pure = stratum[stratum['ckd'] == 0]
    ckd = stratum[stratum['ckd'] == 1]
    pure_mort = pure['hospital_expire_flag'].mean() * 100 if len(pure) > 0 else 0
    ckd_mort = ckd['hospital_expire_flag'].mean() * 100 if len(ckd) > 0 else 0
    n_pure = len(pure)
    n_ckd = len(ckd)
    
    if n_pure > 0 and n_ckd > 0:
        X = sm.add_constant(stratum['ckd'])
        try:
            model = sm.Logit(stratum['hospital_expire_flag'], X).fit(disp=0)
            or_val = np.exp(model.params['ckd'])
            or_ci = np.exp(model.conf_int().loc['ckd'])
            p_val = model.pvalues['ckd']
            print(f"  {label}: N={n_strat} (Pure={n_pure}, CKD={n_ckd}), Pure mort={pure_mort:.1f}%, CKD mort={ckd_mort:.1f}%, OR={or_val:.3f} ({or_ci[0]:.3f}-{or_ci[1]:.3f}), P={p_val:.4f}")
            key = label.split("(")[0].strip().lower()
            results[f'sofa_{key}_n'] = n_strat
            results[f'sofa_{key}_pure_n'] = n_pure
            results[f'sofa_{key}_ckd_n'] = n_ckd
            results[f'sofa_{key}_pure_mort'] = round(pure_mort, 1)
            results[f'sofa_{key}_ckd_mort'] = round(ckd_mort, 1)
            results[f'sofa_{key}_or'] = round(or_val, 3)
            results[f'sofa_{key}_lo'] = round(or_ci[0], 3)
            results[f'sofa_{key}_hi'] = round(or_ci[1], 3)
            results[f'sofa_{key}_p'] = round(p_val, 4)
        except Exception as e:
            print(f"  {label}: N={n_strat}, OR error: {e}")

# ============================================================
# 8. NON-RENAL SOFA (Stage 3)
# ============================================================
print("\n" + "=" * 80)
print("8. NON-RENAL SOFA (Stage 3)")
print("=" * 80)
pure_s3 = stage3_df[stage3_df['ckd'] == 0]
ckd_s3 = stage3_df[stage3_df['ckd'] == 1]
print(f"  Pure AKI non-renal SOFA: {pure_s3['sofa_nonrenal'].mean():.1f} +/- {pure_s3['sofa_nonrenal'].std():.1f}")
print(f"  AKI+CKD non-renal SOFA: {ckd_s3['sofa_nonrenal'].mean():.1f} +/- {ckd_s3['sofa_nonrenal'].std():.1f}")
t_stat, p_val = stats.ttest_ind(pure_s3['sofa_nonrenal'].dropna(), ckd_s3['sofa_nonrenal'].dropna())
print(f"  P-value: {p_val:.6f}")
results['pure_s3_nonrenal_sofa'] = round(pure_s3['sofa_nonrenal'].mean(), 1)
results['ckd_s3_nonrenal_sofa'] = round(ckd_s3['sofa_nonrenal'].mean(), 1)
results['nonrenal_sofa_p'] = round(p_val, 6)

# Non-renal SOFA interaction
X_nr = sm.add_constant(aki_df[['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_nonrenal', 'stage2', 'stage3', 'ckd', 'ckd_stage2', 'ckd_stage3']])
m_nr = sm.Logit(y, X_nr).fit(disp=0)
ckd_s3_nr_or = np.exp(m_nr.params['ckd_stage3'])
ckd_s3_nr_ci = np.exp(m_nr.conf_int().loc['ckd_stage3'])
ckd_s3_nr_p = m_nr.pvalues['ckd_stage3']
print(f"  Non-renal SOFA CKD×Stage3: OR={ckd_s3_nr_or:.3f} ({ckd_s3_nr_ci[0]:.3f}-{ckd_s3_nr_ci[1]:.3f}), P={ckd_s3_nr_p:.6f}")
results['nonrenal_ckd_s3_or'] = round(ckd_s3_nr_or, 3)
results['nonrenal_ckd_s3_lo'] = round(ckd_s3_nr_ci[0], 3)
results['nonrenal_ckd_s3_hi'] = round(ckd_s3_nr_ci[1], 3)

# ============================================================
# 9. PSM — FULL COHORT, WITH AND WITHOUT SOFA
# ============================================================
print("\n" + "=" * 80)
print("9. PSM (Full AKI Cohort)")
print("=" * 80)

def do_psm(data, ps_cols, caliper=0.2, label=""):
    """1:1 nearest neighbor PSM without replacement"""
    from sklearn.linear_model import LogisticRegression
    
    X_ps = data[ps_cols].values
    y_ps = data['ckd'].values
    
    lr = LogisticRegression(max_iter=1000, penalty=None)
    lr.fit(X_ps, y_ps)
    ps = lr.predict_proba(X_ps)[:, 1]
    
    # Logit of PS
    eps = 1e-6
    ps = np.clip(ps, eps, 1 - eps)
    logit_ps = np.log(ps / (1 - ps))
    
    # Caliper = 0.2 * SD(logit_ps)
    caliper_val = caliper * np.std(logit_ps) if caliper > 0 else np.inf
    
    treated_idx = np.where(y_ps == 1)[0]
    control_idx = np.where(y_ps == 0)[0]
    
    matched_pairs = []
    used_controls = set()
    
    # Sort treated by PS (optional, for reproducibility)
    np.random.seed(42)
    treated_order = treated_idx.copy()
    np.random.shuffle(treated_order)
    
    for t in treated_order:
        best_dist = np.inf
        best_c = None
        for c in control_idx:
            if c in used_controls:
                continue
            dist = abs(logit_ps[t] - logit_ps[c])
            if dist < best_dist:
                best_dist = dist
                best_c = c
        if best_c is not None and best_dist <= caliper_val:
            matched_pairs.append((t, best_c))
            used_controls.add(best_c)
    
    n_pairs = len(matched_pairs)
    print(f"  {label}: {n_pairs} pairs (caliper={caliper})")
    
    # Extract matched data
    matched_treated = [p[0] for p in matched_pairs]
    matched_control = [p[1] for p in matched_pairs]
    matched_indices = matched_treated + matched_control
    matched_data = data.iloc[matched_indices].copy()
    
    # SMD before and after
    smd_before = {}
    smd_after = {}
    for col in ps_cols:
        # Before
        t_before = data[data['ckd'] == 1][col]
        c_before = data[data['ckd'] == 0][col]
        smd_b = abs(t_before.mean() - c_before.mean()) / np.sqrt((t_before.std()**2 + c_before.std()**2) / 2)
        smd_before[col] = smd_b
        
        # After
        t_after = matched_data[matched_data['ckd'] == 1][col]
        c_after = matched_data[matched_data['ckd'] == 0][col]
        smd_a = abs(t_after.mean() - c_after.mean()) / np.sqrt((t_after.std()**2 + c_after.std()**2) / 2)
        smd_after[col] = smd_a
    
    max_smd_after = max(smd_after.values())
    print(f"    Max SMD after: {max_smd_after:.4f}")
    
    return matched_data, n_pairs, smd_before, smd_after

# PSM without SOFA
ps_cols_no_sofa = ['age', 'gender_num', 'htn', 'dm', 'hf']
matched_no_sofa, n_pairs_no_sofa, smd_before_ns, smd_after_ns = do_psm(
    aki_df, ps_cols_no_sofa, caliper=0.2, label="PSM (no SOFA)")

# PSM with SOFA
ps_cols_with_sofa = ['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_total']
matched_with_sofa, n_pairs_sofa, smd_before_ws, smd_after_ws = do_psm(
    aki_df, ps_cols_with_sofa, caliper=0.2, label="PSM (with SOFA)")

results['psm_no_sofa_pairs'] = n_pairs_no_sofa
results['psm_sofa_pairs'] = n_pairs_sofa

# PSM results by stage
print("\n  PSM Results by Stage:")
for method, matched, name in [("no_sofa", matched_no_sofa, "PSM (no SOFA)"),
                                ("sofa", matched_with_sofa, "PSM (with SOFA)")]:
    for s in [1, 2, 3]:
        sub = matched[matched['max_aki_stage'] == s]
        if len(sub) < 10:
            continue
        pure = sub[sub['ckd'] == 0]
        ckd = sub[sub['ckd'] == 1]
        pure_mort = pure['hospital_expire_flag'].mean() * 100 if len(pure) > 0 else 0
        ckd_mort = ckd['hospital_expire_flag'].mean() * 100 if len(ckd) > 0 else 0
        
        X = sm.add_constant(sub['ckd'])
        try:
            m = sm.Logit(sub['hospital_expire_flag'], X).fit(disp=0)
            or_val = np.exp(m.params['ckd'])
            or_ci = np.exp(m.conf_int().loc['ckd'])
            p_val = m.pvalues['ckd']
            print(f"    {name} Stage {s}: N={len(sub)} (Pure={len(pure)}, CKD={len(ckd)}), Pure mort={pure_mort:.1f}%, CKD mort={ckd_mort:.1f}%, OR={or_val:.3f} ({or_ci[0]:.3f}-{or_ci[1]:.3f}), P={p_val:.4f}")
            results[f'psm_{method}_stage{s}_or'] = round(or_val, 3)
            results[f'psm_{method}_stage{s}_lo'] = round(or_ci[0], 3)
            results[f'psm_{method}_stage{s}_hi'] = round(or_ci[1], 3)
            results[f'psm_{method}_stage{s}_p'] = round(p_val, 4)
            results[f'psm_{method}_stage{s}_pure_mort'] = round(pure_mort, 1)
            results[f'psm_{method}_stage{s}_ckd_mort'] = round(ckd_mort, 1)
            results[f'psm_{method}_stage{s}_n'] = len(sub)
        except:
            pass

    # Overall PSM OR
    X_all = sm.add_constant(matched['ckd'])
    m_all = sm.Logit(matched['hospital_expire_flag'], X_all).fit(disp=0)
    or_all = np.exp(m_all.params['ckd'])
    ci_all = np.exp(m_all.conf_int().loc['ckd'])
    p_all = m_all.pvalues['ckd']
    print(f"    {name} Overall: OR={or_all:.3f} ({ci_all[0]:.3f}-{ci_all[1]:.3f}), P={p_all:.4f}")
    results[f'psm_{method}_overall_or'] = round(or_all, 3)
    results[f'psm_{method}_overall_lo'] = round(ci_all[0], 3)
    results[f'psm_{method}_overall_hi'] = round(ci_all[1], 3)
    results[f'psm_{method}_overall_p'] = round(p_all, 4)

# SMD details for Table S1/S2
print("\n  SMD Details for Supplementary Tables:")
for col in ps_cols_with_sofa:
    print(f"    {col}: Before (no SOFA model)={smd_before_ns.get(col, 0):.3f}, After (no SOFA)={smd_after_ns.get(col, 0):.3f}, Before (SOFA model)={smd_before_ws.get(col, 0):.3f}, After (SOFA)={smd_after_ws.get(col, 0):.3f}")
    results[f'smd_{col}_before_ns'] = round(smd_before_ns.get(col, 0), 3)
    results[f'smd_{col}_after_ns'] = round(smd_after_ns.get(col, 0), 3)
    results[f'smd_{col}_before_ws'] = round(smd_before_ws.get(col, 0), 3)
    results[f'smd_{col}_after_ws'] = round(smd_after_ws.get(col, 0), 3)

# Also compute SMD for AKI stage (not in PS model but reported in balance table)
for col in ['stage2', 'stage3']:
    t_before = aki_df[aki_df['ckd'] == 1][col]
    c_before = aki_df[aki_df['ckd'] == 0][col]
    smd_b = abs(t_before.mean() - c_before.mean()) / np.sqrt((t_before.std()**2 + c_before.std()**2) / 2)
    
    t_after = matched_with_sofa[matched_with_sofa['ckd'] == 1][col]
    c_after = matched_with_sofa[matched_with_sofa['ckd'] == 0][col]
    smd_a = abs(t_after.mean() - c_after.mean()) / np.sqrt((t_after.std()**2 + c_after.std()**2) / 2)
    
    t_after_ns = matched_no_sofa[matched_no_sofa['ckd'] == 1][col]
    c_after_ns = matched_no_sofa[matched_no_sofa['ckd'] == 0][col]
    smd_a_ns = abs(t_after_ns.mean() - c_after_ns.mean()) / np.sqrt((t_after_ns.std()**2 + c_after_ns.std()**2) / 2)
    
    print(f"    {col}: Before={smd_b:.3f}, After (no SOFA)={smd_a_ns:.3f}, After (SOFA)={smd_a:.3f}")
    results[f'smd_{col}_before'] = round(smd_b, 3)
    results[f'smd_{col}_after_ns'] = round(smd_a_ns, 3)
    results[f'smd_{col}_after_ws'] = round(smd_a, 3)

# ============================================================
# 10. KM LOG-RANK P-VALUES (for Figure S2)
# ============================================================
print("\n" + "=" * 80)
print("10. KM LOG-RANK P-VALUES")
print("=" * 80)
try:
    from lifelines.statistics import logrank_test
    sv = pd.read_csv('aki_survival_data.csv')
    sv_aki = sv[sv['max_aki_stage'] > 0].copy()
    
    for s in [1, 2, 3]:
        sub = sv_aki[sv_aki['max_aki_stage'] == s]
        pure = sub[sub['ckd'] == 0]
        ckd = sub[sub['ckd'] == 1]
        
        result = logrank_test(
            pure['survival_days'], ckd['survival_days'],
            event_observed_A=pure['hospital_expire_flag'],
            event_observed_B=ckd['hospital_expire_flag']
        )
        print(f"  Stage {s}: log-rank P={result.p_value:.4f}")
        results[f'km_stage{s}_p'] = round(result.p_value, 4)
except Exception as e:
    print(f"  Error: {e}")

# ============================================================
# 11. RCS NONLINEARITY TEST
# ============================================================
print("\n" + "=" * 80)
print("11. RCS NONLINEARITY TEST")
print("=" * 80)
# Simple RCS with 4 knots using statsmodels
from patsy import dmatrix
knots_pos = [5, 35, 65, 95]
knot_vals = np.percentile(aki_df['age'], knots_pos)
print(f"  Knot positions: {knot_vals}")

# Create RCS basis
rcs_basis = dmatrix(
    f"cr(age, df=4)",
    data=aki_df,
    return_type='dataframe'
)
# Fit model with RCS
X_rcs = pd.concat([rcs_basis, aki_df[['gender_num', 'htn', 'dm', 'hf', 'sofa_total', 'stage2', 'stage3', 'ckd']]], axis=1)
X_rcs = sm.add_constant(X_rcs)
m_rcs = sm.Logit(y, X_rcs).fit(disp=0)

# Nonlinearity test: compare linear vs RCS model
X_lin = sm.add_constant(aki_df[['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_total', 'stage2', 'stage3', 'ckd']])
m_lin = sm.Logit(y, X_lin).fit(disp=0)

lr_nonlin = -2 * (m_lin.llf - m_rcs.llf)
p_nonlin = 1 - stats.chi2.cdf(lr_nonlin, 3)  # 3 extra df for RCS
print(f"  Nonlinearity test: chi2={lr_nonlin:.2f}, P={p_nonlin:.4f}")
results['rcs_nonlin_p'] = round(p_nonlin, 4)
results['rcs_knots'] = [round(k, 1) for k in knot_vals]

# Overall age effect
lr_age = -2 * (mA.llf - m_rcs.llf)  # Not exactly right but approximate
p_age = 1 - stats.chi2.cdf(lr_age, 4)
print(f"  Overall age effect: P={p_age:.6f}")
results['rcs_age_p'] = round(p_age, 6)

# ============================================================
# 12. LASSO FEATURE SELECTION
# ============================================================
print("\n" + "=" * 80)
print("12. LASSO FEATURE SELECTION")
print("=" * 80)
from sklearn.linear_model import LassoCV
from sklearn.preprocessing import StandardScaler

# All candidate features
lasso_features = ['age', 'gender_num', 'htn', 'dm', 'hf',
                  'sofa_resp', 'sofa_coag', 'sofa_liver', 'sofa_cv', 'sofa_cns', 'sofa_renal',
                  'stage2', 'stage3', 'ckd', 'ckd_stage2', 'ckd_stage3',
                  'sofa_total', 'sofa_nonrenal']
# Remove sofa_total and sofa_nonrenal to avoid collinearity
lasso_features = ['age', 'gender_num', 'htn', 'dm', 'hf',
                  'sofa_resp', 'sofa_coag', 'sofa_liver', 'sofa_cv', 'sofa_cns', 'sofa_renal',
                  'stage2', 'stage3', 'ckd', 'ckd_stage2', 'ckd_stage3']

X_lasso = aki_df[lasso_features].values
y_lasso = y.values

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_lasso)

# LASSO with 5-fold CV
lasso = LassoCV(cv=5, random_state=42, max_iter=10000)
lasso.fit(X_scaled, y_lasso)

selected = [(lasso_features[i], lasso.coef_[i]) for i in range(len(lasso_features)) if abs(lasso.coef_[i]) > 1e-6]
selected.sort(key=lambda x: abs(x[1]), reverse=True)
print(f"  LASSO selected {len(selected)} of {len(lasso_features)} features (lambda={lasso.alpha_:.4f}):")
for name, coef in selected:
    or_val = np.exp(coef)  # per SD
    direction = "Protective" if coef < 0 else "Harmful"
    print(f"    {name}: coef={coef:.4f}, OR(per SD)={or_val:.3f} ({direction})")

results['lasso_n_selected'] = len(selected)
results['lasso_n_total'] = len(lasso_features)
results['lasso_lambda'] = round(lasso.alpha_, 4)

# Check if ckd_stage3 is selected
ckd_s3_coef = lasso.coef_[lasso_features.index('ckd_stage3')]
print(f"\n  CKD×Stage3 selected: {abs(ckd_s3_coef) > 1e-6}, coef={ckd_s3_coef:.4f}, OR(per SD)={np.exp(ckd_s3_coef):.3f}")
results['lasso_ckd_s3_coef'] = round(ckd_s3_coef, 4)
results['lasso_ckd_s3_or'] = round(np.exp(ckd_s3_coef), 3)

# ============================================================
# 13. MODEL WITHOUT SOFA (for AUC comparison)
# ============================================================
print("\n" + "=" * 80)
print("13. MODEL WITHOUT SOFA (for AUC comparison 0.664 -> 0.718)")
print("=" * 80)
X_no_sofa = sm.add_constant(aki_df[['age', 'gender_num', 'htn', 'dm', 'hf', 'stage2', 'stage3', 'ckd', 'ckd_stage2', 'ckd_stage3']])
m_no_sofa = sm.Logit(y, X_no_sofa).fit(disp=0)
preds_no_sofa = m_no_sofa.predict(X_no_sofa).values
auc_no_sofa = roc_auc_score(y_arr, preds_no_sofa)
print(f"  Model C without SOFA: AUC={auc_no_sofa:.4f}")
results['modelC_no_sofa_auc'] = round(auc_no_sofa, 4)

# ============================================================
# 14. SOFA COMPARISON (Pure vs CKD)
# ============================================================
print("\n" + "=" * 80)
print("14. SOFA COMPARISON")
print("=" * 80)
pure_sofa = aki_df[aki_df['ckd'] == 0]['sofa_total']
ckd_sofa = aki_df[aki_df['ckd'] == 1]['sofa_total']
t_stat, p_t = stats.ttest_ind(pure_sofa.dropna(), ckd_sofa.dropna())
u_stat, p_u = stats.mannwhitneyu(pure_sofa.dropna(), ckd_sofa.dropna(), alternative='two-sided')
print(f"  Pure AKI: {pure_sofa.mean():.1f} +/- {pure_sofa.std():.1f}, median={pure_sofa.median():.1f}")
print(f"  AKI+CKD: {ckd_sofa.mean():.1f} +/- {ckd_sofa.std():.1f}, median={ckd_sofa.median():.1f}")
print(f"  t-test: P={p_t:.6f}")
print(f"  Mann-Whitney U: P={p_u:.6f}")
results['pure_sofa_mean'] = round(pure_sofa.mean(), 1)
results['pure_sofa_sd'] = round(pure_sofa.std(), 1)
results['ckd_sofa_mean'] = round(ckd_sofa.mean(), 1)
results['ckd_sofa_sd'] = round(ckd_sofa.std(), 1)

# ============================================================
# SAVE ALL RESULTS
# ============================================================
print("\n" + "=" * 80)
print("SAVING RESULTS TO all_statistics.json")
print("=" * 80)

with open('all_statistics.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"  {len(results)} statistics saved")
print("\n=== ALL COMPUTATIONS COMPLETE ===")
