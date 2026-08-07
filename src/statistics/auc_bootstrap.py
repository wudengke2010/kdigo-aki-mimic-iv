"""
Compute AUC CIs, CV-AUC, Stage 3 non-renal SOFA, and SOFA-stratified analysis.
"""
import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
import warnings
warnings.filterwarnings('ignore')

df = pd.read_csv('aki_paradox_cohort_with_sofa.csv')
aki_df = df[df['max_aki_stage'] > 0].copy()
aki_df['gender_num'] = (aki_df['gender'] == 'M').astype(int)
aki_df['stage2'] = (aki_df['max_aki_stage'] == 2).astype(int)
aki_df['stage3'] = (aki_df['max_aki_stage'] == 3).astype(int)
aki_df['ckd_stage2'] = aki_df['ckd'] * aki_df['stage2']
aki_df['ckd_stage3'] = aki_df['ckd'] * aki_df['stage3']
y = aki_df['hospital_expire_flag']

# ============ AUC WITH 95% CI (Bootstrap) ============
print("=" * 80)
print("AUC WITH 95% CI (Bootstrap, n=2000)")
print("=" * 80)

X_A = sm.add_constant(aki_df[['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_total', 'stage2', 'stage3']])
X_B = sm.add_constant(aki_df[['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_total', 'stage2', 'stage3', 'ckd']])
X_C = sm.add_constant(aki_df[['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_total', 'stage2', 'stage3', 'ckd', 'ckd_stage2', 'ckd_stage3']])
X_D = sm.add_constant(aki_df[['age', 'gender_num', 'htn', 'dm', 'hf', 'sofa_resp', 'sofa_coag', 'sofa_liver', 'sofa_cv', 'sofa_cns', 'sofa_renal', 'stage2', 'stage3', 'ckd', 'ckd_stage2', 'ckd_stage3']])

mA = sm.Logit(y, X_A).fit(disp=0)
mB = sm.Logit(y, X_B).fit(disp=0)
mC = sm.Logit(y, X_C).fit(disp=0)
mD = sm.Logit(y, X_D).fit(disp=0)

np.random.seed(42)
y_arr = y.values
n = len(y_arr)

for name, model, X in [("A", mA, X_A), ("B", mB, X_B), ("C", mC, X_C), ("D", mD, X_D)]:
    preds = model.predict(X).values
    auc_point = roc_auc_score(y_arr, preds)
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
    print(f"  Model {name}: AUC = {auc_point:.4f} ({ci_l:.4f} - {ci_u:.4f})")

# ============ 5-FOLD CV-AUC ============
print("\n" + "=" * 80)
print("5-FOLD CV-AUC")
print("=" * 80)

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
    print(f"  Model {name}: CV-AUC = {np.mean(cv_aucs):.4f} +/- {np.std(cv_aucs):.4f}")

# ============ STAGE 3 NON-RENAL SOFA ============
print("\n" + "=" * 80)
print("STAGE 3 NON-RENAL SOFA (verification)")
print("=" * 80)

stage3 = aki_df[aki_df['max_aki_stage'] == 3]
pure_s3 = stage3[stage3['ckd'] == 0]
ckd_s3 = stage3[stage3['ckd'] == 1]

print(f"  Stage 3 Pure AKI non-renal SOFA: {pure_s3['sofa_nonrenal'].mean():.1f} +/- {pure_s3['sofa_nonrenal'].std():.1f}")
print(f"  Stage 3 AKI+CKD non-renal SOFA: {ckd_s3['sofa_nonrenal'].mean():.1f} +/- {ckd_s3['sofa_nonrenal'].std():.1f}")
t_stat, p_val = stats.ttest_ind(pure_s3['sofa_nonrenal'].dropna(), ckd_s3['sofa_nonrenal'].dropna())
print(f"  P-value: {p_val:.6f}")

# ============ SOFA-STRATIFIED ANALYSIS (Stage 3) ============
print("\n" + "=" * 80)
print("SOFA-STRATIFIED ANALYSIS (Stage 3)")
print("=" * 80)

stage3_df = aki_df[aki_df['max_aki_stage'] == 3].copy()
for label, lo, hi in [("Low (0-5)", 0, 5), ("Medium (6-10)", 6, 10), ("High (>10)", 11, 100)]:
    stratum = stage3_df[(stage3_df['sofa_total'] >= lo) & (stage3_df['sofa_total'] <= hi)]
    n = len(stratum)
    pure = stratum[stratum['ckd'] == 0]
    ckd = stratum[stratum['ckd'] == 1]
    pure_mort = pure['hospital_expire_flag'].mean() * 100 if len(pure) > 0 else 0
    ckd_mort = ckd['hospital_expire_flag'].mean() * 100 if len(ckd) > 0 else 0

    if len(pure) > 0 and len(ckd) > 0:
        X = sm.add_constant(stratum['ckd'])
        try:
            model = sm.Logit(stratum['hospital_expire_flag'], X).fit(disp=0)
            or_val = np.exp(model.params['ckd'])
            or_ci = np.exp(model.conf_int().loc['ckd'])
            p_val = model.pvalues['ckd']
            print(f"  {label}: N={n} (Pure={len(pure)}, CKD={len(ckd)}), Pure mort={pure_mort:.1f}%, CKD mort={ckd_mort:.1f}%, OR={or_val:.3f} ({or_ci[0]:.3f}-{or_ci[1]:.3f}), P={p_val:.4f}")
        except Exception as e:
            print(f"  {label}: N={n}, OR computation error: {e}")

# ============ BRIER SCORES ============
print("\n" + "=" * 80)
print("BRIER SCORES")
print("=" * 80)

from sklearn.metrics import brier_score_loss
for name, model, X in [("A", mA, X_A), ("B", mB, X_B), ("C", mC, X_C), ("D", mD, X_D)]:
    preds = model.predict(X)
    brier = brier_score_loss(y, preds)
    print(f"  Model {name}: Brier = {brier:.4f}")

# ============ SOFA COMPARISON P-VALUE ============
print("\n" + "=" * 80)
print("SOFA COMPARISON: Pure AKI vs AKI+CKD")
print("=" * 80)

pure_sofa = aki_df[aki_df['ckd'] == 0]['sofa_total']
ckd_sofa = aki_df[aki_df['ckd'] == 1]['sofa_total']

# t-test
t_stat, p_t = stats.ttest_ind(pure_sofa.dropna(), ckd_sofa.dropna())
print(f"  t-test: t={t_stat:.2f}, P={p_t:.6f}")

# Mann-Whitney U test
u_stat, p_u = stats.mannwhitneyu(pure_sofa.dropna(), ckd_sofa.dropna(), alternative='two-sided')
print(f"  Mann-Whitney U: P={p_u:.6f}")

print(f"  Pure AKI: {pure_sofa.mean():.1f} +/- {pure_sofa.std():.1f}, median={pure_sofa.median():.1f}")
print(f"  AKI+CKD: {ckd_sofa.mean():.1f} +/- {ckd_sofa.std():.1f}, median={ckd_sofa.median():.1f}")

# ============ NRI/IDI ============
print("\n" + "=" * 80)
print("NRI AND IDI")
print("=" * 80)

def compute_nri(y_true, p_old, p_new):
    n = len(y_true)
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
    idi = idi_events - idi_non
    return idi

preds_A = mA.predict(X_A).values
preds_C = mC.predict(X_C).values

nri, nri_e, nri_ne = compute_nri(y.values, preds_A, preds_C)
idi = compute_idi(y.values, preds_A, preds_C)
print(f"  NRI (C vs A): {nri:.3f} (events: {nri_e:.3f}, non-events: {nri_ne:.3f})")
print(f"  IDI (C vs A): {idi:.4f}")

# NRI for B vs A
preds_B = mB.predict(X_B).values
nri_B, nri_Be, nri_Bne = compute_nri(y.values, preds_A, preds_B)
idi_B = compute_idi(y.values, preds_A, preds_B)
print(f"  NRI (B vs A): {nri_B:.3f}")
print(f"  IDI (B vs A): {idi_B:.4f}")

print("\n" + "=" * 80)
print("ALL COMPUTATIONS COMPLETE")
print("=" * 80)
