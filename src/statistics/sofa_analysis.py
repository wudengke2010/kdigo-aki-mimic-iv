#!/usr/bin/env python3
"""
AKI Paradox: Severity-Adjusted Re-Analysis with SOFA
1. Logistic regression with SOFA (Model comparison)
2. PSM with SOFA
3. CKD x AKI stage interaction with SOFA
4. SOFA-stratified sensitivity analysis
5. Non-renal SOFA analysis
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import chi2_contingency, mannwhitneyu
import statsmodels.api as sm
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 150

OUTPUT_DIR = "C:/Users/admin/WorkBuddy/2026-07-07-19-40-19"
FIG_DIR = os.path.join(OUTPUT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

COLOR_PURE = '#2196F3'
COLOR_CKD = '#FF5722'
COLOR_SOFA = '#4CAF50'

def log(msg):
    print(msg, flush=True)

log("=" * 70)
log("AKI Paradox: Severity-Adjusted Re-Analysis with SOFA")
log("=" * 70)

# ============================================================
# Load data
# ============================================================
df = pd.read_csv(os.path.join(OUTPUT_DIR, "aki_paradox_cohort_with_sofa.csv"))
aki = df[df["max_aki_stage"] >= 1].copy()
aki["male"] = (aki["gender"] == "M").astype(int)
aki["aki_stage2"] = (aki["max_aki_stage"] == 2).astype(int)
aki["aki_stage3"] = (aki["max_aki_stage"] == 3).astype(int)
aki["mortality"] = aki["hospital_expire_flag"]

log(f"Loaded {len(aki):,} AKI patients")
log(f"SOFA: mean={aki['sofa_total'].mean():.1f}, median={aki['sofa_total'].median():.1f}")

# ============================================================
# 1. Logistic Regression: Model Comparison
# ============================================================
log("\n" + "=" * 70)
log("1. Logistic Regression: Model Comparison (with/without SOFA)")
log("=" * 70)

results = {}

# Model 1: Original (age, gender, comorbidities, AKI stage)
X1_cols = ["age", "male", "htn", "dm", "hf", "aki_stage2", "aki_stage3", "ckd"]
X1 = sm.add_constant(aki[X1_cols])
y = aki["mortality"]
m1 = sm.Logit(y, X1).fit(disp=0)
log(f"\nModel 1 (Original): AUC-like pseudo-R2 = {m1.prsquared:.4f}")
log(f"  CKD OR = {np.exp(m1.params['ckd']):.3f} ({np.exp(m1.conf_int().loc['ckd', 0]):.3f}-{np.exp(m1.conf_int().loc['ckd', 1]):.3f}), P = {m1.pvalues['ckd']:.4f}")

# Model 2: + SOFA total
X2_cols = X1_cols + ["sofa_total"]
X2 = sm.add_constant(aki[X2_cols])
m2 = sm.Logit(y, X2).fit(disp=0)
log(f"\nModel 2 (+SOFA): pseudo-R2 = {m2.prsquared:.4f}")
log(f"  CKD OR = {np.exp(m2.params['ckd']):.3f} ({np.exp(m2.conf_int().loc['ckd', 0]):.3f}-{np.exp(m2.conf_int().loc['ckd', 1]):.3f}), P = {m2.pvalues['ckd']:.4f}")
log(f"  SOFA OR = {np.exp(m2.params['sofa_total']):.3f} ({np.exp(m2.conf_int().loc['sofa_total', 0]):.3f}-{np.exp(m2.conf_int().loc['sofa_total', 1]):.3f}), P = {m2.pvalues['sofa_total']:.6f}")

# Model 3: + SOFA + interaction
aki["ckd_stage3"] = aki["ckd"] * aki["aki_stage3"]
aki["ckd_stage2"] = aki["ckd"] * aki["aki_stage2"]
X3 = sm.add_constant(aki[X2_cols + ["ckd_stage2", "ckd_stage3"]])
m3 = sm.Logit(y, X3).fit(disp=0)
log(f"\nModel 3 (+SOFA+Interaction): pseudo-R2 = {m3.prsquared:.4f}")
log(f"  CKD OR = {np.exp(m3.params['ckd']):.3f}, P = {m3.pvalues['ckd']:.4f}")
log(f"  CKD x Stage3 OR = {np.exp(m3.params['ckd_stage3']):.3f} ({np.exp(m3.conf_int().loc['ckd_stage3', 0]):.3f}-{np.exp(m3.conf_int().loc['ckd_stage3', 1]):.3f}), P = {m3.pvalues['ckd_stage3']:.6f}")

# LR test: Model 3 vs Model 2 (interaction significance)
lr_stat = -2 * (m2.llf - m3.llf)
lr_p = stats.chi2.sf(lr_stat, df=2)  # 2 interaction terms
log(f"\n  LR test (interaction): chi2 = {lr_stat:.2f}, df=2, P = {lr_p:.8f}")

# Model 4: + Non-renal SOFA (to remove renal-AKI collinearity)
X4_cols = X1_cols + ["sofa_nonrenal"]
X4 = sm.add_constant(aki[X4_cols + ["ckd_stage2", "ckd_stage3"]])
m4 = sm.Logit(y, X4).fit(disp=0)
log(f"\nModel 4 (+Non-renal SOFA+Interaction): pseudo-R2 = {m4.prsquared:.4f}")
log(f"  CKD x Stage3 OR = {np.exp(m4.params['ckd_stage3']):.3f} ({np.exp(m4.conf_int().loc['ckd_stage3', 0]):.3f}-{np.exp(m4.conf_int().loc['ckd_stage3', 1]):.3f}), P = {m4.pvalues['ckd_stage3']:.6f}")

# AUC comparison
from sklearn.metrics import roc_auc_score
auc1 = roc_auc_score(y, m1.predict(X1))
auc2 = roc_auc_score(y, m2.predict(X2))
auc3 = roc_auc_score(y, m3.predict(X3))
auc4 = roc_auc_score(y, m4.predict(X4))
log(f"\n  AUC Model 1 (original): {auc1:.4f}")
log(f"  AUC Model 2 (+SOFA): {auc2:.4f}")
log(f"  AUC Model 3 (+SOFA+interaction): {auc3:.4f}")
log(f"  AUC Model 4 (+NonrenalSOFA+interaction): {auc4:.4f}")

# Save model comparison table
model_comp = pd.DataFrame({
    "Model": ["M1: Original", "M2: +SOFA", "M3: +SOFA+Interaction", "M4: +NonrenalSOFA+Interaction"],
    "Variables": ["age,sex,comorbidities,AKI stage,CKD", "+SOFA total", "+CKD x Stage interaction", "+Non-renal SOFA + interaction"],
    "CKD_OR": [np.exp(m1.params['ckd']), np.exp(m2.params['ckd']), np.exp(m3.params['ckd']), np.exp(m4.params['ckd'])],
    "CKD_P": [m1.pvalues['ckd'], m2.pvalues['ckd'], m3.pvalues['ckd'], m4.pvalues['ckd']],
    "CKDxStage3_OR": [np.nan, np.nan, np.exp(m3.params['ckd_stage3']), np.exp(m4.params['ckd_stage3'])],
    "CKDxStage3_P": [np.nan, np.nan, m3.pvalues['ckd_stage3'], m4.pvalues['ckd_stage3']],
    "LR_test_P": [np.nan, np.nan, lr_p, np.nan],
    "Pseudo_R2": [m1.prsquared, m2.prsquared, m3.prsquared, m4.prsquared],
    "AUC": [auc1, auc2, auc3, auc4]
})
model_comp.to_csv(os.path.join(OUTPUT_DIR, "table_sofa_model_comparison.csv"), index=False)
log(f"\n  Saved: table_sofa_model_comparison.csv")

# ============================================================
# 2. Subgroup Analysis by Stage (with SOFA adjustment)
# ============================================================
log("\n" + "=" * 70)
log("2. Subgroup Analysis by AKI Stage (SOFA-adjusted)")
log("=" * 70)

subgroup_results = []
for stage in [1, 2, 3]:
    stage_data = aki[aki["max_aki_stage"] == stage].copy()
    n_pure = len(stage_data[stage_data["ckd"] == 0])
    n_ckd = len(stage_data[stage_data["ckd"] == 1])
    mort_pure = stage_data[stage_data["ckd"] == 0]["mortality"].mean() * 100
    mort_ckd = stage_data[stage_data["ckd"] == 1]["mortality"].mean() * 100
    sofa_pure = stage_data[stage_data["ckd"] == 0]["sofa_total"].mean()
    sofa_ckd = stage_data[stage_data["ckd"] == 1]["sofa_total"].mean()
    
    # Unadjusted OR
    a = stage_data[(stage_data["ckd"] == 1) & (stage_data["mortality"] == 1)].shape[0]
    b = stage_data[(stage_data["ckd"] == 1) & (stage_data["mortality"] == 0)].shape[0]
    c = stage_data[(stage_data["ckd"] == 0) & (stage_data["mortality"] == 1)].shape[0]
    d = stage_data[(stage_data["ckd"] == 0) & (stage_data["mortality"] == 0)].shape[0]
    or_unadj = (a * d) / (b * c) if a * d > 0 else np.nan
    
    # SOFA-adjusted OR
    X = sm.add_constant(stage_data[["age", "male", "htn", "dm", "hf", "sofa_total", "ckd"]])
    y_s = stage_data["mortality"]
    try:
        m = sm.Logit(y_s, X).fit(disp=0)
        or_adj = np.exp(m.params['ckd'])
        ci_lo = np.exp(m.conf_int().loc['ckd', 0])
        ci_hi = np.exp(m.conf_int().loc['ckd', 1])
        p_val = m.pvalues['ckd']
    except:
        or_adj = ci_lo = ci_hi = p_val = np.nan
    
    subgroup_results.append({
        "Stage": stage,
        "N_PureAKI": n_pure, "N_AKI_CKD": n_ckd,
        "Mort_PureAKI": f"{mort_pure:.1f}%", "Mort_AKI_CKD": f"{mort_ckd:.1f}%",
        "SOFA_PureAKI": f"{sofa_pure:.1f}", "SOFA_AKI_CKD": f"{sofa_ckd:.1f}",
        "OR_unadjusted": round(or_unadj, 3),
        "OR_SOFA_adjusted": round(or_adj, 3),
        "CI_lower": round(ci_lo, 3), "CI_upper": round(ci_hi, 3),
        "P_value": f"{p_val:.4f}" if not np.isnan(p_val) else "NA"
    })
    
    log(f"\n  Stage {stage} (n={len(stage_data):,}):")
    log(f"    Pure AKI: n={n_pure:,}, mort={mort_pure:.1f}%, SOFA={sofa_pure:.1f}")
    log(f"    AKI+CKD:  n={n_ckd:,}, mort={mort_ckd:.1f}%, SOFA={sofa_ckd:.1f}")
    log(f"    OR unadjusted: {or_unadj:.3f}")
    log(f"    OR SOFA-adjusted: {or_adj:.3f} ({ci_lo:.3f}-{ci_hi:.3f}), P={p_val:.4f}")

subgroup_df = pd.DataFrame(subgroup_results)
subgroup_df.to_csv(os.path.join(OUTPUT_DIR, "table_sofa_subgroup_by_stage.csv"), index=False)
log(f"\n  Saved: table_sofa_subgroup_by_stage.csv")

# ============================================================
# 3. PSM with SOFA
# ============================================================
log("\n" + "=" * 70)
log("3. PSM with SOFA")
log("=" * 70)

psm_covariates = ["age", "male", "htn", "dm", "hf", "sofa_total", 
                   "aki_stage2", "aki_stage3"]

# Standardize
scaler = StandardScaler()
aki_psm = aki.copy()
aki_psm[psm_covariates] = scaler.fit_transform(aki_psm[psm_covariates])

# Propensity score
ps_model = LogisticRegression(max_iter=1000, random_state=42)
ps_model.fit(aki_psm[psm_covariates], aki_psm["ckd"])
aki_psm["ps"] = ps_model.predict_proba(aki_psm[psm_covariates])[:, 1]

# 1:1 nearest neighbor matching
ckd_pts = aki_psm[aki_psm["ckd"] == 1].copy()
pure_pts = aki_psm[aki_psm["ckd"] == 0].copy()

nn = NearestNeighbors(n_neighbors=1, algorithm='ball_tree')
nn.fit(pure_pts[["ps"]].values)
distances, indices = nn.kneighbors(ckd_pts[["ps"]].values)

matched_pure_idx = pure_pts.iloc[indices.flatten()].index
matched_ckd_idx = ckd_pts.index

matched = pd.concat([aki_psm.loc[matched_pure_idx], aki_psm.loc[matched_ckd_idx]])
log(f"  PSM matched pairs: {len(ckd_pts):,} (caliper not applied)")
log(f"  Total matched: {len(matched):,}")

# Check covariate balance
log("\n  Covariate balance after PSM:")
for cov in psm_covariates:
    pure_val = matched.loc[matched_pure_idx, cov].values
    ckd_val = matched.loc[matched_ckd_idx, cov].values
    # SMD (using original values, not standardized)
    orig_pure = aki.loc[matched_pure_idx, cov if cov in aki.columns else cov.replace("aki_stage", "max_aki_stage")].values
    orig_ckd = aki.loc[matched_ckd_idx, cov if cov in aki.columns else cov.replace("aki_stage", "max_aki_stage")].values
    
    if cov.startswith("aki_stage") or cov in ["male", "htn", "dm", "hf"]:
        p1 = orig_pure.mean()
        p2 = orig_ckd.mean()
        smd = abs(p1 - p2) / np.sqrt((p1*(1-p1) + p2*(1-p2)) / 2) if (p1*(1-p1) + p2*(1-p2)) > 0 else 0
    else:
        smd = abs(orig_pure.mean() - orig_ckd.mean()) / np.sqrt((orig_pure.std()**2 + orig_ckd.std()**2) / 2)
    
    log(f"    {cov}: SMD = {smd:.3f} {'OK' if smd < 0.1 else 'UNBALANCED'}")

# PSM results by stage
log("\n  PSM results by stage:")
psm_results = []
for stage in [1, 2, 3]:
    for label, ckd_val in [("Overall", None)]:
        if ckd_val is None:
            stage_matched = matched[matched["max_aki_stage"] == stage]
        pure_m = stage_matched[stage_matched["ckd"] == 0]
        ckd_m = stage_matched[stage_matched["ckd"] == 1]
        
        if len(pure_m) > 0 and len(ckd_m) > 0:
            mort_p = pure_m["mortality"].mean() * 100
            mort_c = ckd_m["mortality"].mean() * 100
            
            a = len(ckd_m[ckd_m["mortality"] == 1])
            b = len(ckd_m[ckd_m["mortality"] == 0])
            c = len(pure_m[pure_m["mortality"] == 1])
            d = len(pure_m[pure_m["mortality"] == 0])
            
            if a * d > 0 and b * c > 0:
                or_val = (a * d) / (b * c)
                se = np.sqrt(1/max(a,0.5) + 1/max(b,0.5) + 1/max(c,0.5) + 1/max(d,0.5))
                ci_lo = np.exp(np.log(or_val) - 1.96 * se)
                ci_hi = np.exp(np.log(or_val) + 1.96 * se)
                z = np.log(or_val) / se
                p = 2 * stats.norm.sf(abs(z))
            else:
                or_val = ci_lo = ci_hi = p = np.nan
            
            log(f"    Stage {stage}: Pure={mort_p:.1f}% vs CKD={mort_c:.1f}%, "
                f"OR={or_val:.3f} ({ci_lo:.3f}-{ci_hi:.3f}), P={p:.4f}, n={len(ckd_m):,} pairs")
            
            psm_results.append({
                "Stage": stage,
                "N_pairs": len(ckd_m),
                "Mort_PureAKI": f"{mort_p:.1f}%",
                "Mort_AKI_CKD": f"{mort_c:.1f}%",
                "OR": round(or_val, 3),
                "CI_lower": round(ci_lo, 3),
                "CI_upper": round(ci_hi, 3),
                "P_value": f"{p:.4f}"
            })

# Overall PSM result
pure_all = matched[matched["ckd"] == 0]
ckd_all = matched[matched["ckd"] == 1]
mort_p = pure_all["mortality"].mean() * 100
mort_c = ckd_all["mortality"].mean() * 100
a = len(ckd_all[ckd_all["mortality"] == 1])
b = len(ckd_all[ckd_all["mortality"] == 0])
c = len(pure_all[pure_all["mortality"] == 1])
d = len(pure_all[pure_all["mortality"] == 0])
or_all = (a * d) / (b * c)
se = np.sqrt(1/a + 1/b + 1/c + 1/d)
ci_lo = np.exp(np.log(or_all) - 1.96 * se)
ci_hi = np.exp(np.log(or_all) + 1.96 * se)
z = np.log(or_all) / se
p_all = 2 * stats.norm.sf(abs(z))
log(f"\n  Overall: Pure={mort_p:.1f}% vs CKD={mort_c:.1f}%, "
    f"OR={or_all:.3f} ({ci_lo:.3f}-{ci_hi:.3f}), P={p_all:.4f}")

psm_df = pd.DataFrame(psm_results)
psm_df.to_csv(os.path.join(OUTPUT_DIR, "table_sofa_psm_results.csv"), index=False)
log(f"\n  Saved: table_sofa_psm_results.csv")

# ============================================================
# 4. SOFA-Stratified Sensitivity Analysis
# ============================================================
log("\n" + "=" * 70)
log("4. SOFA-Stratified Sensitivity Analysis")
log("=" * 70)

# Stratify by SOFA tertiles
sofa_tertiles = pd.qcut(aki["sofa_total"], q=3, labels=["Low (0-5)", "Medium (6-10)", "High (>10)"])
aki["sofa_stratum"] = sofa_tertiles

log("\n  SOFA strata:")
for stratum in ["Low (0-5)", "Medium (6-10)", "High (>10)"]:
    stratum_data = aki[aki["sofa_stratum"] == stratum]
    log(f"    {stratum}: n={len(stratum_data):,}, "
        f"SOFA range={stratum_data['sofa_total'].min()}-{stratum_data['sofa_total'].max()}, "
        f"mortality={stratum_data['mortality'].mean()*100:.1f}%")

log("\n  Stage 3 paradox by SOFA stratum:")
stratified_results = []
for stratum in ["Low (0-5)", "Medium (6-10)", "High (>10)"]:
    stratum_data = aki[(aki["sofa_stratum"] == stratum) & (aki["max_aki_stage"] == 3)]
    if len(stratum_data) < 20:
        log(f"    {stratum} Stage 3: too few patients (n={len(stratum_data)})")
        continue
    
    pure = stratum_data[stratum_data["ckd"] == 0]
    ckd = stratum_data[stratum_data["ckd"] == 1]
    
    if len(pure) > 0 and len(ckd) > 0:
        mort_p = pure["mortality"].mean() * 100
        mort_c = ckd["mortality"].mean() * 100
        
        a = len(ckd[ckd["mortality"] == 1])
        b = len(ckd[ckd["mortality"] == 0])
        c = len(pure[pure["mortality"] == 1])
        d = len(pure[pure["mortality"] == 0])
        
        if a * d > 0 and b * c > 0:
            or_val = (a * d) / (b * c)
            se = np.sqrt(1/max(a,0.5) + 1/max(b,0.5) + 1/max(c,0.5) + 1/max(d,0.5))
            ci_lo = np.exp(np.log(or_val) - 1.96 * se)
            ci_hi = np.exp(np.log(or_val) + 1.96 * se)
            z = np.log(or_val) / se
            p = 2 * stats.norm.sf(abs(z))
        else:
            or_val = ci_lo = ci_hi = p = np.nan
        
        log(f"    {stratum} Stage 3 (n={len(stratum_data):,}): "
            f"Pure={mort_p:.1f}% (n={len(pure)}), CKD={mort_c:.1f}% (n={len(ckd)}), "
            f"OR={or_val:.3f} ({ci_lo:.3f}-{ci_hi:.3f}), P={p:.4f}")
        
        stratified_results.append({
            "SOFA_stratum": stratum,
            "N_total": len(stratum_data),
            "N_PureAKI": len(pure),
            "N_AKI_CKD": len(ckd),
            "Mort_PureAKI": f"{mort_p:.1f}%",
            "Mort_AKI_CKD": f"{mort_c:.1f}%",
            "OR": round(or_val, 3),
            "CI_lower": round(ci_lo, 3),
            "CI_upper": round(ci_hi, 3),
            "P_value": f"{p:.4f}"
        })

strat_df = pd.DataFrame(stratified_results)
strat_df.to_csv(os.path.join(OUTPUT_DIR, "table_sofa_stratified_s3.csv"), index=False)

# Also for all stages
log("\n  All stages by SOFA stratum:")
for stratum in ["Low (0-5)", "Medium (6-10)", "High (>10)"]:
    for stage in [1, 2, 3]:
        stratum_data = aki[(aki["sofa_stratum"] == stratum) & (aki["max_aki_stage"] == stage)]
        if len(stratum_data) < 20:
            continue
        pure = stratum_data[stratum_data["ckd"] == 0]
        ckd = stratum_data[stratum_data["ckd"] == 1]
        if len(pure) > 0 and len(ckd) > 0:
            mort_p = pure["mortality"].mean() * 100
            mort_c = ckd["mortality"].mean() * 100
            log(f"    {stratum} S{stage}: Pure={mort_p:.1f}%(n={len(pure)}), CKD={mort_c:.1f}%(n={len(ckd)})")

log(f"\n  Saved: table_sofa_stratified_s3.csv")

# ============================================================
# 5. Summary
# ============================================================
log("\n" + "=" * 70)
log("SUMMARY: Severity-Adjusted Analysis")
log("=" * 70)

log(f"""
KEY FINDINGS:

1. SOFA distribution:
   - Pure AKI: mean={aki[aki['ckd']==0]['sofa_total'].mean():.1f}
   - AKI+CKD:  mean={aki[aki['ckd']==1]['sofa_total'].mean():.1f}

2. Stage 3 paradox - SOFA-adjusted OR:
   - Unadjusted OR: {subgroup_results[2]['OR_unadjusted']:.3f}
   - SOFA-adjusted OR: {subgroup_results[2]['OR_SOFA_adjusted']:.3f} ({subgroup_results[2]['CI_lower']:.3f}-{subgroup_results[2]['CI_upper']:.3f}), P={subgroup_results[2]['P_value']}

3. Interaction test with SOFA:
   - CKD x Stage3 OR: {np.exp(m3.params['ckd_stage3']):.3f}, P={m3.pvalues['ckd_stage3']:.6f}
   - LR test P: {lr_p:.8f}

4. Model performance:
   - AUC without SOFA: {auc1:.4f}
   - AUC with SOFA: {auc2:.4f}
   - AUC with SOFA + interaction: {auc3:.4f}

5. PSM with SOFA:
   - Overall OR: {or_all:.3f} ({ci_lo:.3f}-{ci_hi:.3f}), P={p_all:.4f}
""")

log("Done! All severity-adjusted analyses complete.")
