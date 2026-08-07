#!/usr/bin/env python3
"""
AKI Paradox: Prediction Modeling & Nomogram
1. Build prediction models for AKI mortality (with/without CKD)
2. Compare AUC, NRI, IDI
3. Generate nomogram
4. Calibration curves
"""

import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
from sklearn.linear_model import LogisticRegression, LassoCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve, brier_score_loss
from sklearn.model_selection import cross_val_score, StratifiedKFold
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 150

OUTPUT_DIR = "C:/Users/admin/WorkBuddy/2026-07-07-19-40-19"
FIG_DIR = os.path.join(OUTPUT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

def log(msg):
    print(msg, flush=True)

log("=" * 70)
log("AKI Paradox: Prediction Modeling & Nomogram")
log("=" * 70)

# ============================================================
# Load data
# ============================================================
df = pd.read_csv(os.path.join(OUTPUT_DIR, "aki_paradox_cohort_with_sofa.csv"))
aki = df[df["max_aki_stage"] >= 1].copy()
aki["male"] = (aki["gender"] == "M").astype(int)
aki["aki_stage2"] = (aki["max_aki_stage"] == 2).astype(int)
aki["aki_stage3"] = (aki["max_aki_stage"] == 3).astype(int)
aki["ckd_stage2"] = aki["ckd"] * aki["aki_stage2"]
aki["ckd_stage3"] = aki["ckd"] * aki["aki_stage3"]
aki["mortality"] = aki["hospital_expire_flag"]

log(f"Loaded {len(aki):,} AKI patients, mortality rate: {aki['mortality'].mean()*100:.1f}%")

# ============================================================
# 1. Model Building Strategy
# ============================================================
log("\n" + "=" * 70)
log("1. Model Building & Comparison")
log("=" * 70)

y = aki["mortality"]

# Model A: Base model (demographics + comorbidities + SOFA + AKI stage)
features_A = ["age", "male", "htn", "dm", "hf", "sofa_total", "aki_stage2", "aki_stage3"]
X_A = sm.add_constant(aki[features_A])
m_A = sm.Logit(y, X_A).fit(disp=0)
pred_A = m_A.predict(X_A)
auc_A = roc_auc_score(y, pred_A)

# Model B: + CKD main effect
features_B = features_A + ["ckd"]
X_B = sm.add_constant(aki[features_B])
m_B = sm.Logit(y, X_B).fit(disp=0)
pred_B = m_B.predict(X_B)
auc_B = roc_auc_score(y, pred_B)

# Model C: + CKD x AKI stage interaction
features_C = features_B + ["ckd_stage2", "ckd_stage3"]
X_C = sm.add_constant(aki[features_C])
m_C = sm.Logit(y, X_C).fit(disp=0)
pred_C = m_C.predict(X_C)
auc_C = roc_auc_score(y, pred_C)

# Model D: + SOFA components (instead of total)
features_D = ["age", "male", "htn", "dm", "hf",
              "sofa_resp", "sofa_coag", "sofa_liver", "sofa_cv", "sofa_cns", "sofa_renal",
              "aki_stage2", "aki_stage3", "ckd", "ckd_stage2", "ckd_stage3"]
X_D = sm.add_constant(aki[features_D])
m_D = sm.Logit(y, X_D).fit(disp=0)
pred_D = m_D.predict(X_D)
auc_D = roc_auc_score(y, pred_D)

log(f"\n  Model A (Base):            AUC = {auc_A:.4f}, pseudo-R2 = {m_A.prsquared:.4f}")
log(f"  Model B (+CKD):            AUC = {auc_B:.4f}, pseudo-R2 = {m_B.prsquared:.4f}")
log(f"  Model C (+CKD interaction): AUC = {auc_C:.4f}, pseudo-R2 = {m_C.prsquared:.4f}")
log(f"  Model D (SOFA components):  AUC = {auc_D:.4f}, pseudo-R2 = {m_D.prsquared:.4f}")

# LR tests
lr_BC = -2 * (m_B.llf - m_C.llf)
p_BC = stats.chi2.sf(lr_BC, df=2)
log(f"\n  LR test (B vs C, interaction): chi2={lr_BC:.2f}, df=2, P={p_BC:.8f}")

lr_AB = -2 * (m_A.llf - m_B.llf)
p_AB = stats.chi2.sf(lr_AB, df=1)
log(f"  LR test (A vs B, CKD main): chi2={lr_AB:.2f}, df=1, P={p_AB:.6f}")

# ============================================================
# 2. NRI & IDI Calculation
# ============================================================
log("\n" + "=" * 70)
log("2. NRI & IDI (Model A -> Model C)")
log("=" * 70)

def calculate_nri_idi(y_true, prob_old, prob_new, n_bins=10):
    """Calculate NRI and IDI"""
    # IDI
    idi = (prob_new[y_true == 1].mean() - prob_old[y_true == 1].mean()) - \
          (prob_new[y_true == 0].mean() - prob_old[y_true == 0].mean())
    
    # NRI
    event_up = 0
    event_down = 0
    nonevent_up = 0
    nonevent_down = 0
    n_event = (y_true == 1).sum()
    n_nonevent = (y_true == 0).sum()
    
    diff = prob_new - prob_old
    
    for i in range(len(y_true)):
        if y_true.iloc[i] == 1:
            if diff.iloc[i] > 0:
                event_up += 1
            elif diff.iloc[i] < 0:
                event_down += 1
        else:
            if diff.iloc[i] > 0:
                nonevent_up += 1
            elif diff.iloc[i] < 0:
                nonevent_down += 1
    
    nri_event = (event_up - event_down) / n_event if n_event > 0 else 0
    nri_nonevent = (nonevent_down - nonevent_up) / n_nonevent if n_nonevent > 0 else 0
    nri = nri_event + nri_nonevent
    
    return nri, nri_event, nri_nonevent, idi

nri, nri_evt, nri_nevt, idi = calculate_nri_idi(y, pred_A, pred_C)
log(f"  NRI (A->C): {nri:.4f} (event: {nri_evt:.4f}, non-event: {nri_nevt:.4f})")
log(f"  IDI (A->C): {idi:.4f}")

nri_BC, nri_evt_BC, nri_nevt_BC, idi_BC = calculate_nri_idi(y, pred_B, pred_C)
log(f"  NRI (B->C): {nri_BC:.4f} (event: {nri_evt_BC:.4f}, non-event: {nri_nevt_BC:.4f})")
log(f"  IDI (B->C): {idi_BC:.4f}")

# ============================================================
# 3. LASSO Feature Selection
# ============================================================
log("\n" + "=" * 70)
log("3. LASSO Feature Selection")
log("=" * 70)

from sklearn.linear_model import LassoCV as LassoCV_sk
from sklearn.preprocessing import StandardScaler

all_features = ["age", "male", "htn", "dm", "hf", "ckd",
                "sofa_total", "sofa_nonrenal",
                "sofa_resp", "sofa_coag", "sofa_liver", "sofa_cv", "sofa_cns", "sofa_renal",
                "aki_stage2", "aki_stage3", "ckd_stage2", "ckd_stage3"]

X_lasso = aki[all_features].copy()
scaler = StandardScaler()
X_lasso_scaled = pd.DataFrame(scaler.fit_transform(X_lasso), columns=all_features)

# Use L1-penalized logistic regression
from sklearn.linear_model import LogisticRegressionCV
lasso_lr = LogisticRegressionCV(
    Cs=np.logspace(-3, 1, 50),
    penalty='l1',
    solver='saga',
    cv=5,
    max_iter=5000,
    random_state=42,
    scoring='roc_auc'
)
lasso_lr.fit(X_lasso_scaled, y)

log(f"  Best C: {lasso_lr.C_[0]:.4f}")
log(f"  LASSO AUC: {roc_auc_score(y, lasso_lr.predict_proba(X_lasso_scaled)[:, 1]):.4f}")

log(f"\n  LASSO selected features (non-zero coefficients):")
lasso_coefs = pd.DataFrame({
    "Feature": all_features,
    "Coefficient": lasso_lr.coef_[0],
    "OR (per SD)": np.exp(lasso_lr.coef_[0])
})
selected = lasso_coefs[lasso_coefs["Coefficient"] != 0].sort_values("OR (per SD)", ascending=False)
for _, row in selected.iterrows():
    log(f"    {row['Feature']}: coef={row['Coefficient']:.4f}, OR(SD)={row['OR (per SD)']:.3f}")

# Check if CKD interaction is selected
ckd_inter_selected = any(lasso_coefs[lasso_coefs["Feature"].isin(["ckd_stage3", "ckd_stage2"])]["Coefficient"] != 0)
log(f"\n  CKD x Stage interaction selected by LASSO: {'Yes' if ckd_inter_selected else 'No'}")

# ============================================================
# 4. Random Forest Comparison
# ============================================================
log("\n" + "=" * 70)
log("4. Random Forest Comparison")
log("=" * 70)

rf_features_A = features_A
rf_features_C = features_C

rf_A = RandomForestClassifier(n_estimators=200, max_depth=8, random_state=42, n_jobs=-1)
rf_A.fit(aki[rf_features_A], y)
pred_rf_A = rf_A.predict_proba(aki[rf_features_A])[:, 1]
auc_rf_A = roc_auc_score(y, pred_rf_A)

rf_C = RandomForestClassifier(n_estimators=200, max_depth=8, random_state=42, n_jobs=-1)
rf_C.fit(aki[rf_features_C], y)
pred_rf_C = rf_C.predict_proba(aki[rf_features_C])[:, 1]
auc_rf_C = roc_auc_score(y, pred_rf_C)

log(f"  RF Model A (base):           AUC = {auc_rf_A:.4f}")
log(f"  RF Model C (with CKD inter): AUC = {auc_rf_C:.4f}")

log(f"\n  Feature importance (RF Model C):")
importance = pd.DataFrame({
    "Feature": rf_features_C,
    "Importance": rf_C.feature_importances_
}).sort_values("Importance", ascending=False)
for _, row in importance.iterrows():
    log(f"    {row['Feature']}: {row['Importance']:.4f}")

# Cross-validation AUC
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
cv_A = cross_val_score(LogisticRegression(max_iter=1000), aki[features_A], y, cv=cv, scoring='roc_auc')
cv_C = cross_val_score(LogisticRegression(max_iter=1000), aki[features_C], y, cv=cv, scoring='roc_auc')
log(f"\n  5-fold CV AUC:")
log(f"    Logistic A (base):     {cv_A.mean():.4f} +/- {cv_A.std():.4f}")
log(f"    Logistic C (+CKD int): {cv_C.mean():.4f} +/- {cv_C.std():.4f}")

# ============================================================
# 5. Nomogram Generation (Model C)
# ============================================================
log("\n" + "=" * 70)
log("5. Nomogram Generation")
log("=" * 70)

fig, ax = plt.subplots(figsize=(14, 8))
ax.set_xlim(0, 100)
ax.set_ylim(0, 12)
ax.axis('off')
ax.set_title('Nomogram for Predicting ICU Mortality in AKI Patients\n(Model C: with CKD x AKI Stage Interaction)',
             fontsize=14, fontweight='bold', pad=20)

# Extract coefficients from Model C
coef = m_C.params
const = coef['const']

# Define variables and their coefficients (excluding constant)
nom_vars = []
for var in features_C:
    if var in coef:
        nom_vars.append({
            'name': var,
            'coef': coef[var],
            'label': {
                'age': 'Age (years)',
                'male': 'Male sex',
                'htn': 'Hypertension',
                'dm': 'Diabetes',
                'hf': 'Heart failure',
                'sofa_total': 'SOFA score',
                'aki_stage2': 'AKI Stage 2',
                'aki_stage3': 'AKI Stage 3',
                'ckd': 'CKD',
                'ckd_stage2': 'CKD x Stage 2',
                'ckd_stage3': 'CKD x Stage 3'
            }.get(var, var),
            'type': 'continuous' if var in ['age', 'sofa_total'] else 'binary',
            'range': (aki[var].min(), aki[var].max()) if var in ['age', 'sofa_total'] else (0, 1)
        })

# Draw each variable row
y_pos = 11
scale_factor = 10  # points per unit of beta * x

for var in nom_vars:
    # Variable label
    ax.text(2, y_pos, var['label'], fontsize=10, va='center', ha='left', fontweight='bold')
    
    if var['type'] == 'continuous':
        # Draw scale
        vmin, vmax = var['range']
        n_ticks = 6
        for i in range(n_ticks + 1):
            val = vmin + (vmax - vmin) * i / n_ticks
            points = abs(var['coef'] * val) * scale_factor
            x_pos = 25 + points
            ax.plot([x_pos, x_pos], [y_pos - 0.15, y_pos + 0.15], 'k-', linewidth=0.5)
            ax.text(x_pos, y_pos - 0.5, f'{val:.0f}', fontsize=8, ha='center', va='top')
        # Axis line
        max_points = abs(var['coef'] * vmax) * scale_factor
        ax.plot([25, 25 + max_points], [y_pos, y_pos], 'k-', linewidth=0.5)
    else:
        # Binary: show 0 and 1
        for val in [0, 1]:
            points = abs(var['coef'] * val) * scale_factor
            x_pos = 25 + points
            ax.plot([x_pos, x_pos], [y_pos - 0.15, y_pos + 0.15], 'k-', linewidth=0.5)
            ax.text(x_pos, y_pos - 0.5, f'{val}', fontsize=8, ha='center', va='top')
        max_points = abs(var['coef'] * 1) * scale_factor
        ax.plot([25, 25 + max_points], [y_pos, y_pos], 'k-', linewidth=0.5)
    
    y_pos -= 1

# Total points scale
ax.text(2, y_pos, 'Total Points', fontsize=10, va='center', ha='left', fontweight='bold')
total_max = sum(abs(v['coef'] * (v['range'][1] if v['type'] == 'continuous' else 1)) * scale_factor for v in nom_vars)
for i in range(11):
    pct = i / 10
    x_pos = 25 + pct * total_max
    ax.plot([x_pos, x_pos], [y_pos - 0.15, y_pos + 0.15], 'k-', linewidth=0.5)
    ax.text(x_pos, y_pos - 0.5, f'{pct * total_max:.0f}', fontsize=7, ha='center', va='top')
ax.plot([25, 25 + total_max], [y_pos, y_pos], 'k-', linewidth=0.5)

# Probability scale
y_pos -= 1.5
ax.text(2, y_pos, 'Predicted Mortality', fontsize=10, va='center', ha='left', fontweight='bold')
for prob in [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80]:
    logit = np.log(prob / (1 - prob))
    points = (logit - const) * scale_factor
    if 25 <= 25 + points <= 25 + total_max:
        x_pos = 25 + points
        ax.plot([x_pos, x_pos], [y_pos - 0.15, y_pos + 0.15], 'k-', linewidth=0.5)
        ax.text(x_pos, y_pos - 0.5, f'{prob:.0%}', fontsize=8, ha='center', va='top')
ax.plot([25, 25 + total_max], [y_pos, y_pos], 'k-', linewidth=0.5)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "figure8_nomogram.png"), dpi=200, bbox_inches='tight')
plt.close()
log(f"  Saved: figure8_nomogram.png")

# ============================================================
# 6. ROC Curves Comparison
# ============================================================
log("\n" + "=" * 70)
log("6. ROC Curves Comparison")
log("=" * 70)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: Logistic regression models
ax = axes[0]
for model_name, pred, auc, color in [
    ("Model A (Base)", pred_A, auc_A, '#2196F3'),
    ("Model B (+CKD)", pred_B, auc_B, '#FF9800'),
    ("Model C (+CKD×Stage)", pred_C, auc_C, '#F44336'),
    ("Model D (SOFA components)", pred_D, auc_D, '#4CAF50')
]:
    fpr, tpr, _ = roc_curve(y, pred)
    ax.plot(fpr, tpr, color=color, linewidth=2, label=f'{model_name} (AUC={auc:.3f})')

ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.5)
ax.set_xlabel('1 - Specificity', fontsize=12)
ax.set_ylabel('Sensitivity', fontsize=12)
ax.set_title('A. Logistic Regression Models', fontsize=13, fontweight='bold')
ax.legend(loc='lower right', fontsize=9)
ax.set_xlim(-0.01, 1.01)
ax.set_ylim(-0.01, 1.01)
ax.grid(True, alpha=0.3)

# Panel B: Feature importance comparison
ax = axes[1]
importance_sorted = importance.sort_values('Importance', ascending=True)
colors = ['#FF5722' if 'ckd' in f.lower() else '#2196F3' for f in importance_sorted['Feature']]
ax.barh(range(len(importance_sorted)), importance_sorted['Importance'], color=colors, edgecolor='white', linewidth=0.5)
ax.set_yticks(range(len(importance_sorted)))
ax.set_yticklabels(importance_sorted['Feature'], fontsize=10)
ax.set_xlabel('Random Forest Importance', fontsize=12)
ax.set_title('B. Feature Importance (RF Model C)', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3, axis='x')

# Add legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#FF5722', label='CKD-related'),
                   Patch(facecolor='#2196F3', label='Other features')]
ax.legend(handles=legend_elements, loc='lower right', fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "figure9_roc_importance.png"), dpi=200, bbox_inches='tight')
plt.close()
log(f"  Saved: figure9_roc_importance.png")

# ============================================================
# 7. Calibration Curves
# ============================================================
log("\n" + "=" * 70)
log("7. Calibration Curves")
log("=" * 70)

fig, ax = plt.subplots(figsize=(8, 8))

from sklearn.calibration import calibration_curve
for model_name, pred, color in [
    ("Model A (Base)", pred_A, '#2196F3'),
    ("Model C (+CKD×Stage)", pred_C, '#F44336')
]:
    frac_pos, mean_pred = calibration_curve(y, pred, n_bins=10, strategy='quantile')
    brier = brier_score_loss(y, pred)
    ax.plot(mean_pred, frac_pos, 'o-', color=color, linewidth=2, markersize=8,
            label=f'{model_name} (Brier={brier:.4f})')

ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.5, label='Perfect calibration')
ax.set_xlabel('Predicted Probability', fontsize=12)
ax.set_ylabel('Observed Mortality Rate', fontsize=12)
ax.set_title('Calibration Curves', fontsize=14, fontweight='bold')
ax.legend(loc='upper left', fontsize=10)
ax.set_xlim(-0.01, 1.01)
ax.set_ylim(-0.01, 1.01)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "figure10_calibration.png"), dpi=200, bbox_inches='tight')
plt.close()
log(f"  Saved: figure10_calibration.png")

# ============================================================
# 8. Save prediction model results
# ============================================================
log("\n" + "=" * 70)
log("8. Save Results")
log("=" * 70)

model_results = pd.DataFrame({
    "Model": ["A: Base", "B: +CKD", "C: +CKD×Stage", "D: SOFA components"],
    "Features": [
        "age,sex,comorbidities,SOFA,AKI stage",
        "+ CKD main effect",
        "+ CKD×Stage2, CKD×Stage3 interaction",
        "SOFA 6 components + CKD + interaction"
    ],
    "N_features": [len(features_A), len(features_B), len(features_C), len(features_D)],
    "AUC": [auc_A, auc_B, auc_C, auc_D],
    "Pseudo_R2": [m_A.prsquared, m_B.prsquared, m_C.prsquared, m_D.prsquared],
    "Brier": [brier_score_loss(y, pred_A), brier_score_loss(y, pred_B),
              brier_score_loss(y, pred_C), brier_score_loss(y, pred_D)],
    "NRI_vs_A": [0, np.nan, nri, np.nan],
    "IDI_vs_A": [0, np.nan, idi, np.nan],
    "LR_test_P_vs_prev": [np.nan, p_AB, p_BC, np.nan]
})
model_results.to_csv(os.path.join(OUTPUT_DIR, "table_prediction_model_comparison.csv"), index=False)
log(f"  Saved: table_prediction_model_comparison.csv")

# Save LASSO results
lasso_coefs.to_csv(os.path.join(OUTPUT_DIR, "table_lasso_coefficients.csv"), index=False)
log(f"  Saved: table_lasso_coefficients.csv")

# Save Model C coefficients for nomogram
nomogram_table = pd.DataFrame({
    "Variable": features_C,
    "Coefficient": [coef.get(v, 0) for v in features_C],
    "OR": [np.exp(coef.get(v, 0)) for v in features_C],
    "OR_CI_lower": [np.exp(m_C.conf_int().loc[v, 0]) if v in m_C.conf_int().index else np.nan for v in features_C],
    "OR_CI_upper": [np.exp(m_C.conf_int().loc[v, 1]) if v in m_C.conf_int().index else np.nan for v in features_C],
    "P_value": [m_C.pvalues.get(v, np.nan) for v in features_C]
})
nomogram_table.to_csv(os.path.join(OUTPUT_DIR, "table_nomogram_coefficients.csv"), index=False)
log(f"  Saved: table_nomogram_coefficients.csv")

# ============================================================
# SUMMARY
# ============================================================
log("\n" + "=" * 70)
log("SUMMARY: Prediction Modeling")
log("=" * 70)
log(f"""
KEY FINDINGS:

1. Model Performance:
   - Model A (base):              AUC = {auc_A:.4f}
   - Model B (+CKD):              AUC = {auc_B:.4f}
   - Model C (+CKD x Stage):      AUC = {auc_C:.4f}
   - Model D (SOFA components):    AUC = {auc_D:.4f}

2. Incremental Value of CKD:
   - Adding CKD main effect:      delta AUC = {auc_B - auc_A:.4f} (P = {p_AB:.4f})
   - Adding CKD x Stage:          delta AUC = {auc_C - auc_B:.4f} (P = {p_BC:.8f})

3. NRI & IDI (A -> C):
   - NRI = {nri:.4f}
   - IDI = {idi:.4f}

4. LASSO selected CKD x Stage3: {ckd_inter_selected}

5. 5-fold CV AUC:
   - Base model:     {cv_A.mean():.4f} +/- {cv_C.std():.4f}
   - With CKD inter: {cv_C.mean():.4f} +/- {cv_C.std():.4f}
""")

log("Done! Prediction modeling complete.")
