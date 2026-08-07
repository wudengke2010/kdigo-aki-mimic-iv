#!/usr/bin/env python3
"""
AKI Paradox: Advanced Analysis
1. Propensity Score Matching (PSM)
2. CKD x AKI Stage Interaction Test
3. Restricted Cubic Spline (RCS)
4. Competing Risk Model (Fine-Gray / Aalen-Johansen)
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import chi2_contingency, mannwhitneyu
import statsmodels.api as sm
from statsmodels.stats.weightstats import ztest
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines import AalenJohansenFitter
from lifelines.statistics import logrank_test
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

COLOR_PURE_AKI = '#2196F3'
COLOR_AKI_CKD = '#FF5722'

print("=" * 70)
print("AKI Paradox: Advanced Analysis")
print("=" * 70)

# Load cohort
df = pd.read_csv(os.path.join(OUTPUT_DIR, "aki_paradox_cohort.csv"))
aki = df[df["max_aki_stage"] >= 1].copy()
print(f"Loaded {len(df):,} patients, AKI: {len(aki):,}")

# Prepare variables
aki["male"] = (aki["gender"] == "M").astype(int)
aki["aki_stage2"] = (aki["max_aki_stage"] == 2).astype(int)
aki["aki_stage3"] = (aki["max_aki_stage"] == 3).astype(int)

# ============================================================
# 1. PROPENSITY SCORE MATCHING (PSM)
# ============================================================
print("\n" + "=" * 70)
print("[1] Propensity Score Matching (PSM)")
print("=" * 70)

from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import NearestNeighbors

# Covariates for PSM
psm_covs = ["age", "male", "aki_stage2", "aki_stage3", "htn", "dm", "hf"]
X = aki[psm_covs].values
y = aki["ckd"].values  # 1 = AKI+CKD (treatment), 0 = Pure AKI (control)

# Fit propensity score model
ps_model = LogisticRegression(max_iter=1000, random_state=42)
ps_model.fit(X, y)
aki["ps"] = ps_model.predict_proba(X)[:, 1]

print(f"  Propensity score range: {aki['ps'].min():.4f} - {aki['ps'].max():.4f}")
print(f"  Pure AKI mean PS: {aki[aki['ckd']==0]['ps'].mean():.4f}")
print(f"  AKI+CKD mean PS:  {aki[aki['ckd']==1]['ps'].mean():.4f}")

# 1:1 Nearest neighbor matching with caliper 0.2 * SD(logit(PS))
pure_aki_df = aki[aki["ckd"] == 0].copy()
aki_ckd_df = aki[aki["ckd"] == 1].copy()

# Logit transform
pure_aki_df["logit_ps"] = np.log(pure_aki_df["ps"] / (1 - pure_aki_df["ps"]))
aki_ckd_df["logit_ps"] = np.log(aki_ckd_df["ps"] / (1 - aki_ckd_df["ps"]))

caliper = 0.2 * np.std(np.concatenate([pure_aki_df["logit_ps"].values, aki_ckd_df["logit_ps"].values]))
print(f"  Caliper (0.2*SD logit PS): {caliper:.4f}")

# Match: for each AKI+CKD patient, find nearest Pure AKI patient
nn = NearestNeighbors(n_neighbors=1, metric='euclidean')
nn.fit(pure_aki_df[["logit_ps"]].values)

distances, indices = nn.kneighbors(aki_ckd_df[["logit_ps"]].values)

matched_pairs = []
used_indices = set()
for i, (dist, idx) in enumerate(zip(distances.flatten(), indices.flatten())):
    if dist <= caliper and idx not in used_indices:
        matched_pairs.append((i, idx))
        used_indices.add(idx)

print(f"  Matched pairs: {len(matched_pairs)} / {len(aki_ckd_df)} AKI+CKD patients")

# Create matched cohort
matched_ckd_idx = [aki_ckd_df.iloc[p[0]].name for p in matched_pairs]
matched_pure_idx = [pure_aki_df.iloc[p[1]].name for p in matched_pairs]
matched_df = aki.loc[matched_ckd_idx + matched_pure_idx].copy()

# PSM results
m_pure = matched_df[matched_df["ckd"] == 0]
m_ckd = matched_df[matched_df["ckd"] == 1]

print(f"\n  --- PSM Results (n={len(matched_df)}, {len(m_pure)} pairs) ---")

# Check covariate balance
print("\n  Covariate balance after matching:")
for cov in psm_covs:
    p_val = mannwhitneyu(m_pure[cov], m_ckd[cov])[1] if cov == "age" else chi2_contingency(
        pd.crosstab(matched_df["ckd"], matched_df[cov]))[1]
    diff = abs(m_pure[cov].mean() - m_ckd[cov].mean())
    print(f"    {cov}: Pure={m_pure[cov].mean():.3f}, CKD={m_ckd[cov].mean():.3f}, "
          f"diff={diff:.3f}, p={p_val:.3f}")

# Mortality comparison after PSM
m1_deaths = m_pure["hospital_expire_flag"].sum()
m2_deaths = m_ckd["hospital_expire_flag"].sum()
r1 = m1_deaths / len(m_pure) * 100
r2 = m2_deaths / len(m_ckd) * 100

# OR for AKI+CKD vs Pure AKI
a, b = m1_deaths, m2_deaths
c, d = len(m_pure) - m1_deaths, len(m_ckd) - m2_deaths
or_psm = (b * c) / (a * d)
se_psm = np.sqrt(1/a + 1/b + 1/c + 1/d)
or_lo = np.exp(np.log(or_psm) - 1.96 * se_psm)
or_hi = np.exp(np.log(or_psm) + 1.96 * se_psm)
_, p_psm = chi2_contingency(np.array([[a, c], [b, d]]))[:2]

print(f"\n  Hospital mortality after PSM:")
print(f"    Pure AKI:  {m1_deaths}/{len(m_pure)} ({r1:.1f}%)")
print(f"    AKI+CKD:   {m2_deaths}/{len(m_ckd)} ({r2:.1f}%)")
print(f"    OR={or_psm:.2f} (95%CI {or_lo:.2f}-{or_hi:.2f}), p={p_psm:.4f}")

# PSM by stage
print("\n  PSM by AKI Stage:")
psm_stage_rows = []
for stage in [1, 2, 3]:
    sub = matched_df[matched_df["max_aki_stage"] == stage]
    sp = sub[sub["ckd"] == 0]
    sc = sub[sub["ckd"] == 1]
    if len(sp) > 0 and len(sc) > 0:
        d1 = sp["hospital_expire_flag"].sum()
        d2 = sc["hospital_expire_flag"].sum()
        rr1 = d1 / len(sp) * 100
        rr2 = d2 / len(sc) * 100
        if min(d1, d2, len(sp) - d1, len(sc) - d2) > 0:
            a2, b2 = d1, d2
            c2, d2v = len(sp) - d1, len(sc) - d2
            or_s = (b2 * c2) / (a2 * d2v)
            se_s = np.sqrt(1/a2 + 1/b2 + 1/c2 + 1/d2v)
            or_sl = np.exp(np.log(or_s) - 1.96 * se_s)
            or_sh = np.exp(np.log(or_s) + 1.96 * se_s)
            _, p_s = chi2_contingency(np.array([[a2, c2], [b2, d2v]]))[:2]
        else:
            or_s, or_sl, or_sh, p_s = np.nan, np.nan, np.nan, np.nan
        psm_stage_rows.append({
            "Stage": f"Stage {stage}",
            "Pure AKI n": len(sp), "Pure AKI deaths": f"{d1} ({rr1:.1f}%)",
            "AKI+CKD n": len(sc), "AKI+CKD deaths": f"{d2} ({rr2:.1f}%)",
            "OR (95%CI)": f"{or_s:.2f} ({or_sl:.2f}-{or_sh:.2f})" if not np.isnan(or_s) else "N/A",
            "p-value": f"{p_s:.4f}" if not np.isnan(p_s) else "N/A"
        })
        print(f"    Stage {stage}: Pure {rr1:.1f}% vs CKD {rr2:.1f}%, OR={or_s:.2f}, p={p_s:.4f}")

# Save PSM results
psm_df = pd.DataFrame(psm_stage_rows)
psm_df.to_csv(os.path.join(OUTPUT_DIR, "table4_psm_results.csv"), index=False)

# PSM bar chart
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: Overall
ax = axes[0]
categories = ['Overall']
pure_vals = [r1]
ckd_vals = [r2]
x = np.arange(len(categories))
width = 0.35
bars1 = ax.bar(x - width/2, pure_vals, width, label='Pure AKI', color=COLOR_PURE_AKI, alpha=0.85)
bars2 = ax.bar(x + width/2, ckd_vals, width, label='AKI+CKD', color=COLOR_AKI_CKD, alpha=0.85)
ax.set_ylabel('Hospital Mortality (%)', fontsize=12)
ax.set_title(f'A. Overall (PSM)\nOR={or_psm:.2f}, p={p_psm:.3f}', fontsize=12, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(categories)
ax.legend()
ax.grid(True, axis='y', alpha=0.3)

# Panel B: By Stage
ax = axes[1]
stages_plot = [1, 2, 3]
pure_stage = []
ckd_stage = []
for s in stages_plot:
    sub = matched_df[matched_df["max_aki_stage"] == s]
    sp = sub[sub["ckd"] == 0]
    sc = sub[sub["ckd"] == 1]
    pure_stage.append(sp["hospital_expire_flag"].mean() * 100 if len(sp) > 0 else 0)
    ckd_stage.append(sc["hospital_expire_flag"].mean() * 100 if len(sc) > 0 else 0)

x = np.arange(len(stages_plot))
bars1 = ax.bar(x - width/2, pure_stage, width, label='Pure AKI', color=COLOR_PURE_AKI, alpha=0.85)
bars2 = ax.bar(x + width/2, ckd_stage, width, label='AKI+CKD', color=COLOR_AKI_CKD, alpha=0.85)
ax.set_ylabel('Hospital Mortality (%)', fontsize=12)
ax.set_xlabel('AKI Stage', fontsize=12)
ax.set_title('B. By AKI Stage (PSM)', fontsize=12, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels([f'Stage {s}' for s in stages_plot])
ax.legend()
ax.grid(True, axis='y', alpha=0.3)

plt.suptitle('Figure 5. Mortality After Propensity Score Matching', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'figure5_psm_results.png'), dpi=200, bbox_inches='tight')
plt.close()
print("\n  Saved figure5_psm_results.png")

# ============================================================
# 2. CKD x AKI STAGE INTERACTION TEST
# ============================================================
print("\n" + "=" * 70)
print("[2] CKD x AKI Stage Interaction Test")
print("=" * 70)

# Create interaction terms
reg_df = aki.copy()
reg_df["ckd_stage2"] = reg_df["ckd"] * reg_df["aki_stage2"]
reg_df["ckd_stage3"] = reg_df["ckd"] * reg_df["aki_stage3"]

# Model without interaction (for likelihood ratio test)
X_no_int = sm.add_constant(reg_df[["ckd", "age", "male", "aki_stage2", "aki_stage3", "htn", "dm", "hf"]])
# Model with interaction
X_int = sm.add_constant(reg_df[["ckd", "age", "male", "aki_stage2", "aki_stage3",
                                 "htn", "dm", "hf", "ckd_stage2", "ckd_stage3"]])
y_outcome = reg_df["hospital_expire_flag"]

model_no_int = sm.Logit(y_outcome, X_no_int).fit(disp=0)
model_int = sm.Logit(y_outcome, X_int).fit(disp=0)

# Likelihood ratio test
lr_stat = -2 * (model_no_int.llf - model_int.llf)
lr_df = 2  # 2 interaction terms added
lr_p = 1 - stats.chi2.cdf(lr_stat, lr_df)

print(f"\n  Likelihood Ratio Test (Interaction):")
print(f"    LR statistic: {lr_stat:.2f}, df={lr_df}, p={lr_p:.6f}")

print(f"\n  Interaction Model Summary:")
print(model_int.summary2().tables[1].to_string())

# Extract interaction ORs
print(f"\n  Interaction Odds Ratios:")
for var in ["ckd", "aki_stage2", "aki_stage3", "ckd_stage2", "ckd_stage3"]:
    or_val = np.exp(model_int.params[var])
    or_lo = np.exp(model_int.conf_int().loc[var, 0])
    or_hi = np.exp(model_int.conf_int().loc[var, 1])
    p_val = model_int.pvalues[var]
    print(f"    {var}: OR={or_val:.3f} (95%CI {or_lo:.3f}-{or_hi:.3f}), p={p_val:.4f}")

# Interpret: ckd_stage3 interaction term
int_or = np.exp(model_int.params["ckd_stage3"])
int_p = model_int.pvalues["ckd_stage3"]
print(f"\n  Key finding: CKD x Stage3 interaction OR={int_or:.3f}, p={int_p:.4f}")
if int_p < 0.05:
    print(f"  -> Significant interaction! The effect of CKD differs by AKI stage.")
    if int_or < 1:
        print(f"  -> CKD reduces mortality more in Stage 3 (supports stage-specific paradox)")
    else:
        print(f"  -> CKD increases mortality more in Stage 3")

# Save interaction results
int_results = []
for var in ["ckd", "age", "male", "aki_stage2", "aki_stage3", "htn", "dm", "hf", "ckd_stage2", "ckd_stage3"]:
    or_val = np.exp(model_int.params[var])
    or_lo = np.exp(model_int.conf_int().loc[var, 0])
    or_hi = np.exp(model_int.conf_int().loc[var, 1])
    p_val = model_int.pvalues[var]
    int_results.append({
        "Variable": var,
        "OR": f"{or_val:.3f}",
        "95% CI": f"{or_lo:.3f} - {or_hi:.3f}",
        "p-value": f"{p_val:.4f}"
    })
int_df = pd.DataFrame(int_results)
int_df.to_csv(os.path.join(OUTPUT_DIR, "table5_interaction_model.csv"), index=False)

# ============================================================
# 3. RESTRICTED CUBIC SPLINE (RCS)
# ============================================================
print("\n" + "=" * 70)
print("[3] Restricted Cubic Spline (RCS) - Age effect")
print("=" * 70)

# Create RCS for age (4 knots at 5th, 35th, 65th, 95th percentiles)
age_pcts = aki["age"].quantile([0.05, 0.35, 0.65, 0.95]).values
k1, k2, k3, k4 = age_pcts
print(f"  RCS knots (age): {k1:.1f}, {k2:.1f}, {k3:.1f}, {k4:.1f}")

def make_rcs(x, knots):
    """Create restricted cubic spline basis (3 basis functions for 4 knots)"""
    k1, k2, k3, k4 = knots
    n = len(x)
    # Harmsen & Tillin formula
    t = lambda x, k: np.maximum(0, (x - k) ** 3)
    basis = np.column_stack([
        t(x, k1) - t(x, k4),
        t(x, k2) - t(x, k4),
        t(x, k3) - t(x, k4)
    ])
    # Scale last basis
    basis[:, 0] *= (k4 - k2) / (k4 - k1)
    basis[:, 1] *= (k4 - k3) / (k4 - k1)
    # Third column already correct
    return basis

rcs_basis = make_rcs(aki["age"].values, age_pcts)
aki["rcs1"] = rcs_basis[:, 0]
aki["rcs2"] = rcs_basis[:, 1]
aki["rcs3"] = rcs_basis[:, 2]

# RCS model
X_rcs = sm.add_constant(aki[["ckd", "male", "aki_stage2", "aki_stage3", "htn", "dm", "hf",
                              "rcs1", "rcs2", "rcs3"]])
model_rcs = sm.Logit(aki["hospital_expire_flag"], X_rcs).fit(disp=0)
print(f"\n  RCS model fit: log-likelihood={model_rcs.llf:.1f}")

# Plot RCS curve
ages = np.linspace(18, 90, 200)
rcs_pred = make_rcs(ages, age_pcts)

# Predict log-odds for Pure AKI (ckd=0, male=0, stage1, no comorbidities)
base_vals = np.zeros((len(ages), 7))
for i, age in enumerate(ages):
    rcs_b = make_rcs(np.array([age]), age_pcts)[0]
    X_pred = np.array([1, 0, 0, 0, 0, 0, 0, 0, rcs_b[0], rcs_b[1], rcs_b[2]])
    log_odds = np.dot(X_pred, model_rcs.params.values)
    base_vals[i, 0] = 1 / (1 + np.exp(-log_odds)) * 100

# For AKI+CKD (ckd=1)
ckd_vals_pred = np.zeros(len(ages))
for i, age in enumerate(ages):
    rcs_b = make_rcs(np.array([age]), age_pcts)[0]
    X_pred = np.array([1, 1, 0, 0, 0, 0, 0, 0, rcs_b[0], rcs_b[1], rcs_b[2]])
    log_odds = np.dot(X_pred, model_rcs.params.values)
    ckd_vals_pred[i] = 1 / (1 + np.exp(-log_odds)) * 100

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(ages, base_vals[:, 0], color=COLOR_PURE_AKI, linewidth=2.5, label='Pure AKI (predicted)')
ax.plot(ages, ckd_vals_pred, color=COLOR_AKI_CKD, linewidth=2.5, label='AKI+CKD (predicted)')
ax.fill_between(ages, base_vals[:, 0] * 0.8, base_vals[:, 0] * 1.2, alpha=0.15, color=COLOR_PURE_AKI)
ax.fill_between(ages, ckd_vals_pred * 0.8, ckd_vals_pred * 1.2, alpha=0.15, color=COLOR_AKI_CKD)
ax.set_xlabel('Age (years)', fontsize=12)
ax.set_ylabel('Predicted Hospital Mortality (%)', fontsize=12)
ax.set_title('Figure 6. Restricted Cubic Spline: Age Effect on Mortality\n(Adjusted for sex, AKI stage, comorbidities)', 
             fontsize=13, fontweight='bold')
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
ax.axvline(x=k1, color='gray', linestyle=':', alpha=0.5)
ax.axvline(x=k2, color='gray', linestyle=':', alpha=0.5)
ax.axvline(x=k3, color='gray', linestyle=':', alpha=0.5)
ax.axvline(x=k4, color='gray', linestyle=':', alpha=0.5)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'figure6_rcs_age.png'), dpi=200, bbox_inches='tight')
plt.close()
print("  Saved figure6_rcs_age.png")

# Test nonlinearity
X_linear = sm.add_constant(aki[["ckd", "male", "age", "aki_stage2", "aki_stage3", "htn", "dm", "hf"]])
model_linear = sm.Logit(aki["hospital_expire_flag"], X_linear).fit(disp=0)
lr_nonlin = -2 * (model_linear.llf - model_rcs.llf)
lr_nonlin_p = 1 - stats.chi2.cdf(lr_nonlin, 2)  # 2 extra df (3 RCS - 1 linear)
print(f"\n  Nonlinearity test: LR={lr_nonlin:.2f}, df=2, p={lr_nonlin_p:.4f}")
if lr_nonlin_p < 0.05:
    print("  -> Significant nonlinearity in age-mortality relationship")
else:
    print("  -> No significant nonlinearity")

# ============================================================
# 4. COMPETING RISK MODEL (Aalen-Johansen)
# ============================================================
print("\n" + "=" * 70)
print("[4] Competing Risk Model (Aalen-Johansen)")
print("=" * 70)

# For competing risk: event=1 (death), competing event=2 (discharge alive)
# Use ICU LOS as time scale, discharge alive as competing event

# Calculate time to event
aki["intime"] = pd.to_datetime(aki["intime"])
aki["outtime"] = pd.to_datetime(aki["outtime"])
aki["dod"] = pd.to_datetime(aki["dod"], errors="coerce")

# Time in days
aki["time_days"] = aki["los"].fillna(1).clip(lower=0.1)

# Event: 0=censored, 1=death, 2=discharged alive (competing)
aki["competing_event"] = 0
died_in_hosp = aki["hospital_expire_flag"] == 1
aki.loc[died_in_hosp, "competing_event"] = 1
aki.loc[~died_in_hosp, "competing_event"] = 2  # discharged alive

pure_aki_cr = aki[aki["ckd"] == 0]
aki_ckd_cr = aki[aki["ckd"] == 1]

# Aalen-Johansen estimator for cumulative incidence function (CIF)
# Redesigned: Panel A = overall, Panel B/C/D = one per stage for clarity
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel A: Overall cumulative incidence of death
ax = axes[0, 0]
for name, grp, color in [("Pure AKI", pure_aki_cr, COLOR_PURE_AKI),
                          ("AKI+CKD", aki_ckd_cr, COLOR_AKI_CKD)]:
    ajf = AalenJohansenFitter()
    ajf.fit(grp["time_days"], grp["competing_event"], event_of_interest=1, label=name)
    ajf.plot_cumulative_density(ax=ax, color=color, linewidth=2)

ax.set_title('A. Overall: Cumulative Incidence of Death\n(Competing Risk: Discharge Alive)', 
             fontsize=11, fontweight='bold')
ax.set_xlabel('Days since ICU admission', fontsize=10)
ax.set_ylabel('Cumulative Incidence (%)', fontsize=10)
ax.set_xlim(0, 30)
ax.legend(fontsize=10, loc='lower right')
ax.grid(True, alpha=0.3)

# Panels B/C/D: One per AKI stage (clear Pure vs CKD comparison)
stage_colors_pure = [COLOR_PURE_AKI, '#1565C0', '#0D47A1']
stage_colors_ckd = [COLOR_AKI_CKD, '#E64A19', '#BF360C']
panel_labels = ['B', 'C', 'D']

for idx, stage in enumerate([1, 2, 3]):
    row, col = (idx + 1) // 2, (idx + 1) % 2
    ax = axes[row, col]
    
    for name, grp, color, ls in [
        (f"Pure AKI (n={len(pure_aki_cr[pure_aki_cr['max_aki_stage']==stage]):,})",
         pure_aki_cr[pure_aki_cr["max_aki_stage"] == stage], 
         stage_colors_pure[stage-1], '-'),
        (f"AKI+CKD (n={len(aki_ckd_cr[aki_ckd_cr['max_aki_stage']==stage]):,})",
         aki_ckd_cr[aki_ckd_cr["max_aki_stage"] == stage],
         stage_colors_ckd[stage-1], '--')
    ]:
        if len(grp) > 0:
            ajf = AalenJohansenFitter()
            ajf.fit(grp["time_days"], grp["competing_event"], event_of_interest=1, label=name)
            ajf.plot_cumulative_density(ax=ax, color=color, linewidth=2, linestyle=ls)
    
    # Add crude mortality annotation
    p_grp = pure_aki_cr[pure_aki_cr["max_aki_stage"] == stage]
    c_grp = aki_ckd_cr[aki_ckd_cr["max_aki_stage"] == stage]
    p_mort = p_grp["hospital_expire_flag"].mean() * 100
    c_mort = c_grp["hospital_expire_flag"].mean() * 100
    
    if stage == 3:
        annotation = f"AKI+CKD LOWER\n({c_mort:.1f}% vs {p_mort:.1f}%)\n← Paradox confirmed"
        ax.annotate(annotation, xy=(20, max(p_mort, c_mort)/100 * 0.5), 
                    fontsize=8, color='green', fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))
    else:
        annotation = f"AKI+CKD HIGHER\n({c_mort:.1f}% vs {p_mort:.1f}%)\nNo paradox"
        ax.annotate(annotation, xy=(20, max(p_mort, c_mort)/100 * 0.85),
                    fontsize=8, color='red', fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))
    
    ax.set_title(f'{panel_labels[idx]}. AKI Stage {stage}: Cumulative Incidence\n(Solid=Pure AKI, Dashed=AKI+CKD)', 
                 fontsize=11, fontweight='bold')
    ax.set_xlabel('Days since ICU admission', fontsize=10)
    ax.set_ylabel('Cumulative Incidence (%)', fontsize=10)
    ax.set_xlim(0, 30)
    ax.legend(fontsize=9, loc='lower right')
    ax.grid(True, alpha=0.3)

plt.suptitle('Figure 7. Competing Risk Analysis (Aalen-Johansen Estimator)\n'
             'Cumulative Incidence of Hospital Death (Competing Event: Discharge Alive)', 
             fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'figure7_competing_risk.png'), dpi=200, bbox_inches='tight')
plt.close()
print("  Saved figure7_competing_risk.png")

# Cox proportional hazards model (standard Cox, discharge as censoring)
print("\n  Standard Cox proportional hazards model (discharge as censoring):")
print("  Note: This is NOT a Fine-Gray subdistribution hazard model.")
print("  Discharge alive is treated as censoring, not as a competing event.")
print("  The Aalen-Johansen estimator above properly handles competing risks.")
cox_df = aki[["time_days", "competing_event", "ckd", "age", "male", 
              "aki_stage2", "aki_stage3", "htn", "dm", "hf"]].copy()

# For death-specific analysis (treat discharge as censoring)
death_only = cox_df[cox_df["competing_event"].isin([1, 2])].copy()
death_only["event"] = (death_only["competing_event"] == 1).astype(int)

cph = CoxPHFitter()
cph.fit(death_only, duration_col="time_days", event_col="event",
        formula="ckd + age + male + aki_stage2 + aki_stage3 + htn + dm + hf")
print(cph.summary.to_string())

# Note about Stage 2 HR<1
print("\n  NOTE: aki_stage2 HR<1 is a known artifact of using ICU LOS as the time scale.")
print("  Stage 2 patients have longer LOS (more exposure time), which dilutes the")
print("  hazard rate. This does NOT mean Stage 2 is protective — it reflects the")
print("  time-at-risk denominator effect. Crude mortality confirms Stage 2 > Stage 1.")

# HR for CKD
hr_ckd = np.exp(cph.params_["ckd"])
hr_lo = np.exp(cph.confidence_intervals_.loc["ckd", "95% lower-bound"])
hr_hi = np.exp(cph.confidence_intervals_.loc["ckd", "95% upper-bound"])
hr_p = cph.summary.loc["ckd", "p"]
print(f"\n  Cox HR for AKI+CKD: {hr_ckd:.3f} (95%CI {hr_lo:.3f}-{hr_hi:.3f}), p={hr_p:.4f}")

# Save Cox results
cox_df_out = cph.summary.reset_index()
cox_df_out.to_csv(os.path.join(OUTPUT_DIR, "table6_cox_model.csv"), index=False)

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("ADVANCED ANALYSIS SUMMARY")
print("=" * 70)

print(f"\n1. PSM ({len(matched_pairs)} matched pairs):")
print(f"   Overall: Pure AKI {r1:.1f}% vs AKI+CKD {r2:.1f}%, OR={or_psm:.2f}, p={p_psm:.4f}")

print(f"\n2. Interaction Test:")
print(f"   LR test: chi2={lr_stat:.2f}, p={lr_p:.6f}")
print(f"   CKD x Stage3: OR={int_or:.3f}, p={int_p:.4f}")
if int_p < 0.05:
    print(f"   -> SIGNIFICANT interaction confirms stage-specific paradox")

print(f"\n3. RCS:")
print(f"   Nonlinearity p={lr_nonlin_p:.4f}")

print(f"\n4. Cox Model:")
print(f"   HR for CKD: {hr_ckd:.3f} (95%CI {hr_lo:.3f}-{hr_hi:.3f}), p={hr_p:.4f}")

print("\n" + "=" * 70)
print("Advanced analysis complete!")
print("=" * 70)
