"""
F7 修复 — 真外部验证 (TRIPOD type 2/3)
======================================================
锁定 MIMIC 模型在 eICU 验证集上的真实外部表现：

  MIMIC 训练集 (mimic_analysis_v2.csv, n=60,506) 拟合 Model C
    CoxPH: 30d in-hospital mortality ~ kdigo_1 + kdigo_2 + kdigo_3 + age
           + sex_male + sofa_nonrenal + vent_24h + ckd_history_prior
           + diabetes + hypertension + heart_failure + copd + liver_disease
    → 提取 β_MIMIC (12 个 log-HR + 1 截距项 lifelines 内部参数化)

  eICU 验证集 (eicu_analysis_v2.csv, n=113,466) 完全不重 fit:
    a) Discrimination
       - C-index (Harrell's) on eICU with MIMIC's LP
       - Time-dependent AUC at 7d / 14d / 30d
       - 对比 MIMIC train C-index (bootstrap CI)

    b) Calibration
       - Calibration intercept (in-the-large): eICU-only refit of LP
       - Calibration slope: eICU-only refit of LP, slope=1 means perfect
       - Predicted vs actual mortality by decile of LP
       - Calibration plot

    c) Decision Curve Analysis
       - Net benefit at thresholds 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35

    d) Cross-cohort parameter consistency
       - HR_kdigo_3 from MIMIC Model C
       - HR_kdigo_3 from eICU-only Model C refit (univariate check)
       - Ratio (consistency test, ideally 0.5-2.0)

    e) KM by LP tertile in eICU

输出: v2_outputs/external_validation_results.json + v2_outputs/external_validation_report.md
"""
import json
import os
import warnings
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index
from sklearn.metrics import roc_auc_score, brier_score_loss
from scipy import stats
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = ['DejaVu Sans']

WORK = r"C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
os.makedirs(os.path.join(WORK, "v2_outputs"), exist_ok=True)

LANDMARK = 24          # hours
TRUNC_H  = 24 * 30     # 30d

# Model C covariates (must match between MIMIC fit and eICU apply)
COVARS = ['kdigo_1','kdigo_2','kdigo_3',
          'age','sex_male','sofa_nonrenal','vent_24h','ckd_history_prior',
          'diabetes','hypertension','heart_failure','copd','liver_disease']

print("=" * 80)
print("F7 EXTERNAL VALIDATION — locked MIMIC Model C on eICU")
print("=" * 80)

# =============================================================================
# 1. Load MIMIC analysis set + fit Model C
# =============================================================================
print("\n[1] Load MIMIC analysis set + fit Model C")
mimic = pd.read_csv(os.path.join(WORK, "mimic_analysis_v2.csv"))
print(f"  MIMIC: {len(mimic):,} stays")

# --- Sanitize event times EXACTLY as survival_mimic_v2.py (consistency) ---
# (a) floor at 0.5 h; (b) decedents with deathtime > dischtime: cap at dischtime+24h
mimic['event_time_h'] = mimic['event_time_h'].clip(lower=0.5)
_coh = pd.read_csv(os.path.join(WORK, "kdigo_cohort_v2.csv"),
                   usecols=['stay_id','intime','dischtime','deathtime'])
_coh['intime']    = pd.to_datetime(_coh['intime'])
_coh['dischtime'] = pd.to_datetime(_coh['dischtime'])
mimic = mimic.merge(_coh[['stay_id','intime','dischtime']], on='stay_id', how='left')
mimic['event_time_h'] = np.where(
    mimic['hospital_expire_flag'] == 1,
    np.minimum(mimic['event_time_h'],
               (mimic['dischtime'] - mimic['intime']).dt.total_seconds()/3600 + 24),
    mimic['event_time_h'])
mimic.loc[mimic['event_time_h'] < 0, 'event_time_h'] = 0.5

# Stage dummies
for s in [1,2,3]:
    mimic[f'kdigo_{s}'] = (mimic['kdigo_stage']==s).astype(int)

# Apply 24h landmark
mimic_surv = mimic[mimic['event_time_h'] >= LANDMARK].copy()
mimic_surv['time_from_lm'] = (mimic_surv['event_time_h'] - LANDMARK).clip(upper=TRUNC_H - LANDMARK)
print(f"  MIMIC after landmark: {len(mimic_surv):,}, events: {mimic_surv.event_ind.sum():,}")

# Fit Model C (identical specification to survival_mimic_v2.py: no penalizer)
print("  Fitting Cox Model C on MIMIC ...")
cox_C = CoxPHFitter(penalizer=0.0)
cox_C.fit(mimic_surv[COVARS + ['time_from_lm','event_ind']],
          duration_col='time_from_lm',
          event_col='event_ind',
          formula=None)
c_mimic = cox_C.concordance_index_
beta_mimic = cox_C.params_             # pd.Series index=COVARS
print(f"  MIMIC Model C: C-index={c_mimic:.4f}, n={cox_C._n_examples:,}")
for c in COVARS:
    hr = np.exp(beta_mimic[c])
    se = cox_C.summary.loc[c, 'se(coef)']
    z = beta_mimic[c] / se
    p = cox_C.summary.loc[c, 'p']
    print(f"    {c:18s} HR={hr:.3f}  95%CI ({np.exp(beta_mimic[c]-1.96*se):.3f}-{np.exp(beta_mimic[c]+1.96*se):.3f})  p={p:.2e}")

# Compute MIMIC linear predictor for bootstrap baseline
mimic_surv['lp'] = mimic_surv[COVARS].values @ beta_mimic.values

# Bootstrap C-index CI on MIMIC train (n=200)
print("  Bootstrap C-index CI on MIMIC (B=200) ...")
rng = np.random.default_rng(42)
boots_c = []
mimic_idx_arr = mimic_surv.index.values
for b in range(200):
    idx = rng.choice(mimic_idx_arr, size=len(mimic_idx_arr), replace=True)
    sub = mimic_surv.loc[idx]
    try:
        c = concordance_index(sub['time_from_lm'], -sub['lp'], sub['event_ind'])
        boots_c.append(c)
    except Exception:
        pass
c_mimic_lo, c_mimic_hi = np.percentile(boots_c, [2.5, 97.5])
print(f"  MIMIC C-index bootstrap 95%CI: {c_mimic:.4f} ({c_mimic_lo:.4f}-{c_mimic_hi:.4f})")

# =============================================================================
# 2. Load eICU analysis set + apply β_MIMIC
# =============================================================================
print("\n[2] Load eICU analysis set + apply MIMIC β")
eicu = pd.read_csv(os.path.join(WORK, "eicu_analysis_v2.csv"))
print(f"  eICU: {len(eicu):,} stays")

for s in [1,2,3]:
    eicu[f'kdigo_{s}'] = (eicu['kdigo_stage']==s).astype(int)

# Apply 24h landmark
eicu_surv = eicu[eicu['event_time_h'] >= LANDMARK].copy()
eicu_surv['time_from_lm'] = (eicu_surv['event_time_h'] - LANDMARK).clip(upper=TRUNC_H - LANDMARK)
print(f"  eICU after landmark: {len(eicu_surv):,}, events: {eicu_surv.event_ind.sum():,}")

# Apply MIMIC β — DO NOT REFIT
eicu_surv['lp'] = eicu_surv[COVARS].values @ beta_mimic.values
print(f"  LP range in eICU: {eicu_surv.lp.min():.2f} to {eicu_surv.lp.max():.2f}")
print(f"  LP mean ± SD: {eicu_surv.lp.mean():.2f} ± {eicu_surv.lp.std():.2f}")

# =============================================================================
# 3. Discrimination on eICU
# =============================================================================
print("\n[3] Discrimination on eICU")
c_eicu = concordance_index(eicu_surv['time_from_lm'], -eicu_surv['lp'], eicu_surv['event_ind'])
print(f"  C-index (MIMIC LP in eICU): {c_eicu:.4f}")
print(f"  Δ C-index (eICU - MIMIC):   {c_eicu - c_mimic:.4f}")
print(f"  MIMIC train C-index bootstrap 95%CI: ({c_mimic_lo:.4f}, {c_mimic_hi:.4f})")
if c_eicu >= c_mimic_lo:
    print(f"  → eICU C-index falls WITHIN MIMIC bootstrap 95%CI → generalisation OK")
else:
    print(f"  → eICU C-index below MIMIC bootstrap 95%CI → some optimism/dataset drift")

# Bootstrap C-index CI for eICU
boots_c_eicu = []
eicu_idx_arr = eicu_surv.index.values
for b in range(200):
    idx = rng.choice(eicu_idx_arr, size=len(eicu_idx_arr), replace=True)
    sub = eicu_surv.loc[idx]
    try:
        c = concordance_index(sub['time_from_lm'], -sub['lp'], sub['event_ind'])
        boots_c_eicu.append(c)
    except Exception:
        pass
c_eicu_lo, c_eicu_hi = np.percentile(boots_c_eicu, [2.5, 97.5])
print(f"  eICU C-index bootstrap 95%CI: {c_eicu:.4f} ({c_eicu_lo:.4f}-{c_eicu_hi:.4f})")

# Time-dependent AUC at 7d / 14d / 28d
def td_auc(df, covars_df, lp, t):
    """AUC for binary event by time t: cases = events by t,
    controls = subjects still at risk past t (time_from_lm > t).
    Used for both MIMIC and eICU with their respective fitted Cox."""
    case = df['event_ind']==1
    case_by_t = case & (df['time_from_lm'] <= t)
    control_t = df['time_from_lm'] > t
    if case_by_t.sum() < 5 or control_t.sum() < 5:
        return float('nan')
    yt = case_by_t.astype(int).values
    ys = lp.values
    try:
        return roc_auc_score(yt, ys)
    except Exception:
        return float('nan')

print("\n  Time-dependent AUC (incident/dynamic approximation):")
auc_7d_mimic  = td_auc(mimic_surv, mimic_surv[COVARS], mimic_surv['lp'], 7*24)
auc_14d_mimic = td_auc(mimic_surv, mimic_surv[COVARS], mimic_surv['lp'], 14*24)
auc_28d_mimic = td_auc(mimic_surv, mimic_surv[COVARS], mimic_surv['lp'], 28*24)
auc_7d_eicu   = td_auc(eicu_surv,  eicu_surv[COVARS],  eicu_surv['lp'],  7*24)
auc_14d_eicu  = td_auc(eicu_surv,  eicu_surv[COVARS],  eicu_surv['lp'],  14*24)
auc_28d_eicu  = td_auc(eicu_surv,  eicu_surv[COVARS],  eicu_surv['lp'],  28*24)
print(f"    7d  MIMIC={auc_7d_mimic:.4f}  eICU={auc_7d_eicu:.4f}  Δ={auc_7d_eicu-auc_7d_mimic:.4f}")
print(f"    14d MIMIC={auc_14d_mimic:.4f}  eICU={auc_14d_eicu:.4f}  Δ={auc_14d_eicu-auc_14d_mimic:.4f}")
print(f"    28d MIMIC={auc_28d_mimic:.4f}  eICU={auc_28d_eicu:.4f}  Δ={auc_28d_eicu-auc_28d_mimic:.4f}")

# =============================================================================
# 4. Calibration on eICU (recalibration regression)
# =============================================================================
print("\n[4] Calibration: refit Cox on eICU using LP alone")
# Cox recalibration: refit LP on eICU with single covariate 'lp'
# output: calibration_slope (should be ~1) + calibration_intercept (should be ~0)
cox_rec = CoxPHFitter(penalizer=0.001)
eicu_surv['lp_mean0'] = eicu_surv['lp'] - eicu_surv['lp'].mean()
cox_rec.fit(eicu_surv[['lp_mean0','time_from_lm','event_ind']],
            duration_col='time_from_lm', event_col='event_ind')
cal_slope = float(cox_rec.params_['lp_mean0'])
cal_intercept = -cox_rec.params_['lp_mean0'] * eicu_surv['lp'].mean()  # = baseline shift
hr_rec = float(np.exp(cox_rec.params_['lp_mean0']))
se_rec = float(cox_rec.summary.loc['lp_mean0','se(coef)'])
p_rec  = float(cox_rec.summary.loc['lp_mean0','p'])
ci_lo_rec = float(np.exp(cox_rec.params_['lp_mean0'] - 1.96*se_rec))
ci_hi_rec = float(np.exp(cox_rec.params_['lp_mean0'] + 1.96*se_rec))
print(f"  Calibration slope (β per unit LP, centered): {cal_slope:.4f}  (ideal=1.0; HR={np.exp(cal_slope):.4f})")
print(f"  95% CI: {ci_lo_rec:.4f}-{ci_hi_rec:.4f} (HR scale)")
print(f"  Calibration intercept (offset):    {cal_intercept:.4f}  (ideal=0)")
print(f"  Cox recalibration C-index: {cox_rec.concordance_index_:.4f}")

# Predicted vs actual mortality by decile of LP (in-time = 30d)
def km_at_t(group, t):
    sub = group.copy()
    sub['t_clip'] = sub['time_from_lm'].clip(upper=t)
    sub['e_clip'] = sub['event_ind'] * (sub['time_from_lm'] <= t)
    kmf = KaplanMeierFitter()
    kmf.fit(sub['t_clip'], sub['e_clip'])
    return float(1 - kmf.predict(t))

print("\n  Decile calibration table (30d mortality):")
eicu_surv['lp_decile'] = pd.qcut(eicu_surv['lp'], 10, labels=False, duplicates='drop') + 1
cal_rows = []
for d in sorted(eicu_surv['lp_decile'].unique()):
    g = eicu_surv[eicu_surv['lp_decile']==d]
    actual = km_at_t(g, TRUNC_H - LANDMARK)
    # Predicted: predict_partial_hazard from Cox_rec
    g_pred_partial = cox_rec.predict_partial_hazard(g[['lp_mean0']])
    base_surv = cox_rec.predict_survival_function(g[['lp_mean0']], times=[TRUNC_H - LANDMARK])
    predicted = float(1 - base_surv.iloc[0].mean())
    cal_rows.append({
        'decile': int(d),
        'n': len(g),
        'events': int(g.event_ind.sum()),
        'actual_30d': round(actual*100, 2),
        'predicted_30d': round(predicted*100, 2),
    })
    print(f"    D{int(d):2d}: n={len(g):>5,}  events={g.event_ind.sum():>4,}  "
          f"actual={actual*100:5.1f}%  predicted={predicted*100:5.1f}%  diff={actual*100-predicted*100:+.1f}pp")

cal_df = pd.DataFrame(cal_rows)

# Hosmer-Lemeshow test on deciles (low df=8)
hl_chi2 = float(((cal_df['actual_30d']/100 - cal_df['predicted_30d']/100)**2 / (cal_df['predicted_30d']/100/(1-cal_df['predicted_30d']/100)/cal_df['n'])).sum())
# This is a rough HL on deciles — the more standard HL uses 10 groups of predicted prob,
# here we use LP deciles as proxy. Higher HL means WORSE fit.
# Use scipy chi2_sf
hl_p = 1 - stats.chi2.cdf(hl_chi2, df=len(cal_df)-2)
print(f"\n  Hosmer-Lemeshow (deciles, χ²={hl_chi2:.2f}, df={len(cal_df)-2}, p={hl_p:.3f})")
if hl_p < 0.05:
    print(f"  → p<0.05 indicates miscalibration")
else:
    print(f"  → p≥0.05 indicates acceptable calibration")

# Brier score at 28d (use 28d to ensure time_from_lm has controls past t)
T_BRIER = 28*24
def brier_at_t(df, t):
    sub = df.copy()
    sub['e_clip'] = sub['event_ind'] * (sub['time_from_lm'] <= t)
    # Predicted prob at t using FULL COVARS (not just LP) — model's perspective
    base_surv = cox_C.predict_survival_function(sub[COVARS], times=[t])
    pred_prob = (1 - base_surv.iloc[0].values)
    y_true = sub['e_clip'].astype(int).values
    return float(brier_score_loss(y_true, pred_prob))

brier_mimic = brier_at_t(mimic_surv, T_BRIER)
brier_eicu  = brier_at_t(eicu_surv,  T_BRIER)
print(f"\n  Brier score at 28d: MIMIC={brier_mimic:.4f}  eICU={brier_eicu:.4f}")

# =============================================================================
# 5. Decision Curve Analysis
# =============================================================================
print("\n[5] Decision Curve Analysis")
def net_benefit(df, threshold, t):
    sub = df.copy()
    sub['e_clip'] = sub['event_ind'] * (sub['time_from_lm'] <= t)
    base_surv = cox_C.predict_survival_function(df[COVARS], times=[t])
    pred_prob = (1 - base_surv.iloc[0].values)
    y = sub['e_clip'].astype(int).values
    # Standard net benefit formula (Vickers 2006):
    # NB = TP/N - FP/N * (pt / (1-pt))
    # where pt = threshold
    n = len(y)
    prevalence = y.mean()
    if threshold <= 0 or threshold >= 1: return float('nan')
    nb = prevalence - (1-prevalence) * (threshold / (1-threshold))
    # Actual model-based NB: for each patient, predict prob > threshold => treat
    treat = pred_prob >= threshold
    tp = (treat & (y==1)).sum()
    fp = (treat & (y==0)).sum()
    nb_model = tp/n - fp/n * (threshold / (1-threshold))
    return float(nb_model), float(nb)

dca_rows = []
print(f"  {'threshold':>10}  {'NB model':>9}  {'NB treat-all':>11}  {'NB treat-none':>13}")
for thr in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]:
    nb_m, nb_a = net_benefit(eicu_surv, thr, TRUNC_H - LANDMARK)
    nb_none = 0.0
    print(f"  {thr:>10.2f}  {nb_m:>9.4f}  {nb_a:>11.4f}  {nb_none:>13.4f}")
    dca_rows.append({'threshold': thr, 'nb_model': nb_m, 'nb_all': nb_a, 'nb_none': nb_none})
dca_df = pd.DataFrame(dca_rows)

# =============================================================================
# 6. Cross-cohort parameter consistency (eICU refit)
# =============================================================================
print("\n[6] Cross-cohort parameter consistency: refit Model C on eICU and compare HRs")
cox_C_eicu = CoxPHFitter(penalizer=0.001)
cox_C_eicu.fit(eicu_surv[COVARS + ['time_from_lm','event_ind']],
               duration_col='time_from_lm', event_col='event_ind')
beta_eicu = cox_C_eicu.params_
print(f"  eICU Model C refit: C-index={cox_C_eicu.concordance_index_:.4f}, "
      f"n={cox_C_eicu._n_examples:,}")
print(f"  {'variable':18s}  {'HR MIMIC':>9}  {'HR eICU':>9}  {'ratio':>6}  {'eICU p':>9}")
hr_consistency = []
for c in COVARS:
    hr_m = float(np.exp(beta_mimic[c]))
    hr_e = float(np.exp(beta_eicu[c]))
    p_e  = float(cox_C_eicu.summary.loc[c,'p'])
    ratio = hr_e / hr_m
    print(f"  {c:18s}  {hr_m:>9.3f}  {hr_e:>9.3f}  {ratio:>6.3f}  {p_e:>9.2e}")
    hr_consistency.append({
        'variable': c,
        'hr_mimic': round(hr_m, 4),
        'hr_mimic_lo': round(float(np.exp(beta_mimic[c] - 1.96*cox_C.summary.loc[c,'se(coef)'])), 4),
        'hr_mimic_hi': round(float(np.exp(beta_mimic[c] + 1.96*cox_C.summary.loc[c,'se(coef)'])), 4),
        'hr_eicu': round(hr_e, 4),
        'hr_eicu_lo': round(float(np.exp(beta_eicu[c] - 1.96*cox_C_eicu.summary.loc[c,'se(coef)'])), 4),
        'hr_eicu_hi': round(float(np.exp(beta_eicu[c] + 1.96*cox_C_eicu.summary.loc[c,'se(coef)'])), 4),
        'hr_ratio': round(ratio, 4),
        'p_eicu': p_e,
    })

# =============================================================================
# 7. KM by LP tertile in eICU
# =============================================================================
print("\n[7] KM by LP tertile in eICU")
eicu_surv['lp_tertile'] = pd.qcut(eicu_surv['lp'], 3, labels=['Low','Med','High'])
kmf_t = {}
km_results_t = []
for grp in ['Low','Med','High']:
    g = eicu_surv[eicu_surv['lp_tertile']==grp]
    kmf = KaplanMeierFitter()
    kmf.fit(g['time_from_lm'], g['event_ind'], label=f'LP_{grp}')
    kmf_t[grp] = kmf
    mort30 = float(1 - kmf.predict(TRUNC_H - LANDMARK))
    km_results_t.append({
        'tertile': grp, 'n': int(len(g)),
        'events': int(g.event_ind.sum()),
        'mort_30d_pct': round(mort30*100, 2)
    })
    print(f"    {grp}: n={len(g):>5,}  events={g.event_ind.sum():>4,}  "
          f"30d mortality={mort30*100:.2f}%")

# Overall logrank
from lifelines.statistics import multivariate_logrank_test
lr = multivariate_logrank_test(eicu_surv['time_from_lm'], eicu_surv['lp_tertile'], eicu_surv['event_ind'])
print(f"  Log-rank (3 tertiles): χ²={lr.test_statistic:.2f}, p={lr.p_value:.2e}")

# =============================================================================
# 8. Calibration plot
# =============================================================================
print("\n[8] Plot calibration & KM by LP tertile")
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

# 8a. Calibration plot
ax = axes[0]
ax.plot([0, 50], [0, 50], 'k--', alpha=0.5, label='Ideal')
ax.scatter(cal_df['predicted_30d'], cal_df['actual_30d'], s=55, color='C0', edgecolor='k', label='eICU deciles')
ax.plot(cal_df['predicted_30d'], cal_df['actual_30d'], 'C0-', alpha=0.6)
ax.set_xlabel('Predicted 30d mortality (%)')
ax.set_ylabel('Observed 30d mortality (%)')
ax.set_title(f'Calibration (eICU)\nslope={cal_slope:.3f} (ideal=1), HL p={hl_p:.3f}')
ax.legend(loc='upper left')
ax.grid(alpha=0.3)

# 8b. KM by LP tertile
ax = axes[1]
colors = {'Low':'C2','Med':'C1','High':'C3'}
for grp in ['Low','Med','High']:
    kmf = kmf_t[grp]
    ax.step(kmf.survival_function_.index, kmf.survival_function_[f'LP_{grp}']*100,
            where='post', color=colors[grp], label=f'LP {grp} (n={sum(eicu_surv.lp_tertile==grp):,})')
ax.set_xlabel('Days from landmark (24h after ICU admission)')
ax.set_ylabel('Survival (%)')
ax.set_title(f'KM by LP tertile (eICU)\nlog-rank p={lr.p_value:.2e}')
ax.set_xlim(0, 30)
ax.set_ylim(0, 100)
ax.legend(loc='lower left', fontsize=8)
ax.grid(alpha=0.3)

# 8c. DCA
ax = axes[2]
ax.plot(dca_df['threshold'], dca_df['nb_model']*100, 'C0-o', label='MIMIC model in eICU')
ax.plot(dca_df['threshold'], dca_df['nb_all']*100,    'C7--', label='Treat all')
ax.axhline(0, color='k', ls=':', label='Treat none')
ax.set_xlabel('Threshold probability (30d mortality)')
ax.set_ylabel('Net benefit (%)')
ax.set_title('Decision Curve (eICU)')
ax.legend(loc='upper right', fontsize=8)
ax.grid(alpha=0.3)

plt.tight_layout()
out_png = os.path.join(WORK, "v2_outputs", "fig_external_validation.png")
plt.savefig(out_png, dpi=150, bbox_inches='tight')
out_pdf = os.path.join(WORK, "v2_outputs", "fig_external_validation.pdf")
plt.savefig(out_pdf, bbox_inches='tight')
plt.close()
print(f"  Saved: {out_png}")
print(f"  Saved: {out_pdf}")

# =============================================================================
# 9. Save JSON + report
# =============================================================================
print("\n[9] Save results")
results = {
    'n_mimic_train': int(len(mimic_surv)),
    'n_eicu_valid':  int(len(eicu_surv)),
    'covars': COVARS,
    'mimic_model_C': {
        'concordance': float(c_mimic),
        'concordance_95ci_lo': float(c_mimic_lo),
        'concordance_95ci_hi': float(c_mimic_hi),
        'hr': {c: round(float(np.exp(beta_mimic[c])), 4) for c in COVARS},
        'hr_lo': {c: round(float(np.exp(beta_mimic[c] - 1.96*cox_C.summary.loc[c,'se(coef)'])), 4) for c in COVARS},
        'hr_hi': {c: round(float(np.exp(beta_mimic[c] + 1.96*cox_C.summary.loc[c,'se(coef)'])), 4) for c in COVARS},
        'p': {c: float(cox_C.summary.loc[c,'p']) for c in COVARS},
    },
    'eicu_external': {
        'concordance_locked_lp': float(c_eicu),
        'concordance_95ci_lo': float(c_eicu_lo),
        'concordance_95ci_hi': float(c_eicu_hi),
        'delta_c_vs_mimic': float(c_eicu - c_mimic),
        'c_falls_in_mimic_ci': bool(c_eicu >= c_mimic_lo),
        'td_auc_7d_mimic':  float(auc_7d_mimic),
        'td_auc_7d_eicu':   float(auc_7d_eicu),
        'td_auc_14d_mimic': float(auc_14d_mimic),
        'td_auc_14d_eicu':  float(auc_14d_eicu),
        'td_auc_28d_mimic': float(auc_28d_mimic),
        'td_auc_28d_eicu':  float(auc_28d_eicu),
        'brier_30d_mimic':  float(brier_mimic),
        'brier_30d_eicu':   float(brier_eicu),
    },
    'calibration': {
        'slope': cal_slope,
        'slope_hr': hr_rec,
        'slope_95ci_lo': ci_lo_rec,
        'slope_95ci_hi': ci_hi_rec,
        'slope_p': p_rec,
        'intercept': cal_intercept,
        'recalibration_C': float(cox_rec.concordance_index_),
        'hosmer_lemeshow_chi2': hl_chi2,
        'hosmer_lemeshow_p':   hl_p,
        'decile_table': cal_rows,
    },
    'dca': dca_rows,
    'cross_cohort_consistency': hr_consistency,
    'km_by_lp_tertile': km_results_t,
    'logrank_chi2': float(lr.test_statistic),
    'logrank_p':     float(lr.p_value),
}
with open(os.path.join(WORK, "v2_outputs", "external_validation_results.json"), 'w') as f:
    json.dump(results, f, indent=2)
print(f"  Saved: v2_outputs/external_validation_results.json")

# Markdown report
report = []
report.append("# F7 真外部验证报告 — 锁定 MIMIC Model C 在 eICU-CRD")
report.append("")
report.append(f"- MIMIC 训练集: n={len(mimic_surv):,} (events={int(mimic_surv.event_ind.sum()):,})")
report.append(f"- eICU 验证集: n={len(eicu_surv):,} (events={int(eicu_surv.event_ind.sum()):,})")
report.append(f"- 协变量集 (n={len(COVARS)}): {', '.join(COVARS)}")
report.append(f"- 时间窗: 24h landmark + 30d 截断")
report.append("")
report.append("## 1. Discrimination")
report.append("")
report.append(f"| 指标 | MIMIC train | eICU valid | Δ |")
report.append(f"|---|---|---|---|")
report.append(f"| C-index (Harrell) | {c_mimic:.4f} ({c_mimic_lo:.4f}-{c_mimic_hi:.4f}) | {c_eicu:.4f} ({c_eicu_lo:.4f}-{c_eicu_hi:.4f}) | {c_eicu-c_mimic:+.4f} |")
report.append(f"| Time-dep AUC 7d | {auc_7d_mimic:.4f} | {auc_7d_eicu:.4f} | {auc_7d_eicu-auc_7d_mimic:+.4f} |")
report.append(f"| Time-dep AUC 14d | {auc_14d_mimic:.4f} | {auc_14d_eicu:.4f} | {auc_14d_eicu-auc_14d_mimic:+.4f} |")
report.append(f"| Time-dep AUC 28d | {auc_28d_mimic:.4f} | {auc_28d_eicu:.4f} | {auc_28d_eicu-auc_28d_mimic:+.4f} |")
report.append(f"| Brier score (28d) | {brier_mimic:.4f} | {brier_eicu:.4f} | {brier_eicu-brier_mimic:+.4f} |")
report.append("")
disc_ok = c_eicu >= c_mimic_lo
report.append(f"**结论 (Discrimination)**: eICU C-index {'落在' if disc_ok else '低于'} MIMIC bootstrap 95%CI → 泛化 {'成功' if disc_ok else '下降'} (Δ = {c_eicu-c_mimic:+.4f})")
report.append("")
report.append("## 2. Calibration")
report.append("")
report.append(f"- 校准斜率 (slope = β per unit LP): **{cal_slope:.4f}** (HR={np.exp(cal_slope):.4f}, 95%CI {ci_lo_rec:.4f}-{ci_hi_rec:.4f}, p={p_rec:.2e})")
report.append(f"- 理想 = 1.0；<1.0 表示模型过度自信，>1.0 表示效应低估")
report.append(f"- Hosmer-Lemeshow (deciles): χ²={hl_chi2:.2f}, p={hl_p:.3f} → {'欠校准' if hl_p<0.05 else '校准良好'}")
report.append("")
report.append("### 十分位校准表")
report.append("")
report.append("| Decile | n | Events | Predicted 30d | Observed 30d | Δ |")
report.append("|---|---|---|---|---|---|")
for r in cal_rows:
    diff = r['actual_30d'] - r['predicted_30d']
    report.append(f"| D{r['decile']} | {r['n']:,} | {r['events']:,} | {r['predicted_30d']:.1f}% | {r['actual_30d']:.1f}% | {diff:+.1f}pp |")
report.append("")
report.append("## 3. Cross-cohort parameter consistency")
report.append("")
report.append("| Variable | HR(MIMIC) | HR(eICU) | Ratio(eICU/MIMIC) | eICU p |")
report.append("|---|---|---|---|---|")
for r in hr_consistency:
    report.append(f"| {r['variable']} | {r['hr_mimic']:.3f} ({r['hr_mimic_lo']:.3f}-{r['hr_mimic_hi']:.3f}) | {r['hr_eicu']:.3f} ({r['hr_eicu_lo']:.3f}-{r['hr_eicu_hi']:.3f}) | {r['hr_ratio']:.3f} | {r['p_eicu']:.2e} |")
report.append("")
report.append("**Ratio ∈ [0.5, 2.0]** 表示参数稳定性可接受。")
report.append("")
report.append("## 4. KM by LP tertile (eICU)")
report.append("")
report.append("| Tertile | n | Events | 30d mortality |")
report.append("|---|---|---|---|")
for r in km_results_t:
    report.append(f"| {r['tertile']} | {r['n']:,} | {r['events']:,} | {r['mort_30d_pct']:.2f}% |")
report.append(f"| Log-rank: χ²={lr.test_statistic:.2f}, p={lr.p_value:.2e}")
report.append("")
report.append("## 5. Decision Curve Analysis")
report.append("")
report.append("| Threshold | Net benefit (model) | Net benefit (treat all) |")
report.append("|---|---|---|")
for r in dca_rows:
    report.append(f"| {r['threshold']:.2f} | {r['nb_model']*100:.2f}% | {r['nb_all']*100:.2f}% |")
report.append("")
report.append("---")
report.append(f"Generated: 2026-09-23 | MIMIC train = {len(mimic_surv):,} | eICU valid = {len(eicu_surv):,}")

with open(os.path.join(WORK, "v2_outputs", "external_validation_report.md"), 'w') as f:
    f.write('\n'.join(report))
print(f"  Saved: v2_outputs/external_validation_report.md")
print(f"\nDone.")