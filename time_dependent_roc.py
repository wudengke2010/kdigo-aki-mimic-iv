#!/usr/bin/env python3
"""
Time-dependent ROC Analysis for KDIGO AKI Staging
==================================================
Computes 7-day, 14-day, and 28-day time-dependent AUCs
for Model B (KDIGO + severity adjusters).

Author: Dengke Wu
Date: 2026-08-08
"""

import pandas as pd
import numpy as np
import json
import os
from sklearn.metrics import roc_auc_score, roc_curve
from lifelines import CoxPHFitter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTPUT_DIR = 'C:/Users/admin/WorkBuddy/2026-07-06-13-22-19'

print("=" * 60)
print("Time-dependent ROC Analysis")
print("=" * 60)

# Load MIMIC-IV analysis data (need to reconstruct from existing data)
# Use the eICU validation data as it has all needed variables
# Also load MIMIC results for comparison

# Check if we have the MIMIC analysis data saved
mimic_data_path = os.path.join(OUTPUT_DIR, 'mimic_analysis_data.csv')
eicu_data_path = os.path.join(OUTPUT_DIR, 'eicu_validation_data.csv')

if os.path.exists(eicu_data_path):
    df = pd.read_csv(eicu_data_path)
    print(f"Loaded eICU data: {len(df):,} patients")
else:
    print("ERROR: No analysis data found")
    exit()

# Prepare variables
df['kdigo_1'] = (df['kdigo_stage'] == 1).astype(int)
df['kdigo_2'] = (df['kdigo_stage'] == 2).astype(int)
df['kdigo_3'] = (df['kdigo_stage'] == 3).astype(int)
df['gender_male'] = (df['gender'] == 'Male').astype(int)

# Ensure duration and event
df['duration'] = df['icu_los_days'].clip(upper=30)
df['event'] = df['mortality']

# Cox model for time-dependent predictions
print("\nFitting Cox model for time-dependent ROC...")
cox_vars = ['kdigo_1', 'kdigo_2', 'kdigo_3', 'apachescore', 'mech_vent', 'age_num', 'gender_male']
cox_data = df[['duration', 'event'] + cox_vars].dropna().astype(float)

cph = CoxPHFitter()
cph.fit(cox_data, duration_col='duration', event_col='event')

# Get risk scores (linear predictor)
risk_scores = cph.predict_partial_hazard(cox_data[cox_vars])

# Time-dependent AUC using cumulative sensitivity / dynamic specificity
# Method: At each time t, subjects are classified as:
#   - Cases: event <= t
#   - Controls: event > t (or censored after t)
# Calculate AUC of risk score for this binary classification

time_points = [7, 14, 28]  # days
td_auc_results = {}

print("\nTime-dependent AUC results:")
for t in time_points:
    # Cases: died before time t
    cases = cox_data[(cox_data['duration'] <= t) & (cox_data['event'] == 1)]
    # Controls: survived past time t
    controls = cox_data[cox_data['duration'] > t]

    if len(cases) > 0 and len(controls) > 0:
        case_scores = risk_scores[cases.index]
        control_scores = risk_scores[controls.index]

        # Combine
        y_binary = np.concatenate([np.ones(len(cases)), np.zeros(len(controls))])
        scores = np.concatenate([case_scores.values.ravel(), control_scores.values.ravel()])

        auc = roc_auc_score(y_binary, scores)
        fpr, tpr, _ = roc_curve(y_binary, scores)

        td_auc_results[f'{t}_day'] = {
            'auc': round(auc, 3),
            'n_cases': int(len(cases)),
            'n_controls': int(len(controls)),
            'fpr': fpr.tolist()[:50],  # Downsample for plotting
            'tpr': tpr.tolist()[:50]
        }
        print(f"  {t}-day AUC: {auc:.3f} (cases={len(cases):,}, controls={len(controls):,})")
    else:
        print(f"  {t}-day: insufficient data (cases={len(cases)}, controls={len(controls)})")

# Also compute KDIGO-only AUC (unadjusted)
print("\nKDIGO-only (unadjusted) time-dependent AUC:")
for t in time_points:
    cases = cox_data[(cox_data['duration'] <= t) & (cox_data['event'] == 1)]
    controls = cox_data[cox_data['duration'] > t]

    if len(cases) > 0 and len(controls) > 0:
        case_kdigo = cox_data.loc[cases.index, ['kdigo_1', 'kdigo_2', 'kdigo_3']].sum(axis=1)
        control_kdigo = cox_data.loc[controls.index, ['kdigo_1', 'kdigo_2', 'kdigo_3']].sum(axis=1)

        y_binary = np.concatenate([np.ones(len(cases)), np.zeros(len(controls))])
        scores = np.concatenate([case_kdigo.values, control_kdigo.values])

        auc = roc_auc_score(y_binary, scores)
        td_auc_results[f'{t}_day_kdigo_only'] = {'auc': round(auc, 3)}
        print(f"  {t}-day AUC (KDIGO only): {auc:.3f}")

# Plot time-dependent ROC curves
print("\nPlotting time-dependent ROC curves...")
fig, ax = plt.subplots(1, 1, figsize=(8, 7))

colors = ['#2166ac', '#ef8a62', '#b2182b']
for i, t in enumerate(time_points):
    key = f'{t}_day'
    if key in td_auc_results and 'fpr' in td_auc_results[key]:
        fpr = np.array(td_auc_results[key]['fpr'])
        tpr = np.array(td_auc_results[key]['tpr'])
        auc_val = td_auc_results[key]['auc']
        ax.plot(fpr, tpr, color=colors[i], linewidth=2.5,
                label=f'{t}-day mortality (AUC = {auc_val:.3f})')

ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.5)
ax.set_xlabel('1 - Specificity', fontsize=13)
ax.set_ylabel('Sensitivity', fontsize=13)
ax.set_title('Time-dependent ROC Curves\n(Model B: KDIGO + APACHE + Mech. Vent + Age + Sex)', fontsize=12)
ax.legend(loc='lower right', fontsize=11, frameon=True, fancybox=True, shadow=True)
ax.set_xlim([-0.01, 1.01])
ax.set_ylim([-0.01, 1.01])
ax.grid(True, alpha=0.3)
ax.set_aspect('equal')

plt.tight_layout()
fig_path = os.path.join(OUTPUT_DIR, 'figures_prism', 'figS6_td_roc.png')
fig.savefig(fig_path, dpi=300, bbox_inches='tight')
fig.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
plt.close()
print(f"  Figure saved: {fig_path}")

# Save results
output_path = os.path.join(OUTPUT_DIR, 'td_roc_results.json')
with open(output_path, 'w', encoding='utf-8') as f:
    json.dump(td_auc_results, f, indent=2)
print(f"  Results saved: {output_path}")

print("\n" + "=" * 60)
print("Time-dependent ROC Analysis COMPLETE")
print("=" * 60)
