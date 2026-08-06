#!/usr/bin/env python3
"""
KDIGO AKI — Supplementary Figures
==================================
Figure S1: RCS nonlinear OR curve
Figure S2: ROC comparison (A) + DCA (B)
Figure S3: Sepsis subgroup forest plot
Figure S4: Calibration curve
"""
import pandas as pd
import numpy as np
import os, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import FancyBboxPatch

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.labelcolor'] = '#333333'
plt.rcParams['xtick.color'] = '#333333'
plt.rcParams['ytick.color'] = '#333333'
plt.rcParams['axes.titleweight'] = 'bold'

IN = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
OUT = os.path.join(IN, "figures_revised")
os.makedirs(OUT, exist_ok=True)

# Color palette (colorblind-friendly, Wong 2011)
C_BLUE = '#0072B2'
C_ORANGE = '#E69F00'
C_RED = '#D55E00'
C_PURPLE = '#7B2D8D'
C_GREEN = '#009E73'
C_GRAY = '#999999'
C_BLACK = '#333333'

# ============================================================
# Figure S1: RCS Nonlinear OR Curve
# ============================================================
print("[Fig S1] RCS nonlinear OR curve...")
rcs = pd.read_csv(os.path.join(IN, 'rcs_results.csv'))

fig, ax = plt.subplots(figsize=(7, 5))

# CI band
ax.fill_between(rcs['cr_ratio'], rcs['OR_low'], rcs['OR_high'],
                alpha=0.15, color=C_BLUE, zorder=2)
# OR curve
ax.plot(rcs['cr_ratio'], rcs['OR'], color=C_BLUE, linewidth=2.2, zorder=3)

# Reference line at OR=1
ax.axhline(y=1, color=C_GRAY, linestyle='--', linewidth=0.8, alpha=0.5, zorder=1)

# KDIGO stage thresholds (staggered vertically to avoid overlap)
thresholds = [
    (1.5, C_ORANGE, 'Stage 1', 5.3),
    (2.0, C_RED, 'Stage 2', 4.5),
    (3.0, C_PURPLE, 'Stage 3', 5.3),
]
for x, color, label, ytext in thresholds:
    ax.axvline(x=x, color=color, linestyle=':', linewidth=1, alpha=0.6, zorder=1)
    ax.text(x, ytext, label, fontsize=8, color=color, fontweight='bold', rotation=0, ha='center')

# Nonlinearity test annotation
with open(os.path.join(IN, 'model_comparison.json')) as f:
    mc = json.load(f)
nl = mc['RCS_nonlinearity']
ax.text(0.97, 0.97,
        f'Restricted cubic spline\n4 knots (Cr ratio: 1.0, 1.5, 2.0, 3.0)\n'
        f'Nonlinearity test:\n$\\chi^2$ = {nl["chi2"]:.1f}, df = {nl["df"]}, p < 0.001',
        transform=ax.transAxes, ha='right', va='top', fontsize=8,
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#cccccc', alpha=0.9))

ax.set_xlabel('Creatinine Ratio (peak / baseline)', fontsize=11, fontweight='bold')
ax.set_ylabel('Adjusted Odds Ratio (vs. Cr ratio = 1.0)', fontsize=11, fontweight='bold')
ax.set_title('Restricted Cubic Spline: Cr Ratio vs Mortality', fontsize=12, fontweight='bold', pad=10)
ax.set_xlim(0.5, 12)
ax.set_ylim(0.5, 6)
ax.grid(alpha=0.25, linestyle='--')

fig.savefig(os.path.join(OUT, 'fig_s1_rcs.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  -> fig_s1_rcs.png")

# ============================================================
# Figure S2: ROC Comparison (A) + DCA (B)
# ============================================================
print("[Fig S2] ROC comparison + DCA...")
roc = pd.read_csv(os.path.join(IN, 'roc_comparison.csv'))
dca = pd.read_csv(os.path.join(IN, 'dca_results.csv'))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: ROC curves
ax1.plot(roc['fpr_a'], roc['tpr_a'], color=C_ORANGE, linewidth=1.8, label=f'Model A (AUC = 0.778)')
ax1.plot(roc['fpr_b'], roc['tpr_b'], color=C_BLUE, linewidth=2.2, label=f'Model B (AUC = 0.828)')
ax1.plot(roc['fpr_c'], roc['tpr_c'], color=C_GREEN, linewidth=1.8, linestyle='--',
         label=f'Model C (AUC = 0.825)')
ax1.plot([0, 1], [0, 1], color=C_GRAY, linestyle=':', linewidth=0.8)

ax1.set_xlabel('1 - Specificity', fontsize=10, fontweight='bold')
ax1.set_ylabel('Sensitivity', fontsize=10, fontweight='bold')
ax1.set_title('A. ROC Curves', fontsize=11, fontweight='bold')
ax1.legend(fontsize=8.5, loc='lower right', framealpha=0.9, edgecolor='#cccccc')
ax1.set_xlim(-0.01, 1.01)
ax1.set_ylim(-0.01, 1.01)
ax1.grid(alpha=0.2, linestyle='--')

# Panel B: Decision Curve Analysis
# Clip to reasonable threshold range (0.01 - 0.50)
mask = (dca['threshold'] >= 0.01) & (dca['threshold'] <= 0.50)
dca_plot = dca[mask]

ax2.plot(dca_plot['threshold'], dca_plot['nb_model_a'], color=C_ORANGE, linewidth=1.8,
         label='Model A')
ax2.plot(dca_plot['threshold'], dca_plot['nb_model_b'], color=C_BLUE, linewidth=2.2,
         label='Model B')
ax2.plot(dca_plot['threshold'], dca_plot['nb_model_c'], color=C_GREEN, linewidth=1.8,
         linestyle='--', label='Model C')
ax2.plot(dca_plot['threshold'], dca_plot['nb_treat_all'], color=C_GRAY, linewidth=1.2,
         linestyle='-.', label='Treat All')
ax2.plot(dca_plot['threshold'], dca_plot['nb_treat_none'], color=C_BLACK, linewidth=1,
         linestyle=':', label='Treat None')

ax2.set_xlabel('Threshold Probability', fontsize=10, fontweight='bold')
ax2.set_ylabel('Net Benefit', fontsize=10, fontweight='bold')
ax2.set_title('B. Decision Curve Analysis', fontsize=11, fontweight='bold')
ax2.legend(fontsize=8, loc='upper right', framealpha=0.9, edgecolor='#cccccc')
ax2.set_xlim(0.01, 0.50)
ax2.set_ylim(-0.05, 0.15)
ax2.grid(alpha=0.2, linestyle='--')

fig.suptitle('Model Discrimination and Clinical Utility', fontsize=13, fontweight='bold', y=1.02)
fig.savefig(os.path.join(OUT, 'fig_s2_roc_dca.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  -> fig_s2_roc_dca.png")

# ============================================================
# Figure S3: Sepsis Subgroup Forest Plot
# ============================================================
print("[Fig S3] Sepsis subgroup forest plot...")
with open(os.path.join(IN, 'sepsis_subgroup.json')) as f:
    sepsis_data = json.load(f)

fig, ax = plt.subplots(figsize=(8, 4.5))

# Data for forest plot
stages = ['Stage 1', 'Stage 2', 'Stage 3']
stage_nums = [1, 2, 3]

sepsis_ors = [sepsis_data['Sepsis']['kdigo_or'][str(s)]['OR'] for s in stage_nums]
sepsis_los = [sepsis_data['Sepsis']['kdigo_or'][str(s)]['CI_low'] for s in stage_nums]
sepsis_his = [sepsis_data['Sepsis']['kdigo_or'][str(s)]['CI_high'] for s in stage_nums]

nonsepsis_ors = [sepsis_data['Non-sepsis']['kdigo_or'][str(s)]['OR'] for s in stage_nums]
nonsepsis_los = [sepsis_data['Non-sepsis']['kdigo_or'][str(s)]['CI_low'] for s in stage_nums]
nonsepsis_his = [sepsis_data['Non-sepsis']['kdigo_or'][str(s)]['CI_high'] for s in stage_nums]

y_pos = np.arange(len(stages))
bar_height = 0.3

# Non-sepsis bars (above)
ax.barh(y_pos + bar_height/2, nonsepsis_ors, height=bar_height, color=C_BLUE, alpha=0.8,
        label=f'Non-sepsis (n={sepsis_data["Non-sepsis"]["n"]:,})', edgecolor='white', linewidth=0.5)
# Sepsis bars (below)
ax.barh(y_pos - bar_height/2, sepsis_ors, height=bar_height, color=C_RED, alpha=0.8,
        label=f'Sepsis (n={sepsis_data["Sepsis"]["n"]:,})', edgecolor='white', linewidth=0.5)

# Error bars
ax.errorbar(nonsepsis_ors, y_pos + bar_height/2,
            xerr=[[o - lo for o, lo in zip(nonsepsis_ors, nonsepsis_los)],
                  [hi - o for o, hi in zip(nonsepsis_ors, nonsepsis_his)]],
            fmt='none', ecolor=C_BLUE, capsize=3, capthick=1, lw=1)
ax.errorbar(sepsis_ors, y_pos - bar_height/2,
            xerr=[[o - lo for o, lo in zip(sepsis_ors, sepsis_los)],
                  [hi - o for o, hi in zip(sepsis_ors, sepsis_his)]],
            fmt='none', ecolor=C_RED, capsize=3, capthick=1, lw=1)

# Labels on bars
for i, (or_val, lo, hi) in enumerate(zip(nonsepsis_ors, nonsepsis_los, nonsepsis_his)):
    ax.text(hi + 0.1, i + bar_height/2, f'{or_val:.2f} ({lo:.2f}-{hi:.2f})',
            va='center', fontsize=7.5, color=C_BLUE)
for i, (or_val, lo, hi) in enumerate(zip(sepsis_ors, sepsis_los, sepsis_his)):
    ax.text(hi + 0.1, i - bar_height/2, f'{or_val:.2f} ({lo:.2f}-{hi:.2f})',
            va='center', fontsize=7.5, color=C_RED)

ax.set_yticks(y_pos)
ax.set_yticklabels(stages, fontsize=10, fontweight='bold')
ax.set_xlabel('Adjusted Odds Ratio (vs. Stage 0)', fontsize=10, fontweight='bold')
ax.set_title('KDIGO Stage ORs: Sepsis vs Non-sepsis Subgroups', fontsize=12, fontweight='bold', pad=10)
ax.axvline(x=1, color=C_GRAY, linestyle='--', linewidth=0.8)
ax.set_xlim(0, max(sepsis_his + nonsepsis_his) * 1.5)
ax.legend(fontsize=9, loc='lower right', framealpha=0.9, edgecolor='#cccccc')
ax.grid(axis='x', alpha=0.2, linestyle='--')

# Interaction p-value annotation
int_data = sepsis_data['interaction']
ax.text(0.97, 0.97,
        f'KDIGO x Sepsis interaction\n$\\chi^2$ = {int_data["chi2"]:.1f}, p < 0.001',
        transform=ax.transAxes, ha='right', va='top', fontsize=8,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#cccccc', alpha=0.9))

fig.savefig(os.path.join(OUT, 'fig_s3_sepsis_subgroup.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  -> fig_s3_sepsis_subgroup.png")

# ============================================================
# Figure S4: Calibration Curve
# ============================================================
print("[Fig S4] Calibration curve...")
cal = pd.read_csv(os.path.join(IN, 'calibration_model_b.csv'))

fig, ax = plt.subplots(figsize=(6, 6))

# Perfect calibration line
ax.plot([0, 1], [0, 1], color=C_GRAY, linestyle='--', linewidth=1, label='Perfect calibration')

# Calibration curve
ax.plot(cal['expected'], cal['observed'], 'o-', color=C_BLUE, linewidth=2, markersize=7,
        markeredgecolor='white', markeredgewidth=1, label='Model B', zorder=3)

# CI for observed
ax.errorbar(cal['expected'], cal['observed'],
            yerr=[cal['observed'] - cal['obs_ci_low'], cal['obs_ci_high'] - cal['observed']],
            fmt='none', ecolor=C_BLUE, capsize=3, capthick=1, lw=1, alpha=0.5)

# Histogram of predicted probabilities (secondary axis)
ax2 = ax.twinx()
ax2.hist(cal['expected'], weights=cal['n'], bins=20, alpha=0.15, color=C_ORANGE, edgecolor='none')
ax2.set_ylabel('Count', fontsize=10, color=C_ORANGE)
ax2.tick_params(axis='y', labelcolor=C_ORANGE)
ax2.set_ylim(0, max(cal['n']) * 3)

# Hosmer-Lemeshow annotation
with open(os.path.join(IN, 'calibration_results.json')) as f:
    cal_results = json.load(f)
hl = cal_results['Model_B']
ax.text(0.05, 0.95,
        f'Hosmer-Lemeshow:\n$\\chi^2$ = {hl["HL_chi2"]:.1f}, p < 0.001\n(large-sample HL test)',
        transform=ax.transAxes, ha='left', va='top', fontsize=8,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#cccccc', alpha=0.9))

ax.set_xlabel('Predicted Probability', fontsize=11, fontweight='bold')
ax.set_ylabel('Observed Mortality Rate', fontsize=11, fontweight='bold')
ax.set_title('Calibration Curve (Model B)', fontsize=12, fontweight='bold', pad=10)
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(-0.02, 1.02)
ax.legend(fontsize=9, loc='lower right', framealpha=0.9, edgecolor='#cccccc')
ax.grid(alpha=0.2, linestyle='--')

fig.savefig(os.path.join(OUT, 'fig_s4_calibration.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  -> fig_s4_calibration.png")

# ============================================================
# Summary
# ============================================================
print("\nAll supplementary figures generated:")
for f in ['fig_s1_rcs.png', 'fig_s2_roc_dca.png', 'fig_s3_sepsis_subgroup.png', 'fig_s4_calibration.png']:
    path = os.path.join(OUT, f)
    size = os.path.getsize(path) / 1024
    print(f"  {f}: {size:.0f} KB")
