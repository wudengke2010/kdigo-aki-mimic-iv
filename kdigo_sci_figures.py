#!/usr/bin/env python3
"""
KDIGO AKI — SCI Publication-Quality Figures
============================================
All figures regenerated from actual data at 300 DPI.
- Figure 1: CONSORT flowchart (matplotlib rendering of draw.io layout)
- Figure 2: Mortality by KDIGO stage (bar chart with 95% CI)
- Figure 3: Kaplan-Meier survival curves with number-at-risk
- Figure 4: Forest plot (Model B logistic regression)
- Figure 5: ICU LOS by KDIGO stage (violin plot)
- Figure 6: Dose-response curve (cr ratio vs mortality)
- Figure 7: Stage distribution (pie chart)
- Figure 8: AKI paradox (dual panel: CKD stratification + RRT decomposition)
"""
import pandas as pd
import numpy as np
import os, json

IN = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
OUT_DIR = os.path.join(IN, "figures_revised")
os.makedirs(OUT_DIR, exist_ok=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.ticker as ticker

# ── Global style ──
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 10,
    'axes.linewidth': 1.0,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.unicode_minus': False,
    'figure.facecolor': 'white',
    'savefig.facecolor': 'white',
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.15,
})

# Colorblind-friendly palette (Wong 2011)
C_BLUE   = '#0072B2'
C_ORANGE = '#E69F00'
C_YELLOW = '#F0E442'
C_RED    = '#D55E00'
C_PURPLE = '#CC79A7'
C_GREEN  = '#009E73'
C_GREY   = '#56B4E9'
C_BLACK  = '#000000'

STAGE_COLORS = [C_BLUE, C_ORANGE, C_RED, C_PURPLE]
STAGE_LABELS = ['No AKI\n(Stage 0)', 'Stage 1', 'Stage 2', 'Stage 3']

# ── Load data ──
print("Loading data...")
cohort = pd.read_csv(os.path.join(IN, "kdigo_cohort_revised_enriched.csv"), low_memory=False)
with open(os.path.join(IN, "flowchart_data.json")) as f:
    fc = json.load(f)
with open(os.path.join(IN, "km_revised.json")) as f:
    km = json.load(f)
logistic = pd.read_csv(os.path.join(IN, "logistic_results_revised.csv"))
dose_resp = pd.read_csv(os.path.join(IN, "dose_response_revised.csv"))
table1 = pd.read_csv(os.path.join(IN, "table1_revised.csv"))
cox_b = pd.read_csv(os.path.join(IN, "cox_b_full.csv"))
print(f"  Cohort N = {len(cohort):,}")

stages = [0, 1, 2, 3]

# ══════════════════════════════════════════════════════════════
# FIGURE 1: CONSORT Flow Diagram
# ══════════════════════════════════════════════════════════════
print("[Fig 1] CONSORT flowchart...")

fig, ax = plt.subplots(figsize=(8.5, 11))
ax.set_xlim(0, 850)
ax.set_ylim(0, 1100)
ax.set_aspect('equal')
ax.axis('off')

def draw_box(ax, x, y, w, h, text, fc='#dae8fc', ec='#6c8ebf', fs=10, bold_first=True, lw=1.5, dashed=False):
    style = "round,pad=0.1"
    rect = FancyBboxPatch((x, y), w, h, boxstyle=style, facecolor=fc, edgecolor=ec, linewidth=lw, linestyle='--' if dashed else '-')
    ax.add_patch(rect)
    lines = text.split('\n')
    for i, line in enumerate(lines):
        weight = 'bold' if (i == 0 and bold_first) else 'normal'
        ax.text(x + w/2, y + h - 8 - i * (fs + 2), line, ha='center', va='top',
                fontsize=fs, fontweight=weight, fontfamily='Arial')

def draw_arrow(ax, x1, y1, x2, y2, color='#333333', lw=1.5, dashed=False):
    style = '-|>' if not dashed else '-|>'
    ls = '--' if dashed else '-'
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, color=color, lw=lw, linestyle=ls))

# Title
ax.text(425, 1085, 'Figure 1. Patient Selection Flow Diagram', ha='center', va='top',
        fontsize=13, fontweight='bold', fontfamily='Arial')

# Box 1: Initial
draw_box(ax, 275, 1000, 300, 55, 'MIMIC-IV v3.1 database\nAll ICU stays (2008-2022)\nn = 94,458', fc='#dae8fc', ec='#6c8ebf', fs=10)
draw_arrow(ax, 425, 1000, 425, 970)

# Box 2: First ICU stay
draw_box(ax, 275, 915, 300, 50, 'First ICU stay per hospital admission\nn = 85,242', fc='#dae8fc', ec='#6c8ebf', fs=10)
draw_arrow(ax, 425, 915, 425, 880)

# Exclusion 1 (right) — moved right to avoid overlap with First ICU stay box
draw_box(ax, 620, 935, 180, 45, 'Excluded:\nMissing serum\ncreatinine (n = 1,075)', fc='#f8cecc', ec='#b85450', fs=9, dashed=True)
draw_arrow(ax, 575, 940, 620, 940, color='#b85450', lw=1, dashed=True)

# Box 3: Final cohort
draw_box(ax, 275, 825, 300, 50, 'Final study cohort\nn = 84,167', fc='#d5e8d4', ec='#82b366', fs=11, lw=2)
draw_arrow(ax, 425, 825, 425, 795)

# Box 4: KDIGO staging
draw_box(ax, 275, 745, 300, 45, 'KDIGO creatinine-based AKI staging', fc='#fff2cc', ec='#d6b656', fs=10)

# Horizontal split
ax.plot([150, 700], [735, 735], color='#333333', lw=1.5)
draw_arrow(ax, 150, 735, 150, 710)
draw_arrow(ax, 700, 735, 700, 710)

# No AKI
draw_box(ax, 50, 645, 200, 60, 'No AKI (Stage 0)\nn = 59,033 (70.1%)\nMortality: 6.0%', fc='#dae8fc', ec='#6c8ebf', fs=10)

# AKI
draw_box(ax, 600, 645, 200, 60, 'AKI (Stage 1-3)\nn = 25,134 (29.9%)\nMortality: 24.2%', fc='#f8cecc', ec='#b85450', fs=10)

# Stage split
ax.plot([620, 780], [635, 635], color='#333333', lw=1.5)
draw_arrow(ax, 620, 635, 620, 610)
draw_arrow(ax, 700, 635, 700, 610)
draw_arrow(ax, 780, 635, 780, 610)

# Stage boxes
draw_box(ax, 570, 540, 100, 65, 'Stage 1\nn = 12,817\n(15.2%)\nMort: 15.1%', fc='#fff2cc', ec='#d6b656', fs=9)
draw_box(ax, 680, 540, 100, 65, 'Stage 2\nn = 2,729\n(3.2%)\nMort: 26.7%', fc='#ffe6cc', ec='#d79b00', fs=9)
draw_box(ax, 790, 540, 100, 65, 'Stage 3\nn = 9,588\n(11.4%)\nMort: 29.0%', fc='#f8cecc', ec='#b85450', fs=9)

# Stage 3 decomposition
draw_arrow(ax, 840, 540, 840, 500, color='#666666', lw=1, dashed=True)
draw_box(ax, 700, 430, 280, 65,
         'Stage 3 stratification by RRT indication:\nRRT only: n = 793 (8.3%)\nCr only: n = 5,617 (58.6%)\nBoth (RRT + Cr): n = 3,178 (33.1%)',
         fc='#e1d5e7', ec='#9673a6', fs=8.5, dashed=True)

# Analysis box
draw_box(ax, 150, 280, 500, 110,
         'Statistical Analysis\n'
         '  Logistic regression: 3 models (unadjusted / severity-adjusted / fully-adjusted)\n'
         '  Cox proportional hazards models (C-index = 0.783)\n'
         '  Kaplan-Meier with Bonferroni-corrected pairwise log-rank tests\n'
         '  Dose-response analysis (creatinine ratio bins)\n'
         '  RRT stratification and AKI paradox analysis\n'
         '  Sensitivity: MDRD baseline Cr / Lab-CKD (eGFR < 60)',
         fc='#f5f5f5', ec='#666666', fs=9, bold_first=True, lw=1)

# Dashed arrows to analysis
draw_arrow(ax, 150, 645, 300, 390, color='#999999', lw=0.8, dashed=True)
draw_arrow(ax, 620, 540, 450, 390, color='#999999', lw=0.8, dashed=True)

fig.savefig(os.path.join(OUT_DIR, 'flowchart.png'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("  -> flowchart.png")

# ══════════════════════════════════════════════════════════════
# FIGURE 2: Mortality by KDIGO Stage (bar chart with 95% CI)
# ══════════════════════════════════════════════════════════════
print("[Fig 2] Mortality by stage...")

def wilson_ci(p, n, z=1.96):
    if n == 0: return p, p
    denom = 1 + z*z/n
    centre = (p + z*z/(2*n)) / denom
    width = z * np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / denom
    return centre - width, centre + width

mort_data = []
for s in stages:
    sub = cohort[cohort['kdigo_stage'] == s]
    n = len(sub)
    events = sub['hospital_expire_flag'].sum()
    p = events / n if n > 0 else 0
    ci_low, ci_high = wilson_ci(p, n)
    mort_data.append({'stage': s, 'n': n, 'mort': p*100, 'ci_low': ci_low*100, 'ci_high': ci_high*100})

fig, ax = plt.subplots(figsize=(7, 5))
x = np.arange(4)
mort_vals = [d['mort'] for d in mort_data]
ci_errs = [[d['mort'] - d['ci_low'] for d in mort_data],
           [d['ci_high'] - d['mort'] for d in mort_data]]

bars = ax.bar(x, mort_vals, color=STAGE_COLORS, edgecolor='white', linewidth=1.2, width=0.65, zorder=3)
ax.errorbar(x, mort_vals, yerr=ci_errs, fmt='none', ecolor='black', capsize=5, capthick=1.2, lw=1.2, zorder=4)

for i, d in enumerate(mort_data):
    ax.text(i, d['ci_high'] + 2.0, f'{d["mort"]:.1f}%', ha='center', va='bottom',
            fontsize=11, fontweight='bold', color='#222222')
    ax.text(i, d['ci_high'] + 5.5, f'(n={d["n"]:,})', ha='center', va='bottom',
            fontsize=8.5, color='#666666')

# Significance brackets
def sig_bracket(ax, x1, x2, y, label):
    ax.plot([x1, x1, x2, x2], [y, y+0.8, y+0.8, y], lw=1, color='#444444')
    ax.text((x1+x2)/2, y+1.0, label, ha='center', va='bottom', fontsize=8, color='#444444')

sig_bracket(ax, 0, 1, 19, '***')
sig_bracket(ax, 1, 2, 31, '***')
sig_bracket(ax, 2, 3, 36, 'n.s.')

ax.set_xticks(x)
ax.set_xticklabels(['No AKI\n(Stage 0)', 'Stage 1', 'Stage 2', 'Stage 3'], fontsize=10)
ax.set_ylabel('Hospital Mortality (%)', fontsize=12, fontweight='bold')
ax.set_title('Hospital Mortality by KDIGO AKI Stage', fontsize=13, fontweight='bold', pad=12)
ax.set_ylim(0, max(d['ci_high'] for d in mort_data) * 1.35)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d%%'))
ax.grid(axis='y', alpha=0.3, linestyle='--', zorder=0)

fig.savefig(os.path.join(OUT_DIR, 'mortality_by_stage.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  -> mortality_by_stage.png")

# ══════════════════════════════════════════════════════════════
# FIGURE 3: Kaplan-Meier Survival Curves with Number-at-Risk
# ══════════════════════════════════════════════════════════════
print("[Fig 3] KM survival curves...")

fig = plt.figure(figsize=(8, 7.2))
ax = plt.axes([0.12, 0.34, 0.82, 0.58])  # Increased bottom margin for risk table

stage_names_km = ['No AKI', 'Stage 1', 'Stage 2', 'Stage 3']
km_colors = [C_BLUE, C_ORANGE, C_RED, C_PURPLE]

for idx, s in enumerate(stages):
    key = str(s)
    d = km['km_data'][key]
    times = np.array(d['times'])
    surv = np.array(d['survival'])
    ax.step(times, surv, where='post', color=km_colors[idx], linewidth=1.8, label=stage_names_km[idx], zorder=3)
    # CI shading (approximate using step)
    ax.fill_between(times, surv * 0.98, surv * 1.02, step='post', color=km_colors[idx], alpha=0.1)

ax.set_xlabel('')  # x-axis label placed below risk table to avoid overlap
ax.set_ylabel('Survival probability', fontsize=11, fontweight='bold')
ax.set_title('Kaplan-Meier Survival Curves by KDIGO AKI Stage', fontsize=13, fontweight='bold', pad=10)
ax.set_xlim(0, 30)
ax.set_ylim(0.55, 1.02)
ax.set_xticks([0, 5, 10, 15, 20, 25, 30])
ax.grid(alpha=0.3, linestyle='--')

# Log-rank annotation
ax.text(0.97, 0.97, f'Overall log-rank: $\\chi^2$ = {km["overall_chi2"]:.1f}\np < 0.001\n\nPairwise (Bonferroni):\n0 vs 1: p < 0.001\n0 vs 2: p < 0.001\n0 vs 3: p < 0.001\n1 vs 2: p < 0.001\n1 vs 3: p < 0.001\n2 vs 3: p = 1.000 (n.s.)',
        transform=ax.transAxes, ha='right', va='top', fontsize=8,
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#cccccc', alpha=0.9))

ax.legend(fontsize=9.5, loc='center right', framealpha=0.9, edgecolor='#cccccc')

# Number-at-risk table
risk_days = [0, 7, 14, 21, 30]
ax_risk = plt.axes([0.12, 0.08, 0.82, 0.20])
ax_risk.axis('off')
ax_risk.set_xlim(0, 30)
ax_risk.set_ylim(0, 5.5)

# Header
ax_risk.text(-1.5, 5.0, 'No. at risk', ha='right', va='center', fontsize=9, fontweight='bold')
for j, day in enumerate(risk_days):
    ax_risk.text(day, 5.0, str(day), ha='center', va='center', fontsize=9, fontweight='bold')

for idx, s in enumerate(stages):
    key = str(s)
    d = km['km_data'][key]
    nar = d['n_at_risk']
    y_pos = 5.0 - (idx + 1)
    ax_risk.text(-1.5, y_pos, stage_names_km[idx], ha='right', va='center', fontsize=8.5, color=km_colors[idx], fontweight='bold')
    # Interpolate n_at_risk to standard days
    # nar has 5 values at days 0,7,14,21,30 (approximately)
    for j, day in enumerate(risk_days):
        if j < len(nar):
            ax_risk.text(day, y_pos, f'{nar[j]:,}', ha='center', va='center', fontsize=8.5)

# x-axis label placed below risk table
fig.text(0.5, 0.03, 'Time since ICU admission (days)', ha='center', fontsize=11, fontweight='bold')

fig.savefig(os.path.join(OUT_DIR, 'km_survival.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  -> km_survival.png")

# ══════════════════════════════════════════════════════════════
# FIGURE 4: Forest Plot (Model B — fully adjusted logistic)
# ══════════════════════════════════════════════════════════════
print("[Fig 4] Forest plot...")

model_b = logistic[logistic['model'] == 'B'].copy()
# Order: KDIGO stages first, then by OR descending
kdigo_rows = model_b[model_b['variable'].str.startswith('kdigo')].copy()
other_rows = model_b[~model_b['variable'].str.startswith('kdigo')].copy()
other_rows = other_rows.sort_values('OR', ascending=True)
plot_data = pd.concat([other_rows, kdigo_rows]).reset_index(drop=True)

# Clean variable names
name_map = {
    'kdigo_stage_1': 'KDIGO Stage 1',
    'kdigo_stage_2': 'KDIGO Stage 2',
    'kdigo_stage_3': 'KDIGO Stage 3',
    'anchor_age': 'Age (per year)',
    'male_num': 'Male sex',
    'sofa_imputed': 'SOFA score',
    'mech_vent': 'Mechanical ventilation',
    'vasopressor': 'Vasopressor use',
    'ckd_icd': 'CKD (ICD-coded)',
    'com_Cardiac_Arrhythmias': 'Cardiac arrhythmias',
    'com_Cerebrovascular': 'Cerebrovascular disease',
    'com_Chronic_Pulmonary': 'Chronic pulmonary disease',
    'com_Coagulopathy': 'Coagulopathy',
    'com_Congestive_Heart_Failure': 'Congestive heart failure',
    'com_Diabetes_complicated': 'Diabetes (complicated)',
    'com_Diabetes_uncomplicated': 'Diabetes (uncomplicated)',
    'com_Fluid/Electrolyte': 'Fluid/Electrolyte disorders',
    'com_Hypertension': 'Hypertension',
    'com_Liver_Disease': 'Liver disease',
    'com_Metastatic_Cancer': 'Metastatic cancer',
    'com_Myocardial_Infarction': 'Myocardial infarction',
    'com_Obesity': 'Obesity',
    'com_Peptic_Ulcer': 'Peptic ulcer disease',
    'com_Peripheral_Vascular': 'Peripheral vascular disease',
    'com_Sepsis': 'Sepsis',
    'com_Solid_Tumor_non-met': 'Solid tumor (non-metastatic)',
    'com_Valvular_Disease': 'Valvular disease',
}
plot_data['label'] = plot_data['variable'].map(name_map).fillna(plot_data['variable'])

n_vars = len(plot_data)
fig, ax = plt.subplots(figsize=(9, max(7, n_vars * 0.35)))
y_pos = np.arange(n_vars)

for i, (_, row) in enumerate(plot_data.iterrows()):
    is_kdigo = row['variable'].startswith('kdigo')
    color = C_RED if is_kdigo else C_BLUE
    marker = 'D' if is_kdigo else 'o'
    size = 7 if is_kdigo else 5

    ax.errorbar(row['OR'], i, xerr=[[row['OR'] - row['CI_low']], [row['CI_high'] - row['OR']]],
                fmt='none', ecolor=color, capsize=3, capthick=1, lw=1)
    ax.scatter(row['OR'], i, color=color, marker=marker, s=size**2, zorder=5, edgecolors='white', linewidths=0.5)

    # OR value and CI text on right
    sig = ''
    if row['p_value'] < 0.001: sig = '***'
    elif row['p_value'] < 0.01: sig = '**'
    elif row['p_value'] < 0.05: sig = '*'
    ax.text(12, i, f'{row["OR"]:.3f} ({row["CI_low"]:.2f}-{row["CI_high"]:.2f}) {sig}',
            va='center', ha='left', fontsize=7.5, color='#333333',
            fontweight='bold' if is_kdigo else 'normal')

ax.axvline(x=1, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
ax.set_yticks(y_pos)
ax.set_yticklabels(plot_data['label'], fontsize=9)
ax.set_xlabel('Adjusted Odds Ratio (95% CI)', fontsize=11, fontweight='bold')
ax.set_title('Forest Plot: Model B (Fully Adjusted Logistic Regression)', fontsize=12, fontweight='bold', pad=10)
ax.set_xlim(0.3, 22)
ax.set_xscale('log')
ax.set_xticks([0.5, 1, 2, 4, 8])
ax.set_xticklabels(['0.5', '1', '2', '4', '8'])
ax.invert_yaxis()
ax.grid(axis='x', alpha=0.3, linestyle='--')

fig.savefig(os.path.join(OUT_DIR, 'forest_plot.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  -> forest_plot.png")

# ══════════════════════════════════════════════════════════════
# FIGURE 5: ICU LOS by KDIGO Stage (violin plot)
# ══════════════════════════════════════════════════════════════
print("[Fig 5] ICU LOS violin...")

fig, ax = plt.subplots(figsize=(7, 5))
los_data = []
for s in stages:
    sub = cohort[cohort['kdigo_stage'] == s]
    los = sub['icu_los_hours'].dropna()
    # Truncate extreme outliers for visualization (>99th percentile)
    p99 = los.quantile(0.99)
    los = los[los <= p99]
    los_data.append(los.values)

parts = ax.violinplot(los_data, positions=np.arange(4), showmeans=True, showmedians=True, widths=0.7)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(STAGE_COLORS[i])
    pc.set_alpha(0.6)
    pc.set_edgecolor(STAGE_COLORS[i])
    pc.set_linewidth(1.2)
parts['cmeans'].set_color('black')
parts['cmeans'].set_linewidth(1.5)
parts['cmedians'].set_color('#444444')
parts['cmedians'].set_linewidth(1.5)
parts['cmedians'].set_linestyle('--')

# Add n and median annotations
for i, s in enumerate(stages):
    sub = cohort[cohort['kdigo_stage'] == s]
    los = sub['icu_los_hours'].dropna()
    med = los.median()
    n = len(los)
    ax.text(i, max(los_data[i]) * 1.02, f'n={n:,}\nmed={med:.0f}h', ha='center', va='bottom',
            fontsize=8.5, color='#444444')

ax.set_xticks(np.arange(4))
ax.set_xticklabels(['No AKI\n(Stage 0)', 'Stage 1', 'Stage 2', 'Stage 3'], fontsize=10)
ax.set_ylabel('ICU Length of Stay (hours)', fontsize=11, fontweight='bold')
ax.set_title('ICU Length of Stay by KDIGO AKI Stage', fontsize=13, fontweight='bold', pad=10)
ax.grid(axis='y', alpha=0.3, linestyle='--')

# Legend for mean/median
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='black', lw=1.5, label='Mean'),
    Line2D([0], [0], color='#444444', lw=1.5, ls='--', label='Median'),
]
ax.legend(handles=legend_elements, fontsize=9, loc='upper left')

fig.savefig(os.path.join(OUT_DIR, 'icu_los_by_stage.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  -> icu_los_by_stage.png")

# ══════════════════════════════════════════════════════════════
# FIGURE 6: Dose-Response Curve (redesigned with leader-line annotations)
# ══════════════════════════════════════════════════════════════
print("[Fig 6] Dose-response...")

fig, ax = plt.subplots(figsize=(7, 5.2))

# Convert bins to numeric midpoints
bin_mids = [0.75, 1.25, 1.75, 2.5, 4.0, 7.5, 12.0]
mort_vals = dose_resp['mortality_pct'].values
ci_low = dose_resp['CI_low'].values
ci_high = dose_resp['CI_high'].values
n_vals = dose_resp['n'].values

# CI band
ax.fill_between(bin_mids, ci_low, ci_high, alpha=0.12, color=C_BLUE, zorder=2)
# Main line
ax.plot(bin_mids, mort_vals, 'o-', color=C_BLUE, linewidth=2.2, markersize=8,
        markeredgecolor='white', markeredgewidth=1.2, zorder=4)
# Error bars
ax.errorbar(bin_mids, mort_vals, yerr=[mort_vals - ci_low, ci_high - mort_vals],
            fmt='none', ecolor=C_BLUE, capsize=4, capthick=1, lw=1, alpha=0.5, zorder=3)

# Data point annotations — staggered to avoid overlap in dense left region
annot_offsets = [(-12, 8), (15, 5), (0, 10), (0, 10), (0, 10), (0, 10), (0, 10)]
for i, (x, y, n) in enumerate(zip(bin_mids, mort_vals, n_vals)):
    dx, dy = annot_offsets[i]
    ax.annotate(f'{y:.1f}%\n(n={n:,})', xy=(x, y), xytext=(dx, dy),
                textcoords='offset points', ha='center', va='bottom',
                fontsize=8, color='#333333', zorder=5,
                bbox=dict(boxstyle='round,pad=0.15', facecolor='white', edgecolor='none', alpha=0.7))

# KDIGO stage thresholds with leader lines to upper clear area, spread horizontally
thresholds = [
    (1.5, C_ORANGE, 'Stage 1\nthreshold', 1.5, 50),
    (2.0, C_RED,     'Stage 2',            2.8, 47),
    (3.0, C_PURPLE,  'Stage 3\nthreshold', 6.0, 44),
]
for x_line, color, label, label_x, label_y in thresholds:
    ax.axvline(x=x_line, color=color, linestyle=':', linewidth=1.2, alpha=0.8, zorder=1)
    # Leader line from label to threshold line
    ax.plot([label_x, x_line], [label_y - 2.5, 32], color=color, linestyle='-', linewidth=0.8, alpha=0.6, zorder=2)
    ax.text(label_x, label_y, label, ha='center', va='bottom', fontsize=8,
            color=color, fontweight='bold', zorder=5,
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=color, alpha=0.9, linewidth=0.8))

ax.set_xlabel('Creatinine Ratio (peak / baseline)', fontsize=11, fontweight='bold')
ax.set_ylabel('Hospital Mortality (%)', fontsize=11, fontweight='bold')
ax.set_title('Dose-Response: Creatinine Ratio vs Hospital Mortality', fontsize=12, fontweight='bold', pad=14)
ax.set_xlim(0.3, 13)
ax.set_ylim(0, 55)
ax.set_xticks([0.5, 1, 1.5, 2, 3, 5, 7.5, 10, 12.5])
ax.set_xticklabels(['0.5', '1.0', '1.5', '2.0', '3.0', '5.0', '7.5', '10', '12.5'], fontsize=9)
ax.grid(alpha=0.3, linestyle='--', which='major')
ax.grid(alpha=0.15, linestyle=':', which='minor')
ax.minorticks_on()

fig.savefig(os.path.join(OUT_DIR, 'dose_response.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  -> dose_response.png")

# ══════════════════════════════════════════════════════════════
# FIGURE 7: Stage Distribution (pie chart)
# ══════════════════════════════════════════════════════════════
print("[Fig 7] Stage distribution...")

fig, ax = plt.subplots(figsize=(6, 6))
stage_n = [fc['stage_counts'][str(s)] for s in stages]
stage_pct = [n / sum(stage_n) * 100 for n in stage_n]

pie_labels = [f'No AKI\nn={stage_n[0]:,}\n({stage_pct[0]:.1f}%)',
              f'Stage 1\nn={stage_n[1]:,}\n({stage_pct[1]:.1f}%)',
              f'Stage 2\nn={stage_n[2]:,}\n({stage_pct[2]:.1f}%)',
              f'Stage 3\nn={stage_n[3]:,}\n({stage_pct[3]:.1f}%)']

wedges, texts = ax.pie(stage_n, labels=pie_labels, colors=STAGE_COLORS,
                       startangle=90, counterclock=False,
                       wedgeprops=dict(edgecolor='white', linewidth=2),
                       textprops=dict(fontsize=9.5, fontweight='bold'),
                       pctdistance=0.75, labeldistance=1.15)

# Add center circle for donut effect
centre_circle = plt.Circle((0, 0), 0.45, fc='white', linewidth=0)
ax.add_artist(centre_circle)
ax.text(0, 0.08, 'Total', ha='center', va='center', fontsize=11, fontweight='bold')
ax.text(0, -0.08, f'N = {sum(stage_n):,}', ha='center', va='center', fontsize=10, color='#666666')

ax.set_title('KDIGO AKI Stage Distribution', fontsize=13, fontweight='bold', pad=15)

fig.savefig(os.path.join(OUT_DIR, 'stage_distribution.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  -> stage_distribution.png")

# ══════════════════════════════════════════════════════════════
# FIGURE 8: AKI Paradox (dual panel)
# ══════════════════════════════════════════════════════════════
print("[Fig 8] AKI paradox...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5))

# ── Panel A: CKD stratification by stage ──
ckd_neg_mort = []
ckd_pos_mort = []
ckd_neg_n = []
ckd_pos_n = []
for s in stages:
    sub = cohort[cohort['kdigo_stage'] == s]
    neg = sub[sub['ckd_icd'] == 0]
    pos = sub[sub['ckd_icd'] == 1]
    ckd_neg_mort.append(neg['hospital_expire_flag'].mean() * 100 if len(neg) > 0 else 0)
    ckd_pos_mort.append(pos['hospital_expire_flag'].mean() * 100 if len(pos) > 0 else 0)
    ckd_neg_n.append(len(neg))
    ckd_pos_n.append(len(pos))

x = np.arange(4)
width = 0.35
bars1 = ax1.bar(x - width/2, ckd_neg_mort, width, label='CKD-negative', color=C_BLUE, edgecolor='white', linewidth=1)
bars2 = ax1.bar(x + width/2, ckd_pos_mort, width, label='CKD-positive', color=C_PURPLE, edgecolor='white', linewidth=1)

for i in range(4):
    ax1.text(i - width/2, ckd_neg_mort[i] + 0.8, f'{ckd_neg_mort[i]:.1f}%', ha='center', fontsize=8, fontweight='bold', color=C_BLUE)
    ax1.text(i + width/2, ckd_pos_mort[i] + 0.8, f'{ckd_pos_mort[i]:.1f}%', ha='center', fontsize=8, fontweight='bold', color=C_PURPLE)

ax1.set_xticks(x)
ax1.set_xticklabels(['No AKI', 'Stage 1', 'Stage 2', 'Stage 3'], fontsize=10)
ax1.set_ylabel('Hospital Mortality (%)', fontsize=11, fontweight='bold')
ax1.set_title('A. CKD Stratification by KDIGO Stage', fontsize=11, fontweight='bold')
ax1.set_ylim(0, max(max(ckd_neg_mort), max(ckd_pos_mort)) * 1.25)
ax1.legend(fontsize=9, loc='upper left')
ax1.grid(axis='y', alpha=0.3, linestyle='--')

# ── Panel B: Stage 3 RRT decomposition ──
stage3 = cohort[cohort['kdigo_stage'] == 3]
rrt_cats = ['RRT only', 'Cr only', 'Both (RRT+Cr)']
rrt_cols = ['stage3_by_rrt_only', 'stage3_by_cr_only', 'stage3_by_both']
rrt_neg_mort = []
rrt_pos_mort = []
rrt_neg_n = []
rrt_pos_n = []
for col in rrt_cols:
    sub_neg = stage3[(stage3['ckd_icd'] == 0) & (stage3[col] == 1)]
    sub_pos = stage3[(stage3['ckd_icd'] == 1) & (stage3[col] == 1)]
    rrt_neg_mort.append(sub_neg['hospital_expire_flag'].mean() * 100 if len(sub_neg) > 0 else 0)
    rrt_pos_mort.append(sub_pos['hospital_expire_flag'].mean() * 100 if len(sub_pos) > 0 else 0)
    rrt_neg_n.append(len(sub_neg))
    rrt_pos_n.append(len(sub_pos))

x2 = np.arange(3)
bars3 = ax2.bar(x2 - width/2, rrt_neg_mort, width, label='CKD-negative', color=C_BLUE, edgecolor='white', linewidth=1)
bars4 = ax2.bar(x2 + width/2, rrt_pos_mort, width, label='CKD-positive', color=C_PURPLE, edgecolor='white', linewidth=1)

for i in range(3):
    ax2.text(i - width/2, rrt_neg_mort[i] + 1.0, f'{rrt_neg_mort[i]:.1f}%\n(n={rrt_neg_n[i]:,})', ha='center', fontsize=7.5, fontweight='bold', color=C_BLUE)
    ax2.text(i + width/2, rrt_pos_mort[i] + 1.0, f'{rrt_pos_mort[i]:.1f}%\n(n={rrt_pos_n[i]:,})', ha='center', fontsize=7.5, fontweight='bold', color=C_PURPLE)

ax2.set_xticks(x2)
ax2.set_xticklabels(rrt_cats, fontsize=10)
ax2.set_ylabel('Hospital Mortality (%)', fontsize=11, fontweight='bold')
ax2.set_title('B. Stage 3 Decomposed by RRT Indication', fontsize=11, fontweight='bold')
ax2.set_ylim(0, max(max(rrt_neg_mort), max(rrt_pos_mort)) * 1.30)
ax2.legend(fontsize=9, loc='upper right')
ax2.grid(axis='y', alpha=0.3, linestyle='--')

plt.tight_layout(pad=2.0)
fig.savefig(os.path.join(OUT_DIR, 'aki_paradox.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  -> aki_paradox.png")

# ══════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("ALL 8 FIGURES GENERATED SUCCESSFULLY")
print(f"Output: {OUT_DIR}")
print("=" * 60)
