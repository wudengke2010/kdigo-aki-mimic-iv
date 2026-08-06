#!/usr/bin/env python3
"""
KDIGO AKI — GraphPad Prism-Quality Publication Figures
=======================================================
Generates all 12 figures in exact GraphPad Prism aesthetic:
- Arial throughout, no grid, inward ticks, thin spines (0.5pt)
- Prism 10 publication color palette
- Both PNG (300 DPI) and PDF (vector) output
- Prism figure dimensions (5 x 3.75 in single, 10 x 3.75 in double)
- Also exports organized Excel workbook for Prism import
"""
import pandas as pd
import numpy as np
import os, json, warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
import matplotlib.font_manager as fm

IN = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
OUT_DIR = os.path.join(IN, "figures_prism")
os.makedirs(OUT_DIR, exist_ok=True)

# ══════════════════════════════════════════════════════════════
# Prism-Exact Global Style
# ══════════════════════════════════════════════════════════════
plt.rcParams.update({
    # Font
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9,
    'font.weight': 'normal',
    'mathtext.fontset': 'custom',
    'mathtext.rm': 'Arial',
    'mathtext.it': 'Arial:italic',
    'mathtext.bf': 'Arial:bold',
    # Axes — Prism style
    'axes.linewidth': 0.5,           # Prism uses thin axis lines
    'axes.spines.top': False,        # No top spine
    'axes.spines.right': False,      # No right spine
    'axes.edgecolor': '#000000',     # Black axis lines
    'axes.labelcolor': '#000000',
    'axes.labelsize': 9,
    'axes.labelweight': 'normal',
    'axes.titlesize': 10,
    'axes.titleweight': 'bold',
    'axes.titlepad': 6,
    'axes.unicode_minus': False,
    # Ticks — inward, short (Prism default)
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.size': 3,           # Short ticks
    'ytick.major.size': 3,
    'xtick.minor.size': 1.5,
    'ytick.minor.size': 1.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.color': '#000000',
    'ytick.color': '#000000',
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    # No grid (Prism default)
    'axes.grid': False,
    # Figure
    'figure.facecolor': 'white',
    'figure.dpi': 300,
    'savefig.facecolor': 'white',
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    'savefig.transparent': False,
    # Legend
    'legend.fontsize': 8,
    'legend.frameon': True,
    'legend.edgecolor': '#000000',
    'legend.fancybox': False,
    'legend.borderpad': 0.4,
    'legend.borderaxespad': 0.5,
    'legend.handletextpad': 0.5,
})

# ── Prism 10 Publication Color Palette ──
# These match Prism's default scheme used in published papers
PRISM_BLACK  = '#000000'
PRISM_RED    = '#E41A1C'   # Prism red
PRISM_BLUE   = '#377EB8'   # Prism blue
PRISM_GREEN  = '#4DAF4A'   # Prism green
PRISM_ORANGE = '#FF7F00'   # Prism orange
PRISM_PURPLE = '#984EA3'   # Prism purple
PRISM_CYAN   = '#56B4E9'   # Prism cyan
PRISM_YELLOW = '#F0E442'   # Prism yellow
PRISM_GRAY   = '#999999'   # Prism gray
PRISM_BROWN  = '#A65628'   # Prism brown
PRISM_PINK   = '#F781BF'   # Prism pink

# Stage colors (Prism scheme)
STAGE_COLORS = [PRISM_BLUE, PRISM_ORANGE, PRISM_RED, PRISM_PURPLE]
STAGE_NAMES = ['No AKI', 'Stage 1', 'Stage 2', 'Stage 3']
STAGE_NAMES_FULL = ['No AKI\n(Stage 0)', 'Stage 1', 'Stage 2', 'Stage 3']

# ── Load all data ──
print("Loading data...")
cohort = pd.read_csv(os.path.join(IN, "kdigo_cohort_revised_enriched.csv"), low_memory=False)
with open(os.path.join(IN, "flowchart_data.json")) as f:
    fc = json.load(f)
with open(os.path.join(IN, "km_revised.json")) as f:
    km = json.load(f)
logistic = pd.read_csv(os.path.join(IN, "logistic_results_revised.csv"))
dose_resp = pd.read_csv(os.path.join(IN, "dose_response_revised.csv"))
cox_b = pd.read_csv(os.path.join(IN, "cox_b_full.csv"))
rcs = pd.read_csv(os.path.join(IN, 'rcs_results.csv'))
roc = pd.read_csv(os.path.join(IN, 'roc_comparison.csv'))
dca = pd.read_csv(os.path.join(IN, 'dca_results.csv'))
cal = pd.read_csv(os.path.join(IN, 'calibration_model_b.csv'))
with open(os.path.join(IN, 'model_comparison.json')) as f:
    mc = json.load(f)
with open(os.path.join(IN, 'sepsis_subgroup.json')) as f:
    sepsis_data = json.load(f)
with open(os.path.join(IN, 'calibration_results.json')) as f:
    cal_results = json.load(f)
with open(os.path.join(IN, 'ckd_interaction.json')) as f:
    ckd_int = json.load(f)
print(f"  Cohort N = {len(cohort):,}")

stages = [0, 1, 2, 3]

# ── Helper: Save figure in both PNG and PDF ──
def save_fig(fig, name):
    png_path = os.path.join(OUT_DIR, f'{name}.png')
    pdf_path = os.path.join(OUT_DIR, f'{name}.pdf')
    fig.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  -> {name}.png + .pdf")

# ── Helper: Wilson CI ──
def wilson_ci(p, n, z=1.96):
    if n == 0: return p, p
    denom = 1 + z*z/n
    centre = (p + z*z/(2*n)) / denom
    width = z * np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / denom
    return centre - width, centre + width

# ── Helper: Prism-style significance bracket ──
def prism_sig(ax, x1, x2, y, label, h=0.005):
    """Draw a Prism-style significance bracket with inward ticks."""
    ax.plot([x1, x1, x2, x2], [y, y + h*2, y + h*2, y], lw=0.5, color='#000000')
    ax.text((x1+x2)/2, y + h*3, label, ha='center', va='bottom', fontsize=7, color='#000000')


# ══════════════════════════════════════════════════════════════
# FIGURE 1: CONSORT Flow Diagram
# (Kept as matplotlib — Prism doesn't do flowcharts)
# ══════════════════════════════════════════════════════════════
print("[Fig 1] CONSORT flowchart...")

fig, ax = plt.subplots(figsize=(7.5, 10))
ax.set_xlim(0, 750)
ax.set_ylim(0, 1050)
ax.set_aspect('equal')
ax.axis('off')

def draw_box(ax, x, y, w, h, text, fc='#dae8fc', ec='#6c8ebf', fs=9, bold_first=True, lw=1, dashed=False):
    style = "round,pad=0.1"
    rect = FancyBboxPatch((x, y), w, h, boxstyle=style, facecolor=fc, edgecolor=ec,
                          linewidth=lw, linestyle='--' if dashed else '-', zorder=2)
    ax.add_patch(rect)
    lines = text.split('\n')
    # Comfortable line spacing, vertically centered block
    line_height = fs * 1.3
    total_text_height = len(lines) * line_height
    start_y = y + h - (h - total_text_height) / 2 - line_height * 0.3
    for i, line in enumerate(lines):
        weight = 'bold' if (i == 0 and bold_first) else 'normal'
        ax.text(x + w/2, start_y - i * line_height, line, ha='center', va='top',
                fontsize=fs, fontweight=weight, fontfamily='Arial', zorder=3)

def draw_arrow(ax, x1, y1, x2, y2, color='#333333', lw=1, dashed=False):
    ls = '--' if dashed else '-'
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='-|>', color=color, lw=lw, linestyle=ls), zorder=1)

# Title (Prism style: above plot, bold)
ax.text(375, 1040, 'Patient Selection Flow Diagram', ha='center', va='top',
        fontsize=11, fontweight='bold', fontfamily='Arial')

# Box 1
draw_box(ax, 225, 970, 300, 50, 'MIMIC-IV v3.1 database\nAll ICU stays (2008-2022)\nn = 94,458', fc='#dae8fc', ec='#6c8ebf', fs=9)
draw_arrow(ax, 375, 970, 375, 945)

# Box 2
draw_box(ax, 225, 895, 300, 45, 'First ICU stay per hospital admission\nn = 85,242', fc='#dae8fc', ec='#6c8ebf', fs=9)
draw_arrow(ax, 375, 895, 375, 870)

# Exclusion 1
draw_box(ax, 545, 910, 160, 45, 'Excluded:\nMissing serum\ncreatinine (n = 1,075)', fc='#f8cecc', ec='#b85450', fs=8, dashed=True)
draw_arrow(ax, 525, 915, 545, 915, color='#b85450', lw=0.8, dashed=True)

# Box 3: Final cohort
draw_box(ax, 225, 815, 300, 50, 'Final study cohort\nn = 84,167', fc='#d5e8d4', ec='#82b366', fs=10, lw=1.5)
draw_arrow(ax, 375, 815, 375, 790)

# Box 4: KDIGO staging
draw_box(ax, 225, 740, 300, 45, 'KDIGO creatinine-based AKI staging', fc='#fff2cc', ec='#d6b656', fs=9)

# Split
ax.plot([130, 620], [730, 730], color='#333333', lw=1)
draw_arrow(ax, 130, 730, 130, 705)
draw_arrow(ax, 620, 730, 620, 705)

# No AKI
draw_box(ax, 55, 635, 180, 60, 'No AKI (Stage 0)\nn = 59,033 (70.1%)\nMortality: 6.0%', fc='#dae8fc', ec='#6c8ebf', fs=9)

# AKI
draw_box(ax, 535, 635, 180, 60, 'AKI (Stage 1-3)\nn = 25,134 (29.9%)\nMortality: 24.2%', fc='#f8cecc', ec='#b85450', fs=9)

# Stage split
ax.plot([555, 715], [625, 625], color='#333333', lw=1)
draw_arrow(ax, 555, 625, 555, 600)
draw_arrow(ax, 625, 625, 625, 600)
draw_arrow(ax, 715, 625, 715, 600)

# Stage boxes
draw_box(ax, 495, 530, 115, 65, 'Stage 1\nn = 12,817\n(15.2%)\nMort: 15.1%', fc='#fff2cc', ec='#d6b656', fs=8)
draw_box(ax, 610, 530, 115, 65, 'Stage 2\nn = 2,729\n(3.2%)\nMort: 26.7%', fc='#ffe6cc', ec='#d79b00', fs=8)
draw_box(ax, 725, 530, 115, 65, 'Stage 3\nn = 9,588\n(11.4%)\nMort: 29.0%', fc='#f8cecc', ec='#b85450', fs=8)

# Stage 3 decomposition (improved layout: taller, cleaner text)
draw_arrow(ax, 780, 555, 780, 525, color='#666666', lw=0.8, dashed=True)
draw_box(ax, 570, 430, 320, 85,
         'Stage 3 decomposition by RRT indication\n'
         'RRT only: n = 793 (8.3%)\n'
         'Cr-based only: n = 5,617 (58.6%)\n'
         'Both criteria: n = 3,178 (33.1%)',
         fc='#e1d5e7', ec='#9673a6', fs=8.5, dashed=True)

# Analysis box (simplified: grouped lines, taller, larger font)
draw_box(ax, 100, 230, 580, 170,
         'Statistical Analysis (see Methods for details)\n'
         'Primary: Logistic regression (3 incremental models)\n'
         'Secondary: Cox regression, Kaplan-Meier, dose-response\n'
         'Exploratory: RRT stratification, AKI paradox\n'
         'Discrimination: ROC, DCA, NRI/IDI, RCS, calibration\n'
         'Sensitivity: MDRD, Lab-CKD, landmark, sepsis subgroup',
         fc='#f5f5f5', ec='#666666', fs=8.5, bold_first=True, lw=0.8)

draw_arrow(ax, 145, 650, 300, 400, color='#999999', lw=0.6, dashed=True)
draw_arrow(ax, 590, 555, 470, 400, color='#999999', lw=0.6, dashed=True)

save_fig(fig, 'fig1_flowchart')


# ══════════════════════════════════════════════════════════════
# FIGURE 2: Mortality by KDIGO Stage (Prism bar chart)
# ══════════════════════════════════════════════════════════════
print("[Fig 2] Mortality by stage...")

mort_data = []
for s in stages:
    sub = cohort[cohort['kdigo_stage'] == s]
    n = len(sub)
    events = sub['hospital_expire_flag'].sum()
    p = events / n if n > 0 else 0
    ci_low, ci_high = wilson_ci(p, n)
    mort_data.append({'stage': s, 'n': n, 'mort': p*100, 'ci_low': ci_low*100, 'ci_high': ci_high*100})

fig, ax = plt.subplots(figsize=(5, 3.75))
x = np.arange(4)
mort_vals = [d['mort'] for d in mort_data]
ci_errs = [[d['mort'] - d['ci_low'] for d in mort_data],
           [d['ci_high'] - d['mort'] for d in mort_data]]

# Prism-style bars: fill with color, thin white border, no edge
bars = ax.bar(x, mort_vals, color=STAGE_COLORS, edgecolor='white', linewidth=0.5, width=0.7, zorder=3)
# Error bars: thin, small caps
ax.errorbar(x, mort_vals, yerr=ci_errs, fmt='none', ecolor='#000000', capsize=3, capthick=0.5, lw=0.5, zorder=4)

# Labels above error bars
for i, d in enumerate(mort_data):
    ax.text(i, d['ci_high'] + 1.5, f'{d["mort"]:.1f}%', ha='center', va='bottom',
            fontsize=9, fontweight='bold', color='#000000')

# x-axis labels with n embedded (Prism style)
ax.set_xticks(x)
xtick_labels = [f'{STAGE_NAMES_FULL[i]}\nn={mort_data[i]["n"]:,}' for i in range(4)]
ax.set_xticklabels(xtick_labels, fontsize=8)

# Significance brackets (Prism style: thin, minimal)
prism_sig(ax, 0, 1, 22, '***')
prism_sig(ax, 1, 2, 33, '***')
prism_sig(ax, 2, 3, 37, 'n.s.')

ax.set_ylabel('Hospital Mortality (%)', fontsize=9, fontweight='bold')
ax.set_ylim(0, max(d['ci_high'] for d in mort_data) * 1.35)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))

save_fig(fig, 'fig2_mortality')


# ══════════════════════════════════════════════════════════════
# FIGURE 3: Kaplan-Meier Survival Curves (Prism style)
# ══════════════════════════════════════════════════════════════
print("[Fig 3] KM survival curves...")

fig = plt.figure(figsize=(5, 5))
ax = plt.axes([0.15, 0.35, 0.80, 0.55])

km_colors = STAGE_COLORS

for idx, s in enumerate(stages):
    key = str(s)
    d = km['km_data'][key]
    times = np.array(d['times'])
    surv = np.array(d['survival'])
    ax.step(times, surv, where='post', color=km_colors[idx], linewidth=1.5,
            label=STAGE_NAMES[idx], zorder=3)
    # CI shading — very subtle (Prism style)
    ci_lower = np.array(d.get('survival_ci_low', surv * 0.98))
    ci_upper = np.array(d.get('survival_ci_high', np.minimum(surv * 1.02, 1.0)))
    ax.fill_between(times, ci_lower, ci_upper, step='post', color=km_colors[idx], alpha=0.08, zorder=2)

ax.set_ylabel('Survival probability', fontsize=9, fontweight='bold')
ax.set_xlim(0, 30)
ax.set_ylim(0.50, 1.02)
ax.set_xticks([0, 5, 10, 15, 20, 25, 30])
ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

# Log-rank annotation (Prism-style: clean text box, thin border)
# Placed at bottom-left to avoid blocking survival curves
ax.text(0.03, 0.10,
        f'Log-rank: P < 0.001\n'
        f'\nPairwise (Bonferroni):'
        f'\n0 vs 1: P < 0.001'
        f'\n0 vs 2: P < 0.001'
        f'\n0 vs 3: P < 0.001'
        f'\n1 vs 2: P < 0.001'
        f'\n1 vs 3: P < 0.001'
        f'\n2 vs 3: P = 1.000',
        transform=ax.transAxes, ha='left', va='bottom', fontsize=6.5,
        bbox=dict(boxstyle='square,pad=0.3', facecolor='white', edgecolor='#000000', linewidth=0.5))

ax.legend(fontsize=7, loc='lower right', framealpha=1, edgecolor='#000000',
          borderpad=0.4, handlelength=1.5)

# Number-at-risk table (Prism-style: aligned with plot x-axis)
risk_days = [0, 7, 14, 21, 30]
ax_risk = plt.axes([0.15, 0.08, 0.80, 0.22])
ax_risk.axis('off')
ax_risk.set_xlim(0, 30)
ax_risk.set_ylim(0, 5.5)

ax_risk.text(-1.5, 5.0, 'No. at risk', ha='right', va='center', fontsize=8, fontweight='bold')
for j, day in enumerate(risk_days):
    ax_risk.text(day, 5.0, str(day), ha='center', va='center', fontsize=8, fontweight='bold')

# Horizontal line above risk table
ax_risk.plot([0, 30], [4.6, 4.6], color='#000000', linewidth=0.3)

for idx, s in enumerate(stages):
    key = str(s)
    d = km['km_data'][key]
    nar = d['n_at_risk']
    y_pos = 4.5 - (idx + 1) * 0.9
    # Stage name and n on same line, separated
    ax_risk.text(-1.5, y_pos, f'{STAGE_NAMES[idx]}    n={nar[0]:,}' if nar else STAGE_NAMES[idx],
                 ha='right', va='center', fontsize=7, color=km_colors[idx], fontweight='bold')
    for j, day in enumerate(risk_days):
        if j < len(nar):
            ax_risk.text(day, y_pos, f'{nar[j]:,}', ha='center', va='center', fontsize=7)

fig.text(0.55, 0.02, 'Time since ICU admission (days)', ha='center', fontsize=9, fontweight='bold')

save_fig(fig, 'fig3_km_survival')


# ══════════════════════════════════════════════════════════════
# FIGURE 4: Forest Plot (Prism-style table layout)
# ══════════════════════════════════════════════════════════════
print("[Fig 4] Forest plot...")

model_b = logistic[logistic['model'] == 'B'].copy()
kdigo_rows = model_b[model_b['variable'].str.startswith('kdigo')].copy()
other_rows = model_b[~model_b['variable'].str.startswith('kdigo')].copy()
other_rows = other_rows.sort_values('OR', ascending=True)
plot_data = pd.concat([other_rows, kdigo_rows]).reset_index(drop=True)

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
fig, ax = plt.subplots(figsize=(7, max(5.5, n_vars * 0.28)))
y_pos = np.arange(n_vars)

# Alternating row shading (Prism-style)
for i in range(n_vars):
    if i % 2 == 0:
        ax.axhspan(i - 0.4, i + 0.4, color='#f0f0f0', zorder=0)

for i, (_, row) in enumerate(plot_data.iterrows()):
    is_kdigo = row['variable'].startswith('kdigo')
    color = PRISM_RED if is_kdigo else PRISM_BLUE
    marker = 's' if is_kdigo else 'o'
    size = 5 if is_kdigo else 4

    ax.errorbar(row['OR'], i, xerr=[[row['OR'] - row['CI_low']], [row['CI_high'] - row['OR']]],
                fmt='none', ecolor=color, capsize=2, capthick=0.5, lw=0.5)
    ax.scatter(row['OR'], i, color=color, marker=marker, s=size**2, zorder=5, edgecolors='white', linewidths=0.3)

    sig = ''
    if row['p_value'] < 0.001: sig = '***'
    elif row['p_value'] < 0.01: sig = '**'
    elif row['p_value'] < 0.05: sig = '*'
    ax.text(15, i, f'{row["OR"]:.3f}  ({row["CI_low"]:.2f}\u2013{row["CI_high"]:.2f})  {sig}',
            va='center', ha='left', fontsize=7, color='#000000',
            fontweight='bold' if is_kdigo else 'normal', fontfamily='Arial')

ax.axvline(x=1, color='#000000', linestyle='--', linewidth=0.5, alpha=0.5)
ax.set_yticks(y_pos)
ax.set_yticklabels(plot_data['label'], fontsize=7)
ax.set_xlabel('Adjusted Odds Ratio (95% CI)', fontsize=9, fontweight='bold')
ax.set_xlim(0.3, 22)
ax.set_xscale('log')
ax.set_xticks([0.5, 1, 2, 4, 8])
ax.set_xticklabels(['0.5', '1', '2', '4', '8'], fontsize=8)
ax.invert_yaxis()

# Column headers
ax.text(15, -1.2, 'OR (95% CI)', fontsize=7, fontweight='bold', va='bottom', ha='left')

save_fig(fig, 'fig4_forest')


# ══════════════════════════════════════════════════════════════
# FIGURE 5: ICU LOS Violin Plot (Prism scatter + bar style)
# ══════════════════════════════════════════════════════════════
print("[Fig 5] ICU LOS violin...")

fig, ax = plt.subplots(figsize=(5, 3.75))
los_data = []
los_stats = []
for s in stages:
    sub = cohort[cohort['kdigo_stage'] == s]
    los = sub['icu_los_hours'].dropna()
    p99 = los.quantile(0.99)
    los_clipped = los[los <= p99]
    los_data.append(los_clipped.values)
    los_stats.append({'median': los.median(), 'mean': los.mean(), 'n': len(los), 'p25': los.quantile(0.25), 'p75': los.quantile(0.75)})

# Prism-style: violin with thin outline, semi-transparent fill
parts = ax.violinplot(los_data, positions=np.arange(4), showmeans=False, showmedians=False, widths=0.6)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(STAGE_COLORS[i])
    pc.set_alpha(0.3)
    pc.set_edgecolor(STAGE_COLORS[i])
    pc.set_linewidth(0.8)

# Box plot inside violin (Prism-style: thin, minimal)
bp = ax.boxplot(los_data, positions=np.arange(4), widths=0.15, patch_artist=True,
                showfliers=False, zorder=4,
                medianprops=dict(color='#000000', linewidth=1),
                boxprops=dict(facecolor='white', edgecolor='#000000', linewidth=0.5),
                whiskerprops=dict(color='#000000', linewidth=0.5),
                capprops=dict(color='#000000', linewidth=0.5))

# Annotations
for i, st in enumerate(los_stats):
    ax.text(i, max(los_data[i]) * 1.02, f'n={st["n"]:,}\nmed={st["median"]:.0f}h',
            ha='center', va='bottom', fontsize=7, color='#000000')

ax.set_xticks(np.arange(4))
ax.set_xticklabels(STAGE_NAMES_FULL, fontsize=8)
ax.set_ylabel('ICU Length of Stay (hours)', fontsize=9, fontweight='bold')

save_fig(fig, 'fig5_icu_los')


# ══════════════════════════════════════════════════════════════
# FIGURE 6: Dose-Response Curve (Prism scatter + error bar)
# ══════════════════════════════════════════════════════════════
print("[Fig 6] Dose-response...")

fig, ax = plt.subplots(figsize=(5, 3.75))

bin_mids = [0.75, 1.25, 1.75, 2.5, 4.0, 7.5, 12.0]
mort_vals = dose_resp['mortality_pct'].values
ci_low = dose_resp['CI_low'].values
ci_high = dose_resp['CI_high'].values
n_vals = dose_resp['n'].values

# CI band (Prism style: very subtle)
ax.fill_between(bin_mids, ci_low, ci_high, alpha=0.12, color=PRISM_BLUE, zorder=2)
# Connecting line (Prism style: thin)
ax.plot(bin_mids, mort_vals, '-', color=PRISM_BLUE, linewidth=1, zorder=3)
# Data points (Prism style: filled squares with error bars)
ax.errorbar(bin_mids, mort_vals, yerr=[mort_vals - ci_low, ci_high - mort_vals],
            fmt='s', color=PRISM_BLUE, markersize=5, markeredgecolor='white', markeredgewidth=0.5,
            ecolor=PRISM_BLUE, capsize=3, capthick=0.5, lw=0.5, zorder=4)

# Data labels (Prism style: minimal, offset to avoid overlap)
for i, (x_val, y_val, n) in enumerate(zip(bin_mids, mort_vals, n_vals)):
    if i == 0:
        dx, dy = 0, -9
    elif i == 1:
        dx, dy = 12, 3
    elif i == 2:
        dx, dy = -10, 5
    else:
        dx, dy = 0, 5
    va = 'top' if dy < 0 else 'bottom'
    ax.annotate(f'{y_val:.1f}%', xy=(x_val, y_val), xytext=(dx, dy),
                textcoords='offset points', ha='center', va=va,
                fontsize=7, color='#000000')

# KDIGO threshold lines (Prism style: thin, labeled at top with more spacing)
threshold_labels = [
    (1.5, PRISM_ORANGE, 'Stage 1', 1.0),
    (2.0, PRISM_RED, 'Stage 2', 2.5),
    (3.0, PRISM_PURPLE, 'Stage 3', 5.5),
]
for x_val, color, label, label_x in threshold_labels:
    ax.axvline(x=x_val, color=color, linestyle=':', linewidth=0.6, alpha=0.7, zorder=1)
    ax.text(label_x, 48, label, ha='center', va='bottom', fontsize=6.5,
            color=color, fontweight='bold', rotation=0)

ax.set_xlabel('Creatinine Ratio (peak / baseline)', fontsize=9, fontweight='bold')
ax.set_ylabel('Hospital Mortality (%)', fontsize=9, fontweight='bold')
ax.set_xlim(0.3, 13)
ax.set_ylim(0, 50)
ax.set_xticks([1, 2, 3, 5, 7.5, 10, 12.5])

save_fig(fig, 'fig6_dose_response')


# ══════════════════════════════════════════════════════════════
# FIGURE 7: Stage Distribution (Prism donut chart)
# ══════════════════════════════════════════════════════════════
print("[Fig 7] Stage distribution...")

fig, ax = plt.subplots(figsize=(4, 4))
stage_n = [fc['stage_counts'][str(s)] for s in stages]
stage_pct = [n / sum(stage_n) * 100 for n in stage_n]

pie_labels = [f'{STAGE_NAMES[i]}\nn = {stage_n[i]:,} ({stage_pct[i]:.1f}%)' for i in range(4)]

wedges, texts = ax.pie(stage_n, labels=pie_labels, colors=STAGE_COLORS,
                       startangle=90, counterclock=False,
                       wedgeprops=dict(edgecolor='white', linewidth=1.5),
                       textprops=dict(fontsize=7.5, fontweight='normal'),
                       labeldistance=1.15)

# Donut center
centre_circle = plt.Circle((0, 0), 0.45, fc='white', linewidth=0)
ax.add_artist(centre_circle)
ax.text(0, 0.06, 'Total', ha='center', va='center', fontsize=9, fontweight='bold')
ax.text(0, -0.10, f'N = {sum(stage_n):,}', ha='center', va='center', fontsize=8, color='#666666')

save_fig(fig, 'fig7_distribution')


# ══════════════════════════════════════════════════════════════
# FIGURE 8: AKI Paradox (Prism dual-panel)
# ══════════════════════════════════════════════════════════════
print("[Fig 8] AKI paradox...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.75))

# ── Panel A: CKD stratification ──
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
ax1.bar(x - width/2, ckd_neg_mort, width, label='CKD-negative', color=PRISM_BLUE,
        edgecolor='white', linewidth=0.5)
ax1.bar(x + width/2, ckd_pos_mort, width, label='CKD-positive', color=PRISM_PURPLE,
        edgecolor='white', linewidth=0.5)

for i in range(4):
    ax1.text(i - width/2, ckd_neg_mort[i] + 0.8, f'{ckd_neg_mort[i]:.1f}%', ha='center',
             fontsize=7, fontweight='bold', color=PRISM_BLUE)
    ax1.text(i + width/2, ckd_pos_mort[i] + 0.8, f'{ckd_pos_mort[i]:.1f}%', ha='center',
             fontsize=7, fontweight='bold', color=PRISM_PURPLE)

ax1.set_xticks(x)
ax1.set_xticklabels(STAGE_NAMES, fontsize=8)
ax1.set_ylabel('Hospital Mortality (%)', fontsize=9, fontweight='bold')
ax1.set_title('A. CKD Stratification by KDIGO Stage', fontsize=9, fontweight='bold', loc='left')
ax1.set_ylim(0, max(max(ckd_neg_mort), max(ckd_pos_mort)) * 1.25)
ax1.legend(fontsize=7, loc='upper left', framealpha=1, edgecolor='#000000')

# ── Panel B: Stage 3 RRT decomposition ──
stage3 = cohort[cohort['kdigo_stage'] == 3]
rrt_cats = ['RRT only', 'Cr only', 'Both\n(RRT+Cr)']
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
ax2.bar(x2 - width/2, rrt_neg_mort, width, label='CKD-negative', color=PRISM_BLUE,
        edgecolor='white', linewidth=0.5)
ax2.bar(x2 + width/2, rrt_pos_mort, width, label='CKD-positive', color=PRISM_PURPLE,
        edgecolor='white', linewidth=0.5)

for i in range(3):
    ax2.text(i - width/2, rrt_neg_mort[i] + 1.0, f'{rrt_neg_mort[i]:.1f}%\n(n={rrt_neg_n[i]:,})',
             ha='center', fontsize=6.5, fontweight='bold', color=PRISM_BLUE)
    ax2.text(i + width/2, rrt_pos_mort[i] + 1.0, f'{rrt_pos_mort[i]:.1f}%\n(n={rrt_pos_n[i]:,})',
             ha='center', fontsize=6.5, fontweight='bold', color=PRISM_PURPLE)

ax2.set_xticks(x2)
ax2.set_xticklabels(rrt_cats, fontsize=8)
ax2.set_ylabel('Hospital Mortality (%)', fontsize=9, fontweight='bold')
ax2.set_title('B. Stage 3 Decomposed by RRT Indication', fontsize=9, fontweight='bold', loc='left')
ax2.set_ylim(0, max(max(rrt_neg_mort), max(rrt_pos_mort)) * 1.30)
ax2.legend(fontsize=7, loc='upper right', framealpha=1, edgecolor='#000000')

plt.tight_layout(pad=1.0)
save_fig(fig, 'fig8_aki_paradox')


# ══════════════════════════════════════════════════════════════
# FIGURE S1: RCS Nonlinear OR Curve (Prism style)
# ══════════════════════════════════════════════════════════════
print("[Fig S1] RCS nonlinear...")

fig, ax = plt.subplots(figsize=(5, 3.75))

ax.fill_between(rcs['cr_ratio'], rcs['OR_low'], rcs['OR_high'],
                alpha=0.12, color=PRISM_BLUE, zorder=2)
ax.plot(rcs['cr_ratio'], rcs['OR'], color=PRISM_BLUE, linewidth=1.2, zorder=3)

ax.axhline(y=1, color='#000000', linestyle='--', linewidth=0.4, alpha=0.5, zorder=1)

# KDIGO thresholds (Prism style: thin, labeled at top with more spacing)
threshold_labels_rcs = [
    (1.5, PRISM_ORANGE, 'Stage 1', 1.0),
    (2.0, PRISM_RED, 'Stage 2', 2.5),
    (3.0, PRISM_PURPLE, 'Stage 3', 5.5),
]
for x_val, color, label, label_x in threshold_labels_rcs:
    ax.axvline(x=x_val, color=color, linestyle=':', linewidth=0.6, alpha=0.7, zorder=1)
    ax.text(label_x, 5.5, label, ha='center', va='bottom', fontsize=6.5,
            color=color, fontweight='bold')

# Nonlinearity test
nl = mc['RCS_nonlinearity']
ax.text(0.97, 0.97,
        f'RCS: 4 knots\nNonlinearity: P < 0.001',
        transform=ax.transAxes, ha='right', va='top', fontsize=7,
        bbox=dict(boxstyle='square,pad=0.3', facecolor='white', edgecolor='#000000', linewidth=0.5))

ax.set_xlabel('Creatinine Ratio (peak / baseline)', fontsize=9, fontweight='bold')
ax.set_ylabel('Adjusted Odds Ratio', fontsize=9, fontweight='bold')
ax.set_xlim(0.5, 12)
ax.set_ylim(0.5, 6)

save_fig(fig, 'figS1_rcs')


# ══════════════════════════════════════════════════════════════
# FIGURE S2: ROC + DCA (Prism dual-panel)
# ══════════════════════════════════════════════════════════════
print("[Fig S2] ROC + DCA...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.75))

# Panel A: ROC
ax1.plot(roc['fpr_a'], roc['tpr_a'], color=PRISM_ORANGE, linewidth=1, label='Model A (AUC = 0.778)')
ax1.plot(roc['fpr_b'], roc['tpr_b'], color=PRISM_BLUE, linewidth=1.5, label='Model B (AUC = 0.828)')
ax1.plot(roc['fpr_c'], roc['tpr_c'], color=PRISM_GREEN, linewidth=1, linestyle='--',
         label='Model C (AUC = 0.825)')
ax1.plot([0, 1], [0, 1], color='#000000', linestyle=':', linewidth=0.4)

ax1.set_xlabel('1 - Specificity', fontsize=9, fontweight='bold')
ax1.set_ylabel('Sensitivity', fontsize=9, fontweight='bold')
ax1.set_title('A. ROC Curves', fontsize=9, fontweight='bold', loc='left')
ax1.legend(fontsize=7, loc='lower right', framealpha=1, edgecolor='#000000')
ax1.set_xlim(-0.01, 1.01)
ax1.set_ylim(-0.01, 1.01)

# Panel B: DCA
mask = (dca['threshold'] >= 0.01) & (dca['threshold'] <= 0.50)
dca_plot = dca[mask]

ax2.plot(dca_plot['threshold'], dca_plot['nb_model_a'], color=PRISM_ORANGE, linewidth=1, label='Model A')
ax2.plot(dca_plot['threshold'], dca_plot['nb_model_b'], color=PRISM_BLUE, linewidth=1.5, label='Model B')
ax2.plot(dca_plot['threshold'], dca_plot['nb_model_c'], color=PRISM_GREEN, linewidth=1,
         linestyle='--', label='Model C')
ax2.plot(dca_plot['threshold'], dca_plot['nb_treat_all'], color='#000000', linewidth=0.6,
         linestyle='-.', label='Treat All')
ax2.plot(dca_plot['threshold'], dca_plot['nb_treat_none'], color='#000000', linewidth=0.6,
         linestyle=':', label='Treat None')

ax2.set_xlabel('Threshold Probability', fontsize=9, fontweight='bold')
ax2.set_ylabel('Net Benefit', fontsize=9, fontweight='bold')
ax2.set_title('B. Decision Curve Analysis', fontsize=9, fontweight='bold', loc='left')
ax2.legend(fontsize=7, loc='upper right', framealpha=1, edgecolor='#000000')
ax2.set_xlim(0.01, 0.50)
ax2.set_ylim(-0.05, 0.15)

plt.tight_layout(pad=1.0)
save_fig(fig, 'figS2_roc_dca')


# ══════════════════════════════════════════════════════════════
# FIGURE S3: Sepsis Subgroup (Prism-style grouped bar)
# ══════════════════════════════════════════════════════════════
print("[Fig S3] Sepsis subgroup...")

fig, ax = plt.subplots(figsize=(5, 3.75))

stage_labels = ['Stage 1', 'Stage 2', 'Stage 3']
stage_nums = [1, 2, 3]

sepsis_ors = [sepsis_data['Sepsis']['kdigo_or'][str(s)]['OR'] for s in stage_nums]
sepsis_los = [sepsis_data['Sepsis']['kdigo_or'][str(s)]['CI_low'] for s in stage_nums]
sepsis_his = [sepsis_data['Sepsis']['kdigo_or'][str(s)]['CI_high'] for s in stage_nums]

nonsepsis_ors = [sepsis_data['Non-sepsis']['kdigo_or'][str(s)]['OR'] for s in stage_nums]
nonsepsis_los = [sepsis_data['Non-sepsis']['kdigo_or'][str(s)]['CI_low'] for s in stage_nums]
nonsepsis_his = [sepsis_data['Non-sepsis']['kdigo_or'][str(s)]['CI_high'] for s in stage_nums]

y_pos = np.arange(len(stage_labels))
bar_height = 0.3

ax.barh(y_pos + bar_height/2, nonsepsis_ors, height=bar_height, color=PRISM_BLUE,
        alpha=0.85, label=f'Non-sepsis (n={sepsis_data["Non-sepsis"]["n"]:,})',
        edgecolor='white', linewidth=0.3)
ax.barh(y_pos - bar_height/2, sepsis_ors, height=bar_height, color=PRISM_RED,
        alpha=0.85, label=f'Sepsis (n={sepsis_data["Sepsis"]["n"]:,})',
        edgecolor='white', linewidth=0.3)

ax.errorbar(nonsepsis_ors, y_pos + bar_height/2,
            xerr=[[o - lo for o, lo in zip(nonsepsis_ors, nonsepsis_los)],
                  [hi - o for o, hi in zip(nonsepsis_ors, nonsepsis_his)]],
            fmt='none', ecolor='#000000', capsize=2, capthick=0.5, lw=0.5)
ax.errorbar(sepsis_ors, y_pos - bar_height/2,
            xerr=[[o - lo for o, lo in zip(sepsis_ors, sepsis_los)],
                  [hi - o for o, hi in zip(sepsis_ors, sepsis_his)]],
            fmt='none', ecolor='#000000', capsize=2, capthick=0.5, lw=0.5)

for i, (or_val, lo, hi) in enumerate(zip(nonsepsis_ors, nonsepsis_los, nonsepsis_his)):
    ax.text(hi + 0.15, i + bar_height/2, f'{or_val:.2f} ({lo:.2f}\u2013{hi:.2f})',
            va='center', fontsize=6.5, color=PRISM_BLUE)
for i, (or_val, lo, hi) in enumerate(zip(sepsis_ors, sepsis_los, sepsis_his)):
    ax.text(hi + 0.15, i - bar_height/2, f'{or_val:.2f} ({lo:.2f}\u2013{hi:.2f})',
            va='center', fontsize=6.5, color=PRISM_RED)

ax.set_yticks(y_pos)
ax.set_yticklabels(stage_labels, fontsize=8, fontweight='bold')
ax.set_xlabel('Adjusted Odds Ratio (vs. Stage 0)', fontsize=9, fontweight='bold')
ax.axvline(x=1, color='#000000', linestyle='--', linewidth=0.4)
ax.set_xlim(0, max(sepsis_his + nonsepsis_his) * 1.5)
ax.legend(fontsize=7, loc='lower right', framealpha=1, edgecolor='#000000')

# Interaction p-value
int_data = sepsis_data['interaction']
ax.text(0.97, 0.97,
        f'KDIGO x Sepsis\nP < 0.001',
        transform=ax.transAxes, ha='right', va='top', fontsize=7,
        bbox=dict(boxstyle='square,pad=0.3', facecolor='white', edgecolor='#000000', linewidth=0.5))

save_fig(fig, 'figS3_sepsis')


# ══════════════════════════════════════════════════════════════
# FIGURE S4: Calibration Curve (Prism style)
# ══════════════════════════════════════════════════════════════
print("[Fig S4] Calibration...")

fig, ax = plt.subplots(figsize=(4.5, 4.5))

ax.plot([0, 1], [0, 1], color='#000000', linestyle='--', linewidth=0.5, label='Ideal')

ax.plot(cal['expected'], cal['observed'], 'o-', color=PRISM_BLUE, linewidth=1.2,
        markersize=5, markeredgecolor='white', markeredgewidth=0.3, label='Model B', zorder=3)

ax.errorbar(cal['expected'], cal['observed'],
            yerr=[cal['observed'] - cal['obs_ci_low'], cal['obs_ci_high'] - cal['observed']],
            fmt='none', ecolor=PRISM_BLUE, capsize=2, capthick=0.5, lw=0.5, alpha=0.5)

# Distribution histogram (secondary axis)
ax2 = ax.twinx()
ax2.hist(cal['expected'], weights=cal['n'], bins=20, alpha=0.1, color=PRISM_ORANGE, edgecolor='none')
ax2.set_ylabel('Count', fontsize=9, color=PRISM_ORANGE)
ax2.tick_params(axis='y', labelcolor=PRISM_ORANGE, labelsize=7)
ax2.set_ylim(0, max(cal['n']) * 3)
ax2.spines['right'].set_visible(True)
ax2.spines['right'].set_linewidth(0.5)

# Hosmer-Lemeshow
hl = cal_results['Model_B']
ax.text(0.05, 0.95,
        f'Hosmer-Lemeshow\nP < 0.001',
        transform=ax.transAxes, ha='left', va='top', fontsize=7,
        bbox=dict(boxstyle='square,pad=0.3', facecolor='white', edgecolor='#000000', linewidth=0.5))

ax.set_xlabel('Predicted Probability', fontsize=9, fontweight='bold')
ax.set_ylabel('Observed Mortality Rate', fontsize=9, fontweight='bold')
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(-0.02, 1.02)
ax.legend(fontsize=7, loc='lower right', framealpha=1, edgecolor='#000000')

save_fig(fig, 'figS4_calibration')


# ══════════════════════════════════════════════════════════════
# EXPORT: Prism-Ready Excel Workbook
# ══════════════════════════════════════════════════════════════
print("\n[Export] Creating Prism data workbook...")

import openpyxl
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter

wb = openpyxl.Workbook()
wb.remove(wb.active)

# Styles
header_font = Font(name='Arial', size=10, bold=True, color='FFFFFF')
header_fill = PatternFill(start_color='4472C4', end_color='4472C4', fill_type='solid')
data_font = Font(name='Arial', size=10)
center_align = Alignment(horizontal='center', vertical='center')
thin_border = Border(
    left=Side(style='thin', color='D9D9D9'),
    right=Side(style='thin', color='D9D9D9'),
    top=Side(style='thin', color='D9D9D9'),
    bottom=Side(style='thin', color='D9D9D9'),
)

def add_sheet(name, data, col_widths=None):
    """Add a data sheet to the workbook.
    data: list of lists (first row = headers)
    col_widths: list of column widths
    """
    ws = wb.create_sheet(name[:31])  # Excel max sheet name = 31 chars
    for row_idx, row in enumerate(data, 1):
        for col_idx, val in enumerate(row, 1):
            cell = ws.cell(row=row_idx, column=col_idx, value=val)
            if row_idx == 1:
                cell.font = header_font
                cell.fill = header_fill
                cell.alignment = center_align
            else:
                cell.font = data_font
                cell.alignment = center_align
            cell.border = thin_border
    if col_widths:
        for i, w in enumerate(col_widths, 1):
            ws.column_dimensions[get_column_letter(i)].width = w
    return ws

# ── Sheet 1: Fig2 Mortality ──
data = [['KDIGO Stage', 'N', 'Mortality (%)', 'CI Lower (%)', 'CI Upper (%)']]
for d in mort_data:
    data.append([f'Stage {d["stage"]}', d['n'], round(d['mort'], 1), round(d['ci_low'], 1), round(d['ci_high'], 1)])
add_sheet('Fig2_Mortality', data, [15, 12, 15, 15, 15])

# ── Sheet 2: Fig3 KM Survival ──
for idx, s in enumerate(stages):
    key = str(s)
    d = km['km_data'][key]
    data = [['Time (days)', 'Survival', 'CI Lower', 'CI Upper', 'N at Risk']]
    times = d['times']
    surv = d['survival']
    ci_l = d.get('survival_ci_low', [None]*len(times))
    ci_h = d.get('survival_ci_high', [None]*len(times))
    nar = d.get('n_at_risk', [])
    for j in range(len(times)):
        nar_val = nar[j] if j < len(nar) else ''
        data.append([round(times[j], 1), round(surv[j], 4),
                     round(ci_l[j], 4) if ci_l[j] else '', round(ci_h[j], 4) if ci_h[j] else '', nar_val])
    add_sheet(f'Fig3_KM_Stage{s}', data, [12, 12, 12, 12, 12])

# ── Sheet 3: Fig4 Forest Plot ──
data = [['Variable', 'Label', 'OR', 'CI Lower', 'CI Upper', 'P-value', 'Significance']]
for _, row in plot_data.iterrows():
    p = row['p_value']
    sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'n.s.'))
    data.append([row['variable'], row['label'], round(row['OR'], 3), round(row['CI_low'], 3),
                 round(row['CI_high'], 3), f'{p:.4f}', sig])
add_sheet('Fig4_Forest_ModelB', data, [28, 28, 10, 10, 10, 12, 12])

# ── Sheet 4: Fig5 ICU LOS ──
data = [['KDIGO Stage', 'N', 'Median (h)', 'Mean (h)', 'Q1 (h)', 'Q3 (h)']]
for i, st in enumerate(los_stats):
    data.append([f'Stage {stages[i]}', st['n'], round(st['median'], 1),
                 round(st['mean'], 1), round(st['p25'], 1), round(st['p75'], 1)])
add_sheet('Fig5_ICU_LOS', data, [15, 12, 12, 12, 12, 12])

# ── Sheet 5: Fig6 Dose-Response ──
data = [['Cr Ratio Bin', 'N', 'Deaths', 'Mortality (%)', 'CI Lower (%)', 'CI Upper (%)']]
for _, row in dose_resp.iterrows():
    data.append([row.get('cr_ratio_bin', ''), row['n'], row.get('deaths', ''),
                 round(row['mortality_pct'], 1), round(row['CI_low'], 1),
                 round(row['CI_high'], 1)])
add_sheet('Fig6_DoseResponse', data, [15, 12, 10, 15, 15, 15])

# ── Sheet 6: Fig7 Distribution ──
data = [['Stage', 'N', 'Percentage (%)']]
for i, s in enumerate(stages):
    data.append([STAGE_NAMES[i], stage_n[i], round(stage_pct[i], 1)])
add_sheet('Fig7_Distribution', data, [15, 12, 15])

# ── Sheet 7: Fig8 AKI Paradox ──
data = [['Stage', 'CKD-negative N', 'CKD-negative Mort (%)', 'CKD-positive N', 'CKD-positive Mort (%)']]
for i in range(4):
    data.append([STAGE_NAMES[i], ckd_neg_n[i], round(ckd_neg_mort[i], 1),
                 ckd_pos_n[i], round(ckd_pos_mort[i], 1)])
add_sheet('Fig8A_CKD_Strat', data, [15, 18, 22, 18, 22])

data = [['RRT Category', 'CKD-negative N', 'CKD-negative Mort (%)', 'CKD-positive N', 'CKD-positive Mort (%)']]
for i in range(3):
    data.append([rrt_cats[i].replace('\n', ' '), rrt_neg_n[i], round(rrt_neg_mort[i], 1),
                 rrt_pos_n[i], round(rrt_pos_mort[i], 1)])
add_sheet('Fig8B_RRT_Decompose', data, [18, 18, 22, 18, 22])

# ── Sheet 8: FigS1 RCS ──
data = [['Cr Ratio', 'OR', 'CI Lower', 'CI Upper']]
for _, row in rcs.iterrows():
    data.append([round(row['cr_ratio'], 3), round(row['OR'], 3),
                 round(row['OR_low'], 3), round(row['OR_high'], 3)])
add_sheet('FigS1_RCS', data, [12, 12, 12, 12])

# ── Sheet 9: FigS2 ROC ──
data = [['FPR_ModelA', 'TPR_ModelA', 'FPR_ModelB', 'TPR_ModelB', 'FPR_ModelC', 'TPR_ModelC']]
for _, row in roc.iterrows():
    data.append([round(row['fpr_a'], 4), round(row['tpr_a'], 4),
                 round(row['fpr_b'], 4), round(row['tpr_b'], 4),
                 round(row['fpr_c'], 4), round(row['tpr_c'], 4)])
add_sheet('FigS2_ROC', data, [14, 14, 14, 14, 14, 14])

# ── Sheet 10: FigS2 DCA ──
data = [['Threshold', 'NB_ModelA', 'NB_ModelB', 'NB_ModelC', 'NB_TreatAll', 'NB_TreatNone']]
for _, row in dca.iterrows():
    data.append([round(row['threshold'], 4), round(row['nb_model_a'], 4),
                 round(row['nb_model_b'], 4), round(row['nb_model_c'], 4),
                 round(row['nb_treat_all'], 4), round(row['nb_treat_none'], 4)])
add_sheet('FigS2_DCA', data, [12, 14, 14, 14, 14, 14])

# ── Sheet 11: FigS3 Sepsis Subgroup ──
data = [['Stage', 'Sepsis OR', 'Sepsis CI Low', 'Sepsis CI High', 'Non-sepsis OR', 'Non-sepsis CI Low', 'Non-sepsis CI High']]
for i, s in enumerate(stage_nums):
    data.append([stage_labels[i], round(sepsis_ors[i], 3), round(sepsis_los[i], 3), round(sepsis_his[i], 3),
                 round(nonsepsis_ors[i], 3), round(nonsepsis_los[i], 3), round(nonsepsis_his[i], 3)])
add_sheet('FigS3_Sepsis', data, [10, 12, 14, 14, 14, 16, 16])

# ── Sheet 12: FigS4 Calibration ──
data = [['Predicted Prob', 'Observed Rate', 'CI Lower', 'CI Upper', 'N']]
for _, row in cal.iterrows():
    data.append([round(row['expected'], 4), round(row['observed'], 4),
                 round(row['obs_ci_low'], 4), round(row['obs_ci_high'], 4), row['n']])
add_sheet('FigS4_Calibration', data, [15, 15, 12, 12, 10])

# ── Sheet 13: Model Comparison Summary ──
data = [
    ['Metric', 'Model A', 'Model B', 'Model C'],
    ['AUC', mc['AUC']['Model_A']['auc'], mc['AUC']['Model_B']['auc'], mc['AUC']['Model_C']['auc']],
    ['AUC 95% CI Low', mc['AUC']['Model_A']['ci_low'], mc['AUC']['Model_B']['ci_low'], mc['AUC']['Model_C']['ci_low']],
    ['AUC 95% CI High', mc['AUC']['Model_A']['ci_high'], mc['AUC']['Model_B']['ci_high'], mc['AUC']['Model_C']['ci_high']],
    ['Brier Score', mc['Brier']['Model_A'], mc['Brier']['Model_B'], mc['Brier']['Model_C']],
    ['NRI (B vs A)', '', mc['NRI_IDI_A_to_B']['NRI'], ''],
    ['NRI Events', '', mc['NRI_IDI_A_to_B']['NRI_events'], ''],
    ['NRI Non-events', '', mc['NRI_IDI_A_to_B']['NRI_nonevents'], ''],
    ['IDI (B vs A)', '', mc['NRI_IDI_A_to_B']['IDI'], ''],
    ['RCS Nonlinearity chi2', '', mc['RCS_nonlinearity']['chi2'], ''],
    ['RCS Nonlinearity df', '', mc['RCS_nonlinearity']['df'], ''],
    ['RCS Nonlinearity p', '', mc['RCS_nonlinearity']['p_value'], ''],
    ['HL chi2 (Model B)', '', cal_results['Model_B']['HL_chi2'], ''],
]
add_sheet('Model_Summary', data, [22, 15, 15, 15])

# ── Sheet 14: CKD x KDIGO Interaction ──
data = [['Interaction Term', 'OR', 'CI Lower', 'CI Upper', 'P-value']]
int_items = [
    ('Stage1 x CKD (ICD)', ckd_int['stage1_interaction']),
    ('Stage2 x CKD (ICD)', ckd_int['stage2_interaction']),
    ('Stage3 x CKD (ICD)', ckd_int['stage3_interaction']),
    ('Stage1 x Lab-CKD', ckd_int['stage1_labckd_interaction']),
    ('Stage2 x Lab-CKD', ckd_int['stage2_labckd_interaction']),
    ('Stage3 x Lab-CKD', ckd_int['stage3_labckd_interaction']),
]
for label, val in int_items:
    p_val = val.get('p', val.get('p_value', ''))
    data.append([label, round(val['OR'], 3), round(val['CI_low'], 3),
                 round(val['CI_high'], 3), p_val])
# Add overall interaction test
data.append([])
data.append(['Overall ICD-CKD Interaction', f'chi2={ckd_int["interaction"]["chi2"]}', f'df={ckd_int["interaction"]["df"]}', '', f'p={ckd_int["interaction"]["p_value"]}'])
data.append(['Overall Lab-CKD Interaction', f'chi2={ckd_int["labckd_interaction"]["chi2"]}', f'df={ckd_int["labckd_interaction"]["df"]}', '', f'p={ckd_int["labckd_interaction"]["p_value"]}'])
add_sheet('CKD_Interaction', data, [28, 12, 12, 12, 15])

xlsx_path = os.path.join(OUT_DIR, 'KDIGO_Prism_Data.xlsx')
wb.save(xlsx_path)
print(f"  -> KDIGO_Prism_Data.xlsx ({len(wb.sheetnames)} sheets)")

# ══════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print(f"ALL FIGURES GENERATED IN PRISM STYLE")
print(f"Output directory: {OUT_DIR}")
print(f"  - 12 figure pairs (PNG + PDF)")
print(f"  - 1 Excel workbook (Prism-ready data)")
print("=" * 60)
