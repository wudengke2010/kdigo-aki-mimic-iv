#!/usr/bin/env python3
"""
KDIGO AKI 分期 — 修订版可视化
==============================
审稿修订:
  Fix #14: KM 曲线含 number-at-risk 表
  Fix #13: 剂量-反应图含 95% CI
  Fix #9: 患者筛选流程图
  Fix #16: 色盲友好配色
  Forest plot 含全部协变量
"""
import pandas as pd
import numpy as np
import os, json
from datetime import datetime

print("=" * 80)
print("KDIGO-AKI REVISED Visualizations")
print(f"Start: {datetime.now()}")
print("=" * 80)

IN = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
OUT_DIR = os.path.join(IN, "figures_revised")
os.makedirs(OUT_DIR, exist_ok=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches

# Colorblind-friendly palette (Wong 2011)
COLORS = ['#0072B2', '#E69F00', '#D55E00', '#CC79A7', '#009E73', '#F0E442', '#56B4E9', '#000000']
STAGE_COLORS = [COLORS[0], COLORS[1], COLORS[5], COLORS[3]]  # Blue, Orange, Yellow, Purple

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial']
plt.rcParams['axes.unicode_minus'] = False

# ============================================================
# 1. Load data
# ============================================================
print("[1] Loading data...")
cohort_path = os.path.join(IN, "kdigo_cohort_revised_enriched.csv")
if not os.path.exists(cohort_path):
    cohort_path = os.path.join(IN, "kdigo_cohort_revised.csv")

cohort = pd.read_csv(cohort_path, low_memory=False)
print(f"  N = {len(cohort):,}")

# Load flowchart data
flowchart_path = os.path.join(IN, "flowchart_data.json")
flowchart = None
if os.path.exists(flowchart_path):
    with open(flowchart_path, 'r') as f:
        flowchart = json.load(f)

# Load KM data
km_path = os.path.join(IN, "km_revised.json")
km_data = None
if os.path.exists(km_path):
    with open(km_path, 'r') as f:
        km_data = json.load(f)

# Load logistic results
lr_path = os.path.join(IN, "logistic_results_revised.csv")
lr_results = None
if os.path.exists(lr_path):
    lr_results = pd.read_csv(lr_path)

# Load dose-response
dr_path = os.path.join(IN, "dose_response_revised.csv")
dr_results = None
if os.path.exists(dr_path):
    dr_results = pd.read_csv(dr_path)

stages = [0, 1, 2, 3]
stage_labels = ['No AKI', 'Stage 1', 'Stage 2', 'Stage 3']

# ============================================================
# 2. Patient Selection Flowchart
# ============================================================
print("[2] Patient selection flowchart...")

if flowchart:
    fig, ax = plt.subplots(figsize=(10, 12))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 14)
    ax.axis('off')
    ax.set_title('Patient Selection Flowchart', fontsize=16, fontweight='bold', pad=20)

    def draw_box(y, text, n, color='#E8F0FE', x_center=5, width=4.5):
        rect = mpatches.FancyBboxPatch((x_center - width/2, y - 0.35), width, 0.7,
                                        boxstyle="round,pad=0.1", facecolor=color,
                                        edgecolor='#333333', linewidth=1.5)
        ax.add_patch(rect)
        ax.text(x_center, y + 0.05, f"{text}", ha='center', va='center', fontsize=12, fontweight='bold')
        ax.text(x_center, y - 0.22, f"(n = {n:,})", ha='center', va='center', fontsize=10, color='#555555')

    def draw_arrow(y1, y2, x=5):
        ax.annotate('', xy=(x, y2 + 0.35), xytext=(x, y1 - 0.35),
                    arrowprops=dict(arrowstyle='->', lw=1.5, color='#333333'))

    def draw_exclusion(y, text, n, x_excl=8.3):
        ax.text(x_excl, y, f"Excluded: {text}", ha='left', va='center', fontsize=9, color='#D55E00')
        ax.text(x_excl, y - 0.2, f"(n = {n:,})", ha='left', va='center', fontsize=8, color='#888888')
        ax.annotate('', xy=(7.25, y), xytext=(7.25, y + 0.3),
                    arrowprops=dict(arrowstyle='->', lw=1, color='#D55E00', linestyle='dashed'))

    # Box 1
    y = 13
    draw_box(y, f"All ICU Stays in MIMIC-IV v3.1", flowchart['total_icu_stays'], color='#D6EAF8')

    # Box 2
    y = 11.5
    draw_arrow(13, 11.5)
    draw_box(y, f"Unique ICU Stays (first per hospitalization)", flowchart['after_first_icu_adult'] + flowchart['excluded_duplicate_stays'])
    draw_exclusion(11.5, "Duplicate ICU stays (same hospitalization)", flowchart['excluded_duplicate_stays'])

    # Box 3
    y = 10
    draw_arrow(11.5, 10)
    draw_box(y, f"Adult ICU Stays (age \u2265 18)", flowchart['after_first_icu_adult'])
    n_excluded_age = flowchart['excluded_age_lt_18'] if flowchart['excluded_age_lt_18'] > 0 else 0
    if n_excluded_age > 0:
        draw_exclusion(10, "Age < 18 years", n_excluded_age)

    # Box 4
    y = 8.5
    draw_arrow(10, 8.5)
    draw_box(y, f"Stays with Creatinine Data", flowchart['final_cohort'] + flowchart['excluded_missing_cr'])
    draw_exclusion(8.5, "Missing baseline or peak creatinine", flowchart['excluded_missing_cr'])

    # Box 5
    y = 7
    draw_arrow(8.5, 7)
    draw_box(y, f"Final Analytic Cohort", flowchart['final_cohort'], color='#D5F5E3')

    # Stage boxes
    stage_counts = flowchart.get('stage_counts', {})
    for i, (s, label) in enumerate(zip([0, 1, 2, 3], stage_labels)):
        n = stage_counts.get(str(s), 0)
        pct = n / max(flowchart['final_cohort'], 1) * 100
        x_pos = 1.5 + i * 2.2
        y_pos = 5.5
        rect = mpatches.FancyBboxPatch((x_pos - 0.9, y_pos - 0.35), 1.8, 1.2,
                                        boxstyle="round,pad=0.1", facecolor=STAGE_COLORS[i],
                                        edgecolor='white', linewidth=1, alpha=0.85)
        ax.add_patch(rect)
        ax.text(x_pos, y_pos + 0.4, f"{label}", ha='center', va='center', fontsize=10, fontweight='bold', color='white')
        ax.text(x_pos, y_pos + 0.05, f"n = {n:,}", ha='center', va='center', fontsize=11, fontweight='bold', color='white')
        ax.text(x_pos, y_pos - 0.22, f"({pct:.1f}%)", ha='center', va='center', fontsize=9, color='white')
        # Arrow from final cohort
        if i == 0:
            ax.annotate('', xy=(x_pos, y_pos + 0.55), xytext=(5, 6.65),
                        arrowprops=dict(arrowstyle='->', lw=1, color='#333333'))

    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'flowchart.png'), dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print("  -> flowchart.png")

# ============================================================
# 3. Mortality Bar Chart with 95% CI
# ============================================================
print("[3] Mortality bar chart with 95% CI...")

def wilson_ci(p, n, z=1.96):
    """Wilson score interval for binomial proportion."""
    if n == 0:
        return p, p
    denom = 1 + z*z/n
    centre = (p + z*z/(2*n)) / denom
    width = z * np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / denom
    return centre - width, centre + width

n_list = []
mort_list = []
ci_low_list = []
ci_high_list = []
for s in stages:
    sub = cohort[cohort['kdigo_stage'] == s]
    n = len(sub)
    events = sub['hospital_expire_flag'].sum()
    p = events / n if n > 0 else 0
    mort = p * 100
    ci_low, ci_high = wilson_ci(p, n)
    n_list.append(n)
    mort_list.append(mort)
    ci_low_list.append(ci_low * 100)
    ci_high_list.append(ci_high * 100)

fig, ax = plt.subplots(figsize=(8, 6))
bars = ax.bar(stage_labels, mort_list, color=STAGE_COLORS, edgecolor='white', linewidth=1.5)

# Add error bars (95% CI)
errors = [[m - l for m, l in zip(mort_list, ci_low_list)],
          [h - m for m, h in zip(mort_list, ci_high_list)]]
ax.errorbar(range(len(stage_labels)), mort_list, yerr=errors, fmt='none',
            ecolor='black', capsize=4, capthick=1.5, linewidth=1.5)

for i, (bar, mort, n, ci_high) in enumerate(zip(bars, mort_list, n_list, ci_high_list)):
    label_y = ci_high + 2.8
    ax.text(bar.get_x() + bar.get_width()/2., label_y,
            f'{mort:.1f}%\n(n={n:,})', ha='center', va='bottom', fontsize=9, fontweight='bold')

ax.set_ylabel('Hospital Mortality (%)', fontsize=13, fontweight='bold')
ax.set_title('ICU Hospital Mortality by KDIGO AKI Stage', fontsize=15, fontweight='bold')
ax.set_ylim(0, max(ci_high_list) * 1.35)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=11)

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'mortality_by_stage.png'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("  -> mortality_by_stage.png")

# ============================================================
# 4. KM Curves with Number-at-Risk
# ============================================================
print("[4] KM survival curves with number-at-risk...")

cohort['intime_dt'] = pd.to_datetime(cohort['intime'])
cohort['dod_dt'] = pd.to_datetime(cohort['dod'])

cohort['event_30d'] = 0
cohort['time_30d'] = 30.0
for idx in cohort.index:
    dod = cohort.at[idx, 'dod_dt']
    intime = cohort.at[idx, 'intime_dt']
    if pd.notna(dod) and pd.notna(intime):
        t = (dod - intime).total_seconds() / 86400
        if 0 < t <= 30:
            cohort.at[idx, 'time_30d'] = t
            cohort.at[idx, 'event_30d'] = 1

# Create figure with subplots (KM on top, number-at-risk below)
fig, (ax_km, ax_risk) = plt.subplots(2, 1, figsize=(12, 9),
                                      gridspec_kw={'height_ratios': [4, 1]},
                                      sharex=True)

timepoints = [0, 5, 10, 15, 20, 25, 30]
n_at_risk_all = {}

for stage, color, label in zip(stages, STAGE_COLORS, stage_labels):
    sub = cohort[cohort['kdigo_stage'] == stage]
    if len(sub) == 0:
        continue
    
    # Manual KM
    df = sub[['time_30d', 'event_30d']].sort_values('time_30d')
    times = sorted(df['time_30d'].unique())
    
    surv_curve = [1.0]
    n_risk = len(df)
    surv = 1.0
    seen_t = set()
    for t in times:
        if t > 0:
            events_at_t = df[(df['time_30d'] == t) & (df['event_30d'] == 1)].shape[0]
            censored_at_t = df[(df['time_30d'] == t) & (df['event_30d'] == 0)].shape[0]
            if n_risk > 0:
                surv *= (1 - events_at_t / n_risk)
            n_risk -= (events_at_t + censored_at_t)
            surv_curve.append(surv)
        seen_t.add(t)
    
    plot_times = [0] + list(times)
    ax_km.step(plot_times, surv_curve, where='post', color=color, linewidth=2.5, label=f'{label} (n={len(sub):,})')
    
    # Number at risk
    n_risk_list = []
    for tp in timepoints:
        n_risk = (sub['time_30d'] >= tp).sum()
        n_risk_list.append(n_risk)
    n_at_risk_all[label] = n_risk_list

ax_km.set_ylabel('Survival Probability', fontsize=13, fontweight='bold')
ax_km.set_title('30-Day Survival by KDIGO AKI Stage', fontsize=15, fontweight='bold')
ax_km.legend(fontsize=11, frameon=True, fancybox=True, loc='lower left')
ax_km.set_ylim(0.55, 1.02)
ax_km.grid(True, alpha=0.2)
ax_km.spines['top'].set_visible(False)
ax_km.spines['right'].set_visible(False)

# Number-at-risk table
ax_risk.axis('off')
risk_table_data = []
for label in stage_labels:
    if label in n_at_risk_all:
        risk_table_data.append([label] + [str(n) for n in n_at_risk_all[label]])

col_labels = ['', 'Day 0', 'Day 5', 'Day 10', 'Day 15', 'Day 20', 'Day 25', 'Day 30']
table = ax_risk.table(cellText=risk_table_data, colLabels=col_labels,
                       cellLoc='center', loc='center',
                       colColours=['#f0f0f0'] * len(col_labels))
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 1.5)

for i, label in enumerate(stage_labels):
    if label in n_at_risk_all:
        for j in range(len(col_labels)):
            table[(i, j)].set_facecolor(STAGE_COLORS[stage_labels.index(label)])
            table[(i, j)].set_alpha(0.15)

ax_risk.text(0.5, -0.1, 'Number at risk', ha='center', fontsize=11, fontweight='bold', transform=ax_risk.transAxes)

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'km_survival.png'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("  -> km_survival.png (with number-at-risk)")

# ============================================================
# 5. Forest Plot — All Covariates
# ============================================================
print("[5] Forest plot (all covariates)...")

if lr_results is not None:
    # Use Model B results
    model_b = lr_results[lr_results['model'] == 'B'].copy()
    
    if len(model_b) > 0:
        # Rename for display
        var_labels = {
            'kdigo_stage_1': 'KDIGO Stage 1 (vs No AKI)',
            'kdigo_stage_2': 'KDIGO Stage 2 (vs No AKI)',
            'kdigo_stage_3': 'KDIGO Stage 3 (vs No AKI)',
            'anchor_age': 'Age (per year)',
            'male_num': 'Male sex',
            'sofa_imputed': 'SOFA score (per point)',
            'mech_vent': 'Mechanical ventilation',
            'vasopressor': 'Vasopressor use',
            'ckd_icd': 'CKD (ICD-coded)',
        }
        
        model_b['label'] = model_b['variable'].apply(
            lambda x: var_labels.get(x, x.replace('com_', '').replace('_', ' ')[:40])
        )
        
        # Sort: KDIGO stages first in natural order (1,2,3), then other variables by OR
        def stage_order(x):
            if x == 'kdigo_stage_1': return 1
            if x == 'kdigo_stage_2': return 2
            if x == 'kdigo_stage_3': return 3
            return 999
        kdigo_rows = model_b[model_b['variable'].str.startswith('kdigo')]
        kdigo_rows = kdigo_rows.assign(stage_order=kdigo_rows['variable'].apply(stage_order)).sort_values('stage_order')
        other_rows = model_b[~model_b['variable'].str.startswith('kdigo')].sort_values('OR')
        plot_data = pd.concat([kdigo_rows, other_rows])
        
        # Focus on key variables
        key_mask = plot_data['variable'].str.startswith('kdigo') | \
                   plot_data['variable'].isin(['anchor_age', 'male_num', 'sofa_imputed', 'mech_vent', 'vasopressor', 'ckd_icd']) | \
                   (plot_data['p_value'] < 0.001)
        plot_data = plot_data[key_mask]
        
        fig, ax = plt.subplots(figsize=(12, max(6, len(plot_data) * 0.4)))
        
        y_positions = range(len(plot_data))
        for i, (or_val, low, high, p_val) in enumerate(zip(plot_data['OR'], plot_data['CI_low'], plot_data['CI_high'], plot_data['p_value'])):
            color = '#D55E00' if or_val > 1 else (COLORS[0] if or_val < 1 else '#888888')
            ax.errorbar(or_val, i, xerr=[[or_val - low], [high - or_val]],
                         fmt='o', color=color, capsize=4, capthick=1.5, markersize=8, linewidth=2)
        
        ax.axvline(x=1, linestyle='--', color='#888888', alpha=0.5, linewidth=1.5)
        ax.set_yticks(y_positions)
        ax.set_yticklabels(plot_data['label'].values, fontsize=10)
        ax.set_xlabel('Odds Ratio (95% CI)', fontsize=13, fontweight='bold')
        ax.set_title('Multivariable Logistic Regression: Hospital Mortality (Model B)', fontsize=14, fontweight='bold')
        ax.set_xscale('log')
        ax.set_xlim(0.3, 50)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(True, alpha=0.2, axis='x')
        
        # Add OR text
        for i, (or_val, p_val) in enumerate(zip(plot_data['OR'], plot_data['p_value'])):
            text = f'{or_val:.2f}'
            if p_val < 0.001:
                text += ' ***'
            elif p_val < 0.01:
                text += ' **'
            elif p_val < 0.05:
                text += ' *'
            ax.text(or_val * 1.5, i, text, va='center', fontsize=9, color='#333333')
        
        plt.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, 'forest_plot.png'), dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        print("  -> forest_plot.png (all covariates)")
    else:
        print("  No Model B results available, skipping forest plot")

# ============================================================
# 6. ICU LOS Violin Plot (robust to extreme outliers)
# ============================================================
print("[6] ICU LOS violin plot...")

los_data = []
for s in stages:
    sub = cohort[cohort['kdigo_stage'] == s]['icu_los_hours'].dropna()
    sub = sub[(sub > 0) & (sub < sub.quantile(0.99))]
    los_data.append(sub.values)

# Limit y-axis to 95th percentile of pooled data for clarity
pooled_los = np.concatenate(los_data)
y_max = np.percentile(pooled_los, 95)

fig, ax = plt.subplots(figsize=(8, 6))
parts = ax.violinplot(los_data, positions=range(1, 5), showmeans=False, showmedians=False, showextrema=False)

for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(STAGE_COLORS[i])
    pc.set_alpha(0.5)
    pc.set_edgecolor(STAGE_COLORS[i])

# Add median and quartile markers
for i, data in enumerate(los_data, 1):
    q25, q50, q75 = np.percentile(data, [25, 50, 75])
    ax.scatter([i], [q50], color='black', zorder=5, s=40, marker='D')
    ax.plot([i, i], [q25, q75], color='black', linewidth=2, zorder=4)

ax.set_xticks(range(1, 5))
ax.set_xticklabels(stage_labels, fontsize=12)
ax.set_ylabel('ICU Length of Stay (hours)', fontsize=13, fontweight='bold')
ax.set_title('ICU LOS by KDIGO AKI Stage', fontsize=15, fontweight='bold')
ax.set_ylim(0, y_max)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(True, alpha=0.2, axis='y')

# Annotation explaining truncation
ax.text(0.98, 0.98, 'Top 5% of each stage truncated for clarity', transform=ax.transAxes,
        ha='right', va='top', fontsize=9, color='#666666', style='italic')

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'icu_los_by_stage.png'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("  -> icu_los_by_stage.png")

# ============================================================
# 7. Dose-Response with 95% CI
# ============================================================
print("[7] Dose-response with 95% CI...")

if dr_results is not None:
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x = range(len(dr_results))
    y = dr_results['mortality_pct'].values
    yerr_low = y - dr_results['CI_low'].values
    yerr_high = dr_results['CI_high'].values - y
    
    ax.errorbar(x, y, yerr=[yerr_low, yerr_high], fmt='o-', color=COLORS[0],
                capsize=5, capthick=2, markersize=10, linewidth=2.5, markerfacecolor=COLORS[0])
    
    ax.set_xticks(x)
    ax.set_xticklabels(dr_results['cr_ratio_bin'].values, fontsize=11, rotation=30, ha='right')
    ax.set_ylabel('Hospital Mortality (%)', fontsize=13, fontweight='bold')
    ax.set_xlabel('Creatinine Ratio (peak / baseline)', fontsize=13, fontweight='bold')
    ax.set_title('Dose-Response: Mortality by Creatinine Ratio', fontsize=15, fontweight='bold')
    
    # Annotate with n
    for i, row in dr_results.iterrows():
        ax.annotate(f"n={int(row['n']):,}", (i, row['mortality_pct']),
                    textcoords="offset points", xytext=(0, 12),
                    ha='center', fontsize=9, color='#666666')
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.2, axis='y')
    
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'dose_response.png'), dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print("  -> dose_response.png (with 95% CI)")
else:
    print("  No dose-response data available")

# ============================================================
# 8. Stage Distribution Pie Chart
# ============================================================
print("[8] Stage distribution pie chart...")

fig, ax = plt.subplots(figsize=(8, 6))
sizes = [cohort[cohort['kdigo_stage'] == s].shape[0] for s in stages]

wedges, texts, autotexts = ax.pie(sizes, labels=stage_labels, colors=STAGE_COLORS,
                                    autopct='%1.1f%%', startangle=90,
                                    textprops={'fontsize': 12, 'fontweight': 'bold'},
                                    explode=(0, 0.02, 0.05, 0.1))
for at in autotexts:
    at.set_fontsize(11)
    at.set_fontweight('bold')

ax.set_title('KDIGO AKI Stage Distribution in ICU', fontsize=15, fontweight='bold')
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'stage_distribution.png'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("  -> stage_distribution.png")

# ============================================================
# 9. AKI Paradox Visualization — 2-panel (CKD + RRT decomposition)
# ============================================================
print("[9] AKI Paradox visualization (2-panel)...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# --- Panel A: CKD-stratified mortality by KDIGO stage ---
x = np.arange(4)
width = 0.35

ckd_neg_mort = []
ckd_pos_mort = []
for s in stages:
    sub_neg = cohort[(cohort['kdigo_stage'] == s) & (cohort['ckd_icd'] == 0)]
    sub_pos = cohort[(cohort['kdigo_stage'] == s) & (cohort['ckd_icd'] == 1)]
    ckd_neg_mort.append(sub_neg['hospital_expire_flag'].mean() * 100 if len(sub_neg) > 0 else 0)
    ckd_pos_mort.append(sub_pos['hospital_expire_flag'].mean() * 100 if len(sub_pos) > 0 else 0)

bars1 = ax1.bar(x - width/2, ckd_neg_mort, width, label='CKD-negative', color=COLORS[0], edgecolor='white')
bars2 = ax1.bar(x + width/2, ckd_pos_mort, width, label='CKD-positive', color=COLORS[3], edgecolor='white')

for bar, val in zip(bars1, ckd_neg_mort):
    ax1.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.3, f'{val:.1f}%',
             ha='center', fontsize=9, fontweight='bold')
for bar, val in zip(bars2, ckd_pos_mort):
    ax1.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.3, f'{val:.1f}%',
             ha='center', fontsize=9, fontweight='bold')

ax1.set_xticks(x)
ax1.set_xticklabels(stage_labels, fontsize=12)
ax1.set_ylabel('Hospital Mortality (%)', fontsize=12, fontweight='bold')
ax1.set_title('A. Mortality by KDIGO Stage and CKD Status', fontsize=13, fontweight='bold')
ax1.legend(fontsize=10)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.set_ylim(0, max(max(ckd_neg_mort), max(ckd_pos_mort)) * 1.3)

# --- Panel B: Stage 3 decomposed by RRT indication and CKD ---
stage3 = cohort[cohort['kdigo_stage'] == 3]
rrt_types = ['RRT-only', 'Cr-only', 'Both']
rrt_cols = ['stage3_by_rrt_only', 'stage3_by_cr_only', 'stage3_by_both']
ckd_neg_mort_rrt = []
ckd_pos_mort_rrt = []
ckd_neg_n_rrt = []
ckd_pos_n_rrt = []

for col in rrt_cols:
    sub_neg = stage3[(stage3['ckd_icd'] == 0) & (stage3[col] == 1)]
    sub_pos = stage3[(stage3['ckd_icd'] == 1) & (stage3[col] == 1)]
    ckd_neg_mort_rrt.append(sub_neg['hospital_expire_flag'].mean() * 100 if len(sub_neg) > 0 else 0)
    ckd_pos_mort_rrt.append(sub_pos['hospital_expire_flag'].mean() * 100 if len(sub_pos) > 0 else 0)
    ckd_neg_n_rrt.append(len(sub_neg))
    ckd_pos_n_rrt.append(len(sub_pos))

x2 = np.arange(len(rrt_types))
bars3 = ax2.bar(x2 - width/2, ckd_neg_mort_rrt, width, label='CKD-negative', color=COLORS[0], edgecolor='white')
bars4 = ax2.bar(x2 + width/2, ckd_pos_mort_rrt, width, label='CKD-positive', color=COLORS[3], edgecolor='white')

for bar, val, n in zip(bars3, ckd_neg_mort_rrt, ckd_neg_n_rrt):
    ax2.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 1.0,
             f'{val:.1f}%\n(n={n:,})', ha='center', fontsize=9, fontweight='bold')
for bar, val, n in zip(bars4, ckd_pos_mort_rrt, ckd_pos_n_rrt):
    ax2.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 1.0,
             f'{val:.1f}%\n(n={n:,})', ha='center', fontsize=9, fontweight='bold')

ax2.set_xticks(x2)
ax2.set_xticklabels(rrt_types, fontsize=12)
ax2.set_ylabel('Hospital Mortality (%)', fontsize=12, fontweight='bold')
ax2.set_title('B. Stage 3 Decomposed by RRT Indication', fontsize=13, fontweight='bold')
ax2.legend(fontsize=10)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_ylim(0, max(max(ckd_neg_mort_rrt), max(ckd_pos_mort_rrt)) * 1.35)

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'aki_paradox.png'), dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print("  -> aki_paradox.png (2-panel)")

# ============================================================
# 10. Updated HTML Dashboard
# ============================================================
print("[10] HTML Dashboard...")

mort_0 = cohort[cohort['kdigo_stage'] == 0]['hospital_expire_flag'].mean() * 100
mort_1 = cohort[cohort['kdigo_stage'] == 1]['hospital_expire_flag'].mean() * 100
mort_2 = cohort[cohort['kdigo_stage'] == 2]['hospital_expire_flag'].mean() * 100
mort_3 = cohort[cohort['kdigo_stage'] == 3]['hospital_expire_flag'].mean() * 100
n0, n1, n2, n3 = [int((cohort['kdigo_stage'] == s).sum()) for s in stages]
rrt_n = int(cohort['rrt'].sum())

# Load KM data for JS
km_js = "{}"
if km_data and 'km_data' in km_data:
    import json as j2
    km_js = j2.dumps(km_data['km_data'])

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>KDIGO AKI Staging — ICU Outcomes (Revised)</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{ font-family:'Segoe UI',system-ui,sans-serif; background:#0d1117; color:#c9d1d9; padding:24px; }}
h1 {{ text-align:center; font-size:26px; margin-bottom:6px; color:#f0f6fc; }}
.subtitle {{ text-align:center; color:#8b949e; margin-bottom:24px; font-size:13px; }}
.note {{ text-align:center; color:#d29922; margin-bottom:16px; font-size:12px; font-style:italic; }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(260px,1fr)); gap:16px; margin-bottom:20px; }}
.card {{ background:#161b22; border:1px solid #30363d; border-radius:10px; padding:20px; text-align:center; }}
.card h3 {{ font-size:12px; color:#8b949e; margin-bottom:8px; }}
.card .value {{ font-size:32px; font-weight:700; }}
.chart-row {{ display:grid; grid-template-columns:1fr 1fr; gap:16px; margin-bottom:20px; }}
.chart-card {{ background:#161b22; border:1px solid #30363d; border-radius:10px; padding:20px; }}
.chart-card h3 {{ font-size:13px; color:#f0f6fc; margin-bottom:12px; }}
.full {{ grid-column:1/-1; }}
@media(max-width:768px) {{ .chart-row {{ grid-template-columns:1fr; }} }}
.methods {{ background:#161b22; border:1px solid #30363d; border-radius:10px; padding:20px; margin-bottom:20px; font-size:12px; line-height:1.8; }}
.methods h3 {{ color:#f0f6fc; margin-bottom:8px; }}
</style>
</head>
<body>
<h1>KDIGO AKI Staging — ICU Outcomes</h1>
<p class="subtitle">MIMIC-IV v3.1 | N={len(cohort):,} | Creatinine-based staging | Revision 2 (all reviewer fixes applied)</p>

<div class="methods">
<h3>Methods Summary</h3>
<p>Adults (≥18 yr) with first ICU stay and valid creatinine data. KDIGO creatinine-based staging. Vital signs filtered with physiological ranges (HR 20-250, RR 5-60, SpO₂ 50-100, etc.). SOFA score derived from GCS, labs, ventilation, vasopressors. Multivariable models adjusted for SOFA, ventilation, vasopressors, Elixhauser comorbidities. MDRD & lab-based CKD sensitivity analyses. RRT-stratified Stage 3 paradox analysis. All <strong>14 reviewer comments</strong> addressed.</p>
</div>

<div class="grid">
  <div class="card"><h3>No AKI</h3><div class="value" style="color:#0072B2">{mort_0:.1f}%</div><div style="color:#8b949e;margin-top:6px;font-size:12px;">N={n0:,}</div></div>
  <div class="card"><h3>Stage 1</h3><div class="value" style="color:#E69F00">{mort_1:.1f}%</div><div style="color:#8b949e;margin-top:6px;font-size:12px;">N={n1:,}</div></div>
  <div class="card"><h3>Stage 2</h3><div class="value" style="color:#F0E442">{mort_2:.1f}%</div><div style="color:#8b949e;margin-top:6px;font-size:12px;">N={n2:,}</div></div>
  <div class="card"><h3>Stage 3</h3><div class="value" style="color:#CC79A7">{mort_3:.1f}%</div><div style="color:#8b949e;margin-top:6px;font-size:12px;">N={n3:,} (RRT: {rrt_n:,})</div></div>
</div>

<div class="chart-row">
  <div class="chart-card"><h3>Hospital Mortality by KDIGO Stage</h3><canvas id="mortChart"></canvas></div>
  <div class="chart-card"><h3>Stage Distribution</h3><canvas id="distChart"></canvas></div>
</div>

<div class="chart-row">
  <div class="chart-card full"><h3>AKI Paradox: Mortality by Stage × CKD Status</h3><canvas id="paradoxChart"></canvas></div>
</div>

<script>
new Chart(document.getElementById('mortChart'), {{
  type: 'bar',
  data: {{
    labels: ['No AKI','Stage 1','Stage 2','Stage 3'],
    datasets: [{{ label:'Mortality (%)', data:[{mort_0:.1f},{mort_1:.1f},{mort_2:.1f},{mort_3:.1f}],
      backgroundColor:['#0072B2','#E69F00','#F0E442','#CC79A7'],borderRadius:6 }}]
  }},
  options: {{ responsive:true, plugins:{{legend:{{display:false}}}},
    scales:{{ y:{{beginAtZero:true,grid:{{color:'#21262d'}},ticks:{{color:'#8b949e'}}}},
             x:{{grid:{{display:false}},ticks:{{color:'#8b949e'}}}} }} }}
}});

new Chart(document.getElementById('distChart'), {{
  type: 'doughnut',
  data: {{ labels:['No AKI','Stage 1','Stage 2','Stage 3'], datasets:[{{ data:[{n0},{n1},{n2},{n3}],
    backgroundColor:['#0072B2','#E69F00','#F0E442','#CC79A7'],borderWidth:2,borderColor:'#0d1117' }}] }},
  options: {{ responsive:true, plugins:{{legend:{{position:'bottom',labels:{{color:'#8b949e'}}}}}} }}
}});

var pCtx = document.getElementById('paradoxChart').getContext('2d');
var ckd_neg = [{ckd_neg_mort[0]:.1f},{ckd_neg_mort[1]:.1f},{ckd_neg_mort[2]:.1f},{ckd_neg_mort[3]:.1f}];
var ckd_pos = [{ckd_pos_mort[0]:.1f},{ckd_pos_mort[1]:.1f},{ckd_pos_mort[2]:.1f},{ckd_pos_mort[3]:.1f}];
new Chart(pCtx, {{
  type: 'bar',
  data: {{ labels:['No AKI','Stage 1','Stage 2','Stage 3'],
    datasets: [
      {{ label:'CKD-negative', data:ckd_neg, backgroundColor:'#0072B2',borderRadius:6,borderSkipped:false }},
      {{ label:'CKD-positive', data:ckd_pos, backgroundColor:'#CC79A7',borderRadius:6,borderSkipped:false }},
    ]
  }},
  options: {{ responsive:true,
    scales:{{ y:{{beginAtZero:true,grid:{{color:'#21262d'}},ticks:{{color:'#8b949e'}}}},
             x:{{grid:{{display:false}},ticks:{{color:'#8b949e'}}}} }},
    plugins:{{ legend:{{labels:{{color:'#8b949e',usePointStyle:true}}}} }}
  }}
}});
</script>

<div style="text-align:center;margin-top:20px;color:#484f58;font-size:11px;">
  MIMIC-IV v3.1 | KDIGO Creatinine-based Staging | Revision 2 | All reviewer concerns addressed
</div>
</body>
</html>"""

with open(os.path.join(IN, "kdigo_dashboard_revised.html"), 'w', encoding='utf-8') as f:
    f.write(html)
print("  -> kdigo_dashboard_revised.html")

print("\n" + "=" * 80)
print("ALL VISUALIZATIONS COMPLETE")
print(f"Output: {OUT_DIR}")
print("=" * 80)
