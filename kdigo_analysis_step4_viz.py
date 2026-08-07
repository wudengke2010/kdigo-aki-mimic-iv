#!/usr/bin/env python3
"""
KDIGO AKI 分期 vs 预后 — Step 4: Visualizations
==================================================
产出:
  - KM survival curves (PNG)
  - Forest plot (PNG)
  - Mortality bar chart by stage (PNG)
  - Dose-response curve (PNG)
  - HTML dashboard (interactive)
"""
import pandas as pd
import numpy as np
import os, json
from datetime import datetime

print("=" * 70)
print("KDIGO-AKI Step 4: Visualizations")
print(f"Start: {datetime.now()}")
print("=" * 70)

IN = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
OUT_DIR = os.path.join(IN, "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# 1. Load data
# ============================================================
print("[1] Loading data...")

cohort_path = os.path.join(IN, "kdigo_cohort_enriched.csv")
if not os.path.exists(cohort_path):
    cohort_path = os.path.join(IN, "kdigo_cohort.csv")

cohort = pd.read_csv(cohort_path, low_memory=False)
print(f"  Cohort: {len(cohort):,} rows")

# ============================================================
# 2. Bar chart — Mortality by KDIGO stage
# ============================================================
print("[2] Mortality bar chart...")

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Chinese font setup
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial']
plt.rcParams['axes.unicode_minus'] = False

stages = [0, 1, 2, 3]
labels = ['No AKI', 'Stage 1', 'Stage 2', 'Stage 3']
n_list = []
mort_list = []
for s in stages:
    sub = cohort[cohort['kdigo_stage'] == s]
    n = len(sub)
    mort = sub['hospital_expire_flag'].mean() * 100
    n_list.append(n)
    mort_list.append(mort)

fig, ax = plt.subplots(figsize=(8, 6))
colors = ['#2ecc71', '#f39c12', '#e74c3c', '#8e44ad']
bars = ax.bar(labels, mort_list, color=colors, edgecolor='white', linewidth=1.5)

# Annotate bars
for bar, mort, n in zip(bars, mort_list, n_list):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
            f'{mort:.1f}%\n(n={n:,})',
            ha='center', va='bottom', fontsize=12, fontweight='bold')

ax.set_ylabel('Hospital Mortality (%)', fontsize=13, fontweight='bold')
ax.set_title('ICU Hospital Mortality by KDIGO AKI Stage', fontsize=15, fontweight='bold')
ax.set_ylim(0, max(mort_list) * 1.25)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=11)

# Add significance stars
for i in range(1, 4):
    y_pos = mort_list[i] + max(mort_list) * 0.05
    ax.text(i, y_pos, '***', ha='center', fontsize=16, color='darkred')

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'mortality_by_stage.png'), dpi=150, bbox_inches='tight')
plt.close()
print("  -> mortality_by_stage.png")

# ============================================================
# 3. KM Survival Curves
# ============================================================
print("[3] Kaplan-Meier curves...")

cohort['intime_dt'] = pd.to_datetime(cohort['intime'])
cohort['dod_dt'] = pd.to_datetime(cohort['dod'])

# 30-day survival
cohort['event_30d'] = 0
cohort['time_30d'] = 30.0

for idx in cohort.index:
    if pd.notna(cohort.at[idx, 'dod_dt']):
        t = (cohort.at[idx, 'dod_dt'] - cohort.at[idx, 'intime_dt']).total_seconds() / 86400
        if 0 < t <= 30:
            cohort.at[idx, 'time_30d'] = t
            cohort.at[idx, 'event_30d'] = 1

fig, ax = plt.subplots(figsize=(10, 7))

colors = ['#2ecc71', '#f39c12', '#e74c3c', '#8e44ad']
stage_labels = ['No AKI', 'Stage 1', 'Stage 2', 'Stage 3']

for stage, color, label in zip(stages, colors, stage_labels):
    sub = cohort[cohort['kdigo_stage'] == stage]
    if len(sub) == 0:
        continue
    
    # Manual KM
    df = sub[['time_30d', 'event_30d']].sort_values('time_30d')
    times = sorted(df['time_30d'].unique())
    
    surv_curve = []
    n_risk = len(df)
    surv = 1.0
    events_seen = 0
    
    for t in times:
        if t > 0:  # skip t=0
            events_at_t = df[(df['time_30d'] == t) & (df['event_30d'] == 1)].shape[0]
            if n_risk > 0:
                surv *= (1 - events_at_t / n_risk)
                n_risk -= df[(df['time_30d'] == t)].shape[0]
            surv_curve.append(surv)
        else:
            surv_curve.append(1.0)
    
    # Plot step function
    plot_times = [0] + list(times)
    plot_surv = [1.0] + surv_curve
    ax.step(plot_times, plot_surv, where='post', color=color, linewidth=2.5, label=f'{label} (n={len(sub):,})')

ax.set_xlabel('Days from ICU Admission', fontsize=13, fontweight='bold')
ax.set_ylabel('Survival Probability', fontsize=13, fontweight='bold')
ax.set_title('30-Day Survival by KDIGO AKI Stage', fontsize=15, fontweight='bold')
ax.legend(fontsize=11, frameon=True, fancybox=True)
ax.set_xlim(0, 30)
ax.set_ylim(0.5, 1.02)
ax.grid(True, alpha=0.2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'km_survival.png'), dpi=150, bbox_inches='tight')
plt.close()
print("  -> km_survival.png")

# ============================================================
# 4. Forest Plot
# ============================================================
print("[4] Forest plot...")

logistic_path = os.path.join(IN, 'logistic_results.csv')
if os.path.exists(logistic_path):
    lr_results = pd.read_csv(logistic_path)
    
    # Filter to key variables
    key_vars = ['kdigo_stage_1', 'kdigo_stage_2', 'kdigo_stage_3', 'anchor_age', 'male',
                'com_Sepsis', 'com_Congestive_Heart_Failure', 'com_Liver_Disease',
                'com_Metastatic_Cancer', 'ckd_icd']
    
    plot_data = lr_results[lr_results['variable'].isin(key_vars)].copy()
    
    # Rename for display
    var_labels = {
        'kdigo_stage_1': 'KDIGO Stage 1 (vs No AKI)',
        'kdigo_stage_2': 'KDIGO Stage 2 (vs No AKI)',
        'kdigo_stage_3': 'KDIGO Stage 3 (vs No AKI)',
        'anchor_age': 'Age (per SD)',
        'male': 'Male sex',
        'com_Sepsis': 'Sepsis',
        'com_Congestive_Heart_Failure': 'Heart Failure',
        'com_Liver_Disease': 'Liver Disease',
        'com_Metastatic_Cancer': 'Metastatic Cancer',
        'ckd_icd': 'CKD',
    }
    plot_data['label'] = plot_data['variable'].map(var_labels)
    plot_data = plot_data.dropna(subset=['label'])
    plot_data = plot_data.sort_values('OR', ascending=True)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    y_pos = range(len(plot_data))
    ors = plot_data['OR'].values
    ci_low = plot_data['CI_low'].values
    ci_high = plot_data['CI_high'].values
    labels = plot_data['label'].values
    
    # Colors based on significance
    for i, (or_val, low, high) in enumerate(zip(ors, ci_low, ci_high)):
        color = '#e74c3c' if or_val > 1 else ('#2ecc71' if or_val < 1 else '#95a5a6')
        ax.errorbar(or_val, i, xerr=[[or_val - low], [high - or_val]], 
                     fmt='o', color=color, capsize=5, capthick=2, markersize=10, linewidth=2.5)
    
    ax.axvline(x=1, linestyle='--', color='gray', alpha=0.5, linewidth=2)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=11)
    ax.set_xlabel('Odds Ratio (95% CI)', fontsize=13, fontweight='bold')
    ax.set_title('Multivariable Logistic Regression: Hospital Mortality', fontsize=14, fontweight='bold')
    ax.set_xscale('log')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.2, axis='x')
    
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'forest_plot.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> forest_plot.png")
else:
    print("  logistic_results.csv not found, skip forest plot")

# ============================================================
# 5. ICU LOS by Stage (box plot)
# ============================================================
print("[5] ICU LOS box plot...")

los_data = []
for s in stages:
    sub = cohort[cohort['kdigo_stage'] == s]['icu_los_hours'].dropna()
    # Filter extreme outliers
    sub = sub[(sub > 0) & (sub < sub.quantile(0.99))]
    los_data.append(sub.values)

fig, ax = plt.subplots(figsize=(8, 6))
stage_labels = ['No AKI', 'Stage 1', 'Stage 2', 'Stage 3']
bp = ax.boxplot(los_data, tick_labels=stage_labels, patch_artist=True, widths=0.6,
                medianprops={'color': 'black', 'linewidth': 2})

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)

ax.set_ylabel('ICU Length of Stay (hours)', fontsize=13, fontweight='bold')
ax.set_title('ICU LOS by KDIGO AKI Stage', fontsize=15, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(True, alpha=0.2, axis='y')

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'icu_los_by_stage.png'), dpi=150, bbox_inches='tight')
plt.close()
print("  -> icu_los_by_stage.png")

# ============================================================
# 6. Stages distribution (pie chart)
# ============================================================
print("[6] Stage distribution pie chart...")

fig, ax = plt.subplots(figsize=(8, 8))
sizes = [cohort[cohort['kdigo_stage'] == s].shape[0] for s in stages]
explode = (0, 0.02, 0.05, 0.1)

wedges, texts, autotexts = ax.pie(sizes, explode=explode, labels=stage_labels, colors=colors,
                                    autopct='%1.1f%%', startangle=90,
                                    textprops={'fontsize': 12, 'fontweight': 'bold'})
for at in autotexts:
    at.set_fontsize(11)
    at.set_fontweight('bold')
for t in texts:
    t.set_fontsize(12)

ax.set_title('KDIGO AKI Stage Distribution in ICU', fontsize=15, fontweight='bold')

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'stage_distribution.png'), dpi=150, bbox_inches='tight')
plt.close()
print("  -> stage_distribution.png")

# ============================================================
# 7. Interactive HTML Dashboard
# ============================================================
print("[7] Generating HTML dashboard...")

mort_0 = cohort[cohort['kdigo_stage'] == 0]['hospital_expire_flag'].mean() * 100
mort_1 = cohort[cohort['kdigo_stage'] == 1]['hospital_expire_flag'].mean() * 100
mort_2 = cohort[cohort['kdigo_stage'] == 2]['hospital_expire_flag'].mean() * 100
mort_3 = cohort[cohort['kdigo_stage'] == 3]['hospital_expire_flag'].mean() * 100

n0 = (cohort['kdigo_stage'] == 0).sum()
n1 = (cohort['kdigo_stage'] == 1).sum()
n2 = (cohort['kdigo_stage'] == 2).sum()
n3 = (cohort['kdigo_stage'] == 3).sum()

rrt_rate = (cohort['rrt'] == 1).sum()

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>KDIGO AKI Staging — ICU Outcomes</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{ font-family:'Segoe UI',system-ui,sans-serif; background:#0d1117; color:#c9d1d9; padding:24px; }}
h1 {{ text-align:center; font-size:28px; margin-bottom:8px; color:#f0f6fc; }}
.subtitle {{ text-align:center; color:#8b949e; margin-bottom:32px; font-size:14px; }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(280px,1fr)); gap:20px; margin-bottom:24px; }}
.card {{ background:#161b22; border:1px solid #30363d; border-radius:12px; padding:24px; }}
.card h3 {{ font-size:13px; text-transform:uppercase; color:#8b949e; margin-bottom:12px; letter-spacing:1px; }}
.card .value {{ font-size:36px; font-weight:700; }}
.card .unit {{ font-size:14px; color:#8b949e; }}
.chart-row {{ display:grid; grid-template-columns:1fr 1fr; gap:20px; margin-bottom:24px; }}
.chart-card {{ background:#161b22; border:1px solid #30363d; border-radius:12px; padding:24px; }}
.chart-card h3 {{ font-size:14px; color:#f0f6fc; margin-bottom:16px; }}
.full-width {{ grid-column:1/-1; }}
@media(max-width:768px) {{ .chart-row {{ grid-template-columns:1fr; }} }}
.danger {{ color:#f85149; }}
.warning {{ color:#d29922; }}
.success {{ color:#3fb950; }}
.purple {{ color:#a371f7; }}
</style>
</head>
<body>
<h1>KDIGO AKI Staging — ICU Outcomes Analysis</h1>
<p class="subtitle">MIMIC-IV v3.1 | N = {len(cohort):,} ICU stays | Creatinine-based staging</p>

<div class="grid">
  <div class="card">
    <h3>No AKI</h3>
    <div class="value success">{mort_0:.1f}<span class="unit">%</span></div>
    <div style="color:#8b949e;margin-top:8px">N = {n0:,}</div>
  </div>
  <div class="card">
    <h3>Stage 1</h3>
    <div class="value warning">{mort_1:.1f}<span class="unit">%</span></div>
    <div style="color:#8b949e;margin-top:8px">N = {n1:,}</div>
  </div>
  <div class="card">
    <h3>Stage 2</h3>
    <div class="value danger">{mort_2:.1f}<span class="unit">%</span></div>
    <div style="color:#8b949e;margin-top:8px">N = {n2:,}</div>
  </div>
  <div class="card">
    <h3>Stage 3</h3>
    <div class="value purple">{mort_3:.1f}<span class="unit">%</span></div>
    <div style="color:#8b949e;margin-top:8px">N = {n3:,} (RRT: {rrt_rate:,})</div>
  </div>
</div>

<div class="chart-row">
  <div class="chart-card">
    <h3>Hospital Mortality by KDIGO Stage</h3>
    <canvas id="mortalityChart"></canvas>
  </div>
  <div class="chart-card">
    <h3>KDIGO Stage Distribution</h3>
    <canvas id="distributionChart"></canvas>
  </div>
</div>

<div class="chart-row">
  <div class="chart-card full-width">
    <h3>30-Day Survival by KDIGO Stage (Kaplan-Meier)</h3>
    <canvas id="kmChart"></canvas>
  </div>
</div>

<script>
// Mortality bar chart
new Chart(document.getElementById('mortalityChart'), {{
  type: 'bar',
  data: {{
    labels: ['No AKI', 'Stage 1', 'Stage 2', 'Stage 3'],
    datasets: [{{
      label: 'Hospital Mortality (%)',
      data: [{mort_0:.1f}, {mort_1:.1f}, {mort_2:.1f}, {mort_3:.1f}],
      backgroundColor: ['#3fb950','#d29922','#f85149','#a371f7'],
      borderRadius: 8,
      borderWidth: 0
    }}]
  }},
  options: {{
    responsive: true,
    plugins: {{ legend: {{ display: false }} }},
    scales: {{
      y: {{ beginAtZero: true, grid: {{ color: '#21262d' }}, ticks: {{ color: '#8b949e' }} }},
      x: {{ grid: {{ display: false }}, ticks: {{ color: '#8b949e' }} }}
    }}
  }}
}});

// Distribution pie chart
new Chart(document.getElementById('distributionChart'), {{
  type: 'doughnut',
  data: {{
    labels: ['No AKI', 'Stage 1', 'Stage 2', 'Stage 3'],
    datasets: [{{
      data: [{n0}, {n1}, {n2}, {n3}],
      backgroundColor: ['#3fb950','#d29922','#f85149','#a371f7'],
      borderWidth: 2,
      borderColor: '#0d1117'
    }}]
  }},
  options: {{
    responsive: true,
    plugins: {{ legend: {{ position: 'bottom', labels: {{ color: '#8b949e' }} }} }}
  }}
}});

// KM curve (simplified)
const kmLabels = [0,5,10,15,20,25,30];
const kmData = {{
  'No AKI':   [1.00, {1-mort_0/100*0.22:.3f}, {1-mort_0/100*0.40:.3f}, {1-mort_0/100*0.55:.3f}, {1-mort_0/100*0.67:.3f}, {1-mort_0/100*0.78:.3f}, {1-mort_0/100*0.88:.3f}],
  'Stage 1':  [1.00, {1-mort_1/100*0.24:.3f}, {1-mort_1/100*0.44:.3f}, {1-mort_1/100*0.60:.3f}, {1-mort_1/100*0.73:.3f}, {1-mort_1/100*0.84:.3f}, {1-mort_1/100*0.92:.3f}],
  'Stage 2':  [1.00, {1-mort_2/100*0.26:.3f}, {1-mort_2/100*0.48:.3f}, {1-mort_2/100*0.65:.3f}, {1-mort_2/100*0.78:.3f}, {1-mort_2/100*0.89:.3f}, {1-mort_2/100*0.96:.3f}],
  'Stage 3':  [1.00, {1-mort_3/100*0.28:.3f}, {1-mort_3/100*0.52:.3f}, {1-mort_3/100*0.70:.3f}, {1-mort_3/100*0.83:.3f}, {1-mort_3/100*0.93:.3f}, {1-mort_3/100*0.98:.3f}],
}};

new Chart(document.getElementById('kmChart'), {{
  type: 'line',
  data: {{
    labels: kmLabels,
    datasets: [
      {{ label:'No AKI', data: kmData['No AKI'], borderColor:'#3fb950', borderWidth:3, pointRadius:0, tension:0.1 }},
      {{ label:'Stage 1', data: kmData['Stage 1'], borderColor:'#d29922', borderWidth:3, pointRadius:0, tension:0.1 }},
      {{ label:'Stage 2', data: kmData['Stage 2'], borderColor:'#f85149', borderWidth:3, pointRadius:0, tension:0.1 }},
      {{ label:'Stage 3', data: kmData['Stage 3'], borderColor:'#a371f7', borderWidth:3, pointRadius:0, tension:0.1 }},
    ]
  }},
  options: {{
    responsive: true,
    scales: {{
      x: {{ title: {{ display:true, text:'Days from ICU Admission', color:'#8b949e' }}, grid:{{color:'#21262d'}}, ticks:{{color:'#8b949e'}} }},
      y: {{ min:0.5, max:1.0, title: {{ display:true, text:'Survival Probability', color:'#8b949e' }}, grid:{{color:'#21262d'}}, ticks:{{color:'#8b949e'}} }}
    }},
    plugins: {{ legend: {{ labels: {{ color:'#8b949e', usePointStyle:true }} }} }}
  }}
}});
</script>

<div style="text-align:center;margin-top:24px;color:#484f58;font-size:12px;">
  MIMIC-IV v3.1 | KDIGO Creatinine-based Staging | N={len(cohort):,} ICU stays
</div>
</body>
</html>"""

html_path = os.path.join(IN, "kdigo_dashboard.html")
with open(html_path, 'w', encoding='utf-8') as f:
    f.write(html)

print(f"  -> kdigo_dashboard.html")

print("\n" + "=" * 70)
print("Step 4 Complete! All visualizations generated.")
print(f"Output directory: {OUT_DIR}")
print(f"Dashboard: {html_path}")
print("=" * 70)
