#!/usr/bin/env python3
"""Regenerate mortality_by_stage.png with fixed label spacing."""
import pandas as pd
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

IN = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
OUT_DIR = os.path.join(IN, "figures_revised")
os.makedirs(OUT_DIR, exist_ok=True)

COLORS = ['#0072B2', '#E69F00', '#D55E00', '#CC79A7', '#009E73', '#F0E442', '#56B4E9', '#000000']
STAGE_COLORS = [COLORS[0], COLORS[1], COLORS[5], COLORS[3]]

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial']
plt.rcParams['axes.unicode_minus'] = False

cohort_path = os.path.join(IN, "kdigo_cohort_revised_enriched.csv")
if not os.path.exists(cohort_path):
    cohort_path = os.path.join(IN, "kdigo_cohort_revised.csv")
cohort = pd.read_csv(cohort_path, low_memory=False)

stages = [0, 1, 2, 3]
stage_labels = ['No AKI', 'Stage 1', 'Stage 2', 'Stage 3']

def wilson_ci(p, n, z=1.96):
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
print(f"Regenerated: {os.path.join(OUT_DIR, 'mortality_by_stage.png')}")
