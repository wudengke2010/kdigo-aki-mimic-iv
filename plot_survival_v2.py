"""
KM 曲线图 — MIMIC vs eICU 并排比较
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

WORK = r"C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
os.makedirs(os.path.join(WORK, "figures_prism"), exist_ok=True)

colors = {0:'#2C3E50', 1:'#3498DB', 2:'#E67E22', 3:'#C0392B'}

fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=False)
for ax, db, km_file, title in [
    (axes[0], 'MIMIC-IV', 'km_curves_mimic.csv', 'MIMIC-IV v3.1 (N=58,554)'),
    (axes[1], 'eICU-CRD', 'km_curves_eicu.csv', 'eICU-CRD v2.0 (N=107,758)')]:
    df = pd.read_csv(os.path.join(WORK, "v2_outputs", km_file))
    for s in [0,1,2,3]:
        sub = df[df.stage == s].sort_values('time_h')
        ax.step(sub['time_h']/24, sub['survival'],
                where='post', label=f'Stage {s} (n={len(sub):,})',
                color=colors[s], linewidth=2)
    ax.set_xlim(0, 30)
    ax.set_ylim(0.3, 1.0)
    ax.set_xlabel('Days from ICU admission (24h landmark)')
    ax.set_ylabel('Survival probability')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='lower left', fontsize=9)

plt.suptitle('Kaplan-Meier survival by KDIGO stage — landmark 24h, 30-day horizon',
             fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(WORK, "figures_prism", "fig_km_mimic_vs_eicu.png"), dpi=200, bbox_inches='tight')
plt.savefig(os.path.join(WORK, "figures_prism", "fig_km_mimic_vs_eicu.pdf"), bbox_inches='tight')
print("Saved fig_km_mimic_vs_eicu.png/pdf")

# Forest plot of time-stratified Cox
import json
import numpy as np
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
for ax, jf, title in [
    (axes[0], 'survival_mimic.json', 'MIMIC-IV v3.1'),
    (axes[1], 'survival_eicu.json', 'eICU-CRD v2.0')]:
    d = json.load(open(os.path.join(WORK, "v2_outputs", jf)))
    ts = d['time_stratified']
    # 3 strata x 3 stages
    for s_idx, (key, lab) in enumerate([('0_7d','0-7d'),('7_14d','7-14d'),('14_30d','14-30d')]):
        for stage_idx, v, lab2 in [(0,'kdigo_1','Stage 1'),(1,'kdigo_2','Stage 2'),(2,'kdigo_3','Stage 3')]:
            r = ts[key]['kdigo'][v]
            y = s_idx*3 + stage_idx
            xerr = [[r['HR'] - r['CI_lo']], [r['CI_hi'] - r['HR']]]
            ax.errorbar(r['HR'], y, xerr=xerr, fmt='o', color=colors[s_idx+1], capsize=3, markersize=7)
            ax.text(r['HR']+0.15, y, f"{r['HR']:.2f}", va='center', fontsize=8)
    ax.axvline(1, color='gray', linestyle='--', alpha=0.6)
    ax.set_yticks(range(9))
    ax.set_yticklabels([f"{s} S{k+1}" for s,lab in enumerate([('0-7d'),('7-14d'),('14-30d')]) for k in range(3)])
    ax.set_xlabel('Hazard ratio (log scale)' if False else 'Hazard ratio')
    ax.set_xscale('log')
    ax.set_xlim(0.5, 10)
    ax.set_title(f'{title}\nTime-stratified Cox (Model C)')
    ax.grid(True, alpha=0.3, axis='x')

plt.suptitle('Time-stratified hazard ratios for KDIGO stage (vs Stage 0)',
             fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(WORK, "figures_prism", "fig_forest_time_stratified.png"), dpi=200, bbox_inches='tight')
plt.savefig(os.path.join(WORK, "figures_prism", "fig_forest_time_stratified.pdf"), bbox_inches='tight')
print("Saved fig_forest_time_stratified.png/pdf")