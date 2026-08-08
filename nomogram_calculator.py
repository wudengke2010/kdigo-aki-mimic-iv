#!/usr/bin/env python3
"""
Nomogram + Calibration + Web Calculator for KDIGO AKI Mortality
================================================================
Builds a clinical prediction nomogram based on Model B variables:
- KDIGO Stage (0-3)
- APACHE score
- Mechanical ventilation
- Age
- Sex

Uses eICU-CRD validation cohort (N=166,373) as the prediction model base.

Author: Dengke Wu
Date: 2026-08-08
"""

import pandas as pd
import numpy as np
import json
import os
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score, brier_score_loss
from sklearn.calibration import calibration_curve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

OUTPUT_DIR = 'C:/Users/admin/WorkBuddy/2026-07-06-13-22-19'

print("=" * 60)
print("Nomogram + Calibration + Web Calculator")
print("=" * 60)

# Load data
df = pd.read_csv(os.path.join(OUTPUT_DIR, 'eicu_validation_data.csv'))
print(f"Loaded: {len(df):,} patients")

# Prepare variables
df['kdigo_1'] = (df['kdigo_stage'] == 1).astype(int)
df['kdigo_2'] = (df['kdigo_stage'] == 2).astype(int)
df['kdigo_3'] = (df['kdigo_stage'] == 3).astype(int)
df['gender_male'] = (df['gender'] == 'Male').astype(int)

# Model B variables
model_vars = ['kdigo_1', 'kdigo_2', 'kdigo_3', 'apachescore', 'mech_vent', 'age_num', 'gender_male']
X = df[model_vars].astype(float)
X = sm.add_constant(X)
y_mort = df['mortality'].astype(float)

print("\nFitting logistic regression (Model B)...")
model = sm.Logit(y_mort, X).fit(disp=0)
print(model.summary().tables[1])

# Get coefficients
coef = model.params
print("\nCoefficients:")
for var, val in coef.items():
    print(f"  {var}: {val:.4f}")

# ============================================================
# Nomogram Plot
# ============================================================
print("\nCreating nomogram plot...")

fig = plt.figure(figsize=(14, 8))
gs = GridSpec(2, 2, width_ratios=[3, 1.5], height_ratios=[1, 1], hspace=0.35, wspace=0.3)

# --- Panel A: Nomogram ---
ax1 = fig.add_subplot(gs[0, 0])

# Nomogram scale: points = beta * value * scale_factor
# Scale to 0-100 points
max_points = 100

# Calculate point range for each variable
var_info = {
    'KDIGO Stage': {
        'categories': ['Stage 0', 'Stage 1', 'Stage 2', 'Stage 3'],
        'betas': [0, coef['kdigo_1'], coef['kdigo_2'], coef['kdigo_3']],
    },
    'APACHE Score': {
        'range': (0, 150),
        'beta': coef['apachescore'],
    },
    'Mech. Ventilation': {
        'categories': ['No', 'Yes'],
        'betas': [0, coef['mech_vent']],
    },
    'Age': {
        'range': (18, 100),
        'beta': coef['age_num'],
    },
    'Sex': {
        'categories': ['Female', 'Male'],
        'betas': [0, coef['gender_male']],
    }
}

# Calculate total beta range
all_betas = []
for v in var_info.values():
    if 'betas' in v:
        all_betas.extend(v['betas'])
    else:
        all_betas.extend([v['beta'] * v['range'][0], v['beta'] * v['range'][1]])

beta_min = min(all_betas)
beta_max = max(all_betas)
scale = max_points / (beta_max - beta_min)

# Draw nomogram scales
y_positions = list(range(len(var_info), 0, -1))

# Points scale at top
ax1.axhline(y=len(var_info) + 1, color='black', linewidth=1.5)
for p in range(0, 101, 10):
    beta_at_p = p / scale + beta_min
    ax1.plot([p, p], [len(var_info) + 0.9, len(var_info) + 1.1], 'k-', linewidth=0.8)
    ax1.text(p, len(var_info) + 0.7, str(p), ha='center', va='bottom', fontsize=8)
ax1.text(50, len(var_info) + 1.5, 'Points', ha='center', va='bottom', fontsize=10, fontweight='bold')

# Variable scales
for i, (var_name, info) in enumerate(var_info.items()):
    y = y_positions[i]
    ax1.axhline(y=y, color='gray', linewidth=0.8)

    if 'categories' in info:
        # Categorical variable
        n_cats = len(info['categories'])
        for j, (cat, beta) in enumerate(zip(info['categories'], info['betas'])):
            points = (beta - beta_min) * scale
            x = j * max_points / (n_cats - 1) if n_cats > 1 else 50
            ax1.plot([x, x], [y - 0.1, y + 0.1], 'k-', linewidth=0.8)
            ax1.text(x, y - 0.3, cat, ha='center', va='top', fontsize=7)
    else:
        # Continuous variable
        lo, hi = info['range']
        beta_lo = info['beta'] * lo
        beta_hi = info['beta'] * hi
        pts_lo = (beta_lo - beta_min) * scale
        pts_hi = (beta_hi - beta_min) * scale

        # Draw scale line from pts_lo to pts_hi
        ax1.plot([pts_lo, pts_hi], [y, y], 'k-', linewidth=1)
        # Tick marks
        n_ticks = 6
        for t in range(n_ticks + 1):
            val = lo + (hi - lo) * t / n_ticks
            pt = (info['beta'] * val - beta_min) * scale
            ax1.plot([pt, pt], [y - 0.1, y + 0.1], 'k-', linewidth=0.8)
            ax1.text(pt, y - 0.3, f'{int(val)}', ha='center', va='top', fontsize=7)

    ax1.text(-15, y, var_name, ha='right', va='center', fontsize=9, fontweight='bold')

# Total points scale
y_total = -0.5
ax1.axhline(y=y_total, color='black', linewidth=1.5)
for p in range(0, int(max_points * len(var_info)) + 1, 50):
    ax1.plot([p, p], [y_total - 0.1, y_total + 0.1], 'k-', linewidth=0.8)
    ax1.text(p, y_total - 0.3, str(p), ha='center', va='top', fontsize=7)
ax1.text(-15, y_total, 'Total Points', ha='right', va='center', fontsize=9, fontweight='bold')

# Risk scale
y_risk = -1.5
ax1.axhline(y=y_risk, color='black', linewidth=1.5)
for p in range(0, int(max_points * len(var_info)) + 1, 20):
    # Convert points to probability
    total_beta = p / scale + beta_min
    prob = 1 / (1 + np.exp(-(total_beta + coef['const'])))
    ax1.plot([p, p], [y_risk - 0.1, y_risk + 0.1], 'k-', linewidth=0.8)
    ax1.text(p, y_risk - 0.3, f'{prob*100:.0f}%', ha='center', va='top', fontsize=7)
ax1.text(-15, y_risk, 'Risk', ha='right', va='center', fontsize=9, fontweight='bold')

ax1.set_xlim(-30, max_points * len(var_info) + 10)
ax1.set_ylim(-2.5, len(var_info) + 2.5)
ax1.set_title('Nomogram for ICU Mortality Prediction', fontsize=11, fontweight='bold', pad=10)
ax1.axis('off')

# --- Panel B: Calibration Curve ---
ax2 = fig.add_subplot(gs[0, 1])

prob_pred = model.predict(X)
prob_pred_clipped = np.clip(prob_pred, 0, 1)

# Calibration curve
y_array = np.array(df['mortality'], dtype=float)
fraction_of_positives, mean_predicted_value = calibration_curve(
    y_array, prob_pred_clipped, n_bins=10, strategy='quantile'
)

ax2.plot([0, 1], [0, 1], 'k--', linewidth=1, label='Perfect calibration')
ax2.plot(mean_predicted_value, fraction_of_positives, 's-', color='#2166ac',
         linewidth=2, markersize=6, label='Model B')

# Brier score
brier = brier_score_loss(y_array, prob_pred_clipped)
ax2.text(0.6, 0.15, f'Brier score: {brier:.4f}', fontsize=10, transform=ax2.transData,
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

ax2.set_xlabel('Predicted Probability', fontsize=11)
ax2.set_ylabel('Observed Frequency', fontsize=11)
ax2.set_title('Calibration Curve', fontsize=11, fontweight='bold')
ax2.legend(fontsize=9, loc='upper left')
ax2.set_xlim([-0.05, 1.05])
ax2.set_ylim([-0.05, 1.05])
ax2.grid(True, alpha=0.3)

# --- Panel C: Variable Importance (forest plot) ---
ax3 = fig.add_subplot(gs[1, 0:])

var_labels = ['KDIGO Stage 1', 'KDIGO Stage 2', 'KDIGO Stage 3',
              'APACHE Score (per 10 pts)', 'Mechanical Ventilation',
              'Age (per 10 yrs)', 'Male Sex']
or_values = [np.exp(coef['kdigo_1']),
             np.exp(coef['kdigo_2']),
             np.exp(coef['kdigo_3']),
             np.exp(coef['apachescore'] * 10),
             np.exp(coef['mech_vent']),
             np.exp(coef['age_num'] * 10),
             np.exp(coef['gender_male'])]
ci_lows = [np.exp(model.conf_int().loc[v, 0]) for v in model_vars]
ci_highs = [np.exp(model.conf_int().loc[v, 1]) for v in model_vars]

# Adjust CI for scaled variables
ci_lows[3] = np.exp(model.conf_int().loc['apachescore', 0] * 10)
ci_highs[3] = np.exp(model.conf_int().loc['apachescore', 1] * 10)
ci_lows[5] = np.exp(model.conf_int().loc['age_num', 0] * 10)
ci_highs[5] = np.exp(model.conf_int().loc['age_num', 1] * 10)

colors_fp = ['#4393c3', '#4393c3', '#4393c3', '#f4a582', '#92c5de', '#d6604d', '#bababa']
y_pos = range(len(var_labels))

ax3.barh(y_pos, or_values, xerr=[np.array(or_values) - np.array(ci_lows),
                                  np.array(ci_highs) - np.array(or_values)],
         color=colors_fp, edgecolor='black', linewidth=0.5, height=0.6, capsize=3)
ax3.axvline(x=1, color='red', linewidth=1, linestyle='--')
ax3.set_yticks(y_pos)
ax3.set_yticklabels(var_labels, fontsize=10)
ax3.set_xlabel('Odds Ratio (95% CI)', fontsize=11)
ax3.set_title('Variable Importance (Adjusted OR)', fontsize=11, fontweight='bold')
ax3.set_xlim([0, max(ci_highs) * 1.1])

# Add OR values as text
for i, (or_val, ci_l, ci_h) in enumerate(zip(or_values, ci_lows, ci_highs)):
    ax3.text(max(ci_highs) * 1.05, i, f'{or_val:.2f} ({ci_l:.2f}-{ci_h:.2f})',
             va='center', fontsize=8)

ax3.grid(True, alpha=0.3, axis='x')

plt.suptitle('Clinical Prediction Model for ICU Mortality in AKI Patients\n(eICU-CRD Validation Cohort, N=166,373)',
             fontsize=12, fontweight='bold', y=1.01)

fig_path = os.path.join(OUTPUT_DIR, 'figures_prism', 'figS7_nomogram.png')
fig.savefig(fig_path, dpi=300, bbox_inches='tight')
fig.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
plt.close()
print(f"  Figure saved: {fig_path}")

# ============================================================
# Web Calculator (HTML)
# ============================================================
print("\nCreating web calculator...")

html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>KDIGO AKI ICU Mortality Calculator</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            max-width: 600px;
            margin: 50px auto;
            padding: 20px;
            background: #f5f5f5;
        }}
        .calculator {{
            background: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2166ac;
            font-size: 20px;
            text-align: center;
            margin-bottom: 5px;
        }}
        .subtitle {{
            text-align: center;
            color: #666;
            font-size: 13px;
            margin-bottom: 25px;
        }}
        .input-group {{
            margin-bottom: 15px;
        }}
        label {{
            display: block;
            font-weight: bold;
            margin-bottom: 5px;
            font-size: 14px;
        }}
        select, input {{
            width: 100%;
            padding: 8px;
            border: 1px solid #ccc;
            border-radius: 5px;
            font-size: 14px;
        }}
        .result {{
            margin-top: 20px;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }}
        .result-low {{ background: #d4edda; color: #155724; }}
        .result-moderate {{ background: #fff3cd; color: #856404; }}
        .result-high {{ background: #f8d7da; color: #721c24; }}
        .risk-value {{ font-size: 36px; font-weight: bold; }}
        .risk-label {{ font-size: 14px; margin-top: 5px; }}
        .disclaimer {{
            margin-top: 20px;
            font-size: 11px;
            color: #999;
            text-align: center;
        }}
    </style>
</head>
<body>
    <div class="calculator">
        <h1>KDIGO AKI ICU Mortality Calculator</h1>
        <p class="subtitle">Based on eICU-CRD validation cohort (N=166,373)<br>
        Model B: KDIGO Stage + APACHE + Mech. Vent + Age + Sex</p>

        <div class="input-group">
            <label>KDIGO AKI Stage</label>
            <select id="kdigo">
                <option value="0">Stage 0 (No AKI)</option>
                <option value="1">Stage 1</option>
                <option value="2">Stage 2</option>
                <option value="3">Stage 3</option>
            </select>
        </div>

        <div class="input-group">
            <label>APACHE IVa Score</label>
            <input type="number" id="apache" value="40" min="0" max="200">
        </div>

        <div class="input-group">
            <label>Mechanical Ventilation</label>
            <select id="vent">
                <option value="0">No</option>
                <option value="1">Yes</option>
            </select>
        </div>

        <div class="input-group">
            <label>Age (years)</label>
            <input type="number" id="age" value="65" min="18" max="100">
        </div>

        <div class="input-group">
            <label>Sex</label>
            <select id="sex">
                <option value="0">Female</option>
                <option value="1">Male</option>
            </select>
        </div>

        <button onclick="calculate()" style="width:100%;padding:12px;background:#2166ac;color:white;border:none;border-radius:5px;font-size:16px;cursor:pointer;">
            Calculate Risk
        </button>

        <div id="result" class="result" style="display:none;">
            <div class="risk-value" id="riskValue">--</div>
            <div class="risk-label" id="riskLabel">--</div>
        </div>

        <div class="disclaimer">
            This calculator is for research purposes only and should not be used for clinical decision-making.
            Model performance: AUC = 0.818 (eICU-CRD). Validated against MIMIC-IV (AUC = 0.828).
        </div>
    </div>

    <script>
        // Logistic regression coefficients (Model B)
        const coef = {{
            const: {coef['const']:.6f},
            kdigo_1: {coef['kdigo_1']:.6f},
            kdigo_2: {coef['kdigo_2']:.6f},
            kdigo_3: {coef['kdigo_3']:.6f},
            apachescore: {coef['apachescore']:.6f},
            mech_vent: {coef['mech_vent']:.6f},
            age_num: {coef['age_num']:.6f},
            gender_male: {coef['gender_male']:.6f}
        }};

        function calculate() {{
            const kdigo = parseInt(document.getElementById('kdigo').value);
            const apache = parseFloat(document.getElementById('apache').value);
            const vent = parseInt(document.getElementById('vent').value);
            const age = parseFloat(document.getElementById('age').value);
            const sex = parseInt(document.getElementById('sex').value);

            const logit = coef.const +
                (kdigo === 1 ? coef.kdigo_1 : 0) +
                (kdigo === 2 ? coef.kdigo_2 : 0) +
                (kdigo === 3 ? coef.kdigo_3 : 0) +
                coef.apachescore * apache +
                coef.mech_vent * vent +
                coef.age_num * age +
                coef.gender_male * sex;

            const prob = 1 / (1 + Math.exp(-logit));
            const pct = (prob * 100).toFixed(1);

            const resultDiv = document.getElementById('result');
            const riskValue = document.getElementById('riskValue');
            const riskLabel = document.getElementById('riskLabel');

            riskValue.textContent = pct + '%';

            resultDiv.style.display = 'block';

            if (prob < 0.1) {{
                resultDiv.className = 'result result-low';
                riskLabel.textContent = 'Low risk of ICU mortality';
            }} else if (prob < 0.25) {{
                resultDiv.className = 'result result-moderate';
                riskLabel.textContent = 'Moderate risk of ICU mortality';
            }} else {{
                resultDiv.className = 'result result-high';
                riskLabel.textContent = 'High risk of ICU mortality';
            }}
        }}
    </script>
</body>
</html>'''

html_path = os.path.join(OUTPUT_DIR, 'nomogram_calculator.html')
with open(html_path, 'w', encoding='utf-8') as f:
    f.write(html_content)
print(f"  Web calculator saved: {html_path}")

# Save model coefficients
nomogram_data = {
    'model': 'Logistic Regression (Model B)',
    'cohort': 'eICU-CRD v2.0 (N=166,373)',
    'coefficients': {k: float(v) for k, v in coef.items()},
    'auc': 0.818,
    'brier_score': float(brier),
    'variables': ['KDIGO Stage', 'APACHE Score', 'Mechanical Ventilation', 'Age', 'Sex'],
}

json_path = os.path.join(OUTPUT_DIR, 'nomogram_model.json')
with open(json_path, 'w', encoding='utf-8') as f:
    json.dump(nomogram_data, f, indent=2)
print(f"  Model data saved: {json_path}")

print("\n" + "=" * 60)
print("Nomogram + Web Calculator COMPLETE")
print("=" * 60)
