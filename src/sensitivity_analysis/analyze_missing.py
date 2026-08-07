#!/usr/bin/env python3
"""Analyze missing data patterns in AKI subset for MICE."""
import pandas as pd
import numpy as np

df = pd.read_csv(r'C:\Users\admin\WorkBuddy\2026-07-07-19-40-19\aki_paradox_cohort_with_sofa.csv')

# AKI subset: max_aki_stage >= 1
aki_df = df[df['max_aki_stage'] >= 1].copy()
print(f"Total cohort: {len(df)}")
print(f"AKI subset (max_aki_stage >= 1): {len(aki_df)}")
print(f"AKI stage distribution: {aki_df['max_aki_stage'].value_counts().sort_index().to_dict()}")
print(f"CKD distribution: {aki_df['ckd'].value_counts().to_dict()}")
print()

# All variables used in models
model_vars = [
    # Demographics
    'age', 'gender',
    # Comorbidities
    'ckd', 'esrd', 'htn', 'dm', 'hf',
    # AKI
    'max_aki_stage',
    # SOFA components (used in Model D)
    'sofa_resp', 'sofa_coag', 'sofa_liver', 'sofa_cv', 'sofa_cns', 'sofa_renal',
    'sofa_total', 'sofa_nonrenal',
    # Outcomes
    'hospital_expire_flag', 'mortality_30d', 'mortality_90d', 'mortality_1y',
    # Underlying labs (source of missingness)
    'bili_max', 'pao2fio2_min', 'platelet_min', 'gcs_min', 'map_min',
    'norepi_max_rate', 'creat_max_24h'
]

print("=== Missing analysis (AKI subset, n={}) ===".format(len(aki_df)))
print(f"{'Variable':<20} {'Missing':>8} {'%':>8}   {'CKD missing %':>14} {'Non-CKD missing %':>18}")
print("-" * 75)

for col in model_vars:
    n_missing = aki_df[col].isna().sum()
    pct = n_missing / len(aki_df) * 100
    ckd_missing = aki_df[aki_df['ckd'] == 1][col].isna().sum() / max(1, (aki_df['ckd'] == 1).sum()) * 100
    nonckd_missing = aki_df[aki_df['ckd'] == 0][col].isna().sum() / max(1, (aki_df['ckd'] == 0).sum()) * 100
    flag = " ***" if pct > 5 else ""
    print(f"{col:<20} {n_missing:>8} {pct:>7.1f}%   {ckd_missing:>13.1f}% {nonckd_missing:>17.1f}%{flag}")

# Check SOFA component zeros - are they "true zero" or "missing imputed as zero"?
print("\n=== SOFA component distributions (AKI subset) ===")
sofa_components = ['sofa_resp', 'sofa_coag', 'sofa_liver', 'sofa_cv', 'sofa_cns', 'sofa_renal']
for comp in sofa_components:
    vals = aki_df[comp].dropna()
    zero_pct = (vals == 0).sum() / len(vals) * 100
    print(f"  {comp}: n={len(vals)}, zero={zero_pct:.1f}%, mean={vals.mean():.2f}, range=[{vals.min()}, {vals.max()}]")

# Cross-tab: missing bili_max vs sofa_liver
print("\n=== Cross-tab: bili_max missing vs sofa_liver ===")
aki_df['bili_missing'] = aki_df['bili_max'].isna()
ct = pd.crosstab(aki_df['bili_missing'], aki_df['sofa_liver'], margins=True)
print(ct)

print("\n=== Cross-tab: pao2fio2_min missing vs sofa_resp ===")
aki_df['pao2_missing'] = aki_df['pao2fio2_min'].isna()
ct2 = pd.crosstab(aki_df['pao2_missing'], aki_df['sofa_resp'], margins=True)
print(ct2)

# Key question: are SOFA components already complete (no NaN)?
print("\n=== SOFA components NaN check ===")
for comp in sofa_components + ['sofa_total', 'sofa_nonrenal']:
    n_nan = aki_df[comp].isna().sum()
    print(f"  {comp}: NaN={n_nan}")

# Stage 3 subset analysis
print("\n=== Stage 3 AKI subset (n={}) ===".format((aki_df['max_aki_stage'] == 3).sum()))
s3 = aki_df[aki_df['max_aki_stage'] == 3]
for col in ['bili_max', 'pao2fio2_min', 'platelet_min', 'gcs_min', 'map_min']:
    n_miss = s3[col].isna().sum()
    pct = n_miss / len(s3) * 100
    ckd_pct = s3[s3['ckd'] == 1][col].isna().sum() / max(1, (s3['ckd'] == 1).sum()) * 100
    pure_pct = s3[s3['ckd'] == 0][col].isna().sum() / max(1, (s3['ckd'] == 0).sum()) * 100
    print(f"  {col}: {pct:.1f}% missing (CKD: {ckd_pct:.1f}%, Pure: {pure_pct:.1f}%)")
