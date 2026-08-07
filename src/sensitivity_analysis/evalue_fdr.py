"""
Compute E-values for key ORs (P1-1) and Benjamini-Hochberg FDR for 13 subgroups (P1-7).
"""
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

# ============================================================
# Part 1: E-value calculation (VanderWeele & Ding, 2017)
# ============================================================
# For OR < 1 (protective): convert to 1/OR, then E = RR + sqrt(RR*(RR-1))
# For OR > 1 (harmful): E = OR + sqrt(OR*(OR-1))
# For CI bound: use the bound closest to 1

def evalue_or(or_val, ci_lo, ci_hi):
    """Compute E-value for odds ratio."""
    # For OR, approximate RR conversion (VanderWeele 2017, Table 2)
    # Since outcome is not rare (~30% mortality), use OR directly as approximation
    # E-value formula works on RR; for OR we use the approximation
    if or_val < 1:
        rr = 1 / or_val
        e = rr + np.sqrt(rr * (rr - 1))
        # CI bound closest to 1 is the upper bound
        if ci_hi < 1:
            rr_ci = 1 / ci_hi
            e_ci = rr_ci + np.sqrt(rr_ci * (rr_ci - 1))
        else:
            e_ci = 1.0  # CI crosses 1, E-value for CI bound = 1
    elif or_val > 1:
        rr = or_val
        e = rr + np.sqrt(rr * (rr - 1))
        # CI bound closest to 1 is the lower bound
        if ci_lo > 1:
            rr_ci = ci_lo
            e_ci = rr_ci + np.sqrt(rr_ci * (rr_ci - 1))
        else:
            e_ci = 1.0
    else:
        e = 1.0
        e_ci = 1.0
    return e, e_ci

print("=" * 70)
print("PART 1: E-VALUE CALCULATION")
print("=" * 70)

key_results = [
    ("Stage 3 unadjusted",          0.755, 0.646, 0.883),
    ("Stage 3 SOFA-adjusted",       0.643, 0.537, 0.769),
    ("Stage 3 PSM+SOFA",            0.587, 0.473, 0.728),
    ("Model C: CKD x Stage3",       0.583, 0.477, 0.712),
    ("Model D: CKD x Stage3",       0.679, 0.554, 0.833),
    ("MICE: CKD x Stage3",          0.605, 0.495, 0.741),
]

print(f"\n{'Finding':<35} {'OR':>6} {'E-value':>8} {'E-value (CI)':>13}")
print("-" * 70)
evalue_results = {}
for name, or_val, ci_lo, ci_hi in key_results:
    e, e_ci = evalue_or(or_val, ci_lo, ci_hi)
    print(f"{name:<35} {or_val:>6.3f} {e:>8.2f} {e_ci:>13.2f}")
    evalue_results[name] = {
        "OR": or_val, "CI_lo": ci_lo, "CI_hi": ci_hi,
        "E_value": round(e, 2), "E_value_CI": round(e_ci, 2)
    }

# ============================================================
# Part 2: FDR correction for 13 subgroups (Benjamini-Hochberg)
# ============================================================
print("\n" + "=" * 70)
print("PART 2: BENJAMINI-HOCHBERG FDR CORRECTION (13 SUBGROUPS)")
print("=" * 70)

# Load data and compute subgroup ORs
print("\nLoading data...")
df = pd.read_csv(r'C:\Users\admin\WorkBuddy\2026-07-07-19-40-19\aki_paradox_cohort_with_sofa.csv')

# Filter Stage 3 AKI patients
s3 = df[df['max_aki_stage'] == 3].copy()
print(f"Stage 3 patients: {len(s3)}")

# Define 13 subgroups (matching R script)
subgroups_def = [
    ("Age < 65",         s3[s3['age'] < 65]),
    ("Age >= 65",        s3[s3['age'] >= 65]),
    ("Male",             s3[s3['gender'] == 'M']),
    ("Female",           s3[s3['gender'] == 'F']),
    ("Hypertension",     s3[s3['htn'] == 1]),
    ("No Hypertension",  s3[s3['htn'] == 0]),
    ("Diabetes",         s3[s3['dm'] == 1]),
    ("No Diabetes",      s3[s3['dm'] == 0]),
    ("Heart Failure",    s3[s3['hf'] == 1]),
    ("No Heart Failure", s3[s3['hf'] == 0]),
    ("SOFA Low (0-5)",   s3[s3['sofa_total'] <= 5]),
    ("SOFA Medium (6-10)", s3[(s3['sofa_total'] > 5) & (s3['sofa_total'] <= 10)]),
    ("SOFA High (>10)",  s3[s3['sofa_total'] > 10]),
]

results = []
for name, sub in subgroups_def:
    n = len(sub)
    if n > 10 and sub['ckd'].nunique() > 1 and min(sub['ckd'].value_counts()) > 3:
        X = sm.add_constant(sub[['ckd']])
        y = sub['hospital_expire_flag']
        try:
            model = sm.Logit(y, X).fit(disp=0)
            coef = model.params['ckd']
            se = model.bse['ckd']
            p = model.pvalues['ckd']
            or_val = np.exp(coef)
            ci_lo = np.exp(coef - 1.96 * se)
            ci_hi = np.exp(coef + 1.96 * se)
            results.append({
                'subgroup': name, 'n': n,
                'OR': or_val, 'CI_lo': ci_lo, 'CI_hi': ci_hi,
                'P_value': p
            })
        except Exception as e:
            print(f"  {name}: Error - {e}")

# Apply Benjamini-Hochberg FDR correction
p_values = np.array([r['P_value'] for r in results])
n_tests = len(p_values)
sorted_idx = np.argsort(p_values)
bh_adjusted = np.zeros(n_tests)

for rank, idx in enumerate(sorted_idx):
    bh_adjusted[idx] = p_values[idx] * n_tests / (rank + 1)

# Enforce monotonicity (BH procedure: make sure adjusted p-values are non-decreasing)
for i in range(n_tests - 2, -1, -1):
    sorted_i = sorted_idx[i]
    sorted_i_plus = sorted_idx[i + 1]
    bh_adjusted[sorted_i] = min(bh_adjusted[sorted_i], bh_adjusted[sorted_i_plus])

bh_adjusted = np.minimum(bh_adjusted, 1.0)

# Print results
print(f"\n{'Subgroup':<25} {'n':>6} {'OR':>6} {'95% CI':>16} {'P':>10} {'BH adj P':>10} {'Sig':>5}")
print("-" * 85)
for i, r in enumerate(results):
    ci_str = f"{r['CI_lo']:.2f}-{r['CI_hi']:.2f}"
    p_str = f"{r['P_value']:.4f}" if r['P_value'] >= 0.001 else "<0.001"
    bh_str = f"{bh_adjusted[i]:.4f}" if bh_adjusted[i] >= 0.001 else "<0.001"
    sig = "***" if bh_adjusted[i] < 0.001 else ("**" if bh_adjusted[i] < 0.01 else ("*" if bh_adjusted[i] < 0.05 else "ns"))
    print(f"{r['subgroup']:<25} {r['n']:>6} {r['OR']:>6.2f} {ci_str:>16} {p_str:>10} {bh_str:>10} {sig:>5}")

# Summary
n_sig_before = sum(1 for p in p_values if p < 0.05)
n_sig_after = sum(1 for p in bh_adjusted if p < 0.05)
print(f"\nSignificant at P<0.05: {n_sig_before}/{n_tests} (raw) -> {n_sig_after}/{n_tests} (BH-adjusted)")

# Save results
output = {
    "evalues": evalue_results,
    "fdr": {
        "n_tests": n_tests,
        "subgroups": [
            {
                "name": r['subgroup'],
                "n": r['n'],
                "OR": round(r['OR'], 3),
                "CI_lo": round(r['CI_lo'], 3),
                "CI_hi": round(r['CI_hi'], 3),
                "P_value": r['P_value'],
                "BH_adjusted_P": bh_adjusted[i],
                "significant_005": bool(bh_adjusted[i] < 0.05)
            }
            for i, r in enumerate(results)
        ],
        "n_significant_raw": n_sig_before,
        "n_significant_bh": n_sig_after
    }
}

with open(r'C:\Users\admin\WorkBuddy\2026-07-07-19-40-19\evalue_fdr_results.json', 'w') as f:
    json.dump(output, f, indent=2)

print("\nResults saved to evalue_fdr_results.json")
