#!/usr/bin/env python3
"""
E-value calculation for KDIGO AKI study (Model B).
E-value quantifies the minimum strength of association that an unmeasured
confounder would need to have with both the exposure and the outcome
(jointly, on the risk ratio scale) to fully explain away an observed
exposure-outcome association.

Reference: VanderWeele TJ, Ding P. Sensitivity Analysis in Observational
Research: Introducing the E-Value. Ann Intern Med. 2017;167(4):268-274.
"""

import json
import math
import csv

# --- Input data from existing analysis results ---
# Mortality in reference group (No AKI): 5.98%
P0 = 0.0598

# Model B logistic regression aORs (from logistic_results_revised.csv / summary_revised.json)
logistic_results = {
    "Stage 1": {"OR": 1.885, "CI_low": 1.765, "CI_high": 2.013},
    "Stage 2": {"OR": 2.978, "CI_low": 2.688, "CI_high": 3.300},
    "Stage 3": {"OR": 2.881, "CI_low": 2.664, "CI_high": 3.115},
}

# Model B Cox regression HRs (from summary_revised.json)
cox_results = {
    "Stage 1": {"HR": 1.624, "CI_low": 1.542, "CI_high": 1.711},
    "Stage 2": {"HR": 2.255, "CI_low": 2.088, "CI_high": 2.435},
    "Stage 3": {"HR": 2.037, "CI_low": 1.917, "CI_high": 2.165},
}

# --- E-value formulas (VanderWeele & Ding 2017) ---

def evalue_for_OR(odds_ratio, p0):
    """
    Convert OR to approximate RR using the baseline risk p0,
    then compute E-value.
    RR_approx = OR / (1 - p0 + p0 * OR)
    E-value = RR + sqrt(RR * (RR - 1))  for RR >= 1
    """
    rr = odds_ratio / (1 - p0 + p0 * odds_ratio)
    if rr < 1:
        # For protective effects, use reciprocal
        rr = 1 / rr
    evalue = rr + math.sqrt(rr * (rr - 1))
    return evalue, rr

def evalue_for_HR(hr):
    """
    E-value for hazard ratio (VanderWeele 2017 Appendix, formula for HR):
    E-value = HR^(1/1.91) ... but the simpler approximation for HR is:
    Convert HR to approximate RR: RR = (1 - 0.5^sqrt(HR)) / (1 - 0.5^sqrt(1/HR))
    Then E-value = RR + sqrt(RR * (RR - 1))

    For simplicity and consistency with the common practice for HRs
    in the E-value calculator (evalue-calculator.com), we use the
    HR directly when the outcome is uncommon, or the OR-to-RR conversion
    when common. Since ICU mortality is ~10% overall, we use the
    direct approximation:
    """
    # Using the square root approximation from VanderWeele 2020
    # For HR > 1:
    if hr >= 1:
        rr = math.sqrt(hr)
    else:
        rr = 1 / math.sqrt(1 / hr)

    if rr < 1:
        rr = 1 / rr
    evalue = rr + math.sqrt(rr * (rr - 1))
    return evalue, rr

# --- Calculate E-values ---

results = []

print("=" * 80)
print("E-VALUE ANALYSIS: KDIGO AKI Study (Model B)")
print("Reference: VanderWeele TJ, Ding P. Ann Intern Med. 2017;167(4):268-274")
print("=" * 80)
print()
print(f"Baseline risk (No AKI mortality): {P0*100:.1f}%")
print()

print("--- Logistic Regression (Model B: aOR) ---")
print(f"{'Stage':<12} {'aOR':>8} {'95% CI':>18} {'Approx RR':>10} {'E-value (point)':>16} {'E-value (CI low)':>18}")
print("-" * 82)

for stage, vals in logistic_results.items():
    or_val = vals["OR"]
    ci_low = vals["CI_low"]
    ci_high = vals["CI_high"]

    ev_point, rr_point = evalue_for_OR(or_val, P0)
    ev_ci_low, rr_ci_low = evalue_for_OR(ci_low, P0)

    print(f"{stage:<12} {or_val:>8.3f} {f'({ci_low:.3f}-{ci_high:.3f})':>18} {rr_point:>10.3f} {ev_point:>16.2f} {ev_ci_low:>18.2f}")

    results.append({
        "model": "Logistic (Model B)",
        "exposure": stage,
        "estimate_type": "aOR",
        "estimate": or_val,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "approx_rr": round(rr_point, 3),
        "evalue_point": round(ev_point, 2),
        "evalue_ci_low": round(ev_ci_low, 2),
    })

print()
print("--- Cox Regression (Model B: HR) ---")
print(f"{'Stage':<12} {'HR':>8} {'95% CI':>18} {'Approx RR':>10} {'E-value (point)':>16} {'E-value (CI low)':>18}")
print("-" * 82)

for stage, vals in cox_results.items():
    hr_val = vals["HR"]
    ci_low = vals["CI_low"]
    ci_high = vals["CI_high"]

    ev_point, rr_point = evalue_for_HR(hr_val)
    ev_ci_low, rr_ci_low = evalue_for_HR(ci_low)

    print(f"{stage:<12} {hr_val:>8.3f} {f'({ci_low:.3f}-{ci_high:.3f})':>18} {rr_point:>10.3f} {ev_point:>16.2f} {ev_ci_low:>18.2f}")

    results.append({
        "model": "Cox (Model B)",
        "exposure": stage,
        "estimate_type": "HR",
        "estimate": hr_val,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "approx_rr": round(rr_point, 3),
        "evalue_point": round(ev_point, 2),
        "evalue_ci_low": round(ev_ci_low, 2),
    })

print()
print("=" * 80)
print("INTERPRETATION")
print("=" * 80)
print()
print("The E-value represents the minimum strength of association (on the")
print("risk ratio scale) that an unmeasured confounder would need to have")
print("with BOTH the exposure (KDIGO stage) and the outcome (ICU mortality)")
print("to fully explain away the observed association.")
print()
print("Conventional benchmarks:")
print("  E-value > 1.5 : Moderate robustness to unmeasured confounding")
print("  E-value > 2.0 : Strong robustness")
print("  E-value > 3.0 : Very strong robustness")
print("  E-value > 4.0 : Extremely strong robustness")
print()
print("The E-value for the CI lower bound indicates the confounding strength")
print("needed to shift the CI to include the null (OR=1 or HR=1).")
print()

# Save to JSON
output_json = {
    "description": "E-value analysis for KDIGO AKI study Model B",
    "reference": "VanderWeele TJ, Ding P. Ann Intern Med. 2017;167(4):268-274",
    "baseline_risk_p0": P0,
    "results": results
}

with open("evalue_results.json", "w") as f:
    json.dump(output_json, f, indent=2)

# Save to CSV
with open("evalue_results.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["model", "exposure", "estimate_type", "estimate", "ci_low", "ci_high", "approx_rr", "evalue_point", "evalue_ci_low"])
    writer.writeheader()
    writer.writerows(results)

print("Results saved to: evalue_results.json + evalue_results.csv")
