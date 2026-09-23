# -*- coding: utf-8 -*-
"""
v3 Phase 3: AKI phenotype (transient vs persistent) — leakage-free landmark.

Fixes 5.5: phenotype is classified ONLY from data up to onset + 72 h;
risk set starts at the landmark; patients dying or being discharged before
onset + 72 h are excluded (landmark attrition, explicitly reported).

Reference for comparison: eligible AKI patients, transient vs persistent
(adjusted Model C covariates); crude 30-d mortality of no-AKI patients
reported for context.

Outputs: v3_outputs/v3_phenotype_results.json + report lines
"""
import pandas as pd
import numpy as np
import os, json
from datetime import datetime
from lifelines import CoxPHFitter

OUT = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
V3 = os.path.join(OUT, "v3_outputs")

COVARS = ["age", "sex_male", "sofa_nonrenal", "vent_24h",
          "ckd_history_prior", "diabetes", "hypertension",
          "heart_failure", "copd", "liver_disease"]

print("=" * 70)
print(f"v3 Phase 3: phenotype landmark analysis  Start: {datetime.now()}")
print("=" * 70)

results = {}

for tag in ["mimic", "eicu"]:
    pat = pd.read_csv(os.path.join(V3, f"v3_{tag}_patients.csv"))
    r = {}

    # no-AKI context
    noaki = pat[pat["onset_any_h"].isna()]
    r["no_aki_n"] = int(len(noaki))
    r["no_aki_crude_death30_pct"] = round(
        100 * (noaki["event_type"] == 1).mean(), 1)

    # eligible AKI: classified + still at risk at landmark
    aki = pat[pat["onset_any_h"].notna()]
    elig = aki[(aki["phenotype"].isin(["transient", "persistent"]))
               & (aki["end_time"] > aki["landmark_h"])].copy()
    r["aki_n"] = int(len(aki))
    r["n_eligible"] = int(len(elig))
    r["n_excluded_died_before_L"] = int(
        ((aki["event_type"] == 1) & (aki["end_time"] <= aki["landmark_h"])).sum())
    r["n_excluded_discharged_before_L"] = int(
        ((aki["event_type"] == 2) & (aki["end_time"] <= aki["landmark_h"])).sum())
    r["n_excluded_unclassifiable"] = int((aki["phenotype"] == "unknown").sum())

    # follow-up from landmark
    elig["futime"] = elig["end_time"] - elig["landmark_h"]
    elig["fuevent"] = (elig["event_type"] == 1).astype(int)
    elig["persist"] = (elig["phenotype"] == "persistent").astype(int)

    for ph in ["transient", "persistent"]:
        sub = elig[elig["phenotype"] == ph]
        r[f"n_{ph}"] = int(len(sub))
        r[f"crude_death30_{ph}_pct"] = round(100 * sub["fuevent"].mean(), 1)
        r[f"median_followup_{ph}_h"] = round(float(sub["futime"].median()), 1)

    # Cox: persistent vs transient (reference), adjusted
    def run_cox(cols, label):
        d = elig[["futime", "fuevent", "persist"] + cols].dropna()
        cph = CoxPHFitter()
        cph.fit(d, duration_col="futime", event_col="fuevent",
                show_progress=False)
        s = cph.summary.loc["persist"]
        out = {"hr": float(s["exp(coef)"]),
               "lo": float(s["exp(coef) lower 95%"]),
               "hi": float(s["exp(coef) upper 95%"]),
               "p": float(s["p"]), "n": int(len(d)),
               "events": int(d["fuevent"].sum())}
        print(f"  [{tag}] {label}: persistent HR={out['hr']:.2f} "
              f"({out['lo']:.2f}-{out['hi']:.2f}) p={out['p']:.1e} "
              f"n={out['n']:,} deaths={out['events']:,}")
        return out

    r["cox_phenotype_unadjusted"] = run_cox([], "unadjusted")
    r["cox_phenotype_adjusted"] = run_cox(COVARS, "adjusted (Model C)")

    # additional: stage at onset as covariate? stage at landmark = severity
    # (uses only pre-landmark data -> legitimate). Sensitivity:
    elig["s_onset"] = np.nan
    # severity proxy at landmark: stage_max within data<=L is captured by
    # phenotype; use sofa_nonrenal (already first-24h). skip extra.

    results[tag] = r

results["generated"] = datetime.now().isoformat()
with open(os.path.join(V3, "v3_phenotype_results.json"), "w") as f:
    json.dump(results, f, indent=2)

print("\nv3 Phase 3 COMPLETE")
