# -*- coding: utf-8 -*-
"""v5 sensitivity refits under the v3 time-varying framework:
(1) MIMIC Model C excluding mechanical ventilation
(2) eICU Model C with APACHE score replacing non-renal SOFA
Also computes locked-logistic sensitivity excluding ventilation (for completeness).
"""
import os
import json
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

BASE = os.path.dirname(os.path.abspath(__file__))
V3 = os.path.join(BASE, "v3_outputs")

COVARS = ["age", "sex_male", "sofa_nonrenal", "vent_24h", "ckd_history_prior",
          "diabetes", "hypertension", "heart_failure", "copd", "liver_disease"]

out = {}

# ---- (1) MIMIC Model C excluding ventilation ----------------------------
tv = pd.read_csv(os.path.join(V3, "v3_mimic_tv.csv.gz"))
pat = pd.read_csv(os.path.join(V3, "v3_mimic_patients.csv"))
tv = tv.merge(pat[["stay_id"] + COVARS], on="stay_id", how="left")
cols = [c for c in ["stage1", "stage2", "stage3"] + COVARS if c != "vent_24h"]
dd = tv[["t_start", "t_stop", "event_death"] + cols].sort_values("t_stop")
cph = CoxPHFitter()
cph.fit(dd, duration_col="t_stop", event_col="event_death",
        entry_col="t_start", show_progress=False)
s3 = cph.summary.loc["stage3"]
out["mimic_novent"] = {
    "s3_hr": round(float(np.exp(s3["coef"])), 2),
    "lo": round(float(s3["exp(coef) lower 95%"]), 2),
    "hi": round(float(s3["exp(coef) upper 95%"]), 2)}
print("MIMIC no-vent S3:", out["mimic_novent"])

# ---- (2) eICU Model C with APACHE replacing non-renal SOFA --------------
ec = pd.read_csv(os.path.join(BASE, "eicu_cohort_v2.csv"),
                 usecols=["patientunitstayid", "apachescore"])
tve = pd.read_csv(os.path.join(V3, "v3_eicu_tv.csv.gz"))
pate = pd.read_csv(os.path.join(V3, "v3_eicu_patients.csv"))
tve = tve.merge(pate[["patientunitstayid"] + COVARS], on="patientunitstayid", how="left")
tve = tve.merge(ec, on="patientunitstayid", how="left")
print("eICU apache non-null patients:", tve["apachescore"].notna().groupby(
    tve["patientunitstayid"]).max().sum(), "/", tve["patientunitstayid"].nunique())
cols2 = [c for c in ["stage1", "stage2", "stage3"] + COVARS
         if c != "sofa_nonrenal"] + ["apachescore"]
dd2 = tve.dropna(subset=["apachescore"])[
    ["t_start", "t_stop", "event_death"] + cols2].sort_values("t_stop")
cph2 = CoxPHFitter()
cph2.fit(dd2, duration_col="t_stop", event_col="event_death",
         entry_col="t_start", show_progress=False)
s3b = cph2.summary.loc["stage3"]
out["eicu_apache"] = {
    "n_rows": int(len(dd2)),
    "s3_hr": round(float(np.exp(s3b["coef"])), 2),
    "lo": round(float(s3b["exp(coef) lower 95%"]), 2),
    "hi": round(float(s3b["exp(coef) upper 95%"]), 2),
    "apache_hr": round(float(np.exp(cph2.params_["apachescore"])), 3)}
print("eICU apache S3:", out["eicu_apache"])

json.dump(out, open(os.path.join(V3, "v5_sensitivity_refits.json"), "w"), indent=1)
print("DONE")
