"""
PH test for eICU only - separate process to avoid state pollution.
"""
import json
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
from scipy.stats import spearmanr
import os

WORK = r"C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"

df = pd.read_csv(os.path.join(WORK, "eicu_cohort_v2.csv"))
df['age'] = df['age_num']
df['sex_male'] = (df['gender'] == 'Male').astype(int)
df['vent_24h'] = df['mech_vent_day1'].fillna(0).astype(int)
df['ckd_history_prior'] = df['ckd_history'].fillna(0).astype(int)
df['kdigo_1'] = (df['kdigo_stage']==1).astype(int)
df['kdigo_2'] = (df['kdigo_stage']==2).astype(int)
df['kdigo_3'] = (df['kdigo_stage']==3).astype(int)
df['event_time_h'] = df['hospital_los_h'].copy()
df['event_ind'] = df['hospital_mortality'].astype(int)
surv = df[df['event_time_h'] >= 24].copy()
surv['time_from_lm'] = (surv['event_time_h'] - 24).clip(upper=696)
df_use = surv[['time_from_lm','event_ind','kdigo_1','kdigo_2','kdigo_3',
               'age','sex_male','apachescore','vent_24h','ckd_history_prior']].dropna()
print(f"eICU n={len(df_use):,}", flush=True)

cph = CoxPHFitter(penalizer=0.0)
cph.fit(df_use, duration_col='time_from_lm', event_col='event_ind')
print("fit done", flush=True)

resid = cph.compute_residuals(df_use, kind='schoenfeld')
print(f"resid shape: {resid.shape}", flush=True)

out = {}
for v in resid.columns:
    if v in cph.params_.index:
        rho, p = spearmanr(resid.index.values, resid[v].values)
        out[v] = {'rho': float(rho), 'p': float(p), 'violates': bool(p < 0.05)}
        print(f"  {v}: rho={rho:+.3f} p={p:.3e}", flush=True)

# Save to JSON for parent process
with open(os.path.join(WORK, "v2_outputs", "_ph_eicu.json"), 'w') as f:
    json.dump(out, f, indent=2)
print("saved", flush=True)