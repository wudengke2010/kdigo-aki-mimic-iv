"""
PH 假设检验 — 手动 Schoenfeld 残差 + 时间秩相关
对 MIMIC 和 eICU 两个数据库的 Cox Model C
(eICU 走独立子进程避免状态污染)
"""
import json
import os
import subprocess
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
from scipy.stats import spearmanr

WORK = r"C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
PY = r"C:/Users/admin/.workbuddy/binaries/python/envs/default/Scripts/python.exe"

def schoenfeld_ph_test(cph, df, time_col, event_col):
    out = {}
    print(f"    [schoenfeld] computing residuals for n_events={int(df[event_col].sum())}", flush=True)
    resid = cph.compute_residuals(df, kind='schoenfeld')
    print(f"    [schoenfeld] residuals computed, shape={resid.shape}", flush=True)
    for v in resid.columns:
        if v in cph.params_.index:
            rho, p = spearmanr(resid.index.values, resid[v].values)
            out[v] = {'rho': float(rho), 'p': float(p), 'violates': bool(p < 0.05)}
    return out

# ====================== MIMIC ======================
print("=" * 70)
print("MIMIC PH test")
print("=" * 70)
df_m = pd.read_csv(os.path.join(WORK, "mimic_analysis_v2.csv"))
df_m['kdigo_1'] = (df_m['kdigo_stage']==1).astype(int)
df_m['kdigo_2'] = (df_m['kdigo_stage']==2).astype(int)
df_m['kdigo_3'] = (df_m['kdigo_stage']==3).astype(int)

df_m['event_time_h'] = df_m['event_time_h'].clip(lower=0.5)
surv = df_m[df_m['event_time_h'] >= 24].copy()
surv['time_from_lm'] = (surv['event_time_h'] - 24).clip(upper=696)
df_m_use = surv[['time_from_lm','event_ind','kdigo_1','kdigo_2','kdigo_3',
                'age','sex_male','sofa_nonrenal','vent_24h','ckd_history_prior',
                'diabetes','hypertension','heart_failure','copd','liver_disease']].dropna()

print(f"  n={len(df_m_use):,}", flush=True)
cph = CoxPHFitter(penalizer=0.0)
cph.fit(df_m_use, duration_col='time_from_lm', event_col='event_ind')
ph_m = schoenfeld_ph_test(cph, df_m_use, 'time_from_lm', 'event_ind')
print(f"  {'Variable':<22s} {'rho':>8s} {'p':>10s}  violates?")
for v, r in ph_m.items():
    flag = '** YES **' if r['violates'] else 'no'
    print(f"  {v:<22s} {r['rho']:>+8.3f} {r['p']:>10.3e}  {flag}")

# ====================== eICU (separate process) ======================
print("\n" + "=" * 70)
print("eICU PH test (separate process)")
print("=" * 70)
r = subprocess.run([PY, os.path.join(WORK, "ph_test_eicu_v2.py")],
                   capture_output=True, text=True, timeout=600)
print(r.stdout, flush=True)
if r.returncode != 0:
    print(f"STDERR: {r.stderr[-1000:]}", flush=True)
    raise RuntimeError(f"eICU PH failed: rc={r.returncode}")
ph_e = json.load(open(os.path.join(WORK, "v2_outputs", "_ph_eicu.json")))

# ====================== Save to JSON ======================
for j_path, ph in [("survival_mimic.json", ph_m), ("survival_eicu.json", ph_e)]:
    p = os.path.join(WORK, "v2_outputs", j_path)
    d = json.load(open(p))
    d['ph_tests'] = {'C_full_adjusted': ph,
                     'method': 'Spearman correlation of Schoenfeld residuals with event time'}
    json.dump(d, open(p, 'w'), indent=2, default=str)
print(f"\n  PH tests saved into both survival_*.json")