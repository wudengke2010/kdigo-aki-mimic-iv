# -*- coding: utf-8 -*-
"""
v3 Phase 0: extract FIRST RRT start time per patient (within first 7 ICU days).
Needed for time-varying Stage 3 (RRT onset time, not just a 7-day flag).

Outputs:
  rrt_times_mimic.csv : stay_id, rrt_time_h   (hours from intime)
  rrt_times_eicu.csv  : patientunitstayid, rrt_time_h (hours from unit admission)
"""
import pandas as pd
import numpy as np
import os, gc
import duckdb
from datetime import datetime

OUT = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
MIMIC_ICU = "E:/mimic-iv/v3.1/physionet.org/files/mimiciv/3.1/icu"
EICU = "E:/mimic-iv/eicu/physionet.org/files/eicu-crd/2.0"
WINDOW_RATIO_DAYS = 7

print("=" * 70)
print(f"v3 Phase 0: RRT first-time extraction  Start: {datetime.now()}")
print("=" * 70)

# ---------------- MIMIC ----------------
print("\n[MIMIC] procedureevents dialysis events ...")
d_items = pd.read_csv(os.path.join(MIMIC_ICU, "d_items.csv.gz"),
                      usecols=["itemid", "label"])
dial_items = d_items[d_items["label"].str.lower().str.contains(
    "dialysis - crrt|dialysis - cvvhd|dialysis - cvvhdf|hemodialysis|dialysis - scuf",
    na=False)]
dial_itemids = set(dial_items["itemid"].tolist())
del d_items; gc.collect()

cohort = pd.read_csv(os.path.join(OUT, "kdigo_cohort_v2.csv"),
                     usecols=["stay_id", "intime"], parse_dates=["intime"])

proc = pd.read_csv(os.path.join(MIMIC_ICU, "procedureevents.csv.gz"),
                   usecols=["subject_id", "hadm_id", "stay_id", "starttime", "itemid"],
                   parse_dates=["starttime"])
proc = proc[proc["itemid"].isin(dial_itemids)]
proc = proc[proc["stay_id"].isin(set(cohort["stay_id"]))]
proc = proc.merge(cohort, on="stay_id", how="left")
proc["t_h"] = (proc["starttime"] - proc["intime"]).dt.total_seconds() / 3600.0
proc = proc[(proc["t_h"] >= 0) & (proc["t_h"] <= WINDOW_RATIO_DAYS * 24)]
rrt_mimic = (proc.groupby("stay_id")["t_h"].min()
             .rename("rrt_time_h").reset_index())
print(f"  MIMIC patients with RRT in [0, 7d]: {len(rrt_mimic):,}")
rrt_mimic.to_csv(os.path.join(OUT, "rrt_times_mimic.csv"), index=False)
del proc, cohort; gc.collect()

# ---------------- eICU ----------------
print("\n[eICU] treatment table dialysis events ...")
CSV_OPTS = ", quote='\"', ignore_errors=true"
cohort_e = pd.read_csv(os.path.join(OUT, "eicu_cohort_v2.csv"),
                       usecols=["patientunitstayid"])
ids = tuple(cohort_e["patientunitstayid"].tolist())
rrt_df = duckdb.query(f"""
    SELECT patientunitstayid, treatmentoffset
    FROM read_csv('{EICU}/treatment.csv.gz'{CSV_OPTS})
    WHERE (treatmentstring ILIKE '%dialysis%' OR treatmentstring ILIKE '%crrt%'
           OR treatmentstring ILIKE '%ultrafiltration%')
      AND treatmentstring NOT ILIKE '%catheter%'
      AND treatmentstring NOT ILIKE '%radiology%'
      AND patientunitstayid IN {ids}
""").df()
rrt_df["treatmentoffset"] = pd.to_numeric(rrt_df["treatmentoffset"], errors="coerce")
rrt_df["t_h"] = rrt_df["treatmentoffset"] / 60.0
rrt_df = rrt_df[(rrt_df["t_h"] >= 0) & (rrt_df["t_h"] <= WINDOW_RATIO_DAYS * 24)]
rrt_eicu = (rrt_df.groupby("patientunitstayid")["t_h"].min()
            .rename("rrt_time_h").reset_index())
print(f"  eICU patients with RRT in [0, 7d]: {len(rrt_eicu):,}")
rrt_eicu.to_csv(os.path.join(OUT, "rrt_times_eicu.csv"), index=False)

print("\nDone.")
