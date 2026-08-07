#!/usr/bin/env python3
"""
AKI Paradox: Fix UO staging and combine results
Run this if aki_paradox_analysis_v2.py crashes at Step 6 (UO staging).
Requires intermediate_creatinine.csv and the cohort demographics from v2.
"""

import pandas as pd
import numpy as np
import os, sys, gc, warnings
warnings.filterwarnings('ignore')

DATA_DIR = "E:/mimic-iv/v3.1/physionet.org/files/mimiciv/3.1"
HOSP_DIR = os.path.join(DATA_DIR, "hosp")
ICU_DIR = os.path.join(DATA_DIR, "icu")
OUTPUT_DIR = "C:/Users/admin/WorkBuddy/2026-07-07-19-40-19"

UO_ITEMIDS = [226559, 226560, 226561, 226584, 226563, 226564, 226565,
              226567, 226557, 226558, 227488, 227489]

def log(msg):
    print(msg, flush=True)

log("=" * 70)
log("AKI Paradox: Fix UO + Combine Results")
log("=" * 70)

# ============================================================
# Reload demographics (same as v2 Steps 1-3)
# ============================================================
log("\n[Reload] Demographics...")
patients = pd.read_csv(os.path.join(HOSP_DIR, "patients.csv.gz"), compression="gzip",
    usecols=["subject_id", "gender", "anchor_age", "anchor_year", "dod"])
admissions = pd.read_csv(os.path.join(HOSP_DIR, "admissions.csv.gz"), compression="gzip",
    usecols=["subject_id", "hadm_id", "admittime", "dischtime", "deathtime",
             "hospital_expire_flag", "admission_type"])
admissions["admittime"] = pd.to_datetime(admissions["admittime"])
admissions["dischtime"] = pd.to_datetime(admissions["dischtime"])
admissions["deathtime"] = pd.to_datetime(admissions["deathtime"])
icustays = pd.read_csv(os.path.join(ICU_DIR, "icustays.csv.gz"), compression="gzip",
    usecols=["subject_id", "hadm_id", "stay_id", "first_careunit", "intime", "outtime", "los"])
icustays["intime"] = pd.to_datetime(icustays["intime"])
icustays["outtime"] = pd.to_datetime(icustays["outtime"])

demo = icustays.merge(admissions, on=["subject_id", "hadm_id"]).merge(patients, on="subject_id")
demo["age"] = demo["anchor_age"] + (demo["admittime"].dt.year - demo["anchor_year"])
demo = demo[demo["age"] >= 18].copy()
demo = demo.sort_values("intime").drop_duplicates("subject_id", keep="first")

# ICD diagnoses
log("[Reload] ICD diagnoses...")
diag = pd.read_csv(os.path.join(HOSP_DIR, "diagnoses_icd.csv.gz"), compression="gzip",
    usecols=["hadm_id", "icd_code", "icd_version"])

def flag_icd(df, prefixes_v9, prefixes_v10):
    mask = ((df["icd_code"].str[:len(prefixes_v9[0])].isin(prefixes_v9)) & (df["icd_version"] == 9)) | \
           ((df["icd_code"].str[:len(prefixes_v10[0])].isin(prefixes_v10)) & (df["icd_version"] == 10))
    return set(df.loc[mask, "hadm_id"].unique())

demo["ckd"] = demo["hadm_id"].isin(flag_icd(diag, ["585"], ["N18"])).astype(int)
demo["esrd"] = demo["hadm_id"].isin(flag_icd(diag, ["V451", "9968"], ["Z992", "T824"])).astype(int)
demo["htn"] = demo["hadm_id"].isin(flag_icd(diag, ["401", "402", "403", "404", "405"], ["I10", "I11", "I12", "I13", "I15"])).astype(int)
demo["dm"] = demo["hadm_id"].isin(flag_icd(diag, ["250"], ["E10", "E11", "E12", "E13", "E14"])).astype(int)
demo["hf"] = demo["hadm_id"].isin(flag_icd(diag, ["4280", "4281", "4282", "4283", "4284", "4289"], ["I50"])).astype(int)
del diag; gc.collect()

# Weights from OMR
log("[Reload] Weights...")
omr = pd.read_csv(os.path.join(HOSP_DIR, "omr.csv.gz"), compression="gzip")
wt = omr[omr["result_name"] == "Weight (Lbs)"].copy()
wt["result_value"] = pd.to_numeric(wt["result_value"], errors="coerce")
wt = wt[(wt["result_value"] > 50) & (wt["result_value"] < 1000)]
wt["weight"] = wt["result_value"] * 0.453592
wt = wt.sort_values("chartdate").drop_duplicates("subject_id", keep="last")[["subject_id", "weight"]]
demo = demo.merge(wt, on="subject_id", how="left")
demo["weight"] = demo["weight"].fillna(70.0)
del omr, wt; gc.collect()

# ============================================================
# Load creatinine staging results
# ============================================================
log("\n[Load] Creatinine data...")
creat_path = os.path.join(OUTPUT_DIR, "intermediate_creatinine.csv")
if os.path.exists(creat_path):
    creat_df = pd.read_csv(creat_path)
    creat_df["charttime"] = pd.to_datetime(creat_df["charttime"])
    log(f"  Loaded {len(creat_df):,} creatinine values for {creat_df['stay_id'].nunique():,} stays")
else:
    log("  ERROR: intermediate_creatinine.csv not found!")
    sys.exit(1)

# ============================================================
# KDIGO Creatinine Staging (using groupby for speed)
# ============================================================
log("\n[Step 5] KDIGO creatinine staging (groupby)...")
stay_ids_list = creat_df["stay_id"].unique()
log(f"  Processing {len(stay_ids_list):,} stays...")

aki_results = []
processed = 0
for stay_id, grp in creat_df.groupby("stay_id"):
    grp = grp.sort_values("charttime")
    n = len(grp)
    if n == 0:
        continue
    times = grp["charttime"].values
    creats = grp["creatinine"].values
    max_stage = 0

    for i in range(n):
        t_prev_7d = times[i] - np.timedelta64(7, 'D')
        t_prev_48h = times[i] - np.timedelta64(48, 'h')
        mask_7d = (times[:i] >= t_prev_7d) & (times[:i] < times[i])
        mask_48h = (times[:i] >= t_prev_48h) & (times[:i] < times[i])

        cl7 = creats[:i][mask_7d].min() if mask_7d.any() else creats[i]
        cl48 = creats[:i][mask_48h].min() if mask_48h.any() else creats[i]
        ratio = creats[i] / cl7 if cl7 > 0 else 1.0

        stage = 0
        if creats[i] >= cl7 * 3.0:
            stage = 3
        elif creats[i] >= 4.0 and (cl48 <= 3.7 or ratio >= 1.5):
            stage = 3
        elif creats[i] >= cl7 * 2.0:
            stage = 2
        elif creats[i] >= cl48 + 0.3:
            stage = 1
        elif ratio >= 1.5:
            stage = 1

        if stage > max_stage:
            max_stage = stage

    aki_results.append({"stay_id": stay_id, "max_aki_stage_creat": int(max_stage)})
    processed += 1
    if processed % 10000 == 0:
        log(f"  Processed {processed:,}/{len(stay_ids_list):,} stays...")

max_aki_creat = pd.DataFrame(aki_results)
log(f"  Done. Stages: " + ", ".join(
    f"S{s}={(max_aki_creat['max_aki_stage_creat']==s).sum()}" for s in range(4)))

del creat_df; gc.collect()

# ============================================================
# Urine Output - FIXED version
# ============================================================
log("\n[Step 6] UO extraction + staging (fixed)...")
oe = pd.read_csv(os.path.join(ICU_DIR, "outputevents.csv.gz"), compression="gzip",
    usecols=["stay_id", "charttime", "itemid", "value"],
    dtype={"stay_id": "float64", "itemid": "int32", "value": "float64"})
stay_ids_set = set(demo["stay_id"].values)
uo = oe[oe["itemid"].isin(UO_ITEMIDS) & oe["stay_id"].isin(stay_ids_set)].copy()
uo["charttime"] = pd.to_datetime(uo["charttime"])
uo.loc[uo["itemid"] == 227488, "value"] = -uo.loc[uo["itemid"] == 227488, "value"].abs()
uo = uo.groupby(["stay_id", "charttime"], as_index=False)["value"].sum()
uo = uo.rename(columns={"value": "uo"})
del oe; gc.collect()

wt_lookup = demo[["stay_id", "weight", "intime"]].copy()
uo = uo.merge(wt_lookup, on="stay_id", how="inner")
log(f"  UO records: {len(uo):,}, stays: {uo['stay_id'].nunique():,}")

aki_uo_results = []
processed = 0
for stay_id, grp in uo.groupby("stay_id"):
    grp = grp.sort_values("charttime")
    intime = grp["intime"].iloc[0]
    weight = max(grp["weight"].iloc[0], 1.0)
    times = grp["charttime"].values
    vals = grp["uo"].values
    n = len(times)
    max_stage = 0

    for i in range(n):
        # Skip first 6 hours
        seconds_since = (times[i] - np.datetime64(intime)).astype('timedelta64[s]').astype(np.float64)
        if seconds_since < 6 * 3600:
            continue

        # 6-hour window
        t6 = times[i] - np.timedelta64(6, 'h')
        m6 = (times[:i+1] > t6) & (times[:i+1] <= times[i])
        h6 = m6.sum()
        if h6 >= 6:
            rate6 = vals[:i+1][m6].sum() / weight / 6
            if rate6 < 0.5:
                max_stage = max(max_stage, 1)

        # 12-hour window
        t12 = times[i] - np.timedelta64(12, 'h')
        m12 = (times[:i+1] > t12) & (times[:i+1] <= times[i])
        h12 = m12.sum()
        if h12 >= 12:
            rate12 = vals[:i+1][m12].sum() / weight / 12
            if rate12 == 0:
                max_stage = max(max_stage, 3)
            elif rate12 < 0.5:
                max_stage = max(max_stage, 2)

        # 24-hour window
        t24 = times[i] - np.timedelta64(24, 'h')
        m24 = (times[:i+1] > t24) & (times[:i+1] <= times[i])
        h24 = m24.sum()
        if h24 >= 24:
            rate24 = vals[:i+1][m24].sum() / weight / 24
            if rate24 < 0.3:
                max_stage = max(max_stage, 3)

    aki_uo_results.append({"stay_id": stay_id, "max_aki_stage_uo": int(max_stage)})
    processed += 1
    if processed % 10000 == 0:
        log(f"  UO: Processed {processed:,} stays...")

max_aki_uo = pd.DataFrame(aki_uo_results)
log(f"  UO staging done. Stages: " + ", ".join(
    f"S{s}={(max_aki_uo['max_aki_stage_uo']==s).sum()}" for s in range(4)))
del uo; gc.collect()

# ============================================================
# Step 7: Combine
# ============================================================
log("\n[Step 7] Combining results...")
cohort = demo[["stay_id", "hadm_id", "subject_id", "age", "gender",
               "hospital_expire_flag", "ckd", "esrd", "htn", "dm", "hf",
               "los", "first_careunit", "admission_type", "weight",
               "intime", "outtime"]].copy()

cohort = cohort.merge(max_aki_creat, on="stay_id", how="left")
cohort["max_aki_stage_creat"] = cohort["max_aki_stage_creat"].fillna(0).astype(int)
cohort = cohort.merge(max_aki_uo, on="stay_id", how="left")
cohort["max_aki_stage_uo"] = cohort["max_aki_stage_uo"].fillna(0).astype(int)
cohort["max_aki_stage"] = cohort[["max_aki_stage_creat", "max_aki_stage_uo"]].max(axis=1)

# Exclude ESRD
cohort = cohort[cohort["esrd"] == 0].copy()

# Mortality
dod = patients[["subject_id", "dod"]].copy()
dod["dod"] = pd.to_datetime(dod["dod"])
cohort = cohort.merge(dod, on="subject_id", how="left")
cohort["intime"] = pd.to_datetime(cohort["intime"])
cohort["mortality_30d"] = ((cohort["dod"].notna()) &
    (cohort["dod"] <= cohort["intime"] + pd.Timedelta(days=30))).astype(int)
cohort["mortality_90d"] = ((cohort["dod"].notna()) &
    (cohort["dod"] <= cohort["intime"] + pd.Timedelta(days=90))).astype(int)
cohort["mortality_1y"] = ((cohort["dod"].notna()) &
    (cohort["dod"] <= cohort["intime"] + pd.Timedelta(days=365))).astype(int)

# Groups
aki_pts = cohort[cohort["max_aki_stage"] >= 1]
pure_aki = aki_pts[aki_pts["ckd"] == 0]
aki_ckd = aki_pts[aki_pts["ckd"] == 1]

log(f"\n  Total cohort (excl ESRD): {len(cohort):,}")
log(f"  AKI patients: {len(aki_pts):,}")
log(f"  Pure AKI: {len(pure_aki):,}")
log(f"  AKI+CKD: {len(aki_ckd):,}")

log(f"\n  AKI stages:")
for s in range(4):
    n = (cohort["max_aki_stage"] == s).sum()
    log(f"    Stage {s}: {n:,} ({n/len(cohort)*100:.1f}%)")

log(f"\n  Hospital mortality:")
if len(pure_aki) > 0:
    log(f"    Pure AKI: {pure_aki['hospital_expire_flag'].sum()}/{len(pure_aki)} ({pure_aki['hospital_expire_flag'].mean()*100:.1f}%)")
if len(aki_ckd) > 0:
    log(f"    AKI+CKD:  {aki_ckd['hospital_expire_flag'].sum()}/{len(aki_ckd)} ({aki_ckd['hospital_expire_flag'].mean()*100:.1f}%)")

cohort.to_csv(os.path.join(OUTPUT_DIR, "aki_paradox_cohort.csv"), index=False)
log(f"\n  Saved aki_paradox_cohort.csv ({len(cohort):,} rows)")

log("\n" + "=" * 70)
log("Pipeline complete!")
log("=" * 70)
