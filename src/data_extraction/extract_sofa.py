#!/usr/bin/env python3
"""
SOFA Score Extraction from MIMIC-IV v3.1
Extracts first 24h SOFA components for AKI patients.
Optimized with vectorized groupby operations.
"""

import duckdb
import pandas as pd
import numpy as np
import os
import time

OUTPUT_DIR = "C:/Users/admin/WorkBuddy/2026-07-07-19-40-19"
DATA_DIR = "E:/mimic-iv/v3.1/physionet.org/files/mimiciv/3.1"
ICU_DIR = f"{DATA_DIR}/icu"
HOSP_DIR = f"{DATA_DIR}/hosp"

def log(msg):
    t = time.strftime("%H:%M:%S")
    print(f"[{t}] {msg}", flush=True)

log("=" * 70)
log("SOFA Score Extraction from MIMIC-IV v3.1")
log("=" * 70)

# ============================================================
# Load cohort
# ============================================================
log("Loading cohort...")
df = pd.read_csv(os.path.join(OUTPUT_DIR, "aki_paradox_cohort.csv"))
aki = df[df["max_aki_stage"] >= 1].copy()
aki["intime"] = pd.to_datetime(aki["intime"])
log(f"  Total cohort: {len(df):,}, AKI patients: {len(aki):,}")

aki_subset = aki[["stay_id", "hadm_id", "intime"]].copy()
aki_subset["intime_str"] = aki_subset["intime"].dt.strftime("%Y-%m-%d %H:%M:%S")
aki_subset[["stay_id", "hadm_id", "intime_str"]].to_csv(
    os.path.join(OUTPUT_DIR, "aki_subset_for_sofa.csv"), index=False
)

con = duckdb.connect()
con.execute(f"""
    CREATE TABLE cohort AS 
    SELECT stay_id, hadm_id, CAST(intime_str AS TIMESTAMP) as intime
    FROM read_csv_auto('{OUTPUT_DIR}/aki_subset_for_sofa.csv')
""")
log(f"  Registered {len(aki_subset):,} AKI patients in DuckDB")

# ============================================================
# Step 1: Query chartevents (GCS, MAP, FiO2)
# ============================================================
log("\nStep 1: Querying chartevents (3.3GB)...")
t0 = time.time()

chartevents = con.execute(f"""
    SELECT ce.stay_id, ce.itemid, ce.charttime, ce.valuenum
    FROM read_csv_auto('{ICU_DIR}/chartevents.csv.gz') ce
    WHERE ce.itemid IN (
        220739, 223900, 223901,
        184, 454, 467, 654,
        220052, 220181, 225312,
        223835
      )
      AND ce.stay_id IN (SELECT stay_id FROM cohort)
""").fetchdf()

log(f"  Chartevents: {len(chartevents):,} rows in {time.time()-t0:.1f}s")

chartevents["charttime"] = pd.to_datetime(chartevents["charttime"])
intime_map = aki.set_index("stay_id")["intime"].to_dict()
chartevents["intime"] = chartevents["stay_id"].map(intime_map)
chartevents_24h = chartevents[
    (chartevents["charttime"] >= chartevents["intime"]) &
    (chartevents["charttime"] < chartevents["intime"] + pd.Timedelta(hours=24))
].copy()
log(f"  After 24h filter: {len(chartevents_24h):,} rows")

# ============================================================
# Step 2: Query labevents (platelets, bilirubin, PaO2, BUN)
# ============================================================
log("\nStep 2: Querying labevents (2.5GB)...")
t0 = time.time()

labevents = con.execute(f"""
    SELECT le.hadm_id, le.itemid, le.charttime, le.valuenum
    FROM read_csv_auto('{HOSP_DIR}/labevents.csv.gz') le
    WHERE le.itemid IN (
        51265, 50885, 50821, 50816, 51222
      )
      AND le.hadm_id IN (SELECT hadm_id FROM cohort)
""").fetchdf()

log(f"  Labevents: {len(labevents):,} rows in {time.time()-t0:.1f}s")

labevents["charttime"] = pd.to_datetime(labevents["charttime"])
hadm_to_stay = aki.set_index("hadm_id")["stay_id"].to_dict()
hadm_to_intime = aki.set_index("hadm_id")["intime"].to_dict()
labevents["stay_id"] = labevents["hadm_id"].map(hadm_to_stay)
labevents["intime"] = labevents["hadm_id"].map(hadm_to_intime)
labevents_24h = labevents[
    (labevents["charttime"] >= labevents["intime"]) &
    (labevents["charttime"] < labevents["intime"] + pd.Timedelta(hours=24))
].copy()
log(f"  After 24h filter: {len(labevents_24h):,} rows")

# ============================================================
# Step 3: Query inputevents (vasopressors)
# ============================================================
log("\nStep 3: Querying inputevents (383MB)...")
t0 = time.time()

inputevents = con.execute(f"""
    SELECT ie.stay_id, ie.itemid, ie.starttime, ie.endtime,
           ie.amount, ie.rate, ie.rateuom, ie.patientweight
    FROM read_csv_auto('{ICU_DIR}/inputevents.csv.gz') ie
    WHERE ie.itemid IN (
        221906, 221289, 221662, 221653, 222315
      )
      AND ie.stay_id IN (SELECT stay_id FROM cohort)
""").fetchdf()

log(f"  Inputevents: {len(inputevents):,} rows in {time.time()-t0:.1f}s")

inputevents["starttime"] = pd.to_datetime(inputevents["starttime"])
inputevents["intime"] = inputevents["stay_id"].map(intime_map)
inputevents_24h = inputevents[
    (inputevents["starttime"] >= inputevents["intime"]) &
    (inputevents["starttime"] < inputevents["intime"] + pd.Timedelta(hours=24))
].copy()
log(f"  After 24h filter: {len(inputevents_24h):,} rows")

# ============================================================
# Step 4: Vectorized SOFA component calculation
# ============================================================
log("\nStep 4: Calculating SOFA components (vectorized)...")

all_stay_ids = aki["stay_id"].unique()
sofa_df = pd.DataFrame({"stay_id": all_stay_ids})

# --- 4a: GCS ---
log("  GCS...")
gcs_comp = chartevents_24h[chartevents_24h["itemid"].isin(
    [220739, 223900, 223901, 454, 467, 654])].copy()
gcs_legacy = chartevents_24h[chartevents_24h["itemid"] == 184]

if len(gcs_comp) > 0:
    gcs_comp["component"] = gcs_comp["itemid"].map({
        220739: "eye", 454: "eye",
        223900: "verbal", 467: "verbal",
        223901: "motor", 654: "motor"
    })
    gcs_pivot = gcs_comp.pivot_table(
        index=["stay_id", "charttime"],
        columns="component",
        values="valuenum",
        aggfunc="min"
    )
    if all(c in gcs_pivot.columns for c in ["eye", "verbal", "motor"]):
        gcs_pivot["gcs_total"] = gcs_pivot[["eye", "verbal", "motor"]].sum(axis=1, min_count=1)
        gcs_min_computed = gcs_pivot.groupby("stay_id")["gcs_total"].min()
    else:
        gcs_min_computed = pd.Series(dtype=float)
else:
    gcs_min_computed = pd.Series(dtype=float)

gcs_min_legacy = gcs_legacy.groupby("stay_id")["valuenum"].min() if len(gcs_legacy) > 0 else pd.Series(dtype=float)

gcs_min = gcs_min_computed.to_dict()
for stay_id, val in gcs_min_legacy.items():
    if stay_id not in gcs_min or pd.isna(gcs_min[stay_id]):
        gcs_min[stay_id] = val

sofa_df["gcs_min"] = sofa_df["stay_id"].map(gcs_min)

# --- 4b: MAP ---
log("  MAP...")
map_data = chartevents_24h[chartevents_24h["itemid"].isin([220052, 220181, 225312])]
if len(map_data) > 0:
    map_min = map_data.groupby("stay_id")["valuenum"].min()
    sofa_df["map_min"] = sofa_df["stay_id"].map(map_min)
else:
    sofa_df["map_min"] = np.nan

# --- 4c: FiO2 ---
log("  FiO2...")
fio2_data = chartevents_24h[chartevents_24h["itemid"] == 223835]
if len(fio2_data) > 0:
    fio2_max = fio2_data.groupby("stay_id")["valuenum"].max()
    sofa_df["fio2_max"] = sofa_df["stay_id"].map(fio2_max)
else:
    sofa_df["fio2_max"] = np.nan

# --- 4d: PaO2 ---
log("  PaO2...")
pao2_data = labevents_24h[labevents_24h["itemid"] == 50821]
if len(pao2_data) > 0:
    pao2_min = pao2_data.groupby("stay_id")["valuenum"].min()
    sofa_df["pao2_min"] = sofa_df["stay_id"].map(pao2_min)
else:
    sofa_df["pao2_min"] = np.nan

# --- 4e: PaO2/FiO2 ratio ---
log("  PaO2/FiO2 ratio...")
def calc_pafi(row):
    pao2 = row["pao2_min"]
    fio2 = row["fio2_max"]
    if pd.isna(pao2) or pd.isna(fio2) or fio2 == 0:
        return np.nan
    fio2_pct = fio2 if fio2 > 1 else fio2 * 100
    return pao2 / fio2_pct * 100

sofa_df["pao2fio2_min"] = sofa_df.apply(calc_pafi, axis=1)

# --- 4f: Platelets ---
log("  Platelets...")
plt_data = labevents_24h[labevents_24h["itemid"] == 51265]
if len(plt_data) > 0:
    plt_min = plt_data.groupby("stay_id")["valuenum"].min()
    sofa_df["platelet_min"] = sofa_df["stay_id"].map(plt_min)
else:
    sofa_df["platelet_min"] = np.nan

# --- 4g: Bilirubin ---
log("  Bilirubin...")
bili_data = labevents_24h[labevents_24h["itemid"] == 50885]
if len(bili_data) > 0:
    bili_max = bili_data.groupby("stay_id")["valuenum"].max()
    sofa_df["bili_max"] = sofa_df["stay_id"].map(bili_max)
else:
    sofa_df["bili_max"] = np.nan

# --- 4h: Vasopressors ---
log("  Vasopressors...")
vaso_rates = {}
for stay_id in all_stay_ids:
    vaso_rates[stay_id] = {"norepi": 0, "epi": 0, "dopa": 0, "doba": 0}

if len(inputevents_24h) > 0:
    vaso_pivot = inputevents_24h.groupby(["stay_id", "itemid"])["rate"].max()
    
    vaso_itemid_map = {
        221906: "norepi",
        221289: "epi",
        221662: "dopa",
        221653: "doba"
    }
    
    for (stay_id, itemid), rate in vaso_pivot.items():
        if itemid in vaso_itemid_map and not pd.isna(rate):
            comp = vaso_itemid_map[itemid]
            if stay_id in vaso_rates:
                vaso_rates[stay_id][comp] = max(vaso_rates[stay_id][comp], rate)

sofa_df["norepi_max_rate"] = sofa_df["stay_id"].map(lambda s: vaso_rates.get(s, {}).get("norepi", 0))
sofa_df["epi_max_rate"] = sofa_df["stay_id"].map(lambda s: vaso_rates.get(s, {}).get("epi", 0))
sofa_df["dopa_max_rate"] = sofa_df["stay_id"].map(lambda s: vaso_rates.get(s, {}).get("dopa", 0))
sofa_df["doba_max_rate"] = sofa_df["stay_id"].map(lambda s: vaso_rates.get(s, {}).get("doba", 0))

# --- 4i: Creatinine max in 24h ---
log("  Creatinine (24h max)...")
creat_df = pd.read_csv(os.path.join(OUTPUT_DIR, "intermediate_creatinine.csv"))
creat_df["charttime"] = pd.to_datetime(creat_df["charttime"])
creat_df = creat_df.merge(aki[["stay_id", "intime"]], on="stay_id")
creat_24h = creat_df[
    (creat_df["charttime"] >= creat_df["intime"]) &
    (creat_df["charttime"] < creat_df["intime"] + pd.Timedelta(hours=24))
]
creat_max_24h = creat_24h.groupby("stay_id")["creatinine"].max()
sofa_df["creat_max_24h"] = sofa_df["stay_id"].map(creat_max_24h)

# ============================================================
# Step 5: Calculate SOFA scores
# ============================================================
log("\nStep 5: Calculating SOFA scores...")

def sofa_respiratory(pao2fio2):
    if pd.isna(pao2fio2):
        return 0
    if pao2fio2 >= 400: return 0
    elif pao2fio2 >= 300: return 1
    elif pao2fio2 >= 200: return 2
    elif pao2fio2 >= 100: return 3
    else: return 4

def sofa_coagulation(plt):
    if pd.isna(plt): return 0
    if plt >= 150: return 0
    elif plt >= 100: return 1
    elif plt >= 50: return 2
    elif plt >= 20: return 3
    else: return 4

def sofa_liver(bili):
    if pd.isna(bili): return 0
    if bili < 1.2: return 0
    elif bili < 2.0: return 1
    elif bili < 6.0: return 2
    elif bili < 12.0: return 3
    else: return 4

def sofa_cv(row):
    norepi = row["norepi_max_rate"]
    epi = row["epi_max_rate"]
    dopa = row["dopa_max_rate"]
    doba = row["doba_max_rate"]
    map_min = row["map_min"]
    
    if dopa > 15 or epi > 0.1 or norepi > 0.1:
        return 4
    elif dopa > 5 or (0 < epi <= 0.1) or (0 < norepi <= 0.1):
        return 3
    elif dopa > 0 or doba > 0:
        return 2
    elif not pd.isna(map_min) and map_min < 70:
        return 1
    else:
        return 0

def sofa_cns(gcs):
    if pd.isna(gcs): return 0
    if gcs >= 15: return 0
    elif gcs >= 13: return 1
    elif gcs >= 10: return 2
    elif gcs >= 6: return 3
    else: return 4

def sofa_renal(creat):
    if pd.isna(creat): return 0
    if creat >= 5.0: return 4
    elif creat >= 3.5: return 3
    elif creat >= 2.0: return 2
    elif creat >= 1.2: return 1
    else: return 0

sofa_df["sofa_resp"] = sofa_df["pao2fio2_min"].apply(sofa_respiratory)
sofa_df["sofa_coag"] = sofa_df["platelet_min"].apply(sofa_coagulation)
sofa_df["sofa_liver"] = sofa_df["bili_max"].apply(sofa_liver)
sofa_df["sofa_cv"] = sofa_df.apply(sofa_cv, axis=1)
sofa_df["sofa_cns"] = sofa_df["gcs_min"].apply(sofa_cns)
sofa_df["sofa_renal"] = sofa_df["creat_max_24h"].apply(sofa_renal)

sofa_df["sofa_total"] = (sofa_df["sofa_resp"] + sofa_df["sofa_coag"] +
                          sofa_df["sofa_liver"] + sofa_df["sofa_cv"] +
                          sofa_df["sofa_cns"] + sofa_df["sofa_renal"])

# Also calculate SOFA without renal component (non-renal SOFA)
sofa_df["sofa_nonrenal"] = (sofa_df["sofa_resp"] + sofa_df["sofa_coag"] +
                             sofa_df["sofa_liver"] + sofa_df["sofa_cv"] +
                             sofa_df["sofa_cns"])

log(f"  SOFA calculated for {len(sofa_df):,} patients")

# ============================================================
# Step 6: Merge with full cohort and save
# ============================================================
log("\nStep 6: Merging with cohort...")

sofa_cols = ["stay_id", "sofa_resp", "sofa_coag", "sofa_liver", "sofa_cv",
             "sofa_cns", "sofa_renal", "sofa_total", "sofa_nonrenal",
             "gcs_min", "map_min", "platelet_min", "bili_max",
             "pao2fio2_min", "norepi_max_rate", "epi_max_rate",
             "dopa_max_rate", "doba_max_rate", "creat_max_24h"]

sofa_final = sofa_df[sofa_cols].copy()

cohort_with_sofa = df.merge(sofa_final, on="stay_id", how="left")
for col in ["sofa_resp", "sofa_coag", "sofa_liver", "sofa_cv", "sofa_cns",
            "sofa_renal", "sofa_total", "sofa_nonrenal"]:
    cohort_with_sofa[col] = cohort_with_sofa[col].fillna(0)

cohort_with_sofa.to_csv(os.path.join(OUTPUT_DIR, "aki_paradox_cohort_with_sofa.csv"), index=False)
log(f"  Saved: {len(cohort_with_sofa):,} patients")

# ============================================================
# Summary statistics
# ============================================================
log("\n" + "=" * 70)
log("SOFA Score Summary (AKI patients)")
log("=" * 70)

aki_sofa = cohort_with_sofa[cohort_with_sofa["max_aki_stage"] >= 1]
log(f"\nAKI patients: {len(aki_sofa):,}")
log(f"SOFA Total: mean={aki_sofa['sofa_total'].mean():.1f}, "
    f"median={aki_sofa['sofa_total'].median():.1f}, "
    f"IQR=[{aki_sofa['sofa_total'].quantile(0.25):.0f}-"
    f"{aki_sofa['sofa_total'].quantile(0.75):.0f}]")

log("\nSOFA by AKI stage:")
for stage in [1, 2, 3]:
    d = aki_sofa[aki_sofa["max_aki_stage"] == stage]
    log(f"  Stage {stage} (n={len(d):,}): SOFA={d['sofa_total'].mean():.1f}+/-{d['sofa_total'].std():.1f}")

log("\nSOFA by group:")
for name, ckd in [("Pure AKI", 0), ("AKI+CKD", 1)]:
    d = aki_sofa[aki_sofa["ckd"] == ckd]
    log(f"  {name} (n={len(d):,}): SOFA={d['sofa_total'].mean():.1f}+/-{d['sofa_total'].std():.1f}")

log("\nSOFA by group x stage:")
for stage in [1, 2, 3]:
    for name, ckd in [("PureAKI", 0), ("AKI+CKD", 1)]:
        d = aki_sofa[(aki_sofa["max_aki_stage"] == stage) & (aki_sofa["ckd"] == ckd)]
        if len(d) > 0:
            log(f"  S{stage} {name} (n={len(d):,}): SOFA={d['sofa_total'].mean():.1f}+/-{d['sofa_total'].std():.1f}, "
                f"mort={d['hospital_expire_flag'].mean()*100:.1f}%")

log("\nNon-renal SOFA by group x stage (excludes renal component):")
for stage in [1, 2, 3]:
    for name, ckd in [("PureAKI", 0), ("AKI+CKD", 1)]:
        d = aki_sofa[(aki_sofa["max_aki_stage"] == stage) & (aki_sofa["ckd"] == ckd)]
        if len(d) > 0:
            log(f"  S{stage} {name}: nonrenal_SOFA={d['sofa_nonrenal'].mean():.1f}+/-{d['sofa_nonrenal'].std():.1f}")

log("\nData completeness (AKI patients):")
for col in ["gcs_min", "map_min", "platelet_min", "bili_max", "pao2fio2_min",
            "norepi_max_rate", "creat_max_24h"]:
    nn = aki_sofa[col].notna().sum()
    nz = (aki_sofa[col] > 0).sum() if col != "creat_max_24h" else nn
    log(f"  {col}: {nn:,}/{len(aki_sofa):,} ({nn/len(aki_sofa)*100:.1f}%)")

log("\nSOFA component distribution:")
for comp in ["sofa_resp", "sofa_coag", "sofa_liver", "sofa_cv", "sofa_cns", "sofa_renal"]:
    dist = aki_sofa[comp].value_counts().sort_index()
    log(f"  {comp}: {dict(dist)}")

log("\nDone! Output: aki_paradox_cohort_with_sofa.csv")
