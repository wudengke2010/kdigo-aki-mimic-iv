#!/usr/bin/env python3
"""
MIMIC-IV first-24h SOFA (mimic-code style) — F3 修复
=====================================================
按 mimic-code 验证 SQL 的语义重建 first-day SOFA，从 MIMIC-IV v3.1 原始表直接计算。
对照评审报告 F3 (kdigo_revised_data.py:637-682) 的"自制近似SOFA"问题：

  - v1 评分：任意通气=3, 任意升压药=3, 肾脏项直接用 peak_cr（与暴露共享信息）,
            缺失默认 0
  - v2 评分：呼吸 = PaO2/FiO2 比值 + 通气交互（<100=4, <200=3, <300=2, <400=1）
            凝血 = min(血小板)
            肝脏 = max(总胆红素)
            心血管 = min(MAP) + 升压药分级（norepi/epi>0=4, dopa>0=3, dobu>0=2）
            神经 = min(GCS总)
            肾脏 = max(肌酐) + 24h 尿量（仅 full SOFA 计入）

主分析使用 **non-renal SOFA** (5 个分量，最大 20 分)，
灵敏度分析使用 full SOFA (6 个分量，最大 24 分)。

时间窗：[intime - 6h, intime + 24h]，与 mimic-code 一致。

每个分量对缺失值记为 NULL，最终仅当全部分量可用时求和；若任一分量缺失
则在最终记录 sofa_n_missing 计数（最差为 6，最好为 0），下游统计层据此判断
是否以 median/0 填充。
"""
import duckdb
import pandas as pd
import numpy as np
import os, json, gc
from datetime import datetime

print("=" * 80)
print("MIMIC-IV first-24h SOFA (mimic-code style) — F3 修复")
print(f"Start: {datetime.now()}")
print("=" * 80)

HOSP = "E:/mimic-iv/v3.1/physionet.org/files/mimiciv/3.1/hosp"
ICU  = "E:/mimic-iv/v3.1/physionet.org/files/mimiciv/3.1/icu"
OUT = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
CSV_OPTS = ", quote='\"', ignore_errors=true"

def q(sql):
    return duckdb.query(sql).df()

# ============================================================
# Itemid mapping (verified against d_items / d_labitems 2026-09-23)
# ============================================================
ITEMID = {
    # Vasopressors (inputevents)
    # norepi/epi in mcg/kg/min (mimic-code thresholds: >0.1 = high = score 4)
    'norepinephrine': 221906,
    'epinephrine':   221289,
    'epinephrine_alt': 229617,   # 'Epinephrine.' alternate itemid
    'dopamine':      221662,
    'dobutamine':    221653,
    # MAP (chartevents) — 220052 is the standard metavision "Arterial BP Mean"
    # (3M+ rows / 36k+ stays), 225312 is ART-line specific, 226766/227023 are
    # APACHE-recomputed (sparse, used as fallback only)
    'map_primary':    220052,
    'map_art_line':   225312,
    'map_apache':     226766,    # MapApacheIIValue (validated, sparse)
    'map_apacheiv':   227023,    # MAP_ApacheIV
    # GCS (chartevents)
    'gcs_eye':       220739,
    'gcs_motor':     223901,
    'gcs_verbal':    223900,
    'gcs_total_apache':   226755,   # GcsApacheIIScore
    'gcs_total_apacheiv': 227013,   # GcsScore_ApacheIV
    # FiO2 (chartevents)
    'fio2_apache':   226754,
    'fio2_set':      223835,    # Inspired O2 Fraction
    # PaO2 (labevents)
    'po2_blood':     50821,
    # Bilirubin, Platelets, Creatinine (labevents)
    'bilirubin_total': 50885,
    'platelet':      51265,
    'creatinine':    50912,     # same as cohort script
    # Urine output (outputevents)
    'urine_foley':   226559,
    'urine_void':    226560,
    'urine_irrigant':227489,
    'urine_gu_irrigant_out': 226566,
    # Ventilation (procedureevents)
    'invasive_vent': 225792,
}

# ============================================================
# Load cohort (MIMIC v2)
# ============================================================
print("\n[Load cohort] kdigo_cohort_v2.csv ...")
cohort = pd.read_csv(os.path.join(OUT, "kdigo_cohort_v2.csv"))
print(f"  Cohort: {len(cohort):,} stays")

stay_ids = set(cohort['stay_id'])
intime_lookup = cohort.set_index('stay_id')['intime'].to_dict()

# ============================================================
# Build per-stay time window in DuckDB
# ============================================================
print("\n[Time-window joiner] ...")
duckdb.query("""
    CREATE OR REPLACE TEMP TABLE cohort_window AS
    SELECT stay_id, subject_id, hadm_id,
           intime,
           intime - INTERVAL '6' HOUR AS win_start,
           intime + INTERVAL '24' HOUR AS win_end
    FROM read_csv_auto('""" + OUT + """/kdigo_cohort_v2.csv')
""")

# ============================================================
# COMPONENT 1: Cardiovascular (MAP + vasopressors)
# ============================================================
print("\n[Cardiovascular] MAP_min + norepi/epi/dopa/dobu ...")

# MAP: prefer MapApacheIIValue (validated)
map_df = q(f"""
    SELECT cw.stay_id, MIN(ce.valuenum) AS map_min
    FROM cohort_window cw
    JOIN read_csv('{ICU}/chartevents.csv.gz'{CSV_OPTS}) ce
      ON ce.stay_id = cw.stay_id
     AND ce.itemid IN ({ITEMID['map_primary']}, {ITEMID['map_art_line']},
                       {ITEMID['map_apache']}, {ITEMID['map_apacheiv']})
     AND ce.charttime BETWEEN cw.win_start AND cw.win_end
     AND ce.valuenum BETWEEN 20 AND 250
    GROUP BY cw.stay_id
""")
print(f"  MAP records: {len(map_df):,} stays")

# Vasopressors (inputevents) — per-stay max rate in mcg/kg/min
vaso = q(f"""
    SELECT stay_id, itemid,
           MAX(CASE WHEN LOWER(rateuom) = 'mcg/kg/min' THEN rate
                                WHEN LOWER(rateuom) = 'mg/kg/min' THEN rate * 1000
                                ELSE NULL END) AS max_rate_mcgkgmin
    FROM read_csv('{ICU}/inputevents.csv.gz'{CSV_OPTS})
    WHERE itemid IN ({ITEMID['norepinephrine']}, {ITEMID['epinephrine']},
                     {ITEMID['epinephrine_alt']},
                     {ITEMID['dopamine']}, {ITEMID['dobutamine']})
    GROUP BY stay_id, itemid
""")
# Per-stay aggregates
vaso['norepi_rate'] = vaso.loc[vaso['itemid'].isin(
    [ITEMID['norepinephrine']]), 'max_rate_mcgkgmin']
vaso['epi_rate']    = vaso.loc[vaso['itemid'].isin(
    [ITEMID['epinephrine'], ITEMID['epinephrine_alt']]), 'max_rate_mcgkgmin']
vaso['dopa_rate']   = vaso.loc[vaso['itemid'].isin(
    [ITEMID['dopamine']]), 'max_rate_mcgkgmin']
vaso['dobu_rate']   = vaso.loc[vaso['itemid'].isin(
    [ITEMID['dobutamine']]), 'max_rate_mcgkgmin']
vaso_agg = vaso.groupby('stay_id')[['norepi_rate','epi_rate','dopa_rate','dobu_rate']].max().reset_index()
print(f"  Vaso stays: {len(vaso_agg):,}")

cardio = map_df.merge(vaso_agg, on='stay_id', how='outer')

# Score with dose tiers (mimic-code style)
def score_cardio(row):
    norepi = row.get('norepi_rate')
    epi    = row.get('epi_rate')
    dopa   = row.get('dopa_rate')
    dobu   = row.get('dobu_rate')
    # norepi/epi > 0.1 mcg/kg/min = high dose -> score 4
    if pd.notna(norepi) and norepi > 0.1: return 4
    if pd.notna(epi)    and epi > 0.1:    return 4
    # norepi/epi 0 < r <= 0.1 = low dose -> score 3
    if pd.notna(norepi) and norepi > 0:   return 3
    if pd.notna(epi)    and epi > 0:      return 3
    # dopamine tiers: >15 = 4, 5-15 = 3, 0-5 = 2
    if pd.notna(dopa) and dopa > 15: return 4
    if pd.notna(dopa) and dopa > 5:  return 3
    if pd.notna(dopa) and dopa > 0:  return 2
    # dobutamine any dose -> score 2
    if pd.notna(dobu) and dobu > 0:  return 2
    # MAP < 70 -> score 1
    if pd.notna(row.get('map_min')) and row['map_min'] < 70: return 1
    if pd.notna(row.get('map_min')): return 0
    return np.nan  # all missing

cardio['cardiovascular'] = cardio.apply(score_cardio, axis=1)
cardio = cardio[['stay_id','cardiovascular']]
print(f"  CV score distribution: {cardio['cardiovascular'].value_counts(dropna=False).sort_index().to_dict()}")

# ============================================================
# COMPONENT 2: Coagulation (platelets_min)
# ============================================================
print("\n[Coagulation] platelets_min ...")
plt = q(f"""
    SELECT cw.stay_id, MIN(le.valuenum) AS platelet_min
    FROM cohort_window cw
    JOIN read_csv('{HOSP}/labevents.csv.gz'{CSV_OPTS}) le
      ON le.hadm_id = cw.hadm_id
     AND le.itemid = {ITEMID['platelet']}
     AND le.charttime BETWEEN cw.win_start AND cw.win_end
     AND le.valuenum BETWEEN 0 AND 1000
    GROUP BY cw.stay_id
""")

def score_plt(p):
    if pd.isna(p): return np.nan
    if p < 20: return 4
    if p < 50: return 3
    if p < 100: return 2
    if p < 150: return 1
    return 0

plt['coagulation'] = plt['platelet_min'].apply(score_plt)
plt = plt[['stay_id','coagulation']]
print(f"  Platelet stays: {len(plt):,}, score: {plt['coagulation'].value_counts(dropna=False).sort_index().to_dict()}")

# ============================================================
# COMPONENT 3: Liver (bilirubin_total_max)
# ============================================================
print("\n[Liver] bilirubin_total_max ...")
bili = q(f"""
    SELECT cw.stay_id, MAX(le.valuenum) AS bili_max
    FROM cohort_window cw
    JOIN read_csv('{HOSP}/labevents.csv.gz'{CSV_OPTS}) le
      ON le.hadm_id = cw.hadm_id
     AND le.itemid = {ITEMID['bilirubin_total']}
     AND le.charttime BETWEEN cw.win_start AND cw.win_end
     AND le.valuenum BETWEEN 0 AND 50
    GROUP BY cw.stay_id
""")

def score_bili(b):
    if pd.isna(b): return np.nan
    if b >= 12.0: return 4
    if b >= 6.0: return 3
    if b >= 2.0: return 2
    if b >= 1.2: return 1
    return 0

bili['liver'] = bili['bili_max'].apply(score_bili)
bili = bili[['stay_id','liver']]
print(f"  Bili stays: {len(bili):,}, score: {bili['liver'].value_counts(dropna=False).sort_index().to_dict()}")

# ============================================================
# COMPONENT 4: CNS (GCS_min via 3 components summed)
# ============================================================
print("\n[CNS] GCS_min (Eye+Motor+Verbal summed per row) ...")
gcs = q(f"""
    WITH gcs_raw AS (
        SELECT cw.stay_id, ce.charttime,
               MAX(CASE WHEN ce.itemid = {ITEMID['gcs_eye']}    THEN ce.valuenum END) AS eye,
               MAX(CASE WHEN ce.itemid = {ITEMID['gcs_motor']}  THEN ce.valuenum END) AS motor,
               MAX(CASE WHEN ce.itemid = {ITEMID['gcs_verbal']} THEN ce.valuenum END) AS verbal
        FROM cohort_window cw
        JOIN read_csv('{ICU}/chartevents.csv.gz'{CSV_OPTS}) ce
          ON ce.stay_id = cw.stay_id
         AND ce.itemid IN ({ITEMID['gcs_eye']}, {ITEMID['gcs_motor']}, {ITEMID['gcs_verbal']})
         AND ce.charttime BETWEEN cw.win_start AND cw.win_end
        GROUP BY cw.stay_id, ce.charttime
    )
    SELECT stay_id,
           MIN(eye + motor + verbal) AS gcs_min
    FROM gcs_raw
    WHERE eye BETWEEN 1 AND 4
      AND motor BETWEEN 1 AND 6
      AND verbal BETWEEN 1 AND 5
    GROUP BY stay_id
""")

def score_gcs(g):
    if pd.isna(g): return np.nan
    if g < 6: return 4
    if g <= 9: return 3
    if g <= 12: return 2
    if g <= 14: return 1
    return 0

gcs['cns'] = gcs['gcs_min'].apply(score_gcs)
gcs = gcs[['stay_id','cns']]
print(f"  GCS stays: {len(gcs):,}, score: {gcs['cns'].value_counts(dropna=False).sort_index().to_dict()}")

# ============================================================
# COMPONENT 5: Respiration (PaO2/FiO2 + vent interaction)
# ============================================================
print("\n[Respiration] PaO2/FiO2 + vent ...")

# Define vent = any invasive vent procedure in window
vented = q(f"""
    SELECT DISTINCT cw.stay_id
    FROM cohort_window cw
    JOIN read_csv('{ICU}/procedureevents.csv.gz'{CSV_OPTS}) pe
      ON pe.stay_id = cw.stay_id
     AND pe.itemid = {ITEMID['invasive_vent']}
     AND pe.starttime BETWEEN cw.win_start AND cw.win_end
""")
vented_set = set(vented['stay_id'])
print(f"  Ventilated (invasive vent proc in 24h): {len(vented_set):,}")

# Min PaO2 (any pO2 in blood) over window
pao2 = q(f"""
    SELECT cw.stay_id, MIN(le.valuenum) AS pao2_min
    FROM cohort_window cw
    JOIN read_csv('{HOSP}/labevents.csv.gz'{CSV_OPTS}) le
      ON le.hadm_id = cw.hadm_id
     AND le.itemid = {ITEMID['po2_blood']}
     AND le.charttime BETWEEN cw.win_start AND cw.win_end
     AND le.valuenum BETWEEN 20 AND 600
    GROUP BY cw.stay_id
""")

# Max FiO2 (use FiO2ApacheIIValue as primary; if missing use Inspired O2 Fraction)
fio2 = q(f"""
    SELECT cw.stay_id,
           MAX(CASE WHEN ce.itemid = {ITEMID['fio2_apache']}
                    THEN ce.valuenum END) AS fio2_apache_max,
           MAX(CASE WHEN ce.itemid = {ITEMID['fio2_set']}
                    THEN ce.valuenum END) AS fio2_set_max
    FROM cohort_window cw
    JOIN read_csv('{ICU}/chartevents.csv.gz'{CSV_OPTS}) ce
      ON ce.stay_id = cw.stay_id
     AND ce.itemid IN ({ITEMID['fio2_apache']}, {ITEMID['fio2_set']})
     AND ce.charttime BETWEEN cw.win_start AND cw.win_end
     AND ce.valuenum BETWEEN 0.21 AND 1.0
    GROUP BY cw.stay_id
""")

resp = pao2.merge(fio2, on='stay_id', how='outer')
resp['fio2_max'] = resp[['fio2_apache_max','fio2_set_max']].max(axis=1)
resp['fio2_max'] = resp['fio2_max'].clip(lower=0.21, upper=1.0)
resp['isvent'] = resp['stay_id'].isin(vented_set).astype(int)

def score_resp(row):
    pao2 = row.get('pao2_min')
    fio2 = row.get('fio2_max')
    if pd.isna(pao2): return np.nan
    ratio = pao2 / fio2 if pd.notna(fio2) and fio2 > 0 else pao2
    if row['isvent'] == 1:
        if ratio < 100: return 4
        if ratio < 200: return 3
        return 2  # ventilated with ratio >= 200 -> 2 (mimic-code default for vent)
    else:
        if ratio < 300: return 2
        if ratio < 400: return 1
        return 0

resp['respiration'] = resp.apply(score_resp, axis=1)
resp = resp[['stay_id','respiration']]
print(f"  Resp stays: {len(resp):,}, score: {resp['respiration'].value_counts(dropna=False).sort_index().to_dict()}")

# ============================================================
# COMPONENT 6: Renal (creatinine_max + 24h urine output) — full SOFA only
# ============================================================
print("\n[Renal] creatinine_max + 24h urine output ...")
cr_max = q(f"""
    SELECT cw.stay_id, MAX(le.valuenum) AS cr_max
    FROM cohort_window cw
    JOIN read_csv('{HOSP}/labevents.csv.gz'{CSV_OPTS}) le
      ON le.hadm_id = cw.hadm_id
     AND le.itemid = {ITEMID['creatinine']}
     AND le.charttime BETWEEN cw.win_start AND cw.win_end
     AND le.valuenum BETWEEN 0.1 AND 30
    GROUP BY cw.stay_id
""")

urine = q(f"""
    SELECT cw.stay_id, SUM(oe.value) AS urine_24h
    FROM cohort_window cw
    JOIN read_csv('{ICU}/outputevents.csv.gz'{CSV_OPTS}) oe
      ON oe.stay_id = cw.stay_id
     AND oe.itemid IN ({ITEMID['urine_foley']}, {ITEMID['urine_void']},
                       {ITEMID['urine_irrigant']}, {ITEMID['urine_gu_irrigant_out']})
     AND oe.charttime BETWEEN cw.win_start AND cw.win_end
    GROUP BY cw.stay_id
""")

renal = cr_max.merge(urine, on='stay_id', how='outer')

def score_renal(row):
    cr = row.get('cr_max')
    uo = row.get('urine_24h')
    # Urine criteria take priority per mimic-code (worst-of)
    if pd.notna(uo):
        if uo < 200: return 4
        if uo < 500: return 3
    if pd.notna(cr):
        if cr >= 5.0: return 4
        if cr >= 3.5: return 3
        if cr >= 2.0: return 2
        if cr >= 1.2: return 1
        return 0
    return np.nan

renal['renal'] = renal.apply(score_renal, axis=1)
renal = renal[['stay_id','renal']]
print(f"  Renal stays: {len(renal):,}, score: {renal['renal'].value_counts(dropna=False).sort_index().to_dict()}")

# ============================================================
# Merge all 6 components + compute sums
# ============================================================
print("\n[Merge] 6 components + non-renal + full SOFA ...")
sofa = cohort[['stay_id']].copy()
for d in [cardio, plt, bili, gcs, resp, renal]:
    sofa = sofa.merge(d, on='stay_id', how='left')

component_cols = ['cardiovascular','coagulation','liver','cns','respiration','renal']
sofa['sofa_n_missing'] = sofa[component_cols].isna().sum(axis=1)
sofa['sofa_nonrenal'] = sofa[['cardiovascular','coagulation','liver','cns','respiration']].sum(axis=1, min_count=1)
sofa['sofa_full']     = sofa[component_cols].sum(axis=1, min_count=1)

# Available-cohort coverage
print(f"  Non-renal SOFA computed: {sofa['sofa_nonrenal'].notna().sum():,} "
      f"(mean {sofa['sofa_nonrenal'].mean():.2f}, std {sofa['sofa_nonrenal'].std():.2f})")
print(f"  Full SOFA computed:     {sofa['sofa_full'].notna().sum():,}")
print(f"  Missing any component:  {(sofa['sofa_n_missing']>0).sum():,}")
print(f"  All 6 missing:          {(sofa['sofa_n_missing']==6).sum():,}")

# ============================================================
# Update cohort CSV in-place
# ============================================================
print("\n[Save] merging SOFA into kdigo_cohort_v2.csv ...")
out_cols = component_cols + ['sofa_nonrenal','sofa_full','sofa_n_missing']
final = cohort.merge(sofa[['stay_id'] + out_cols], on='stay_id', how='left')
final.to_csv(os.path.join(OUT, "kdigo_cohort_v2.csv"), index=False)

flow = {
    'cardio_records':   int(len(cardio)),
    'platelet_records': int(len(plt)),
    'bili_records':     int(len(bili)),
    'gcs_records':      int(len(gcs)),
    'resp_records':     int(len(resp)),
    'renal_records':    int(len(renal)),
    'sofa_nonrenal_available': int(sofa['sofa_nonrenal'].notna().sum()),
    'sofa_full_available':     int(sofa['sofa_full'].notna().sum()),
    'sofa_missing_any':        int((sofa['sofa_n_missing']>0).sum()),
    'sofa_nonrenal_mean':      float(sofa['sofa_nonrenal'].mean()),
    'sofa_nonrenal_std':       float(sofa['sofa_nonrenal'].std()),
    'sofa_full_mean':          float(sofa['sofa_full'].mean()),
    'sofa_full_std':           float(sofa['sofa_full'].std()),
}
with open(os.path.join(OUT, "sofa_v2_summary.json"), 'w') as f:
    json.dump(flow, f, indent=2)

print(f"\nDone: {datetime.now()}")