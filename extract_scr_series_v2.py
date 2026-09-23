#!/usr/bin/env python3
"""
Extract & cache serum creatinine TIME SERIES for sensitivity analyses (v2).

Needed for:
  (2) eGFR-imputed baseline re-staging (CKD-EPI 2021, eGFR=75)
  (4) Transient vs persistent AKI phenotype (resolution within 72h of onset)

MIMIC: labevents itemid 50912, restricted to [admittime-365d, intime+7d]
       per patient (cohort = kdigo_cohort_v2.csv, patient-level)
eICU : lab table labname='creatinine', labresultoffset in
       [min(hospitaladmitoffset, -10080), +10080] min per stay
       (cohort = eicu_cohort_v2.csv, patient-level)

Outputs (cache):
  mimic_scr_series.csv.gz  (subject_id, charttime, valuenum)
  eicu_scr_series.csv.gz   (patientunitstayid, labresultoffset, labresult)
"""
import pandas as pd
import numpy as np
import duckdb
import os, gc
from datetime import datetime

WORK = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
MIMIC_HOSP = "E:/mimic-iv/v3.1/physionet.org/files/mimiciv/3.1/hosp"
EICU = "E:/mimic-iv/eicu/physionet.org/files/eicu-crd/2.0"
CREATININE_ITEMID = 50912
CR_MIN, CR_MAX = 0.1, 30.0

print("=" * 70)
print("SCr series extraction for sensitivity analyses")
print(f"Start: {datetime.now()}")
print("=" * 70)

# ============================================================
# MIMIC-IV
# ============================================================
out_mimic = os.path.join(WORK, "mimic_scr_series.csv.gz")
if os.path.exists(out_mimic):
    print(f"\n[MIMIC] cache exists, skip: {out_mimic}")
else:
    print("\n[MIMIC] scanning labevents (chunked) ...")
    coh = pd.read_csv(os.path.join(WORK, "kdigo_cohort_v2.csv"),
                      usecols=['subject_id', 'admittime', 'intime'])
    coh['admittime'] = pd.to_datetime(coh['admittime'])
    coh['intime'] = pd.to_datetime(coh['intime'])
    adm_map = coh.set_index('subject_id')['admittime']
    icu_map = coh.set_index('subject_id')['intime']
    subjects = set(coh['subject_id'])
    print(f"  cohort subjects: {len(subjects):,}")

    parts = []
    chunk_size = 5_000_000
    n_chunks = 0
    for chunk in pd.read_csv(os.path.join(MIMIC_HOSP, "labevents.csv.gz"),
                             compression='gzip', chunksize=chunk_size,
                             usecols=['subject_id', 'itemid', 'charttime', 'valuenum'],
                             dtype={'subject_id': 'int64', 'itemid': 'int64',
                                    'valuenum': 'float64'},
                             low_memory=False):
        n_chunks += 1
        m = chunk[(chunk['itemid'] == CREATININE_ITEMID)
                  & (chunk['subject_id'].isin(subjects))]
        if len(m) > 0:
            parts.append(m)
        if n_chunks % 5 == 0:
            print(f"    ... scanned {n_chunks * chunk_size / 1e6:.0f}M rows "
                  f"({datetime.now().strftime('%H:%M:%S')})")
    cr = pd.concat(parts, ignore_index=True)
    del parts, chunk; gc.collect()
    cr['charttime'] = pd.to_datetime(cr['charttime'])
    cr = cr[(cr['valuenum'] >= CR_MIN) & (cr['valuenum'] <= CR_MAX)]
    cr['admittime'] = cr['subject_id'].map(adm_map)
    cr['intime'] = cr['subject_id'].map(icu_map)
    cr = cr[(cr['charttime'] >= cr['admittime'] - pd.Timedelta(days=365))
            & (cr['charttime'] <= cr['intime'] + pd.Timedelta(days=7))]
    cr = cr.sort_values(['subject_id', 'charttime'])
    cr[['subject_id', 'charttime', 'valuenum']].to_csv(out_mimic,
                                                       index=False, compression='gzip')
    print(f"  saved: {out_mimic} ({len(cr):,} rows, "
          f"{cr['subject_id'].nunique():,} patients)")

# ============================================================
# eICU-CRD
# ============================================================
out_eicu = os.path.join(WORK, "eicu_scr_series.csv.gz")
if os.path.exists(out_eicu):
    print(f"\n[eICU] cache exists, skip: {out_eicu}")
else:
    print("\n[eICU] querying lab table (DuckDB) ...")
    CSV_OPTS = ", quote='\"', ignore_errors=true"
    coh = pd.read_csv(os.path.join(WORK, "eicu_cohort_v2.csv"),
                      usecols=['patientunitstayid'])
    # hospitaladmitoffset lives in patient.csv.gz (not in cohort CSV)
    hosp_off = duckdb.query(f"""
        SELECT patientunitstayid, hospitaladmitoffset
        FROM read_csv('{EICU}/patient.csv.gz'{CSV_OPTS})
    """).df()
    coh = coh.merge(hosp_off, on='patientunitstayid', how='left')
    coh['hospitaladmitoffset'] = pd.to_numeric(coh['hospitaladmitoffset'],
                                               errors='coerce').fillna(-10080)
    # per-stay lower bound: max(hospital admit offset, -7d) - small margin
    lo = np.maximum(coh['hospitaladmitoffset'], -10080)
    lo_map = dict(zip(coh['patientunitstayid'], lo))
    # build VALUES list for join (113k stays -> fine as temp DF)
    lo_df = pd.DataFrame({'patientunitstayid': coh['patientunitstayid'],
                          'lo_off': lo})

    cr = duckdb.query(f"""
        SELECT l.patientunitstayid, l.labresultoffset, l.labresult
        FROM read_csv('{EICU}/lab.csv.gz'{CSV_OPTS}) l
        WHERE l.labname = 'creatinine'
          AND l.labresult >= {CR_MIN}
          AND l.labresult <= {CR_MAX}
          AND l.labresultoffset IS NOT NULL
    """).df()
    print(f"  all creatinine rows: {len(cr):,}")
    cr = cr.merge(lo_df, on='patientunitstayid', how='inner')
    cr['labresultoffset'] = cr['labresultoffset'].astype(float)
    cr = cr[(cr['labresultoffset'] >= cr['lo_off'] - 1440)
            & (cr['labresultoffset'] <= 10080)]
    cr = cr.drop(columns='lo_off')
    cr = cr.sort_values(['patientunitstayid', 'labresultoffset'])
    cr.to_csv(out_eicu, index=False, compression='gzip')
    print(f"  saved: {out_eicu} ({len(cr):,} rows, "
          f"{cr['patientunitstayid'].nunique():,} stays)")

print(f"\nDone: {datetime.now()}")
