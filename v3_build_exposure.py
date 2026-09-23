# -*- coding: utf-8 -*-
"""
v3 Phase 1: time-varying KDIGO exposure + competing-event follow-up.

Fixes (vs v2):
  5.1 immortal time : stage is a monotone step function computed from data up to
                      time t only; exposure enters Cox as time-varying covariate
                      in counting-process (start, stop] format from ICU admission.
  5.2 30-d truncation: event = death AND death time <= 720 h; deaths after 720 h
                      are censored at 720 h (NOT counted as events).
  5.3 competing risk: discharge alive retained as event type 2 (CIF analysis);
                      for cause-specific Cox it censors (assumption stated).
  5.5 leakage       : phenotype (transient/persistent) classified ONLY from data
                      up to AKI onset + 72 h; deaths/discharges before the
                      landmark are excluded (landmark attrition, reported).
  flow reconciliation: every count asserted non-negative and summed.

Outputs (v3_outputs/):
  v3_mimic_patients.csv / v3_eicu_patients.csv   (patient level)
  v3_mimic_tv.csv.gz    / v3_eicu_tv.csv.gz      (counting process long)
  v3_flow.json                                   (auto reconciled)
"""
import pandas as pd
import numpy as np
import os, json, gc
from datetime import datetime

OUT = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
V3 = os.path.join(OUT, "v3_outputs")
os.makedirs(V3, exist_ok=True)

WINDOW_H = 168.0      # 7-day staging window
RISE_H = 48.0         # rolling 0.3 mg/dL window
ABS_CUT = 4.0         # Stage 3 absolute SCr threshold
TAU = 720.0           # 30-day horizon
SPLIT1, SPLIT2 = 168.0, 336.0   # period boundaries for time-stratified HRs

COVARS = ["age", "sex_male", "sofa_nonrenal", "vent_24h",
          "ckd_history_prior", "diabetes", "hypertension",
          "heart_failure", "copd", "liver_disease"]

print("=" * 70)
print(f"v3 Phase 1: time-varying exposure build  Start: {datetime.now()}")
print("=" * 70)


# ----------------------------------------------------------------------
def stage_measurements(t, v, base, rrt_h):
    """Per-timepoint achieved stage (mirrors v2 criteria) + inserted RRT point.
    t: sorted hours from ICU admission; v: SCr values. Returns (t2, s2)."""
    n = len(v)
    s = np.zeros(n, dtype=np.int8)
    if n:
        ratio = v / base if base > 0 else np.zeros(n)
        s_ratio = np.where(ratio >= 3.0, 3, np.where(ratio >= 2.0, 2,
                            np.where(ratio >= 1.5, 1, 0)))
        # rolling 48h rise: v[j] - min(v[i]) over t[i] >= t[j]-48
        s_rise = np.zeros(n, dtype=np.int8)
        for j in range(n):
            lo = np.searchsorted(t, t[j] - RISE_H, side="left")
            s_rise[j] = 1 if (v[j] - v[lo:j + 1].min()) >= 0.3 else 0
        # absolute stage 3 with acute change
        s_abs = np.where((v >= ABS_CUT) & ((s_ratio >= 1) | (s_rise >= 1)), 3, 0)
        s = np.maximum(np.maximum(s_ratio, s_rise), s_abs).astype(np.int8)
    if rrt_h is not None and 0.0 <= rrt_h <= WINDOW_H:
        # insert RRT timepoint with stage 3
        k = np.searchsorted(t, rrt_h)
        t = np.insert(t, k, rrt_h)
        s = np.insert(s, k, 3).astype(np.int8)
    return t, s


def build_patient_record(t, s, end_time):
    """Stage step function. Returns (change_times, change_stages,
    stage_at_24h, stage_max, onset_any_h)."""
    if len(s) == 0:
        return [], [], 0, 0, np.nan
    cm = np.maximum.accumulate(s)
    stage_at_24 = int(cm[np.searchsorted(t, 24.0, side="right") - 1]) \
        if np.searchsorted(t, 24.0, side="right") > 0 else 0
    stage_max = int(cm[-1])
    idx_any = np.argmax(s >= 1) if (s >= 1).any() else None
    onset = float(t[idx_any]) if idx_any is not None else np.nan
    # step changes (within follow-up only; later handled by caller)
    ch_t, ch_s = [], []
    cur = 0
    for ti, si in zip(t, cm):
        if si > cur:
            ch_t.append(float(ti)); ch_s.append(int(si)); cur = int(si)
    return ch_t, ch_s, stage_at_24, stage_max, onset


def classify_phenotype(t, s, onset, rrt_h):
    """Transient vs persistent at landmark L = onset + 72 h, using data <= L only.
    transient = last measurement in (onset, L] shows NO active KDIGO criterion
    (and no RRT by L)."""
    L = onset + 72.0
    if rrt_h is not None and rrt_h <= L:
        return "persistent", L
    m = (t > onset) & (t <= L)
    if not m.any():
        return "unknown", L
    last_stage = int(s[np.where(m)[0][-1]])
    return ("transient" if last_stage == 0 else "persistent"), L


# ----------------------------------------------------------------------
def build_cohort(tag, cohort_file, cov_file, scr_file, rrt_file,
                 id_col, scr_key, time_mode):
    """tag: 'mimic'/'eicu'; time_mode: 'abs' (MIMIC datetimes) or 'offset' (eICU)."""
    print(f"\n{'='*60}\n[{tag.upper()}] building\n{'='*60}")
    flow = {}

    cohort = pd.read_csv(cohort_file)
    cov = pd.read_csv(cov_file)
    scr = pd.read_csv(scr_file)
    rrt = pd.read_csv(rrt_file)

    # --- merge covariates (drop v2 leftover dup columns from cohort first) ---
    cov_ids = cov["stay_id"] if "stay_id" in cov.columns else cov[id_col]
    cov = cov.copy()
    cov[id_col] = cov_ids
    drop_cols = [c for c in COVARS if c in cohort.columns]
    if drop_cols:
        cohort = cohort.drop(columns=drop_cols)
    cohort = cohort.merge(cov[[id_col] + COVARS], on=id_col, how="left")

    # --- event times (hours from ICU admission) ---
    if time_mode == "abs":
        cohort["intime"] = pd.to_datetime(cohort["intime"])
        for c in ["outtime", "dischtime", "deathtime"]:
            cohort[c] = pd.to_datetime(cohort[c], errors="coerce")
        died = cohort["hospital_expire_flag"] == 1
        t_death = np.where(died,
                           (np.minimum(cohort["deathtime"], cohort["dischtime"])
                            - cohort["intime"]).dt.total_seconds() / 3600.0,
                           np.inf)
        t_disc = (cohort["dischtime"] - cohort["intime"]).dt.total_seconds() / 3600.0
    else:
        died = cohort["hospitaldischargestatus"] == "Expired"
        t_hosp = pd.to_numeric(cohort["hospitaldischargeoffset"],
                               errors="coerce") / 60.0
        t_death = np.where(died, t_hosp, np.inf)
        t_disc = t_hosp

    t_death = np.clip(np.nan_to_num(t_death, nan=np.inf), 0, None)
    t_disc = np.clip(np.nan_to_num(t_disc, nan=np.inf), 0, None)

    event_type = np.zeros(len(cohort), dtype=np.int8)   # 0 censored, 1 death, 2 discharge
    end_time = np.full(len(cohort), TAU)
    death_le = died.values & (t_death <= TAU)
    event_type[death_le] = 1
    end_time[death_le] = t_death[death_le]
    disc_le = (~died.values | (t_death > TAU)) & (t_disc <= TAU)
    event_type[disc_le] = 2
    end_time[disc_le] = t_disc[disc_le]
    end_time = np.clip(end_time, 0.1, None)

    cohort["end_time"] = end_time
    cohort["event_type"] = event_type

    flow["n_cohort"] = int(len(cohort))
    flow["n_death_le_30d"] = int((event_type == 1).sum())
    flow["n_discharge_le_30d"] = int((event_type == 2).sum())
    flow["n_censored_30d"] = int((event_type == 0).sum())
    flow["n_death_after_30d_not_event"] = int((died.values & (t_death > TAU)).sum())

    # --- SCr series -> hours from admission, window [0, 168] ---
    if time_mode == "abs":
        scr = scr.merge(cohort[["subject_id", "intime"]], on="subject_id",
                        how="inner")
        scr["t_h"] = (pd.to_datetime(scr["charttime"]) - scr["intime"]
                      ).dt.total_seconds() / 3600.0
    else:
        scr = scr.merge(cohort[["patientunitstayid"]], on="patientunitstayid",
                        how="inner")
        scr["t_h"] = pd.to_numeric(scr["labresultoffset"], errors="coerce") / 60.0
        scr["valuenum"] = scr["labresult"]
    scr = scr[(scr["t_h"] >= 0) & (scr["t_h"] <= WINDOW_H)]
    scr = scr.dropna(subset=["valuenum"])
    scr = scr.sort_values([scr_key, "t_h"])
    print(f"  SCr measurements in window: {len(scr):,} "
          f"({scr[scr_key].nunique():,} patients)")

    # --- RRT times ---
    rrt_map = rrt.set_index(id_col if id_col in rrt.columns else rrt.columns[0])
    rrt_map = rrt_map["rrt_time_h"].to_dict()

    # --- per-patient loop ---
    grouped = dict(tuple(scr.groupby(scr_key)))
    base_map = cohort.set_index(id_col)["baseline_cr"].to_dict()
    ids = cohort[id_col].values
    keys = cohort[scr_key].values  # subject_id (MIMIC) / stay id (eICU)
    end_map = dict(zip(ids, end_time))
    ev_map = dict(zip(ids, event_type))

    pat_rows, tv_rows = [], []
    n_stage_up = 0
    for i, (pid, key) in enumerate(zip(ids, keys)):
        gr = grouped.get(key)
        if gr is None:
            t, v = np.array([]), np.array([])
        else:
            t = gr["t_h"].values.astype(float)
            v = gr["valuenum"].values.astype(float)
        rrt_h = rrt_map.get(pid, None)
        t2, s2 = stage_measurements(t, v, base_map.get(pid, np.nan), rrt_h)
        end = end_map[pid]
        ch_t, ch_s, s24, smax, onset = build_patient_record(t2, s2, end)

        # phenotype (among patients with AKI onset, before attrition screening)
        phen, L = (None, None)
        if not np.isnan(onset):
            phen, L = classify_phenotype(t2, s2, onset, rrt_h)

        n_stage_up += len(ch_t)

        # counting-process intervals: breakpoints = change times + period splits
        bps = sorted({0.0} | {x for x in ch_t if 0 < x < end}
                     | {x for x in (SPLIT1, SPLIT2) if 0 < x < end} | {end})
        cur_stage = 0
        ci = 0
        for a, b in zip(bps[:-1], bps[1:]):
            while ci < len(ch_t) and ch_t[ci] <= a:
                cur_stage = ch_s[ci]; ci += 1
            is_last = (b == end)
            tv_rows.append((pid, a, b, cur_stage,
                            int(ev_map[pid]) if is_last else 0,
                            int(ev_map[pid] == 1) if is_last else 0,
                            int(ev_map[pid] == 2) if is_last else 0))
        if len(bps) == 1:  # end == 0 edge (should not happen, end>=0.1)
            tv_rows.append((pid, 0.0, end, 0, int(ev_map[pid]),
                            int(ev_map[pid] == 1), int(ev_map[pid] == 2)))

        pat_rows.append({
            id_col: pid, "stage_at_24h": s24, "stage_max_7d": smax,
            "onset_any_h": onset, "phenotype": phen, "landmark_h": L,
            "end_time": end, "event_type": int(ev_map[pid]),
            "rrt_time_h": rrt_h})
        if (i + 1) % 20000 == 0:
            print(f"    ... {i+1:,} patients processed")

    pat = pd.DataFrame(pat_rows)
    tv = pd.DataFrame(tv_rows, columns=[id_col, "t_start", "t_stop", "stage",
                                        "event_type", "event_death", "event_disc"])
    tv["stage1"] = (tv["stage"] == 1).astype(int)
    tv["stage2"] = (tv["stage"] == 2).astype(int)
    tv["stage3"] = (tv["stage"] == 3).astype(int)
    tv["period"] = np.where(tv["t_start"] < SPLIT1, 1,
                    np.where(tv["t_start"] < SPLIT2, 2, 3))

    # --- merge covariates + event info into patient file ---
    pat = pat.merge(cohort[[id_col] + COVARS + ["end_time", "event_type"]],
                    on=id_col, how="left", suffixes=("", "_dup"))
    pat = pat.drop(columns=[c for c in pat.columns if c.endswith("_dup")])

    # --- phenotype attrition flow ---
    aki = pat[pat["onset_any_h"].notna()]
    L_ok = aki[(aki["end_time"] > aki["landmark_h"]) | (aki["landmark_h"].isna())]
    # eligible: still at risk at landmark (no death/discharge before L)
    flow["n_aki_any"] = int(len(aki))
    flow["n_phen_died_before_landmark"] = int(
        ((aki["event_type"] == 1) & (aki["end_time"] <= aki["landmark_h"])).sum())
    flow["n_phen_discharged_before_landmark"] = int(
        ((aki["event_type"] == 2) & (aki["end_time"] <= aki["landmark_h"])).sum())
    flow["n_phen_unknown"] = int((aki["phenotype"] == "unknown").sum())
    phen_elig = aki[(aki["phenotype"].isin(["transient", "persistent"]))
                    & (aki["end_time"] > aki["landmark_h"])]
    flow["n_phen_eligible"] = int(len(phen_elig))
    flow["n_transient"] = int((phen_elig["phenotype"] == "transient").sum())
    flow["n_persistent"] = int((phen_elig["phenotype"] == "persistent").sum())
    flow["n_tv_rows"] = int(len(tv))
    flow["n_stage_changes"] = int(n_stage_up)

    # --- reconciliation assertions ---
    assert flow["n_death_le_30d"] + flow["n_discharge_le_30d"] \
        + flow["n_censored_30d"] == flow["n_cohort"]
    assert flow["n_transient"] + flow["n_persistent"] == flow["n_phen_eligible"]
    assert all(v >= 0 for k, v in flow.items() if isinstance(v, int))
    for k, v in flow.items():
        print(f"  {k}: {v:,}")

    pat.to_csv(os.path.join(V3, f"v3_{tag}_patients.csv"), index=False)
    tv.to_csv(os.path.join(V3, f"v3_{tag}_tv.csv.gz"), index=False, compression="gzip")
    del grouped, scr
    gc.collect()
    return flow


# ----------------------------------------------------------------------
flow_all = {}
flow_all["mimic"] = build_cohort(
    "mimic",
    os.path.join(OUT, "kdigo_cohort_v2.csv"),
    os.path.join(OUT, "mimic_analysis_v2.csv"),
    os.path.join(OUT, "mimic_scr_series.csv.gz"),
    os.path.join(OUT, "rrt_times_mimic.csv"),
    "stay_id", "subject_id", "abs")

flow_all["eicu"] = build_cohort(
    "eicu",
    os.path.join(OUT, "eicu_cohort_v2.csv"),
    os.path.join(OUT, "eicu_analysis_v2.csv"),
    os.path.join(OUT, "eicu_scr_series.csv.gz"),
    os.path.join(OUT, "rrt_times_eicu.csv"),
    "patientunitstayid", "patientunitstayid", "offset")

flow_all["generated"] = datetime.now().isoformat()
with open(os.path.join(V3, "v3_flow.json"), "w") as f:
    json.dump(flow_all, f, indent=2)

print("\n" + "=" * 70)
print("v3 Phase 1 COMPLETE — all flow assertions passed")
print("=" * 70)
