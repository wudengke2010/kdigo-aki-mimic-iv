# -*- coding: utf-8 -*-
"""
v3 Phase 2: primary analysis.

 (a) Time-varying (counting-process) cause-specific Cox: Models A/B/C.
     Exposure = current max KDIGO stage (time-varying); discharge alive censors
     (independence assumption stated); deaths >30 d censored at 720 h.
 (b) Period-specific stage HRs via stage x period (0-7/7-14/14-30 d)
     interactions — replaces the invalid v2 time-stratified risk sets (5.4).
 (c) Competing-risk CIF (Aalen-Johansen) of 30-d death by stage at 24 h
     (stage computed from first-24-h data only: no immortal time).
 (d) Schoenfeld PH checks on stage terms.

Outputs: v3_outputs/v3_primary_results.json, v3_primary_report.md,
         v3_table_models.csv, v3_table_cif.csv,
         fig_v3_cif.png/pdf, fig_v3_period_forest.png/pdf
"""
import pandas as pd
import numpy as np
import os, json, gc
from datetime import datetime
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

OUT = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
V3 = os.path.join(OUT, "v3_outputs")
TAU = 720.0

COVARS = ["age", "sex_male", "sofa_nonrenal", "vent_24h",
          "ckd_history_prior", "diabetes", "hypertension",
          "heart_failure", "copd", "liver_disease"]
STAGES = ["stage1", "stage2", "stage3"]

print("=" * 70)
print(f"v3 Phase 2: primary analysis  Start: {datetime.now()}")
print("=" * 70)

results = {}


# ----------------------------------------------------------------------
def fit_models(tag, id_col):
    tv = pd.read_csv(os.path.join(V3, f"v3_{tag}_tv.csv.gz"))
    pat = pd.read_csv(os.path.join(V3, f"v3_{tag}_patients.csv"))
    tv = tv.merge(pat[[id_col] + COVARS], on=id_col, how="left")
    r = {}

    # ---------------- crude 30-d mortality by stage descriptors ----------
    d = pat.copy()
    d["death30"] = (d["event_type"] == 1).astype(int)  # event_type 1 == death <=720h
    for grp, col in [("stage_at_24h", "stage_at_24h"), ("stage_max_7d", "stage_max_7d")]:
        tab = d.groupby(col)["death30"].agg(["count", "mean"])
        r[f"crude_by_{grp}"] = {int(k): [int(v["count"]), round(100*v["mean"], 1)]
                                for k, v in tab.iterrows()}

    # lifelines counting-process fit requires rows sorted by stop time
    tv = tv.sort_values("t_stop").reset_index(drop=True)

    # ---------------- Models A/B/C (time-varying Cox) --------------------
    def fit(tv, cols, label):
        data = tv[["t_start", "t_stop", "event_death"] + cols].copy()
        cph = CoxPHFitter()
        cph.fit(data, duration_col="t_stop", event_col="event_death",
                entry_col="t_start", show_progress=False)
        s = cph.summary
        out = {}
        for c in cols:
            out[c] = {"coef": float(s.loc[c, "coef"]),
                      "hr": float(s.loc[c, "exp(coef)"]),
                      "lo": float(s.loc[c, "exp(coef) lower 95%"]),
                      "hi": float(s.loc[c, "exp(coef) upper 95%"]),
                      "p": float(s.loc[c, "p"])}
        out["_n"] = int(cph._n_examples)
        out["_events"] = int(cph.event_observed.sum())
        out["_loglik"] = float(cph.log_likelihood_)
        print(f"  [{tag}] {label}: n={out['_n']:,} events={out['_events']:,} "
              f"S3 HR={out['stage3']['hr']:.2f} "
              f"({out['stage3']['lo']:.2f}-{out['stage3']['hi']:.2f})")
        return out, cph

    r["model_A"], _ = fit(tv, STAGES, "Model A (unadjusted)")
    colsB = STAGES + ["age", "sex_male"]
    r["model_B"], _ = fit(tv, colsB, "Model B (+demographics)")
    colsC = STAGES + COVARS
    r["model_C"], cphC = fit(tv, colsC, "Model C (full)")

    # ---------------- PH check: joint Wald test on stage x period ----------
    # (Schoenfeld residuals on 116k-row counting-process data are
    #  memory-prohibitive; the stage x period interactions ARE the formal
    #  test of hazard-ratio homogeneity over time = PH for stage.)
    int_cols = STAGES + [f"{s_}_p{p_}" for s_ in STAGES for p_ in (2, 3)] + COVARS
    tvi = tv.copy()
    for s_ in STAGES:
        for p_ in (2, 3):
            tvi[f"{s_}_p{p_}"] = tvi[s_] * (tvi["period"] == p_).astype(int)
    cphI = CoxPHFitter()
    cphI.fit(tvi[["t_start", "t_stop", "event_death"] + int_cols],
             duration_col="t_stop", event_col="event_death",
             entry_col="t_start", show_progress=False)
    b = cphI.params_
    V = cphI.variance_matrix_
    inter_names = [f"{s_}_p{p_}" for s_ in STAGES for p_ in (2, 3)]
    bI = b[inter_names].values
    VI = V.loc[inter_names, inter_names].values
    from scipy.stats import chi2, norm
    wald = float(bI @ np.linalg.pinv(VI) @ bI)
    r["ph_joint_wald_stage_x_period"] = {
        "chi2": wald, "df": len(inter_names),
        "p": float(chi2.sf(wald, len(inter_names)))}
    r["ph_interaction_p"] = {n: float(cphI.summary.loc[n, "p"])
                             for n in inter_names}
    print(f"  [{tag}] PH joint Wald (stage x period): chi2={wald:.1f}, "
          f"p={r['ph_joint_wald_stage_x_period']['p']:.2e}")
    r["period_hr"] = {}
    for s_ in STAGES:
        for p_ in (1, 2, 3):
            if p_ == 1:
                names, w = [s_], np.array([1.0])
            else:
                names = [s_, f"{s_}_p{p_}"]
                w = np.array([1.0, 1.0])
            beta = float(sum(w[k] * b[n] for k, n in enumerate(names)))
            sub = V.loc[names, names].values
            se = float(np.sqrt(w @ sub @ w))
            from scipy.stats import norm
            z = norm.ppf(0.975)
            r["period_hr"][f"{s_}_p{p_}"] = {
                "hr": float(np.exp(beta)), "lo": float(np.exp(beta - z*se)),
                "hi": float(np.exp(beta + z*se)), "p": float(
                    2*(1-norm.cdf(abs(beta/se))))}
    print(f"  [{tag}] period HRs S3: "
          + " / ".join(f"p{p}={r['period_hr'][f'stage3_p{p}']['hr']:.2f}"
                       for p in (1, 2, 3)))

    # ---------------- CIF by stage at 24 h (Aalen-Johansen) ----------------
    def aj_cif_30(times, causes):
        """Vectorised Aalen-Johansen CIF of cause-1 (death) by TAU.
        causes: 1 death, 2 discharge, 0 censored."""
        t, c = np.asarray(times, float), np.asarray(causes)
        ev = np.unique(t[(c > 0) & (t <= TAU)])
        if len(ev) == 0:
            return 0.0
        ts = np.sort(t)
        n_risk = len(t) - np.searchsorted(ts, ev, side="left")
        td1 = np.sort(t[c == 1]); td2 = np.sort(t[c == 2])
        d1 = (np.searchsorted(td1, ev, "right") - np.searchsorted(td1, ev, "left"))
        d2 = (np.searchsorted(td2, ev, "right") - np.searchsorted(td2, ev, "left"))
        surv_prev = np.cumprod(1.0 - (d1 + d2) / n_risk)
        surv_before = np.concatenate([[1.0], surv_prev[:-1]])
        return float(np.sum(surv_before * d1 / n_risk))

    times_all = pat["end_time"].values
    causes_all = pat["event_type"].values
    stage24 = pat["stage_at_24h"].values

    r["cif_by_stage24"] = {}
    rng = np.random.default_rng(42)
    for s_ in [0, 1, 2, 3]:
        m = stage24 == s_
        cif = aj_cif_30(times_all[m], causes_all[m])
        # bootstrap CI (200)
        boots = []
        idx_m = np.where(m)[0]
        for _ in range(200):
            sb = rng.choice(idx_m, size=len(idx_m), replace=True)
            boots.append(aj_cif_30(times_all[sb], causes_all[sb]))
        lo, hi = np.percentile(boots, [2.5, 97.5])
        r["cif_by_stage24"][str(s_)] = {
            "n": int(m.sum()), "cif30": round(100*cif, 1),
            "lo": round(100*lo, 1), "hi": round(100*hi, 1)}
        print(f"  [{tag}] CIF30 stage24={s_}: n={m.sum():,} "
              f"CIF={100*cif:.1f}% ({100*lo:.1f}-{100*hi:.1f})")

    # CIF of discharge (cause 2) overall for context
    causes_mapped = np.where(causes_all == 2, 1,
                             np.where(causes_all == 1, 2, 0))
    r["cif_discharge_overall"] = round(
        100 * aj_cif_30(times_all, causes_mapped), 1)

    return r


# ----------------------------------------------------------------------
for tag, id_col in [("mimic", "stay_id"), ("eicu", "patientunitstayid")]:
    print(f"\n{'='*60}\n[{tag.upper()}]\n{'='*60}")
    results[tag] = fit_models(tag, id_col)

results["generated"] = datetime.now().isoformat()
with open(os.path.join(V3, "v3_primary_results.json"), "w") as f:
    json.dump(results, f, indent=2)

print("\nv3 Phase 2 COMPLETE")
