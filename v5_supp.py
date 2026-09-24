# -*- coding: utf-8 -*-
"""v5 supplementary: Table 1 baseline characteristics by max-7d stage,
subgroup x stage3 interaction analysis under the v3 time-varying framework,
and supplementary forest figures (subgroup + phenotype).

Outputs: v3_outputs/v5_table1.json, v3_outputs/v5_subgroup.json,
         v3_outputs/fig_v5_subgroup_forest.png/pdf,
         v3_outputs/fig_v5_phenotype_forest.png/pdf
"""
import os
import json
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
V3 = os.path.join(BASE, "v3_outputs")
PY = "C:/Users/admin/.workbuddy/binaries/python/versions/3.13.12/python.exe"

COVARS = ["age", "sex_male", "sofa_nonrenal", "vent_24h", "ckd_history_prior",
          "diabetes", "hypertension", "heart_failure", "copd", "liver_disease"]

MPL = {"font.size": 8, "font.family": "Arial",
       "axes.edgecolor": "#444444", "axes.linewidth": 0.8}
matplotlib.rcParams.update(MPL)

# ----------------------------------------------------------------- Table 1
def table1(tag):
    pat = pd.read_csv(os.path.join(V3, f"v3_{tag}_patients.csv"))
    out = {}
    for s in [0, 1, 2, 3]:
        d = pat[pat["stage_max_7d"] == s]
        n = len(d)
        death30 = ((d["event_type"] == 1)).mean() * 100  # event_type 1 = death<=30d
        row = {
            "n": n,
            "age": round(d["age"].mean(), 1),
            "male": round(d["sex_male"].mean() * 100, 1),
            "sofa": round(d["sofa_nonrenal"].mean(), 2),
            "vent": round(d["vent_24h"].mean() * 100, 1),
            "ckd": round(d["ckd_history_prior"].mean() * 100, 1),
            "dm": round(d["diabetes"].mean() * 100, 1),
            "htn": round(d["hypertension"].mean() * 100, 1),
            "hf": round(d["heart_failure"].mean() * 100, 1),
            "copd": round(d["copd"].mean() * 100, 1),
            "liver": round(d["liver_disease"].mean() * 100, 1),
            "death30": round(death30, 1),
        }
        out[s] = row
    return out

t1 = {tag: table1(tag) for tag in ["mimic", "eicu"]}
json.dump(t1, open(os.path.join(V3, "v5_table1.json"), "w"), indent=1)
print("Table 1 computed")

# ------------------------------------------------------- subgroup analysis
def subgroups(tag):
    tv = pd.read_csv(os.path.join(V3, f"v3_{tag}_tv.csv.gz"))
    pat = pd.read_csv(os.path.join(V3, f"v3_{tag}_patients.csv"))
    idc = "stay_id" if tag == "mimic" else "patientunitstayid"
    tv = tv.merge(pat[[idc] + COVARS], on=idc, how="left")

    sofa_t1, sofa_t2 = pat["sofa_nonrenal"].quantile([1/3, 2/3])

    subs = {
        "Age >= 65": pat["age"] >= 65,
        "Age < 65": pat["age"] < 65,
        "Male": pat["sex_male"] == 1,
        "Female": pat["sex_male"] == 0,
        "Prior CKD": pat["ckd_history_prior"] == 1,
        "No prior CKD": pat["ckd_history_prior"] == 0,
        "Diabetes": pat["diabetes"] == 1,
        "No diabetes": pat["diabetes"] == 0,
        "Ventilated (24h)": pat["vent_24h"] == 1,
        "Not ventilated": pat["vent_24h"] == 0,
        "Non-renal SOFA high tertile": pat["sofa_nonrenal"] > sofa_t2,
        "Non-renal SOFA low/mid": pat["sofa_nonrenal"] <= sofa_t2,
    }
    # subgroup membership map (per patient)
    submap = {k: dict(zip(pat[idc], v)) for k, v in subs.items()}

    res = {}
    # interaction pairs: (labelA, labelB, variable)
    pairs = [("Age >= 65", "Age < 65", "age65"),
             ("Male", "Female", "male"),
             ("Prior CKD", "No prior CKD", "ckd"),
             ("Diabetes", "No diabetes", "dm"),
             ("Ventilated (24h)", "Not ventilated", "vent"),
             ("Non-renal SOFA high tertile", "Non-renal SOFA low/mid", "sofa_hi")]
    vardef = {
        "age65": (pat["age"] >= 65).astype(int),
        "male": pat["sex_male"],
        "ckd": pat["ckd_history_prior"],
        "dm": pat["diabetes"],
        "vent": pat["vent_24h"],
        "sofa_hi": (pat["sofa_nonrenal"] > sofa_t2).astype(int),
    }
    varmap = {k: dict(zip(pat[idc], v)) for k, v in vardef.items()}

    for la, lb, var in pairs:
        # subgroup-specific S3 HR: fit Model C within each subgroup
        hr = {}
        for lab in (la, lb):
            keep = set(pat.loc[subs[lab], idc])
            dd = tv[tv[idc].isin(keep)].copy()
            cols = ["stage1", "stage2", "stage3"] + COVARS
            # drop zero-variance covariates within subgroup (avoids collinearity)
            vary = [c for c in cols if dd[c].std() > 0]
            dd = dd[["t_start", "t_stop", "event_death"] + vary].sort_values("t_stop")
            cph = CoxPHFitter()
            cph.fit(dd, duration_col="t_stop", event_col="event_death",
                    entry_col="t_start", show_progress=False)
            s3 = cph.summary.loc["stage3"]
            n_pat = len(keep)
            hr[lab] = {"n": n_pat, "hr": round(float(np.exp(s3["coef"])), 2),
                       "lo": round(float(s3["exp(coef) lower 95%"]), 2),
                       "hi": round(float(s3["exp(coef) upper 95%"]), 2)}
        # interaction test: stage3 x var (+ all stage x var for completeness)
        dd = tv.copy()
        dd["_v"] = dd[idc].map(varmap[var])
        for s_ in ["stage1", "stage2", "stage3"]:
            dd[f"{s_}x"] = dd[s_] * dd["_v"]
        # remove the corresponding main covariate to avoid collinearity
        main_col = {"age65": "age", "male": "sex_male", "ckd": "ckd_history_prior",
                    "dm": "diabetes", "vent": "vent_24h", "sofa_hi": "sofa_nonrenal"}[var]
        cols = [c for c in ["stage1", "stage2", "stage3"] + COVARS if c != main_col]
        cols = cols + ["_v", "stage1x", "stage2x", "stage3x"]
        dd = dd[["t_start", "t_stop", "event_death"] + cols].sort_values("t_stop")
        cph = CoxPHFitter()
        cph.fit(dd, duration_col="t_stop", event_col="event_death",
                entry_col="t_start", show_progress=False)
        b = cph.params_
        V = cph.variance_matrix_
        # joint Wald for stage3 interaction alone and all three
        def wald(names):
            bb = np.array([b[n] for n in names])
            VV = V.loc[names, names].values
            chi2 = float(bb @ np.linalg.pinv(VV) @ bb)
            return chi2
        chi2_3 = wald(["stage3x"])
        from scipy import stats as st
        p3 = float(1 - st.chi2.cdf(chi2_3, 1))
        res[var] = {"subgroups": hr, "interaction_p_s3": round(p3, 4)}
        print(f"  [{tag}] {var}: {la} S3={hr[la]['hr']} vs {lb} S3={hr[lb]['hr']}, p_int={p3:.4f}")
    return res

sg = {tag: subgroups(tag) for tag in ["mimic", "eicu"]}
json.dump(sg, open(os.path.join(V3, "v5_subgroup.json"), "w"), indent=1)
print("Subgroups computed")

# ------------------------------------------------------ subgroup forest fig
fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.2), sharey=True)
for ai, (tag, ax) in enumerate(zip(["mimic", "eicu"], axes)):
    labels, hrs, los, his = [], [], [], []
    for var, pair in [("age65", ("Age >= 65", "Age < 65")),
                      ("male", ("Male", "Female")),
                      ("ckd", ("Prior CKD", "No prior CKD")),
                      ("dm", ("Diabetes", "No diabetes")),
                      ("vent", ("Ventilated (24h)", "Not ventilated")),
                      ("sofa_hi", ("SOFA high tertile", "SOFA low/mid"))]:
        for lab in pair:
            h = sg[tag][var]["subgroups"][lab]
            labels.append(lab)
            hrs.append(h["hr"]); los.append(h["lo"]); his.append(h["hi"])
    y = np.arange(len(labels))[::-1]
    ax.errorbar(hrs, y, xerr=[np.array(hrs) - np.array(los),
                              np.array(his) - np.array(hrs)],
                fmt="s", ms=3, color="#2F5597", ecolor="#2F5597",
                elinewidth=0.9, capsize=2)
    ax.axvline(1, color="#999999", lw=0.7, ls="--")
    ax.set_yticks(y); ax.set_yticklabels(labels if ai == 0 else [""] * len(labels))
    ax.set_xlabel("Adjusted Stage 3 HR (95% CI)", fontsize=8)
    ax.set_title("MIMIC-IV" if tag == "mimic" else "eICU-CRD", fontsize=9)
    ax.tick_params(length=2.5, direction="in")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(V3, "fig_v5_subgroup_forest.png"), dpi=600)
fig.savefig(os.path.join(V3, "fig_v5_subgroup_forest.pdf"))
plt.close(fig)
print("subgroup forest saved")

# ----------------------------------------------------- phenotype forest fig
pheno = json.load(open(os.path.join(V3, "v3_phenotype_results.json")))
fig, ax = plt.subplots(figsize=(4.6, 2.6))
rows = [("MIMIC-IV\npersistent vs transient", pheno["mimic"]),
        ("eICU-CRD\npersistent vs transient", pheno["eicu"])]
labels, hrs, los, his = [], [], [], []
for lab, d in rows:
    a = d["cox_phenotype_adjusted"]
    labels.append(lab); hrs.append(a["hr"]); los.append(a["lo"]); his.append(a["hi"])
y = np.arange(len(labels))[::-1]
ax.errorbar(hrs, y, xerr=[np.array(hrs) - np.array(los),
                          np.array(his) - np.array(hrs)],
            fmt="s", ms=4, color="#C00000", ecolor="#C00000",
            elinewidth=1, capsize=2.5)
ax.axvline(1, color="#999999", lw=0.7, ls="--")
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
ax.set_xlabel("Adjusted HR (95% CI), onset+72h landmark", fontsize=8)
ax.tick_params(length=2.5, direction="in")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(V3, "fig_v5_phenotype_forest.png"), dpi=600)
fig.savefig(os.path.join(V3, "fig_v5_phenotype_forest.pdf"))
plt.close(fig)
print("phenotype forest saved")
print("ALL DONE")
