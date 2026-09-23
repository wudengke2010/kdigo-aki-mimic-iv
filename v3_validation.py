# -*- coding: utf-8 -*-
"""
v3 Phase 4: locked external validation (TRIPOD type 3).

Target: 30-day IN-HOSPITAL death (CIF estimand, consistent with the
competing-risk primary analysis):
    event     = death <= 720 h in hospital
    non-event = discharge alive <= 720 h
    censored  = still hospitalised at 720 h (excluded from the core analysis;
                IPCW sensitivity included — censoring is mild: 1.5-3.7%)

Model: logistic regression, stage at 24 h (dummies; first-24-h data only ->
no immortal time) + baseline covariates. LOCKED on MIMIC; applied to eICU
with ZERO refitting.

Metrics: AUC, Brier score, calibration slope (logistic on logit(p)),
calibration deciles, E-axis O:E ratio; MIMIC 5-fold CV for optimism.

Outputs: v3_outputs/v3_validation_results.json, fig_v3_calibration.png/pdf
"""
import pandas as pd
import numpy as np
import os, json
from datetime import datetime
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, brier_score_loss

OUT = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
V3 = os.path.join(OUT, "v3_outputs")

COVARS = ["age", "sex_male", "sofa_nonrenal", "vent_24h",
          "ckd_history_prior", "diabetes", "hypertension",
          "heart_failure", "copd", "liver_disease"]
PRED = ["s1_24", "s2_24", "s3_24"] + COVARS

print("=" * 70)
print(f"v3 Phase 4: locked external validation  Start: {datetime.now()}")
print("=" * 70)


def prep(tag):
    pat = pd.read_csv(os.path.join(V3, f"v3_{tag}_patients.csv"))
    pat["s1_24"] = (pat["stage_at_24h"] == 1).astype(int)
    pat["s2_24"] = (pat["stage_at_24h"] == 2).astype(int)
    pat["s3_24"] = (pat["stage_at_24h"] == 3).astype(int)
    # In-hospital death within 30 d: FULLY OBSERVED binary outcome.
    #  death <= 720 h            -> y = 1
    #  discharged alive <= 720 h -> y = 0 (no in-hospital death)
    #  still hospitalised at 720 -> y = 0 (alive at day 30 by definition)
    pat["y"] = (pat["event_type"] == 1).astype(int)
    return pat


mimic = prep("mimic")
eicu = prep("eicu")

res = {}
for tag, d in [("mimic", mimic), ("eicu", eicu)]:
    res[f"{tag}_n"] = int(len(d))
    res[f"{tag}_events"] = int(d["y"].sum())
    res[f"{tag}_nonevents"] = int((d["y"] == 0).sum())
    res[f"{tag}_still_hosp_720"] = int((d["event_type"] == 0).sum())

# --- LOCKED model: fit on MIMIC (full cohort, fully observed outcome) ----
tr = mimic.dropna(subset=PRED)
model = LogisticRegression(penalty=None, max_iter=2000)
model.fit(tr[PRED].values, tr["y"].values)
coef = model.coef_[0]
intercept = float(model.intercept_[0])
res["locked_coefs"] = {p: float(c) for p, c in zip(PRED, coef)}
res["locked_intercept"] = intercept
print(f"Locked logistic: n={len(tr):,}, events={int(tr['y'].sum()):,}")
for p, c in zip(PRED, coef):
    print(f"    {p:>18s}: {c:+.4f}")
print(f"    {'intercept':>18s}: {intercept:+.4f}")


def evaluate(d, tag, label):
    dd = d.dropna(subset=PRED)
    X = dd[PRED].values
    y = dd["y"].values
    p = model.predict_proba(X)[:, 1]   # LOCKED coefficients, zero refit
    auc = roc_auc_score(y, p)
    brier = brier_score_loss(y, p)
    # calibration slope: logistic of y on logit(p) (free slope+intercept)
    from scipy.special import logit
    lp = logit(np.clip(p, 1e-6, 1 - 1e-6))
    import statsmodels.api as sm
    fitc = sm.GLM(y, sm.add_constant(lp), family=sm.families.Binomial()).fit()
    slope = float(fitc.params[1])
    slope_ci = [float(x) for x in fitc.conf_int()[1]]
    oe = y.sum() / p.sum()
    out = {"n": int(len(dd)), "events": int(y.sum()), "auc": round(auc, 3),
           "brier": round(brier, 4), "calib_slope": round(slope, 3),
           "calib_slope_ci": [round(x, 3) for x in slope_ci],
           "oe_ratio": round(oe, 3)}
    # deciles
    q = pd.qcut(p, 10, duplicates="drop")
    dec = pd.DataFrame({"p": p, "y": y, "q": q}).groupby(
        "q", observed=True).agg(p_mean=("p", "mean"), y_mean=("y", "mean"),
                                n=("p", "size"))
    out["deciles"] = [[round(r.p_mean, 4), round(r.y_mean, 4), int(r.n)]
                      for r in dec.itertuples()]
    print(f"  [{label}] n={out['n']:,} events={out['events']:,} "
          f"AUC={out['auc']:.3f} Brier={out['brier']:.4f} "
          f"slope={out['calib_slope']} O:E={out['oe_ratio']}")
    return out, p, y, dd


res["mimic_apparent"], pm, ym, _ = evaluate(mimic, "mimic", "MIMIC apparent")
res["eicu_locked"], pe, ye, _ = evaluate(eicu, "eicu", "eICU locked")

# --- MIMIC internal optimism: 5-fold CV AUC ---
aucs = []
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
Xall, yall = tr[PRED].values, tr["y"].values
for tri, tei in skf.split(Xall, yall):
    m = LogisticRegression(penalty=None, max_iter=2000)
    m.fit(Xall[tri], yall[tri])
    aucs.append(roc_auc_score(yall[tei], m.predict_proba(Xall[tei])[:, 1]))
res["mimic_cv5_auc"] = round(float(np.mean(aucs)), 3)
res["mimic_cv5_auc_sd"] = round(float(np.std(aucs)), 4)
print(f"  MIMIC 5-fold CV AUC: {res['mimic_cv5_auc']} (sd {res['mimic_cv5_auc_sd']})")

# --- sensitivity: exclude patients still hospitalised at 720 h ----------
def exclude_still_hosp(d, label):
    dd = d[d["event_type"] != 0].dropna(subset=PRED)
    p = model.predict_proba(dd[PRED].values)[:, 1]
    y = dd["y"].values
    auc = roc_auc_score(y, p)
    print(f"  [{label}] excl-still-hospitalised: n={len(dd):,} AUC={auc:.3f}")
    return {"n": int(len(dd)), "auc": round(auc, 3)}

res["eicu_sens_excl_still_hosp"] = exclude_still_hosp(eicu, "eICU")
res["mimic_sens_excl_still_hosp"] = exclude_still_hosp(mimic, "MIMIC")

res["generated"] = datetime.now().isoformat()
with open(os.path.join(V3, "v3_validation_results.json"), "w") as f:
    json.dump(res, f, indent=2)

# --- calibration figure --------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
for ax, (p, y, title) in zip(axes, [
        (pm, ym, f"MIMIC-IV (derivation, apparent)"),
        (pe, ye, "eICU-CRD (validation, LOCKED model)")]):
    q = pd.qcut(p, 10, duplicates="drop")
    dec = pd.DataFrame({"p": p, "y": y, "q": q}).groupby(
        "q", observed=True).agg(pm=("p", "mean"), ym=("y", "mean"))
    ax.plot([0, 0.6], [0, 0.6], "k--", lw=1, label="ideal")
    ax.plot(dec.pm, dec.ym, "o-", color="#d62728", ms=5, label="observed")
    ax.set_xlabel("Predicted 30-day in-hospital death risk")
    ax.set_ylabel("Observed risk")
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=9)
    ax.set_xlim(0, 0.6); ax.set_ylim(0, 0.6)
fig.suptitle("Calibration: locked MIMIC model (stage at 24 h + covariates)",
             fontsize=12)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(V3, f"fig_v3_calibration.{ext}"), dpi=300,
                bbox_inches="tight")

print("\nv3 Phase 4 COMPLETE")
