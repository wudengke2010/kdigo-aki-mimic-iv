# -*- coding: utf-8 -*-
"""
Figure 1 (v2 cohort flow diagram) + Table S1 (baseline SCr source distribution)
v2.1 reanalysis pipeline — all numbers verified against cohort flow audit
(2026-09-23): MIMIC window-SCr exclusion = 473 (not 155); eICU = 2,031 (not -1,806).

Outputs:
  v2_outputs/fig1_flowchart_v2.png / .pdf
  v2_outputs/table_S1.md
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

WORK = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(WORK, "v2_outputs")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- Figure 1 --
# Verified exclusion chain (audit_report_v4.md + cohort flow JSONs)
# ('main', text) = cohort box on the main chain; ('excl', text) = exclusion box
MIMIC = [
    ("main", "ICU stays screened\n94,458 stays / 65,366 patients"),
    ("excl", "Excluded:\n3,991 stays with ICU LOS < 12 h"),
    ("main", "First eligible ICU stay per patient\n63,487 patients"),
    ("excl", "Excluded:\n26,980 repeat ICU stays"),
    ("excl", "Excluded:\n2,083 ESRD / chronic dialysis"),
    ("excl", "Excluded:\n898 without usable SCr\n(425 no baseline;\n473 no in-window SCr)"),
    ("main", "Final cohort\n60,506 patients\nStage 0/1/2/3 =\n43,011/11,882/2,249/3,364"),
    ("excl", "Excluded:\n1,952 died or discharged < 24 h"),
    ("main", "24-h landmark analysis set\n58,554 patients"),
]
EICU = [
    ("main", "Unit stays screened\n200,859 stays / 139,367 patients\n(age > 89 y recoded to 90, n = 7,081)"),
    ("excl", "Excluded:\n625 aged < 18 y or age missing"),
    ("excl", "Excluded:\n27,842 stays with unit LOS < 12 h"),
    ("main", "First eligible unit stay per patient\n130,361 patients"),
    ("excl", "Excluded:\n42,031 repeat unit stays"),
    ("excl", "Excluded:\n6,523 ESRD / chronic dialysis\n(ICD n = 2,801; pastHistory n = 6,219;\noverlap n = 2,497)"),
    ("excl", "Excluded:\n10,372 without usable SCr\n(8,341 no baseline;\n2,031 no in-window SCr)"),
    ("main", "Final cohort\n113,466 patients\nStage 0/1/2/3 =\n95,203/11,999/1,874/4,390"),
    ("excl", "Excluded:\n5,708 died or discharged < 24 h"),
    ("main", "24-h landmark analysis set\n107,758 patients"),
]

ROW_H = 0.92
X_MIMIC, X_EXCL_M = 3.1, 5.25          # MIMIC main chain x, exclusion box left x
X_EICU, X_EXCL_E = 10.6, 12.75
COL_HALF_W = 1.75                      # main box half width
BOX_MAIN = dict(boxstyle="round,pad=0.40", fc="#FFFFFF", ec="#1a3d5c", lw=1.4)
BOX_EXCL = dict(boxstyle="round,pad=0.30", fc="#F5F7FA", ec="#8a97a5", lw=1.0, ls="--")


def draw_column(ax, items, x_main, x_excl):
    ys = {}
    mains = []
    for i, (kind, text) in enumerate(items):
        y = -i * ROW_H
        ys[i] = y
        if kind == "main":
            bold = i == 0 or i == len(items) - 1 or text.startswith("Final")
            ax.text(x_main, y, text, ha="center", va="center", fontsize=8.4,
                    linespacing=1.35, fontweight="bold" if bold else "normal",
                    bbox=BOX_MAIN)
            mains.append(i)
        else:
            ax.text(x_excl, y, text, ha="left", va="center", fontsize=7.4,
                    linespacing=1.3, color="#4a5568", bbox=BOX_EXCL)
            # horizontal arrow from the previous main box to this exclusion box
            prev_m = max(m for m in mains)
            ax.add_patch(FancyArrowPatch((x_main + COL_HALF_W - 0.02, ys[prev_m]),
                                         (x_excl - 0.10, y),
                                         arrowstyle="-|>", mutation_scale=9,
                                         color="#8a97a5", lw=0.9))
    # vertical arrows between consecutive main boxes
    for a, b in zip(mains[:-1], mains[1:]):
        ax.add_patch(FancyArrowPatch((x_main, ys[a] - 0.34), (x_main, ys[b] + 0.34),
                                     arrowstyle="-|>", mutation_scale=12,
                                     color="#1a3d5c", lw=1.3))
    return min(ys.values())


n_rows = max(len(MIMIC), len(EICU))
fig, ax = plt.subplots(figsize=(14.5, n_rows * ROW_H + 2.4))
ax.set_xlim(0, 16.2)
ax.set_ylim(-(n_rows - 1) * ROW_H - 1.1, 1.15)
ax.axis("off")

ax.text(X_MIMIC, 0.85, "MIMIC-IV v3.1 (derivation)", ha="center", fontsize=12.5,
        fontweight="bold", color="#1a3d5c")
ax.text(X_EICU, 0.85, "eICU-CRD v2.0 (validation)", ha="center", fontsize=12.5,
        fontweight="bold", color="#1a3d5c")

draw_column(ax, MIMIC, X_MIMIC, X_EXCL_M)
draw_column(ax, EICU, X_EICU, X_EXCL_E)

fig.suptitle("Figure 1. Patient-level cohort construction and 24-hour landmark exclusions",
             fontsize=13.5, fontweight="bold", y=0.99)
fig.text(0.5, 0.012,
         "Patient-level cohorts: only the first eligible ICU stay per patient was retained. "
         "Dashed boxes denote exclusions.\nKDIGO = Kidney Disease: Improving Global Outcomes; "
         "LOS = length of stay; ESRD = end-stage renal disease; SCr = serum creatinine.",
         ha="center", fontsize=8.2, color="#4a5568")

plt.savefig(os.path.join(OUT, "fig1_flowchart_v2.png"), dpi=350, bbox_inches="tight")
plt.savefig(os.path.join(OUT, "fig1_flowchart_v2.pdf"), bbox_inches="tight")
plt.close()
print("[Fig1] saved fig1_flowchart_v2.png/.pdf")

# --------------------------------------------------------------- Table S1 ---
mimic = pd.read_csv(os.path.join(WORK, "mimic_analysis_v2.csv"))
eicu = pd.read_csv(os.path.join(WORK, "eicu_analysis_v2.csv"))

MIMIC_SRC = {
    "prior_0_7d": "Pre-admission SCr, 0–7 d before admission",
    "prior_any_7_365d": "Pre-admission SCr (any setting), 7–365 d before admission",
    "prior_outpt_7_365d": "Outpatient SCr, 7–365 d before admission",
    "hosp_first24h": "First in-hospital SCr (first 24 h)",
}
EICU_SRC = {
    "preICU_0_7d": "Pre-ICU SCr, 0–7 d before ICU admission",
    "icu_first24h": "First SCr in ICU (first 24 h)",
}

def table_block(df, src_map, cohort_name, pre_prefixes):
    lines = [f"**{cohort_name} (N = {len(df):,})**"]
    lines.append("| Baseline SCr source | n (%) | Stage 0 | Stage 1 | Stage 2 | Stage 3 | AKI % | Hospital mortality % |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for key, label in src_map.items():
        sub = df[df.baseline_source == key]
        n = len(sub)
        cells = []
        for s in [0, 1, 2, 3]:
            ns = (sub.kdigo_stage == s).sum()
            cells.append(f"{ns:,} ({100 * ns / n:.1f}%)")
        aki = 100 * (sub.kdigo_stage > 0).mean()
        mort = 100 * sub.hospital_expire_flag.mean()
        lines.append(f"| {label} | {n:,} ({100 * n / len(df):.1f}%) | " +
                     " | ".join(cells) + f" | {aki:.1f} | {mort:.1f} |")
    pre = df[df.baseline_source.isin(pre_prefixes)]
    lines.append(f"| *Any pre-admission baseline (pooled)* | *{len(pre):,} ({100*len(pre)/len(df):.1f}%)* | " +
                 " | ".join(f"*{(pre.kdigo_stage==s).sum():,}*"
                            for s in [0, 1, 2, 3]) +
                 f" | *{100*(pre.kdigo_stage>0).mean():.1f}* | *{100*pre.hospital_expire_flag.mean():.1f}* |")
    lines.append("")
    return lines

doc = ["# Supplementary Table S1",
       "",
       "**Table S1. Baseline serum creatinine (SCr) source distribution by cohort and KDIGO stage.**",
       "",
       "Values are n (% within baseline-source row). AKI = KDIGO Stage 1–3. Hospital mortality is the "
       "in-hospital mortality of the full (pre-landmark) cohort. The pre-admission vs first-ICU-period "
       "asymmetry between cohorts (MIMIC-IV 80.2% vs eICU-CRD 39.4% pre-admission) is the principal "
       "explanation for the lower AKI prevalence in eICU-CRD, because a first-ICU creatinine baseline "
       "inflates the denominator and suppresses ratio-defined AKI.",
       ""]
doc += table_block(mimic, MIMIC_SRC, "MIMIC-IV v3.1 (derivation cohort)",
                   ["prior_0_7d", "prior_any_7_365d", "prior_outpt_7_365d"])
doc += table_block(eicu, EICU_SRC, "eICU-CRD v2.0 (validation cohort)",
                   ["preICU_0_7d"])

with open(os.path.join(OUT, "table_S1.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(doc))
print("[TableS1] saved table_S1.md")
