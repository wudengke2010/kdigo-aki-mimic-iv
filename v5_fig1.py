# -*- coding: utf-8 -*-
"""Figure 1 for v5 manuscript: cohort flow WITHOUT the 24-hour landmark
(the v3 primary analysis follows patients from ICU admission).
Stage counts updated to the v3 time-varying implementation.
Output: v3_outputs/fig_v5_flowchart.png / .pdf
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

WORK = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(WORK, "v3_outputs")

MIMIC = [
    ("main", "ICU stays screened\n94,458 stays / 65,366 patients"),
    ("excl", "Excluded:\n3,991 stays with ICU LOS < 12 h"),
    ("main", "First eligible ICU stay per patient\n63,487 patients"),
    ("excl", "Excluded:\n26,980 repeat ICU stays"),
    ("excl", "Excluded:\n2,083 ESRD / chronic dialysis"),
    ("excl", "Excluded:\n898 without usable SCr\n(425 no baseline;\n473 no in-window SCr)"),
    ("main", "Final cohort — followed from\nICU admission (30-day follow-up)\n60,506 patients\nStage 0/1/2/3 =\n43,014/11,879/2,249/3,364"),
]
EICU = [
    ("main", "Unit stays screened\n200,859 stays / 139,367 patients\n(age > 89 y recoded to 90, n = 7,081)"),
    ("excl", "Excluded:\n625 aged < 18 y or age missing"),
    ("excl", "Excluded:\n27,842 stays with unit LOS < 12 h"),
    ("main", "First eligible unit stay per patient\n130,361 patients"),
    ("excl", "Excluded:\n42,031 repeat unit stays"),
    ("excl", "Excluded:\n6,523 ESRD / chronic dialysis\n(ICD n = 2,801; pastHistory n = 6,219;\noverlap n = 2,497)"),
    ("excl", "Excluded:\n10,372 without usable SCr\n(8,341 no baseline;\n2,031 no in-window SCr)"),
    ("main", "Final cohort — followed from\nICU admission (30-day follow-up)\n113,466 patients\nStage 0/1/2/3 =\n95,220/11,983/1,874/4,389"),
]

ROW_H = 0.92
X_MIMIC, X_EXCL_M = 3.1, 5.25
X_EICU, X_EXCL_E = 10.6, 12.75
COL_HALF_W = 1.75
BOX_MAIN = dict(boxstyle="round,pad=0.40", fc="#FFFFFF", ec="#1a3d5c", lw=1.4)
BOX_EXCL = dict(boxstyle="round,pad=0.30", fc="#F5F7FA", ec="#8a97a5", lw=1.0, ls="--")


def draw_column(ax, items, x_main, x_excl):
    ys = {}
    mains = []
    for i, (kind, text) in enumerate(items):
        y = -i * ROW_H
        ys[i] = y
        if kind == "main":
            bold = i == 0 or text.startswith("Final")
            ax.text(x_main, y, text, ha="center", va="center", fontsize=8.4,
                    linespacing=1.35, fontweight="bold" if bold else "normal",
                    bbox=BOX_MAIN)
            mains.append(i)
        else:
            ax.text(x_excl, y, text, ha="left", va="center", fontsize=7.4,
                    linespacing=1.3, color="#4a5568", bbox=BOX_EXCL)
            prev_m = max(m for m in mains)
            ax.add_patch(FancyArrowPatch((x_main + COL_HALF_W - 0.02, ys[prev_m]),
                                         (x_excl - 0.10, y),
                                         arrowstyle="-|>", mutation_scale=9,
                                         color="#8a97a5", lw=0.9))
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

fig.suptitle("Figure 1. Patient-level cohort construction (follow-up from ICU admission; no landmark exclusion)",
             fontsize=13.5, fontweight="bold", y=0.99)
fig.text(0.5, 0.012,
         "Patient-level cohorts: only the first eligible ICU stay per patient was retained. "
         "Dashed boxes denote exclusions.\nKDIGO = Kidney Disease: Improving Global Outcomes; "
         "LOS = length of stay; ESRD = end-stage renal disease; SCr = serum creatinine.",
         ha="center", fontsize=8.2, color="#4a5568")

fig.savefig(os.path.join(OUT, "fig_v5_flowchart.png"), dpi=350, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig_v5_flowchart.pdf"), bbox_inches="tight")
plt.close(fig)
print("[Fig1 v5] saved")
