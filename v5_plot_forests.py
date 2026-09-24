# -*- coding: utf-8 -*-
"""Plot-only: v5 subgroup + phenotype forest figures from saved JSONs."""
import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
V3 = os.path.join(BASE, "v3_outputs")
matplotlib.rcParams.update({"font.size": 8, "font.family": "Arial",
                            "axes.edgecolor": "#444444", "axes.linewidth": 0.8})

sg = json.load(open(os.path.join(V3, "v5_subgroup.json")))
pheno = json.load(open(os.path.join(V3, "v3_phenotype_results.json")))

# ------------------------------------------------------ subgroup forest fig
fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.4), sharey=True)
pairs = [("age65", ("Age >= 65", "Age < 65")),
         ("male", ("Male", "Female")),
         ("ckd", ("Prior CKD", "No prior CKD")),
         ("dm", ("Diabetes", "No diabetes")),
         ("vent", ("Ventilated (24h)", "Not ventilated")),
         ("sofa_hi", ("Non-renal SOFA high tertile", "Non-renal SOFA low/mid"))]
for ai, (tag, ax) in enumerate(zip(["mimic", "eicu"], axes)):
    labels, hrs, los, his = [], [], [], []
    for var, (la, lb) in pairs:
        for lab in (la, lb):
            h = sg[tag][var]["subgroups"][lab]
            labels.append(lab)
            hrs.append(h["hr"]); los.append(h["lo"]); his.append(h["hi"])
    y = np.arange(len(labels))[::-1]
    ax.errorbar(hrs, y, xerr=[np.array(hrs) - np.array(los),
                              np.array(his) - np.array(hrs)],
                fmt="s", ms=3, color="#2F5597", ecolor="#2F5597",
                elinewidth=0.9, capsize=2)
    ax.axvline(1, color="#999999", lw=0.7, ls="--")
    ax.set_yticks(y); ax.set_yticklabels(labels if ai == 0 else [""] * len(labels),
                                         fontsize=7.2)
    ax.set_xlabel("Adjusted Stage 3 HR (95% CI)", fontsize=8)
    ax.set_title("MIMIC-IV" if tag == "mimic" else "eICU-CRD", fontsize=9)
    ax.tick_params(length=2.5, direction="in")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.set_xlim(0.5, 4.5)
fig.tight_layout()
fig.savefig(os.path.join(V3, "fig_v5_subgroup_forest.png"), dpi=600)
fig.savefig(os.path.join(V3, "fig_v5_subgroup_forest.pdf"))
plt.close(fig)
print("subgroup forest saved")

# ----------------------------------------------------- phenotype forest fig
fig, ax = plt.subplots(figsize=(4.6, 2.4))
rows = [("MIMIC-IV", pheno["mimic"]), ("eICU-CRD", pheno["eicu"])]
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
ax.set_xlabel("Persistent vs transient AKI: adjusted HR (95% CI)", fontsize=8)
ax.tick_params(length=2.5, direction="in")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(V3, "fig_v5_phenotype_forest.png"), dpi=600)
fig.savefig(os.path.join(V3, "fig_v5_phenotype_forest.pdf"))
plt.close(fig)
print("phenotype forest saved")

# dump full subgroup numbers for manuscript
for tag in ["mimic", "eicu"]:
    print(f"=== {tag} ===")
    for var, (la, lb) in pairs:
        d = sg[tag][var]
        print(f"{var}: {la} {d['subgroups'][la]['hr']} ({d['subgroups'][la]['lo']}-{d['subgroups'][la]['hi']}, n={d['subgroups'][la]['n']}) vs "
              f"{lb} {d['subgroups'][lb]['hr']} ({d['subgroups'][lb]['lo']}-{d['subgroups'][lb]['hi']}, n={d['subgroups'][lb]['n']}), p_int={d['interaction_p_s3']}")
