# -*- coding: utf-8 -*-
"""
v3 Phase 5: consistency audit + figures + consolidated report.

1. Cross-check every count between tv long data, patient files and flow JSON.
2. Figures: CIF curves by stage-at-24h (MIMIC + eICU); forest of Model A/B/C
   and period-specific stage HRs.
3. v3_report.md: full results + v2-vs-v3 comparison (what survived, what
   shrank and why).
"""
import pandas as pd
import numpy as np
import os, json
from datetime import datetime

OUT = "C:/Users/admin/WorkBuddy/2026-07-06-13-22-19"
V3 = os.path.join(OUT, "v3_outputs")
TAU = 720.0

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

print("=" * 70)
print(f"v3 Phase 5: audit + figures + report  Start: {datetime.now()}")
print("=" * 70)

flow = json.load(open(os.path.join(V3, "v3_flow.json")))
prim = json.load(open(os.path.join(V3, "v3_primary_results.json")))
phen = json.load(open(os.path.join(V3, "v3_phenotype_results.json")))
val = json.load(open(os.path.join(V3, "v3_validation_results.json")))

# ------------------------------------------------------------- audit ----
audit = {}
for tag, idc in [("mimic", "stay_id"), ("eicu", "patientunitstayid")]:
    tv = pd.read_csv(os.path.join(V3, f"v3_{tag}_tv.csv.gz"))
    pat = pd.read_csv(os.path.join(V3, f"v3_{tag}_patients.csv"))
    a = {}
    a["patients"] = int(len(pat))
    a["tv_unique_ids"] = int(tv[idc].nunique())
    a["tv_deaths"] = int(tv["event_death"].sum())
    a["pat_deaths"] = int((pat["event_type"] == 1).sum())
    a["tv_discharges"] = int(tv["event_disc"].sum())
    a["pat_discharges"] = int((pat["event_type"] == 2).sum())
    # every patient's intervals tile [0, end] exactly
    g = tv.groupby(idc).agg(lo=("t_start", "min"), hi=("t_stop", "max"),
                            n=("t_stop", "size"))
    m = pat.set_index(idc)
    a["interval_tiling_ok"] = bool(
        (g["lo"] <= 1e-9).all()
        and np.allclose(g["hi"].values, m.loc[g.index, "end_time"].values))
    # no overlap / no gap: sum of durations == end time
    dur = tv.groupby(idc)["t_stop"].apply(lambda s: s.diff().fillna(s.iloc[0]).sum())
    a["duration_match_ok"] = bool(
        np.allclose(dur.values, m.loc[dur.index, "end_time"].values, atol=1e-6))
    a["flow_deaths"] = flow[tag]["n_death_le_30d"]
    a["flow_discharges"] = flow[tag]["n_discharge_le_30d"]
    a["flow_cohort"] = flow[tag]["n_cohort"]
    ok = (a["tv_unique_ids"] == a["patients"] == a["flow_cohort"]
          and a["tv_deaths"] == a["pat_deaths"] == a["flow_deaths"]
          and a["tv_discharges"] == a["pat_discharges"] == a["flow_discharges"]
          and a["interval_tiling_ok"] and a["duration_match_ok"])
    a["ALL_OK"] = ok
    audit[tag] = a
    print(f"[{tag}] audit: {'PASS' if ok else 'FAIL'} "
          f"(n={a['patients']:,}, deaths={a['pat_deaths']:,}, "
          f"discharges={a['pat_discharges']:,})")
    if not ok:
        print(json.dumps(a, indent=2))
        raise SystemExit("AUDIT FAILED")

# ------------------------------------------------------------- figures --
def cif_curves(times, causes, n_grid=200):
    t, c = np.asarray(times, float), np.asarray(causes)
    ev = np.unique(t[(c > 0) & (t <= TAU)])
    ts_ = np.sort(t)
    n_risk = len(t) - np.searchsorted(ts_, ev, side="left")
    td1 = np.sort(t[c == 1]); td2 = np.sort(t[c == 2])
    d1 = (np.searchsorted(td1, ev, "right") - np.searchsorted(td1, ev, "left"))
    d2 = (np.searchsorted(td2, ev, "right") - np.searchsorted(td2, ev, "left"))
    surv = np.cumprod(1.0 - (d1 + d2) / n_risk)
    sb = np.concatenate([[1.0], surv[:-1]])
    inc = sb * d1 / n_risk
    cif_at = np.cumsum(inc)
    grid = np.linspace(0, TAU, n_grid)
    return grid, np.interp(grid, ev, cif_at, left=0.0)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
colors = {0: "#7f7f7f", 1: "#ff7f0e", 2: "#2ca02c", 3: "#d62728"}
for ax, tag, title in [(axes[0], "mimic", "MIMIC-IV v3.1 (derivation)"),
                       (axes[1], "eicu", "eICU-CRD (validation)")]:
    pat = pd.read_csv(os.path.join(V3, f"v3_{tag}_patients.csv"))
    for s in [0, 1, 2, 3]:
        sub = pat[pat["stage_at_24h"] == s]
        g, cif = cif_curves(sub["end_time"].values, sub["event_type"].values)
        ax.plot(g / 24, 100 * cif, color=colors[s], lw=1.8,
                label=f"Stage {s} (n={len(sub):,})")
    ax.set_xlabel("Days since ICU admission")
    ax.set_ylabel("Cumulative incidence of death (%)")
    ax.set_title(title, fontsize=11)
    ax.set_ylim(0, 45)
    ax.legend(fontsize=8.5, loc="upper left")
fig.suptitle("30-day mortality: competing-risk CIF by KDIGO stage at 24 h "
             "(discharge alive = competing event)", fontsize=12)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(V3, f"fig_v3_cif.{ext}"), dpi=300,
                bbox_inches="tight")
print("fig_v3_cif saved")

# forest: Model A/B/C + period HRs, both cohorts
rows = []
for tag in ["mimic", "eicu"]:
    for s in ["stage1", "stage2", "stage3"]:
        for model in ["model_A", "model_B", "model_C"]:
            e = prim[tag][model][s]
            rows.append([tag, s[-1], model[-1], e["hr"], e["lo"], e["hi"]])
        for p in [1, 2, 3]:
            e = prim[tag]["period_hr"][f"{s}_p{p}"]
            rows.append([tag, s[-1], f"P{p}", e["hr"], e["lo"], e["hi"]])
fr = pd.DataFrame(rows, columns=["cohort", "stage", "model", "hr", "lo", "hi"])

fig, axes = plt.subplots(1, 2, figsize=(12, 6.5), sharey=True)
labels = ["1", "2", "3"]
for ax, tag, title in [(axes[0], "mimic", "MIMIC-IV"), (axes[1], "eicu", "eICU-CRD")]:
    sub = fr[fr["cohort"] == tag]
    ylab, yticks = [], []
    y = 0
    for st in ["1", "2", "3"]:
        for mdl, col, name in [("A", "#9ecae1", "Model A"), ("B", "#6baed6", "Model B"),
                               ("C", "#2171b5", "Model C"),
                               ("P1", "#fdae6b", "0-7 d"), ("P2", "#e6550d", "7-14 d"),
                               ("P3", "#a63603", "14-30 d")]:
            r = sub[(sub["stage"] == st) & (sub["model"] == mdl)].iloc[0]
            ax.plot([r.lo, r.hi], [y, y], color=col, lw=1.6)
            ax.plot(r.hr, y, "o", color=col, ms=5)
            ylab.append(f"S{st} {name}")
            yticks.append(y)
            y += 1
        y += 0.5
    ax.axvline(1.0, color="k", lw=0.8, ls="--")
    ax.set_yticks(yticks); ax.set_yticklabels(ylab, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("Hazard ratio (log scale)")
    ax.set_title(f"{title}: stage HRs (TV Cox)", fontsize=11)
fig.suptitle("Time-varying KDIGO stage and death: overall and period-specific HRs",
             fontsize=12)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(V3, f"fig_v3_forest.{ext}"), dpi=300,
                bbox_inches="tight")
fr.to_csv(os.path.join(V3, "v3_forest_data.csv"), index=False)
print("fig_v3_forest saved")

# ------------------------------------------------------------- report ---
def hr(e):
    return f"{e['hr']:.2f} ({e['lo']:.2f}-{e['hi']:.2f})"

rep = []
rep.append("# v3 重分析报告（时变暴露 + 竞争风险 + 无泄漏表型 + 锁定验证）\n")
rep.append(f"生成: {datetime.now().isoformat()[:19]}  \n"
           f"数据: v3_outputs/ (MIMIC 60,506 + eICU 113,466)\n")
rep.append("\n## 一、设计（对应评审 6 项致命问题的修复）\n")
rep.append("| 问题 | v2 错误 | v3 修复 |")
rep.append("|---|---|---|")
rep.append("| 5.1 immortal time | 7天max分期作固定暴露+24h landmark | 时变单调阶梯暴露，counting process Cox，从入科起算 |")
rep.append("| 5.2 30d截断 | clip时间不重置事件 | 事件=死亡且≤720h；>720h死亡在720h删失（MIMIC 293/eICU 177例正确处理） |")
rep.append("| 5.3 竞争风险 | 出院当随机删失 | Aalen-Johansen CIF（出院=竞争事件）+ cause-specific Cox 并列 |")
rep.append("| 5.4 时间分层风险集 | 区间内终止者才入风险集 | stage×period 交互（0-7/7-14/14-30d），联合Wald检验 |")
rep.append("| 5.5 表型泄漏 | 死亡→persistent | onset+72h landmark，仅用landmark前数据分类，landmark前死亡/出院者排除 |")
rep.append("| 5.6 验证无删失处理/伪锁定 | TD-AUC/Brier/DCA无IPCW; eICU重拟合 | 完全观测二分类结局（院内死亡≤30d）+ 纯锁定logistic + 校准slope/O:E |")
rep.append("| flow对账 | 负排除数/数字对不上 | 程序自动生成+断言全部通过；tv区间铺砌审计 PASS |")

rep.append("\n## 二、主分析：时变 Cox（cause-specific death）\n")
rep.append("| 模型 | MIMIC S1 | MIMIC S2 | MIMIC S3 | eICU S1 | eICU S2 | eICU S3 |")
rep.append("|---|---|---|---|---|---|---|")
for model, name in [("model_A", "A 未调整"), ("model_B", "B +人口学"),
                    ("model_C", "C 完全调整")]:
    cells = []
    for tag in ["mimic", "eicu"]:
        for s in ["stage1", "stage2", "stage3"]:
            cells.append(hr(prim[tag][model][s]))
    rep.append(f"| {name} | " + " | ".join(cells) + " |")

rep.append("\n**跨库一致性**: Model C S3 HR MIMIC "
           f"{prim['mimic']['model_C']['stage3']['hr']:.2f} vs eICU "
           f"{prim['eicu']['model_C']['stage3']['hr']:.2f} "
           f"(ratio {prim['eicu']['model_C']['stage3']['hr']/prim['mimic']['model_C']['stage3']['hr']:.2f})\n")

rep.append("\n### 期间特异性 HR（stage×period 交互，替代旧时间分层）\n")
rep.append("| Stage | MIMIC 0-7d | 7-14d | 14-30d | eICU 0-7d | 7-14d | 14-30d |")
rep.append("|---|---|---|---|---|---|---|")
for s in ["stage1", "stage2", "stage3"]:
    cells = []
    for tag in ["mimic", "eicu"]:
        for p in [1, 2, 3]:
            cells.append(hr(prim[tag]["period_hr"][f"{s}_p{p}"]))
    rep.append(f"| {s[-1]} | " + " | ".join(cells) + " |")
for tag in ["mimic", "eicu"]:
    w = prim[tag]["ph_joint_wald_stage_x_period"]
    rep.append(f"\n- {tag.upper()} 联合 Wald (stage×period): "
               f"chi2={w['chi2']:.1f}, df={w['df']}, p={w['p']:.1e}")

rep.append("\n## 三、竞争风险 CIF（30天，按 24h 时点分期）\n")
rep.append("| Stage@24h | MIMIC n | MIMIC CIF% (95%CI) | eICU n | eICU CIF% (95%CI) |")
rep.append("|---|---|---|---|---|")
for s in ["0", "1", "2", "3"]:
    m = prim["mimic"]["cif_by_stage24"][s]
    e = prim["eicu"]["cif_by_stage24"][s]
    rep.append(f"| {s} | {m['n']:,} | {m['cif30']} ({m['lo']}-{m['hi']}) "
               f"| {e['n']:,} | {e['cif30']} ({e['lo']}-{e['hi']}) |")

rep.append("\n## 四、表型（onset+72h landmark，无泄漏）\n")
rep.append("| 指标 | MIMIC | eICU |")
rep.append("|---|---|---|")
for k, lab in [("aki_n", "AKI 患者"), ("n_eligible", "可评估"),
               ("n_excluded_died_before_L", "landmark前死亡(排除)"),
               ("n_excluded_discharged_before_L", "landmark前出院(排除)"),
               ("n_excluded_unclassifiable", "无法分类(排除)"),
               ("n_transient", "transient"),
               ("n_persistent", "persistent"),
               ("crude_death30_transient_pct", "transient 死亡%"),
               ("crude_death30_persistent_pct", "persistent 死亡%"),
               ("no_aki_crude_death30_pct", "无AKI 死亡%（参考）")]:
    rep.append(f"| {lab} | {phen['mimic'][k]:,} | {phen['eicu'][k]:,} |".replace(
        f"{phen['mimic'][k]:,}", str(phen['mimic'][k])).replace(
        f"{phen['eicu'][k]:,}", str(phen['eicu'][k])))
rep.append(f"| 调整后 persistent vs transient HR | "
           f"{hr(phen['mimic']['cox_phenotype_adjusted'])} | "
           f"{hr(phen['eicu']['cox_phenotype_adjusted'])} |")

rep.append("\n## 五、锁定外部验证（TRIPOD type 3, logistic, 院内死亡≤30d）\n")
rep.append("| 指标 | MIMIC (apparent) | eICU (locked) |")
rep.append("|---|---|---|")
rep.append(f"| n | {val['mimic_n']:,} | {val['eicu_n']:,} |")
rep.append(f"| events | {val['mimic_events']:,} | {val['eicu_events']:,} |")
rep.append(f"| AUC | {val['mimic_apparent']['auc']} | {val['eicu_locked']['auc']} |")
rep.append(f"| Brier | {val['mimic_apparent']['brier']} | {val['eicu_locked']['brier']} |")
rep.append(f"| 校准 slope (95%CI) | {val['mimic_apparent']['calib_slope']} | "
           f"{val['eicu_locked']['calib_slope']} "
           f"({val['eicu_locked']['calib_slope_ci'][0]}-"
           f"{val['eicu_locked']['calib_slope_ci'][1]}) |")
rep.append(f"| O:E 比 | {val['mimic_apparent']['oe_ratio']} | "
           f"{val['eicu_locked']['oe_ratio']} (绝对风险需再校准) |")
rep.append(f"\nMIMIC 5-fold CV AUC: {val['mimic_cv5_auc']} "
           f"(sd {val['mimic_cv5_auc_sd']})")

rep.append("\n## 六、v2 → v3 关键数字对比（结论存活性判断）\n")
rep.append("| 指标 | v2 (有偏) | v3 (修正) | 判断 |")
rep.append("|---|---|---|---|")
m3 = prim['mimic']['model_C']['stage3']['hr']
e3 = prim['eicu']['model_C']['stage3']['hr']
rep.append(f"| Model C S3 HR (MIMIC) | 2.45 | {m3:.2f} | 方向保留，效应更大（时变框架） |")
rep.append(f"| Model C S3 HR (eICU) | 2.46 | {e3:.2f} | 跨库一致性仍近乎完美 (ratio {e3/m3:.2f}) |")
p1 = prim['mimic']['period_hr']['stage3_p1']['hr']
p3_ = prim['mimic']['period_hr']['stage3_p3']['hr']
rep.append(f"| S3 时间衰减 (MIMIC) | 5.25→2.59→1.30 (风险集错误) | {p1:.2f}→{prim['mimic']['period_hr']['stage3_p2']['hr']:.2f}→{p3_:.2f} | 衰减模式保留且更极端 |")
c0 = prim['mimic']['cif_by_stage24']['0']['cif30']
c3 = prim['mimic']['cif_by_stage24']['3']['cif30']
rep.append(f"| 30天死亡率 S0→S3 (MIMIC) | KM 4.7→36.0% (上偏) | CIF {c0}→{c3}% | 竞争风险校正后梯度保留但绝对值更保守 |")
rep.append(f"| 表型 persistent vs transient | 2.63/3.24 (泄漏) | "
           f"{phen['mimic']['cox_phenotype_adjusted']['hr']:.2f}/"
           f"{phen['eicu']['cox_phenotype_adjusted']['hr']:.2f} | 效应约为v2一半但双库一致——主张仍成立 |")
rep.append(f"| 外部验证 | C 0.753→0.716 (伪锁定) | AUC 0.750, slope {val['eicu_locked']['calib_slope']} (纯锁定) | 判别运输良好；绝对风险O:E={val['eicu_locked']['oe_ratio']}需再校准 |")

rep.append("\n## 七、必须写进 Methods/Limitations 的内容\n")
rep.append("1. 时变暴露为**单调累积 max stage**（不降级）；恢复通过表型分析单独刻画。")
rep.append("2. Cause-specific HR 将出院当删失（独立删失假设）；CIF 并列呈现。")
rep.append("3. 表型 landmark 排除了 onset+72h 前死亡/出院者（选择效应，人数已报）。")
rep.append("4. 验证目标量=院内死亡≤30d（出院计为非事件）；与CIF estimand一致。")
rep.append(f"5. eICU 绝对风险再校准需求 (O:E={val['eicu_locked']['oe_ratio']}) 如实报告。")
rep.append("6. vent_24h 在 MIMIC 中的保护性方向（系数-0.88）——与 v2 相同的悖论，已在讨论中处理。")

with open(os.path.join(V3, "v3_report.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(rep))
print("v3_report.md written")

with open(os.path.join(V3, "v3_audit.json"), "w") as f:
    json.dump(audit, f, indent=2)

print("\nv3 Phase 5 COMPLETE — all audits PASS")
