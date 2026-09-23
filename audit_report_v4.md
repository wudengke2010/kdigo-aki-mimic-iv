# manuscript_revised_v4.md 逐字逐句审核报告

**审核日期:** 2026-09-23 | **稿件版本:** v4 → v4.1（审核后修订）
**核对方法:** 全部数字逐一对照源文件（cohort_v2_flow.json / eicu_cohort_v2_flow.json / sofa_v2_summary.json / mimic·eicu_analysis_v2_meta.json / survival_mimic.json / survival_eicu.json / external_validation_results.json / _ph_mimic.pkl / 两个分析CSV），全部 Methods 声明逐句对照脚本源码（kdigo_cohort_v2.py / eicu_cohort_v2.py / build_mimic_analysis_v2.py / eicu_full_covariates.py / survival_*_v2.py / external_validation_v2.py）。

---

## 一、发现并已修复的问题

### P0-1｜MIMIC vent_24b 无首24h时间过滤（数据级错误，已修复+全链路重跑）
- **位置:** build_mimic_analysis_v2.py:37-48（原SQL仅有itemid过滤，无starttime窗口）
- **影响:** 稿件全部"first 24 h ventilation"表述失实；Table 1、Model C、锁定模型、跨库一致性、Discussion vent段落全部受污染。与导致v1被拒的F4（跨时间窗污染）同类。
- **修复:** 加 `[intime, intime+24h]` 过滤 → vent 19,720 (32.6%) vs 原 22,859 (全stay)；重跑 survival_mimic_v2 + PH检验 + external_validation_v2 + 森林图。
- **新数字:** Model C vent HR 0.597 (0.555–0.642)；S3 HR 2.45 (2.26–2.66)；C-index 0.753；time-stratified S3 5.25/2.59/1.30；eICU锁定 C 0.716 (Δ −0.037)；kdigo_3 跨库 ratio **1.004**（更完美）；vent ratio 2.46（反转更显著，讨论已重写并补充去vent敏感性 S3 2.61）。

### P0-2｜基线肌酐层级描述与实现不符（Methods失实，已改写）
- **实际实现（kdigo_cohort_v2.py:221-239）:** 四层级、outpatient优先：(1) outpatient 7–365d **median** → (2) any-setting 7–365d median → (3) 0–7d pre-admission median → (4) 住院首24h **first**。eICU: pre-ICU median（窗口=max(住院入科, −7d)）→ ICU首24h first。
- **原稿错误:** 写成三层级、0-7d优先、用"lowest"（实际median）。
- **修复:** Methods按实现改写为四层级。

### P0-3｜MIMIC 院前基线比例错误（Results数字错误，已改）
- 原稿"42.2%"= prior_0_7d 单一来源占比；**院前基线合计 = (25,615+21,330+1,594)/60,506 = 80.2%**。eICU 39.4%。不对称方向表述已纠正。

### P0-4｜流程排除数与flow.json不可调和（已修正）
- MIMIC "(580 excluded)" → 实际 **898**（425无基线 + 473有基线但无窗口SCr）。flow.json的"excluded_no_window_scr=155"是脚本计数口径错误（按全队列而非基线可得子集），不影响CSV数据。
- eICU "baseline creatinine requirements" → 实际 **10,372**（8,341 + 2,031）。flow.json的"-1,806"负数为同一计数bug。
- 稿件已按真实排除数改写；Figure 1流程图重绘时须用新数字。

### P0-5｜351例异常死亡时间处理描述失实（已改写）
- 实际（survival_mimic_v2.py:50-55）: decedent事件时间**封顶于 dischtime+24h**（仍计为事件），非"censored at discharge"；另有 **5例** 无deathtime死亡者被landmark排除。Methods已准确改写。

### P1-1｜eICU ESRD双源数字有重叠（已澄清）
- ICD 2,801 + pastHistory 6,219 = 9,020 ≠ 并集 6,523（重叠2,497）；pastHistory实际含三条路径（含"not currently dialyzed"）。已加"with overlap"及路径说明。

### P1-2｜表名/时间基准/统计细节（已修正）
- `procedureevents_mv` → v3.1 实为 `procedureevents`；RRT条目（CRRT/CVVHD/CVVHDF/HD/SCUF；eICU排除catheter/radiology）已写明。
- TD-AUC "7/14/28 days" → 明确为 **landmark后** 7/14/28天；Brier "30d" → 实为 landmark后28天（随访终点），已改。
- HL → 明确为"linear-predictor十分位的Hosmer-Lemeshow-type统计量"（脚本注释自述rough HL）。
- MIMIC time-stratified 用 λ=0.01 轻ridge（Table 3脚注已披露）；eICU ts协变量集已在脚注写明。
- MIMIC合并症 → 由"coding algorithms [13,14]"改为如实列出Elixhauser-style ICD前缀实现。
- MIMIC ESRD/CKD → 补充"限定index admission或更早/严格早于index admission"（与代码一致，已验证 kdigo_cohort_v2.py:143-160）。
- eICU事件时间 hospitaldischargeoffset/hospitaldischargestatus（已验证）。

### P1-3｜两处 Model C 系数微差（已消除）
- survival_mimic_v2（封顶351例、penalizer=0）与 external_validation_v2（未封顶、penalizer=0.001）拟合的"同一"Model C 系数不一致（2.453 vs 2.455）。**已将external_validation_v2预处理改为与survival脚本完全一致并重跑** → 两处现在完全一致（S3 HR 2.453）。

---

## 二、审核确认无误的项目（抽样列举）

| 项目 | 验证结果 |
|---|---|
| 队列数 60,506 / 113,466；landmark 58,554 / 107,758（死亡5,318/8,616） | ✅ CSV+JSON一致 |
| KDIGO分布 43,011/11,882/2,249/3,364；95,203/11,999/1,874/4,390 | ✅ |
| AKI 28.9% vs 16.1% | ✅ |
| 30d死亡率 MIMIC 4.7/11.8/24.2/36.0%；eICU 5.0/18.0/31.7/30.2% | ✅ survival json |
| 粗死亡率、SOFA by stage（4.23/6.82/7.71/9.05；2.87/4.66/5.62/5.68） | ✅ CSV重算 |
| KM log-rank χ² 1,826.3 / 2,675.1；S2vsS3 Bonferroni p<0.001 / 0.040 | ✅ |
| Model A/B 全部HR与CI（两队列） | ✅ |
| eICU同协变量refit 1.86/2.67/2.46；APACHE敏感性 1.72 (n=91,530, C 0.805) | ✅ |
| PH: MIMIC全过（KDIGO p≥0.29, copd p=0.046）；eICU仅age (ρ−0.038, p=8.4e-4) | ✅ pickle+json |
| 外部验证全套（C 0.716/0.753；AUC；Brier 0.0710/0.0802；slope 0.957；HL 867；tertile 19.9/32.1/52.4；DCA 0.034vs0.032；13个HR ratio） | ✅ 重跑后JSON |
| E-value: HR 2.45 → 2.45+√(2.45×1.45)=4.33≈4.3 | ✅ 重算 |
| 去vent敏感性 S3 2.61 (2.40–2.83), C 0.749 | ✅ 新增 |
| 173,972总数；351/5例处理；KDIGO分期逻辑（48h滚动窗/7d比值/≥4.0需急性变化/RRT限7d） | ✅ 代码逐行 |
| 34篇引用、基金、伦理、ORCID | ✅ 与v3一致 |

## 三、遗留待办
1. **Figure 1 流程图**须按新排除数重绘（MIMIC 898；eICU 10,372）
2. Table S1（基线来源分布）正式生成
3. 选刊 → 格式压缩 → DOCX → GitHub推送 v2.1 代码（vent修复+外部验证一致性修复）
