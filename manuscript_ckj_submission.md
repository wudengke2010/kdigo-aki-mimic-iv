# Association of KDIGO Creatinine-Based AKI Staging with ICU Mortality: A Dual-Cohort Retrospective Study Using MIMIC-IV v3.1 and eICU-CRD

---

## Authors

**Jiqiang Liu**<sup>1,2</sup>, Dengke Wu<sup>1,2*</sup>

<sup>1</sup> Department of Emergency Medicine, The Second Xiangya Hospital of Central South University, Changsha 410011, Hunan, China

<sup>2</sup> Emergency Medicine and Difficult Diseases Institute, The Second Xiangya Hospital of Central South University, Changsha 410011, Hunan, China

**\*Corresponding author:** Dengke Wu, Department of Emergency Medicine, Second Xiangya Hospital, Central South University, Changsha, Hunan, China; Emergency Medicine and Difficult Diseases Institute, The Second Xiangya Hospital of Central South University, Changsha 410011, Hunan, China. Electronic address: wudk2010@csu.edu.cn.

**Running head:** KDIGO AKI staging and ICU mortality

---

## Abstract

**Background:** KDIGO creatinine-based AKI staging is widely used, but illness severity confounding, Stage 2–3 mortality convergence, and the "AKI paradox" (CKD patients tolerating severe AKI better) remain incompletely characterised. No study has simultaneously quantified these with reclassification metrics, time-stratified analysis, external validation, and a prediction tool.

**Methods:** This STROBE-compliant retrospective cohort study used two independent databases: MIMIC-IV v3.1 (84,167 ICU stays; derivation) and eICU-CRD v2.0 (166,373 ICU stays; external validation; total N = 250,540). AKI was staged by KDIGO creatinine criteria. Three incremental logistic models were fitted with severity adjustment (SOFA in MIMIC; APACHE IVa in eICU). Cox regression with time-stratified analysis, Kaplan-Meier with Bonferroni correction, restricted cubic splines, 24-hour landmark analysis, time-dependent ROC (7/14/28-day), and a clinical nomogram with web calculator were developed.

**Results:** AKI occurred in 29.9% (MIMIC) and 31.8% (eICU) of ICU stays. Mortality was 6.0%/15.1%/26.7%/29.0% (MIMIC Stages 0–3) and 4.6%/12.5%/18.0%/19.2% (eICU). Severity adjustment attenuated KDIGO ORs by 32–59% (NRI = 0.625; AUC: 0.778→0.828 MIMIC; 0.672→0.818 eICU). The AKI paradox was replicated in both cohorts (CKD × Stage 3 interaction OR = 0.476, MIMIC; 0.726, eICU; both p < 0.001). Time-dependent AUCs were 0.643/0.639/0.650 at 7/14/28 days. A nomogram (KDIGO, APACHE, ventilation, age, sex) achieved AUC = 0.818 (Brier = 0.062).

**Conclusions:** KDIGO staging is independently associated with ICU mortality across two independent cohorts (N = 250,540), but illness severity accounts for a substantial portion of this association. The Stage 3-specific AKI paradox is reproducible across databases, supporting its biological rather than artefactual nature.

---

## Key Learning Points

- **What was known:** KDIGO creatinine-based AKI staging is associated with ICU mortality, but the degree of illness severity confounding, the Stage 2–3 mortality convergence, and the "AKI paradox" (better tolerance of severe AKI in CKD patients) have not been simultaneously quantified with reclassification metrics or validated across independent databases.

- **This study adds:** In 250,540 ICU stays from two databases, severity accounted for 32–59% of the KDIGO–mortality association (NRI = 0.625), the AKI paradox was replicated cross-database (interaction OR 0.476 and 0.726), and a nomogram with web calculator (AUC = 0.818) was deployed for bedside use.

- **Potential impact:** These findings support continued use of KDIGO staging while advocating for CKD-aware, severity-adjusted, and temporally nuanced risk stratification. The web calculator enables real-time mortality risk estimation using routinely available ICU variables.

---

## Keywords

acute kidney injury; intensive care unit; KDIGO; mortality; nomogram

---

## Introduction

Acute kidney injury (AKI) affects 30–60% of critically ill patients [6], with the landmark AKI-EPI multinational cohort reporting an incidence of 57.7% and a mortality gradient from 16.1% (no AKI) to 50.0% (Stage 3) [1,2]. The KDIGO 2012 criteria standardised AKI classification using serum creatinine (SCr), urine output, and renal replacement therapy (RRT) [3], and several studies have validated its prognostic value [4–7]. However, critical knowledge gaps persist.

**First**, the degree to which the KDIGO–mortality association is confounded by overall illness severity has not been quantified with formal reclassification metrics (NRI/IDI) in large cohorts. Girling et al. [24] applied Bradford Hill criteria to AKI-associated mortality and identified severity confounding as the principal obstacle to causal interpretation.

**Second**, the incremental prognostic distinction between Stage 2 and Stage 3 remains debated, as Stage 3 is enriched with chronic kidney disease (CKD) patients exhibiting different mortality trajectories [8,9]. Whether the observed Stage 2–3 mortality convergence represents a biological ceiling or a staging artefact has not been resolved through time-stratified analysis.

**Third**, the "AKI paradox"—where CKD patients appear to tolerate severe AKI better than non-CKD patients [10,11]—has not been systematically decomposed in a cohort large enough for formal interaction testing, nor validated across independent databases.

**Fourth**, no study has translated KDIGO-based risk stratification into a practical, deployable clinical prediction tool validated in an external cohort. While machine-learning approaches have advanced AKI prediction in specific subgroups [26,27], none has integrated KDIGO staging into an interpretable nomogram with cross-database validation.

MIMIC-IV v3.1 (84,167 ICU stays, 2008–2022) [12] and eICU-CRD v2.0 (200,859 ICU stays from 208 hospitals) [23] provide independent cohorts for derivation and validation. Lin et al. [25] established the precedent for this dual-database design in AKI subphenotyping. This study aimed to: (1) quantify illness severity confounding using incremental regression with NRI/IDI; (2) characterise survival trajectories with time-stratified Cox analysis; (3) decompose the AKI paradox by RRT indication; (4) externally validate all findings in eICU-CRD; (5) assess time-dependent predictive performance; and (6) develop a deployable nomogram with web calculator.

---

## Methods

### Study Design and Data Source

This STROBE-compliant retrospective dual-cohort study used MIMIC-IV v3.1 (derivation; single-centre, BIDMC, 2008–2022) [12] and eICU-CRD v2.0 (external validation; 208 US hospitals, 2014–2015) [23].

### Study Population

Adult patients (≥18 years) with at least one ICU stay and at least one SCr measurement were included. For patients with multiple ICU stays, only the first was retained. After exclusions (Figure 1), 84,167 MIMIC-IV and 166,373 eICU-CRD stays were analysed.

### Exposure: KDIGO Creatinine-Based AKI Staging

AKI was staged by KDIGO 2012 creatinine criteria: Stage 1 (≥1.5× baseline or ≥0.3 mg/dL increase), Stage 2 (≥2.0× baseline), Stage 3 (≥3.0× baseline, SCr ≥4.0 mg/dL, or RRT). Baseline SCr was the first hospital measurement; peak SCr was the highest ICU value. This implementation was verified against the original KDIGO criteria [20]. Urine output was not used due to incomplete recording.

### Outcomes and Covariates

The primary outcome was in-hospital mortality. Covariates included age, sex, Elixhauser comorbidities [13,14], SOFA score (MIMIC) or APACHE IVa (eICU), mechanical ventilation, vasopressor use, and ICD-defined CKD. Laboratory-based CKD (CKD-EPI eGFR <60) was assessed in sensitivity analysis. Full covariate details are in Supplementary Material.

### Statistical Analysis

Three incremental logistic models were fitted: Model A (KDIGO + age + sex + comorbidities), Model B (Model A + SOFA/APACHE + ventilation + vasopressors + CKD), and Model C (continuous creatinine ratio replacing KDIGO stages). Severity confounding was quantified by OR attenuation (Model A→B) and NRI/IDI reclassification [21]. VIF diagnostics assessed collinearity [17].

Kaplan-Meier survival curves were compared with Bonferroni-corrected log-rank tests. Cox regression with Schoenfeld residual testing was performed; time-stratified Cox models (0–7, 7–14, 14–30 days) addressed proportional hazards (PH) violations, following the approach of Singh et al. [28]. Restricted mean survival time (RMST) was computed as a PH-independent metric.

Restricted cubic splines (4-knot) modelled the nonlinear creatinine ratio–mortality relationship. Model discrimination (AUC), calibration (Hosmer-Lemeshow), and decision curve analysis (DCA) were assessed. Sepsis and CKD × KDIGO interaction terms were tested.

Five sensitivity analyses were performed: (1) MDRD-based baseline SCr; (2) laboratory-based CKD; (3) complete-case SOFA; (4) 24-hour landmark analysis; (5) temporal trend analysis (2008–2022).

The eICU-CRD analysis replicated the primary pipeline with APACHE IVa replacing SOFA. Time-dependent ROC curves (7/14/28-day, cumulative sensitivity/dynamic specificity framework) and a clinical nomogram with web calculator were developed from the eICU Model B. E-values assessed unmeasured confounding [22]. A directed acyclic graph (DAG) defined the causal structure (Figure S5). Analyses used Python 3.13 (pandas 2.3.3, statsmodels 0.14.6, lifelines 0.30.3).

---

## Results

### Derivation Cohort (MIMIC-IV)

From 94,458 ICU stays, 84,167 were included (Figure 1). AKI (Stage 1–3) occurred in 25,134 (29.9%): Stage 1 12,817 (15.2%), Stage 2 2,729 (3.2%), Stage 3 9,588 (11.4%). Among Stage 3, 3,971 (41.4%) received RRT.

Baseline characteristics are presented in **Table 1**. The mean age was 63.0 ± 16.8 years; 55.8% were male. Illness severity increased markedly: SOFA 4.3→9.9 (Stages 0→3), mechanical ventilation 43.0%→63.4%, vasopressors 33.8%→63.6% (all p < 0.001). CKD (ICD) increased from 19.4% (Stage 0) to 70.7% (Stage 3). Notably, Stage 2 patients had the lowest baseline SCr (0.90 ± 0.36 mg/dL) and highest eGFR (84.3 mL/min/1.73 m²), reflecting the fold-change nature of KDIGO staging: patients with low baseline SCr require smaller absolute increases to reach the ≥2.0-fold threshold.

### Primary Outcome: In-Hospital Mortality

Mortality was 6.0% (Stage 0), 15.1% (Stage 1), 26.7% (Stage 2), 29.0% (Stage 3); overall 10.7%. The Stage 2→3 increment was modest (26.7% vs. 29.0%).

### Multivariable Logistic Regression

**Model A** ORs: Stage 1 2.762, Stage 2 5.261, Stage 3 7.029 (all p < 0.001). **Model B** (severity-adjusted) ORs: Stage 1 1.885 (95% CI 1.765–2.013), Stage 2 2.978 (2.688–3.300), Stage 3 2.881 (2.664–3.115); all p < 0.001 (**Table 2**). Severity adjustment attenuated ORs by 31.8% (Stage 1), 43.4% (Stage 2), 59.0% (Stage 3), with NRI = 0.625, IDI = 0.058, and AUC improvement from 0.778 to 0.828. In Model B, SOFA was the strongest predictor (aOR 1.333 per point). E-values were 2.98–4.62, indicating robustness to unmeasured confounding. VIF values were all <5 (SOFA 4.73, ventilation 2.68). **Model C** (continuous creatinine ratio) aOR 2.135 confirmed the independent prognostic value of continuous creatinine change (AUC 0.825).

### Survival Analysis and Time-Stratified Cox

Kaplan-Meier analysis confirmed significant stage separation (log-rank p < 0.001), except Stage 2 vs. 3 (Bonferroni p = 1.000). RMST converged at 23.5 days for both Stage 2 and 3.

Schoenfeld testing identified PH violations for Stage 3 (p < 0.001). Time-stratified Cox revealed Stage 3 HR progressively increased: 1.41 (0–7d) → 1.61 (7–14d) → 1.86 (14–30d), while Stage 2 effects were front-loaded (HR 1.54, 0–7d). This temporal divergence explains the Stage 2–3 convergence: Stage 3's mortality signal emerges later and persists longer.

### AKI Paradox

CKD-negative Stage 3 mortality was 49.4% versus CKD-positive 20.5% (2.4-fold difference). The formal interaction was Stage 3-specific: CKD × Stage 3 OR = 0.476 (p < 0.001), with no interaction at Stages 1–2. RRT-stratified decomposition showed the paradox persisted in the creatinine-only subgroup (42.3% vs. 18.5%), indicating it is not solely attributable to differential RRT initiation. Laboratory-based CKD (eGFR <60) confirmed the paradox (44.4% vs. 25.8%), while ICD-CKD showed a protective association (aOR 0.620) and lab-CKD showed increased mortality (aOR 1.385), suggesting these definitions capture different clinical phenotypes [30,32].

### Sensitivity Analyses

The 24-hour landmark analysis (excluding 1,658 early deaths, 18.5% of all mortality) yielded 12–29% higher ORs (Stage 1: 2.127, Stage 2: 3.644, Stage 3: 3.341), indicating the primary estimates are conservative—immortal time bias attenuated the full-cohort ORs toward the null. MDRD-based staging showed 68.9% agreement with primary staging, with MDRD Stage 3 aOR higher (4.750 vs. 2.881). Complete-case SOFA analysis (N = 17,285) showed 11–23% lower ORs, confirming SOFA imputation overestimates the independent KDIGO effect. Temporal trends showed AKI incidence rising from 28.4% (2008–2010) to 33.1% (2020–2022).

### External Validation (eICU-CRD)

The eICU-CRD cohort (N = 166,373) showed consistent results: AKI prevalence 31.8%, mortality 4.6%/12.5%/18.0%/19.2% (Stages 0–3), Model B ORs 1.849/2.811/2.581 (Stages 1–3), and AUC improvement from 0.672 to 0.818 (Model A→B). The severity attenuation closely paralleled MIMIC: Stage 1 by 38.3% (MIMIC: 31.8%), Stage 2 by 38.6% (43.4%), Stage 3 by 48.0% (59.0%). The AKI paradox was replicated: CKD × Stage 3 interaction OR = 0.726 (95% CI 0.603–0.872, p < 0.001), again Stage 3-specific with no interaction at Stages 1–2. **Table 3** presents the cross-cohort comparison.

### Time-Dependent ROC

The full Model B maintained stable discrimination: 7-day AUC 0.643, 14-day 0.639, 28-day 0.650. KDIGO stage alone showed declining performance: 0.484 → 0.452 → 0.404, confirming that static staging cannot capture evolving risk.

### Nomogram and Web Calculator

A five-variable nomogram (KDIGO stage, APACHE score, mechanical ventilation, age, sex) achieved AUC = 0.818 (Brier = 0.062) with good calibration. APACHE score contributed the most points (0–100), followed by KDIGO Stage 3 (0–48), ventilation (0–30), age (0–18), and sex (0–1). A standalone web calculator enables real-time bedside risk estimation (Supplementary Material).

---

## Discussion

In this dual-cohort study of 250,540 ICU stays, we comprehensively evaluated KDIGO creatinine-based AKI staging. Six key findings emerge.

Our mortality rates (MIMIC: 6.0%/15.1%/26.7%/29.0%; eICU: 4.6%/12.5%/18.0%/19.2% across Stages 0–3) are consistent with prior cohorts, though absolute rates vary with case mix and AKI definition. The AKI-EPI study reported 16.1%–50.0% using both creatinine and urine output criteria [2], while creatinine-only studies generally report lower mortality [19]. Hoste et al. [33] reported RIFLE Stage 1–3 mortality of 8.8%, 11.4%, and 26.3% in 5,383 ICU patients, and Nisula et al. [34] reported KDIGO Stage 1–3 mortality of 29.3%, 34.1%, and 39.0% in the FINNAKI study of 2,901 patients. The lower eICU mortality likely reflects its multi-centre community hospital case mix.

**First**, illness severity accounted for 32–59% of the KDIGO–mortality association (NRI = 0.625), providing direct empirical evidence for the causal-versus-correlational question raised by Girling et al. [24]. The landmark analysis confirmed that primary ORs were conservative (12–29% lower than landmark estimates).

**Second**, RCS analysis confirmed a significant nonlinear creatinine ratio–mortality relationship (χ² = 161.29, p < 0.001), with the OR peaking at ratio ~2.6 and plateauing at higher values, consistent with the Stage 2–3 convergence and the RMST convergence at 23.5 days. The plateau at extreme ratios supports a biological ceiling effect, and may reflect CKD enrichment at high ratios. Long et al. [31] showed that KDIGO's absolute (≥26.5 µmol/L) and relative (≥1.5× baseline) Stage 1 criteria identify different populations with different long-term outcomes stratified by CKD status, supporting that criterion type—not just stage number—drives prognosis.

**Third**, time-stratified Cox analysis revealed that Stage 3's hazard progressively increased over follow-up (HR 1.41→1.86), explaining the PH violation and convergence with Stage 2. The time-dependent ROC confirmed this: KDIGO-only AUC declined over time (0.484→0.404), while the full model maintained performance (0.643→0.650). Lin et al. [25] similarly demonstrated that trajectory-based AKI subphenotyping captures within-stage heterogeneity that static KDIGO categories cannot resolve.

**Fourth**, the AKI paradox was formally confirmed as Stage 3-specific and **replicated in the independent eICU-CRD cohort** (interaction OR 0.726, p < 0.001, vs. MIMIC 0.476). This cross-database replication—across different hospital systems, severity scores, and data platforms—provides the strongest evidence that the paradox is biological rather than artefactual. Gleeson et al. [29] demonstrated that muscle mass superseded creatinine as a survival predictor in ICU AKI, and that higher creatinine at RRT initiation paradoxically predicted better survival. Creatinine is a muscle-derived marker; elevated creatinine in non-CKD patients with preserved muscle mass may reflect less catastrophic injury than the same level in sarcopenic patients who reach Stage 3 at deceptively low values. This mechanism could explain both our AKI paradox and the Stage 2 low-baseline-creatinine anomaly, and is consistent with evidence that KDIGO criteria may be inadequate for CKD patients [18,20].

**Fifth**, sepsis significantly modified the KDIGO–mortality association (ORs 1.3–1.8-fold higher in sepsis), suggesting that AKI in sepsis carries a distinct pathophysiological signature [5,15].

**Sixth**, the nomogram (AUC = 0.818, Brier = 0.062) provides practical clinical translation. While CKJ studies have deployed ensemble ML models with SHAP interpretability for sepsis-associated AKI [26] and pancreatitis-associated AKI [27], our logistic nomogram prioritises bedside deployability—transportable coefficients, no computational infrastructure, and real-time web calculation—while achieving comparable discrimination (Luo et al. [26] reported AUC 0.732–0.778 across external cohorts). The choice of five routinely available variables (KDIGO, APACHE, ventilation, age, sex) maximises generalisability across ICU settings.

### What Distinguishes This Study

This study addresses four literature gaps through: (1) dual-cohort validation (total N = 250,540); (2) comprehensive confounding quantification (NRI/IDI, DCA, E-values); (3) temporal dynamics analysis (time-stratified Cox, time-dependent ROC); and (4) clinical translation (nomogram with web calculator).

### Strengths and Limitations

Strengths include the large combined sample, external validation, three-model comparison with reclassification metrics, Bonferroni correction, landmark analysis, time-stratified Cox, time-dependent ROC, formal CKD × KDIGO interaction tests, E-values, and a deployable nomogram.

Limitations include: (1) creatinine-only staging without urine output [16,25,28]; (2) single-centre derivation cohort mitigated by multi-centre validation; (3) ICD-based comorbidity coding bias [30,32]; (4) SOFA imputation bias (11–23% OR overestimation by complete-case analysis); (5) APACHE IVa and SOFA are not directly interchangeable; (6) serum creatinine is confounded by muscle mass and sarcopenia [29]; (7) observational design precludes causal inference. The DAG (Figure S5) formally identified colliders, and E-values ≥2.77 suggest robustness to unmeasured confounding.

### Clinical Implications

AKI staging should be interpreted in the context of illness severity. Non-CKD patients developing Stage 3 AKI represent an extremely high-risk phenotype (49.4% mortality). The nomogram and web calculator enable real-time risk estimation using five routinely available variables.

### Conclusion

KDIGO creatinine-based AKI staging demonstrated a robust, independent association with ICU mortality across 250,540 patients from two independent databases. Illness severity accounted for a substantial portion of this association. The AKI paradox—validated as a Stage 3-specific, cross-database phenomenon—likely reflects muscle mass–creatinine confounding. A clinically deployable nomogram and web calculator provide practical risk stratification, supporting continued KDIGO use while advocating for CKD-aware, severity-adjusted refinement.

---

## Acknowledgements

The authors thank the PhysioNet team and the Laboratory for Computational Physiology at MIT for creating and maintaining the MIMIC-IV and eICU-CRD databases.

---

## References

1. Susantitaphong P, Cruz DN, Cerda J, et al. World incidence of AKI: a meta-analysis. Clin J Am Soc Nephrol. 2013;8(9):1482–1493. doi:10.2215/CJN.00710113
2. Hoste EAJ, Bagshaw SM, Bellomo R, et al. Epidemiology of acute kidney injury in critically ill patients: the multinational AKI-EPI study. Intensive Care Med. 2015;41(8):1411–1423. doi:10.1007/s00134-015-3934-7
3. Kidney Disease: Improving Global Outcomes (KDIGO) Acute Kidney Injury Work Group. KDIGO Clinical Practice Guideline for Acute Kidney Injury. Kidney Int Suppl. 2012;2(1):1–138. doi:10.1038/kisup.2012.1
4. Lassnigg A, Schmidlin D, Mouhieddine M, et al. Minimal changes of serum creatinine predict prognosis in patients after cardiothoracic surgery. J Am Soc Nephrol. 2004;15(6):1597–1605. doi:10.1097/01.ASN.0000130340.93930.DD
5. Bagshaw SM, Uchino S, Bellomo R, et al. Septic acute kidney injury in critically ill patients. Clin J Am Soc Nephrol. 2007;2(3):431–439. doi:10.2215/CJN.03681106
6. Uchino S, Kellum JA, Bellomo R, et al. Acute renal failure in critically ill patients. JAMA. 2005;294(7):813–818. doi:10.1001/jama.294.7.813
7. Luo X, Jiang L, Du B, et al. A comparison of different diagnostic criteria of acute kidney injury in critically ill patients. Crit Care. 2014;18(4):R144. doi:10.1186/cc13977
8. Thomas ME, Blaine C, Dawnay A, et al. The definition of acute kidney injury and its use in practice. Kidney Int. 2015;87(1):62–73. doi:10.1038/ki.2014.328
9. Pannu N, James M, Hemmelgarn B, Klarenbach S. Association between AKI, recovery of renal function, and long-term outcomes. Clin J Am Soc Nephrol. 2013;8(2):194–202. doi:10.2215/CJN.06480612
10. James MT, Ghali WA, Knudtson ML, et al. Associations between AKI and cardiovascular and renal outcomes after coronary angiography. Circulation. 2011;123(4):409–416. doi:10.1161/CIRCULATIONAHA.110.970160
11. Chertow GM, Burdick E, Honour M, et al. Acute kidney injury, mortality, length of stay, and costs. J Am Soc Nephrol. 2005;16(11):3365–3370. doi:10.1681/ASN.2004090740
12. Johnson AEW, Bulgarelli L, Shen L, et al. MIMIC-IV, a freely accessible electronic health record dataset. Sci Data. 2023;10(1):1. doi:10.1038/s41597-022-01899-x
13. Elixhauser A, Steiner C, Harris DR, Coffey RM. Comorbidity measures for use with administrative data. Med Care. 1998;36(1):8–27. doi:10.1097/00005650-199801000-00004
14. Quan H, Sundararajan V, Halfon P, et al. Coding algorithms for defining comorbidities in ICD-9-CM and ICD-10. Med Care. 2005;43(11):1130–1139. doi:10.1097/01.mlr.0000182534.19832.83
15. Prowle JR, Echeverri JE, Ligabo EV, et al. Fluid balance and acute kidney injury. Nat Rev Nephrol. 2010;6(2):107–115. doi:10.1038/nrneph.2009.213
16. Kellum JA, Sileanu FE, Murugan R, et al. Classifying AKI by urine output versus serum creatinine level. J Am Soc Nephrol. 2015;26(9):2231–2238. doi:10.1681/ASN.2014070724
17. Harrell FE Jr, Lee KL, Mark DB. Multivariable prognostic models. Stat Med. 1996;15(4):361–387. doi:10.1002/(SICI)1097-0258(19960229)15:4<361::AID-SIM168>3.0.CO;2-4
18. Hong D, Ren Q, Zhang J, et al. A new criteria for acute on preexisting kidney dysfunction in critically ill patients. Ren Fail. 2023;45(1):2173498. doi:10.1080/0886022X.2023.2173498
19. Birkelo BC, Pannu N, Siew ED. Overview of diagnostic criteria and epidemiology of acute kidney injury. Clin J Am Soc Nephrol. 2022;17(5):717–735. doi:10.2215/CJN.14181021
20. Yasrebi-de Kom IAR, Dongelmans DA, Abu-Hanna A, et al. Incorrect application of the KDIGO acute kidney injury staging criteria. Clin Kidney J. 2022;15(5):937–941. doi:10.1093/ckj/sfab256
21. Dong GY, Qin JP, An Y, et al. Utilizing reclassification to explore characteristics and prognosis of KDIGO SCr AKI subgroups. Ren Fail. 2021;43(1):1569–1576. doi:10.1080/0886022X.2021.1997761
22. VanderWeele TJ, Ding P. Sensitivity Analysis in Observational Research: Introducing the E-Value. Ann Intern Med. 2017;167(4):268–274. doi:10.7326/M16-2607
23. Pollard TJ, Johnson AEW, Raffa JD, et al. The eICU Collaborative Research Database. Sci Data. 2018;5:180178. doi:10.1038/sdata.2018.178
24. Girling BJ, Channon SW, Haines RW, Prowle JR. Acute kidney injury and adverse outcomes of critical illness: correlation or causation? Clin Kidney J. 2020;13(2):133–141. doi:10.1093/ckj/sfz158
25. Lin J, Liu L, Zhu S, et al. Machine learning-derived multivariate renal function trajectories in acute kidney injury in critically ill patients. Clin Kidney J. 2025;18(6):sfaf142. doi:10.1093/ckj/sfaf142
26. Luo S, Lai J, Mo L, et al. Prediction of hospital mortality in sepsis-associated AKI using machine learning. Clin Kidney J. 2026;19(1):sfaf372. doi:10.1093/ckj/sfaf372
27. Liu Y, Zhu X, Xue J, et al. Machine learning models for mortality prediction in acute pancreatitis-associated AKI. Clin Kidney J. 2024;17(10):sfae284. doi:10.1093/ckj/sfae284
28. Singh S, Andonovic M, Traynor JP, et al. Short- and long-term outcomes in oliguric and non-oliguric AKI in intensive care. Clin Kidney J. 2025;18(6):sfaf170. doi:10.1093/ckj/sfaf170
29. Gleeson PJ, Crippa IA, Sannier A, et al. Critically ill patients with acute kidney injury: clinical determinants and post-mortem histology. Clin Kidney J. 2023;16(10):1664–1673. doi:10.1093/ckj/sfad113
30. Høyer S, Heide-Jørgensen U, Jensen SK, et al. Sensitivity and positive predictive value of diagnosis codes for acute kidney injury in Denmark. Clin Kidney J. 2026;19(3):sfag019. doi:10.1093/ckj/sfag019
31. Long TE, Helgason D, Helgadottir S, et al. Mild Stage 1 post-operative AKI: association with CKD and long-term survival. Clin Kidney J. 2021;14(1):237–244. doi:10.1093/ckj/sfz197
32. Esposito P, Cappadona F, Marengo M, et al. Recognition patterns of acute kidney injury in hospitalized patients. Clin Kidney J. 2024;17(8):sfae231. doi:10.1093/ckj/sfae231
33. Hoste EAJ, Clermont G, Kersten A, et al. RIFLE criteria for acute kidney injury are associated with hospital mortality. Crit Care. 2006;10(3):R73. doi:10.1186/cc4915
34. Nisula S, Kaukonen KM, Vaara ST, et al. Incidence, risk factors and 90-day mortality of AKI in Finnish ICUs: the FINNAKI study. Intensive Care Med. 2013;39(3):420–428. doi:10.1007/s00134-012-2796-5

---

## Tables

### Table 1. Baseline Characteristics of the Derivation Cohort Stratified by KDIGO AKI Stage (MIMIC-IV, N = 84,167)

| Variable | No AKI (n = 59,033) | Stage 1 (n = 12,817) | Stage 2 (n = 2,729) | Stage 3 (n = 9,588) | Overall (N = 84,167) |
|:---|:---:|:---:|:---:|:---:|:---:|
| **Demographics** | | | | | |
| Age (years) | 62.3 ± 17.3 | 67.1 ± 14.9 | 64.7 ± 15.3 | 61.5 ± 15.3 | 63.0 ± 16.8 |
| Male, n (%) | 31,924 (54.1%) | 7,747 (60.4%) | 1,424 (52.2%) | 5,879 (61.3%) | 46,974 (55.8%) |
| Race — White, n (%) | 40,031 (67.8%) | 8,677 (67.7%) | 1,790 (65.6%) | 5,362 (55.9%) | 55,860 (66.4%) |
| Race — Black, n (%) | 5,697 (9.7%) | 1,358 (10.6%) | 270 (9.9%) | 1,942 (20.3%) | 9,267 (11.0%) |
| **Severity Scores** | | | | | |
| SOFA score | 4.3 ± 3.9 | 7.0 ± 4.0 | 8.0 ± 4.0 | 9.9 ± 4.0 | 5.5 ± 4.4 |
| GCS | 6.4 ± 3.6 | 5.8 ± 3.6 | 6.4 ± 3.8 | 7.0 ± 3.9 | 6.4 ± 3.6 |
| Mechanical ventilation, n (%) | 25,404 (43.0%) | 7,714 (60.2%) | 1,754 (64.3%) | 6,082 (63.4%) | 40,954 (48.7%) |
| Vasopressor use, n (%) | 19,925 (33.8%) | 6,776 (52.9%) | 1,646 (60.3%) | 6,096 (63.6%) | 34,443 (40.9%) |
| **Renal Function** | | | | | |
| Baseline Cr (mg/dL) | 1.03 ± 0.54 | 1.22 ± 0.57 | 0.90 ± 0.36 | 4.41 ± 3.46 | 1.44 ± 1.66 |
| Peak Cr (mg/dL) | 1.10 ± 0.54 | 1.76 ± 0.73 | 2.07 ± 0.85 | 6.25 ± 3.30 | 1.81 ± 2.04 |
| Cr ratio (peak/baseline) | 1.08 ± 0.11 | 1.48 ± 0.20 | 2.30 ± 0.28 | 2.13 ± 1.76 | 1.30 ± 0.72 |
| CKD (ICD), n (%) | 11,477 (19.4%) | 5,219 (40.7%) | 737 (27.0%) | 6,778 (70.7%) | 24,211 (28.8%) |
| CKD (lab, eGFR<60), n (%) | 14,449 (24.5%) | 5,563 (43.4%) | 615 (22.5%) | 7,954 (83.0%) | 28,581 (34.0%) |
| eGFR (CKD-EPI) | 81.3 ± 28.9 | 68.6 ± 29.1 | 84.3 ± 27.3 | 30.3 ± 32.6 | 73.7 ± 33.5 |
| **Comorbidities** | | | | | |
| Sepsis, n (%) | 12,808 (21.7%) | 3,752 (29.3%) | 1,084 (39.7%) | 5,004 (52.2%) | 22,648 (26.9%) |
| Congestive Heart Failure, n (%) | 17,151 (29.1%) | 6,270 (48.9%) | 1,159 (42.5%) | 5,461 (57.0%) | 30,041 (35.7%) |
| Liver Disease, n (%) | 8,499 (14.4%) | 2,538 (19.8%) | 726 (26.6%) | 3,181 (33.2%) | 14,944 (17.8%) |
| Coagulopathy, n (%) | 12,423 (21.0%) | 4,086 (31.9%) | 1,041 (38.1%) | 4,454 (46.5%) | 22,004 (26.1%) |
| **Outcomes** | | | | | |
| In-hospital mortality, n (%) | 3,530 (6.0%) | 1,932 (15.1%) | 730 (26.7%) | 2,778 (29.0%) | 8,970 (10.7%) |
| ICU LOS (hours) | 67.7 ± 89.2 | 100.5 ± 124.0 | 153.9 ± 203.5 | 150.5 ± 216.4 | 84.9 ± 124.5 |
| Hospital LOS (days) | 8.1 ± 8.4 | 13.2 ± 12.8 | 18.7 ± 18.1 | 17.1 ± 19.5 | 10.2 ± 11.9 |
| RRT, n (%) | 0 (0.0%) | 0 (0.0%) | 0 (0.0%) | 3,971 (41.4%) | 3,971 (4.7%) |

*Continuous variables: mean ± SD. Categorical: n (%). All across-stage comparisons: p < 0.001 (ANOVA or chi-square).*

### Table 2. Multivariable Logistic Regression: Comparison of Three Incremental Models (MIMIC-IV)

| Variable | Model A OR (95% CI) | *p* | Model B OR (95% CI) | *p* | Model C OR (95% CI) | *p* |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| **KDIGO Stage 1** | 2.762 (2.595–2.940) | <0.001 | 1.885 (1.765–2.013) | <0.001 | — | — |
| **KDIGO Stage 2** | 5.261 (4.775–5.797) | <0.001 | 2.978 (2.688–3.300) | <0.001 | — | — |
| **KDIGO Stage 3** | 7.029 (6.589–7.499) | <0.001 | 2.881 (2.664–3.115) | <0.001 | — | — |
| Log-Cr ratio (continuous) | — | — | — | — | 2.135 (2.003–2.276) | <0.001 |
| Age (per year) | 1.036 (1.034–1.038) | <0.001 | 1.038 (1.036–1.040) | <0.001 | 1.036 (1.034–1.038) | <0.001 |
| Male | 0.894 (0.852–0.938) | <0.001 | 0.765 (0.727–0.804) | <0.001 | 0.876 (0.836–0.918) | <0.001 |
| SOFA score | — | — | 1.333 (1.318–1.349) | <0.001 | 1.325 (1.310–1.340) | <0.001 |
| Mechanical ventilation | — | — | 0.441 (0.407–0.478) | <0.001 | 0.471 (0.435–0.509) | <0.001 |
| Vasopressor use | — | — | 0.660 (0.614–0.709) | <0.001 | 0.688 (0.641–0.739) | <0.001 |
| CKD (ICD) | — | — | 0.620 (0.582–0.660) | <0.001 | 0.626 (0.588–0.667) | <0.001 |
| Sepsis | 1.805 (1.713–1.902) | <0.001 | 1.598 (1.512–1.690) | <0.001 | 1.671 (1.583–1.764) | <0.001 |
| Metastatic Cancer | 1.547 (1.415–1.691) | <0.001 | 1.947 (1.773–2.139) | <0.001 | 1.601 (1.464–1.751) | <0.001 |
| Liver Disease | 1.590 (1.496–1.691) | <0.001 | 1.164 (1.089–1.244) | <0.001 | 1.220 (1.144–1.301) | <0.001 |
| **AIC** | 48,347 | | 44,716 | | 45,132 | |
| **AUC (95% CI)** | 0.778 (0.773–0.784) | | 0.828 (0.823–0.833) | | 0.825 (0.820–0.831) | |
| **Brier Score** | 0.0829 | | 0.0779 | | 0.0788 | |

*Reference: KDIGO Stage 0. Full covariate results (27 variables) in Table S2. VIF for all covariates <5.*

### Table 3. External Validation: KDIGO Distribution, Mortality, and Model Performance in eICU-CRD (N = 166,373)

| Metric | eICU-CRD | MIMIC-IV (comparison) |
|:---|:---:|:---:|
| **KDIGO Distribution** | | |
| Stage 0, n (%) | 113,389 (68.2%) | 59,033 (70.1%) |
| Stage 1, n (%) | 28,129 (16.9%) | 12,817 (15.2%) |
| Stage 2, n (%) | 6,523 (3.9%) | 2,729 (3.2%) |
| Stage 3, n (%) | 18,332 (11.0%) | 9,588 (11.4%) |
| AKI prevalence | 31.8% | 29.9% |
| **Mortality by Stage** | | |
| Stage 0 | 4.6% | 6.0% |
| Stage 1 | 12.5% | 15.1% |
| Stage 2 | 18.0% | 26.7% |
| Stage 3 | 19.2% | 29.0% |
| Overall | 8.1% | 10.7% |
| **Model B Adjusted ORs** | | |
| Stage 1 | 1.849 (1.761–1.942) | 1.885 (1.765–2.013) |
| Stage 2 | 2.811 (2.603–3.035) | 2.978 (2.688–3.300) |
| Stage 3 | 2.581 (2.450–2.719) | 2.881 (2.664–3.115) |
| **CKD × Stage 3 Interaction** | | |
| Interaction OR | 0.726 (0.603–0.872) | 0.476 |
| p-value | <0.001 | <0.001 |
| **Model Discrimination** | | |
| Model A AUC | 0.672 | 0.778 |
| Model B AUC | 0.818 | 0.828 |

---

## Figure Legends

- **Figure 1.** CONSORT-style patient selection flowchart. From 94,458 total ICU stays, 84,167 were included after excluding non-adult, duplicate, and 1,075 stays without serum creatinine data.
  - *Alt text:* Flowchart showing patient selection from 94,458 ICU stays to 84,167 included patients, with exclusion boxes for non-adult, duplicate, and missing creatinine data.
- **Figure 2.** In-hospital mortality rate by KDIGO AKI stage (bar chart with 95% Wilson CI).
  - *Alt text:* Bar chart of in-hospital mortality percentages across KDIGO Stages 0–3, with error bars showing 95% Wilson confidence intervals.
- **Figure 3.** Kaplan-Meier 30-day survival curves stratified by KDIGO stage, with number-at-risk table and Bonferroni-corrected pairwise log-rank test results.
  - *Alt text:* Kaplan-Meier survival curves for KDIGO Stages 0–3 over 30 days, with number-at-risk table below and pairwise log-rank p-values.
- **Figure 4.** Forest plot of adjusted odds ratios from Model B (full 27 covariate set), showing all variables with 95% CI.
  - *Alt text:* Forest plot of adjusted odds ratios from the multivariable logistic regression Model B, with 95% confidence intervals for all 27 covariates.
- **Figure 5.** Dose-response curve: creatinine ratio versus mortality rate with 95% CI.
  - *Alt text:* Line graph showing the nonlinear relationship between creatinine ratio and mortality rate, with shaded 95% confidence interval band.
- **Figure 6.** AKI paradox. Panel A: in-hospital mortality by CKD status within each KDIGO stage. Panel B: Stage 3 decomposed by RRT indication within each CKD stratum.
  - *Alt text:* Two-panel figure: Panel A shows grouped bar chart of mortality by CKD status across KDIGO stages; Panel B shows Stage 3 mortality stratified by RRT and CKD status.
- **Figure 7.** Time-dependent ROC analysis. Panel A: ROC curves at 7, 14, and 28 days for full Model B (AUC = 0.643, 0.639, 0.650). Panel B: KDIGO stage alone (AUC = 0.484, 0.452, 0.404) showing declining performance.
  - *Alt text:* Two-panel figure: Panel A shows ROC curves for the full model at 7, 14, and 28 days; Panel B shows ROC curves for KDIGO-only model at the same time points.
- **Figure 8.** Nomogram for ICU mortality prediction. Panel A: Nomogram points. Panel B: Calibration curve (Brier = 0.062). Based on eICU-CRD Model B (N = 166,373, AUC = 0.818).
  - *Alt text:* Two-panel figure: Panel A shows a nomogram with point scales for KDIGO stage, APACHE score, ventilation, age, and sex; Panel B shows calibration plot comparing predicted vs observed probabilities.

---

## Supplementary Material

### Supplementary Figures

- **Figure S1.** Restricted cubic spline analysis of the nonlinear creatinine ratio–mortality relationship.
  - *Alt text:* Line graph showing restricted cubic spline fitted curve of creatinine ratio versus log-odds of mortality, with 95% CI and reference line.
- **Figure S2.** Model discrimination and clinical utility. Panel A: ROC curves for Models A–C. Panel B: DCA showing net benefit of Model B.
  - *Alt text:* Two-panel figure: Panel A shows overlapping ROC curves for three logistic models; Panel B shows decision curve analysis with net benefit across threshold probabilities.
- **Figure S3.** Sepsis subgroup analysis. Forest plot comparing adjusted ORs between sepsis and non-sepsis subgroups.
  - *Alt text:* Forest plot comparing adjusted odds ratios for KDIGO stages between sepsis and non-sepsis subgroups with 95% CIs.
- **Figure S4.** Calibration plot for Model B.
  - *Alt text:* Calibration plot showing predicted versus observed mortality probabilities for Model B, with a reference diagonal line.
- **Figure S5.** Directed acyclic graph (DAG) representing the assumed causal structure.
  - *Alt text:* Directed acyclic graph with nodes for KDIGO stage, illness severity, CKD, mortality, and confounders, with arrows showing assumed causal directions.
- **Figure S6.** Variable importance forest plot from the nomogram model.
  - *Alt text:* Horizontal bar chart showing variable importance scores for the five nomogram predictors (KDIGO, APACHE, ventilation, age, sex).
- **Supplementary File.** Web-based mortality risk calculator (standalone HTML file).
- **STROBE Checklist.** Completed STROBE checklist for observational cohort studies.

### Supplementary Tables

- **Table S1.** Systematic missing data report stratified by KDIGO stage.
- **Table S2.** Full logistic regression results for Model B (all 27 covariates).
- **Table S3.** Time-stratified Cox proportional hazards regression (0–7, 7–14, 14–30 days).
- **Table S4.** Temporal trend analysis (2008–2022) by 3-year periods.
- **Table S5.** 24-hour landmark sensitivity analysis.
- **Table S6.** Dose-response analysis: creatinine ratio binned mortality rates.
- **Table S7.** AKI paradox: CKD-stratified mortality, RRT-stratified decomposition, and formal CKD × KDIGO interaction test.
- **Table S8.** E-value sensitivity analysis for Model B associations.

### Supplementary Methods

Detailed descriptions of: (1) SOFA score calculation and missing component handling; (2) MDRD back-calculation of baseline creatinine; (3) CKD-EPI equation for laboratory-based CKD; (4) Elixhauser comorbidity derivation from ICD-9/10 codes; (5) complete statistical analysis plan including VIF diagnostics, dose-response binning, RCS knot placement, DCA methodology, and time-dependent ROC framework; (6) eICU-CRD KDIGO staging implementation and APACHE IVa extraction; (7) nomogram construction methodology and web calculator deployment; and (8) sensitivity analysis details (MDRD, lab-CKD, complete-case SOFA, landmark, temporal trends).

---

## Funding

This work was supported by the Chronic Disease Management Research Project of National Health Commission Capacity Building and Continuing Education Center (Grant No. GWJJMB202510024181), Changsha Science and Technology Bureau Project (Grant No. kq2014242), and Natural Science Foundation of Hunan Province (Grant No. 2021JJ30959).

---

## Data Availability

MIMIC-IV and eICU-CRD are publicly available at PhysioNet. All analysis code is at GitHub (https://github.com/wudengke2010/kdigo-aki-mimic-iv).

---

## Ethics Approval

Access to MIMIC-IV and eICU-CRD was obtained through PhysioNet after CITI Data or Specimens Only Research course completion. Use of de-identified, publicly available data did not require additional IRB approval.

---

## Consent for Publication

Not applicable. This study used de-identified, publicly available data.

---

## Patient and Public Involvement

Patients and/or the public were not involved in the design, conduct, reporting, or dissemination of this research.

---

## Author Contributions

**Jiqiang Liu:** Conceptualisation, Data curation, Formal analysis, Investigation, Methodology, Writing – original draft. **Dengke Wu:** Conceptualisation, Funding acquisition, Project administration, Supervision, Writing – review & editing.

---

## Disclosure of Interest

The authors report no competing interests.

---

## ORCID

Dengke Wu https://orcid.org/0009-0008-1363-9621

---

## Declaration of Generative AI in Scientific Writing

During preparation, the authors used LLMs for language editing. After using these tools, the authors reviewed and edited the content and take full responsibility for the publication.
