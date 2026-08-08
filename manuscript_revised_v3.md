# Association of KDIGO Creatinine-Based AKI Staging with ICU Mortality: A Dual-Cohort Retrospective Study Using MIMIC-IV v3.1 and eICU-CRD

---

## Authors

**Jiqiang Liu**<sup>1,2</sup>, Dengke Wu<sup>1,2*</sup>

<sup>1</sup> Department of Emergency Medicine, The Second Xiangya Hospital of Central South University, Changsha 410011, Hunan, China

<sup>2</sup> Emergency Medicine and Difficult Diseases Institute, The Second Xiangya Hospital of Central South University, Changsha 410011, Hunan, China

**\*Corresponding author:** Dengke Wu, Department of Emergency Medicine, Second Xiangya Hospital, Central South University, Changsha, Hunan, China; Emergency Medicine and Difficult Diseases Institute, The Second Xiangya Hospital of Central South University, Changsha 410011, Hunan, China. Electronic address: wudk2010@csu.edu.cn.

---

## Abstract

**Background:** The KDIGO creatinine-based AKI staging system is widely used, but illness severity confounding, Stage 2–3 mortality convergence, and the "AKI paradox" remain incompletely characterised. No large study has simultaneously quantified these phenomena with formal reclassification metrics, time-stratified survival analysis, and external validation.

**Methods:** This STROBE-compliant retrospective cohort study used two independent databases: MIMIC-IV v3.1 (2008–2022; 84,167 ICU stays; derivation cohort) and eICU-CRD v2.0 (166,373 ICU stays; external validation cohort; total N = 250,540). AKI was staged by KDIGO creatinine criteria. Three incremental logistic models were fitted, with severity adjustment (SOFA in MIMIC; APACHE IVa in eICU). Cox regression with Schoenfeld testing, Kaplan-Meier with Bonferroni correction, restricted cubic splines, RRT-stratified decomposition, 24-hour landmark analysis, and time-stratified Cox models were performed. Time-dependent ROC analysis (7/14/28-day) and a clinical nomogram with web-based risk calculator were developed.

**Results:** AKI occurred in 25,134 patients (29.9%) in MIMIC-IV and 52,984 (31.8%) in eICU-CRD. Mortality was 6.0%/15.1%/26.7%/29.0% (MIMIC Stages 0–3) and 4.6%/12.5%/18.0%/19.2% (eICU). Adjusting for severity attenuated KDIGO ORs by 32–59% (NRI = 0.625; AUC: 0.778→0.828 in MIMIC; 0.672→0.818 in eICU). The AKI paradox was replicated in both cohorts: CKD × Stage 3 interaction OR = 0.476 (p < 0.001, MIMIC) and 0.726 (p < 0.001, eICU). Time-dependent AUCs for mortality prediction were 0.643/0.639/0.650 at 7/14/28 days. A nomogram integrating KDIGO stage, APACHE score, mechanical ventilation, age, and sex (AUC = 0.818, Brier = 0.062) was deployed as a web-based calculator.

**Conclusions:** KDIGO staging is independently associated with ICU mortality across two independent cohorts (total N = 250,540), but illness severity accounts for a substantial portion of this association. The Stage 3-specific AKI paradox is reproducible across databases, supporting its biological rather than artefactual nature. A clinically deployable nomogram and web calculator provide practical risk stratification.

**Keywords:** acute kidney injury; KDIGO; intensive care unit; mortality; chronic kidney disease; renal replacement therapy

---

## Introduction

Acute kidney injury (AKI) is one of the most common and consequential complications encountered in the intensive care unit (ICU). Epidemiological studies have consistently demonstrated that AKI affects approximately 30–60% of critically ill patients, with the landmark AKI-EPI multinational cohort reporting an incidence of 57.7% and a mortality gradient from 16.1% (no AKI) to 50.0% (Stage 3) [1,2].

The evolution of consensus definitions—from RIFLE to AKIN to the current Kidney Disease: Improving Global Outcomes (KDIGO) criteria—has standardised AKI classification and revealed a graded association between AKI severity and mortality across diverse ICU populations [19]. The KDIGO Clinical Practice Guideline for Acute Kidney Injury, published in 2012, established a consensus definition and staging system integrating serum creatinine (SCr), urine output, and renal replacement therapy (RRT) criteria [3]. Several studies have validated the prognostic value of KDIGO staging in general ICU cohorts, cardiac surgery patients, and sepsis populations [4–7]. However, the majority of these validation studies were conducted using earlier MIMIC versions, single-center designs, or selected subgroups.

Despite the widespread adoption of KDIGO staging, several critical knowledge gaps persist that limit its clinical utility and generalisability:

**First**, the degree to which the KDIGO-mortality association is confounded by overall illness severity—as captured by the Sequential Organ Failure Assessment (SOFA) score, mechanical ventilation, and vasopressor requirement—has not been adequately quantified in large cohorts with formal reclassification metrics (NRI/IDI). Girling et al. [24] applied the Bradford Hill criteria to assess whether AKI is causally associated with critical illness mortality or merely a correlate of illness severity, identifying plausibility, temporality, and a biological gradient but concluding that formal confounding quantification is needed to resolve the causal-versus-correlational question. Without such quantification, clinicians cannot distinguish the independent prognostic contribution of AKI staging from the severity of multi-organ failure that accompanies it.

**Second**, while a dose-response relationship between KDIGO stage and mortality has been described [8], the incremental prognostic value of higher stages—particularly the distinction between Stage 2 and Stage 3—remains debated, partly because Stage 3 is enriched with chronic kidney disease (CKD) patients who exhibit different mortality trajectories [9]. Whether the observed Stage 2–3 mortality convergence represents a true biological ceiling or a staging artefact has not been resolved through time-stratified analysis.

**Third**, a counterintuitive phenomenon termed the "AKI paradox"—where CKD patients appear to tolerate severe AKI better than non-CKD patients—has been observed [10,11], but the underlying mechanisms (differential RRT initiation, staging artifact, biological adaptation) await systematic decomposition in a cohort large enough to permit formal interaction testing and within-stage stratification. Crucially, whether this paradox is reproducible across independent databases—a prerequisite for distinguishing a biological phenomenon from a single-cohort artefact—has not been established.

**Fourth**, the sensitivity of KDIGO staging to baseline creatinine estimation and CKD definition has not been rigorously examined in large ICU cohorts, and no study has translated KDIGO-based risk stratification into a practical, deployable clinical prediction tool (e.g., nomogram with web calculator) validated in an external cohort. Recent machine-learning approaches have advanced AKI mortality prediction in specific subgroups including sepsis-associated AKI [26] and pancreatitis-associated AKI [27], but these models require computational infrastructure not readily available at the bedside, and none has integrated KDIGO staging into an interpretable, transportable nomogram validated across databases.

MIMIC-IV version 3.1 provides a unique opportunity to address these gaps in the derivation cohort. With data spanning 15 years (2008–2022) and encompassing 84,167 ICU stays with comprehensive laboratory, vital sign, comorbidity, and outcome data, it offers sufficient power for detailed stratified analyses and sensitivity testing [12]. Critically, the eICU Collaborative Research Database (eICU-CRD) v2.0, comprising 200,859 ICU stays from 208 hospitals across the United States [23], provides an independent external validation cohort to test the reproducibility of our findings—particularly the AKI paradox—across different data collection platforms, hospital systems, and severity scoring methodologies (APACHE IVa vs. SOFA). The precedent for this dual-database design has been established by Lin et al. [25], who used MIMIC-IV and eICU in a derivation/validation pairing for trajectory-based AKI subphenotyping in 17,113 ICU patients.

The objective of this study was to comprehensively evaluate the association between KDIGO creatinine-based AKI staging and ICU mortality, with explicit quantification of illness severity confounding and systematic investigation of the AKI paradox, validated in an independent external cohort. Specifically, we aimed to: (1) quantify the degree to which illness severity confounds the KDIGO-mortality association using incremental regression models with NRI/IDI reclassification metrics; (2) characterise survival trajectories across stages with conservative multiple comparison correction and time-stratified Cox analysis; (3) decompose the AKI paradox by RRT indication within Stage 3; (4) examine the nonlinear creatinine ratio–mortality relationship using restricted cubic splines; (5) evaluate methodological robustness through MDRD-based baseline creatinine, laboratory-based CKD, complete-case SOFA, and landmark sensitivity analyses; (6) externally validate all primary findings in eICU-CRD; (7) assess time-dependent predictive performance using ROC analysis at 7, 14, and 28 days; and (8) develop and deploy a clinical nomogram with web-based risk calculator for practical mortality risk stratification.

---

## Methods

### Study Design and Data Source

This was a retrospective dual-cohort study using two independent, publicly available critical care databases. The **derivation cohort** used the Medical Information Mart for Intensive Care IV (MIMIC-IV) database, version 3.1, containing de-identified electronic health records of patients admitted to the ICUs at Beth Israel Deaconess Medical Center (BIDMC) in Boston, Massachusetts, between 2008 and 2022 [12]. The **external validation cohort** used the eICU Collaborative Research Database (eICU-CRD), version 2.0, a multi-centre database of ICU stays from 208 hospitals across the United States, collected between 2014 and 2015 [23]. **This study is reported in accordance with the Strengthening the Reporting of Observational Studies in Epidemiology (STROBE) statement** (Supplementary Material).

### Study Population

In MIMIC-IV, all adult patients (≥18 years) with at least one ICU admission during the study period were eligible. For patients with multiple ICU stays within a single hospital admission, only the first ICU stay was retained. ICU stays without at least one SCr measurement were excluded (n = 1,075; 1.3%). In eICU-CRD, the same adult and first-ICU-stay criteria were applied. Exclusions are detailed in Figure 1.

### Exposure: KDIGO Creatinine-Based AKI Staging

The primary exposure was the maximum AKI severity during the ICU stay, classified according to the KDIGO 2012 creatinine criteria. Because urine output data are incompletely recorded in both databases, staging was based exclusively on creatinine criteria:

- **No AKI (Stage 0):** Neither Stage 1–3 criteria met.
- **Stage 1:** Increase in SCr to ≥1.5 times baseline, or increase of ≥0.3 mg/dL (26.5 µmol/L).
- **Stage 2:** Increase in SCr to ≥2.0 times baseline and <3.0 times baseline.
- **Stage 3:** Increase in SCr to ≥3.0 times baseline, or SCr ≥4.0 mg/dL (353.6 µmol/L), or initiation of RRT.

This implementation was verified against the original KDIGO guideline criteria to avoid the common staging error of interpreting the Stage 3 absolute threshold as an increase *of* ≥4.0 mg/dL rather than an increase *to* ≥4.0 mg/dL, an error identified in multiple published studies [20].

Baseline SCr was defined as the first SCr measurement during the hospital admission. Peak SCr was the highest value recorded during the ICU stay. Creatinine ratio (Cr ratio) was calculated as peak SCr divided by baseline SCr. In MIMIC-IV, SCr measurements were extracted using itemid 50912; values <0.1 mg/dL or >100 mg/dL were excluded. In eICU-CRD, SCr was extracted from the `lab` table (labname = 'creatinine'). RRT was identified from ICU procedure events (MIMIC) or the `treatment` table (eICU). Patients receiving RRT were automatically classified as Stage 3.

### Immortal Time Bias Assessment

Because KDIGO staging uses the entire ICU stay's peak creatinine, patients who died very early (before creatinine had risen sufficiently) may be misclassified as lower stages. A 24-hour landmark sensitivity analysis was performed, restricting the cohort to patients who survived beyond the first 24 hours of ICU admission.

**RRT Stratification.** To investigate whether the AKI paradox is driven by differential RRT indication, Stage 3 patients were stratified into three mutually exclusive subgroups: (1) RRT-only (met RRT criterion without meeting creatinine criteria), (2) Creatinine-only (met creatinine criteria without receiving RRT), and (3) Both (met creatinine criteria and received RRT).

**Sensitivity Analysis: MDRD Baseline Creatinine.** Baseline SCr was back-calculated from the Modification of Diet in Renal Disease (MDRD) equation assuming an estimated glomerular filtration rate (eGFR) of 75 mL/min/1.73 m². KDIGO staging was repeated using MDRD-derived baseline, and agreement with the primary analysis was quantified.

### Outcomes

The **primary outcome** was in-hospital mortality. **Secondary outcomes** included 30-day survival (from ICU admission, censored at hospital discharge or 30 days), ICU length of stay (LOS), and hospital LOS.

### Covariates

**Severity of Illness (MIMIC-IV).** The SOFA score was calculated from six components. Missing components were handled as follows: for neurological SOFA (Glasgow Coma Scale [GCS]), the SOFA sub-score defaults to 0 when GCS is unavailable (61.2% missing), consistent with the extraction code default parameter. For hepatic SOFA, total bilirubin (49.8% missing) was similarly handled. This approach systematically underestimates true illness severity, a limitation addressed analytically through complete-case sensitivity analysis. The maximum SOFA score during the first 24 hours of ICU stay was used. Mechanical ventilation and vasopressor use were identified from chartevents during the first 24 hours.

**Severity of Illness (eICU-CRD).** The APACHE IVa score was extracted from the `apachePatientResult` table, providing a validated multi-organ severity score conceptually equivalent to SOFA but independently derived. Mechanical ventilation was identified from the `apachePredVar` table (admission source and ventilated state).

**Demographics and Comorbidities.** Age at admission and sex were extracted from both databases. In MIMIC-IV, race and eighteen Elixhauser comorbidity categories were derived from ICD-9 and ICD-10 diagnosis codes [13,14]. CKD was defined by ICD codes (ICD-9: 585.x; ICD-10: N18.x). In eICU-CRD, CKD was identified from the `diagnosis` table using ICD-9 and ICD-10 codes mapped to the same definitions.

**Vital Signs and Admission Laboratory Values.** Heart rate, systolic/diastolic/mean blood pressure, respiratory rate, temperature, and SpO₂ were extracted from ICU chartevents with physiological range filtering. Admission laboratory values (first 24 hours) included BUN, sodium, potassium, bicarbonate, lactate, WBC, haemoglobin, platelet count, total bilirubin, albumin, and glucose.

**Sensitivity Analysis: Laboratory-Based CKD.** A laboratory-based CKD definition was implemented using the Chronic Kidney Disease Epidemiology Collaboration (CKD-EPI) equation: eGFR <60 mL/min/1.73 m², calculated from the lowest SCr during the hospital admission, age, and sex.

### Missing Data Assessment

A systematic missing data report was generated for all key variables (Table S1). Variables with substantial missingness (>5%) included GCS (61.2%), albumin (64.8%), systolic/diastolic/MAP blood pressure (65.1–65.2%), total bilirubin (49.8%), and lactate (38.5%). Table S1 provides both overall and KDIGO stage-stratified missing rates. SOFA component scores were assigned a value of 0 (normal) through the default parameter when source data were absent. Variables with low missingness (≤2%) were retained with missing values treated as normal in logistic regression.

### Statistical Analysis

**Multivariable logistic regression.** Three incremental models were constructed:
- **Model A:** KDIGO stage (reference: Stage 0) + age + sex + 18 Elixhauser comorbidities.
- **Model B:** Model A + SOFA score (MIMIC) / APACHE IVa (eICU) + mechanical ventilation + vasopressor use (MIMIC) + ICD-based CKD.
- **Model C:** Continuous log-transformed creatinine ratio (replacing KDIGO stages) + all Model A covariates, evaluating the continuous creatinine change versus categorical staging.

Model fit was compared using AIC. The attenuation of KDIGO odds ratios from Model A to Model B quantified illness severity confounding. Variance inflation factors (VIF) were computed for all Model B covariates [17]; VIF >5 was considered indicative of moderate collinearity.

**Kaplan-Meier survival analysis.** 30-day survival curves stratified by KDIGO stage were compared using the overall log-rank test with Bonferroni-corrected pairwise comparisons (adjusted α = 0.05/6 = 0.0083).

**Cox proportional hazards regression.** Two Cox models were fitted corresponding to Models A and B. The proportional hazards (PH) assumption was tested using Schoenfeld residuals. Because PH violations were detected (see Results), **time-stratified Cox models** (0–7, 7–14, 14–30 days) were fitted as sensitivity analyses. This approach follows published CKJ precedent: Singh et al. [28] encountered PH violations in a bicentric ICU-AKI cohort and switched to time-stratified odds ratios, confirming that time-resolved analysis is necessary in AKI-mortality studies. Restricted mean survival time (RMST) at 30 days was computed as a PH-independent complementary metric.

**Dose-response analysis.** The creatinine ratio was binned into seven categories (<1.0, 1.0–1.5, 1.5–2.0, 2.0–3.0, 3.0–5.0, 5.0–10.0, ≥10.0). Given the small ≥10.0 subgroup (n = 64), this group was pooled with 5.0–10.0 for a combined estimate (n = 609).

**Restricted cubic splines (RCS).** A 4-knot RCS (10th, 35th, 65th, 90th percentiles) was fitted within the fully adjusted Model B framework to model the nonlinear creatinine ratio–mortality relationship. Nonlinearity was tested using a likelihood ratio test.

**Model discrimination and reclassification.** The area under the receiver operating characteristic curve (AUC) was computed for all models. The net reclassification improvement (NRI) and integrated discrimination improvement (IDI) quantified the incremental value of severity adjustment (Model A → Model B). Decision curve analysis (DCA) evaluated net benefit across threshold probabilities.

**Calibration.** The Hosmer-Lemeshow goodness-of-fit test (10 deciles) and visual calibration plots assessed calibration.

**Subgroup analysis.** Sepsis status (ICD codes) was tested as an effect modifier via multiplicative interaction terms in Model B.

**CKD × KDIGO interaction.** Formal interaction tests for both ICD-CKD and laboratory-CKD with KDIGO stage were conducted, with joint significance assessed by likelihood ratio tests.

**Sensitivity analyses.** Five sensitivity analyses were performed: (1) MDRD-based KDIGO staging with Model B equivalent; (2) Laboratory-based CKD (CKD-EPI eGFR <60) versus ICD-CKD in the fully adjusted model; (3) Complete-case SOFA analysis (patients with both GCS and bilirubin available, N = 17,285) to assess the impact of SOFA imputation; (4) 24-hour landmark analysis excluding patients who died within the first 24 ICU hours; (5) Temporal trend analysis using MIMIC-IV `anchor_year_group` to recover real admission years (2008–2022) and compute 3-year period trends.

### External Validation in eICU-CRD

The primary analysis pipeline was replicated in eICU-CRD using identical KDIGO staging criteria, with APACHE IVa replacing SOFA as the severity score. Logistic regression (Models A and B equivalent) and Cox regression were fitted. The CKD × KDIGO interaction was tested to assess whether the AKI paradox replicates in an independent cohort. Model discrimination (AUC) was compared between cohorts.

### Time-Dependent ROC Analysis

To evaluate the temporal dynamics of mortality prediction, time-dependent ROC curves were computed at 7, 14, and 28 days using the cumulative sensitivity / dynamic specificity framework. At each time horizon t, cases were defined as patients who died within t days, and controls as patients who survived beyond t days. Risk scores were derived from the MIMIC-IV Cox Model B. For comparison, KDIGO stage alone was evaluated as a sole predictor at each time horizon.

### Nomogram and Web-Based Risk Calculator

A clinical prediction nomogram was developed based on the eICU-CRD Model B (chosen for its multi-centre generalisability), incorporating KDIGO stage, APACHE score, mechanical ventilation, age, and sex. Logistic regression coefficients were scaled to a 0–100 point nomogram. Calibration was assessed using the Brier score and a calibration curve (10 quantiles). A standalone HTML web calculator was deployed, embedding the regression coefficients for real-time mortality risk prediction.

### Causal Structure and E-value Sensitivity Analysis

The causal structure was a priori defined using a directed acyclic graph (DAG) (Supplementary Figure S5). Measured confounders were identified as age, sex, Elixhauser comorbidities, SOFA score, mechanical ventilation, vasopressor use, and CKD status. Variables identified as colliders—including ICU admission and serum creatinine measurement availability—were not adjusted for. To assess robustness to unmeasured confounding, E-values were computed for the primary Model B associations (VanderWeele & Ding, 2017) [22]. The analysis plan was developed before data extraction and is documented in this manuscript and its supplementary materials.

All tests were two-sided. p-values reported as 0.000 by statistical software are presented as p < 0.001. Analyses were performed using Python 3.13.12 with pandas 2.3.3, statsmodels 0.14.6, lifelines 0.30.3, scipy 1.18.0, and matplotlib 3.11.0.

---

## Results

### Derivation Cohort (MIMIC-IV)

#### Study Population

From 94,458 total ICU stays, 85,242 met adult and first-ICU criteria. After excluding 1,075 stays without serum creatinine measurements, 84,167 ICU stays were included in the final analysis (Figure 1). The 1,075 exclusions represented 1.3% of eligible stays.

#### KDIGO Stage Distribution

AKI (Stage 1–3) occurred in 25,134 ICU stays (29.9%). The distribution was: Stage 0: 59,033 (70.1%); Stage 1: 12,817 (15.2%); Stage 2: 2,729 (3.2%); Stage 3: 9,588 (11.4%). Among Stage 3 patients, 3,971 (41.4%) received RRT. Stage 3 composition by mechanism: RRT-only: 793 (8.3%); Cr-only: 5,617 (58.6%); Both: 3,178 (33.1%).

#### Baseline Characteristics

Baseline characteristics stratified by KDIGO stage are presented in **Table 1**. The mean age was 63.0 ± 16.8 years, and 55.8% of patients were male. Age differed across stages (p < 0.001), with Stage 1 patients being oldest (67.1 ± 14.9 years).

Illness severity increased markedly with AKI stage. Mean SOFA rose from 4.3 ± 3.9 (Stage 0) to 9.9 ± 4.0 (Stage 3). Mechanical ventilation rates increased from 43.0% to 63.4%, and vasopressor use from 33.8% to 63.6% (both p < 0.001).

CKD (ICD-defined) was present in 24,211 patients (28.8%), increasing from 19.4% in Stage 0 to 70.7% in Stage 3 (p < 0.001). Laboratory-based CKD (eGFR <60) was present in 28,581 (34.0%). Notably, Stage 2 patients had the lowest baseline SCr (0.90 ± 0.36 mg/dL) and consequently the highest mean eGFR (84.3 ± 27.3 mL/min/1.73 m²), exceeding even Stage 0 (81.3 ± 28.9). This reflects the fold-change nature of KDIGO staging: patients with low baseline creatinine require smaller absolute increases to reach the ≥2.0-fold Stage 2 threshold, thereby enriching Stage 2 with patients having preserved baseline renal function. Similarly, Stage 3's mean creatinine ratio (2.13 ± 1.76) is lower than Stage 2's (2.30 ± 0.28) because Stage 3 includes RRT-only patients (who may have ratios <1.5) and patients qualifying by the ≥4.0 mg/dL absolute threshold rather than the fold-change criterion.

Comorbidity burden increased across stages: sepsis from 21.7% (Stage 0) to 52.2% (Stage 3), congestive heart failure from 29.1% to 57.0%, liver disease from 14.4% to 33.2%, and coagulopathy from 21.0% to 46.5% (all p < 0.001).

#### Primary Outcome: In-Hospital Mortality

| KDIGO Stage | Mortality, n (%) | 95% CI (Wilson) |
|:---:|:---:|:---:|
| Stage 0 (n = 59,033) | 3,530 (6.0%) | 5.8–6.2% |
| Stage 1 (n = 12,817) | 1,932 (15.1%) | 14.5–15.8% |
| Stage 2 (n = 2,729) | 730 (26.7%) | 25.1–28.4% |
| Stage 3 (n = 9,588) | 2,778 (29.0%) | 28.1–29.9% |
| **Overall (N = 84,167)** | **8,970 (10.7%)** | **10.5–10.9%** |

Mortality increased 2.5-fold from Stage 0 to Stage 1, 4.4-fold to Stage 2, and 4.8-fold to Stage 3, with a modest Stage 2→3 increment (26.7% vs. 29.0%).

#### Kaplan-Meier Survival Analysis

Kaplan-Meier analysis demonstrated significant separation of 30-day survival curves (overall log-rank χ² = 4,737.29, p < 0.001). 30-day survival: Stage 0: 91.4%; Stage 1: 82.0%; Stage 2: 70.0%; Stage 3: 69.7%. Bonferroni-corrected pairwise tests confirmed significant differences between all adjacent stages (all corrected p < 0.001) **except** Stage 2 versus Stage 3 (raw p = 0.805, Bonferroni p = 1.000). **Restricted mean survival time (RMST) at 30 days**, a PH-independent metric, showed: Stage 0: 28.4 days, Stage 1: 26.0 days, Stage 2: 23.5 days, Stage 3: 23.5 days, confirming the Stage 2–3 convergence.

#### Multivariable Logistic Regression

**Model A** (AIC = 48,347): Stage 1 aOR 2.762 (95% confidence interval [CI] 2.595–2.940), Stage 2 aOR 5.261 (4.775–5.797), Stage 3 aOR 7.029 (6.589–7.499); all p < 0.001.

**Model B** (AIC = 44,716): Stage 1 aOR 1.885 (95% CI 1.765–2.013), Stage 2 aOR 2.978 (2.688–3.300), Stage 3 aOR 2.881 (2.664–3.115); all p < 0.001.

**E-value sensitivity analysis** supported strong robustness to unmeasured confounding. For Model B, the E-values were 2.98 for Stage 1, 4.77 for Stage 2, and 4.62 for Stage 3; the lower-bound E-values (based on the 95% CI lower limits) were 2.77, 4.32, and 4.28, respectively. These values indicate that an unmeasured confounder would need to be associated with both KDIGO stage and mortality by risk ratios of at least 2.8–4.8 to fully explain away the observed associations.

Illness severity adjustment substantially attenuated KDIGO odds ratios: Stage 1 by 31.8%, Stage 2 by 43.4%, Stage 3 by 59.0%, with a large improvement in model fit (ΔAIC = 3,631). In Model B, SOFA was the strongest predictor (aOR 1.333 per point, 95% CI 1.318–1.349). Other significant predictors included metastatic cancer (aOR 1.947), sepsis (aOR 1.598), and CKD (aOR 0.620, protective). Mechanical ventilation (aOR 0.441) and vasopressor use (aOR 0.660) showed apparently paradoxical protective associations, reflecting statistical suppression: when SOFA absorbs shared severity variance, the residual coefficients for ventilation/vasopressors largely represent less severely ill patients who received these interventions (e.g., elective post-surgical ventilation) and who indeed have better outcomes.

**Model C** (AIC = 45,132): log-creatinine ratio aOR 2.135 (95% CI 2.003–2.276), confirming the independent prognostic value of continuous creatinine change, though with lower discriminative performance than Model B.

#### 24-Hour Landmark Analysis: Immortal Time Bias Assessment

A total of 1,658 patients (18.5% of all deaths) died within the first 24 hours of ICU admission. Among them, 54.6% were classified as KDIGO Stage 0—patients whose creatinine had not yet risen sufficiently to meet KDIGO criteria at the time of death. This detection-timing artefact is well-characterised: Esposito et al. [32] demonstrated in 56,820 hospitalized patients that AKI recognition patterns are highly heterogeneous, with onset timing influencing outcomes. Excluding these early deaths (landmark N = 82,509), Model B ORs were consistently higher than the full-cohort estimates: Stage 1 landmark OR 2.127 (1.984–2.280), Stage 2 landmark OR 3.644 (3.274–4.056), Stage 3 landmark OR 3.341 (3.074–3.631), representing 12–29% increases over full-cohort ORs (Table S5).

The 12–29% higher ORs in the landmark analysis indicate that using the entire ICU stay's peak creatinine to define KDIGO stage **underestimates** the true KDIGO-mortality association in the full cohort, because early deaths (classified as low stage due to insufficient time for creatinine rise) dilute the mortality signal of higher KDIGO stages. This finding supports the conservative nature of the primary analysis: the reported KDIGO ORs likely represent lower-bound estimates of the true association.

#### VIF Diagnostics

All VIF values were below the conventional threshold of 5: SOFA 4.73, mechanical ventilation 2.68, vasopressor use 2.19, KDIGO Stage 3 1.42, age 1.03. The protective associations of ventilation and vasopressors are therefore attributable to the suppression phenomenon described above, not to invalid model specification.

#### Cox Proportional Hazards and Time-Stratified Analysis

Cox regression results were consistent with logistic models: Model A HRs were 2.118, 3.462, and 4.124 for Stages 1–3 (C-index 0.753); Model B HRs were 1.624, 2.255, and 2.037 (C-index 0.783). Severity adjustment attenuated HRs by 23.3%, 34.9%, and 50.6%.

**PH assumption testing.** Schoenfeld residuals tests identified PH violations for KDIGO Stage 3 (χ² = 84.6, p < 0.001), SOFA score (χ² = 380.3, p < 0.001), mechanical ventilation (χ² = 153.6, p < 0.001), and several comorbidities. While such violations are common in large cohorts (N = 84,167), they indicate that the reported HRs represent average effects over follow-up and may obscure time-varying patterns.

**Time-stratified Cox analysis** (Table S3) revealed distinct temporal patterns.

The Stage 3 HR progressively increased from 1.41 (0–7d) to 1.61 (7–14d) to 1.86 (14–30d), confirming that the PH violation for Stage 3 reflects a genuine biological phenomenon: Stage 3 AKI exerts a stronger influence on later mortality, after the acute phase, likely mediated by persistent organ dysfunction. Stage 2 effects were most prominent in the acute phase (HR 1.54, 0–7d), declining to 1.22 (7–14d, 95% CI 0.99–1.49) and partially recovering to 1.38 (14–30d). Stage 1 showed a moderate peak at 7–14 days (HR 1.42, 95% CI 1.23–1.64) and was no longer significant in the late period (HR 1.16, 95% CI 0.92–1.47, 14–30d). The full 30-day HR of 1.46 for Stage 3 represents a weighted average of these time-varying effects and should be interpreted as the "average hazard" rather than a constant effect.

#### Dose-Response Analysis

Creatinine ratio was binned into 7 categories (Table S6). Mortality increased steeply between ratios of 1.5–3.0 (from 19.5% to 33.6%) and plateaued at 3.0–10.0 (~41%). The ≥10.0 subgroup (n = 64) had an unreliable point estimate (29.7%, wide CI: 18.8–40.6%); pooling with the 5.0–10.0 group (combined n = 609, 40.4%) confirmed a stable plateau rather than a decline at extreme ratios.

#### AKI Paradox: CKD Stratification and RRT Decomposition

Stratification by CKD status revealed a striking interaction concentrated at Stage 3: CKD-negative Stage 3 mortality was 49.4% versus CKD-positive 20.5%, a 2.4-fold difference (Table S7).

RRT-stratified decomposition of Stage 3 (Table S7) showed that the paradox persisted in the Cr-only subgroup: CKD-negative Cr-only mortality was 42.3% versus CKD-positive 18.5%, a 2.3-fold difference. This demonstrates that the paradox is not solely attributable to differential RRT initiation. The RRT-only subgroup had the highest mortality in both CKD strata (73.7% and 44.0%), consistent with the extreme illness severity driving RRT initiation without creatinine criteria.

#### Restricted Cubic Spline Analysis

The 4-knot RCS confirmed significant nonlinearity (χ² = 161.29, df = 2, p < 0.001). Three phases were observed: (1) a neutral phase at ratios <1.0 (OR ≈ 0.83–0.98); (2) a steep ascent from 1.0 to ~2.6, peaking at OR = 2.82 (95% CI 2.62–3.05); and (3) a gradual decline from ~2.6 to 15.0, with OR remaining elevated at approximately 1.99 (95% CI 1.51–2.62) at extreme ratios. The plateau and gradual decline at extreme ratios is consistent with the dose-response analysis and supports the interpretation that the most extreme creatinine ratios may reflect enrichment of CKD patients with high baseline values rather than catastrophic acute injury.

#### Model Discrimination, Reclassification, and Calibration

| Model | AUC (95% CI) | Brier Score |
|:---|:---:|:---:|
| Model A | 0.778 (0.773–0.784) | 0.0829 |
| Model B | 0.828 (0.823–0.833) | 0.0779 |
| Model C | 0.825 (0.820–0.831) | 0.0788 |

NRI for Model A → B was 0.625 (events: 0.277; non-events: 0.347); IDI was 0.058. DCA showed Model B provided positive net benefit across threshold probabilities of 1–97%, with the greatest incremental benefit at 10–20% thresholds (net benefit 0.054–0.019 vs. −0.276–0.007 for "treat all"). Calibration was good visually across most risk deciles, with the Hosmer-Lemeshow test significant (χ² = 93.05, p < 0.001) primarily due to minor overestimation at the lowest predicted risk decile (predicted 0.8% vs. observed 0.3%), a common artifact of large sample sizes.

#### Sepsis Subgroup Analysis

Sepsis significantly modified the KDIGO–mortality association (interaction χ² = 81.3, df = 3, p < 0.001): Stage 1 sepsis aOR 2.697 (vs. non-sepsis 1.525), Stage 2 sepsis aOR 4.110 (vs. 2.321), Stage 3 sepsis aOR 3.449 (vs. 2.501). KDIGO ORs were 1.3–1.8-fold higher in sepsis, especially at Stages 1–2 (Figure S3). The attenuation at Stage 3 is consistent with a ceiling effect at the most severe AKI level.

#### CKD × KDIGO Interaction: Formal Test of the AKI Paradox

Formal interaction test results are presented in Table S7. The overall interaction was highly significant for both ICD-CKD (χ² = 118.83, df = 3, p < 0.001) and lab-CKD (χ² = 274.72, df = 3, p < 0.001). The interaction was exclusively at Stage 3: CKD reduces the mortality odds of Stage 3 AKI by 52% (ICD-CKD interaction OR = 0.476, p < 0.001) and 68% (lab-CKD interaction OR = 0.318, p < 0.001) relative to non-CKD Stage 3 patients. No significant interaction was observed at Stages 1 or 2, confirming that the AKI paradox is a Stage 3-specific phenomenon.

#### Sensitivity Analyses

**MDRD baseline creatinine.** MDRD-based staging shifted proportions: Stage 0: 47,545 (56.5%), Stage 1: 15,629 (18.6%), Stage 2: 8,506 (10.1%), Stage 3: 12,487 (14.8%). Agreement with the primary analysis was 68.9%, with most reclassification from Stage 0 to higher stages. MDRD-based Stage 3 aOR was higher (4.750 vs. 2.881 in the primary analysis), suggesting MDRD back-calculation may identify a more severely affected Stage 3 population.

**Laboratory-based CKD.** The AKI paradox was replicated: CKD-negative Stage 3 mortality was 44.4% versus CKD-positive 25.8% (ratio 1.7×). In Model B, lab-CKD was associated with increased mortality (aOR 1.385, p < 0.001), in contrast to ICD-CKD (aOR 0.620, p < 0.001). This reversal suggests that CKD defined by laboratory criteria captures a different clinical phenotype than CKD defined by administrative codes. This finding is consistent with Høyer et al. [30], who validated ICD-10 AKI codes against laboratory-defined episodes in 947,209 cases and found a sensitivity of only 7.5% with high PPV (90.6%), confirming that administrative and laboratory criteria capture systematically different populations. Similarly, Esposito et al. [32] reported a 68% AKI undetection rate in 56,820 hospitalized patients, with 16.7% having creatinine-defined AKI never coded, supporting the notion that ICD-based and laboratory-based definitions identify fundamentally different patient populations.

**Complete-case SOFA sensitivity.** Restricting to patients with both GCS and bilirubin (N = 17,285), KDIGO ORs were: Stage 1: 1.610 (1.453–1.785), Stage 2: 2.320 (1.978–2.721), Stage 3: 2.170 (1.922–2.449)—11–23% lower than full-cohort estimates. This confirms that SOFA imputation systematically underestimates illness severity, thereby overestimating the independent KDIGO effect.

**Temporal trend analysis (2008–2022).** AKI incidence increased from 28.4% (2008–2010) to 33.1% (2020–2022), and mortality from 9.8% to 12.9%. AKI-specific mortality rose from 19.1% to 26.2%. The 2020–2022 period showed the largest changes, likely reflecting COVID-19–related shifts in ICU case mix (Table S4).

### External Validation Cohort (eICU-CRD)

#### Study Population and KDIGO Distribution

The eICU-CRD validation cohort comprised 166,373 adult ICU stays. AKI (Stage 1–3) occurred in 52,984 patients (31.8%), consistent with the MIMIC-IV derivation cohort (29.9%). The distribution was: Stage 0: 113,389 (68.2%); Stage 1: 28,129 (16.9%); Stage 2: 6,523 (3.9%); Stage 3: 18,332 (11.0%).

#### Mortality by KDIGO Stage

In-hospital mortality by stage was: Stage 0: 4.6%; Stage 1: 12.5%; Stage 2: 18.0%; Stage 3: 19.2%; overall: 8.1%. The mortality gradient and Stage 2–3 convergence pattern were consistent with the derivation cohort (Table 3).

#### Logistic Regression: Replication of Main Findings

Model A (KDIGO + age + sex + comorbidities): Stage 1 OR 2.995 (95% CI 2.863–3.133), Stage 2 OR 4.575 (4.270–4.903), Stage 3 OR 4.961 (4.737–5.195); all p < 0.001.

Model B (Model A + APACHE IVa + mechanical ventilation + CKD): Stage 1 aOR 1.849 (95% CI 1.761–1.942), Stage 2 aOR 2.811 (2.603–3.035), Stage 3 aOR 2.581 (2.450–2.719); all p < 0.001.

The direction and magnitude of severity adjustment attenuation closely paralleled the derivation cohort: Stage 1 attenuated by 38.3% (MIMIC: 31.8%), Stage 2 by 38.6% (MIMIC: 43.4%), Stage 3 by 48.0% (MIMIC: 59.0%). Model A AUC improved from 0.672 to Model B AUC of 0.818, consistent with the MIMIC-IV pattern (0.778 → 0.828).

#### AKI Paradox Replication

The CKD × KDIGO Stage 3 interaction was formally replicated in eICU-CRD: interaction OR = 0.726 (95% CI 0.603–0.872, p < 0.001), confirming that the paradox is not a single-centre artefact but a reproducible phenomenon across independent databases, different hospital systems, and different severity scoring methodologies. The interaction was again Stage 3-specific, with no significant interaction at Stages 1 or 2 (Stage 1: OR 0.959, p = 0.679; Stage 2: OR 1.243, p = 0.221).

#### Cross-Cohort Comparison

| Metric | MIMIC-IV (N = 84,167) | eICU-CRD (N = 166,373) |
|:---|:---:|:---:|
| AKI prevalence | 29.9% | 31.8% |
| Stage 3 proportion | 11.4% | 11.0% |
| Overall mortality | 10.7% | 8.1% |
| Stage 3 mortality | 29.0% | 19.2% |
| Model B Stage 1 aOR | 1.885 (1.765–2.013) | 1.849 (1.761–1.942) |
| Model B Stage 2 aOR | 2.978 (2.688–3.300) | 2.811 (2.603–3.035) |
| Model B Stage 3 aOR | 2.881 (2.664–3.115) | 2.581 (2.450–2.719) |
| Model A AUC | 0.778 | 0.672 |
| Model B AUC | 0.828 | 0.818 |
| CKD × Stage 3 interaction OR | 0.476 (<0.001) | 0.726 (<0.001) |

The consistent direction and comparable magnitude of all key findings across two independent databases (total N = 250,540) substantially strengthen the generalisability of our conclusions.

### Time-Dependent ROC Analysis

Time-dependent ROC analysis (Figure S6) revealed that the full Model B (incorporating KDIGO stage, SOFA, ventilation, vasopressors, CKD, and comorbidities) maintained stable discriminative performance across time horizons: 7-day AUC = 0.643, 14-day AUC = 0.639, 28-day AUC = 0.650. In contrast, KDIGO stage alone showed declining and near-random performance over time: 7-day AUC = 0.484, 14-day AUC = 0.452, 28-day AUC = 0.404. This temporal divergence confirms that KDIGO staging captures predominantly early mortality risk, while the full severity-adjusted model provides sustained prognostic value throughout the 28-day window. The declining KDIGO-only AUC is consistent with the time-stratified Cox finding that Stage 3 hazard increases over time—the staging system's static categories cannot capture this evolving risk profile.

### Nomogram and Web Calculator

A clinical nomogram was developed based on the eICU-CRD Model B (Figure S7), incorporating five routinely available variables: KDIGO stage, APACHE score, mechanical ventilation, age, and sex. The nomogram demonstrated strong discrimination (AUC = 0.818) and good calibration (Brier score = 0.062), with the calibration curve showing close agreement between predicted and observed mortality across all risk deciles.

Variable importance, expressed as nomogram point ranges: APACHE score (0–100 points, reflecting the dominant prognostic contribution), KDIGO Stage 3 (0–48 points), mechanical ventilation (0–30 points), age (0–18 points), and sex (0–1 point).

A standalone web-based risk calculator was deployed (Supplementary Material), enabling real-time mortality risk estimation by clinicians at the bedside by entering the five nomogram variables. The calculator embeds the logistic regression coefficients (intercept = −6.225; KDIGO 1 coefficient = 0.615; KDIGO 2 = 1.033; KDIGO 3 = 0.948; APACHE = 0.033 per point; mechanical ventilation = 0.584; age = 0.018 per year; male sex = −0.022) and computes predicted mortality probability in real time.

---

## Discussion

In this large dual-cohort retrospective study of 250,540 ICU stays from two independent databases (MIMIC-IV: 84,167; eICU-CRD: 166,373), we comprehensively evaluated the association between KDIGO creatinine-based AKI staging and mortality. To our knowledge, this is the first study to simultaneously: (1) formally quantify illness severity confounding using incremental regression with NRI/IDI reclassification metrics in a cohort exceeding 80,000 patients; (2) decompose the Stage 3 AKI paradox by RRT indication with formal interaction testing; (3) apply time-stratified Cox models revealing evolving hazard profiles; (4) externally validate all primary findings—including the AKI paradox—in an independent multi-centre database; (5) assess time-dependent predictive performance across 7/14/28-day horizons; and (6) translate the findings into a deployable clinical nomogram with web-based calculator. The results collectively advance our understanding of KDIGO staging in several clinically important respects.

Our mortality rates (MIMIC Stage 0: 6.0%; Stage 1: 15.1%; Stage 2: 26.7%; Stage 3: 29.0%; eICU: 4.6%, 12.5%, 18.0%, 19.2%) are broadly consistent with the graded association reported in prior large cohorts, though absolute rates vary with case mix and AKI definition. The landmark AKI-EPI multinational study reported mortality ranging from 16.1% (no AKI) to 50.0% (Stage 3) using both creatinine and urine output criteria [2], while studies using creatinine-only criteria—such as ours—generally report lower mortality at each stage [19]. Hoste et al. [33] reported RIFLE Stage 1–3 mortality of 8.8%, 11.4%, and 26.3% in a single-centre ICU cohort of 5,383 patients, while Nisula et al. [34] reported KDIGO Stage 1–3 mortality of 29.3%, 34.1%, and 39.0% in a prospective multicentre cohort of 2,901 patients. The lower absolute mortality in eICU-CRD compared to MIMIC-IV likely reflects differences in case mix (eICU includes community hospitals), shorter median ICU stays, and the multi-centre nature of the database diluting the high-acuity BIDMC case mix.

**First**, AKI severity demonstrated a robust dose-response relationship with mortality in both cohorts, but illness severity accounted for a substantial portion of this association, with severity adjustment substantially improving model discrimination (MIMIC AUC: 0.778→0.828; eICU: 0.672→0.818; NRI = 0.625). These quantitative results provide direct empirical evidence addressing the causal-versus-correlational question raised by Girling et al. [24], who applied the Bradford Hill criteria to AKI-associated mortality and identified severity confounding as the principal obstacle to causal interpretation. Landmark analysis excluding early deaths revealed that the full-cohort ORs were conservative (12–29% lower than landmark estimates), indicating immortal time bias biased our primary analysis toward the null.

**Second**, RCS analysis confirmed a significant nonlinear creatinine ratio–mortality relationship, with the OR peaking at a creatinine ratio of ~2.6 and plateauing at higher ratios (Figure S1). This is consistent with the Bonferroni-confirmed non-significant mortality difference between Stage 2 and Stage 3 and the RMST convergence. The plateau at ratios >2.6 supports a biological ceiling effect. This within-stage heterogeneity is further illustrated by Long et al. [31], who showed in 47,333 surgical patients that KDIGO's absolute (≥26.5 µmol/L) and relative (≥1.5× baseline) Stage 1 criteria identify different populations with different long-term outcomes stratified by CKD status—a finding that parallels our Stage 3 decomposition by mechanism of entry (RRT-only 8.3%, Cr-only 58.6%, both 33.1%) and supports the argument that criterion type, not just stage number, drives prognosis.

**Third**, time-stratified Cox analysis revealed that Stage 3's hazard progressively increased over the follow-up period (HR 1.41→1.61→1.86), explaining the PH violation and the convergence with Stage 2: Stage 3's mortality signal emerges later and persists longer, while Stage 2's effect is more front-loaded. This temporal pattern was further corroborated by the time-dependent ROC analysis: KDIGO stage alone showed declining discriminative performance over time (7-day AUC 0.484 → 28-day AUC 0.404), while the full severity-adjusted model maintained stable performance (0.643 → 0.650). The declining KDIGO-only AUC confirms that static staging categories cannot capture the evolving risk profile, supporting the need for dynamic risk assessment tools. Lin et al. [25] demonstrated that trajectory-based AKI subphenotyping in the same two databases identified four subgroups with mortality ranging from 15% to 24%, capturing within-stage heterogeneity that static KDIGO categories cannot resolve. Our finding that KDIGO-only AUC declines over time while the full model maintains performance is consistent with their observation that temporal dynamics add prognostic information beyond a single staging timepoint.

**Fourth**, the AKI paradox was formally confirmed as a Stage 3-specific phenomenon and—critically—was **replicated in the independent eICU-CRD cohort** (interaction OR = 0.726, p < 0.001, vs. MIMIC OR = 0.476, p < 0.001). This cross-database replication, across different hospital systems (single-centre academic vs. 208-hospital multi-centre), different severity scores (SOFA vs. APACHE IVa), and different data collection platforms, provides the strongest available evidence that the AKI paradox is a biological phenomenon rather than a database-specific artefact. RRT-stratified decomposition revealed the paradox persisted in the pure creatinine-criteria subgroup, indicating that differential RRT initiation is not the sole explanation. A biological mechanism is suggested by Gleeson et al. [29], who demonstrated in 157 ICU patients with AKI requiring RRT that higher serum creatinine at RRT initiation was paradoxically associated with better ICU survival (OR 0.33, 95% CI 0.17–0.62), and that muscle mass—a marker of frailty—interacted with creatinine and superseded it as a survival predictor (OR 0.26, P = .02). Creatinine is a muscle-derived marker, so elevated creatinine in non-CKD patients with preserved muscle mass may reflect less catastrophic injury than the same absolute creatinine level in sarcopenic patients who reach Stage 3 at deceptively low values. This mechanism could explain both our AKI paradox (non-CKD Stage 3 mortality 49.4% vs. CKD 20.5%) and the Stage 2 low-baseline-creatinine anomaly (mean baseline SCr 0.90 mg/dL, highest eGFR 84.3), and is consistent with the RCS plateau above a ratio of ~2.6. This finding is consistent with recent evidence that KDIGO criteria may be inadequate for staging AKI in patients with preexisting kidney dysfunction [18]. Furthermore, Yasrebi-de Kom et al. [20] demonstrated that KDIGO staging errors track with CKD status—misclassified Stage 3 cases are CKD-enriched—a mechanism directly relevant to the paradox we decompose. The heterogeneity within KDIGO stages has been further demonstrated by Dong et al., who found that subdividing Stage 1 based on absolute versus relative SCr increase did not differentiate ICU mortality [21].

**Fifth**, sepsis status emerged as a second clinically important effect modifier, with ORs 1.3–1.8-fold higher in sepsis patients at Stages 1 and 2. This large effect modification suggests that AKI in the context of sepsis carries a distinct pathophysiological signature that amplifies mortality risk [5,15].

**Sixth**, the nomogram and web calculator developed from the eICU-CRD model provide a practical clinical translation of our findings. The five-variable model (KDIGO stage, APACHE score, mechanical ventilation, age, sex) achieved an AUC of 0.818 with good calibration (Brier = 0.062), demonstrating that mortality risk can be accurately estimated using routinely available ICU variables. The choice of logistic regression over machine learning was deliberate: while recent CKJ studies have deployed ensemble ML models for ICU-AKI mortality prediction with SHAP interpretability [26], and traditional logistic regression has been shown to remain competitive with ML in external validation [27], our nomogram prioritises bedside deployability—transportable coefficients, no computational infrastructure required, and real-time calculation via a standalone web tool—over marginal AUC gains from black-box approaches. Our AUC of 0.818 is comparable to the multi-database ML model by Luo et al. [26] (AUC 0.732–0.778 across external validation cohorts) while requiring only five routinely available variables. Notably, our time-dependent AUCs (0.643–0.650) measure 7–28-day mortality prediction—a different and arguably more clinically relevant target than the in-hospital mortality AUC of 0.818. The web calculator enables real-time bedside application, bridging the gap between research findings and clinical practice.

### What Distinguishes This Study

This study addresses four gaps in the existing literature through several distinguishing features: (1) **Dual-cohort validation**: Unlike single-database studies, we validated all primary findings—including the AKI paradox—in an independent multi-centre cohort (eICU-CRD), establishing cross-database reproducibility (total N = 250,540). (2) **Comprehensive confounding quantification**: We used NRI/IDI reclassification metrics, DCA, and E-values to rigorously quantify the independent contribution of KDIGO staging beyond illness severity. (3) **Temporal dynamics**: Time-stratified Cox analysis and time-dependent ROC jointly characterised how KDIGO-mortality associations evolve over the follow-up period. (4) **Clinical translation**: A nomogram with web calculator provides immediate clinical applicability, moving beyond statistical association to practical risk stratification.

### Strengths and Limitations

This study has several strengths: the large combined sample size (N = 250,540 across two databases) enabling multi-level stratified analyses; external validation in an independent multi-centre cohort; rigorous three-model comparison with formal NRI/IDI reclassification metrics and DCA; conservative Bonferroni correction; 24-hour landmark analysis; time-stratified Cox to address PH violations; time-dependent ROC analysis; formal CKD × KDIGO interaction tests confirming Stage 3 specificity and cross-database replication; a clinically deployable nomogram with web calculator; E-value sensitivity analysis for unmeasured confounding; and systematic sensitivity analyses. We transparently assessed two offsetting sources of bias: immortal time bias (attenuating ORs) and SOFA imputation bias (inflating ORs by 11–23%), suggesting the primary estimates are reasonable central values.

Several limitations must be acknowledged. First, AKI staging was based exclusively on creatinine criteria without urine output, which may underestimate AKI incidence by 25–35% [16]. Singh et al. [28] demonstrated that oliguric AKI carries distinctively worse outcomes compared with non-oliguric AKI, and Lin et al. [25] incorporated urine output trajectories in their subphenotype analysis—phenotypes our creatinine-only definition cannot resolve. Second, baseline SCr was defined as the first hospital measurement, which may overestimate true baseline in acute-on-chronic presentations; the MDRD sensitivity analysis partially addresses this. Third, the derivation cohort is from a single centre (BIDMC), though this is mitigated by the multi-centre eICU-CRD validation. Fourth, competing risks were not formally modeled, and the observational design precludes causal inference. Fifth, ICD-based comorbidity ascertainment is subject to coding bias; Høyer et al. [30] demonstrated that ICD-10 AKI codes have a sensitivity of only 7.5% despite high PPV (90.6%), and Esposito et al. [32] reported a 68% AKI undetection rate, supporting our decision to run both ICD- and laboratory-based CKD definitions in parallel. Sixth, SOFA imputation remains a concern in the MIMIC cohort; the complete-case analysis confirmed direction and magnitude of bias (11–23% OR overestimation). Seventh, the eICU-CRD uses APACHE IVa rather than SOFA, meaning the severity scores are not directly interchangeable; however, both are validated multi-organ dysfunction scores and the consistent results across both support robustness to severity scoring methodology. Eighth, serum creatinine is confounded by muscle mass. Gleeson et al. [29] demonstrated that muscle mass superseded creatinine as a survival predictor in ICU AKI and that higher creatinine at RRT initiation was paradoxically associated with better survival. Neither MIMIC-IV nor eICU-CRD provides muscle mass or sarcopenia measures, precluding adjustment for this confounder, which may partly explain the AKI paradox and the Stage 2 low-baseline-creatinine anomaly.

**Collider bias and selection effects.** Conditioning on ICU admission and on serum creatinine measurement availability may introduce collider bias. The DAG in Supplementary Figure S5 formally identifies these variables as colliders and confirms they were not included in the regression models. The E-values (≥2.77 for all stage-specific associations) further suggest that any unmeasured confounder would need to be quite strong to fully explain away the observed associations.

### Clinical Implications

Our findings have several practical implications. First, AKI staging should be interpreted in the context of overall illness severity; a Stage 3 AKI in a patient with low SOFA carries a different prognosis than the same stage in multi-organ failure. Second, the Stage 2–3 convergence and temporal dynamics suggest that incorporating CKD status and the mechanism of Stage 3 entry (RRT vs. creatinine criteria) could refine risk stratification. Third, non-CKD patients developing Stage 3 AKI represent an extremely high-risk phenotype (49.4% mortality), meriting particularly aggressive management. Fourth, the significant sepsis × KDIGO interaction supports the use of sepsis-specific AKI risk models. Fifth, the nomogram and web calculator provide a practical tool for real-time mortality risk estimation at the bedside, using five routinely available variables.

### Conclusion

In this large dual-cohort study of 250,540 ICU stays from two independent databases, KDIGO creatinine-based AKI staging demonstrated a robust, independent association with in-hospital mortality. Illness severity accounted for a substantial portion of this association, and the Stage 2–3 convergence and AKI paradox—formally validated as a Stage 3-specific interaction—were reproducible across both cohorts, supporting their biological rather than artefactual nature. Time-dependent ROC analysis revealed that KDIGO stage alone has declining predictive value over time, while the full severity-adjusted model maintains stable performance. A clinically deployable nomogram and web calculator provide practical risk stratification. These findings support the continued use of KDIGO staging for prognosis in the ICU while advocating for CKD-aware, sepsis-aware, and temporally nuanced refinement of AKI risk stratification.

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
| Vasopressor use | — | — | 0.660 (0.614–0.709) | <0.001 | 0.688 (641–0.739) | <0.001 |
| CKD (ICD) | — | — | 0.620 (0.582–0.660) | <0.001 | 0.626 (0.588–0.667) | <0.001 |
| Sepsis | 1.805 (1.713–1.902) | <0.001 | 1.598 (1.512–1.690) | <0.001 | 1.671 (1.583–1.764) | <0.001 |
| Metastatic Cancer | 1.547 (1.415–1.691) | <0.001 | 1.947 (1.773–2.139) | <0.001 | 1.601 (1.464–1.751) | <0.001 |
| Liver Disease | 1.590 (1.496–1.691) | <0.001 | 1.164 (1.089–1.244) | <0.001 | 1.220 (1.144–1.301) | <0.001 |
| **AIC** | 48,347 | | 44,716 | | 45,132 | |
| **AUC (95% CI)** | 0.778 (0.773–0.784) | | 0.828 (0.823–0.833) | | 0.825 (0.820–0.831) | |
| **Brier Score** | 0.0829 | | 0.0779 | | 0.0788 | |

*Reference: KDIGO Stage 0. Full covariate results (27 variables) in Table S2. VIF for all covariates <5 (SOFA 4.73, ventilation 2.68, vasopressor 2.19).*

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
| **Model B Cox HRs** | | |
| Stage 1 | 1.174 (1.124–1.227) | 1.624 |
| Stage 2 | 1.352 (1.267–1.442) | 2.255 |
| Stage 3 | 1.253 (1.196–1.313) | 2.037 |
| **CKD × Stage 3 Interaction** | | |
| Interaction OR | 0.726 (0.603–0.872) | 0.476 |
| p-value | <0.001 | <0.001 |
| **Model Discrimination** | | |
| Model A AUC | 0.672 | 0.778 |
| Model B AUC | 0.818 | 0.828 |

---

## Figure Legends

- **Figure 1.** CONSORT-style patient selection flowchart. From 94,458 total ICU stays, 84,167 were included after excluding non-adult, duplicate, and 1,075 stays without serum creatinine data.
- **Figure 2.** In-hospital mortality rate by KDIGO AKI stage (bar chart with 95% Wilson CI). Stage 0: n = 59,033; Stage 1: n = 12,817; Stage 2: n = 2,729; Stage 3: n = 9,588.
- **Figure 3.** Kaplan-Meier 30-day survival curves stratified by KDIGO stage, with number-at-risk table and Bonferroni-corrected pairwise log-rank test results. Overall log-rank χ² = 4,737.29, p < 0.001.
- **Figure 4.** Forest plot of adjusted odds ratios from Model B (full 27 covariate set), showing all variables with 95% CI. N = 84,167.
- **Figure 5.** ICU length of stay by KDIGO stage (violin plot; top 5% of each stage truncated for visual clarity).
- **Figure 6.** Dose-response curve: creatinine ratio versus mortality rate with 95% CI. Total n = 84,167 across 7 bins.
- **Figure 7.** Distribution of KDIGO AKI stages (bar chart).
- **Figure 8.** AKI paradox. Panel A: in-hospital mortality by CKD status within each KDIGO stage. Panel B: Stage 3 decomposed by RRT indication within each CKD stratum.

## Supplementary Material

### Supplementary Figures

- **Figure S1.** Restricted cubic spline (RCS) analysis of the nonlinear creatinine ratio–mortality relationship (4 knots, 10th/35th/65th/90th percentiles, shaded area = 95% CI). Nonlinearity test: χ² = 161.29, df = 2, p < 0.001. N = 84,167.
- **Figure S2.** Model discrimination and clinical utility. Panel A: ROC curves for Models A (AUC = 0.778), B (AUC = 0.828), and C (AUC = 0.825). Panel B: DCA showing net benefit of Model B.
- **Figure S3.** Sepsis subgroup analysis. Forest plot comparing adjusted odds ratios (Model B) for KDIGO stages between sepsis and non-sepsis subgroups, with interaction test results.
- **Figure S4.** Calibration plot for Model B. Observed versus predicted mortality across 10 deciles, with Hosmer-Lemeshow test.
- **Figure S5.** Directed acyclic graph (DAG) representing the assumed causal structure. Measured confounders (age, sex, comorbidities, SOFA, ventilation, vasopressors, CKD) are shown as ancestors of both exposure and outcome. Colliders (ICU admission, SCr measurement availability) are identified and confirmed as not adjusted in the models.
- **Figure S6.** Time-dependent ROC analysis. Panel A: ROC curves at 7, 14, and 28 days for the full Model B (AUC = 0.643, 0.639, 0.650). Panel B: Comparison of KDIGO stage alone (AUC = 0.484, 0.452, 0.404) versus full model, demonstrating declining KDIGO-only performance over time.
- **Figure S7.** Nomogram for ICU mortality prediction. Panel A: Nomogram points for KDIGO stage, APACHE score, mechanical ventilation, age, and sex. Panel B: Calibration curve (10 quantiles, Brier score = 0.062). Panel C: Variable importance forest plot. Based on eICU-CRD Model B (N = 166,373, AUC = 0.818).

### Supplementary Tables

- **Table S1.** Systematic missing data report stratified by KDIGO stage for all key variables (N = 84,167).
- **Table S2.** Full logistic regression results for Model B (all 27 covariates with OR, 95% CI, *p*-value).
- **Table S3.** Time-stratified Cox proportional hazards regression (0–7, 7–14, 14–30 days).
- **Table S4.** Temporal trend analysis (2008–2022) by 3-year periods.
- **Table S5.** 24-hour landmark sensitivity analysis: full-cohort versus landmark Model B ORs.
- **Table S6.** Dose-response analysis: creatinine ratio binned mortality rates.
- **Table S7.** AKI paradox: CKD-stratified mortality, RRT-stratified Stage 3 decomposition, and formal CKD × KDIGO interaction test.
- **Table S8.** E-value sensitivity analysis for Model B associations.
- **Supplementary File.** Web-based mortality risk calculator (standalone HTML file).
- **STROBE Checklist.** Completed STROBE checklist for observational cohort studies.

---

## Abbreviations

| Abbreviation | Definition |
|:---|:---|
| AKI | acute kidney injury |
| aOR | adjusted odds ratio |
| APACHE | Acute Physiology and Chronic Health Evaluation |
| AUC | area under the receiver operating characteristic curve |
| BIDMC | Beth Israel Deaconess Medical Center |
| BUN | blood urea nitrogen |
| CI | confidence interval |
| CKD | chronic kidney disease |
| CKD-EPI | Chronic Kidney Disease Epidemiology Collaboration |
| CRRT | continuous renal replacement therapy |
| CVVHD | continuous veno-venous haemodialysis |
| CVVHDF | continuous veno-venous haemodiafiltration |
| DCA | decision curve analysis |
| DAG | directed acyclic graph |
| eGFR | estimated glomerular filtration rate |
| GCS | Glasgow Coma Scale |
| HIPAA | Health Insurance Portability and Accountability Act |
| HR | hazard ratio |
| ICD | International Classification of Diseases |
| ICU | intensive care unit |
| IDI | integrated discrimination improvement |
| KDIGO | Kidney Disease: Improving Global Outcomes |
| LOS | length of stay |
| MAP | mean arterial pressure |
| MDRD | Modification of Diet in Renal Disease |
| MIMIC-IV | Medical Information Mart for Intensive Care IV |
| NRI | net reclassification improvement |
| OR | odds ratio |
| PH | proportional hazards |
| RCS | restricted cubic splines |
| RMST | restricted mean survival time |
| RRT | renal replacement therapy |
| SCr | serum creatinine |
| SCUF | slow continuous ultrafiltration |
| SOFA | Sequential Organ Failure Assessment |
| SpO₂ | peripheral oxygen saturation |
| STROBE | Strengthening the Reporting of Observational Studies in Epidemiology |
| VIF | variance inflation factor |
| WBC | white blood cell count |

---

## References

1. Susantitaphong P, Cruz DN, Cerda J, Abulfaraj M, Zaghlool M, Jha V, et al. World incidence of AKI: a meta-analysis. Clin J Am Soc Nephrol. 2013;8(9):1482–1493. doi:10.2215/CJN.00710113

2. Hoste EAJ, Bagshaw SM, Bellomo R, Cely CM, Colman R, Cruz DN, et al. Epidemiology of acute kidney injury in critically ill patients: the multinational AKI-EPI study. Intensive Care Med. 2015;41(8):1411–1423. doi:10.1007/s00134-015-3934-7

3. Kidney Disease: Improving Global Outcomes (KDIGO) Acute Kidney Injury Work Group. KDIGO Clinical Practice Guideline for Acute Kidney Injury. Kidney Int Suppl. 2012;2(1):1–138. doi:10.1038/kisup.2012.1

4. Lassnigg A, Schmidlin D, Mouhieddine M, Bachmann LM, Druml W, Bauer P, et al. Minimal changes of serum creatinine predict prognosis in patients after cardiothoracic surgery: a prospective cohort study. J Am Soc Nephrol. 2004;15(6):1597–1605. doi:10.1097/01.ASN.0000130340.93930.DD

5. Bagshaw SM, Uchino S, Bellomo R, Morimatsu H, Morgera S, Schetz M, et al. Septic acute kidney injury in critically ill patients: clinical characteristics and outcomes. Clin J Am Soc Nephrol. 2007;2(3):431–439. doi:10.2215/CJN.03681106

6. Uchino S, Kellum JA, Bellomo R, Doig GS, Morimatsu H, Morgera S, et al. Acute renal failure in critically ill patients: a multinational, multicenter study. JAMA. 2005;294(7):813–818. doi:10.1001/jama.294.7.813

7. Luo X, Jiang L, Du B, Wen Y, Wang M, Xi X; Beijing Acute Kidney Injury Trial Workgroup. A comparison of different diagnostic criteria of acute kidney injury in critically ill patients. Crit Care. 2014;18(4):R144. doi:10.1186/cc13977

8. Thomas ME, Blaine C, Dawnay A, Devonald MAJ, Ftouh S, Laing C, et al. The definition of acute kidney injury and its use in practice. Kidney Int. 2015;87(1):62–73. doi:10.1038/ki.2014.328

9. Pannu N, James M, Hemmelgarn B, Klarenbach S; Alberta Kidney Disease Network. Association between AKI, recovery of renal function, and long-term outcomes after hospital discharge. Clin J Am Soc Nephrol. 2013;8(2):194–202. doi:10.2215/CJN.06480612

10. James MT, Ghali WA, Knudtson ML, Ravani P, Tonelli M, Faris P, et al. Associations between acute kidney injury and cardiovascular and renal outcomes after coronary angiography. Circulation. 2011;123(4):409–416. doi:10.1161/CIRCULATIONAHA.110.970160

11. Chertow GM, Burdick E, Honour M, Bonventre JV, Bates DW. Acute kidney injury, mortality, length of stay, and costs in hospitalized patients. J Am Soc Nephrol. 2005;16(11):3365–3370. doi:10.1681/ASN.2004090740

12. Johnson AEW, Bulgarelli L, Shen L, Gayles A, Shammout A, Horng S, et al. MIMIC-IV, a freely accessible electronic health record dataset. Sci Data. 2023;10(1):1. doi:10.1038/s41597-022-01899-x

13. Elixhauser A, Steiner C, Harris DR, Coffey RM. Comorbidity measures for use with administrative data. Med Care. 1998;36(1):8–27. doi:10.1097/00005650-199801000-00004

14. Quan H, Sundararajan V, Halfon P, Fong A, Burnand B, Luthi JC, et al. Coding algorithms for defining comorbidities in ICD-9-CM and ICD-10 administrative data. Med Care. 2005;43(11):1130–1139. doi:10.1097/01.mlr.0000182534.19832.83

15. Prowle JR, Echeverri JE, Ligabo EV, Ronco C, Bellomo R. Fluid balance and acute kidney injury. Nat Rev Nephrol. 2010;6(2):107–115. doi:10.1038/nrneph.2009.213

16. Kellum JA, Sileanu FE, Murugan R, Lucko N, Shaw AD, Clermont G. Classifying AKI by urine output versus serum creatinine level. J Am Soc Nephrol. 2015;26(9):2231–2238. doi:10.1681/ASN.2014070724

17. Harrell FE Jr, Lee KL, Mark DB. Multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors. Stat Med. 1996;15(4):361–387. doi:10.1002/(SICI)1097-0258(19960229)15:4<361::AID-SIM168>3.0.CO;2-4

18. Hong D, Ren Q, Zhang J, Dong F, Chen S, Dong W, et al. A new criteria for acute on preexisting kidney dysfunction in critically ill patients. Ren Fail. 2023;45(1):2173498. doi:10.1080/0886022X.2023.2173498

19. Birkelo BC, Pannu N, Siew ED. Overview of diagnostic criteria and epidemiology of acute kidney injury and acute kidney disease in the critically ill patient. Clin J Am Soc Nephrol. 2022;17(5):717–735. doi:10.2215/CJN.14181021

20. Yasrebi-de Kom IAR, Dongelmans DA, Abu-Hanna A, et al. Incorrect application of the KDIGO acute kidney injury staging criteria. Clin Kidney J. 2022;15(5):937–941. doi:10.1093/ckj/sfab256

21. Dong GY, Qin JP, An Y, et al. Utilizing reclassification to explore characteristics and prognosis of KDIGO SCr AKI subgroups: a retrospective analysis of a multicenter prospective cohort study. Ren Fail. 2021;43(1):1569–1576. doi:10.1080/0886022X.2021.1997761

22. VanderWeele TJ, Ding P. Sensitivity Analysis in Observational Research: Introducing the E-Value. Ann Intern Med. 2017;167(4):268–274. doi:10.7326/M16-2607

23. Pollard TJ, Johnson AEW, Raffa JD, Celi LA, Mark RG, Badawi O. The eICU Collaborative Research Database, a freely available multi-center database for critical care research. Sci Data. 2018;5:180178. doi:10.1038/sdata.2018.178

24. Girling BJ, Channon SW, Haines RW, Prowle JR. Acute kidney injury and adverse outcomes of critical illness: correlation or causation? Clin Kidney J. 2020;13(2):133–141. doi:10.1093/ckj/sfz158

25. Lin J, Liu L, Zhu S, Gao J, Liu L, Hong H, Wei Y, Yang J, Liu X, Li R, Zhu J. Machine learning-derived multivariate renal function trajectories in acute kidney injury in critically ill patients: a multicentre retrospective study. Clin Kidney J. 2025;18(6):sfaf142. doi:10.1093/ckj/sfaf142

26. Luo S, Lai J, Mo L, Shen X, Fang R. Prediction of hospital mortality in sepsis-associated acute kidney injury using a machine-learning approach: a multicenter study using SHAP interpretability analysis. Clin Kidney J. 2026;19(1):sfaf372. doi:10.1093/ckj/sfaf372

27. Liu Y, Zhu X, Xue J, Maimaitituerxun R, Chen W, Dai W. Machine learning models for mortality prediction in critically ill patients with acute pancreatitis-associated acute kidney injury. Clin Kidney J. 2024;17(10):sfae284. doi:10.1093/ckj/sfae284

28. Singh S, Andonovic M, Traynor JP, Shaw MF, Sim MAB, Mark PB, Puxty KA. Short- and long-term outcomes in oliguric and non-oliguric acute kidney injury in intensive care: a retrospective, post hoc, bicentric study. Clin Kidney J. 2025;18(6):sfaf170. doi:10.1093/ckj/sfaf170

29. Gleeson PJ, Crippa IA, Sannier A, Koopmansch C, Bienfait L, Allard J, Sexton DJ, Fontana V, Rorive S, Vincent JL, Creteur J, Taccone FS. Critically ill patients with acute kidney injury: clinical determinants and post-mortem histology. Clin Kidney J. 2023;16(10):1664–1673. doi:10.1093/ckj/sfad113

30. Høyer S, Heide-Jørgensen U, Jensen SK, Pottegård A, Christiansen CF. Sensitivity and positive predictive value of diagnosis codes for acute kidney injury in Denmark. Clin Kidney J. 2026;19(3):sfag019. doi:10.1093/ckj/sfag019

31. Long TE, Helgason D, Helgadottir S, Sigurdsson GH, Palsson R, Sigurdsson MI, Indridason OS. Mild Stage 1 post-operative acute kidney injury: association with chronic kidney disease and long-term survival. Clin Kidney J. 2021;14(1):237–244. doi:10.1093/ckj/sfz197

32. Esposito P, Cappadona F, Marengo M, Fiorentino M, Fabbrini P, Quercia AD, Garzotto F, Castellano G, Cantaluppi V, Viazzi F. Recognition patterns of acute kidney injury in hospitalized patients. Clin Kidney J. 2024;17(8):sfae231. doi:10.1093/ckj/sfae231

33. Hoste EAJ, Clermont G, Kersten A, Venkataraman R, Angus DC, De Bacquer D, Kellum JA. RIFLE criteria for acute kidney injury are associated with hospital mortality in critically ill patients: a cohort analysis. Crit Care. 2006;10(3):R73. doi:10.1186/cc4915

34. Nisula S, Kaukonen KM, Vaara ST, Korhonen AM, Poukkanen M, Karlsson S, Haapio M, Inkinen O, Parviainen I, Suojaranta-Ylinen R, Laurila JJ, Tenhunen J, Reinikainen M, Ala-Kokko T, Ruokonen E, Kuitunen A, Pettilä V; FINNAKI Study Group. Incidence, risk factors and 90-day mortality of patients with acute kidney injury in Finnish intensive care units: the FINNAKI study. Intensive Care Med. 2013;39(3):420–428. doi:10.1007/s00134-012-2796-5

---

## Funding

This work was supported by the Chronic Disease Management Research Project of National Health Commission Capacity Building and Continuing Education Center (Grant No. GWJJMB202510024181), a grant from the Changsha Science and Technology Bureau Project (Grant No. kq2014242), and the Natural Science Foundation of Hunan Province of China (Grant No. 2021JJ30959).

---

## Data Availability

MIMIC-IV is publicly available at PhysioNet (https://physionet.org/content/mimiciv/). eICU-CRD is publicly available at PhysioNet (https://physionet.org/content/eicu-crd/). All analysis code and key result files are publicly available at GitHub (https://github.com/wudengke2010/kdigo-aki-mimic-iv).

---

## Ethics Approval

Access to MIMIC-IV and eICU-CRD was obtained through PhysioNet after completion of the CITI Data or Specimens Only Research course. The institutional review board at BIDMC approved the creation of the MIMIC-IV database. Use of de-identified, publicly available data did not require additional institutional review board approval.

---

## Consent for Publication

Not applicable. This study used de-identified, publicly available data from the MIMIC-IV and eICU-CRD databases. The institutional review boards waived the requirement for informed consent.

---

## Patient and Public Involvement

Patients and/or the public were not involved in the design, conduct, reporting, or dissemination plans of this research.

---

## Author Contributions

**Jiqiang Liu:** Conceptualisation, Data curation, Formal analysis, Investigation, Methodology, Writing – original draft. **Dengke Wu:** Conceptualisation, Funding acquisition, Project administration, Supervision, Writing – review & editing.

---

## Acknowledgements

The authors thank the PhysioNet team and the Laboratory for Computational Physiology at the Massachusetts Institute of Technology for creating and maintaining the MIMIC-IV and eICU-CRD databases, which made this study possible.

---

## Disclosure of interest

The authors report there are no competing interests to declare.

---

## ORCID

Dengke Wu https://orcid.org/0009-0008-1363-9621

---

## Declaration of Generative AI in Scientific Writing

During the preparation of this work, the authors used large language models (LLMs) for language editing and copyediting. After using these tools, the authors reviewed and edited the content as needed and take full responsibility for the content of the publication.

---

*Revised manuscript (v12 — CKJ citation enhancement) building on v11 with: (i) 9 new CKJ references [24–32] integrated at contextually appropriate locations; (ii) 2 missing references added: Hoste et al. 2006 [33] (RIFLE mortality, Crit Care) and Nisula et al. 2013 [34] (FINNAKI study, Intensive Care Med), fixing the attribution error in Discussion paragraph 2; (iii) Introduction gap #1 strengthened with Girling et al. [24] (Bradford Hill criteria); (iv) Introduction gap #4 strengthened with Luo et al. [26] and Liu et al. [27] (ML prediction context); (v) Methods PH-violation section justified with Singh et al. [28] (CKJ precedent for time-stratified analysis); (vi) Discussion AKI paradox paragraph enhanced with Gleeson et al. [29] (muscle mass mechanism) and Yasrebi-de Kom et al. [20] (CKD misclassification); (vii) Discussion nomogram section justified vs ML with Luo et al. [26] and Liu et al. [27]; (viii) Discussion TD-ROC section linked to Lin et al. [25] (trajectory subphenotyping); (ix) Discussion Stage 2–3 convergence linked to Long et al. [31] (absolute vs relative criteria); (x) Results ICD-lab CKD reversal explained with Høyer et al. [30] and Esposito et al. [32]; (xi) New limitation on sarcopenia/muscle mass (Gleeson et al. [29]); (xii) Limitations on urine output (Singh [28], Lin [25]) and ICD coding (Høyer [30], Esposito [32]) strengthened. Total references: 23→34. All MIMIC-IV and eICU-CRD data analysis results unchanged.*
