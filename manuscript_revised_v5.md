# Time-Resolved KDIGO Creatinine Staging and 30-Day Mortality in Critically Ill Patients: A Dual-Cohort Study with Time-Varying Exposure Modelling, Competing-Risk Analysis, and Locked-Model External Validation (MIMIC-IV v3.1 and eICU-CRD)

---

## Authors

**Jiqiang Liu**<sup>1,2</sup>, Dengke Wu<sup>1,2*</sup>

<sup>1</sup> Department of Emergency Medicine, The Second Xiangya Hospital of Central South University, Changsha 410011, Hunan, China

<sup>2</sup> Emergency Medicine and Difficult Diseases Institute, The Second Xiangya Hospital of Central South University, Changsha 410011, Hunan, China

**\*Corresponding author:** Dengke Wu, Department of Emergency Medicine, Second Xiangya Hospital, Central South University, Changsha, Hunan, China; Emergency Medicine and Difficult Diseases Institute, The Second Xiangya Hospital of Central South University, Changsha 410011, Hunan, China. Electronic address: wudk2010@csu.edu.cn.

---

## Abstract

**Background:** The KDIGO creatinine-based AKI staging system is widely used, but many large validation studies applied the criteria without their mandated temporal windows, treated maximum stage as a fixed baseline exposure (immortal-time bias), censored live discharge as random, and adjusted for severity scores embedding renal function. The independent, time-resolved association between strictly implemented KDIGO stage and mortality, and genuine locked-model transportability across databases, remain incompletely characterised.

**Methods:** This STROBE- and TRIPOD-compliant retrospective dual-cohort study used MIMIC-IV v3.1 (derivation) and eICU-CRD v2.0 (external validation). Patient-level cohorts were staged by KDIGO creatinine criteria with explicit temporal windows (≥0.3 mg/dL within 48 h; ≥1.5×/2×/3× baseline within 7 days) against a hierarchical pre-admission-prioritised baseline. Stage entered models as a **time-varying, monotonically accumulating exposure** (counting-process Cox models, cause-specific death, live discharge censored); absolute risks used Aalen-Johansen cumulative incidence functions with live discharge as a competing event. Stage-by-period interactions, an onset+72-hour AKI-duration landmark, and a locked logistic model for in-hospital death within 30 days (TRIPOD type 3) completed the analysis.

**Results:** Among 60,506 MIMIC-IV patients (AKI 28.9%) and 113,466 eICU-CRD patients (AKI 16.1%), 5,697 and 9,373 died within 30 days. Fully adjusted Stage 3 hazard ratios were 3.35 (95% CI 3.08–3.63) and 3.33 (3.12–3.55) (cross-cohort ratio 1.00). Hazards were strongly time-dependent: Stage 3 declined from 5.41 (days 0–7) to 1.11 (days 14–30) in MIMIC-IV, replicated in eICU-CRD. The 30-day cumulative incidence of death rose from 7.1% (Stage 0 at 24 h) to 33.6% (Stage 3) in MIMIC-IV; in eICU-CRD, Stage 3 (26.9%) did not exceed Stage 2 (37.5%), reflecting RRT-driven stage assignment. AKI persisting beyond 72 hours carried roughly twice the adjusted hazard of transient AKI (HR 1.95, 1.77–2.16; 1.79, 1.62–1.97). The locked model's AUC in eICU-CRD was 0.750, with calibration slope 0.946 but systematic underestimation of absolute risk (observed-to-expected ratio 1.55).

**Conclusions:** Strictly implemented, time-varying KDIGO creatinine staging retains a graded, cross-cohort-consistent association with 30-day mortality after conservative severity adjustment; the excess hazard concentrates in the first week, AKI duration adds prognostic information beyond peak stage, and transported models require absolute-risk recalibration.

**Keywords:** acute kidney injury; KDIGO; intensive care unit; mortality; competing risks; external validation

---

## Introduction

Acute kidney injury (AKI) is one of the most common and consequential complications encountered in the intensive care unit (ICU). Epidemiological studies have consistently demonstrated that AKI affects approximately 30–60% of critically ill patients, with the landmark AKI-EPI multinational cohort reporting an incidence of 57.7% and a mortality gradient from 16.1% (no AKI) to 50.0% (Stage 3) [1,2].

The evolution of consensus definitions—from RIFLE to AKIN to the current Kidney Disease: Improving Global Outcomes (KDIGO) criteria—has standardised AKI classification and revealed a graded association between AKI severity and mortality across diverse ICU populations [19]. The KDIGO 2012 guideline established a staging system integrating serum creatinine (SCr), urine output, and renal replacement therapy (RRT) criteria [3], validated across general ICU, cardiac surgery, and sepsis populations [4–7]. However, the methodological quality of large database studies implementing these criteria varies substantially, and several threats to validity recur.

**First**, the KDIGO creatinine criteria are explicitly time-bound: the absolute criterion (increase ≥0.3 mg/dL) applies within a 48-hour window, and ratio criteria (≥1.5× baseline) apply within a 7-day window. Many retrospective studies stage AKI using the peak-to-baseline creatinine ratio over the entire admission, ignoring these windows entirely. Yasrebi-de Kom et al. documented widespread incorrect application of KDIGO staging criteria, including the misinterpretation of the Stage 3 absolute threshold (an increase *to* ≥4.0 mg/dL, not *of* ≥4.0 mg/dL) [20]. Window-free staging inflates apparent AKI severity by attributing slow creatinine drifts over weeks to "acute" injury, and inflates mortality associations through reverse causation (dying patients accumulate more measurements). Whether the KDIGO-mortality association survives *strict* window-based implementation at scale is not established.

**Second**, most large studies treat the maximum stage attained over the ascertainment window as a fixed, baseline exposure, assigning patients to their eventual stage from time zero. Because higher stages are, by construction, attained only after surviving long enough to accumulate creatinine rises, this design embeds immortal time and overstates hazards. AKI stage is a **time-dependent exposure**: it should enter the risk set only once the staging criteria are actually met, in a counting-process (start–stop) framework [17]. Similarly, live discharge is a competing event for in-hospital death—treating it as random censoring overestimates absolute mortality when many patients are discharged alive—yet competing-risk estimates remain rare in AKI database studies.

**Third**, illness severity confounding has not been adequately resolved. Girling et al. applied the Bradford Hill criteria to assess whether AKI is causally associated with critical illness mortality or merely a correlate of illness severity, concluding that formal confounding quantification is needed [24]. A specific methodological trap is that commonly used severity scores (SOFA, and derived APACHE components) embed a renal term that is a mathematical function of the exposure itself; adjusting for such scores overadjusts and biases KDIGO effect estimates toward the null [19]. Studies quantifying the KDIGO-mortality association while adjusting for a **non-renal** severity score in large cohorts are lacking.

**Fourth**, claims of "external validation" in dual-database studies frequently involve refitting all model coefficients in the validation cohort, which validates the *covariate structure* but not the *model*. Transportability requires applying locked derivation coefficients to the validation cohort, in the spirit of the TRIPOD framework. The precedent for MIMIC-IV/eICU-CRD derivation-validation pairing was established by Lin et al. for AKI trajectory subphenotyping [25], but locked-model transport of a KDIGO-based mortality model has not, to our knowledge, been reported.

We therefore conducted a dual-cohort study in MIMIC-IV v3.1 (derivation) and eICU-CRD v2.0 (external validation) with the following objectives: (1) implement KDIGO creatinine staging with its mandated temporal windows and a hierarchical pre-admission-referenced baseline creatinine; (2) quantify the independent association between KDIGO stage and 30-day in-hospital mortality using a **time-varying exposure** in counting-process Cox models, adjusting for a non-renal first-24-h SOFA score; (3) estimate absolute risks with cumulative incidence functions treating live discharge as a competing event, characterise hazard time-dependence through stage-by-period interactions, and examine whether AKI duration (transient vs persistent) adds prognostic information beyond peak stage; and (4) externally validate a locked derivation model in eICU-CRD without any refitting, quantifying discrimination and calibration under a fully observed outcome definition.

---

## Methods

### Study Design and Data Sources

This was a retrospective dual-cohort study. The **derivation cohort** used the Medical Information Mart for Intensive Care IV (MIMIC-IV) database, version 3.1, containing de-identified electronic health records of patients admitted to ICUs at Beth Israel Deaconess Medical Center (Boston, MA, USA) between 2008 and 2022 [12]. The **external validation cohort** used the eICU Collaborative Research Database (eICU-CRD), version 2.0, a multi-centre database of ICU stays from 208 hospitals across the United States, collected 2014–2015 [23]. The study is reported in accordance with the STROBE statement for cohort studies; the external validation component follows the TRIPOD framework (development plus independent external validation of a locked model; TRIPOD type 3).

### Study Population

In both databases, we constructed **patient-level** cohorts: for each patient, only the first ICU stay meeting all eligibility criteria was retained, eliminating within-person correlation from repeated stays. Eligibility required age ≥18 years and ICU length of stay ≥12 hours. Patients with end-stage renal disease (ESRD) or chronic dialysis were excluded, identified in MIMIC-IV from ICD-9/ICD-10 diagnosis codes (585.6, N18.6, V45.1, Z99.2; codes stored without decimal points in MIMIC-IV; restricted to diagnoses dated on or before the index admission) and in eICU-CRD from both ICD codes and structured `pastHistory` entries (chronic haemodialysis, peritoneal dialysis, or non-dialysed renal failure; dual-source identification with overlap). Patients without a usable baseline creatinine or without any serum creatinine within the staging window were excluded. Full cohort flow is shown in Figure 1.

### Exposure: KDIGO Creatinine-Based AKI Staging with Temporal Windows

AKI severity was staged according to the KDIGO 2012 **creatinine** criteria with their mandated temporal windows (urine output criteria were not applied because hourly urine records are incompletely captured in both databases):

- **Stage 1:** increase in SCr ≥0.3 mg/dL (26.5 µmol/L) within any **48-hour rolling window**, or increase to ≥1.5× baseline within **7 days** of ICU admission;
- **Stage 2:** increase to ≥2.0× and <3.0× baseline within 7 days;
- **Stage 3:** increase to ≥3.0× baseline within 7 days, **or** SCr increase to ≥4.0 mg/dL (353.6 µmol/L) with an acute change, or initiation of RRT within 7 days.

The Stage 3 absolute threshold was implemented as an increase *to* ≥4.0 mg/dL requiring concurrent acute change (ratio ≥1.5× baseline or rise ≥0.3 mg/dL within 48 h), per Yasrebi-de Kom et al. [20]. In MIMIC-IV, SCr was extracted from `labevents` (itemid 50912) with physiological range filtering (0.1–100 mg/dL); in eICU-CRD, from the `lab` table (labname = 'creatinine'). RRT was identified from dialysis procedure events (`procedureevents`; CRRT, CVVHD, CVVHDF, haemodialysis, SCUF) in MIMIC-IV and from the `treatment` table (dialysis, CRRT, or ultrafiltration entries, excluding catheter insertion and radiology procedures) in eICU-CRD, in both cases restricted to the first 7 ICU days of the index stay.

**Time-varying exposure implementation.** For each patient, KDIGO stage was modelled as a **monotonically non-decreasing step function of time since ICU admission**: at any time *t*, the patient's stage equals the highest stage whose criteria have been met using creatinine and RRT information observed up to and including *t* only. RRT initiation immediately advanced the stage to 3. Stages were never down-graded; renal recovery is a distinct process addressed separately by the phenotype analysis below. The resulting exposure process yielded 21,312 stage transitions in MIMIC-IV and 21,741 in eICU-CRD, expanded into counting-process (start–stop) intervals (116,206 and 180,829 rows, respectively) covering follow-up from ICU admission. Because a patient contributes person-time at the stage actually attained at each moment, this implementation eliminates the immortal-time bias inherent in assigning the eventual maximum stage from time zero. For descriptive stratification at a fixed time point (cumulative incidence analysis and the prediction model), we used the **stage attained at 24 hours** after ICU admission, which is determined by 24-hour data only.

**Baseline SCr hierarchy.** Baseline creatinine was defined hierarchically, prioritising pre-admission kidney function. In MIMIC-IV: (1) median outpatient SCr 7–365 days before the index admission; (2) if unavailable, median SCr from any care setting 7–365 days before admission; (3) if unavailable, median SCr 0–7 days before admission; (4) if unavailable, first SCr within 24 h of hospital admission. In eICU-CRD: (1) median SCr between hospital admission (or 7 days before ICU admission, whichever is later) and ICU admission; (2) if unavailable, first SCr within the first 24 ICU hours. The source of each patient's baseline was recorded and reported.

### Outcomes and Follow-up

The primary outcome was **in-hospital death within 30 days of ICU admission**. Event times were derived from recorded death and discharge timestamps, not from length-of-stay surrogates: in MIMIC-IV, time-to-event was `deathtime − intime` for decedents and `dischtime − intime` (censored) for survivors; 351 decedents whose recorded death time exceeded their discharge time by more than 2 hours (a known data anomaly) had event times capped at discharge time + 24 hours. In eICU-CRD, `hospitaldischargeoffset` (minutes, converted to hours) defined follow-up with `hospitaldischargestatus` defining the event. Follow-up was truncated at 30 days (720 hours): deaths occurring after 30 days were censored at day 30 (MIMIC-IV n = 293; eICU-CRD n = 177), and patients still hospitalised and alive at day 30 were administratively censored (MIMIC-IV n = 2,263, 3.7%; eICU-CRD n = 1,679, 1.5%). Live discharge before day 30 was treated as **censoring** in cause-specific Cox models (with the independent-censoring assumption stated and its implications discussed) and as a **competing event** in cumulative incidence analyses.

### Covariates

**Illness severity: non-renal first-24-h SOFA.** A first-24-hour SOFA score was computed for every patient following the mimic-code first-day SOFA specification [17], comprising five non-renal components: respiration (PaO₂/FiO₂ with mechanical ventilation), coagulation (platelets), liver (total bilirubin), cardiovascular (mean arterial pressure; norepinephrine, epinephrine, dopamine, dobutamine, and vasopressin doses with dose-tiered scoring), and central nervous system (Glasgow Coma Scale). **The renal component was deliberately excluded** to avoid adjusting for a mathematical function of the exposure; this non-renal SOFA is the primary severity adjuster in both cohorts, ensuring covariate harmonisation. In MIMIC-IV, components used charted itemids (e.g., MAP itemid 220052; GCS items 220739/223901/223900); in eICU-CRD, corresponding `lab`, `nurseCharting`, `respiratoryCare`, `infusionDrug`, and `medication` fields. Components with no data in the first 24 h were scored 0 (consistent with mimic-code defaults).

**Mechanical ventilation** was defined in MIMIC-IV as any invasive ventilation, intubation, or ventilator-care procedure event *starting within the first 24 ICU hours* of the index stay (`procedureevents`), and in eICU-CRD as ventilated status on ICU day 1 (`apachePredVar` ventday1).

**Comorbidities** (diabetes, hypertension, heart failure, COPD, liver disease) were identified in MIMIC-IV from index-admission ICD-9/ICD-10 codes using Elixhauser-style definitions (diabetes 250.x/E10–E11; hypertension 401–405.x/I10–I15; heart failure 428.x/I50.x; COPD 490–496.x/J44.x; liver disease 570–572.x/K70–K74.x) [13,14] and in eICU-CRD from the structured `pastHistory` table using equivalent diagnostic categories. **Prior chronic kidney disease (CKD)** was identified from diagnoses dated strictly *before* the index admission (ICD-9 585.3–585.5 / ICD-10 N18.3–N18.5 in MIMIC-IV; ICD codes plus pastHistory creatinine entries in eICU-CRD) to avoid conditioning on a downstream consequence of AKI. Age and sex were included as demographics. All covariates were fixed at baseline (measured within the first 24 h or before admission).

### Statistical Analysis

**Time-varying Cox models.** Three incremental cause-specific Cox models were fitted in each cohort on the counting-process data: **Model A** (time-varying KDIGO stage, reference Stage 0); **Model B** (+ age, sex); **Model C** (+ non-renal SOFA, first-24-h ventilation, five comorbidities, prior CKD). These models estimate cause-specific hazard ratios for death, censored at live discharge; they quantify the relative rate of death among patients still hospitalised and require the assumption that discharge is non-informative for the death hazard conditional on covariates.

**Competing-risk absolute risks.** Thirty-day cumulative incidence functions (CIFs) of death were estimated by the Aalen-Johansen estimator, treating live discharge as a competing event, stratified by the stage attained at 24 hours; 95% confidence intervals used 200 bootstrap replicates resampling patients. The overall 30-day CIF of live discharge is reported for context.

**Hazard time-dependence.** Because the counting-process framework makes Schoenfeld-residual diagnostics for time-varying exposures unstable at this scale, time-dependence was characterised directly through **stage-by-period interactions**: Model C was extended with interactions between each stage indicator and follow-up period (days 0–7, 7–14, 14–30), yielding period-specific hazard ratios and a joint Wald test of proportional hazards for the stage effects.

**AKI phenotype (duration) analysis.** Among patients meeting any KDIGO creatinine criterion (AKI onset = time of first qualifying criterion; median onset 28.1 h [IQR 10.1–52.1] in MIMIC-IV and 33.7 h [15.6–59.6] in eICU-CRD), a **landmark analysis at onset + 72 hours** classified AKI as **transient** (SCr below 1.5× baseline and below baseline + 0.3 mg/dL at any measurement up to the landmark) or **persistent** (otherwise, including RRT before the landmark). Classification used only information available up to the landmark; patients who died (MIMIC-IV n = 1,540; eICU-CRD n = 2,421) or were discharged alive (n = 2,240; n = 3,131) before the landmark were excluded, as were those still hospitalised but without classifying measurements (n = 296; n = 655); these exclusions are reported and their selection implications discussed. From the landmark, a Cox model for 30-day mortality (from ICU admission) compared persistent with transient AKI, adjusted for the Model C covariates. This design classifies exposure strictly before follow-up begins, eliminating outcome-dependent exposure assignment.

**External validation (locked model).** The validation target was the fully observed binary outcome *in-hospital death within 30 days* (death = event; live discharge before day 30 or still hospitalised alive at day 30 = non-event), which requires no censoring assumptions and aligns with the cumulative-incidence estimand. A logistic regression on this outcome with 13 predictors (stage at 24 h, age, sex, non-renal SOFA, first-24-h ventilation, five comorbidities, prior CKD) was fitted in MIMIC-IV and **locked** (all coefficients and the intercept frozen). The locked model was applied to the entire eICU-CRD cohort without any refitting: discrimination was quantified by the AUC and the Brier score, calibration by the logistic calibration slope and the observed-to-expected (O:E) event ratio across linear-predictor deciles. Internal optimism was quantified by 5-fold cross-validation in MIMIC-IV, and a sensitivity analysis excluding patients still hospitalised at day 30 was pre-specified. **Cross-cohort parameter consistency** was assessed by refitting the identical time-varying Model C in eICU-CRD and computing the ratio of eICU to MIMIC hazard ratios (ratio ≈ 1 indicates parameter transportability).

**Sensitivity analyses.** Subgroup-specific Stage 3 hazard ratios (age ≥ 65, sex, prior CKD, diabetes, first-24-h ventilation, non-renal SOFA tertiles) with formal interaction tests were computed in both cohorts under the time-varying framework (Figure S2, Table S2).

All tests were two-sided. p-values below 0.001 are reported as p < 0.001. Analyses used Python 3.13.12 with pandas 2.3.3, lifelines 0.30.3, statsmodels 0.14.6, scipy 1.18.0, scikit-learn, and matplotlib 3.11.0. All code is publicly available (see Data Availability).

---

## Results

### Cohort Assembly and Follow-up

In MIMIC-IV, 94,458 ICU stays among 65,366 patients were screened. After excluding 3,991 stays with ICU length of stay <12 h, retaining the first eligible stay per patient (26,980 repeat stays excluded), excluding 2,083 patients with ESRD/chronic dialysis, and requiring a baseline creatinine and at least one SCr within the staging window (898 excluded: 425 without baseline, 473 without in-window SCr), **60,506 patients** formed the derivation cohort. In eICU-CRD, 200,859 unit stays among 139,367 patients were screened; after age exclusions (625), length-of-stay exclusion (27,842), patient-level deduplication (42,031), dual-source ESRD exclusion (6,523; ICD-identified n = 2,801, pastHistory-identified n = 6,219, with overlap), and requiring baseline and in-window creatinine (10,372 excluded), **113,466 patients** formed the validation cohort (Figure 1).

AKI prevalence differed substantially between databases: 28.9% in MIMIC-IV (maximum 7-day Stage 1/2/3 = 11,879/2,249/3,364) versus 16.1% in eICU-CRD (11,983/1,874/4,389). This asymmetry was largely attributable to baseline creatinine availability: 80.2% of MIMIC-IV patients had a pre-admission baseline creatinine, whereas only 39.4% of eICU-CRD patients had a pre-ICU baseline — 60.6% relied on the first ICU-period creatinine, which inflates the baseline and suppresses ratio-defined AKI (Table S1).

During 30-day follow-up, 5,697 MIMIC-IV patients (9.4%) and 9,373 eICU-CRD patients (8.3%) died in hospital; 52,546 (86.8%) and 102,414 (90.3%) were discharged alive; and 2,263 (3.7%) and 1,679 (1.5%) remained hospitalised at day 30. The high proportion discharged alive underscores the relevance of the competing-risk framework.

### Baseline Characteristics by Maximum KDIGO Stage (Table 1)

In both cohorts, higher KDIGO stage was associated with older age (except Stage 3), more mechanical ventilation, and higher non-renal SOFA (MIMIC-IV: 4.23/6.82/7.71/9.05 across Stages 0–3; eICU-CRD: 2.87/4.66/5.62/5.68). Liver disease in MIMIC-IV rose steeply with stage (6.4% → 39.6%), consistent with hepatorenal contributions to creatinine rise. Prior CKD was most frequent at Stage 3 (MIMIC-IV 8.8%; eICU-CRD 7.2%). Crude 30-day mortality increased monotonically with maximum stage in MIMIC-IV (5.4%/12.8%/26.5%/37.6%) but not in eICU-CRD (5.3%/19.2%/33.0%/31.1%), where Stage 3 did not exceed Stage 2 (see below).

### Time-Varying Cox Models (Table 2)

In MIMIC-IV, the unadjusted Stage 3 cause-specific hazard ratio (Model A: HR 5.71, 95% CI 5.31–6.13) attenuated by 41% after full adjustment (Model C: HR 3.35, 3.08–3.63). Non-renal SOFA was the dominant confounder (Model C HR per point 1.158, 1.148–1.168). The monotonic exposure-response was preserved after full adjustment (Stage 1: 1.81, 1.68–1.94; Stage 2: 2.76, 2.51–3.04; Stage 3: 3.35, 3.08–3.63).

In eICU-CRD, the fully adjusted model yielded Stage 3 HR 3.33 (3.12–3.55) — virtually identical to MIMIC-IV (cross-cohort ratio 1.00). Stage 1 and 2 hazards were higher in eICU-CRD (2.48, 2.36–2.61 and 3.54, 3.24–3.85) than in MIMIC-IV, consistent with the baseline-creatinine asymmetry concentrating more severely injured patients into eICU Stage 1–2 (see Discussion).

Other Model C covariates behaved as expected and were directionally consistent across cohorts (age HR 1.034 and 1.029 per year; liver disease 1.36 and 1.37); prior CKD was not independently associated with mortality in either cohort (HR 1.02, p = 0.71; 0.99, p = 0.89). First-24-h mechanical ventilation showed an apparently protective adjusted association in MIMIC-IV (HR 0.58, 0.54–0.62) but a harmful one in eICU-CRD (HR 1.47, 1.41–1.54); this heterogeneity is examined in the Discussion.

### Hazard Time-Dependence (Table 3, Figure 3)

Stage-by-period interactions revealed a strong, replicated time structure. In MIMIC-IV, the Stage 3 hazard ratio declined from **5.41 (4.91–5.96)** during days 0–7, to **2.10 (1.81–2.43)** during days 7–14, to **1.11 (0.92–1.34)** during days 14–30 (joint Wald χ² = 332.7, df = 6, p = 8 × 10⁻⁶⁹). The same pattern held in eICU-CRD: Stage 3 HR 5.07 (4.70–5.47) at 0–7 d, 1.87 (1.65–2.13) at 7–14 d, and 1.09 (0.90–1.30) at 14–30 d (joint Wald χ² = 545.2, p = 2 × 10⁻¹¹⁴). The pattern was graded across stages: by days 14–30, no stage retained a statistically detectable excess hazard in either cohort. Thus the excess mortality risk of AKI stages is concentrated in the first week after criteria are met; among 2-week survivors, stage-specific excess hazards are near null.

### Cumulative Incidence of Death under Competing Risks (Figure 2)

Stratifying by the stage attained at 24 hours, the 30-day Aalen-Johansen cumulative incidence of death rose monotonically in MIMIC-IV: 7.1% (95% CI 6.9–7.3) for Stage 0, 20.7% (19.6–21.7) for Stage 1, 29.6% (27.4–31.8) for Stage 2, and 33.6% (31.6–36.2) for Stage 3. In eICU-CRD the gradient was monotonic only up to Stage 2: 7.0% (6.8–7.1), 29.1% (27.6–30.4), 37.5% (33.5–40.7), and 26.9% (25.0–28.9) — Stage 3 did not exceed Stage 2. This non-monotonicity was attributable to the composition of Stage 3 at 24 hours in eICU-CRD: 62.2% of these patients had reached Stage 3 via RRT initiation (versus 39.0% in MIMIC-IV), and RRT-initiated eICU patients had *lower* 30-day mortality than creatinine-only Stage 3 patients (25.0% vs 30.0%), whereas the reverse held in MIMIC-IV (42.8% vs 27.8%). The overall 30-day CIF of live discharge was 86.8% (MIMIC-IV) and 90.3% (eICU-CRD).

### AKI Duration: Transient versus Persistent Phenotype (Figure S3, Table S3)

Among 17,492 MIMIC-IV AKI patients, 13,416 (76.7%) were evaluable at the onset+72-hour landmark (8,611 transient, 64.2%; 4,805 persistent); 1,540 died and 2,240 were discharged before the landmark, and 296 lacked classifying measurements. In eICU-CRD, 12,039 of 18,246 AKI patients were evaluable (7,077 transient, 58.8%; 4,962 persistent; 2,421 pre-landmark deaths, 3,131 pre-landmark discharges, 655 unclassifiable). Post-landmark 30-day mortality was 7.6% (transient) versus 23.8% (persistent) in MIMIC-IV, against 5.4% in patients without AKI; in eICU-CRD, 9.9% versus 21.6%, against 5.3%. In fully adjusted Cox models from the landmark, **persistent AKI carried roughly twice the hazard of transient AKI (HR 1.95, 95% CI 1.77–2.16 in MIMIC-IV; 1.79, 1.62–1.97 in eICU-CRD)** — a consistent, cross-cohort gradient that is smaller than crude contrasts (2.22 and 1.76 before adjustment) because stage, severity, and comorbidity partly explain the duration-mortality association. Transient AKI nevertheless carried a real excess mortality over no-AKI (7.6% vs 5.4%; 9.9% vs 5.3%).

### External Validation of the Locked MIMIC-IV Model (Table 4, Figure 4)

The 13-predictor logistic model for in-hospital death within 30 days was fitted in MIMIC-IV (n = 60,506; 5,697 events) and locked.

**Discrimination.** The apparent AUC in MIMIC-IV was 0.800 (5-fold cross-validated AUC 0.799, SD 0.003); the locked model achieved an AUC of 0.750 in eICU-CRD (n = 113,466; 9,373 events). Excluding patients still hospitalised at day 30 changed the AUCs negligibly (MIMIC-IV 0.806, n = 58,243; eICU-CRD 0.751, n = 111,787).

**Calibration.** The calibration slope in eICU-CRD was 0.946 (95% CI 0.924–0.969) — close to unity, indicating that the relative weighting of predictors transported well. However, the observed-to-expected event ratio was 1.55: the locked model systematically **underestimated** absolute risk in eICU-CRD (e.g., top predicted-risk decile: predicted 20.3%, observed 28.1%), reflecting the lower acuity mix (mean non-renal SOFA 3.2 vs 5.1) but similar crude mortality of the eICU population. The Brier score was 0.0723 (MIMIC-IV) and 0.0703 (eICU-CRD).

**Cross-cohort parameter consistency.** Hazard-ratio ratios (eICU/MIMIC) from the identical time-varying Model C were 1.00 for Stage 3, 1.00 for non-renal SOFA, 1.00 for age, and 0.95–1.10 for all comorbidities — indicating excellent transportability of the core parameter set. Stage 1 and 2 hazards were 37–28% higher in eICU-CRD (ratios 1.37 and 1.28), and first-24-h mechanical ventilation reversed direction (0.58 vs 1.47; ratio 2.55).

### Subgroup Analyses (Figure S2, Table S2)

The adjusted Stage 3 hazard ratio remained between 1.95 and 3.75 in every pre-specified subgroup of both cohorts, with no direction reversal. Two interactions were statistically significant in both databases: **sex** (higher Stage 3 hazard in males: MIMIC-IV 3.69, 3.32–4.11 vs 2.94, 2.60–3.34, p(interaction) = 0.003; eICU-CRD 3.57, 3.28–3.88 vs 3.03, 2.74–3.35, p = 0.056) and **first-24-h ventilation**, though the latter ran in opposite directions (MIMIC-IV: ventilated 3.65 vs not ventilated 3.33, p = 0.0008; eICU-CRD: ventilated 3.19 vs not ventilated 3.45, p < 0.0001) — a further manifestation of the ventilation heterogeneity discussed above. A notable cohort-specific finding was **prior CKD in eICU-CRD**, where the Stage 3 hazard was substantially attenuated among CKD patients (1.95, 1.44–2.65 vs 3.40, 3.19–3.64 without CKD, p = 0.003) but not in MIMIC-IV (2.90 vs 3.38, p = 0.16); plausibly, ratio-based staging against an already-elevated CKD baseline misclassifies chronic worsening as severe AKI less discriminatively. Age, diabetes, and non-renal SOFA tertiles showed no consistent interactions across cohorts.

---

## Discussion

### Principal Findings

In two large, independent, patient-level ICU cohorts totalling 173,972 patients, strictly windowed KDIGO creatinine staging — modelled as a time-varying exposure to eliminate immortal-time bias — showed a robust, monotonic association with 30-day in-hospital mortality that survived conservative severity adjustment using a deliberately non-renal SOFA score, and that was virtually identical across databases (Stage 3 HR 3.35 vs 3.33; cross-cohort ratio 1.00). The association was strongly time-dependent, concentrated in the first week after criteria are met and near null among 2-week survivors. AKI duration added prognostic information: persisting AKI carried roughly twice the hazard of transient AKI in both cohorts. A locked prediction model transported with preserved discrimination (AUC 0.800 → 0.750) and relative calibration (slope 0.946) but systematically underestimated absolute risk in the validation cohort (O:E 1.55), a transparent illustration of the distinction between relative and absolute transportability.

### Strength of the Association and Severity Confounding

Our fully adjusted Stage 3 hazard ratio of 3.35 (MIMIC-IV) and 3.33 (eICU-CRD) is more conservative than the crude 5.71, quantifying the confounding contribution of non-renal organ failure at roughly two-fifths of the crude effect — directly addressing the correlation-versus-causation question raised by Girling et al. [24]. Importantly, we adjusted for a **non-renal** SOFA score. An additional eICU sensitivity model replacing non-renal SOFA with the APACHE acute physiology score — which embeds renal function — attenuated the Stage 3 hazard ratio from 3.33 to 2.29 (2.13–2.45), a useful empirical demonstration that severity scores embedding renal function absorb exposure information and should not be used unmodified when quantifying the independent prognostic value of AKI staging.

Notably, the time-varying Stage 3 hazard (3.35) exceeds the estimate obtained under the previous fixed-exposure landmark implementation of the same cohort (2.45). This is expected: assigning the eventual maximum stage from time zero forces patients who attain Stage 3 late (and must survive to do so) to contribute person-time at Stage 3 before their criteria are met, diluting the stage contrast; the counting-process implementation removes this dilution as well as the opposite immortal-time distortion. Both biases operate in fixed-exposure designs, in directions that partially mask each other — reinforcing the case for time-varying exposure modelling as the default in AKI database studies [17].

### Time-Dependence of AKI Hazards

The marked attenuation of the Stage 3 hazard from 5.41 (days 0–7) to 1.11 (days 14–30), replicated almost exactly in eICU-CRD, has two implications. First, single averaged hazard ratios reported by most registry studies mix a large early excess with a near-null late phase and therefore overstate the residual risk carried by survivors of severe AKI; a clinician interpreting a Stage 3 diagnosis in a patient alive at day 14 should weigh a much smaller excess hazard. Second, this structure justifies explicit hazard-period modelling of the kind advocated by Singh et al. after proportional hazards violations [28] — in our design, the stage-by-period joint Wald tests (p < 10⁻⁶⁸ in both cohorts) formalise what global PH tests on averaged effects would understate.

### Absolute Risks and the Competing-Risk Framework

Because 87–90% of patients were discharged alive within 30 days, Kaplan-Meier-type estimates that treat live discharge as random censoring conflate "discharged alive" with "event-free at risk," and the CIF estimates we report are the appropriate absolute-risk quantities. Under competing risks, the 30-day probability of death rose from 7.1% (Stage 0 at 24 h) to 33.6% (Stage 3 at 24 h) in MIMIC-IV — a flatter gradient than Kaplan-Meier contrasts would suggest, and one that reflects the real clinical competition between recovery-with-discharge and death.

The eICU-CRD non-monotonicity (Stage 3 at 24 h, 26.9%, below Stage 2, 37.5%) warrants explicit attention. It is an artefact of exposure composition rather than evidence against staging: 62% of eICU patients classified as Stage 3 at 24 hours reached that stage through RRT initiation, versus 39% in MIMIC-IV, and RRT-initiated eICU patients had lower 30-day mortality (25.0%) than creatinine-only Stage 3 patients (30.0%). Early RRT in eICU-CRD appears to select or rescue a group whose mortality does not follow the creatinine gradient; in MIMIC-IV, where RRT-received Stage 3 patients had higher mortality (42.8%), the monotonicity is preserved. RRT-based stage assignment therefore mixes an illness-severity signal with a treatment-decision signal whose direction varies between practice environments — a previously underappreciated threat to the interpretability of pooled KDIGO gradients.

### AKI Duration Adds Prognostic Information

At the onset+72-hour landmark — where phenotype classification strictly precedes follow-up, eliminating outcome-dependent exposure assignment — persistent AKI carried roughly twice the adjusted hazard of transient AKI (1.95 and 1.79), with cross-cohort consistency. Two caveats temper this finding. First, the landmark excluded 23% of MIMIC-IV and 34% of eICU-CRD AKI patients (pre-landmark deaths and discharges); the comparison is conditional on surviving and remaining hospitalised to 72 hours after onset, so absolute risks must not be extrapolated to all AKI patients. Second, adjustment attenuated the crude contrast (2.22 → 1.95 in MIMIC-IV), indicating that part of the duration-mortality association reflects severity and comorbidity rather than duration per se. Nevertheless, the direction and magnitude replicated across databases, and transient AKI itself carried a modest excess over no-AKI (7.6–9.9% vs 5.3–5.4% crude). Duration is bedside-available at day 3 of AKI and meaningfully refines the peak-stage snapshot — consistent with trajectory-based subphenotyping reports [25] — but it complements rather than replaces staging.

### External Validation: What Transports and What Does Not

To our knowledge, this is among the first KDIGO staging studies to perform genuine locked-model external validation across independent critical care databases in the TRIPOD spirit, rather than refitting in the validation cohort. The validation target — in-hospital death within 30 days, with live discharge and continued hospitalisation counted as non-events — is fully observed, requiring no censoring assumptions and aligning with the cumulative-incidence estimand of the primary analysis. Discrimination transported acceptably (AUC 0.800 → 0.750; cross-validated 0.799), and the calibration slope of 0.946 shows that the *relative* weighting of predictors was stable. However, the O:E ratio of 1.55 means the locked MIMIC-IV intercept understates eICU-CRD absolute risk by roughly a third: eICU-CRD patients have lower measured acuity (non-renal SOFA 3.2 vs 5.1) yet similar crude mortality, so part of their risk is not captured by the transported covariates. The practical implication is direct: a transported model's discrimination and slope may validate while its absolute predictions fail, and recalibration (intercept adjustment) is mandatory before any clinical use. Most published "external validations" that refit the intercept or entire model would mask exactly this failure mode.

At the parameter level, Stage 3 (ratio 1.00), age, non-renal SOFA, and all comorbidities transported essentially perfectly. Stage 1 and 2 hazards were 28–37% higher in eICU-CRD, plausibly reflecting the baseline creatinine availability asymmetry: a higher baseline suppresses ratio-based staging, so eICU Stage 1–2 patients represent more truly injured patients relative to MIMIC-IV (supported by their higher stage-specific mortality). The **mechanical ventilation reversal** (HR 0.58 apparently protective in MIMIC-IV vs 1.47 harmful in eICU-CRD, both adjusted for non-renal SOFA) persists under the time-varying framework and warrants explicit caution. The two databases operationalise day-1 ventilation differently (delivered invasive ventilation procedures in MIMIC-IV versus a ventilated-status flag in eICU-CRD), differ in hospital mix and ventilation practice, and early ventilation is strongly conditioned on surviving initial support; residual confounding by indication in either direction cannot be excluded. Crucially, the KDIGO estimates are not affected by this instability: excluding ventilation from Model C increased the MIMIC-IV Stage 3 hazard ratio only from 3.35 to 3.59 (3.32–3.90), confirming that the KDIGO coefficients are not driven by the ventilation term. We retain ventilation in the locked model for covariate-set fidelity, flag its coefficient as non-transportable, and recommend cohort-specific recalibration of any transported model that includes it [25,28].

### Comparison with Prior Literature

The AKI-EPI study reported Stage 3 mortality of 50.0% [2]; Hoste et al. reported RIFLE Failure mortality of 26.3% [33]; FINNAKI reported KDIGO Stage 3 90-day mortality of 39.0% [34]. Our 30-day Stage 3 CIF of 33.6% (MIMIC-IV) sits within this range, as expected for creatinine-only staging (urine output criteria add substantial, mostly earlier-onset AKI [16,28]). Our lower AKI prevalence versus AKI-EPI reflects both the creatinine-only approach and our exclusion of ESRD, and the MIMIC-IV/eICU-CRD difference (28.9% vs 16.1%) is a database property — chiefly baseline creatinine availability — not a biological finding, and itself illustrates how implementation choices drive the apparent epidemiology of AKI [30,32].

### Limitations

First, staging was creatinine-based only; urine output criteria were not applied, so our stages represent the creatinine-defined subset of AKI, generally of later onset [16,28]. Second, baseline creatinine was unmeasurable for a majority of eICU-CRD patients; our hierarchical definition mitigates but cannot eliminate misclassification, and prior CKD identification relies on history/coding with limited sensitivity [30,32]. Third, the time-varying stage was modelled as monotonically non-decreasing; real creatinine trajectories fluctuate, and our design deliberately separates "peak severity attained by time t" (primary exposure) from "recovery" (phenotype analysis), at the cost of not modelling stage resolution as a time-varying process. Fourth, cause-specific Cox models censored live discharge under an independent-censoring assumption; patients discharged alive may have systematically different post-discharge hazards, and our CIF analysis addresses the absolute-risk consequences of discharge without modelling post-discharge death (no vital-status follow-up exists in either database). Fifth, the phenotype landmark excluded pre-landmark deaths and discharges, conditioning the duration-mortality comparison on 72-hour survival and hospitalisation. Sixth, the validation outcome was in-hospital death within 30 days — a composite of timing and discharge-alive status chosen for its complete observability; studies with post-discharge vital status may find different transportability. Seventh, the two databases differ in era (2008–2022 vs 2014–2015), hospital mix, and ventilation/severity ascertainment; we quantified rather than eliminated these differences. Finally, although we adjusted for measured confounders and used a non-renal severity score, unmeasured confounding (e.g., indication for RRT, nephrotoxin exposure) persists; the E-value for the fully adjusted Stage 3 hazard ratio (HR 3.35) is approximately 6.2, meaning an unmeasured confounder would need to be associated with both Stage 3 AKI and death by a risk ratio of at least this magnitude, above and beyond the measured covariates, to explain away the finding [22].

### Conclusions

KDIGO creatinine staging, implemented with its mandated temporal windows in patient-level cohorts, modelled as a time-varying exposure, and adjusted conservatively for non-renal illness severity, retains a graded association with 30-day ICU mortality that is highly reproducible across two independent databases — with the caveats that the excess hazard of severe AKI is concentrated in the first week, that RRT-based Stage 3 assignment can distort absolute-risk gradients in a practice-dependent manner, and that AKI persisting beyond 72 hours carries roughly twice the hazard of transient AKI. A locked MIMIC-IV model transported to a multi-centre cohort with preserved discrimination and relative calibration but required intercept recalibration for absolute risk, and individual parameter transportability was heterogeneous, with mechanical ventilation showing cohort-dependent direction. Future dual-database studies should model AKI stage as time-varying, treat live discharge as a competing event, incorporate AKI duration alongside peak stage, and validate locked rather than refitted models.

---

## Tables

### Table 1. Baseline characteristics by maximum 7-day KDIGO stage, derivation (MIMIC-IV) and validation (eICU-CRD) cohorts

| Characteristic | MIMIC-IV Stage 0 (n=43,014) | Stage 1 (n=11,879) | Stage 2 (n=2,249) | Stage 3 (n=3,364) | eICU Stage 0 (n=95,220) | Stage 1 (n=11,983) | Stage 2 (n=1,874) | Stage 3 (n=4,389) |
|---|---|---|---|---|---|---|---|---|
| Age, years, mean | 63.5 | 68.8 | 66.6 | 63.5 | 63.0 | 67.1 | 65.7 | 62.7 |
| Male sex, % | 54.4 | 61.6 | 53.8 | 61.4 | 53.7 | 56.5 | 49.4 | 61.1 |
| Non-renal SOFA (first 24 h), mean | 4.23 | 6.82 | 7.71 | 9.05 | 2.87 | 4.66 | 5.62 | 5.68 |
| Mechanical ventilation (first 24 h), % | 27.6 | 45.7 | 41.9 | 44.0 | 20.2 | 37.8 | 43.1 | 41.1 |
| Diabetes, % | 23.7 | 35.6 | 35.0 | 35.8 | 25.1 | 33.2 | 28.3 | 35.6 |
| Hypertension, % | 59.0 | 73.3 | 67.4 | 65.1 | 48.1 | 58.3 | 54.3 | 54.2 |
| Heart failure, % | 17.3 | 35.3 | 32.3 | 32.9 | 11.2 | 19.9 | 16.3 | 16.7 |
| COPD, % | 15.8 | 18.2 | 18.1 | 16.6 | 13.3 | 16.0 | 16.1 | 12.1 |
| Liver disease, % | 6.4 | 11.5 | 24.0 | 39.6 | 2.1 | 3.3 | 3.5 | 5.4 |
| Prior CKD, % | 2.0 | 5.5 | 3.9 | 8.8 | 1.5 | 3.6 | 1.6 | 7.2 |
| 30-day in-hospital mortality, % | 5.4 | 12.8 | 26.5 | 37.6 | 5.3 | 19.2 | 33.0 | 31.1 |

CKD, chronic kidney disease; COPD, chronic obstructive pulmonary disease; SOFA, Sequential Organ Failure Assessment. Stage = maximum KDIGO creatinine stage attained within 7 days of ICU admission.

### Table 2. Time-varying cause-specific Cox models for 30-day in-hospital mortality by cohort

| Variable | MIMIC-IV Model A HR (95% CI) | Model B HR (95% CI) | Model C HR (95% CI) | eICU-CRD Model A HR (95% CI) | Model B HR (95% CI) | Model C HR (95% CI) |
|---|---|---|---|---|---|---|
| KDIGO Stage 1 (time-varying) | 2.51 (2.34–2.68) | 2.34 (2.19–2.51) | 1.81 (1.68–1.94) | 3.50 (3.32–3.68) | 3.34 (3.17–3.51) | 2.48 (2.36–2.61) |
| KDIGO Stage 2 (time-varying) | 4.28 (3.90–4.69) | 4.18 (3.81–4.59) | 2.76 (2.51–3.04) | 5.36 (4.92–5.84) | 5.29 (4.86–5.76) | 3.54 (3.24–3.85) |
| KDIGO Stage 3 (time-varying) | 5.71 (5.31–6.13) | 6.08 (5.66–6.54) | 3.35 (3.08–3.63) | 4.48 (4.21–4.77) | 4.82 (4.52–5.13) | 3.33 (3.12–3.55) |
| Age (per year) | — | 1.025 (1.023–1.027) | 1.034 (1.032–1.036) | — | 1.024 (1.022–1.025) | 1.029 (1.027–1.030) |
| Male sex | — | 0.889 (0.844–0.937) | 0.856 (0.812–0.902) | — | 0.962 (0.923–1.002) | 0.915 (0.878–0.953) |
| Non-renal SOFA (per point) | — | — | 1.158 (1.148–1.168) | — | — | 1.158 (1.151–1.165) |
| Ventilation (first 24 h) | — | — | 0.576 (0.537–0.618) | — | — | 1.473 (1.408–1.542) |
| Prior CKD | — | — | 1.024 (0.904–1.160) | — | — | 0.991 (0.881–1.115) |
| Diabetes | — | — | 0.880 (0.829–0.933) | — | — | 0.876 (0.836–0.918) |
| Hypertension | — | — | 0.798 (0.752–0.846) | — | — | 0.884 (0.847–0.922) |
| Heart failure | — | — | 0.959 (0.904–1.018) | — | — | 1.088 (1.030–1.150) |
| COPD | — | — | 1.123 (1.052–1.200) | — | — | 1.231 (1.168–1.298) |
| Liver disease | — | — | 1.357 (1.266–1.455) | — | — | 1.372 (1.246–1.511) |
| Events / person-time rows | 5,697 / 116,206 | 5,697 / 116,206 | 5,697 / 116,206 | 9,373 / 180,829 | 9,373 / 180,829 | 9,373 / 180,829 |

Model A: time-varying KDIGO stage only. Model B: + age, sex. Model C: + non-renal SOFA, ventilation, comorbidities, prior CKD. Live discharge censored (cause-specific hazards). All KDIGO coefficients p < 0.001. An eICU-CRD sensitivity model substituting the APACHE acute physiology score for non-renal SOFA (n = 96,150) yielded Stage 3 HR 2.29 (2.13–2.45); a MIMIC-IV sensitivity model excluding ventilation yielded Stage 3 HR 3.59 (3.32–3.90).

### Table 3. Period-specific adjusted hazard ratios from stage-by-period interaction models (Model C extended)

| Interval | MIMIC-IV Stage 1 HR (95% CI) | Stage 2 HR (95% CI) | Stage 3 HR (95% CI) | eICU-CRD Stage 1 HR (95% CI) | Stage 2 HR (95% CI) | Stage 3 HR (95% CI) |
|---|---|---|---|---|---|---|
| Days 0–7 | 2.17 (1.98–2.37) | 3.99 (3.55–4.48) | 5.41 (4.91–5.96) | 3.18 (2.99–3.38) | 5.38 (4.87–5.95) | 5.07 (4.70–5.47) |
| Days 7–14 | 1.31 (1.15–1.50) | 1.73 (1.43–2.10) | 2.10 (1.81–2.43) | 1.56 (1.40–1.74) | 1.61 (1.32–1.96) | 1.87 (1.65–2.13) |
| Days 14–30 | 1.09 (0.93–1.29) | 1.08 (0.84–1.39) | 1.11 (0.92–1.34) | 1.08 (0.91–1.27) | 1.24 (0.95–1.63) | 1.09 (0.90–1.30) |

Joint Wald test of stage-by-period interactions (df = 6): MIMIC-IV χ² = 332.7, p = 8 × 10⁻⁶⁹; eICU-CRD χ² = 545.2, p = 2 × 10⁻¹¹⁴.

### Table 4. External validation of the locked MIMIC-IV logistic model for in-hospital death within 30 days in eICU-CRD (no refitting; TRIPOD type 3)

| Metric | MIMIC-IV (derivation, apparent) | eICU-CRD (validation, locked) |
|---|---|---|
| N (events) | 60,506 (5,697) | 113,466 (9,373) |
| AUC | 0.800 | 0.750 |
| AUC, 5-fold cross-validation | 0.799 (SD 0.003) | — |
| Brier score | 0.0723 | 0.0703 |
| Calibration slope (95% CI) | 1.0 (by construction) | 0.946 (0.924–0.969) |
| Observed-to-expected ratio | 1.0 | 1.55 |
| Highest predicted-risk decile, predicted vs observed | 37.9% vs 38.6% | 20.3% vs 28.1% |
| Sensitivity excluding still-hospitalised at day 30 | 0.806 (n = 58,243) | 0.751 (n = 111,787) |

Predictors (13): stage at 24 h (indicators), age, sex, non-renal SOFA, first-24-h ventilation, prior CKD, diabetes, hypertension, heart failure, COPD, liver disease. Outcome: in-hospital death ≤30 days (live discharge or continued hospitalisation = non-event). Cross-cohort hazard-ratio ratios from the identical time-varying Model C: Stage 3 1.00; non-renal SOFA 1.00; age 1.00; Stage 1 1.37; Stage 2 1.28; ventilation 2.55 (direction reversal).

---

## Figure Legends

**Figure 1.** Cohort flow diagram for MIMIC-IV (derivation) and eICU-CRD (validation), showing patient-level deduplication, end-stage renal disease exclusion, and baseline/window creatinine requirements. Follow-up begins at ICU admission; no landmark exclusion applies to the primary analysis.

**Figure 2.** Aalen-Johansen cumulative incidence of 30-day in-hospital death by KDIGO stage attained at 24 hours, with live discharge as a competing event. (A) MIMIC-IV. (B) eICU-CRD. Shaded areas, bootstrap 95% confidence intervals.

**Figure 3.** Forest plot of period-specific adjusted hazard ratios for KDIGO stages (days 0–7, 7–14, 14–30) from stage-by-period interaction models in both cohorts, demonstrating early hazard concentration and subsequent attenuation.

**Figure 4.** External validation of the locked MIMIC-IV logistic model in eICU-CRD. (A) ROC curves, derivation (apparent) and validation (locked). (B) Calibration by predicted-risk decile (observed vs predicted 30-day in-hospital mortality), illustrating preserved relative calibration (slope 0.946) with systematic underestimation of absolute risk (O:E 1.55).

**Figure S1 (Supplementary).** Directed acyclic graph of the assumed causal structure, indicating non-renal SOFA, ventilation, comorbidities, age, and sex as measured confounders and the renal SOFA component as part of the exposure-outcome pathway.

**Figure S2 (Supplementary).** Forest plot of fully adjusted, time-varying Stage 3 vs Stage 0 hazard ratios by pre-specified subgroup (age, sex, prior CKD, diabetes, ventilation, non-renal SOFA tertiles) in MIMIC-IV and eICU-CRD, with interaction p-values.

**Figure S3 (Supplementary).** Adjusted 30-day mortality hazard ratios for persistent versus transient AKI (72-hour resolution criterion, onset+72-hour landmark) in both cohorts.

**Table S1 (Supplementary).** Baseline creatinine source distribution by cohort and KDIGO stage.

**Table S2 (Supplementary).** Subgroup-specific Stage 3 hazard ratios and interaction tests in both cohorts (time-varying framework).

**Table S3 (Supplementary).** AKI phenotype (transient vs persistent) attrition, 30-day mortality, and fully adjusted Cox models at the onset+72-hour landmark in both cohorts.

---

## Abbreviations

AKI, acute kidney injury; APACHE, Acute Physiology and Chronic Health Evaluation; AUC, area under the receiver operating characteristic curve; CIF, cumulative incidence function; CKD, chronic kidney disease; COPD, chronic obstructive pulmonary disease; DAG, directed acyclic graph; eICU-CRD, eICU Collaborative Research Database; ESRD, end-stage renal disease; GCS, Glasgow Coma Scale; HR, hazard ratio; ICU, intensive care unit; IQR, interquartile range; KDIGO, Kidney Disease: Improving Global Outcomes; MIMIC-IV, Medical Information Mart for Intensive Care IV; PH, proportional hazards; RRT, renal replacement therapy; SCr, serum creatinine; SOFA, Sequential Organ Failure Assessment; STROBE, Strengthening the Reporting of Observational Studies in Epidemiology; TRIPOD, Transparent Reporting of a multivariable prediction model for Individual Prognosis Or Diagnosis.

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

MIMIC-IV is publicly available at PhysioNet (https://physionet.org/content/mimiciv/). eICU-CRD is publicly available at PhysioNet (https://physionet.org/content/eicu-crd/). All analysis code (cohort construction, SOFA computation, time-varying exposure construction, competing-risk and phenotype analyses, and locked-model external validation) and key result files are publicly available at GitHub (https://github.com/wudengke2010/kdigo-aki-mimic-iv).

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

The authors thank the PhysioNet team and the Laboratory for Computational Physiology at the Massachusetts Institute of Technology for creating and maintaining the MIMIC-IV and eICU-CRD databases, which made this study possible. The authors also acknowledge the MIT-LCP `mimic-code` repository, whose validated first-day SOFA SQL specification informed the severity score implementation.

---

## Disclosure of interest

The authors report there are no competing interests to declare.

---

## ORCID

Dengke Wu https://orcid.org/0000-0003-4101-8461

---

## Declaration of Generative AI in Scientific Writing

During the preparation of this work, the authors used large language models (LLMs) for language editing and copyediting. After using these tools, the authors reviewed and edited the content as needed and take full responsibility for the content of the publication.

---

*Version 5 (v3 reanalysis). This manuscript replaces all analyses in prior versions following a second independent methodological review that identified time-structure and estimand errors in the previous landmark implementation (immortal time from assigning the 7-day maximum stage as a fixed exposure, 30-day truncation that retained post-truncation deaths as events, live discharge treated as random censoring, invalid time-stratified risk sets, outcome-dependent phenotype classification, and pseudo-locked external validation). All results derive from the v3 reanalysis pipeline (v3_build_exposure, v3_primary_analysis, v3_phenotype, v3_validation, v5_supp), with cohort flow and interval-tiling audits passing automatically. The eGFR-imputed and pre-admission-baseline re-staging sensitivities from prior versions were discontinued because they were implemented under the invalidated fixed-exposure framework.*
