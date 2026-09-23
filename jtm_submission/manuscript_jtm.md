# KDIGO Creatinine Staging and Transient Versus Persistent AKI Phenotypes Are Independently Associated with 30-Day Mortality in Critically Ill Patients: A Dual-Cohort Study with Locked-Model External Validation (MIMIC-IV v3.1 and eICU-CRD)

---

## Authors

**Jiqiang Liu**<sup>1,2</sup>, Dengke Wu<sup>1,2*</sup>

<sup>1</sup> Department of Emergency Medicine, The Second Xiangya Hospital of Central South University, Changsha 410011, Hunan, China

<sup>2</sup> Emergency Medicine and Difficult Diseases Institute, The Second Xiangya Hospital of Central South University, Changsha 410011, Hunan, China

**\*Corresponding author:** Dengke Wu, Department of Emergency Medicine, Second Xiangya Hospital, Central South University, Changsha, Hunan, China; Emergency Medicine and Difficult Diseases Institute, The Second Xiangya Hospital of Central South University, Changsha 410011, Hunan, China. Electronic address: wudk2010@csu.edu.cn.

---

## Abstract

**Background:** The KDIGO creatinine-based AKI staging system is widely used, but many large validation studies applied the criteria without their required temporal windows, staged AKI using entire-admission creatinine values, and adjusted for severity scores that embed renal function. The independent, time-resolved association between strictly implemented KDIGO stage and mortality, and the prognostic value of AKI duration relative to peak stage, remain incompletely characterised.

**Methods:** In this STROBE- and TRIPOD-compliant retrospective dual-cohort study, patient-level cohorts (one ICU stay per patient) were constructed from MIMIC-IV v3.1 (derivation) and eICU-CRD v2.0 (external validation). AKI was staged according to KDIGO creatinine criteria with explicit temporal windows (≥0.3 mg/dL within 48 h; ≥1.5×/≥2.0×/≥3.0× baseline within 7 days), a hierarchical baseline creatinine definition prioritising pre-admission values, and end-stage renal disease exclusion. Illness severity was adjusted using a first-24-h non-renal SOFA score. A 24-hour landmark analysis with real event times assessed 30-day mortality: Kaplan-Meier estimates, incremental Cox models, proportional hazards testing, time-stratified Cox models, and a transient (resolution within 72 h) versus persistent AKI phenotype analysis. The fully adjusted MIMIC-IV model was locked and applied, without any refitting, to eICU-CRD.

**Results:** Of 60,506 MIMIC-IV (AKI 28.9%) and 113,466 eICU-CRD patients (AKI 16.1%), 58,554 and 107,758 formed the landmark cohorts. Thirty-day mortality rose monotonically across stages in both cohorts. Adjusted Stage 3 hazard ratios were 2.45 (95% CI 2.26–2.66) in MIMIC-IV and 2.46 (2.31–2.63) in eICU-CRD (cross-cohort ratio 1.00), declining from 5.25 (days 0–7) to 1.30 (days 14–30) in MIMIC-IV. Locked-model external validation yielded a C-index of 0.716 (0.708–0.722) versus 0.753 (0.746–0.760) in derivation, calibration slope 0.96, and tertile 30-day mortality of 19.9%/32.1%/52.4%. Transient AKI carried a near-null adjusted hazard (HR 1.07, 0.99–1.15 in MIMIC-IV; 1.24, 1.16–1.32 in eICU-CRD), whereas persistent AKI carried 2.63- and 3.24-fold hazards; the phenotype model discriminated better than peak stage alone (C-index 0.768 vs 0.753; 0.789 vs 0.716).

**Conclusions:** Strictly implemented KDIGO creatinine staging retains a graded, highly transportable association with 30-day mortality, with the excess hazard concentrated in the first week and, decisively, in AKI persisting beyond 72 hours. Creatinine trajectory outperforms peak stage and should complement staging at the bedside.

**Keywords:** acute kidney injury; KDIGO; intensive care unit; mortality; SOFA; external validation

---

## Introduction

Acute kidney injury (AKI) is one of the most common and consequential complications encountered in the intensive care unit (ICU). Epidemiological studies have consistently demonstrated that AKI affects approximately 30–60% of critically ill patients, with the landmark AKI-EPI multinational cohort reporting an incidence of 57.7% and a mortality gradient from 16.1% (no AKI) to 50.0% (Stage 3) [1,2].

The evolution of consensus definitions—from RIFLE to AKIN to the current Kidney Disease: Improving Global Outcomes (KDIGO) criteria—has standardised AKI classification and revealed a graded association between AKI severity and mortality across diverse ICU populations [19]. The KDIGO 2012 guideline established a staging system integrating serum creatinine (SCr), urine output, and renal replacement therapy (RRT) criteria [3], validated across general ICU, cardiac surgery, and sepsis populations [4–7]. However, the methodological quality of large database studies implementing these criteria varies substantially, and several threats to validity recur.

**First**, the KDIGO creatinine criteria are explicitly time-bound: the absolute criterion (increase ≥0.3 mg/dL) applies within a 48-hour window, and ratio criteria (≥1.5× baseline) apply within a 7-day window. Many retrospective studies stage AKI using the peak-to-baseline creatinine ratio over the entire admission, ignoring these windows entirely. Yasrebi-de Kom et al. documented widespread incorrect application of KDIGO staging criteria, including the misinterpretation of the Stage 3 absolute threshold (an increase *to* ≥4.0 mg/dL, not *of* ≥4.0 mg/dL) [20]. Window-free staging inflates apparent AKI severity by attributing slow creatinine drifts over weeks to "acute" injury, and inflates mortality associations through reverse causation (dying patients accumulate more measurements). Whether the KDIGO-mortality association survives *strict* window-based implementation at scale is not established.

**Second**, illness severity confounding has not been adequately resolved. Girling et al. applied the Bradford Hill criteria to assess whether AKI is causally associated with critical illness mortality or merely a correlate of illness severity, concluding that formal confounding quantification is needed [24]. A specific methodological trap is that commonly used severity scores (SOFA, and derived APACHE components) embed a renal term that is a mathematical function of the exposure itself; adjusting for such scores overadjusts and biases KDIGO effect estimates toward the null [19]. Studies quantifying the KDIGO-mortality association while adjusting for a **non-renal** severity score in large cohorts are lacking.

**Third**, the temporal structure of the AKI-mortality association is rarely modelled. AKI stage is a time-dependent exposure, and its hazard is plausibly non-proportional—early deaths after AKI onset are frequent, whereas survivors of severe AKI who leave the ICU face competing hazards. Singh et al. encountered proportional hazards violations in a bicentric ICU-AKI cohort and resolved them with time-stratified estimates [28], but this approach remains uncommon in large database studies, which typically report a single average hazard ratio.

**Fourth**, claims of "external validation" in dual-database studies frequently involve refitting all model coefficients in the validation cohort, which validates the *covariate structure* but not the *model*. Transportability requires applying locked derivation coefficients to the validation cohort, in the spirit of the TRIPOD framework. The precedent for MIMIC-IV/eICU-CRD derivation-validation pairing was established by Lin et al. for AKI trajectory subphenotyping [25], but locked-model transport of a KDIGO-based mortality model has not, to our knowledge, been reported.

We therefore conducted a dual-cohort study in MIMIC-IV v3.1 (derivation) and eICU-CRD v2.0 (external validation) with the following objectives: (1) implement KDIGO creatinine staging with its mandated temporal windows and a hierarchical pre-admission-referenced baseline creatinine; (2) quantify the independent association between KDIGO stage and 30-day mortality using a landmark design, real event times, and a non-renal first-24-h SOFA score for severity adjustment; (3) characterise the time-dependence of AKI hazards through proportional hazards testing and time-stratified Cox models; and (4) externally validate the locked, fully adjusted derivation model in eICU-CRD without any refitting, quantifying discrimination, calibration, and cross-cohort parameter consistency.

---

## Methods

### Study Design and Data Sources

This was a retrospective dual-cohort study. The **derivation cohort** used the Medical Information Mart for Intensive Care IV (MIMIC-IV) database, version 3.1, containing de-identified electronic health records of patients admitted to ICUs at Beth Israel Deaconess Medical Center (Boston, MA, USA) between 2008 and 2022 [12]. The **external validation cohort** used the eICU Collaborative Research Database (eICU-CRD), version 2.0, a multi-centre database of ICU stays from 208 hospitals across the United States, collected 2014–2015 [23]. The study is reported in accordance with the STROBE statement for cohort studies; the external validation component follows the TRIPOD framework (development plus independent external validation of a locked model).

### Study Population

In both databases, we constructed **patient-level** cohorts: for each patient, only the first ICU stay meeting all eligibility criteria was retained, eliminating within-person correlation from repeated stays. Eligibility required age ≥18 years and ICU length of stay ≥12 hours. Patients with end-stage renal disease (ESRD) or chronic dialysis were excluded, identified in MIMIC-IV from ICD-9/ICD-10 diagnosis codes (585.6, N18.6, V45.1, Z99.2; codes stored without decimal points in MIMIC-IV; restricted to diagnoses dated on or before the index admission) and in eICU-CRD from both ICD codes and structured `pastHistory` entries (chronic haemodialysis, peritoneal dialysis, or non-dialysed renal failure; dual-source identification with overlap). Patients without a usable baseline creatinine or without any serum creatinine within the staging window were excluded. Full cohort flow is shown in Figure 1.

### Exposure: KDIGO Creatinine-Based AKI Staging with Temporal Windows

The primary exposure was maximum AKI severity within 7 days of ICU admission, staged according to the KDIGO 2012 **creatinine** criteria with their mandated temporal windows (urine output criteria were not applied because hourly urine records are incompletely captured in both databases):

- **Stage 1:** increase in SCr ≥0.3 mg/dL (26.5 µmol/L) within any **48-hour rolling window**, or increase to ≥1.5× baseline within **7 days** of ICU admission;
- **Stage 2:** increase to ≥2.0× and <3.0× baseline within 7 days;
- **Stage 3:** increase to ≥3.0× baseline within 7 days, **or** SCr increase to ≥4.0 mg/dL (353.6 µmol/L) with an acute change, or initiation of RRT within 7 days.

The Stage 3 absolute threshold was implemented as an increase *to* ≥4.0 mg/dL requiring concurrent acute change (ratio ≥1.5× baseline or rise ≥0.3 mg/dL within 48 h), per Yasrebi-de Kom et al. [20]. In MIMIC-IV, SCr was extracted from `labevents` (itemid 50912) with physiological range filtering (0.1–100 mg/dL); in eICU-CRD, from the `lab` table (labname = 'creatinine'). RRT was identified from dialysis procedure events (`procedureevents`; CRRT, CVVHD, CVVHDF, haemodialysis, SCUF) in MIMIC-IV and from the `treatment` table (dialysis, CRRT, or ultrafiltration entries, excluding catheter insertion and radiology procedures) in eICU-CRD, in both cases restricted to the first 7 ICU days of the index stay.

**Baseline SCr hierarchy.** Baseline creatinine was defined hierarchically, prioritising pre-admission kidney function. In MIMIC-IV: (1) median outpatient SCr 7–365 days before the index admission; (2) if unavailable, median SCr from any care setting 7–365 days before admission; (3) if unavailable, median SCr 0–7 days before admission; (4) if unavailable, first SCr within 24 h of hospital admission. In eICU-CRD: (1) median SCr between hospital admission (or 7 days before ICU admission, whichever is later) and ICU admission; (2) if unavailable, first SCr within the first 24 ICU hours. The source of each patient's baseline was recorded and reported.

### Outcome: 30-Day Mortality with Real Event Times

The primary outcome was death within 30 days of ICU admission. **Event times were derived from recorded death and discharge timestamps**, not from length-of-stay surrogates: in MIMIC-IV, time-to-event was `deathtime − intime` for decedents and `dischtime − intime` (censored) for survivors; 351 decedents whose recorded death time exceeded their discharge time by more than 2 hours (a known data anomaly) had event times capped at discharge time + 24 hours, and five decedents lacking any death timestamp were excluded by the landmark. In eICU-CRD, `hospitaldischargeoffset` (minutes, converted to hours) defined follow-up with `hospitaldischargestatus` defining the event. Follow-up was truncated at 30 days.

**Landmark design.** Because KDIGO stage is ascertained over the first 7 days, patients dying before staging can be assessed are misclassified as lower stages (immortal time bias). A **24-hour landmark analysis** was therefore the primary design: patients who died or were discharged within the first 24 hours were excluded, follow-up started at hour 24, and staging used creatinine values from ICU admission onward. The 24-hour landmark was chosen as a compromise: it removes the early deaths that are most severely misclassified (patients dying within 24 hours generally cannot yet meet the 48-hour or 7-day criteria) while retaining the at-risk population during the acute staging window. Residual staging-related misclassification among patients dying between 24 hours and day 7 is inherent to window-based ascertainment and is additionally quantified by the time-stratified analyses.

### Covariates

**Illness severity: non-renal first-24-h SOFA.** A first-24-hour SOFA score was computed for every patient following the mimic-code first-day SOFA specification [17], comprising five non-renal components: respiration (PaO₂/FiO₂ with mechanical ventilation), coagulation (platelets), liver (total bilirubin), cardiovascular (mean arterial pressure; norepinephrine, epinephrine, dopamine, dobutamine, and vasopressin doses with dose-tiered scoring), and central nervous system (Glasgow Coma Scale). **The renal component was deliberately excluded** to avoid adjusting for a mathematical function of the exposure; this non-renal SOFA is the primary severity adjuster in both cohorts, ensuring covariate harmonisation. In MIMIC-IV, components used charted itemids (e.g., MAP itemid 220052; GCS items 220739/223901/223900); in eICU-CRD, corresponding `lab`, `nurseCharting`, `respiratoryCare`, `infusionDrug`, and `medication` fields. Components with no data in the first 24 h were scored 0 (consistent with mimic-code defaults).

**Mechanical ventilation** was defined in MIMIC-IV as any invasive ventilation, intubation, or ventilator-care procedure event *starting within the first 24 ICU hours* of the index stay (`procedureevents`), and in eICU-CRD as ventilated status on ICU day 1 (`apachePredVar` ventday1).

**Comorbidities** (diabetes, hypertension, heart failure, COPD, liver disease) were identified in MIMIC-IV from index-admission ICD-9/ICD-10 codes using Elixhauser-style definitions (diabetes 250.x/E10–E11; hypertension 401–405.x/I10–I15; heart failure 428.x/I50.x; COPD 490–496.x/J44.x; liver disease 570–572.x/K70–K74.x) [13,14] and in eICU-CRD from the structured `pastHistory` table using equivalent diagnostic categories. **Prior chronic kidney disease (CKD)** was identified from diagnoses dated strictly *before* the index admission (ICD-9 585.3–585.5 / ICD-10 N18.3–N18.5 in MIMIC-IV; ICD codes plus pastHistory creatinine entries in eICU-CRD) to avoid conditioning on a downstream consequence of AKI. Age and sex were included as demographics.

### Statistical Analysis

**Kaplan-Meier analysis.** Thirty-day survival from the 24-hour landmark was estimated by KDIGO stage in each cohort, with overall and Bonferroni-corrected pairwise log-rank tests (adjusted α = 0.05/6 = 0.0083).

**Cox proportional hazards models.** Three incremental models were fitted in each cohort: **Model A** (KDIGO stage, reference Stage 0); **Model B** (+ age, sex); **Model C** (+ non-renal SOFA, first-24-h ventilation, five comorbidities, prior CKD).

**Proportional hazards testing.** The PH assumption for Model C was tested using Spearman rank correlations between scaled Schoenfeld residuals and event time, with p < 0.001 (conservative threshold given the very large samples) flagging a violation.

**Time-stratified Cox models.** To characterise hazard time-dependence, Model C was refitted within three pre-specified intervals: 0–7, 7–14, and 14–30 days after the landmark, with at-risk individuals entering each interval.

**External validation (locked model).** The MIMIC-IV Model C was locked (all coefficients and the baseline hazard frozen, with event-time sanitisation identical to the derivation analysis) and applied to the eICU-CRD landmark cohort **without any refitting**. The linear predictor was transported to eICU-CRD and evaluated for: Harrell's C-index with bootstrap 95% CI (B = 200); time-dependent AUC at 7, 14, and 28 days after the landmark; Brier score at 28 days after the landmark (end of follow-up); calibration slope and intercept (Cox recalibration framework); a Hosmer-Lemeshow-type statistic across linear-predictor deciles; Kaplan-Meier mortality by linear-predictor tertiles; and decision curve analysis. **Cross-cohort parameter consistency** was assessed by refitting Model C in eICU-CRD with the identical covariate set and computing the ratio of eICU to MIMIC hazard ratios for each covariate (ratio ≈ 1 indicates parameter transportability). As a pre-specified sensitivity analysis, Model C was also refitted excluding mechanical ventilation (the covariate with the greatest cross-cohort heterogeneity) to confirm the stability of KDIGO estimates.

**Sensitivity and supplementary analyses.** Four pre-specified analyses were performed in both cohorts (Supplementary report, Tables S2–S4, Figures S2–S3): (i) restricting to patients with a pre-admission baseline creatinine; (ii) re-staging after imputing the baseline creatinine from an assumed eGFR of 75 mL/min/1.73 m² (CKD-EPI 2021 race-free equation) in patients lacking a pre-admission baseline; (iii) subgroup analyses of the Stage 3 hazard (age ≥ 65, sex, prior CKD, diabetes, ventilation, non-renal SOFA tertiles) with formal interaction tests; and (iv) an AKI phenotype analysis classifying each Stage ≥ 1 patient as **transient** (SCr falling below 1.5× baseline and below baseline + 0.3 mg/dL within 72 h of AKI onset) or **persistent** (otherwise, including RRT within 7 days or death before resolution), modelled in place of KDIGO stage with the same covariate set.

All tests were two-sided. p-values below 0.001 are reported as p < 0.001. Analyses used Python 3.13.12 with pandas 2.3.3, lifelines 0.30.3, statsmodels 0.14.6, scipy 1.18.0, scikit-learn, and matplotlib 3.11.0. All code is publicly available (see Data Availability). During manuscript preparation, large language models were used for language editing and copyediting only; all scientific content, analyses, and conclusions were generated and verified by the authors, who take full responsibility for the manuscript.

---

## Results

### Cohort Assembly

In MIMIC-IV, 94,458 ICU stays among 65,366 patients were screened. After excluding 3,991 stays with ICU length of stay <12 h, retaining the first eligible stay per patient (26,980 repeat stays excluded), excluding 2,083 patients with ESRD/chronic dialysis, and requiring a baseline creatinine and at least one SCr within the staging window (898 excluded: 425 without baseline, 473 without in-window SCr), **60,506 patients** formed the derivation cohort. In eICU-CRD, 200,859 unit stays among 139,367 patients were screened; after age exclusions (625), length-of-stay exclusion (27,842), patient-level deduplication (42,031), dual-source ESRD exclusion (6,523; ICD-identified n = 2,801, pastHistory-identified n = 6,219, with overlap), and requiring baseline and in-window creatinine (10,372 excluded), **113,466 patients** formed the validation cohort (Figure 1).

AKI prevalence differed substantially between databases: 28.9% in MIMIC-IV (Stage 1/2/3 = 11,882/2,249/3,364) versus 16.1% in eICU-CRD (11,999/1,874/4,390). This asymmetry was largely attributable to baseline creatinine availability: 80.2% of MIMIC-IV patients had a pre-admission baseline creatinine (outpatient or any-setting values within 365 days, or values within 7 days before admission), whereas only 39.4% of eICU-CRD patients had a pre-ICU baseline — 60.6% relied on the first ICU-period creatinine, which inflates the baseline and suppresses ratio-defined AKI (Table S1).

After the 24-hour landmark, **58,554 MIMIC-IV patients (5,318 deaths within 30 days)** and **107,758 eICU-CRD patients (8,616 deaths)** were analysed.

### Baseline Characteristics by KDIGO Stage (Table 1)

In both cohorts, higher KDIGO stage was associated with older age (except Stage 3), more mechanical ventilation, and higher non-renal SOFA (MIMIC-IV: 4.23/6.82/7.71/9.05 across Stages 0–3; eICU-CRD: 2.87/4.66/5.62/5.68). Liver disease in MIMIC-IV rose steeply with stage (6.4% → 39.6%), consistent with hepatorenal contributions to creatinine rise. Prior CKD was most frequent at Stage 3 (MIMIC-IV 8.8%; eICU 7.2%).

### Crude Mortality and Kaplan-Meier Analysis (Figure 2)

Crude hospital mortality increased monotonically with KDIGO stage in MIMIC-IV (5.7%/13.5%/27.6%/39.7%) and in eICU-CRD (5.4%/19.6%/33.9%/31.7%); in eICU-CRD, Stage 3 hospital mortality did not exceed Stage 2. In the landmark cohorts, 30-day mortality was 4.7%/11.8%/24.2%/36.0% (MIMIC-IV) and 5.0%/18.0%/31.7%/30.2% (eICU-CRD). Kaplan-Meier curves separated cleanly and progressively across all four stages in both cohorts (overall log-rank χ² = 1,826.3 and 2,675.1, both p < 0.001); all six pairwise comparisons remained significant after Bonferroni correction in both cohorts, including Stage 2 vs Stage 3 (MIMIC-IV p < 0.001; eICU-CRD adjusted p = 0.040).

### Cox Proportional Hazards Models (Table 2)

In MIMIC-IV, the unadjusted Stage 3 hazard ratio (Model A: HR 4.09, 95% CI 3.81–4.39) attenuated by 40% after full adjustment (Model C: HR 2.45, 2.26–2.66), with intermediate models confirming that demographic adjustment alone accounted for little of this (Model B: HR 4.29). Non-renal SOFA was the dominant confounder (Model C HR per point 1.136, 1.126–1.146). The monotonic exposure-response was preserved after full adjustment (Stage 1: 1.34, 1.25–1.43; Stage 2: 2.03, 1.84–2.24; Stage 3: 2.45, 2.26–2.66). Model C discrimination was 0.753 in MIMIC-IV. In a sensitivity model excluding mechanical ventilation, KDIGO estimates were essentially unchanged (Stage 3 HR 2.61, 2.40–2.83; C-index 0.749), confirming that the KDIGO coefficients are not driven by the ventilation term.

In eICU-CRD, Model A hazard ratios were 2.48 (2.35–2.61), 3.85 (3.53–4.19), and 3.25 (3.06–3.46) for Stages 1–3; notably, crude Stage 3 hazards did not exceed Stage 2 in this cohort. With the covariate set harmonised to MIMIC-IV (including non-renal SOFA), the fully adjusted eICU model yielded Stage 3 HR 2.46 (2.31–2.63) — nearly identical to MIMIC-IV. An eICU sensitivity model substituting APACHE IVa for non-renal SOFA (available for 91,530 patients) further attenuated the Stage 3 HR to 1.72 (1.60–1.84) with C-index 0.805, consistent with APACHE's renal component partially absorbing the AKI signal.

Other Model C covariates in MIMIC-IV behaved as expected (age HR 1.034 per year; liver disease HR 1.50; diabetes HR 0.89; hypertension HR 0.81); prior CKD was not independently associated with mortality (HR 1.06, p = 0.42). First-24-h mechanical ventilation showed an apparently protective adjusted association in MIMIC-IV (HR 0.60, 0.56–0.64) but a harmful one in eICU-CRD (HR 1.47, 1.40–1.54); this heterogeneity is examined below and in the Discussion.

### Proportional Hazards and Time-Stratified Analysis (Table 3, Figure 3)

In MIMIC-IV Model C, no covariate violated the proportional hazards assumption at p < 0.001 (all KDIGO stage p ≥ 0.29). In eICU-CRD, only age showed a mild violation (ρ = −0.038, p = 0.0008); all KDIGO terms passed (p ≥ 0.068).

Nevertheless, time-stratified models revealed a strong time-structure in the KDIGO-mortality association. In MIMIC-IV, the Stage 3 hazard ratio declined from **5.25 (4.73–5.84)** during days 0–7, to **2.59 (2.21–3.05)** during days 7–14, to **1.30 (1.09–1.54)** during days 14–30. In eICU-CRD the same pattern held: Stage 3 HR 2.67 (2.45–2.91) at 0–7 d, 2.22 (1.91–2.59) at 7–14 d, and 1.10 (0.90–1.33, p = 0.35) at 14–30 d. Thus the excess mortality risk of KDIGO Stage 3 is concentrated in the first week after staging; among 2-week survivors, Stage 3 retains only a modest (MIMIC-IV) or no statistically detectable (eICU-CRD) excess hazard.

### External Validation of the Locked MIMIC-IV Model (Table 4, Figure 4)

The locked 13-covariate MIMIC-IV Model C was applied to the eICU-CRD landmark cohort without refitting:

**Discrimination.** The C-index was 0.716 (95% CI 0.708–0.722) in eICU-CRD versus 0.753 (0.746–0.760) in MIMIC-IV (Δ = −0.037). Time-dependent AUCs at 7/14/28 days after the landmark were 0.796/0.800/0.799 in derivation and 0.750/0.748/0.746 in validation.

**Calibration.** The calibration slope was 0.957 (slightly below unity; eICU effects marginally weaker than the transported weights implied) with intercept −2.22. The Hosmer-Lemeshow-type statistic was χ² = 867 (p < 0.001; expected given n > 100,000); across linear-predictor deciles, observed 30-day mortality tracked predictions within roughly 4–8 percentage points (e.g., lowest decile 12.7% observed vs 9.1% predicted; highest decile 61.8% observed vs 65.7% predicted). The Brier score at 28 days after the landmark was 0.0710 (eICU-CRD) versus 0.0802 (MIMIC-IV).

**Risk stratification.** Kaplan-Meier 30-day mortality by linear-predictor tertiles in eICU-CRD was **19.9% / 32.1% / 52.4%** (log-rank χ² = 2,921, p < 0.001).

**Decision curve analysis.** At a low treatment threshold (5%), the locked model yielded higher net benefit than treating all patients (0.034 vs 0.032); at thresholds ≥10%, neither the model nor treating all patients outperformed treating none, given the overall low 30-day mortality.

**Cross-cohort parameter consistency.** Hazard ratio ratios (eICU/MIMIC) were 1.00 for Stage 3, 0.99 for age, 1.01 for non-renal SOFA, and 0.95–1.11 for all comorbidities — indicating excellent transportability of the core parameter set. Two exceptions were notable: Stage 1 and 2 hazards were 31–40% higher in eICU-CRD (ratios 1.40 and 1.31), and **first-24-h mechanical ventilation reversed direction** (MIMIC-IV HR 0.60, 0.56–0.64 vs eICU-CRD HR 1.47, 1.40–1.54; ratio 2.46). This reversal is examined in the Discussion.

### Sensitivity, Subgroup, and AKI Phenotype Analyses (Supplementary)

Restricting both cohorts to patients with a pre-admission baseline creatinine (46,901 MIMIC-IV; 42,868 eICU-CRD) left the Stage 3 hazard materially unchanged (2.38, 2.17–2.62 and 2.70, 2.46–2.97, respectively). Imputing the baseline from an assumed eGFR of 75 mL/min/1.73 m² (mean imputed baseline 1.00 mg/dL) and re-staging raised AKI prevalence from 28.9% to 30.3% (MIMIC-IV) and from 16.1% to 25.1% (eICU-CRD), and the Stage 3 hazard became virtually identical across cohorts (2.52, 2.32–2.73 vs 2.50, 2.36–2.66) — directly supporting the interpretation that the prevalence asymmetry between databases reflects baseline creatinine availability rather than true differences in AKI incidence.

In subgroup analyses, the Stage 3 hazard ratio remained between 1.96 and 2.91 in every pre-specified subgroup of both cohorts, with no direction reversal. The only interaction replicated in both databases was sex (higher Stage 3 hazard in males: MIMIC-IV 2.70 vs 2.15, p(interaction) = 0.0005; eICU-CRD 2.62 vs 2.29, p = 0.030); eICU-CRD interactions with ventilation, diabetes, and high non-renal SOFA were statistically significant but modest, and did not replicate as significant in MIMIC-IV (Figure S2, Table S3).

In the phenotype analysis, 56.0% of MIMIC-IV and 44.6% of eICU-CRD AKI patients had transient AKI (SCr resolution within 72 h of onset). Thirty-day mortality was 9.4% (transient) versus 32.1% (persistent) in MIMIC-IV and 12.7% versus 35.8% in eICU-CRD. In fully adjusted Cox models replacing stage with phenotype, **transient AKI carried a near-null adjusted hazard in MIMIC-IV (HR 1.07, 0.99–1.15) and a small hazard in eICU-CRD (1.24, 1.16–1.32), whereas persistent AKI carried a 2.63-fold (2.45–2.81) and 3.24-fold (3.08–3.40) hazard, respectively**; the phenotype model discriminated better than creatinine stage alone (C-index 0.768 vs 0.753 in MIMIC-IV; 0.789 vs 0.716 in eICU-CRD) (Figure S3, Table S4).

---

## Discussion

### Principal Findings

In two large, independent, patient-level ICU cohorts totalling 173,972 patients, strictly windowed KDIGO creatinine staging showed a robust, monotonic association with 30-day mortality that survived conservative severity adjustment using a deliberately non-renal SOFA score. The locked derivation model transported to a multi-centre database with a modest 0.037 decrement in C-index and virtually identical Stage 3 hazard ratios (2.45 vs 2.46). The association was, however, strongly time-dependent: the Stage 3 excess hazard was concentrated in the first week and largely dissipated among 2-week survivors.

### Strength of the Association and Severity Confounding

Our fully adjusted Stage 3 hazard ratio of 2.45 (MIMIC-IV) and 2.46 (eICU-CRD) is more conservative than the crude 4.09, quantifying the confounding contribution of non-renal organ failure at roughly two-fifths of the crude effect — directly addressing the correlation-versus-causation question raised by Girling et al. [24]. Importantly, we adjusted for a **non-renal** SOFA score. The eICU sensitivity model using APACHE IVa — which contains renal and creatinine-based elements — attenuated the Stage 3 HR to 1.72, a useful empirical demonstration that severity scores embedding renal function absorb exposure information and should not be used unmodified when quantifying the independent prognostic value of AKI staging.

Our Stage 2–3 separation contrasts with earlier reports of stage convergence [21]. Two design choices likely contribute: strict temporal windows prevent slow drifts (common in prolonged stays) from accumulating into Stage 2–3, and the 24-hour landmark removes early deaths that would otherwise be misclassified as Stage 0. Long et al. showed that mild, ratio-defined Stage 1 AKI behaves differently from absolute-criterion AKI [31]; our window-based implementation respects precisely this distinction.

### Time-Dependence of AKI Hazards

The marked attenuation of the Stage 3 hazard from 5.25 (days 0–7) to 1.30 (days 14–30) in MIMIC-IV — replicated in direction in eICU-CRD — has two practical implications. First, single averaged hazard ratios reported by most registry studies overstate the late-phase risk of AKI stages; clinicians interpreting a Stage 3 diagnosis in a patient surviving to day 14 should weigh a much smaller residual excess hazard. Second, this structure justifies the time-stratified analytic approach advocated by Singh et al. after PH violations [28]. In our landmark design, formal PH tests passed in both cohorts, yet the time-stratified models still revealed substantive effect modification by time — a reminder that passing a global PH test at conservative thresholds does not preclude clinically meaningful hazard decay.

### AKI Duration Matters More Than Peak Stage

The phenotype analysis provides the study's most clinically actionable refinement: among otherwise comparable patients, **transient AKI — resolution of creatinine criteria within 72 h of onset — carried essentially no excess 30-day mortality hazard in MIMIC-IV (HR 1.07) and only a small hazard in eICU-CRD (1.24), whereas AKI persisting beyond 72 h carried a 2.6–3.2-fold hazard, and the phenotype model discriminated better than creatinine stage itself**. This replicates, within a strictly windowed KDIGO framework and across two independent databases, the transient/persistent distinction reported in smaller single-centre cohorts, and refracts the staging–mortality association in a form closer to what clinicians actually need: a Stage 1 value that resolves by day 3 is prognostically closer to no AKI than to Stage 3. Serial creatinine trajectory, not the peak stage snapshot, carries the prognostic information.

### External Validation: What Transports and What Does Not

To our knowledge, this is among the first KDIGO staging studies to perform genuine locked-model external validation across independent critical care databases in the TRIPOD spirit, rather than refitting in the validation cohort. Discrimination transported acceptably (C-index 0.753 → 0.716; time-dependent AUCs within 0.05), the calibration slope was close to unity (0.96), and — most tellingly — the Stage 3 hazard ratio was virtually identical across cohorts (2.45 vs 2.46, ratio 1.00), as were age and non-renal SOFA weights. This parameter-level stability is stronger evidence of transportability than any single summary statistic, and supports the biological gradient interpretation of the KDIGO-mortality association [24].

Two parameters did not transport. Stage 1 and 2 hazards were 31–40% higher in eICU-CRD, plausibly reflecting the baseline creatinine availability asymmetry (60.6% of eICU baselines from the first ICU period): a higher baseline suppresses ratio-based staging, so eICU Stage 1–2 patients represent more truly injured patients relative to MIMIC-IV. This interpretation is supported by the corresponding mortality data (eICU Stage 1 mortality 18.0% vs MIMIC-IV 11.8%).

The **mechanical ventilation reversal** (HR 0.60 apparently protective in MIMIC-IV vs 1.47 harmful in eICU-CRD, among 24-hour landmark survivors, both adjusted for non-renal SOFA) warrants explicit caution. Although both variables capture ventilation on ICU day 1, they operationalise it differently — delivered invasive ventilation or intubation procedures in MIMIC-IV versus a ventilated-status flag in eICU-CRD — and the databases differ in hospital mix and ventilation practice. Within a 24-hour landmark cohort, early ventilation is strongly conditioned on surviving the initial period of support; patients ventilated on day 1 who survive to the landmark are a selected, successfully supported group, and residual confounding by indication in either direction cannot be excluded. Crucially, the KDIGO estimates are not affected by this instability: excluding ventilation from Model C shifted the Stage 3 HR only from 2.45 to 2.61. We retain ventilation in the locked model for covariate-set fidelity, flag its coefficient as non-transportable, and recommend cohort-specific recalibration of any transported model that includes it [25,28].

### Comparison with Prior Literature

The AKI-EPI study reported Stage 3 mortality of 50.0% [2]; Hoste et al. reported RIFLE Failure mortality of 26.3% [33]; FINNAKI reported KDIGO Stage 3 90-day mortality of 39.0% [34]. Our 30-day Stage 3 mortality (36.0% MIMIC-IV; 30.2% eICU-CRD) and adjusted HRs sit within this range, as expected for creatinine-only staging (urine output criteria add substantial, mostly earlier-onset, AKI [16,28]). Our lower AKI prevalence versus AKI-EPI reflects both the creatinine-only approach and our exclusion of ESRD, and the MIMIC-IV/eICU-CRD difference (28.9% vs 16.1%) is a database property — chiefly baseline creatinine availability — not a biological finding, and itself illustrates how implementation choices drive apparent epidemiology of AKI [30,32].

### Limitations

First, staging was creatinine-based only; urine output criteria were not applied, so our stages represent the creatinine-defined subset of AKI, generally of later onset [16,28]. Second, baseline creatinine was unmeasurable for a majority of eICU-CRD patients; our hierarchical definition mitigates but cannot eliminate misclassification, and prior CKD identification relies on history/coding with limited sensitivity [30,32]. Third, creatinine is affected by muscle mass, fluid balance, and haemodilution; sarcopenic patients may be understaged, and fluid resuscitation may mask injury — residual confounding by these mechanisms cannot be excluded [29]. Fourth, the two databases differ in era (2008–2022 vs 2014–2015), hospital mix, and ventilation/severity ascertainment, as discussed above; we quantified rather than eliminated these differences. Fifth, we lacked post-discharge vital status, restricting inference to 30-day in-hospital mortality. Sixth, although we adjusted for measured confounders and used a non-renal severity score, unmeasured confounding (e.g., indication for RRT, nephrotoxin exposure) persists; the E-value for the fully adjusted Stage 3 hazard ratio (HR 2.45) is approximately 4.3, meaning an unmeasured confounder would need to be associated with both Stage 3 AKI and death by a risk ratio of at least this magnitude, above and beyond the measured covariates, to explain away the finding [22].

### Conclusions

KDIGO creatinine staging, implemented with its mandated temporal windows in patient-level cohorts and adjusted conservatively for non-renal illness severity, retains a graded, monotonic association with 30-day ICU mortality that is highly reproducible across two independent databases — with the caveat that the excess hazard of severe AKI is concentrated in the first week and, more decisively, in AKI that persists beyond 72 h: transient AKI carries little if any excess mortality risk, and creatinine trajectory outperforms peak stage prognostically. A locked MIMIC-IV model transported to a multi-centre cohort with acceptable discrimination and calibration, but individual parameter transportability was heterogeneous, with mechanical ventilation showing cohort-dependent direction. Future dual-database studies should implement staging criteria exactly as specified, use non-renal severity adjustment, incorporate AKI duration alongside peak stage, and validate locked rather than refitted models.

---

## Tables

### Table 1. Baseline characteristics by KDIGO stage, derivation (MIMIC-IV) and validation (eICU-CRD) cohorts

| Characteristic | MIMIC-IV Stage 0 (n=43,011) | Stage 1 (n=11,882) | Stage 2 (n=2,249) | Stage 3 (n=3,364) | eICU Stage 0 (n=95,203) | Stage 1 (n=11,999) | Stage 2 (n=1,874) | Stage 3 (n=4,390) |
|---|---|---|---|---|---|---|---|---|
| Age, years, mean | 63.5 | 68.8 | 66.6 | 63.5 | 63.0 | 67.1 | 65.7 | 62.7 |
| Male sex, % | 54.4 | 61.6 | 53.8 | 61.4 | 53.7 | 56.5 | 49.4 | 61.1 |
| Non-renal SOFA (first 24 h), mean | 4.23 | 6.82 | 7.71 | 9.05 | 2.87 | 4.66 | 5.62 | 5.68 |
| Mechanical ventilation (first 24 h), % | 27.6 | 45.7 | 41.9 | 44.0 | 20.2 | 37.8 | 43.1 | 41.1 |
| Diabetes, % | 23.7 | 35.6 | 35.0 | 35.8 | 25.1 | 33.3 | 28.3 | 35.6 |
| Hypertension, % | 59.0 | 73.3 | 67.4 | 65.1 | 48.1 | 58.3 | 54.3 | 54.2 |
| Heart failure, % | 17.3 | 35.3 | 32.3 | 32.9 | 11.2 | 19.9 | 16.3 | 16.8 |
| COPD, % | 15.8 | 18.2 | 18.1 | 16.6 | 13.3 | 16.1 | 16.1 | 12.1 |
| Liver disease, % | 6.4 | 11.5 | 24.0 | 39.6 | 2.1 | 3.3 | 3.5 | 5.4 |
| Prior CKD, % | 2.0 | 5.5 | 3.9 | 8.8 | 1.5 | 3.6 | 1.6 | 7.2 |
| Hospital mortality, % | 5.7 | 13.5 | 27.6 | 39.7 | 5.4 | 19.6 | 33.9 | 31.7 |
| 30-day mortality (landmark cohort), % | 4.7 | 11.8 | 24.2 | 36.0 | 5.0 | 18.0 | 31.7 | 30.2 |

CKD, chronic kidney disease; COPD, chronic obstructive pulmonary disease; SOFA, Sequential Organ Failure Assessment. 30-day mortality computed in the 24-hour landmark cohorts (MIMIC-IV n = 58,554; eICU-CRD n = 107,758).

### Table 2. Cox proportional hazards models for 30-day mortality by cohort

| Variable | MIMIC-IV Model A HR (95% CI) | Model B HR (95% CI) | Model C HR (95% CI) | eICU-CRD Model A HR (95% CI) | Model B HR (95% CI) |
|---|---|---|---|---|---|
| KDIGO Stage 1 | 1.77 (1.65–1.89) | 1.63 (1.53–1.75) | 1.34 (1.25–1.43) | 2.48 (2.35–2.61) | 2.35 (2.23–2.47) |
| KDIGO Stage 2 | 3.01 (2.74–3.31) | 2.90 (2.64–3.19) | 2.03 (1.84–2.24) | 3.85 (3.53–4.19) | 3.76 (3.45–4.09) |
| KDIGO Stage 3 | 4.09 (3.81–4.39) | 4.29 (3.99–4.61) | 2.45 (2.26–2.66) | 3.25 (3.06–3.46) | 3.45 (3.24–3.68) |
| Age (per year) | — | 1.026 (1.024–1.028) | 1.034 (1.032–1.036) | — | 1.024 (1.022–1.025) |
| Male sex | — | 0.884 (0.837–0.933) | 0.856 (0.811–0.904) | — | 0.974 (0.933–1.016) |
| Non-renal SOFA (per point) | — | — | 1.136 (1.126–1.146) | — | — |
| Ventilation (first 24 h) | — | — | 0.597 (0.555–0.642) | — | — |
| Prior CKD | — | — | 1.055 (0.928–1.199) | — | — |
| Diabetes | — | — | 0.888 (0.836–0.945) | — | — |
| Hypertension | — | — | 0.814 (0.766–0.865) | — | — |
| Heart failure | — | — | 1.006 (0.946–1.069) | — | — |
| COPD | — | — | 1.121 (1.047–1.200) | — | — |
| Liver disease | — | — | 1.495 (1.392–1.606) | — | — |
| C-index | 0.653 | 0.705 | 0.753 | 0.637 | 0.690 |

Model A: KDIGO stage only. Model B: + age, sex. Model C: + non-renal SOFA, ventilation, comorbidities, prior CKD (MIMIC-IV n = 58,554). eICU-CRD Models A/B: n = 107,758. The eICU-CRD fully adjusted model with the identical covariate set (presented within the external validation framework) yielded HR 1.86 (1.77–1.96), 2.67 (2.45–2.92), and 2.46 (2.31–2.63) for Stages 1–3; an APACHE IVa-based sensitivity model (n = 91,530) yielded Stage 3 HR 1.72 (1.60–1.84), C-index 0.805. A MIMIC-IV sensitivity model excluding ventilation yielded Stage 3 HR 2.61 (2.40–2.83), C-index 0.749. All KDIGO coefficients p < 0.001.

### Table 3. Time-stratified Cox models (fully adjusted), hazard ratios by interval after 24-hour landmark

| Interval | MIMIC-IV n (events) | Stage 1 HR (95% CI) | Stage 2 HR (95% CI) | Stage 3 HR (95% CI) | eICU-CRD n (events) | Stage 1 HR (95% CI) | Stage 2 HR (95% CI) | Stage 3 HR (95% CI) |
|---|---|---|---|---|---|---|---|---|
| Days 0–7 | 37,612 (3,131) | 1.59 (1.46–1.74) | 3.89 (3.44–4.40) | 5.25 (4.73–5.84) | 67,042 (5,459) | 1.97 (1.84–2.10) | 3.42 (3.08–3.81) | 2.67 (2.45–2.91) |
| Days 7–14 | 12,223 (1,176) | 1.33 (1.16–1.52) | 1.81 (1.48–2.22) | 2.59 (2.21–3.05) | 16,630 (1,411) | 1.56 (1.39–1.77) | 2.20 (1.76–2.75) | 2.22 (1.91–2.59) |
| Days 14–30 | 8,719 (1,011) | 1.15 (1.00–1.33) | 1.14 (0.92–1.43) | 1.30 (1.09–1.54) | 7,846 (813) | 1.25 (1.06–1.47) | 1.56 (1.19–2.03) | 1.10 (0.90–1.33) |

MIMIC-IV Model C proportional hazards tests: all covariates p ≥ 0.046 (Schoenfeld residual Spearman correlations; no violation at the pre-specified p < 0.001 threshold). eICU-CRD: only age violated (ρ = −0.038, p = 0.0008); all KDIGO terms p ≥ 0.068. MIMIC-IV interval models include the full 13-covariate set (with light L2 penalisation, λ = 0.01, for interval stability); eICU-CRD interval models include age, sex, APACHE score, ventilation, and prior CKD.

### Table 4. External validation of the locked MIMIC-IV Model C in eICU-CRD (no refitting)

| Metric | MIMIC-IV (derivation) | eICU-CRD (validation) |
|---|---|---|
| N (deaths) | 58,554 (5,318) | 107,758 (8,616) |
| C-index (95% CI) | 0.753 (0.746–0.760) | 0.716 (0.708–0.722) |
| Time-dependent AUC, 7 d after landmark | 0.796 | 0.750 |
| Time-dependent AUC, 14 d after landmark | 0.800 | 0.748 |
| Time-dependent AUC, 28 d after landmark | 0.799 | 0.746 |
| Brier score, 28 d after landmark | 0.0802 | 0.0710 |
| Calibration slope | 1.0 (by construction) | 0.957 |
| Calibration intercept | 0 | −2.22 |
| Hosmer-Lemeshow-type χ² (LP deciles) | — | 867 (p < 0.001) |
| 30-d mortality by risk tertile | — | 19.9% / 32.1% / 52.4% (log-rank p < 0.001) |

Cross-cohort hazard-ratio ratios (eICU/MIMIC): Stage 3 1.00; age 0.99; non-renal SOFA 1.01; liver disease 0.95; Stage 1 1.40; Stage 2 1.31; ventilation 2.46 (direction reversal).

---

## Figure Legends

**Figure 1.** Cohort flow diagram for MIMIC-IV (derivation) and eICU-CRD (validation), showing patient-level deduplication, end-stage renal disease exclusion, baseline creatinine requirements, and 24-hour landmark exclusions.

**Figure 2.** Kaplan-Meier 30-day survival by KDIGO stage (24-hour landmark). (A) MIMIC-IV (overall log-rank χ² = 1,826.3, p < 0.001). (B) eICU-CRD (χ² = 2,675.1, p < 0.001). All pairwise comparisons significant after Bonferroni correction.

**Figure 3.** Forest plot of time-stratified, fully adjusted Cox hazard ratios for KDIGO stages by interval (days 0–7, 7–14, 14–30) in both cohorts, demonstrating early hazard concentration and subsequent attenuation.

**Figure 4.** External validation of the locked MIMIC-IV model in eICU-CRD. (A) C-index derivation vs validation. (B) Calibration by predicted-risk decile (observed vs predicted 30-day mortality). (C) Kaplan-Meier curves by linear-predictor tertile. (D) Decision curve analysis.

**Figure S1 (Supplementary).** Directed acyclic graph of the assumed causal structure, indicating non-renal SOFA, ventilation, comorbidities, age, and sex as measured confounders and the renal SOFA component as part of the exposure-outcome pathway.

**Figure S2 (Supplementary).** Forest plot of fully adjusted Stage 3 vs Stage 0 hazard ratios by pre-specified subgroup (age, sex, prior CKD, diabetes, ventilation, non-renal SOFA tertiles) in MIMIC-IV and eICU-CRD, with interaction p-values.

**Figure S3 (Supplementary).** Adjusted 30-day mortality hazard ratios for transient versus persistent AKI (72 h resolution criterion), each versus no AKI, in both cohorts.

**Table S1 (Supplementary).** Baseline creatinine source distribution by cohort and KDIGO stage.

**Table S2 (Supplementary).** Sensitivity analyses: pre-admission baseline restriction and eGFR-imputed (CKD-EPI 2021, 75 mL/min/1.73 m²) baseline re-staging, with KDIGO hazard ratios and AKI prevalence in both cohorts.

**Table S3 (Supplementary).** Subgroup-specific Stage 3 hazard ratios and interaction tests in both cohorts.

**Table S4 (Supplementary).** AKI phenotype distribution, 30-day mortality, and fully adjusted Cox models for transient versus persistent AKI.

---

## Abbreviations

AKI, acute kidney injury; APACHE, Acute Physiology and Chronic Health Evaluation; CKD, chronic kidney disease; COPD, chronic obstructive pulmonary disease; DAG, directed acyclic graph; eICU-CRD, eICU Collaborative Research Database; ESRD, end-stage renal disease; GCS, Glasgow Coma Scale; HR, hazard ratio; ICU, intensive care unit; IQR, interquartile range; KDIGO, Kidney Disease: Improving Global Outcomes; MIMIC-IV, Medical Information Mart for Intensive Care IV; PH, proportional hazards; RRT, renal replacement therapy; SCr, serum creatinine; SOFA, Sequential Organ Failure Assessment; STROBE, Strengthening the Reporting of Observational Studies in Epidemiology; TRIPOD, Transparent Reporting of a multivariable prediction model for Individual Prognosis Or Diagnosis.

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

## Ethics approval and consent to participate

Access to MIMIC-IV and eICU-CRD was obtained through PhysioNet after completion of the CITI Data or Specimens Only Research course. The institutional review board at BIDMC approved the creation of the MIMIC-IV database. Use of de-identified, publicly available data did not require additional institutional review board approval, and the requirement for informed consent was waived.

---

## Consent for publication

Not applicable. This study used de-identified, publicly available data from the MIMIC-IV and eICU-CRD databases.

---

## Availability of data and materials

MIMIC-IV is publicly available at PhysioNet (https://physionet.org/content/mimiciv/). eICU-CRD is publicly available at PhysioNet (https://physionet.org/content/eicu-crd/). All analysis code (cohort construction, SOFA computation, survival analysis, and locked-model external validation) and key result files are publicly available at GitHub (https://github.com/wudengke2010/kdigo-aki-mimic-iv).

---

## Competing interests

The authors declare that they have no competing interests.

---

## Funding

This work was supported by the Chronic Disease Management Research Project of National Health Commission Capacity Building and Continuing Education Center (Grant No. GWJJMB202510024181), a grant from the Changsha Science and Technology Bureau Project (Grant No. kq2014242), and the Natural Science Foundation of Hunan Province of China (Grant No. 2021JJ30959). The funders had no role in study design, data collection and analysis, decision to publish, or preparation of the manuscript.

---

## Authors' contributions

**Jiqiang Liu:** Conceptualisation, Data curation, Formal analysis, Investigation, Methodology, Writing – original draft. **Dengke Wu:** Conceptualisation, Funding acquisition, Project administration, Supervision, Writing – review & editing. All authors read and approved the final manuscript.

---

## Acknowledgements

The authors thank the PhysioNet team and the Laboratory for Computational Physiology at the Massachusetts Institute of Technology for creating and maintaining the MIMIC-IV and eICU-CRD databases, which made this study possible. The authors also acknowledge the MIT-LCP `mimic-code` repository, whose validated first-day SOFA SQL specification informed the severity score implementation.

---

## Authors' information

Dengke Wu https://orcid.org/0000-0003-4101-8461

---

## Declaration of Generative AI in Scientific Writing

During the preparation of this work, the authors used large language models (LLMs) for language editing and copyediting. After using these tools, the authors reviewed and edited the content as needed and take full responsibility for the content of the publication.
