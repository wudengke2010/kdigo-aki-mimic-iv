# KDIGO AKI Paradox Analysis Using MIMIC-IV

This repository contains the analysis code for the study:

**"The AKI Paradox Persists After Severity Adjustment: A SOFA-Adjusted Propensity-Matched Analysis with Prediction Modeling Using MIMIC-IV"**

## Study Overview

This study investigates the "AKI paradox" — the counterintuitive finding that chronic kidney disease (CKD) appears protective in patients with severe (Stage 3) acute kidney injury (AKI), using the MIMIC-IV v3.1 database.

### Key Findings

- In Stage 3 AKI, CKD patients had lower mortality than Pure AKI patients (35.4% vs 42.0%)
- The CKD x Stage 3 interaction remained robust after SOFA adjustment (OR 0.643), PSM (OR 0.587), and multiple imputation (OR 0.605)
- The paradox was most pronounced in the sickest patients (SOFA > 10: OR 0.679)
- The protective effect was temporally specific to acute hospitalization and not significant by 90 days

## Data Source

The MIMIC-IV v3.1 database is publicly available at https://physionet.org/content/mimiciv/3.1/ after completion of CITI Data or Specimens Only Research training and signing of the data use agreement.

Due to the data use agreement, individual patient data cannot be shared. This repository provides all analysis code for reproducibility.

## Repository Structure

```
kdigo-aki-mimic-iv/
|-- src/
|   |-- data_extraction/        # Cohort extraction and SOFA scoring
|   |   |-- extract_cohort.py   # MIMIC-IV cohort extraction (DuckDB)
|   |   |-- extract_sofa.py     # SOFA score calculation
|   |-- statistics/             # Primary and secondary analyses
|   |   |-- compute_all_stats.py    # Unified statistics computation
|   |   |-- primary_analysis.py     # Logistic regression Models A-D
|   |   |-- sofa_analysis.py        # SOFA-stratified and nonrenal SOFA analysis
|   |   |-- temporal_analysis.py    # Temporal dynamics (in-hospital to 1-year)
|   |   |-- auc_bootstrap.py        # AUC with bootstrap 95% CI
|   |-- sensitivity_analysis/   # Robustness checks
|   |   |-- mice_imputation.py      # MICE multiple imputation (m=20)
|   |   |-- evalue_fdr.py           # E-values and Benjamini-Hochberg FDR
|   |   |-- create_dag.py           # Directed Acyclic Graph (DAG)
|   |   |-- analyze_missing.py      # Missing data pattern analysis
|   |-- prediction_model/       # Nomogram and LASSO prediction
|   |   |-- nomogram_lasso.py       # LASSO feature selection + nomogram
|   |   |-- advanced_models.py      # NRI/IDI and calibration
|   |-- figures/                # Figure generation
|       |-- generate_all_figures.R  # R figures (forest, SOFA, ROC, calibration, etc.)
|       |-- generate_flow_diagram.R # CONSORT flow diagram (R)
|       |-- generate_flow_diagram.py # CONSORT flow diagram (Python)
|-- README.md
|-- requirements.txt
|-- .gitignore
```

## Software Requirements

### Python 3.10+
```
pandas>=2.0
numpy>=1.24
scipy>=1.10
scikit-learn>=1.3
statsmodels>=0.14
lifelines>=0.27
duckdb>=0.9
matplotlib>=3.7
```

### R 4.6+
```
ggplot2 (>= 4.0)
survminer (>= 0.5)
survival (>= 3.8)
rms (>= 8.1)
forestplot (>= 3.2)
cmprsk (>= 2.2)
ggrepel (>= 0.9)
patchwork (>= 1.3)
cowplot (>= 1.2)
scales (>= 1.4)
```

## How to Reproduce

1. **Obtain MIMIC-IV access**: Complete CITI training and request access at https://physionet.org/content/mimiciv/3.1/

2. **Extract the cohort**: Modify `src/data_extraction/extract_cohort.py` to point to your MIMIC-IV data path, then run:
   ```bash
   python src/data_extraction/extract_cohort.py
   ```

3. **Compute SOFA scores**: Run `src/data_extraction/extract_sofa.py` on the extracted cohort.

4. **Primary analysis**: Run `src/statistics/compute_all_stats.py` to generate all statistics, followed by individual analysis scripts as needed.

5. **Sensitivity analyses**: Execute scripts in `src/sensitivity_analysis/` for MICE, E-values, FDR, and DAG.

6. **Generate figures**: Run `src/figures/generate_all_figures.R` in R.

## Citation

If you use this code, please cite:

> Liu J, Wu D. The AKI Paradox Persists After Severity Adjustment: A SOFA-Adjusted Propensity-Matched Analysis with Prediction Modeling Using MIMIC-IV. *Kidney International Reports* (2026) [Under review].

## License

This analysis code is released under the MIT License. The MIMIC-IV database is governed by its own data use agreement and is not included here.

## Authors

- **Jiqiang Liu** (ORCID: 0009-0000-9884-3089) - Conceptualization, Methodology, Data curation, Formal analysis, Software, Visualization, Writing - original draft
- **Dengke Wu** (ORCID: 0009-0008-1363-9621) - Conceptualization, Methodology, Supervision, Funding acquisition, Writing - review and editing

## Acknowledgments

We thank the MIMIC-IV program investigators and the PhysioNet team for making the database publicly available.
