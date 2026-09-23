# -*- coding: utf-8 -*-
"""Generate STROBE checklist DOCX for JTM submission (manuscript_jtm.md)."""
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH

# (item, recommendation, location in manuscript)
ITEMS = [
    ("Title and Abstract", "", ""),
    ("1", "Indicate the study's design with a commonly used term in the title or in the abstract.", "Title; Abstract (Design: retrospective dual-cohort study)"),
    ("Introduction", "", ""),
    ("2", "Explain the scientific background and rationale for the investigation being reported.", "Introduction, paragraphs 1-3"),
    ("3", "State specific objectives, including any prespecified hypotheses.", "Introduction, final paragraph"),
    ("Methods", "", ""),
    ("4", "Present key elements of study design early in the paper.", "Methods, Study Design and Data Sources"),
    ("5", "Describe the setting, locations, and relevant dates, including periods of recruitment, exposure, follow-up, and data collection.", "Methods, Study Design and Data Sources (MIMIC-IV 2008-2022; eICU-CRD 2014-2015)"),
    ("6", "(a) Cohorts: Give the eligibility criteria, and the sources and methods of selection of participants.", "Methods, Cohort Construction; Results, Cohort Assembly (Figure 1)"),
    ("7", "Clearly define all outcomes, exposures, predictors, potential confounders, and effect modifiers.", "Methods, Outcomes; KDIGO Staging; Covariates; Statistical Analysis (subgroups/interactions)"),
    ("8", "For each variable of interest, give sources of data and details of methods of assessment (measurement).", "Methods, Data Sources; KDIGO Staging; Non-renal SOFA; Covariates (laboratory itemids, procedureevents, ICD codes)"),
    ("9", "Describe any efforts to address potential sources of bias.", "Methods (patient-level cohorts, temporal windows, non-renal SOFA, hierarchical baseline, 24-h landmark); Discussion, Strengths"),
    ("10", "Explain how missing data were addressed.", "Methods, Baseline Creatinine (hierarchy); Results, Cohort Assembly (exclusions for missing baseline/window creatinine); Table S1"),
    ("11", "Explain how quantitative variables were handled in the analyses. (If applicable, describe which groupings were chosen and why.)", "Methods (age per decade; non-renal SOFA continuous; SOFA tertiles in subgroup analysis)"),
    ("12", "(a) Describe all statistical methods, including those used to control for confounding; (b) Describe any methods used to examine subgroups and interactions; (c) Explain how missing data were addressed; (d) If applicable, describe analytical methods taking account of sampling strategy; (e) Describe any sensitivity analyses.", "Methods, Statistical Analysis (Cox models A/B/C, Schoenfeld PH tests, time-stratified Cox, locked-model external validation, subgroup/interaction tests, four sensitivity analyses: pre-admission baseline, eGFR-75 back-calculation, subgroups, transient vs persistent phenotype)"),
    ("Results", "", ""),
    ("13", "(a) Report numbers of individuals at each stage of study. Give reasons for non-participation at each stage. (b) Consider use of a flow diagram.", "Results, Cohort Assembly; Figure 1"),
    ("14", "(a) Give characteristics of study participants and information on exposures and potential confounders. (b) Indicate number of participants with missing data for each variable of interest.", "Table 1; Table S1"),
    ("15", "Report numbers of outcome events or summary measures over time.", "Results, Crude Mortality (stage-specific 30-day mortality); Figure 2"),
    ("16", "(a) Give unadjusted estimates and, if applicable, confounder-adjusted estimates and their precision. (b) Report category boundaries when continuous variables were categorized. (c) If relevant, consider translating estimates of relative risk into absolute risk for a meaningful time period.", "Table 2 (Models A/B/C HRs with 95% CIs); Table 3 (time-stratified HRs); Results (tertile absolute mortality 19.9/32.1/52.4%)"),
    ("17", "Report other analyses done - e.g., analyses of subgroups and interactions, and sensitivity analyses.", "Results, Sensitivity Analyses; Figure S2; Tables S2-S4; Figure S3"),
    ("Discussion", "", ""),
    ("18", "Summarize key results with reference to study objectives.", "Discussion, first paragraph (Key Findings)"),
    ("19", "Discuss limitations of the study, taking into account sources of potential bias or imprecision.", "Discussion, Strengths and Limitations"),
    ("20", "Give a cautious overall interpretation of results considering objectives, limitations, multiplicity of analyses, results from similar studies, and other relevant evidence.", "Discussion, paragraphs 2-7 (comparison with AKI-EPI, FINNAKI, Lin et al. CKJ 2025)"),
    ("21", "Discuss the generalisability (external validity) of the study results.", "Discussion, locked-model external validation paragraph (C-index transport, cross-cohort HR ratio 1.00)"),
    ("Other information", "", ""),
    ("22", "Give the source of funding and the role of the funders for the present study.", "Declarations, Funding"),
    ("23", "Describe the role of the study funders, if any, in study design; collection, analysis, and interpretation of data; writing of the report; and the decision to submit the report for publication.", "Declarations, Funding (funders had no role)"),
]

doc = Document()
# Title
h = doc.add_heading("STROBE Statement — Checklist of Items for Reports of Observational Studies", level=1)
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = p.add_run("Manuscript: Association of KDIGO creatinine-based AKI staging with 30-day mortality in critically ill patients: a dual-cohort retrospective study with phenotype analysis and locked-model external validation (MIMIC-IV v3.1 and eICU-CRD)\n")
r.italic = True
r.font.size = Pt(10)
p2 = doc.add_paragraph()
p2.alignment = WD_ALIGN_PARAGRAPH.CENTER
r2 = p2.add_run("Completed by: Jiqiang Liu, Dengke Wu   |   Date: 2026-09-23")
r2.font.size = Pt(9)

doc.add_paragraph("Note: This study additionally complies with TRIPOD guidance for the locked-model external validation component (Model C transported without refitting), reported in Results (External Validation) and Table 4.").italic = True

table = doc.add_table(rows=1, cols=4)
table.style = "Table Grid"
hdr = table.rows[0].cells
for i, t in enumerate(["Item", "STROBE Recommendation", "Location in manuscript (section/paragraph/table/figure)", "Page"]):
    hdr[i].text = ""
    rr = hdr[i].paragraphs[0].add_run(t)
    rr.bold = True
    rr.font.size = Pt(9)

for item, rec, loc in ITEMS:
    row = table.add_row().cells
    if rec == "":  # section header
        row[0].text = ""
        rr = row[0].paragraphs[0].add_run(item)
        rr.bold = True
        rr.font.size = Pt(9)
        for c in row[1:]:
            c.text = ""
    else:
        row[0].text = item
        row[1].text = rec
        row[2].text = loc
        row[3].text = ""
        for c in row:
            for par in c.paragraphs:
                for rr in par.runs:
                    rr.font.size = Pt(8)

table.columns[0].width = Inches(0.5)
table.columns[1].width = Inches(2.8)
table.columns[2].width = Inches(2.8)
table.columns[3].width = Inches(0.5)
for row in table.rows:
    row.cells[0].width = Inches(0.5)
    row.cells[1].width = Inches(2.8)
    row.cells[2].width = Inches(2.8)
    row.cells[3].width = Inches(0.5)

doc.save("jtm_submission/strobe_checklist_jtm.docx")
print("Saved jtm_submission/strobe_checklist_jtm.docx")
