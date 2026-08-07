#!/usr/bin/env python3
"""
Generate Directed Acyclic Graph (DAG) for KDIGO AKI study.
Output: Supplementary Figure S5 (PNG + PDF)
Tool: graphviz Python package
Reference format: DAGitty conventions
"""

from graphviz import Digraph
import os

# Create DAG
dag = Digraph('KDIGO_AKI_DAG', format='png')
dag.attr(rankdir='LR', size='10,7', dpi='300', bgcolor='white',
         fontname='Helvetica', labelloc='t', fontsize='14',
         label='Supplementary Figure S5. Directed Acyclic Graph (DAG) for KDIGO AKI Study')
dag.attr('node', fontname='Helvetica', fontsize='10', style='filled')
dag.attr('edge', fontname='Helvetica', fontsize='9', color='#333333')

# ============================
# EXPOSURE (red)
# ============================
dag.node('KDIGO', 'KDIGO Stage\n(Exposure)',
         shape='box', fillcolor='#FFB3B3', color='#CC0000', penwidth='2')

# ============================
# OUTCOME (blue)
# ============================
dag.node('Mortality', 'ICU Mortality\n(Outcome)',
         shape='box', fillcolor='#B3D9FF', color='#0066CC', penwidth='2')

# ============================
# MEASURED CONFOUNDERS (green) - adjusted in Model B
# ============================
dag.node('Age', 'Age', shape='ellipse', fillcolor='#B3FFB3', color='#009900')
dag.node('Sex', 'Sex', shape='ellipse', fillcolor='#B3FFB3', color='#009900')
dag.node('Comorb', 'Comorbidities\n(CHF, CKD, Liver,\nCancer, Sepsis, etc.)',
         shape='ellipse', fillcolor='#B3FFB3', color='#009900')
dag.node('SOFA', 'SOFA Score\n(Illness Severity)',
         shape='ellipse', fillcolor='#B3FFB3', color='#009900')
dag.node('Vent', 'Mechanical\nVentilation',
         shape='ellipse', fillcolor='#B3FFB3', color='#009900')
dag.node('Vaso', 'Vasopressor\nUse',
         shape='ellipse', fillcolor='#B3FFB3', color='#009900')

# Confounders -> Exposure
dag.edge('Age', 'KDIGO')
dag.edge('Sex', 'KDIGO')
dag.edge('Comorb', 'KDIGO')
dag.edge('SOFA', 'KDIGO')
dag.edge('Vent', 'KDIGO')
dag.edge('Vaso', 'KDIGO')

# Confounders -> Outcome
dag.edge('Age', 'Mortality')
dag.edge('Sex', 'Mortality')
dag.edge('Comorb', 'Mortality')
dag.edge('SOFA', 'Mortality')
dag.edge('Vent', 'Mortality')
dag.edge('Vaso', 'Mortality')

# ============================
# UNMEASURED CONFOUNDERS (orange, dashed)
# ============================
dag.node('Unmeasured', 'Unmeasured\nConfounders\n(genetics, SES,\npre-ICU health)',
         shape='ellipse', fillcolor='#FFE0B3', color='#FF8800', style='filled,dashed')
dag.edge('Unmeasured', 'KDIGO', style='dashed', color='#FF8800')
dag.edge('Unmeasured', 'Mortality', style='dashed', color='#FF8800')

# ============================
# COLLIDERS (yellow, with warning) - NOT adjusted
# ============================
dag.node('ICUAdmit', 'ICU Admission\n(Collider)\n[conditioning inherent\nto MIMIC-IV design]',
         shape='diamond', fillcolor='#FFFFB3', color='#CCCC00', penwidth='1.5')
dag.node('CrMeas', 'SCr Measurement\nAvailability\n(Collider)\n[1,075 excluded\nwithout SCr]',
         shape='diamond', fillcolor='#FFFFB3', color='#CCCC00', penwidth='1.5')

# Collider paths
dag.edge('Unmeasured', 'ICUAdmit', style='dashed', color='#FF8800')
dag.edge('Comorb', 'ICUAdmit', color='#009900')
dag.edge('SOFA', 'ICUAdmit', color='#009900')
dag.edge('ICUAdmit', 'KDIGO', color='#CCCC00')
dag.edge('ICUAdmit', 'Mortality', color='#CCCC00')

dag.edge('Unmeasured', 'CrMeas', style='dashed', color='#FF8800')
dag.edge('SOFA', 'CrMeas', color='#009900')
dag.edge('CrMeas', 'KDIGO', color='#CCCC00')

# ============================
# MEDIATOR (purple) - RRT
# ============================
dag.node('RRT', 'RRT Initiation\n(Mediator)',
         shape='ellipse', fillcolor='#E6B3FF', color='#9900CC')
dag.edge('KDIGO', 'RRT', color='#9900CC')
dag.edge('RRT', 'Mortality', color='#9900CC')
dag.edge('SOFA', 'RRT', color='#009900')

# ============================
# LEGEND
# ============================
dag.attr(label='''Supplementary Figure S5. Directed Acyclic Graph (DAG) for KDIGO AKI Study

Legend:
  Red box     = Exposure (KDIGO creatinine-based AKI staging)
  Blue box    = Outcome (ICU in-hospital mortality)
  Green ellipse = Measured confounders (adjusted in Model B)
  Orange ellipse (dashed) = Unmeasured confounders
  Yellow diamond = Collider variables (NOT adjusted; conditioning acknowledged in Limitations)
  Purple ellipse = Mediator (RRT; adjusted in Model B as part of severity)

Note: ICU admission and SCr measurement availability are structural colliders inherent to
the MIMIC-IV design. Conditioning on ICU admission is unavoidable in ICU-based cohorts.
SCr measurement availability affected 1,075 (1.3%) of 85,242 ICU stays. Neither was
included as a covariate in the regression models. RRT initiation is both a component
of KDIGO Stage 3 definition and a mediator on the AKI → mortality pathway; it was
adjusted in Model B as a severity marker and further decomposed in stratified analyses.''')

# ============================
# Render
# ============================
out_dir = os.path.join(os.path.dirname(__file__), 'figures_prism')
os.makedirs(out_dir, exist_ok=True)

# PNG
dag.render(filename='figS5_dag', directory=out_dir, format='png', cleanup=True)
print(f"-> {out_dir}/figS5_dag.png")

# PDF
dag.render(filename='figS5_dag', directory=out_dir, format='pdf', cleanup=True)
print(f"-> {out_dir}/figS5_dag.pdf")

print("\nDAG generated successfully.")
print("Nodes: 1 exposure, 1 outcome, 6 measured confounders, 1 unmeasured confounder, 2 colliders, 1 mediator")
