"""
Create Supplementary Figure S5: Directed Acyclic Graph (DAG)
Shows causal structure: CKD (exposure) -> Mortality (outcome)
with confounders, mediators, and colliders marked.
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(12, 9))
ax.set_xlim(0, 12)
ax.set_ylim(0, 9)
ax.axis('off')

# Color scheme
C_EXPOSURE = '#1565C0'   # Blue for exposure
C_OUTCOME = '#C62828'    # Red for outcome
C_CONFOUNDER = '#2E7D32' # Green for confounders
C_MEDIATOR = '#F57F17'   # Orange for mediator
C_COLLIDER = '#6A1B9A'   # Purple for collider (not adjusted)
C_ARROW_CONF = '#2E7D32'
C_ARROW_MED = '#F57F17'
C_ARROW_COLL = '#6A1B9A'

# Node positions
nodes = {
    'CKD':           (2.0, 5.0),
    'Mortality':     (10.0, 5.0),
    'Age':           (4.0, 8.0),
    'Sex':           (2.5, 8.0),
    'HTN':           (5.5, 8.0),
    'DM':            (7.0, 8.0),
    'HF':            (8.5, 8.0),
    'SOFA':          (6.0, 5.0),
    'AKI Stage':     (6.0, 2.5),
    'ICU LOS':       (4.0, 1.5),
    'AKI Timing':    (8.0, 1.5),
}

# Draw nodes
def draw_node(ax, x, y, label, color, fontsize=10, boxstyle='round,pad=0.3'):
    box = mpatches.FancyBboxPatch(
        (x - 0.65, y - 0.3), 1.3, 0.6,
        boxstyle=boxstyle,
        facecolor=color, edgecolor='white',
        alpha=0.15, linewidth=1.5
    )
    ax.add_patch(box)
    ax.text(x, y, label, ha='center', va='center',
            fontsize=fontsize, fontweight='bold', color=color)

# Draw arrows
def draw_arrow(ax, x1, y1, x2, y2, color='black', style='-', lw=1.2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color,
                               linestyle=style, lw=lw,
                               connectionstyle='arc3,rad=0'))

# Exposure and outcome
draw_node(ax, *nodes['CKD'], 'CKD\n(Exposure)', C_EXPOSURE, fontsize=11)
draw_node(ax, *nodes['Mortality'], 'Hospital\nMortality\n(Outcome)', C_OUTCOME, fontsize=11)

# Confounders (adjusted)
for name in ['Age', 'Sex', 'HTN', 'DM', 'HF']:
    draw_node(ax, *nodes[name], name, C_CONFOUNDER, fontsize=9)

# Mediator (adjusted)
draw_node(ax, *nodes['SOFA'], 'SOFA\n(Mediator)', C_MEDIATOR, fontsize=9)
draw_node(ax, *nodes['AKI Stage'], 'AKI Stage\n(Effect Modifier)', C_MEDIATOR, fontsize=9)

# Colliders (NOT adjusted)
draw_node(ax, *nodes['ICU LOS'], 'ICU LOS\n(Collider)', C_COLLIDER, fontsize=8)
draw_node(ax, *nodes['AKI Timing'], 'AKI Diagnostic\nTiming (Collider)', C_COLLIDER, fontsize=8)

# Arrows: Confounders -> CKD and -> Mortality
for name in ['Age', 'Sex', 'HTN', 'DM', 'HF']:
    x, y = nodes[name]
    # To CKD
    draw_arrow(ax, x-0.3, y-0.3, nodes['CKD'][0]+0.3, nodes['CKD'][1]+0.3,
               color=C_ARROW_CONF, lw=1.0)
    # To Mortality
    draw_arrow(ax, x+0.3, y-0.3, nodes['Mortality'][0]-0.3, nodes['Mortality'][1]+0.3,
               color=C_ARROW_CONF, lw=1.0)

# CKD -> Mortality (main effect, thick)
draw_arrow(ax, nodes['CKD'][0]+0.65, nodes['CKD'][1],
           nodes['Mortality'][0]-0.65, nodes['Mortality'][1],
           color=C_EXPOSURE, lw=2.5)

# CKD -> SOFA -> Mortality (mediation pathway)
draw_arrow(ax, nodes['CKD'][0]+0.3, nodes['CKD'][1]+0.3,
           nodes['SOFA'][0]-0.65, nodes['SOFA'][1],
           color=C_ARROW_MED, lw=1.2)
draw_arrow(ax, nodes['SOFA'][0]+0.65, nodes['SOFA'][1],
           nodes['Mortality'][0]-0.3, nodes['Mortality'][1]-0.3,
           color=C_ARROW_MED, lw=1.2)

# AKI Stage -> Mortality
draw_arrow(ax, nodes['AKI Stage'][0]+0.3, nodes['AKI Stage'][1]+0.3,
           nodes['Mortality'][0]-0.3, nodes['Mortality'][1]-0.3,
           color=C_ARROW_MED, lw=1.2)

# CKD -> AKI Stage (CKD patients may have different AKI severity)
draw_arrow(ax, nodes['CKD'][0]+0.3, nodes['CKD'][1]-0.3,
           nodes['AKI Stage'][0]-0.3, nodes['AKI Stage'][1]+0.3,
           color=C_ARROW_MED, lw=1.0)

# Colliders: ICU LOS and AKI Timing
# ICU LOS <- CKD, ICU LOS <- Mortality (collider)
draw_arrow(ax, nodes['CKD'][0]+0.3, nodes['CKD'][1]-0.3,
           nodes['ICU LOS'][0]-0.3, nodes['ICU LOS'][1]+0.3,
           color=C_ARROW_COLL, lw=0.8, style='--')
draw_arrow(ax, nodes['Mortality'][0]-0.5, nodes['Mortality'][1]-0.3,
           nodes['ICU LOS'][0]+0.3, nodes['ICU LOS'][1]+0.3,
           color=C_ARROW_COLL, lw=0.8, style='--')

# AKI Timing <- CKD, AKI Timing <- Mortality (collider)
draw_arrow(ax, nodes['AKI Stage'][0]+0.3, nodes['AKI Stage'][1]-0.3,
           nodes['AKI Timing'][0]-0.3, nodes['AKI Timing'][1]+0.3,
           color=C_ARROW_COLL, lw=0.8, style='--')
draw_arrow(ax, nodes['Mortality'][0]-0.3, nodes['Mortality'][1]-0.3,
           nodes['AKI Timing'][0]+0.3, nodes['AKI Timing'][1]+0.3,
           color=C_ARROW_COLL, lw=0.8, style='--')

# Legend
legend_items = [
    (C_EXPOSURE, 'Exposure'),
    (C_OUTCOME, 'Outcome'),
    (C_CONFOUNDER, 'Confounder (adjusted)'),
    (C_MEDIATOR, 'Mediator / Effect Modifier (adjusted)'),
    (C_COLLIDER, 'Collider (NOT adjusted)'),
]

for i, (color, label) in enumerate(legend_items):
    y_leg = 0.5 + i * 0.35
    ax.add_patch(mpatches.FancyBboxPatch(
        (0.5, y_leg - 0.12), 0.35, 0.24,
        boxstyle='round,pad=0.1',
        facecolor=color, edgecolor='white', alpha=0.3, linewidth=1
    ))
    ax.text(1.0, y_leg, label, fontsize=9, va='center', color=color, fontweight='bold')

# Arrow legend
ax.annotate('', xy=(3.5, 0.5), xytext=(2.8, 0.5),
            arrowprops=dict(arrowstyle='->', color='black', lw=1.2))
ax.text(3.6, 0.5, 'Causal path', fontsize=9, va='center')

ax.annotate('', xy=(6.0, 0.5), xytext=(5.3, 0.5),
            arrowprops=dict(arrowstyle='->', color=C_COLLIDER, lw=1.0, linestyle='--'))
ax.text(6.1, 0.5, 'Collider path (dashed)', fontsize=9, va='center', color=C_COLLIDER)

# Title
ax.text(6.0, 8.8, 'Supplementary Figure S5. Directed Acyclic Graph (DAG)',
        ha='center', va='center', fontsize=13, fontweight='bold')
ax.text(6.0, 8.45, 'Causal structure of CKD effect on hospital mortality in AKI patients',
        ha='center', va='center', fontsize=10, style='italic', color='gray')

plt.tight_layout()
out_dir = r'C:\Users\admin\WorkBuddy\2026-07-07-19-40-19\submission_figures_R'
import os
os.makedirs(out_dir, exist_ok=True)

plt.savefig(os.path.join(out_dir, 'Figure_S5_DAG.png'), dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig(os.path.join(out_dir, 'Figure_S5_DAG.pdf'), bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()

print("DAG saved to submission_figures_R/Figure_S5_DAG.png and .pdf")
