"""
Figure 1: Study Population Flow Diagram
Rendered with matplotlib to match draw.io professional quality.
300 DPI, SCI publication standard.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.font_manager as fm
import numpy as np

# ===== Global style =====
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 11,
    'axes.linewidth': 0.8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
})

# ===== Colors =====
C_BLUE = '#dae8fc'
C_BLUE_EDGE = '#6c8ebf'
C_PURPLE = '#e1d5e7'
C_PURPLE_EDGE = '#9673a6'
C_YELLOW = '#fff2cc'
C_YELLOW_EDGE = '#d6b656'
C_GREEN = '#d5e8d4'
C_GREEN_EDGE = '#82b366'
C_RED = '#f8cecc'
C_RED_EDGE = '#b85450'
C_GREY = '#f5f5f5'
C_GREY_EDGE = '#666666'
C_DARK = '#333333'

def draw_box(ax, x, y, w, h, text, fc=C_BLUE, ec=C_BLUE_EDGE, fs=11, fw='normal',
             rounded=True, lw=1.5, text_color=C_DARK, style='round'):
    """Draw a rounded rectangle with centered text."""
    if rounded:
        box = FancyBboxPatch((x, y), w, h,
                             boxstyle=f"round,pad=0.1,rounding_size=0.3",
                             facecolor=fc, edgecolor=ec, linewidth=lw, zorder=2)
    else:
        box = FancyBboxPatch((x, y), w, h,
                             boxstyle="square,pad=0.1",
                             facecolor=fc, edgecolor=ec, linewidth=lw, zorder=2)
    ax.add_patch(box)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center',
            fontsize=fs, fontweight=fw, color=text_color, zorder=3,
            linespacing=1.4)
    return (x, y, w, h)

def draw_arrow(ax, x1, y1, x2, y2, color='#4A4A4A', lw=2, style='->', connectionstyle='arc3,rad=0'):
    """Draw an arrow from (x1,y1) to (x2,y2)."""
    arrow = FancyArrowPatch((x1, y1), (x2, y2),
                            arrowstyle=style, color=color, linewidth=lw,
                            connectionstyle=connectionstyle,
                            mutation_scale=15, zorder=1)
    ax.add_patch(arrow)

def draw_arrow_with_points(ax, points, color='#4A4A4A', lw=2):
    """Draw an arrow through multiple points (orthogonal routing)."""
    arrow = FancyArrowPatch(points[0], points[-1],
                            arrowstyle='->', color=color, linewidth=lw,
                            connectionstyle=f"bar,fraction=-0.1,angle=0",
                            mutation_scale=15, zorder=1)
    # For orthogonal routing, draw line segments manually
    for i in range(len(points) - 1):
        x1, y1 = points[i]
        x2, y2 = points[i + 1]
        if i == len(points) - 2:
            # Last segment: draw arrow
            arrow = FancyArrowPatch((x1, y1), (x2, y2),
                                    arrowstyle='->', color=color, linewidth=lw,
                                    mutation_scale=15, zorder=1)
            ax.add_patch(arrow)
        else:
            ax.plot([x1, x2], [y1, y2], color=color, linewidth=lw, zorder=1, solid_capstyle='round')

# ===== Create figure =====
fig, ax = plt.subplots(figsize=(14, 12))
ax.set_xlim(0, 14)
ax.set_ylim(0, 12)
ax.axis('off')

# Title
ax.text(7, 11.7, 'Figure 1. Study Population Flow Diagram',
        ha='center', va='center', fontsize=16, fontweight='bold', color=C_DARK)

# --- Box 1: MIMIC-IV Total ---
b1 = draw_box(ax, 5, 10.4, 4, 0.7,
              'MIMIC-IV v3.1 ICU Admissions (2008–2022)   n = 63,798',
              fc=C_BLUE, ec=C_BLUE_EDGE, fs=12, fw='bold')

# ===== Layout coordinates =====
# Canvas: 14 x 12 units
# Center column: x=5 to x=9 (width=4)
# Left column: x=0.5 to x=3 (width=2.5)
# Right column: x=11 to x=13.5

# Arrow 1->2
draw_arrow(ax, 7, 10.4, 7, 9.8, color=C_BLUE_EDGE)

# --- Box 2: Inclusion criteria ---
b2 = draw_box(ax, 5, 9.5, 4, 0.6,
              'Inclusion Criteria\nAge ≥ 18; First ICU admission; ICU LOS ≥ 24h',
              fc=C_PURPLE, ec=C_PURPLE_EDGE, fs=10)

# Arrow 2->3
draw_arrow(ax, 7, 9.5, 7, 9.0, color=C_PURPLE_EDGE)

# --- Box 3: KDIGO Assessment ---
b3 = draw_box(ax, 5, 8.4, 4, 0.5,
              'KDIGO AKI Assessment  (Creatinine + Urine Output criteria)',
              fc=C_YELLOW, ec=C_YELLOW_EDGE, fs=11, fw='bold')

# --- Split: No AKI (left) vs AKI (center) ---

# No AKI box (left)
b_noaki = draw_box(ax, 0.5, 6.6, 3, 0.7,
                   'No AKI (KDIGO Stage 0)\nn = 39,594  (excluded)',
                   fc=C_RED, ec=C_RED_EDGE, fs=10, lw=1.5)

# Arrow from KDIGO to No AKI (orthogonal: left then down)
draw_arrow_with_points(ax, [(5, 8.65), (3.5, 8.65), (3.5, 7.5), (2, 7.5)],
                       color=C_RED_EDGE, lw=2)

# AKI cohort box (center)
b_aki = draw_box(ax, 5, 6.6, 4, 0.7,
                 'AKI Cohort (KDIGO Stage 1–3)   n = 24,204',
                 fc=C_GREEN, ec=C_GREEN_EDGE, fs=14, fw='bold', lw=2)

# Arrow from KDIGO to AKI
draw_arrow(ax, 7, 8.4, 7, 7.3, color=C_GREEN_EDGE, lw=2.5)

# ===== Three AKI Stage Branches =====
# Stage 1 (left-center), Stage 2 (center), Stage 3 (right)

b_s1 = draw_box(ax, 1.0, 5.0, 2.5, 0.7,
                'AKI Stage 1   n = 13,066',
                fc=C_GREEN, ec=C_GREEN_EDGE, fs=12, fw='bold')

b_s2 = draw_box(ax, 5.75, 5.0, 2.5, 0.7,
                'AKI Stage 2   n = 8,120',
                fc=C_GREEN, ec=C_GREEN_EDGE, fs=12, fw='bold')

b_s3 = draw_box(ax, 10.5, 5.0, 2.5, 0.7,
                'AKI Stage 3   n = 3,018',
                fc=C_GREEN, ec=C_GREEN_EDGE, fs=12, fw='bold', lw=2.5)

# Arrows from AKI to 3 stages
# To Stage 1: down-left
draw_arrow_with_points(ax, [(6, 6.6), (6, 6.0), (2.25, 6.0), (2.25, 5.7)],
                       color=C_GREEN_EDGE, lw=2)
# To Stage 2: straight down
draw_arrow(ax, 7, 6.6, 7, 5.7, color=C_GREEN_EDGE, lw=2)
# To Stage 3: down-right
draw_arrow_with_points(ax, [(8, 6.6), (8, 6.0), (11.75, 6.0), (11.75, 5.7)],
                       color=C_GREEN_EDGE, lw=2)

# ===== Stage 1: Pure AKI vs AKI+CKD =====
b_s1p = draw_box(ax, 0.3, 3.5, 2.2, 0.6,
                 'Pure AKI  n = 10,070\nMortality: 13.1%',
                 fc=C_GREY, ec=C_GREY_EDGE, fs=9.5)

b_s1c = draw_box(ax, 2.8, 3.5, 2.2, 0.6,
                 'AKI + CKD  n = 2,996\nMortality: 16.7%',
                 fc=C_GREY, ec=C_GREY_EDGE, fs=9.5)

# Arrows Stage 1 -> Pure / CKD
draw_arrow_with_points(ax, [(1.8, 5.0), (1.8, 4.5), (1.4, 4.5), (1.4, 4.1)],
                       color=C_GREY_EDGE, lw=1.5)
draw_arrow_with_points(ax, [(2.2, 5.0), (2.2, 4.5), (3.9, 4.5), (3.9, 4.1)],
                       color=C_GREY_EDGE, lw=1.5)

# ===== Stage 2: Pure AKI vs AKI+CKD =====
b_s2p = draw_box(ax, 5.0, 3.5, 2.2, 0.6,
                 'Pure AKI  n = 6,764\nMortality: 16.5%',
                 fc=C_GREY, ec=C_GREY_EDGE, fs=9.5)

b_s2c = draw_box(ax, 7.5, 3.5, 2.2, 0.6,
                 'AKI + CKD  n = 1,356\nMortality: 19.5%',
                 fc=C_GREY, ec=C_GREY_EDGE, fs=9.5)

# Arrows Stage 2 -> Pure / CKD
draw_arrow_with_points(ax, [(6.5, 5.0), (6.5, 4.5), (6.1, 4.5), (6.1, 4.1)],
                       color=C_GREY_EDGE, lw=1.5)
draw_arrow_with_points(ax, [(6.9, 5.0), (6.9, 4.5), (8.6, 4.5), (8.6, 4.1)],
                       color=C_GREY_EDGE, lw=1.5)

# ===== Stage 3: Pure AKI vs AKI+CKD (KEY FINDING) =====
b_s3p = draw_box(ax, 9.8, 3.5, 2.2, 0.6,
                 'Pure AKI  n = 1,998\nMortality: 42.0%',
                 fc='#ffe6e6', ec=C_RED_EDGE, fs=9.5, lw=2)

b_s3c = draw_box(ax, 12.2, 3.5, 1.6, 0.6,
                 'AKI + CKD  n = 1,020\nMortality: 35.4%',
                 fc='#e6ffe6', ec=C_GREEN_EDGE, fs=9.5, lw=2)

# Arrows Stage 3 -> Pure / CKD
draw_arrow_with_points(ax, [(11.3, 5.0), (11.3, 4.5), (10.9, 4.5), (10.9, 4.1)],
                       color=C_RED_EDGE, lw=1.5)
draw_arrow_with_points(ax, [(11.7, 5.0), (11.7, 4.5), (13.0, 4.5), (13.0, 4.1)],
                       color=C_GREEN_EDGE, lw=1.5)

# ===== Key Finding Callout (moved to far right, separate column) =====
# Yellow callout box at right (separate from CKD note which is center)
finding_box = FancyBboxPatch((9.5, 0.5), 4.3, 1.4,
                              boxstyle="round,pad=0.15,rounding_size=0.2",
                              facecolor=C_YELLOW, edgecolor=C_YELLOW_EDGE,
                              linewidth=2, zorder=2)
ax.add_patch(finding_box)
ax.text(11.65, 1.2,
        'Key Finding: The AKI Paradox (Stage 3)\n\n'
        'AKI+CKD mortality LOWER than Pure AKI\n'
        '(35.4% vs 42.0%, OR = 0.643, P < 0.0001)\n'
        'Paradox STRENGTHENS after SOFA adjustment',
        ha='center', va='center', fontsize=9.5, fontweight='normal',
        color=C_DARK, zorder=3, linespacing=1.5)

# Dashed arrow from Stage 3 CKD to finding
ax.plot([13.0, 13.0], [3.5, 2.4], color=C_YELLOW_EDGE, linewidth=2,
        linestyle='--', zorder=1)
ax.plot([13.0, 12.0], [2.4, 2.0], color=C_YELLOW_EDGE, linewidth=2,
        linestyle='--', zorder=1)
arrow_to_finding = FancyArrowPatch((12.0, 2.0), (11.7, 1.9),
                                    arrowstyle='->', color=C_YELLOW_EDGE,
                                    linewidth=2, mutation_scale=12, zorder=1)
ax.add_patch(arrow_to_finding)

# ===== Analysis box (bottom left) =====
analysis_box = FancyBboxPatch((0.5, 0.5), 5.0, 1.4,
                               boxstyle="round,pad=0.15,rounding_size=0.2",
                               facecolor=C_PURPLE, edgecolor=C_PURPLE_EDGE,
                               linewidth=1.5, zorder=2)
ax.add_patch(analysis_box)
ax.text(3.0, 1.2,
        'Statistical Analysis\n\n'
        '• Multivariable logistic regression (Models A–D)\n'
        '• PSM (Propensity Score Matching) ± SOFA\n'
        '• CKD × AKI stage interaction testing\n'
        '• SOFA-stratified sensitivity analysis\n'
        '• Prediction modeling + Nomogram + LASSO',
        ha='center', va='center', fontsize=9,
        color=C_DARK, zorder=3, linespacing=1.4)

# ===== CKD note (bottom center) =====
ckd_box = FancyBboxPatch((5.8, 0.5), 3.5, 0.7,
                          boxstyle="round,pad=0.1,rounding_size=0.15",
                          facecolor=C_GREY, edgecolor='#999999',
                          linewidth=1, zorder=2)
ax.add_patch(ckd_box)
ax.text(7.55, 0.85,
        'CKD Identification\n'
        'ICD-9: 585.x; ICD-10: N18.x\n'
        'n = 5,372 (22.2% of AKI)',
        ha='center', va='center', fontsize=8.5,
        color='#666666', zorder=3, style='italic', linespacing=1.3)

# ===== Save =====
import os
os.makedirs('submission_figures', exist_ok=True)
plt.tight_layout(pad=0.5)
fig.savefig('submission_figures/Figure_1_Flow_Diagram.png',
            dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.3)
fig.savefig('submission_figures/Figure_1_Flow_Diagram.pdf',
            bbox_inches='tight', facecolor='white', pad_inches=0.3)
plt.close(fig)
print('Figure 1 saved (matplotlib flow diagram, 300 DPI)')
