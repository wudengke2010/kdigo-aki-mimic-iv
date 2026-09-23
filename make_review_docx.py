#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build internal-review DOCX: manuscript with figures embedded after legends
and supplementary Tables S1-S4 appended. Output: jtm_submission/manuscript_jtm_review.docx"""

import re
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING

from md_to_docx_ckj import (
    add_page_number, set_run_format, add_formatted_text,
    parse_table, add_table_to_doc,
)

BASE = r'C:\Users\admin\WorkBuddy\2026-07-06-13-22-19'
FIGDIR = BASE + r'\jtm_submission\figures'

FIG_IMAGES = {
    'Figure 1.': ('Figure1_flowchart.png', 5.8),
    'Figure 2.': ('Figure2_KM.png', 6.5),
    'Figure 3.': ('Figure3_forest_time_stratified.png', 6.5),
    'Figure 4.': ('Figure4_external_validation.png', 6.5),
    'Figure S1': ('FigureS1_DAG.png', 4.6),
    'Figure S2': ('FigureS2_subgroup_forest.png', 6.5),
    'Figure S3': ('FigureS3_phenotype_forest.png', 6.5),
}


def insert_image(doc, filename, width_in):
    import os
    path = os.path.join(FIGDIR, filename)
    if not os.path.exists(path):
        p = doc.add_paragraph()
        r = p.add_run(f'[MISSING IMAGE: {filename}]')
        r.font.color.rgb = None
        return
    para = doc.add_paragraph()
    para.alignment = WD_ALIGN_PARAGRAPH.CENTER
    para.paragraph_format.space_before = Pt(6)
    para.paragraph_format.space_after = Pt(12)
    run = para.add_run()
    run.add_picture(path, width=Inches(width_in))


def convert_lines(doc, lines, embed_figs=False):
    i = 0
    in_references = False
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if not stripped or stripped == '---':
            i += 1
            continue
        # Title
        if stripped.startswith('# ') and not stripped.startswith('## '):
            para = doc.add_paragraph()
            para.alignment = WD_ALIGN_PARAGRAPH.CENTER
            para.paragraph_format.space_after = Pt(24)
            para.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            run = para.add_run(stripped[2:])
            set_run_format(run, font_size=14, bold=True)
            i += 1
            continue
        if stripped.startswith('### '):
            para = doc.add_heading(level=2)
            para.paragraph_format.space_before = Pt(12)
            para.paragraph_format.space_after = Pt(6)
            para.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            add_formatted_text(para, stripped[4:], font_size=12)
            for run in para.runs:
                run.font.bold = True
                run.font.name = 'Times New Roman'
                run.font.size = Pt(12)
            i += 1
            continue
        if stripped.startswith('## '):
            heading_text = stripped[3:]
            in_references = 'References' in heading_text
            para = doc.add_heading(level=1)
            para.paragraph_format.space_before = Pt(18)
            para.paragraph_format.space_after = Pt(12)
            para.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            add_formatted_text(para, heading_text, font_size=14)
            for run in para.runs:
                run.font.bold = True
                run.font.name = 'Times New Roman'
                run.font.size = Pt(14)
            i += 1
            continue
        if stripped.startswith('|'):
            rows, next_i = parse_table(lines, i)
            if rows:
                fs = 9 if (len(rows) > 2 and len(rows[0]) > 3) else 10
                add_table_to_doc(doc, rows, font_size=fs)
            i = next_i
            continue
        if stripped.startswith('- '):
            para = doc.add_paragraph(style='List Bullet')
            para.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            para.paragraph_format.space_after = Pt(3)
            add_formatted_text(para, stripped[2:], font_size=12)
            i += 1
            continue
        if in_references and re.match(r'^\d+\.\s', stripped):
            ref_text = re.sub(r'^\d+\.\s', '', stripped)
            para = doc.add_paragraph()
            para.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            para.paragraph_format.space_after = Pt(6)
            para.paragraph_format.left_indent = Inches(0.3)
            para.paragraph_format.first_line_indent = Inches(-0.3)
            add_formatted_text(para, ref_text, font_size=10)
            i += 1
            continue
        # Regular paragraph
        para = doc.add_paragraph()
        para.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
        para.paragraph_format.first_line_indent = Inches(0.5)
        add_formatted_text(para, stripped, font_size=12)
        # Embed figure image after its legend
        if embed_figs:
            for key, (fn, w) in FIG_IMAGES.items():
                if stripped.startswith('**' + key):
                    insert_image(doc, fn, w)
                    break
        i += 1


def main():
    doc = Document()
    style = doc.styles['Normal']
    style.font.name = 'Times New Roman'
    style.font.size = Pt(12)
    style.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    style.paragraph_format.space_after = Pt(0)
    for section in doc.sections:
        section.top_margin = Inches(1)
        section.bottom_margin = Inches(1)
        section.left_margin = Inches(1)
        section.right_margin = Inches(1)
    add_page_number(doc)

    # 1. Main manuscript with embedded figures
    with open(BASE + r'\jtm_submission\manuscript_jtm.md', encoding='utf-8') as f:
        ms_lines = f.read().split('\n')
    convert_lines(doc, ms_lines, embed_figs=True)

    # 2. Supplementary tables for review
    doc.add_page_break()
    h = doc.add_heading('Supplementary Tables (embedded for internal review)', level=1)
    for run in h.runs:
        run.font.name = 'Times New Roman'
        run.font.size = Pt(14)
    note = doc.add_paragraph()
    nr = note.add_run('Note: In the final submission, Tables S1-S4 are uploaded as separate Additional files; '
                      'they are embedded here for convenience of internal review.')
    nr.italic = True
    nr.font.size = Pt(10)

    # Table S1 from table_S1.md (skip its own title line, keep cohort sub-headers)
    with open(BASE + r'\v2_outputs\table_S1.md', encoding='utf-8') as f:
        s1_lines = [l for l in f.read().split('\n') if not l.startswith('# ')]
    convert_lines(doc, s1_lines, embed_figs=False)

    # Sensitivity report S1-S4 sections -> Tables S2-S5 style, keep section text+tables
    with open(BASE + r'\v2_outputs\sensitivity_report.md', encoding='utf-8') as f:
        rep = f.read().split('\n')
    # drop the top title, keep from first "## S1." onward
    start = next(idx for idx, l in enumerate(rep) if l.startswith('## S1.'))
    convert_lines(doc, rep[start:], embed_figs=False)

    out = BASE + r'\jtm_submission\manuscript_jtm_review.docx'
    doc.save(out)
    import os
    print(f'Saved {out} ({os.path.getsize(out)/1024:.0f} KB)')


if __name__ == '__main__':
    main()
