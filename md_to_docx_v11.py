#!/usr/bin/env python3
"""
Convert manuscript_revised_v3.md (v11 dual-cohort) to submission-formatted DOCX.
Output: manuscript_dual_cohort_v11.docx

Formatting standards:
- Font: Times New Roman, 12pt body / 14pt headings
- Line spacing: Double (1.5 for tables)
- Margins: 1 inch (2.54 cm)
- Title: Bold, centered, 14pt
- Tables: Grid style with borders
- Page numbers: Bottom center
"""

import re
import os
from docx import Document
from docx.shared import Pt, Inches, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_ALIGN_VERTICAL
from docx.enum.section import WD_ORIENT
from docx.oxml.ns import qn, nsdecls
from docx.oxml import parse_xml

def add_page_number(doc):
    """Add page number to footer (centered)."""
    section = doc.sections[0]
    footer = section.footer
    footer_para = footer.paragraphs[0] if footer.paragraphs else footer.add_paragraph()
    footer_para.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = footer_para.add_run()
    fldChar1 = parse_xml(f'<w:fldChar {nsdecls("w")} w:fldCharType="begin"/>')
    run._element.append(fldChar1)
    run2 = footer_para.add_run()
    instrText = parse_xml(f'<w:instrText {nsdecls("w")} xml:space="preserve">PAGE</w:instrText>')
    run2._element.append(instrText)
    run3 = footer_para.add_run()
    fldChar2 = parse_xml(f'<w:fldChar {nsdecls("w")} w:fldCharType="end"/>')
    run3._element.append(fldChar2)

def set_cell_border(cell, **kwargs):
    """Set cell border."""
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    tcBorders = parse_xml(f'<w:tcBorders {nsdecls("w")}>'
                          '<w:top w:val="single" w:sz="4" w:space="0" w:color="000000"/>'
                          '<w:left w:val="single" w:sz="4" w:space="0" w:color="000000"/>'
                          '<w:bottom w:val="single" w:sz="4" w:space="0" w:color="000000"/>'
                          '<w:right w:val="single" w:sz="4" w:space="0" w:color="000000"/>'
                          '</w:tcBorders>')
    tcPr.append(tcBorders)

def set_run_format(run, font_name='Times New Roman', font_size=12, bold=False, italic=False, superscript=False):
    """Set run formatting."""
    run.font.name = font_name
    run.font.size = Pt(font_size)
    run.font.bold = bold
    run.font.italic = italic
    if superscript:
        run.font.superscript = True
    # Set East Asian font
    rPr = run._element.get_or_add_rPr()
    rFonts = rPr.find(qn('w:rFonts'))
    if rFonts is None:
        rFonts = parse_xml(f'<w:rFonts {nsdecls("w")}/>')
        rPr.append(rFonts)
    rFonts.set(qn('w:eastAsia'), font_name)

def add_formatted_text(paragraph, text, font_size=12, bold_marker='**', italic_marker='*'):
    """Add text with inline formatting (bold, italic, superscript)."""
    # Process superscripts: <sup>text</sup>
    parts = re.split(r'(<sup>.*?</sup>)', text)
    for part in parts:
        if part.startswith('<sup>'):
            content = part[5:-6]  # Strip <sup></sup>
            run = paragraph.add_run(content)
            set_run_format(run, font_size=font_size, superscript=True)
        else:
            # Process bold and italic within this part
            add_bold_italic_text(paragraph, part, font_size)

def add_bold_italic_text(paragraph, text, font_size=12):
    """Process **bold** and *italic* markers."""
    # Split on bold first
    bold_parts = re.split(r'(\*\*.*?\*\*)', text)
    for bp in bold_parts:
        if bp.startswith('**') and bp.endswith('**') and len(bp) > 4:
            content = bp[2:-2]
            # Check for italic within bold
            italic_parts = re.split(r'(\*.*?\*)', content)
            for ip in italic_parts:
                if ip.startswith('*') and ip.endswith('*') and len(ip) > 2:
                    run = paragraph.add_run(ip[1:-1])
                    set_run_format(run, font_size=font_size, bold=True, italic=True)
                elif ip:
                    run = paragraph.add_run(ip)
                    set_run_format(run, font_size=font_size, bold=True)
        else:
            # Check for italic in non-bold
            italic_parts = re.split(r'(\*.*?\*)', bp)
            for ip in italic_parts:
                if ip.startswith('*') and ip.endswith('*') and len(ip) > 2:
                    content = ip[1:-1]
                    # Check if it's actually an italic marker vs a multiplication sign
                    # Only treat as italic if preceded by non-word char or start
                    run = paragraph.add_run(content)
                    set_run_format(run, font_size=font_size, italic=True)
                elif ip:
                    run = paragraph.add_run(ip)
                    set_run_format(run, font_size=font_size)

def parse_table(lines, start_idx):
    """Parse a markdown table starting at start_idx. Returns (rows, next_idx)."""
    rows = []
    i = start_idx
    while i < len(lines) and '|' in lines[i].strip():
        line = lines[i].strip()
        if line.startswith('|'):
            line = line[1:]
        if line.endswith('|'):
            line = line[:-1]
        cells = [c.strip() for c in line.split('|')]
        rows.append(cells)
        i += 1
    # Remove separator row (---)
    if len(rows) > 1 and all(re.match(r'^[-:]+$', c.replace(' ', '')) for c in rows[1]):
        rows.pop(1)
    return rows, i

def add_table_to_doc(doc, rows, font_size=10):
    """Add a formatted table to the document."""
    if not rows:
        return
    n_cols = len(rows[0])
    n_rows = len(rows)
    table = doc.add_table(rows=n_rows, cols=n_cols)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.style = 'Table Grid'

    for i, row in enumerate(rows):
        for j, cell_text in enumerate(row):
            if j >= n_cols:
                break
            cell = table.rows[i].cells[j]
            cell.text = ''  # Clear default
            cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
            para = cell.paragraphs[0]
            para.alignment = WD_ALIGN_PARAGRAPH.CENTER if i == 0 else WD_ALIGN_PARAGRAPH.LEFT
            para.paragraph_format.space_before = Pt(2)
            para.paragraph_format.space_after = Pt(2)
            para.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE

            is_header = (i == 0)
            add_formatted_text(para, cell_text, font_size=font_size)
            # Make header bold
            if is_header:
                for run in para.runs:
                    run.font.bold = True

    # Add spacing after table
    doc.add_paragraph()

def process_inline_elements(text):
    """Convert markdown inline elements to simpler forms for docx."""
    # Convert <sup> tags - will be handled in add_formatted_text
    return text

def convert_md_to_docx(md_path, docx_path):
    """Convert markdown to formatted DOCX."""
    with open(md_path, 'r', encoding='utf-8') as f:
        content = f.read()

    lines = content.split('\n')

    doc = Document()

    # Set default style
    style = doc.styles['Normal']
    style.font.name = 'Times New Roman'
    style.font.size = Pt(12)
    style.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    style.paragraph_format.space_after = Pt(0)

    # Set margins
    for section in doc.sections:
        section.top_margin = Inches(1)
        section.bottom_margin = Inches(1)
        section.left_margin = Inches(1)
        section.right_margin = Inches(1)

    add_page_number(doc)

    i = 0
    in_references = False

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # Skip empty lines
        if not stripped:
            i += 1
            continue

        # Skip horizontal rules
        if stripped == '---':
            i += 1
            continue

        # Skip the version note at the end
        if stripped.startswith('*Revised manuscript (v1'):
            i += 1
            continue

        # Title (first # heading)
        if stripped.startswith('# ') and not stripped.startswith('## '):
            title_text = stripped[2:]
            para = doc.add_paragraph()
            para.alignment = WD_ALIGN_PARAGRAPH.CENTER
            para.paragraph_format.space_before = Pt(0)
            para.paragraph_format.space_after = Pt(24)
            para.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            run = para.add_run(title_text)
            set_run_format(run, font_size=14, bold=True)
            i += 1
            continue

        # Headings
        if stripped.startswith('### '):
            heading_text = stripped[4:]
            para = doc.add_heading(level=2)
            para.paragraph_format.space_before = Pt(12)
            para.paragraph_format.space_after = Pt(6)
            para.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            add_formatted_text(para, heading_text, font_size=12, bold_marker='**')
            for run in para.runs:
                run.font.bold = True
                run.font.italic = False
                run.font.name = 'Times New Roman'
                run.font.size = Pt(12)
            i += 1
            continue

        if stripped.startswith('## '):
            heading_text = stripped[3:]
            # Check for References section
            if 'References' in heading_text:
                in_references = True
            else:
                in_references = False

            para = doc.add_heading(level=1)
            para.paragraph_format.space_before = Pt(18)
            para.paragraph_format.space_after = Pt(12)
            para.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            add_formatted_text(para, heading_text, font_size=14, bold_marker='**')
            for run in para.runs:
                run.font.bold = True
                run.font.italic = False
                run.font.name = 'Times New Roman'
                run.font.size = Pt(14)
            i += 1
            continue

        if stripped.startswith('#### '):
            heading_text = stripped[5:]
            para = doc.add_heading(level=3)
            para.paragraph_format.space_before = Pt(6)
            para.paragraph_format.space_after = Pt(3)
            para.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            add_formatted_text(para, heading_text, font_size=12, bold_marker='**')
            for run in para.runs:
                run.font.bold = True
                run.font.italic = True
                run.font.name = 'Times New Roman'
                run.font.size = Pt(12)
            i += 1
            continue

        # Tables
        if stripped.startswith('|'):
            rows, next_i = parse_table(lines, i)
            if rows:
                # Check if it's a large table (Tables 1-3) vs inline
                if len(rows) > 2 and len(rows[0]) > 3:
                    # Landscape for wide tables
                    add_table_to_doc(doc, rows, font_size=9)
                else:
                    add_table_to_doc(doc, rows, font_size=10)
            i = next_idx = next_i
            continue

        # Bullet lists
        if stripped.startswith('- '):
            bullet_text = stripped[2:]
            # Check if it's a sub-bullet
            indent_level = 0
            if line.startswith('  - '):
                indent_level = 1
                bullet_text = stripped[2:]

            para = doc.add_paragraph(style='List Bullet')
            para.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            para.paragraph_format.space_before = Pt(0)
            para.paragraph_format.space_after = Pt(3)
            if indent_level > 0:
                para.paragraph_format.left_indent = Inches(0.5 + 0.25 * indent_level)
            add_formatted_text(para, bullet_text, font_size=12)
            i += 1
            continue

        # References (numbered list with special formatting)
        if in_references and re.match(r'^\d+\.\s', stripped):
            ref_text = re.sub(r'^\d+\.\s', '', stripped)
            para = doc.add_paragraph()
            para.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            para.paragraph_format.space_before = Pt(0)
            para.paragraph_format.space_after = Pt(6)
            para.paragraph_format.left_indent = Inches(0.3)
            para.paragraph_format.first_line_indent = Inches(-0.3)  # Hanging indent
            add_formatted_text(para, ref_text, font_size=10)
            i += 1
            continue

        # Regular paragraph
        para = doc.add_paragraph()
        para.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
        para.paragraph_format.space_before = Pt(0)
        para.paragraph_format.space_after = Pt(0)
        para.paragraph_format.first_line_indent = Inches(0.5)
        add_formatted_text(para, stripped, font_size=12)
        i += 1

    # Save
    doc.save(docx_path)
    print(f"DOCX saved: {docx_path}")
    print(f"File size: {os.path.getsize(docx_path) / 1024:.1f} KB")

if __name__ == '__main__':
    md_path = r'C:\Users\admin\WorkBuddy\2026-07-06-13-22-19\manuscript_revised_v3.md'
    docx_path = r'C:\Users\admin\WorkBuddy\2026-07-06-13-22-19\manuscript_dual_cohort_v12.docx'
    convert_md_to_docx(md_path, docx_path)
