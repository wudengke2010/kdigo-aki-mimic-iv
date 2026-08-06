"""
Convert manuscript Markdown files to DOCX for Renal Failure submission.
Handles: headings, bold/italic, tables, references, figure legends.
"""
import re
import sys
from pathlib import Path
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT

def add_formatted_runs(paragraph, text):
    """Parse markdown bold/italic and add runs to paragraph."""
    # Handle bold+italic, bold, italic patterns
    patterns = [
        (r'\*\*\*(.+?)\*\*\*', 'bold_italic'),
        (r'\*\*(.+?)\*\*', 'bold'),
        (r'\*(.+?)\*', 'italic'),
    ]
    
    remaining = text
    runs_to_add = []
    
    while remaining:
        earliest_pos = len(remaining)
        earliest_match = None
        earliest_type = None
        
        for pattern, fmt_type in patterns:
            match = re.search(pattern, remaining)
            if match and match.start() < earliest_pos:
                earliest_pos = match.start()
                earliest_match = match
                earliest_type = fmt_type
        
        if earliest_match:
            # Add text before the match
            if earliest_pos > 0:
                runs_to_add.append(('normal', remaining[:earliest_pos]))
            
            # Add the matched text
            inner_text = earliest_match.group(1)
            runs_to_add.append((earliest_type, inner_text))
            
            remaining = remaining[earliest_match.end():]
        else:
            if remaining:
                runs_to_add.append(('normal', remaining))
            remaining = ''
    
    for fmt, text in runs_to_add:
        run = paragraph.add_run(text)
        if fmt == 'bold' or fmt == 'bold_italic':
            run.bold = True
        if fmt == 'italic' or fmt == 'bold_italic':
            run.italic = True

def parse_table(lines, start_idx):
    """Parse a markdown pipe table starting at start_idx. Returns (rows, end_idx)."""
    rows = []
    i = start_idx
    while i < len(lines) and '|' in lines[i] and lines[i].strip():
        line = lines[i].strip()
        if line.startswith('|'):
            line = line[1:]
        if line.endswith('|'):
            line = line[:-1]
        cells = [c.strip() for c in line.split('|')]
        rows.append(cells)
        i += 1
    # Skip separator row (---|---|---)
    if len(rows) >= 2 and all(re.match(r'^[-:]+$', c) for c in rows[1] if c):
        rows = [rows[0]] + rows[2:]
    return rows, i

def add_table_to_doc(doc, rows):
    """Add a formatted table to the document."""
    if not rows:
        return
    num_cols = len(rows[0])
    table = doc.add_table(rows=len(rows), cols=num_cols)
    table.style = 'Table Grid'
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    
    for i, row in enumerate(rows):
        for j, cell_text in enumerate(row):
            if j < num_cols:
                cell = table.cell(i, j)
                cell.text = ''
                para = cell.paragraphs[0]
                add_formatted_runs(para, cell_text)
                para.style.font.size = Pt(9)
                if i == 0:
                    for run in para.runs:
                        run.bold = True

def convert_md_to_docx(md_path, docx_path, exclude_title_page=False):
    """Convert a Markdown file to DOCX."""
    with open(md_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    lines = content.split('\n')
    doc = Document()
    
    # Set default font
    style = doc.styles['Normal']
    font = style.font
    font.name = 'Times New Roman'
    font.size = Pt(12)
    
    # Set margins
    for section in doc.sections:
        section.top_margin = Inches(1)
        section.bottom_margin = Inches(1)
        section.left_margin = Inches(1)
        section.right_margin = Inches(1)
    
    i = 0
    in_title_page = False
    skip_until_next_h2 = False
    
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        
        # Skip horizontal rules
        if stripped == '---':
            i += 1
            continue
        
        # Skip empty lines
        if not stripped:
            i += 1
            continue
        
        # Skip metadata/version note (italic lines at end)
        if stripped.startswith('*Revised manuscript'):
            i += 1
            continue
        
        # Handle headings
        if stripped.startswith('### '):
            heading_text = stripped[4:]
            heading = doc.add_heading(level=3)
            add_formatted_runs(heading, heading_text)
            i += 1
            continue
        
        if stripped.startswith('## '):
            heading_text = stripped[3:]
            
            # If excluding title page, skip everything from "## Authors" to next "## ---" or "## Abstract"
            if exclude_title_page and heading_text in ['Authors']:
                skip_until_next_h2 = True
                i += 1
                continue
            
            if skip_until_next_h2:
                skip_until_next_h2 = False
            
            heading = doc.add_heading(level=2)
            add_formatted_runs(heading, heading_text)
            i += 1
            continue
        
        if stripped.startswith('# '):
            heading_text = stripped[2:]
            heading = doc.add_heading(level=1)
            add_formatted_runs(heading, heading_text)
            i += 1
            continue
        
        # Skip if in title page section
        if skip_until_next_h2:
            i += 1
            continue
        
        # Handle tables
        if '|' in stripped and stripped.startswith('|'):
            rows, end_idx = parse_table(lines, i)
            if rows:
                add_table_to_doc(doc, rows)
            i = end_idx
            continue
        
        # Handle numbered references (e.g., "1. Author...")
        ref_match = re.match(r'^(\d+)\.\s+(.+)', stripped)
        if ref_match:
            ref_text = ref_match.group(2)
            para = doc.add_paragraph()
            para.paragraph_format.left_indent = Inches(0.3)
            para.paragraph_format.first_line_indent = Inches(-0.3)
            run = para.add_run(f"{ref_match.group(1)}. ")
            add_formatted_runs(para, ref_text)
            para.style.font.size = Pt(10)
            i += 1
            continue
        
        # Handle bullet points
        if stripped.startswith('- '):
            bullet_text = stripped[2:]
            para = doc.add_paragraph(style='List Bullet')
            add_formatted_runs(para, bullet_text)
            i += 1
            continue
        
        # Regular paragraph
        para = doc.add_paragraph()
        add_formatted_runs(para, stripped)
        i += 1
    
    doc.save(docx_path)
    print(f"Converted: {md_path} -> {docx_path}")

if __name__ == '__main__':
    base = Path(r'C:\Users\admin\WorkBuddy\2026-07-06-13-22-19')
    
    # 1. Manuscript (exclude title page section)
    convert_md_to_docx(
        base / 'manuscript_revised_v2.md',
        base / 'manuscript_renal_failure_v9.docx',
        exclude_title_page=True
    )
    
    # 2. Title page
    convert_md_to_docx(
        base / 'title_page_renal_failure.md',
        base / 'title_page_renal_failure.docx'
    )
    
    # 3. Cover letter
    convert_md_to_docx(
        base / 'cover_letter_renal_failure.md',
        base / 'cover_letter_renal_failure.docx'
    )
    
    # 4. Supplementary materials index
    convert_md_to_docx(
        base / 'supplementary_materials_index.md',
        base / 'supplementary_materials_index.docx'
    )
    
    print("\nAll DOCX files generated successfully.")
