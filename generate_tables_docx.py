"""Extract Tables 1-2 from manuscript and generate standalone DOCX files."""
import re
from pathlib import Path
from docx import Document
from docx.shared import Inches, Pt, Cm, RGBColor
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

BASE = Path(r"C:\Users\admin\WorkBuddy\2026-07-06-13-22-19")
MANUSCRIPT = BASE / "manuscript_revised_v2.md"

def set_cell_shading(cell, color):
    """Set cell background color."""
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    shading = OxmlElement("w:shd")
    shading.set(qn("w:fill"), color)
    shading.set(qn("w:val"), "clear")
    tcPr.append(shading)

def set_cell_border(cell, **kwargs):
    """Set cell borders."""
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    tcBorders = OxmlElement("w:tcBorders")
    for edge, val in kwargs.items():
        element = OxmlElement(f"w:{edge}")
        element.set(qn("w:val"), val.get("val", "single"))
        element.set(qn("w:sz"), val.get("sz", "4"))
        element.set(qn("w:color"), val.get("color", "000000"))
        element.set(qn("w:space"), "0")
        tcBorders.append(element)
    tcPr.append(tcBorders)

def parse_markdown_table(lines):
    """Parse a markdown table into headers and rows."""
    # Find separator row
    header_idx = None
    sep_idx = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("|") and stripped.endswith("|"):
            if re.match(r'^\|[\s\-:|]+\|$', stripped):
                if header_idx is not None and sep_idx is None:
                    sep_idx = i
            elif header_idx is None:
                header_idx = i
    
    if header_idx is None or sep_idx is None:
        return None, None, []
    
    # Parse headers
    header_line = lines[header_idx].strip()
    headers = [h.strip() for h in header_line.split("|")[1:-1]]
    
    # Parse data rows
    rows = []
    for i in range(sep_idx + 1, len(lines)):
        line = lines[i].strip()
        if not line.startswith("|"):
            break
        cells = [c.strip() for c in line.split("|")[1:-1]]
        rows.append(cells)
    
    return headers, rows, lines[sep_idx + len(rows) + 1:]

def extract_tables(md_path):
    """Extract Table 1 and Table 2 from the markdown manuscript."""
    with open(md_path, "r", encoding="utf-8") as f:
        content = f.read()
    
    lines = content.split("\n")
    
    tables = {}
    
    # Find table sections
    in_table = False
    table_lines = []
    table_title = ""
    table_footnotes = []  # Will collect footnotes after table
    
    for i, line in enumerate(lines):
        if line.startswith("### Table 1"):
            in_table = True
            table_title = "Table 1. Baseline Characteristics of the Study Population Stratified by KDIGO AKI Stage (N = 84,167)"
            table_lines = []
            continue
        elif line.startswith("### Table 2"):
            if in_table and table_lines:
                tables["Table 1"] = (table_title, table_lines)
            in_table = True
            table_title = "Table 2. Multivariable Logistic Regression: Comparison of Three Incremental Models"
            table_lines = []
            continue
        
        if in_table:
            if line.strip().startswith("### ") or line.strip().startswith("## "):
                if table_lines:
                    tables[table_title.split(".")[0]] = (table_title, table_lines)
                in_table = False
                table_lines = []
            else:
                table_lines.append(line)
    
    # Don't forget last table
    if in_table and table_lines:
        tables[table_title.split(".")[0]] = (table_title, table_lines)
    
    return tables

def clean_cell(text):
    """Clean markdown bold markers from cell text."""
    text = re.sub(r'\*\*(.*?)\*\*', r'\1', text)
    return text

def create_table_docx(table_name, table_title, table_lines, output_path, col_widths=None):
    """Create a standalone DOCX with a single formatted table."""
    doc = Document()
    
    # Page margins
    for section in doc.sections:
        section.top_margin = Cm(2.54)
        section.bottom_margin = Cm(2.54)
        section.left_margin = Cm(2.54)
        section.right_margin = Cm(2.54)
    
    # Title
    title_p = doc.add_paragraph()
    title_p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    title_run = title_p.add_run(table_title)
    title_run.bold = True
    title_run.font.size = Pt(11)
    title_run.font.name = "Times New Roman"
    
    doc.add_paragraph()  # spacing
    
    # Find table data in lines
    headers = None
    rows = []
    footnotes = []
    
    # Find table rows (lines starting with |)
    found_table_start = False
    for line in table_lines:
        stripped = line.strip()
        if stripped.startswith("|"):
            if not found_table_start:
                found_table_start = True
                # This is header or separator, skip separator
            if re.match(r'^\|[\s\-:|]+\|$', stripped):
                continue  # skip separator
            cells = [clean_cell(c.strip()) for c in stripped.split("|")[1:-1]]
            if headers is None:
                headers = cells
            else:
                rows.append(cells)
        elif found_table_start and stripped and not stripped.startswith("|"):
            # Footnote or trailing text
            footnotes.append(stripped)
    
    if headers is None:
        print(f"  WARNING: Could not parse table {table_name}")
        return
    
    print(f"  {table_name}: {len(headers)} cols, {len(rows)} rows")
    
    # Default col widths based on table
    if col_widths is None:
        if table_name == "Table 1":
            col_widths = [Inches(1.5)] + [Inches(1.0)] * 5
        else:
            col_widths = [Inches(1.6)] + [Inches(1.0)] * 6
    
    # Create table
    num_cols = len(headers)
    table = doc.add_table(rows=len(rows) + 1, cols=num_cols)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = True
    
    # Header row
    header_color = "D9E2F3"
    for j, header in enumerate(headers):
        cell = table.rows[0].cells[j]
        cell.text = ""
        p = cell.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run = p.add_run(header)
        run.bold = True
        run.font.size = Pt(9)
        run.font.name = "Times New Roman"
        set_cell_shading(cell, header_color)
    
    # Data rows
    for i, row_data in enumerate(rows):
        for j, cell_text in enumerate(row_data):
            cell = table.rows[i + 1].cells[j]
            cell.text = ""
            p = cell.paragraphs[0]
            # Align: left for first col, center for numeric cols
            if j == 0:
                p.alignment = WD_ALIGN_PARAGRAPH.LEFT
                # Check if bold (sub-header like **Demographics**)
                if cell_text and not cell_text.strip().startswith("*"):
                    run = p.add_run(cell_text)
                    run.font.size = Pt(9)
                    run.font.name = "Times New Roman"
                    if row_data[0].startswith("**") and row_data[0].endswith("**"):
                        run.bold = True
                else:
                    run = p.add_run(cell_text)
                    run.font.size = Pt(9)
                    run.font.name = "Times New Roman"
            else:
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                run = p.add_run(cell_text)
                run.font.size = Pt(9)
                run.font.name = "Times New Roman"
    
    # Add footnotes
    if footnotes:
        doc.add_paragraph()
        for fn in footnotes:
            if fn:
                p = doc.add_paragraph()
                p.paragraph_format.space_after = Pt(2)
                p.paragraph_format.space_before = Pt(0)
                run = p.add_run(fn)
                run.font.size = Pt(8)
                run.font.name = "Times New Roman"
                run.italic = True
    
    doc.save(str(output_path))
    print(f"  Saved: {output_path.name}")

def main():
    tables = extract_tables(MANUSCRIPT)
    
    for table_name in ["Table 1", "Table 2"]:
        if table_name in tables:
            title, lines = tables[table_name]
            output_path = BASE / f"table_{table_name.split()[-1]}_renal_failure.docx"
            create_table_docx(table_name, title, lines, output_path)
        else:
            print(f"  WARNING: {table_name} not found in manuscript")
    
    print("\nDone.")

if __name__ == "__main__":
    main()
