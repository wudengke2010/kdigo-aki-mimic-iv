"""Extract text from 10 non-renal journal PDFs for benchmark analysis."""
import os
import PyPDF2

PDF_DIR = r"I:\6篇mimic renal\非renal杂志"
OUT_DIR = r"C:\Users\admin\WorkBuddy\2026-07-06-13-22-19\bench_nonrenal"
os.makedirs(OUT_DIR, exist_ok=True)

pdfs = [
    "IRNF_42_1788581.pdf",
    "IRNF_47_2562445.pdf",
    "IRNF_46_2313172.pdf",
    "IRNF_45_2282708.pdf",
    "13054_2024_Article_4935.pdf",
    "134_2023_Article_7138.pdf",
    "CJN.14181021.pdf",
    "IRNF_45_2173498.pdf",
    "IRNF_43_1997761.pdf",
    "sfab256.pdf",
]

for pdf_name in pdfs:
    pdf_path = os.path.join(PDF_DIR, pdf_name)
    out_path = os.path.join(OUT_DIR, pdf_name.replace(".pdf", ".txt"))
    print(f"\n{'='*80}\nProcessing: {pdf_name}")
    try:
        with open(pdf_path, "rb") as f:
            reader = PyPDF2.PdfReader(f)
            n_pages = len(reader.pages)
            text_parts = []
            for i, page in enumerate(reader.pages):
                try:
                    txt = page.extract_text() or ""
                except Exception as e:
                    txt = f"[PAGE {i+1} EXTRACT ERROR: {e}]"
                text_parts.append(f"\n--- PAGE {i+1}/{n_pages} ---\n{txt}")
            full_text = "".join(text_parts)
        with open(out_path, "w", encoding="utf-8") as out:
            out.write(full_text)
        # Print first 500 chars to identify the paper
        print(f"  Pages: {n_pages}, Chars: {len(full_text)}")
        print(f"  Saved: {out_path}")
        # Print first 800 chars for quick identification
        preview = full_text[:800].replace("\n", " ")
        print(f"  Preview: {preview}")
    except Exception as e:
        print(f"  ERROR: {e}")

print(f"\n{'='*80}\nDone. All text files in: {OUT_DIR}")
