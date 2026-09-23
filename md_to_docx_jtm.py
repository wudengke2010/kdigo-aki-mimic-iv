#!/usr/bin/env python3
"""Convert manuscript_jtm.md to JTM (BMC) submission DOCX. Reuses CKJ converter."""
from md_to_docx_ckj import convert_md_to_docx

if __name__ == '__main__':
    convert_md_to_docx(
        r'C:\Users\admin\WorkBuddy\2026-07-06-13-22-19\jtm_submission\manuscript_jtm.md',
        r'C:\Users\admin\WorkBuddy\2026-07-06-13-22-19\jtm_submission\manuscript_jtm.docx'
    )
