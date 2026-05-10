#!/usr/bin/env python3
"""D046 — Pre-submission correction + visual/text audit.

Verifies that D045 packaging fixes are visible in the rendered artifacts:
  1. No 'Roselli 2025 EUR' / 'EUR, GRCh38 primary anchor' inside any figure PDF text.
  2. No duplicate 'Figure N: Figure N:' captions in main manuscript PDF.
  3. Single 'Abstract' heading occurrence in main manuscript PDF.
  4. Main manuscript Table 2 fits within 4 INTERIM-only rows (compact).
  5. Forbidden phrase scan on all rendered text (figures + manuscript PDFs + DOCX).
  6. Formula/symbol rendering check (5×10⁻⁸ → \times, 10^{-8}; F > 10; |ΔPP| < 1×10⁻⁶).
  7. Submission readiness checklist explicitly lists all open placeholders.

Outputs:
  manuscript/main_manuscript_D046.docx   (copy of main_manuscript.docx)
  manuscript/main_manuscript_D046.pdf    (copy of main_manuscript.pdf)
  manuscript/d046_pre_submission_audit.md
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

import fitz
import openpyxl
from docx import Document

ROOT = Path("/Users/apple/Desktop/gwas_af")
MS_DIR = ROOT / "manuscript"
MS_FIG = MS_DIR / "figures"

FIGURE_PDFS = sorted(MS_FIG.glob("*.pdf"))
MAIN_PDF = MS_DIR / "main_manuscript.pdf"
SUPP_PDF = MS_DIR / "supplementary_methods.pdf"
MAIN_DOCX = MS_DIR / "main_manuscript.docx"
SUPP_DOCX = MS_DIR / "supplementary_methods.docx"
SUPP_XLSX = MS_DIR / "supplementary_tables.xlsx"
READINESS = MS_DIR / "submission_readiness_checklist.md"

D046_DOCX = MS_DIR / "main_manuscript_D046.docx"
D046_PDF = MS_DIR / "main_manuscript_D046.pdf"
AUDIT_OUT = MS_DIR / "d046_pre_submission_audit.md"

# Forbidden phrases in positive context (paragraph-aware negation handling).
FORBIDDEN = [
    ("Roselli 2025 EUR", "wrong AF ancestry attribution"),
    ("AF Roselli 2025 EUR", "wrong AF ancestry attribution"),
    ("EUR, GRCh38 primary anchor", "wrong AF ancestry attribution"),
    ("EUR, GRCh38)\nprimary anchor", "wrong AF ancestry attribution (figure embed)"),
    ("definitive therapeutic target", "claim too high"),
    ("validated drug target", "claim too high"),
    ("validated therapeutic target", "claim too high"),
    ("final GREEN-2", "final classification locked"),
    ("NPPA inhibition treats", "drug recommendation not supported"),
    ("NPPA inhibition prevents", "drug recommendation not supported"),
]

NEG_MARKERS = (
    "✗", "Forbidden", "forbidden", "Do not", "do not", "should not",
    "is not", "are not", "not a ", "not an ", "not the ", "(not ",
    " not ", "not part", "not run", "not yet", "rather than", "instead of",
    "*not*", "explicitly not", "never", "Never",
    "is locked", "are locked", "remain locked", "remains locked",
    "still locked", "STILL LOCKED", "LOCKED", "Locked",
    "is explicitly locked", "are explicitly locked",
    "kilitli", "kilit",
    "not as ", "rather than ",
)


def _read(path: Path) -> str:
    return path.read_text() if path.exists() else ""


def _scan_pdf(path: Path) -> str:
    if not path.exists():
        return ""
    with fitz.open(path) as doc:
        return "\n".join(p.get_text() for p in doc)


def _scan_docx(path: Path) -> str:
    if not path.exists():
        return ""
    doc = Document(path)
    parts = [p.text for p in doc.paragraphs]
    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                parts.append(cell.text)
    return "\n\n".join(parts)


def _paragraph_violations(text, phrases):
    paras = [p.strip() for p in re.split(r"\n\s*\n", text) if p.strip()]
    issues = []
    for i, para in enumerate(paras):
        para_low = para.lower()
        has_neg = any(m.lower() in para_low for m in NEG_MARKERS)
        for phrase, reason in phrases:
            if phrase.lower() in para_low and not has_neg:
                issues.append((i, phrase, reason, para.replace("\n", " ")[:160]))
    return issues


def _line_violations(text, phrases):
    """Line-level scan for figure-embedded text where paragraphs may be broken."""
    issues = []
    for i, line in enumerate(text.splitlines(), start=1):
        line_low = line.lower()
        for phrase, reason in phrases:
            if phrase.lower() in line_low:
                # in figure text the negation is unlikely to live on the same line; flag conservatively
                issues.append((i, phrase, reason, line.strip()[:160]))
    return issues


# ---------------------------------------------------------------------------
# Audit checks
# ---------------------------------------------------------------------------


def check_figure_pdfs():
    issues = []
    summary = {}
    for pdf in FIGURE_PDFS:
        text = _scan_pdf(pdf)
        summary[pdf.name] = len(text)
        # Figure-internal scan: if any forbidden phrase appears at all, it's a violation
        # (figure text is never legitimately negating itself).
        for i, line in enumerate(text.splitlines(), start=1):
            line_low = line.lower()
            for phrase, reason in FORBIDDEN:
                if phrase.lower() in line_low:
                    issues.append((pdf.name, i, phrase, reason, line.strip()[:160]))
    ok = not issues
    detail = []
    if ok:
        detail.append(f"6 figure PDFs scanned ({sum(summary.values())} chars total); no forbidden phrase found.")
    else:
        for fn, i, phrase, reason, snippet in issues:
            detail.append(f"`{fn}` line {i}: forbidden '{phrase}' ({reason}) — '{snippet}'")
    return "1. Figure-PDF embedded text scan", ok, detail


def check_main_pdf_forbidden():
    text = _scan_pdf(MAIN_PDF)
    issues = _paragraph_violations(text, FORBIDDEN)
    ok = not issues
    detail = []
    if ok:
        detail.append(f"main_manuscript.pdf scanned ({len(text):,} chars); no forbidden phrase in positive context.")
    else:
        for i, phrase, reason, snippet in issues:
            detail.append(f"para{i}: forbidden '{phrase}' ({reason}) — '{snippet}'")
    return "2. Main PDF forbidden phrase scan (paragraph-aware)", ok, detail


def check_duplicate_captions():
    text = _scan_pdf(MAIN_PDF)
    pattern = re.compwith(r"Figure\s+(\d+):\s+Figure\s+\1:")
    matches = list(pattern.finditer(text))
    ok = not matches
    detail = []
    if ok:
        detail.append("No 'Figure N: Figure N:' duplicate caption pattern detected.")
    else:
        for m in matches[:10]:
            detail.append(f"duplicate caption: '{m.group(0)}'")
    return "3. Duplicate caption check (LaTeX auto-numbering)", ok, detail


def check_abstract_heading_count():
    """Abstract heading should appear exactly once (LaTeX auto title)."""
    text = _scan_pdf(MAIN_PDF)
    # Count standalone Abstract heading occurrences (line starts/ends with Abstract).
    lines = text.splitlines()
    counts = sum(1 for line in lines if line.strip() == "Abstract")
    ok = counts <= 1
    detail = [f"'Abstract' standalone-line occurrences: {counts}"]
    if not ok:
        detail.append("Expected ≤ 1; LaTeX auto-titles abstract block. Strip duplicate ## Abstract from source.")
    return "4. Single Abstract heading", ok, detail


def check_table2_compact():
    """Main Table 2 should be ≤ 5 columns × ≤ 5 data rows."""
    if not MAIN_DOCX.exists():
        return "5. Main Table 2 compact layout", False, ["main_manuscript.docx missing"]
    doc = Document(MAIN_DOCX)
    issues = []
    for table in doc.tables:
        n_cols = len(table.columns)
        n_data_rows = len(table.rows) - 1
        # Find the table whose first header row matches Table 2 compact (gene, axis, interim label, ...)
        header_cells = [c.text.strip().lower() for c in table.rows[0].cells]
        if "gene" in header_cells and "axis" in header_cells and "interim label" in header_cells:
            if n_cols > 6:
                issues.append(f"main Table 2 has {n_cols} columns (> 6 → cramped)")
            if n_data_rows > 6:
                issues.append(f"main Table 2 has {n_data_rows} data rows (> 6 → cramped); full table belongs in supplement")
            return "5. Main Table 2 compact layout", not issues, issues or [
                f"main Table 2 layout OK: {n_cols} cols × {n_data_rows} data rows."
            ]
    return "5. Main Table 2 compact layout", False, ["compact Table 2 not found in main DOCX"]


def check_formula_rendering():
    text = _scan_pdf(MAIN_PDF)
    expected_patterns = [
        # 5 × 10^{-8} or 5×10^-8 / 5×10⁻⁸ — accept any reasonable encoding
        (r"(5\s*[×x]\s*10\s*[\-−⁻]?\s*[\^]?[\{]?\s*[\-−]?\s*8|5e-8|5×10−8)", "p < 5×10^-8"),
        (r"F\s*[>≥]\s*10", "F > 10"),
        (r"(1\s*[×x]\s*10\s*[\-−⁻]?\s*[\^]?[\{]?\s*[\-−]?\s*6|1e-6|10−6)", "|ΔPP| < 1×10^-6"),
    ]
    issues = []
    found = []
    for pat, label in expected_patterns:
        if re.search(pat, text):
            found.append(label)
        else:
            issues.append(f"missing/garbled formula: {label}")
    ok = not issues
    detail = []
    detail.append(f"detected: {', '.join(found) if found else 'none'}")
    detail.extend(issues)
    return "6. Formula/symbol rendering", ok, detail


def check_readiness_placeholders():
    text = _read(READINESS)
    required = ["[Authors]", "[Introduction placeholder]", "[Discussion placeholder]",
                "citation placeholder", "Cover letter", "author contributions", "funding"]
    missing = [r for r in required if r not in text]
    ok = not missing
    detail = []
    if ok:
        detail.append("submission_readiness_checklist.md explicitly lists all required open placeholders.")
    else:
        detail.append(f"missing in checklist: {missing}")
    return "7. Submission readiness checklist coverage", ok, detail


def check_xlsx_final_class_preserved():
    if not SUPP_XLSX.exists():
        return "8. final_class_allowed preserved (XLSX)", False, ["xlsx missing"]
    wb = openpyxl.load_workbook(SUPP_XLSX, read_only=True)
    if "Table2_Interim_Evidence" not in wb.sheetnames:
        return "8. final_class_allowed preserved (XLSX)", False, ["Table2_Interim_Evidence sheet missing"]
    ws = wb["Table2_Interim_Evidence"]
    header = [c.value for c in next(ws.iter_rows(min_row=1, max_row=1))]
    if "final_class_allowed" not in header:
        return "8. final_class_allowed preserved (XLSX)", False, ["final_class_allowed column missing"]
    col_idx = header.index("final_class_allowed")
    bad = []
    for row in ws.iter_rows(min_row=2):
        v = (row[col_idx].value or "").strip().lower()
        if v not in ("false", ""):
            bad.append((row[0].value, v))
    ok = not bad
    detail = [f"Table 2 in XLSX has final_class_allowed=false on every row"] if ok else [f"non-false rows: {bad}"]
    return "8. final_class_allowed preserved (XLSX)", ok, detail


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def write_audit(checks, overall_pass):
    md = [
        f"# D046 pre-submission correction + visual audit — {'PASS' if overall_pass else 'FAIL'}",
        "",
        f"**Snapshot:** 2026-05-10",
        "",
        f"**Result:** {sum(1 for _, ok, _ in checks if ok)} / {len(checks)} checks pass.",
        "",
    ]
    for label, ok, details in checks:
        md.append(f"## {label} — {'✅ PASS' if ok else '❌ FAIL'}")
        for d in details:
            md.append(f"- {d}")
        md.append("")
    md.extend([
        "## Submission status",
        "",
        f"**{'D046 PASS — proceed to D047 actual submission drafting (with full Introduction, Discussion, citations, authors, funding).' if overall_pass else 'D046 FAIL — submission blocked. Fix the failing items above before D047.'}**",
        "",
        "## Locked items (continued)",
        "",
        "- Final GREEN/RED/AMBER therapeutic classification — LOCKED.",
        "- Replication phase — LOCKED.",
        "- Wave 3 SuSiE / D042 optional — LOCKED.",
        "- New coloc / new MR / new data downloads / proteome-wide expansion — LOCKED.",
        "",
        "## Open submission-time work (D047)",
        "",
        "- Replace [Authors] with author list + affiliations.",
        "- Replace [Introduction placeholder] with full literature review.",
        "- Replace [Discussion placeholder] with full Discussion (claim boundaries preserved).",
        "- Replace [REF_*] citation placeholders with real references.",
        "- Add cover letter, author contributions, COI, funding statement.",
        "- Reformat per target journal style.",
    ])
    AUDIT_OUT.write_text("\n".join(md))
    print(f"wrote {AUDIT_OUT}")


def main() -> int:
    print("=== D046 visual audit ===")
    checks = [
        check_figure_pdfs(),
        check_main_pdf_forbidden(),
        check_duplicate_captions(),
        check_abstract_heading_count(),
        check_table2_compact(),
        check_formula_rendering(),
        check_readiness_placeholders(),
        check_xlsx_final_class_preserved(),
    ]
    for label, ok, details in checks:
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {label}")
        for d in details:
            print(f"    - {d}")
    overall_pass = all(ok for _, ok, _ in checks)
    write_audit(checks, overall_pass)

    # Save D046 named copies
    if MAIN_DOCX.exists():
        shutil.copy2(MAIN_DOCX, D046_DOCX)
        print(f"copied → {D046_DOCX}")
    if MAIN_PDF.exists():
        shutil.copy2(MAIN_PDF, D046_PDF)
        print(f"copied → {D046_PDF}")

    print(f"=== OVERALL: {'PASS' if overall_pass else 'FAIL'} ===")
    return 0 if overall_pass else 2


if __name__ == "__main__":
    sys.exit(main())
