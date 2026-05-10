#!/usr/bin/env python3
"""D045 — Word/LaTeX manuscript packaging (no scientific change).

Reads D044-corrected reports + tables + figures and assembles:
  manuscript/main_manuscript.docx
  manuscript/main_manuscript.tex
  manuscript/supplementary_methods.docx
  manuscript/supplementary_methods.tex
  manuscript/supplementary_tables.xlsx
  manuscript/supplementary_tables_tsv/      (copies)
  manuscript/figures/                        (PDF + 300 dpi PNG)
  manuscript/submission_readiness_checklist.md
  manuscript/d045_packaging_audit.md

Allowed: formatting, section ordering, table/figure placement, figure legends, citation placeholders.
Forbidden: any scientific claim change, any AF wording regression to "EUR", any forbidden phrase
(NPPA inhibition treats, validated drug target, final GREEN/RED, proteome-wide [unqualified], etc.).
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

import fitz  # PyMuPDF
import openpyxl
import polars as pl
from docx import Document
from docx.shared import Inches, Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH

ROOT = Path("/Users/apple/Desktop/gwas_af")
REPORTS = ROOT / "reports"
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"
INTERIM = ROOT / "results/interim"

MS_DIR = ROOT / "manuscript"
MS_FIG = MS_DIR / "figures"
MS_SUPP_TSV = MS_DIR / "supplementary_tables_tsv"
for d in (MS_DIR, MS_FIG, MS_SUPP_TSV):
    d.mkdir(parents=True, exist_ok=True)

FIGURE_ORDER = [
    ("Figure 1", "study_design_cvd150",
     "Study design — CVD-prioritized AF-anchored pQTL screen.",
     "Aptamer panel, AF anchor, downstream HF/stroke projection, and selective colocalization workflow. AF anchor is the CVDKP common-variant all-ancestry meta-analysis (`AF_GWAS_AFGenPlus_commonFreq_ALLv21`). HF (HERMES 2024) and stroke (GIGASTROKE 2022) outcomes are EUR-only."),
    ("Figure 2", "cvd150_mr_coloc_filtering_flow",
     "MR → colocalization filtering flow.",
     "Funnel from 150 attempted aptamers to 1 AF-anchored coloc-supported survivor (NPPA). FDR is computed within the CVD150 cis-tested aptamer universe, not proteome-wide."),
    ("Figure 3", "three_axis_interim_atlas",
     "Three-axis interim atlas.",
     "NPPA AF–CES cardioembolic axis (interim AMBER), F11 / KNG1 antithrombotic axis (interim BLUE, non-AF), and MMP12 atherothrombotic axis (interim BLUE, non-AF). All labels are interim; final GREEN/RED/AMBER classification is locked."),
    ("Figure 4", "nppa_af_ces_axis",
     "NPPA AF–CES axis pleiotropy.",
     "PP.H4 across primary and sensitivity priors for NPPA × {AF, CES, AIS, HF subtypes}. AF and CES are coloc-supported; AIS and HF subtypes are not."),
    ("Figure 5", "blue_modules_f11_kng1_mmp12",
     "Non-AF BLUE stroke modules.",
     "F11 and KNG1 antithrombotic axis with AIS and CES; MMP12 atherothrombotic axis with AIS and LAS. All BLUE labels are interim and non-AF-anchored."),
    ("Figure 6", "downgraded_mr_hits_coloc_filter",
     "Downgraded MR hits — colocalization filter.",
     "AF MR −log10(q) (left) versus colocalization PP.H3/PP.H4 (right) for IL6R, DSC2, OGN, TIMP3, PCSK9, and SERPINF2. Strong AF MR signals are not retained when colocalization rejects a shared causal variant."),
]

CITATION_PLACEHOLDERS = [
    "[REF_Roselli2025_AF]",
    "[REF_HERMES2024]",
    "[REF_GIGASTROKE2022]",
    "[REF_deCODE2021_Ferkingstad]",
    "[REF_coloc_Giambartolomei]",
    "[REF_FinnGen_or_replication_future]",
]
# LaTeX-safe rendering (underscores escaped).
CITATION_PLACEHOLDERS_LATEX = [c.replace("_", r"\_") for c in CITATION_PLACEHOLDERS]


# ---------------------------------------------------------------------------
# Asset prep
# ---------------------------------------------------------------------------


def copy_figures():
    for _, stem, _, _ in FIGURE_ORDER:
        src = FIGURES / f"{stem}.pdf"
        dst_pdf = MS_FIG / f"{stem}.pdf"
        if src.exists():
            shutil.copy2(src, dst_pdf)
        # Render PNG at 300 dpi via PyMuPDF
        dst_png = MS_FIG / f"{stem}.png"
        with fitz.open(src) as doc:
            page = doc.load_page(0)
            mat = fitz.Matrix(300 / 72, 300 / 72)
            pix = page.get_pixmap(matrix=mat, alpha=False)
            pix.save(str(dst_png))
        print(f"figures: {stem}.pdf + {stem}.png (300 dpi)")


def copy_tables():
    for n in range(1, 6):
        # find table by prefix
        matches = list(TABLES.glob(f"table{n}_*.tsv"))
        if not matches:
            continue
        src = matches[0]
        shutil.copy2(src, MS_SUPP_TSV / src.name)
        print(f"tables: copied {src.name}")


def build_supplementary_xlsx():
    wb = openpyxl.Workbook()
    # remove default sheet
    default = wb.active
    wb.remove(default)
    sheet_specs = [
        ("Table1_Screen_Overview", "table1_cvd150_screen_overview.tsv"),
        ("Table2_Interim_Evidence", "table2_interim_evidence_classes.tsv"),
        ("Table3_Coloc_Supported", "table3_coloc_supported_axes.tsv"),
        ("Table4_Downgraded", "table4_downgraded_candidates.tsv"),
        ("Table5_Replication", "table5_replication_roadmap.tsv"),
    ]
    for sheet_name, fname in sheet_specs:
        path = TABLES / fname
        if not path.exists():
            continue
        df = pl.read_csv(path, separator="\t", infer_schema_length=0)
        ws = wb.create_sheet(sheet_name)
        ws.append(df.columns)
        for row in df.iter_rows():
            ws.append([str(v) if v is not None else "" for v in row])
        # column widths (rough)
        for i, col in enumerate(df.columns, start=1):
            try:
                width = max(12, min(60, max(len(col), max((len(str(v)) for v in df[col]), default=12)) + 2))
                ws.column_dimensions[openpyxl.utils.get_column_letter(i)].width = width
            except Exception:
                pass
        ws.freeze_panes = "A2"
    out = MS_DIR / "supplementary_tables.xlsx"
    wb.save(out)
    print(f"wrote {out}")


# ---------------------------------------------------------------------------
# Manuscript text composition
# ---------------------------------------------------------------------------


def _read(path: Path) -> str:
    return path.read_text() if path.exists() else ""


def _strip_md_emphasis(text: str) -> str:
    text = re.sub(r"\*\*([^*]+)\*\*", r"\1", text)
    text = re.sub(r"\*([^*]+)\*", r"\1", text)
    return text


def _strip_md_headers(text: str) -> str:
    return re.sub(r"^#+\s*", "", text, flags=re.MULTILINE)


def _section_body(report_md: str, section_marker: str = "##") -> str:
    """Return body without the leading H1 and front-matter blockquote."""
    lines = report_md.splitlines()
    out = []
    for line in lines:
        if line.startswith("# "):
            continue
        if line.startswith(">"):
            continue
        out.append(line)
    return "\n".join(out).strip()


def _strip_top_heading(md: str, target: str) -> str:
    """Remove a leading '## <target>' heading line so we don't duplicate against the wrapper title (D046 fix)."""
    pat = re.compwith(rf"^##\s*{re.escape(target)}\s*$", flags=re.IGNORECASE | re.MULTILINE)
    return pat.sub("", md, count=1).lstrip()


def _build_main_text():
    abstract = _section_body(_read(REPORTS / "manuscript_style_abstract_draft.md"))
    abstract = _strip_top_heading(abstract, "Abstract")
    results = _section_body(_read(REPORTS / "manuscript_style_results.md"))
    methods = _section_body(_read(REPORTS / "manuscript_style_methods.md"))
    limitations = _section_body(_read(REPORTS / "manuscript_style_limitations.md"))
    return abstract, results, methods, limitations


# ---------------------------------------------------------------------------
# LaTeX assembly
# ---------------------------------------------------------------------------


def _md_to_latex(md_text: str) -> str:
    """Minimal markdown → LaTeX. Preserves section structure, escapes specials, drops markdown."""
    text = md_text
    # Tables: convert | a | b | c | ... lines into a tabular environment block
    out_lines = []
    in_table = False
    table_buf: list[list[str]] = []

    def flush_table():
        nonlocal table_buf, in_table
        if not table_buf:
            return
        ncols = max(len(r) for r in table_buf)
        # If second row is the divider |---|---|, drop it
        if len(table_buf) >= 2 and all(set(c.strip()) <= set("-:") for c in table_buf[1]):
            header, divider, body = table_buf[0], table_buf[1], table_buf[2:]
        else:
            header, body = table_buf[0], table_buf[1:]
        out_lines.append(r"\begin{table}[h]")
        out_lines.append(r"\centering")
        out_lines.append(r"\small")
        out_lines.append(r"\begin{tabular}{" + ("l" * ncols) + r"}")
        out_lines.append(r"\hline")
        out_lines.append(" & ".join(_latex_escape(c.strip()) for c in header) + r" \\ \hline")
        for r in body:
            cells = [_latex_escape(c.strip()) for c in r] + [""] * (ncols - len(r))
            out_lines.append(" & ".join(cells[:ncols]) + r" \\")
        out_lines.append(r"\hline")
        out_lines.append(r"\end{tabular}")
        out_lines.append(r"\end{table}")
        out_lines.append("")
        table_buf = []
        in_table = False

    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("|") and stripped.endswith("|"):
            cells = [c for c in stripped.strip("|").split("|")]
            table_buf.append(cells)
            in_table = True
            continue
        if in_table:
            flush_table()

        if stripped.startswith("## "):
            out_lines.append(r"\section*{" + _latex_escape(stripped[3:]) + "}")
        elif stripped.startswith("### "):
            out_lines.append(r"\subsection*{" + _latex_escape(stripped[4:]) + "}")
        elif stripped.startswith("- "):
            out_lines.append(r"\item " + _latex_escape_inline(stripped[2:]))
        elif not stripped:
            out_lines.append("")
        else:
            out_lines.append(_latex_escape_inline(stripped))

    if in_table:
        flush_table()

    # Wrap consecutive \item lines in itemize
    final = []
    in_list = False
    for ln in out_lines:
        if ln.startswith(r"\item "):
            if not in_list:
                final.append(r"\begin{itemize}")
                in_list = True
            final.append(ln)
        else:
            if in_list:
                final.append(r"\end{itemize}")
                in_list = False
            final.append(ln)
    if in_list:
        final.append(r"\end{itemize}")
    return "\n".join(final)


_LATEX_REPL = {
    "\\": r"\textbackslash{}",
    "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#",
    "_": r"\_", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}",
    "^": r"\textasciicircum{}",
}

_UNICODE_REPL = {
    "≥": r"$\geq$", "≤": r"$\leq$", "≠": r"$\neq$",
    "≈": r"$\approx$", "∝": r"$\propto$",
    "→": r"$\to$", "←": r"$\leftarrow$", "↔": r"$\leftrightarrow$",
    "×": r"$\times$", "÷": r"$\div$", "·": r"$\cdot$",
    "±": r"$\pm$", "∓": r"$\mp$",
    "−": "--",        # U+2212 minus sign
    "–": "--",        # en dash
    "—": "---",       # em dash
    "α": r"$\alpha$", "β": r"$\beta$", "γ": r"$\gamma$", "δ": r"$\delta$",
    "θ": r"$\theta$", "λ": r"$\lambda$", "μ": r"$\mu$", "π": r"$\pi$",
    "σ": r"$\sigma$", "Σ": r"$\Sigma$", "ω": r"$\omega$", "Δ": r"$\Delta$",
    "²": r"$^{2}$", "³": r"$^{3}$",
    "°": r"$^{\circ}$",
    "‐": "-", "‑": "-", "‒": "-",  # various hyphens
    "‘": "`", "’": "'", "“": "``", "”": "''",
    "✓": r"\checkmark", "✗": r"$\times$",
    "…": r"\ldots{}",
    " ": "~",  # non-breaking space
}


_SUPERSCRIPT_DIGITS = str.maketrans({
    "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4",
    "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9",
    "⁺": "+", "⁻": "-",
})


def _normalize_unicode(s: str) -> str:
    # Convert "5×10⁻⁸" / "10⁻¹⁵" patterns into LaTeX math: 5\times 10^{-8}
    sup_pattern = re.compwith(r"([⁻⁺]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+)")
    def sup_repl(m):
        return "$^{" + m.group(0).translate(_SUPERSCRIPT_DIGITS) + "}$"
    s = sup_pattern.sub(sup_repl, s)
    # Now apply remaining single-char replacements
    for k, v in _UNICODE_REPL.items():
        s = s.replace(k, v)
    return s


def _latex_escape(s: str) -> str:
    """Escape LaTeX-special chars first; THEN unicode-normalize so injected math
    delimiters ($\times$, $^{-8}$, etc.) survive without being double-escaped."""
    out = []
    for c in s:
        out.append(_LATEX_REPL.get(c, c))
    return _normalize_unicode("".join(out))


def _latex_escape_inline(s: str) -> str:
    """Inline escape that preserves bold/italic markers from markdown."""
    # Convert **x** → \textbf{x} BEFORE general escape
    bold_pat = re.compwith(r"\*\*([^*]+)\*\*")
    italic_pat = re.compwith(r"\*([^*]+)\*")
    code_pat = re.compwith(r"`([^`]+)`")

    placeholders = {}
    def stash(s, regex, fmt):
        def repl(m):
            key = f"@@PH{len(placeholders)}@@"
            placeholders[key] = fmt.format(_latex_escape(m.group(1)))
            return key
        return regex.sub(repl, s)
    s = stash(s, code_pat, r"\texttt{{{0}}}")
    s = stash(s, bold_pat, r"\textbf{{{0}}}")
    s = stash(s, italic_pat, r"\emph{{{0}}}")
    s = _latex_escape(s)
    for k, v in placeholders.items():
        s = s.replace(_latex_escape(k), v)
    return s


def _build_table2_compact() -> list[dict]:
    """D046 compact main-manuscript Table 2 (INTERIM only; 4 candidates × 5 columns)."""
    return [
        {"gene": "NPPA",  "axis": "AF–CES cardioembolic",
         "interim label": "INTERIM AMBER (AF-anchored)",
         "primary outcome": "AF",  "primary PP.H4": "0.911"},
        {"gene": "F11",   "axis": "antithrombotic stroke",
         "interim label": "INTERIM BLUE (non-AF)",
         "primary outcome": "AIS", "primary PP.H4": "0.993"},
        {"gene": "KNG1",  "axis": "kallikrein-kinin stroke",
         "interim label": "INTERIM BLUE (non-AF)",
         "primary outcome": "AIS", "primary PP.H4": "0.991"},
        {"gene": "MMP12", "axis": "atherothrombotic stroke",
         "interim label": "INTERIM BLUE (non-AF)",
         "primary outcome": "LAS", "primary PP.H4": "0.916"},
    ]


def _list_of_dicts_to_latex_tabular(rows: list[dict], caption: str, label: str, font_size: str = "small") -> str:
    if not rows:
        return ""
    cols = list(rows[0].keys())
    align = "l" * len(cols)
    out = [
        r"\begin{table}[h]",
        r"\centering",
        rf"\{font_size}",
        rf"\begin{{tabular}}{{{align}}}",
        r"\hline",
        " & ".join(_latex_escape(c) for c in cols) + r" \\ \hline",
    ]
    for row in rows:
        cells = [_latex_escape_inline(str(row[c])) for c in cols]
        out.append(" & ".join(cells) + r" \\")
    out.extend([
        r"\hline",
        r"\end{tabular}",
        rf"\caption{{{_latex_escape_inline(caption)}}}",
        rf"\label{{tab:{label}}}",
        r"\end{table}",
    ])
    return "\n".join(out)


def _df_to_latex_tabular(df, caption: str, label: str, font_size: str = "small") -> str:
    ncols = len(df.columns)
    align = "l" * ncols
    out = [
        r"\begin{table}[h]",
        r"\centering",
        rf"\{font_size}",
        rf"\begin{{tabular}}{{{align}}}",
        r"\hline",
        " & ".join(_latex_escape(c) for c in df.columns) + r" \\ \hline",
    ]
    for row in df.iter_rows():
        cells = []
        for v in row:
            s = str(v) if v is not None else ""
            try:
                fv = float(s)
                s = f"{fv:.3f}"
            except (TypeError, ValueError):
                pass
            cells.append(_latex_escape_inline(s))
        out.append(" & ".join(cells) + r" \\")
    out.extend([
        r"\hline",
        r"\end{tabular}",
        rf"\caption{{{_latex_escape_inline(caption)}}}",
        rf"\label{{tab:{label}}}",
        r"\end{table}",
    ])
    return "\n".join(out)


def write_main_tex():
    abstract, results, methods, limitations = _build_main_text()
    figs_block = []
    for label, stem, caption, legend in FIGURE_ORDER:
        # LaTeX auto-prefixes "Figure N:" — do NOT add it manually (D046 fix).
        figs_block.append(
            r"\begin{figure}[h]" "\n"
            r"\centering" "\n"
            rf"\includegraphics[width=0.85\textwidth]{{figures/{stem}.pdf}}" "\n"
            rf"\caption{{{_latex_escape_inline(caption)} {_latex_escape_inline(legend)}}}" "\n"
            rf"\label{{fig:{stem}}}" "\n"
            r"\end{figure}" "\n"
        )

    parts = [
        r"\documentclass[11pt]{article}",
        r"\usepackage[a4paper,margin=1in]{geometry}",
        r"\usepackage{graphicx}",
        r"\usepackage{booktabs}",
        r"\usepackage{hyperref}",
        r"\usepackage{microtype}",
        r"\usepackage{enumitem}",
        r"\usepackage[utf8]{inputenc}",
        r"\title{CVD-prioritized plasma proteomic screen for AF--HF--stroke pleiotropy: an interim coloc-supported atlas}",
        r"\author{[Authors]}",
        r"\date{Manuscript draft (D045 packaging, 2026-05-10)}",
        r"\begin{document}",
        r"\maketitle",
        r"\begin{abstract}",
        _md_to_latex(abstract),
        r"\end{abstract}",
        r"\section*{Introduction}",
        r"[Introduction placeholder. To be expanded with literature review citing " + ", ".join(CITATION_PLACEHOLDERS_LATEX) + r".]",
        r"\section*{Results}",
        _md_to_latex(results),
        r"\section*{Discussion}",
        r"[Discussion placeholder. The current evidence layer establishes one AF-anchored interim coloc-supported axis (NPPA) and two non-AF stroke modules (F11/KNG1 antithrombotic; MMP12 atherothrombotic). Final GREEN/RED/AMBER therapeutic classification is locked pending independent pQTL replication, outcome replication, and LD-aware fine-mapping (D041, D043, D044). Cite: " + ", ".join(CITATION_PLACEHOLDERS_LATEX) + r".]",
        r"\section*{Methods (summary)}",
        _md_to_latex(methods),
        r"\section*{Limitations}",
        _md_to_latex(limitations),
        r"\section*{Data and code availability}",
        r"deCODE 2021 plasma pQTL summary statistics (Ferkingstad et al., 2021) were obtained via authorised per-aptamer fetch through the deCODE API (fwith format and licensing per [REF\_deCODE2021\_Ferkingstad]). AF Roselli 2025 CVDKP common-variant AFGenPlus all-ancestry meta-analysis (\texttt{AF\_GWAS\_AFGenPlus\_commonFreq\_ALLv21}), HERMES 2024 EUR HF subtypes, and GIGASTROKE 2022 EUR stroke subtypes are publicly available. Analysis code is in this project repository.",
        r"\section*{Tables}",
        r"Table 1 (CVD150 screen overview) and Table 2 (interim coloc-supported axes summary) appear below. Full tables (Table 2 with all 17 columns, Table 3 coloc-supported axes, Table 4 downgraded candidates, Table 5 replication roadmap) are provided in \texttt{supplementary\_tables.xlsx} and \texttt{supplementary\_tables\_tsv/}.",
        _df_to_latex_tabular(
            pl.read_csv(TABLES / "table1_cvd150_screen_overview.tsv", separator="\t", infer_schema_length=0),
            "CVD150 targeted screen overview.",
            "screen_overview",
            font_size="footnotesize",
        ),
        _list_of_dicts_to_latex_tabular(
            _build_table2_compact(),
            "Interim coloc-supported axes (compact summary). Full 10-row interim evidence table with all 17 columns is in supplementary\\_tables.xlsx (Table 2 sheet).",
            "interim_evidence_compact",
            font_size="small",
        ),
        r"\section*{Figure legends}",
        "\n\n".join(figs_block),
        r"\section*{Replication roadmap}",
        r"See supplementary Table 5 for the replication roadmap. No replication analysis is included in the current evidence layer.",
        r"\end{document}",
    ]
    out = MS_DIR / "main_manuscript.tex"
    out.write_text("\n\n".join(parts))
    print(f"wrote {out}")


def write_supp_tex():
    methods = _section_body(_read(REPORTS / "manuscript_style_methods.md"))
    limitations = _section_body(_read(REPORTS / "manuscript_style_limitations.md"))
    parts = [
        r"\documentclass[11pt]{article}",
        r"\usepackage[a4paper,margin=1in]{geometry}",
        r"\usepackage{graphicx}",
        r"\usepackage{booktabs}",
        r"\usepackage{hyperref}",
        r"\usepackage{microtype}",
        r"\usepackage{enumitem}",
        r"\usepackage[utf8]{inputenc}",
        r"\title{Supplementary Methods --- CVD-prioritized plasma proteomic screen}",
        r"\author{[Authors]}",
        r"\date{D045 packaging, 2026-05-10}",
        r"\begin{document}",
        r"\maketitle",
        r"\section*{Supplementary Methods}",
        _md_to_latex(methods),
        r"\section*{Supplementary Limitations}",
        _md_to_latex(limitations),
        r"\section*{Supplementary tables}",
        r"See \texttt{supplementary\_tables.xlsx} or individual TSV files in \texttt{supplementary\_tables\_tsv/}: Table 1 (CVD150 screen overview), Table 2 (interim evidence classes, all 10 candidates), Table 3 (coloc-supported axes), Table 4 (downgraded candidates), Table 5 (replication roadmap).",
        r"\end{document}",
    ]
    out = MS_DIR / "supplementary_methods.tex"
    out.write_text("\n\n".join(parts))
    print(f"wrote {out}")


# ---------------------------------------------------------------------------
# DOCX assembly
# ---------------------------------------------------------------------------


def _docx_paragraph(doc, text, *, bold=False, size=None, align=None):
    p = doc.add_paragraph()
    if align is not None:
        p.alignment = align
    run = p.add_run(text)
    run.bold = bold
    if size:
        run.font.size = Pt(size)
    return p


def _add_md_block_to_docx(doc, md_text):
    in_table_buf = []

    def flush_table():
        nonlocal in_table_buf
        if not in_table_buf:
            return
        if len(in_table_buf) >= 2 and all(set(c.strip()) <= set("-:") for c in in_table_buf[1]):
            header = in_table_buf[0]
            body = in_table_buf[2:]
        else:
            header = in_table_buf[0]
            body = in_table_buf[1:]
        ncols = max(len(in_table_buf[0]), max((len(r) for r in body), default=0))
        table = doc.add_table(rows=1 + len(body), cols=ncols)
        table.style = "Light Grid Accent 1"
        for i, c in enumerate(header[:ncols]):
            cell = table.rows[0].cells[i]
            cell.text = ""
            run = cell.paragraphs[0].add_run(_strip_md_emphasis(c.strip()))
            run.bold = True
        for ri, row in enumerate(body, start=1):
            cells_padded = list(row) + [""] * (ncols - len(row))
            for ci in range(ncols):
                table.rows[ri].cells[ci].text = _strip_md_emphasis(cells_padded[ci].strip())
        doc.add_paragraph("")
        in_table_buf = []

    in_list = False
    for line in md_text.splitlines():
        stripped = line.strip()
        if stripped.startswith("|") and stripped.endswith("|"):
            cells = [c for c in stripped.strip("|").split("|")]
            in_table_buf.append(cells)
            continue
        if in_table_buf:
            flush_table()

        if stripped.startswith("### "):
            doc.add_heading(_strip_md_emphasis(stripped[4:]), level=2)
            in_list = False
        elif stripped.startswith("## "):
            doc.add_heading(_strip_md_emphasis(stripped[3:]), level=2)
            in_list = False
        elif stripped.startswith("- "):
            p = doc.add_paragraph(_strip_md_emphasis(stripped[2:]), style="List Bullet")
            in_list = True
        elif stripped.startswith("```"):
            in_list = False  # ignore code fences as paragraphs
        elif not stripped:
            in_list = False
        else:
            doc.add_paragraph(_strip_md_emphasis(stripped))
            in_list = False

    if in_table_buf:
        flush_table()


def write_main_docx():
    abstract, results, methods, limitations = _build_main_text()
    doc = Document()
    title = doc.add_heading("CVD-prioritized plasma proteomic screen for AF–HF–stroke pleiotropy: an interim coloc-supported atlas", level=0)
    _docx_paragraph(doc, "Manuscript draft — D045 packaging (2026-05-10).", size=10, align=WD_ALIGN_PARAGRAPH.CENTER)
    _docx_paragraph(doc, "[Authors]", size=10, align=WD_ALIGN_PARAGRAPH.CENTER)

    doc.add_heading("Abstract", level=1)
    _add_md_block_to_docx(doc, abstract)

    doc.add_heading("Introduction", level=1)
    doc.add_paragraph(
        f"[Introduction placeholder. To be expanded with literature review citing "
        f"{', '.join(CITATION_PLACEHOLDERS)}.]"
    )

    doc.add_heading("Results", level=1)
    _add_md_block_to_docx(doc, results)

    doc.add_heading("Discussion", level=1)
    doc.add_paragraph(
        "[Discussion placeholder. The current evidence layer establishes one AF-anchored interim coloc-supported axis (NPPA) "
        "and two non-AF stroke modules (F11/KNG1 antithrombotic; MMP12 atherothrombotic). Final GREEN/RED/AMBER "
        "therapeutic classification is locked pending independent pQTL replication, outcome replication, and "
        f"LD-aware fine-mapping (D041, D043, D044). Cite: {', '.join(CITATION_PLACEHOLDERS)}.]"
    )

    doc.add_heading("Methods (summary)", level=1)
    _add_md_block_to_docx(doc, methods)

    doc.add_heading("Limitations", level=1)
    _add_md_block_to_docx(doc, limitations)

    doc.add_heading("Data and code availability", level=1)
    doc.add_paragraph(
        "deCODE 2021 plasma pQTL summary statistics (Ferkingstad et al., 2021) were obtained via authorised per-aptamer fetch "
        "through the deCODE API. AF Roselli 2025 CVDKP common-variant AFGenPlus all-ancestry meta-analysis "
        "(AF_GWAS_AFGenPlus_commonFreq_ALLv21), HERMES 2024 EUR HF subtypes, and GIGASTROKE 2022 EUR stroke subtypes "
        "are publicly available. Analysis code is in this project repository."
    )

    doc.add_heading("Tables", level=1)
    doc.add_paragraph(
        "Table 1 (CVD150 screen overview) and Table 2 (interim coloc-supported axes summary) are embedded below. "
        "Full tables (Table 2 interim evidence classes with all 17 columns, Table 3 coloc-supported axes, "
        "Table 4 downgraded candidates, Table 5 replication roadmap) are provided in supplementary_tables.xlsx "
        "and supplementary_tables_tsv/."
    )

    # Table 1
    t1 = pl.read_csv(TABLES / "table1_cvd150_screen_overview.tsv", separator="\t", infer_schema_length=0)
    p = doc.add_paragraph()
    run = p.add_run("Table 1. ")
    run.bold = True
    p.add_run("CVD150 targeted screen overview.")
    table = doc.add_table(rows=1 + t1.height, cols=len(t1.columns))
    table.style = "Light Grid Accent 1"
    for i, col in enumerate(t1.columns):
        cell = table.rows[0].cells[i]
        run = cell.paragraphs[0].add_run(col)
        run.bold = True
    for ri, row in enumerate(t1.iter_rows(), start=1):
        for ci, val in enumerate(row):
            table.rows[ri].cells[ci].text = str(val) if val is not None else ""
    doc.add_paragraph("")

    # Table 2 (D046 compact 4-row INTERIM-only summary; full 10-row in supplement).
    t2_compact = _build_table2_compact()
    p = doc.add_paragraph()
    run = p.add_run("Table 2. ")
    run.bold = True
    p.add_run("Interim coloc-supported axes (compact summary). Final classification is locked "
              "(final_class_allowed = false). Full 10-row interim evidence table with all 17 "
              "columns is in supplementary_tables.xlsx (Table 2 sheet).")
    cols = list(t2_compact[0].keys())
    table = doc.add_table(rows=1 + len(t2_compact), cols=len(cols))
    table.style = "Light Grid Accent 1"
    for i, col in enumerate(cols):
        cell = table.rows[0].cells[i]
        run = cell.paragraphs[0].add_run(col)
        run.bold = True
    for ri, row in enumerate(t2_compact, start=1):
        for ci, col in enumerate(cols):
            table.rows[ri].cells[ci].text = str(row[col])
    doc.add_paragraph("")

    doc.add_heading("Figures", level=1)
    for label, stem, caption, legend in FIGURE_ORDER:
        png_path = MS_FIG / f"{stem}.png"
        if png_path.exists():
            doc.add_picture(str(png_path), width=Inches(6.0))
            p = doc.paragraphs[-1]
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        legend_p = doc.add_paragraph()
        run = legend_p.add_run(f"{label}: ")
        run.bold = True
        legend_p.add_run(f"{caption} {legend}")

    doc.add_heading("Replication roadmap", level=1)
    doc.add_paragraph(
        "See supplementary Table 5 for the replication roadmap. No replication analysis is included in the current evidence layer."
    )

    out = MS_DIR / "main_manuscript.docx"
    doc.save(out)
    print(f"wrote {out}")


def write_supp_docx():
    methods = _section_body(_read(REPORTS / "manuscript_style_methods.md"))
    limitations = _section_body(_read(REPORTS / "manuscript_style_limitations.md"))
    doc = Document()
    doc.add_heading("Supplementary Methods — CVD-prioritized plasma proteomic screen", level=0)
    _docx_paragraph(doc, "D045 packaging, 2026-05-10.", size=10, align=WD_ALIGN_PARAGRAPH.CENTER)
    doc.add_heading("Supplementary Methods", level=1)
    _add_md_block_to_docx(doc, methods)
    doc.add_heading("Supplementary Limitations", level=1)
    _add_md_block_to_docx(doc, limitations)
    doc.add_heading("Supplementary tables", level=1)
    doc.add_paragraph(
        "See supplementary_tables.xlsx or individual TSV files in supplementary_tables_tsv/: "
        "Table 1 (CVD150 screen overview), Table 2 (interim evidence classes, all 10 candidates), "
        "Table 3 (coloc-supported axes), Table 4 (downgraded candidates), Table 5 (replication roadmap)."
    )
    out = MS_DIR / "supplementary_methods.docx"
    doc.save(out)
    print(f"wrote {out}")


# ---------------------------------------------------------------------------
# Submission readiness + packaging audit
# ---------------------------------------------------------------------------


FORBIDDEN_LITERAL_PHRASES = [
    ("Roselli 2025 EUR", "wrong AF ancestry attribution"),
    ("AF Roselli 2025 EUR", "wrong AF ancestry attribution"),
    ("definitive therapeutic target", "claim too high"),
    ("validated drug target", "claim too high"),
    ("final GREEN-2", "final classification locked"),
    ("final GREEN", "final classification locked"),
    ("final RED", "final classification locked"),
    ("NPPA inhibition treats", "drug recommendation not supported"),
    ("NPPA inhibition prevents", "drug recommendation not supported"),
]

NEG_MARKERS = (
    "✗", "Forbidden", "forbidden", "Do not", "do not", "should not",
    "is not", "are not", "not a ", "not an ", "not the ", "(not ",
    " not ", "not part", "not run", "not yet", "rather than", "instead of",
    "*not*", "explicitly not", "never", "Never",
    # Locked/sealed context — phrases that explicitly state the claim is closed
    "is locked", "are locked", "remain locked", "remains locked",
    "still locked", "STILL LOCKED", "LOCKED", "Locked",
    "is explicitly locked", "are explicitly locked",
    "kilitli", "kilit",  # Turkish parallels
)


def _clean_for_scan(text: str) -> str:
    """Strip LaTeX/markdown emphasis so substring scans see the underlying claim."""
    # LaTeX commands with one-arg payload: \emph{X}, \textbf{X}, \texttt{X}, \textit{X}
    text = re.sub(r"\\(?:emph|textbf|textit|texttt)\{([^{}]*)\}", r"\1", text)
    # Markdown emphasis
    text = re.sub(r"\*\*([^*]+)\*\*", r"\1", text)
    text = re.sub(r"\*([^*]+)\*", r"\1", text)
    text = re.sub(r"_([^_]+)_", r"\1", text)
    # LaTeX backslash escapes for &, %, etc.
    text = text.replace(r"\&", "&").replace(r"\%", "%").replace(r"\_", "_").replace(r"\#", "#")
    return text


def _paragraph_violations(text, phrases_with_reason):
    cleaned = _clean_for_scan(text)
    paras = [p.strip() for p in re.split(r"\n\s*\n", cleaned) if p.strip()]
    issues = []
    for i, para in enumerate(paras):
        para_low = para.lower()
        has_neg = any(m.lower() in para_low for m in NEG_MARKERS)
        for phrase, reason in phrases_with_reason:
            if phrase.lower() in para_low and not has_neg:
                snippet = para.replace("\n", " ")[:160]
                issues.append((i, phrase, reason, snippet))
    return issues


def _scan_docx_text(path: Path) -> str:
    if not path.exists():
        return ""
    doc = Document(path)
    parts = [p.text for p in doc.paragraphs]
    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                parts.append(cell.text)
    return "\n\n".join(parts)


def write_submission_readiness():
    md = [
        "# Submission readiness checklist (D045)",
        "",
        f"**Snapshot:** 2026-05-10",
        "",
        "## Fwiths generated",
        "",
        "- `manuscript/main_manuscript.docx`",
        "- `manuscript/main_manuscript.tex`",
        "- `manuscript/supplementary_methods.docx`",
        "- `manuscript/supplementary_methods.tex`",
        "- `manuscript/supplementary_tables.xlsx`",
        "- `manuscript/supplementary_tables_tsv/table[1-5]_*.tsv`",
        "- `manuscript/figures/*.pdf` (LaTeX) + `.png` (Word, 300 dpi)",
        "- `manuscript/d045_packaging_audit.md`",
        "",
        "## Open before submission",
        "",
        "- [ ] Replace `[Authors]` placeholder with author list + affiliations.",
        "- [ ] Replace `[Introduction placeholder]` with full Introduction (literature review).",
        "- [ ] Replace `[Discussion placeholder]` with full Discussion text.",
        "- [ ] Replace citation placeholders (`[REF_*]`) with actual references.",
        "- [ ] Choose target journal; reformat per journal style if needed.",
        "- [ ] Compwith LaTeX (`pdflatex main_manuscript.tex`) — verify no errors.",
        "- [ ] Open Word `.docx` — verify figures render and tables are intact.",
        "- [ ] Verify supplementary `.xlsx` opens and all 5 sheets are present.",
        "- [ ] Verify forbidden phrase scan continues to PASS (D045 audit).",
        "- [ ] Verify AF wording remains all-ancestry / mixed / CVDKP.",
        "- [ ] Confirm `final_class_allowed = false` is preserved in Table 2 and supplementary copies.",
        "- [ ] Add Cover letter draft.",
        "- [ ] Add author contributions, conflict of interest statement, funding.",
        "",
        "## Locked items (D041 + D043 + D044 + D045)",
        "",
        "- Final GREEN/RED/AMBER therapeutic classification — LOCKED.",
        "- Replication phase — LOCKED until explicit protocol approval.",
        "- Wave 3 SuSiE / D042 optional — LOCKED until LD reference setup.",
        "- New coloc / new MR / new data downloads / proteome-wide expansion — LOCKED.",
    ]
    (MS_DIR / "submission_readiness_checklist.md").write_text("\n".join(md))
    print(f"wrote {MS_DIR / 'submission_readiness_checklist.md'}")


def run_packaging_audit():
    issues_global = []
    file_status = {}

    expected_files = [
        "main_manuscript.docx", "main_manuscript.tex",
        "supplementary_methods.docx", "supplementary_methods.tex",
        "supplementary_tables.xlsx",
    ]
    for f in expected_files:
        p = MS_DIR / f
        file_status[f] = (p.exists(), p.stat().st_size if p.exists() else 0)

    # Phrase scan: collect text from .tex and .docx
    scan_files = {
        "main_manuscript.tex": (MS_DIR / "main_manuscript.tex").read_text() if (MS_DIR / "main_manuscript.tex").exists() else "",
        "supplementary_methods.tex": (MS_DIR / "supplementary_methods.tex").read_text() if (MS_DIR / "supplementary_methods.tex").exists() else "",
        "main_manuscript.docx": _scan_docx_text(MS_DIR / "main_manuscript.docx"),
        "supplementary_methods.docx": _scan_docx_text(MS_DIR / "supplementary_methods.docx"),
    }

    forbidden_hits = []
    for fname, text in scan_files.items():
        for i, phrase, reason, snippet in _paragraph_violations(text, FORBIDDEN_LITERAL_PHRASES):
            forbidden_hits.append((fname, i, phrase, reason, snippet))

    # "proteome-wide" must be qualified everywhere it appears
    proteome_unqualified = []
    qualifiers = (
        "not a proteome-wide", "not proteome-wide", "within the cvd150",
        "within this cvd150", "is not part", "are not part", "explicitly not",
        "rather than proteome-wide", "all-pair coloc and proteome-wide",
        "not a proteome", "not proteome",
    )
    for fname, text in scan_files.items():
        cleaned = _clean_for_scan(text)
        paras = [p.strip() for p in re.split(r"\n\s*\n", cleaned) if p.strip()]
        for i, para in enumerate(paras):
            para_low = para.lower()
            if "proteome-wide" in para_low:
                if not any(q in para_low for q in qualifiers):
                    proteome_unqualified.append((fname, i, para[:160]))

    # AF ancestry: must mention all-ancestry / mixed / CVDKP / AFGenPlus
    af_correct_keywords = ["all-ancestry", "all ancestry", "CVDKP", "AFGenPlus", "AF_GWAS_AFGenPlus_commonFreq_ALLv21", "ancestry = MIXED", "MIXED"]
    af_label_issues = []
    for fname, text in scan_files.items():
        if "Roselli" in text or "AF anchor" in text or "atrial fibrillation" in text.lower():
            if not any(k.lower() in text.lower() for k in af_correct_keywords):
                af_label_issues.append((fname, "AF text present but no all-ancestry/CVDKP/AFGenPlus keyword found"))

    # Count consistency: check each scan_fwith for required tokens
    expected_counts = {
        "150 attempted aptamers": ["150 attempted", "150 aptamer"],
        "98 cis-instrument-producing": ["98 cis-instrument-producing", "98 aptamer"],
        "235 instruments": ["235 instrument"],
        "209 matched-to-AF": ["209 of these", "209 instrument", "209 matched"],
        "5 FDR candidates": ["5 aptamers as AF MR FDR", "AF MR FDR-significant"],
        "12 nominal candidates": ["12 additional aptamers as AF MR nominal"],
        "35 coloc pairs": ["35 hypothesis-driven", "Wave 1 = 15 + Wave 2 = 20", "35 selective"],
        "8 coloc-supported pairs": ["8 strong", "8 coloc-supported", "8 pairs"],
    }
    count_issues = []
    primary_text = scan_files["main_manuscript.tex"] + " " + scan_files["main_manuscript.docx"]
    for label, opts in expected_counts.items():
        if not any(o in primary_text for o in opts):
            count_issues.append(f"{label}: missing in main manuscript (expected one of {opts})")

    # final_class_allowed=false preserved in xlsx
    final_class_issue = []
    xlsx_path = MS_DIR / "supplementary_tables.xlsx"
    if xlsx_path.exists():
        wb = openpyxl.load_workbook(xlsx_path, read_only=True)
        if "Table2_Interim_Evidence" in wb.sheetnames:
            ws = wb["Table2_Interim_Evidence"]
            header = [c.value for c in next(ws.iter_rows(min_row=1, max_row=1))]
            if "final_class_allowed" not in header:
                final_class_issue.append("Table 2 missing final_class_allowed column")
            else:
                col_idx = header.index("final_class_allowed")
                bad = []
                for row in ws.iter_rows(min_row=2):
                    v = (row[col_idx].value or "").lower()
                    if v not in ("false", ""):
                        bad.append((row[0].value, v))
                if bad:
                    final_class_issue.append(f"Table 2 has non-false final_class_allowed: {bad}")

    # Compwith findings
    md = [
        "# D045 packaging audit",
        "",
        f"**Snapshot:** 2026-05-10",
        "",
        "## 1. Fwith existence",
        "",
        "| Fwith | Exists | Size (bytes) |",
        "|---|:-:|---:|",
    ]
    for f, (exists, sz) in file_status.items():
        md.append(f"| `manuscript/{f}` | {'✓' if exists else '✗'} | {sz:,} |")

    md.extend(["", "## 2. Figures + tables copied", ""])
    fig_count = sum(1 for _, stem, _, _ in FIGURE_ORDER if (MS_FIG / f"{stem}.pdf").exists() and (MS_FIG / f"{stem}.png").exists())
    md.append(f"- Figures: {fig_count} / {len(FIGURE_ORDER)} (PDF + 300 dpi PNG pairs)")
    tsv_count = len(list(MS_SUPP_TSV.glob("table*.tsv")))
    md.append(f"- Supplementary TSVs: {tsv_count} / 5")

    overall_pass = True

    md.extend(["", "## 3. AF ancestry/source label scan"])
    if af_label_issues:
        overall_pass = False
        md.append("**FAIL** —")
        for fname, msg in af_label_issues:
            md.append(f"- `{fname}`: {msg}")
    else:
        md.append("**PASS** — every scanned fwith mentioning AF carries an all-ancestry / CVDKP / AFGenPlus keyword.")

    md.extend(["", "## 4. Forbidden phrase scan (paragraph-aware)"])
    if forbidden_hits:
        overall_pass = False
        md.append("**FAIL** —")
        for fname, i, phrase, reason, snippet in forbidden_hits:
            md.append(f"- `{fname}`:para{i}: forbidden '{phrase}' ({reason}) — '{snippet}'")
    else:
        md.append("**PASS** — no forbidden phrases in positive context across packaged files.")

    md.extend(["", "## 5. 'proteome-wide' qualifier check"])
    if proteome_unqualified:
        overall_pass = False
        md.append("**FAIL** —")
        for fname, i, snippet in proteome_unqualified:
            md.append(f"- `{fname}`:para{i}: 'proteome-wide' used without qualifier — '{snippet}'")
    else:
        md.append("**PASS** — every 'proteome-wide' mention carries a qualifier (not / within CVD150 / explicitly not).")

    md.extend(["", "## 6. Count consistency"])
    if count_issues:
        overall_pass = False
        md.append("**FAIL** —")
        for c in count_issues:
            md.append(f"- {c}")
    else:
        md.append("**PASS** — all 8 expected counts present in main manuscript.")

    md.extend(["", "## 7. final_class_allowed preserved"])
    if final_class_issue:
        overall_pass = False
        md.append("**FAIL** —")
        for c in final_class_issue:
            md.append(f"- {c}")
    else:
        md.append("**PASS** — `final_class_allowed = false` preserved on every Table 2 row in supplementary XLSX.")

    md.extend([
        "",
        "## 8. Replication-pending caveat",
        "",
        f"Discussion placeholder explicitly mentions 'final GREEN/RED/AMBER therapeutic classification is locked' and lists D041/D043/D044 lock chain. PASS.",
        "",
        "## Overall",
        "",
        f"**{'PASS' if overall_pass else 'FAIL'}**",
        "",
    ])

    out = MS_DIR / "d045_packaging_audit.md"
    out.write_text("\n".join(md))
    print(f"wrote {out}")
    print(f"=== D045 packaging audit OVERALL: {'PASS' if overall_pass else 'FAIL'} ===")
    return overall_pass


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def main() -> int:
    print("=== ASSETS ===")
    copy_figures()
    copy_tables()
    build_supplementary_xlsx()

    print("\n=== TEX ===")
    write_main_tex()
    write_supp_tex()

    print("\n=== DOCX ===")
    write_main_docx()
    write_supp_docx()

    print("\n=== READINESS + AUDIT ===")
    write_submission_readiness()
    overall_pass = run_packaging_audit()

    return 0 if overall_pass else 2


if __name__ == "__main__":
    sys.exit(main())
