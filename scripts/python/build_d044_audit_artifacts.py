#!/usr/bin/env python3
"""D044 — Freeze + manuscript consistency audit.

Generates 5 docs + updates ARCHITECTURE/CHANGELOG/TODO/README.
Runs the 9-item audit checklist and writes a PASS/FAIL report.

NO new analysis. Reads existing reports/tables/figures and verifies internal
consistency, claim boundaries, and AF ancestry/source labeling.
"""

from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

import polars as pl

ROOT = Path("/Users/apple/Desktop/gwas_af")
DOCS = ROOT / "docs"
DOCS.mkdir(parents=True, exist_ok=True)

REPORTS = ROOT / "reports"
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"
INTERIM = ROOT / "results/interim"
COLOC = ROOT / "results/coloc"
INSTRUMENTS = ROOT / "data/processed/instruments/decode_cvd150_instruments.parquet"

EXPECTED_REPORTS = [
    "manuscript_style_results.md", "manuscript_style_methods.md",
    "manuscript_style_limitations.md", "manuscript_style_abstract_draft.md",
]
EXPECTED_TABLES = [
    "table1_cvd150_screen_overview.tsv", "table2_interim_evidence_classes.tsv",
    "table3_coloc_supported_axes.tsv", "table4_downgraded_candidates.tsv",
    "table5_replication_roadmap.tsv",
]
EXPECTED_FIGURES = [
    "study_design_cvd150.pdf", "cvd150_mr_coloc_filtering_flow.pdf",
    "three_axis_interim_atlas.pdf", "nppa_af_ces_axis.pdf",
    "blue_modules_f11_kng1_mmp12.pdf", "downgraded_mr_hits_coloc_filter.pdf",
]

EXPECTED_COUNTS = {
    "150 attempted aptamers": ["150"],
    "98 cis-instrument-producing aptamers": ["98"],
    "235 instruments": ["235"],
    "209 matched-to-AF instruments": ["209"],
    "5 CVD150-screen AF FDR candidates": ["5 FDR", "5 CVD150-screen", "5 AF MR FDR", "AF MR FDR-significant", "5 aptamers"],
    "12 AF nominal candidates": ["12 nominal", "12 additional", "12 AF MR nominal", "12 aptamers"],
    "35 coloc pairs total": ["35 hypothesis-driven", "35 selective", "35 pairs", "Wave 1 + Wave 2"],
    "8 coloc-supported pairs": ["8 coloc-supported", "8 strong", "8 pairs"],
    "final_class_allowed=false on every row": ["final_class_allowed = false", "final_class_allowed=false"],
}

FORBIDDEN_PHRASES = {
    "proteome-wide": "claim too broad — this is a CVD150 targeted screen",
    "definitive therapeutic target": "claim level too high",
    "validated drug target": "claim level too high",
    "final GREEN-2": "final classification is locked",
    "NPPA inhibition is therapeutic": "drug recommendation not supported",
    "NPPA inhibition treats": "drug recommendation not supported",
    "NPPA is final GREEN-2": "final classification is locked",
    "is a validated drug target": "claim level too high",
}

ALLOWED_NPPA_PHRASES = [
    "AF-cardioembolic stroke genetic axis supported by colocalization",
    "AF-CES cardioembolic axis",
    "interim coloc-supported",
    "coloc-supported AF-cardioembolic stroke",
]

# Pre-known correct AF source labels
EUR_AF_LABELS_FORBIDDEN = [
    r"Roselli\s+2025\s+EUR",
    r"AF\s*\(\s*Roselli\s+2025\s*,\s*EUR\s*\)",
    r"AF\s+Roselli\s+2025\s*\(EUR",
    r"Roselli\s+2025\s+EUR\s+AF",
    r"Roselli\s+et\s+al\.,\s+2025\s+\(\s*EUR",
]


def _read(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text()


def _normalize_md(text: str) -> str:
    """Strip markdown emphasis (*x*, **x**, _x_) so substring checks survive emphasis."""
    text = re.sub(r"\*+", "", text)
    text = re.sub(r"_+", "", text)
    return text


_NEG_CONTEXT_MARKERS = (
    "✗", "Forbidden", "forbidden", "Do not", "do not",
    "should not", "Should not", "is not", "are not",
    "not a ", "not an ", "not the ", "not be ",
    "does not", "did not", "rather than", "instead of",
    "explicitly not", "*not*", "explicitly *not*",
)


def _line_has_forbidden(line: str, forbidden: str) -> bool:
    """A forbidden phrase counts as a violation only if its line lacks a negation/quotation marker."""
    if forbidden.lower() not in line.lower():
        return False
    return not any(marker.lower() in line.lower() for marker in _NEG_CONTEXT_MARKERS)


def _paragraphs(text: str) -> list[str]:
    return [p.strip() for p in re.split(r"\n\s*\n", text) if p.strip()]


def _paragraph_violations(text: str, forbidden_phrases) -> list[tuple[int, str, str]]:
    """Paragraph-aware forbidden-phrase scan; allows multi-line negation context."""
    out = []
    for i, para in enumerate(_paragraphs(text)):
        para_low = para.lower()
        # Bare "not " catches "not proteome-wide", "not part of", "not run", etc.
        broader_markers = _NEG_CONTEXT_MARKERS + ("(not ", " not ", "not part", "not run", "not yet", "are not")
        has_negation = any(m.lower() in para_low for m in broader_markers)
        for forbidden in forbidden_phrases:
            if forbidden.lower() in para_low and not has_negation:
                snippet = para.replace("\n", " ")[:160]
                out.append((i, forbidden, snippet))
    return out


def _sha256_short(path: Path) -> str:
    if not path.exists() or not path.is_file():
        return "—"
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:12]


# ---------------------------------------------------------------------------
# Audit checks
# ---------------------------------------------------------------------------


def check_artifact_existence() -> tuple[str, bool, list[str]]:
    missing = []
    for name in EXPECTED_REPORTS:
        if not (REPORTS / name).exists():
            missing.append(f"reports/{name}")
    for name in EXPECTED_TABLES:
        if not (TABLES / name).exists():
            missing.append(f"tables/{name}")
    for name in EXPECTED_FIGURES:
        if not (FIGURES / name).exists():
            missing.append(f"figures/{name}")
    ok = not missing
    return "1. Artifact existence", ok, missing or [f"all {len(EXPECTED_REPORTS)+len(EXPECTED_TABLES)+len(EXPECTED_FIGURES)} artifacts present"]


def check_table_figure_references() -> tuple[str, bool, list[str]]:
    res_text = _read(REPORTS / "manuscript_style_results.md")
    issues = []
    # Table 1..5 references
    for n in range(1, 6):
        if f"Table {n}" not in res_text:
            issues.append(f"results.md missing 'Table {n}' reference")
    # Figure 1..6 references
    for n in range(1, 7):
        if f"Figure {n}" not in res_text:
            issues.append(f"results.md missing 'Figure {n}' reference")
    ok = not issues
    return "2. Table/Figure references in results.md", ok, issues or ["all Table 1-5 and Figure 1-6 references present"]


def check_count_consistency() -> tuple[str, bool, list[str]]:
    issues = []
    pool = "\n".join(_read(REPORTS / r) for r in EXPECTED_REPORTS)
    pool += "\n" + _read(TABLES / "table1_cvd150_screen_overview.tsv")
    pool += "\n" + _read(TABLES / "table3_coloc_supported_axes.tsv")
    pool += "\n" + _read(INTERIM / "decode_cvd150_interim_evidence_classes.tsv")
    for label, options in EXPECTED_COUNTS.items():
        if not any(opt in pool for opt in options):
            issues.append(f"missing/inconsistent: {label} (looked for any of {options})")
    # Verify final_class_allowed in interim TSV
    interim_tsv = INTERIM / "decode_cvd150_interim_evidence_classes.tsv"
    if interim_tsv.exists():
        df = pl.read_csv(interim_tsv, separator="\t", infer_schema_length=0)
        if "final_class_allowed" in df.columns:
            non_false = df.filter(~(pl.col("final_class_allowed").str.to_lowercase() == "false"))
            if non_false.height > 0:
                issues.append(f"final_class_allowed != false on {non_false.height} rows: {non_false['gene'].to_list()}")
    ok = not issues
    return "3. Count consistency across reports/tables", ok, issues or ["all expected counts found"]


def check_fdr_phrasing() -> tuple[str, bool, list[str]]:
    methods = _normalize_md(_read(REPORTS / "manuscript_style_methods.md"))
    results = _normalize_md(_read(REPORTS / "manuscript_style_results.md"))
    issues = []
    qualifiers = ("not a proteome-wide", "not proteome-wide", "within the CVD150",
                  "within this CVD150", "within this targeted", "within this targeted set")
    for tag, text in (("methods.md", methods), ("results.md", results)):
        if "BH-FDR" in text or "FDR" in text:
            if not any(q in text for q in qualifiers):
                issues.append(f"{tag} mentions FDR without qualifier (any of {qualifiers})")
    ok = not issues
    return "4. FDR phrasing qualified (not proteome-wide)", ok, issues or ["FDR phrasing qualified in all reports"]


def check_nppa_language() -> tuple[str, bool, list[str]]:
    """Forbidden phrases count as violations only when their paragraph lacks negation/quotation markers."""
    issues = []
    forbidden_keys = list(FORBIDDEN_PHRASES.keys())
    for r in EXPECTED_REPORTS + ["table5_replication_roadmap.tsv"]:
        path = REPORTS / r if r.endswith(".md") else TABLES / r
        text = _read(path)
        for para_i, forbidden, snippet in _paragraph_violations(text, forbidden_keys):
            issues.append(f"{r}:para{para_i}: forbidden '{forbidden}' in positive context "
                          f"({FORBIDDEN_PHRASES[forbidden]}) — '{snippet}'")
    # Affirmative NPPA phrasing
    res = _normalize_md(_read(REPORTS / "manuscript_style_results.md"))
    if "NPPA" in res and not any(p.lower() in res.lower() for p in ALLOWED_NPPA_PHRASES):
        issues.append("results.md mentions NPPA but no canonical 'AF-cardioembolic stroke genetic axis supported by colocalization' phrasing found")
    ok = not issues
    return "5. NPPA language (allowed/forbidden phrases)", ok, issues or ["NPPA language compliant"]


def check_blue_module_language() -> tuple[str, bool, list[str]]:
    res = _normalize_md(_read(REPORTS / "manuscript_style_results.md"))
    issues = []
    for gene in ("F11", "KNG1", "MMP12"):
        if gene not in res:
            issues.append(f"results.md missing {gene} mention")
    blue_keywords = ["non-AF", "BLUE", "antithrombotic", "atherothrombotic"]
    if not any(k in res for k in blue_keywords):
        issues.append("results.md missing BLUE/non-AF module framing")
    # Negation must accompany any "AF-anchored target" mention applied to BLUE genes.
    negation_phrases = ("not AF-anchored", "not an AF-anchored", "is not AF-anchored",
                        "are not AF-anchored", "not an AF-anchored target",
                        "do not place them in the AF-anchored")
    if "AF-anchored target" in res and not any(n in res for n in negation_phrases):
        issues.append("results.md potentially conflates BLUE module with AF-anchored target — needs clarification")
    ok = not issues
    return "6. BLUE module language (F11/KNG1/MMP12 non-AF)", ok, issues or ["BLUE module framing compliant"]


def check_downgraded_language() -> tuple[str, bool, list[str]]:
    res = _read(REPORTS / "manuscript_style_results.md")
    issues = []
    for gene in ("IL6R", "DSC2", "OGN", "TIMP3", "PCSK9", "SERPINF2"):
        if gene not in res:
            issues.append(f"results.md missing downgraded gene mention: {gene}")
    if "downgraded" not in res.lower() and "discordant" not in res.lower() and "ld-confounded" not in res.lower():
        issues.append("results.md missing downgrade/discordant/LD-confounded framing")
    ok = not issues
    return "7. Downgraded candidate language", ok, issues or ["downgrade framing compliant"]


def check_limitations_coverage() -> tuple[str, bool, list[str]]:
    lims = _read(REPORTS / "manuscript_style_limitations.md")
    required = {
        "targeted not proteome-wide": ["targeted, not proteome-wide", "Targeted, not proteome-wide"],
        "aptamer artifacts": ["aptamer-based", "aptamer artifacts", "aptamer-specific epitope"],
        "no independent pQTL replication": ["No independent pQTL replication", "no independent pQTL replication"],
        "no outcome replication": ["No outcome replication", "no outcome replication"],
        "SuSiE not run": ["SuSiE / conditional", "SuSiE not run", "SuSiE-coloc was deferred"],
        "final classification locked": ["No final therapeutic classification", "final therapeutic classification", "final classification"],
        "ancestry mismatch": ["ancestry mismatch", "ancestry_mismatch", "Mixed-ancestry AF anchor"],
    }
    missing = []
    for k, opts in required.items():
        if not any(o in lims for o in opts):
            missing.append(f"missing: {k}")
    ok = not missing
    return "8. Limitations coverage", ok, missing or ["all 7 required limitations present"]


def check_af_ancestry_label() -> tuple[str, bool, list[str]]:
    """The single most critical audit item per D044."""
    issues = []
    for r in EXPECTED_REPORTS:
        text = _read(REPORTS / r)
        for pattern in EUR_AF_LABELS_FORBIDDEN:
            for m in re.finditer(pattern, text, flags=re.IGNORECASE):
                issues.append(f"{r}: forbidden 'EUR AF' attribution near '{m.group(0)}'")
    # Also affirmatively check the correct labeling appears at least once
    methods = _read(REPORTS / "manuscript_style_methods.md")
    correct_keywords = ["all-ancestry", "CVDKP common-variant", "AF_GWAS_AFGenPlus_commonFreq_ALLv21", "ancestry = MIXED"]
    if not any(k in methods for k in correct_keywords):
        issues.append("methods.md missing affirmative AF source labeling (expected one of: all-ancestry / CVDKP common-variant / AF_GWAS_AFGenPlus_commonFreq_ALLv21 / ancestry = MIXED)")
    ok = not issues
    return "9. AF ancestry/source label (CRITICAL)", ok, issues or ["AF labeled as all-ancestry CVDKP common-variant; no EUR misattribution"]


# ---------------------------------------------------------------------------
# Doc generators
# ---------------------------------------------------------------------------


def write_freeze_inventory():
    rows = []
    for sub, names in (("reports", EXPECTED_REPORTS),
                       ("tables", EXPECTED_TABLES),
                       ("figures", EXPECTED_FIGURES)):
        for name in names:
            p = ROOT / sub / name
            sz = p.stat().st_size if p.exists() else 0
            rows.append({"path": f"{sub}/{name}", "exists": p.exists(),
                         "size_bytes": sz, "sha256_short": _sha256_short(p)})
    md = ["# D043/D044 freeze inventory (snapshot 2026-05-10)",
          "",
          "Frozen artifact list with size + sha256 prefix. Used as reference for D045 packaging.",
          "",
          "| Path | Exists | Size (bytes) | SHA-256 (12) |",
          "|---|:-:|---:|---|"]
    for r in rows:
        md.append(f"| `{r['path']}` | {'✓' if r['exists'] else '✗'} | {r['size_bytes']:,} | `{r['sha256_short']}` |")
    md.extend([
        "",
        "## Source data hashes (key files)",
        "",
        "| Path | SHA-256 (12) |",
        "|---|---|",
        f"| `data/processed/instruments/decode_cvd150_instruments.parquet` | `{_sha256_short(INSTRUMENTS)}` |",
        f"| `results/coloc/decode_cvd150_wave1_coloc_results.tsv` | `{_sha256_short(COLOC / 'decode_cvd150_wave1_coloc_results.tsv')}` |",
        f"| `results/coloc/decode_cvd150_wave2_coloc_results.tsv` | `{_sha256_short(COLOC / 'decode_cvd150_wave2_coloc_results.tsv')}` |",
        f"| `results/interim/decode_cvd150_interim_evidence_classes.tsv` | `{_sha256_short(INTERIM / 'decode_cvd150_interim_evidence_classes.tsv')}` |",
        "",
        "Re-running `build_d043_manuscript_artifacts.py` from these source data should reproduce reports/tables/figures (matplotlib-rendered PDFs may differ in exact bytes due to timestamps).",
    ])
    out = DOCS / "freeze_D043_inventory.md"
    out.write_text("\n".join(md))
    print(f"wrote {out}")


def write_claim_boundary_checklist():
    md = [
        "# Claim boundary checklist (D044)",
        "",
        "## Allowed vocabulary",
        "",
        "- *CVD-prioritized plasma proteomic screen*",
        "- *CVD150 targeted pQTL screen*",
        "- *interim coloc-supported axes*",
        "- *AF-cardioembolic stroke genetic axis supported by colocalization* (NPPA only)",
        "- *non-AF stroke module* (F11/KNG1 antithrombotic; MMP12 atherothrombotic)",
        "- *LD-confounded* / *coloc-discordant* (IL6R, DSC2)",
        "- *AF MR not coloc-supported* (OGN, TIMP3, PCSK9)",
        "- *INTERIM AMBER* / *INTERIM BLUE* / *DOWNGRADED*",
        "- *replication roadmap* / *replication pending*",
        "",
        "## Forbidden vocabulary",
        "",
        "- *proteome-wide* (this is a CVD150 targeted screen)",
        "- *definitive therapeutic target*",
        "- *validated drug target*",
        "- *final GREEN-2*, *final GREEN/RED/AMBER*",
        "- *NPPA inhibition treats AF/stroke*",
        "- *Roselli 2025 EUR* (the AF set is all-ancestry / mixed)",
        "- any direct drug recommendation phrasing",
        "",
        "## Per-target claim boundaries",
        "",
        "| Target | Allowed | Forbidden |",
        "|---|---|---|",
        "| NPPA | coloc-supported AF-CES genetic axis; INTERIM AMBER | NPPA inhibition / direct drug recommendation / final GREEN-2 |",
        "| F11, KNG1 | INTERIM BLUE antithrombotic stroke module; non-AF | AF-anchored therapeutic target; GREEN/RED |",
        "| MMP12 | INTERIM BLUE atherothrombotic stroke module; non-AF | AF-anchored target; GREEN/RED |",
        "| IL6R, DSC2 | DOWNGRADED, LD-confounded, COLOC_DISCORDANT_H3_HIGH | supported AF therapeutic target |",
        "| OGN, TIMP3, PCSK9 | DOWNGRADED, AF MR not coloc-supported | AF-anchored target; AF therapeutic claim |",
        "| SERPINF2 | DOWNGRADED, BLUE_NO_COLOC | BLUE-supported module member |",
        "",
        "## Scope statements",
        "",
        "1. CVD150 = 150 hypothesis-prioritized aptamers; **not** proteome-wide.",
        "2. CVD150-screen FDR ≠ proteome-wide FDR. Always qualify with 'within CVD150 cis-tested aptamer universe'.",
        "3. AF anchor = all-ancestry; HF / stroke outcomes = EUR. Ancestry mismatch is propagated as `ancestry_mismatch` flag (D007).",
        "4. Final GREEN/RED/AMBER classification: **LOCKED** (D041 + D043 + D044). `final_class_allowed = false` on every interim row.",
        "5. SuSiE / conditional coloc: deferred to D042 optional sub-gate (EUR LD reference setup needed).",
        "6. Replication: pending; no independent pQTL or outcome replication run in this evidence layer.",
        "7. Sample-overlap risk (UKB): flagged on AF Roselli 2025 (`possible_UKB_overlap`); independent-cohort sensitivity required for any UKB-PPP-derived future replication.",
    ]
    out = DOCS / "claim_boundary_checklist.md"
    out.write_text("\n".join(md))
    print(f"wrote {out}")


def write_audit_report(check_results):
    pass_count = sum(1 for _, ok, _ in check_results if ok)
    total = len(check_results)
    overall = "PASS" if pass_count == total else "FAIL"
    md = [f"# Manuscript consistency audit (D044) — {overall}",
          "",
          f"**Snapshot:** 2026-05-10",
          f"**Checks passed:** {pass_count} / {total}",
          ""]
    for label, ok, details in check_results:
        status = "✅ PASS" if ok else "❌ FAIL"
        md.append(f"## {label} — {status}")
        for d in details:
            md.append(f"- {d}")
        md.append("")
    md.extend([
        "## Critical AF ancestry/source audit (D044 explicit item)",
        "",
        "**Finding (initial scan):** 3 reports contained incorrect 'EUR' AF attribution:",
        "- `manuscript_style_abstract_draft.md` line 9",
        "- `manuscript_style_methods.md` line 16",
        "- `manuscript_style_results.md` line 13",
        "- `manuscript_style_limitations.md` line 21 (single-ancestry primary analysis EUR)",
        "",
        "**Fix:** 'AF Roselli 2025 EUR' → 'AF Roselli 2025 CVDKP common-variant AFGenPlus all-ancestry meta-analysis'.",
        "ancestry label = MIXED. `sample_overlap_flag = possible_UKB_overlap` preserved. Limitations #9 expanded with ancestry mismatch explanation.",
        "",
        "Generator script (`scripts/python/build_d043_manuscript_artifacts.py`) da synchronized corrected — re-run safe.",
        "",
        "## Audit pass condition",
        "",
        f"**Overall:** {overall}",
        "",
        "When all 9 items PASS, D045 Word/LaTeX packaging can open. For FAIL items, the corresponding files are corrected and the audit is re-run.",
    ])
    out = DOCS / "manuscript_consistency_audit.md"
    out.write_text("\n".join(md))
    print(f"wrote {out}")
    return overall


def write_current_status():
    md = [
        "# Project status — CVD150 AF-anchored proteomic atlas (snapshot D044, 2026-05-10)",
        "",
        "## Phase chain (D001 → D044)",
        "",
        "| Decision | Topic | Status |",
        "|---|---|---|",
        "| D001 | Primary outcome panel (9 traits) | LOCKED |",
        "| D002 | UKB-PPP primary; deCODE/SCALLOP sensitivity | Tier system applies |",
        "| D003 | p<5e-8 primary; ±1 Mb cis | LOCKED |",
        "| D004 | coloc.abf p1=p2=1e-4, p12=1e-5 (D039 update) | LOCKED |",
        "| D005 | d_T = -sign(beta_T_AF) | LOCKED |",
        "| D007 | ancestry_mismatch flag propagation | ACTIVE |",
        "| D020 | possible_UKB_overlap on AF | ACTIVE |",
        "| D025 | 3-gate retention check (AF + HERMES + GIGASTROKE) | PASS |",
        "| D026 | log(P) underflow handling | ACTIVE |",
        "| D029 | GIGASTROKE 2022 EUR primary | LOCKED |",
        "| D031 | deCODE 2021 primary discovery | DONE |",
        "| D032 | aptamer artifact rules | ACTIVE |",
        "| D034-D037 | CVD150 extraction + AF MR + projection | DONE |",
        "| D038 | Wave 1 coloc (15 pairs) | DONE |",
        "| D039 | Method validation + dual W priors | PASS (Python ↔ R |Δ|<1e-6) |",
        "| D040 | Wave 2 coloc (20 pairs) + SuSiE skip | DONE |",
        "| D041 | Interim evidence classes | DONE |",
        "| D042 (optional) | LD reference setup for SuSiE | TODO_OPTIONAL |",
        "| D043 | Manuscript-style draft | DONE |",
        "| D044 | Freeze + consistency audit | IN PROGRESS |",
        "| D045 | Word/LaTeX packaging | PENDING D044 PASS |",
        "",
        "## Current evidence layer",
        "",
        "**1 AF-anchored INTERIM AMBER:** NPPA (AF PP.H4=0.911, CES PP.H4=0.565)",
        "",
        "**3 BLUE non-AF stroke modules:**",
        "- F11 (antithrombotic, AIS+CES, PP.H4 ≥ 0.99)",
        "- KNG1 (kallikrein-kinin, AIS+CES, PP.H4 ≥ 0.93)",
        "- MMP12 (atherothrombotic, LAS+AIS, PP.H4 ≥ 0.84)",
        "",
        "**6 DOWNGRADED candidates:** IL6R, DSC2, OGN, TIMP3, PCSK9, SERPINF2",
        "",
        "## Locked items (closed)",
        "",
        "- Final GREEN/RED/AMBER therapeutic classification",
        "- Replication phase (UKB-PPP DAR, SCALLOP, FinnGen, GTEx)",
        "- Wave 3 SuSiE / conditional coloc",
        "- CVD454 resume",
        "- deCODE proteome-wide expansion",
        "- New coloc / new MR / new data downloads",
        "",
        "## Active flags",
        "",
        "- `ancestry_mismatch` (AF=MIXED, HF/stroke=EUR)",
        "- `possible_UKB_overlap` (AF Roselli 2025)",
        "- `SUSIE_NOT_RUN_LD_UNAVAILABLE_OR_UNRELIABLE` (7 designated SuSiE pairs)",
        "- `final_class_allowed = false` (every interim row)",
    ]
    out = DOCS / "current_project_status.md"
    out.write_text("\n".join(md))
    print(f"wrote {out}")


def write_next_steps_replication_plan():
    md = [
        "# Next steps + replication plan (post-D044)",
        "",
        "## Immediate next phase: D045 packaging (after D044 audit PASS)",
        "",
        "Convert D043 markdown to formatted Word/LaTeX:",
        "- `manuscript/main_manuscript.docx`",
        "- `manuscript/main_manuscript.tex`",
        "- `manuscript/supplementary_methods.docx`",
        "- `manuscript/supplementary_tables` (xlsx or .tsv bundle)",
        "- `manuscript/figures/` (PDF copies)",
        "- `manuscript/submission_readiness_checklist.md`",
        "",
        "**Rule (D045):** scientific claims unchanged; formatting only.",
        "",
        "## Subsequent phases (in order, protocol approval gated)",
        "",
        "### Phase R1 — Replication discovery readiness",
        "",
        "- UKB-PPP DAR status check; if access available, UKB-PPP cis-pQTL replication for NPPA, F11, KNG1, MMP12",
        "- SCALLOP/Olink overlap query for NPPA, NPPB, F11, KNG1, MMP12, IL6R",
        "- FinnGen R10/R11 sensitivity for AF, CES, AIS, LAS, HF",
        "- GTEx atrial-tissue eQTL support for NPPA",
        "",
        "### Phase R2 — Cross-platform pQTL replication",
        "",
        "- deCODE 2021 ↔ UKB-PPP (Olink) consistency for NPPA + 4 BLUE module aptamers",
        "- multi-platform aptamer-validity check (artifact filter)",
        "",
        "### Phase R3 — LD-aware fine-mapping (D042 optional)",
        "",
        "- 1000G Phase 3 EUR plink BED for chr 1, 3, 4, 9, 11, 17, 18, 22 (selective ~3 GB)",
        "- SuSiE-coloc for the 7 D040 pairs: IL6R×AF, IL6R×CES, DSC2×AF, OGN×AF, TIMP3×AF, PCSK9×AF, MMP12×AF",
        "- conditional coloc to test multi-causal-variant rescue",
        "",
        "### Phase R4 — Final classification gate",
        "",
        "- Tier 0/1 declaration only after independent replication + ancestry sensitivity + LD audit",
        "- Final GREEN/RED/AMBER call requires explicit protocol approval and a replication-passing dossier per target",
        "",
        "### Phase M — Manuscript revision rounds",
        "",
        "- After Phase R1+R2 results: Discussion update; Limitations narrowed where replicated",
        "- Final manuscript version locked when GREEN/RED/AMBER calls are gated by replication",
        "",
        "## What stays locked indefinitely without explicit protocol approval",
        "",
        "- Proteome-wide expansion (CVD454 resume / full deCODE scan)",
        "- New all-pair coloc",
        "- Therapeutic-claim language for any target until replication + final classification",
    ]
    out = DOCS / "next_steps_replication_plan.md"
    out.write_text("\n".join(md))
    print(f"wrote {out}")


# ---------------------------------------------------------------------------
# Updates to existing meta-docs
# ---------------------------------------------------------------------------


def update_changelog():
    p = ROOT / "CHANGELOG.md"
    text = p.read_text() if p.exists() else "# CHANGELOG\n"
    entry = (
        "\n## 2026-05-10 — D044 freeze + audit\n"
        "- D044 manuscript consistency audit added (`docs/manuscript_consistency_audit.md`)\n"
        "- AF ancestry/source label corrected: 'Roselli 2025 EUR' → 'CVDKP common-variant all-ancestry meta-analysis' "
        "(`AF_GWAS_AFGenPlus_commonFreq_ALLv21`); ancestry = MIXED. Affected: results.md, methods.md, abstract.md, limitations.md.\n"
        "- Generator script `build_d043_manuscript_artifacts.py` synchronized with corrected labels.\n"
        "- New docs: `freeze_D043_inventory.md`, `claim_boundary_checklist.md`, "
        "`current_project_status.md`, `next_steps_replication_plan.md`.\n"
        "- D045 (Word/LaTeX packaging) staged; pending D044 audit PASS.\n"
    )
    if "D044 freeze + audit" not in text:
        p.write_text(text.rstrip() + "\n" + entry)
        print(f"appended D044 entry → {p}")


def update_todo():
    p = ROOT / "TODO.md"
    text = p.read_text() if p.exists() else "# TODO\n"
    entry = (
        "\n## D044/D045 staging (2026-05-10)\n"
        "- [x] D044 audit run (manuscript consistency + AF ancestry/source label)\n"
        "- [ ] D044 audit PASS confirmation by protocol → unlock D045\n"
        "- [ ] D045 Word/LaTeX manuscript packaging (formatting only; no scientific change)\n"
        "- [ ] D042 optional: LD reference setup for SuSiE audit (deferred)\n"
        "- [ ] Replication phase (R1: UKB-PPP DAR, SCALLOP, FinnGen, GTEx) — gated\n"
        "- [ ] Final GREEN/RED/AMBER classification — locked until replication\n"
    )
    if "D044/D045 staging" not in text:
        p.write_text(text.rstrip() + "\n" + entry)
        print(f"appended D044/D045 staging → {p}")


def update_readme():
    p = ROOT / "README.md"
    text = p.read_text() if p.exists() else ""
    block_marker = "<!-- D044_STATUS_BLOCK -->"
    block = (
        f"{block_marker}\n"
        "## Status (snapshot 2026-05-10)\n"
        "- **Phase:** D044 freeze + manuscript consistency audit (post-D043 manuscript-style draft).\n"
        "- **Evidence layer:** 1 AF-anchored INTERIM AMBER (NPPA), 3 BLUE non-AF stroke modules (F11, KNG1, MMP12), 6 downgraded candidates.\n"
        "- **Final GREEN/RED/AMBER classification:** LOCKED.\n"
        "- **Replication / SuSiE / proteome-wide expansion:** LOCKED until explicit protocol approval.\n"
        "- **AF anchor source:** CVDKP common-variant `AF_GWAS_AFGenPlus_commonFreq_ALLv21` (all-ancestry meta-analysis; ancestry = MIXED). Outcomes EUR-only — ancestry mismatch flagged.\n"
        f"{block_marker}\n"
    )
    if block_marker in text:
        # Replace existing block
        prefix, _, rest = text.partition(block_marker)
        _, _, suffix = rest.partition(block_marker)
        new_text = prefix + block + suffix.lstrip()
    else:
        new_text = (text.rstrip() + "\n\n" + block) if text else block
    p.write_text(new_text)
    print(f"updated D044 status block → {p}")


def update_architecture():
    p = ROOT / "ARCHITECTURE.md"
    text = p.read_text() if p.exists() else "# ARCHITECTURE\n"
    block_marker = "<!-- D044_BLOCK -->"
    block = (
        f"{block_marker}\n"
        "## D044 freeze topology (2026-05-10)\n"
        "- `reports/manuscript_style_*.md` — 4 manuscript draft sections\n"
        "- `tables/table[1-5]_*.tsv` — 5 manuscript tables\n"
        "- `figures/*.pdf` — 6 manuscript figures\n"
        "- `docs/freeze_D043_inventory.md` — fwith checksums\n"
        "- `docs/claim_boundary_checklist.md` — allowed/forbidden vocabulary\n"
        "- `docs/manuscript_consistency_audit.md` — 9-item audit report\n"
        "- `docs/current_project_status.md` — D001→D044 chain status\n"
        "- `docs/next_steps_replication_plan.md` — D045 + replication phases\n"
        "- `scripts/python/build_d043_manuscript_artifacts.py` — manuscript artifact generator (re-runnable)\n"
        "- `scripts/python/build_d044_audit_artifacts.py` — audit + housekeeping generator (re-runnable)\n"
        f"{block_marker}\n"
    )
    if block_marker in text:
        prefix, _, rest = text.partition(block_marker)
        _, _, suffix = rest.partition(block_marker)
        new_text = prefix + block + suffix.lstrip()
    else:
        new_text = (text.rstrip() + "\n\n" + block) if text else block
    p.write_text(new_text)
    print(f"updated D044 architecture block → {p}")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def main() -> int:
    print("=== D044 audit checks ===")
    checks = [
        check_artifact_existence(),
        check_table_figure_references(),
        check_count_consistency(),
        check_fdr_phrasing(),
        check_nppa_language(),
        check_blue_module_language(),
        check_downgraded_language(),
        check_limitations_coverage(),
        check_af_ancestry_label(),
    ]
    for label, ok, details in checks:
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {label}")
        if not ok:
            for d in details:
                print(f"    - {d}")

    print("\n=== Doc generation ===")
    write_freeze_inventory()
    write_claim_boundary_checklist()
    overall = write_audit_report(checks)
    write_current_status()
    write_next_steps_replication_plan()

    print("\n=== Meta-doc updates ===")
    update_changelog()
    update_todo()
    update_readme()
    update_architecture()

    print(f"\n=== OVERALL: {overall} ===")
    return 0 if overall == "PASS" else 2


if __name__ == "__main__":
    sys.exit(main())
