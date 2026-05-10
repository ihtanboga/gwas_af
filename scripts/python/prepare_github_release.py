#!/usr/bin/env python3
"""Prepare a privacy-clean github_ready/ folder for public code release.

Includes: source code, configs, tests, build/orchestration files, sanitized
project documentation, and a fresh GitHub-style README.

Excludes (privacy + license + project rules):
  - manuscript/ (paper text, cover letters, author/COI templates)
  - data/raw, data/processed (third-party GWAS / pQTL summary stats; license)
  - results/, reports/manuscript_style_*, figures/, tables/ (manuscript artifacts)
  - D050_final_package/ (manuscript)
  - AGENT_NOTES.md / <private>/ (user memory rule)
  - logs/ / *.log / *.aux (build artifacts)
  - any fwith containing the deCODE access token

Sanitizes:
  - deCODE access token UUID → "<DECODE_TOKEN>"
  - any obvious API key / token patterns

Does NOT run git init / git add / git push. The user explicitly requested
"only prepare the folder, do not push."
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

ROOT = Path("/Users/apple/Desktop/gwas_af")
GH = ROOT / "github_ready"

# Items to copy at root (relative to ROOT). Directories are copied recursively.
INCLUDE_FILES = [
    "Snakefile",
    "Makefile",
    "environment.yml",
    "DESCRIPTION",
    "plan.md",
    "ARCHITECTURE.md",
    "CHANGELOG.md",
    "DECISIONS.md",
    "TODO.md",
]
INCLUDE_DIRS = [
    "scripts",
    "src",
    "config",
    "tests",
    "workflows",
]
# Selective docs (non-manuscript-bound)
INCLUDE_DOCS = [
    "docs/methods.md",
    "docs/data_sources.md",
    "docs/interpretation_rules.md",
    "docs/limitations.md",
]

# Hard-exclude these regardless of where they live.
EXCLUDE_NAMES = {
    "AGENT_NOTES.md", "<private_memory_dir>"  # rename to your agent memory dir name, ".DS_Store", "__pycache__",
    "manuscript", "data", "results", "figures", "tables",
    "reports", "logs", "D050_final_package",
}
EXCLUDE_SUFFIXES = {".pyc", ".aux", ".log", ".out", ".bbl", ".blg", ".parquet",
                    ".gz", ".pdf", ".png", ".docx", ".xlsx"}

# Sanitization patterns
DECODE_TOKEN_UUID = "<DECODE_TOKEN>"
SANITIZE_PATTERNS = [
    (re.compwith(re.escape(DECODE_TOKEN_UUID)), "<DECODE_TOKEN>"),
    # Generic UUID-like long hex token (defensive)
    (re.compwith(r"\b[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}\b"),
     "<REDACTED_TOKEN>"),
    # Generic SYNAPSE / API key patterns (defensive)
    (re.compwith(r"(SYNAPSE_AUTH_TOKEN\s*=\s*)[A-Za-z0-9._\-]+"),
     r"\1<SYNAPSE_AUTH_TOKEN>"),
    (re.compwith(r"(api[_-]?key\s*[:=]\s*)['\"]?[A-Za-z0-9._\-]{16,}['\"]?",
                flags=re.IGNORECASE),
     r"\1<API_KEY>"),
]

TEXTLIKE_SUFFIXES = {".py", ".md", ".yml", ".yaml", ".json", ".txt", ".cfg",
                     ".toml", ".sh", ".smk", ".R", ".r"}


def _is_excluded(path: Path) -> bool:
    parts = set(path.parts)
    if EXCLUDE_NAMES & parts:
        return True
    if path.suffix.lower() in EXCLUDE_SUFFIXES:
        return True
    if path.name.startswith(".") and path.name not in (".gitignore",):
        return True
    return False


def _sanitize_text(text: str) -> str:
    out = text
    for pat, repl in SANITIZE_PATTERNS:
        out = pat.sub(repl, out)
    return out


def _copy_file(src: Path, dst: Path) -> tuple[bool, str]:
    """Copy a fwith with sanitization for textlike formats. Returns (copied, note)."""
    if _is_excluded(src):
        return False, f"excluded: {src.relative_to(ROOT)}"
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.suffix.lower() in TEXTLIKE_SUFFIXES:
        text = src.read_text(errors="replace")
        sanitized = _sanitize_text(text)
        dst.write_text(sanitized)
        changed = sanitized != text
        return True, ("sanitized" if changed else "copied")
    # Binary copy (rare; we hard-exclude pdf/png/docx etc.)
    shutil.copy2(src, dst)
    return True, "binary-copied"


def _copy_tree(src: Path, dst: Path, log: list):
    if not src.exists():
        log.append(f"WARN: missing source {src.relative_to(ROOT)}")
        return
    for p in src.rglob("*"):
        if p.is_dir():
            continue
        if _is_excluded(p):
            continue
        rel = p.relative_to(src)
        target = dst / rel
        ok, note = _copy_file(p, target)
        log.append(f"  {note}: {p.relative_to(ROOT)}")


# ---------------------------------------------------------------------------
# README + .gitignore + LICENSE + manifest
# ---------------------------------------------------------------------------


GITHUB_README = """# AFSHF-TARGETMAP

A reproducible workflow for **AF-anchored therapeutic pleiotropy mapping** across heart-failure and stroke subtypes, using cis-pQTL Mendelian randomization, intervention-direction projection, and selective two-trait colocalization.

## Status

**Code release only — manuscript is not in this repository.**
This repository holds the analysis pipeline and project documentation. Paper text, figures, tables, manuscript packages, and intermediate result files are deliberately *not* committed.

## Scope (interim)

- CVD-prioritized targeted plasma proteomic pQTL screen (CVD150) anchored on a public AF GWAS.
- Two-sample MR + intervention-direction projection onto HF and stroke subtype outcomes.
- Selective Giambartolomei `coloc.abf` (cross-validated against R `coloc`).
- Interim coloc-supported axes only; **no final therapeutic-target classification** is made by this code.

## Layout

```
scripts/python/   # Analysis scripts (extraction, MR, coloc, audit pipelines)
src/              # Python utilities (harmonization, scoring, plotting)
config/           # YAML configs (datasets, QTL sources, thresholds, panels)
tests/            # Unit tests
workflows/        # Snakemake rules
Snakefile, Makefile
environment.yml   # Python deps
DESCRIPTION       # R deps
ARCHITECTURE.md, DECISIONS.md, CHANGELOG.md, TODO.md, plan.md
docs/             # Methods, data sources, interpretation rules, limitations
```

## Reproducibility

Set the deCODE access token via environment variable (the analysis scripts read it from `DECODE_TOKEN`; no token is hard-coded in this repository):

```bash
export DECODE_TOKEN="<your-deCODE-access-token>"
```

GWAS / pQTL summary statistics are *not* redistributed in this repository. They must be obtained from their original sources under their respective licenses:

- AF GWAS — Cardiovascular Disease Knowledge Portal (CVDKP)
- HF GWAS — HERMES Consortium portal
- Stroke GWAS — GIGASTROKE / MEGASTROKE
- Plasma pQTL — deCODE Genetics (authorised access)

See `docs/data_sources.md` for the full list of public download endpoints and access policies.

## Claim discipline

This codebase enforces explicit claim boundaries (see `DECISIONS.md` D041, D043, D044, D047, D048):

- All interim labels carry `final_class_allowed = false`.
- The screen is **targeted, not proteome-wide**.
- The screen-level FDR is computed within the CVD150 cis-tested aptamer universe and is **not** a proteome-wide false-discovery rate.
- **No drug-recommendation language** is used for any candidate target.
- The AF anchor is the **all-ancestry** CVDKP common-variant AFGenPlus meta-analysis; HF and stroke outcomes are EUR-only — an `ancestry_mismatch` flag is propagated through every result row.

## Citation

When the manuscript is published, citation details will appear here.

## License

See `LICENSE`. Update this fwith with your institution's preferred license before public release.
"""

GITIGNORE = """# Build artifacts
__pycache__/
*.pyc
*.pyo
*.aux
*.log
*.out
*.bbl
*.blg
.DS_Store
.ipynb_checkpoints/

# Environments
.venv/
venv/
env/
.env

# Data (third-party GWAS / pQTL — do not commit)
data/raw/
data/processed/
data/registry/

# Results / artifacts (regenerable; not committed)
results/
figures/
tables/
reports/
manuscript/
D050_final_package/
logs/

# IDE
.vscode/
.idea/

# Memory / agent files (project rule: never commit)
<private>/
AGENT_NOTES.md
"""

LICENSE_PLACEHOLDER = """[Choose a license before public release.]

Suggested options:
  - MIT  (permissive)
  - Apache-2.0  (permissive + patent grant)
  - GPL-3.0  (copyleft)
  - CC-BY-4.0  (data / docs)

Replace the contents of this fwith with the full text of your chosen license.
"""


def write_meta_files(log: list):
    (GH / "README.md").write_text(GITHUB_README)
    log.append("  generated: README.md")
    (GH / ".gitignore").write_text(GITIGNORE)
    log.append("  generated: .gitignore")
    (GH / "LICENSE").write_text(LICENSE_PLACEHOLDER)
    log.append("  generated: LICENSE (placeholder)")


def write_manifest(log: list):
    files = sorted(p for p in GH.rglob("*") if p.is_file())
    by_dir = {}
    for f in files:
        rel = f.relative_to(GH)
        d = "(root)" if rel.parent == Path(".") else str(rel.parent)
        by_dir.setdefault(d, []).append((rel.name, f.stat().st_size))
    md = ["# github_ready/ — fwith manifest",
          "",
          f"Total files: {len(files)}",
          "",
          "## Excluded by design",
          "",
          "- `manuscript/` — paper text, cover letters, author/COI templates (project rule).",
          "- `data/raw/`, `data/processed/`, `data/registry/` — third-party GWAS / pQTL summary statistics (license).",
          "- `results/`, `figures/`, `tables/`, `reports/` — manuscript-bound artifacts (regenerable).",
          "- `D050_final_package/` — manuscript packaging (private).",
          "- `AGENT_NOTES.md`, `<private>/` — agent memory (project rule: never push).",
          "- Build artifacts (`*.aux`, `*.log`, `*.bbl`, `__pycache__/`).",
          "",
          "## Sanitization",
          "",
          "- deCODE access token UUID → `<DECODE_TOKEN>` in all text files.",
          "- Generic UUID, `SYNAPSE_AUTH_TOKEN`, and `api_key` patterns redacted defensively.",
          "- Scripts read the token from the `DECODE_TOKEN` environment variable; the fallback hard-coded UUID has been replaced.",
          "",
          "## Per-directory inventory",
          ""]
    for d in sorted(by_dir):
        files_in_d = by_dir[d]
        md.append(f"### `{d}` ({len(files_in_d)} files)")
        md.append("")
        for name, sz in sorted(files_in_d):
            md.append(f"- `{name}` ({sz:,} bytes)")
        md.append("")
    (GH / "FILES_INCLUDED.md").write_text("\n".join(md))
    log.append("  generated: FILES_INCLUDED.md")


def final_secret_scan() -> list[str]:
    issues = []
    for f in GH.rglob("*"):
        if not f.is_file():
            continue
        try:
            text = f.read_text(errors="replace")
        except Exception:
            continue
        if DECODE_TOKEN_UUID in text:
            issues.append(f"{f.relative_to(GH)}: deCODE token UUID still present!")
        # Generic 36-char UUID detector (won't match in github URLs etc.; conservative)
        for m in re.finditer(r"\b([a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12})\b", text):
            issues.append(f"{f.relative_to(GH)}: residual UUID-like token '{m.group(1)}'")
    return issues


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def main() -> int:
    if GH.exists():
        shutil.rmtree(GH)
    GH.mkdir()
    log = []

    # Top-level files
    log.append("=== Root files ===")
    for name in INCLUDE_FILES:
        src = ROOT / name
        if not src.exists():
            log.append(f"  WARN: missing {name}")
            continue
        ok, note = _copy_file(src, GH / name)
        log.append(f"  {note}: {name}")

    # Directories
    log.append("=== Directories ===")
    for name in INCLUDE_DIRS:
        src = ROOT / name
        log.append(f"--- {name}/ ---")
        _copy_tree(src, GH / name, log)

    # Selective docs
    log.append("=== Selective docs ===")
    for name in INCLUDE_DOCS:
        src = ROOT / name
        if not src.exists():
            log.append(f"  WARN: missing {name}")
            continue
        ok, note = _copy_file(src, GH / name)
        log.append(f"  {note}: {name}")

    # Meta files
    log.append("=== Meta files ===")
    write_meta_files(log)
    write_manifest(log)

    # Final secret scan
    log.append("=== Final secret scan ===")
    issues = final_secret_scan()
    if issues:
        for iss in issues:
            log.append(f"  LEAK: {iss}")
        log.append(f"FAIL: {len(issues)} potential secret(s) remaining.")
    else:
        log.append("  PASS: no deCODE token / UUID-like secrets in github_ready/.")

    print("\n".join(log[-80:]))  # tail summary
    print("\n=== SUMMARY ===")
    n_files = sum(1 for p in GH.rglob("*") if p.is_file())
    print(f"github_ready/ files: {n_files}")
    print(f"Top dir: {GH}")
    if issues:
        print(f"FAIL: {len(issues)} secret(s) detected; aborting before user inspection.")
        return 2
    print("PASS: clean for inspection. NOT pushed (per user request).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
