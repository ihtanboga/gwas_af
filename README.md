# AFSHF-TARGETMAP

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
scripts/python/   Analysis scripts (extraction, MR, coloc, audit pipelines)
src/              Python utilities (harmonization, scoring, plotting) + R bridge
config/           YAML configs (datasets, QTL sources, thresholds, panels)
tests/            Unit tests
workflows/rules/  Snakemake rules
docs/             Methods, data sources, interpretation rules, limitations
Snakefile, Makefile, environment.yml, DESCRIPTION, ARCHITECTURE.md
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

This codebase enforces explicit claim boundaries:

- All interim labels carry `final_class_allowed = false`.
- The screen is **targeted, not proteome-wide**.
- The screen-level FDR is computed within the CVD150 cis-tested aptamer universe and is **not** a proteome-wide false-discovery rate.
- **No drug-recommendation language** is used for any candidate target.
- The AF anchor is the **all-ancestry** CVDKP common-variant AFGenPlus meta-analysis; HF and stroke outcomes are EUR-only — an `ancestry_mismatch` flag is propagated through every result row.

## Citation

When the manuscript is published, citation details will appear here.

## License

See `LICENSE`. Update this file with your institution's preferred license before public release.
