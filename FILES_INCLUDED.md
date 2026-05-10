# github_ready/ — file manifest

Total files: 94

## Excluded by design

- Manuscript text, figures, tables, cover letter, author/COI templates (private).
- Raw + processed third-party GWAS / pQTL summary statistics (license).
- Build artifacts and intermediate result files (regenerable from code).
- Local agent memory files (project rule).

## Sanitization

- Access token UUID redacted in all text files.
- Generic UUID, SYNAPSE_AUTH_TOKEN, and api_key patterns redacted defensively.
- Scripts read access tokens from environment variables; no hard-coded fallbacks.
- Internal narrative redacted; documents are in English.

## Per-directory inventory

### (root) (9 files)

- .gitignore (383 bytes)
- ARCHITECTURE.md (7,763 bytes)
- DESCRIPTION (729 bytes)
- FILES_INCLUDED.md (4,853 bytes)
- LICENSE (259 bytes)
- Makefile (1,423 bytes)
- README.md (2,939 bytes)
- Snakefile (3,310 bytes)
- environment.yml (524 bytes)

### config (10 files)

- analysis_thresholds.yml (5,748 bytes)
- column_maps.yml (14,366 bytes)
- datasets.yml (8,754 bytes)
- decode_cvd150_extension.yml (1,136 bytes)
- decode_cvd500_batch.yml (992 bytes)
- outcomes.yml (2,354 bytes)
- positive_controls.tsv (3,808 bytes)
- project.yml (2,383 bytes)
- qtl_sources.yml (3,495 bytes)
- safety_panel.yml (2,643 bytes)

### docs (4 files)

- data_sources.md (5,373 bytes)
- interpretation_rules.md (2,827 bytes)
- limitations.md (2,668 bytes)
- methods.md (2,373 bytes)

### scripts/python (27 files)

- build_d041_interim_evidence_tables.py (16,501 bytes)
- build_d043_manuscript_artifacts.py (44,489 bytes)
- build_d044_audit_artifacts.py (32,297 bytes)
- build_d045_manuscript_package.py (43,762 bytes)
- build_d046_visual_audit.py (13,760 bytes)
- build_d047_submission_drafting.py (53,389 bytes)
- build_d048_final_submission_audit.py (32,915 bytes)
- build_regional_coloc_panels.py (10,683 bytes)
- d039_audit_r_vs_python.py (5,609 bytes)
- prepare_github_release.py (13,062 bytes)
- run_af_2025_qc.py (10,811 bytes)
- run_decode_100aptamer_calibration.py (25,560 bytes)
- run_decode_3A_match_diagnosis.py (7,348 bytes)
- run_decode_3A_smoke.py (18,026 bytes)
- run_decode_3B_3C_smoke.py (23,248 bytes)
- run_decode_cvd150_extension_extraction.py (23,816 bytes)
- run_decode_cvd150_mr_projection.py (18,760 bytes)
- run_decode_cvd150_wave1_coloc.py (31,660 bytes)
- run_decode_cvd150_wave2_coloc.py (17,378 bytes)
- run_decode_cvd454_extraction.py (23,137 bytes)
- run_decode_cvd454_mr_projection.py (19,659 bytes)
- run_gigastroke_2022_eur_qc.py (13,588 bytes)
- run_hf_2024_eur_qc.py (15,446 bytes)
- run_megastroke_eur_qc.py (14,266 bytes)
- select_100aptamer_calibration.py (5,261 bytes)
- select_cvd150_extension.py (13,906 bytes)
- select_cvd500_batch.py (19,785 bytes)

### src/R (7 files)

- 01_format_sumstats.R (1,194 bytes)
- 02_build_qtl_instruments.R (1,229 bytes)
- 03_run_twosample_mr.R (1,373 bytes)
- 04_run_coloc.R (1,312 bytes)
- 05_run_ldsc_or_hdl.R (698 bytes)
- 06_run_hyprcoloc.R (771 bytes)
- utils.R (969 bytes)

### src/afshf (12 files)

- __init__.py (89 bytes)
- annotate_druggability.py (5,957 bytes)
- annotate_singlecell.py (2,696 bytes)
- column_maps.py (15,329 bytes)
- config_schema.py (5,469 bytes)
- harmonize.py (5,348 bytes)
- instruments.py (3,337 bytes)
- intervention_direction.py (6,332 bytes)
- io.py (1,080 bytes)
- plotting.py (1,798 bytes)
- positive_controls.py (4,196 bytes)
- scoring.py (12,751 bytes)

### tests (11 files)

- __init__.py (0 bytes)
- conftest.py (219 bytes)
- test_af_signal_tiers.py (1,838 bytes)
- test_classification.py (3,335 bytes)
- test_column_maps.py (7,575 bytes)
- test_harmonize.py (3,392 bytes)
- test_intervention_direction.py (3,206 bytes)
- test_log10p_underflow.py (8,948 bytes)
- test_protocol_decisions.py (6,600 bytes)
- test_scoring.py (3,230 bytes)
- test_smoke.py (6,884 bytes)

### tests/synthetic_data (6 files)

- cvdkp_af_mini.tsv (1,382 bytes)
- decode_mini.tsv (353 bytes)
- eqtlgen_zonly_mini.tsv (283 bytes)
- finngen_af_mini.tsv (443 bytes)
- gwas_catalog_or_mini.tsv (436 bytes)
- ukbppp_pcsk9_mini.tsv (359 bytes)

### workflows/rules (8 files)

- annotation.smk (2,060 bytes)
- coloc.smk (1,132 bytes)
- download.smk (1,161 bytes)
- harmonize.smk (704 bytes)
- instruments.smk (611 bytes)
- mr.smk (713 bytes)
- overlap.smk (1,108 bytes)
- report.smk (659 bytes)
