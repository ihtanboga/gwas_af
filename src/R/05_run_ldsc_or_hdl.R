#!/usr/bin/env Rscript
# 05_run_ldsc_or_hdl.R — outcome-outcome genetic correlation matrix.
#
# CLI:
#   Rscript src/R/05_run_ldsc_or_hdl.R \
#       --inputs "<sumstats_1.parquet> <sumstats_2.parquet> ..." \
#       --output results/overlap/genetic_correlations.tsv
#
# SKELETON. Phase 2 implementation:
#   - Convert each parquet sumstat to LDSC munged format (or HDL format).
#   - Run LDSC pairwif then (or HDL) for the AF/HF/stroke matrix.
#   - Write tidy TSV: trait1, trait2, rg, se, p, intercept, h2_obs.

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
})

source("src/R/utils.R")
afshf_log("05_run_ldsc_or_hdl.R skeleton — implement Phase 2.")
quit(status = 0)
