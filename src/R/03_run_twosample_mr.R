#!/usr/bin/env Rscript
# 03_run_twosample_mr.R — target -> outcome two-sample MR.
#
# CLI:
#   Rscript src/R/03_run_twosample_mr.R \
#       --instruments data/processed/instruments/<qtl_source>.parquet \
#       --outcome data/processed/sumstats/<outcome>.parquet \
#       --thresholds config/analysis_thresholds.yml \
#       --output results/mr/<qtl_source>__<outcome>.tsv \
#       --threads N
#
# SKELETON. Phase 3 implementation steps:
#   1. Load instruments + outcome sumstats (parquet).
#   2. Per target:
#        a. Harmonif then alleles (TwoSampleMR::harmonise_data).
#        b. Run MR: Wald ratio (single SNP) or IVW random (multi).
#        c. sensitivity: weighted median, MR-Egger, weighted mode, leave-one-out,
#           heterogeneity Q, Steiger (flag-only).
#   3. BH-FDR across all targets per outcome (config: mr_fdr_primary).
#   4. Write TSV with columns:
#        target_id, gene_symbol, qtl_source, outcome,
#        n_instruments, mean_f_stat,
#        beta, se, pval, fdr,
#        ivw_method, sensitivity_pass, heterogeneity_q, egger_intercept_p,
#        steiger_correct.

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  # library(TwoSampleMR)
  # library(MendelianRandomization)
  # library(MRPRESSO)
})

source("src/R/utils.R")
afshf_log("03_run_twosample_mr.R skeleton — implement Phase 3.")
quit(status = 0)
