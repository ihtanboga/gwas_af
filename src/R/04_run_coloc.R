#!/usr/bin/env Rscript
# 04_run_coloc.R — QTL <-> outcome colocalization.
#
# CLI sub-commands:
#   Rscript src/R/04_run_coloc.R prepare \
#       --mr <mr.tsv> --outcome <outcome.parquet> \
#       --instruments <instruments.parquet> \
#       --output <regions.parquet>
#
#   Rscript src/R/04_run_coloc.R run \
#       --regions <regions.parquet> \
#       --output <coloc.tsv>
#
# SKELETON. Phase 4 implementation steps:
#   prepare:
#     - For every (target, outcome) pair where MR p < exploratory threshold,
#       extract a coloc region (lead +/- region_kb) including ALL SNPs
#       (no significance filtering — spec critical rule).
#     - Persist as parquet for fast re-use.
#   run:
#     - For each region, run coloc::coloc.abf (priors from thresholds.yml).
#     - When LD reference is available, also run coloc::coloc.susie.
#     - Output: target, outcome, region, PP.H0..PP.H4,
#               H4_over_H3, p12_used, ld_diagnostic_flag, susie_credible_set_id.

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  # library(coloc)
  # library(susieR)
})

source("src/R/utils.R")

args <- commandArgs(trailingOnly = TRUE)
sub <- if (length(args) >= 1L) args[1] else "help"
afshf_log(sprintf("04_run_coloc.R skeleton invoked with sub-command '%s'.", sub))
quit(status = 0)
