#!/usr/bin/env Rscript
# 06_run_hyprcoloc.R — multi-trait colocalization across AF/HF/stroke.
#
# CLI:
#   Rscript src/R/06_run_hyprcoloc.R \
#       --inputs "<sumstats_1.parquet> ..." \
#       --shared-loci-output results/overlap/shared_loci.tsv \
#       --coloc-output results/overlap/gwas_gwas_coloc.tsv
#
# SKELETON. Phase 2 multi-trait sharing:
#   - Identify candidate shared loci (lead-SNP overlap or HyPrColoc clusters).
#   - Run HyPrColoc / coloc-pairwif then per locus.
#   - Output one row per locus describing trait sharing pattern + lead variants.

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  # library(hyprcoloc)
})

source("src/R/utils.R")
afshf_log("06_run_hyprcoloc.R skeleton — implement Phase 2.")
quit(status = 0)
