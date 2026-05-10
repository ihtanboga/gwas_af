#!/usr/bin/env Rscript
# 02_build_qtl_instruments.R — extract cis instruments per target from a QTL source.
#
# CLI:
#   Rscript src/R/02_build_qtl_instruments.R \
#       --qtl-source <ukbppp|decode|gtex_v8|...> \
#       --registry config/qtl_sources.yml \
#       --thresholds config/analysis_thresholds.yml \
#       --output data/processed/instruments/<qtl_source>.parquet
#
# SKELETON. Implementation steps:
#   1. Load gene/protein annotation (Ensembl coordinates) for the assay panel.
#   2. For each target gene/protein:
#        a. Subset QTL associations to the cis window (config: cis_window_kb).
#        b. Apply primary p-threshold (config: qtl_p_primary).
#        c. Compute F-stat, drop F < min_f_stat.
#        d. LD clump (PLINK) at clump_r2, clump_kb against the EUR reference.
#   3. Write parquet with columns:
#        target_id, gene_symbol, ensembl_gene_id, uniprot_id (if pQTL),
#        instrument_snp, chr, pos, effect_allele, other_allele, beta, se, pval,
#        eaf, f_stat, qtl_source, ancestry, build.

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
})

source("src/R/utils.R")
afshf_log("02_build_qtl_instruments.R skeleton — implement Phase 1B.")
quit(status = 0)
