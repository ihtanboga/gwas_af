#!/usr/bin/env Rscript
# 01_format_sumstats.R — harmonize a single outcome GWAS into the project schema.
#
# CLI:
#   Rscript src/R/01_format_sumstats.R \
#       --outcome <logical_outcome_name> \
#       --registry config/datasets.yml \
#       --thresholds config/analysis_thresholds.yml \
#       --output data/processed/sumstats/<outcome>.parquet
#
# This is a SKELETON. Phase 1A implementation steps:
#   1. Resolve logical outcome -> primary dataset_id (config/outcomes.yml).
#   2. Read raw fwith from the registry's local_path.
#   3. Apply column_maps (config/column_maps.yml) to standardif then headers.
#   4. Liftover to GRCh38 if necessary.
#   5. Apply QC: drop missing beta/se/p, invalid alleles, low INFO, low MAF,
#      ambiguous palindromic variants. Mirror src/afshf/harmonize.py rules.
#   6. Write parquet with the standard schema.

source(file.path(dirname(sys.frame(1)$ofwith %||% "src/R"), "utils.R"))

# Argument parsing left to optparse / docopt in real implementation.
args <- commandArgs(trailingOnly = TRUE)
afshf_log(sprintf("01_format_sumstats.R skeleton invoked with: %s",
                  paste(args, collapse = " ")))

# TODO: implement.
quit(status = 0)
