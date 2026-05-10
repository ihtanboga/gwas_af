# utils.R — shared helpers for the AFSHF-TARGETMAP R pipeline.

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
})

afshf_load_config <- function(path) {
  yaml::read_yaml(path)
}

afshf_die <- function(msg) {
  stop(paste0("[afshf] ", msg), call. = FALSE)
}

afshf_log <- function(msg) {
  cat(sprintf("[afshf %s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

# Standard summary-stat schema (mirror src/afshf/harmonize.py STANDARD_COLUMNS).
afshf_standard_columns <- c(
  "trait", "source", "variant_id", "rsid", "chr", "pos", "build",
  "effect_allele", "other_allele", "beta", "se", "or", "pval", "eaf",
  "n", "n_cases", "n_controls", "ancestry"
)

afshf_check_standard_columns <- function(dt) {
  missing_cols <- setdiff(afshf_standard_columns, names(dt))
  if (length(missing_cols) > 0L) {
    afshf_die(sprintf("Missing required columns: %s",
                      paste(missing_cols, collapse = ", ")))
  }
  invisible(dt)
}
