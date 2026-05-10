"""AFSHF-TARGETMAP top-level Snakefile.

Phase-gated rule chain:

    Phase 0  download metadata + dataset registry
    Phase 1  harmonize sumstats + build cis instruments
    Phase 2  outcome-outcome overlap (LDSC / HDL / HyPrColoc)
    Phase 3  target -> outcome two-sample MR
    Phase 4  QTL <-> outcome colocalization
    Phase 5  intervention direction + therapeutic balance scoring + classification
    Phase 6  druggability + single-cell annotation + Quarto report

Invoke specific phases via dedicated targets, e.g.

    snakemake -j8 phase1
    snakemake -j8 results/final/target_pleiotropy_master.tsv

Heavy steps (downloads, MR, coloc) are written as wildcards over the logical
outcome panel and the QTL source panel, both defined in `config/`.
"""

from pathlib import Path

import yaml

configfile: "config/project.yml"


# ---- load auxiliary configs ------------------------------------------------

with open("config/outcomes.yml") as f:
    OUTCOMES_CFG = yaml.safe_load(f)
with open("config/qtl_sources.yml") as f:
    QTL_CFG = yaml.safe_load(f)


PRIMARY_OUTCOMES = OUTCOMES_CFG["panels"]["primary"]
SAFETY_OUTCOMES = OUTCOMES_CFG["panels"]["safety"]
ALL_OUTCOMES = PRIMARY_OUTCOMES + SAFETY_OUTCOMES
QTL_SOURCES = list(QTL_CFG["qtl_sources"].keys())
PRIMARY_QTL_SOURCES = [k for k, v in QTL_CFG["qtl_sources"].items() if v["role"] == "discovery"]


# ---- include rule modules --------------------------------------------------

include: "workflows/rules/download.smk"
include: "workflows/rules/harmonize.smk"
include: "workflows/rules/instruments.smk"
include: "workflows/rules/overlap.smk"
include: "workflows/rules/mr.smk"
include: "workflows/rules/coloc.smk"
include: "workflows/rules/annotation.smk"
include: "workflows/rules/report.smk"


# ---- top-level targets -----------------------------------------------------

rule all:
    input:
        "results/final/target_pleiotropy_master.tsv",
        "results/final/red_flag_targets.tsv",
        "results/final/green_shared_benefit_targets.tsv",
        "results/final/shared_loci_af_hf_stroke.tsv",
        "results/final/druggability_annotation.tsv",
        "reports/afshf_targetmap.html",


rule phase0_metadata:
    input:
        "data/registry/dataset_metadata.tsv",
        "data/registry/qtl_source_metadata.tsv",


rule phase1_harmonize:
    input:
        expand(
            "data/processed/sumstats/{outcome}.parquet",
            outcome=ALL_OUTCOMES,
        ),


rule phase1_instruments:
    input:
        expand(
            "data/processed/instruments/{qtl_source}.parquet",
            qtl_source=PRIMARY_QTL_SOURCES,
        ),


rule phase2_overlap:
    input:
        "results/overlap/genetic_correlations.tsv",
        "results/overlap/shared_loci.tsv",


rule phase3_mr:
    input:
        expand(
            "results/mr/{qtl_source}__{outcome}.tsv",
            qtl_source=PRIMARY_QTL_SOURCES,
            outcome=PRIMARY_OUTCOMES,
        ),


rule phase4_coloc:
    input:
        expand(
            "results/coloc/{qtl_source}__{outcome}.tsv",
            qtl_source=PRIMARY_QTL_SOURCES,
            outcome=PRIMARY_OUTCOMES,
        ),


rule phase5_classify:
    input:
        "results/final/target_pleiotropy_master.tsv",


rule phase6_report:
    input:
        "reports/afshf_targetmap.html",
