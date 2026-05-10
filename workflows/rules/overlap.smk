"""Phase 2 — outcome-outcome genetic overlap (LDSC / HDL / HyPrColoc)."""


rule genetic_correlations:
    input:
        sumstats=expand(
            "data/processed/sumstats/{outcome}.parquet",
            outcome=PRIMARY_OUTCOMES,
        ),
    output:
        "results/overlap/genetic_correlations.tsv",
    log:
        "logs/phase2/genetic_correlations.log",
    shell:
        r"""
        Rscript src/R/05_run_ldsc_or_hdl.R \
            --inputs "{input.sumstats}" \
            --output {output} \
            > {log} 2>&1
        """


rule shared_loci_and_hyprcoloc:
    input:
        sumstats=expand(
            "data/processed/sumstats/{outcome}.parquet",
            outcome=PRIMARY_OUTCOMES,
        ),
    output:
        "results/overlap/shared_loci.tsv",
        "results/overlap/gwas_gwas_coloc.tsv",
    log:
        "logs/phase2/shared_loci.log",
    shell:
        r"""
        Rscript src/R/06_run_hyprcoloc.R \
            --inputs "{input.sumstats}" \
            --shared-loci-output {output[0]} \
            --coloc-output {output[1]} \
            > {log} 2>&1
        """
