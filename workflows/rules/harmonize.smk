"""Phase 1A — harmonize raw outcome / QTL summary statistics."""


rule harmonize_outcome:
    input:
        registry="config/datasets.yml",
        thresholds="config/analysis_thresholds.yml",
        column_maps="config/column_maps.yml" if Path("config/column_maps.yml").exists() else "config/datasets.yml",
    output:
        "data/processed/sumstats/{outcome}.parquet",
    log:
        "logs/phase1/harmonize_{outcome}.log",
    shell:
        r"""
        Rscript src/R/01_format_sumstats.R \
            --outcome {wildcards.outcome} \
            --registry {input.registry} \
            --thresholds {input.thresholds} \
            --output {output} \
            > {log} 2>&1
        """
