"""Phase 1B — build cis-pQTL / cis-eQTL instrument tables."""


rule build_instruments:
    input:
        registry="config/qtl_sources.yml",
        thresholds="config/analysis_thresholds.yml",
    output:
        "data/processed/instruments/{qtl_source}.parquet",
    log:
        "logs/phase1/instruments_{qtl_source}.log",
    shell:
        r"""
        Rscript src/R/02_build_qtl_instruments.R \
            --qtl-source {wildcards.qtl_source} \
            --registry {input.registry} \
            --thresholds {input.thresholds} \
            --output {output} \
            > {log} 2>&1
        """
