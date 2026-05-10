"""Phase 3 — target -> outcome two-sample MR."""


rule run_target_mr:
    input:
        instruments="data/processed/instruments/{qtl_source}.parquet",
        outcome="data/processed/sumstats/{outcome}.parquet",
        thresholds="config/analysis_thresholds.yml",
    output:
        "results/mr/{qtl_source}__{outcome}.tsv",
    log:
        "logs/phase3/{qtl_source}__{outcome}.log",
    threads: 4
    shell:
        r"""
        Rscript src/R/03_run_twosample_mr.R \
            --instruments {input.instruments} \
            --outcome {input.outcome} \
            --thresholds {input.thresholds} \
            --output {output} \
            --threads {threads} \
            > {log} 2>&1
        """
