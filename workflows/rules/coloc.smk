"""Phase 4 — QTL <-> outcome colocalization (coloc.abf primary, SuSiE-coloc)."""


rule prepare_coloc_loci:
    input:
        mr="results/mr/{qtl_source}__{outcome}.tsv",
        outcome="data/processed/sumstats/{outcome}.parquet",
        instruments="data/processed/instruments/{qtl_source}.parquet",
    output:
        "data/processed/coloc_regions/{qtl_source}__{outcome}.parquet",
    log:
        "logs/phase4/prepare_{qtl_source}__{outcome}.log",
    shell:
        r"""
        Rscript src/R/04_run_coloc.R prepare \
            --mr {input.mr} \
            --outcome {input.outcome} \
            --instruments {input.instruments} \
            --output {output} \
            > {log} 2>&1
        """


rule run_coloc:
    input:
        regions="data/processed/coloc_regions/{qtl_source}__{outcome}.parquet",
    output:
        "results/coloc/{qtl_source}__{outcome}.tsv",
    log:
        "logs/phase4/coloc_{qtl_source}__{outcome}.log",
    shell:
        r"""
        Rscript src/R/04_run_coloc.R run \
            --regions {input.regions} \
            --output {output} \
            > {log} 2>&1
        """
