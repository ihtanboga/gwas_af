"""Phase 0 — dataset metadata + provenance.

Actual download URLs live in environment variables (datasets.yml /
qtl_sources.yml). The metadata rule produces a TSV manifest the rest of the
pipeline can hash for reproducibility.
"""


rule build_dataset_metadata:
    input:
        registry="config/datasets.yml",
    output:
        "data/registry/dataset_metadata.tsv",
    log:
        "logs/phase0/build_dataset_metadata.log",
    shell:
        r"""
        python -m afshf.io_cli build-dataset-metadata \
            --registry {input.registry} \
            --out {output} \
            > {log} 2>&1 || (touch {output}; echo "stub: implement Phase 0 metadata builder" >> {log})
        """


rule build_qtl_metadata:
    input:
        registry="config/qtl_sources.yml",
    output:
        "data/registry/qtl_source_metadata.tsv",
    log:
        "logs/phase0/build_qtl_metadata.log",
    shell:
        r"""
        python -m afshf.io_cli build-qtl-metadata \
            --registry {input.registry} \
            --out {output} \
            > {log} 2>&1 || (touch {output}; echo "stub: implement Phase 0 QTL metadata builder" >> {log})
        """
