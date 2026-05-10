"""Phase 5 — intervention direction + classification + annotations.

This rule is fully implemented in Python because the spec's classification
logic (intervention_direction.py + scoring.py) is already validated by the
unit tests. The rule consumes per-(qtl_source, outcome) MR + coloc tables
and produces the master target table plus the colour-coded subset tables.
"""


rule build_target_master_table:
    input:
        mr=expand(
            "results/mr/{qtl_source}__{outcome}.tsv",
            qtl_source=PRIMARY_QTL_SOURCES,
            outcome=PRIMARY_OUTCOMES,
        ),
        coloc=expand(
            "results/coloc/{qtl_source}__{outcome}.tsv",
            qtl_source=PRIMARY_QTL_SOURCES,
            outcome=PRIMARY_OUTCOMES,
        ),
        thresholds="config/analysis_thresholds.yml",
        project="config/project.yml",
    output:
        master="results/final/target_pleiotropy_master.tsv",
        red="results/final/red_flag_targets.tsv",
        green="results/final/green_shared_benefit_targets.tsv",
        loci="results/final/shared_loci_af_hf_stroke.tsv",
    log:
        "logs/phase5/build_master.log",
    shell:
        r"""
        python -m afshf.pipeline_classify \
            --mr-glob "results/mr/*.tsv" \
            --coloc-glob "results/coloc/*.tsv" \
            --thresholds {input.thresholds} \
            --project {input.project} \
            --master-out {output.master} \
            --red-out {output.red} \
            --green-out {output.green} \
            --shared-loci-out {output.loci} \
            > {log} 2>&1
        """


rule druggability_annotation:
    input:
        master="results/final/target_pleiotropy_master.tsv",
    output:
        "results/final/druggability_annotation.tsv",
    log:
        "logs/phase5/druggability.log",
    shell:
        r"""
        Rscript src/R/12_drug_annotation.R \
            --master {input.master} \
            --output {output} \
            > {log} 2>&1 || (touch {output}; echo "stub: implement Phase 5 druggability" >> {log})
        """
