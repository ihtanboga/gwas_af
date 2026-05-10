"""Phase 6 — Quarto report + figures."""


rule render_report:
    input:
        master="results/final/target_pleiotropy_master.tsv",
        red="results/final/red_flag_targets.tsv",
        green="results/final/green_shared_benefit_targets.tsv",
        loci="results/final/shared_loci_af_hf_stroke.tsv",
        druggability="results/final/druggability_annotation.tsv",
    output:
        "reports/afshf_targetmap.html",
    log:
        "logs/phase6/render_report.log",
    shell:
        r"""
        quarto render reports/afshf_targetmap.qmd \
            --to html \
            --output afshf_targetmap.html \
            > {log} 2>&1
        """
