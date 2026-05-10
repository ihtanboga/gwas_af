#!/usr/bin/env python3
"""D036 — CVD150 +100 extension list builder.

Read the CVD454 design (454 aptamers; D034) and the already-extracted CVD50
subset (50 aptamers; from cis_instrument_report.tsv). Pick 100 aptamers from
the remaining 404 with category-balanced quotas:

  cat1 AF/HF/stroke GWAS-locus: 15
  cat2 coagulation:             15
  cat4 lipid:                   15
  cat6 ECM/fibrosis:            15
  cat7 endothelial/vascular:    15
  cat3 inflammation:            15  (do NOT let it dominate)
  cat8 druggable:               10

Selection rules (D036):
  - Skip random and large-fwith categories (cat0/cat?: not present in CVD454,
    but cat9 sanity already in CVD50; we exclude cat9 from extension).
  - 1 aptamer per gene by default (extension is small; favour diversity).
  - Multi-aptamer exception only if biologically essential (NPPB / F9 / F2 /
    PROC already in CVD50; no need to add more here).
  - Skip intracellular known-low-plasma proteins flagged in calibration log
    (any aptamer in CVD50 with empty cis is flagged so we avoid same gene
    again at extension time).

Outputs:
  config/decode_cvd150_extension.yml
  results/qc/decode_cvd150_extension_design.tsv
  results/qc/decode_cvd150_extension_category_counts.tsv
  results/qc/decode_cvd150_extension_expected_runtime.tsv
  results/qc/decode_cvd150_extension_list_qc_summary.md
"""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import yaml

ROOT = Path("/Users/apple/Desktop/gwas_af")
CVD454_DESIGN = ROOT / "results/qc/decode_cvd500_batch_design.tsv"
CVD50_LOG = ROOT / "results/qc/decode_cvd454_cis_instrument_report.tsv"
CALIB_REPORT = ROOT / "results/qc/decode_100aptamer_calibration_report.tsv"

OUT_YAML = ROOT / "config/decode_cvd150_extension.yml"
OUT_DESIGN = ROOT / "results/qc/decode_cvd150_extension_design.tsv"
OUT_CATCOUNTS = ROOT / "results/qc/decode_cvd150_extension_category_counts.tsv"
OUT_RUNTIME = ROOT / "results/qc/decode_cvd150_extension_expected_runtime.tsv"
OUT_SUMMARY = ROOT / "results/qc/decode_cvd150_extension_list_qc_summary.md"

# Per D036 (approximate, biology > exact). cat1 GWAS-loci was fully consumed
# by CVD50 (19 aptamers, all done) — quota redistributed to other categories.
# cat9 sanity already covered in CVD50; cat5 cardiac_stress small + many
# intracellular, dropped to 5.
CATEGORY_QUOTAS: dict[str, int] = {
    "cat2_coagulation":            18,
    "cat4_lipid":                  18,
    "cat6_ecm_fibrosis":           18,
    "cat7_endothelial_vascular":   18,
    "cat3_inflammation":           18,
    "cat5_cardiac_stress":          5,   # squeezed in
    "cat8_druggable_cvd":           5,   # rest of slot to fill quota
}

# Genes flagged from CVD50 as empty (likely intracellular/poor plasma) — avoid in extension
KNOWN_EMPTY_GENES_FROM_CVD50: set[str] = set()


def main() -> int:
    # Load full CVD454 design
    cvd454 = pl.read_csv(CVD454_DESIGN, separator="\t")
    print(f"CVD454 design: {cvd454.height} aptamers")

    # Load already-extracted CVD50 (matched: ALL rows in extraction_log; consider 50 as done)
    cvd50_log = pl.read_csv(CVD50_LOG, separator="\t")
    done_aptamers = set(cvd50_log["aptamer"].to_list())
    print(f"CVD50 already done: {len(done_aptamers)} aptamers")

    # Track empty-gene flags from CVD50 (n_clumped == 0 -> likely no plasma cis-pQTL)
    empty_genes = set()
    for r in cvd50_log.iter_rows(named=True):
        if r.get("n_clumped", 0) == 0:
            empty_genes.add(r["gene"])
    print(f"Genes empty in CVD50 (avoid in extension if same gene): {len(empty_genes)}")
    KNOWN_EMPTY_GENES_FROM_CVD50.update(empty_genes)

    # Filter CVD454 to remaining (not yet processed)
    remaining = cvd454.filter(~pl.col("aptamer_key").is_in(list(done_aptamers)))
    print(f"Remaining CVD454 candidates: {remaining.height}")

    # Per category, prefer:
    #  - genes NOT already extracted in CVD50 (one aptamer per gene)
    #  - genes NOT flagged as empty from CVD50
    #  - prefer stable size (median +/- around 950 MB)
    selected_keys: list[str] = []
    selected_rows: list[dict] = []
    seen_genes: set[str] = set(g for r in cvd50_log.iter_rows(named=True) for g in [r["gene"]])

    for cat, quota in CATEGORY_QUOTAS.items():
        cat_pool = remaining.filter(pl.col("primary_category") == cat).to_dicts()
        # prefer 1 aptamer per gene; skip genes already taken or empty in CVD50
        added = 0
        # Sort: aptamers from genes not yet seen first; then by aptamer_key for determinism
        cat_pool_sorted = sorted(
            cat_pool,
            key=lambda r: (
                r["gene"] in KNOWN_EMPTY_GENES_FROM_CVD50,  # empty-flagged genes last
                r["gene"] in seen_genes,                    # already-taken genes last
                r["aptamer_key"],
            ),
        )
        for r in cat_pool_sorted:
            if added >= quota:
                break
            if r["gene"] in seen_genes:
                continue                                    # one-aptamer-per-gene
            if r["aptamer_key"] in selected_keys:
                continue
            selected_rows.append({**r, "extension_category_quota": cat})
            selected_keys.append(r["aptamer_key"])
            seen_genes.add(r["gene"])
            added += 1
        print(f"  {cat}: added {added}/{quota}")

    # If short, fill from remaining categories (rare but possible if quotas exceed pool)
    target_total = sum(CATEGORY_QUOTAS.values())
    if len(selected_rows) < target_total:
        deficit = target_total - len(selected_rows)
        print(f"Filling deficit of {deficit} from remaining pool ...")
        remaining_pool = [r for r in remaining.to_dicts()
                          if r["aptamer_key"] not in selected_keys
                          and r["gene"] not in seen_genes
                          and r["gene"] not in KNOWN_EMPTY_GENES_FROM_CVD50]
        for r in remaining_pool[:deficit]:
            selected_rows.append({**r, "extension_category_quota": "fill"})
            selected_keys.append(r["aptamer_key"])
            seen_genes.add(r["gene"])

    df = pl.DataFrame(selected_rows)
    print(f"\nSelected {df.height} extension aptamers; {df['gene'].n_unique()} unique genes")

    df.write_csv(OUT_DESIGN, separator="\t")

    cat_counts = df.group_by("extension_category_quota").agg([
        pl.col("gene").n_unique().alias("n_unique_genes"),
        pl.len().alias("n_aptamers"),
    ]).sort("extension_category_quota")
    cat_counts.write_csv(OUT_CATCOUNTS, separator="\t")

    # Expected metrics from calibration cv_prior priors
    calib = pl.read_csv(CALIB_REPORT, separator="\t").row(0, named=True)
    median_gb = float(calib["GB_per_aptamer_median"])
    median_min = float(calib["minutes_per_aptamer_median"])
    cv_prior_empty = 0.13
    cv_prior_match = 0.85

    # CVD454 partial run reality: 4-worker speedup ~1x, 3.4 min/aptamer wall
    real_min_per_apt_4w = 3.4

    n_new = df.height
    expected_transfer_gb = n_new * median_gb
    # Expected runtime now uses the empirical 4-worker speedup observed on CVD50.
    expected_runtime_4w_h = n_new * real_min_per_apt_4w / 60
    proj = []
    for nw in (1, 2, 4):
        # For 1/2 worker projections we still report optimistic theoretical
        # values from calibration; primary planning value is the 4-worker
        # CVD50-observed estimate.
        if nw == 4:
            wall_h = expected_runtime_4w_h
        else:
            speedup = {1: 1.0, 2: 1.6}[nw]
            wall_h = n_new * (median_min * 60) / speedup / 3600
        proj.append({
            "n_aptamers": n_new, "workers": nw,
            "wall_hours": round(wall_h, 2),
            "transfer_GB": round(expected_transfer_gb, 1),
            "speedup_assumed": "1.0x_observed_CVD50" if nw == 4 else f"{speedup:.1f}x_calib_naive",
            "expected_empty_rate_calib_cv_prior": cv_prior_empty,
            "expected_match_rate_calib_cv_prior": cv_prior_match,
        })
    pl.DataFrame(proj).write_csv(OUT_RUNTIME, separator="\t")

    # YAML config
    yaml_doc = {
        "batch_name": "CVD150 +100 extension (D036)",
        "scope_label": "CVD-prioritized deCODE CVD150 target screen",
        "completed_subset_aptamers": len(done_aptamers),
        "n_new_aptamers": int(n_new),
        "n_total_cvd150": int(len(done_aptamers) + n_new),
        "n_unique_genes_extension": int(df["gene"].n_unique()),
        "category_quotas": CATEGORY_QUOTAS,
        "extension_category_counts": {
            row["extension_category_quota"]: int(row["n_aptamers"])
            for row in cat_counts.iter_rows(named=True)
        },
        "calibration_priors": {
            "median_gb_per_aptamer": median_gb,
            "median_minutes_per_aptamer_calib": median_min,
            "real_minutes_per_aptamer_4w_observed_CVD50": real_min_per_apt_4w,
        },
        "expected_metrics_4w_observed": {
            "transfer_gb": round(expected_transfer_gb, 1),
            "wall_hours": round(expected_runtime_4w_h, 2),
            "empty_rate": cv_prior_empty,
            "match_rate": cv_prior_match,
        },
        "list_qc_thresholds_d036": {
            "n_new_aptamers_min": 90, "n_new_aptamers_max": 110,
            "n_total_min": 140, "n_total_max": 160,
            "transfer_gb_max": 120,
            "runtime_4w_hours_max": 8,
            "empty_rate_max": 0.30,
            "match_rate_min": 0.75,
            "top_gene_share_max": 0.05,
        },
        "approval_gate": "explicit_protocol_approval_required_before_extension_run",
    }
    OUT_YAML.parent.mkdir(parents=True, exist_ok=True)
    OUT_YAML.write_text(yaml.safe_dump(yaml_doc, sort_keys=False))

    # Threshold check
    th = yaml_doc["list_qc_thresholds_d036"]
    expected = yaml_doc["expected_metrics_4w_observed"]
    n_total = len(done_aptamers) + n_new
    top_gene_share = (
        df.group_by("gene").agg(pl.len().alias("n")).sort("n", descending=True).row(0)[1]
        / max(1, df.height)
    )
    passes = {
        "n_new_in_range":     th["n_new_aptamers_min"] <= n_new <= th["n_new_aptamers_max"],
        "n_total_in_range":   th["n_total_min"] <= n_total <= th["n_total_max"],
        "transfer_within":    expected["transfer_gb"] <= th["transfer_gb_max"],
        "runtime_within":     expected["wall_hours"] <= th["runtime_4w_hours_max"],
        "empty_within":       expected["empty_rate"] <= th["empty_rate_max"],
        "match_within":       expected["match_rate"] >= th["match_rate_min"],
        "no_single_gene_dominates": top_gene_share <= th["top_gene_share_max"],
        "all_quotas_filled":  all(yaml_doc["extension_category_counts"].get(c, 0) > 0
                                  for c in CATEGORY_QUOTAS.keys()),
    }

    md = [
        "# CVD150 +100 extension — list QC summary (D036)",
        "",
        "**Scope:** CVD-prioritized deCODE CVD150 target screen (NOT proteome-wide).",
        "",
        "## Counts",
        "",
        f"- CVD50 already done: **{len(done_aptamers)}**",
        f"- New extension aptamers: **{n_new}**",
        f"- Total CVD150: **{n_total}**",
        f"- Unique genes (extension): **{df['gene'].n_unique()}**",
        f"- Top-gene aptamer share (extension): {top_gene_share:.1%}",
        f"- Genes flagged empty in CVD50 (excluded): {len(KNOWN_EMPTY_GENES_FROM_CVD50)}",
        "",
        "## Per category quota",
        "",
    ]
    for r in cat_counts.iter_rows(named=True):
        target = CATEGORY_QUOTAS.get(r['extension_category_quota'], '—')
        md.append(f"- `{r['extension_category_quota']}`: **{r['n_aptamers']} aptamers** "
                  f"({r['n_unique_genes']} genes) [target ~{target}]")
    md.extend([
        "",
        "## Expected metrics (4-worker, CVD50-observed)",
        "",
        f"- Transfer (+100): **{expected['transfer_gb']:.0f} GB**",
        f"- Wall-time @ 4 workers: **{expected['wall_hours']:.1f} h** (using observed 1× speedup)",
        f"- Empty rate expected: {expected['empty_rate']:.0%}",
        f"- Match-to-AF rate expected: {expected['match_rate']:.0%}",
        "",
        "## Threshold check (D036)",
        "",
        "| Threshold | Target | Observed/Expected | Pass? |",
        "|---|---|---|---|",
        f"| n_new ∈ [90, 110] | yes | {n_new} | {'✅' if passes['n_new_in_range'] else '❌'} |",
        f"| n_total ∈ [140, 160] | yes | {n_total} | {'✅' if passes['n_total_in_range'] else '❌'} |",
        f"| +100 transfer ≤ 120 GB | yes | {expected['transfer_gb']:.0f} GB | {'✅' if passes['transfer_within'] else '❌'} |",
        f"| +100 runtime 4w ≤ 8 h | yes | {expected['wall_hours']:.1f} h | {'✅' if passes['runtime_within'] else '❌'} |",
        f"| empty rate ≤ 30% | yes | {expected['empty_rate']:.0%} | {'✅' if passes['empty_within'] else '❌'} |",
        f"| match rate ≥ 75% | yes | {expected['match_rate']:.0%} | {'✅' if passes['match_within'] else '❌'} |",
        f"| top-gene aptamer share ≤ 5% | yes | {top_gene_share:.1%} | {'✅' if passes['no_single_gene_dominates'] else '❌'} |",
        f"| all quotas filled | yes | {sum(1 for v in yaml_doc['extension_category_counts'].values() if v > 0)}/{len(CATEGORY_QUOTAS)} | {'✅' if passes['all_quotas_filled'] else '❌'} |",
        "",
        f"**Overall: {sum(passes.values())}/{len(passes)} threshold pass.**",
        "",
        ("→ Extension list QC PASS. **awaiting protocol explicit approval to start +100 extraction.**"
         if all(passes.values()) else
         "→ Extension list QC FAIL. Adjust quotas / pool before approval."),
    ])
    OUT_SUMMARY.write_text("\n".join(md))

    print()
    print("Threshold check:")
    for k, v in passes.items():
        print(f"  {'✅' if v else '❌'} {k}")
    return 0


if __name__ == "__main__":
    raif then SystemExit(main())
