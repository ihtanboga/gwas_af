#!/usr/bin/env python3
"""D033 Step 1 — diagnose why 32/49 smoke instruments did NOT match AF.

For each instrument:
  1. Try chr:pos + allele match in AF (existing path).
  2. If miss: try rsid fallback in AF.
  3. Classify unmatched reasons:
     - rare_pqtl_below_AF_maf (deCODE MAF < 0.01)
     - indel_or_multibase  (one of effect/other allele is not single-base ACGT)
     - allele_set_mismatch (chr:pos in AF but alleles differ entirely)
     - rsid_only_match     (recovered via rsid fallback; chr:pos diverges)
     - missing_truly       (neither chr:pos nor rsid in AF)
     - other

Outputs:
  results/qc/decode_smoke_match_diagnosis.tsv
  results/qc/decode_smoke_match_diagnosis_summary.md
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import polars as pl

ROOT = Path("/Users/apple/Desktop/gwas_af")
INSTR = ROOT / "data/processed/instruments/decode_2021_smoke.parquet"
AF = ROOT / "data/processed/sumstats/AF.parquet"
OUT_TSV = ROOT / "results/qc/decode_smoke_match_diagnosis.tsv"
OUT_MD = ROOT / "results/qc/decode_smoke_match_diagnosis_summary.md"

VALID_ALLELES = {"A", "C", "G", "T"}


def _log(msg: str) -> None:
    print(f"[match-diag {time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> int:
    instr = pl.read_parquet(INSTR)
    _log(f"loaded {instr.height} instruments")

    # Build AF lookup datasets: by (chr, pos), by rsid.
    af = pl.read_parquet(AF).select([
        "chr", "pos", "rsid", "effect_allele", "other_allele",
        "beta", "se", "pval", "eaf",
    ])
    _log(f"AF parquet: {af.height:,} rows")

    # Index by (chr,pos) — use polars join for efficiency on the smoke set
    keys_chrpos = pl.DataFrame({
        "chr": [str(c) for c in instr["chrom"].to_list()],
        "pos": [int(p) for p in instr["pos"].to_list()],
    })
    af_chrpos_hit = af.join(keys_chrpos, on=["chr", "pos"], how="inner")
    af_chrpos_dict = {(str(r["chr"]), int(r["pos"])): r
                       for r in af_chrpos_hit.iter_rows(named=True)}

    # rsid fallback: keep only instruments with a non-empty, non-rs???? string starting with 'rs'
    instr_rsids = [r for r in instr["rsid"].to_list() if r and str(r).startswith("rs")]
    af_rsid_hit = af.filter(pl.col("rsid").is_in(instr_rsids))
    af_rsid_dict = {str(r["rsid"]): r for r in af_rsid_hit.iter_rows(named=True)}
    _log(f"AF chr:pos hits {len(af_chrpos_dict)}, rsid hits {len(af_rsid_dict)}")

    rows = []
    for r in instr.iter_rows(named=True):
        chrom = str(r["chrom"])
        pos = int(r["pos"])
        rsid = str(r.get("rsid") or "")
        ea = (r["effect_allele"] or "").upper()
        oa = (r["other_allele"] or "").upper()
        maf = r.get("imp_maf")

        single_base = (ea in VALID_ALLELES) and (oa in VALID_ALLELES)

        chrpos_hit = af_chrpos_dict.get((chrom, pos))
        rsid_hit = af_rsid_dict.get(rsid) if rsid.startswith("rs") else None

        path = ""
        reason = ""
        af_ea = af_oa = ""
        af_pos = None

        if chrpos_hit is not None:
            af_ea = (chrpos_hit["effect_allele"] or "").upper()
            af_oa = (chrpos_hit["other_allele"] or "").upper()
            af_pos = int(chrpos_hit["pos"])
            if {ea, oa} == {af_ea, af_oa}:
                path = "matched_chrpos"
                reason = "ok"
            else:
                path = "chrpos_present_alleles_differ"
                reason = "allele_set_mismatch"
        elif rsid_hit is not None:
            af_ea = (rsid_hit["effect_allele"] or "").upper()
            af_oa = (rsid_hit["other_allele"] or "").upper()
            af_pos = int(rsid_hit["pos"])
            if af_pos != pos:
                path = "rsid_only_match"
                reason = "rsid_recoverable_pos_drift"
            elif {ea, oa} == {af_ea, af_oa}:
                path = "rsid_only_match"
                reason = "ok_rsid_fallback"
            else:
                path = "rsid_only_match"
                reason = "rsid_match_alleles_differ"
        else:
            # truly missing — diagnose by feature
            if not single_base:
                reason = "indel_or_multibase"
            elif maf is not None and float(maf) < 0.01:
                reason = "rare_pqtl_below_AF_maf"
            else:
                reason = "missing_truly"
            path = "no_match"

        rows.append({
            "gene": r["gene"], "aptamer": r["aptamer"],
            "chrom": chrom, "pos": pos, "rsid": rsid,
            "ea_decode": ea, "oa_decode": oa,
            "imp_maf_decode": maf,
            "single_base": single_base,
            "af_pos": af_pos, "af_ea": af_ea, "af_oa": af_oa,
            "match_path": path, "reason": reason,
        })

    df = pl.DataFrame(rows)
    df.write_csv(OUT_TSV, separator="\t")
    _log(f"wrote {OUT_TSV}")

    # Summary by reason
    counts = df.group_by("reason").agg(pl.len().alias("n")).sort("n", descending=True)
    summary_lines = [
        "# deCODE smoke match diagnosis (D033 Step 1)\n",
        f"Total instruments: {df.height}",
        "",
        "## Match outcome breakdown",
        "",
        "| reason | n |",
        "|---|---:|",
    ]
    for r in counts.iter_rows(named=True):
        summary_lines.append(f"| {r['reason']} | {r['n']} |")
    summary_lines.append("")

    # Recoverable via rsid fallback
    rsid_recoverable = df.filter(pl.col("reason").is_in([
        "ok_rsid_fallback", "rsid_recoverable_pos_drift",
    ])).height
    summary_lines.append(f"**rsid fallback recoverable:** {rsid_recoverable} additional instruments.")
    summary_lines.append("")

    # Per-gene summary of ok vs missing
    per_gene = df.group_by("gene").agg([
        pl.len().alias("n_total"),
        (pl.col("reason") == "ok").sum().alias("n_chrpos_ok"),
        pl.col("reason").is_in(["ok_rsid_fallback", "rsid_recoverable_pos_drift"]).sum().alias("n_rsid_ok"),
        pl.col("reason").is_in(["allele_set_mismatch"]).sum().alias("n_allele_mismatch"),
        pl.col("reason").is_in(["rare_pqtl_below_AF_maf"]).sum().alias("n_rare"),
        pl.col("reason").is_in(["indel_or_multibase"]).sum().alias("n_indel"),
        pl.col("reason").is_in(["missing_truly"]).sum().alias("n_missing"),
    ]).sort("gene")
    summary_lines.append("## Per gene")
    summary_lines.append("")
    summary_lines.append("| gene | total | chrpos_ok | rsid_ok | allele_mismatch | rare | indel | truly_missing |")
    summary_lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for r in per_gene.iter_rows(named=True):
        summary_lines.append(
            f"| {r['gene']} | {r['n_total']} | {r['n_chrpos_ok']} | {r['n_rsid_ok']} | "
            f"{r['n_allele_mismatch']} | {r['n_rare']} | {r['n_indel']} | {r['n_missing']} |"
        )
    OUT_MD.write_text("\n".join(summary_lines))
    _log(f"wrote {OUT_MD}")

    matched_chrpos = (df["reason"] == "ok").sum()
    matched_rsid = df["reason"].is_in(["ok_rsid_fallback"]).sum()
    rate_chrpos = matched_chrpos / df.height
    rate_with_rsid = (matched_chrpos + matched_rsid) / df.height
    _log(f"current match rate (chrpos only): {matched_chrpos}/{df.height} = {rate_chrpos:.1%}")
    _log(f"with rsid fallback: {matched_chrpos + matched_rsid}/{df.height} = {rate_with_rsid:.1%}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
