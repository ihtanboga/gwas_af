"""Positive-control locus registry + smoke-test helpers.

Loads ``config/positive_controls.tsv`` and exposes two helpers:

* ``load_positive_controls`` — return a DataFrame of canonical loci.
* ``check_positive_control_retention`` — given a harmonized sumstat
  DataFrame and a trait label, verify that the cis window of every
  positive control for that trait still contains at least one variant.

These are *sanity* gates. They do NOT verify that the variants are
significant — that comes later in the MR pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd


@dataclass(frozen=True)
class PositiveControl:
    gene_symbol: str
    chr: str
    gene_start: int
    gene_end: int
    build: str
    trait: str
    role: str
    expected_direction: str
    lead_variant_rsid_approx: str
    reference: str


def load_positive_controls(path: str | Path) -> pd.DataFrame:
    """Read the positive-controls TSV (comment lines start with ``#``)."""
    df = pd.read_csv(path, sep="\t", comment="#")
    df.columns = [c.strip() for c in df.columns]
    df["chr"] = df["chr"].astype(str)
    df["gene_start"] = df["gene_start"].astype(int)
    df["gene_end"] = df["gene_end"].astype(int)
    return df


def positive_controls_for_trait(controls: pd.DataFrame, trait: str) -> pd.DataFrame:
    return controls.loc[controls["trait"] == trait].copy()


def cis_bounds(row: pd.Series, window_kb: int) -> tuple[str, int, int]:
    pad = int(window_kb) * 1000
    chrom = str(row["chr"])
    start = max(0, int(row["gene_start"]) - pad)
    end = int(row["gene_end"]) + pad
    return chrom, start, end


@dataclass
class RetentionResult:
    trait: str
    window_kb: int
    total_loci: int
    retained: int
    missing: list[str]

    @property
    def all_retained(self) -> bool:
        return self.retained == self.total_loci


def check_positive_control_retention(sumstats: pd.DataFrame,
                                     controls: pd.DataFrame,
                                     trait: str,
                                     *,
                                     window_kb: int = 1000,
                                     chr_col: str = "chr",
                                     pos_col: str = "pos") -> RetentionResult:
    """Confirm the cis window of every positive control retains >=1 variant.

    Both ``sumstats[chr_col]`` and ``controls['chr']`` are coerced to string
    so ``"4"`` and ``"chr4"`` mismatches surface explicitly.
    """
    if sumstats.empty:
        return RetentionResult(trait=trait, window_kb=window_kb,
                               total_loci=0, retained=0,
                               missing=["<empty sumstats>"])

    pc = positive_controls_for_trait(controls, trait)
    if pc.empty:
        return RetentionResult(trait=trait, window_kb=window_kb,
                               total_loci=0, retained=0, missing=[])

    sumstats = sumstats.copy()
    sumstats[chr_col] = sumstats[chr_col].astype(str)
    sumstats[pos_col] = pd.to_numeric(sumstats[pos_col], errors="coerce").astype("Int64")

    missing: list[str] = []
    retained = 0
    for _, row in pc.iterrows():
        chrom, start, end = cis_bounds(row, window_kb)
        same_chr = sumstats[chr_col] == chrom
        in_range = (sumstats[pos_col] >= start) & (sumstats[pos_col] <= end)
        if (same_chr & in_range).any():
            retained += 1
        else:
            missing.append(str(row["gene_symbol"]))

    return RetentionResult(
        trait=trait,
        window_kb=window_kb,
        total_loci=len(pc),
        retained=retained,
        missing=missing,
    )


def summarise_retention(results: Iterable[RetentionResult]) -> pd.DataFrame:
    rows = []
    for r in results:
        rows.append({
            "trait": r.trait,
            "window_kb": r.window_kb,
            "total_loci": r.total_loci,
            "retained": r.retained,
            "missing_count": len(r.missing),
            "missing": ",".join(r.missing) if r.missing else "",
            "all_retained": r.all_retained,
        })
    return pd.DataFrame(rows)
