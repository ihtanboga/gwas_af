"""Instrument selection helpers (cis-pQTL / cis-eQTL).

Heavy lifting (LD clumping, SuSiE) is delegated to R scripts. This module only
provides the pure-Python pieces: cis-window filtering, F-statistic calculation
and a couple of convenience helpers used in tests.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Optional

import pandas as pd


@dataclass(frozen=True)
class CisWindow:
    chr: str
    start: int
    end: int

    def contains(self, chrom: str, pos: int) -> bool:
        return str(chrom) == str(self.chr) and self.start <= int(pos) <= self.end


def make_cis_window(gene_chr: str,
                    gene_start: int,
                    gene_end: int,
                    window_kb: int) -> CisWindow:
    pad = int(window_kb) * 1000
    return CisWindow(
        chr=str(gene_chr),
        start=max(0, int(gene_start) - pad),
        end=int(gene_end) + pad,
    )


def filter_to_cis(df: pd.DataFrame,
                  window: CisWindow,
                  chr_col: str = "chr",
                  pos_col: str = "pos") -> pd.DataFrame:
    same_chr = df[chr_col].astype(str) == window.chr
    in_range = (df[pos_col].astype(int) >= window.start) & (df[pos_col].astype(int) <= window.end)
    return df.loc[same_chr & in_range].copy()


def compute_f_statistic(beta: float, se: float) -> float:
    """Per-SNP F-statistic = (beta/se)^2 — the conventional MR proxy."""
    if not (isfinite(beta) and isfinite(se)) or se == 0.0:
        return 0.0
    return (beta / se) ** 2


def annotate_f_statistic(df: pd.DataFrame,
                         beta_col: str = "beta",
                         se_col: str = "se",
                         out_col: str = "f_stat") -> pd.DataFrame:
    out = df.copy()
    out[out_col] = [
        compute_f_statistic(b, s) for b, s in zip(df[beta_col], df[se_col])
    ]
    return out


def filter_by_pvalue(df: pd.DataFrame,
                     p_threshold: float,
                     pval_col: str = "pval") -> pd.DataFrame:
    return df.loc[df[pval_col].astype(float) < float(p_threshold)].copy()


def filter_by_fstat(df: pd.DataFrame,
                    min_f_stat: float,
                    f_col: str = "f_stat") -> pd.DataFrame:
    if f_col not in df.columns:
        df = annotate_f_statistic(df)
    return df.loc[df[f_col] >= min_f_stat].copy()


def select_instruments_for_gene(qtl_df: pd.DataFrame,
                                gene_chr: str,
                                gene_start: int,
                                gene_end: int,
                                *,
                                window_kb: int = 1000,
                                p_threshold: float = 5e-8,
                                min_f_stat: float = 10.0) -> pd.DataFrame:
    """End-to-end cis instrument selection short of LD clumping.

    LD clumping must be performed downstream (R / PLINK). This function returns
    every cis-significant SNP that also passes the F-statistic threshold; the
    R clumping step then prunes correlated variants.
    """
    window = make_cis_window(gene_chr, gene_start, gene_end, window_kb)
    cis = filter_to_cis(qtl_df, window)
    sig = filter_by_pvalue(cis, p_threshold)
    sig = annotate_f_statistic(sig)
    return filter_by_fstat(sig, min_f_stat)
