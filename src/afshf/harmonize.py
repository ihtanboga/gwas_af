"""GWAS / QTL summary statistics harmonization (spec section 2).

Implements the standard schema and palindromic / MAF / build-aware QC rules.
The functions are pure (DataFrame in, DataFrame out) so they can be unit
tested with synthetic inputs.

Standard schema columns (in order):

    trait, source, variant_id, rsid, chr, pos, build,
    effect_allele, other_allele, beta, se, or, pval, eaf,
    n, n_cases, n_controls, ancestry
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import pandas as pd


STANDARD_COLUMNS: tuple[str, ...] = (
    "trait", "source", "variant_id", "rsid", "chr", "pos", "build",
    "effect_allele", "other_allele", "beta", "se", "or", "pval", "eaf",
    "n", "n_cases", "n_controls", "ancestry",
)

VALID_ALLELES = frozenset({"A", "C", "G", "T"})
PALINDROMIC_PAIRS = frozenset({("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")})


@dataclass(frozen=True)
class HarmonizationReport:
    n_in: int
    n_out: int
    dropped_missing_betas: int = 0
    dropped_invalid_alleles: int = 0
    dropped_low_info: int = 0
    dropped_low_maf: int = 0
    dropped_palindromic_unresolved: int = 0


def _coerce_str_alleles(s: pd.Series) -> pd.Series:
    return s.astype("string").str.upper().str.strip()


def is_valid_allele_pair(ea: str, oa: str) -> bool:
    return ea in VALID_ALLELES and oa in VALID_ALLELES and ea != oa


def is_palindromic(ea: str, oa: str) -> bool:
    return (ea, oa) in PALINDROMIC_PAIRS


def palindromic_decision(ea: str,
                         oa: str,
                         eaf: float | None,
                         drop_if_eaf_missing: bool = True,
                         eaf_ambiguous_window: float = 0.05) -> str:
    """Return one of {'keep', 'drop_missing_eaf', 'drop_ambiguous'}."""
    if not is_palindromic(ea, oa):
        return "keep"
    if eaf is None or pd.isna(eaf):
        return "drop_missing_eaf" if drop_if_eaf_missing else "keep"
    if abs(float(eaf) - 0.5) < eaf_ambiguous_window:
        return "drop_ambiguous"
    return "keep"


def harmonize_sumstats(df: pd.DataFrame,
                       *,
                       drop_info_below: float = 0.8,
                       drop_maf_below: float = 0.01,
                       palindromic_drop_if_missing_eaf: bool = True,
                       palindromic_eaf_window: float = 0.05,
                       keep_low_freq_separately: bool = False
                       ) -> tuple[pd.DataFrame, HarmonizationReport]:
    """Apply the QC pipeline described in spec section 2.

    The function never mutates the input DataFrame.
    """

    work = df.copy()
    n_in = len(work)

    if "effect_allele" in work.columns:
        work["effect_allele"] = _coerce_str_alleles(work["effect_allele"])
    if "other_allele" in work.columns:
        work["other_allele"] = _coerce_str_alleles(work["other_allele"])

    drop_missing = work["beta"].isna() | work["se"].isna() | work["pval"].isna()
    n_missing = int(drop_missing.sum())
    work = work.loc[~drop_missing]

    valid_alleles = work.apply(
        lambda r: is_valid_allele_pair(r["effect_allele"], r["other_allele"]), axis=1
    )
    n_invalid = int((~valid_alleles).sum())
    work = work.loc[valid_alleles]

    n_low_info = 0
    if "info" in work.columns:
        info_ok = work["info"].fillna(1.0) >= drop_info_below
        n_low_info = int((~info_ok).sum())
        work = work.loc[info_ok]

    n_low_maf = 0
    low_freq_subset = None
    if "eaf" in work.columns:
        eaf = work["eaf"].astype(float)
        maf = eaf.where(eaf <= 0.5, 1.0 - eaf)
        # Only enforce the MAF filter when EAF is observed; missing EAF is
        # routed through the palindromic step further down.
        eaf_observed = ~eaf.isna()
        maf_ok = (~eaf_observed) | (maf >= drop_maf_below)
        n_low_maf = int((~maf_ok).sum())
        if keep_low_freq_separately:
            low_freq_subset = work.loc[~maf_ok].copy()
        work = work.loc[maf_ok]

    pal_decisions = work.apply(
        lambda r: palindromic_decision(
            r["effect_allele"],
            r["other_allele"],
            r.get("eaf"),
            drop_if_eaf_missing=palindromic_drop_if_missing_eaf,
            eaf_ambiguous_window=palindromic_eaf_window,
        ),
        axis=1,
    )
    keep = pal_decisions == "keep"
    n_pal = int((~keep).sum())
    work = work.loc[keep]

    work = work.reset_index(drop=True)
    report = HarmonizationReport(
        n_in=n_in,
        n_out=len(work),
        dropped_missing_betas=n_missing,
        dropped_invalid_alleles=n_invalid,
        dropped_low_info=n_low_info,
        dropped_low_maf=n_low_maf,
        dropped_palindromic_unresolved=n_pal,
    )
    if keep_low_freq_separately and low_freq_subset is not None:
        # Attach as attribute for downstream rules to pick up.
        work.attrs["low_freq_variants"] = low_freq_subset
    return work, report


def ensure_standard_columns(df: pd.DataFrame,
                            extra_required: Iterable[str] = ()) -> None:
    missing = [c for c in STANDARD_COLUMNS if c not in df.columns]
    extra = [c for c in extra_required if c not in df.columns]
    if missing or extra:
        raif then ValueError(
            f"Input is missing required columns: standard={missing}, extra={extra}"
        )
