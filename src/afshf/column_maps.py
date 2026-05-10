"""Source-specific column-alias resolver and effect-scale derivations.

Loads ``config/column_maps.yml`` and provides a single ``apply_column_map``
function that turns a raw GWAS/QTL DataFrame into the project's standard
schema (``afshf.harmonize.STANDARD_COLUMNS``) prior to QC.

Supported derivations:

* ``beta_from_or``  : if ``or`` is present and ``beta`` is missing,
  set ``beta = log(or)``.
* ``beta_from_z``   : eQTLGen-style ``z + N + EAF`` -> ``beta + se`` using
  the standard z-to-beta approximation:

      var_g = 2 * MAF * (1 - MAF) * (N + z**2)
      se    = 1 / sqrt(var_g)
      beta  = z * se

* ``beta_from_log10p`` : UKB-PPP ships ``LOG10P``; we recover
  ``p = 10 ** (-log10p)`` so downstream FDR / coloc rules work.
* ``se_from_or``    : when ``se`` is reported on the OR scale (rare),
  divide by OR to obtain the log-OR-scale SE.

The function never mutates the input DataFrame.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import log, sqrt
from pathlib import Path
from typing import Iterable, Mapping, Optional

import numpy as np
import pandas as pd
import yaml


STANDARD_FIELDS: tuple[str, ...] = (
    "variant_id", "rsid", "chr", "pos",
    "effect_allele", "other_allele",
    "beta", "or", "se", "pval", "log10p", "log_p", "eaf",
    "n", "n_cases", "n_controls", "info", "z",
    # D026 — underflow-safe significance columns:
    "log10p_signed", "neg_log10p", "pval_capped_for_tools", "pval_underflow_flag",
)


# D026 — pval=0 floor for downstream R/TwoSampleMR/coloc compatibility.
# Floor chosen so the smallest representable double (~1e-308) never appears.
PVAL_CAPPED_FLOOR = 1e-300


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ColumnMap:
    name: str
    description: str
    column_aliases: dict[str, list[str]]
    derive: dict[str, bool]
    metadata: dict[str, str]


@dataclass(frozen=True)
class ColumnMapRegistry:
    defaults: dict
    sources: dict[str, ColumnMap]

    def get(self, source: str) -> ColumnMap:
        if source not in self.sources:
            raif then KeyError(
                f"Unknown source '{source}'. Available: {sorted(self.sources)}"
            )
        return self.sources[source]


def load_column_maps(path: str | Path) -> ColumnMapRegistry:
    with open(path, "r", encoding="utf-8") as f:
        raw = yaml.safe_load(f)
    sources = {}
    for key, body in raw.get("sources", {}).items():
        sources[key] = ColumnMap(
            name=key,
            description=str(body.get("description", "")),
            column_aliases={k: list(v or []) for k, v in body.get("column_aliases", {}).items()},
            derive=dict(body.get("derive", {}) or {}),
            metadata=dict(body.get("metadata", {}) or {}),
        )
    return ColumnMapRegistry(defaults=raw.get("defaults", {}) or {}, sources=sources)


# ---------------------------------------------------------------------------
# Alias resolution
# ---------------------------------------------------------------------------


def _resolve_alias(columns: Iterable[str],
                   candidates: Iterable[str],
                   case_sensitive_first: bool = True) -> Optional[str]:
    cols = list(columns)
    cols_lower = {c.lower(): c for c in cols}
    for cand in candidates:
        if not cand:
            continue
        if case_sensitive_first and cand in cols:
            return cand
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    return None


def resolve_columns(df: pd.DataFrame,
                    column_map: ColumnMap,
                    case_sensitive_first: bool = True) -> dict[str, Optional[str]]:
    """Return a {standard_field: source_column_or_None} dict for the DataFrame."""
    resolved: dict[str, Optional[str]] = {}
    for std_field in STANDARD_FIELDS:
        candidates = column_map.column_aliases.get(std_field, []) or []
        resolved[std_field] = _resolve_alias(df.columns, candidates, case_sensitive_first)
    return resolved


# ---------------------------------------------------------------------------
# Effect-scale derivations
# ---------------------------------------------------------------------------


def derive_beta_from_or(df: pd.DataFrame,
                        beta_col: str = "beta",
                        or_col: str = "or") -> pd.DataFrame:
    """Set ``beta = log(or)`` for rows where beta is missing but OR is present."""
    out = df.copy()
    has_or = or_col in out.columns
    if not has_or:
        return out
    if beta_col not in out.columns:
        out[beta_col] = np.nan
    or_vals = pd.to_numeric(out[or_col], errors="coerce")
    needs = out[beta_col].isna() & or_vals.gt(0.0)
    if needs.any():
        out.loc[needs, beta_col] = np.log(or_vals.loc[needs])
    return out


def derive_or_from_beta(df: pd.DataFrame,
                        beta_col: str = "beta",
                        or_col: str = "or") -> pd.DataFrame:
    """Set ``or = exp(beta)`` when missing — useful for downstream reporting."""
    out = df.copy()
    if beta_col not in out.columns:
        return out
    if or_col not in out.columns:
        out[or_col] = np.nan
    needs = out[or_col].isna() & out[beta_col].notna()
    if needs.any():
        out.loc[needs, or_col] = np.exp(out.loc[needs, beta_col])
    return out


def derive_se_from_or_scale(df: pd.DataFrame,
                            se_col: str = "se",
                            or_col: str = "or") -> pd.DataFrame:
    """Convert SE from OR-scale to log-OR-scale: ``se_log_or = se_or / or``."""
    out = df.copy()
    if se_col not in out.columns or or_col not in out.columns:
        return out
    or_vals = pd.to_numeric(out[or_col], errors="coerce")
    se_vals = pd.to_numeric(out[se_col], errors="coerce")
    valid = or_vals.gt(0.0) & se_vals.notna()
    if valid.any():
        out.loc[valid, se_col] = se_vals.loc[valid] / or_vals.loc[valid]
    return out


def derive_beta_from_z(df: pd.DataFrame,
                       z_col: str = "z",
                       eaf_col: str = "eaf",
                       n_col: str = "n",
                       beta_col: str = "beta",
                       se_col: str = "se") -> pd.DataFrame:
    """eQTLGen-style Z + EAF + N -> beta + SE.

    Uses the standard approximation:
        var_g = 2 * MAF * (1 - MAF) * (N + z**2)
        se    = 1 / sqrt(var_g)
        beta  = z * se
    """
    out = df.copy()
    needed = {z_col, eaf_col, n_col}
    if not needed.issubset(out.columns):
        return out

    z = pd.to_numeric(out[z_col], errors="coerce")
    eaf = pd.to_numeric(out[eaf_col], errors="coerce")
    n = pd.to_numeric(out[n_col], errors="coerce")
    maf = eaf.where(eaf <= 0.5, 1.0 - eaf)

    valid = z.notna() & maf.gt(0.0) & maf.lt(0.5) & n.gt(0.0)
    var_g = 2.0 * maf * (1.0 - maf) * (n + z ** 2)
    se = 1.0 / np.sqrt(var_g.where(valid))
    beta = z * se

    if beta_col not in out.columns:
        out[beta_col] = np.nan
    if se_col not in out.columns:
        out[se_col] = np.nan

    needs_beta = out[beta_col].isna() & valid
    needs_se = out[se_col].isna() & valid
    out.loc[needs_beta, beta_col] = beta.loc[needs_beta]
    out.loc[needs_se, se_col] = se.loc[needs_se]
    return out


def derive_pval_from_log10p(df: pd.DataFrame,
                            log10p_col: str = "log10p",
                            pval_col: str = "pval") -> pd.DataFrame:
    """UKB-PPP releases ship ``LOG10P``; recover ``p = 10**(-log10p)``."""
    out = df.copy()
    if log10p_col not in out.columns:
        return out
    log10p = pd.to_numeric(out[log10p_col], errors="coerce")
    if pval_col not in out.columns:
        out[pval_col] = np.nan
    needs = out[pval_col].isna() & log10p.notna()
    if needs.any():
        out.loc[needs, pval_col] = np.power(10.0, -log10p.loc[needs])
    return out


def derive_log10p_columns(df: pd.DataFrame,
                          *,
                          source_signed_col: str | None = "log_p",
                          source_unsigned_col: str | None = "log10p",
                          base_is_natural: bool = False,
                          floor: float = PVAL_CAPPED_FLOOR
                          ) -> pd.DataFrame:
    """D026 — derive underflow-safe significance columns.

    Adds (when at least one source column is present):

    * ``log10p_signed``           : signed log10(P) (negative for significant).
    * ``neg_log10p``              : -log10(P) (positive for significant).
    * ``pval_capped_for_tools``   : ``max(10**log10p_signed, floor)``.
    * ``pval_underflow_flag``     : True when raw 10**log10p_signed < floor.

    Source priority:
        1. ``source_signed_col`` (e.g. Roselli 2025 ``log_p``) — already signed.
        2. ``source_unsigned_col`` (e.g. UKB-PPP ``log10p``) — flipped to signed.

    If ``base_is_natural`` then the source signed column is divided by ln(10)
    to convert to log10. Unsigned UKB-PPP-style log10p is always base-10.
    """
    out = df.copy()

    log10p_signed: pd.Series | None = None
    if source_signed_col and source_signed_col in out.columns:
        s = pd.to_numeric(out[source_signed_col], errors="coerce")
        if base_is_natural:
            s = s / np.log(10.0)
        log10p_signed = s
    elif source_unsigned_col and source_unsigned_col in out.columns:
        u = pd.to_numeric(out[source_unsigned_col], errors="coerce")
        log10p_signed = -u
    else:
        return out  # nothing to do

    out["log10p_signed"] = log10p_signed
    out["neg_log10p"] = -log10p_signed

    # Cap below the floor: any log10p_signed < log10(floor) is set to log10(floor).
    log_floor = float(np.log10(floor))           # e.g. log10(1e-300) = -300
    capped = np.clip(log10p_signed, log_floor, None)
    out["pval_capped_for_tools"] = np.power(10.0, capped)
    out["pval_underflow_flag"] = (log10p_signed < log_floor).fillna(False)
    return out


def derive_pval_from_log_p(df: pd.DataFrame,
                           log_p_col: str = "log_p",
                           pval_col: str = "pval",
                           base: float | str = 10,
                           sign: str = "auto") -> pd.DataFrame:
    """Recover P from ``log(P)``. Supports natural log and base-10, signed or
    unsigned logs.

    Parameters
    ----------
    base : 10 or "e"
        log base used by the source.
    sign : "auto" | "signed" | "unsigned"
        - "signed" : ``log_p`` is negative for significant variants
                     -> ``pval = base ** log_p``.
        - "unsigned": ``log_p`` is positive for significant variants
                     -> ``pval = base ** (-log_p)``.
        - "auto"   : pick "signed" if any value is negative, else "unsigned".

    Roselli 2025 ships "log(P)" — base + sign are inferred empirically.
    """
    out = df.copy()
    if log_p_col not in out.columns:
        return out
    log_p = pd.to_numeric(out[log_p_col], errors="coerce")
    if pval_col not in out.columns:
        out[pval_col] = np.nan
    needs = out[pval_col].isna() & log_p.notna()
    if not needs.any():
        return out

    if sign == "auto":
        sign = "signed" if (log_p.dropna() < 0).any() else "unsigned"

    base_val = np.e if (isinstance(base, str) and base.lower() == "e") else float(base)
    exponent = log_p if sign == "signed" else -log_p
    out.loc[needs, pval_col] = np.power(base_val, exponent.loc[needs])
    return out


# ---------------------------------------------------------------------------
# Top-level mapper
# ---------------------------------------------------------------------------


def apply_column_map(df: pd.DataFrame,
                     column_map: ColumnMap,
                     *,
                     trait: Optional[str] = None,
                     source: Optional[str] = None,
                     overrides: Mapping[str, str] | None = None,
                     case_sensitive_first: bool = True
                     ) -> pd.DataFrame:
    """Rename columns, derive missing scales, attach metadata.

    Returns a new DataFrame with at least the project's standard schema columns
    (extras pass through unchanged).
    """
    overrides = dict(overrides or {})

    resolved = resolve_columns(df, column_map, case_sensitive_first=case_sensitive_first)

    rename_map: dict[str, str] = {}
    for std_field, src_col in resolved.items():
        if src_col is None:
            continue
        if std_field in overrides:
            continue
        if src_col != std_field:
            rename_map[src_col] = std_field

    out = df.rename(columns=rename_map).copy()

    # Apply caller-side overrides (rare).
    for std_field, src_col in overrides.items():
        if src_col in out.columns and src_col != std_field:
            out[std_field] = out[src_col]

    derive = column_map.derive
    if "log10p" in out.columns:
        # log10p is now a first-class standard field; recover pval if missing.
        out = derive_pval_from_log10p(out)

    if derive.get("pval_from_log_p", False) and "log_p" in out.columns:
        out = derive_pval_from_log_p(
            out,
            base=derive.get("log_p_base", 10),
            sign=derive.get("log_p_sign", "auto"),
        )

    # D026 — always populate underflow-safe significance columns when either
    # signed (`log_p`) or unsigned (`log10p`) source is available.
    if "log_p" in out.columns or "log10p" in out.columns:
        base_is_natural = (str(derive.get("log_p_base", 10)).lower() == "e")
        out = derive_log10p_columns(out, base_is_natural=base_is_natural)

    if derive.get("beta_from_or", False):
        out = derive_beta_from_or(out)

    if derive.get("beta_from_z", False):
        out = derive_beta_from_z(out)

    if derive.get("se_from_or", False):
        out = derive_se_from_or_scale(out)

    out = derive_or_from_beta(out)

    metadata = {**column_map.metadata, **(overrides or {})}
    if "build" in metadata and "build" not in out.columns:
        out["build"] = metadata["build"]
    if "ancestry" in metadata and "ancestry" not in out.columns:
        out["ancestry"] = metadata["ancestry"]
    if trait and "trait" not in out.columns:
        out["trait"] = trait
    if source and "source" not in out.columns:
        out["source"] = source

    return out


# ---------------------------------------------------------------------------
# Convenience: full pipeline (alias map -> harmonize)
# ---------------------------------------------------------------------------


def to_standard_schema(df: pd.DataFrame,
                       registry: ColumnMapRegistry,
                       source: str,
                       *,
                       trait: Optional[str] = None,
                       overrides: Mapping[str, str] | None = None
                       ) -> pd.DataFrame:
    column_map = registry.get(source)
    case_sensitive_first = bool(registry.defaults.get("case_sensitive_first", True))
    return apply_column_map(
        df,
        column_map,
        trait=trait,
        source=source,
        overrides=overrides,
        case_sensitive_first=case_sensitive_first,
    )
