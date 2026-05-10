#!/usr/bin/env python3
"""Phase 3B + 3C smoke MR + projection (deCODE 2021 × AF + HF + Stroke).

per protocol decision 2026-05-09 (D030 + D032):
  3B: deCODE pQTL → AF (anchor) MR — aptamer-level, never collapse genes yet.
  3C: AF-anchored projection onto HF (3 primary + HFpEF exploratory) +
      stroke (4 primary). Theta = d_T * beta_T_outcome.

Strict smoke constraints:
  - All result rows tagged SMOKE_ONLY = TRUE.
  - No GREEN/RED/AMBER classification.
  - No biological claims from 20-protein run.
  - HFpEF exploratory; never drives primary class logic (D001 REVISED).
  - Pval thresholds use neg_log10p (D026); pval_capped_for_tools only as
    R/coloc placeholder.
  - Multi-aptamer rows kept separate (D032).

Build handling:
  AF Roselli 2025         GRCh38   (matches deCODE)
  GIGASTROKE 2022 EUR     GRCh38   (matches deCODE)
  HERMES 2024 EUR         GRCh37   --> liftover deCODE instruments to b37
                                       per-outcome via pyliftover (chain hg38ToHg19).

Outputs:
  results/qc/decode_20protein_to_AF_harmonization_qc.tsv
  results/smoke/decode_20protein_to_AF_MR.tsv
  results/smoke/decode_20protein_AF_direction_table.tsv
  results/smoke/decode_20protein_hf_stroke_projection.tsv
  results/smoke/decode_20protein_direction_table.tsv
  results/smoke/decode_20protein_mr_summary.md
"""

from __future__ import annotations

import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import polars as pl

ROOT = Path("/Users/apple/Desktop/gwas_af")
INSTR_PARQUET = ROOT / "data/processed/instruments/decode_2021_smoke.parquet"

OUTCOME_PARQUETS: dict[str, dict] = {
    # outcome_label: (parquet, build, role)
    "AF":               {"path": ROOT / "data/processed/sumstats/AF.parquet",                          "build": "GRCh38", "role": "anchor"},
    "HF_overall":       {"path": ROOT / "data/processed/sumstats/HF_2024_EUR_HF_overall.parquet",     "build": "GRCh37", "role": "primary"},
    "HF_nonischemic":   {"path": ROOT / "data/processed/sumstats/HF_2024_EUR_HF_nonischemic.parquet", "build": "GRCh37", "role": "primary"},
    "HF_ni_HFrEF":      {"path": ROOT / "data/processed/sumstats/HF_2024_EUR_HF_ni_HFrEF.parquet",    "build": "GRCh37", "role": "primary"},
    "HF_ni_HFpEF":      {"path": ROOT / "data/processed/sumstats/HF_2024_EUR_HF_ni_HFpEF.parquet",    "build": "GRCh37", "role": "exploratory"},
    "ischemic_stroke":      {"path": ROOT / "data/processed/sumstats/STROKE_AIS_EUR.parquet", "build": "GRCh38", "role": "primary"},
    "cardioembolic_stroke": {"path": ROOT / "data/processed/sumstats/STROKE_CES_EUR.parquet", "build": "GRCh38", "role": "primary"},
    "large_artery_stroke":  {"path": ROOT / "data/processed/sumstats/STROKE_LAS_EUR.parquet", "build": "GRCh38", "role": "primary"},
    "small_vessel_stroke":  {"path": ROOT / "data/processed/sumstats/STROKE_SVS_EUR.parquet", "build": "GRCh38", "role": "primary"},
}

CHAIN_PATH = ROOT / "data/raw/annotations/hg38ToHg19.over.chain.gz"

QC_AF        = ROOT / "results/qc/decode_20protein_to_AF_harmonization_qc.tsv"
MR_AF        = ROOT / "results/smoke/decode_20protein_to_AF_MR.tsv"
DIR_AF       = ROOT / "results/smoke/decode_20protein_AF_direction_table.tsv"
PROJ_TBL     = ROOT / "results/smoke/decode_20protein_hf_stroke_projection.tsv"
DIR_TBL      = ROOT / "results/smoke/decode_20protein_direction_table.tsv"
SUMMARY_MD   = ROOT / "results/smoke/decode_20protein_mr_summary.md"

PALINDROMIC = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
PALINDROMIC_AMBIG_WINDOW = 0.05


def _log(msg: str) -> None:
    print(f"[3b-3c {time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Liftover for HF (b38 -> b37) — applies only to instrument coordinates
# ---------------------------------------------------------------------------


def liftover_instruments_b37(instr: pl.DataFrame) -> pl.DataFrame:
    from pyliftover import LiftOver
    lo = LiftOver(str(CHAIN_PATH))
    pos_b37 = []
    ok = []
    for chrom, pos in zip(instr["chrom"].to_list(), instr["pos"].to_list()):
        r = lo.convert_coordinate("chr" + str(chrom), int(pos))
        if r:
            pos_b37.append(r[0][1])
            ok.append(True)
        else:
            pos_b37.append(None)
            ok.append(False)
    return instr.with_columns(
        pl.Series("pos_b37", pos_b37, dtype=pl.Int64),
        pl.Series("liftover_ok", ok, dtype=pl.Boolean),
    )


# ---------------------------------------------------------------------------
# Per-instrument harmonization against an outcome parquet
# ---------------------------------------------------------------------------


@dataclass
class HarmRow:
    aptamer: str
    gene: str
    rsid: str
    chrom: str
    pos: int
    pos_outcome: int
    ea_exposure: str
    oa_exposure: str
    ea_outcome: str
    oa_outcome: str
    beta_exposure: float
    se_exposure: float
    pval_exposure: float
    beta_outcome: float
    se_outcome: float
    pval_outcome: float
    neg_log10p_outcome: float | None
    eaf_outcome: float | None
    sign_flip: bool
    drop_reason: str            # "" if kept


def _is_palindromic(ea: str, oa: str) -> bool:
    return (ea, oa) in PALINDROMIC


def harmonize_for_outcome(instruments: pl.DataFrame,
                          outcome_df: pl.DataFrame,
                          *, build: str) -> tuple[list[HarmRow], dict]:
    """Match each instrument row to outcome by (chr, pos) — pos column for the
    relevant build. Apply allele harmonization + palindromic drop. Returns the
    list of kept rows + harmonization counters.
    """
    pos_col = "pos_b37" if build == "GRCh37" else "pos"
    counters = {
        "n_instruments_input": instruments.height,
        "n_instruments_matched": 0,
        "n_instruments_dropped_allele_mismatch": 0,
        "n_instruments_dropped_palindromic": 0,
        "n_instruments_missing_in_outcome": 0,
        "n_instruments_liftover_fawithd": 0,
    }

    # Drop rows that fawithd liftover (b37 path).
    if build == "GRCh37":
        before = instruments.height
        instruments = instruments.filter(pl.col("liftover_ok"))
        counters["n_instruments_liftover_fawithd"] = before - instruments.height

    if instruments.height == 0:
        return [], counters

    # Build lookup keyset (chr, pos)
    lookup_pos = instruments[pos_col].to_list()
    lookup_chr = [str(c) for c in instruments["chrom"].to_list()]

    # Filter outcome to candidate keys (efficient anti-cross-product join)
    keys = pl.DataFrame({
        "chr": lookup_chr,
        "pos": [int(p) for p in lookup_pos],
    })
    outcome_subset = outcome_df.join(keys, on=["chr", "pos"], how="inner")

    # Build outcome lookup map (chr,pos) -> first match (deCODE GIGASTROKE
    # files include only one row per variant; AF/HF rare duplicates will pick first).
    by_key = {}
    for r in outcome_subset.iter_rows(named=True):
        k = (str(r["chr"]), int(r["pos"]))
        if k not in by_key:
            by_key[k] = r

    kept: list[HarmRow] = []
    for r in instruments.iter_rows(named=True):
        key = (str(r["chrom"]), int(r[pos_col]))
        oc = by_key.get(key)
        if oc is None:
            counters["n_instruments_missing_in_outcome"] += 1
            kept.append(HarmRow(
                aptamer=r["aptamer"], gene=r["gene"], rsid=r.get("rsid", "") or "",
                chrom=str(r["chrom"]), pos=int(r["pos"]),
                pos_outcome=int(r[pos_col]),
                ea_exposure=r["effect_allele"], oa_exposure=r["other_allele"],
                ea_outcome="", oa_outcome="",
                beta_exposure=float(r["beta"]), se_exposure=float(r["se"]),
                pval_exposure=float(r["pval"]),
                beta_outcome=math.nan, se_outcome=math.nan, pval_outcome=math.nan,
                neg_log10p_outcome=None, eaf_outcome=None,
                sign_flip=False, drop_reason="missing_in_outcome",
            ))
            continue
        ea_e = (r["effect_allele"] or "").upper()
        oa_e = (r["other_allele"] or "").upper()
        ea_o = (oc.get("effect_allele") or "").upper()
        oa_o = (oc.get("other_allele") or "").upper()

        # Allele compatibility: same alleles, possibly flipped.
        sign_flip = False
        if {ea_e, oa_e} != {ea_o, oa_o}:
            counters["n_instruments_dropped_allele_mismatch"] += 1
            kept.append(HarmRow(
                aptamer=r["aptamer"], gene=r["gene"], rsid=r.get("rsid", "") or "",
                chrom=str(r["chrom"]), pos=int(r["pos"]),
                pos_outcome=int(r[pos_col]),
                ea_exposure=ea_e, oa_exposure=oa_e,
                ea_outcome=ea_o, oa_outcome=oa_o,
                beta_exposure=float(r["beta"]), se_exposure=float(r["se"]),
                pval_exposure=float(r["pval"]),
                beta_outcome=float(oc["beta"]), se_outcome=float(oc["se"]),
                pval_outcome=float(oc["pval"]),
                neg_log10p_outcome=oc.get("neg_log10p"),
                eaf_outcome=oc.get("eaf"),
                sign_flip=False, drop_reason="allele_mismatch",
            ))
            continue
        if ea_e == oa_o and oa_e == ea_o:
            sign_flip = True

        # Palindromic check
        if _is_palindromic(ea_e, oa_e):
            eaf_o = oc.get("eaf")
            if eaf_o is None or math.isnan(float(eaf_o)) or abs(float(eaf_o) - 0.5) < PALINDROMIC_AMBIG_WINDOW:
                counters["n_instruments_dropped_palindromic"] += 1
                kept.append(HarmRow(
                    aptamer=r["aptamer"], gene=r["gene"], rsid=r.get("rsid", "") or "",
                    chrom=str(r["chrom"]), pos=int(r["pos"]),
                    pos_outcome=int(r[pos_col]),
                    ea_exposure=ea_e, oa_exposure=oa_e,
                    ea_outcome=ea_o, oa_outcome=oa_o,
                    beta_exposure=float(r["beta"]), se_exposure=float(r["se"]),
                    pval_exposure=float(r["pval"]),
                    beta_outcome=float(oc["beta"]), se_outcome=float(oc["se"]),
                    pval_outcome=float(oc["pval"]),
                    neg_log10p_outcome=oc.get("neg_log10p"),
                    eaf_outcome=eaf_o,
                    sign_flip=sign_flip, drop_reason="palindromic_unresolved",
                ))
                continue

        beta_o = float(oc["beta"]) * (-1.0 if sign_flip else 1.0)
        kept.append(HarmRow(
            aptamer=r["aptamer"], gene=r["gene"], rsid=r.get("rsid", "") or "",
            chrom=str(r["chrom"]), pos=int(r["pos"]),
            pos_outcome=int(r[pos_col]),
            ea_exposure=ea_e, oa_exposure=oa_e,
            ea_outcome=ea_o, oa_outcome=oa_o,
            beta_exposure=float(r["beta"]), se_exposure=float(r["se"]),
            pval_exposure=float(r["pval"]),
            beta_outcome=beta_o, se_outcome=float(oc["se"]),
            pval_outcome=float(oc["pval"]),
            neg_log10p_outcome=oc.get("neg_log10p"),
            eaf_outcome=oc.get("eaf"),
            sign_flip=sign_flip, drop_reason="",
        ))
        counters["n_instruments_matched"] += 1
    return kept, counters


# ---------------------------------------------------------------------------
# MR estimators
# ---------------------------------------------------------------------------


def mr_wald(beta_e: float, se_e: float,
            beta_o: float, se_o: float) -> tuple[float, float, float]:
    if beta_e == 0 or se_e <= 0:
        return math.nan, math.nan, math.nan
    beta_mr = beta_o / beta_e
    se_mr = abs(se_o / beta_e)
    z = beta_mr / se_mr if se_mr > 0 else math.nan
    p = 2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(z) / math.sqrt(2.0)))) if math.isfinite(z) else math.nan
    return beta_mr, se_mr, p


def mr_ivw(rows: list[HarmRow]) -> tuple[float, float, float, float, int, float]:
    """Inverse-variance weighted MR. Returns (beta, se, p, Q, n, Q_p)."""
    n = len(rows)
    if n == 0:
        return (math.nan,) * 4 + (0, math.nan)
    ws = []
    bxs = []
    bys = []
    for r in rows:
        if r.se_outcome <= 0:
            continue
        w = 1.0 / (r.se_outcome ** 2)
        ws.append(w)
        bxs.append(r.beta_exposure)
        bys.append(r.beta_outcome)
    if not ws:
        return (math.nan,) * 4 + (0, math.nan)
    num = sum(w * x * y for w, x, y in zip(ws, bxs, bys))
    den = sum(w * x * x for w, x in zip(ws, bxs))
    if den <= 0:
        return (math.nan,) * 4 + (0, math.nan)
    beta_ivw = num / den
    se_ivw = math.sqrt(1.0 / den)
    z = beta_ivw / se_ivw
    p = 2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(z) / math.sqrt(2.0))))

    # Heterogeneity Q
    Q = sum(w * (y - beta_ivw * x) ** 2 for w, x, y in zip(ws, bxs, bys))
    df = max(0, len(ws) - 1)
    if df > 0 and math.isfinite(Q):
        from math import lgamma, exp
        # Survival function of chi-square via regularised gamma — fall back to
        # a coarse approximation: report Q value, df, and let downstream stats
        # tools compute exact p. (No scipy in venv yet.)
        Q_p = math.nan
    else:
        Q_p = math.nan
    return beta_ivw, se_ivw, p, Q, len(ws), Q_p


def mr_weighted_median(rows: list[HarmRow]) -> tuple[float, float]:
    """Bowden et al. 2016 weighted median estimator (point estimate only;
    bootstrap SE deferred to full-scale R/TwoSampleMR run)."""
    n = len(rows)
    if n < 3:
        return math.nan, math.nan
    ratios = []
    weights = []
    for r in rows:
        if r.beta_exposure == 0 or r.se_outcome <= 0:
            continue
        ratios.append(r.beta_outcome / r.beta_exposure)
        weights.append((r.beta_exposure / r.se_outcome) ** 2)
    if len(ratios) < 3:
        return math.nan, math.nan
    order = sorted(range(len(ratios)), key=lambda i: ratios[i])
    ratios_s = [ratios[i] for i in order]
    weights_s = [weights[i] for i in order]
    total = sum(weights_s)
    cumulative = 0.0
    median = math.nan
    for r, w in zip(ratios_s, weights_s):
        cumulative += w
        if cumulative >= total / 2.0:
            median = r
            break
    return median, math.nan  # SE bootstrapping deferred


# ---------------------------------------------------------------------------
# Per-aptamer / per-outcome MR
# ---------------------------------------------------------------------------


def run_mr_per_aptamer(harm_rows: list[HarmRow]) -> dict:
    """Aggregate kept HarmRows for a single aptamer × outcome pair into MR result."""
    kept = [r for r in harm_rows if not r.drop_reason]
    n = len(kept)
    if n == 0:
        return {
            "n_instruments": 0, "method": "no_signal",
            "beta": math.nan, "se": math.nan, "p": math.nan,
            "neg_log10p": None, "Q": math.nan, "n_used": 0,
            "weighted_median_beta": math.nan,
        }
    if n == 1:
        b, s, p = mr_wald(kept[0].beta_exposure, kept[0].se_exposure,
                          kept[0].beta_outcome, kept[0].se_outcome)
        nlp = -math.log10(p) if (math.isfinite(p) and p > 0) else None
        return {
            "n_instruments": 1, "method": "wald_ratio",
            "beta": b, "se": s, "p": p, "neg_log10p": nlp,
            "Q": math.nan, "n_used": 1,
            "weighted_median_beta": math.nan,
        }
    beta_ivw, se_ivw, p_ivw, Q, n_used, _ = mr_ivw(kept)
    nlp = -math.log10(p_ivw) if (math.isfinite(p_ivw) and p_ivw > 0) else None
    wm_b, _ = mr_weighted_median(kept)
    return {
        "n_instruments": n, "method": "ivw_random",
        "beta": beta_ivw, "se": se_ivw, "p": p_ivw, "neg_log10p": nlp,
        "Q": Q, "n_used": n_used,
        "weighted_median_beta": wm_b,
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def main() -> int:
    QC_AF.parent.mkdir(parents=True, exist_ok=True)
    MR_AF.parent.mkdir(parents=True, exist_ok=True)

    instr = pl.read_parquet(INSTR_PARQUET)
    _log(f"loaded {instr.height} instruments / {instr['aptamer'].n_unique()} aptamers")

    # Liftover for HF outcomes
    instr_b37 = liftover_instruments_b37(instr)
    n_lift_fawithd = (~instr_b37["liftover_ok"]).sum()
    _log(f"liftover b38->b37 ok={instr_b37['liftover_ok'].sum()} fawithd={n_lift_fawithd}")

    # ----- Per-outcome harmonization + per-aptamer MR -----
    harmonization_qc_rows = []
    mr_rows = []          # one row per (aptamer, outcome)
    direction_rows_per_outcome = []
    af_dir = {}            # aptamer -> dict with d_T, beta_T_AF, etc.

    for outcome_label, meta in OUTCOME_PARQUETS.items():
        path = meta["path"]
        if not path.exists():
            _log(f"{outcome_label}: parquet missing; skip")
            continue
        outcome_df = pl.read_parquet(path)
        _log(f"{outcome_label} (build {meta['build']}, role {meta['role']}, "
             f"{outcome_df.height:,} rows)")

        instr_for_outcome = instr_b37 if meta["build"] == "GRCh37" else instr
        kept, counters = harmonize_for_outcome(
            instr_for_outcome, outcome_df, build=meta["build"]
        )
        counters["outcome"] = outcome_label
        counters["outcome_role"] = meta["role"]
        counters["outcome_build"] = meta["build"]
        harmonization_qc_rows.append(counters)

        # Group by aptamer
        by_aptamer: dict[str, list[HarmRow]] = {}
        for hr in kept:
            by_aptamer.setdefault(hr.aptamer, []).append(hr)

        for aptamer, hrs in by_aptamer.items():
            gene = hrs[0].gene
            mr = run_mr_per_aptamer(hrs)
            row = {
                "aptamer": aptamer, "gene": gene,
                "outcome": outcome_label, "outcome_role": meta["role"],
                "outcome_build": meta["build"],
                **mr,
                "smoke_only": True,
            }
            mr_rows.append(row)

            # Direction per outcome (anchor or projection)
            if outcome_label == "AF":
                if math.isfinite(row["beta"]) and row["beta"] != 0.0:
                    d_T = -1 if row["beta"] > 0 else 1
                else:
                    d_T = 0
                af_dir[aptamer] = {
                    "gene": gene,
                    "beta_T_AF": row["beta"], "se_T_AF": row["se"],
                    "p_T_AF": row["p"], "neg_log10p_T_AF": row["neg_log10p"],
                    "method_AF": row["method"], "n_instruments_AF": row["n_instruments"],
                    "d_T": d_T,
                }

    # ----- Build direction tables -----
    af_direction_rows = []
    for aptamer, d in af_dir.items():
        af_direction_rows.append({
            "aptamer": aptamer, "gene": d["gene"],
            "beta_T_AF": d["beta_T_AF"], "se_T_AF": d["se_T_AF"],
            "p_T_AF": d["p_T_AF"], "neg_log10p_T_AF": d["neg_log10p_T_AF"],
            "method_AF": d["method_AF"], "n_instruments_AF": d["n_instruments_AF"],
            "d_T": d["d_T"],
            "provisional_label": (
                "AF_MR_FDR_candidate" if d["p_T_AF"] is not None and math.isfinite(d["p_T_AF"]) and d["p_T_AF"] < 0.05 / max(1, len(af_dir)) else
                "AF_MR_nominal_candidate" if d["p_T_AF"] is not None and math.isfinite(d["p_T_AF"]) and d["p_T_AF"] < 0.05 else
                "AF_MR_no_signal"
            ),
            "smoke_only": True,
        })

    # Projection rows (HF + stroke) using per-aptamer d_T
    projection_rows = []
    for r in mr_rows:
        if r["outcome"] == "AF":
            continue
        d = af_dir.get(r["aptamer"])
        if d is None:
            continue
        d_T = d["d_T"]
        beta_o = r["beta"]
        theta = (d_T * beta_o) if (d_T != 0 and math.isfinite(beta_o)) else math.nan
        OR_intervention = math.exp(theta) if math.isfinite(theta) else math.nan
        # Provisional direction labels
        labels = []
        if d_T != 0 and math.isfinite(theta) and r["p"] is not None and math.isfinite(r["p"]) and r["p"] < 0.05:
            if theta < 0:
                labels.append("possible_shared_benefit_direction")
            else:
                if r["outcome"].startswith("HF"):
                    labels.append("possible_HF_harm_direction")
                elif "stroke" in r["outcome"]:
                    labels.append("possible_stroke_harm_direction")
        projection_rows.append({
            "aptamer": r["aptamer"], "gene": r["gene"],
            "outcome": r["outcome"], "outcome_role": r["outcome_role"],
            "outcome_build": r["outcome_build"],
            "d_T": d_T,
            "beta_T_outcome": beta_o, "se_T_outcome": r["se"], "p_T_outcome": r["p"],
            "neg_log10p_T_outcome": r["neg_log10p"],
            "theta": theta, "OR_intervention": OR_intervention,
            "method": r["method"], "n_instruments": r["n_instruments"],
            "weighted_median_beta": r.get("weighted_median_beta"),
            "Q_heterogeneity": r["Q"],
            "provisional_labels": "|".join(labels),
            "smoke_only": True,
        })

    # ----- Persist -----
    pl.DataFrame(harmonization_qc_rows).write_csv(QC_AF, separator="\t")
    pl.DataFrame(mr_rows).write_csv(MR_AF, separator="\t")
    pl.DataFrame(af_direction_rows).write_csv(DIR_AF, separator="\t")
    pl.DataFrame(projection_rows).write_csv(PROJ_TBL, separator="\t")
    pl.DataFrame(projection_rows).write_csv(DIR_TBL, separator="\t")

    # MD summary
    md = ["# deCODE 20-protein smoke MR + projection (Phase 3B + 3C, SMOKE_ONLY)\n",
          "All rows tagged `smoke_only = True`. **No final classification.**",
          "",
          "## AF anchor MR (Phase 3B)",
          "",
          "| Aptamer | Gene | n_instr | method | beta_T_AF | p | neg_log10p | d_T | label |",
          "|---|---|---:|---|---:|---:|---:|---:|---|",
    ]
    for r in af_direction_rows:
        if r['p_T_AF'] is not None and math.isfinite(r['p_T_AF']):
            nlp = f"{r['neg_log10p_T_AF']:.2f}" if r['neg_log10p_T_AF'] is not None else "—"
            md.append(
                f"| {r['aptamer'][:35]} | {r['gene']} | {r['n_instruments_AF']} | {r['method_AF']} | "
                f"{r['beta_T_AF']:.3f} | {r['p_T_AF']:.2e} | {nlp} | "
                f"{r['d_T']:+d} | {r['provisional_label']} |"
            )
        else:
            md.append(
                f"| {r['aptamer'][:35]} | {r['gene']} | 0 | no_signal | — | — | — | 0 | AF_MR_no_signal |"
            )

    md.append("\n## HF/stroke projection (Phase 3C, theta = d_T * beta_outcome)\n")
    md.append("Provisional direction labels only. HFpEF exploratory.\n")
    n_with_label = sum(1 for r in projection_rows if r["provisional_labels"])
    md.append(f"- aptamer × outcome rows: {len(projection_rows)}")
    md.append(f"- rows with provisional_label fired: {n_with_label}")
    SUMMARY_MD.write_text("\n".join(md))

    _log(f"wrote {MR_AF}")
    _log(f"wrote {DIR_AF}")
    _log(f"wrote {PROJ_TBL}")
    _log(f"wrote {SUMMARY_MD}")
    _log(f"smoke MR + projection complete: {len(af_direction_rows)} AF MR rows, "
         f"{len(projection_rows)} projection rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
