#!/usr/bin/env python3
"""D037 — CVD150 AF MR + candidate-only HF/stroke projection.

Reads `data/processed/instruments/decode_cvd150_instruments.parquet`
(combined CVD50 + extension), runs aptamer-level MR against AF, applies
BH-FDR within the CVD150 universe, builds candidate shortlist with
CVD150_AF_* labels, and projects HF/stroke for candidates only.

Coloc + final classification remain LOCKED (D037).
"""

from __future__ import annotations

import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import polars as pl

ROOT = Path("/Users/apple/Desktop/gwas_af")
INSTR_PARQUET = ROOT / "data/processed/instruments/decode_cvd150_instruments.parquet"
CHAIN_PATH = ROOT / "data/raw/annotations/hg38ToHg19.over.chain.gz"

OUTCOMES = {
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

OUT_AF_HARM = ROOT / "results/qc/decode_cvd150_to_AF_harmonization_qc.tsv"
OUT_AF_MR = ROOT / "results/mr/decode_cvd150_pqtl_to_AF.tsv"
OUT_CANDIDATES = ROOT / "results/mr/decode_cvd150_AF_candidate_shortlist.tsv"
OUT_PROJ = ROOT / "results/mr/decode_cvd150_candidate_hf_stroke_projection.tsv"
OUT_DIRECTION = ROOT / "results/mr/decode_cvd150_candidate_direction_table.tsv"
OUT_SUMMARY = ROOT / "results/mr/decode_cvd150_candidate_summary.md"

SANITY_TARGETS = {"IL6R", "NPPA", "NPPB", "PCSK9", "F11", "F2", "F9",
                  "ADAMTS13", "MMP12", "DSC2", "APOE", "F13B"}

PALINDROMIC = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
PALINDROMIC_AMBIG_WINDOW = 0.05


def _log(msg: str) -> None:
    print(f"[cvd150-mr {time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Liftover (deCODE b38 -> b37 for HERMES HF outcomes)
# ---------------------------------------------------------------------------


def liftover_b38_to_b37(instr: pl.DataFrame) -> pl.DataFrame:
    from pyliftover import LiftOver
    lo = LiftOver(str(CHAIN_PATH))
    pos_b37, ok = [], []
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
# Harmonization vs outcome (mirrors CVD454 MR script)
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
    beta_exposure: float
    se_exposure: float
    pval_exposure: float
    beta_outcome: float
    se_outcome: float
    pval_outcome: float
    neg_log10p_outcome: float | None
    eaf_outcome: float | None
    sign_flip: bool


def harmonize(instruments: pl.DataFrame, outcome_df: pl.DataFrame,
              build: str) -> tuple[list[HarmRow], dict]:
    pos_col = "pos_b37" if build == "GRCh37" else "pos"
    counters = {
        "n_instruments_input": instruments.height,
        "n_instruments_matched": 0,
        "n_instruments_dropped_allele_mismatch": 0,
        "n_instruments_dropped_palindromic": 0,
        "n_instruments_missing_in_outcome": 0,
        "n_instruments_liftover_fawithd": 0,
    }
    if build == "GRCh37":
        before = instruments.height
        instruments = instruments.filter(pl.col("liftover_ok"))
        counters["n_instruments_liftover_fawithd"] = before - instruments.height
    if instruments.height == 0:
        return [], counters

    keys = pl.DataFrame({
        "chr": [str(c) for c in instruments["chrom"].to_list()],
        "pos": [int(p) for p in instruments[pos_col].to_list()],
    })
    outcome_subset = outcome_df.join(keys, on=["chr", "pos"], how="inner")
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
            continue
        ea_e = (r["effect_allele"] or "").upper()
        oa_e = (r["other_allele"] or "").upper()
        ea_o = (oc.get("effect_allele") or "").upper()
        oa_o = (oc.get("other_allele") or "").upper()
        sign_flip = False
        if {ea_e, oa_e} != {ea_o, oa_o}:
            counters["n_instruments_dropped_allele_mismatch"] += 1
            continue
        if ea_e == oa_o and oa_e == ea_o:
            sign_flip = True
        if (ea_e, oa_e) in PALINDROMIC:
            eaf_o = oc.get("eaf")
            if eaf_o is None or eaf_o != eaf_o or abs(float(eaf_o) - 0.5) < PALINDROMIC_AMBIG_WINDOW:
                counters["n_instruments_dropped_palindromic"] += 1
                continue
        beta_o = float(oc["beta"]) * (-1.0 if sign_flip else 1.0)
        kept.append(HarmRow(
            aptamer=r["aptamer"], gene=r["gene"], rsid=r.get("rsid", "") or "",
            chrom=str(r["chrom"]), pos=int(r["pos"]),
            pos_outcome=int(r[pos_col]),
            ea_exposure=ea_e, oa_exposure=oa_e,
            beta_exposure=float(r["beta"]), se_exposure=float(r["se"]),
            pval_exposure=float(r["pval"]),
            beta_outcome=beta_o, se_outcome=float(oc["se"]),
            pval_outcome=float(oc["pval"]),
            neg_log10p_outcome=oc.get("neg_log10p"),
            eaf_outcome=oc.get("eaf"),
            sign_flip=sign_flip,
        ))
        counters["n_instruments_matched"] += 1
    return kept, counters


# ---------------------------------------------------------------------------
# MR estimators
# ---------------------------------------------------------------------------


def _normal_cdf(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def mr_wald(beta_e, se_e, beta_o, se_o):
    if beta_e == 0 or se_e <= 0 or se_o <= 0:
        return math.nan, math.nan, math.nan
    beta_mr = beta_o / beta_e
    se_mr = abs(se_o / beta_e)
    z = beta_mr / se_mr
    p = 2.0 * (1.0 - _normal_cdf(abs(z)))
    return beta_mr, se_mr, p


def mr_ivw(rows):
    if not rows:
        return {"beta": math.nan, "se": math.nan, "p": math.nan, "n": 0, "Q": math.nan}
    ws, bxs, bys = [], [], []
    for r in rows:
        if r.se_outcome <= 0:
            continue
        ws.append(1.0 / (r.se_outcome ** 2))
        bxs.append(r.beta_exposure); bys.append(r.beta_outcome)
    if not ws:
        return {"beta": math.nan, "se": math.nan, "p": math.nan, "n": 0, "Q": math.nan}
    num = sum(w * x * y for w, x, y in zip(ws, bxs, bys))
    den = sum(w * x * x for w, x in zip(ws, bxs))
    if den <= 0:
        return {"beta": math.nan, "se": math.nan, "p": math.nan, "n": len(ws), "Q": math.nan}
    beta = num / den
    se = math.sqrt(1.0 / den)
    z = beta / se
    p = 2.0 * (1.0 - _normal_cdf(abs(z)))
    Q = sum(w * (y - beta * x) ** 2 for w, x, y in zip(ws, bxs, bys))
    return {"beta": beta, "se": se, "p": p, "n": len(ws), "Q": Q}


def mr_weighted_median(rows):
    if len(rows) < 3:
        return math.nan
    ratios, weights = [], []
    for r in rows:
        if r.beta_exposure == 0 or r.se_outcome <= 0:
            continue
        ratios.append(r.beta_outcome / r.beta_exposure)
        weights.append((r.beta_exposure / r.se_outcome) ** 2)
    if len(ratios) < 3:
        return math.nan
    order = sorted(range(len(ratios)), key=lambda i: ratios[i])
    ratios_s = [ratios[i] for i in order]
    weights_s = [weights[i] for i in order]
    total = sum(weights_s)
    cumulative = 0.0
    for r, w in zip(ratios_s, weights_s):
        cumulative += w
        if cumulative >= total / 2.0:
            return r
    return math.nan


def run_mr(rows):
    if not rows:
        return {"n_instruments": 0, "method": "no_signal", "beta": math.nan,
                "se": math.nan, "p": math.nan, "Q": math.nan,
                "weighted_median_beta": math.nan, "neg_log10p": None}
    if len(rows) == 1:
        b, s, p = mr_wald(rows[0].beta_exposure, rows[0].se_exposure,
                          rows[0].beta_outcome, rows[0].se_outcome)
        return {"n_instruments": 1, "method": "wald_ratio",
                "beta": b, "se": s, "p": p, "Q": math.nan,
                "weighted_median_beta": math.nan,
                "neg_log10p": (-math.log10(p) if (math.isfinite(p) and p > 0) else None)}
    ivw = mr_ivw(rows)
    return {"n_instruments": len(rows), "method": "ivw_random",
            "beta": ivw["beta"], "se": ivw["se"], "p": ivw["p"],
            "Q": ivw["Q"], "weighted_median_beta": mr_weighted_median(rows),
            "neg_log10p": (-math.log10(ivw["p"]) if (math.isfinite(ivw["p"]) and ivw["p"] > 0) else None)}


def benjamini_hochberg(pvals):
    n = len(pvals)
    if n == 0:
        return []
    indexed = [(p, i) for i, p in enumerate(pvals)]
    indexed_sorted = sorted(indexed, key=lambda x: (math.inf if (x[0] is None or not math.isfinite(x[0])) else x[0]))
    q = [None] * n
    prev_q = math.inf
    for rank, (p, orig_i) in enumerate(reversed(indexed_sorted), start=1):
        k = n - rank + 1
        if p is None or not math.isfinite(p):
            q[orig_i] = math.nan
            continue
        adj = p * n / k
        prev_q = min(prev_q, adj)
        q[orig_i] = min(prev_q, 1.0)
    return q


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def main() -> int:
    if not INSTR_PARQUET.exists():
        _log(f"FATAL: instruments parquet missing: {INSTR_PARQUET}")
        return 2
    OUT_AF_HARM.parent.mkdir(parents=True, exist_ok=True)
    OUT_AF_MR.parent.mkdir(parents=True, exist_ok=True)

    instr = pl.read_parquet(INSTR_PARQUET)
    matched = instr.filter(pl.col("af_match") == "matched")
    _log(f"loaded {instr.height} instruments / {instr['aptamer'].n_unique()} aptamers; "
         f"matched-to-AF: {matched.height} / {matched['aptamer'].n_unique()}")

    instr_b37 = liftover_b38_to_b37(matched)

    af_df = pl.read_parquet(OUTCOMES["AF"]["path"]).select([
        "chr", "pos", "rsid", "effect_allele", "other_allele",
        "beta", "se", "pval", "neg_log10p", "eaf",
    ])
    _log(f"AF parquet ({af_df.height:,})")

    # ----- AF MR per aptamer -----
    mr_rows = []
    af_harm_qc = []
    by_apt = matched.group_by("aptamer").agg(pl.first("gene").alias("gene"))
    for r in by_apt.iter_rows(named=True):
        apt, gene = r["aptamer"], r["gene"]
        sub = matched.filter(pl.col("aptamer") == apt)
        kept, counters = harmonize(sub, af_df, build="GRCh38")
        counters["aptamer"] = apt; counters["gene"] = gene
        af_harm_qc.append(counters)
        m = run_mr(kept)
        mr_rows.append({
            "aptamer": apt, "gene": gene, **m,
            "scope": "CVD150_screen", "fdr_label": "CVD150_screen_FDR_NOT_proteome_wide",
        })
    pl.DataFrame(af_harm_qc).write_csv(OUT_AF_HARM, separator="\t")
    _log(f"wrote {OUT_AF_HARM}")

    # BH-FDR within CVD150 universe
    pvals = [r["p"] for r in mr_rows]
    fdrs = benjamini_hochberg(pvals)
    for r, q in zip(mr_rows, fdrs):
        r["fdr_bh_within_cvd150"] = q
    pl.DataFrame(mr_rows).write_csv(OUT_AF_MR, separator="\t")
    _log(f"wrote {OUT_AF_MR} (BH-FDR within CVD150 universe; NOT proteome-wide)")

    # ----- Candidate shortlist -----
    candidates = []
    for r in mr_rows:
        p = r["p"]; q = r["fdr_bh_within_cvd150"]
        if p is None or not math.isfinite(p):
            label = "CVD150_AF_no_signal"
        elif q is not None and math.isfinite(q) and q < 0.05:
            label = "CVD150_AF_FDR_candidate"
        elif p < 0.05:
            label = "CVD150_AF_nominal_candidate"
        else:
            label = "CVD150_AF_no_signal"
        sanity = r["gene"] in SANITY_TARGETS
        carry = (label in ("CVD150_AF_FDR_candidate", "CVD150_AF_nominal_candidate") or sanity)
        d_T = 0
        if r.get("beta") is not None and math.isfinite(r["beta"]) and r["beta"] != 0:
            d_T = -1 if r["beta"] > 0 else 1
        candidates.append({
            "aptamer": r["aptamer"], "gene": r["gene"],
            "method_AF": r["method"], "n_instruments_AF": r["n_instruments"],
            "beta_T_AF": r["beta"], "se_T_AF": r["se"],
            "p_T_AF": r["p"], "fdr_bh_within_cvd150": r["fdr_bh_within_cvd150"],
            "neg_log10p_T_AF": r["neg_log10p"],
            "weighted_median_beta_AF": r.get("weighted_median_beta"),
            "Q_heterogeneity_AF": r.get("Q"),
            "d_T": d_T,
            "candidate_label": label,
            "is_sanity_target": sanity,
            "carry_to_projection": carry,
            "scope": "CVD150_screen_NOT_proteome_wide",
            "fdr_label": "CVD150_screen_FDR_NOT_proteome_wide",
        })
    pl.DataFrame(candidates).write_csv(OUT_CANDIDATES, separator="\t")
    _log(f"wrote {OUT_CANDIDATES}")

    n_fdr = sum(1 for c in candidates if c["candidate_label"] == "CVD150_AF_FDR_candidate")
    n_nom = sum(1 for c in candidates if c["candidate_label"] == "CVD150_AF_nominal_candidate")
    n_no = sum(1 for c in candidates if c["candidate_label"] == "CVD150_AF_no_signal")
    n_carry = sum(1 for c in candidates if c["carry_to_projection"])
    _log(f"CVD150 AF MR: FDR={n_fdr} nominal={n_nom} no_signal={n_no} carry={n_carry}")

    # ----- Projection (candidates only) -----
    carry_aptamers = {c["aptamer"]: c for c in candidates if c["carry_to_projection"]}
    projection_rows = []

    for outcome_label, meta in OUTCOMES.items():
        if outcome_label == "AF":
            continue
        if not meta["path"].exists():
            continue
        outcome_df = pl.read_parquet(meta["path"])
        _log(f"{outcome_label} (build {meta['build']}, {outcome_df.height:,} rows)")
        instr_for = instr_b37 if meta["build"] == "GRCh37" else matched
        sub = instr_for.filter(pl.col("aptamer").is_in(list(carry_aptamers.keys())))

        for apt, cand in carry_aptamers.items():
            apt_sub = sub.filter(pl.col("aptamer") == apt)
            if apt_sub.height == 0:
                continue
            kept, _ = harmonize(apt_sub, outcome_df, build=meta["build"])
            m = run_mr(kept)
            d_T = cand["d_T"]
            theta = (d_T * m["beta"]) if (d_T != 0 and math.isfinite(m["beta"])) else math.nan
            OR_int = math.exp(theta) if math.isfinite(theta) else math.nan
            labels = []
            if d_T != 0 and math.isfinite(theta) and m["p"] is not None and math.isfinite(m["p"]) and m["p"] < 0.05:
                if theta < 0:
                    labels.append("possible_shared_benefit_direction")
                else:
                    if outcome_label.startswith("HF"):
                        labels.append("possible_HF_harm_direction")
                    elif "stroke" in outcome_label:
                        labels.append("possible_stroke_harm_direction")
            projection_rows.append({
                "aptamer": apt, "gene": cand["gene"],
                "outcome": outcome_label, "outcome_role": meta["role"],
                "outcome_build": meta["build"],
                "candidate_label": cand["candidate_label"],
                "d_T": d_T,
                "beta_outcome": m["beta"], "se_outcome": m["se"], "p_outcome": m["p"],
                "neg_log10p_outcome": m["neg_log10p"],
                "n_instruments": m["n_instruments"], "method": m["method"],
                "theta": theta, "OR_intervention": OR_int,
                "weighted_median_beta_outcome": m.get("weighted_median_beta"),
                "Q_heterogeneity_outcome": m.get("Q"),
                "provisional_labels": "|".join(labels),
                "scope": "CVD150_screen",
            })
    pl.DataFrame(projection_rows).write_csv(OUT_PROJ, separator="\t")
    pl.DataFrame(projection_rows).write_csv(OUT_DIRECTION, separator="\t")
    _log(f"wrote {OUT_PROJ}")

    # Summary md
    n_with_label = sum(1 for r in projection_rows if r["provisional_labels"])
    md = [
        "# CVD150 staged discovery — candidate summary (D037, NOT proteome-wide)",
        "",
        f"**Scope:** CVD-prioritized deCODE CVD150 target screen.",
        f"**FDR:** Benjamini-Hochberg WITHIN the CVD150 aptamer universe (n={len(mr_rows)}).",
        f"  → labelled `CVD150_screen_FDR_NOT_proteome_wide`.",
        "",
        "## AF anchor MR (Phase 3B)",
        "",
        f"- CVD150_AF_FDR_candidate: **{n_fdr}**",
        f"- CVD150_AF_nominal_candidate: **{n_nom}**",
        f"- CVD150_AF_no_signal: **{n_no}**",
        f"- carry_to_projection (FDR + nominal + sanity): **{n_carry}**",
        "",
        "## Projection (Phase 3C; candidates only)",
        "",
        f"- aptamer × outcome rows: {len(projection_rows)}",
        f"- rows with provisional_label fired: {n_with_label}",
        "",
        "## Coloc + final classification",
        "",
        "**STILL LOCKED (D037)** — protocol explicit approval required for Phase 3D coloc.",
    ]
    OUT_SUMMARY.write_text("\n".join(md))
    _log(f"wrote {OUT_SUMMARY}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
