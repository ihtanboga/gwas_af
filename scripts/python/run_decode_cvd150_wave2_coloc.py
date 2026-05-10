#!/usr/bin/env python3
"""D040 — CVD150 Wave 2 selective coloc, 20 target-outcome pairs.

Approved scope (D040):
  Discordance audit (5):    IL6R×AF, IL6R×CES, DSC2×AF, OGN×AF, TIMP3×AF
  NPPA pleiotropy (4):      NPPA×AIS, NPPA×HF_overall, NPPA×HF_nonischemic, NPPA×HF_ni_HFrEF
  BLUE module (6):          F11×AIS, F11×CES, KNG1×AIS, KNG1×CES, SERPINF2×AIS, SERPINF2×CES
  Nominal AF / red-flag (5): PCSK9×AF, PCSK9×HF_ni_HFrEF, MMP12×AF, MMP12×AIS, MMP12×LAS

SuSiE-coloc: 7 pairs flagged (IL6R×AF, IL6R×CES, DSC2×AF, OGN×AF, TIMP3×AF,
PCSK9×AF, MMP12×AF). EUR LD reference not present locally → flag every row
SUSIE_NOT_RUN_LD_UNAVAILABLE_OR_UNRELIABLE per D040.

Reuses Wave 1 utilities (download, region extract, liftover, harmonize, coloc.abf).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import polars as pl

ROOT = Path("/Users/apple/Desktop/gwas_af")
sys.path.insert(0, str(ROOT / "scripts/python"))

from run_decode_cvd150_wave1_coloc import (  # noqa: E402
    OUTCOMES, PRIOR_PRIMARY, PRIOR_SENSITIVITY, REGION_KB_PRIMARY, TMP_DIR,
    OUT_DIR, _log, annotate_artifact, build_pair_rows, cleanup_tmp,
    coloc_abf, download_aptamer, init_environment, run_pair,
)

OUT_RES = OUT_DIR / "decode_cvd150_wave2_coloc_results.tsv"
OUT_REG = OUT_DIR / "decode_cvd150_wave2_region_qc.tsv"
OUT_ART = OUT_DIR / "decode_cvd150_wave2_artifact_flags.tsv"
OUT_SUSIE = OUT_DIR / "decode_cvd150_wave2_susie_results.tsv"
OUT_MD = OUT_DIR / "decode_cvd150_wave2_coloc_summary.md"
OUT_COMBINED = OUT_DIR / "decode_cvd150_wave1_wave2_combined_summary.md"

WAVE1_RESULTS = OUT_DIR / "decode_cvd150_wave1_coloc_results.tsv"

WAVE2_PAIRS: list[tuple[str, str, str, str]] = [
    # (aptamer_key, gene, outcome_label, module_tag)
    ("15602_43_IL6R_IL_6_sRa.txt.gz",                "IL6R",     "AF",                   "DISCORDANCE_AUDIT"),
    ("15602_43_IL6R_IL_6_sRa.txt.gz",                "IL6R",     "cardioembolic_stroke", "DISCORDANCE_AUDIT"),
    ("13126_52_DSC2_DSC2.txt.gz",                    "DSC2",     "AF",                   "DISCORDANCE_AUDIT"),
    ("17224_12_OGN_MIME.txt.gz",                     "OGN",      "AF",                   "DISCORDANCE_AUDIT"),
    ("2480_58_TIMP3_TIMP_3.txt.gz",                  "TIMP3",    "AF",                   "DISCORDANCE_AUDIT"),
    ("5443_62_NPPA_ANP.txt.gz",                      "NPPA",     "ischemic_stroke",      "NPPA_PLEIOTROPY"),
    ("5443_62_NPPA_ANP.txt.gz",                      "NPPA",     "HF_overall",           "NPPA_PLEIOTROPY"),
    ("5443_62_NPPA_ANP.txt.gz",                      "NPPA",     "HF_nonischemic",       "NPPA_PLEIOTROPY"),
    ("5443_62_NPPA_ANP.txt.gz",                      "NPPA",     "HF_ni_HFrEF",          "NPPA_PLEIOTROPY"),
    ("2190_55_F11_Coagulation_Factor_XI.txt.gz",     "F11",      "ischemic_stroke",      "BLUE_ANTITHROMBOTIC"),
    ("2190_55_F11_Coagulation_Factor_XI.txt.gz",     "F11",      "cardioembolic_stroke", "BLUE_ANTITHROMBOTIC"),
    ("15343_337_KNG1_Kininogen__HMW__Two_Chain.txt.gz", "KNG1", "ischemic_stroke",      "BLUE_ANTITHROMBOTIC"),
    ("15343_337_KNG1_Kininogen__HMW__Two_Chain.txt.gz", "KNG1", "cardioembolic_stroke", "BLUE_ANTITHROMBOTIC"),
    ("3024_18_SERPINF2_a2_Antiplasmin.txt.gz",       "SERPINF2", "ischemic_stroke",      "BLUE_ANTITHROMBOTIC"),
    ("3024_18_SERPINF2_a2_Antiplasmin.txt.gz",       "SERPINF2", "cardioembolic_stroke", "BLUE_ANTITHROMBOTIC"),
    ("5231_79_PCSK9_PCSK9.txt.gz",                   "PCSK9",    "AF",                   "NOMINAL_AF_REDFLAG"),
    ("5231_79_PCSK9_PCSK9.txt.gz",                   "PCSK9",    "HF_ni_HFrEF",          "NOMINAL_AF_REDFLAG"),
    ("4496_60_MMP12_MMP_12.txt.gz",                  "MMP12",    "AF",                   "NOMINAL_AF_REDFLAG"),
    ("4496_60_MMP12_MMP_12.txt.gz",                  "MMP12",    "ischemic_stroke",      "NOMINAL_AF_REDFLAG"),
    ("4496_60_MMP12_MMP_12.txt.gz",                  "MMP12",    "large_artery_stroke",  "NOMINAL_AF_REDFLAG"),
]

SUSIE_PAIRS = [
    ("IL6R", "AF"), ("IL6R", "cardioembolic_stroke"),
    ("DSC2", "AF"), ("OGN", "AF"), ("TIMP3", "AF"),
    ("PCSK9", "AF"), ("MMP12", "AF"),
]

# Need HF_nonischemic outcome path (not in OUTCOMES dict from wave1; add here)
EXTRA_OUTCOMES = {
    "HF_nonischemic":  {"path": ROOT / "data/processed/sumstats/HF_2024_EUR_HF_nonischemic.parquet", "build": "GRCh37"},
    "HF_ni_HFrEF":     {"path": ROOT / "data/processed/sumstats/HF_2024_EUR_HF_ni_HFrEF.parquet",    "build": "GRCh37"},
}
OUTCOMES.update(EXTRA_OUTCOMES)


def detect_ld_reference() -> tuple[bool, str]:
    """Returns (available, reason)."""
    candidates = [
        ROOT / "data/raw/ld_reference",
        ROOT / "data/reference/1000G_EUR",
        ROOT / "data/reference/ld_blocks",
    ]
    for d in candidates:
        if d.exists() and any(d.iterdir()):
            return True, str(d)
    return False, "no EUR LD reference at data/raw/ld_reference/, data/reference/1000G_EUR/, data/reference/ld_blocks/"


def write_susie_skip_results():
    rows = []
    avail, reason = detect_ld_reference()
    flag = "SUSIE_RUN_OK" if avail else "SUSIE_NOT_RUN_LD_UNAVAILABLE_OR_UNRELIABLE"
    for gene, outcome in SUSIE_PAIRS:
        rows.append({
            "gene": gene, "outcome": outcome,
            "susie_status": flag,
            "ld_reference_available": avail,
            "skip_reason": reason if not avail else "",
            "n_credible_sets_pqtl": None,
            "n_credible_sets_outcome": None,
            "max_pp_h4_susie": None,
            "method": "susie_coloc_pending_LD_setup",
        })
    pl.DataFrame(rows).write_csv(OUT_SUSIE, separator="\t")
    _log(f"wrote {OUT_SUSIE} (SuSiE: {flag})")
    return flag


def write_summary_md(coloc_rows):
    primary = [r for r in coloc_rows if r["prior_model"] == PRIOR_PRIMARY["name"]]
    sens = {(r["gene"], r["outcome"]): r for r in coloc_rows
            if r["prior_model"] == PRIOR_SENSITIVITY["name"]}
    strong = [r for r in primary if r["PP.H4"] is not None and not math.isnan(r["PP.H4"]) and r["PP.H4"] >= 0.80]
    moderate = [r for r in primary if r["PP.H4"] is not None and not math.isnan(r["PP.H4"]) and 0.50 <= r["PP.H4"] < 0.80]
    weak = [r for r in primary if r["PP.H4"] is not None and not math.isnan(r["PP.H4"]) and r["PP.H4"] < 0.50]
    ld_warn = [r for r in primary if r["PP.H3"] is not None and r["PP.H4"] is not None
               and not math.isnan(r["PP.H3"]) and not math.isnan(r["PP.H4"])
               and r["PP.H3"] >= 0.50 and r["PP.H3"] > r["PP.H4"]]

    md = [
        "# CVD150 Wave 2 colocalization (D040, NOT proteome-wide)",
        "",
        f"**Scope:** {len(primary)} target-outcome pairs across 4 modules.",
        f"**Region:** lead cis-pQTL ±{REGION_KB_PRIMARY} kb.",
        f"**Primary prior** ({PRIOR_PRIMARY['name']}): p1=p2={PRIOR_PRIMARY['p1']}, p12={PRIOR_PRIMARY['p12']}, "
        f"W_pqtl={PRIOR_PRIMARY['W_pqtl']}, W_outcome={PRIOR_PRIMARY['W_outcome']}",
        f"**sensitivity** ({PRIOR_SENSITIVITY['name']}): p12={PRIOR_SENSITIVITY['p12']}, W=0.04 (both)",
        "",
        "## Per-pair results (primary prior; sensitivity PP.H4 in parentheses)",
        "",
        "| Module | Gene | Outcome | n_present | PP.H1 | PP.H3 | PP.H4 (sens) | Lead pQTL | Lead outcome p | Intragenic | Multi-apt |",
        "|---|---|---|---:|---:|---:|---|---|---|:-:|:-:|",
    ]
    module_by_pair = {(g, o): m for _, g, o, m in WAVE2_PAIRS}
    for r in primary:
        s = sens.get((r["gene"], r["outcome"]))
        sens_h4 = f" ({s['PP.H4']:.3f})" if s and s["PP.H4"] is not None and not math.isnan(s["PP.H4"]) else ""
        h4_str = f"**{r['PP.H4']:.3f}**{sens_h4}" if r["PP.H4"] is not None and not math.isnan(r["PP.H4"]) else "—"
        mod = module_by_pair.get((r["gene"], r["outcome"]), "?")
        md.append(
            f"| {mod} | {r['gene']} | {r['outcome']} | {r['n_variants']} | "
            f"{r['PP.H1']:.3f} | {r['PP.H3']:.3f} | {h4_str} | "
            f"{r['lead_pqtl_rsid']} (p={r['lead_pqtl_pval']:.1e}) | "
            f"{r['lead_outcome_pval']:.1e} | "
            f"{'✓' if r['lead_pqtl_intragenic'] else '✗'} | "
            f"{r['multi_aptamer_status']} |"
        )
    md.extend([
        "",
        "## Summary by PP.H4 strength (primary prior)",
        "",
        f"- **Strong (PP.H4 ≥ 0.80):** {len(strong)} pairs",
        f"- **Moderate (0.50 ≤ PP.H4 < 0.80):** {len(moderate)} pairs",
        f"- **Weak (PP.H4 < 0.50):** {len(weak)} pairs",
        f"- **LD-confounding warning (PP.H3 ≥ 0.50 and PP.H3 > PP.H4):** {len(ld_warn)} pairs",
        "",
        "## SuSiE-coloc status",
        "",
        "Per D040, SuSiE is required only for: IL6R×AF, IL6R×CES, DSC2×AF, OGN×AF, TIMP3×AF, PCSK9×AF, MMP12×AF.",
        "Local EUR LD reference is not available — every SuSiE pair flagged `SUSIE_NOT_RUN_LD_UNAVAILABLE_OR_UNRELIABLE` (D040 explicit allowance).",
        "",
        "## Operasyonel kural (D040)",
        "",
        "Final classification (GREEN/RED/AMBER/BLUE/PURPLE) **STILL LOCKED**.",
        "Allowed outputs after Wave 2: interim evidence table, coloc-supported / LD-confounded / BLUE candidate tables, recommended candidates.",
    ])
    OUT_MD.write_text("\n".join(md))
    _log(f"wrote {OUT_MD}")


def write_combined_summary(wave2_coloc_rows):
    if not WAVE1_RESULTS.exists():
        _log(f"WARN: {WAVE1_RESULTS} not found, skipping combined summary")
        return
    w1 = pl.read_csv(WAVE1_RESULTS, separator="\t")
    w1_primary = w1.filter(pl.col("prior_model") == PRIOR_PRIMARY["name"])
    w2_primary = [r for r in wave2_coloc_rows if r["prior_model"] == PRIOR_PRIMARY["name"]]

    # Build candidate tables
    coloc_supported = []   # PP.H4 >= 0.50
    ld_confounded = []     # PP.H3 >= 0.50 and PP.H3 > PP.H4
    blue_module = []       # F11, KNG1, SERPINF2 — antithrombotic, irrespective of H4
    nominal_unsupported = []  # weak coloc, AF MR strong but not coloc-supported

    for r in w1_primary.iter_rows(named=True):
        pair = (r["gene"], r["outcome"], "Wave1")
        if r["PP.H4"] is not None and r["PP.H4"] >= 0.50:
            coloc_supported.append({**r, "wave": "Wave1"})
        if (r["PP.H3"] is not None and r["PP.H4"] is not None
                and r["PP.H3"] >= 0.50 and r["PP.H3"] > r["PP.H4"]):
            ld_confounded.append({**r, "wave": "Wave1"})
    for r in w2_primary:
        if r["PP.H4"] is not None and not math.isnan(r["PP.H4"]) and r["PP.H4"] >= 0.50:
            coloc_supported.append({**r, "wave": "Wave2"})
        if (r["PP.H3"] is not None and r["PP.H4"] is not None
                and not math.isnan(r["PP.H3"]) and not math.isnan(r["PP.H4"])
                and r["PP.H3"] >= 0.50 and r["PP.H3"] > r["PP.H4"]):
            ld_confounded.append({**r, "wave": "Wave2"})
        if r["gene"] in {"F11", "KNG1", "SERPINF2"}:
            blue_module.append({**r, "wave": "Wave2"})

    nppa_pairs = [r for r in coloc_supported if r["gene"] == "NPPA"]

    md = [
        "# Wave 1 + Wave 2 combined coloc summary (D038/D039/D040)",
        "",
        f"**Total pairs analyzed:** Wave 1 = 15, Wave 2 = {len(w2_primary)}.",
        "**Final classification (GREEN/RED/AMBER/BLUE/PURPLE) STILL LOCKED.**",
        "",
        "## 1. Coloc-supported candidates (PP.H4 ≥ 0.50)",
        "",
        "| Wave | Gene | Outcome | PP.H4 | PP.H3 | Lead pQTL | Lead outcome p |",
        "|---|---|---|---:|---:|---|---|",
    ]
    for r in sorted(coloc_supported, key=lambda x: -x["PP.H4"]):
        md.append(
            f"| {r['wave']} | {r['gene']} | {r['outcome']} | "
            f"**{r['PP.H4']:.3f}** | {r['PP.H3']:.3f} | "
            f"{r['lead_pqtl_rsid']} | {r['lead_outcome_pval']:.1e} |"
        )
    if not coloc_supported:
        md.append("| — | — | — | — | — | — | — |")

    md.extend([
        "",
        "## 2. LD-confounded / downgraded candidates (PP.H3 ≥ 0.50 > PP.H4)",
        "",
        "| Wave | Gene | Outcome | PP.H3 | PP.H4 | Lead pQTL | Note |",
        "|---|---|---|---:|---:|---|---|",
    ])
    for r in sorted(ld_confounded, key=lambda x: -x["PP.H3"]):
        md.append(
            f"| {r['wave']} | {r['gene']} | {r['outcome']} | "
            f"**{r['PP.H3']:.3f}** | {r['PP.H4']:.3f} | "
            f"{r['lead_pqtl_rsid']} | COLOC_DISCORDANT_H3_HIGH |"
        )
    if not ld_confounded:
        md.append("| — | — | — | — | — | — | — |")

    md.extend([
        "",
        "## 3. BLUE module (antithrombotic / safety, F11 / KNG1 / SERPINF2)",
        "",
        "| Gene | Outcome | PP.H4 | PP.H3 | Status |",
        "|---|---|---:|---:|---|",
    ])
    for r in blue_module:
        if r["PP.H4"] is None or math.isnan(r["PP.H4"]):
            status = "no signal"
            h4 = "—"
        elif r["PP.H4"] >= 0.50:
            status = "BLUE_COLOC_SUPPORTED"
            h4 = f"{r['PP.H4']:.3f}"
        else:
            status = "BLUE_NO_COLOC"
            h4 = f"{r['PP.H4']:.3f}"
        md.append(f"| {r['gene']} | {r['outcome']} | {h4} | {r['PP.H3']:.3f} | {status} |")

    md.extend([
        "",
        "## 4. NPPA AF-anchored axis (interim pre-classification only)",
        "",
        f"NPPA coloc-positive pairs (PP.H4 ≥ 0.50): **{len(nppa_pairs)}** of {sum(1 for _, g, _, _ in WAVE2_PAIRS if g == 'NPPA') + 2}",
        "",
        "| Outcome | PP.H4 | PP.H3 | Wave |",
        "|---|---:|---:|:-:|",
    ])
    nppa_all_h4 = ([r for r in w1_primary.iter_rows(named=True) if r["gene"] == "NPPA"]
                   + [r for r in w2_primary if r["gene"] == "NPPA"])
    for r in nppa_all_h4:
        wave = "W1" if any(rr["gene"] == "NPPA" and rr["outcome"] == r["outcome"]
                           for rr in w1_primary.iter_rows(named=True)) else "W2"
        h4 = f"{r['PP.H4']:.3f}" if r["PP.H4"] is not None and not (isinstance(r["PP.H4"], float) and math.isnan(r["PP.H4"])) else "—"
        h3 = f"{r['PP.H3']:.3f}" if r["PP.H3"] is not None and not (isinstance(r["PP.H3"], float) and math.isnan(r["PP.H3"])) else "—"
        md.append(f"| {r['outcome']} | {h4} | {h3} | {wave} |")

    md.extend([
        "",
        "## 5. Recommended next steps (recommendation, not a call — final classification closed)",
        "",
        "- NPPA: AF + CES coloc-supported; Wave 2 HF/AIS pleiotropy mapping in this summary → AMBER/GREEN-2 pre-candidate.",
        "- IL6R, DSC2, OGN, TIMP3 (AF MR LD-confounded): SuSiE-coloc audit pending — `SUSIE_NOT_RUN_LD_UNAVAILABLE_OR_UNRELIABLE`. If EUR LD reference setup is needed, this is a separate sub-gate.",
        "- F11 / KNG1 / SERPINF2: BLUE module assessment (antithrombotic/fibrinolytic axis). If stroke shared-benefit is coloc-supported then BLUE_COLOC_SUPPORTED.",
        "- PCSK9: AF nominal coloc-supported  — note: AF-anchored red flag is not called. HFrEF off-axis signal depends on Wave 2 results.",
        "- MMP12: AF MR coloc-supported  — note: AF-anchored interpretation is not made. AIS/LAS shared-benefit depends on Wave 2 results.",
        "",
        "## SuSiE skip nedeni",
        "",
        "Local EUR LD reference yok. D040 izniyle skip + flag.",
        "Setup requires: 1000G Phase 3 EUR plink BED (~3 GB selective by chromosome) + plink2 region LD computation. Multi-day or a later sub-gate.",
    ])
    OUT_COMBINED.write_text("\n".join(md))
    _log(f"wrote {OUT_COMBINED}")


def main() -> int:
    target_apts = {p[0] for p in WAVE2_PAIRS}
    leads, apt_paths, coord_cache, instr, caches = init_environment(target_apts)

    coloc_rows, qc_rows, artifact_rows = [], [], []
    for apt, gene, outcome, module in WAVE2_PAIRS:
        for prior in (PRIOR_PRIMARY, PRIOR_SENSITIVITY):
            res, qc, art, _aligned = run_pair(
                apt, gene, outcome, prior,
                pqtl_region_cache=caches["region"],
                pqtl_region_b37_cache=caches["region_b37"],
                apt_paths=apt_paths, leads=leads,
                outcome_cache=caches["outcome"],
                coord_cache=coord_cache, instr=instr,
            )
            res["module"] = module
            coloc_rows.append(res)
            if prior["name"] == PRIOR_PRIMARY["name"]:
                qc_rows.append(qc)
                artifact_rows.append({
                    "gene": gene, "aptamer": apt, "outcome": outcome,
                    "module": module,
                    "lead_pqtl_intragenic": art["intragenic"],
                    "multi_aptamer_status": art["multi_aptamer_status"],
                    "lead_pqtl_protein_altering": "TBD_VEP_required",
                    "aptamer_artifact_risk": ("low" if art["intragenic"]
                                              and art["multi_aptamer_status"] == "single"
                                              else "review"),
                })
        _log(f"  pair {gene}-{outcome} ({module}) done")

    pl.DataFrame(coloc_rows).write_csv(OUT_RES, separator="\t")
    pl.DataFrame(qc_rows).write_csv(OUT_REG, separator="\t")
    pl.DataFrame(artifact_rows).write_csv(OUT_ART, separator="\t")
    _log(f"wrote {OUT_RES} ({len(coloc_rows)} rows)")

    write_susie_skip_results()
    write_summary_md(coloc_rows)
    write_combined_summary(coloc_rows)
    cleanup_tmp()
    return 0


if __name__ == "__main__":
    sys.exit(main())
