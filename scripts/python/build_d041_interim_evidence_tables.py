#!/usr/bin/env python3
"""D041 — Build interim evidence classes (Wave 1 + Wave 2 coloc + AF MR).

Outputs (per D041):
  results/interim/decode_cvd150_interim_evidence_classes.tsv
  results/interim/decode_cvd150_interim_candidate_summary.md
  results/interim/decode_cvd150_downgraded_candidates.tsv
  results/interim/decode_cvd150_blue_module_candidates.tsv
  results/interim/decode_cvd150_nppa_af_ces_axis_summary.md

NOTE: final_class_allowed = false on every row (D041 explicit).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import polars as pl

ROOT = Path("/Users/apple/Desktop/gwas_af")
INTERIM_DIR = ROOT / "results/interim"
COLOC_DIR = ROOT / "results/coloc"
MR_DIR = ROOT / "results/mr"

W1 = COLOC_DIR / "decode_cvd150_wave1_coloc_results.tsv"
W2 = COLOC_DIR / "decode_cvd150_wave2_coloc_results.tsv"
MR_SHORTLIST = MR_DIR / "decode_cvd150_AF_candidate_shortlist.tsv"
MR_FULL = MR_DIR / "decode_cvd150_pqtl_to_AF.tsv"

# D041 specs (per gene; primary outcome chosen by protocol interpretation):
#   gene → (aptamer, primary_outcome, module, interim_label, AF_anchored, blue_module, downgrade_reason, replication_required, notes)
GENE_SPECS = {
    "NPPA":     ("5443_62_NPPA_ANP.txt.gz",
                 "AF",                          # primary anchor outcome
                 "AF_CES_CARDIOEMBOLIC_AXIS",
                 "INTERIM_AMBER_AF_CES_SHARED_BENEFIT_CANDIDATE",
                 True, False, "",
                 "independent pQTL replication (UKB-PPP/SCALLOP) + replicated CES support",
                 "Strong AF coloc (PP.H4=0.911); moderate CES coloc (0.565); HF/AIS coloc-not-supported. AF→CES specific axis."),
    "F11":      ("2190_55_F11_Coagulation_Factor_XI.txt.gz",
                 "ischemic_stroke",
                 "ANTITHROMBOTIC_STROKE_AXIS",
                 "INTERIM_BLUE_COLOC_SUPPORTED_ANTITHROMBOTIC_STROKE_AXIS",
                 False, True, "",
                 "stroke outcome replication (FinnGen, MEGASTROKE)",
                 "F11×AIS PP.H4=0.993; F11×CES PP.H4=0.989; AF MR no_signal. Antithrombotic anchor; outside AF-anchored atlas."),
    "KNG1":     ("15343_337_KNG1_Kininogen__HMW__Two_Chain.txt.gz",
                 "ischemic_stroke",
                 "KALLIKREIN_KININ_STROKE_AXIS",
                 "INTERIM_BLUE_COLOC_SUPPORTED_KALLIKREIN_KININ_STROKE_AXIS",
                 False, True, "",
                 "stroke outcome replication",
                 "KNG1×AIS PP.H4=0.991; KNG1×CES PP.H4=0.933; AF MR no_signal. Kallikrein-kinin axis."),
    "MMP12":    ("4496_60_MMP12_MMP_12.txt.gz",
                 "large_artery_stroke",
                 "ATHEROTHROMBOTIC_STROKE_AXIS",
                 "INTERIM_BLUE_COLOC_SUPPORTED_ATHEROTHROMBOTIC_STROKE_AXIS",
                 False, True, "",
                 "stroke outcome replication (MEGASTROKE LAS)",
                 "MMP12×LAS PP.H4=0.916; MMP12×AIS PP.H4=0.838; AF coloc not supported. Atherosclerosis/elastin axis."),
    "IL6R":     ("15602_43_IL6R_IL_6_sRa.txt.gz",
                 "AF",
                 "DOWNGRADED",
                 "DOWNGRADED_COLOC_DISCORDANT_H3_HIGH",
                 False, False, "AF coloc PP.H3=1.0, PP.H4=0; cis-pQTL not shared with AF causal variant. SuSiE/conditional rescue not run.",
                 "EUR LD reference + SuSiE-coloc audit (D042 optional)",
                 "AF MR strong (q=4.1e-13) but coloc rejects shared variant. LD-confounded under coloc.abf."),
    "DSC2":     ("13126_52_DSC2_DSC2.txt.gz",
                 "AF",
                 "DOWNGRADED",
                 "DOWNGRADED_COLOC_DISCORDANT_H3_HIGH",
                 False, False, "AF coloc PP.H3=1.0, PP.H4=0.",
                 "EUR LD reference + SuSiE-coloc audit (D042 optional)",
                 "AF MR strong (q=9e-6) but coloc rejects shared variant. LD-confounded under coloc.abf."),
    "OGN":      ("17224_12_OGN_MIME.txt.gz",
                 "AF",
                 "DOWNGRADED",
                 "DOWNGRADED_AF_MR_NOT_COLOC_SUPPORTED",
                 False, False, "AF coloc PP.H4=0.009, PP.H1=0.568, PP.H3=0.424; borderline distinct/non-coloc.",
                 "SuSiE-coloc audit (D042 optional)",
                 "AF MR nominal but coloc not supported."),
    "TIMP3":    ("2480_58_TIMP3_TIMP_3.txt.gz",
                 "AF",
                 "DOWNGRADED",
                 "DOWNGRADED_AF_MR_NOT_COLOC_SUPPORTED",
                 False, False, "AF coloc PP.H4=0.003, PP.H1=0.932; outcome region without independent signal.",
                 "SuSiE-coloc audit (D042 optional)",
                 "AF MR nominal but coloc not supported (PP.H1 dominant)."),
    "SERPINF2": ("3024_18_SERPINF2_a2_Antiplasmin.txt.gz",
                 "ischemic_stroke",
                 "DOWNGRADED",
                 "DOWNGRADED_BLUE_NO_COLOC",
                 False, False, "AIS PP.H3=0.703 (LD-confounded); CES PP.H4=0.096 weak.",
                 "stroke outcome replication + LD audit",
                 "Excluded from BLUE supported module."),
    "PCSK9":    ("5231_79_PCSK9_PCSK9.txt.gz",
                 "AF",
                 "DOWNGRADED",
                 "DOWNGRADED_AF_MR_NOT_COLOC_SUPPORTED",
                 False, False, "AF PP.H4=0.094 (H1=0.852); HFrEF PP.H4=0.201 (H1=0.611). Neither AF nor HFrEF coloc-supported.",
                 "SuSiE-coloc audit (D042 optional)",
                 "Smoke-level nominal AF/HFrEF signals not rescued by region-level coloc. No AF-anchored red flag."),
}

PRIMARY_OUTCOME_LABEL_MAP = {
    "AF": "AF",
    "ischemic_stroke": "AIS",
    "cardioembolic_stroke": "CES",
    "large_artery_stroke": "LAS",
    "small_vessel_stroke": "SVS",
    "HF_overall": "HF_overall",
    "HF_ni_HFpEF": "HF_ni_HFpEF",
    "HF_ni_HFrEF": "HF_ni_HFrEF",
    "HF_nonischemic": "HF_nonischemic",
}


def _load_coloc():
    rows = []
    for path in (W1, W2):
        if path.exists():
            df = pl.read_csv(path, separator="\t")
            df = df.filter(pl.col("prior_model") == "default_quant_cc")
            for r in df.iter_rows(named=True):
                rows.append({"gene": r["gene"], "outcome": r["outcome"],
                             "PP.H4": r["PP.H4"], "PP.H3": r["PP.H3"],
                             "wave": "Wave1" if path == W1 else "Wave2"})
    sens_rows = {}
    for path in (W1, W2):
        if path.exists():
            df = pl.read_csv(path, separator="\t")
            df = df.filter(pl.col("prior_model") == "flat_W04")
            for r in df.iter_rows(named=True):
                sens_rows[(r["gene"], r["outcome"])] = r["PP.H4"]
    by_gene_outcome = {}
    for r in rows:
        by_gene_outcome[(r["gene"], r["outcome"])] = r
    return by_gene_outcome, sens_rows


def _load_mr():
    """Per (gene, aptamer): raw p + FDR q from CVD150 AF MR full table."""
    af_mr = {}
    if MR_FULL.exists():
        df = pl.read_csv(MR_FULL, separator="\t", infer_schema_length=0)
        for r in df.iter_rows(named=True):
            try:
                p = float(r["p"]) if r["p"] not in (None, "") else None
            except (ValueError, TypeError):
                p = None
            try:
                q = float(r["fdr_bh_within_cvd150"]) if r["fdr_bh_within_cvd150"] not in (None, "") else None
            except (ValueError, TypeError):
                q = None
            af_mr[(r["gene"], r["aptamer"])] = {"p": p, "q": q}
    return af_mr


def evidence_strength(coloc_h4, h3) -> str:
    if coloc_h4 is None or (isinstance(coloc_h4, float) and math.isnan(coloc_h4)):
        return "no_signal"
    if coloc_h4 >= 0.80:
        return "strong"
    if coloc_h4 >= 0.50:
        return "moderate"
    if h3 is not None and h3 >= 0.50 and h3 > coloc_h4:
        return "ld_confounded"
    return "weak"


def main() -> int:
    INTERIM_DIR.mkdir(parents=True, exist_ok=True)
    by_pair, sens = _load_coloc()
    af_mr = _load_mr()

    rows = []
    for gene, spec in GENE_SPECS.items():
        (apt, primary_outcome, module, label,
         af_anchored, blue, dn_reason, repl, notes) = spec
        af_coloc = by_pair.get((gene, "AF"), {})
        primary_coloc = by_pair.get((gene, primary_outcome), {})
        af_h4 = af_coloc.get("PP.H4")
        primary_h4 = primary_coloc.get("PP.H4")
        primary_h3 = primary_coloc.get("PP.H3")
        sens_h4 = sens.get((gene, primary_outcome))
        mr = af_mr.get((gene, apt), {})
        rows.append({
            "gene": gene,
            "aptamer_id": apt,
            "module": module,
            "interim_label": label,
            "AF_MR_p": mr.get("p"),
            "AF_MR_q_CVD150": mr.get("q"),
            "AF_coloc_PPH4": af_h4,
            "primary_outcome": PRIMARY_OUTCOME_LABEL_MAP.get(primary_outcome, primary_outcome),
            "primary_outcome_coloc_PPH4": primary_h4,
            "sensitivity_PPH4": sens_h4,
            "evidence_strength": evidence_strength(primary_h4, primary_h3),
            "AF_anchored": af_anchored,
            "blue_module": blue,
            "downgrade_reason": dn_reason,
            "replication_required": repl,
            "final_class_allowed": False,
            "notes": notes,
        })

    df = pl.DataFrame(rows)
    out_main = INTERIM_DIR / "decode_cvd150_interim_evidence_classes.tsv"
    df.write_csv(out_main, separator="\t")
    print(f"wrote {out_main} ({df.height} rows)")

    # Sub-tables
    blue = df.filter(pl.col("blue_module"))
    blue_path = INTERIM_DIR / "decode_cvd150_blue_module_candidates.tsv"
    blue.write_csv(blue_path, separator="\t")
    print(f"wrote {blue_path} ({blue.height} rows)")

    downgraded = df.filter(pl.col("interim_label").str.starts_with("DOWNGRADED"))
    downgraded_path = INTERIM_DIR / "decode_cvd150_downgraded_candidates.tsv"
    downgraded.write_csv(downgraded_path, separator="\t")
    print(f"wrote {downgraded_path} ({downgraded.height} rows)")

    # Markdown summary
    md = [
        "# CVD150 interim evidence classes (D041, AF-anchored therapeutic atlas)",
        "",
        "**Final GREEN/RED/AMBER classification STILL LOCKED.** D041 opened only interim classes.",
        "",
        f"**Total candidates:** {df.height} ({blue.height} BLUE, {downgraded.height} downgraded, {df.height - blue.height - downgraded.height} AF-anchored interim).",
        "",
        "## 1. AF-anchored interim candidates (NPPA only)",
        "",
        "| Gene | Module | Interim label | AF coloc PP.H4 | Primary outcome | Primary PP.H4 | Sens PP.H4 | Replication required |",
        "|---|---|---|---:|---|---:|---:|---|",
    ]
    for r in df.filter(pl.col("AF_anchored")).iter_rows(named=True):
        af_h4_str = f"{r['AF_coloc_PPH4']:.3f}" if r["AF_coloc_PPH4"] is not None else "—"
        p_h4 = f"{r['primary_outcome_coloc_PPH4']:.3f}" if r["primary_outcome_coloc_PPH4"] is not None else "—"
        s_h4 = f"{r['sensitivity_PPH4']:.3f}" if r["sensitivity_PPH4"] is not None else "—"
        md.append(f"| {r['gene']} | {r['module']} | {r['interim_label']} | "
                  f"**{af_h4_str}** | {r['primary_outcome']} | **{p_h4}** | {s_h4} | "
                  f"{r['replication_required']} |")

    md.extend([
        "",
        "## 2. BLUE module (non-AF-anchored, antithrombotic/atherothrombotic)",
        "",
        "| Gene | Module | Interim label | Primary outcome | Primary PP.H4 | Sens PP.H4 |",
        "|---|---|---|---|---:|---:|",
    ])
    for r in blue.iter_rows(named=True):
        p_h4 = f"{r['primary_outcome_coloc_PPH4']:.3f}" if r["primary_outcome_coloc_PPH4"] is not None else "—"
        s_h4 = f"{r['sensitivity_PPH4']:.3f}" if r["sensitivity_PPH4"] is not None else "—"
        md.append(f"| {r['gene']} | {r['module']} | {r['interim_label']} | "
                  f"{r['primary_outcome']} | **{p_h4}** | {s_h4} |")

    md.extend([
        "",
        "## 3. Downgraded / discordant candidates",
        "",
        "| Gene | Interim label | Downgrade reason | AF MR q | AF coloc PP.H4 |",
        "|---|---|---|---:|---:|",
    ])
    for r in downgraded.iter_rows(named=True):
        af_q = f"{r['AF_MR_q_CVD150']:.2e}" if r['AF_MR_q_CVD150'] is not None else "—"
        af_h4 = f"{r['AF_coloc_PPH4']:.3f}" if r['AF_coloc_PPH4'] is not None else "—"
        md.append(f"| {r['gene']} | {r['interim_label']} | {r['downgrade_reason'][:80]}... | "
                  f"{af_q} | {af_h4} |")

    md.extend([
        "",
        "## 4. Three-axis story (atlas headline)",
        "",
        "```",
        "NPPA       = AF-CES cardioembolic axis           (AF-anchored, INTERIM AMBER)",
        "F11/KNG1   = antithrombotic AIS/CES axis         (non-AF, BLUE)",
        "MMP12      = atherothrombotic AIS/LAS axis       (non-AF, BLUE)",
        "```",
        "",
        "## 5. final_class_allowed = false (D041 explicit)",
        "",
        "These interim classes are NOT final GREEN-2/RED/AMBER calls. Final classification requires:",
        "- Independent pQTL replication (UKB-PPP DAR, SCALLOP)",
        "- Outcome replication (FinnGen, MEGASTROKE)",
        "- LD-aware fine-mapping (D042 optional SuSiE) — IL6R/DSC2 rescue for",
        "",
        "## 6. Bilimsel yorumlama",
        "",
        "1. CVD150 targeted pQTL screen produced apparent AF MR candidates.",
        "2. Colocalization applied a sharp filter.",
        "3. **Only NPPA AF-anchored coloc-supported axis olarak survived.**",
        "4. NPPA is specific to AF + cardioembolic stroke; no pleiotropy on HF or general stroke.",
        "5. Outside the AF anchor **two strong stroke modules** were identified:",
        "   - F11/KNG1 antithrombotic AIS/CES axis",
        "   - MMP12 atherothrombotic AIS/LAS axis",
        "6. Attractive MR hits were dropped at colocalization: IL6R, DSC2, OGN, TIMP3, PCSK9, SERPINF2.",
        "",
        "**Atlas message:** The project's message is not 'MR found everything'; MR + colocalization eliminated false targets and separated genuine mechanism axes.",
    ])
    summary_path = INTERIM_DIR / "decode_cvd150_interim_candidate_summary.md"
    summary_path.write_text("\n".join(md))
    print(f"wrote {summary_path}")

    # NPPA dedicated summary
    nppa_md = [
        "# NPPA AF-CES cardioembolic axis — interim summary (D041)",
        "",
        f"**Aptamer:** 5443_62_NPPA_ANP",
        f"**Lead cis-pQTL:** rs145488887 (chr1:11767739, p=4.2e-12)",
        f"**Interim label:** INTERIM_AMBER_AF_CES_SHARED_BENEFIT_CANDIDATE",
        f"**Final classification:** LOCKED (D041)",
        "",
        "## Pleiotropy mapping (Wave 1 + Wave 2)",
        "",
        "| Outcome | PP.H4 | Sens PP.H4 | PP.H3 | Status |",
        "|---|---:|---:|---:|---|",
    ]
    nppa_outcomes = [
        ("AF", "Wave1"), ("cardioembolic_stroke", "Wave1"),
        ("ischemic_stroke", "Wave2"), ("HF_overall", "Wave2"),
        ("HF_nonischemic", "Wave2"), ("HF_ni_HFrEF", "Wave2"),
    ]
    for outcome, _wave in nppa_outcomes:
        c = by_pair.get(("NPPA", outcome), {})
        s = sens.get(("NPPA", outcome))
        h4 = c.get("PP.H4")
        h3 = c.get("PP.H3")
        if h4 is None:
            status = "missing"
            h4_s = "—"
        elif h4 >= 0.80:
            status = "STRONG_COLOC"
            h4_s = f"**{h4:.3f}**"
        elif h4 >= 0.50:
            status = "MODERATE_COLOC"
            h4_s = f"**{h4:.3f}**"
        else:
            status = "no_coloc_support"
            h4_s = f"{h4:.3f}"
        s_s = f"{s:.3f}" if s is not None else "—"
        h3_s = f"{h3:.3f}" if h3 is not None else "—"
        nppa_md.append(f"| {outcome} | {h4_s} | {s_s} | {h3_s} | {status} |")

    nppa_md.extend([
        "",
        "## Interpretation",
        "",
        "- **AF coloc strong** (PP.H4=0.911, sens 0.829) — primary AF-anchored signal confirmed.",
        "- **CES coloc moderate** (PP.H4=0.565, sens 0.395) — sensitivity prior conservative; AF→CES axis supported.",
        "- **AIS, HF_overall, HF_nonischemic, HF_ni_HFrEF: no coloc support** — NPPA AF-anchored signal does NOT generalize to other HF/AIS subtypes.",
        "",
        "**Conclusion:** NPPA is interim AF-CES cardioembolic axis candidate. Final GREEN-2 requires independent pQTL replication (UKB-PPP/SCALLOP) + replicated CES support.",
    ])
    nppa_path = INTERIM_DIR / "decode_cvd150_nppa_af_ces_axis_summary.md"
    nppa_path.write_text("\n".join(nppa_md))
    print(f"wrote {nppa_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
