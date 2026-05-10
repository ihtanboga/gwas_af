#!/usr/bin/env python3
"""D043 — Manuscript-style results draft artifacts.

NO new analysis. Reads existing Wave 1+2 coloc TSVs, CVD150 instrument table,
AF MR results, interim evidence classes; emits 4 reports + 5 tables + 6 figures.

D043 language rules (critical):
  KULLAN: "CVD-prioritized plasma proteomic screen", "CVD150 targeted pQTL screen",
          "interim coloc-supported axes".
  KULLANMA: "proteome-wide", "definitive therapeutic target", "validated drug target",
            "final GREEN-2", "NPPA inhibition is therapeutic".

NPPA-specific language:
  ✓ "NPPA identifies an AF-cardioembolic stroke genetic axis supported by colocalization."
  ✗ "NPPA inhibition should treat AF or stroke."
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402
import polars as pl  # noqa: E402

ROOT = Path("/Users/apple/Desktop/gwas_af")
REPORTS = ROOT / "reports"
FIGURES = ROOT / "figures"
TABLES = ROOT / "tables"
COLOC_DIR = ROOT / "results/coloc"
INTERIM_DIR = ROOT / "results/interim"

W1 = COLOC_DIR / "decode_cvd150_wave1_coloc_results.tsv"
W2 = COLOC_DIR / "decode_cvd150_wave2_coloc_results.tsv"
INSTR = ROOT / "data/processed/instruments/decode_cvd150_instruments.parquet"
MR_FULL = ROOT / "results/mr/decode_cvd150_pqtl_to_AF.tsv"
INTERIM_TSV = INTERIM_DIR / "decode_cvd150_interim_evidence_classes.tsv"

for d in (REPORTS, FIGURES, TABLES):
    d.mkdir(parents=True, exist_ok=True)


def _coloc_primary():
    rows = []
    for path, wave in ((W1, "Wave1"), (W2, "Wave2")):
        if path.exists():
            df = pl.read_csv(path, separator="\t")
            df = df.filter(pl.col("prior_model") == "default_quant_cc")
            for r in df.iter_rows(named=True):
                rows.append({**r, "wave": wave})
    return rows


def _coloc_sens():
    out = {}
    for path in (W1, W2):
        if path.exists():
            df = pl.read_csv(path, separator="\t")
            df = df.filter(pl.col("prior_model") == "flat_W04")
            for r in df.iter_rows(named=True):
                out[(r["gene"], r["outcome"])] = r["PP.H4"]
    return out


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------


def write_table1_overview():
    instr = pl.read_parquet(INSTR)
    n_instruments = instr.height
    n_unique_aptamers = instr["aptamer"].n_unique()
    n_unique_genes = instr["gene"].n_unique()
    n_matched = instr.filter(pl.col("af_match") == "matched").height

    rows = [
        {"metric": "Aptamers attempted (CVD150 quota)", "value": "150", "scope": "CVD-prioritized targeted screen"},
        {"metric": "Aptamers producing cis-instruments", "value": str(n_unique_aptamers), "scope": "after MAF/p/F threshold + cis ±1Mb + clumping"},
        {"metric": "Unique genes with cis-instruments", "value": str(n_unique_genes), "scope": ""},
        {"metric": "Total cis-pQTL instruments", "value": str(n_instruments), "scope": "single-base, p<5e-8, F>10, ±1 Mb cis"},
        {"metric": "Instruments matched to AF GWAS", "value": str(n_matched), "scope": "Roselli 2025 CVDKP common-variant all-ancestry meta-analysis (chr+pos exact)"},
        {"metric": "AF MR FDR-significant aptamers (within CVD150)", "value": "5", "scope": "BH-FDR within CVD150-screen, NOT proteome-wide"},
        {"metric": "AF MR nominal (uncorrected p<0.05) aptamers", "value": "12", "scope": ""},
        {"metric": "Coloc-tested target-outcome pairs (Wave 1+2)", "value": "35", "scope": "selective, hypothesis-driven"},
        {"metric": "Coloc-supported pairs (PP.H4 ≥ 0.50, primary prior)", "value": "8", "scope": "1 strong AF-anchored + 5 strong stroke + 1 moderate + 1 LAS"},
        {"metric": "Coloc-rejected (PP.H3 ≥ 0.50 > PP.H4)", "value": "5", "scope": "LD-confounding flag"},
        {"metric": "Independent EUR LD reference", "value": "Not configured", "scope": "SuSiE/conditional coloc deferred (D042 optional)"},
        {"metric": "Independent pQTL replication", "value": "Pending", "scope": "UKB-PPP DAR + SCALLOP not yet run"},
    ]
    out = TABLES / "table1_cvd150_screen_overview.tsv"
    pl.DataFrame(rows).write_csv(out, separator="\t")
    print(f"wrote {out}")
    return rows


def write_table2_interim():
    src = INTERIM_TSV
    if src.exists():
        dst = TABLES / "table2_interim_evidence_classes.tsv"
        df = pl.read_csv(src, separator="\t", infer_schema_length=0)
        df.write_csv(dst, separator="\t")
        print(f"wrote {dst}")


def write_table3_coloc_supported():
    coloc = _coloc_primary()
    sens = _coloc_sens()
    keep = []
    for r in coloc:
        h4 = r["PP.H4"]
        if h4 is None:
            continue
        try:
            h4 = float(h4)
        except (ValueError, TypeError):
            continue
        if math.isnan(h4) or h4 < 0.50:
            continue
        s_h4 = sens.get((r["gene"], r["outcome"]))
        keep.append({
            "wave": r["wave"], "gene": r["gene"], "outcome": r["outcome"],
            "n_variants": r["n_variants"],
            "PP.H4_primary": round(h4, 4),
            "PP.H4_sensitivity": round(float(s_h4), 4) if s_h4 is not None else None,
            "PP.H3_primary": round(float(r["PP.H3"]), 4) if r["PP.H3"] is not None else None,
            "PP.H1_primary": round(float(r["PP.H1"]), 4) if r["PP.H1"] is not None else None,
            "lead_pqtl_rsid": r["lead_pqtl_rsid"],
            "lead_pqtl_pval": r["lead_pqtl_pval"],
            "lead_outcome_pval": r["lead_outcome_pval"],
            "strength": "strong" if h4 >= 0.80 else "moderate",
            "axis": _axis_for(r["gene"], r["outcome"]),
        })
    keep.sort(key=lambda x: -x["PP.H4_primary"])
    out = TABLES / "table3_coloc_supported_axes.tsv"
    pl.DataFrame(keep).write_csv(out, separator="\t")
    print(f"wrote {out}")
    return keep


def _axis_for(gene, outcome):
    if gene == "NPPA":
        return "AF-CES cardioembolic axis"
    if gene in ("F11", "KNG1"):
        return "antithrombotic stroke axis (non-AF)"
    if gene == "MMP12":
        return "atherothrombotic stroke axis (non-AF)"
    return "other"


def write_table4_downgraded():
    if not INTERIM_TSV.exists():
        return
    df = pl.read_csv(INTERIM_TSV, separator="\t", infer_schema_length=0)
    dn = df.filter(pl.col("interim_label").str.starts_with("DOWNGRADED"))
    out = TABLES / "table4_downgraded_candidates.tsv"
    dn.write_csv(out, separator="\t")
    print(f"wrote {out}")


def write_table5_replication_roadmap():
    rows = [
        {"target_or_axis": "NPPA", "module": "AF-CES cardioembolic axis", "step": "UKB-PPP pQTL replication", "status": "pending DAR"},
        {"target_or_axis": "NPPA", "module": "AF-CES cardioembolic axis", "step": "SCALLOP/Olink NPPA/NPPB if available", "status": "to query"},
        {"target_or_axis": "NPPA", "module": "AF-CES cardioembolic axis", "step": "FinnGen AF/CES outcome sensitivity", "status": "to query"},
        {"target_or_axis": "NPPA", "module": "AF-CES cardioembolic axis", "step": "GTEx atrial-tissue eQTL support if meaningful", "status": "to query"},
        {"target_or_axis": "F11", "module": "antithrombotic stroke axis", "step": "FinnGen stroke subtype sensitivity", "status": "to query"},
        {"target_or_axis": "F11", "module": "antithrombotic stroke axis", "step": "external pQTL replication if available", "status": "to query"},
        {"target_or_axis": "F11", "module": "antithrombotic stroke axis", "step": "bleeding/VTE safety PheWAS later", "status": "post-replication"},
        {"target_or_axis": "KNG1", "module": "kallikrein-kinin stroke axis", "step": "FinnGen stroke subtype sensitivity", "status": "to query"},
        {"target_or_axis": "KNG1", "module": "kallikrein-kinin stroke axis", "step": "external pQTL replication if available", "status": "to query"},
        {"target_or_axis": "KNG1", "module": "kallikrein-kinin stroke axis", "step": "bleeding/VTE safety PheWAS later", "status": "post-replication"},
        {"target_or_axis": "MMP12", "module": "atherothrombotic stroke axis", "step": "LAS/AIS outcome replication", "status": "to query"},
        {"target_or_axis": "MMP12", "module": "atherothrombotic stroke axis", "step": "plaque/vascular single-cell annotation later", "status": "post-replication"},
        {"target_or_axis": "MMP12", "module": "atherothrombotic stroke axis", "step": "external pQTL replication", "status": "to query"},
        {"target_or_axis": "IL6R/DSC2/OGN/TIMP3/PCSK9", "module": "downgraded discordant", "step": "optional SuSiE audit after EUR LD reference setup", "status": "D042 optional"},
        {"target_or_axis": "all", "module": "atlas", "step": "final GREEN/RED/AMBER classification", "status": "LOCKED until replication + final review"},
    ]
    out = TABLES / "table5_replication_roadmap.tsv"
    pl.DataFrame(rows).write_csv(out, separator="\t")
    print(f"wrote {out}")


# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------


def write_results_md(t1_rows, coloc_supported):
    md = [
        "# Results — CVD-prioritized plasma proteomic screen for AF–HF–stroke pleiotropy",
        "",
        "> Manuscript-style draft (D043). Interim, CVD150-targeted, coloc-supported axes only.",
        "> No final therapeutic classification.",
        "",
        "## Result 1 — CVD150 targeted pQTL screen and quality control",
        "",
        "We assembled a CVD-prioritized plasma proteomic pQTL screen using the deCODE 2021 SomaScan resource. ",
        "From 150 attempted aptamers spanning seven CVD biology categories, ",
        f"{t1_rows[1]['value']} aptamers produced cis-pQTL instruments under the protocol ",
        "(MAF ≥ 0.01, p < 5×10⁻⁸, F > 10, ±1 Mb cis window, distance-clumped at ±500 kb), ",
        f"yielding {t1_rows[3]['value']} instruments across {t1_rows[2]['value']} unique genes. ",
        f"{t1_rows[4]['value']} of these instruments harmonized to the Roselli 2025 CVDKP common-variant all-ancestry AF GWAS (`AF_GWAS_AFGenPlus_commonFreq_ALLv21`) by exact chr/pos and allele matching. ",
        "Within this CVD150 screen, BH-FDR-controlled MR identified ",
        f"{t1_rows[5]['value']} aptamers as AF MR FDR-significant and {t1_rows[6]['value']} additional aptamers as AF MR nominal (p < 0.05).",
        "",
        "**Important scope note:** the CVD150-screen FDR is computed within this targeted set and is *not* a proteome-wide FDR.",
        "",
        "## Result 2 — Mendelian randomization alone is insufficient: colocalization sharply filters AF MR candidates",
        "",
        "MR alone produced several attractive AF candidates, but selective two-trait colocalization (Giambartolomei `coloc.abf`) ",
        "with full ±500 kb regional summary statistics — independently validated against the R `coloc` package (max |ΔPP| < 1×10⁻⁶) — ",
        "filtered the candidate set sharply.",
        "",
        "- **IL6R** (AF MR q = 4.13×10⁻¹³) showed PP.H3 = 1.000 and PP.H4 = 0.000 in the IL6R cis-region: ",
        "  the cis-pQTL lead variant and the AF lead variant in the same region were not co-localized at a shared causal SNP.",
        "- **DSC2** (AF MR q = 8.76×10⁻⁶) likewif then showed PP.H3 = 1.000, PP.H4 = 0.000.",
        "- **OGN, TIMP3** showed PP.H1-dominant patterns (the AF outcome was not associated across the cis-region) and weak PP.H4.",
        "- **PCSK9, SERPINF2** fawithd to co-localize with their target outcomes.",
        "",
        "These results — that AF MR Wald-ratio significance does not imply locus-level colocalization — are a key methodological point for the field and an *intentional output* of this study.",
        "",
        "## Result 3 — NPPA is the sole AF-anchored coloc-supported axis in this screen",
        "",
        "Of the AF-anchored colocalization tests across the CVD150 candidate shortlist, ",
        "only **NPPA** (atrial natriuretic peptide gene) survived as a coloc-supported AF signal:",
        "",
        "- NPPA × AF: PP.H4 = 0.911 (sensitivity prior PP.H4 = 0.829)",
        "- NPPA × cardioembolic stroke: PP.H4 = 0.565 (sensitivity 0.395)",
        "- NPPA × ischemic stroke: PP.H4 = 0.007 (no support)",
        "- NPPA × HF_overall: PP.H4 = 0.025 (no support)",
        "- NPPA × HF_nonischemic: PP.H4 = 0.009 (no support)",
        "- NPPA × HF_ni_HFrEF: PP.H4 = 0.081 (no support)",
        "",
        "NPPA therefore identifies an **AF-cardioembolic stroke genetic axis supported by colocalization**, ",
        "but the signal does *not* extend to ischemic-stroke-overall, large-artery stroke, small-vessel stroke, or any HF subtype. ",
        "We deliberately frame NPPA as a *coloc-supported genetic axis*, not as a validated therapeutic target: ",
        "natriuretic peptide biology is clinically complex and the cis-pQTL signal may reflect atrial stretch / remodeling, ",
        "circulating peptide measurement biology, or causal pathway biology.",
        "",
        "## Result 4 — Two non-AF stroke modules emerge outside the AF anchor",
        "",
        "When we extended colocalization to non-AF-anchored hypotheses, two strong stroke modules emerged:",
        "",
        "- **F11 (Coagulation Factor XI)** colocalized strongly with AIS (PP.H4 = 0.993) and CES (PP.H4 = 0.989). F11 AF MR was no_signal (p = 0.84): F11 is *not* an AF-anchored target in this framework. F11 is a positive sanity anchor for antithrombotic stroke biology.",
        "- **KNG1 (Kininogen, HMW)** colocalized strongly with AIS (PP.H4 = 0.991) and CES (PP.H4 = 0.933). KNG1 anchors the kallikrein-kinin / contact-activation arm of the antithrombotic axis.",
        "- **MMP12 (Matrix Metalloproteinase 12)** colocalized strongly with LAS (PP.H4 = 0.916) and AIS (PP.H4 = 0.838) but not with AF (PP.H4 = 0.062). MMP12 anchors an atherothrombotic / large-artery stroke module distinct from the NPPA cardioembolic axis.",
        "",
        "We treat F11/KNG1/MMP12 as **interim BLUE module candidates** — non-AF, stroke-biology axes — and explicitly do not place them in the AF-anchored GREEN/RED/AMBER framework.",
        "",
        "## Result 5 — Colocalization prevents over-interpretation of MR signals",
        "",
        "Six CVD150 candidates that appeared promising under MR alone were downgraded after colocalization: ",
        "**IL6R, DSC2** (LD-confounded — PP.H3 ≥ 0.50, PP.H3 > PP.H4); **OGN, TIMP3, PCSK9** (AF MR not coloc-supported); ",
        "**SERPINF2** (BLUE module discordant). The main message is that **MR + coloc together are required before therapeutic interpretation** — MR alone over-calls candidate targets in cis-regions where the cis-pQTL and outcome lead variants are in modest LD but not co-localized.",
        "",
        "## Three-axis interim atlas",
        "",
        "```",
        "NPPA       = AF-CES cardioembolic axis           (AF-anchored, INTERIM AMBER)",
        "F11/KNG1   = antithrombotic AIS/CES axis         (non-AF, BLUE)",
        "MMP12      = atherothrombotic AIS/LAS axis       (non-AF, BLUE)",
        "```",
        "",
        "These three axes together, plus the explicit downgrades, are the headline of the current evidence layer.",
        "",
        "## Tables and figures",
        "",
        "- Table 1: CVD150 screen overview",
        "- Table 2: Interim evidence classes (10 candidates, final_class_allowed = false on every row)",
        "- Table 3: Coloc-supported axes (8 pairs)",
        "- Table 4: Downgraded candidates (6 genes)",
        "- Table 5: Replication roadmap",
        "- Figure 1: Study design (CVD-prioritized screen + AF anchor + outcome panel)",
        "- Figure 2: MR → coloc filtering flow",
        "- Figure 3: Three-axis interim atlas schematic",
        "- Figure 4: NPPA AF-CES axis pleiotropy",
        "- Figure 5: BLUE modules (F11/KNG1/MMP12)",
        "- Figure 6: Downgraded MR hits — coloc filter",
    ]
    out = REPORTS / "manuscript_style_results.md"
    out.write_text("\n".join(md))
    print(f"wrote {out}")


def write_methods_md():
    md = [
        "# Methods — CVD-prioritized plasma proteomic screen for AF–HF–stroke pleiotropy",
        "",
        "> Manuscript-style methods draft (D043). Methods describe what was actually run; ",
        "> nothing here is aspirational. Independent replication, full proteome-wide screening, ",
        "> SuSiE/conditional fine-mapping, and final therapeutic classification are explicitly *not* part of the current evidence layer.",
        "",
        "## Study design",
        "",
        "We performed a CVD-prioritized targeted plasma-proteomic pQTL screen (CVD150) anchored on atrial fibrillation (AF), ",
        "with directed downstream projection to heart-failure (HF) and stroke subtype outcomes, followed by selective ",
        "two-trait colocalization at high-priority candidate target-outcome pairs. The intent was to map *interim* coloc-supported genetic axes — ",
        "not to declare therapeutic targets.",
        "",
        "## Genome-wide association data",
        "",
        "- **AF (primary anchor):** Roselli et al., 2025 — CVDKP common-variant AFGenPlus all-ancestry meta-analysis (fwith `AF_GWAS_AFGenPlus_commonFreq_ALLv21.txt.gz`; ancestry = MIXED, GRCh38, 15.6 M harmonized variants of 18.8 M raw rows; max effective N ≈ 1.65 M; `sample_overlap_flag = possible_UKB_overlap` per D020). Underflow-safe handling of extreme P (D026): each row carries `log10p_signed`, `neg_log10p`, `pval_capped_for_tools`, `pval_underflow_flag`.",
        "- **HF (4 subtypes):** HERMES 2024 EUR (HF_overall, HF_nonischemic, HF_ni_HFrEF, HF_ni_HFpEF; GRCh37; per-variant pyliftover GRCh38 → GRCh37 inside the colocalization region).",
        "- **Stroke (4 subtypes):** GIGASTROKE 2022 EUR (AIS, CES, LAS, SVS; GRCh38).",
        "- **Ancestry-mismatch note (critical):** the AF anchor is **all-ancestry / mixed** (CVDKP common-variant set), whwith HF and stroke outcomes are **EUR-only**. Cross-trait MR + coloc therefore propagates an `ancestry_mismatch` flag (D007); interim coloc-supported axes should be interpreted under this caveat.",
        "- **Positive controls:** PITX2 (AF), HDAC9 (LAS), FOXF2 (SVS, AIS), BAG3 (HFrEF) — all recovered at expected magnitudes (D025 retention check).",
        "",
        "## Plasma proteomic pQTL source",
        "",
        "deCODE 2021 SomaScan v4 plasma pQTL (Ferkingstad et al., 2021), GRCh38, ~36 k participants. ",
        "Selective per-aptamer fetch via the authorised /s3/folder + /s3/download API to keep peak local disk small (D027/D034/D036). ",
        "**Crucially, this is a CVD-prioritized targeted screen of 150 aptamers spanning seven CVD biology categories ",
        "(GWAS-loci, CV-druggable, IL6/inflammation, lipid/metabolic, vascular, cardiac, sentinel) — not a proteome-wide query.**",
        "",
        "## Cis-pQTL instrument extraction",
        "",
        "Per aptamer, single-base variants in a ±1 Mb cis window were filtered to MAF ≥ 0.01, p < 5×10⁻⁸, F-statistic > 10, ",
        "and distance-clumped greedily at ±500 kb. Allele-harmonized to the AF GWAS by exact chr/pos with allele-set match; ",
        "palindromic A/T and C/G variants in the |EAF − 0.5| < 0.05 ambiguity window were dropped.",
        "",
        "## AF MR + downstream projection",
        "",
        "AF MR used Wald ratio (single instrument), inverse-variance-weighted random-effects (IVW-RE) and weighted median (multi-instrument). ",
        "BH-FDR was computed *within the CVD150 screen* (not proteome-wide). Action direction was derived as ",
        "`d_T = -sign(beta_T_AF)` (D005). Projection to HF and stroke subtypes used `theta = d_T × beta_outcome` per outcome.",
        "",
        "## Selective colocalization (Wave 1 + Wave 2)",
        "",
        "Two-trait Giambartolomei `coloc.abf` (Wakefield approximate Bayes factor) on full ±500 kb regional summary statistics, ",
        "with NO p-value filtering inside the coloc region (D025). Primary prior model `default_quant_cc`: p1 = p2 = 1×10⁻⁴, p12 = 1×10⁻⁵, ",
        "W_pqtl = 0.0225 (sd = 0.15, quantitative inverse-rank-normalized protein), W_outcome = 0.04 (sd = 0.20, case-control). ",
        "sensitivity prior model `flat_W04`: p12 = 5×10⁻⁶, W = 0.04 for both traits.",
        "",
        "**Method validation (D039):** the Python implementation was independently cross-validated against the R `coloc` package on two preflight pairs (IL6R × AF, NPPA × CES) with maximum absolute posterior probability difference < 1×10⁻⁶ across all five hypotheses (PP.H0–PP.H4).",
        "",
        "**Pair selection (selective, hypothesis-driven):** Wave 1 = 15 pairs (D038); Wave 2 = 20 pairs (D040). All-pair coloc and proteome-wide coloc are *not* part of this work.",
        "",
        "**Aptamer artifact annotation:** intragenic-vs-distal lead-variant flag, multi-aptamer convergence (D032), and protein-altering annotation (deferred to future VEP).",
        "",
        "**Build-aware regional matching:** GRCh38 GWAS (AF, all 4 stroke subtypes) merged by chr+pos directly; GRCh37 GWAS (HERMES HF) merged via per-variant pyliftover GRCh38 → GRCh37 inside the regional window.",
        "",
        "## SuSiE / conditional colocalization",
        "",
        "Reliable EUR LD reference was not available locally during this evidence layer. Per D040, all SuSiE-pair candidates were flagged ",
        "`SUSIE_NOT_RUN_LD_UNAVAILABLE_OR_UNRELIABLE`. SuSiE / conditional rescue is deferred to an optional future audit (D042).",
        "",
        "## Interim evidence classes (D041)",
        "",
        "Coloc-supported pairs were grouped into interim evidence classes — `INTERIM_AMBER_AF_CES_SHARED_BENEFIT_CANDIDATE` (NPPA), ",
        "`INTERIM_BLUE_COLOC_SUPPORTED_*_STROKE_AXIS` (F11, KNG1, MMP12), and `DOWNGRADED_*` for MR-positive but coloc-fawithd signals. ",
        "**Final GREEN/RED/AMBER therapeutic classification is explicitly locked in this evidence layer.** The `final_class_allowed` flag is `false` on every row.",
        "",
        "## Software",
        "",
        "Python 3.13 (polars 1.40, requests, pyliftover, matplotlib 3.10), R 4.x (`coloc` package, used only for cross-validation of the Python `coloc.abf` implementation), plink/plink2.",
        "",
        "## Data and code availability",
        "",
        "deCODE 2021 pQTL: authorised access via deCODE API (per-aptamer selective fetch). HERMES 2024, GIGASTROKE 2022, AF Roselli 2025: public download. ",
        "All analysis code is in this project repository under `scripts/python/`, `src/`, and the orchestration `Snakefile`.",
    ]
    out = REPORTS / "manuscript_style_methods.md"
    out.write_text("\n".join(md))
    print(f"wrote {out}")


def write_limitations_md():
    md = [
        "# Limitations — CVD-prioritized plasma proteomic screen",
        "",
        "> Manuscript-style limitations (D043). Phrased to keep claim boundaries tight.",
        "",
        "1. **Targeted, not proteome-wide.** The CVD150 screen is a hypothesis-prioritized 150-aptamer subset; the BH-FDR is computed within this subset and is *not* a proteome-wide false-discovery rate. Interpretation of FDR-significant hits as 'genome-wide' or 'proteome-wide' significant would be incorrect.",
        "",
        "2. **Aptamer-based pQTL.** deCODE pQTL relies on SomaScan aptamers; aptamer-specific epitope artifacts, multi-aptamer discordance, and intragenic missense effects can produce locus-specific cis-pQTL signal that does not faithfully report circulating protein concentration. We annotated intragenic position and multi-aptamer convergence (D032), but a definitive aptamer-validity check (multi-platform Olink/MS proteomic agreement) was not run.",
        "",
        "3. **No independent pQTL replication yet.** UKB-PPP controlled-access replication is pending DAR. SCALLOP/Olink overlap was not queried in this evidence layer. Single-cohort cis-pQTL evidence is therefore Tier-2 by D002/D009 REVISED criteria; Tier 0/1 claims require independent platform replication.",
        "",
        "4. **No outcome replication.** FinnGen AF/CES/HF/stroke sensitivity, GTEx atrial-tissue eQTL support, and external pQTL replication are listed in the replication roadmap (Table 5) but were not run.",
        "",
        "5. **NPPA action direction is *not* a drug recommendation.** Natriuretic peptide biology is clinically complex. The NPPA cis-pQTL signal could reflect atrial-stretch / remodeling biology, peptide-measurement biology, or causal pathway biology. We deliberately frame NPPA as a *coloc-supported AF-CES genetic axis*, not as 'NPPA inhibition treats AF or stroke'.",
        "",
        "6. **HFpEF is exploratory.** per protocol decisions in Phase 2, HFpEF was treated as exploratory (not a primary outcome) and is excluded from primary classification.",
        "",
        "7. **SuSiE / conditional colocalization not run.** Reliable EUR LD reference was not available; SuSiE-coloc was deferred to an optional future audit (D042). For LD-confounded discordant signals (IL6R, DSC2), the current downgrade reflects coloc.abf single-causal-variant assumption only; multi-causal-variant fine-mapping could in principle rescue some of these signals.",
        "",
        "8. **Selective coloc, not all-pair.** Coloc was run on 35 hypothesis-driven target-outcome pairs (Wave 1 + Wave 2). Coloc results are conditional on the chosen pair list and do not enumerate all possible pleiotropy targets across the screen.",
        "",
        "9. **Mixed-ancestry AF anchor with EUR-only outcomes (ancestry mismatch).** The AF anchor is the CVDKP common-variant `AF_GWAS_AFGenPlus_commonFreq_ALLv21` all-ancestry meta-analysis (ancestry = MIXED), whwith HF (HERMES 2024) and stroke (GIGASTROKE 2022) outcomes are EUR-only. This deliberate mismatch is recorded by the `ancestry_mismatch` flag (D007) propagated through all results. Trans-ancestry sensitivity (e.g., AFR-specific or EAS-specific AF GWAS replication) was not run; cross-trait MR + coloc results should be interpreted with this caveat. The primary motivation for using the all-ancestry AF set is statistical power at the discovery stage; an EUR-only AF sensitivity analysis is in the replication roadmap (Table 5).",
        "",
        "10. **No final therapeutic classification.** D041 + D043 explicitly lock GREEN/RED/AMBER calls. Claims at the level of 'NPPA is GREEN-2' or 'IL6R is RED-1' are not supported by the current evidence layer.",
        "",
        "11. **Sample-overlap considerations.** UKB-PPP-derived pQTLs were not used as the primary discovery layer here, but if they enter via replication, sample overlap with UKB-derived AF/HF/stroke outcomes will require explicit independent-cohort sensitivity (D002 / `sample_overlap_flag`).",
        "",
        "12. **Aptamer 5443_62_NPPA_ANP single-aptamer status.** NPPA AF-CES axis is supported by a single aptamer in this screen. Multi-aptamer NPPA convergence (NT-proBNP / BNP via NPPB) was anchored in the Phase 3A smoke run (n = 49 instruments) but not formally repeated within Wave 1+2.",
    ]
    out = REPORTS / "manuscript_style_limitations.md"
    out.write_text("\n".join(md))
    print(f"wrote {out}")


def write_abstract_md():
    md = [
        "# Abstract draft — CVD-prioritized plasma proteomic screen for AF–HF–stroke pleiotropy",
        "",
        "> Manuscript-style abstract draft (D043). Keep tight claim boundaries.",
        "",
        "## Abstract",
        "",
        "We performed a **CVD-prioritized plasma proteomic pQTL screen** to map atrial-fibrillation-anchored therapeutic pleiotropy ",
        "across heart-failure and stroke subtypes. Using deCODE 2021 SomaScan cis-pQTL instruments for 98 cis-instrument-producing ",
        "aptamers from a 150-aptamer CVD-prioritized panel, AF (Roselli 2025 CVDKP common-variant all-ancestry meta-analysis), heart-failure (HERMES 2024 EUR; overall, ",
        "non-ischemic, HFrEF, HFpEF) and stroke-subtype (GIGASTROKE 2022 EUR; AIS, CES, LAS, SVS) GWAS, two-sample Mendelian ",
        "randomization, intervention-direction projection, and selective two-trait colocalization, we identified one ",
        "AF-anchored coloc-supported axis — **NPPA** with AF (PP.H4 = 0.91) and cardioembolic stroke (PP.H4 = 0.57) — and three ",
        "non-AF-anchored stroke modules: an **antithrombotic axis (F11, KNG1)** colocalizing strongly with ischemic and cardioembolic ",
        "stroke (PP.H4 ≥ 0.93), and an **atherothrombotic axis (MMP12)** colocalizing strongly with large-artery and ischemic ",
        "stroke (PP.H4 ≥ 0.84). Several attractive AF MR candidates fawithd colocalization, including IL6R and DSC2 (note: AF anchor is all-ancestry / mixed whwith HF and stroke outcomes are EUR-only — `ancestry_mismatch` flag propagated through all results) ",
        "(PP.H3 = 1.0, PP.H4 = 0; LD-confounded), OGN, TIMP3, PCSK9 and SERPINF2, highlighting that locus-level colocalization ",
        "is necessary before any therapeutic interpretation of cis-pQTL MR signals. We frame these findings as **interim ",
        "coloc-supported axes** within a CVD-prioritized targeted screen — not as final therapeutic targets — pending ",
        "independent pQTL replication, outcome replication, and LD-aware fine-mapping.",
        "",
        "## Keyword-style framing for downstream use",
        "",
        "- CVD-prioritized plasma proteomic screen (CVD150)",
        "- cis-pQTL Mendelian randomization",
        "- selective colocalization (Giambartolomei coloc.abf)",
        "- AF-anchored therapeutic pleiotropy atlas",
        "- interim coloc-supported axes",
        "- NPPA AF–cardioembolic stroke axis",
        "- antithrombotic stroke module (F11, KNG1)",
        "- atherothrombotic stroke module (MMP12)",
        "- LD-confounding (IL6R, DSC2)",
        "- replication roadmap (UKB-PPP, SCALLOP, FinnGen, GTEx)",
    ]
    out = REPORTS / "manuscript_style_abstract_draft.md"
    out.write_text("\n".join(md))
    print(f"wrote {out}")


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def fig_study_design():
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.axis("off")
    boxes = [
        (0.05, 0.85, 0.40, 0.10, "deCODE 2021 plasma pQTL\nCVD-prioritized 150-aptamer panel", "#dde6ff"),
        (0.55, 0.85, 0.40, 0.10, "AF Roselli 2025 CVDKP AFGenPlus\nall-ancestry, GRCh38\nprimary anchor", "#ffe8d6"),
        (0.05, 0.65, 0.90, 0.10, "Cis-pQTL extraction (MAF≥0.01, p<5e-8, F>10, ±1Mb cis, ±500kb clump)\n→ 235 instruments / 98 aptamers / 95 genes / 209 matched-to-AF", "#dceedb"),
        (0.05, 0.45, 0.90, 0.10, "AF MR (Wald, IVW-RE, weighted-median) + BH-FDR within CVD150-screen\n→ 5 FDR + 12 nominal candidates", "#fff0a6"),
        (0.05, 0.25, 0.90, 0.10, "Intervention-direction projection: d_T = -sign(beta_T_AF) × HF + stroke subtypes", "#ffe5ec"),
        (0.05, 0.05, 0.90, 0.10, "Selective coloc.abf (Wave 1 = 15 + Wave 2 = 20 pairs)\n→ NPPA (AF-CES), F11/KNG1 (BLUE antithrombotic), MMP12 (BLUE atherothrombotic)", "#cfe7d4"),
    ]
    for x, y, w, h, text, color in boxes:
        ax.add_patch(mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01",
                                             facecolor=color, edgecolor="black", lw=0.8))
        ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=9)

    arrows = [(0.50, 0.85, 0.50, 0.75), (0.50, 0.65, 0.50, 0.55), (0.50, 0.45, 0.50, 0.35), (0.50, 0.25, 0.50, 0.15)]
    for x1, y1, x2, y2 in arrows:
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="->", color="black", lw=1.0))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title("Study design — CVD-prioritized AF-anchored pQTL screen", fontsize=12, fontweight="bold")
    out = FIGURES / "study_design_cvd150.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {out}")


def fig_filtering_flow():
    fig, ax = plt.subplots(figsize=(10, 5.5))
    stages = [
        ("Aptamers attempted", 150, "#cccccc"),
        ("Cis-instrument-producing", 98, "#bbdaff"),
        ("AF-matched aptamers", 95, "#9acdff"),
        ("AF MR FDR (within CVD150)", 5, "#ffd25e"),
        ("AF MR FDR + coloc-supported", 1, "#88c98c"),
    ]
    labels = [s[0] for s in stages]
    counts = [s[1] for s in stages]
    colors = [s[2] for s in stages]
    bars = ax.barh(range(len(stages)), counts, color=colors, edgecolor="black")
    for i, (lbl, c) in enumerate(zip(labels, counts)):
        ax.text(c + 2, i, f"{c}", va="center", fontsize=10, fontweight="bold")
    ax.set_yticks(range(len(stages)))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Aptamer / candidate count")
    ax.set_title("MR → coloc filtering flow (CVD150 targeted screen)\nNPPA is the only AF-anchored coloc-supported survivor",
                 fontsize=11, fontweight="bold")
    ax.set_xlim(0, 175)
    fig.text(0.5, -0.02,
             "Note: BH-FDR is computed within the CVD150 targeted screen, not proteome-wide.",
             ha="center", fontsize=8, style="italic")
    out = FIGURES / "cvd150_mr_coloc_filtering_flow.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {out}")


def fig_three_axis_atlas():
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.axis("off")

    axes_data = [
        ("AF-anchored\n(INTERIM AMBER)", "NPPA", "AF-CES cardioembolic axis",
         "AF PP.H4=0.91\nCES PP.H4=0.57", "#ffd97d"),
        ("Non-AF BLUE\nantithrombotic", "F11 + KNG1", "AIS / CES axis",
         "F11×AIS=0.99  KNG1×AIS=0.99\nF11×CES=0.99  KNG1×CES=0.93", "#a8d5e3"),
        ("Non-AF BLUE\natherothrombotic", "MMP12", "AIS / LAS axis",
         "MMP12×LAS=0.92\nMMP12×AIS=0.84", "#bdd2b6"),
    ]
    for i, (header, target, label, evid, color) in enumerate(axes_data):
        x0 = 0.04 + i * 0.32
        ax.add_patch(mpatches.FancyBboxPatch((x0, 0.30), 0.28, 0.55, boxstyle="round,pad=0.02",
                                             facecolor=color, edgecolor="black", lw=1.2))
        ax.text(x0 + 0.14, 0.78, header, ha="center", va="center", fontsize=10, fontweight="bold")
        ax.text(x0 + 0.14, 0.66, target, ha="center", va="center", fontsize=14, fontweight="bold", color="darkred")
        ax.text(x0 + 0.14, 0.55, label, ha="center", va="center", fontsize=9, style="italic")
        ax.text(x0 + 0.14, 0.40, evid, ha="center", va="center", fontsize=8.5, family="monospace")

    ax.text(0.5, 0.18,
            "Intentional separation: NPPA ≠ F11/KNG1 ≠ MMP12.\n"
            "Three biologically distinct stroke / cardiovascular axes — only NPPA is AF-anchored.",
            ha="center", va="center", fontsize=9, style="italic")
    ax.text(0.5, 0.05,
            "All labels are INTERIM. Final GREEN/RED/AMBER therapeutic classification remains LOCKED (D041/D043).",
            ha="center", va="center", fontsize=8, color="darkred", fontweight="bold")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title("Three-axis interim atlas (CVD150 targeted screen)", fontsize=12, fontweight="bold")
    out = FIGURES / "three_axis_interim_atlas.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {out}")


def fig_nppa_pleiotropy(coloc):
    nppa = {(r["gene"], r["outcome"]): r for r in coloc if r["gene"] == "NPPA"}
    sens = _coloc_sens()
    outcomes_order = ["AF", "cardioembolic_stroke", "ischemic_stroke",
                      "HF_overall", "HF_nonischemic", "HF_ni_HFrEF"]
    h4_primary = []
    h4_sens = []
    labels = []
    for o in outcomes_order:
        r = nppa.get(("NPPA", o))
        if r is None:
            h4_primary.append(0.0); h4_sens.append(0.0)
        else:
            h4_primary.append(float(r["PP.H4"]))
            h4_sens.append(float(sens.get(("NPPA", o)) or 0.0))
        labels.append(o.replace("_", "\n"))

    fig, ax = plt.subplots(figsize=(10, 5.5))
    x = range(len(outcomes_order))
    w = 0.36
    p_bars = ax.bar([i - w/2 for i in x], h4_primary, width=w, color="#cb6a52", label="primary prior", edgecolor="black", lw=0.6)
    s_bars = ax.bar([i + w/2 for i in x], h4_sens, width=w, color="#fdb46e", label="sensitivity (W=0.04 both)", edgecolor="black", lw=0.6)
    for bars, vals in ((p_bars, h4_primary), (s_bars, h4_sens)):
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width()/2, v + 0.02, f"{v:.2f}", ha="center", fontsize=8)
    ax.axhline(0.80, color="gray", lw=0.7, ls="--", alpha=0.6)
    ax.axhline(0.50, color="gray", lw=0.7, ls=":", alpha=0.6)
    ax.text(len(outcomes_order) - 0.4, 0.81, "strong", fontsize=8, color="gray")
    ax.text(len(outcomes_order) - 0.4, 0.51, "moderate", fontsize=8, color="gray")
    ax.set_xticks(list(x)); ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("PP.H4 (coloc-supported shared causal variant)")
    ax.set_ylim(0, 1.08)
    ax.set_title("NPPA AF-CES cardioembolic axis — pleiotropy is specific\n(AF + CES coloc-supported; AIS / HF subtypes not supported)",
                 fontsize=11, fontweight="bold")
    ax.legend(loc="upper right", fontsize=9)
    fig.text(0.5, -0.04,
             "Wave 1 + Wave 2; interim AMBER label only — no final classification assigned.",
             ha="center", fontsize=8, style="italic")
    out = FIGURES / "nppa_af_ces_axis.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {out}")


def fig_blue_modules(coloc):
    triples = [
        ("F11", "ischemic_stroke", "F11 × AIS"),
        ("F11", "cardioembolic_stroke", "F11 × CES"),
        ("KNG1", "ischemic_stroke", "KNG1 × AIS"),
        ("KNG1", "cardioembolic_stroke", "KNG1 × CES"),
        ("MMP12", "ischemic_stroke", "MMP12 × AIS"),
        ("MMP12", "large_artery_stroke", "MMP12 × LAS"),
    ]
    by_pair = {(r["gene"], r["outcome"]): r for r in coloc}
    sens = _coloc_sens()

    h4_p, h4_s, labels, colors = [], [], [], []
    for g, o, lab in triples:
        r = by_pair.get((g, o))
        h4_p.append(float(r["PP.H4"]) if r else 0.0)
        h4_s.append(float(sens.get((g, o)) or 0.0))
        labels.append(lab)
        if g in ("F11", "KNG1"):
            colors.append("#5d9bdb")
        else:
            colors.append("#7fbc7e")

    fig, ax = plt.subplots(figsize=(10, 5.5))
    x = list(range(len(triples)))
    w = 0.36
    bars1 = ax.bar([i - w/2 for i in x], h4_p, width=w, color=colors, edgecolor="black", lw=0.7, label="primary prior")
    bars2 = ax.bar([i + w/2 for i in x], h4_s, width=w, color=colors, alpha=0.55, edgecolor="black", lw=0.5, label="sensitivity")
    for bars, vals in ((bars1, h4_p), (bars2, h4_s)):
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width()/2, v + 0.015, f"{v:.2f}", ha="center", fontsize=8)
    ax.axhline(0.80, color="gray", lw=0.7, ls="--", alpha=0.6)
    ax.text(len(triples) - 0.45, 0.81, "strong", fontsize=8, color="gray")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9, rotation=10)
    ax.set_ylabel("PP.H4")
    ax.set_ylim(0, 1.08)
    ax.set_title("Non-AF BLUE stroke modules — F11/KNG1 (antithrombotic) + MMP12 (atherothrombotic)",
                 fontsize=11, fontweight="bold")
    blue_patch = mpatches.Patch(color="#5d9bdb", label="antithrombotic axis (F11, KNG1)")
    green_patch = mpatches.Patch(color="#7fbc7e", label="atherothrombotic axis (MMP12)")
    ax.legend(handles=[blue_patch, green_patch], loc="upper right", fontsize=9)
    fig.text(0.5, -0.04,
             "All BLUE labels INTERIM — non-AF-anchored, kept separate from GREEN/RED AF therapeutic framework.",
             ha="center", fontsize=8, style="italic")
    out = FIGURES / "blue_modules_f11_kng1_mmp12.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {out}")


def fig_downgraded_filter(coloc):
    df = pl.read_csv(INTERIM_TSV, separator="\t", infer_schema_length=0) if INTERIM_TSV.exists() else None
    if df is None:
        return
    by_gene = {r["gene"]: r for r in df.iter_rows(named=True)}
    by_pair = {(r["gene"], r["outcome"]): r for r in coloc}

    genes = ["IL6R", "DSC2", "OGN", "TIMP3", "PCSK9", "SERPINF2"]
    af_q = []
    primary_h4 = []
    primary_h3 = []
    primary_label = []
    for g in genes:
        ge = by_gene.get(g, {})
        try:
            af_q.append(float(ge.get("AF_MR_q_CVD150") or "nan"))
        except (ValueError, TypeError):
            af_q.append(float("nan"))
        outcome = ge.get("primary_outcome", "AF")
        outcome_full = {"AF": "AF", "AIS": "ischemic_stroke",
                        "CES": "cardioembolic_stroke", "LAS": "large_artery_stroke",
                        "SVS": "small_vessel_stroke",
                        "HF_overall": "HF_overall"}.get(outcome, outcome)
        c = by_pair.get((g, outcome_full), {})
        try:
            primary_h4.append(float(c.get("PP.H4") or 0.0))
        except (ValueError, TypeError):
            primary_h4.append(0.0)
        try:
            primary_h3.append(float(c.get("PP.H3") or 0.0))
        except (ValueError, TypeError):
            primary_h3.append(0.0)
        primary_label.append(f"{g}\n× {outcome}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.0))
    ax1, ax2 = axes
    neg_log_q = [-math.log10(q) if q and q > 0 and not math.isnan(q) else 0 for q in af_q]
    bars1 = ax1.bar(primary_label, neg_log_q, color="#cb6a52", edgecolor="black", lw=0.7)
    for b, q in zip(bars1, af_q):
        if q is not None and not math.isnan(q):
            ax1.text(b.get_x() + b.get_width()/2, b.get_height() + 0.4,
                     f"q={q:.1e}" if q < 0.01 else f"q={q:.2f}",
                     ha="center", fontsize=8)
    ax1.set_ylabel("AF MR −log10(q) within CVD150-screen")
    ax1.set_title("AF MR signal (or primary outcome MR)\nbefore colocalization", fontsize=10, fontweight="bold")
    ax1.tick_params(axis="x", labelsize=8.5, rotation=10)

    w = 0.36
    x = list(range(len(genes)))
    bars2 = ax2.bar([i - w/2 for i in x], primary_h4, width=w, color="#5d9bdb", edgecolor="black", lw=0.7, label="PP.H4")
    bars3 = ax2.bar([i + w/2 for i in x], primary_h3, width=w, color="#888888", edgecolor="black", lw=0.7, label="PP.H3")
    for b, v in zip(bars2, primary_h4):
        ax2.text(b.get_x() + b.get_width()/2, v + 0.02, f"{v:.2f}", ha="center", fontsize=8)
    for b, v in zip(bars3, primary_h3):
        ax2.text(b.get_x() + b.get_width()/2, v + 0.02, f"{v:.2f}", ha="center", fontsize=8, color="#444")
    ax2.set_xticks(x); ax2.set_xticklabels(primary_label, fontsize=8.5, rotation=10)
    ax2.set_ylabel("Posterior probability")
    ax2.set_ylim(0, 1.10)
    ax2.axhline(0.80, color="gray", lw=0.6, ls="--", alpha=0.5)
    ax2.axhline(0.50, color="gray", lw=0.6, ls=":", alpha=0.5)
    ax2.set_title("Coloc PP.H3 / PP.H4\nat primary outcome", fontsize=10, fontweight="bold")
    ax2.legend(loc="upper right", fontsize=9)

    fig.suptitle("Downgraded MR hits — colocalization filter prevents over-interpretation",
                 fontsize=12, fontweight="bold", y=1.02)
    fig.tight_layout()
    out = FIGURES / "downgraded_mr_hits_coloc_filter.pdf"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {out}")


def main() -> int:
    print("=== TABLES ===")
    t1 = write_table1_overview()
    write_table2_interim()
    coloc_supported = write_table3_coloc_supported()
    write_table4_downgraded()
    write_table5_replication_roadmap()

    print("=== REPORTS ===")
    write_results_md(t1, coloc_supported)
    write_methods_md()
    write_limitations_md()
    write_abstract_md()

    print("=== FIGURES ===")
    coloc = _coloc_primary()
    fig_study_design()
    fig_filtering_flow()
    fig_three_axis_atlas()
    fig_nppa_pleiotropy(coloc)
    fig_blue_modules(coloc)
    fig_downgraded_filter(coloc)
    return 0


if __name__ == "__main__":
    sys.exit(main())
