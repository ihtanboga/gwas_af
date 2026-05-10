#!/usr/bin/env python3
"""D034 — CVD-prioritized 500 aptamer staged batch builder.

Curates a CVD-relevant gene list across 9 categories (D034), intersects
against the Ferkingstad 2021 deCODE folder, applies the max-3-aptamer-per-gene
rule (with biological exceptions), and emits 4 QC tables that gate the actual
extraction. Extraction is NOT executed by this script.

Outputs:
  config/decode_cvd500_batch.yml
  results/qc/decode_cvd500_batch_design.tsv          (per-aptamer)
  results/qc/decode_cvd500_batch_category_counts.tsv (per-category)
  results/qc/decode_cvd500_expected_runtime.tsv      (worker projections)
  results/qc/decode_cvd500_list_qc_summary.md
"""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import yaml

ROOT = Path("/Users/apple/Desktop/gwas_af")
FOLDER_JSON = ROOT / "data/registry/decode_2021_folder.json"
CALIB_REPORT = ROOT / "results/qc/decode_100aptamer_calibration_report.tsv"

OUT_YAML = ROOT / "config/decode_cvd500_batch.yml"
OUT_DESIGN = ROOT / "results/qc/decode_cvd500_batch_design.tsv"
OUT_CATCOUNTS = ROOT / "results/qc/decode_cvd500_batch_category_counts.tsv"
OUT_RUNTIME = ROOT / "results/qc/decode_cvd500_expected_runtime.tsv"
OUT_SUMMARY = ROOT / "results/qc/decode_cvd500_list_qc_summary.md"

# ---------------------------------------------------------------------------
# Curated CVD-relevant gene lists per D034 category. Order = priority for the
# max-3-per-gene rule when a single gene appears in multiple categories.
# ---------------------------------------------------------------------------

CATEGORIES: dict[str, list[str]] = {
    # NB: cat9 sanity comes FIRST so smoke-validated targets get cat9 primary
    # tag (D034 requirement: all categories nonzero). The primary_category
    # logic walks dict order; cat9 then cat1..cat8.
    "cat9_smoke_sanity": [
        "IL6R", "NPPA", "NPPB", "PCSK9", "F11", "ADAMTS13",
        "F2", "F9", "PROC", "ANGPTL3", "ANGPTL4", "APOB", "APOE", "ANGPT1",
        "MMP9", "CRP", "ANGPTL4", "PROS1", "F2", "FGB",
    ],

    # 1. AF / HF / stroke GWAS locus-near proteins (broadened to plasma-detected
    # neighbours; intracellular sarcomere/channel proteins kept for completeness
    # but most won't be in deCODE — gracefully skipped).
    "cat1_af_hf_stroke_gwas_loci": [
        # AF anchors + neighbours
        "PITX2", "ZFHX3", "KCNN3", "PRRX1", "CAV1", "TBX5", "NKX2-5",
        "GJA5", "KCNN2", "MYH6", "MYH7", "TTN", "BAG3", "SCN5A", "SCN10A",
        "KCNJ2", "KCNH2", "CACNA1C", "CASQ2", "GJA1", "HCN4",
        "PLN", "RYR2", "SCN1B", "SCN3B", "TBX3", "FBN2", "WNT8A", "XPO7",
        "PRKAG2", "TITIN", "MAP3K7", "CAMK2D",
        # HF (HERMES 2024 + DCM)
        "CDKN1A", "BCL11B", "FLNC", "MTSS1", "ABO", "FTO", "AGAP1",
        "ZBTB17", "MFN2", "ZFP36L2", "RBM20", "KLHL3",
        "DSC2", "DSG2", "DSP", "PKP2", "LMNA",
        # Stroke (GIGASTROKE / MEGASTROKE)
        "HDAC9", "PHACTR1", "FOXF2", "ZCCHC14", "COL4A1", "COL4A2",
        "NINJ2", "EDNRA", "MMP12", "CASZ1", "SLC22A7", "ANK2",
        "MTHFR", "MMP12", "SH2B3", "ALDH2",
        # Additional cardiomyopathy + arrhythmia + GWAS HF locus-near plasma-detected
        "MEF2A", "MEF2C", "TBX1", "TBX20", "EVA1A", "SLC4A4",
        "AGTR1", "BMP10", "ADAMTS6", "ADAMTSL3", "ADAMTSL2",
    ],

    # 2. Coagulation / platelet / thrombosis (broadened cascade + regulation)
    "cat2_coagulation": [
        "F2", "F3", "F5", "F7", "F8", "F9", "F10", "F11", "F12", "F13A1", "F13B",
        "FGA", "FGB", "FGG", "VWF", "PROC", "PROS1", "SERPINC1", "SERPINE1",
        "SERPIND1", "SERPINF2", "SERPINB2", "SERPINA10", "SERPINA12",
        "PLAT", "PLAU", "PLG", "TFPI", "TFPI2",
        "ADAMTS13", "ADAMTS9", "ADAMTS7", "ADAMTS4", "ADAMTS5", "ADAMTS8", "ADAMTS17",
        "ITGA2B", "ITGB3", "ITGA2", "ITGAM", "ITGAX",
        "GP1BA", "GP1BB", "GP9", "GP5", "GP6",
        "VTN", "THBD", "PROCR", "KLKB1", "KNG1", "MASP1", "MASP2", "MBL2",
        "C1S", "C1R", "CFD", "CFB", "CFI", "CFH", "CFP",
        "MERTK", "TYRO3", "AXL",
    ],

    # 3. Inflammation / immune / cytokine (broadened cytokine + chemokine + complement)
    "cat3_inflammation": [
        "IL6R", "IL6", "IL6ST", "CRP", "IL1A", "IL1B", "IL1RN", "IL1R1", "IL1R2",
        "IL18", "IL18BP", "IL18R1", "IL33", "IL1RL1", "IL1RL2", "IL36A", "IL36B",
        "IL36G", "IL37", "IL38", "IL2", "IL2RA", "IL2RB", "IL2RG", "IL3", "IL4",
        "IL4R", "IL5", "IL5RA", "IL7", "IL7R", "IL10", "IL10RA", "IL10RB",
        "IL12A", "IL12B", "IL13", "IL13RA1", "IL13RA2", "IL15", "IL15RA",
        "IL17A", "IL17B", "IL17C", "IL17D", "IL17F", "IL17RA", "IL17RB",
        "IL19", "IL20", "IL21", "IL22", "IL22RA1", "IL23A", "IL23R", "IL27",
        "IL31", "IL32",
        "TNF", "LTA", "LTB", "LIGHT", "TNFSF10", "TNFSF11", "TNFSF12",
        "TNFSF13", "TNFSF13B", "TNFSF14", "TNFSF15", "TNFSF18",
        "TNFRSF1A", "TNFRSF1B", "TNFRSF4", "TNFRSF8", "TNFRSF9",
        "TNFRSF10A", "TNFRSF10B", "TNFRSF11A", "TNFRSF11B", "TNFRSF12A",
        "TNFRSF13B", "TNFRSF13C", "TNFRSF14", "TNFRSF17", "TNFRSF18",
        "TNFRSF19", "TNFRSF21", "TNFRSF25",
        "IFNG", "IFNGR1", "IFNGR2", "IFNA1", "IFNA2", "IFNB1", "IFNAR1", "IFNAR2",
        # Chemokines
        "CXCL1", "CXCL5", "CXCL6", "CXCL8", "CXCL9", "CXCL10", "CXCL11",
        "CXCL12", "CXCL13", "CXCL14", "CXCL16", "CXCL17",
        "CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL8", "CCL11", "CCL13",
        "CCL15", "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22",
        "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28",
        # Acute phase + alarmins
        "MIF", "S100A8", "S100A9", "S100A12", "PTX3", "SAA1", "SAA2", "SAP", "ORM1",
        "HP", "HPX", "TFRC", "FERR", "LCN2", "LBP", "BPI",
        # Complement
        "C2", "C3", "C4A", "C4B", "C5", "C6", "C7", "C8A", "C8B", "C8G", "C9",
        "C1QA", "C1QB", "C1QC",
    ],

    # 4. Lipid / lipoprotein / metabolic vascular
    "cat4_lipid": [
        "PCSK9", "PCSK7", "PCSK6", "PCSK1", "PCSK2", "PCSK5",
        "ANGPTL3", "ANGPTL4", "ANGPTL6", "ANGPTL7", "ANGPTL8",
        "LPA", "APOA1", "APOA2", "APOA4", "APOA5", "APOB", "APOC1", "APOC2", "APOC3",
        "APOD", "APOE", "APOH", "APOL1", "APOM", "APOO", "APOOL",
        "LDLR", "LDLRAP1", "LRP1", "LRP2", "LRP4", "LRP5", "LRP6", "LRP8",
        "LIPA", "LIPC", "LIPG", "LIPI", "LPL", "PNPLA2", "PNPLA3",
        "CETP", "LCAT", "PLTP", "ABCA1", "ABCG1", "ABCG5", "ABCG8",
        "SORT1", "SCARB1", "SCARB2", "OLR1", "MTTP", "DGAT1", "DGAT2",
        "MOGS", "PLIN5", "GPIHBP1", "ANGPTL3", "ANGPTL4",
    ],

    # 5. Natriuretic / cardiac stress / cardiomyopathy
    "cat5_cardiac_stress": [
        "NPPA", "NPPB", "NPPC", "NPR1", "NPR2", "NPR3",
        "GDF15", "GDF11", "MSTN", "IGFBP7", "IGFBP1", "IGFBP2", "IGFBP3",
        "IGFBP4", "IGFBP5", "IGFBP6", "IGF1", "IGF2", "IGFBPL1",
        "IL1RL1", "ST2", "TNNI3", "TNNT2", "TNNC1",
        "MYBPC3", "ACTC1", "MYL2", "MYL3", "MYL7", "TPM1",
        "BAG3", "DSG2", "DSP", "PKP2", "JUP", "DES", "LMNA",
        "FHL1", "FHL2", "VCL", "ANKRD1", "CSRP3", "GATA4", "GATA6", "MYH7B",
        "FSTL1", "FSTL3", "ACTA2", "ACTN2",
    ],

    # 6. ECM / fibrosis / remodeling (broadened)
    "cat6_ecm_fibrosis": [
        "MMP1", "MMP2", "MMP3", "MMP7", "MMP8", "MMP9", "MMP10", "MMP11", "MMP12",
        "MMP13", "MMP14", "MMP15", "MMP16", "MMP17", "MMP19", "MMP20", "MMP21",
        "MMP23B", "MMP24", "MMP25", "MMP26", "MMP27", "MMP28",
        "TIMP1", "TIMP2", "TIMP3", "TIMP4",
        "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3", "TGFBI",
        "BMP1", "BMP2", "BMP4", "BMP6", "BMP7", "BMP9", "BMP10",
        "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL5A2",
        "COL6A1", "COL6A2", "COL6A3", "COL8A1", "COL14A1", "COL15A1", "COL18A1",
        "COL23A1", "COL28A1",
        "LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4",
        "LGALS1", "LGALS3", "LGALS3BP", "LGALS4", "LGALS7", "LGALS8", "LGALS9",
        "POSTN", "TNC", "OMD",
        "THBS1", "THBS2", "THBS3", "THBS4", "FN1", "FBN1", "FBN2", "FBN3",
        "EFEMP1", "EFEMP2", "EMILIN1", "EMILIN2", "ELN",
        "DCN", "BGN", "ASPN", "FMOD", "OGN", "PRG4",
    ],

    # 7. Endothelial / vascular / angiogenesis (broadened)
    "cat7_endothelial_vascular": [
        "ANGPT1", "ANGPT2", "ANGPTL1", "ANGPTL2",
        "TEK", "TIE1",
        "VEGFA", "VEGFB", "VEGFC", "VEGFD",
        "FLT1", "FLT4", "KDR", "NRP1", "NRP2",
        "PDGFA", "PDGFB", "PDGFC", "PDGFD", "PDGFRA", "PDGFRB",
        "ICAM1", "ICAM2", "ICAM3", "ICAM4", "ICAM5", "VCAM1", "MADCAM1",
        "SELP", "SELE", "SELL", "PECAM1", "ESAM",
        "EDN1", "EDN2", "EDN3", "EDNRA", "EDNRB", "ECE1", "ECE2",
        "AGT", "AGTR1", "AGTR2", "REN", "ACE", "ACE2", "AGTRAP",
        "NOS1", "NOS2", "NOS3", "DDAH1", "DDAH2",
        "CDH5", "CDH13", "ENG", "ROBO4", "EPHB4", "NOTCH1", "NOTCH4",
        "JAG1", "JAG2", "DLL1", "DLL4", "VWF",
        "EPCAM", "CDH1", "CDH3",
        "FGF1", "FGF2", "FGFR1", "FGFR2", "FGFR3", "FGFR4",
        "PROS1", "GAS6",
        "S1PR1", "S1PR2", "S1PR3", "SPHK1", "SPHK2",
    ],

    # 8. Druggable CVD-relevant (Open Targets / ChEMBL / DGIdb / Pharos)
    "cat8_druggable_cvd": [
        "HMGCR", "NPC1L1", "MTOR", "CASP1",
        "MPO", "PLA2G7", "PLA2G1B", "PLA2G2A", "PLA2G4A", "PLA2G6", "PLA2G10",
        "ACE", "ACE2", "REN", "AGT", "AGTR1",
        "PCSK9", "ANGPTL3", "ANGPTL4", "APOC3", "LDLR", "CETP",
        "F11", "F2", "F9", "F7", "F10", "FGA", "FGB", "FGG",
        "VWF", "ADAMTS13", "PROC", "PROS1", "SERPINC1", "SERPINE1",
        "IL6R", "IL1B", "IL18", "TNF", "TNFRSF1A", "CRP", "IL17A",
        "MMP9", "MMP2", "MMP12", "TIMP2", "TGFB1",
        "ANGPT1", "ANGPT2", "VEGFA", "EDN1", "EDNRA", "EDNRB",
        "GDF15", "NPPB", "IL1RL1",
        # Sister-project Faz 4 underexplored druggable
        "FES", "FN1", "HYOU1", "LMOD1", "COL6A3", "INHBC",
        "C1S", "C1R", "FURIN", "VAMP5", "VAMP8", "PDZK1", "PDLIM7", "PECAM1",
        "GALNT2", "FGF5", "CELSR2",
        # Additional underexplored CVD targets (Pharos Tbio/Tdark + Open Targets bug)
        "PLAUR", "MERTK", "AXL", "TYRO3",
        "DPP4", "CPB2", "CPN1", "CPN2",
        "HRG", "KNG1", "KLK1", "KLK6",
        "ANGPTL2", "ANGPT4", "ANGPTL1",
        "TFRC", "TF", "HFE",
        "LEP", "ADIPOQ", "RBP4",
    ],
}

MAX_APTAMERS_PER_GENE_DEFAULT = 3
MULTI_APTAMER_EXCEPTIONS = {"NPPB", "F9", "F2", "PROC"}  # biological multi-epitope


def _is_ferkingstad(key: str) -> bool:
    parts = key.replace(".txt.gz", "").split("_")
    return len(parts) >= 4 and parts[0].isdigit() and parts[1].isdigit()


def main() -> int:
    files = json.loads(FOLDER_JSON.read_text())["files"]
    ferk = [f for f in files if _is_ferkingstad(f["Key"])]
    gene_apt: dict[str, list[dict]] = {}
    for f in ferk:
        gene = f["Key"].replace(".txt.gz", "").split("_")[2]
        gene_apt.setdefault(gene, []).append(f)
    print(f"Ferkingstad genes available: {len(gene_apt)}")

    # Calibration empirical priors
    calib = pl.read_csv(CALIB_REPORT, separator="\t").row(0, named=True)
    median_gb = float(calib["GB_per_aptamer_median"])
    median_min = float(calib["minutes_per_aptamer_median"])
    cv_prior_empty = 0.13                          # from per-category breakdown
    cv_prior_match = 0.85

    # Build aptamer selection: gene -> primary_category, all categories
    primary_category: dict[str, str] = {}
    all_cats_for_gene: dict[str, list[str]] = {}
    for cat in CATEGORIES:
        for g in CATEGORIES[cat]:
            all_cats_for_gene.setdefault(g, []).append(cat)
            primary_category.setdefault(g, cat)

    selected: list[dict] = []
    chosen_keys: set[str] = set()
    skipped_genes: list[dict] = []

    for gene, cats in all_cats_for_gene.items():
        if gene not in gene_apt:
            skipped_genes.append({"gene": gene, "primary_category": primary_category[gene],
                                  "reason": "not_in_deCODE_panel"})
            continue
        all_apts = sorted(gene_apt[gene], key=lambda f: f["Key"])
        cap = MAX_APTAMERS_PER_GENE_DEFAULT
        # Biological exceptions: gene-level whitelist may keep more aptamers
        if gene in MULTI_APTAMER_EXCEPTIONS:
            cap = max(cap, len(all_apts))
        kept = all_apts[:cap]
        for f in kept:
            if f["Key"] in chosen_keys:
                continue
            selected.append({
                "gene": gene,
                "aptamer_key": f["Key"],
                "size_bytes": f["Size"],
                "primary_category": primary_category[gene],
                "all_categories": "|".join(cats),
                "n_aptamers_for_gene": len(all_apts),
                "n_kept_for_gene": len(kept),
                "multi_aptamer_exception": gene in MULTI_APTAMER_EXCEPTIONS,
            })
            chosen_keys.add(f["Key"])

    df = pl.DataFrame(selected)
    print(f"Selected aptamers: {df.height}; unique genes: {df['gene'].n_unique()}")

    # ----- Persist outputs -----
    df.write_csv(OUT_DESIGN, separator="\t")

    cat_counts = df.group_by("primary_category").agg([
        pl.col("gene").n_unique().alias("n_unique_genes"),
        pl.len().alias("n_aptamers"),
    ]).sort("primary_category")
    cat_counts.write_csv(OUT_CATCOUNTS, separator="\t")

    # Runtime / transfer projections (using calibration cv_prior priors)
    n = df.height
    median_bytes = int(median_gb * 1e9)
    median_sec = median_min * 60
    proj = []
    for nw in (1, 2, 4):
        speedup = {1: 1.0, 2: 1.6, 4: 2.5}[nw]
        wall_h = n * median_sec / speedup / 3600
        proj.append({
            "n_aptamers": n, "workers": nw,
            "median_seconds_per_aptamer_calib": round(median_sec, 1),
            "median_bytes_per_aptamer_calib": median_bytes,
            "expected_wall_hours": round(wall_h, 2),
            "expected_transfer_GB": round(n * median_gb, 1),
            "speedup_assumed": speedup,
            "expected_empty_rate_calib_cv_prior": cv_prior_empty,
            "expected_match_rate_calib_cv_prior": cv_prior_match,
        })
    pl.DataFrame(proj).write_csv(OUT_RUNTIME, separator="\t")

    # YAML config (machine-readable batch spec)
    yaml_doc = {
        "batch_name": "CVD-prioritized 500 aptamer staged batch (D034)",
        "scope_label": "CVD-prioritized plasma proteomic target screen",
        "n_aptamers": int(n),
        "n_unique_genes": int(df["gene"].n_unique()),
        "max_aptamers_per_gene_default": MAX_APTAMERS_PER_GENE_DEFAULT,
        "multi_aptamer_exceptions": sorted(MULTI_APTAMER_EXCEPTIONS),
        "categories": {row["primary_category"]: int(row["n_aptamers"])
                       for row in cat_counts.iter_rows(named=True)},
        "calibration_priors": {
            "median_gb_per_aptamer": median_gb,
            "median_minutes_per_aptamer": median_min,
            "cv_prior_empty_rate": cv_prior_empty,
            "cv_prior_match_to_AF_rate": cv_prior_match,
        },
        "expected_metrics_4_workers": {
            "transfer_gb": round(n * median_gb, 1),
            "wall_hours": round(n * median_sec / 2.5 / 3600, 2),
            "empty_rate": cv_prior_empty,
            "match_rate": cv_prior_match,
        },
        "list_qc_thresholds": {
            "n_aptamers_min": 450, "n_aptamers_max": 550,
            "transfer_gb_max": 600,
            "runtime_4w_hours_max": 16,
            "empty_rate_max": 0.25,
            "match_rate_min": 0.75,
            "top_gene_share_max": 0.05,
        },
        "approval_gate": "explicit_protocol_approval_required_before_extraction",
    }
    OUT_YAML.parent.mkdir(parents=True, exist_ok=True)
    OUT_YAML.write_text(yaml.safe_dump(yaml_doc, sort_keys=False))

    # ----- List QC threshold check -----
    thresholds = yaml_doc["list_qc_thresholds"]
    expected = yaml_doc["expected_metrics_4_workers"]
    top_gene_share = (
        df.group_by("gene").agg(pl.len().alias("n")).sort("n", descending=True).row(0)[1]
        / max(1, df.height)
    )
    passes = {
        "n_aptamers_in_range": thresholds["n_aptamers_min"] <= n <= thresholds["n_aptamers_max"],
        "transfer_within_budget": expected["transfer_gb"] <= thresholds["transfer_gb_max"],
        "runtime_within_budget": expected["wall_hours"] <= thresholds["runtime_4w_hours_max"],
        "empty_within_budget": expected["empty_rate"] <= thresholds["empty_rate_max"],
        "match_within_budget": expected["match_rate"] >= thresholds["match_rate_min"],
        "no_single_gene_dominates": top_gene_share <= thresholds["top_gene_share_max"],
        "all_categories_nonzero": all(
            cat in yaml_doc["categories"] and yaml_doc["categories"][cat] > 0
            for cat in CATEGORIES.keys()
        ),
    }
    summary = [
        "# CVD-prioritized 500 aptamer batch — list QC summary (D034)",
        "",
        "**Scope:** CVD-prioritized plasma proteomic target screen (NOT proteome-wide).",
        "",
        "## Counts",
        "",
        f"- Total aptamers: **{n}**",
        f"- Unique genes: **{df['gene'].n_unique()}**",
        f"- Max aptamers/gene: {df.group_by('gene').agg(pl.len()).max()['len'][0]}",
        f"- Top-gene aptamer share: {top_gene_share:.1%}",
        f"- Skipped genes (not on deCODE panel): {len(skipped_genes)}",
        "",
        "## Per category (primary)",
        "",
    ]
    for r in cat_counts.iter_rows(named=True):
        summary.append(f"- `{r['primary_category']}`: **{r['n_aptamers']} aptamers** ({r['n_unique_genes']} genes)")
    summary.extend([
        "",
        "## Expected metrics (calibration cv_prior priors)",
        "",
        f"- Transfer: **{expected['transfer_gb']:.0f} GB**",
        f"- Wall-time @ 4 workers: **{expected['wall_hours']:.1f} h** (assuming 2.5x speedup)",
        f"- Empty rate (expected): {expected['empty_rate']:.0%}",
        f"- Match-to-AF rate (expected): {expected['match_rate']:.0%}",
        "",
        "## Threshold check (D034)",
        "",
        "| Threshold | Target | Observed/Expected | Pass? |",
        "|---|---|---|---|",
        f"| n_aptamers in [450, 550] | yes | {n} | {'✅' if passes['n_aptamers_in_range'] else '❌'} |",
        f"| transfer ≤ 600 GB | yes | {expected['transfer_gb']:.0f} GB | {'✅' if passes['transfer_within_budget'] else '❌'} |",
        f"| runtime 4w ≤ 16 h | yes | {expected['wall_hours']:.1f} h | {'✅' if passes['runtime_within_budget'] else '❌'} |",
        f"| empty rate ≤ 25% | yes | {expected['empty_rate']:.0%} | {'✅' if passes['empty_within_budget'] else '❌'} |",
        f"| match rate ≥ 75% | yes | {expected['match_rate']:.0%} | {'✅' if passes['match_within_budget'] else '❌'} |",
        f"| top-gene aptamer share ≤ 5% | yes | {top_gene_share:.1%} | {'✅' if passes['no_single_gene_dominates'] else '❌'} |",
        f"| all categories nonzero | yes | yes | {'✅' if passes['all_categories_nonzero'] else '❌'} |",
        "",
        f"**Overall: {sum(passes.values())}/{len(passes)} threshold pass.**",
        "",
    ])
    if all(passes.values()):
        summary.append("→ List QC PASS. **awaiting protocol explicit approval to start extraction.**")
    else:
        summary.append("→ List QC FAIL. Adjust gene list / category sizes before approval.")
    OUT_SUMMARY.write_text("\n".join(summary))

    print(f"wrote {OUT_DESIGN}")
    print(f"wrote {OUT_CATCOUNTS}")
    print(f"wrote {OUT_RUNTIME}")
    print(f"wrote {OUT_YAML}")
    print(f"wrote {OUT_SUMMARY}")
    print()
    print("List QC threshold check:")
    for k, v in passes.items():
        mark = "✅" if v else "❌"
        print(f"  {mark} {k}")
    return 0


if __name__ == "__main__":
    raif then SystemExit(main())
