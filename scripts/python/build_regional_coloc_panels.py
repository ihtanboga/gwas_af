#!/usr/bin/env python3
"""Regional colocalization panel figures (D048 supplemental request).

Produces 4 regional overlay panels (top: pQTL -log10 p; bottom: outcome -log10 p
across the same ±500 kb window) and saves to figures/ + manuscript/figures/:

  fig_regional_coloc_NPPA_AF.{pdf,png}            — strong shared (PP.H4=0.911)
  fig_regional_coloc_NPPA_CES.{pdf,png}           — moderate shared (PP.H4=0.565)
  fig_regional_coloc_IL6R_AF_H3_discordance.{pdf,png}  — H3=1.000, distinct leads
  fig_regional_coloc_MMP12_LAS.{pdf,png}          — strong non-AF shared (PP.H4=0.916)

Reuses existing aligned audit TSVs for IL6R×AF and NPPA×CES; regenerates aligned
data for NPPA×AF and MMP12×LAS via the Wave 1/Wave 2 download + extract path.
"""

from __future__ import annotations

import math
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import polars as pl  # noqa: E402

ROOT = Path("/Users/apple/Desktop/gwas_af")
sys.path.insert(0, str(ROOT / "scripts/python"))

from run_decode_cvd150_wave1_coloc import (  # noqa: E402
    OUTCOMES, REGION_KB_PRIMARY, TMP_DIR, _log,
    download_aptamer, extract_region, harmonize_region, init_environment,
)

FIG_DIR = ROOT / "figures"
MS_FIG_DIR = ROOT / "manuscript/figures"
AUDIT_DIR = ROOT / "results/coloc/preflight_audit"
AUDIT_DIR.mkdir(parents=True, exist_ok=True)

PANELS = [
    {
        "tag": "NPPA_AF",
        "gene": "NPPA",
        "aptamer": "5443_62_NPPA_ANP.txt.gz",
        "outcome": "AF",
        "title": "NPPA pQTL × AF — strong coloc (PP.H4 = 0.911)",
        "subtitle": "Wave 1 primary; AF-anchored axis. Lead pQTL rs145488887.",
        "outcome_label": "AF",
    },
    {
        "tag": "NPPA_CES",
        "gene": "NPPA",
        "aptamer": "5443_62_NPPA_ANP.txt.gz",
        "outcome": "cardioembolic_stroke",
        "title": "NPPA pQTL × cardioembolic stroke — moderate coloc (PP.H4 = 0.565)",
        "subtitle": "Wave 1 primary; AF-anchored downstream axis (sensitivity PP.H4 = 0.395).",
        "outcome_label": "CES",
    },
    {
        "tag": "IL6R_AF_H3_discordance",
        "gene": "IL6R",
        "aptamer": "15602_43_IL6R_IL_6_sRa.txt.gz",
        "outcome": "AF",
        "title": "IL6R pQTL × AF — LD-confounded discordance (PP.H3 = 1.000, PP.H4 = 0.000)",
        "subtitle": "Wave 1 primary; cis-pQTL lead variant ≠ AF lead variant. AF MR q = 4.1e-13 but coloc rejects shared causal variant.",
        "outcome_label": "AF",
    },
    {
        "tag": "MMP12_LAS",
        "gene": "MMP12",
        "aptamer": "4496_60_MMP12_MMP_12.txt.gz",
        "outcome": "large_artery_stroke",
        "title": "MMP12 pQTL × large-artery stroke — strong coloc (PP.H4 = 0.916)",
        "subtitle": "Wave 2 primary; non-AF atherothrombotic axis. Lead pQTL rs470530.",
        "outcome_label": "LAS",
    },
]


# ---------------------------------------------------------------------------
# Step 1: ensure each panel has an aligned TSV
# ---------------------------------------------------------------------------


def aligned_tsv_path(panel: dict) -> Path:
    return AUDIT_DIR / f"aligned_{panel['gene']}_{panel['outcome']}.tsv"


def ensure_aligned_tsvs():
    """Regenerate any missing aligned TSV by re-downloading the aptamer + extracting region."""
    needed = []
    for panel in PANELS:
        if not aligned_tsv_path(panel).exists():
            needed.append(panel)
    if not needed:
        _log("all 4 aligned TSVs already exist; skipping deCODE fetch")
        return

    _log(f"need to regenerate {len(needed)} aligned TSV(s)")
    target_apts = {p["aptamer"] for p in needed}
    leads, apt_paths, _coord_cache, _instr, caches = init_environment(target_apts)
    pqtl_region_cache = caches["region"]

    for panel in needed:
        apt = panel["aptamer"]
        if apt not in apt_paths or apt not in leads:
            _log(f"SKIP {panel['tag']}: download/lead missing for {apt}")
            continue
        if apt not in pqtl_region_cache:
            lead = leads[apt]
            start = max(0, lead["pos"] - REGION_KB_PRIMARY * 1000)
            end = lead["pos"] + REGION_KB_PRIMARY * 1000
            pqtl_region_cache[apt] = extract_region(apt_paths[apt], lead["chrom"], start, end)
            _log(f"  region {apt}: {len(pqtl_region_cache[apt])} variants")
        outcome_meta = OUTCOMES[panel["outcome"]]
        outcome_df = pl.read_parquet(outcome_meta["path"]).select([
            "chr", "pos", "rsid", "effect_allele", "other_allele",
            "beta", "se", "pval", "eaf",
        ])
        # Build-aware: if outcome is GRCh37, liftover region first
        if outcome_meta["build"] == "GRCh37":
            from pyliftover import LiftOver
            lo = LiftOver(str(ROOT / "data/raw/annotations/hg38ToHg19.over.chain.gz"))
            pqtl_for_outcome = []
            for r in pqtl_region_cache[apt]:
                cv = lo.convert_coordinate("chr" + str(r["chrom"]), int(r["pos"]))
                if cv:
                    r2 = dict(r)
                    r2["pos_b37"] = int(cv[0][1])
                    pqtl_for_outcome.append(r2)
        else:
            pqtl_for_outcome = pqtl_region_cache[apt]
        aligned, counters = harmonize_region(pqtl_for_outcome, outcome_df, build=outcome_meta["build"])
        if not aligned:
            _log(f"WARN {panel['tag']}: no aligned variants ({counters})")
            continue
        out = aligned_tsv_path(panel)
        pl.DataFrame(aligned).write_csv(out, separator="\t")
        _log(f"  wrote {out} ({len(aligned)} variants)")

    # Cleanup downloaded .gz to keep disk lean
    for f in TMP_DIR.iterdir():
        try:
            f.unlink()
        except OSError:
            pass


# ---------------------------------------------------------------------------
# Step 2: render the 4 regional overlay panels
# ---------------------------------------------------------------------------


def _safe_neg_log10(p):
    try:
        v = float(p)
        if v <= 0 or v != v:
            return float("nan")
        return -math.log10(v)
    except (TypeError, ValueError):
        return float("nan")


def render_panel(panel: dict):
    tsv = aligned_tsv_path(panel)
    if not tsv.exists():
        _log(f"SKIP render {panel['tag']}: aligned TSV missing")
        return
    df = pl.read_csv(tsv, separator="\t", infer_schema_length=0)
    pos = df["pos"].cast(pl.Int64).to_list()
    pval_pqtl = [float(p) for p in df["pval"].to_list()]
    pval_outcome = [float(p) for p in df["pval_outcome"].to_list()]
    nlp_pqtl = [_safe_neg_log10(p) for p in pval_pqtl]
    nlp_outcome = [_safe_neg_log10(p) for p in pval_outcome]

    # Lead variants
    lead_pqtl_idx = min(range(len(pval_pqtl)), key=lambda i: pval_pqtl[i])
    lead_outcome_idx = min(range(len(pval_outcome)), key=lambda i: pval_outcome[i])
    lead_pqtl_pos = pos[lead_pqtl_idx]
    lead_outcome_pos = pos[lead_outcome_idx]
    lead_pqtl_rsid = df["rsid"][lead_pqtl_idx] or df["name"][lead_pqtl_idx]
    lead_outcome_rsid = df["rsid"][lead_outcome_idx] or df["name"][lead_outcome_idx]

    # Coloring: nearest-lead (pQTL lead) by distance bucket
    max_dist = max(abs(p - lead_pqtl_pos) for p in pos) or 1
    distances = [abs(p - lead_pqtl_pos) / max_dist for p in pos]

    chrom = df["chrom"][0]

    fig, axes = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True,
                              gridspec_kw={"height_ratios": [1, 1]})
    ax1, ax2 = axes

    # pQTL track
    ax1.scatter(pos, nlp_pqtl, c=distances, cmap="viridis_r",
                s=14, edgecolor="black", lw=0.2, alpha=0.85)
    ax1.axvline(lead_pqtl_pos, color="#cb6a52", lw=1.0, ls="--", alpha=0.9)
    ax1.scatter([lead_pqtl_pos], [nlp_pqtl[lead_pqtl_idx]],
                s=80, marker="D", facecolor="#cb6a52", edgecolor="black", lw=1.0,
                label=f"pQTL lead {lead_pqtl_rsid}", zorder=5)
    ax1.set_ylabel("pQTL  −log₁₀(p)")
    ax1.legend(loc="upper right", fontsize=8)
    ax1.grid(alpha=0.2)
    ax1.axhline(-math.log10(5e-8), color="gray", lw=0.5, ls=":", alpha=0.5)

    # Outcome track
    ax2.scatter(pos, nlp_outcome, c=distances, cmap="viridis_r",
                s=14, edgecolor="black", lw=0.2, alpha=0.85)
    ax2.axvline(lead_pqtl_pos, color="#cb6a52", lw=1.0, ls="--", alpha=0.6,
                label="pQTL lead")
    ax2.axvline(lead_outcome_pos, color="#5d9bdb", lw=1.0, ls="--", alpha=0.9)
    ax2.scatter([lead_outcome_pos], [nlp_outcome[lead_outcome_idx]],
                s=80, marker="o", facecolor="#5d9bdb", edgecolor="black", lw=1.0,
                label=f"{panel['outcome_label']} lead {lead_outcome_rsid}", zorder=5)
    ax2.set_ylabel(f"{panel['outcome_label']}  −log₁₀(p)")
    ax2.set_xlabel(f"chr{chrom} position (bp, GRCh38)")
    ax2.axhline(-math.log10(5e-8), color="gray", lw=0.5, ls=":", alpha=0.5)
    ax2.legend(loc="upper right", fontsize=8)
    ax2.grid(alpha=0.2)

    fig.suptitle(panel["title"], fontsize=11, fontweight="bold", y=0.995)
    fig.text(0.5, 0.945, panel["subtitle"], ha="center", fontsize=9, style="italic")

    # Distance from lead colorbar
    sm = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(0, 1), cmap="viridis_r")
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation="vertical", fraction=0.025, pad=0.02)
    cbar.set_label("distance from pQTL lead (normalized)", fontsize=8)

    # Footer with claim discipline
    fig.text(0.5, 0.005,
             "Aligned variants in ±500 kb window after allele harmonization. "
             "All labels INTERIM (D041 / D043 / D044 / D047 / D048).",
             ha="center", fontsize=7, style="italic", color="gray")

    out_pdf = FIG_DIR / f"fig_regional_coloc_{panel['tag']}.pdf"
    out_png = FIG_DIR / f"fig_regional_coloc_{panel['tag']}.png"
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.close(fig)
    _log(f"wrote {out_pdf.name} + {out_png.name}")

    # Also copy into manuscript/figures/ so the staged packages get them
    shutil.copy2(out_pdf, MS_FIG_DIR / out_pdf.name)
    shutil.copy2(out_png, MS_FIG_DIR / out_png.name)


def main() -> int:
    ensure_aligned_tsvs()
    for panel in PANELS:
        render_panel(panel)
    print("\n=== Regional coloc panels (4) generated ===")
    for panel in PANELS:
        for ext in ("pdf", "png"):
            f = FIG_DIR / f"fig_regional_coloc_{panel['tag']}.{ext}"
            print(f"  {f}: {'OK' if f.exists() else 'MISSING'} ({f.stat().st_size if f.exists() else 0:,} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
