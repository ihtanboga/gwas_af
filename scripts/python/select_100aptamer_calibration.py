#!/usr/bin/env python3
"""D033 Step 3 — pick 100 aptamer files for the calibration batch.

Categories (per protocol):
  30 cardiovascular-prior  (curated list, present in deCODE folder)
  20 multi-aptamer ambiguous (genes with >=2 aptamers)
  20 large-fwith / high-risk  (largest Ferkingstad files not already picked)
  30 random                  (deterministic seed; remaining genes)

Output:
  data/registry/decode_100aptamer_calibration_targets.tsv
"""

from __future__ import annotations

import json
import random
from pathlib import Path

import polars as pl

ROOT = Path("/Users/apple/Desktop/gwas_af")
FOLDER = ROOT / "data/registry/decode_2021_folder.json"
OUT = ROOT / "data/registry/decode_100aptamer_calibration_targets.tsv"

# 30 CV-prior candidates (will be intersected with the panel; first 30 hits).
CV_PRIOR = [
    # natriuretic / cardiac stress
    "NPPA", "NPPB", "GDF15", "IGFBP7",
    # lipid metabolism
    "PCSK9", "APOB", "APOE", "APOA1", "APOA4", "ANGPTL3", "ANGPTL4", "LDLR", "CETP",
    # coagulation / platelet
    "F2", "F11", "F9", "F7", "FGB", "VWF", "PROC", "PROS1", "SERPINC1", "ADAMTS13",
    "PLAT", "PLG", "SERPINE1", "ITGA2B",
    # inflammation
    "IL6R", "IL1B", "IL18", "CRP", "TNFRSF1A", "TNFRSF1B", "VCAM1", "ICAM1", "SELE", "SELP",
    # ECM / fibrosis / remodeling
    "MMP1", "MMP2", "MMP3", "MMP9", "TIMP1", "TIMP2", "COL1A1", "TGFB1", "TGFB2",
    # endothelial / vascular
    "ANGPT1", "ANGPT2", "VEGFA", "VEGFC", "PDGFB", "EDN1",
    # cardiomyopathy
    "BAG3", "MYBPC3", "ACTC1",
    # additional anchors / druggable
    "LPA",  # may be absent
    "REN", "AGT",  # RAAS
]

# Curated random seed for reproducibility
RANDOM_SEED = 1742_05_09


def _is_ferkingstad(key: str) -> bool:
    """Return True for SomaScan_v4 Ferkingstad-pattern keys."""
    parts = key.replace(".txt.gz", "").split("_")
    return len(parts) >= 4 and parts[0].isdigit() and parts[1].isdigit()


def main() -> int:
    files = json.loads(FOLDER.read_text())["files"]
    ferk = [f for f in files if _is_ferkingstad(f["Key"])]
    print(f"Ferkingstad files: {len(ferk)}")

    # Build gene -> [aptamer files]
    gene_apt: dict[str, list[dict]] = {}
    for f in ferk:
        gene = f["Key"].replace(".txt.gz", "").split("_")[2]
        gene_apt.setdefault(gene, []).append(f)
    print(f"Unique genes (Ferkingstad subset): {len(gene_apt)}")

    selected: list[dict] = []
    used_keys: set[str] = set()

    # 1) CV-prior — 30 aptamers (one per gene; first available aptamer)
    cv_added = 0
    cv_skipped = []
    for gene in CV_PRIOR:
        if cv_added >= 30:
            break
        if gene not in gene_apt:
            cv_skipped.append(gene)
            continue
        # pick first aptamer (sorted by SeqId for determinism)
        apts = sorted(gene_apt[gene], key=lambda f: f["Key"])
        f = apts[0]
        if f["Key"] in used_keys:
            continue
        selected.append({
            "category": "cv_prior", "gene": gene,
            "aptamer_key": f["Key"], "size_bytes": f["Size"],
        })
        used_keys.add(f["Key"])
        cv_added += 1
    print(f"CV-prior selected: {cv_added} (skipped not on panel: {cv_skipped})")

    # 2) Multi-aptamer ambiguous — 20 aptamers (pick from genes with >=2 aptamers,
    # take first 2 aptamers per gene until 20 reached)
    multi_genes = sorted(g for g, a in gene_apt.items() if len(a) >= 2)
    rng = random.Random(RANDOM_SEED)
    rng.shuffle(multi_genes)
    multi_added = 0
    for gene in multi_genes:
        if multi_added >= 20:
            break
        for f in sorted(gene_apt[gene], key=lambda f: f["Key"]):
            if f["Key"] in used_keys:
                continue
            selected.append({
                "category": "multi_aptamer", "gene": gene,
                "aptamer_key": f["Key"], "size_bytes": f["Size"],
            })
            used_keys.add(f["Key"])
            multi_added += 1
            if multi_added >= 20:
                break
    print(f"Multi-aptamer selected: {multi_added}")

    # 3) Large-file/high-risk — 20 largest remaining
    remaining = [f for f in ferk if f["Key"] not in used_keys]
    remaining_sorted = sorted(remaining, key=lambda f: -f["Size"])
    for f in remaining_sorted[:20]:
        gene = f["Key"].replace(".txt.gz", "").split("_")[2]
        selected.append({
            "category": "large_file", "gene": gene,
            "aptamer_key": f["Key"], "size_bytes": f["Size"],
        })
        used_keys.add(f["Key"])
    print(f"Large-fwith selected: 20")

    # 4) Random — 30 aptamers from remaining (deterministic shuffle)
    remaining = [f for f in ferk if f["Key"] not in used_keys]
    rng2 = random.Random(RANDOM_SEED + 1)
    rng2.shuffle(remaining)
    for f in remaining[:30]:
        gene = f["Key"].replace(".txt.gz", "").split("_")[2]
        selected.append({
            "category": "random", "gene": gene,
            "aptamer_key": f["Key"], "size_bytes": f["Size"],
        })
        used_keys.add(f["Key"])
    print(f"Random selected: 30")

    print(f"Total: {len(selected)}")
    df = pl.DataFrame(selected)
    df.write_csv(OUT, separator="\t")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raif then SystemExit(main())
