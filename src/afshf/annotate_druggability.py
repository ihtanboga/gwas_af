"""Druggability + disease-prior + clinical-phase annotation (spec §12, D024).

Only thin client logic lives here. Network calls are gated behind feature
flags so unit tests can run offline.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional

import pandas as pd


# Disease-prior knowledge tags (D024)
KNOWN_AF_TARGET = "known_AF_target"
KNOWN_HF_TARGET = "known_HF_target"
KNOWN_STROKE_TARGET = "known_stroke_target"
KNOWN_CVD_TARGET = "known_CVD_target"
EMERGING_TARGET = "emerging_target"
DARK_BUT_DRUGGABLE = "dark_but_druggable"

VALID_DISEASE_PRIOR_TAGS = frozenset({
    KNOWN_AF_TARGET, KNOWN_HF_TARGET, KNOWN_STROKE_TARGET, KNOWN_CVD_TARGET,
    EMERGING_TARGET, DARK_BUT_DRUGGABLE,
})


@dataclass
class DruggabilityRecord:
    target_id: str
    gene_symbol: str
    open_targets_score: Optional[float] = None
    chembl_id: Optional[str] = None
    max_clinical_phase: Optional[int] = None
    known_drugs_count: int = 0
    drug_modality: Optional[str] = None
    pharos_tdl: Optional[str] = None
    novelty_class: Optional[str] = None  # established / emerging / dark / etc.
    repurposing_ready: bool = False
    # D024 — disease-prior knowledge tags
    disease_prior_tags: list[str] = field(default_factory=list)
    known_drug_target: bool = False


def classify_novelty(record: DruggabilityRecord) -> str:
    """Classify a target into the spec section-12 novelty buckets."""
    if record.repurposing_ready or (record.max_clinical_phase or 0) >= 3:
        return "Repurposing-ready"
    if record.pharos_tdl == "Tclin":
        return "Established"
    if record.pharos_tdl == "Tchem":
        return "Emerging"
    if record.pharos_tdl == "Tbio":
        return "Emerging"
    if record.pharos_tdl == "Tdark":
        return "Dark but druggable"
    return "Undruggable current"


# ---------------------------------------------------------------------------
# Disease-prior tag derivation from positive controls + Open Targets scores
# ---------------------------------------------------------------------------


_POSITIVE_CONTROL_TRAIT_TO_TAG = {
    "AF": KNOWN_AF_TARGET,
    "HF_overall": KNOWN_HF_TARGET,
    "HF_nonischemic": KNOWN_HF_TARGET,
    "HF_ni_HFrEF": KNOWN_HF_TARGET,
    "HF_ni_HFpEF": KNOWN_HF_TARGET,
    "ischemic_stroke": KNOWN_STROKE_TARGET,
    "cardioembolic_stroke": KNOWN_STROKE_TARGET,
    "large_artery_stroke": KNOWN_STROKE_TARGET,
    "small_vessel_stroke": KNOWN_STROKE_TARGET,
}


def disease_prior_tags_from_controls(gene_symbol: str,
                                     positive_controls: pd.DataFrame
                                     ) -> list[str]:
    """Derive disease-prior tags from positive_controls.tsv membership."""
    if positive_controls is None or positive_controls.empty:
        return []
    rows = positive_controls.loc[positive_controls["gene_symbol"] == gene_symbol]
    tags: list[str] = []
    for trait in rows["trait"].unique():
        tag = _POSITIVE_CONTROL_TRAIT_TO_TAG.get(str(trait))
        if tag and tag not in tags:
            tags.append(tag)
    if tags and KNOWN_CVD_TARGET not in tags:
        tags.append(KNOWN_CVD_TARGET)
    return tags


def annotate_disease_priors(record: DruggabilityRecord,
                            positive_controls: pd.DataFrame | None = None,
                            opentargets_disease_score: float | None = None,
                            opentargets_known_threshold: float = 0.5) -> DruggabilityRecord:
    """Populate ``disease_prior_tags`` from positive controls + OT score.

    Mutates a *copy* of the record and returns the new record.
    """
    tags = list(record.disease_prior_tags)
    if positive_controls is not None:
        for tag in disease_prior_tags_from_controls(record.gene_symbol, positive_controls):
            if tag not in tags:
                tags.append(tag)
    if (opentargets_disease_score is not None
            and opentargets_disease_score >= opentargets_known_threshold
            and KNOWN_CVD_TARGET not in tags):
        tags.append(KNOWN_CVD_TARGET)
    if not tags:
        # No prior cardiac evidence -> tag according to druggability frontier.
        if record.pharos_tdl == "Tdark":
            tags.append(DARK_BUT_DRUGGABLE)
        elif record.pharos_tdl in {"Tbio", "Tchem"}:
            tags.append(EMERGING_TARGET)

    return DruggabilityRecord(
        target_id=record.target_id,
        gene_symbol=record.gene_symbol,
        open_targets_score=record.open_targets_score,
        chembl_id=record.chembl_id,
        max_clinical_phase=record.max_clinical_phase,
        known_drugs_count=record.known_drugs_count,
        drug_modality=record.drug_modality,
        pharos_tdl=record.pharos_tdl,
        novelty_class=record.novelty_class,
        repurposing_ready=record.repurposing_ready,
        disease_prior_tags=tags,
        known_drug_target=bool(record.known_drugs_count > 0
                               or (record.max_clinical_phase or 0) >= 1),
    )


def to_table(records: Iterable[DruggabilityRecord]) -> pd.DataFrame:
    rows = []
    for r in records:
        rows.append({
            "target_id": r.target_id,
            "gene_symbol": r.gene_symbol,
            "open_targets_score": r.open_targets_score,
            "chembl_id": r.chembl_id,
            "max_clinical_phase": r.max_clinical_phase,
            "known_drugs_count": r.known_drugs_count,
            "drug_modality": r.drug_modality,
            "pharos_tdl": r.pharos_tdl,
            "novelty_class": r.novelty_class or classify_novelty(r),
            "repurposing_ready": r.repurposing_ready,
            "disease_prior_tags": "|".join(r.disease_prior_tags),
            "known_drug_target": r.known_drug_target,
        })
    return pd.DataFrame(rows)


# Real Open Targets / ChEMBL / DGIdb / Pharos clients live in a separate
# module added during Phase 4 once API tokens are configured.
