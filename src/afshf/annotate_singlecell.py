"""Single-cell / spatial expression annotation (spec section 11).

Skeleton: actual data ingestion (Heart Cell Atlas v2 .h5ad, plaque atlases)
arrives in Phase 5. This module exposes the scoring math so it can be
unit-tested independently of the data files.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable


# Cell type families relevant to the project's safety reasoning.
ATRIAL_CARDIOMYOCYTE = "atrial_cardiomyocyte"
VENTRICULAR_CARDIOMYOCYTE = "ventricular_cardiomyocyte"
CONDUCTION_SYSTEM = "conduction_system"
FIBROBLAST = "fibroblast"
ENDOTHELIAL = "endothelial"
SMOOTH_MUSCLE = "vascular_smooth_muscle"
MACROPHAGE = "macrophage"


@dataclass
class CellExpressionProfile:
    target_id: str
    gene_symbol: str
    expression_specificity: dict[str, float] = field(default_factory=dict)
    disease_state_enrichment: dict[str, float] = field(default_factory=dict)
    tissue_relevance_bonus: float = 0.0


def cell_relevance_score(profile: CellExpressionProfile,
                         cells_of_interest: Iterable[str] | None = None) -> float:
    """Spec §11 scoring: specificity + disease enrichment + tissue bonus."""
    cells = list(cells_of_interest) if cells_of_interest is not None else None
    if cells is None:
        spec = sum(profile.expression_specificity.values())
        enrich = sum(profile.disease_state_enrichment.values())
    else:
        spec = sum(profile.expression_specificity.get(c, 0.0) for c in cells)
        enrich = sum(profile.disease_state_enrichment.get(c, 0.0) for c in cells)
    return spec + enrich + profile.tissue_relevance_bonus


def safety_flags_from_expression(profile: CellExpressionProfile,
                                 specificity_threshold: float = 0.3) -> list[str]:
    """Return interpretation flags driven by the cell-type expression context.

    Implements the rule matrix from spec §11.
    """
    flags: list[str] = []
    spec = profile.expression_specificity

    if spec.get(CONDUCTION_SYSTEM, 0.0) >= specificity_threshold:
        flags.append("conduction_safety_panel_required")
    if spec.get(VENTRICULAR_CARDIOMYOCYTE, 0.0) >= specificity_threshold:
        flags.append("hf_cardiomyopathy_safety_panel_required")
    if spec.get(ENDOTHELIAL, 0.0) >= specificity_threshold or spec.get(
        SMOOTH_MUSCLE, 0.0
    ) >= specificity_threshold:
        flags.append("stroke_subtype_safety_panel_required")
    if spec.get(FIBROBLAST, 0.0) >= specificity_threshold:
        flags.append("hfpef_fibrosis_interpretation")
    if spec.get(MACROPHAGE, 0.0) >= specificity_threshold:
        flags.append("inflammation_safety_annotation")
    return flags
