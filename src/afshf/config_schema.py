"""Pydantic schemas for the AFSHF-TARGETMAP configuration files.

These schemas validate `config/*.yml` at pipeline start. They are intentionally
permissive about optional fields so the registry can be populated incrementally
during Phase 0 (dataset discovery / download).
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal, Optional

import yaml
from pydantic import BaseModel, ConfigDict, Field


# ---------------------------------------------------------------------------
# project.yml
# ---------------------------------------------------------------------------


class ProjectMeta(BaseModel):
    name: str
    short_name: str
    description: str
    genome_build_preferred: Literal["GRCh37", "GRCh38"]
    primary_ancestry: Literal["EUR", "MIXED", "EAS", "AFR", "AMR"]
    allow_mixed_ancestry_secondary: bool = True


class Thresholds(BaseModel):
    qtl_p_primary: float
    qtl_p_exploratory: float
    min_f_stat: float
    clump_r2: float
    clump_kb: int
    cis_window_kb: int
    coloc_strong_pph4: float
    coloc_moderate_pph4: float
    coloc_lcl_h3_warn: float
    mr_fdr_primary: float
    mr_p_replication: float
    maf_min_primary: float


class InterventionConfig(BaseModel):
    define_action_from_af: bool
    action_rule: str
    report_or_scale: str
    af_signal_min_p: float
    af_signal_min_pph4: float


class OutcomePanels(BaseModel):
    primary: list[str]
    safety: list[str]
    exploratory: list[str] = Field(default_factory=list)


class ScoringConfig(BaseModel):
    harm_penalty_lambda: float
    max_logp_component: float
    severity_weights: dict[str, float]
    coloc_weight_pph4_strong: float = 1.0
    coloc_weight_pph4_moderate: float = 0.6
    coloc_weight_pph4_weak: float = 0.2
    replication_weight_replicated: float = 1.0
    replication_weight_partial: float = 0.6
    replication_weight_none: float = 0.3


class PathsConfig(BaseModel):
    data_raw: str
    data_processed: str
    results: str
    figures: str
    reports: str
    logs: str


class ProjectConfig(BaseModel):
    """Top-level config matching ``config/project.yml``."""

    project: ProjectMeta
    thresholds: Thresholds
    intervention: InterventionConfig
    outcomes: OutcomePanels
    scoring: ScoringConfig
    paths: PathsConfig

    model_config = ConfigDict(extra="forbid")

    @classmethod
    def load(cls, path: str | Path) -> "ProjectConfig":
        with open(path, "r", encoding="utf-8") as f:
            return cls.model_validate(yaml.safe_load(f))


# ---------------------------------------------------------------------------
# datasets.yml
# ---------------------------------------------------------------------------


SampleOverlapFlag = Literal[
    "verified_UKB_free",
    "verified_includes_UKB",
    "possible_UKB_overlap",
    "possible_internal_ISGC_overlap",        # D028 — ISGC consortium internal cohort reuse
    "unknown",
]


class OutcomeDataset(BaseModel):
    trait: str
    trait_group: Literal["AF", "HF", "STROKE", "SAFETY"]
    source: str
    dataset_label: str
    ancestry: Literal["EUR", "MIXED", "EAS", "AFR", "AMR"]
    genome_build: Literal["GRCh37", "GRCh38"]
    n_cases: Optional[int] = None
    n_controls: Optional[int] = None
    sample_overlap_notes: str = ""
    sample_overlap_flag: SampleOverlapFlag = "unknown"  # D022
    download_url_env: Optional[str] = None
    local_path: str
    license: Optional[str] = None
    primary_or_replication: Literal["primary", "replication", "sensitivity"]


class DatasetRegistry(BaseModel):
    outcomes: dict[str, OutcomeDataset]

    @classmethod
    def load(cls, path: str | Path) -> "DatasetRegistry":
        with open(path, "r", encoding="utf-8") as f:
            return cls.model_validate(yaml.safe_load(f))


# ---------------------------------------------------------------------------
# qtl_sources.yml
# ---------------------------------------------------------------------------


class QtlSource(BaseModel):
    type: Literal["pQTL", "eQTL"]
    assay: str
    role: Literal["discovery", "replication", "secondary", "validation"]    # D031
    ancestry: Literal["EUR", "MIXED", "EAS", "AFR", "AMR"]
    genome_build: Literal["GRCh37", "GRCh38"]
    cis_window_kb: int
    p_primary: float
    p_exploratory: float
    download_url_env: Optional[str] = None
    local_path: str
    license: Optional[str] = None
    sample_overlap_warning: str = ""

    model_config = ConfigDict(extra="allow")


class QtlRegistry(BaseModel):
    qtl_sources: dict[str, QtlSource]

    @classmethod
    def load(cls, path: str | Path) -> "QtlRegistry":
        with open(path, "r", encoding="utf-8") as f:
            return cls.model_validate(yaml.safe_load(f))


# ---------------------------------------------------------------------------
# outcomes.yml (logical -> concrete dataset mapping)
# ---------------------------------------------------------------------------


class OutcomeMapping(BaseModel):
    primary: str
    sensitivity: list[str] = Field(default_factory=list)
    replication: list[str] = Field(default_factory=list)


class OutcomePanel(BaseModel):
    panels: dict[str, list[str]]
    mappings: dict[str, OutcomeMapping]
    af_to_outcome_paths: list[str] = Field(default_factory=list)

    @classmethod
    def load(cls, path: str | Path) -> "OutcomePanel":
        with open(path, "r", encoding="utf-8") as f:
            return cls.model_validate(yaml.safe_load(f))
