"""Lightweight I/O helpers used across the pipeline."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import yaml


def read_yaml(path: str | Path) -> Any:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def write_yaml(obj: Any, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(obj, f, sort_keys=False)


def read_table(path: str | Path) -> pd.DataFrame:
    """Read a tab-separated or parquet table inferred by extension."""
    p = Path(path)
    if p.suffix in {".parquet", ".pq"}:
        return pd.read_parquet(p)
    return pd.read_csv(p, sep="\t", compression="infer", low_memory=False)


def write_table(df: pd.DataFrame, path: str | Path) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.suffix in {".parquet", ".pq"}:
        df.to_parquet(p, index=False)
    else:
        df.to_csv(p, sep="\t", index=False, compression="infer")
