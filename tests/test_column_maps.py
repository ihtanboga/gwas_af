"""Column-map alias resolution + effect-scale derivations."""

from math import isclose, log
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from afshf.column_maps import (
    apply_column_map,
    derive_beta_from_or,
    derive_beta_from_z,
    derive_or_from_beta,
    derive_pval_from_log10p,
    derive_se_from_or_scale,
    load_column_maps,
    resolve_columns,
)


CONFIG = Path(__file__).resolve().parents[1] / "config" / "column_maps.yml"
FIXTURES = Path(__file__).resolve().parent / "synthetic_data"


@pytest.fixture(scope="module")
def registry():
    return load_column_maps(CONFIG)


# ---------------------------------------------------------------------------
# Alias resolution per source
# ---------------------------------------------------------------------------


def test_load_column_maps_has_all_sources(registry):
    expected = {
        "CVDKP_AF", "CVDKP_HERMES", "CVDKP_MEGASTROKE",
        "FinnGen", "GWAS_Catalog", "UKBPPP", "deCODE", "GTEx", "eQTLGen",
    }
    assert expected.issubset(set(registry.sources))


def test_resolve_cvdkp_af_aliases(registry):
    df = pd.read_csv(FIXTURES / "cvdkp_af_mini.tsv", sep="\t")
    resolved = resolve_columns(df, registry.get("CVDKP_AF"))
    assert resolved["chr"] == "CHR"
    assert resolved["pos"] == "BP"
    assert resolved["effect_allele"] == "Allele1"
    assert resolved["other_allele"] == "Allele2"
    assert resolved["beta"] == "Effect"
    assert resolved["se"] == "StdErr"
    assert resolved["pval"] == "P-value"
    assert resolved["eaf"] == "Freq1"
    assert resolved["rsid"] == "rsID"


def test_resolve_finngen_aliases(registry):
    df = pd.read_csv(FIXTURES / "finngen_af_mini.tsv", sep="\t")
    resolved = resolve_columns(df, registry.get("FinnGen"))
    assert resolved["chr"] == "#chrom"
    assert resolved["effect_allele"] == "alt"
    assert resolved["other_allele"] == "ref"
    assert resolved["beta"] == "beta"
    assert resolved["se"] == "sebeta"
    assert resolved["eaf"] == "af_alt"


def test_resolve_ukbppp_aliases(registry):
    df = pd.read_csv(FIXTURES / "ukbppp_pcsk9_mini.tsv", sep="\t")
    resolved = resolve_columns(df, registry.get("UKBPPP"))
    assert resolved["chr"] == "CHROM"
    assert resolved["pos"] == "GENPOS"
    assert resolved["effect_allele"] == "ALLELE1"
    assert resolved["other_allele"] == "ALLELE0"
    assert resolved["beta"] == "BETA"
    assert resolved["se"] == "SE"
    assert resolved["eaf"] == "A1FREQ"
    # UKB-PPP ships LOG10P, not P; log10p resolves separately and pval is derived.
    assert resolved["log10p"] == "LOG10P"
    assert resolved["pval"] is None


def test_resolve_decode_aliases(registry):
    df = pd.read_csv(FIXTURES / "decode_mini.tsv", sep="\t")
    resolved = resolve_columns(df, registry.get("deCODE"))
    assert resolved["chr"] == "Chrom"
    assert resolved["pos"] == "Pos"
    assert resolved["effect_allele"] == "effectAllele"
    assert resolved["beta"] == "Beta"


def test_resolve_eqtlgen_aliases(registry):
    df = pd.read_csv(FIXTURES / "eqtlgen_zonly_mini.tsv", sep="\t")
    resolved = resolve_columns(df, registry.get("eQTLGen"))
    assert resolved["chr"] == "SNPChr"
    assert resolved["pos"] == "SNPPos"
    assert resolved["z"] == "Zscore"
    assert resolved["beta"] is None  # not provided
    assert resolved["se"] is None    # not provided
    assert resolved["eaf"] == "AssessedAlleleFreq"
    assert resolved["n"] == "NrSamples"


# ---------------------------------------------------------------------------
# OR -> beta + GWAS Catalog end-to-end
# ---------------------------------------------------------------------------


def test_derive_beta_from_or_keeps_existing_beta():
    df = pd.DataFrame({"or": [1.5, 2.0], "beta": [0.5, np.nan]})
    out = derive_beta_from_or(df)
    assert out.loc[0, "beta"] == 0.5            # untouched
    assert isclose(out.loc[1, "beta"], log(2.0))


def test_derive_or_from_beta():
    df = pd.DataFrame({"beta": [log(2.0), 0.0]})
    out = derive_or_from_beta(df)
    assert isclose(out.loc[0, "or"], 2.0)
    assert isclose(out.loc[1, "or"], 1.0)


def test_apply_column_map_gwas_catalog_or_to_beta(registry):
    raw = pd.read_csv(FIXTURES / "gwas_catalog_or_mini.tsv", sep="\t")
    out = apply_column_map(raw, registry.get("GWAS_Catalog"),
                            trait="atrial_fibrillation",
                            source="GWAS_Catalog")
    assert "beta" in out.columns
    assert "or" in out.columns
    expected_beta = float(np.log(1.154))
    matched = out.loc[out["rsid"] == "rs2200733", "beta"].iloc[0]
    assert isclose(matched, expected_beta, abs_tol=1e-6)
    assert "build" in out.columns and out["build"].iloc[0] == "GRCh38"
    assert "ancestry" in out.columns


def test_derive_se_from_or_scale_only_when_requested():
    df = pd.DataFrame({"or": [2.0, 1.5], "se": [0.4, 0.3]})
    out = derive_se_from_or_scale(df)
    assert isclose(out.loc[0, "se"], 0.2)        # 0.4 / 2.0
    assert isclose(out.loc[1, "se"], 0.2, abs_tol=1e-9)  # 0.3 / 1.5


# ---------------------------------------------------------------------------
# Z -> beta (eQTLGen)
# ---------------------------------------------------------------------------


def test_derive_beta_from_z_eqtlgen_formula():
    """beta = z * se where var_g = 2*MAF*(1-MAF)*(N + z^2)."""
    df = pd.DataFrame({
        "z": [-12.0, 0.5, np.nan],
        "eaf": [0.015, 0.30, 0.10],
        "n": [31684, 31684, 31684],
    })
    out = derive_beta_from_z(df)
    # Row 0: large negative effect -> negative beta with small SE
    z, maf, n = -12.0, 0.015, 31684
    var_g = 2.0 * maf * (1 - maf) * (n + z ** 2)
    expected_se = 1.0 / (var_g ** 0.5)
    expected_beta = z * expected_se
    assert isclose(out.loc[0, "beta"], expected_beta, abs_tol=1e-9)
    assert isclose(out.loc[0, "se"], expected_se, abs_tol=1e-9)
    # Row 2 has missing z -> beta/se stay NaN
    assert pd.isna(out.loc[2, "beta"])
    assert pd.isna(out.loc[2, "se"])


def test_apply_column_map_eqtlgen_end_to_end(registry):
    raw = pd.read_csv(FIXTURES / "eqtlgen_zonly_mini.tsv", sep="\t")
    out = apply_column_map(raw, registry.get("eQTLGen"),
                            source="eQTLGen", trait="PCSK9_eQTL")
    assert {"beta", "se", "z"}.issubset(out.columns)
    rs = out.set_index("rsid")
    assert rs.loc["rs11591147", "z"] == -12.0
    # The strongest Z should yield the largest |beta|
    assert abs(rs.loc["rs11591147", "beta"]) > abs(rs.loc["rs562556", "beta"])
    # Build/ancestry from metadata
    assert out["build"].iloc[0] == "GRCh37"
    assert out["ancestry"].iloc[0] == "EUR"


# ---------------------------------------------------------------------------
# LOG10P -> P (UKB-PPP)
# ---------------------------------------------------------------------------


def test_derive_pval_from_log10p_basic():
    df = pd.DataFrame({"log10p": [10.0, 0.0]})
    out = derive_pval_from_log10p(df)
    assert isclose(out.loc[0, "pval"], 1e-10)
    assert isclose(out.loc[1, "pval"], 1.0)


def test_apply_column_map_ukbppp_recovers_pval(registry):
    raw = pd.read_csv(FIXTURES / "ukbppp_pcsk9_mini.tsv", sep="\t")
    out = apply_column_map(raw, registry.get("UKBPPP"),
                            source="UKBPPP", trait="PCSK9_pQTL")
    assert "pval" in out.columns
    rs = out.set_index("rsid")
    assert rs.loc["rs11591147", "pval"] < 1e-50      # very strong signal
    assert isclose(rs.loc["rs562556", "pval"], 10 ** -0.5, abs_tol=1e-12)
    assert "beta" in out.columns
    assert rs.loc["rs11591147", "beta"] == -1.234
