"""Spec section 18 — Test 3 (palindromic variant handling) plus core QC."""

import pandas as pd

from afshf.harmonize import (
    STANDARD_COLUMNS,
    HarmonizationReport,
    ensure_standard_columns,
    harmonize_sumstats,
    is_palindromic,
    palindromic_decision,
)


def _row(**overrides) -> dict:
    base = {
        "trait": "AF", "source": "TEST", "variant_id": "rs1", "rsid": "rs1",
        "chr": "1", "pos": 100000, "build": "GRCh38",
        "effect_allele": "A", "other_allele": "G",
        "beta": 0.1, "se": 0.02, "or": None, "pval": 1e-5, "eaf": 0.3,
        "n": 10000, "n_cases": 1000, "n_controls": 9000, "ancestry": "EUR",
    }
    base.update(overrides)
    return base


def test_palindromic_detector():
    assert is_palindromic("A", "T")
    assert is_palindromic("C", "G")
    assert not is_palindromic("A", "G")
    assert not is_palindromic("a", "t")  # uppercase required


def test_palindromic_decision_missing_eaf():
    assert palindromic_decision("A", "T", eaf=None) == "drop_missing_eaf"


def test_palindromic_decision_keep_when_eaf_far_from_half():
    assert palindromic_decision("A", "T", eaf=0.12) == "keep"


def test_palindromic_decision_drop_when_eaf_near_half():
    assert palindromic_decision("A", "T", eaf=0.49) == "drop_ambiguous"


def test_harmonize_drops_invalid_alleles_and_palindromic_unresolved():
    rows = [
        _row(rsid="rs_keep_1"),
        _row(rsid="rs_invalid", effect_allele="N", other_allele="G"),
        _row(rsid="rs_pal_drop", effect_allele="A", other_allele="T", eaf=None),
        _row(rsid="rs_pal_drop_ambig", effect_allele="C", other_allele="G", eaf=0.49),
        _row(rsid="rs_pal_keep", effect_allele="A", other_allele="T", eaf=0.12),
        _row(rsid="rs_missing_beta", beta=None),
    ]
    df = pd.DataFrame(rows)
    out, report = harmonize_sumstats(df)

    rsid_kept = set(out["rsid"])
    assert "rs_keep_1" in rsid_kept
    assert "rs_pal_keep" in rsid_kept
    assert "rs_invalid" not in rsid_kept
    assert "rs_pal_drop" not in rsid_kept
    assert "rs_pal_drop_ambig" not in rsid_kept
    assert "rs_missing_beta" not in rsid_kept

    assert isinstance(report, HarmonizationReport)
    assert report.n_in == 6
    assert report.dropped_invalid_alleles == 1
    assert report.dropped_missing_betas == 1
    assert report.dropped_palindromic_unresolved == 2


def test_harmonize_drops_low_maf():
    rows = [
        _row(rsid="rs_common", eaf=0.3),
        _row(rsid="rs_rare", eaf=0.005),
    ]
    df = pd.DataFrame(rows)
    out, report = harmonize_sumstats(df, drop_maf_below=0.01)
    assert "rs_rare" not in set(out["rsid"])
    assert "rs_common" in set(out["rsid"])
    assert report.dropped_low_maf == 1


def test_ensure_standard_columns_raises_on_missing():
    df = pd.DataFrame({"variant_id": ["rs1"], "beta": [0.1]})
    try:
        ensure_standard_columns(df)
    except ValueError as e:
        assert "missing required columns" in str(e)
    else:
        raif then AssertionError("expected ValueError")


def test_standard_columns_is_complete():
    # Sanity: spec section 2 lists exactly these columns.
    expected = {
        "trait", "source", "variant_id", "rsid", "chr", "pos", "build",
        "effect_allele", "other_allele", "beta", "se", "or", "pval", "eaf",
        "n", "n_cases", "n_controls", "ancestry",
    }
    assert set(STANDARD_COLUMNS) == expected
