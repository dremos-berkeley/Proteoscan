import filecmp
import os

import numpy as np
import pandas as pd
import pytest

from conftest import SAMPLES, require_ref
from peaksparity import peaks_export as pe
from peaksparity.peaks_export import (
    TableLayout,
    canonical_columns,
    empty_table,
    format_value,
    java_float_to_string,
    read_export,
    read_table,
    write_export,
    write_table,
)

CSV_TABLES = list(pe.TABLES)


def _ref_path(d, exp, table):
    return os.path.join(d, pe.TABLES[table].format(sample=exp.sample_name))


@pytest.mark.parametrize("k", SAMPLES)
def test_round_trip_byte_identical(k, tmp_path):
    d = require_ref(k)
    exp = read_export(d)
    assert set(exp.tables()) == set(CSV_TABLES)
    assert exp.sample_name == f"Sample {k}"
    paths = write_export(exp, tmp_path)
    assert set(paths) == set(CSV_TABLES)
    for table, p in paths.items():
        ref = _ref_path(d, exp, table)
        assert os.path.basename(p) == os.path.basename(ref)
        assert filecmp.cmp(ref, p, shallow=False), f"{table} differs from {ref}"
    # read -> write -> read gives identical DataFrames and layouts
    exp2 = read_export(tmp_path)
    for table in CSV_TABLES:
        a, b = getattr(exp, table), getattr(exp2, table)
        pd.testing.assert_frame_equal(a, b, check_exact=True)
        assert exp.layouts[table] == exp2.layouts[table]


@pytest.mark.parametrize("k", SAMPLES)
def test_static_spec_matches_files(k, tmp_path):
    """The built-in PEAKS layout alone (no captured layout) reproduces every file."""
    d = require_ref(k)
    exp = read_export(d)
    for table, layout in exp.layouts.items():
        assert not layout.notes, (table, layout.notes)
        spec = TableLayout.for_columns(table, layout.columns)
        assert spec.kinds == layout.kinds, table
        assert spec.header_quoted == layout.header_quoted, table
        for c, q_file, q_spec in zip(layout.columns, layout.quotes, spec.quotes):
            assert q_file == q_spec, (table, c)
    exp.layouts = {}
    for table, p in write_export(exp, tmp_path).items():
        assert filecmp.cmp(_ref_path(d, exp, table), p, shallow=False), table


def test_canonical_headers_match_reference():
    d = require_ref(1)
    for table in CSV_TABLES:
        path = os.path.join(d, pe.TABLES[table].format(sample="Sample 1"))
        cols = list(read_table(path, table)[0].columns)
        assert canonical_columns(table, "Sample 1") == cols, table


def test_empty_table_header_line(tmp_path):
    d = require_ref(1)
    for table in CSV_TABLES:
        ref = os.path.join(d, pe.TABLES[table].format(sample="Sample 1"))
        df = empty_table(table, "Sample 1")
        p = tmp_path / f"{table}.csv"
        write_table(df, p, table)
        with open(ref) as fh:
            first = fh.readline()
        assert p.read_text() == first, table


def test_pandas_read_csv_frames_write_identically(tmp_path):
    """A pipeline DataFrame with ordinary pandas dtypes (NaN for missing text,
    float for int columns with gaps) is written in exact PEAKS format."""
    d = require_ref(1)
    exp = read_export(d, tables=[])
    for table in CSV_TABLES:
        ref = os.path.join(d, pe.TABLES[table].format(sample=exp.sample_name))
        df = pd.read_csv(ref)
        p = tmp_path / os.path.basename(ref)
        write_table(df, p, table)
        assert filecmp.cmp(ref, p, shallow=False), table


def test_dtypes_sample1():
    d = require_ref(1)
    exp = read_export(d)
    psms = exp.psms
    assert psms["Scan"].dtype == np.int64
    assert psms["z"].dtype == np.int64
    assert psms["-10LgP"].dtype == np.float64
    assert psms["Area"].isna().sum() > 0  # empty Area where no feature
    assert (psms.loc[psms["Area"].isna(), "Feature Id"] == 0).all()
    assert exp.layouts["psms"].entry("MS2 Correlation")[2] == "NaN"  # missing -> "NaN" in file
    assert psms["PTM"].map(type).eq(str).all()  # missing text is "", not NaN
    assert exp.proteins["Top"].dtype == bool
    assert psms["Accession"].str.contains(";").any()
    assert not psms["Accession"].str.contains(":").any()
    assert exp.peptides["Accession"].str.contains(":").any()
    assert exp.ptm_profile.loc[exp.ptm_profile["Name"] == "Deamidation (NQ)", "Mass (Da)"].iloc[0] == 0.98


def test_params_sample1():
    d = require_ref(1)
    p = read_export(d, tables=[]).params
    assert p["enzyme"] == "Trypsin"
    assert p["digest_mode"] == "Semi-Specific"
    assert p["precursor_tolerance"] == 15.0 and p["precursor_tolerance_unit"] == "ppm"
    assert p["fragment_tolerance"] == 0.5 and p["fragment_tolerance_unit"] == "Da"
    assert p["max_missed_cleavages"] == 2
    assert (p["peptide_length_min"], p["peptide_length_max"]) == (6, 45)
    assert p["fixed_modifications"] == [{"name": "Carbamidomethylation", "delta": 57.02}]
    assert p["variable_modifications"] == [
        {"name": "Deamidation (NQ)", "delta": 0.98},
        {"name": "Oxidation (M)", "delta": 15.99},
    ]
    assert p["max_variable_ptm_per_peptide"] == 3
    assert p["database"] == "Homo_sapiens"
    assert p["searched_entries"] == 20597
    assert p["deep_learning_boost"] is True
    assert p["fdr_estimation"] is True
    assert p["peptide_min_10lgp"] == 28.36
    assert p["protein_min_10lgp"] == 15
    assert p["protein_min_unique_peptides"] == 1
    assert p["fdr_psm_pct"] == 0.0
    assert p["fdr_peptide_pct"] == 1.0
    assert p["fdr_protein_group_pct"] == 0.0
    assert p["denovo_only_spectra"] == 35082
    assert p["statistics"]["Sample 1"]["psms"] == 794


def test_ms2_correlation_nan():
    exp = read_export(require_ref(4), tables=["psms"])
    assert exp.psms["MS2 Correlation"].isna().sum() > 0  # written as "NaN", never empty


def test_nullable_int_columns():
    exp = read_export(require_ref(4), tables=["modified_peptides"])
    mp = exp.modified_peptides
    assert mp["Length"].dtype == "Int64"  # empty when no unmodified counterpart
    assert mp["Length"].isna().sum() > 0
    assert (mp.loc[mp["Length"].isna(), "Peptide"] == "").all()


@pytest.mark.parametrize("k,thr,psm_fdr", [(2, 26.79, 0.1), (3, 24.08, 0.2), (4, 23.02, 0.5)])
def test_params_other_samples(k, thr, psm_fdr):
    p = read_export(require_ref(k), tables=[]).params
    assert p["peptide_min_10lgp"] == thr
    assert p["fdr_psm_pct"] == psm_fdr
    # de novo only count in the summary equals the rows of the denovoOnly export
    exp = read_export(require_ref(k), tables=["denovo_only"])
    assert p["denovo_only_spectra"] == len(exp.denovo_only)


@pytest.mark.parametrize("x,s", [
    (2.95233664e8, "2.95233664E8"),   # JDK<=18 exact-integer path, not shortest
    (3.0000496e7, "3.0000496E7"),
    (129882.984, "129882.984"),
    (4531940.0, "4531940.0"),
    (114.257935, "114.257935"),
    (491.0, "491.0"),
    (0.0005, "5.0E-4"),
    (0.01, "0.01"),
])
def test_java_float(x, s):
    assert java_float_to_string(x) == s


@pytest.mark.parametrize("kind,x,s", [
    ("sci", 6.58e9, "6.58E9"), ("sci", 1.1e5, "1.1E5"), ("sci", 3e7, "3E7"),
    ("sci", 0.5, "5E-1"), ("sci", 99850.0, "9.98E4"), ("sci", 99950.0, "1E5"),  # HALF_EVEN
    ("dec2", 0.984016, ".98"), ("dec2", 57.021464, "57.02"), ("dec2", -17.026549, "-17.03"),
    ("dec2", 57.0, "57"),
    ("fixed2", 0.125, "0.13"), ("fixed1", 0.15, "0.2"), ("fixed4", 14.4835, "14.4835"),
    ("fixed1", -0.04, "-0.0"),
    ("int", 12571, "12571"), ("bool", True, "true"),
])
def test_format_value(kind, x, s):
    assert format_value(kind, x) == s


def test_absent_files_are_none(tmp_path):
    d = require_ref(1)
    exp = read_export(d, tables=["psms", "proteins"])
    write_export(exp, tmp_path)
    exp2 = read_export(tmp_path)
    assert exp2.psms is not None and exp2.proteins is not None
    assert exp2.peptides is None and exp2.features is None and exp2.denovo_only is None
    assert exp2.params == {}
