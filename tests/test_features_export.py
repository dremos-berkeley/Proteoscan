"""features_export (PEAKS <Sample>.db.features.csv) and the feature-derived peptide columns."""
import os

import numpy as np
import pandas as pd
import pytest

from conftest import require_ref
from peaksparity.export_builder import build_peptides
from peaksparity.features_export import PROTON, build_features_table, feature_concordance
from peaksparity.peaks_export import canonical_columns, read_export, write_table

SRC = "run1.raw"


def _features():
    return pd.DataFrame({
        "feature_id": [1, 2, 3, 4, 5],
        "mz": [400.2, 500.3, 600.4, 700.5, 800.6],
        "charge": [2, 2, 3, 2, 2],
        "rt": [10.0, 20.0, 30.0, 40.0, 50.0],
        "rt_begin": [9.8, 19.5, 29.9, 39.0, 49.0],
        "rt_end": [10.3, 20.6, 30.2, 41.0, 51.0],
        "area": [1.0e5, 2.5e6, 3.3e4, 7.7e7, 1.0e3],
        "quality": [0.9, 0.8, 0.7, 0.6, 0.5],
    })


def _psms():
    rows = [
        # feature 2: two scans, best 45.1; I/L variants at the best score (split by protein mapping)
        dict(Peptide="TLVDLAK", **{"-10LgP": 45.1}, Mass=744.4432, Scan=100, Accession="B|B_HUMAN", fid=2),
        dict(Peptide="TLVDIAK", **{"-10LgP": 45.1}, Mass=744.4432, Scan=100, Accession="A|A_HUMAN;C|C_HUMAN", fid=2),
        dict(Peptide="TLVDIAK", **{"-10LgP": 30.0}, Mass=744.4432, Scan=101, Accession="A|A_HUMAN;C|C_HUMAN", fid=2),
        # feature 4: one PSM
        dict(Peptide="PEPTIDEK", **{"-10LgP": 25.5}, Mass=927.4549, Scan=200, Accession="D|D_HUMAN", fid=4),
        # featureless PSM
        dict(Peptide="SAMPLER", **{"-10LgP": 33.0}, Mass=802.4122, Scan=300, Accession="E|E_HUMAN", fid=0),
    ]
    df = pd.DataFrame(rows).rename(columns={"fid": "Feature Id"})
    df["Source File"] = SRC
    return df


def test_features_table_db_rows():
    out = build_features_table(_features(), _psms(), sample_name="S")
    assert list(out.columns) == canonical_columns("features")
    assert out["Feature Id"].tolist() == [2, 4]                   # only features with a PSM, sorted
    r = out.set_index("Feature Id")
    assert r.loc[2, "DB Peptide"] == "TLVDIAK;TLVDLAK"             # I/L variants, ';'-joined
    assert r.loc[2, "-10LgP"] == pytest.approx(45.1)                # max over the feature's PSMs
    assert r.loc[2, "Accession"] == "A|A_HUMAN:C|C_HUMAN"           # first variant's proteins, ':'
    assert r.loc[4, "Accession"] == "D|D_HUMAN"
    assert (out["Quality"] == 0).all() and (out["ALC (%)"] == 0).all()
    assert (out["Denovo Peptide"] == "").all()
    assert r.loc[2, "RT Begin"] == 19.5 and r.loc[2, "RT End"] == 20.6
    assert r.loc[2, "Area"] == pytest.approx(2.5e6)
    assert (out["Source File"] == SRC).all()
    # ppm: uncalibrated feature m/z against the exact peptide mass
    from pyteomics import mass as pm
    m = pm.fast_mass("PEPTIDEK")
    assert r.loc[4, "ppm"] == pytest.approx(((700.5 - PROTON) * 2 - m) / m * 1e6, abs=1e-6)


def test_features_table_denovo_rows():
    dn = pd.DataFrame({
        "Feature Id": [2, 1, 1, 3, 0],
        "Peptide": ["XXXX", "LVYATGK", "LVYATGR", "AEQFKK", "NOFEAT"],
        "ALC (%)": [99.0, 60.0, 81.0, 45.0, 90.0],
        "Mass": [700.0, 736.4, 764.4, 750.0, 800.0],
        "ppm": [1.0, 2.0, 3.5, 4.0, 5.0],
        "Accession": ["", "", "P1|X_HUMAN;P2|Y_HUMAN", "", ""],
    })
    out = build_features_table(_features(), _psms(), dn, sample_name="S")
    assert out["Feature Id"].tolist() == [1, 2, 4]   # 3 is below ALC 50; 2 keeps its DB peptide
    r = out.set_index("Feature Id")
    assert r.loc[1, "Denovo Peptide"] == "LVYATGR" and r.loc[1, "ALC (%)"] == 81.0   # best ALC
    assert r.loc[1, "DB Peptide"] == "" and r.loc[1, "-10LgP"] == 0.0
    assert r.loc[1, "ppm"] == 3.5
    assert r.loc[1, "Accession"] == "P1|X_HUMAN:P2|Y_HUMAN"
    assert r.loc[2, "Denovo Peptide"] == ""


def test_features_table_written_like_peaks(tmp_path):
    out = build_features_table(_features(), _psms(), sample_name="S")
    p = tmp_path / "S.db.features.csv"
    write_table(out, str(p), "features")
    lines = p.read_text().splitlines()
    assert lines[0].startswith('"Source File","Feature Id","m/z"')
    assert lines[1].split(",")[:9] == ['"run1.raw"', "2", "500.30000", "2", "20.0000", "19.5000", "20.6000",
                                      "2500000.0", "0.00"]
    assert lines[1].split(",")[9:13] == ["TLVDIAK;TLVDLAK", "45.10", "", "0.0"]


def test_features_table_empty():
    ps = _psms()
    ps["Feature Id"] = 0
    out = build_features_table(_features(), ps)
    assert len(out) == 0 and list(out.columns) == canonical_columns("features")


def test_peptide_area_sums_distinct_features():
    ps = pd.DataFrame({
        "Peptide": ["AAK", "AAK", "AAK", "CCK"], "-10LgP": [40.0, 30.0, 20.0, 25.0],
        "Mass": [1.0] * 4, "Tag Length": 0, "CAA (%)": 0.0, "Length": 3, "Delta RT": 0.0,
        "MS2 Correlation": 0.9, "ppm": 1.0, "m/z": 1.0, "z": 2, "RT": 1.0,
        "Area": [1e6, 1e6, 3e5, np.nan], "Scan": [1, 2, 3, 4], "Source File": SRC,
        "Feature Id": [7, 7, 9, 0], "Accession": "X", "PTM": "", "AScore": "",
        "Positional Confidence": "", "tag(>=0.0%)": "",
    })
    pep = build_peptides(ps, "S").set_index("Peptide")
    assert pep.loc["AAK", "Area S"] == pytest.approx(1.3e6)    # feature 7 counted once
    assert pep.loc["AAK", "#Feature"] == 2
    assert np.isnan(pep.loc["CCK", "Area S"]) and pep.loc["CCK", "#Feature"] == 0


def test_feature_concordance_identity():
    t = build_features_table(_features(), _psms())
    res = feature_concordance(t, t)
    assert res["recall"] == 1.0 and res["precision"] == 1.0
    assert res["median_log10_ratio"] == 0.0 and res["peptide_agree_IL"] == 1.0
    t2 = t.copy()
    t2["Area"] = t2["Area"] * 2
    assert feature_concordance(t, t2)["median_log10_ratio"] == pytest.approx(np.log10(2))


@pytest.mark.parametrize("k", [1, 4])
def test_reproduces_peaks_features_from_its_psms(k):
    d = require_ref(k)
    ex = read_export(d, tables=["psms", "peptides", "features", "denovo_only"], verify=False)
    if ex.features is None:
        pytest.skip("features.csv absent")
    out = build_features_table(ex.features, ex.psms, ex.denovo_only, ex.sample_name)
    assert len(out) == len(ex.features)
    assert out["Feature Id"].tolist() == ex.features["Feature Id"].tolist()
    m = ex.features.merge(out, on="Feature Id", suffixes=("_r", "_t"))
    for c in ("-10LgP", "Area", "RT Begin", "RT End", "ALC (%)"):
        assert np.allclose(m[c + "_r"].astype(float), m[c + "_t"].astype(float))
    sets = [set(str(a).split(";")) == set(str(b).split(";"))
            for a, b in zip(m["DB Peptide_r"].fillna(""), m["DB Peptide_t"].fillna(""))]
    assert np.mean(sets) >= 0.999
    assert (m["Accession_r"].fillna("") == m["Accession_t"].fillna("")).mean() >= 0.999
    assert (m["Denovo Peptide_r"].fillna("") == m["Denovo Peptide_t"].fillna("")).mean() == 1.0
    assert np.mean(np.abs(m["ppm_r"].astype(float) - m["ppm_t"].astype(float)) <= 0.15) >= 0.998
    # the peptide table's Area / #Feature come from the distinct features of its PSMs
    pep = build_peptides(ex.psms, ex.sample_name).set_index("Peptide")
    ref = ex.peptides.set_index("Peptide")
    a = f"Area {ex.sample_name}"
    j = ref.join(pep[[a, "#Feature"]], rsuffix="_t")
    assert (j["#Feature"] == j["#Feature_t"]).mean() == 1.0
    ok = np.isclose(j[a], j[a + "_t"], rtol=6e-3) | (j[a].isna() & j[a + "_t"].isna())
    assert ok.mean() >= 0.995
