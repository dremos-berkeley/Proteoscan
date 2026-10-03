import numpy as np
import pytest

from proteoscan.fragpred import DEFAULT_FRAG_TYPES, annotate, fragment_mz, parse_peptide, valid_mask
from proteoscan.rescore_features import parse_peptdeep_spec


def test_parse_peptide_formats_agree():
    a = parse_peptide("PEPC[+57.0215]M[+15.9949]N[+0.9840]K")
    b = parse_peptide("PEPC(+57.02)M(+15.99)N(+0.98)K")
    assert a[0] == b[0] == "PEPCMNK"
    assert np.allclose(a[1], b[1], atol=1e-4)  # PEAKS's rounded deltas snap to the exact masses
    assert a[2] == "Carbamidomethyl@C;Oxidation@M;Deamidated@N" and a[3] == "4;5;6" and a[4]


def test_unknown_mod_flagged():
    seq, res, mods, sites, ok = parse_peptide("PEPS[+79.9663]K")
    assert not ok and abs(res[3] - (87.032028 + 79.9663)) < 1e-6


def test_fragment_mz_matches_alphabase():
    alphabase = pytest.importorskip("alphabase.peptide.fragment")
    import pandas as pd
    s, res, mods, sites, _ = parse_peptide("[+42.0106]-AC[+57.0215]DEM[+15.9949]K")
    df = pd.DataFrame({"sequence": [s], "mods": [mods], "mod_sites": [sites], "charge": [2], "nAA": [len(s)]})
    ref = alphabase.create_fragment_mz_dataframe(df, DEFAULT_FRAG_TYPES)[DEFAULT_FRAG_TYPES].to_numpy()
    assert np.abs(ref - fragment_mz(res, DEFAULT_FRAG_TYPES)).max() < 1e-3


def test_annotate_and_mask():
    theo = np.array([[100.0, 200.0], [300.0, 400.0]])
    mz = np.array([100.3, 199.0, 400.1]); it = np.array([5.0, 7.0, 9.0])
    obs = annotate(theo, mz, it, 0.5)
    assert obs.tolist() == [[5.0, 0.0], [0.0, 9.0]]
    m = valid_mask(theo, ["b_z1", "b_z2"], 1, 150, 1000)
    assert m.tolist() == [[False, False], [True, False]]


def test_spec_parsing():
    assert parse_peptdeep_spec("peptdeep:released") == ({"*": "released"}, "by")
    assert parse_peptdeep_spec("peptdeep:hcd=/a.pth,cid=/b.pth#all") == ({"hcd": "/a.pth", "cid": "/b.pth"}, "all")
