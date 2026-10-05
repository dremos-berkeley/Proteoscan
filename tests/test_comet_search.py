"""Tests for proteoscan.comet_search: static/variable modification handling in comet.params and
explicit fixed mods in read_comet output (docs/STATUS.md "Corrections", P1)."""
import re
import shutil
from pathlib import Path

import pandas as pd
import pytest

import proteoscan.comet_search as cs
from proteoscan.config import Modification, SearchConfig

TEMPLATE = Path(__file__).parent / "data" / "comet.params.template"  # `comet -p`, Comet 2026.01 rev. 1


@pytest.fixture
def template(monkeypatch):
    text = TEMPLATE.read_text()
    monkeypatch.setattr(cs, "_default_params", lambda: text)
    return text


def _params(text):
    out = {}
    for line in text.splitlines():
        m = re.match(r"^(\w+) = (.*?)\s*(#.*)?$", line)
        if m:
            out[m.group(1)] = m.group(2)
    return out


def _cfg(**kw):
    kw.setdefault("database", "/x/db.fasta")
    return SearchConfig(**kw)


PROP_VARIABLE = [Modification("Carbamidomethylation", "C"), Modification("Deamidation", "NQ"),
                 Modification("Oxidation", "M"), Modification("Propionamide", "C")]


def test_template_has_nonzero_defaults(template):
    # the reason for the fix: comet -p ships a static +57 on C and a variable oxidation
    p = _params(template)
    assert float(p["add_C_cysteine"]) == pytest.approx(57.021464)
    assert p["variable_mod01"].startswith("15.9949 M")


def test_no_fixed_mods_zeroes_every_static_mod(template):
    p = _params(cs.build_comet_params(_cfg(fixed_modifications=[], variable_modifications=PROP_VARIABLE)))
    adds = {k: v for k, v in p.items() if k.startswith("add_")}
    assert len(adds) == 30  # 26 residue letters + 4 termini
    assert all(float(v) == 0.0 for v in adds.values()), {k: v for k, v in adds.items() if float(v) != 0}


def test_prop_config_variable_slots(template):
    p = _params(cs.build_comet_params(_cfg(fixed_modifications=[], variable_modifications=PROP_VARIABLE)))
    var = {k: v for k, v in p.items() if k.startswith("variable_mod")}
    used = [v for v in var.values() if not v.startswith("0.0 X")]
    assert sorted(float(v.split()[0]) for v in used) == pytest.approx([0.984016, 15.994915, 57.021464, 71.037114])
    assert any(v.startswith("71.037114 C ") for v in used)  # propionamide searched as variable on C
    assert any(v.startswith("57.021464 C ") for v in used)
    unused = [v for v in var.values() if v.startswith("0.0 X")]
    assert len(var) == 5 and len(unused) == 1 and unused[0] == cs.EMPTY_VARIABLE_MOD


def test_no_variable_mods_drops_default_oxidation(template):
    p = _params(cs.build_comet_params(_cfg(variable_modifications=[])))
    assert all(v == cs.EMPTY_VARIABLE_MOD for k, v in p.items() if k.startswith("variable_mod"))


def test_fixed_cam_only_sets_cysteine(template):
    p = _params(cs.build_comet_params(_cfg()))  # default config: fixed CAM, variable deamidation + oxidation
    nonzero = {k: v for k, v in p.items() if k.startswith("add_") and float(v) != 0.0}
    assert nonzero == {"add_C_cysteine": "57.021464"}


def test_fixed_terminal_and_rare_residues(template):
    cfg = _cfg(fixed_modifications=[Modification("TMT", "N-term", 229.162932), Modification("TMT", "K", 229.162932),
                                    Modification("x", "U", 1.5)], variable_modifications=[])
    assert cs.fixed_mod_deltas(cfg) == {"N-term": 229.162932, "K": 229.162932, "U": 1.5}
    p = _params(cs.build_comet_params(cfg))
    assert float(p["add_Nterm_peptide"]) == pytest.approx(229.162932)
    assert float(p["add_K_lysine"]) == pytest.approx(229.162932)
    assert float(p["add_U_selenocysteine"]) == 1.5
    assert float(p["add_C_cysteine"]) == 0.0


def test_unmappable_fixed_residue_raises():
    with pytest.raises(ValueError):
        cs.fixed_mod_deltas(_cfg(fixed_modifications=[Modification("weird", "1", 1.0)]))


def test_fixed_and_variable_on_same_residue_raises():
    with pytest.raises(ValueError, match="overlaps fixed"):
        _cfg(fixed_modifications=[Modification("Carbamidomethylation", "C")],
             variable_modifications=[Modification("Carbamidomethylation", "C")])
    with pytest.raises(ValueError):
        _cfg(fixed_modifications=[Modification("Carbamidomethylation", "C")],
             variable_modifications=[Modification("Propionamide", "C")])
    # the propionamide config (no fixed mods, CAM and propionamide variable on C) is valid
    _cfg(fixed_modifications=[], variable_modifications=PROP_VARIABLE)


def _comet_txt(tmp_path):
    cols = ["scan", "num", "charge", "exp_neutral_mass", "calc_neutral_mass", "e-value", "xcorr", "delta_cn",
            "sp_score", "ions_matched", "ions_total", "plain_peptide", "modified_peptide", "prev_aa", "next_aa",
            "protein", "retention_time_sec", "sp_rank"]
    rows = [[10, 1, 2, 1000.0, 1000.0, 1e-5, 3.0, 0.2, 500, 10, 20, "ACDCK", "K.ACDCK.L", "K", "L", "sp|P1|X", 600, 1],
            [11, 1, 2, 900.0, 900.0, 1e-3, 2.0, 0.1, 300, 8, 18, "PEPMK", "R.PEPM[15.9949]K.A", "R", "A", "DECOY_sp|P2|Y", 610, 1]]
    p = tmp_path / "x.txt"
    with open(p, "w") as fh:
        fh.write("CometVersion 2026.01\n")
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    return p


def test_read_comet_fixed_mods_explicit(tmp_path):
    txt = _comet_txt(tmp_path)
    cfg = _cfg()
    df = cs.read_comet(txt, fixed_mods=cs.fixed_mod_deltas(cfg))
    assert df["peptide_sage"].tolist() == ["AC[+57.0215]DC[+57.0215]K", "PEPM[+15.9949]K"]
    assert df["is_decoy"].tolist() == [False, True]


def test_read_comet_no_fixed_mods_bare(tmp_path):
    df = cs.read_comet(_comet_txt(tmp_path), fixed_mods=cs.fixed_mod_deltas(_cfg(fixed_modifications=[])))
    assert df["peptide_sage"].tolist() == ["ACDCK", "PEPM[+15.9949]K"]


def test_read_comet_terminal_fixed_mod(tmp_path):
    df = cs.read_comet(_comet_txt(tmp_path), fixed_mods={"N-term": 229.162932, "C-term": 1.0})
    assert df["peptide_sage"].tolist() == ["[+229.1629]-ACDCK-[+1.0000]", "[+229.1629]-PEPM[+15.9949]K-[+1.0000]"]
    assert cs._add_terminal("[+42.0106]-PEPK", "N-term", 1.0) == "[+43.0106]-PEPK"


@pytest.mark.skipif(shutil.which("comet") is None, reason="comet not on PATH")
def test_real_comet_has_every_key():
    real = _params(cs._default_params())
    tmpl = _params(TEMPLATE.read_text())
    keys = {k for k in tmpl if k.startswith(("add_", "variable_mod"))}
    assert keys <= set(real)
    text = cs.build_comet_params(_cfg(fixed_modifications=[], variable_modifications=PROP_VARIABLE))
    assert all(float(v) == 0.0 for k, v in _params(text).items() if k.startswith("add_"))
