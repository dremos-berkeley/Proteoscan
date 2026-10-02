import pandas as pd
import pytest

from conftest import SAMPLES, require_ref
from proteoscan import peaks_export as pe
from proteoscan.seqnorm import (
    canonical,
    equal_modseq,
    format_delta,
    il_collapse,
    mod_delta,
    mods_of,
    parse,
    strip_mods,
)


@pytest.mark.parametrize("inp,out", [
    # PEAKS
    ("KVPQ(+0.98)VSTPTLVEVSR", "KVPQ(+0.98)VSTPTLVEVSR"),
    ("N(+.98)GK", "N(+0.98)GK"),
    ("YGM(+15.99)KNALLVR", "YGM(+15.99)KNALLVR"),
    ("FSSC(+57.02)GGGGGSFGAGGGFGSR", "FSSC(+57.02)GGGGGSFGAGGGFGSR"),
    ("(+42.01)PEPTIDE", "(+42.01)PEPTIDE"),
    # ProForma names
    ("PEPTM[Oxidation]IDE", "PEPTM(+15.99)IDE"),
    ("PEPTC[Carbamidomethyl]IDE", "PEPTC(+57.02)IDE"),
    ("PEPN[Deamidated]IDE", "PEPN(+0.98)IDE"),
    ("[Acetyl]-PEPTIDE", "(+42.01)PEPTIDE"),
    ("PEPTC[Propionamide]IDE", "PEPTC(+71.04)IDE"),
    ("[Gln->pyro-Glu]-QPEPTIDE", "(-17.03)QPEPTIDE"),
    ("[Glu->pyro-Glu]-EPEPTIDE", "(-18.01)EPEPTIDE"),
    ("Q[Pyro-glu from Q]PEPTIDE", "Q(-17.03)PEPTIDE"),
    ("E[Pyro-glu from E]PEPTIDE", "E(-18.01)PEPTIDE"),
    ("[Carbamyl]-PEPTIDEK", "(+43.01)PEPTIDEK"),
    ("PEPTIDEK[Carbamyl]", "PEPTIDEK(+43.01)"),
    # ProForma masses and Unimod accessions
    ("PEPTM[+15.9949]IDE", "PEPTM(+15.99)IDE"),
    ("PEPTM[+15.994915]IDE", "PEPTM(+15.99)IDE"),
    ("PEPN[+0.984016]IDE", "PEPN(+0.98)IDE"),
    ("[+42.010565]-PEPTIDE", "(+42.01)PEPTIDE"),
    ("[-17.026549]-QPEPTIDE", "(-17.03)QPEPTIDE"),
    ("PEPTM[UNIMOD:35]IDE", "PEPTM(+15.99)IDE"),
    ("PEPTIDE-[Amidated]", "PEPTIDE-(-0.98)"),
    # other dialects
    ("PEPTM(Oxidation)IDE", "PEPTM(+15.99)IDE"),
    (".(Acetyl)PEPTM(Oxidation)IDE", "(+42.01)PEPTM(+15.99)IDE"),
    ("n[42.0106]PEPTC[57.0215]K", "(+42.01)PEPTC(+57.02)K"),
    ("_(Acetyl (Protein N-term))PEPTM(Oxidation (M))IDE_", "(+42.01)PEPTM(+15.99)IDE"),
    ("K.PEPTM(+15.99)IDE.R", "PEPTM(+15.99)IDE"),
    ("-.PEPTIDE.R", "PEPTIDE"),
    # two mods on one site are summed
    ("M[Oxidation][Acetyl]K", "M(+58.01)K"),
])
def test_canonical(inp, out):
    assert canonical(inp) == out


def test_il_distinct_and_site_matters():
    assert not equal_modseq("PEPTIDE", "PEPTLDE")
    assert not equal_modseq("VEIIAN(+0.98)DQGNR", "VEILAN(+0.98)DQGNR")
    assert not equal_modseq("PQ(+0.98)VSTPTLVEVSN", "PQVSTPTLVEVSN(+0.98)")
    assert not equal_modseq("(+42.01)SPEPTIDE", "S(+42.01)PEPTIDE")
    assert not equal_modseq("PEPTM(+15.99)IDE", "PEPTMIDE")
    assert equal_modseq("PEPTM[Oxidation]IDE", "PEPTM(+15.99)IDE")
    assert equal_modseq("PEPTM[+15.9949]IDE", "PEPTM(+15.99)IDE")
    assert equal_modseq("[Acetyl]-PEPTIDE", "(+42.01)PEPTIDE")


def test_strip_and_mods():
    assert strip_mods("(+42.01)PEPTM(+15.99)IDE") == "PEPTMIDE"
    assert strip_mods("[Acetyl]-PEPTM[Oxidation]IDE-[Amidated]") == "PEPTMIDE"
    assert mods_of("(+42.01)PEPTM(+15.99)IDE") == [(0, "N-term", "+42.01"), (5, "M", "+15.99")]
    assert mods_of("PEPTIDE-[Amidated]") == [(8, "C-term", "-0.98")]
    assert mods_of("PEPTIDE") == []
    assert parse("N(+0.98)K") == ("NK", ((1, 0.98),))


def test_format_delta_and_names():
    assert format_delta(15.9949) == "+15.99"
    assert format_delta(0.984016) == "+0.98"
    assert format_delta(-17.026549) == "-17.03"
    assert format_delta(71.037114) == "+71.04"
    assert format_delta(43.005814) == "+43.01"
    assert format_delta(0.0) == "+0.00"
    for name, d in [("Carbamidomethyl", "+57.02"), ("Oxidation", "+15.99"), ("Deamidated", "+0.98"),
                    ("Acetyl", "+42.01"), ("Propionamide", "+71.04"), ("Pyro-glu from Q", "-17.03"),
                    ("Pyro-glu from E", "-18.01"), ("Carbamyl", "+43.01"), ("Deamidation (NQ)", "+0.98"),
                    ("Oxidation (M)", "+15.99"), ("Carbamidomethylation", "+57.02")]:
        assert format_delta(mod_delta(name)) == d, name


def test_fixed_and_nterm_options():
    assert canonical("PEPTCIDE", fixed={"C": "Carbamidomethyl"}) == "PEPTC(+57.02)IDE"
    assert canonical("PEPTC(+57.02)IDE", fixed={"C": 57.021464}) == "PEPTC(+57.02)IDE"
    assert canonical("(+42.01)MPEP", nterm_to_first_residue=True) == "M(+42.01)PEP"


@pytest.mark.parametrize("bad", ["PEPTX[Foo]IDE", "PEP(+15.99", "", "pepTIDE", "PEPT[UNIMOD:999999]"])
def test_invalid(bad):
    with pytest.raises(ValueError):
        canonical(bad)


def test_il_collapse():
    assert il_collapse("LIFAGK") == "LLFAGK"


@pytest.mark.parametrize("k", SAMPLES)
def test_peaks_strings_are_already_canonical(k):
    """Every peptide string PEAKS exported is a fixed point of canonical()."""
    d = require_ref(k)
    exp = pe.read_export(d, tables=["psms", "peptides", "denovo_only", "modified_peptides"], verify=False)
    peps = set(exp.psms["Peptide"]) | set(exp.peptides["Peptide"]) | set(exp.denovo_only["Peptide"])
    peps |= set(exp.modified_peptides["Peptide Mod"])
    peps |= {p for p in exp.modified_peptides["Peptide"] if p}
    bad = [p for p in peps if canonical(p) != p]
    assert not bad, bad[:10]
    # PEAKS "tag" column = peptide with the same mod annotation
    assert (exp.psms["tag(>=0.0%)"] == exp.psms["Peptide"]).mean() > 0.5
