"""Tests for peaksparity.proteins (PEAKS protein inference imitation)."""
import importlib.util
import os

import numpy as np
import pandas as pd
import pytest

from conftest import REF_ROOT, REPO, SAMPLES, require_ref
from peaksparity.peaks_export import canonical_columns, format_value, java_float_to_string
from peaksparity.proteins import (
    Enzyme,
    average_mass,
    digest_positions,
    extended_peptide,
    infer_proteins,
    map_peptides,
    protein_peptide_mass,
    protein_score,
    read_protein_db,
)

# ---------------------------------------------------------------------------
# toy database
# ---------------------------------------------------------------------------
P1, P2, P3, P4 = "AAAPEPTK", "LLSMAYK", "GGFFYTPR", "DLIELNK"
TOY = {
    "sp|P00001|A_HUMAN Protein A OS=Homo sapiens OX=9606 GN=GA PE=1 SV=1": "MK" + P1 + P2 + P3 + "WWW",
    "sp|P00002|B_HUMAN Protein B OS=Homo sapiens OX=9606 GN=GB PE=1 SV=1": "MR" + P1 + P2 + "HHH",
    "sp|P00003|C_HUMAN Protein C OS=Homo sapiens OX=9606 PE=1 SV=1": "MK" + P3 + P2 + P1 + "CCC",
    "sp|P00004|D_HUMAN Protein D OS=Homo sapiens OX=9606 GN=GD PE=1 SV=1": "MR" + P3 + P4 + "YYY",
    "sp|P00005|E_HUMAN Protein E OS=Homo sapiens OX=9606 GN=GE PE=1 SV=1": "MK" + P4.replace("I", "L") + "FFF",
    "sp|P00006|F_HUMAN Protein F OS=Homo sapiens OX=9606 GN=GF PE=1 SV=1": "MR" + P1 + "QQQ" + P3,
    "sp|P00007|G_HUMAN Protein G OS=Homo sapiens OX=9606 GN=GG PE=1 SV=1": "MAQ" + P2 + "PW",
    "sp|P00008|H_HUMAN Protein H OS=Homo sapiens OX=9606 GN=GH PE=1 SV=1": "MK" + P3 + "NNN",
    "#DECOY#sp|P00001|A_HUMAN decoy": "KPTPEAAAKYAMSLLRPTYFFGGWWWM",
}


@pytest.fixture
def toy_fasta(tmp_path):
    p = tmp_path / "toy.fasta"
    p.write_text("".join(f">{h}\n{s}\n" for h, s in TOY.items()))
    return str(p)


def toy_peptides():
    return pd.DataFrame({
        "Peptide": [P1, P2, "LLSM(+15.99)AYK", P3, P4],
        "-10LgP": [50.0, 40.0, 35.0, 30.0, 20.0],
        "Area Sample 1": [1e6, 2e6, 3e5, 4e5, np.nan],
        "#Feature": [1, 1, 1, 1, 0],
        "PTM": ["", "", "Oxidation (M)", "", ""],
        "Scan": [1, 2, 3, 4, 5],
    })


def test_read_protein_db(toy_fasta):
    db = read_protein_db(toy_fasta)
    assert len(db) == 8  # decoy skipped
    i = db.index["P00001|A_HUMAN"]
    assert db.gene[i] == "GA"
    assert db.description[i] == "Protein A OS=Homo sapiens OX=9606 GN=GA PE=1 SV=1"
    assert db.gene[db.index["P00003|C_HUMAN"]] == ""


def test_enzyme_rules():
    e = Enzyme()
    seq = "MAQ" + P2 + "PW"
    assert not e.valid(seq, 3, len(P2))  # Q| ... K|P : no enzymatic terminus
    assert e.valid("MR" + P1 + "QQQ", 2, len(P1))  # R|AAAPEPTK|Q : N-term tryptic
    assert e.valid("M" + P1 + "Q", 1, len(P1))  # initiator Met removed counts as N-term
    assert not Enzyme(specificity="specific").valid("MR" + P1 + "QQQ", 3, len(P1) - 1)
    # K before P is not a cleavage site (also not a missed cleavage)
    assert not e._site("AKPA", 2)
    assert Enzyme(max_missed=0).valid("KAAKAAKAA", 1, 3) is True
    assert Enzyme(max_missed=0).valid("KAAKAAKAA", 1, 6) is False


def test_map_and_digest_positions(toy_fasta):
    db = read_protein_db(toy_fasta)
    hits = map_peptides(db, {P2, P4})
    accs = {db.accession[i] for i, _, _ in hits[P2]}
    assert "P00007|G_HUMAN" not in accs  # non-enzymatic occurrence
    assert accs == {"P00001|A_HUMAN", "P00002|B_HUMAN", "P00003|C_HUMAN"}
    il = {db.accession[i]: exact for i, _, exact in hits[P4]}
    assert il == {"P00004|D_HUMAN": True, "P00005|E_HUMAN": False}
    dp = digest_positions(db, ["LLSM(+15.99)AYK"])
    row = dp[dp.Accession == "P00001|A_HUMAN"].iloc[0]
    assert (row.Start, row.End) == (11, 17)
    assert row["Extended Peptide"] == "K.LLSM(+15.99)AYK.G"


def test_extended_peptide_termini():
    assert extended_peptide("MAGSGAGVR", "MAGSGAGVRC", 0, 9) == "MAGSGAGVR.C"
    assert extended_peptide("TLYGFGG", "RTLYGFGG", 1, 7) == "R.TLYGFGG"


def test_protein_score():
    assert protein_score([30.0]) == pytest.approx(30.0)
    assert protein_score([30.0, 50.0, 40.0]) == pytest.approx(50 + 20 + 10, abs=1e-5)
    # float32 accumulation (reference value from the Jurkat export, P14678 = 41.878494)
    assert java_float_to_string(protein_score([41.878494])) == "41.878494"


def test_masses():
    ccl21 = ("MAQSLALSLLILVLAFGIPRTQGSDGGAQDCCLKYSQRKIPAKVVRSYRKQEPSLGCSIPAILFLPRKRSQAELCADPKE"
             "LWVQQLMQHLDKTPSPQKPAQGCRKDRGASKTGKKGKGSKGCKRTERSQTPKGP")
    assert java_float_to_string(average_mass(ccl21)) == "14646.14"
    assert format_value("fixed4", protein_peptide_mass("AILFLPR")) == "828.5221"
    assert format_value("fixed4", protein_peptide_mass("LVTDLTK")) == "788.4643"  # HALF_DOWN tie
    assert format_value("fixed4", protein_peptide_mass("GENLVSM(+15.99)TVEGPPPK")) == "1569.7708"


def test_infer_toy(toy_fasta):
    pr, pp = infer_proteins(toy_peptides(), toy_fasta)
    assert list(pr.columns) == canonical_columns("proteins", ("Sample 1",))
    assert list(pp.columns) == canonical_columns("protein_peptides", ("Sample 1",))
    r = pr.set_index("Accession")
    # H only has the shared GGFFYTPR -> #Unique 0 -> filtered; G has no enzymatic hit
    assert set(r.index) == {f"P0000{i}|{c}_HUMAN" for i, c in zip("123456", "ABCDEF")}
    g1 = {"P00001|A_HUMAN", "P00002|B_HUMAN", "P00003|C_HUMAN", "P00006|F_HUMAN"}
    assert set(r.index[r["Protein Group"] == 1]) == g1
    assert set(r.index[r["Protein Group"] == 2]) == {"P00004|D_HUMAN", "P00005|E_HUMAN"}
    assert set(r.index[r.Top]) == {"P00001|A_HUMAN", "P00003|C_HUMAN", "P00004|D_HUMAN"}
    a = r.loc["P00001|A_HUMAN"]
    assert a["-10LgP"] == pytest.approx(80.0, abs=1e-4)  # PTM form of LLSMAYK counts once
    assert (a["#Peptides"], a["#Unique"]) == (3, 2)
    assert a["Coverage(%)"] == pytest.approx(100 * 23 / 28)
    assert a["Area Sample 1"] == pytest.approx(3.3e6, rel=1e-6)  # unique rows only
    assert a["#Spec Sample 1"] == 4.0  # one per peptide row without psms_df
    assert a["PTM"] == "Oxidation (M)"
    assert a["Gene"] == "GA"
    d = r.loc["P00004|D_HUMAN"]
    assert (d["#Peptides"], d["#Unique"]) == (2, 1)
    assert d["-10LgP"] == pytest.approx(40.0, abs=1e-4)
    e = r.loc["P00005|E_HUMAN"]
    assert e["-10LgP"] == pytest.approx(20.0) and not e.Top
    u = pp.set_index(["Accession", "Peptide"])["Unique"]
    assert not u[("P00001|A_HUMAN", P3)] and u[("P00001|A_HUMAN", P1)]
    # E contains the peptide only with I->L; listed with its input sequence at the I/L position
    er = pp[pp.Accession == "P00005|E_HUMAN"].iloc[0]
    assert (er.Peptide, er.Start, er.End) == (P4, 3, 9)


def test_infer_score_threshold(toy_fasta):
    pr, _ = infer_proteins(toy_peptides(), toy_fasta, score_threshold=25.0)
    assert "P00005|E_HUMAN" not in set(pr.Accession)  # -10lgP 20 < 25


def test_infer_accession_mapping(toy_fasta):
    pep = toy_peptides()
    pep["Accession"] = ["P00001|A_HUMAN:P00002|B_HUMAN"] + [""] * 4
    pr, pp = infer_proteins(pep, toy_fasta, mapping="accession")
    assert set(pr.Accession) == {"P00001|A_HUMAN", "P00002|B_HUMAN"}
    assert pr.Top.sum() == 2  # same family set


# ---------------------------------------------------------------------------
# reference exports: exact reproduction from PEAKS's own peptides
# ---------------------------------------------------------------------------
LOCAL_FASTA = "/drive_a/drew/openms-data/260921-DrewR_Adam_91526/fasta/human_sp.fasta"


def _eval_module():
    spec = importlib.util.spec_from_file_location("eval_proteins", os.path.join(REPO, "scripts", "eval_proteins.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def merged_db(tmp_path_factory):
    if not os.path.exists(LOCAL_FASTA) or not os.path.isdir(REF_ROOT):
        pytest.skip("reference FASTA / exports not available")
    mod = _eval_module()
    path = tmp_path_factory.mktemp("fasta") / "merged.fasta"
    mod.build_merged_fasta(REF_ROOT, LOCAL_FASTA, str(path))
    return mod, read_protein_db(str(path))


EXACT_PROTEIN_FIELDS = ["Top", "-10LgP", "-10LgP (2 dp)", "#Peptides", "#Unique", "Coverage(%)", "#Spec",
                        "Area", "Average Mass", "PTM", "Gene", "Description"]
EXACT_PP_FIELDS = ["Top", "Unique", "Extended Peptide", "-10LgP", "Mass", "Start", "End", "Area",
                   "#Feature", "PTM", "AScore", "Scan"]


@pytest.mark.parametrize("k", SAMPLES)
def test_reference_accession_mapping_exact(merged_db, k):
    require_ref(k)
    mod, db = merged_db
    res, pr, pp = mod.run_sample(REF_ROOT, k, db, "accession", True)
    g = res["groups"]
    assert g["identical_membership"] == g["ref"] == g["test"]
    assert res["proteins"]["ref_only"] == [] and res["proteins"]["test_only"] == []
    for f in EXACT_PROTEIN_FIELDS:
        assert res["protein_fields"][f] == 1.0, (f, res["protein_mismatches"].get(f))
    ppr = res["protein_peptides"]
    assert ppr["ref"] == ppr["test"] == ppr["shared"]
    for f in EXACT_PP_FIELDS:
        assert res["protein_peptide_fields"][f] == 1.0, (f, res["protein_peptide_mismatches"].get(f))
    # group numbers: PEAKS's key (rows of the first protein, descending) is never violated;
    # the order inside ties is not reproducible
    assert res["numbering_rule"]["violations"] == 0
    assert res["group_numbering_pair_order"] > 0.8


@pytest.mark.parametrize("k", SAMPLES)
def test_reference_fasta_mapping(merged_db, k):
    """Mapping by sequence: identical except where our (newer) UniProt release has entries
    or sequences PEAKS's database lacked (IGK_HUMAN; NACA_HUMAN and LAP2B_HUMAN)."""
    require_ref(k)
    mod, db = merged_db
    res, pr, pp = mod.run_sample(REF_ROOT, k, db, "fasta", True)
    allowed = {"P0DOX7|IGK_HUMAN", "Q13765|NACA_HUMAN", "P42167|LAP2B_HUMAN"}
    assert res["proteins"]["ref_only"] == []
    assert set(res["proteins"]["test_only"]) <= allowed
    g = res["groups"]
    assert g["ref"] - g["identical_membership"] <= len(res["proteins"]["test_only"])
    for f in ("Top", "-10LgP", "#Peptides", "#Unique", "Coverage(%)", "#Spec", "Area", "Average Mass"):
        assert res["protein_fields"][f] == 1.0, (f, res["protein_mismatches"].get(f))
