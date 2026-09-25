"""Tests for peaksparity.ptm_tables (db.ptmProfile.csv and db.modifiedPeptide.csv).

Reference tests feed PEAKS's own psms / peptides / proteins / protein-peptides tables back in and
compare the written CSV text with the reference files.  See docs/PTM_TABLES.md.
"""
import csv
import os

import numpy as np
import pandas as pd
import pytest

from conftest import SAMPLES, ref_dir, require_ref
from peaksparity.peaks_export import read_export, write_table
from peaksparity.ptm_tables import (
    build_modified_peptides,
    build_ptm_profile,
    peaks_reference_scores,
    read_isolation_targets,
)

FASTA = os.environ.get("PEAKSPARITY_MERGED_FASTA",
                       "/drive_a/drew/peaks-parity-work/fasta/human_sp_peaks_merged.fasta")
MZML_ROOT = os.environ.get("PEAKSPARITY_MZML_ROOT",
                           "/drive_a/drew/openms-data/260921-DrewR_Adam_91526/mzml")
CACHE = "/drive_a/drew/peaks-parity-work/ptm_tables/cache"

#: RT shift cells that differ by 0.0001 because the CSV RTs carry 4 decimals and PEAKS
#: subtracts the unrounded feature RTs (docs/PTM_TABLES.md, "RT shift").
RT_SHIFT_RESIDUAL = {1: 0, 2: 3, 3: 0, 4: 5}


# ---------------------------------------------------------------------------
# toy data
# ---------------------------------------------------------------------------

def _psm(pep, score, scan, ptm="", ascore="", pc=None, fid=0, rt=10.0, z=2, mz=500.0):
    n = len(pep.replace("(+15.99)", "").replace("(+0.98)", "").replace("(+57.02)", ""))
    return {"Peptide": pep, "-10LgP": score, "Mass": 1000.0, "Length": n, "ppm": 1.0, "m/z": mz,
            "raw m/z": 0.0, "z": z, "RT": rt, "Scan": scan, "Source File": "a.raw", "Feature Id": fid,
            "Accession": "P1|A_HUMAN", "PTM": ptm, "AScore": ascore,
            "Positional Confidence": " ".join(str(x) for x in (pc or [5] * n))}


def _peptides_from(psms):
    rows = []
    for pep, g in psms.groupby("Peptide", sort=True):
        b = g.sort_values("-10LgP", ascending=False, kind="stable").iloc[0]
        rows.append({**b.to_dict(), "Area Sample 1": 1e6, "#Feature": 1, "#Feature Sample 1": 1,
                     "Accession": "P1|A_HUMAN", "Found By": "DB Search"})
    return pd.DataFrame(rows)


@pytest.fixture
def toy():
    psms = pd.DataFrame([
        _psm("PEPM(+15.99)K", 40.0, 10, "Oxidation (M)", "M4:Oxidation (M):1001.00", [5, 5, 5, 5, 5], fid=1),
        _psm("PEPM(+15.99)K", 40.0, 11, "Oxidation (M)", "M4:Oxidation (M):1001.00", [5, 5, 5, 1, 5], fid=1),
        _psm("PEPMK", 35.0, 12, fid=0, rt=9.0),
        _psm("PEPMK", 30.0, 13, fid=2, rt=9.5),
        _psm("AN(+0.98)C(+57.02)K", 30.0, 14, "Deamidation (NQ);Carbamidomethylation",
             "N2:Deamidation (NQ):14.00", [9, 9, 9, 9]),
        _psm("AN(+0.98)C(+57.02)K", 20.0, 15, "Deamidation (NQ);Carbamidomethylation",
             "N2:Deamidation (NQ):30.00", [9, 9, 9, 9]),
        _psm("ANC(+57.02)K", 30.0, 16, "Carbamidomethylation", "", [9, 9, 9, 9]),
    ])
    peptides = _peptides_from(psms)
    pp = pd.DataFrame({
        "Protein Group": [1, 1, 1, 1, 2],
        "Accession": ["P1|A_HUMAN"] * 4 + ["P2|B_HUMAN"],
        "Top": [True, True, True, True, True],
        "Unique": [True, True, True, True, False],
        "Peptide": ["PEPM(+15.99)K", "PEPMK", "AN(+0.98)C(+57.02)K", "ANC(+57.02)K", "PEPM(+15.99)K"],
        "Start": [5, 5, 20, 20, 1],
        "PTM": ["Oxidation (M)", "", "Deamidation (NQ); Carbamidomethylation", "Carbamidomethylation",
                "Oxidation (M)"],
    })
    return psms, peptides, pp


def test_modified_peptides_toy(toy):
    psms, peptides, pp = toy
    mp = build_modified_peptides(psms, peptides, None, pp, "Sample 1", isolation_mz={10: 1.0, 11: 2.0, 12: 3.0})
    assert list(mp["Peptide Mod"]) == ["PEPM(+15.99)K", "AN(+0.98)C(+57.02)K"]
    ox = mp.iloc[0]
    # the "Mod" representative is the LAST tied PSM (scan 11); its PC at M4 is 1 -> not confident
    assert ox["Scan Mod"] == 11 and ox["m/z Mod"] == 2.0
    assert ox["PTM"] == "M4:Oxidation (M):false"
    assert ox["AScore"] == "M4:Oxidation (M):1000.00"
    # unmodified counterpart from peptides; #Feature counts Feature Ids {0, 2} including 0
    assert ox["Peptide"] == "PEPMK" and ox["Scan"] == 12 and ox["#Feature Sample 1"] == 2
    assert ox["m/z"] == 3.0 and ox["RT shift"] == pytest.approx(1.0)
    dm = mp.iloc[1]
    assert dm["PTM"] == "N2:Deamidation (NQ):false;C3:Carbamidomethylation:true"  # AScore 14 < 15
    assert dm["Peptide"] == "ANC(+57.02)K"


def test_ptm_profile_toy(toy):
    psms, peptides, pp = toy
    prof = build_ptm_profile(psms, peptides, None, pp, None, "Sample 1").set_index("Name")
    # all three have one peptide sequence: ties put variable PTMs first, then by name
    assert list(prof.index) == ["Deamidation (NQ)", "Oxidation (M)", "Carbamidomethylation"]
    assert prof.loc["Oxidation (M)", "# peptides"] == 1
    assert prof.loc["Oxidation (M)", "# proteins"] == 1          # B_HUMAN row is not Unique
    assert prof.loc["Oxidation (M)", "# modification sites"] == 1  # A_HUMAN M8
    assert prof.loc["Oxidation (M)", "# confident modification sites"] == 1  # scan 10 is confident
    # deamidation: one PSM has AScore 14 (no), the other 30 with PC 9 (yes) -> confident
    assert prof.loc["Deamidation (NQ)", "# confident modification sites"] == 1
    assert prof.loc["Deamidation (NQ)", "Position"] == "NQ"
    assert prof.loc["Carbamidomethylation", "# peptides"] == 1  # ANCK, counted once
    assert prof.loc["Carbamidomethylation", "Position"] == "C"


# ---------------------------------------------------------------------------
# reference exports
# ---------------------------------------------------------------------------

def _load(k):
    d = require_ref(k)
    if not os.path.exists(os.path.join(d, "db.ptmProfile.csv")):
        pytest.skip("reference ptmProfile not available")
    return d, read_export(d, tables=["psms", "peptides", "proteins", "protein_peptides"])


def _cells(path):
    with open(path, newline="") as fh:
        return list(csv.reader(fh))


@pytest.mark.parametrize("k", SAMPLES)
def test_reference_ptm_profile(k, tmp_path):
    d, e = _load(k)
    fasta = FASTA if os.path.exists(FASTA) else None
    prof = build_ptm_profile(e.psms, e.peptides, e.proteins, e.protein_peptides, fasta, e.sample_name)
    out = tmp_path / "db.ptmProfile.csv"
    write_table(prof, out, "ptm_profile")
    assert out.read_bytes() == open(os.path.join(d, "db.ptmProfile.csv"), "rb").read()


def _cache_dir():
    try:
        os.makedirs(CACHE, exist_ok=True)
        return CACHE if os.access(CACHE, os.W_OK) else None
    except OSError:
        return None


@pytest.mark.parametrize("k", SAMPLES)
def test_reference_modified_peptides(k, tmp_path):
    d, e = _load(k)
    mzml = os.path.join(MZML_ROOT, f"DrewR_Adam_91526_DDA_OT_{k}.mzML")
    if not os.path.exists(mzml):
        pytest.skip(f"mzML not available: {mzml}")
    if not os.path.exists(os.path.join(d, "db.mzidentml_1_2_0.xml")):
        pytest.skip("mzIdentML not available")
    scores = peaks_reference_scores(d, e.psms, cache_dir=_cache_dir())
    iso = read_isolation_targets(mzml)
    mp = build_modified_peptides(e.psms, e.peptides, e.proteins, e.protein_peptides, e.sample_name,
                                 scores=scores, isolation_mz=iso)
    out = tmp_path / "db.modifiedPeptide.csv"
    write_table(mp, out, "modified_peptides")
    ours, ref = _cells(out), _cells(os.path.join(d, "db.modifiedPeptide.csv"))
    assert ours[0] == ref[0] and len(ours) == len(ref)
    shift = ref[0].index("RT shift")
    residual = 0
    for a, b in zip(ours[1:], ref[1:]):
        for j, (x, y) in enumerate(zip(a, b)):
            if x == y:
                continue
            assert j == shift, (ref[0][j], x, y, b[0])
            assert abs(float(x) - float(y)) <= 1.0001e-4
            residual += 1
    assert residual <= RT_SHIFT_RESIDUAL[k]
    # every line except the RT-shift residual rows is byte-identical
    ol = out.read_text().split("\n")
    rl = open(os.path.join(d, "db.modifiedPeptide.csv")).read().split("\n")
    assert sum(x != y for x, y in zip(ol, rl)) == residual


@pytest.mark.parametrize("k", SAMPLES)
def test_reference_modified_peptides_without_extras(k):
    """Without full scores / isolation targets the row SET and the non-positional columns still
    agree; only the tie-dependent choices (representative PSM, row order) and m/z can differ."""
    d, e = _load(k)
    mp = build_modified_peptides(e.psms, e.peptides, e.proteins, e.protein_peptides, e.sample_name)
    ref = read_export(d, tables=["modified_peptides"]).modified_peptides
    assert sorted(mp["Peptide Mod"]) == sorted(ref["Peptide Mod"])
    a = mp.set_index("Peptide Mod").sort_index()
    b = ref.set_index("Peptide Mod").sort_index()
    for c in ("-10LgP Mod", "Mass Mod", "Accession", "Peptide", "-10LgP", "Scan"):
        x, y = a[c].astype(str).replace({"<NA>": "", "nan": ""}), b[c].astype(str).replace({"<NA>": "", "nan": ""})
        assert (x == y).all(), c
