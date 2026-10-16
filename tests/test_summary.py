"""Tests for proteoscan.summary (PEAKS-style HTML export pages) on a small synthetic export."""
import os
import re
import struct

import numpy as np
import pandas as pd
import pytest

from proteoscan import summary
from proteoscan.config import SearchConfig
from proteoscan.peaks_export import parse_summary_html

from test_ui import FASTA, PEPTIDES, _fragments, _synthetic_export, _write_mzml

pytest.importorskip("matplotlib")


def _scored(path, n_decoy=30):
    """psms_scored.parquet: the 8 exported target PSMs plus low-scoring targets and decoys."""
    rows = []
    for pep, score, scan, z, accs in PEPTIDES:
        stripped = re.sub(r"\([^)]*\)", "", pep)
        rows.append(dict(source_file="RUN_1.mzML", scan=scan, precursor_idx=0, charge=float(z), peptide_sage=pep,
                         stripped=stripped, proteins=";".join("sp|" + a for a in accs), is_decoy=False,
                         minus10lgp=score, precursor_ppm=1.0))
    rng = np.random.default_rng(1)
    for i in range(n_decoy):
        seq = "".join(rng.choice(list("ACDEFGHKLMNPQRSTVWY"), 9))
        rows.append(dict(source_file="RUN_1.mzML", scan=1000 + i, precursor_idx=0, charge=2.0, peptide_sage=seq,
                         stripped=seq, proteins=f"#DECOY#sp|D{i:05d}|DEC_HUMAN", is_decoy=True,
                         minus10lgp=float(rng.uniform(0, 25)), precursor_ppm=float(rng.uniform(-15, 15))))
        rows.append(dict(source_file="RUN_1.mzML", scan=2000 + i, precursor_idx=0, charge=2.0, peptide_sage=seq[::-1],
                         stripped=seq[::-1], proteins="sp|P00009|X_HUMAN", is_decoy=False,
                         minus10lgp=float(rng.uniform(0, 25)), precursor_ppm=float(rng.uniform(-15, 15))))
    pd.DataFrame(rows).to_parquet(path)
    return str(path)


@pytest.fixture
def run(tmp_path):
    """Pipeline run-directory layout: <run>/export, psms_scored.parquet, <stem>/precursors|ms1_features."""
    rd = tmp_path / "run"
    exp = _synthetic_export(str(rd / "export"))
    _scored(rd / "psms_scored.parquet")
    (rd / "RUN_1").mkdir()
    # scans 2 and 3 have two precursor entries (chimeric), 4 and 5 one
    pd.DataFrame({"scan": [2, 2, 3, 3, 4, 5], "precursor_idx": [0, 1, 0, 1, 0, 0]}).to_parquet(rd / "RUN_1" / "precursors.parquet")
    pd.DataFrame({"feature_id": range(1, 11), "mz": np.linspace(400, 900, 10), "rt": np.linspace(5, 50, 10),
                  "area": np.logspace(3, 7, 10)}).to_parquet(rd / "RUN_1" / "ms1_features.parquet")
    fa = tmp_path / "db.fasta"
    fa.write_text("".join(f">{h}\n{s}\n" for h, s in FASTA.items()))
    cfg = SearchConfig(database=str(fa))
    cfg.to_yaml(rd / "search_config.yml")
    return {"dir": str(rd), "export": exp, "fasta": str(fa), "tmp": tmp_path}


def _text(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def _png_size(path):
    with open(path, "rb") as fh:
        b = fh.read(24)
    assert b[:8] == b"\x89PNG\r\n\x1a\n"
    return struct.unpack(">II", b[16:24])


def test_full_summary(run, tmp_path):
    mzml = _write_mzml(tmp_path / "RUN_1.mzML")
    out = tmp_path / "html"
    info = summary.build_summary(run["export"], str(out), mzml=[mzml], processes=1)
    for f in ("db.summary.html", "db.proteins.html", "db.peptides.html", "export.css"):
        assert (out / f).exists()
    s = _text(out / "db.summary.html")
    # PEAKS section headings, figure and table numbering, in order
    heads = ["1. Notes", "2. Result Statistics", "3. Experiment Control", "4. Other Information"]
    assert [h for h in re.findall(r"<h1>(.*?)</h1>", s)] == heads
    figs = re.findall(r"<b>Figure (\d)\.</b>", s)
    tabs = re.findall(r"<b>Table (\d)\.</b>", s)
    assert figs == ["1", "2", "3", "4", "5", "6"] and tabs == ["1", "2", "3", "4", "5", "6"]
    assert s.count("Not applicable to only one sample") == 3
    assert "False discovery rate (FDR) curve" in s and "Statistics of data." in s
    assert "Parent Mass Error Tolerance: 15.0 ppm" in s and "&nbsp; Oxidation (M) (+15.99)" in s
    # every referenced image exists; no external resources
    for src in re.findall(r'src="([^"]+)"', s):
        assert (out / src).exists(), src
    assert "http" not in re.sub(r"https?://www\.w3\.org", "", s.split("<body>")[0])
    # Table 1 from mzML (1 MS1, 5 MS2), precursor list (2 chimeric), feature table (10)
    t1 = info["table1"]
    assert (t1["ms1_scans"], t1["ms2_scans"], t1["chimera_scans"], t1["features"]) == (1, 5, 2, 10)
    assert (t1["psms"], t1["identified_scans"], t1["peptides"], t1["sequences"]) == (8, 8, 4, 4)
    assert (t1["protein_groups"], t1["proteins_all"], t1["proteins_top"]) == (2, 2, 2)
    # the page round-trips through the PEAKS summary parser
    p = parse_summary_html(str(out / "db.summary.html"))
    assert p["statistics"]["Total"]["chimera_scans"] == 2 and p["statistics"]["Sample 1"]["peptides"] == 4
    assert p["precursor_tolerance"] == 15.0 and p["fragment_tolerance"] == 0.5
    assert p["max_missed_cleavages"] == 2 and p["digest_mode"] == "Semi-Specific"
    assert [m["name"] for m in p["variable_modifications"]] == ["Deamidation (NQ)", "Oxidation (M)"]
    assert p["missed_cleavages"]["Sample 1"]["0"] == 4
    assert p["protein_min_10lgp"] == 15 and isinstance(p["peptide_min_10lgp"], float)
    assert isinstance(p["fdr_peptide_pct"], float)
    # threshold and FDRs from the scored decoys
    assert info["peptide_threshold"] is not None and info["table3"]["seq"].endswith("%")
    for k in ("fdr_curve", "score_histogram", "score_vs_ppm", "feature_mz", "feature_rt", "feature_abundance",
              "precursor_error_histogram", "precursor_error_vs_mz"):
        assert (out / "img" / f"{k}.png").exists(), k
    # one coverage map per protein, one spectrum per peptide
    assert info["coverage_maps"] == 2 and info["spectra"]["drawn"] == 4
    pr = _text(out / "db.proteins.html")
    assert pr.count('class="bar"') == 3 and "total 2 proteins" in pr and "Best Unique PSM (" in pr
    for src in re.findall(r'src="([^"]+)"', pr):
        assert (out / src).exists()
    pe = _text(out / "db.peptides.html")
    assert "total 4 peptides" in pe and len(re.findall(r'<img src="img/spec_', pe)) == 4
    w, h = _png_size(out / "img" / "spec_00001.png")
    assert w == 1200 and h > 200


def test_spectrum_cap_and_degraded_inputs(run, tmp_path):
    """No scored decoys, no run dir, no mzML for MS1: the pages still build and say what is missing."""
    mzml = _write_mzml(tmp_path / "RUN_1.mzML")
    out = tmp_path / "html2"
    exp = _synthetic_export(str(tmp_path / "bare_export"))  # not inside a run directory
    info = summary.build_summary(exp, str(out), mzml=[mzml], max_spectra=2, processes=1)
    s = _text(out / "db.summary.html")
    assert "Not available: needs scored target and decoy PSMs" in s
    assert not (out / "img" / "fdr_curve.png").exists() and (out / "img" / "score_histogram.png").exists()
    assert info["table1"]["features"] is None and "#Features (needs ms1_features.parquet)" in s
    assert "N/A" in info["table3"]["psm"]
    pe = _text(out / "db.peptides.html")
    # 2 top peptides + best unique PSMs of the proteins (AAAPEPTK is shared, so both proteins add a unique one)
    assert info["spectra"]["beyond_cap"] == 2 and "highest-scoring peptides only" in pe
    assert "No FASTA (--fasta): protein coverage maps are not drawn." in s


def test_no_inputs_at_all(tmp_path):
    exp = _synthetic_export(str(tmp_path / "e"))
    info = summary.build_summary(exp, str(tmp_path / "o"), processes=1, images=False)
    assert info["table1"]["ms1_scans"] is None and info["n_images"] >= 3
    s = _text(tmp_path / "o" / "db.summary.html")
    assert "No mzML given" in s and "Table 6." in s


def test_annotate_fragments_finds_all_ions():
    b, y = _fragments("AAAPEPTK")
    mz = np.array(sorted(b + y + [333.33, 777.77]))
    ann = summary.annotate_fragments("AAAPEPTK", [], 2, mz, np.ones(len(mz)), tol=0.02, unit="Da")
    got = {a["ion"] for a in ann}
    assert got == {f"b{i}" for i in range(1, 8)} | {f"y{i}" for i in range(1, 8)}
    assert max(abs(a["err"]) for a in ann) < 1e-3
    # a modification shifts the ions after its position
    b2, _y2 = _fragments("LLSMAYK", extra={4: 15.994915})
    ann = summary.annotate_fragments("LLSMAYK", [(4, 15.994915)], 2, np.array(b2), np.ones(len(b2)), tol=0.01)
    assert {a["ion"] for a in ann} >= {"b4", "b5", "b6"}


def test_cli_subcommand(run, tmp_path):
    from proteoscan.cli import main as cli_main
    out = tmp_path / "cli_out"
    cli_main(["summary", "--export", run["export"], "--out", str(out), "--no-images", "--processes", "1"])
    s = _text(out / "db.summary.html")
    # run-directory inputs were found automatically (search_config.yml, psms_scored, precursor list)
    assert "Search parameters from the pipeline configuration" in s
    assert os.path.exists(out / "img" / "fdr_curve.png")
    assert "Searched Entries: 2" in s
