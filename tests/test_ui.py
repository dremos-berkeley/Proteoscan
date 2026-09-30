"""Tests for peaksparity.ui (self-contained HTML results viewer) on a small synthetic export."""
import json
import os
import re
import shutil
import subprocess
from html.parser import HTMLParser

import numpy as np
import pandas as pd
import pytest

from peaksparity import ui
from peaksparity.peaks_export import PeaksExport, canonical_columns, write_export

SRC = "RUN_1.raw"
PROT_A = "MKAAAPEPTKLLSMAYKGGFFYTPRWWW"
PROT_B = "MRAAAPEPTKDLIELNKHHH"
FASTA = {
    "sp|P00001|A_HUMAN Protein A OS=Homo sapiens OX=9606 GN=GA PE=1 SV=1": PROT_A,
    "sp|P00002|B_HUMAN Protein B OS=Homo sapiens OX=9606 GN=GB PE=1 SV=1": PROT_B,
}
# peptide, score, scan, z, accessions
PEPTIDES = [
    ("AAAPEPTK", 50.0, 2, 2, ["P00001|A_HUMAN", "P00002|B_HUMAN"]),
    ("LLSM(+15.99)AYK", 40.0, 3, 2, ["P00001|A_HUMAN"]),
    ("GGFFYTPR", 35.0, 4, 3, ["P00001|A_HUMAN"]),
    ("DLIELN(+0.98)K", 30.0, 5, 2, ["P00002|B_HUMAN"]),
]
MONO = {"A": 71.037114, "K": 128.094963, "P": 97.052764, "E": 129.042593, "T": 101.047679,
        "L": 113.084064, "S": 87.032028, "M": 131.040485, "Y": 163.063329, "G": 57.021464,
        "F": 147.068414, "R": 156.101111, "D": 115.026943, "I": 113.084064, "N": 114.042927}
PROTON, H2O = 1.007276467, 18.010565


def _fill(table, rows):
    """DataFrame with the PEAKS columns of ``table``; unspecified columns get neutral values."""
    cols = canonical_columns(table, ("Sample 1",))
    out = {}
    for c in cols:
        vals = [r.get(c) for r in rows]
        if all(v is None for v in vals):
            vals = [np.nan] * len(rows)
        out[c] = vals
    return pd.DataFrame(out, columns=cols)


def _mass(stripped, extra=0.0):
    return sum(MONO[a] for a in stripped) + H2O + extra


def _synthetic_export(d):
    psm_rows, pep_rows, pp_rows, prot_rows = [], [], [], []
    for pep, score, scan, z, accs in PEPTIDES:
        stripped = re.sub(r"\([^)]*\)", "", pep)
        extra = 15.994915 if "+15.99" in pep else (0.984016 if "+0.98" in pep else 0.0)
        m = _mass(stripped, extra)
        ptm = "Oxidation (M)" if "+15.99" in pep else ("Deamidation (NQ)" if "+0.98" in pep else "")
        common = {"Peptide": pep, "-10LgP": score, "Mass": m, "Length": len(stripped), "ppm": 1.0,
                  "m/z": (m + z * PROTON) / z, "z": z, "RT": 10.0 + scan, "Scan": scan, "PTM": ptm,
                  "Tag Length": 0, "CAA (%)": 0.0, "Delta RT": 0.0, "MS2 Correlation": 0.9}
        for extra_scan in (0, 100):  # a second, lower-scoring PSM of each peptide
            psm_rows.append(dict(common, **{"-10LgP": score - extra_scan / 10, "Scan": scan + extra_scan,
                                            "Source File": SRC, "Feature Id": scan, "raw m/z": 0.0,
                                            "raw RT": 10.0, "Accession": ";".join(accs), "AScore": ""}))
        pep_rows.append(dict(common, **{"Source File": SRC, "Accession": ":".join(accs), "#Feature": 1,
                                        "#Feature Sample 1": 1, "Found By": "DB Search", "AScore": ""}))
        for acc in accs:
            seq = PROT_A if acc.startswith("P00001") else PROT_B
            st = seq.find(stripped)
            pp_rows.append(dict(common, **{"Protein Group": 1 if acc.startswith("P00001") else 2,
                                           "Accession": acc, "Gene": "G" + acc[5], "Top": True,
                                           "Extended Peptide": pep, "Unique": len(accs) == 1,
                                           "Source File": SRC, "#Feature": 1, "Start": st + 1,
                                           "End": st + len(stripped), "Found By": "DB Search", "AScore": ""}))
    for g, acc, gene in ((1, "P00001|A_HUMAN", "GA"), (2, "P00002|B_HUMAN", "GB")):
        prot_rows.append({"Protein Group": g, "Top": True, "Accession": acc, "Gene": gene, "-10LgP": 80.0 - g,
                          "Coverage(%)": 50.0, "Coverage(%) Sample 1": 50.0, "#Peptides": 2, "#Unique": 1,
                          "#Spec Sample 1": 4.0, "PTM": "", "Average Mass": 3000.0,
                          "Description": f"Protein {gene[-1]}, test OS=Homo sapiens GN={gene}"})
    exp = PeaksExport(psms=_fill("psms", psm_rows), peptides=_fill("peptides", pep_rows),
                      proteins=_fill("proteins", prot_rows), protein_peptides=_fill("protein_peptides", pp_rows),
                      sample_name="Sample 1")
    for name, df in exp.tables().items():
        for c in df.columns:
            if df[c].dtype == object and df[c].isna().any():
                df[c] = df[c].fillna("")
    write_export(exp, d)
    return d


def _fragments(stripped, z=1, extra=None):
    """b and y 1+ m/z (extra: {1-based position: delta})."""
    extra = extra or {}
    r = [MONO[a] + extra.get(i + 1, 0.0) for i, a in enumerate(stripped)]
    b = [sum(r[:i]) + PROTON for i in range(1, len(r))]
    y = [sum(r[len(r) - i:]) + H2O + PROTON for i in range(1, len(r))]
    return b, y


def _write_mzml(path):
    oms = pytest.importorskip("pyopenms")
    exp = oms.MSExperiment()
    for scan in range(1, 7):
        s = oms.MSSpectrum()
        s.setNativeID(f"controllerType=0 controllerNumber=1 scan={scan}")
        s.setRT(60.0 * scan)
        if scan == 1:
            s.setMSLevel(1)
            mz = np.array([400.0, 500.0, 600.0])
        else:
            s.setMSLevel(2)
            p = oms.Precursor()
            p.setMZ(500.0)
            p.setCharge(2)
            s.setPrecursors([p])
            if scan == 2:  # b/y ions of AAAPEPTK plus noise, 200 peaks in total
                b, y = _fragments("AAAPEPTK")
                noise = list(np.linspace(150.0, 900.0, 200 - len(b) - len(y)) + 0.37)
                mz = np.array(sorted(b + y + noise))
            else:
                mz = np.linspace(200.0, 800.0, 50)
        s.set_peaks((mz, np.linspace(1000.0, 5000.0, len(mz))))
        exp.addSpectrum(s)
    oms.MzMLFile().store(str(path), exp)
    return str(path)


class _Scripts(HTMLParser):
    def __init__(self):
        super().__init__()
        self.scripts, self.cur, self.external = [], None, []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == "script":
            if a.get("src"):
                self.external.append(a["src"])
            self.cur = [a, []]
        if tag == "link" and a.get("href"):
            self.external.append(a["href"])

    def handle_data(self, data):
        if self.cur is not None:
            self.cur[1].append(data)

    def handle_endtag(self, tag):
        if tag == "script" and self.cur is not None:
            self.scripts.append((self.cur[0], "".join(self.cur[1])))
            self.cur = None


def _parse(path):
    p = _Scripts()
    with open(path, encoding="utf-8") as fh:
        p.feed(fh.read())
    by_id = {a.get("id"): body for a, body in p.scripts if a.get("id")}
    code = [body for a, body in p.scripts if not a.get("id")]
    return p, json.loads(by_id["pp-data"]), json.loads(by_id["pp-spectra"]), code


@pytest.fixture
def export_dir(tmp_path):
    return _synthetic_export(str(tmp_path / "export"))


@pytest.fixture
def fasta(tmp_path):
    p = tmp_path / "db.fasta"
    p.write_text("".join(f">{h}\n{s}\n" for h, s in FASTA.items()))
    return str(p)


def test_report_without_spectra(export_dir, tmp_path):
    out = tmp_path / "r.html"
    info = ui.build_report(export_dir, str(out))
    assert out.exists() and info["bytes"] == os.path.getsize(out)
    parser, data, spectra, code = _parse(out)
    assert parser.external == []  # nothing loaded from the network
    assert spectra == {}
    assert len(code) == 1 and "DataTable" in code[0]
    t = data["tables"]
    assert set(t) == {"psms", "peptides", "proteins", "protein_peptides"}
    assert len(t["psms"]["rows"]) == 8 and len(t["peptides"]["rows"]) == 4 and len(t["proteins"]["rows"]) == 2
    assert t["psms"]["columns"] == canonical_columns("psms")
    assert all(len(r) == len(t["psms"]["columns"]) for r in t["psms"]["rows"])
    assert data["counts"]["peptide_rows"] == 4 and data["counts"]["protein_groups"] == 2
    assert data["frag_tol"] == {"value": 0.5, "unit": "Da", "source": "default (no fragment tolerance in the export)"}
    assert data["sequences"] == {}  # no FASTA
    assert data["missed_cleavages"]["computed"]["0"] == 4
    assert data["ptm_profile"]["rows"]  # computed (no db.ptmProfile.csv)
    # exact modification masses for the spectrum annotation
    assert data["pepinfo"]["LLSM(+15.99)AYK"] == ["LLSMAYK", [[4, 15.994915]]]
    assert data["pepinfo"]["DLIELN(+0.98)K"][1] == [[6, 0.984016]]
    assert data["concordance"] is None


def test_report_with_fasta_mzml_ref(export_dir, fasta, tmp_path):
    mzml = _write_mzml(tmp_path / "RUN_1.mzML")
    out = tmp_path / "full.html"
    ui.build_report(export_dir, str(out), mzml=[mzml], fasta=fasta, ref=export_dir, top_n=150)
    _p, data, spectra, _code = _parse(out)
    # best PSM of each peptide = scans 2..5; the second PSMs (scan + 100) are not embedded
    assert sorted(spectra) == [f"RUN_1|{s}" for s in (2, 3, 4, 5)]
    e = spectra["RUN_1|2"]
    assert len(e["mz"]) == 150 and e["mz"] == sorted(e["mz"]) and max(e["i"]) == 100.0
    assert e["pz"] == 2 and e["id"].endswith("scan=2")
    assert data["spectra_info"]["embedded"] == 4
    assert data["sequences"] == {"P00001|A_HUMAN": PROT_A, "P00002|B_HUMAN": PROT_B}
    m = data["concordance"]["metrics"]["levels"]
    assert m["psm"]["recall"] == 1.0 and m["peptide"]["precision"] == 1.0 and m["protein_group"]["recall"] == 1.0
    assert data["tables"]["psms"]["columns"][0] == "In ref"
    assert all(r[0] is True for r in data["tables"]["psms"]["rows"])


def test_spectrum_size_cap(export_dir, tmp_path):
    mzml = _write_mzml(tmp_path / "RUN_1.mzML")
    out = tmp_path / "cap.html"
    ui.build_report(export_dir, str(out), mzml=[mzml], max_spectra_mb=0.004)
    _p, data, spectra, _code = _parse(out)
    info = data["spectra_info"]
    assert info["embedded"] == len(spectra) < 4
    assert info["skipped_size_cap"] == 4 - len(spectra)
    assert "RUN_1|3" not in spectra or "RUN_1|2" in spectra  # highest-scoring peptide first


def test_cli_subcommand(export_dir, tmp_path):
    from peaksparity.cli import main
    out = tmp_path / "cli.html"
    main(["ui", "--export", export_dir, "--out", str(out), "--frag-tol", "20", "--frag-unit", "ppm"])
    _p, data, _s, _c = _parse(out)
    assert data["frag_tol"]["value"] == 20.0 and data["frag_tol"]["unit"] == "ppm"


def test_embedded_json_is_script_safe(tmp_path, export_dir):
    s = ui._dumps({"x": "</script><!--  "})
    assert "</" not in s and "<!--" not in s and json.loads(s)["x"] == "</script><!--  "


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_javascript_runs_in_node(export_dir, tmp_path):
    """The page script parses, and its fragment calculator matches b/y masses computed here."""
    out = tmp_path / "r.html"
    ui.build_report(export_dir, str(out))
    _p, data, _s, code = _parse(out)
    js = code[0]
    start = js.index("function theoreticalIons")
    end = js.index("function ionLabel")
    prog = (
        "const PROTON = 1.007276467, H2O = 18.010565, NH3 = 17.026549;\n"
        + js[js.index("const AA = {"): js.index("};", js.index("const AA = {")) + 2] + "\n"
        + js[start:end]
        + "\nconst info = " + json.dumps(data["pepinfo"]["LLSM(+15.99)AYK"]) + ";\n"
        + "const r = theoreticalIons(info, 2, {b: true, y: true, two: false, loss: false});\n"
        + "console.log(JSON.stringify(r.ions.map(i => [i.ser, i.k, i.mz])));\n"
    )
    script = tmp_path / "frag.js"
    script.write_text(prog)
    res = subprocess.run(["node", str(script)], capture_output=True, text=True, check=True)
    ions = json.loads(res.stdout)
    b, y = _fragments("LLSMAYK", extra={4: 15.994915})
    got_b = [mz for s, _k, mz in ions if s == "b"]
    got_y = [mz for s, _k, mz in ions if s == "y"]
    assert np.allclose(got_b, b, atol=1e-6) and np.allclose(got_y, y, atol=1e-6)
    check = tmp_path / "page.js"
    check.write_text(js)
    subprocess.run(["node", "--check", str(check)], check=True)
