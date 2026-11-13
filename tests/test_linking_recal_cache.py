"""Precursor entry classification, feature_linking, automatic precursor recalibration, cache fingerprints and
entrapment twins in the protein-support feature."""
import json

import numpy as np
import pandas as pd
import pytest

from proteoscan import pipeline as pl
from proteoscan import twopass
from proteoscan.config import ISOTOPE_ERROR_OPTIONS, SearchConfig
from proteoscan.precursors import build_precursor_list
from proteoscan.rescore_features import FEATURE_COLUMNS

C13 = 1.0033548


def _ms2(mz=500.0, z=2):
    return pd.DataFrame([dict(index=1, scan=10, rt_min=5.0, precursor_mz=mz, target_mz=mz, half_width=0.8,
                              charge=z, trailer_mono=mz, master_ms1=None)])


def _features():
    rows = [
        (1, 500.0005, 2, 4e6),              # the selected ion itself
        (2, 500.0 - C13 / 2, 2, 1e6),       # its monoisotope one isotope peak below: still the instrument precursor
        (3, 500.30, 2, 2e6),                # another species, same charge
        (4, 500.0003, 3, 3e6),              # same m/z, other charge
    ]
    return pd.DataFrame([dict(feature_id=i, mz=m, charge=z, rt=5.0, rt_begin=4.8, rt_end=5.2, area=a)
                         for i, m, z, a in rows])


def test_entries_carry_primary_flag_and_area_share():
    e = build_precursor_list(None, _features(), _ms2()).set_index("feature_id")
    assert e.loc[1, "primary"] == 1 and e.loc[2, "primary"] == 1
    assert e.loc[3, "primary"] == 0 and e.loc[4, "primary"] == 0
    assert e.loc[1, "area_frac"] == pytest.approx(0.4)
    assert e.loc[[1, 2, 3, 4], "area_frac"].sum() == pytest.approx(1.0)
    assert (e["n_entries"] == len(e)).all()
    featureless = e[e["source"] == "featureless"]
    assert (featureless["primary"] == 1).all() and (featureless["area_frac"] == 0).all()


def test_feature_linking_primary_drops_coisolated_species():
    allm = build_precursor_list(None, _features(), _ms2(), feature_linking="all")
    prim = build_precursor_list(None, _features(), _ms2(), feature_linking="primary")
    assert set(allm["feature_id"]) >= {1, 2, 3, 4}
    assert set(prim["feature_id"]) - {0} == {1, 2}
    assert list(prim["precursor_idx"]) == list(range(len(prim)))
    with pytest.raises(ValueError):
        build_precursor_list(None, _features(), _ms2(), feature_linking="some")


def test_entry_features_reach_the_rescorer():
    assert {"entry_primary", "entry_area_frac", "entry_n"} <= set(FEATURE_COLUMNS)
    entries = pd.DataFrame(dict(query_id=[100, 101], scan=[10, 10], precursor_idx=[0, 1], feature_id=[7, 0],
                                area=[1e6, np.nan], rt=[5.0, 5.0], mz=[500.01, 500.0], mz_uncal=[500.0, 499.99],
                                source=["feature", "featureless"], primary=[0, 1], area_frac=[1.0, 0.0],
                                n_entries=[2, 2]))
    cands = pd.DataFrame(dict(scan=[100, 101], rt=[5.0, 5.0]))
    c = pl.attach_entries(cands, entries)
    assert list(c["scan"]) == [10, 10] and list(c["query_id"]) == [100, 101]
    assert list(c["entry_primary"]) == [0.0, 1.0] and list(c["entry_n"]) == [2.0, 2.0]
    assert list(c["entry_mz"]) == [500.0, 499.99]      # the measured m/z, not the recalibrated one


def test_config_validation():
    for iso in ISOTOPE_ERROR_OPTIONS:
        SearchConfig(isotope_errors=iso)
    with pytest.raises(ValueError, match="isotope_errors"):
        SearchConfig(isotope_errors=(-2, 2))
    with pytest.raises(ValueError):
        SearchConfig(feature_linking="most")
    with pytest.raises(ValueError):
        SearchConfig(precursor_recalibration="always")
    assert SearchConfig(precursor_recalibration="offset").precursor_recalibration == "offset"
    assert SearchConfig(precursor_recalibration=False).precursor_recalibration == "off"   # YAML `off`
    assert SearchConfig(precursor_recalibration=True).precursor_recalibration == "auto"
    c = SearchConfig()
    assert c.precursor_recalibration == "auto" and c.feature_linking == "all" and c.clip_nterm_methionine


def _cands(n, ppm, evalue=1e-4, seed=0):
    r = np.random.default_rng(seed)
    return pd.DataFrame(dict(rank=1, is_decoy=False, isotope_error=0, comet_evalue=evalue,
                             hyperscore=r.uniform(1, 5, n), precursor_ppm=ppm + r.normal(0, 1.5, n)))


def test_offset_estimate_uses_confident_targets_only():
    good = _cands(1000, 9.0)
    noise = _cands(4000, -20.0, evalue=5.0, seed=1)
    decoys = _cands(500, 30.0, seed=2).assign(is_decoy=True)
    iso = _cands(500, 30.0, seed=3).assign(isotope_error=1)
    est = pl.estimate_precursor_offset(pd.concat([good, noise, decoys, iso], ignore_index=True))
    assert est["offset_ppm"] == pytest.approx(9.0, abs=0.2) and est["n"] == 1000
    assert est["robust_sd_ppm"] == pytest.approx(1.5, abs=0.3)


def test_offset_estimate_needs_enough_psms():
    est = pl.estimate_precursor_offset(_cands(50, 9.0))
    assert est["offset_ppm"] is None and est["n"] < pl.RECAL_MIN_PSMS


class _FakeSearch:
    """Stands in for stage_precursors + Comet: the measured error is true offset minus the applied offset."""

    def __init__(self, true_ppm, n=1000):
        self.true_ppm, self.n, self.calls = true_ppm, n, []

    def install(self, monkeypatch):
        state = {}

        def stage_precursors(mzml, rdir, features_path, cfg, offset_ppm=None):
            state["offset"] = offset_ppm
            return pd.DataFrame(dict(query_id=[1])), rdir / "x.expanded.mzML"

        def run_comet(cfg, mzmls, out, num_output_lines=10):
            self.calls.append((out.name, state["offset"], cfg.parent_mass_error_tolerance_ppm, cfg.digest_mode,
                               num_output_lines))
            return [out / "x.txt"]

        def read_comet(txt, **kw):
            c = _cands(self.n, self.true_ppm - state["offset"])
            return c.assign(proteins="sp|P1|A", scan=1, rt=1.0)

        def attach_entries(c, entries):
            return c.assign(raw_rt=c["rt"])

        for name, fn in (("stage_precursors", stage_precursors), ("run_comet", run_comet),
                         ("read_comet", read_comet), ("attach_entries", attach_entries),
                         ("canonicalize_clipped_nterm", lambda c, fasta: c)):
            monkeypatch.setattr(pl, name, fn)
        return self


def _search(tmp_path, cfg):
    return pl.search_run(tmp_path / "run.mzML", tmp_path, None, cfg, {}, 10, 0.0)


def test_auto_recalibration_shifts_and_narrows_the_search(tmp_path, monkeypatch):
    fake = _FakeSearch(9.0).install(monkeypatch)
    cands, rec = _search(tmp_path, SearchConfig(parent_mass_error_tolerance_ppm=15.0, digest_mode="Semi-Specific"))
    assert rec["source"] == "auto" and rec["offset_ppm"] == pytest.approx(9.0, abs=0.3)
    est, final = fake.calls
    # quick estimation search: widened tolerance, fully specific, one candidate
    assert est == ("comet_estimate", 0.0, pl.RECAL_SEARCH_TOL_PPM, "Specific", 1)
    # real search: offset applied, tolerance 10 robust SDs (1.5 ppm here) but within [5 ppm, configured]
    assert final[0] == "comet_recal" and final[1] == pytest.approx(9.0, abs=0.3)
    assert final[2] == pytest.approx(15.0) and rec["tolerance_ppm"] == pytest.approx(15.0)
    assert final[3:] == ("Semi-Specific", 10)
    assert abs(cands["precursor_ppm"].median()) < 0.3


def test_calibrated_tolerance_bounds():
    assert pl.calibrated_tolerance(15.0, 0.3) == pl.RECAL_MIN_TOL_PPM       # tight instrument: the floor
    assert pl.calibrated_tolerance(15.0, 0.8) == pytest.approx(8.0)          # 10 robust SDs
    assert pl.calibrated_tolerance(15.0, 4.0) == 15.0                        # never wider than configured
    assert pl.calibrated_tolerance(3.0, 0.3) == 3.0
    assert pl.calibrated_tolerance(15.0, None) == 15.0                       # no estimate: as configured


def test_offset_mode_keeps_the_configured_tolerance(tmp_path, monkeypatch):
    fake = _FakeSearch(9.0).install(monkeypatch)
    monkeypatch.setattr(pl, "estimate_precursor_offset", lambda c: dict(offset_ppm=9.0, n=1000, robust_sd_ppm=0.3))
    _, rec = _search(tmp_path, SearchConfig(parent_mass_error_tolerance_ppm=15.0, precursor_recalibration="offset"))
    assert rec["tolerance_ppm"] == 15.0 and fake.calls[1][:3] == ("comet_recal", 9.0, 15.0)
    fake = _FakeSearch(9.0).install(monkeypatch)
    _, rec = _search(tmp_path, SearchConfig(parent_mass_error_tolerance_ppm=15.0))        # auto narrows
    assert rec["tolerance_ppm"] == 5.0 and fake.calls[1][:3] == ("comet_recal", 9.0, 5.0)


def test_auto_recalibration_does_not_shift_small_offsets(tmp_path, monkeypatch):
    fake = _FakeSearch(0.3).install(monkeypatch)
    _, rec = _search(tmp_path, SearchConfig(parent_mass_error_tolerance_ppm=25.0))
    assert rec["offset_ppm"] == 0.0 and rec["estimate"]["offset_ppm"] == pytest.approx(0.3, abs=0.2)
    assert [c[:2] for c in fake.calls] == [("comet_estimate", 0.0), ("comet", 0.0)]
    assert fake.calls[0][2] == 25.0 and fake.calls[1][2] == pytest.approx(15.0, abs=2.0)   # 10 x 1.5 ppm robust SD
    assert rec["tolerance_ppm"] == fake.calls[1][2]


def test_recalibration_off_or_configured_offset_searches_once(tmp_path, monkeypatch):
    fake = _FakeSearch(9.0).install(monkeypatch)
    _, rec = _search(tmp_path, SearchConfig(precursor_recalibration="off"))
    assert rec == dict(offset_ppm=0.0, source="none", tolerance_ppm=15.0)
    assert [c[:3] for c in fake.calls] == [("comet", 0.0, 15.0)]
    fake = _FakeSearch(9.0).install(monkeypatch)
    _, rec = _search(tmp_path, SearchConfig(precursor_mass_offset_ppm=8.5))
    assert rec == dict(offset_ppm=8.5, source="config", tolerance_ppm=15.0)
    assert [c[:2] for c in fake.calls] == [("comet", 8.5)]


def test_too_few_psms_falls_back_to_uncalibrated(tmp_path, monkeypatch):
    fake = _FakeSearch(9.0, n=40).install(monkeypatch)
    _, rec = _search(tmp_path, SearchConfig(parent_mass_error_tolerance_ppm=15.0))
    assert rec["offset_ppm"] == 0.0 and rec["estimate"]["offset_ppm"] is None and rec["tolerance_ppm"] == 15.0
    assert [c[:3] for c in fake.calls] == [("comet_estimate", 0.0, 20.0), ("comet", 0.0, 15.0)]


def test_second_pass_uses_the_first_pass_tolerance(tmp_path):
    cfg = SearchConfig(parent_mass_error_tolerance_ppm=15.0)
    assert pl.run_tolerance(tmp_path, "a.mzML", cfg) == 15.0
    (tmp_path / "precursor_recalibration.json").write_text(json.dumps({"a.mzML": {"tolerance_ppm": 5.0}}))
    assert pl.run_tolerance(tmp_path, "a.mzML", cfg) == 5.0 and pl.run_tolerance(tmp_path, "b.mzML", cfg) == 15.0


def test_clipped_protein_nterm_is_reported_the_same_for_targets_and_decoys(tmp_path):
    from proteoscan.comet_search import canonicalize_clipped_nterm
    fa = tmp_path / "db.fasta"
    fa.write_text(">sp|P1|A_HUMAN x\nMSTPEPTIDEKAAAR\n>sp|P2|B_HUMAN\nAMGGGGK\n>sp|P3|C_HUMAN\nMPLLLK\n")
    c = pd.DataFrame(dict(
        stripped=["STPEPTIDEK", "EDITPEPTSK", "GGGGK", "STPEPT", "PLLLK", "AAAR"],
        prev_aa=["M", "M", "M", "M", "M", "K"], next_aa=["A", "A", "-", "I", "-", "-"],
        proteins=["sp|P1|A_HUMAN", "DECOY_sp|P1|A_HUMAN", "sp|P2|B_HUMAN", "sp|P9|Z_HUMAN;sp|P1|A_HUMAN",
                  "sp|P3|C_HUMAN", "sp|P1|A_HUMAN"],
        semi_enzymatic=[1, 1, 1, 1, 1, 0]))
    out = canonicalize_clipped_nterm(c.copy(), fa)
    # target after Met-1, its decoy (reversed, C-terminal residue kept), internal Met, semi-specific C-terminus,
    # Pro after the clipped Met (not an enzymatic N-terminus), untouched row
    assert list(out["prev_aa"]) == ["-", "-", "M", "-", "-", "K"]
    assert list(out["semi_enzymatic"]) == [0, 0, 1, 1, 1, 0]


def test_rescoring_is_independent_of_row_order():
    from proteoscan.rescore import rescore_tda
    r = np.random.default_rng(3)
    n = 6000
    decoy = r.random(n) < 0.5
    good = (~decoy) & (r.random(n) < 0.5)
    df = pd.DataFrame(dict(source_file="a.mzML", scan=np.arange(n), precursor_idx=0, rank=1, charge=2,
                           peptide_sage=[f"PEP{i}K" for i in range(n)], proteins=np.where(decoy, "#DECOY#P", "P"),
                           is_decoy=decoy, hyperscore=r.normal(1, 1, n) + 3 * good,
                           delta_next=r.normal(0, 1, n) + 2 * good, peptide_len=r.integers(7, 30, n)))
    cols = ["hyperscore", "delta_next", "peptide_len"]
    a = rescore_tda(df, feature_cols=cols).set_index("scan")["score"]
    b = rescore_tda(df.sample(frac=1.0, random_state=1), feature_cols=cols).set_index("scan")["score"]
    assert np.array_equal(a.sort_index().to_numpy(), b.sort_index().to_numpy())


def _fp(tmp_path, **kw):
    mz = tmp_path / "a.mzML"
    if not mz.exists():
        mz.write_text("x")
    return pl.cache_fingerprint(SearchConfig(**kw), [mz], None, 10, True, "ms2pip")


def test_cache_fingerprint_ignores_threads_and_tracks_settings(tmp_path):
    assert _fp(tmp_path, threads=4)["sha1"] == _fp(tmp_path, threads=32)["sha1"]
    assert _fp(tmp_path)["sha1"] != _fp(tmp_path, parent_mass_error_tolerance_ppm=10.0)["sha1"]
    assert _fp(tmp_path)["sha1"] != _fp(tmp_path, feature_linking="primary")["sha1"]
    a = _fp(tmp_path)
    (tmp_path / "a.mzML").write_text("another file of another size")
    assert a["sha1"] != _fp(tmp_path)["sha1"]


def test_check_cache_refuses_stale_stages_and_overwrite_clears_them(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    fp = _fp(tmp_path)
    pl.check_cache(out, fp)                                   # fresh directory: fingerprint written
    assert json.loads((out / "cache_fingerprint.json").read_text())["sha1"] == fp["sha1"]
    rd = out / "a"
    (rd / "comet").mkdir(parents=True)
    for n in ("precursors.parquet", "precursors_raw.parquet", "a.expanded.mzML", "ms1_features.parquet"):
        (rd / n).write_text("x")
    (out / "candidates_features.parquet").write_text("x")
    pl.check_cache(out, fp)                                   # same settings: cache kept
    assert (out / "candidates_features.parquet").exists()
    other = _fp(tmp_path, parent_mass_error_tolerance_ppm=10.0)
    with pytest.raises(RuntimeError, match="parent_mass_error_tolerance_ppm"):
        pl.check_cache(out, other)
    assert (out / "candidates_features.parquet").exists()
    pl.check_cache(out, other, overwrite=True)
    assert not (out / "candidates_features.parquet").exists() and not (rd / "comet").exists()
    assert not (rd / "precursors.parquet").exists() and not (rd / "a.expanded.mzML").exists()
    assert (rd / "ms1_features.parquet").exists()             # same input file: detected features are still valid
    assert json.loads((out / "cache_fingerprint.json").read_text())["sha1"] == other["sha1"]


def test_check_cache_refuses_directories_without_fingerprint(tmp_path):
    out = tmp_path / "old"
    out.mkdir()
    (out / "psms_scored.parquet").write_text("x")
    with pytest.raises(RuntimeError, match="no fingerprint"):
        pl.check_cache(out, _fp(tmp_path))


def test_entrapment_twin_inherits_protein_support():
    acc = pd.DataFrame(dict(stripped=["PEPTIDEK", "PEPTLDEK", "ANOTHERK", "TWINONLYK"],
                            proteins=["sp|P1|A_HUMAN", "sp|P1|A_HUMAN", "sp|P1|A_HUMAN;sp|P2|B_HUMAN",
                                      "ENTRAP_sp|P2|B_HUMAN"]))
    support, seqs = twopass.protein_support(acc, twin_prefix="ENTRAP_")
    assert support["P1"] == 2 and support["P2"] == 1          # I=L collapses; the twin's own hit does not count
    assert support["ENTRAP_P1"] == 2 and support["ENTRAP_P2"] == 1
    assert "ENTRAP_P2" not in seqs
    plain, _ = twopass.protein_support(acc)
    assert "ENTRAP_P1" not in plain
    assert twopass._acc("DECOY_ENTRAP_sp|P1|A_HUMAN", "ENTRAP_") == ("ENTRAP_P1", True)
