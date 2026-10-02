#!/usr/bin/env python
"""Reproduce PEAKS protein inference from PEAKS's own peptides and diff against its export.

For each reference sample, feed db.peptides.csv (+ db.psms.csv for #Spec and features,
+ <sample>.db.features.csv for unrounded feature areas, + full-precision PSM scores from
db.mzidentml_1_2_0.xml) into ``proteoscan.proteins.infer_proteins`` and compare with
db.proteins.csv and db.protein-peptides.csv.

Usage::

    python scripts/eval_proteins.py [--ref ROOT] [--fasta FASTA] [--samples 1 2 3 4]
        [--mapping fasta|accession|both] [--no-mzid] [--out DIR]

Every comparison is an exact string match of the values as PEAKS would print them
(``proteoscan.peaks_export.format_value``), except "-10LgP (2 dp)" which compares the
protein score rounded to 2 decimals.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from proteoscan.peaks_export import format_value, read_export  # noqa: E402
from proteoscan.proteins import infer_proteins, read_mzid_psm_scores, read_protein_db  # noqa: E402

REF = "/drive_a/drew/vendor-exports/DrewR_Adam_91526_DDA_OT"
FASTA = "/drive_a/drew/proteoscan-work/proteins/human_sp_peaks_merged.fasta"
LOCAL_FASTA = "/drive_a/drew/openms-data/260921-DrewR_Adam_91526/fasta/human_sp.fasta"
OUT = "/drive_a/drew/proteoscan-work/proteins/eval"


def build_merged_fasta(ref_root, local_fasta, out_path):
    """Local SwissProt FASTA with entries replaced/added from PEAKS's own db.proteins.fasta
    files (PEAKS's DB is an older UniProt release: renamed entries such as GRP75 -> HSPA9,
    changed descriptions/genes, 8 changed sequences, 2 entries missing locally)."""
    from proteoscan.proteins import _read_fasta

    peaks = {}
    for d in sorted(os.listdir(ref_root)):
        f = os.path.join(ref_root, d, "db.proteins.fasta")
        if os.path.exists(f):
            for h, s in _read_fasta(f):
                peaks[h[1:].split("|")[1]] = (h, s)
    seen = set()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as out:
        for h, s in _read_fasta(local_fasta):
            a = h[1:].split("|")[1]
            if a in peaks:
                h, s = peaks[a]
            seen.add(a)
            out.write(f"{h}\n{s}\n")
        added = [a for a in peaks if a not in seen]
        for a in added:
            h, s = peaks[a]
            out.write(f"{h}\n{s}\n")
    return added


def _fmt(kind, v):
    if v is None or (isinstance(v, float) and np.isnan(v)) or v is pd.NA:
        return ""
    if isinstance(v, str):
        return v
    return format_value(kind, v)


PROTEIN_FIELDS = [
    ("Protein Group", "int"),
    ("Top", "bool"),
    ("-10LgP", "jfloat"),
    ("-10LgP (2 dp)", "fixed2"),
    ("#Peptides", "int"),
    ("#Unique", "int"),
    ("Coverage(%)", "fixed2"),
    ("#Spec", "jfloat"),
    ("Area", "sci"),
    ("Average Mass", "jfloat"),
    ("PTM", "str"),
    ("Gene", "str"),
    ("Description", "str"),
]
PP_FIELDS = [
    ("Protein Group", "int"), ("Top", "bool"), ("Unique", "bool"), ("Extended Peptide", "str"),
    ("-10LgP", "fixed2"), ("Mass", "fixed4"), ("Start", "int"), ("End", "int"),
    ("Area", "sci"), ("#Feature", "int"), ("PTM", "str"), ("AScore", "str"), ("Scan", "int"),
]


def _col(df, name, sample):
    if name in ("#Spec", "Area"):
        return f"{name} {sample}"
    if name == "-10LgP (2 dp)":
        return "-10LgP"
    return name


def partition_agreement(ref, test):
    """Groups as sets of accessions: exact set matches."""
    rg = {frozenset(g["Accession"]) for _, g in ref.groupby("Protein Group")}
    tg = {frozenset(g["Accession"]) for _, g in test.groupby("Protein Group")}
    return len(rg), len(tg), len(rg & tg), sorted(sorted(x) for x in rg - tg), sorted(sorted(x) for x in tg - rg)


def order_agreement(ref, test):
    """Protein Group numbering: fraction of pairs of shared groups in the same order."""
    r = ref[ref.Top].groupby("Protein Group").Accession.apply(frozenset)
    t = test[test.Top].groupby("Protein Group").Accession.apply(frozenset)
    rmap = {v: k for k, v in r.items()}
    tmap = {v: k for k, v in t.items()}
    shared = [s for s in rmap if s in tmap]
    rv = np.array([rmap[s] for s in shared])
    tv = np.array([tmap[s] for s in shared])
    if len(shared) < 2:
        return 1.0, len(shared)
    o = np.argsort(rv)
    tv = tv[o]
    conc = disc = 0
    for i in range(len(tv)):
        conc += int((tv[i + 1:] > tv[i]).sum())
        disc += int((tv[i + 1:] < tv[i]).sum())
    return conc / max(conc + disc, 1), len(shared)


def numbering_key_check(ref_pr, ref_pp):
    """PEAKS numbers groups by the number of peptide rows of the group's first protein,
    descending.  Returns (#adjacent violations, #tied adjacent pairs) over the reference
    Top proteins (a violation would falsify the rule; ties are unresolved)."""
    nrows = ref_pp.groupby("Accession").size()
    top = ref_pr[ref_pr.Top].groupby("Protein Group").Accession.first().sort_index()
    v = nrows.reindex(top.values).to_numpy()
    d = np.diff(v)
    return int((d > 0).sum()), int((d == 0).sum())


def compare(ref_pr, ref_pp, pr, pp, sample):
    res = {}
    viol, ties = numbering_key_check(ref_pr, ref_pp)
    res["numbering_rule"] = {"violations": viol, "adjacent_ties": ties}
    n_ref, n_test, n_same, ref_only_g, test_only_g = partition_agreement(ref_pr, pr)
    res["groups"] = {"ref": n_ref, "test": n_test, "identical_membership": n_same,
                     "rate": n_same / max(n_ref, 1), "ref_only": ref_only_g[:20],
                     "test_only": test_only_g[:20]}
    conc, nshared = order_agreement(ref_pr, pr)
    res["group_numbering_pair_order"] = conc
    ra, ta = set(ref_pr.Accession), set(pr.Accession)
    res["proteins"] = {"ref": len(ra), "test": len(ta), "shared": len(ra & ta),
                       "ref_only": sorted(ra - ta), "test_only": sorted(ta - ra)}
    m = ref_pr.merge(pr, on="Accession", suffixes=("_ref", "_test"))
    fields = {}
    mism = {}
    for name, kind in PROTEIN_FIELDS:
        c = _col(ref_pr, name, sample)
        a = [_fmt(kind, v) for v in m[c + "_ref"]]
        b = [_fmt(kind, v) for v in m[c + "_test"]]
        ok = [x == y for x, y in zip(a, b)]
        fields[name] = sum(ok) / max(len(ok), 1)
        bad = [(acc, x, y) for acc, x, y, o in zip(m.Accession, a, b, ok) if not o]
        if bad:
            mism[name] = bad[:15]
    res["protein_fields"] = fields
    res["protein_mismatches"] = mism
    key = ["Accession", "Peptide"]
    rk = set(map(tuple, ref_pp[key].values))
    tk = set(map(tuple, pp[key].values))
    res["protein_peptides"] = {"ref": len(rk), "test": len(tk), "shared": len(rk & tk),
                               "ref_only": sorted(rk - tk)[:15], "test_only": sorted(tk - rk)[:15]}
    mm = ref_pp.merge(pp, on=key, suffixes=("_ref", "_test"))
    f2, mis2 = {}, {}
    for name, kind in PP_FIELDS:
        c = f"Area {sample}" if name == "Area" else name
        a = [_fmt(kind, v) for v in mm[c + "_ref"]]
        b = [_fmt(kind, v) for v in mm[c + "_test"]]
        ok = [x == y for x, y in zip(a, b)]
        f2[name] = sum(ok) / max(len(ok), 1)
        bad = [(acc, p, x, y) for acc, p, x, y, o in zip(mm.Accession, mm.Peptide, a, b, ok) if not o]
        if bad:
            mis2[name] = bad[:10]
    res["protein_peptide_fields"] = f2
    res["protein_peptide_mismatches"] = mis2
    return res


def run_sample(ref_root, k, fasta, mapping, use_mzid):
    d = os.path.join(ref_root, f"Analysis 1 - Sample {k}")
    exp = read_export(d, tables=["psms", "peptides", "proteins", "protein_peptides", "features"],
                      verify=False)
    sample = exp.sample_name or f"Sample {k}"
    pep = exp.peptides.copy()
    if use_mzid and os.path.exists(os.path.join(d, "db.mzidentml_1_2_0.xml")):
        mz = read_mzid_psm_scores(os.path.join(d, "db.mzidentml_1_2_0.xml"))
        best = mz.groupby("Peptide").score.max()
        pep["score"] = pep["Peptide"].map(best)
        miss = pep["score"].isna().sum()
        pep["score"] = pep["score"].fillna(pep["-10LgP"])
        # sanity: the rounded mzid score equals the CSV value
        bad = (np.round(pep["score"] + 1e-9, 2) - pep["-10LgP"]).abs() > 0.001
        info = {"mzid_missing": int(miss), "mzid_disagree": int(bad.sum())}
    else:
        info = {}
    pr, pp = infer_proteins(pep, fasta, psms_df=exp.psms, features_df=exp.features,
                            sample=sample, mapping=mapping)
    res = compare(exp.proteins, exp.protein_peptides, pr, pp, sample)
    res["input"] = info
    return res, pr, pp


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", default=REF)
    ap.add_argument("--fasta", default=FASTA)
    ap.add_argument("--local-fasta", default=LOCAL_FASTA)
    ap.add_argument("--samples", type=int, nargs="+", default=[1, 2, 3, 4])
    ap.add_argument("--mapping", default="both", choices=["fasta", "accession", "both"])
    ap.add_argument("--no-mzid", action="store_true")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args(argv)
    if not os.path.exists(a.fasta):
        added = build_merged_fasta(a.ref, a.local_fasta, a.fasta)
        print(f"built {a.fasta} (+{len(added)} entries)")
    db = read_protein_db(a.fasta)
    mappings = ["fasta", "accession"] if a.mapping == "both" else [a.mapping]
    os.makedirs(a.out, exist_ok=True)
    allres = {}
    for mapping in mappings:
        for k in a.samples:
            res, pr, pp = run_sample(a.ref, k, db, mapping, not a.no_mzid)
            allres[f"{mapping}/sample{k}"] = res
            tag = f"{mapping}_sample{k}"
            pr.to_csv(os.path.join(a.out, f"{tag}.proteins.csv"), index=False)
            pp.to_csv(os.path.join(a.out, f"{tag}.protein-peptides.csv"), index=False)
            g = res["groups"]
            print(f"\n== mapping={mapping} sample {k}: groups {g['identical_membership']}/{g['ref']} identical "
                  f"(test {g['test']}), proteins ref {res['proteins']['ref']} test {res['proteins']['test']}, "
                  f"numbering pair-order {res['group_numbering_pair_order']:.4f} "
                  f"(rule violations {res['numbering_rule']['violations']}, ties {res['numbering_rule']['adjacent_ties']}) {res['input']}")
            print("   proteins: " + ", ".join(f"{n} {v:.4f}" for n, v in res["protein_fields"].items()))
            ppr = res["protein_peptides"]
            print(f"   protein-peptides rows ref {ppr['ref']} test {ppr['test']} shared {ppr['shared']}: "
                  + ", ".join(f"{n} {v:.4f}" for n, v in res["protein_peptide_fields"].items()))
            for n, v in res["protein_mismatches"].items():
                print(f"   x {n}: {v[:4]}")
            if res["groups"]["ref_only"]:
                print(f"   ref-only groups: {res['groups']['ref_only'][:4]}")
                print(f"   test-only groups: {res['groups']['test_only'][:4]}")
    with open(os.path.join(a.out, "eval_proteins.json"), "w") as fh:
        json.dump(allres, fh, indent=1, default=str)
    print(f"\nwrote {a.out}/eval_proteins.json")


if __name__ == "__main__":
    main()
