"""``casanovo`` CLI with fast result writing: ``python -m proteoscan.casanovo_fastio sequence ...``

Casanovo 5.2.1 parses every predicted peptide with spectrum_utils' ProForma (lark) parser after
prediction, for the run report and the mzTab ``sequence`` / ``modifications`` columns: ~3.5 ms
per PSM, i.e. ~50 min for 165k spectra x 5 beams (as long as the GPU prediction itself). Here the
plain sequence is obtained with a regex and the mzTab ``modifications`` column is written as
``null``; the ProForma column (``opt_global_cv_MS:1003169_proforma_peptidoform_sequence``),
scores and per-residue scores, which are what proteoscan reads, are unchanged.
"""
import re
import sys

_MOD_RE = re.compile(r"\[[^\]]*\]-?")


def _patch():
    from casanovo.data import psm as _psm
    _psm.PepSpecMatch.aa_sequence = property(lambda self: _MOD_RE.sub("", self.sequence))
    _psm.PepSpecMatch.modifications = property(lambda self: "null")


def main():
    _patch()
    from casanovo.casanovo import main as cli
    sys.argv[0] = "casanovo"
    cli()


if __name__ == "__main__":
    main()
