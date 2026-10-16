"""Command-line entry point: python -m proteoscan.cli <command> ..."""
from __future__ import annotations

import argparse
from pathlib import Path


def cmd_convert(a):
    from .convert import raw_to_mzml
    for r in a.raw:
        print(raw_to_mzml(r, a.out))


def cmd_search(a):
    from .config import SearchConfig
    from .sage_search import run_sage
    cfg = SearchConfig.from_yaml(a.config)
    if a.threads:
        cfg.threads = a.threads
    if a.chunks > 1:
        from .sage_search import run_sage_chunked
        print(run_sage_chunked(cfg, a.mzml, a.out, n_chunks=a.chunks, chimera=a.chimera, report_psms=a.report_psms))
    else:
        print(run_sage(cfg, a.mzml, a.out, chimera=a.chimera, report_psms=a.report_psms))


def cmd_concordance(a):
    from .concordance import main as conc_main
    conc_main(["--ref", a.ref, "--test", a.test, "--out", a.out])


def cmd_ui(a, rest):
    from .ui import main as ui_main
    ui_main(rest)


def cmd_summary(a, rest):
    from .summary import main as summary_main
    return summary_main(rest)


def cmd_run(a, rest):
    from .pipeline import main as pipeline_main
    pipeline_main(rest)


def main(argv=None):
    p = argparse.ArgumentParser(prog="proteoscan")
    sub = p.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("convert", help="Thermo raw -> mzML"); c.add_argument("raw", nargs="+"); c.add_argument("--out", required=True); c.set_defaults(fn=cmd_convert)
    s = sub.add_parser("search", help="Sage search from a PEAKS-style config"); s.add_argument("mzml", nargs="+"); s.add_argument("--config", required=True); s.add_argument("--out", required=True); s.add_argument("--threads", type=int); s.add_argument("--chimera", action="store_true"); s.add_argument("--report-psms", type=int, default=1); s.add_argument("--chunks", type=int, default=1, help="split the FASTA into N chunks to bound memory"); s.set_defaults(fn=cmd_search)
    k = sub.add_parser("concordance", help="compare an export dir to a PEAKS reference"); k.add_argument("--ref", required=True); k.add_argument("--test", required=True); k.add_argument("--out", required=True); k.set_defaults(fn=cmd_concordance)
    r = sub.add_parser("run", help="end-to-end pipeline, mzML/raw -> PEAKS-format export (options: proteoscan run -h)", add_help=False); r.set_defaults(fn=None)
    u = sub.add_parser("ui", help="self-contained HTML results viewer (options: python -m proteoscan.ui -h)", add_help=False); u.set_defaults(fn=None)
    m = sub.add_parser("summary", help="PEAKS-style HTML export pages: db.summary/proteins/peptides.html + img/ (options: python -m proteoscan.summary -h)", add_help=False); m.set_defaults(fn=None)
    a, rest = p.parse_known_args(argv)
    if a.cmd == "run":
        return cmd_run(a, rest)
    if a.cmd == "ui":
        return cmd_ui(a, rest)
    if a.cmd == "summary":
        return cmd_summary(a, rest)
    if rest:
        p.error("unrecognized arguments: " + " ".join(rest))
    a.fn(a)


if __name__ == "__main__":
    main()
