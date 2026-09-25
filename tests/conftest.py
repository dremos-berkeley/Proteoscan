"""Shared fixtures.  Reference-dependent tests skip when the PEAKS exports are absent.

Set PEAKSPARITY_REF_ROOT to point at another copy of the DrewR_Adam_91526_DDA_OT exports.
"""
import os
import sys

# Keep BLAS thread pools small on the shared server (must precede numpy import).
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "4")

import pytest  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

REF_ROOT = os.environ.get(
    "PEAKSPARITY_REF_ROOT", "/drive_a/drew/vendor-exports/DrewR_Adam_91526_DDA_OT"
)
SAMPLES = (1, 2, 3, 4)


def ref_dir(k):
    return os.path.join(REF_ROOT, f"Analysis 1 - Sample {k}")


def require_ref(k):
    d = ref_dir(k)
    if not os.path.exists(os.path.join(d, "db.psms.csv")):
        pytest.skip(f"reference export not available: {d}")
    return d
