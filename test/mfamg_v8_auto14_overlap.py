#!/usr/bin/env python3
"""Read-only, post-hoc geometric audit of generated14 versus offline18 Ritz.

The six offline Ritz directions are used ONLY as a diagnostic control, never
as generator inputs. Subspace overlap is neither a convergence predictor nor
proof of why the 14-mode fixed-A solve failed. No A/b/EXE inputs, file writes,
paths, vectors or input digests appear in the report.
"""
import argparse
import hashlib
import math
from pathlib import Path
import sys

from mfamg_mode_audit import dot, read_modes
from mfamg_mode_ritz_selection import orthonormalize


def sha256_input(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for chunk in iter(lambda: source.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.digest()


def analyze(base_path, generated_path, control_path):
    paths = (base_path, generated_path, control_path)
    before = tuple(sha256_input(p) for p in paths)
    n0, base = read_modes(base_path)
    ng, generated = read_modes(generated_path)
    no, control = read_modes(control_path)
    if n0 != ng or n0 != no or (len(base), len(generated), len(control)) != (12, 14, 18):
        raise ValueError("expected common row count and 12/14/18 mode columns")
    if generated[:12] != base or control[:12] != base:
        raise ValueError("baseline prefix differs")
    qbase = orthonormalize(base, tolerance=1e-10)
    qgen = orthonormalize(generated[12:], prior=qbase, tolerance=1e-10)
    qoffline = orthonormalize(control[12:], prior=qbase, tolerance=1e-10)
    # The 2x6 overlap matrix measures the geometry of *extra* subspaces
    # AFTER projecting out the same 12-dimensional baseline space.
    overlap = [[dot(g, o) for o in qoffline] for g in qgen]
    a = math.fsum(x * x for x in overlap[0])
    d = math.fsum(x * x for x in overlap[1])
    b = math.fsum(x * y for x, y in zip(overlap[0], overlap[1]))
    halfgap = 0.5 * math.hypot(a - d, 2 * b)
    mean = 0.5 * (a + d)
    eigenvalues = (mean - halfgap, mean + halfgap)
    if not all(math.isfinite(x) and -1e-8 <= x <= 1 + 1e-8 for x in eigenvalues):
        raise ValueError("overlap Gram matrix is numerically inconsistent")
    principal_cosines = [math.sqrt(min(1.0, max(0.0, x))) for x in eigenvalues]
    if tuple(sha256_input(p) for p in paths) != before:
        raise ValueError("mode inputs changed during post-hoc overlap audit")
    return {
        "n": n0,
        "mean_squared_projection": min(1.0, max(0.0, mean)),
        "principal_cosine_min": principal_cosines[0],
        "principal_cosine_max": principal_cosines[1],
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--baseline", required=True)
    p.add_argument("--generated", required=True)
    p.add_argument("--offline-control", required=True)
    args = p.parse_args()
    try:
        summary = analyze(args.baseline, args.generated, args.offline_control)
    except (ValueError, OSError, OverflowError, MemoryError) as exc:
        print("V8_AUTO14_OVERLAP=BLOCKED", file=sys.stderr)
        print("V8_AUTO14_OVERLAP_REASON=" + type(exc).__name__, file=sys.stderr)
        return 2
    print("V8_AUTO14_OVERLAP=OK")
    print("V8_AUTO14_OVERLAP_BASE_MODES=12")
    print("V8_AUTO14_OVERLAP_GENERATED_EXTRA=2")
    print("V8_AUTO14_OVERLAP_OFFLINE_EXTRA=6")
    print("V8_AUTO14_OVERLAP_MEAN_SQUARED_PROJECTION=%.12g" % summary["mean_squared_projection"])
    print("V8_AUTO14_OVERLAP_PRINCIPAL_COSINE_MIN=%.12g" % summary["principal_cosine_min"])
    print("V8_AUTO14_OVERLAP_PRINCIPAL_COSINE_MAX=%.12g" % summary["principal_cosine_max"])
    print("V8_AUTO14_OVERLAP_INTERPRETATION=GEOMETRY_ONLY_NOT_CAUSAL_OR_CONVERGENCE_PROOF")
    print("V8_AUTO14_OVERLAP_PRIVATE_INPUTS=READ_ONLY")
    return 0


if __name__ == "__main__":
    sys.exit(main())
