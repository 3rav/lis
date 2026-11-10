#!/usr/bin/env python3
"""Read-only structural preflight for USER-DECLARED consecutive Newton snapshots.

The CSV is authoritative only for the ORDER ITS AUTHOR DECLARED.
No filenames, timestamps, A01/A02 suffixes, or file creation dates prove
physical Newton chronology. This tool cannot certify provenance.
It reads no private data into GitHub or an output directory, writes nothing,
and prints aggregate classifications only. Use before the real triplet runner.
"""
import argparse
import math
import os
from pathlib import Path
import subprocess
import sys

from mfamg_mode_audit import read_modes
from mfamg_v8_snapshot_sequence import read_steps, inspect_preflight


def mm_shape(path, kind):
    """Validate declared Matrix Market A or RHS shape from its small header."""
    with path.open("r", encoding="ascii") as f:
        header = f.readline().strip().lower().split()
        if len(header) != 5 or header[:2] != ["%%matrixmarket", "matrix"]:
            raise ValueError("not a supported Matrix Market matrix header")
        if kind == "rhs":
            if header[2:] != ["array", "real", "general"]:
                raise ValueError("RHS must be Matrix Market real general array")
        elif header[2] != "coordinate" or header[3] not in ("real", "integer"):
            raise ValueError("A must be a real/integer coordinate Matrix Market matrix")
        for line in f:
            if line.strip() and not line.lstrip().startswith("%"):
                try:
                    dims = [int(x) for x in line.split()]
                except ValueError as exc:
                    raise ValueError("invalid Matrix Market dimensions") from exc
                expected_count = 3 if kind == "matrix" else 2
                if len(dims) != expected_count or any(v <= 0 for v in dims):
                    raise ValueError("invalid Matrix Market shape")
                return dims[:2]
    raise ValueError("missing Matrix Market shape")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe", type=Path, required=True)
    p.add_argument("--modes", type=Path, required=True)
    p.add_argument("--snapshots", type=Path, required=True,
                   help="Explicit curated CSV: step,matrix,rhs; chronology USER-declared")
    p.add_argument("--max-windows", type=int, default=1,
                   help="1=first three rows; or 0=all declared transitions")
    args = p.parse_args()
    if args.max_windows < 0 or args.max_windows > 48:
        p.error("--max-windows must be 0..48")
    try:
        steps = read_steps(args.snapshots.resolve())
        if len(steps) < 3:
            raise ValueError("need at least three user-declared consecutive rows")
        paths = [args.exe.resolve(), args.modes.resolve(), args.snapshots.resolve()]
        for row in steps:
            paths.extend([row["matrix"], row["rhs"]])
        if len(set(paths)) != len(paths):
            raise ValueError("duplicated input file or EXE path")
        if any(not x.is_file() or not x.stat().st_size for x in paths):
            raise ValueError("missing/empty declared input or solver")
        n, modes = read_modes(args.modes)
        if n <= 0 or len(modes) < 1 or len(modes) > 256:
            raise ValueError("invalid mode array")
        for row in steps:
            r, c = mm_shape(row["matrix"], "matrix")
            v, cols = mm_shape(row["rhs"], "rhs")
            if r != n or c != n or v != n or cols != 1:
                raise ValueError("matrix/RHS/mode shape mismatch for one declared step")

        windows = len(steps) - 2 if args.max_windows == 0 else min(
            len(steps) - 2, args.max_windows)
        tested_pairs = windows + 1
        base_env = {key: value for key, value in os.environ.items()
                    if not key.startswith("LIS_")}
        states = []
        for i in range(tested_pairs):
            left, right = steps[i], steps[i + 1]
            p_run = subprocess.run(
                [str(args.exe.resolve()), "--compare-csr",
                 str(left["matrix"]), str(right["matrix"])],
                cwd=str(left["matrix"].parent), env=base_env,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, errors="replace", timeout=300, check=False)
            graph, changed, n_detected = inspect_preflight(
                p_run.returncode, p_run.stdout)
            if n_detected is not None and n_detected != n:
                raise ValueError("CSR preflight matrix size differs from modes")
            states.append((graph, changed))
        changed = sum(v[0] != "unchanged" for v in states)
        if not changed:
            raise ValueError("no actual operator changes in this declared sequence")
        print("V8_NEWTON_PREFLIGHT=USER_DECLARED_SEQUENCE_READY")
        print("V8_NEWTON_PREFLIGHT_STEPS=%d" % len(steps))
        print("V8_NEWTON_PREFLIGHT_MODES=%d" % len(modes))
        print("V8_NEWTON_PREFLIGHT_WINDOWS=%d" % windows)
        print("V8_NEWTON_PREFLIGHT_GRAPHS=" +
              ",".join(graph for graph, count in states))
        print("V8_NEWTON_CHRONOLOGY=CSV_ORDER_DECLARED_NOT_INDEPENDENTLY_PROVEN")
        print("V8_NEWTON_DATA_PRIVACY=LOCAL_READ_ONLY")
        return 0
    except (ValueError, OSError, subprocess.TimeoutExpired, UnicodeError) as exc:
        # Deliberately suppress private filenames and raw importer output.
        print("V8_NEWTON_PREFLIGHT=BLOCKED", file=sys.stderr)
        print("V8_NEWTON_PREFLIGHT_REASON=" + type(exc).__name__, file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
