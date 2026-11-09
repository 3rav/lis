#!/usr/bin/env python3
"""Research-only automatic candidate vectors from A plus baseline modes.

Deterministic projected gradient filtering on A^T A, without offline
enriched-mode input. These are NOT certified eigenvectors or proven
MFAMG contact near-nullspace corrections. No private data is uploaded.
"""
import argparse
from array import array
import math
import os
import random
import sys

from mfamg_mode_audit import header_and_sizes, read_modes, dot, norm, remove_span
from mfamg_mode_ritz_selection import orthonormalize


def read_sparse(path, n, max_nnz):
    with open(path, encoding="ascii") as stream:
        symmetry, (rows, cols, nnz) = header_and_sizes(stream, "coordinate")
        if rows != n or cols != n or not (0 <= nnz <= max_nnz):
            raise ValueError("matrix dimensions/entries invalid or cap exceeded")
        row_idx, col_idx, vals = array("I"), array("I"), array("d")
        row_abs, col_abs = [0.0] * n, [0.0] * n
        for k in range(nnz):
            line = stream.readline()
            while line and (not line.strip() or line.lstrip().startswith("%")):
                line = stream.readline()
            parts = line.split()
            if len(parts) != 3:
                raise ValueError("malformed matrix entry %d" % (k + 1))
            try:
                i, j = int(parts[0]) - 1, int(parts[1]) - 1
                value = float(parts[2])
            except ValueError as exc:
                raise ValueError("invalid matrix entry %d" % (k + 1)) from exc
            if not (0 <= i < n and 0 <= j < n and math.isfinite(value)):
                raise ValueError("nonfinite/out-of-range matrix entry %d" % (k + 1))
            if symmetry == "skew-symmetric" and i == j and value != 0:
                raise ValueError("nonzero diagonal in skew-symmetric matrix")
            row_idx.append(i)
            col_idx.append(j)
            vals.append(value)
            ab = abs(value)
            row_abs[i] += ab
            col_abs[j] += ab
            if symmetry != "general" and i != j:
                row_idx.append(j)
                col_idx.append(i)
                vals.append(-value if symmetry == "skew-symmetric" else value)
                row_abs[j] += ab
                col_abs[i] += ab
        if any(line.strip() and not line.lstrip().startswith("%") for line in stream):
            raise ValueError("extra matrix entries")
    bound = max(row_abs) * max(col_abs)
    if not math.isfinite(bound) or bound <= 0:
        raise ValueError("nonfinite/zero matrix norm bound")
    return row_idx, col_idx, vals, bound


def multiply(rows, cols, vals, vector, transpose=False):
    result = [0.0] * len(vector)
    if transpose:
        for i, j, val in zip(rows, cols, vals):
            result[j] += val * vector[i]
    else:
        for i, j, val in zip(rows, cols, vals):
            result[i] += val * vector[j]
    if not all(math.isfinite(x) for x in result):
        raise ValueError("nonfinite sparse multiplication")
    return result


def normalized(vec, basis, tolerance):
    projected = remove_span(vec, basis)
    size = norm(projected)
    if size <= tolerance * max(norm(vec), 1e-300):
        raise ValueError("dependent candidate vector")
    return [x / size for x in projected]


def generate(args):
    if os.path.lexists(args.output):
        raise ValueError("output already exists")
    n, base = read_modes(args.baseline)
    if n > 250000 or len(base) + args.count > 256:
        raise ValueError("matrix size or mode count exceeds prototype limits")
    basis = orthonormalize(base, tolerance=args.rank_tol)
    if n - len(basis) < args.count:
        raise ValueError("insufficient free dimensions for candidates")
    rows, cols, vals, bound = read_sparse(args.matrix, n, args.max_nnz)
    step = args.step_scale / bound
    rng = random.Random(args.seed)
    candidates = []
    for j in range(args.count):
        vec = normalized([rng.uniform(-1.0, 1.0) for _ in range(n)],
                         basis, args.rank_tol)
        before = norm(multiply(rows, cols, vals, vec))
        for _ in range(args.steps):
            av = multiply(rows, cols, vals, vec)
            grad = multiply(rows, cols, vals, av, transpose=True)
            vec = normalized([x - step * g for x, g in zip(vec, grad)],
                             basis, args.rank_tol)
        after = norm(multiply(rows, cols, vals, vec))
        if not math.isfinite(after) or after > before * (1 + 1e-10):
            raise ValueError("candidate %d failed A-norm reduction" % (j + 1))
        candidates.append(vec)
        basis.append(vec)
        print("AUTO_MODE_%02d_A_NORM_BEFORE=%.12g" % (j + 1, before))
        print("AUTO_MODE_%02d_A_NORM_AFTER=%.12g" % (j + 1, after))
    with open(args.output, "x", encoding="ascii", newline="\n") as stream:
        stream.write("%%MatrixMarket matrix array real general\n")
        stream.write("% Experimental automatic A^T*A filtered candidate directions\n")
        stream.write("% Not eigencertified; real contact validation required\n")
        stream.write("%d %d\n" % (n, len(base) + len(candidates)))
        for column in base + candidates:
            for value in column:
                stream.write("%.17g\n" % value)
    print("AUTO_MODE_OUTPUT=" + args.output)
    print("AUTO_MODE_EXTRA_COUNT=%d" % len(candidates))
    print("AUTO_MODE_NOTE=automatic_candidates_not_validated_for_contact")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("matrix")
    p.add_argument("baseline")
    p.add_argument("--output", required=True)
    p.add_argument("--count", type=int, default=6)
    p.add_argument("--steps", type=int, default=60)
    p.add_argument("--seed", type=int, default=20261109)
    p.add_argument("--step-scale", type=float, default=0.8)
    p.add_argument("--rank-tol", type=float, default=1e-10)
    p.add_argument("--max-nnz", type=int, default=10000000)
    args = p.parse_args()
    if not (1 <= args.count <= 32 and 1 <= args.steps <= 500 and
            math.isfinite(args.step_scale) and 0 < args.step_scale <= 1 and
            math.isfinite(args.rank_tol) and 0 < args.rank_tol < 1 and
            args.max_nnz > 0):
        p.error("invalid count/steps/filter tolerances")
    try:
        generate(args)
        return 0
    except (OSError, ValueError, OverflowError, MemoryError) as exc:
        print("AUTO_MODE_ERROR=" + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
