#!/usr/bin/env python3
"""Read-only, standard-library Matrix Market snapshot checks for FE contact.

This is an offline operator audit, not a solver or a CalculiX exporter.
"""
import argparse
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path


def read_matrix_market(path):
    """Return (nrow, ncol, sparse entries) with 0-based coordinates.

    Supports real/integer/pattern coordinate matrices, real/integer general
    arrays and real symmetric coordinate matrices. Duplicate coordinates
    are summed, as required by the Matrix Market coordinate convention.
    """
    with open(path, "r", encoding="utf-8") as stream:
        header = stream.readline().strip().split()
        if len(header) != 5 or [v.lower() for v in header[:2]] != [
            "%%matrixmarket", "matrix"
        ]:
            raise ValueError("invalid Matrix Market matrix header")
        layout, field, symmetry = [v.lower() for v in header[2:]]
        if layout not in ("coordinate", "array"):
            raise ValueError("unsupported Matrix Market layout")
        if field not in ("real", "integer", "pattern"):
            raise ValueError("only real/integer/pattern matrices are supported")
        if symmetry not in ("general", "symmetric"):
            raise ValueError("only general and real symmetric matrices are supported")
        if layout == "array" and (field == "pattern" or symmetry != "general"):
            raise ValueError("only general real/integer arrays are supported")

        def records():
            for line in stream:
                value = line.strip()
                if value and not value.startswith("%"):
                    yield value.split()

        lines = records()
        try:
            dimensions = next(lines)
        except StopIteration as exc:
            raise ValueError("missing Matrix Market dimensions") from exc
        if len(dimensions) != (3 if layout == "coordinate" else 2):
            raise ValueError("invalid Matrix Market dimensions")
        try:
            dims = [int(x) for x in dimensions]
        except ValueError as exc:
            raise ValueError("non-integer matrix dimensions") from exc
        rows, cols = dims[:2]
        if rows <= 0 or cols <= 0:
            raise ValueError("matrix dimensions must be positive")
        if symmetry == "symmetric" and rows != cols:
            raise ValueError("symmetric Matrix Market matrix must be square")

        data = defaultdict(float)
        expected = dims[2] if layout == "coordinate" else rows * cols
        if expected < 0:
            raise ValueError("negative entry count")
        count = 0
        for tokens in lines:
            if count >= expected:
                raise ValueError("extra Matrix Market data")
            if layout == "coordinate":
                if len(tokens) != (2 if field == "pattern" else 3):
                    raise ValueError("invalid coordinate entry")
                try:
                    i, j = int(tokens[0]) - 1, int(tokens[1]) - 1
                    v = 1.0 if field == "pattern" else float(tokens[2])
                except ValueError as exc:
                    raise ValueError("invalid coordinate value") from exc
                if not (0 <= i < rows and 0 <= j < cols):
                    raise ValueError("coordinate index outside matrix")
            else:
                if len(tokens) != 1:
                    raise ValueError("invalid array entry")
                i, j = count % rows, count // rows
                try:
                    v = float(tokens[0])
                except ValueError as exc:
                    raise ValueError("invalid array value") from exc
            if not math.isfinite(v):
                raise ValueError("non-finite matrix coefficient")
            data[(i, j)] += v
            if symmetry == "symmetric" and i != j:
                data[(j, i)] += v
            count += 1
        if count != expected:
            raise ValueError("incomplete Matrix Market file")
        result = dict(data)
        if not all(math.isfinite(v) for v in result.values()):
            raise ValueError("non-finite accumulated matrix coefficient")
        return rows, cols, result


def read_vector(path, n):
    rows, cols, entries = read_matrix_market(path)
    if rows != n or cols != 1:
        raise ValueError("vector dimensions do not match the operator")
    result = [0.0] * n
    for (i, j), value in entries.items():
        result[i] += value
    return result


def matvec(n, matrix, x):
    result = [0.0] * n
    for (i, j), value in matrix.items():
        result[i] += value * x[j]
    return result


def norm(x):
    return math.sqrt(math.fsum(v * v for v in x))


def relative_residual(n, matrix, b, x):
    ax = matvec(n, matrix, x)
    denominator = norm(b)
    return norm([bi - ai for bi, ai in zip(b, ax)]) / (
        denominator if denominator > 0.0 else 1.0
    )


def read_dof_map(path, n):
    """Validate complete 0-based equation,node,dof,body CSV metadata."""
    seen_equations = set()
    seen_node_dofs = set()
    node_dofs = defaultdict(set)
    with open(path, newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or not {"eq", "node", "dof", "body"}.issubset(
            set(reader.fieldnames)
        ):
            raise ValueError("DOF CSV requires eq,node,dof,body columns")
        for entry in reader:
            try:
                eq = int(entry["eq"])
                node = int(entry["node"])
                dof = int(entry["dof"])
            except (ValueError, TypeError) as exc:
                raise ValueError("invalid DOF map integer") from exc
            body = entry["body"].strip()
            if not (0 <= eq < n) or dof not in (1, 2, 3) or node < 0 or not body:
                raise ValueError("invalid DOF map entry")
            if eq in seen_equations or (body, node, dof) in seen_node_dofs:
                raise ValueError("duplicate equation or physical DOF")
            seen_equations.add(eq)
            seen_node_dofs.add((body, node, dof))
            node_dofs[(body, node)].add(dof)
    if len(seen_equations) != n:
        raise ValueError("DOF map must include every active equation exactly once")
    return {
        "equations": len(seen_equations),
        "physical_nodes": len(node_dofs),
        "nodes_with_incomplete_dofs": sum(len(v) < 3 for v in node_dofs.values()),
        "structural_bodies": len({body for body, node in node_dofs}),
    }


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def audit(matrix_path, rhs_path=None, initial_path=None, solution_path=None,
          modes_path=None, dof_map_path=None):
    n, cols, matrix = read_matrix_market(matrix_path)
    if n != cols:
        raise ValueError("contact stiffness matrix must be square")

    row_norm = [0.0] * n
    diag = [0.0] * n
    numerator = 0.0
    denominator = 0.0
    for (i, j), value in matrix.items():
        row_norm[i] += abs(value)
        denominator += value * value
        if i == j:
            diag[i] = value
        difference = value - matrix.get((j, i), 0.0)
        numerator += difference * difference
    if not math.isfinite(numerator) or not math.isfinite(denominator):
        raise ValueError("matrix norms overflowed")
    inf_norm = max(row_norm)
    nonzero_diag = [abs(v) for v in diag if v != 0.0]
    result = {
        "schema_version": 1,
        "operator": {
            "n": n,
            "nnz": len(matrix),
            "file_sha256": digest(matrix_path),
            "frobenius_relative_asymmetry": (
                math.sqrt(numerator / denominator) if denominator else 0.0
            ),
            "infinity_norm": inf_norm,
            "zero_diagonal_count": n - len(nonzero_diag),
            "min_nonzero_abs_diagonal": min(nonzero_diag) if nonzero_diag else None,
            "max_abs_diagonal": max(abs(v) for v in diag),
        },
    }
    if rhs_path is None and (initial_path or solution_path):
        raise ValueError("--initial and --solution require --rhs")
    if rhs_path:
        b = read_vector(rhs_path, n)
        result["rhs"] = {"norm2": norm(b), "file_sha256": digest(rhs_path)}
        x0 = read_vector(initial_path, n) if initial_path else [0.0] * n
        result["initial_relative_residual"] = relative_residual(n, matrix, b, x0)
        if solution_path:
            x = read_vector(solution_path, n)
            result["solution_relative_residual"] = relative_residual(n, matrix, b, x)
    if modes_path:
        mode_rows, num_modes, modes = read_matrix_market(modes_path)
        if mode_rows != n:
            raise ValueError("near-nullspace mode rows do not match matrix")
        reports = []
        for col in range(num_modes):
            mode = [modes.get((i, col), 0.0) for i in range(n)]
            energy = math.fsum(u * v for u, v in zip(mode, matvec(n, matrix, mode)))
            norm2 = math.fsum(v * v for v in mode)
            avnorm = norm(matvec(n, matrix, mode))
            vnorm = math.sqrt(norm2)
            reports.append({
                "index": col,
                "rayleigh_quotient": energy / norm2 if norm2 else None,
                "relative_operator_action": (
                    avnorm / (inf_norm * vnorm)
                    if inf_norm > 0.0 and vnorm > 0.0 else None
                ),
            })
        result["near_nullspace_modes"] = reports
    if dof_map_path:
        result["dof_map"] = read_dof_map(dof_map_path, n)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--rhs", type=Path)
    parser.add_argument("--initial", type=Path)
    parser.add_argument("--solution", type=Path)
    parser.add_argument("--modes", type=Path)
    parser.add_argument("--dof-map", type=Path)
    args = parser.parse_args(argv)
    try:
        report = audit(args.matrix, args.rhs, args.initial, args.solution,
                       args.modes, args.dof_map)
    except (ValueError, OSError, OverflowError, ZeroDivisionError) as exc:
        print("contact-matrix-audit: " + str(exc), file=sys.stderr)
        return 2
    json.dump(report, sys.stdout, sort_keys=True, indent=2, allow_nan=False)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
