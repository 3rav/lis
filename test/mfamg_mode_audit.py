#!/usr/bin/env python3
"""Offline MFAMG basis audit, using only Python's standard library.

Read-only: streams the sparse matrix once; never exports private input data.
The spectral diagnostics concern the symmetric part of the *projected*
operator, not a certified eigendecomposition of the full matrix.
"""
import argparse
import csv
import math
import sys


def fail(message):
    raise ValueError(message)


def header_and_sizes(source, expected_format):
    header = source.readline().strip().lower().split()
    if (len(header) != 5 or header[:3] != ['%%matrixmarket', 'matrix', expected_format]
            or header[3] not in ('real', 'integer')):
        fail('expected Matrix Market matrix %s real/integer' % expected_format)
    symmetry = header[4]
    if expected_format == 'array' and symmetry != 'general':
        fail('mode arrays must be Matrix Market general')
    if expected_format == 'coordinate' and symmetry not in ('general', 'symmetric', 'skew-symmetric'):
        fail('unsupported matrix symmetry: %s' % symmetry)
    for line in source:
        if line.strip() and not line.lstrip().startswith('%'):
            try:
                dims = list(map(int, line.split()))
            except ValueError:
                fail('invalid Matrix Market dimensions')
            if len(dims) != (3 if expected_format == 'coordinate' else 2):
                fail('invalid number of Matrix Market dimensions')
            return symmetry, dims
    fail('missing Matrix Market dimensions')


def read_modes(path):
    with open(path, 'r', encoding='ascii') as source:
        _, (n, k) = header_and_sizes(source, 'array')
        if not (n > 0 and 0 < k <= 256):
            fail('invalid modes shape')
        columns = []
        for j in range(k):
            col = []
            for i in range(n):
                line = source.readline()
                try:
                    value = float(line.strip())
                except ValueError:
                    fail('invalid mode at column %d row %d' % (j + 1, i + 1))
                if not math.isfinite(value):
                    fail('non-finite mode at column %d row %d' % (j + 1, i + 1))
                col.append(value)
            columns.append(col)
        for line in source:
            if line.strip() and not line.lstrip().startswith('%'):
                fail('extra values in modes array')
        return n, columns


def dot(a, b):
    return math.fsum(x * y for x, y in zip(a, b))


def norm(a):
    return math.sqrt(max(0.0, dot(a, a)))


def axpy_in_place(x, a, y):
    for i in range(len(x)):
        x[i] += a * y[i]


def matmul_sparse(path, n, cols):
    outputs = [[0.0] * n for _ in cols]
    with open(path, 'r', encoding='ascii') as source:
        symmetry, (rows, width, entries) = header_and_sizes(source, 'coordinate')
        if rows != n or width != n:
            fail('A dimension %dx%d does not match modes rows=%d' % (rows, width, n))
        if entries < 0:
            fail('invalid matrix entry count')
        for line_idx in range(entries):
            line = source.readline()
            while line and (not line.strip() or line.lstrip().startswith('%')):
                line = source.readline()
            tokens = line.split()
            if len(tokens) != 3:
                fail('invalid coordinate at entry %d' % (line_idx + 1))
            try:
                i, j = int(tokens[0]) - 1, int(tokens[1]) - 1
                v = float(tokens[2])
            except ValueError:
                fail('invalid numeric coordinate at entry %d' % (line_idx + 1))
            if not (0 <= i < n and 0 <= j < n and math.isfinite(v)):
                fail('invalid matrix coordinate at entry %d' % (line_idx + 1))
            if symmetry == 'skew-symmetric' and i == j and v != 0:
                fail('nonzero diagonal in skew-symmetric coordinate matrix')
            for x, y in zip(cols, outputs):
                y[i] += v * x[j]
                if symmetry == 'symmetric' and i != j:
                    y[j] += v * x[i]
                elif symmetry == 'skew-symmetric' and i != j:
                    y[j] -= v * x[i]
        for line in source:
            if line.strip() and not line.lstrip().startswith('%'):
                fail('more matrix entries than declared')
    return outputs, symmetry, entries


def remove_span(v, basis):
    w = v[:]
    for _ in range(2):
        for q in basis:
            axpy_in_place(w, -dot(q, w), q)
    return w


def jacobi_eigenvalues_symmetric(matrix):
    a = [row[:] for row in matrix]
    n = len(a)
    if n < 2:
        return [a[0][0]] if n else []
    for _ in range(80 * n * n):
        p, q = max(((i, j) for i in range(n) for j in range(i + 1, n)),
                   key=lambda t: abs(a[t[0]][t[1]]))
        if abs(a[p][q]) <= 1e-13 * max(1.0, max(abs(a[i][i]) for i in range(n))):
            break
        theta = 0.5 * math.atan2(2 * a[p][q], a[q][q] - a[p][p])
        c, s = math.cos(theta), math.sin(theta)
        app, aqq, apq = a[p][p], a[q][q], a[p][q]
        a[p][p] = c*c*app - 2*s*c*apq + s*s*aqq
        a[q][q] = s*s*app + 2*s*c*apq + c*c*aqq
        a[p][q] = a[q][p] = 0.0
        for k in range(n):
            if k != p and k != q:
                akp, akq = a[k][p], a[k][q]
                a[k][p] = a[p][k] = c*akp - s*akq
                a[k][q] = a[q][k] = s*akp + c*akq
    return sorted(a[i][i] for i in range(n))


def audit(matrix, baseline, enriched, csv_path, rank_tolerance, match_tolerance):
    n, base = read_modes(baseline)
    n2, extra_modes = read_modes(enriched)
    if n != n2 or len(extra_modes) <= len(base):
        fail('enriched basis must match row count and contain more columns')
    b = len(base)
    max_relative_mismatch = 0.0
    for j, ref in enumerate(base):
        delta = [x - y for x, y in zip(ref, extra_modes[j])]
        relative = norm(delta) / max(norm(ref), 1e-300)
        if not math.isfinite(relative):
            fail('nonfinite baseline comparison')
        max_relative_mismatch = max(max_relative_mismatch, relative)
    if max_relative_mismatch > match_tolerance:
        fail('enriched first %d modes differ from baseline (max relative %.4g)' %
             (b, max_relative_mismatch))
    av, symmetry, entries = matmul_sparse(matrix, n, extra_modes)
    basis, abasis = [], []
    extra_q, extra_aq = [], []
    rows = []
    min_extra_distance = 1.0
    base_rank = 0
    for j, (v, avj) in enumerate(zip(extra_modes, av)):
        size = norm(v)
        if size <= 0:
            fail('zero-length mode %d' % (j + 1))
        rayleigh = dot(v, avj) / (size * size)
        eigen_error = norm([a - rayleigh * z for a, z in zip(avj, v)]) / max(
            norm(avj), abs(rayleigh) * size, 1e-300)
        distance = None
        if j >= b:
            distance = norm(remove_span(v, basis[:base_rank])) / size
            min_extra_distance = min(distance, min_extra_distance)
        w, aw = v[:], avj[:]
        for _ in range(2):
            for q, aq in zip(basis, abasis):
                alpha = dot(q, w)
                axpy_in_place(w, -alpha, q)
                axpy_in_place(aw, -alpha, aq)
        remaining = norm(w)
        independent = int(remaining > rank_tolerance * size)
        if independent:
            inv = 1.0 / remaining
            q = [z * inv for z in w]
            aq = [z * inv for z in aw]
            basis.append(q)
            abasis.append(aq)
            if j < b:
                base_rank += 1
            else:
                extra_q.append(q)
                extra_aq.append(aq)
        rows.append((j + 1, 'base' if j < b else 'extra', size, rayleigh,
                     eigen_error, distance, independent))
    projected = [[dot(q, aq) for aq in extra_aq] for q in extra_q]
    asym = max((abs(projected[i][j] - projected[j][i])
                for i in range(len(projected)) for j in range(len(projected))), default=0.0)
    sympart = [[0.5 * (projected[i][j] + projected[j][i])
                for j in range(len(projected))] for i in range(len(projected))]
    eigs = jacobi_eigenvalues_symmetric(sympart)
    max_orthogonality = max((abs(dot(q, p) - (1.0 if i == j else 0.0))
                             for i, q in enumerate(basis)
                             for j, p in enumerate(basis)), default=0.0)
    with open(csv_path, 'w', encoding='utf-8', newline='') as target:
        writer = csv.writer(target)
        writer.writerow(['mode_index', 'group', 'norm', 'rayleigh', 'relative_eigen_residual',
                         'distance_from_base_span', 'independent_after_previous'])
        for index, group, size, lam, err, distance, independent in rows:
            writer.writerow([index, group, '%.12g' % size, '%.12g' % lam,
                             '%.12g' % err, 'NA' if distance is None else '%.12g' % distance,
                             independent])
    print('MODE_AUDIT_N=%d' % n)
    print('MODE_AUDIT_A_NNZ_FILE=%d' % entries)
    print('MODE_AUDIT_A_MM_SYMMETRY=%s' % symmetry)
    print('MODE_AUDIT_BASE_COUNT=%d' % b)
    print('MODE_AUDIT_EXTRA_COUNT=%d' % (len(extra_modes) - b))
    print('MODE_AUDIT_FIRST_BASE_MAX_REL_DIFF=%.12g' % max_relative_mismatch)
    print('MODE_AUDIT_BASE_RANK=%d' % base_rank)
    print('MODE_AUDIT_EXTRA_NEW_RANK=%d' % len(extra_q))
    print('MODE_AUDIT_EXTRA_MIN_DISTANCE_FROM_BASE=%.12g' % min_extra_distance)
    print('MODE_AUDIT_MAX_ORTH_ERROR=%.12g' % max_orthogonality)
    print('MODE_AUDIT_PROJECTED_ASYMMETRY_MAX=%.12g' % asym)
    print('MODE_AUDIT_EXTRA_SYMPART_NEGATIVE_COUNT=%d' % sum(z < -1e-10 for z in eigs))
    print('MODE_AUDIT_EXTRA_SYMPART_EIGENVALUES=' + ','.join('%.12g' % z for z in eigs))
    print('MODE_AUDIT_CSV=%s' % csv_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('matrix', help='A.mtx, sparse Matrix Market coordinate')
    parser.add_argument('baseline', help='modes04.mtx, dense Matrix Market array')
    parser.add_argument('enriched', help='modes04_12plus6negative.mtx')
    parser.add_argument('--csv-out', required=True)
    parser.add_argument('--rank-tol', type=float, default=1e-10)
    parser.add_argument('--match-tol', type=float, default=1e-10)
    args = parser.parse_args()
    if not (0 < args.rank_tol < 1 and 0 <= args.match_tol < 1):
        parser.error('invalid tolerances')
    try:
        audit(args.matrix, args.baseline, args.enriched, args.csv_out,
              args.rank_tol, args.match_tol)
    except (ValueError, OSError, OverflowError) as exc:
        print('MODE_AUDIT_ERROR=%s' % exc, file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
