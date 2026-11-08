#!/usr/bin/env python3
"""Offline basis experiment: rotate provided extra modes into Ritz directions.

This REQUIRES the extra modes already prepared offline. It does not
discover additional modes from A and the original 12 modes alone.
The symmetric projected operator is 0.5*(Q^T A Q + Q^T A^T Q).
"""
import argparse
import math
import os
import sys

from mfamg_mode_audit import (
    axpy_in_place, dot, header_and_sizes, matmul_sparse, norm, read_modes,
)


def fail(message):
    raise ValueError(message)


def orthonormalize(columns, prior=None, tolerance=1e-10):
    vectors = list(prior) if prior is not None else []
    out = []
    for j, column in enumerate(columns):
        orig = norm(column)
        if orig == 0:
            fail('zero mode %d' % (j + 1))
        w = column[:]
        for _ in range(2):
            for q in vectors:
                axpy_in_place(w, -dot(q, w), q)
        length = norm(w)
        if length <= tolerance * orig:
            fail('dependent mode %d (relative residual %.5g)' %
                 (j + 1, length / orig))
        q = [v / length for v in w]
        vectors.append(q)
        out.append(q)
    return out


def symmetric_eigensystem(matrix):
    """Jacobi diagonalization: sorted eigenvalues and eigenvector columns."""
    n = len(matrix)
    a = [row[:] for row in matrix]
    eigenvectors = [[float(i == j) for j in range(n)] for i in range(n)]
    for _ in range(100 * n * n):
        if n < 2:
            break
        p, q = max(((i, j) for i in range(n) for j in range(i + 1, n)),
                   key=lambda pair: abs(a[pair[0]][pair[1]]))
        scale = max(1.0, max(abs(a[i][i]) for i in range(n)))
        if abs(a[p][q]) < 1e-14 * scale:
            break
        theta = 0.5 * math.atan2(2.0 * a[p][q], a[q][q] - a[p][p])
        c, s = math.cos(theta), math.sin(theta)
        app, aqq, apq = a[p][p], a[q][q], a[p][q]
        a[p][p] = c*c*app - 2*c*s*apq + s*s*aqq
        a[q][q] = s*s*app + 2*c*s*apq + c*c*aqq
        a[p][q] = a[q][p] = 0.0
        for k in range(n):
            if k != p and k != q:
                akp, akq = a[k][p], a[k][q]
                a[k][p] = a[p][k] = c*akp - s*akq
                a[k][q] = a[q][k] = s*akp + c*akq
        for k in range(n):
            vp, vq = eigenvectors[k][p], eigenvectors[k][q]
            eigenvectors[k][p] = c*vp - s*vq
            eigenvectors[k][q] = s*vp + c*vq
    order = sorted(range(n), key=lambda j: a[j][j])
    values = [a[j][j] for j in order]
    weights = [[eigenvectors[i][j] for i in range(n)] for j in order]
    if not all(math.isfinite(z) for z in values):
        fail('nonfinite projected eigenvalues')
    return values, weights


def ritz_combination(columns, weights):
    n = len(columns[0])
    result = [0.0] * n
    for column, w in zip(columns, weights):
        axpy_in_place(result, w, column)
    return result


def write_array(path, columns):
    n = len(columns[0])
    with open(path, 'x', encoding='ascii', newline='\n') as stream:
        stream.write('%%MatrixMarket matrix array real general\n')
        stream.write('% Offline projected Ritz basis; not an online solver mode.\n')
        stream.write('% Generated for local experimental replay only.\n')
        stream.write('%d %d\n' % (n, len(columns)))
        for column in columns:
            if len(column) != n:
                fail('incompatible output column dimension')
            for v in column:
                if not math.isfinite(v):
                    fail('nonfinite output column')
                stream.write('%.17g\n' % v)


def select(args):
    n, base = read_modes(args.baseline)
    n2, enriched = read_modes(args.enriched)
    if n != n2 or len(enriched) <= len(base):
        fail('enriched modes must have more columns at the same matrix dimension')
    m = len(enriched) - len(base)
    if m > 64:
        fail('too many extra modes; this prototype supports up to 64')
    for j, original in enumerate(base):
        delta = norm([x - y for x, y in zip(original, enriched[j])])
        if delta / max(norm(original), 1e-300) > args.match_tol:
            fail('enriched first %d modes differ from baseline' % len(base))
    qb = orthonormalize(base, tolerance=args.rank_tol)
    qe = orthonormalize(enriched[len(base):], prior=qb, tolerance=args.rank_tol)
    aq, symmetry, _ = matmul_sparse(args.matrix, n, qe)
    raw = [[dot(q, aqj) for aqj in aq] for q in qe]
    asym = max((abs(raw[i][j] - raw[j][i]) for i in range(m)
                for j in range(m)), default=0.0)
    sym = [[0.5 * (raw[i][j] + raw[j][i]) for j in range(m)]
           for i in range(m)]
    eigenvalues, eigenweights = symmetric_eigensystem(sym)
    scale = max(1.0, max(abs(z) for z in eigenvalues))
    threshold = args.sign_tol * scale
    if any(abs(z) <= threshold for z in eigenvalues):
        fail('ambiguous near-zero Ritz value; adjust sign_tol or review')
    negatives = [j for j, val in enumerate(eigenvalues) if val < -threshold]
    if not negatives:
        fail('no negative Ritz directions were identified')
    ritz = [ritz_combination(qe, weights) for weights in eigenweights]
    ritz_A = [ritz_combination(aq, weights) for weights in eigenweights]
    print('RITZ_SELECTION_DIM=%d' % n)
    print('RITZ_SELECTION_ORIGINAL_MODES=%d' % len(base))
    print('RITZ_SELECTION_EXTRA_MODES=%d' % m)
    print('RITZ_SELECTION_NEGATIVE_COUNT=%d' % len(negatives))
    print('RITZ_SELECTION_POSITIVE_COUNT=%d' % (m - len(negatives)))
    print('RITZ_SELECTION_PROJECTED_EIGENVALUES=' +
          ','.join('%.12g' % x for x in eigenvalues))
    print('RITZ_SELECTION_PROJECTED_ASYMMETRY=%.12g' % asym)
    print('RITZ_SELECTION_MATRIX_SYMMETRY=%s' % symmetry)
    for j, (val, v, av) in enumerate(zip(eigenvalues, ritz, ritz_A)):
        residual = norm([a - val * x for a, x in zip(av, v)])
        relative = residual / max(norm(av), abs(val) * norm(v), 1e-300)
        print('RITZ_SELECTION_DIRECTION_%d_EIGENRESIDUAL=%.12g' %
              (j + 1, relative))
    folder = args.output_dir
    if not os.path.isdir(folder):
        fail('output directory must already exist: %s' % folder)
    negfile = os.path.join(folder, 'modes_negative_ritz.mtx')
    allfile = os.path.join(folder, 'modes_all_ritz.mtx')
    if os.path.exists(negfile) or os.path.exists(allfile):
        fail('output files exist; choose a new empty output directory')
    if len(base) + m > 256:
        fail('too many columns for replay')
    # Original first 12 columns remain unaltered. Full-space control
    # spans the same 18-dimensional space as the original enriched data.
    write_array(negfile, base + [ritz[j] for j in negatives])
    write_array(allfile, base + ritz)
    print('RITZ_SELECTION_NEGATIVE_OUTPUT=%s' % negfile)
    print('RITZ_SELECTION_ALL_OUTPUT=%s' % allfile)
    print('RITZ_SELECTION_NOTE=offline_rotation_not_online_mode_discovery')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('matrix')
    parser.add_argument('baseline')
    parser.add_argument('enriched')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--rank-tol', type=float, default=1e-10)
    parser.add_argument('--match-tol', type=float, default=1e-10)
    parser.add_argument('--sign-tol', type=float, default=1e-10)
    args = parser.parse_args()
    if not (0 < args.rank_tol < 1 and 0 <= args.match_tol < 1
            and 0 < args.sign_tol < 1):
        parser.error('invalid tolerances')
    try:
        select(args)
    except (ValueError, OSError, OverflowError) as exc:
        print('RITZ_SELECTION_ERROR=%s' % exc, file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
