/* Native restarted symmetric Krylov eigensolver for LIS.
 * BSD-3-Clause: Copyright (c) 2026 LIS contributors.
 * This self-contained API does not depend on an external eigensolver.
 */
#ifndef LIS_KRYLOV_H
#define LIS_KRYLOV_H
#ifdef __cplusplus
extern "C" {
#endif

typedef int (*lis_krylov_matvec)(void *ctx, int n, const double *x, double *y);

/* A x = lambda B x; B is the identity in STANDARD mode.
 *  STANDARD: symmetric A, metric I, operator A.
 *  SHIFT_INVERT_B: SPD B, operator (A-sigma*B)^-1 B, metric B.
 *  SHIFT_INVERT_A: SPD A, symmetric B, operator (A-sigma*B)^-1 A,
 *                  metric A. sigma must be nonzero.
 * 'solve' supplies the inverse action of a previously factorized shifted
 * matrix. The implementation never assembles or factorizes that matrix.
 * All callbacks must accept nonaliasing input/output arrays.
 */
enum {
    LIS_KRYLOV_STANDARD = 0,
    LIS_KRYLOV_SHIFT_INVERT_B = 1,
    LIS_KRYLOV_SHIFT_INVERT_A = 2
};
typedef struct {
    int n;
    int transform;
    double sigma;
    void *ctx;
    lis_krylov_matvec apply_a;
    lis_krylov_matvec apply_b;
    lis_krylov_matvec solve;
} lis_krylov_problem;

typedef struct {
    int nev;               /* 1 <= nev < n */
    int ncv;               /* 0: automatic; otherwise nev < ncv <= n */
    int maxiter;           /* 0: 300 Rayleigh-Ritz expansion cycles */
    double tol;            /* 0: 1e-10 relative original residual */
    const char *which;     /* LM (default), SM, LA, SA; transformed spectrum */
} lis_krylov_options;

typedef struct {
    int nconv;
    int iterations;
    int restarts;
    int op_count;
    int metric_count;
    int solve_count;
    double max_residual;
} lis_krylov_stats;

enum {
    LIS_KRYLOV_SUCCESS = 0,
    LIS_KRYLOV_BAD_ARGUMENT = 1,
    LIS_KRYLOV_NO_MEMORY = 2,
    LIS_KRYLOV_CALLBACK_ERROR = 3,
    LIS_KRYLOV_BREAKDOWN = 4,
    LIS_KRYLOV_NO_CONVERGENCE = 5
};

/* Sought eigenpairs ordered by 'which'. Arrays belong to caller:
 * eigenvalues[nev], eigenvectors[n * nev] (column-major), residuals[nev].
 * Only [0,stats->nconv) entries are certified on nonzero return.
 * Original-equation residual = ||A*x - lambda*B*x|| /
 * (||A*x|| + |lambda|*||B*x||); B=I for STANDARD.
 */
int lis_krylov_symmetric(const lis_krylov_problem *problem,
                         const lis_krylov_options *options,
                         double *eigenvalues, double *eigenvectors,
                         double *residuals, lis_krylov_stats *stats);
#ifdef __cplusplus
}
#endif
#endif
