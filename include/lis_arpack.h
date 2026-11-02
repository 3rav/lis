/* Optional ARPACK-NG bridge for real symmetric generalized FE eigenproblems.
 * Copyright (c) 2026 LIS contributors. BSD-3-Clause (same as LIS).
 * This API uses serial, local double vectors, independently of LIS_SCALAR.
 */
#ifndef LIS_ARPACK_H
#define LIS_ARPACK_H

#ifdef __cplusplus
extern "C" {
#endif

typedef int (*lis_arpack_apply_fn)(void *context, int n,
                                    const double *x, double *y);

/* Callback 'solve' solves (K - sigma * M) y = x with a FACTORIZATION
 * prepared and owned by the caller. The bridge never re-factorizes.
 * The callback must return zero on success. In buckling mode M=KG.
 * No global callbacks/state are used (independent problems are reentrant).
 */
typedef struct {
    int n;
    int mode;                /* 3: frequency, 4: buckling */
    double sigma;
    void *context;
    lis_arpack_apply_fn stiffness;    /* K*x */
    lis_arpack_apply_fn second;       /* M*x (mode 3), KG*x (mode 4) */
    lis_arpack_apply_fn solve;        /* (K-sigma*M)^-1*x */
} lis_arpack_problem;

typedef struct {
    int nev;                 /* 1 <= nev < n */
    int ncv;                 /* 0 = automatic, otherwise nev < ncv <= n */
    int maxiter;             /* 0 = 300 implicit restarts */
    double tol;              /* 0 = machine precision */
    const char *which;       /* NULL = "LM" (nearest shift) */
} lis_arpack_options;

typedef struct {
    int nconv;
    int iterations;
    int op_count;
    int b_count;
    int solve_count;
    int arpack_info;          /* raw ARPACK status, zero on success */
} lis_arpack_stats;

enum {
    LIS_ARPACK_SUCCESS = 0,
    LIS_ARPACK_BAD_ARGUMENT = 1,
    LIS_ARPACK_NO_MEMORY = 2,
    LIS_ARPACK_CALLBACK_ERROR = 3,
    LIS_ARPACK_BACKEND_ERROR = 4,
    LIS_ARPACK_NO_CONVERGENCE = 5
};

/* Output buffers: eigenvalues[nev], eigenvectors[n*nev] (column-major),
 * relative_residuals[nev]. Only the first stats->nconv entries are valid.
 * Residual = ||K*x - lambda*M*x|| / (||K*x|| + |lambda|*||M*x||).
 * All buffers belong to the caller. No implicit matrix storage/conversion.
 */
int lis_arpack_symmetric(const lis_arpack_problem *problem,
                         const lis_arpack_options *options,
                         double *eigenvalues, double *eigenvectors,
                         double *relative_residuals, lis_arpack_stats *stats);

#ifdef __cplusplus
}
#endif
#endif
