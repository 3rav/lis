/* ARPACK-NG symmetric reverse-communication bridge for CalculiX-style FE.
 * Copyright (c) 2026 LIS contributors. BSD-3-Clause.
 */
#include "lis_arpack.h"
#include <arpack.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int valid_which(const char *which)
{
    return strcmp(which, "LM") == 0 || strcmp(which, "SM") == 0 ||
           strcmp(which, "LA") == 0 || strcmp(which, "SA") == 0 ||
           strcmp(which, "BE") == 0;
}

static double norm2(int n, const double *x)
{
    double scale = 0.0, ssq = 1.0;
    int i;
    for (i = 0; i < n; ++i) {
        double a = fabs(x[i]);
        if (a != 0.0) {
            if (scale < a) { double t = scale/a; ssq = 1.0 + ssq*t*t; scale = a; }
            else { double t = a/scale; ssq += t*t; }
        }
    }
    return scale == 0.0 ? 0.0 : scale*sqrt(ssq);
}

int lis_arpack_symmetric(const lis_arpack_problem *p,
                         const lis_arpack_options *o,
                         double *values, double *vectors,
                         double *residuals, lis_arpack_stats *stats)
{
    a_int ido = 0, info = 0, iparam[11] = {0}, ipntr[11] = {0};
    a_int n, nev, ncv, lworkl, *select = NULL;
    double *resid = NULL, *v = NULL, *workd = NULL, *workl = NULL;
    double *d = NULL, *rhs = NULL, *ax = NULL, *bx = NULL;
    const char *which;
    size_t wn, wv, wl;
    int status = LIS_ARPACK_SUCCESS, i, j;

    if (stats) memset(stats, 0, sizeof(*stats));
    if (!p || !o || !values || !vectors || !residuals || !stats ||
        !p->stiffness || !p->second || !p->solve ||
        p->n < 2 || o->nev < 1 || o->nev >= p->n ||
        (p->mode != 3 && p->mode != 4) ||
        (p->mode == 4 && p->sigma == 0.0) ||
        !isfinite(p->sigma) ||
        !isfinite(o->tol) || o->tol < 0.0 || o->maxiter < 0 ||
        p->n > INT_MAX/3)
        return LIS_ARPACK_BAD_ARGUMENT;
    which = o->which ? o->which : "LM";
    if (!valid_which(which)) return LIS_ARPACK_BAD_ARGUMENT;
    n = (a_int)p->n;
    nev = (a_int)o->nev;
    if (o->ncv) ncv = (a_int)o->ncv;
    else {
        int target = o->nev > (INT_MAX-8)/2 ? p->n : 2*o->nev + 8;
        ncv = (a_int)(target < p->n ? target : p->n);
    }
    if (ncv <= nev || ncv > n || ncv > INT_MAX/11) return LIS_ARPACK_BAD_ARGUMENT;
    wn = (size_t)n;
    wv = wn*(size_t)ncv;
    wl = (size_t)ncv*((size_t)ncv + 8u);
    if (wv > SIZE_MAX/sizeof(double) || wl > SIZE_MAX/sizeof(double) ||
        3u*wn > SIZE_MAX/sizeof(double) || wl > INT_MAX)
        return LIS_ARPACK_BAD_ARGUMENT;
    lworkl = (a_int)wl;
    resid = (double *)calloc(wn, sizeof(double));
    v = (double *)calloc(wv, sizeof(double));
    workd = (double *)calloc(3u*wn, sizeof(double));
    workl = (double *)calloc(wl, sizeof(double));
    d = (double *)calloc((size_t)nev+1u, sizeof(double));
    rhs = (double *)malloc(wn*sizeof(double));
    ax = (double *)malloc(wn*sizeof(double));
    bx = (double *)malloc(wn*sizeof(double));
    select = (a_int *)calloc((size_t)ncv, sizeof(a_int));
    if (!resid || !v || !workd || !workl || !d || !rhs || !ax || !bx || !select) {
        status = LIS_ARPACK_NO_MEMORY;
        goto cleanup;
    }
    iparam[0] = 1; /* ARPACK exact implicit shifts */
    iparam[2] = o->maxiter ? (a_int)o->maxiter : 300;
    iparam[3] = 1;
    iparam[6] = (a_int)p->mode;

    for (;;) {
        double *y, *x;
        lis_arpack_apply_fn b = p->mode == 3 ? p->second : p->stiffness;
        dsaupd_c(&ido, "G", n, which, nev, o->tol,
                 resid, ncv, v, n, iparam, ipntr,
                 workd, workl, lworkl, &info);
        if (ido == 99) break;
        if (ido != -1 && ido != 1 && ido != 2) {
            status = LIS_ARPACK_BACKEND_ERROR;
            goto cleanup;
        }
        if (ipntr[0] < 1 || ipntr[0] > 2*n+1 ||
            ipntr[1] < 1 || ipntr[1] > 2*n+1) {
            status = LIS_ARPACK_BACKEND_ERROR;
            goto cleanup;
        }
        x = workd + ipntr[0] - 1;
        y = workd + ipntr[1] - 1;
        if (ido == 2) {
            if (b(p->context, p->n, x, y)) {
                status = LIS_ARPACK_CALLBACK_ERROR;
                goto cleanup;
            }
        } else {
            if (ido == 1) {
                if (ipntr[2] < 1 || ipntr[2] > 2*n+1) {
                    status = LIS_ARPACK_BACKEND_ERROR;
                    goto cleanup;
                }
                memcpy(rhs, workd + ipntr[2] - 1, wn*sizeof(double));
            } else if (b(p->context, p->n, x, rhs)) {
                status = LIS_ARPACK_CALLBACK_ERROR;
                goto cleanup;
            }
            if (p->solve(p->context, p->n, rhs, y)) {
                status = LIS_ARPACK_CALLBACK_ERROR;
                goto cleanup;
            }
            stats->solve_count++;
        }
    }
    stats->arpack_info = (int)info;
    stats->nconv = (int)iparam[4];
    stats->iterations = (int)iparam[2];
    stats->op_count = (int)iparam[8];
    stats->b_count = (int)iparam[9];
    if (info < 0) {
        status = LIS_ARPACK_BACKEND_ERROR;
        goto cleanup;
    }
    if (stats->nconv > 0) {
        a_int ierr = 0;
        dseupd_c(1, "A", select, d, vectors, n, p->sigma,
                 "G", n, which, nev, o->tol, resid, ncv, v, n,
                 iparam, ipntr, workd, workl, lworkl, &ierr);
        if (ierr != 0) {
            stats->arpack_info = (int)ierr;
            stats->nconv = 0;
            status = LIS_ARPACK_BACKEND_ERROR;
            goto cleanup;
        }
        for (j = 0; j < stats->nconv && j < o->nev; ++j) {
            double na, nb, nr;
            double lambda = d[j];
            double *z = vectors + (size_t)j*wn;
            if (!isfinite(lambda) ||
                p->stiffness(p->context, p->n, z, ax) ||
                p->second(p->context, p->n, z, bx)) {
                stats->nconv = j;
                status = LIS_ARPACK_CALLBACK_ERROR;
                goto cleanup;
            }
            na = norm2(p->n, ax);
            nb = norm2(p->n, bx);
            for (i = 0; i < p->n; ++i) ax[i] -= lambda*bx[i];
            nr = norm2(p->n, ax);
            values[j] = lambda;
            residuals[j] = nr/(na + fabs(lambda)*nb + 1e-300);
        }
    }
    if (info != 0 || stats->nconv < o->nev)
        status = LIS_ARPACK_NO_CONVERGENCE;
cleanup:
    free(resid); free(v); free(workd); free(workl); free(d);
    free(rhs); free(ax); free(bx); free(select);
    return status;
}
