/* Experimental application-side reference driver; not an installed LIS API. */
#ifndef LIS_ANCHOR_PILOT_H
#define LIS_ANCHOR_PILOT_H
#include "lis.h"

typedef struct lis_anchor_pilot lis_anchor_pilot;
typedef struct {
    LIS_INT raw_status, pilot_iterations, continuation_iterations, fresh_iterations;
    int accepted, continued, fallback, factorizations;
    double relative_residual, ratio, setup_seconds;
} lis_anchor_pilot_report;

/* Serial real-double, unscaled assembled CSR, one active thread.
 * P is borrowed and must stay alive and unchanged until destroy. A/b/x
 * share its DOF ordering. The driver owns the solver and factors, starts
 * from zero for each RHS, and never modifies x unless physically accepted.
 * Native numeric statuses are recorded separately from certified acceptance.
 * API/argument failures are not retried. Recreate after a hard solve error. */
LIS_INT lis_anchor_pilot_create(LIS_MATRIX P, lis_anchor_pilot **out);
LIS_INT lis_anchor_pilot_solve(lis_anchor_pilot *ctx, LIS_MATRIX A,
    LIS_VECTOR b, LIS_VECTOR x, lis_anchor_pilot_report *report);
/* Matched fresh-factorization reference, same options and certification. */
LIS_INT lis_anchor_fresh_solve(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x,
    lis_anchor_pilot_report *report);
void lis_anchor_pilot_destroy(lis_anchor_pilot *ctx);
#endif
