# Native LIS SuperLU fallback evaluation on CCX contact matrices

## Scope

This is an **optional external LIS preconditioner** using the SuperLU
system library. It does **not** alter default LIS/ILUT(P) behavior, add a
mandatory dependency, or enable automatic CalculiX fallback by itself.

The tests use the *original* CCX Newton `A.mtx` and `b.mtx` from
two private 60,265-DOF contact states. The original input files remain
private; this report publishes numerical aggregate results only.

## Native LIS + SuperLU (this contribution)

With real-double CSR, FGMRES, one thread, `-scale none`,
`-p superlu`, `restart=10`, `maxiter=20`, target `tol=1e-10`:

| Contact | FGMRES iterations | Factor CPU seconds | Solve CPU seconds | Original `||Ax-b||/||b||` |
|---|---:|---:|---:|---:|
| **256 springs** | **1** | **6.468** | **0.068** | **5.12e-13** |
| **187 springs** | **1** | **4.368** | **0.048** | **1.73e-14** |

For 256 springs a complete command timed with GNU time used approximately
**8.63 seconds wall clock** and **645,980 KiB maximum process RSS**
(~631 MiB). Full LU requires substantially more fill and may use much
more memory on larger or differently ordered models.

This is a genuinely native **LIS + C SuperLU** integration test. No
SciPy routine solves the matrix in this section. The plugin is
registered as `-p superlu`, factors once, and reuses its factors in
the FGMRES preconditioner callback.

## Robustness tests

- A strongly row-scaled nonsymmetric 6x6 matrix: two repeated forward
  and two transpose solves, maximum equation residual < 4e-16.
- AddressSanitizer / UndefinedBehaviorSanitizer / LeakSanitizer:
  no errors on the forward/transpose/reuse test.
- Rank-deficient 3x3 matrix: SuperLU reports a singular pivot;
  the callback returns `LIS_BREAKDOWN` and releases partial
  factorization state under sanitizer instrumentation.
- Continuous integration compiles and runs the contributed test
  on Linux with explicit `libsuperlu-dev` + BLAS. Other LIS
  configurations remain independent of SuperLU.

## Earlier iterative baseline (separate LIS experiments)

For 187 springs, ILUTP `tau=0.5` converged in 76 iterations,
true relative residual 9.53e-6. For 256 springs, ILUT/ILUTP with
the same contact profile suffered BREAKDOWN after 1–2 iterations
and produced a nonfinite Arnoldi process. All stored incomplete
factors were finite, but triangular applications amplified the
preconditioned vector to magnitudes of about 1e230.

A diagonal-shifted preconditioner of
`M=A+alpha*diag(max_j |A_ij|)`, still solving original `Ax=b`,
stabilized FGMRES but did not converge at 1e-5. For alpha=3.0
the physical residual was 0.01071 after 1000 iterations; for
alpha=10 it was 0.00986.

A weighted maximum-product row matching changed 114 out of 60,265
positions in the 256-spring case and did not prevent BREAKDOWN.
Spectral diagnostics (independent SciPy only) found 30 negative
and 34 positive eigenvalues among 64 closest to zero in that case,
versus one negative and 63 positive for 187 springs. Deflating
64 near-zero modes did not materially change MINRES convergence.

## Integration boundary

The contributed module must be **linked into** the target application
and `lis_superlu_direct_register()` called after `lis_initialize()`.
Only then does `-p superlu` work. Existing CCX binaries do
not automatically gain this feature by setting an environment variable.

A future **separately tested** CCX fallback should detect a real
numerical ILUT(P) failure, rebuild a fresh direct preconditioner,
keep the exact Newton operator and RHS, check the final physical
residual, and enforce explicit memory/failure controls. The upstream
ILUTP PR should remain a distinct, experimental change.

Do not confuse success in solving an individual Newton linearization
with acceptance or convergence of the complete nonlinear CCX increment.
