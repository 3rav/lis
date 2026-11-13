# Optional SuperLU Direct preconditioner for LIS

This is a **separate, opt-in extension**, not part of LIS's default build
and not part of the experimental ILUTP upstream PR #229.

It registers a user-defined LIS preconditioner named **superlu** using
the existing public \`lis_precon_register_ex\` API. One complete sparse
SuperLU factorization is created when LIS builds its preconditioner. All
FGMRES \`M^-1 b\` and adjoint \`M^-H b\` applications reuse the same factors.

This is an *exact sparse direct solve* used through LIS's preconditioner
interface, **not** an incomplete-LU variant. Its memory consumption may be
larger than ILUT/ILUTP, especially with fill-heavy contact systems.

## Build (Linux)

Build LIS as normal, with serial real-double CSR support. Install system
SuperLU and BLAS (for example \`libsuperlu-dev libblas-dev\` on Debian).
From the LIS source root:

\`\`\`sh
cc -std=c99 -O2 -I./include -I./contrib/superlu_direct -DHAVE_CONFIG_H \
  contrib/superlu_direct/testsuperlu_direct.c \
  contrib/superlu_direct/lis_superlu_direct.c \
  src/.libs/liblis.a $(pkg-config --cflags --libs superlu) -lblas -lm \
  -o testsuperlu_direct
./testsuperlu_direct
\`\`\`

For direct validation of the original CCX snapshot pair (no numerical
data is published), compile the standalone validator from the LIS source
root:

```sh
cc -std=c99 -O2 -I./include -I./contrib/superlu_direct -DHAVE_CONFIG_H \
  contrib/superlu_direct/benchsuperlu.c \
  contrib/superlu_direct/lis_superlu_direct.c \
  src/.libs/liblis.a $(pkg-config --cflags --libs superlu) -lblas -lm \
  -o benchsuperlu
OMP_NUM_THREADS=1 ./benchsuperlu A.mtx b.mtx
```

It computes the **physical** `||A*x-b||/||b||` and reports the
FGMRES iteration count. The RHS must be the original CCX Matrix Market
`matrix array real general` vector.

A separate `testsuperlu_singular.c` validates controlled
`LIS_BREAKDOWN` and cleanup for a rank-deficient input; run it under
ASan/UBSan with leak detection where available.

For other applications, compile/link \`lis_superlu_direct.c\` and the
required SuperLU libraries into the application, then register it
**after** \`lis_initialize()\` and **before** setting LIS options:

\`\`\`c
lis_initialize(&argc, &argv);
lis_superlu_direct_register();
lis_solver_set_option("-i fgmres -p superlu -scale none "
                      "-omp_num_threads 1", solver);
\`\`\`

Do not set \`-p superlu\` without registering the extension. Existing LIS
and CCX executables cannot enable this plugin merely by changing an
environment variable: the binary must be explicitly linked and must call
\`lis_superlu_direct_register()\`.

The registration name \`superlu\` is distinct from the built-in ILUT
preconditioner and does not affect existing default options.

## Numerical method

- Converts LIS CSR into privately owned SuperLU CSC.
- Applies SuperLU's \`dgsequ/dlaqgs\` equilibrations.
- Uses COLAMD column ordering and threshold partial **row** pivoting
  (\`DiagPivotThresh=1.0\`).
- Builds full L/U once using \`dgstrf\`.
- Reuses the factorization through \`dgstrs\`, with correct forward and
  transpose row/column scaling.
- Destroys all factors and CSC buffers through LIS's registered
  \`user_destroy\` callback, including partially constructed state.

## Supported and unsupported configurations

Currently supported: serial, single-process **real double CSR**, indices
within the 32-bit SuperLU API range. Unsupported configurations return an
error rather than falling back silently. MPI, complex, FP16, and long-double
support are not implemented. FGMRES remains a serial solve when this
preconditioner is active.

## Validation

Using a local build of LIS and system SuperLU 7.0:

- A nonsymmetric, strongly row-scaled 6x6 matrix: repeated forward and
  transpose solves, maximum *scaled equation residual* 3.9e-16,
  no ASan/UBSan/LeakSanitizer errors.
- A private CCX Newton matrix with 256 contact springs, 60,265 DOF:
  **one FGMRES iteration**, original \`||Ax-b||/||b||=5.12e-13\`,
  one factorization 6.47 CPU seconds, peak process RSS ~631 MiB.
- A private CCX Newton matrix with 187 contact springs, 60,265 DOF:
  **one FGMRES iteration**, original \`||Ax-b||/||b||=1.73e-14\`,
  factorization 4.37 CPU seconds.

The original CCX snapshots are **not stored in the repository**. These
measurements validate linear solves only; nonlinear increment acceptance in
CalculiX has not been tested with this extension.

## Follow-up

A CCX adapter may implement an explicit opt-in fallback on LIS
\`BREAKDOWN\`/non-finite residual, using a freshly built \`-p superlu\`
preconditioner. That fallback is not enabled by this source file alone.
Any such integration must preserve the original Newton A,b, verify the
physical residual, and handle factor memory, repeated Newton solves,
transposition, and sequence cache invalidation explicitly.

Additional optional regression checks: run `LIS_SUPERLU_TEST_BOTH_SCALES=1 ./testsuperlu_direct` to exercise **both row and column** equilibration (`equil=B`) for ordinary and transposed reuse. `testsuperlu_singular.c` checks zero-row and rank-deficient setup rejection plus 100 complete create/destroy cycles under ASan/UBSan/LSan. These checks are included in the serial Linux CI smoke step; they do not change the default LIS build.
