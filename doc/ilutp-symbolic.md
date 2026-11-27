# Opt-in symbolic complete-fill LU

This proposal extracts the recovered RC32 scalar symbolic factorization,
exact minimum-degree ordering, direct permuted CSR assembly and numerical
refresh. It does not import the historical benchmark archive or GLI patches.

Configuration (all new options default to 0, accept only 0/1):

    -i fgmres -p ilut -ilutp 1 -ilutp_full 1 -ilutp_threshold 0
    -ilutp_symbolic 1 -ilutp_md 1 -ilutp_reuse 1

Symbolic factors require a sorted, duplicate-free, structurally symmetric
CSR matrix in a non-MPI real-double int32 build, no scaling, half storage,
RAS, Schur, or additive Schwarz. Numerical factorization uses the original
scalar row kernel and does not dynamically pivot. Explicit symbolic
selection rejects unsupported inputs instead of silently selecting a
different factorization. OpenMP outer Krylov operations may use threads;
triangular substitutions remain serial. Minimum degree is optional;
reuse requires minimum degree and frozen-pattern sequence refresh policy.

The cache belongs to the preconditioner, not application user_data.
A factor validity flag is set only on successful construction/refresh;
failed refresh invalidates the factor and the sequence caller destroys it
before rebuilding. Solver sequence configuration snapshots include all new
options. No environment variable controls the symbolic path.

The new fields are appended to the public solver/preconditioner structs:
existing field offsets are preserved but sizeof changes. This requires
ABI review and rebuilding clients, not a claim of binary compatibility.

The grouped 2/3/4/8 kernels, DOT4, RUNS and parallel numerical variants
remain in the historical RC32 branch and are not included here.
No performance claim for the private CCX fixture is made.
The focused CI compares known solutions, fresh factors and sequence
solutions, checks state and user-data isolation, and exercises failure
recovery. Broader structural patterns and measured performance are still
needed before production acceptance.
