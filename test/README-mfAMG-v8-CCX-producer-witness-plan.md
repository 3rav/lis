# CCX producer-side Newton witness: final A/b capture integration plan

Scope: **research-only isolated CalculiX/LIS integration**. No changes to the
production `ccx.exe`, no PR merge, no automatic edit to local private sources.
FGMRES remains production default. This document is a **patch contract, not
a patch claiming compatibility with the installed binary**.

## Evidence that motivates the patch (2026-11-11)

Local-only, lexical inspection of 1,311 C/Fortran sources found exactly one
definition of `lisccx_snapshot_capture` in one source file, and two static
callsite candidates:
* final solver input: `solve_A, ccx_lis_b`, literal `"newton"` tag;
* separate state input: `state->A, state->b`, `state->profile_label` tag.

Neither capture-call argument list nor capture function definition receives
`istep`, `iinc`, cutback count, `iit` or a producer-side solve ordinal.
The exporter has no manifest or SHA-256 markers. There are 10 complete
numbered local A/b pairs but only 5 recorded LIS linear solves, all logged
with the same Newton number. **Static calls are not runtime call counts,
and 10 pairs are not 10 proven Newton iterations.**

The exact relationship between the inspected source and the executed
research `ccx.exe` is **NOT ESTABLISHED**. Until a compatible source and
isolated build are verified, the integration below is **BLOCKED**, and no
existing export can be reclassified as a certified Newton chronology.

## Minimal producer-side integration — acceptance contract

1. In the **actual nonlinear Newton driver**, assign one monotonically
   increasing `linear_solve_seq` **at entry to every attempted linear solve**.
   Obtain `istep`, `iinc`, actual Newton iteration `iit`, and cutback /
   retry generation from the producer's authoritative state. If the driver
   lacks a cutback number, instrument the code that actually requests each
   retry. **Do not use file numbers, file times, or guessed counters.**
2. Thread these values through an explicit `producer_context` struct or
   parameter list to the **final** adapter call. A separate source-tag
   constant must identify `FINAL_SOLVE_A_B`. Use the final `solve_A`
   and `ccx_lis_b` that will be passed to `lis_solve()`, after all actual
   assembly/transformations that the solver consumes. Include a stable
   DOF-ordering/map fingerprint (the same n alone is insufficient).
3. Make precisely **one authoritative export per attempted linear solve**.
   The `state->A,state->b` path may remain as a distinctly tagged *debug*
   export, but it must never share names or ordinals with authoritative
   Newton captures; in the research producer it should normally be disabled.
   Do not delete it without reviewing runtime control flow and any users
   of the state-snapshot behavior.
4. For every `linear_solve_seq`, record a private producer-side `BEGIN`
   event **before** capture and a terminal `END` or `ABORT` event after
   capture/solve. This preserves all failed attempts, cutbacks, repeated
   solves and skipped writes. Do not suppress failures to manufacture
   adjacent successful Newton iterations.
5. When both Matrix Market files have been successfully written and closed,
   record their exact paths and **SHA-256 measured immediately by the exporter**
   in the `END` event. Use one unambiguous run identifier, monotonically
   increasing `export_seq`, and an append-only manifest. Never write a
   success event for a partially written pair. On I/O error, record `ABORT`
   and fail closed. Ensure one event line cannot be interleaved by concurrent
   writers (single owner/writer, lock or equivalent) and flush it.
6. A minimal private event format (CSV with proper quoting or JSONL) must
   preserve, without reconstruction:

   `run_id,linear_solve_seq,event,istep,iinc,icutb,iit,export_seq,tag,n,nnz,dof_map_id,A_path,b_path,A_sha256,b_sha256,solve_status`

   Every field is from the producer at the indicated event. Fields not
   yet known at `BEGIN` must be null, not fabricated. Include active spring
   count when the nonlinear driver actually supplies it; absence is UNKNOWN.
   Record solver return code and step/cutback outcome, including nonzero
   or abandoned attempts.
7. Make capture and manifest **explicit opt-in** for a private experimental
   executable. Do not upload matrices, RHS, model, modes, manifest, file
   hashes, source tree or raw logs to GitHub or chat. Only share aggregate
   status classifications. Never assume `LIS_MFAMG_REPLAY_METHOD` controls
   the live CalculiX executable; it is an *offline replay* variable.

## Validation sequence: minimum evidence, no unnecessary runs

**Gate A — source correspondence.** Establish the matching LIS-enabled
adapter/driver source and reproducible isolated build before changing any
source. Prior 779-file staged trees had zero LIS markers at key dispatch
sites; a candidate 44-file patch failed clean `git apply --check`.
Those checks are not an established build base. **No force, fuzzy patch,
`--reject`, or production overwrite.**

**Gate B — synthetic integration.** In a *private isolated copy* only,
compile the new producer witness; run a small known converging model with
snapshot capture enabled. Require a 1:1 correspondence between *actual*
linear solve entries and unique authoritative final A/b manifests, matching
on-disk hashes, no collisions, terminal events and no dropped retries.
Validate a failure/cutback case to ensure aborted solves are visible.
Do **not** claim that simple synthetic iterations prove contact chronology.

**Gate C — one real, converging nonlinear run.** Use a genuinely converging
live solver in the research copy (e.g. an actually verified MFAMG profile,
or a properly integrated direct solver). The earlier ILU(2) configuration
failed and led to repeated `newton=2` entries; it is not a suitable
provenance producer. Preserve the full run and manifest locally. For a
candidate triple, require: one run, one physical increment and cutback
generation, `iit` values increasing by exactly one, consecutive
`linear_solve_seq` **with no intervening events**, same DOF mapping, both
files present and unchanged since export, and confirmed terminal results.
Verify any contact-induced CSR changes separately.

**Gate D — research cross-system replay.** Only after Gate C, use a single
persistent solver on those *proven* systems with opt-in FGCRO-DR. Compare
fresh and recycled solves, true relative residuals, recycle
KEEP/TRANSPORT/REJECT, graph policy and fully converged wall time.
The currently enforced CSR structural-change REJECT policy must not be
silently disabled for contact. Passing a three-row user-created CSV or
`MFAMG_V8_NEWTON_CONSECUTIVE_CONFIRMED=1` is not provenance proof.

Status as of this document: single-matrix experimental normwise gate PASS
(`EXIT=0`, 2138 iters, true residual `9.860966704249e-06`);
**Gate A, B, C, D remain NOT CERTIFIED**. PR #47 is Draft / NO MERGE.
