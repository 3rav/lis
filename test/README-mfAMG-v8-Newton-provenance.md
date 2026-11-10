# MFAMG v8 / PR #47 — provenance for real consecutive Newton systems

**Research protocol only; not permission to merge or edit production CalculiX.**
The existing three-step runner can validate a *user-declared* matrix sequence,
but **cannot certify where that sequence came from**. A successful CSV
preflight or three-solve replay is never proof of Newton chronology.

## Local source/model compatibility probe (read-only)

Before attempting any instrumented export, test the exact local CalculiX
source tree and the original input deck **without copying or uploading them**.
Use MSYS2 UCRT64 and substitute your PRIVATE paths only in the terminal:

```bash
CCX_ROOT="/local/private/CalculiX-source"
CCX_INP="/local/private/model.inp"
(
  set -euo pipefail
  cd "$HOME/lis-v8"
  BRANCH=experiment/mfamg-fgcrodr-replay-v8-20261108
  git fetch --no-tags origin "refs/heads/$BRANCH:refs/remotes/origin/$BRANCH"
  git show "refs/remotes/origin/$BRANCH:test/mfamg_v8_ccx_capture_probe.py" |
    python3 - --source-root "$(cygpath -m "$CCX_ROOT")"               --input-file "$(cygpath -m "$CCX_INP")"
)
```

The probe prints **only** aggregate `CCX_CAPTURE_*` values. It checks that
`src/nonlingeo.c` (or `nonlingeo.c`) exists, contains the expected
assembly/Newton-loop signatures, and reports coarse source solver signatures
and the **top-level only** `*STEP`, `*STATIC`, `*DYNAMIC`, contact,
spring and `*INCLUDE` input directives. Included files and runtime solver
selection are **not** resolved. It verifies both files stayed byte-for-byte
unchanged, makes no outputs and never prints their paths or hashes.

`CCX_CAPTURE_PROBE=SOURCE_AND_MODEL_READABLE` is NOT certification of
the correct interception point, physical chronology, a valid .inp, or that
three Newton linear systems actually occur. The exact solver branch and
Matrix Market serialization must be reviewed against that LOCAL version
**before** editing an isolated experimental source copy. In particular,
the stiffness assembly is not automatically the final operator and RHS
consumed by the solver. Never instrument or rebuild the production tree.

## Local isolated source-only staging (no build / no export)

The local 31-spring installation checks yielded two candidate source
directories with 779 matching immediate source/build files and available
GCC, GFortran and Make. These results justify an **isolated staging copy**,
but do not prove that either directory was used to build the installed EXE.

`test/mfamg_v8_ccx_isolated_stage.py` copies allowed C/Fortran source
files, headers, include files and root Makefiles from the selected
candidate to a fresh PRIVATE directory. It checks both selected and peer
source files are byte-identical, hashes inputs before and after copying,
never overwrites an existing output, and does not copy the model,
installed EXE or MKL libraries.

Run in the original MSYS2 UCRT64 terminal using the SAME candidate
enumeration used for the preceding local probe:

```bash
(
set -euo pipefail
CCX_DIR="$(cygpath -u 'C:\Program Files (x86)\bConverged\CalculiX')"
CCX_EXE="$CCX_DIR/ccx/ccx.exe"
CCX_MKL="$CCX_DIR/MKL_PARDISO"
DOCS="$(cygpath -u 'C:\Users\rafal\Documents')"
mapfile -t SOURCES < <(
  find "$CCX_DIR" "$HOME" "$DOCS" -maxdepth 8 -type f \
    -name 'nonlingeo.c' -print 2>/dev/null | sort -u
)
if (( ${#SOURCES[@]} != 3 )); then
  echo CCX_STAGE=BLOCKED_SOURCE_COUNT_CHANGED
  exit 2
fi
cd "$HOME/lis-v8"
BRANCH=experiment/mfamg-fgcrodr-replay-v8-20261108
git fetch --no-tags origin "refs/heads/$BRANCH:refs/remotes/origin/$BRANCH"
git show "refs/remotes/origin/$BRANCH:test/mfamg_v8_ccx_isolated_stage.py" |
  python3 - \
    --source-dir "$(cygpath -m "${SOURCES[1]%/nonlingeo.c}")" \
    --peer-dir "$(cygpath -m "${SOURCES[2]%/nonlingeo.c}")" \
    --install-dir "$(cygpath -m "$CCX_DIR")" \
    --installed-exe "$(cygpath -m "$CCX_EXE")" \
    --mkl-dir "$(cygpath -m "$CCX_MKL")" \
    --out-dir "$(cygpath -m "$HOME/mfamg_replay/ccx_newton_research_stage_20261110")"
)
```

Print and share only `CCX_STAGE_*` aggregate lines.
This stage **does write** a new local private source-only snapshot; it
does not modify anything in the bConverged installation. Do not run
`make` on it until the Makefiles' toolchain and MKL/PARDISO settings
have been reviewed, and never use its mere existence as evidence of
runtime solver identity or a genuine Newton trajectory.

## Source of truth: the ORIGINAL nonlinear solve

When an authentic sequence is available (or when a future research-only
instrumented CalculiX run exports it), retain a **local, private producer-side
witness** that binds each exported A/b pair to its position in the actual
Newton/contact computation. Record, at the moment of export:

- Identity of one original nonlinear run and one load/time increment
  (for example a run UUID and increment identifier), established by the
  exporting process rather than reconstructed from filenames afterward.
- Actual consecutively increasing **linear-system assembly/solve ordinals**
  within that run, the corresponding Newton iteration numbers, and whether
  any additional linear system was assembled/solved between them. If the
  producer performs retries or repeated linear solves in one Newton iteration,
  record these explicitly; do not omit them to manufacture adjacency.
- An original producer log entry linking each ordinal to **both** exported
  Matrix Market A and RHS b. Capture SHA-256 of each file *at export time*,
  not solely later during analysis. Preserve the original log privately.
- The exact model/configuration, export point, DOF mapping/permutation and
  any rejected/failed Newton attempts relevant to the three selected systems.
  Only use the same topology/DOF indexing if supported; an actual changed CSR
  graph is permitted, but the solver's safety policy must handle it.

**Never** infer order from `A01/A02/A03`, file timestamps, number of rows,
CSR graph compatibility, residual history, or the fact that files have
similar numerical values. These are not evidence of physical succession.

## Private hand-off to the existing read-only preflight

After inspecting the **original producer-side witness**, choose three
*actually adjacent* linear solves from that same nonlinear run/increment.
Create a PRIVATE UTF-8 CSV with exactly this header and three rows:

```csv
step,matrix,rhs
```

Each `step` is a unique safe label (letters, digits, `_`, `-`);
each `matrix` and `rhs` is the path to the exact exported pair corresponding
to that ordinal (absolute paths or paths relative to the CSV file).
The runner does not read a separate witness file; keep the mapping from
CSV rows to original producer ordinals and SHA-256 in a separate private
record and recheck hashes before the replay.

If the witness is missing, contradictory, contains gaps or cannot bind both
A and b to each ordinal, **STOP**:
`V8_NEWTON_GATE=NEEDS_EXPLICIT_CURATED_STEPS_CSV` remains an unmet
research-data gate. A manually typed CSV by itself is **not** certification.

Run `test/mfamg_v8_newton_triplet_one_shot.sh` in default (preflight-only)
mode first. It checks existence, dimensions, CSR transitions and
unchanged local inputs, then prints
`V8_NEWTON_CHRONOLOGY=CSV_ORDER_DECLARED_NOT_INDEPENDENTLY_PROVEN`.
Only after independent confirmation of original provenance should one
consider opt-in execution with BOTH
`MFAMG_V8_NEWTON_EXECUTE=1` and
`MFAMG_V8_NEWTON_CONSECUTIVE_CONFIRMED=1`; this flag records the
operator's explicit confirmation, it does **not** prove it.

Three consecutive systems in one replay `LIS_SOLVER` are still an offline
approximation, not a full nonlinear CalculiX run or evidence of speedup.
Compare actual residuals, `AZ=C`, KEEP/TRANSPORT/REJECT and the fresh
third solve before any separate in-application validation.

## Privacy and scope

Do NOT upload the witness, contact matrices, RHS, modes, file hashes,
original raw logs or locally identifying paths to GitHub or chat. Share
only aggregate `V8_NEWTON_GATE`, `V8_NEWTON_PREFLIGHT`, and eventual
solver-result classifications. No automatic harvesting, filename guessing,
production source modification, merging or reset/clean of local checkout.
If the original nonlinear-run export is unavailable, leave this gate
**NOT TESTED**, and the PR **Draft / NO MERGE**.
