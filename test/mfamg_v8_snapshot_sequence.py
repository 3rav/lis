#!/usr/bin/env python3
"""LOCAL pairwise real-contact sequence probe. No private uploads.

Each transition (A_i,b_i) -> (A_{i+1},b_{i+1}) is tested in one
LIS_SOLVER *within that pair*. The solver is new at the beginning of
each pair; this is NOT a persistent multi-step Newton trajectory.

The binary itself validates true residuals and A*Z=C. Here we verify
the exact CSR preflight, policy counters, diagnostics, provenance, and
the optional independent fresh-A2 control. No speedup claim follows
from one pass.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys

from mfamg_mode_audit import read_modes
from mfamg_replay_v8_repeats import ENV

FIELDS = ("pair", "from_step", "to_step", "graph", "mode", "changed_nnz",
          "exit", "first_status", "first_true_relative", "second_status",
          "second_true_relative", "second_iter", "second_wall_seconds",
          "cross_uses_delta", "keeps_delta", "transports_delta", "rejects_delta",
          "orth_error", "relation_error", "fresh_exit", "fresh_status",
          "fresh_true_relative", "fresh_iter", "fresh_wall_seconds",
          "pair_pass", "error", "preflight_log", "replay_log", "fresh_log")

STRUCTURE_ERRORS = (
    "replace-values CSR row graph mismatch",
    "replace-values CSR column graph mismatch",
    "replace-values requires a genuinely changed operator",
)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return {"sha256": h.hexdigest(), "bytes": path.stat().st_size}


def snapshot_hashes(paths):
    return {name: sha256(path) for name, path in paths.items()}


def check_unchanged(paths, expected):
    for name, path in paths.items():
        if sha256(path) != expected[name]:
            raise ValueError("input/EXE changed during pair: " + name)


def read_steps(path):
    with path.open("r", newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ["step", "matrix", "rhs"]:
            raise ValueError("snapshot CSV header must be exactly step,matrix,rhs")
        rows = list(reader)
    if not (2 <= len(rows) <= 50):
        raise ValueError("need 2..50 consecutive snapshot rows")
    seen = set()
    steps = []
    for row in rows:
        label = row.get("step", "").strip()
        if not re.fullmatch(r"[A-Za-z0-9_-]{1,48}", label) or label in seen:
            raise ValueError("missing/duplicate/unsafe step identifier")
        seen.add(label)
        files = {}
        for key in ("matrix", "rhs"):
            raw = row.get(key)
            if raw is None or not raw.strip():
                raise ValueError("missing " + key + " for step " + label)
            p = Path(raw.strip()).expanduser()
            if not p.is_absolute():
                p = path.parent / p
            files[key] = p.resolve()
        steps.append({"step": label, **files})
    return steps


def one_line(log, marker):
    matched = [s for s in log.splitlines() if s.startswith(marker + " ")]
    if len(matched) != 1:
        raise ValueError("expected exactly one " + marker + " line")
    result = {}
    for part in matched[0].split()[1:]:
        if "=" in part:
            key, value = part.split("=", 1)
            if key in result:
                raise ValueError("duplicate " + marker + " field " + key)
            result[key] = value
    return result


def integer(fields, key):
    v = fields.get(key, "")
    if not re.fullmatch(r"-?[0-9]+", v):
        raise ValueError("missing/invalid integer " + key)
    return int(v)


def real(fields, key):
    try:
        v = float(fields[key])
    except (KeyError, ValueError) as exc:
        raise ValueError("missing/invalid finite number " + key) from exc
    if not math.isfinite(v):
        raise ValueError("nonfinite " + key)
    return v


def parse_solve(log, marker, tol):
    v = one_line(log, marker)
    rc, status = integer(v, "solve_rc"), integer(v, "status")
    iter_count = integer(v, "iter")
    residual = real(v, "true_relative")
    seconds = real(v, "wall_seconds")
    if iter_count < 0 or residual < 0 or seconds < 0:
        raise ValueError("invalid solver numerical metrics")
    return {
        "status": status, "true_relative": residual,
        "iter": iter_count, "wall_seconds": seconds,
        "ok": rc == 0 and status == 0 and residual <= tol,
    }


def parse_check(log, marker):
    v = one_line(log, marker)
    return {
        "orth_error": real(v, "orth_error"),
        "relation_error": real(v, "relation_error"),
        "valid": integer(v, "rc") == 0
        and 0 <= real(v, "orth_error") <= 1e-8
        and 0 <= real(v, "relation_error") <= 1e-8,
    }


def invoke(argv, cwd, env, log):
    with log.open("x", encoding="utf-8") as stream:
        result = subprocess.run(
            argv, cwd=str(cwd), env=env, stdout=stream,
            stderr=subprocess.STDOUT, check=False)
    return result.returncode


def read_log(path):
    return path.read_text(encoding="utf-8", errors="replace")


def inspect_preflight(exit_code, log):
    if exit_code == 0:
        fields = one_line(log, "REPLAY_CSR_COMPARE")
        if integer(fields, "compatible") != 1:
            raise ValueError("zero preflight exit without compatible=1")
        n, nnz = integer(fields, "n"), integer(fields, "nnz")
        changed = integer(fields, "changed_nnz")
        if n <= 0 or nnz <= 0 or not 1 <= changed <= nnz:
            raise ValueError("invalid CSR preflight counts")
        return "same_csr", changed, n
    fields = one_line(log, "REPLAY_CSR_COMPARE")
    if integer(fields, "compatible") != 0:
        raise ValueError("failed preflight without compatible=0")
    if "replace-values requires a genuinely changed operator" in log:
        return "unchanged", 0, None
    if any(item in log for item in STRUCTURE_ERRORS[:2]):
        return "changed_graph", 0, None
    # Bad matrix input, invalid numerics or format: do not hide by rebuild.
    raise ValueError("CSR preflight failed for unsupported reason")


def validate_header(fields, n, mode_count, tol):
    if integer(fields, "n") != n or integer(fields, "modes") != mode_count:
        raise ValueError("replay header row count or mode count mismatch")
    if (fields.get("method") != "fgcrodr" or
            integer(fields, "recycle_k") != 8 or
            integer(fields, "restart") != 180 or
            integer(fields, "maxiter") != 3000 or
            not math.isclose(real(fields, "tol"), tol, rel_tol=1e-12)):
        raise ValueError("replay solver configuration mismatch")


def evaluate_replay(log, exit_code, graph, tol, n, mode_count,
                    expected_changed):
    header = one_line(log, "REPLAY")
    validate_header(header, n, mode_count, tol)
    step = one_line(log, "REPLAY_NEWTON_STEP")
    if (integer(step, "step") != 2 or integer(step, "n") != n or
            integer(step, "modes") != mode_count or
            step.get("method") != "fgcrodr" or
            integer(step, "recycle_k") != 8):
        raise ValueError("second Newton step diagnostics mismatch")
    first = parse_solve(log, "REPLAY_RESULT", tol)
    second = parse_solve(log, "REPLAY_NEWTON_RESULT", tol)
    first_check = parse_check(log, "REPLAY_RECYCLE_CHECK")
    check = parse_check(log, "REPLAY_NEWTON_RECYCLE_CHECK")
    policy = one_line(log, "REPLAY_NEWTON_POLICY")
    matrix = one_line(log, "REPLAY_NEWTON_MATRIX")
    mode = matrix.get("mode", "")
    expected_mode = ("inplace_replace_values" if graph == "same_csr"
                     else "separate_object")
    if mode != expected_mode:
        raise ValueError("replay second matrix mode differs from preflight")
    p = {key: integer(policy, key) for key in
         ("cross_uses_delta", "keeps_delta",
          "transports_delta", "rejects_delta")}
    if any(x < 0 for x in p.values()):
        raise ValueError("negative policy deltas")
    if graph != "same_csr" and (p["cross_uses_delta"] != 0 or
                                p["keeps_delta"] != 0 or
                                p["transports_delta"] != 0):
        raise ValueError("stale basis reused across a distinct matrix object")
    if graph == "same_csr" and integer(matrix, "changed_nnz") != expected_changed:
        raise ValueError("in-place changed_nnz disagrees with CSR preflight")
    # A policy may legitimately reject an unstable transported basis.
    return first, second, check, p, (
        exit_code == 0 and first["ok"] and second["ok"]
        and first_check["valid"] and check["valid"])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe", type=Path, required=True)
    p.add_argument("--modes", type=Path, required=True)
    p.add_argument("--snapshots", type=Path, required=True,
                   help="Private LOCAL CSV with step,matrix,rhs; paths relative to CSV")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--tol", type=float, default=1e-5)
    p.add_argument("--fresh-control", action="store_true",
                   help="One additional independent A2 solve per pair")
    p.add_argument("--continue-on-failure", action="store_true")
    args = p.parse_args()
    if not math.isfinite(args.tol) or args.tol != 1e-5:
        p.error("current replay supports only its hardcoded 1e-5 tolerance")
    try:
        steps = read_steps(args.snapshots.resolve())
        n, modes = read_modes(args.modes)
        paths = {"exe": args.exe.resolve(), "modes": args.modes.resolve(),
                 "snapshots_csv": args.snapshots.resolve()}
        for step in steps:
            paths[step["step"] + "_A"] = step["matrix"]
            paths[step["step"] + "_b"] = step["rhs"]
        if len({v for v in paths.values()}) != len(paths):
            raise ValueError("input paths must be distinct (no silent duplicated snapshots)")
        for name, path in paths.items():
            if not path.is_file() or path.stat().st_size <= 0:
                raise ValueError("missing/empty local input: " + name)
        baseline_hashes = snapshot_hashes(paths)
        out = args.out_dir.resolve()
        if out.exists():
            raise ValueError("output directory already exists; choose a NEW local directory")
        out.mkdir(parents=True, mode=0o700)
        os.chmod(out, 0o700)
        manifest = {
            "schema": "mfamg-v8-pairwise-contact-sequence-v1",
            "step_ids": [s["step"] for s in steps],
            "hashes": baseline_hashes,
            "mode_count": len(modes),
            "n_from_modes": n,
            "fixed_solver_environment": ENV,
            "fresh_control": args.fresh_control,
            "warning": "Fresh solver per adjacent PAIR. NOT a continuous Newton trajectory. Private LOCAL only.",
        }
        (out / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        env = {k: v for k, v in os.environ.items() if not k.startswith("LIS_")}
        env.update(ENV)
        entries = []
        with (out / "transitions.csv").open(
                "x", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, FIELDS)
            writer.writeheader()
            stream.flush()
            for i, (left, right) in enumerate(zip(steps, steps[1:]), 1):
                prefix = "pair%02d" % i
                logs = [prefix + "_csr.log", prefix + "_replay.log",
                        prefix + "_fresh.log" if args.fresh_control else ""]
                row = dict.fromkeys(FIELDS, "")
                row.update(pair=i, from_step=left["step"], to_step=right["step"],
                           preflight_log=logs[0], replay_log=logs[1],
                           fresh_log=logs[2], pair_pass=0)
                active_paths = {
                    "exe": paths["exe"], "modes": paths["modes"],
                    "snapshots_csv": paths["snapshots_csv"],
                    left["step"] + "_A": left["matrix"],
                    left["step"] + "_b": left["rhs"],
                    right["step"] + "_A": right["matrix"],
                    right["step"] + "_b": right["rhs"],
                }
                try:
                    check_unchanged(active_paths, baseline_hashes)
                    preflight_cmd = [str(paths["exe"]), "--compare-csr",
                                     str(left["matrix"]), str(right["matrix"])]
                    preflight_exit = invoke(preflight_cmd, out, env, out / logs[0])
                    graph, changed, n_pre = inspect_preflight(
                        preflight_exit, read_log(out / logs[0]))
                    if n_pre is not None and n_pre != n:
                        raise ValueError("CSR preflight rows disagree with modes")
                    row.update(graph=graph, changed_nnz=changed)
                    if graph == "same_csr":
                        new_mode = "inplace_replace_values"
                    else:
                        new_mode = "separate_object"
                    row["mode"] = new_mode
                    pair_env = env.copy()
                    pair_env["LIS_MFAMG_REPLAY_NEXT_A"] = str(right["matrix"])
                    pair_env["LIS_MFAMG_REPLAY_NEXT_B"] = str(right["rhs"])
                    if graph == "same_csr":
                        pair_env["LIS_MFAMG_REPLAY_NEXT_REPLACE_VALUES"] = "1"
                    argv = [str(paths["exe"]), str(left["matrix"]),
                            str(left["rhs"]), str(paths["modes"])]
                    code = invoke(argv, left["matrix"].parent, pair_env, out / logs[1])
                    row["exit"] = code
                    check_unchanged(active_paths, baseline_hashes)
                    first, second, check, policy, ok = evaluate_replay(
                        read_log(out / logs[1]), code, graph, args.tol, n,
                        len(modes), changed)
                    row.update(first_status=first["status"],
                               first_true_relative=first["true_relative"],
                               second_status=second["status"],
                               second_true_relative=second["true_relative"],
                               second_iter=second["iter"],
                               second_wall_seconds=second["wall_seconds"],
                               orth_error=check["orth_error"],
                               relation_error=check["relation_error"], **policy)
                    if args.fresh_control:
                        fresh_env = env.copy()
                        code2 = invoke(
                            [str(paths["exe"]), str(right["matrix"]),
                             str(right["rhs"]), str(paths["modes"])],
                            right["matrix"].parent, fresh_env, out / logs[2])
                        row["fresh_exit"] = code2
                        check_unchanged(active_paths, baseline_hashes)
                        fresh_header = one_line(read_log(out / logs[2]), "REPLAY")
                        validate_header(fresh_header, n, len(modes), args.tol)
                        fresh_check = parse_check(
                            read_log(out / logs[2]), "REPLAY_RECYCLE_CHECK")
                        fresh = parse_solve(read_log(out / logs[2]), "REPLAY_RESULT", args.tol)
                        row.update(fresh_status=fresh["status"],
                                   fresh_true_relative=fresh["true_relative"],
                                   fresh_iter=fresh["iter"],
                                   fresh_wall_seconds=fresh["wall_seconds"])
                        ok = ok and code2 == 0 and fresh["ok"] and fresh_check["valid"]
                    row["pair_pass"] = int(ok)
                    if not ok:
                        row["error"] = "one or more solves/algebraic checks did not converge"
                except (OSError, ValueError) as exc:
                    row["error"] = str(exc)
                entries.append(row)
                writer.writerow(row)
                stream.flush()
                status = "PASS" if row["pair_pass"] else "FAIL"
                print("V8_SEQUENCE_PAIR=%s from=%s to=%s graph=%s result=%s" %
                      (i, left["step"], right["step"], row["graph"], status),
                      flush=True)
                summary = {
                    "assessment": "INCOMPLETE",
                    "completed_pairs": len(entries),
                    "expected_pairs": len(steps) - 1,
                    "passed_pairs": sum(int(z["pair_pass"]) for z in entries),
                    "warning": manifest["warning"],
                    "transitions_csv": "transitions.csv",
                }
                (out / "summary.json").write_text(
                    json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
                if not row["pair_pass"] and not args.continue_on_failure:
                    break
        total = len(steps) - 1
        passed = sum(int(z["pair_pass"]) for z in entries)
        assessment = ("PAIRWISE_PASS" if len(entries) == total and passed == total
                      else "INCONCLUSIVE_OR_FAILED")
        final = {
            "assessment": assessment, "completed_pairs": len(entries),
            "expected_pairs": total, "passed_pairs": passed,
            "warning": manifest["warning"],
            "transitions_csv": "transitions.csv",
        }
        (out / "summary.json").write_text(
            json.dumps(final, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print("V8_SEQUENCE_ASSESSMENT=" + assessment)
        print("V8_SEQUENCE_SUMMARY=" + str(out / "summary.json"))
        return 0 if assessment == "PAIRWISE_PASS" else 1
    except (OSError, ValueError) as exc:
        print("V8_SEQUENCE_ERROR=" + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
