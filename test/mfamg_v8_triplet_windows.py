#!/usr/bin/env python3
"""LOCAL sliding three-step FGCRO-DR windows, one LIS_SOLVER per triplet.

A1 -> A2 -> A3 truly share a solver WITHIN each window. For a longer
CSV the next window starts in a fresh process; this is not a production
Newton trajectory. All real contact matrices/RHS/logs remain private.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import sys

from mfamg_mode_audit import read_modes
from mfamg_replay_v8_repeats import ENV
from mfamg_v8_snapshot_sequence import (
    read_steps, snapshot_hashes, check_unchanged, invoke, read_log,
    inspect_preflight, one_line, integer, real, parse_solve,
    parse_check, validate_header, evaluate_replay,
)

FIELDS = ("window", "from_step", "mid_step", "to_step", "graph12", "graph23",
          "changed12", "changed23", "mode2", "mode3", "exit",
          "first_status", "second_status", "third_status",
          "second_true_relative", "third_true_relative", "third_iter",
          "third_wall_seconds", "third_cross_uses_delta",
          "third_keeps_delta", "third_transports_delta", "third_rejects_delta",
          "third_orth_error", "third_relation_error", "fresh_exit",
          "fresh_status", "fresh_true_relative", "pass", "error",
          "csr12_log", "csr23_log", "replay_log", "fresh_log")


def check_triplet(log, exit_code, graph12, changed12, graph23, changed23,
                  n, mode_count, tolerance):
    first, second, check2, policy2, first_two = evaluate_replay(
        log, 0, graph12, tolerance, n, mode_count, changed12)
    header = one_line(log, "REPLAY_SEQUENCE_STEP")
    if (integer(header, "step") != 3 or integer(header, "n") != n or
            integer(header, "modes") != mode_count or
            integer(header, "recycle_k") != 8 or
            header.get("method") != "fgcrodr"):
        raise ValueError("third step configuration disagrees with replay input")
    matrix = one_line(log, "REPLAY_SEQUENCE_MATRIX")
    expected_mode = ("inplace_replace_values" if graph23 == "same_csr"
                     else "separate_object")
    if integer(matrix, "step") != 3 or matrix.get("mode") != expected_mode:
        raise ValueError("third-step matrix mode disagrees with CSR preflight")
    if graph23 == "same_csr" and integer(matrix, "changed_nnz") != changed23:
        raise ValueError("third-step CSR changed_nnz differs from preflight")
    status = parse_solve(log, "REPLAY_SEQUENCE_RESULT", tolerance)
    check = parse_check(log, "REPLAY_SEQUENCE_RECYCLE_CHECK")
    policy = one_line(log, "REPLAY_SEQUENCE_POLICY")
    if integer(policy, "step") != 3:
        raise ValueError("third policy step marker missing")
    values = {name: integer(policy, name) for name in
              ("cross_uses_delta", "keeps_delta",
               "transports_delta", "rejects_delta")}
    if any(v < 0 for v in values.values()):
        raise ValueError("negative third policy counters")
    if graph23 != "same_csr" and any(values[key] != 0 for key in
                                      ("cross_uses_delta", "keeps_delta",
                                       "transports_delta")):
        raise ValueError("third-step stale recycling accepted for new graph")
    ok = exit_code == 0 and first_two and status["ok"] and check["valid"]
    return first, second, status, check, values, ok


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe", type=Path, required=True)
    p.add_argument("--modes", type=Path, required=True)
    p.add_argument("--snapshots", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--fresh-third", action="store_true",
                   help="Compare each A3 with an independent fresh-solver control")
    p.add_argument("--continue-on-failure", action="store_true")
    p.add_argument("--max-windows", type=int, default=0,
                   help="0=all; otherwise run only first N adjacent triplets")
    p.add_argument("--fine-agg-factor", type=int, choices=(2, 4), default=4,
                   help="Common local MFAMG factor (2=116 springs, 4=31 springs)")
    args = p.parse_args()
    if not 0 <= args.max_windows <= 48:
        p.error("max-windows must be 0..48")
    try:
        steps = read_steps(args.snapshots.resolve())
        if len(steps) < 3:
            raise ValueError("three-step probe needs at least 3 consecutive snapshots")
        n, modes = read_modes(args.modes)
        paths = {"exe": args.exe.resolve(), "modes": args.modes.resolve(),
                 "snapshots_csv": args.snapshots.resolve()}
        for step in steps:
            paths[step["step"] + "_A"] = step["matrix"]
            paths[step["step"] + "_b"] = step["rhs"]
        if len(set(paths.values())) != len(paths):
            raise ValueError("input paths must be unique")
        for label, path in paths.items():
            if not path.is_file() or path.stat().st_size <= 0:
                raise ValueError("missing/empty local input: " + label)
        expected = snapshot_hashes(paths)
        out = args.out_dir.resolve()
        if out.exists():
            raise ValueError("output directory already exists; choose a NEW local directory")
        # POSIX 0700 is not a reliable permissions model on Windows/MSYS2.
        # Respect Windows ACLs instead of forcing chmod on the local result dir.
        out.mkdir(parents=True, mode=0o700)
        if os.name != "nt":
            os.chmod(out, 0o700)
        max_windows = len(steps) - 2
        if args.max_windows:
            max_windows = min(max_windows, args.max_windows)
        fixed = dict(ENV)
        fixed["LIS_MFAMG_FINE_AGG_FACTOR"] = str(args.fine_agg_factor)
        manifest = {
            "schema": "mfamg-v8-persistent-triplet-v1",
            "step_ids": [s["step"] for s in steps],
            "hashes": expected,
            "modes": len(modes),
            "n": n,
            "windows_requested": max_windows,
            "fresh_third": args.fresh_third,
            "solver_environment": fixed,
            "limitation": "One LIS_SOLVER persists for EXACTLY 3 consecutive linear systems per window. Each window is a NEW process. Synthetic PASS does not certify contact Newton.",
        }
        (out / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        env = {k: v for k, v in os.environ.items() if not k.startswith("LIS_")}
        env.update(fixed)
        rows = []
        with (out / "triplets.csv").open(
                "x", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, FIELDS)
            writer.writeheader()
            stream.flush()
            for j in range(max_windows):
                first, second, third = steps[j:j + 3]
                prefix = "window%02d" % (j + 1)
                logs = [prefix + "_csr12.log", prefix + "_csr23.log",
                        prefix + "_replay.log",
                        prefix + "_fresh.log" if args.fresh_third else ""]
                row = dict.fromkeys(FIELDS, "")
                row.update(window=j + 1, from_step=first["step"],
                           mid_step=second["step"], to_step=third["step"],
                           csr12_log=logs[0], csr23_log=logs[1],
                           replay_log=logs[2], fresh_log=logs[3],
                           **{"pass": 0})
                active = {"exe": paths["exe"], "modes": paths["modes"],
                          "snapshots_csv": paths["snapshots_csv"]}
                for step in (first, second, third):
                    active[step["step"] + "_A"] = step["matrix"]
                    active[step["step"] + "_b"] = step["rhs"]
                try:
                    check_unchanged(active, expected)
                    graph = []
                    changed = []
                    for k, (left, right) in enumerate(
                            ((first, second), (second, third))):
                        cmd = [str(paths["exe"]), "--compare-csr",
                               str(left["matrix"]), str(right["matrix"])]
                        code = invoke(cmd, out, env, out / logs[k])
                        g, nchanged, n_detected = inspect_preflight(
                            code, read_log(out / logs[k]))
                        if n_detected is not None and n_detected != n:
                            raise ValueError("CSR rows disagree with mode rows")
                        graph.append(g)
                        changed.append(nchanged)
                    row.update(graph12=graph[0], graph23=graph[1],
                               changed12=changed[0], changed23=changed[1])
                    pair_env = env.copy()
                    pair_env.update({
                        "LIS_MFAMG_REPLAY_NEXT_A": str(second["matrix"]),
                        "LIS_MFAMG_REPLAY_NEXT_B": str(second["rhs"]),
                        "LIS_MFAMG_REPLAY_THIRD_A": str(third["matrix"]),
                        "LIS_MFAMG_REPLAY_THIRD_B": str(third["rhs"]),
                    })
                    if graph[0] == "same_csr":
                        pair_env["LIS_MFAMG_REPLAY_NEXT_REPLACE_VALUES"] = "1"
                    if graph[1] == "same_csr":
                        pair_env["LIS_MFAMG_REPLAY_THIRD_REPLACE_VALUES"] = "1"
                    row.update(
                        mode2=("inplace_replace_values" if graph[0] == "same_csr"
                               else "separate_object"),
                        mode3=("inplace_replace_values" if graph[1] == "same_csr"
                               else "separate_object"))
                    exit_code = invoke(
                        [str(paths["exe"]), str(first["matrix"]),
                         str(first["rhs"]), str(paths["modes"])],
                        first["matrix"].parent, pair_env, out / logs[2])
                    row["exit"] = exit_code
                    check_unchanged(active, expected)
                    a, b, c, diag3, policy, ok = check_triplet(
                        read_log(out / logs[2]), exit_code,
                        graph[0], changed[0], graph[1], changed[1],
                        n, len(modes), 1e-5)
                    row.update(
                        first_status=a["status"], second_status=b["status"],
                        third_status=c["status"],
                        second_true_relative=b["true_relative"],
                        third_true_relative=c["true_relative"],
                        third_iter=c["iter"],
                        third_wall_seconds=c["wall_seconds"],
                        third_cross_uses_delta=policy["cross_uses_delta"],
                        third_keeps_delta=policy["keeps_delta"],
                        third_transports_delta=policy["transports_delta"],
                        third_rejects_delta=policy["rejects_delta"],
                        third_orth_error=diag3["orth_error"],
                        third_relation_error=diag3["relation_error"])
                    if args.fresh_third:
                        fresh_code = invoke(
                            [str(paths["exe"]), str(third["matrix"]),
                             str(third["rhs"]), str(paths["modes"])],
                            third["matrix"].parent, env, out / logs[3])
                        row["fresh_exit"] = fresh_code
                        check_unchanged(active, expected)
                        fresh_log = read_log(out / logs[3])
                        validate_header(one_line(fresh_log, "REPLAY"),
                                        n, len(modes), 1e-5)
                        fresh = parse_solve(fresh_log, "REPLAY_RESULT", 1e-5)
                        valid = parse_check(fresh_log, "REPLAY_RECYCLE_CHECK")
                        row.update(fresh_status=fresh["status"],
                                   fresh_true_relative=fresh["true_relative"])
                        ok = ok and fresh_code == 0 and fresh["ok"] and valid["valid"]
                    row["pass"] = int(ok)
                    if not ok:
                        row["error"] = "convergence or recycle validation failed"
                except (OSError, ValueError) as exc:
                    row["error"] = str(exc)
                rows.append(row)
                writer.writerow(row)
                stream.flush()
                print("V8_TRIPLET window=%d %s/%s/%s graph12=%s graph23=%s pass=%d" %
                      (j + 1, first["step"], second["step"], third["step"],
                       row["graph12"], row["graph23"], row["pass"]), flush=True)
                interim = {
                    "assessment": "INCOMPLETE",
                    "completed_windows": len(rows),
                    "requested_windows": max_windows,
                    "passed_windows": sum(int(r["pass"]) for r in rows),
                    "limitation": manifest["limitation"],
                }
                (out / "summary.json").write_text(
                    json.dumps(interim, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")
                if not row["pass"] and not args.continue_on_failure:
                    break
        succeeded = sum(int(row["pass"]) for row in rows)
        result = "PERSISTENT_TRIPLETS_PASS" if (
            len(rows) == max_windows and succeeded == max_windows) else "INCONCLUSIVE_OR_FAILED"
        summary = {
            "assessment": result,
            "requested_windows": max_windows,
            "completed_windows": len(rows),
            "passed_windows": succeeded,
            "limitation": manifest["limitation"],
            "triplets_csv": "triplets.csv",
        }
        (out / "summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8")
        print("V8_TRIPLET_ASSESSMENT=" + result)
        print("V8_TRIPLET_SUMMARY=" + str(out / "summary.json"))
        return 0 if result == "PERSISTENT_TRIPLETS_PASS" else 1
    except (OSError, ValueError) as exc:
        print("V8_TRIPLET_ERROR=" + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
