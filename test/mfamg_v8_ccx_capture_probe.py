#!/usr/bin/env python3
"""Private local read-only readiness probe for a CalculiX Newton capture.

Checks source signatures and top-level INP directives only. Does not locate
the precise runtime solver dispatch, instrument CalculiX, infer chronology,
solve a model, or publish contents/paths/hashes. Source and model stay local.
"""
import argparse
import hashlib
from pathlib import Path
import sys


def digest(path):
    d = hashlib.sha256()
    with path.open("rb") as stream:
        for data in iter(lambda: stream.read(1024 * 1024), b""):
            d.update(data)
    return d.digest()


def inspect(root, model):
    root = Path(root).expanduser().resolve()
    model = Path(model).expanduser().resolve()
    candidates = (root / "src" / "nonlingeo.c", root / "nonlingeo.c")
    found = [p for p in candidates if p.is_file() and p.stat().st_size > 0]
    if len(found) != 1 or not model.is_file() or not model.stat().st_size:
        raise ValueError("missing or ambiguous source/model")
    source = found[0]
    if source == model:
        raise ValueError("source and model are not distinct")
    before = (digest(source), digest(model))
    code = source.read_text(encoding="utf-8", errors="replace")
    signatures = {
        "assembly": "mafillsmmain(" in code,
        "newton_counters": "iinc" in code and "iit" in code,
        "direct_solver": ("spooles(" in code or "pardiso_main(" in code or
                          "pastix_main(" in code),
        "iterative_solver": "preiter(" in code or "predgmres_struct" in code,
        "matrixstorage": "matrixstorage(" in code,
    }
    if not signatures["assembly"] or not signatures["newton_counters"]:
        raise ValueError("source variant requires manual inspection")
    flags = {
        "step": False, "static": False, "dynamic": False,
        "contact": False, "spring": False, "include": False,
        "solver_keyword": False,  # legacy flag; not the *STATIC solver setting
    }
    static_count = 0
    explicit_static_solvers = []
    with model.open("r", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            stripped = line.lstrip()
            if not stripped.startswith("*") or stripped.startswith("**"):
                continue
            directive = stripped.split(",", 1)[0].strip().upper()
            if directive == "*STEP":
                flags["step"] = True
            elif directive == "*STATIC":
                flags["static"] = True
                static_count += 1
                for parameter in stripped.split(",")[1:]:
                    key, sep, value = parameter.partition("=")
                    if sep and key.strip().upper() == "SOLVER":
                        # Only fixed public solver categories are emitted.
                        solver = value.strip().upper()
                        explicit_static_solvers.append(
                            solver if solver in (
                                "PARDISO", "SPOOLES", "PASTIX",
                                "ITERATIVESCALING", "ITERATIVECHOLESKY",
                                "SGI", "TAUCS"
                            ) else "OTHER"
                        )
            elif directive == "*DYNAMIC":
                flags["dynamic"] = True
            elif directive.startswith("*CONTACT") or directive == "*SURFACE INTERACTION":
                flags["contact"] = True
            elif directive.startswith("*SPRING"):
                flags["spring"] = True
            elif directive == "*INCLUDE":
                flags["include"] = True
            elif directive == "*SOLVER":
                flags["solver_keyword"] = True
    if not static_count:
        static_solver = "NO_TOP_LEVEL_STATIC"
    elif not explicit_static_solvers:
        static_solver = "UNSPECIFIED_TOP_LEVEL"
    elif (len(explicit_static_solvers) != static_count or
          len(set(explicit_static_solvers)) != 1):
        static_solver = "MIXED_OR_PARTIAL_TOP_LEVEL"
    else:
        static_solver = explicit_static_solvers[0]
    if (digest(source), digest(model)) != before:
        raise ValueError("source/model changed during probe")
    return signatures, flags, static_solver


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", required=True)
    parser.add_argument("--input-file", required=True)
    args = parser.parse_args()
    try:
        signatures, flags, static_solver = inspect(args.source_root, args.input_file)
    except (ValueError, OSError, UnicodeError, OverflowError) as exc:
        print("CCX_CAPTURE_PROBE=BLOCKED", file=sys.stderr)
        print("CCX_CAPTURE_PROBE_REASON=" + type(exc).__name__, file=sys.stderr)
        return 2
    print("CCX_CAPTURE_PROBE=SOURCE_AND_MODEL_READABLE")
    for name in ("assembly", "newton_counters", "direct_solver",
                 "iterative_solver", "matrixstorage"):
        print("CCX_CAPTURE_SOURCE_%s=%d" %
              (name.upper(), int(signatures[name])))
    for name in ("step", "static", "dynamic", "contact", "spring",
                 "include", "solver_keyword"):
        print("CCX_CAPTURE_INP_%s=%d" %
              (name.upper(), int(flags[name])))
    print("CCX_CAPTURE_INP_STATIC_SOLVER=" + static_solver)
    print("CCX_CAPTURE_INP_SCOPE=TOP_LEVEL_ONLY")
    print("CCX_CAPTURE_INP_STATIC_SOLVER_IS_RUNTIME_PROOF=0")
    print("CCX_CAPTURE_SOLVER_HOOK=NEEDS_BRANCH_SPECIFIC_REVIEW")
    print("CCX_CAPTURE_EXPORT_ENABLED=0")
    print("CCX_CAPTURE_LOCAL_INPUTS_UNCHANGED=1")
    print("CCX_CAPTURE_PROVENANCE=NOT_YET_ESTABLISHED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
