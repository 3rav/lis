#!/usr/bin/env python3
"""Stage a PRIVATE, isolated source-only CalculiX research copy.

This does not rebuild or instrument CalculiX. The original source trees,
installed EXE, MKL directory, INP and any private solver matrices remain
untouched. No source bytes, paths, hashes or private logs are printed.
"""
import argparse
import hashlib
from pathlib import Path
import shutil
import sys


EXTENSIONS = {".c", ".h", ".f", ".for", ".f90", ".inc"}


def fail(reason):
    raise ValueError(reason)


def allowed_file(path):
    return path.suffix.lower() in EXTENSIONS or path.name.lower().startswith("makefile")


def collect(folder):
    if not folder.is_dir() or folder.is_symlink():
        fail("invalid source directory")
    result = {}
    for entry in folder.iterdir():
        if entry.is_symlink():
            # Preserve neither links nor a possible escape outside the tree.
            continue
        if entry.is_file() and allowed_file(entry):
            result[entry.name] = entry
    if not result or "nonlingeo.c" not in result or "CalculiX.h" not in result:
        fail("missing required CalculiX sources")
    if not any(k in result for k in ("ccx.c", "ccx_2.22.c", "ccx_2.23.c")) \
            and not any(x.startswith("ccx_") and x.endswith(".c") for x in result):
        fail("missing main CalculiX source")
    if not any(x in result for x in ("mafillsmmain.c", "mafillsmmain.f")):
        fail("missing assembly bridge source")
    return result


def digest(file):
    d = hashlib.sha256()
    with file.open("rb") as stream:
        for piece in iter(lambda: stream.read(1024 * 1024), b""):
            d.update(piece)
    return d.digest()


def source_parent(dirpath):
    return dirpath.parent if dirpath.name.lower() == "src" else dirpath


def makefiles_in_root(source_dir):
    parent = source_parent(source_dir)
    if parent == source_dir:
        return {}
    return {p.name: p for p in parent.iterdir()
            if p.is_file() and not p.is_symlink()
            and p.name.lower().startswith("makefile")}


def stage(args):
    selected = Path(args.source_dir).expanduser().resolve()
    peer = Path(args.peer_dir).expanduser().resolve()
    out = Path(args.out_dir).expanduser().resolve()
    install = Path(args.install_dir).expanduser().resolve()
    exe = Path(args.installed_exe).expanduser().resolve()
    mkl = Path(args.mkl_dir).expanduser().resolve()
    if selected == peer or not exe.is_file() or not exe.stat().st_size:
        fail("source pair or installed binary missing")
    if not install.is_dir() or not mkl.is_dir():
        fail("missing install or MKL dependency directory")
    if not selected.is_relative_to(install) or not peer.is_relative_to(install):
        fail("candidate sources must be under the specified installation")
    if not exe.is_relative_to(install) or not mkl.is_relative_to(install):
        fail("binary/MKL must be inside the specified installation")
    if out.exists() or out.is_symlink():
        fail("output already exists")
    if out == Path(out.anchor) or out.is_relative_to(install):
        fail("output must be separate from installation")
    if not out.parent.is_dir():
        fail("output parent must exist, no implicit directory creation")
    if out.is_relative_to(selected) or out.is_relative_to(peer):
        fail("output cannot be under source")
    left, right = collect(selected), collect(peer)
    if set(left) != set(right):
        fail("source sets differ")
    # Hash every file before copy; never log private digests or paths.
    hashes = {name: digest(path) for name, path in left.items()}
    if any(digest(right[name]) != hashes[name] for name in left):
        fail("candidate source contents differ")
    if len(left) < 50:
        fail("source tree is suspiciously small")
    configs = makefiles_in_root(selected)
    # Prevent recursive copying and refuse unexpected output state.
    out.mkdir()
    try:
        (out / "src").mkdir()
        for name, src in left.items():
            dest = out / "src" / name
            with dest.open("xb") as f_out, src.open("rb") as f_in:
                shutil.copyfileobj(f_in, f_out, length=1024 * 1024)
            if digest(dest) != hashes[name]:
                fail("staged source integrity mismatch")
        for name, src in configs.items():
            if (out / name).exists():
                fail("duplicate root makefile")
            with (out / name).open("xb") as f_out, src.open("rb") as f_in:
                shutil.copyfileobj(f_in, f_out, length=1024 * 1024)
        if (set(collect(selected)) != set(left) or
                any(digest(left[name]) != hashes[name] for name in left)):
            fail("source changed during staging")
        if any(digest(right[name]) != hashes[name] for name in right):
            fail("peer source changed during staging")
    except Exception:
        # Remove only the new directory created by this invocation; never
        # remove user-provided pre-existing directories or the source itself.
        shutil.rmtree(out)
        raise
    print("CCX_STAGE=ISOLATED_SOURCE_COPY_READY")
    print("CCX_STAGE_SOURCE_CANDIDATES_MATCH=1")
    print("CCX_STAGE_SOURCE_FILES=%d" % len(left))
    print("CCX_STAGE_ROOT_MAKEFILES=%d" % len(configs))
    print("CCX_STAGE_ORIGINAL_SOURCES_UNCHANGED=1")
    print("CCX_STAGE_INSTALLED_EXE_UNCHANGED=1")
    print("CCX_STAGE_MKL_NOT_COPIED=1")
    print("CCX_STAGE_MODEL_NOT_COPIED=1")
    print("CCX_STAGE_BUILD=NOT_ATTEMPTED")
    print("CCX_STAGE_EXPORT=NOT_INSTRUMENTED")
    print("CCX_STAGE_BINARY_CORRESPONDENCE=NOT_PROVEN")
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", required=True)
    parser.add_argument("--peer-dir", required=True)
    parser.add_argument("--install-dir", required=True)
    parser.add_argument("--installed-exe", required=True)
    parser.add_argument("--mkl-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    try:
        return stage(args)
    except (ValueError, OSError, OverflowError, UnicodeError) as exc:
        print("CCX_STAGE=BLOCKED", file=sys.stderr)
        print("CCX_STAGE_REASON=" + type(exc).__name__, file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
