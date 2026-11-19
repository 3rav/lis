#!/usr/bin/env python3
"""Independently compare byte-exact ILUTP factor snapshots from pinned LIS builds.

This does NOT modify the production source in either commit. Both checkouts
receive exactly the same diagnostic-only edits to test/benchilutp.c, then build
that benchmark against their own unchanged LIS libraries.
"""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys

from bench_ilutp_topk_public import FIELDS, VARIANTS, generate_inputs

C_HELPER = r'''
/* Optional post-setup binary snapshot. Never included in setup_cpu_s. */
static int snapshot_write(FILE *f,const void *ptr,size_t size,size_t count)
{
    return count==0 || fwrite(ptr,size,count,f)==count;
}
static int snapshot_factors(LIS_PRECON pc,LIS_INT n,const char *variant)
{
    const char *prefix=getenv("LIS_ILUTP_FACTOR_PREFIX");
    static const char magic[]="LIS-factor-raw-v1";
    char path[4096];
    FILE *f;
    LIS_INT i,m;
    int kind,has_perm;
    if(!prefix || !*prefix) return 0;
    if(snprintf(path,sizeof(path),"%s.%s.bin",prefix,variant)<0 ||
       strlen(path)>=sizeof(path)-1) return 1;
    f=fopen(path,"wb");
    if(!f)return 1;
    if(!snapshot_write(f,magic,1,sizeof(magic)) ||
       !snapshot_write(f,&n,sizeof(n),1))goto fail;
    for(kind=0;kind<2;kind++)
    {
        LIS_MATRIX_ILU factor=kind ? pc->U : pc->L;
        if(!factor)goto fail;
        for(i=0;i<n;i++)
        {
            m=factor->nnz[i];
            if(m<0 || (m>0 && (!factor->index[i] || !factor->value[i])))
                goto fail;
            if(!snapshot_write(f,&m,sizeof(m),1) ||
               !snapshot_write(f,factor->index[i],sizeof(LIS_INT),(size_t)m) ||
               !snapshot_write(f,factor->value[i],sizeof(LIS_SCALAR),(size_t)m))
                goto fail;
        }
    }
    if(!pc->D || !pc->D->value ||
       !snapshot_write(f,pc->D->value,sizeof(LIS_SCALAR),(size_t)n))
        goto fail;
    has_perm=(pc->ilutp_colperm!=NULL);
    if(!snapshot_write(f,&has_perm,sizeof(has_perm),1) ||
       (has_perm && !snapshot_write(f,pc->ilutp_colperm,
                                    sizeof(LIS_INT),(size_t)n)) ||
       !snapshot_write(f,&pc->ilutp_swap_count,
                       sizeof(pc->ilutp_swap_count),1))
        goto fail;
    return fclose(f)!=0;
fail:
    fclose(f);
    return 1;
}
'''

MARKER = "static int run_case("
CALL = "    err=lis_vector_set_all((LIS_SCALAR)0.0,x);"
INJECT_CALL = '''    if(snapshot_factors(pc,A->n,name))
    {
        fprintf(stderr,"Cannot serialize factor snapshot for %s\\n",name);
        exit(92);
    }
'''

CASES = (
    ("banded", 1200, 0, False),
    ("selection", 1200, 8, False),
    ("pivoted", 1200, 8, True),
)
RATES = (1, 2, 4, 8)
REFS = {
    "baseline": "3d027903ab079158d99bfc91f8ab6867c212dc20",
    "candidate": "78e694a82fe15575b22df11977b7277617e9b039",
}


def instrument(root, source):
    """Make the same controlled diagnostic change in each clean checkout."""
    text = source.decode("utf-8")
    if text.count(MARKER) != 1 or text.count(CALL) != 1:
        raise RuntimeError("native benchmark source no longer matches patch anchors")
    patched = text.replace(MARKER, C_HELPER + "\n" + MARKER, 1)
    patched = patched.replace(CALL, INJECT_CALL + CALL, 1)
    (root / "test" / "benchilutp.c").write_text(patched, encoding="utf-8")


def run_one(root, a, b, rate, prefix):
    binpath = root / "test" / "benchilutp"
    args = [str(binpath), str(a), str(b), "--rate", str(rate),
            "--maxiter", "220", "--restart", "50", "--rtol", "1e-8",
            "--drop", "0", "--pivot-tol", "1e-9", "--rcm", "0",
            "--row-scale", "0", "--equil", "0"]
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE",
               OPENBLAS_NUM_THREADS="1", LC_ALL="C",
               LIS_ILUTP_FACTOR_PREFIX=str(prefix))
    env["LD_LIBRARY_PATH"] = (str(root / "src" / ".libs") + os.pathsep +
                              env.get("LD_LIBRARY_PATH", ""))
    proc = subprocess.run(args, capture_output=True, text=True, env=env)
    prefix.with_suffix(".stdout.log").write_text(proc.stdout, encoding="utf-8")
    prefix.with_suffix(".stderr.log").write_text(proc.stderr, encoding="utf-8")
    if proc.returncode:
        raise RuntimeError(f"{prefix}: native benchmark exited {proc.returncode}")
    lines = proc.stdout.splitlines(keepends=True)
    try:
        i = next(k for k, line in enumerate(lines)
                 if line.strip() == ",".join(FIELDS))
    except StopIteration:
        raise RuntimeError(f"{prefix}: missing native CSV header")
    reader = csv.DictReader(io.StringIO("".join(lines[i:])))
    rows = list(reader)
    if reader.fieldnames != list(FIELDS) or [r["variant"] for r in rows] != list(VARIANTS):
        raise RuntimeError(f"{prefix}: native benchmark CSV schema mismatch")
    return {r["variant"]: r for r in rows}


def file_sha256(path):
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for buf in iter(lambda: stream.read(1 << 20), b""):
            sha.update(buf)
    return sha.hexdigest()


def equal_bytes(a, b):
    if a.stat().st_size != b.stat().st_size:
        return False
    with a.open("rb") as left, b.open("rb") as right:
        while True:
            x, y = left.read(1 << 20), right.read(1 << 20)
            if x != y:
                return False
            if not x:
                return True


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--baseline-root", type=Path, required=True)
    ap.add_argument("--candidate-root", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    roots = {"baseline": args.baseline_root.resolve(),
             "candidate": args.candidate_root.resolve()}
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    native = {k: (root / "test" / "benchilutp.c").read_bytes()
              for k, root in roots.items()}
    if native["baseline"] != native["candidate"]:
        ap.error("pinned revisions have different native benchmark sources")
    native_sha = hashlib.sha256(native["baseline"]).hexdigest()
    for label, root in roots.items():
        instrument(root, native[label])
        subprocess.run(["make", "-C", str(root / "test"), "benchilutp"],
                       check=True, stdout=(out / f"{label}-build.stdout.log").open("w"),
                       stderr=(out / f"{label}-build.stderr.log").open("w"))
    report, issues, inputs = [], [], {}
    for case, n, far, swapped in CASES:
        folder = out / case
        folder.mkdir(parents=True, exist_ok=True)
        a, b, input_manifest = generate_inputs(n, far, 20261119,
                                               folder / "inputs", swapped)
        inputs[case] = input_manifest
        for rate in RATES:
            snapshots = {}
            results = {}
            for label, root in roots.items():
                prefix = folder / f"rate{rate}-{label}"
                results[label] = run_one(root, a, b, rate, prefix)
                snapshots[label] = prefix
            for variant in VARIANTS:
                left, right = (results[name][variant] for name in
                               ("baseline", "candidate"))
                al = Path(str(snapshots["baseline"]) + f".{variant}.bin")
                bl = Path(str(snapshots["candidate"]) + f".{variant}.bin")
                if not al.is_file() or not bl.is_file():
                    issues.append(f"{case} rate{rate} {variant}: missing snapshot")
                    continue
                identical = equal_bytes(al, bl)
                digest_a = file_sha256(al)
                digest_b = file_sha256(bl)
                row = dict(case=case, rate=rate, variant=variant,
                           factor_bytes_identical=int(identical),
                           baseline_sha256=digest_a, candidate_sha256=digest_b,
                           baseline_bytes=al.stat().st_size,
                           candidate_bytes=bl.stat().st_size,
                           baseline_status=left["status"],
                           candidate_status=right["status"],
                           baseline_iterations=left["iterations"],
                           candidate_iterations=right["iterations"],
                           baseline_moved=left["moved_columns"],
                           candidate_moved=right["moved_columns"],
                           baseline_residual=left["physical_relres"],
                           candidate_residual=right["physical_relres"])
                report.append(row)
                if not identical:
                    issues.append(f"{case} rate{rate} {variant}: factor bytes differ")
                for key in ("status", "iterations", "nnz_L", "nnz_U",
                            "moved_columns", "solver_status", "api_error"):
                    if left[key] != right[key]:
                        issues.append(f"{case} rate{rate} {variant}: {key} differs")
                if variant != "legacy_ilut":
                    for label, rec in (("baseline", left), ("candidate", right)):
                        try:
                            resid = float(rec["physical_relres"])
                        except ValueError:
                            resid = float("inf")
                        if (rec["status"] != "converged" or
                                rec["api_error"] != "0" or
                                not 0 <= resid <= 1e-7):
                            issues.append(f"{case} rate{rate} {variant} {label}: "
                                          "ILUTP did not meet residual/convergence gate")
                    if swapped and int(left["moved_columns"]) == 0:
                        issues.append(f"{case} rate{rate} {variant}: no actual pivot")
    if report:
        with (out / "factor_identity.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(report[0]))
            writer.writeheader()
            writer.writerows(report)
    manifest = dict(baseline_ref=REFS["baseline"],
                    candidate_ref=REFS["candidate"],
                    common_original_native_benchmark_sha256=native_sha,
                    snapshots="native LIS_INT and LIS_SCALAR bytes, same runner ABI",
                    inputs=inputs, rates=RATES, variants=VARIANTS,
                    cases=[c[0] for c in CASES])
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    result = ("FAIL:\n" + "\n".join(issues) + "\n" if issues else
              f"PASS: {len(report)} factor snapshots byte-identical; "
              "including full L/U index+value arrays, D and permutation\n")
    (out / "validation.txt").write_text(result)
    print(result)
    return 1 if issues else 0


if __name__ == "__main__":
    sys.exit(main())
