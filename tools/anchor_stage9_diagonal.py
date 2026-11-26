#!/usr/bin/env python3
"""Native diagonal-fixture benchmark, not CalculiX or production timing evidence."""
import hashlib
import json
import math
import os
import re
import statistics
import subprocess
import sys
from decimal import Decimal, localcontext
from pathlib import Path

MODES = ("pilot", "pilot-restart-matched", "pilot-single-pass")
INTS = "step accepted api raw_status pilot_iter continuation_iter fresh_iter continued fallback factors".split()
FLOATS = "ratio residual setup_s total_s".split()
N = 80
STEPS = 15

def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)

def digest(p):
    with p.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()

def array(p):
    ls = p.read_text().splitlines()
    need(ls[0] == "%%MatrixMarket matrix array real general", str(p))
    body = [l for l in ls[1:] if not l.startswith("%")]
    need(body[0].split() == [str(N), "1"] and len(body) == N+1, str(p))
    vals = [float(v) for v in body[1:]]
    need(all(map(math.isfinite, vals)), str(p))
    return [Decimal.from_float(v) for v in vals]

def diagonal(p):
    ls = p.read_text().splitlines()
    need(ls[0] == "%%MatrixMarket matrix coordinate real general", str(p))
    need(ls[1].split() == [str(N)]*3 and len(ls)==N+2, str(p))
    result = []
    for i, line in enumerate(ls[2:], 1):
        a,b,v = line.split()
        need(a==b==str(i) and math.isfinite(float(v)), str(p))
        result.append(Decimal.from_float(float(v)))
    return result

def physical(a,b,x):
    with localcontext() as c:
        c.prec = 70
        nr = sum((d*y-z)**2 for d,z,y in zip(a,b,x)).sqrt()
        nb = sum(z*z for z in b).sqrt()
        return nr/nb if nb else nr

def generate(root):
    root.mkdir()
    args, inputs, dofs = [], [], []
    for k in range(1,STEPS+1):
        stem = root/f"pivot1e9_test2_{k:020d}_newton"
        a,b,d = (Path(str(stem)+s) for s in ("_A.mtx","_b.mtx","_dofs.csv"))
        kind = (k-1)%3
        vals = [1.0 if kind==0 else float(i+1) if kind==1 else 10**(-8+8*i/(N-1)) for i in range(N)]
        a.write_text("%%MatrixMarket matrix coordinate real general\n"+f"{N} {N} {N}\n"+"".join(f"{i+1} {i+1} {v:.17g}\n" for i,v in enumerate(vals)))
        b.write_text("%%MatrixMarket matrix array real general\n"+f"{N} 1\n"+"1\n"*N)
        d.write_text("dof\n"+"".join(f"{i+1}\n" for i in range(N)))
        args.extend([str(a),str(b)])
        inputs.append((a,b))
        dofs.append(d)
    need(len({digest(d) for d in dofs})==1,"DOF mismatch")
    return args,inputs,dofs

def kv(line):
    return dict(v.split("=",1) for v in line.split()[1:])

def parse(log,mode):
    rows=[]; summaries=[]
    for line in log.read_text().splitlines():
        if line.startswith("STEP "):
            v=kv(line)
            r={k:int(v[k]) for k in INTS}
            r.update({k:float(v[k]) for k in FLOATS})
            rows.append(r)
        elif line.startswith("SUMMARY "):
            summaries.append(kv(line))
    need(len(rows)==STEPS and len(summaries)==1,f"incomplete {log}")
    for i,r in enumerate(rows):
        need(r["step"]==i+1 and r["accepted"]==1 and r["api"]==0 and 0<=r["residual"]<=1e-8,f"acceptance {log} {i+1}")
        need(r["continued"] in (0,1) and r["fallback"] in (0,1),f"flags {log}")
        need(r["factors"]==int(i==0)+r["fallback"],f"factors {log}")
        need(bool(r["continuation_iter"])==bool(r["continued"]),f"continuation {log}")
        need(r["fallback"] or not r["fresh_iter"],f"unexpected fresh {log}")
        caps=(21,81,101) if mode=="pilot" else (20,80,101)
        need(all(0<=r[k]<=cap for k,cap in zip(("pilot_iter","continuation_iter","fresh_iter"),caps)),f"budget {log}")
        if r["continued"]:
            need(min(r["pilot_iter"],20)==20 and 0<=r["ratio"]<=0.25,f"checkpoint {log}")
        need(0<=r["setup_s"]<=r["total_s"]+1e-6 and math.isfinite(r["total_s"]),f"time {log}")
    s=summaries[0]
    need(s["accepted"]==f"{STEPS}/{STEPS}" and int(s["factors"])==sum(r["factors"] for r in rows),f"summary {log}")
    for k in ("setup_s","total_s"):
        need(abs(float(s[k])-sum(r[k] for r in rows))<=(STEPS+1)*1e-6,f"summary timing {log}")
    return rows,s

def signature(rows):
    return [tuple(r[k] for k in INTS) for r in rows]

def main():
    bench=Path(sys.argv[1]).resolve()
    root=Path(sys.argv[2]).resolve()
    root.mkdir()
    args,inputs,dofs=generate(root/"input")
    tracked=[bench,*[p for pair in inputs for p in pair],*dofs]
    hashes={str(p):digest(p) for p in tracked}
    (root/"inputs.sha256").write_text("".join(f"{h}  {p}\n" for p,h in hashes.items()))
    env=os.environ.copy()
    env.update(OMP_NUM_THREADS="1",OMP_DYNAMIC="FALSE",LC_ALL="C")
    env.pop("LIS_ANCHOR_AUDIT_DIR",None)
    cpu=min(os.sched_getaffinity(0))
    signatures={}; timed={}; measurements=[]; certs=[]
    for kind in ("timing","audit"):
        for rep in range(3):
            for offset in range(3):
                mode=MODES[(rep+offset)%3]
                stem=root/f"{kind}-rep{rep+1}-{mode}"
                log=stem.with_suffix(".log")
                tim=stem.with_suffix(".time")
                audit=Path(str(stem)+"-solutions")
                runenv=env.copy()
                if kind=="audit":
                    audit.mkdir()
                    runenv["LIS_ANCHOR_AUDIT_DIR"]=str(audit)
                with log.open("w") as f:
                    p=subprocess.run(["/usr/bin/time","-v","-o",str(tim),"taskset","-c",str(cpu),str(bench),mode,*args],env=runenv,stdout=f,stderr=subprocess.STDOUT)
                stem.with_suffix(".exit").write_text(str(p.returncode)+"\n")
                need(p.returncode==0,f"exit {p.returncode}: {log}")
                rows,s=parse(log,mode)
                sig=signature(rows)
                need(signatures.setdefault(mode,sig)==sig,f"repeat mismatch {stem}")
                if kind=="timing":
                    timed[rep,mode]=rows
                    mem=re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)",tim.read_text())
                    need(mem is not None,"RSS missing")
                    measurements.append(dict(rep=rep+1,mode=mode,total_s=float(s["total_s"]),setup_s=float(s["setup_s"]),
                        anchor_iterations=sum(min(r["pilot_iter"],20)+min(r["continuation_iter"],80) for r in rows),
                        fresh_iterations_capped=sum(min(r["fresh_iter"],100) for r in rows),fallbacks=sum(r["fallback"] for r in rows),
                        factors=int(s["factors"]),rss_kib=int(mem[1])))
                else:
                    t=timed[rep,mode]
                    need(signature(t)==sig,f"paired mismatch {stem}")
                    expected={audit/f"x-{i:03d}.mtx" for i in range(1,STEPS+1)}
                    need(set(audit.iterdir())==expected,f"solution set {stem}")
                    for i,(r,tr,(a,b)) in enumerate(zip(rows,t,inputs),1):
                        if r["continued"]:
                            need(r["ratio"]==tr["ratio"],f"ratio mismatch {stem}")
                        x=audit/f"x-{i:03d}.mtx"
                        h=digest(x)
                        res=physical(diagonal(a),array(b),array(x))
                        need(res<=Decimal("1e-8"),f"physical residual {res} {x}")
                        need(h==digest(x),f"changed solution {x}")
                        certs.append(dict(rep=rep+1,mode=mode,step=i,residual=str(res),reported=r["residual"],sha256=h,path=str(x)))
                print(f"RUN {kind} rep={rep+1} mode={mode} accepted={s['accepted']} factors={s['factors']}",flush=True)
    need(len(certs)==135 and len(measurements)==9,"incomplete benchmark")
    need(all(digest(Path(p))==h for p,h in hashes.items()),"input/binary changed")
    med={m:{k:statistics.median(r[k] for r in measurements if r["mode"]==m) for k in measurements[0] if k not in ("rep","mode")} for m in MODES}
    report=dict(status="PASS",fixture="synthetic diagonal 80x80; not CalculiX",cpu=cpu,paired_cases=len(certs),native_solves=270,
                max_independent_residual=str(max(Decimal(c["residual"]) for c in certs)),medians=med,measurements=measurements,certifications=certs,inputs_sha256=hashes)
    tmp=root/"report.tmp"
    tmp.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    tmp.replace(root/"report.json")
    print("DIAGONAL_PAIRED_GATE=PASS",json.dumps({k:report[k] for k in ("paired_cases","max_independent_residual","medians")}),flush=True)

if __name__=="__main__":
    main()
