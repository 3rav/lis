#!/usr/bin/env python3
"""Contact-like sparse sequences: spring network + penalty constraints + nonsymmetric coupling.
Not a CCX simulation. Each strategy sees the same fixed DOF ordering.
"""
import argparse
import hashlib
import json
import math
import os
import re
import statistics
import subprocess
from pathlib import Path

MODES=("fresh","pilot","pilot-low-memory","pilot-rolling","pilot-restart-matched","pilot-single-pass")
INTS="step accepted api raw_status pilot_iter continuation_iter fresh_iter continued fallback factors".split()
FLOATS="ratio residual setup_s total_s".split()
# active fraction, normal penalty, nonsymmetric coupling; deliberate opening and reclosure.
CASES=[(0,0,0),(.15,2,0),(.35,5,0),(.65,10,0),(1,20,0),(1,60,0),
       (1,150,0),(.8,80,.12),(.6,40,.25),(.2,20,.25),(0,0,0),
       (.5,25,0),(.9,80,.15),(1,150,.25),(.4,10,0)]

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def sha(p):
    with p.open("rb") as f: return hashlib.file_digest(f,"sha256").hexdigest()

def array(p,n):
    with p.open() as f:
        need(f.readline().strip()=="%%MatrixMarket matrix array real general",str(p))
        shape=f.readline().split()
        need(shape==[str(n),"1"],str(p))
        x=[float(v) for v in f.read().split()]
    need(len(x)==n and all(map(math.isfinite,x)),str(p))
    return x

def generate(root,w,h):
    root.mkdir()
    n=w*h
    base=[{} for _ in range(n)]
    def add(rows,i,j,v):
        rows[i][j]=rows[i].get(j,0.0)+v
    def spring(rows,i,j,k):
        add(rows,i,i,k); add(rows,j,j,k)
        add(rows,i,j,-k); add(rows,j,i,-k)
    for y in range(h):
        for x in range(w):
            i=y*w+x
            add(base,i,i,.025+(.5 if y==0 or y==h-1 else 0))
            if x+1<w: spring(base,i,i+1,1+.2*math.sin(i*.17))
            if y+1<h: spring(base,i,i+w,1+.2*math.cos(i*.11))
    # Multiple interface rows, with nonlocal contacts across a gap in the grid.
    pairs=[((y-1)*w+x,(y+2)*w+x) for y in range(2,h-2, max(3,h//8)) for x in range(w)]
    pairs.sort(key=lambda ij: ((ij[0]*2654435761)&0xffffffff,ij[0]))
    args=[]; records=[]
    dof=root/"dofs.csv"
    dof.write_text("dof\n"+"".join(f"{i}\n" for i in range(1,n+1)))
    for k,(fraction,penalty,friction) in enumerate(CASES,1):
        rows=[r.copy() for r in base]
        active=pairs[:round(fraction*len(pairs))]
        for i,j in active:
            p=penalty*(1+.15*math.sin(i*.019+k*.13))
            spring(rows,i,j,p)
            # Nonsymmetric rank-one term p*mu*(e_i-e_j)*(e_t-e_j)^T.
            # This is an algebraic friction surrogate, not a constitutive law.
            t=i+1 if i%w<w-1 else i-1
            if friction:
                q=p*friction
                add(rows,i,t,q); add(rows,i,j,-q)
                add(rows,j,t,-q); add(rows,j,j,q)
        rows=[{j:v for j,v in r.items() if v!=0} for r in rows]
        xstar=[math.sin(.037*i+.07*k)+.2*math.cos(.091*i-.03*k) for i in range(n)]
        rhs=[math.fsum(v*xstar[j] for j,v in sorted(r.items())) for r in rows]
        a=root/f"A-{k:03}.mtx"; b=root/f"b-{k:03}.mtx"
        nnz=sum(map(len,rows))
        with a.open("w") as f:
            f.write(f"%%MatrixMarket matrix coordinate real general\n{n} {n} {nnz}\n")
            for i,r in enumerate(rows):
                f.writelines(f"{i+1} {j+1} {v:.17g}\n" for j,v in sorted(r.items()))
        b.write_text(f"%%MatrixMarket matrix array real general\n{n} 1\n"+"".join(f"{v:.17g}\n" for v in rhs))
        args.extend([str(a),str(b)])
        records.append(dict(step=k,active=len(active),penalty=penalty,friction=friction,nnz=nnz,a=str(a),b=str(b)))
    return args,records,dof

def certify(a,b,x,n):
    # Independent streaming read of stored Matrix Market entries, no LIS APIs.
    # math.fsum accumulates each row and norms accurately in binary64.
    rhs=array(b,n); sol=array(x,n)
    residual=[]
    with a.open() as f:
        need(f.readline().strip()=="%%MatrixMarket matrix coordinate real general",str(a))
        nr,nc,nnz=map(int,f.readline().split())
        need(nr==nc==n,str(a))
        row=0; terms=[]; lastj=-1; seen=0
        for line in f:
            i,j,v=line.split(); i=int(i)-1; j=int(j)-1; v=float(v)
            need(0<=i<n and 0<=j<n and math.isfinite(v),str(a))
            if i!=row:
                need(i==row+1 and terms,str(a))
                residual.append(math.fsum(terms+[-rhs[row]]))
                row=i; terms=[]; lastj=-1
            need(j>lastj,str(a)); lastj=j; seen+=1
            terms.append(v*sol[j])
        need(seen==nnz and row==n-1 and terms,str(a))
        residual.append(math.fsum(terms+[-rhs[row]]))
    rn=math.sqrt(math.fsum(v*v for v in residual))
    bn=math.sqrt(math.fsum(v*v for v in rhs))
    return rn/bn if bn else rn

def kv(line):
    return dict(v.split("=",1) for v in line.split()[1:])

def parse(log,mode,gate=.25,steps=15):
    rows=[]; summary=[]
    for line in log.read_text().splitlines():
        if line.startswith("STEP "):
            d=kv(line); r={k:int(d[k]) for k in INTS}
            r.update({k:float(d[k]) for k in FLOATS}); rows.append(r)
        elif line.startswith("SUMMARY "): summary.append(kv(line))
    need(len(rows)==steps and len(summary)==1,f"incomplete {log}")
    for i,r in enumerate(rows):
        need(r["step"]==i+1 and r["accepted"]==1 and r["api"]==0 and 0<=r["residual"]<=1e-8,f"acceptance {log} {i+1}")
        need(r["continued"] in (0,1) and r["fallback"] in (0,1),f"flags {log}")
        need(0<=r["setup_s"]<=r["total_s"]+1e-6 and math.isfinite(r["total_s"]),f"time {log}")
        need(0<=r["pilot_iter"]<=21 and 0<=r["continuation_iter"]<=81 and 0<=r["fresh_iter"]<=101,f"budget {log}")
        if mode=="fresh":
            need(r["factors"]==1 and not any(r[k] for k in ("pilot_iter","continuation_iter","continued","fallback")),f"fresh {log}")
        else:
            need(r["factors"]>=int(i==0)+r["fallback"],f"factors {log}")
            if r["continued"]:
                need(r["pilot_iter"] in (20,21) and 0<=r["ratio"]<=gate,f"gate {log}")
            need(r["fallback"] or r["fresh_iter"]==0,f"unexpected fresh {log}")
    s=summary[0]
    need(s["accepted"]==f"{steps}/{steps}" and int(s["factors"])==sum(r["factors"] for r in rows),f"summary {log}")
    for k in ("setup_s","total_s"):
        need(abs(float(s[k])-sum(r[k] for r in rows))<=(steps+2)*1e-6,f"summary time {log}")
    return rows,s

def signature(rows):
    # Exact raw counters, including legacy maxiter+1; no rounding in pairing.
    return [tuple(r[k] for k in INTS)+tuple(None if math.isnan(r[k]) else r[k] for k in ("ratio","residual")) for r in rows]

def main():
    p=argparse.ArgumentParser()
    p.add_argument("bench",type=Path); p.add_argument("root",type=Path)
    p.add_argument("--width",type=int,required=True); p.add_argument("--height",type=int,required=True)
    p.add_argument("--reps",type=int,default=3)
    a=p.parse_args(); n=a.width*a.height
    need(a.width>=4 and a.height>=8 and a.reps>=1,"bad dimensions/reps")
    bench=a.bench.resolve(); root=a.root.resolve(); root.mkdir()
    args,records,dof=generate(root/"input",a.width,a.height)
    hashes={str(p):sha(p) for p in [bench,dof,*map(Path,args)]}
    (root/"manifest.json").write_text(json.dumps(dict(n=n,records=records,hashes=hashes),indent=2))
    env=os.environ.copy(); env.update(OMP_NUM_THREADS="1",OMP_DYNAMIC="FALSE",LC_ALL="C")
    env.pop("LIS_ANCHOR_AUDIT_DIR",None); cpu=min(os.sched_getaffinity(0))
    signatures={}; timings={}; measurements=[]; certs=[]; failures=[]
    for kind in ("timing","audit"):
        for rep in range(a.reps):
            for offset in range(len(MODES)):
                mode=MODES[(rep+offset)%len(MODES)]
                stem=root/f"{kind}-rep{rep+1}-{mode}"
                log=stem.with_suffix(".log"); tim=stem.with_suffix(".time")
                audit=Path(str(stem)+"-solutions"); runenv=env.copy()
                if kind=="audit":
                    audit.mkdir(); runenv["LIS_ANCHOR_AUDIT_DIR"]=str(audit)
                try:
                    with log.open("w") as f:
                        r=subprocess.run(["/usr/bin/time","-v","-o",str(tim),"taskset","-c",str(cpu),str(bench),mode,*args],env=runenv,stdout=f,stderr=subprocess.STDOUT,timeout=420)
                    stem.with_suffix(".exit").write_text(str(r.returncode)+"\n")
                    need(r.returncode==0,f"exit={r.returncode}")
                    rows,s=parse(log,mode); sig=signature(rows)
                    need(signatures.setdefault(mode,sig)==sig,"repeat mismatch")
                    if kind=="timing":
                        timings[rep,mode]=sig
                        mem=re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)",tim.read_text())
                        need(mem is not None,"RSS missing")
                        measurements.append(dict(rep=rep+1,mode=mode,total_s=float(s["total_s"]),setup_s=float(s["setup_s"]),
                            anchor_iterations_capped=sum(min(r["pilot_iter"],20)+min(r["continuation_iter"],80) for r in rows),
                            fresh_iterations_capped=sum(min(r["fresh_iter"],100) for r in rows),fallbacks=sum(r["fallback"] for r in rows),
                            factors=int(s["factors"]),rss_kib=int(mem[1])))
                    else:
                        need(timings.get((rep,mode))==sig,"paired mismatch")
                        need(set(audit.iterdir())=={audit/f"x-{k:03}.mtx" for k in range(1,16)},"solution set")
                        for row,rec in zip(rows,records):
                            x=audit/f"x-{row['step']:03}.mtx"; h=sha(x)
                            res=certify(Path(rec["a"]),Path(rec["b"]),x,n)
                            need(math.isfinite(res) and res<=1e-8,f"independent residual={res} {x}")
                            need(h==sha(x),"solution changed")
                            certs.append(dict(rep=rep+1,mode=mode,step=row["step"],residual=res,reported=row["residual"],sha256=h))
                        # Raw solution files remain archived for reproducibility.
                    print(f"RUN n={n} {kind} rep={rep+1} mode={mode} accepted={s['accepted']} factors={s['factors']}",flush=True)
                except (RuntimeError,ValueError,KeyError,OSError,subprocess.TimeoutExpired) as e:
                    failures.append(dict(kind=kind,rep=rep+1,mode=mode,error=str(e)))
                    print("RUN_FAILED",failures[-1],flush=True)
    need(all(sha(Path(p))==h for p,h in hashes.items()),"input/binary changed")
    med={}
    for m in MODES:
        rs=[r for r in measurements if r["mode"]==m]
        if len(rs)==a.reps:
            med[m]={k:statistics.median(r[k] for r in rs) for k in rs[0] if k not in ("rep","mode")}
    ok=not failures and len(certs)==15*len(MODES)*a.reps and len(measurements)==len(MODES)*a.reps
    report=dict(status="PASS" if ok else "FAIL",n=n,reps=a.reps,cpu=cpu,paired_cases=len(certs),
                max_independent_residual=max((c["residual"] for c in certs),default=None),medians=med,
                measurements=measurements,certifications=certs,failures=failures,fixture="synthetic spring/contact surrogate; not CCX")
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("CONTACT_GATE="+report["status"],json.dumps({k:report[k] for k in ("n","reps","paired_cases","max_independent_residual","medians","failures")}),flush=True)
    return 0 if ok else 1

if __name__=="__main__":
    raise SystemExit(main())
