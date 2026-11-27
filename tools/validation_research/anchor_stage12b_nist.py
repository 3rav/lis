#!/usr/bin/env python3
"""Public stiffness matrices, symmetric diagonal equilibration, manufactured RHS and synthetic contact perturbations."""
import argparse
import hashlib
import json
import math
import os
import statistics
import subprocess
from pathlib import Path
import numpy as np
from scipy import io, sparse
import anchor_stage10_contact as contact

STATES=[(0,0,0),(.15,.1,0),(.5,1,0),(1,10,.2),(0,0,0),(.7,3,.1)]
EXPECTED={"bcsstk16":(4884,290378),"bcsstk17":(10974,428650),"bcsstk18":(11948,149090)}

def generate(source,root,name):
    root.mkdir()
    raw=sparse.csr_matrix(io.mmread(source),dtype=np.float64)
    raw.sum_duplicates(); raw.eliminate_zeros(); raw.sort_indices()
    n,nnz=EXPECTED[name]
    contact.need(raw.shape==(n,n) and raw.nnz==nnz,"unexpected input shape/nnz")
    d=raw.diagonal()
    contact.need(np.all(np.isfinite(raw.data)) and np.all(d>0),"nonfinite input or nonpositive diagonal")
    scale=sparse.diags(1/np.sqrt(d))
    base=(scale@raw@scale).tocsr()
    base.sum_duplicates(); base.eliminate_zeros(); base.sort_indices()
    rng=np.random.default_rng(20261126)
    pairs=rng.permutation(n)[:2*(n//10)].reshape(-1,2)
    args=[]; records=[]
    dof=root/"dofs.csv"; dof.write_text("dof\n"+"".join(f"{i}\n" for i in range(1,n+1)))
    for k,(fraction,penalty,mu) in enumerate(STATES,1):
        m=round(fraction*len(pairs)); i,j=pairs[:m].T
        p=penalty*(1+.1*np.sin(i*.019+k*.13)); t=(i+1)%n
        rr=np.concatenate([i,j,i,j,i,i,j,j])
        cc=np.concatenate([i,j,j,i,t,j,t,j])
        vv=np.concatenate([p,p,-p,-p,p*mu,-p*mu,-p*mu,p*mu])
        a=(base+sparse.coo_matrix((vv,(rr,cc)),shape=(n,n))).tocsr()
        a.sum_duplicates(); a.eliminate_zeros(); a.sort_indices()
        u=np.sin(.037*np.arange(n)+.07*k)+.2*np.cos(.091*np.arange(n)-.03*k)
        b=a@u
        ap=root/f"A-{k:03}.mtx"; bp=root/f"b-{k:03}.mtx"
        with ap.open("w") as f:
            f.write(f"%%MatrixMarket matrix coordinate real general\n{n} {n} {a.nnz}\n")
            for row in range(n):
                lo,hi=a.indptr[row:row+2]
                f.writelines(f"{row+1} {int(jj)+1} {v:.17g}\n" for jj,v in zip(a.indices[lo:hi],a.data[lo:hi]))
        bp.write_text(f"%%MatrixMarket matrix array real general\n{n} 1\n"+"".join(f"{v:.17g}\n" for v in b))
        args.extend([str(ap),str(bp)])
        records.append(dict(step=k,a=str(ap),b=str(bp),active=m,penalty=penalty,friction=mu,nnz=int(a.nnz)))
    return args,records,dof,n

def main():
    p=argparse.ArgumentParser()
    p.add_argument("bench",type=Path); p.add_argument("source",type=Path); p.add_argument("root",type=Path)
    p.add_argument("--name",choices=EXPECTED,required=True); p.add_argument("--reps",type=int,required=True)
    a=p.parse_args()
    root=a.root.resolve(); root.mkdir()
    bench=a.bench.resolve(); source=a.source.resolve()
    args,records,dof,n=generate(source,root/"input",a.name)
    files=[bench,source,dof,*map(Path,args)]
    hashes={str(v):contact.sha(v) for v in files}
    (root/"manifest.json").write_text(json.dumps(dict(name=a.name,n=n,source="NIST Matrix Market / Harwell-Boeing BCSSTRUC2",scaled=True,
        perturbations="synthetic, not original contact history",records=records,hashes=hashes),indent=2))
    env=os.environ.copy(); env.update(OMP_NUM_THREADS="1",OMP_DYNAMIC="FALSE",LC_ALL="C",OPENBLAS_NUM_THREADS="1")
    env.pop("LIS_ANCHOR_AUDIT_DIR",None); cpu=min(os.sched_getaffinity(0))
    signatures={}; timings={}; measurements=[]; certs=[]; failures=[]
    modes=contact.MODES
    for kind in ("timing","audit"):
        for rep in range(a.reps):
            for offset in range(len(modes)):
                mode=modes[(rep+offset)%len(modes)]
                stem=root/f"{kind}-{rep+1}-{mode}"
                log=stem.with_suffix(".log"); tim=stem.with_suffix(".time")
                audit=Path(str(stem)+"-solutions"); runenv=env.copy()
                if kind=="audit":
                    audit.mkdir(); runenv["LIS_ANCHOR_AUDIT_DIR"]=str(audit)
                try:
                    with log.open("w") as f:
                        q=subprocess.run(["/usr/bin/time","-v","-o",str(tim),"taskset","-c",str(cpu),str(bench),mode,*args],
                            stdout=f,stderr=subprocess.STDOUT,env=runenv,timeout=600)
                    stem.with_suffix(".exit").write_text(str(q.returncode))
                    contact.need(q.returncode==0,f"exit={q.returncode}")
                    rows,s=contact.parse(log,mode,steps=len(STATES)); sig=contact.signature(rows)
                    contact.need(signatures.setdefault(mode,sig)==sig,"repeat mismatch")
                    if kind=="timing":
                        timings[rep,mode]=sig
                        measurements.append(dict(mode=mode,rep=rep+1,total_s=float(s["total_s"]),setup_s=float(s["setup_s"]),factors=int(s["factors"]),fallbacks=sum(v["fallback"] for v in rows),
                            iterations=sum(min(v["pilot_iter"],20)+min(v["continuation_iter"],80)+min(v["fresh_iter"],100) for v in rows)))
                    else:
                        contact.need(timings.get((rep,mode))==sig,"paired mismatch")
                        contact.need(set(audit.iterdir())=={audit/f"x-{k:03}.mtx" for k in range(1,len(STATES)+1)},"solution set")
                        batch=[]
                        for row,rec in zip(rows,records):
                            x=audit/f"x-{row['step']:03}.mtx"
                            h=contact.sha(x); res=contact.certify(Path(rec["a"]),Path(rec["b"]),x,n)
                            contact.need(math.isfinite(res) and res<=1e-8 and h==contact.sha(x),f"physical residual={res}")
                            batch.append(dict(mode=mode,rep=rep+1,step=row["step"],residual=res,sha256=h))
                        certs.extend(batch)
                    print("PUBLIC_RUN",a.name,kind,rep+1,mode,s["accepted"],flush=True)
                except (RuntimeError,ValueError,KeyError,OSError,subprocess.TimeoutExpired) as e:
                    failures.append(dict(mode=mode,rep=rep+1,kind=kind,error=str(e)))
                    print("PUBLIC_FAILURE",failures[-1],flush=True)
    contact.need(all(contact.sha(Path(v))==h for v,h in hashes.items()),"input/binary changed")
    med={}
    for m in modes:
        ms=[v for v in measurements if v["mode"]==m]; cs=[v for v in certs if v["mode"]==m]
        if len(ms)==a.reps and len(cs)==len(STATES)*a.reps and not any(v["mode"]==m for v in failures):
            med[m]={k:statistics.median(v[k] for v in ms) for k in ms[0] if k not in ("mode","rep")}
            med[m]["max_residual"]=max(v["residual"] for v in cs)
    ok=len(med)==len(modes)
    report=dict(status="PASS" if ok else "FAIL",name=a.name,n=n,reps=a.reps,eligible=med,failures=failures,certifications=certs,measurements=measurements,
        note="Equilibrated public stiffness matrix with synthetic perturbations; not CCX")
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("PUBLIC_GATE="+report["status"],json.dumps({k:report[k] for k in ("name","n","reps","eligible","failures")}),flush=True)
    return 0 if ok else 1

if __name__=="__main__":
    raise SystemExit(main())
