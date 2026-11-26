#!/usr/bin/env python3
"""Fixed-budget gate/restart screening; no change to default production settings."""
import argparse
import json
import math
import os
import statistics
import subprocess
from pathlib import Path
import anchor_stage10_contact as contact

CONFIGS=((.10,60),(.25,60),(.50,60),(.75,60),(.25,40),(.25,100))
STRESS=[(0,0,0),(1,300,0),(.05,2,0),(1,1000,0),(.5,100,.1),
        (0,0,0),(1,1500,.2),(.15,10,0),(1,500,.3),(.25,50,0),
        (0,0,0),(.9,800,.15),(1,1500,.3),(.1,5,0),(0,0,0)]

def main():
    p=argparse.ArgumentParser()
    p.add_argument("root",type=Path)
    p.add_argument("--width",type=int,required=True)
    p.add_argument("--height",type=int,required=True)
    p.add_argument("--reps",type=int,default=3)
    p.add_argument("--stress",action="store_true")
    a=p.parse_args()
    root=a.root.resolve(); root.mkdir()
    n=a.width*a.height
    if a.stress: contact.CASES=STRESS
    args,records,dof=contact.generate(root/"input",a.width,a.height)
    source=Path("tools/lis_anchor_pilot.c").read_text()
    contact.need(source.count("<=0.25")==2 and source.count("-restart 60 ")==1,"source anchors changed")
    bins={}
    for gate,restart in CONFIGS:
        tag=f"g{int(gate*100):02}-r{restart}"
        src=root/f"{tag}.c"; exe=root/f"{tag}.bin"
        src.write_text(source.replace("<=0.25",f"<={gate:.2f}").replace("-restart 60 ",f"-restart {restart} "))
        subprocess.run(["cc","-O2","-DNDEBUG","-fopenmp","-DHAVE_CONFIG_H","-Iinclude","-Itools",str(src),"tools/bench_anchor_pilot.c","src/.libs/liblis.a","-lm","-o",str(exe)],check=True)
        bins[gate,restart]=exe
    hashes={str(v):contact.sha(v) for v in [dof,*map(Path,args),*bins.values()]}
    (root/"manifest.json").write_text(json.dumps(dict(n=n,stress=a.stress,records=records,hashes=hashes),indent=2))
    candidates=[(.25,60,"fresh")]+[(g,r,m) for g,r in CONFIGS for m in ("pilot-restart-matched","pilot-single-pass")]
    env=os.environ.copy(); env.update(OMP_NUM_THREADS="1",OMP_DYNAMIC="FALSE",LC_ALL="C")
    env.pop("LIS_ANCHOR_AUDIT_DIR",None); cpu=min(os.sched_getaffinity(0))
    measurements=[]; certs=[]; failed=[]; sigs={}; paired={}
    for kind in ("timing","audit"):
        for rep in range(a.reps):
            order=candidates[rep%len(candidates):]+candidates[:rep%len(candidates)]
            if rep%2: order=order[::-1]
            for g,r,m in order:
                tag=f"g{int(g*100):02}-r{r}-{m}"
                stem=root/f"{kind}-{rep+1}-{tag}"
                log=stem.with_suffix(".log"); audit=Path(str(stem)+"-solutions")
                runenv=env.copy()
                if kind=="audit":
                    audit.mkdir(); runenv["LIS_ANCHOR_AUDIT_DIR"]=str(audit)
                try:
                    with log.open("w") as f:
                        q=subprocess.run(["taskset","-c",str(cpu),str(bins[g,r]),m,*args],stdout=f,stderr=subprocess.STDOUT,env=runenv,timeout=420)
                    stem.with_suffix(".exit").write_text(str(q.returncode))
                    contact.need(q.returncode==0,f"exit={q.returncode}")
                    rows,s=contact.parse(log,m,gate=g); sig=contact.signature(rows)
                    contact.need(sigs.setdefault(tag,sig)==sig,"repeat mismatch")
                    if kind=="timing":
                        paired[rep,tag]=sig
                        measurements.append(dict(tag=tag,rep=rep+1,total_s=float(s["total_s"]),setup_s=float(s["setup_s"]),factors=int(s["factors"]),fallbacks=sum(v["fallback"] for v in rows),
                            iterations=sum(min(v["pilot_iter"],20)+min(v["continuation_iter"],80)+min(v["fresh_iter"],100) for v in rows)))
                    else:
                        contact.need(paired.get((rep,tag))==sig,"paired mismatch")
                        contact.need(set(audit.iterdir())=={audit/f"x-{k:03}.mtx" for k in range(1,16)},"solution set")
                        batch=[]
                        for row,rec in zip(rows,records):
                            x=audit/f"x-{row['step']:03}.mtx"
                            h=contact.sha(x); res=contact.certify(Path(rec["a"]),Path(rec["b"]),x,n)
                            contact.need(math.isfinite(res) and res<=1e-8 and h==contact.sha(x),f"physical residual={res}")
                            batch.append(dict(tag=tag,rep=rep+1,step=row["step"],residual=res,sha256=h))
                        certs.extend(batch)
                    print("SWEEP_RUN",kind,rep+1,tag,s["accepted"],flush=True)
                except (RuntimeError,ValueError,KeyError,OSError,subprocess.TimeoutExpired) as e:
                    failed.append(dict(tag=tag,rep=rep+1,kind=kind,error=str(e)))
                    print("SWEEP_FAILURE",failed[-1],flush=True)
    contact.need(all(contact.sha(Path(v))==h for v,h in hashes.items()),"input/binary changed")
    med={}
    for g,r,m in candidates:
        tag=f"g{int(g*100):02}-r{r}-{m}"
        ms=[v for v in measurements if v["tag"]==tag]; cs=[v for v in certs if v["tag"]==tag]
        if len(ms)==a.reps and len(cs)==15*a.reps and not any(v["tag"]==tag for v in failed):
            med[tag]={k:statistics.median(v[k] for v in ms) for k in ms[0] if k not in ("tag","rep")}
            med[tag]["max_residual"]=max(v["residual"] for v in cs)
    baseline="g25-r60-fresh"
    controls=[baseline,"g25-r60-pilot-restart-matched","g25-r60-pilot-single-pass"]
    ok=all(t in med for t in controls)
    report=dict(status="PASS" if ok else "FAIL",n=n,stress=a.stress,reps=a.reps,eligible=med,failures=failed,certifications=certs,measurements=measurements,
                note="Screening only. Noneligible variants excluded. No automatic change of defaults.")
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("SWEEP_GATE="+report["status"],json.dumps({k:report[k] for k in ("n","stress","reps","eligible","failures")}),flush=True)
    return 0 if ok else 1

if __name__=="__main__":
    raise SystemExit(main())
