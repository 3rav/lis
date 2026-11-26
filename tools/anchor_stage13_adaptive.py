#!/usr/bin/env python3
"""Holdout screening of early stagnation and conservative predicted continuation."""
import argparse
import json
import math
import os
import statistics
import subprocess
from pathlib import Path
import anchor_stage10_contact as contact

CONFIGS=("fixed","early","predict")
STRESS=[(0,0,0),(.9,200,.1),(.1,5,0),(1,700,.25),(.3,30,.1),
        (1,1200,0),(.5,150,.3),(0,0,0),(1,900,.15),(.05,3,0),
        (.8,400,.25),(1,1500,.2),(0,0,0),(.2,20,0),(1,600,.3)]

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
    old_callback="static LIS_INT checkpoint(LIS_INT iter,LIS_REAL res,void *data)\n{\n    pilot_history *h=data;\n    h->last=iter;\n    if(iter==10)h->r10=res;\n    if(iter!=20)return LIS_FALSE;\n    h->r20=res;\n    return !(isfinite(h->r10) && h->r10>0 && isfinite(h->r20) &&\n             h->r20>=0 && h->r20/h->r10<=0.25);\n}\n"
    new_callback="static LIS_INT checkpoint(LIS_INT iter,LIS_REAL res,void *data)\n{\n    pilot_history *h=data;\n    double ratio;\n    h->last=iter;\n    if(iter==4)h->r4=res;\n    if(iter==8)h->r8=res;\n    if(iter==10)h->r10=res;\n    /* Two disjoint four-step windows must both stagnate.\n       A near-tolerance projected residual is never rejected early. */\n    if(iter==12 && isfinite(h->r4) && isfinite(h->r8) && isfinite(res) &&\n       h->r4>0 && h->r8>0 && res>1e-4 &&\n       h->r8/h->r4>=0.97 && res/h->r8>=0.97)return LIS_TRUE;\n    if(iter!=20)return LIS_FALSE;\n    h->r20=res;\n    if(!isfinite(h->r10) || h->r10<=0 || !isfinite(res) || res<0)return LIS_TRUE;\n    ratio=res/h->r10;\n    if(ratio<=0.25)return LIS_FALSE;\n    /* Experimental only: continue if a conservative geometric extrapolation\n       reaches the projected target before the first restart at iteration 60.\n       Physical acceptance and the hard budget remain unchanged. */\n    PREDICT_BRANCH\n    return LIS_TRUE;\n}\n"
    contact.need(source.count(old_callback)==1 and source.count("double r10,r20; LIS_INT last;")==1,"source anchors changed")
    bins={}
    for tag in CONFIGS:
        src=root/f"{tag}.c"; exe=root/f"{tag}.bin"
        body=source
        if tag!="fixed":
            prediction="if(ratio<1.0 && res>1e-8 && 10.0*log(1e-8/res)/log(ratio)<=35.0)return LIS_FALSE;" if tag=="predict" else ""
            body=source.replace("double r10,r20; LIS_INT last;","double r10,r20; LIS_INT last; double r4,r8;").replace(old_callback,new_callback.replace("PREDICT_BRANCH",prediction))
        src.write_text(body)
        subprocess.run(["cc","-O2","-DNDEBUG","-fopenmp","-DHAVE_CONFIG_H","-Iinclude","-Itools",str(src),"tools/bench_anchor_pilot.c","src/.libs/liblis.a","-lm","-o",str(exe)],check=True)
        bins[tag]=exe
    hashes={str(v):contact.sha(v) for v in [dof,*map(Path,args),*bins.values()]}
    (root/"manifest.json").write_text(json.dumps(dict(n=n,stress=a.stress,records=records,hashes=hashes),indent=2))
    candidates=[("fixed","fresh"),("fixed","pilot-restart-matched")]+[(tag,"pilot-single-pass") for tag in CONFIGS]
    env=os.environ.copy(); env.update(OMP_NUM_THREADS="1",OMP_DYNAMIC="FALSE",LC_ALL="C")
    env.pop("LIS_ANCHOR_AUDIT_DIR",None); cpu=min(os.sched_getaffinity(0))
    measurements=[]; certs=[]; failed=[]; sigs={}; paired={}
    for kind in ("timing","audit"):
        for rep in range(a.reps):
            order=candidates[rep%len(candidates):]+candidates[:rep%len(candidates)]
            if rep%2: order=order[::-1]
            for variant,m in order:
                tag=f"{variant}-{m}"
                stem=root/f"{kind}-{rep+1}-{tag}"
                log=stem.with_suffix(".log"); audit=Path(str(stem)+"-solutions")
                runenv=env.copy()
                if kind=="audit":
                    audit.mkdir(); runenv["LIS_ANCHOR_AUDIT_DIR"]=str(audit)
                try:
                    with log.open("w") as f:
                        q=subprocess.run(["taskset","-c",str(cpu),str(bins[variant]),m,*args],stdout=f,stderr=subprocess.STDOUT,env=runenv,timeout=420)
                    stem.with_suffix(".exit").write_text(str(q.returncode))
                    contact.need(q.returncode==0,f"exit={q.returncode}")
                    rows,s=contact.parse(log,m,gate=1.0 if variant=="predict" else .25); sig=contact.signature(rows)
                    contact.need(sigs.setdefault(tag,sig)==sig,"repeat mismatch")
                    if kind=="timing":
                        paired[rep,tag]=sig
                        measurements.append(dict(tag=tag,rep=rep+1,total_s=float(s["total_s"]),setup_s=float(s["setup_s"]),factors=int(s["factors"]),fallbacks=sum(v["fallback"] for v in rows),fallback_before20=sum(v["fallback"] and v["pilot_iter"]<20 for v in rows),relaxed_continuations=sum(v["continued"] and v["ratio"]>.25 for v in rows),
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
                    print("ADAPTIVE_RUN",kind,rep+1,tag,s["accepted"],flush=True)
                except (RuntimeError,ValueError,KeyError,OSError,subprocess.TimeoutExpired) as e:
                    failed.append(dict(tag=tag,rep=rep+1,kind=kind,error=str(e)))
                    print("ADAPTIVE_FAILURE",failed[-1],flush=True)
    contact.need(all(contact.sha(Path(v))==h for v,h in hashes.items()),"input/binary changed")
    med={}
    for variant,m in candidates:
        tag=f"{variant}-{m}"
        ms=[v for v in measurements if v["tag"]==tag]; cs=[v for v in certs if v["tag"]==tag]
        if len(ms)==a.reps and len(cs)==15*a.reps and not any(v["tag"]==tag for v in failed):
            med[tag]={k:statistics.median(v[k] for v in ms) for k in ms[0] if k not in ("tag","rep")}
            med[tag]["max_residual"]=max(v["residual"] for v in cs)
    baseline="fixed-fresh"
    controls=[baseline,"fixed-pilot-restart-matched","fixed-pilot-single-pass"]
    ok=all(t in med for t in controls)
    report=dict(status="PASS" if ok else "FAIL",n=n,stress=a.stress,reps=a.reps,eligible=med,failures=failed,certifications=certs,measurements=measurements,
                note="Holdout screening only; no automatic change of defaults.")
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("ADAPTIVE_GATE="+report["status"],json.dumps({k:report[k] for k in ("n","stress","reps","eligible","failures")}),flush=True)
    return 0 if ok else 1

if __name__=="__main__":
    raise SystemExit(main())
