#!/usr/bin/env python3
"""Online cost-aware reuse, rolling-only ablation, fresh and fixed controls.
Every timed sequence exports accepted solutions outside per-step timed regions.
The adaptive decision can vary with wall time; certify each actual repetition.
"""
import argparse,json,math,os,re,statistics,subprocess
from pathlib import Path
import anchor_stage10_contact as contact
import anchor_stage12b_nist as nist

PROFILES=("full-ilutp","incomplete-ilutp","ilut","ilu0","ilu1")
MODES=("fresh","fixed","rolling","cost-aware")
HOLDOUT=[(0,0,0),(.1,1,0),(.9,200,.3),(.95,180,.2),(.85,150,.1),
         (0,0,0),(.05,3,0),(.75,300,.4),(.7,220,.3),(.65,160,.1),
         (.1,2,0),(0,0,0),(.6,30,0),(.8,100,.2),(.2,20,.1)]
def parse(log,steps,mode):
    rows=[];summary=[]
    for line in log.read_text().splitlines():
        if line.startswith("COST "):
            d=contact.kv(line)
            r={k:int(d[k]) for k in contact.INTS+["skip","reason","promote","attempted"]}
            r.update({k:float(d[k]) for k in contact.FLOATS})
            rows.append(r)
        elif line.startswith("COST_SUMMARY "):summary.append(contact.kv(line))
    contact.need(len(rows)==steps and len(summary)==1,"incomplete sequence")
    for i,r in enumerate(rows):
        contact.need(r["step"]==i+1 and r["accepted"] in (0,1) and r["api"]==0,"row")
        contact.need(math.isfinite(r["total_s"]) and 0<=r["setup_s"]<=r["total_s"]+1e-8,"time")
        contact.need(0<=r["pilot_iter"]<=20 and 0<=r["continuation_iter"]<=80 and 0<=r["fresh_iter"]<=500,"budget")
        contact.need(all(r[k] in (0,1) for k in ("skip","promote","attempted","fallback","continued")),"flags")
        contact.need(not(r["skip"] and r["attempted"]) and (not r["fallback"] or r["attempted"]),"decision")
        if r["accepted"]:contact.need(math.isfinite(r["residual"]) and 0<=r["residual"]<=1e-8,"physical acceptance")
        if mode=="fresh":contact.need(r["factors"]==1 and not r["attempted"] and not r["skip"],"fresh policy")
        if mode=="fixed":contact.need(r["factors"]==int(i==0)+r["fallback"] and r["attempted"],"fixed policy")
        if mode in ("rolling","cost-aware"):
            contact.need(r["factors"]==int(not r["attempted"] or r["fallback"]),"rolling factor count")
            contact.need(r["promote"]==int(r["accepted"] and r["factors"]>0),"promotion")
        if mode!="cost-aware":contact.need(not r["skip"],"unexpected skip")
    s=summary[0];accepted=sum(r["accepted"] for r in rows);factors=sum(r["factors"] for r in rows)
    contact.need(s["accepted"]==f"{accepted}/{steps}" and int(s["factors"])==factors,"summary")
    for k in ("setup_s","total_s"):
        contact.need(abs(float(s[k])-sum(r[k] for r in rows))<=(steps+2)*1e-9,"summary time")
    fills=[int(l.split()[1]) for l in log.read_text().splitlines() if l.startswith("FACTOR_NNZ ")]
    contact.need(len(fills)==factors,"fill records")
    metrics=dict(accepted=accepted,factors=factors,setup_s=float(s["setup_s"]),total_s=float(s["total_s"]),
        fallbacks=sum(r["fallback"] for r in rows),skips=sum(r["skip"] for r in rows),
        promotions=sum(r["promote"] for r in rows),anchor_iterations=sum(r["pilot_iter"]+r["continuation_iter"] for r in rows),
        fresh_iterations=sum(r["fresh_iter"] for r in rows),peak_factor_entries=max(fills) if fills and all(f>=0 for f in fills) else None)
    return rows,metrics
def main():
    p=argparse.ArgumentParser()
    p.add_argument("bench",type=Path);p.add_argument("root",type=Path)
    p.add_argument("--case",choices=("grid32","grid100","holdout48","bcsstk16","bcsstk17","bcsstk18"),required=True)
    p.add_argument("--source",type=Path);p.add_argument("--reps",type=int,default=3)
    a=p.parse_args();contact.need(a.reps>=3,"minimum three repetitions")
    root=a.root.resolve();root.mkdir();bench=a.bench.resolve()
    if a.case.startswith("bcsstk"):
        args,records,dof,n=nist.generate(a.source.resolve(),root/"input",a.case)
    else:
        w={"grid32":32,"grid100":100,"holdout48":48}[a.case]
        if a.case=="holdout48":contact.CASES=HOLDOUT
        args,records,dof=contact.generate(root/"input",w,w);n=w*w
    profiles=PROFILES if a.case.startswith("grid") or a.case=="holdout48" else PROFILES[1:]
    modes=[(p,m) for p in profiles for m in MODES]
    hashes={str(v):contact.sha(v) for v in [bench,dof,*map(Path,args)]}
    if a.case.startswith("bcsstk"):hashes[str(a.source.resolve())]=contact.sha(a.source.resolve())
    env=os.environ.copy();env.update(OMP_NUM_THREADS="1",OMP_DYNAMIC="FALSE",LC_ALL="C",OPENBLAS_NUM_THREADS="1")
    cpu=min(os.sched_getaffinity(0));steps=len(records)
    (root/"manifest.json").write_text(json.dumps(dict(case=a.case,n=n,steps=steps,reps=a.reps,cpu=cpu,hashes=hashes,records=records,profiles=profiles,
        cost_policy="last observed fresh total; 80% budget; projected slope every five after ten; 1.25 safety factor; two failed pilots -> two fresh steps",
        timing="per-step timer excludes input read and solution export; all timed solutions audited; sequential rotated order"),indent=2))
    measurements=[];certs=[];failures=[]
    for rep in range(a.reps):
        for offset in range(len(modes)):
            profile,mode=modes[(rep+offset)%len(modes)];key=profile+":"+mode
            stem=root/f"{rep+1}-{profile}-{mode}";log=stem.with_suffix(".log");tim=stem.with_suffix(".time")
            audit=Path(str(stem)+"-solutions");audit.mkdir()
            runenv=env.copy();runenv.update(LIS_ANCHOR_PROFILE=profile,LIS_ANCHOR_AUDIT_DIR=str(audit))
            try:
                with log.open("w") as f:
                    q=subprocess.run(["/usr/bin/time","-v","-o",str(tim),"taskset","-c",str(cpu),str(bench),mode,*args],stdout=f,stderr=subprocess.STDOUT,env=runenv,timeout=900)
                stem.with_suffix(".exit").write_text(str(q.returncode))
                rows,m=parse(log,steps,mode)
                contact.need(q.returncode==int(m["accepted"]!=steps),"exit mismatch")
                wanted={audit/f"x-{r['step']:03}.mtx" for r in rows if r["accepted"]}
                contact.need(set(audit.iterdir())==wanted,"audit set")
                for r,rec in zip(rows,records):
                    if not r["accepted"]:continue
                    x=audit/f"x-{r['step']:03}.mtx";h=contact.sha(x)
                    res=contact.certify(Path(rec["a"]),Path(rec["b"]),x,n)
                    contact.need(math.isfinite(res) and res<=1e-8 and contact.sha(x)==h,"independent residual")
                    certs.append(dict(mode=key,rep=rep+1,step=r["step"],residual=res,sha256=h))
                mem=re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)",tim.read_text())
                contact.need(mem is not None,"RSS")
                measurements.append(dict(mode=key,rep=rep+1,rss_kib=int(mem[1]),**m))
                if m["accepted"]!=steps:failures.append(dict(mode=key,rep=rep+1,error=f"accepted {m['accepted']}/{steps}",rejected=[r["step"] for r in rows if not r["accepted"]]))
                print("COST_RUN",a.case,key,rep+1,m["accepted"],flush=True)
            except (RuntimeError,ValueError,KeyError,OSError,subprocess.TimeoutExpired) as e:
                failures.append(dict(mode=key,rep=rep+1,error=str(e)))
                print("COST_FAILURE",failures[-1],flush=True)
    contact.need(all(contact.sha(Path(v))==h for v,h in hashes.items()),"input/binary changed")
    result={}
    for profile,mode in modes:
        key=profile+":"+mode;ms=[v for v in measurements if v["mode"]==key];cs=[v for v in certs if v["mode"]==key]
        ok=len(ms)==a.reps and len(cs)==steps*a.reps and not any(v["mode"]==key for v in failures)
        metrics={k:statistics.median(v[k] for v in ms) for k in ms[0] if k not in ("mode","rep") and all(v[k] is not None for v in ms)} if ms else {}
        if ms:metrics.update(time_min_s=min(v["total_s"] for v in ms),time_max_s=max(v["total_s"] for v in ms))
        result[key]=dict(eligible=ok,completed_runs=len(ms),audits=len(cs),max_residual=max((v["residual"] for v in cs),default=None),**metrics)
    pairs={}
    for profile in profiles:
        f=result[profile+":fresh"]
        if f["eligible"]:
            pairs[profile]={mode:f["total_s"]/result[profile+":"+mode]["total_s"] for mode in MODES[1:] if result[profile+":"+mode]["eligible"]}
    report=dict(case=a.case,result=result,pairs=pairs,failures=failures,measurements=measurements,certifications=certs)
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("COST_RESULT",json.dumps(dict(case=a.case,result=result,pairs=pairs,failures=failures),allow_nan=False),flush=True)
    ok=all(v["eligible"] for v in result.values())
    print("COST_GATE="+("PASS" if ok else "FAIL"),flush=True)
    return 0 if ok else 1
if __name__=="__main__":raise SystemExit(main())
