#!/usr/bin/env python3
"""Same-runner three-repeat benchmark; retain rejected methods and raw data.
UMFPACK reuses symbolic only for exact observed CSC patterns, no future input.
"""
import argparse,json,math,os,re,statistics,subprocess
from pathlib import Path
import anchor_stage10_contact as contact
import anchor_stage12b_nist as nist

PROFILES=("full-ilutp","incomplete-ilutp","ilut","ilu0","ilu1")
POLICIES=("fresh","pilot-single-pass")

def parse(log,mode,steps):
    profile,policy=mode.split(":")
    rows=[];summary=[]
    for line in log.read_text().splitlines():
        if line.startswith("STEP "):
            d=contact.kv(line);r={k:int(d[k]) for k in contact.INTS}
            r.update({k:float(d[k]) for k in contact.FLOATS});rows.append(r)
        elif line.startswith("SUMMARY "):summary.append(contact.kv(line))
    contact.need(len(rows)==steps and len(summary)==1,f"incomplete run ({len(rows)}/{steps}); see raw log")
    for i,r in enumerate(rows):
        contact.need(r["step"]==i+1 and r["accepted"] in (0,1),"row order/acceptance")
        if r["accepted"]:
            contact.need(r["api"]==0 and math.isfinite(r["residual"]) and 0<=r["residual"]<=1e-8,"physical acceptance")
        contact.need(math.isfinite(r["total_s"]) and 0<=r["setup_s"]<=r["total_s"]+1e-6,"time")
        contact.need(r["fallback"] in (0,1) and r["continued"] in (0,1),"flags")
        contact.need(0<=r["pilot_iter"]<=20 and 0<=r["continuation_iter"]<=80 and 0<=r["fresh_iter"]<=500,"budget")
        if policy=="fresh":
            contact.need(r["factors"]==1 and not any(r[k] for k in ("pilot_iter","continuation_iter","continued","fallback")),"fresh factor count")
        else:
            contact.need(r["factors"]==int(i==0)+r["fallback"],"retained factor count")
            if r["continued"]:contact.need(r["pilot_iter"]==20 and 0<=r["ratio"]<=.25,"gate")
            contact.need(r["fallback"] or r["fresh_iter"]==0,"unexpected fresh solve")
    s=summary[0];accepted=sum(r["accepted"] for r in rows);factors=sum(r["factors"] for r in rows)
    contact.need(s["accepted"]==f"{accepted}/{steps}" and int(s["factors"])==factors,"summary")
    for k in ("total_s","setup_s"):
        contact.need(abs(float(s[k])-sum(r[k] for r in rows))<=(steps+2)*1e-6,"summary time")
    fills=[int(l.split()[1]) for l in log.read_text().splitlines() if l.startswith("FACTOR_NNZ ")]
    contact.need(len(fills)==factors,"fill records")
    sig=[tuple(r[k] for k in contact.INTS)+tuple(r[k] if math.isfinite(r[k]) else None for k in ("ratio","residual")) for r in rows]
    metrics=dict(accepted=accepted,setup_s=float(s["setup_s"]),total_s=float(s["total_s"]),factors=factors,
        fallbacks=sum(r["fallback"] for r in rows),
        anchor_iterations=sum(r["pilot_iter"]+r["continuation_iter"] for r in rows),
        fresh_iterations=sum(r["fresh_iter"] for r in rows),
        peak_factor_entries=max(fills) if fills and all(f>=0 for f in fills) else None)
    return rows,sig,metrics

def main():
    p=argparse.ArgumentParser()
    p.add_argument("anchor",type=Path);p.add_argument("root",type=Path)
    p.add_argument("--case",choices=("grid32","grid100","bcsstk16","bcsstk17","bcsstk18"),required=True)
    p.add_argument("--source",type=Path);p.add_argument("--reps",type=int,default=3)
    a=p.parse_args();contact.need(a.reps>=3,"minimum three repetitions")
    root=a.root.resolve();root.mkdir()
    if a.case.startswith("bcsstk"):
        source=a.source.resolve();args,records,dof,n=nist.generate(source,root/"input",a.case)
    else:
        w=32 if a.case=="grid32" else 100
        args,records,dof=contact.generate(root/"input",w,w);n=w*w
    steps=len(records)
    exe=a.anchor.resolve()
    files=[exe,dof,*map(Path,args)]
    if a.case.startswith("bcsstk"):files.append(source)
    # Full-fill control is already expensive on these two matrices; compare
    # only the four cheaper profiles on bcsstk17 and bcsstk18.
    profiles=PROFILES[1:] if a.case in ("bcsstk17","bcsstk18") else PROFILES
    MODES=tuple(f"{p}:{m}" for p in profiles for m in POLICIES)
    hashes={str(v):contact.sha(v) for v in files}
    env=os.environ.copy();env.update(OMP_NUM_THREADS="1",OMP_DYNAMIC="FALSE",LC_ALL="C",OPENBLAS_NUM_THREADS="1",BLIS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    env.pop("LIS_ANCHOR_AUDIT_DIR",None);cpu=min(os.sched_getaffinity(0))
    manifest=dict(case=a.case,n=n,steps=steps,reps=a.reps,cpu=cpu,hashes=hashes,records=records,
        native_fresh_budget=500,anchor_budget=100,anchor_fallback_budget=500,
        profiles=profiles,policies=POLICIES,
        fill_definition="maximum entries in one factor; L+U+one diagonal, null if unavailable",
        timing="input load excluded; process RSS includes input and library; same runner, sequential, rotated order")
    (root/"manifest.json").write_text(json.dumps(manifest,indent=2))
    measurements=[];certs=[];failures=[];signatures={}
    for rep in range(a.reps):
        for kind in ("timing","audit"):
            for offset in range(len(MODES)):
                mode=MODES[(rep+offset)%len(MODES)]
                stem=root/f"{kind}-{rep+1}-{mode.replace(chr(58),chr(45))}"
                log=stem.with_suffix(".log");tim=stem.with_suffix(".time")
                audit=Path(str(stem)+"-solutions");runenv=env.copy()
                profile,policy=mode.split(":")
                runenv["LIS_ANCHOR_PROFILE"]=profile
                if kind=="audit":audit.mkdir();runenv["LIS_ANCHOR_AUDIT_DIR"]=str(audit)
                try:
                    with log.open("w") as f:
                        q=subprocess.run(["/usr/bin/time","-v","-o",str(tim),"taskset","-c",str(cpu),str(exe),policy,*args],
                            stdout=f,stderr=subprocess.STDOUT,env=runenv,timeout=900)
                    stem.with_suffix(".exit").write_text(str(q.returncode))
                    rows,sig,m=parse(log,mode,steps)
                    contact.need(q.returncode==(0 if m["accepted"]==steps else 1),"exit/acceptance mismatch")
                    if kind=="timing":
                        contact.need(signatures.setdefault(("repeat",mode),sig)==sig,"repeat decision mismatch")
                        signatures[rep,mode]=sig
                        mem=re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)",tim.read_text())
                        contact.need(mem is not None,"missing RSS")
                        measurements.append(dict(mode=mode,rep=rep+1,rss_kib=int(mem[1]),**m))
                    else:
                        contact.need(signatures.get((rep,mode))==sig,"paired decision mismatch")
                        wanted={audit/f"x-{r['step']:03}.mtx" for r in rows if r["accepted"]}
                        contact.need(set(audit.iterdir())==wanted,"wrong audit vector set")
                        for r,rec in zip(rows,records):
                            if not r["accepted"]:continue
                            x=audit/f"x-{r['step']:03}.mtx";h=contact.sha(x)
                            res=contact.certify(Path(rec["a"]),Path(rec["b"]),x,n)
                            contact.need(math.isfinite(res) and res<=1e-8 and contact.sha(x)==h,"physical residual rejected")
                            certs.append(dict(mode=mode,rep=rep+1,step=r["step"],residual=res,sha256=h))
                    if m["accepted"]!=steps:failures.append(dict(mode=mode,rep=rep+1,kind=kind,error=f"accepted {m['accepted']}/{steps}",rejected=[r["step"] for r in rows if not r["accepted"]]))
                    print("REUSE_RUN",a.case,rep+1,kind,mode,m["accepted"],flush=True)
                except (RuntimeError,ValueError,KeyError,OSError,subprocess.TimeoutExpired) as e:
                    failures.append(dict(mode=mode,rep=rep+1,kind=kind,error=str(e)))
                    print("REUSE_FAILURE",failures[-1],flush=True)
    contact.need(all(contact.sha(Path(v))==h for v,h in hashes.items()),"input/binary changed")
    result={}
    for mode in MODES:
        ms=[v for v in measurements if v["mode"]==mode];cs=[v for v in certs if v["mode"]==mode]
        ok=len(ms)==a.reps and len(cs)==steps*a.reps and not any(v["mode"]==mode for v in failures)
        metrics={}
        if ms:
            for k in ms[0]:
                if k not in ("mode","rep") and all(v[k] is not None for v in ms):
                    metrics[k]=statistics.median(v[k] for v in ms)
            metrics["time_min_s"]=min(v["total_s"] for v in ms);metrics["time_max_s"]=max(v["total_s"] for v in ms)
        result[mode]=dict(eligible=ok,completed_timing_runs=len(ms),audits=len(cs),max_residual=max((v["residual"] for v in cs),default=None),**metrics)
    report=dict(case=a.case,result=result,failures=failures,measurements=measurements,certifications=certs)
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    pairs={}
    for profile in profiles:
        f=result[profile+":fresh"];r=result[profile+":pilot-single-pass"]
        if f["eligible"] and r["eligible"]:
            pairs[profile]=dict(speedup=f["total_s"]/r["total_s"],fresh_s=f["total_s"],reuse_s=r["total_s"],
                fresh_factors=f["factors"],reuse_factors=r["factors"],fallbacks=r["fallbacks"])
    report["pairs"]=pairs
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("REUSE_RESULT",json.dumps(dict(case=a.case,result=result,pairs=pairs,failures=failures),allow_nan=False),flush=True)
    # An unqualified profile is never promoted; keep evidence of all failures.
    ok=all(result[m]["eligible"] for m in MODES)
    print("REUSE_GATE="+("PASS" if ok else "FAIL"),flush=True)
    return 0 if ok else 1

if __name__=="__main__":raise SystemExit(main())
