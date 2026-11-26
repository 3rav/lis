#!/usr/bin/env python3
"""Same-runner three-repeat benchmark; retain rejected methods and raw data.
UMFPACK reuses symbolic only for exact observed CSC patterns, no future input.
"""
import argparse,json,math,os,re,statistics,subprocess
from pathlib import Path
import anchor_stage10_contact as contact
import anchor_stage12b_nist as nist

ANCHOR=("fresh","pilot","pilot-restart-matched","pilot-single-pass")
OTHER=("fresh-ilut","fresh-ilutp-incomplete","fresh-ilu0","fresh-ilu1","fresh-mfamg","umfpack-fresh","umfpack-symbolic-cache")
MODES=ANCHOR+OTHER

def parse(log,mode,steps):
    if mode in ANCHOR:
        rows,s=contact.parse(log,mode,steps=steps)
        fills=[int(l.split()[1]) for l in log.read_text().splitlines() if l.startswith("FACTOR_NNZ ")]
        contact.need(len(fills)==int(s["factors"]),"missing factor fill record")
        sig=contact.signature(rows)
        metrics=dict(accepted=steps,setup_s=float(s["setup_s"]),total_s=float(s["total_s"]),
            factors=int(s["factors"]),symbolic=None,fallbacks=sum(r["fallback"] for r in rows),
            iterations=sum(min(r["pilot_iter"],20)+min(r["continuation_iter"],80)+min(r["fresh_iter"],100) for r in rows),
            peak_factor_entries=max(fills) if all(f>=0 for f in fills) else None)
        return rows,sig,metrics
    rows=[contact.kv(l) for l in log.read_text().splitlines() if l.startswith("BASE ")]
    contact.need(len(rows)==steps,"incomplete baseline")
    for i,r in enumerate(rows):
        for k in ("step","accepted","api","status","iterations","factors","symbolic","fill"):r[k]=int(r[k])
        for k in ("residual","setup_s","total_s"):r[k]=float(r[k])
        contact.need(r["step"]==i+1 and r["accepted"] in (0,1),"bad baseline row")
        contact.need(0<=r["setup_s"]<=r["total_s"]+1e-8 and math.isfinite(r["total_s"]),"bad timing")
        if r["accepted"]:
            contact.need(r["api"]==0 and math.isfinite(r["residual"]) and 0<=r["residual"]<=1e-8,"bad acceptance")
    sig=[tuple(r[k] for k in ("step","accepted","api","status","iterations","factors","symbolic","fill"))+
         (r["residual"] if math.isfinite(r["residual"]) else None,) for r in rows]
    fills=[r["fill"] for r in rows if r["factors"]]
    metrics=dict(accepted=sum(r["accepted"] for r in rows),setup_s=sum(r["setup_s"] for r in rows),
        total_s=sum(r["total_s"] for r in rows),factors=sum(r["factors"] for r in rows),
        symbolic=sum(r["symbolic"] for r in rows) if mode.startswith("umfpack") else None,
        fallbacks=0,iterations=sum(r["iterations"] for r in rows) if not mode.startswith("umfpack") else None,
        peak_factor_entries=max(fills) if fills and all(f>=0 for f in fills) else None)
    return rows,sig,metrics

def main():
    p=argparse.ArgumentParser()
    p.add_argument("anchor",type=Path);p.add_argument("baseline",type=Path);p.add_argument("root",type=Path)
    p.add_argument("--case",choices=("grid32","grid100","bcsstk16"),required=True)
    p.add_argument("--source",type=Path);p.add_argument("--reps",type=int,default=3)
    a=p.parse_args();contact.need(a.reps>=3,"minimum three repetitions")
    root=a.root.resolve();root.mkdir()
    if a.case=="bcsstk16":
        source=a.source.resolve();args,records,dof,n=nist.generate(source,root/"input",a.case)
    else:
        w=32 if a.case=="grid32" else 100
        args,records,dof=contact.generate(root/"input",w,w);n=w*w
    steps=len(records)
    bins=[a.anchor.resolve(),a.baseline.resolve()]
    files=[*bins,dof,*map(Path,args)]
    if a.case=="bcsstk16":files.append(source)
    hashes={str(v):contact.sha(v) for v in files}
    env=os.environ.copy();env.update(OMP_NUM_THREADS="1",OMP_DYNAMIC="FALSE",LC_ALL="C",OPENBLAS_NUM_THREADS="1",BLIS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    env.pop("LIS_ANCHOR_AUDIT_DIR",None);cpu=min(os.sched_getaffinity(0))
    manifest=dict(case=a.case,n=n,steps=steps,reps=a.reps,cpu=cpu,hashes=hashes,records=records,
        native_fresh_budget=500,anchor_budget=100,anchor_fallback_budget=100,
        direct_symbolic="exact-pattern online cache, numeric factorization always repeated",
        fill_definition="maximum entries in one factor; L+U+one diagonal, null if unavailable",
        timing="input load excluded; process RSS includes input and library; same runner, sequential, rotated order")
    (root/"manifest.json").write_text(json.dumps(manifest,indent=2))
    measurements=[];certs=[];failures=[];signatures={}
    for rep in range(a.reps):
        for kind in ("timing","audit"):
            for offset in range(len(MODES)):
                mode=MODES[(rep+offset)%len(MODES)]
                stem=root/f"{kind}-{rep+1}-{mode}"
                log=stem.with_suffix(".log");tim=stem.with_suffix(".time")
                audit=Path(str(stem)+"-solutions");runenv=env.copy()
                if kind=="audit":audit.mkdir();runenv["LIS_ANCHOR_AUDIT_DIR"]=str(audit)
                try:
                    exe=bins[0 if mode in ANCHOR else 1]
                    with log.open("w") as f:
                        q=subprocess.run(["/usr/bin/time","-v","-o",str(tim),"taskset","-c",str(cpu),str(exe),mode,*args],
                            stdout=f,stderr=subprocess.STDOUT,env=runenv,timeout=900)
                    stem.with_suffix(".exit").write_text(str(q.returncode))
                    rows,sig,m=parse(log,mode,steps)
                    contact.need(q.returncode==(0 if m["accepted"]==steps else 1),"exit/acceptance mismatch")
                    if kind=="timing":
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
                    print("FAIR_RUN",a.case,rep+1,kind,mode,m["accepted"],flush=True)
                except (RuntimeError,ValueError,KeyError,OSError,subprocess.TimeoutExpired) as e:
                    failures.append(dict(mode=mode,rep=rep+1,kind=kind,error=str(e)))
                    print("FAIR_FAILURE",failures[-1],flush=True)
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
    print("FAIR_RESULT",json.dumps(dict(case=a.case,result=result,failures=failures),allow_nan=False),flush=True)
    # A rejected competing method is reported, not silently discarded and not
    # itself an infrastructure failure. Require every anchor and direct mode.
    ok=all(result[m]["eligible"] for m in (*ANCHOR,"umfpack-fresh","umfpack-symbolic-cache"))
    print("FAIR_GATE="+("PASS" if ok else "FAIL"),flush=True)
    return 0 if ok else 1

if __name__=="__main__":raise SystemExit(main())
