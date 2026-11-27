#!/usr/bin/env python3
"""Paired native audit; instrumented timings are diagnostic only."""
import argparse, csv, json, math, os, subprocess
from pathlib import Path
import anchor_stage10_contact as contact
import anchor_stage12b_nist as nist

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("baseline",type=Path); p.add_argument("instrumented",type=Path)
    p.add_argument("source",type=Path); p.add_argument("root",type=Path)
    p.add_argument("--name",choices=nist.EXPECTED,required=True)
    a=p.parse_args(); root=a.root.resolve(); root.mkdir()
    args,records,dof,n=nist.generate(a.source.resolve(),root/"input",a.name)
    inputs=[Path(x) for x in args]+[dof,a.source.resolve()]
    hashes={str(x):contact.sha(x) for x in inputs}
    (root/"manifest.json").write_text(json.dumps({"name":a.name,"records":records,"hashes":hashes},indent=2))
    env=os.environ.copy(); env.update(OMP_NUM_THREADS="1",OMP_DYNAMIC="FALSE",LC_ALL="C",OPENBLAS_NUM_THREADS="1")
    env.pop("LIS_ANCHOR_FACTOR_DIAGNOSTICS",None); cpu=min(os.sched_getaffinity(0))
    results=[]; factors=[]
    for mode in ("fresh","pilot-restart-matched","pilot-single-pass"):
        previous=None
        for label,exe in (("baseline",a.baseline),("disabled",a.instrumented),("enabled",a.instrumented)):
            stem=root/(mode+"-"+label); audit=Path(str(stem)+"-solutions"); audit.mkdir()
            runenv=dict(env,LIS_ANCHOR_AUDIT_DIR=str(audit))
            if label=="enabled": runenv["LIS_ANCHOR_FACTOR_DIAGNOSTICS"]="1"
            log=stem.with_suffix(".log")
            with log.open("w") as f:
                subprocess.run(["/usr/bin/time","-v","-o",str(stem.with_suffix(".time")),"taskset","-c",str(cpu),str(exe.resolve()),mode,*args],stdout=f,stderr=subprocess.STDOUT,env=runenv,check=True,timeout=900)
            rows,summary=contact.parse(log,mode,steps=len(records))
            signature=contact.signature(rows); cert=[]
            pending=[]; step=0
            for line in log.read_text().splitlines():
                if line.startswith("FACTOR "):
                    pending.append(dict(item.split("=",1) for item in line.split()[1:]))
                elif line.startswith("STEP "):
                    step+=1; expected=rows[step-1]["factors"] if label=="enabled" else 0
                    need(len(pending)==expected,"factor count mismatch")
                    for f in pending:
                        need(int(f["rc"])==0 and int(f["stats_rc"])==0,"factor diagnostic error")
                        for k in ("n","nnz_a","nnz_l","nnz_u","swaps","rss_peak_kib"):
                            f[k]=int(f[k]); need(f[k]>=0,k)
                        for k in ("min_pivot","max_pivot","fill_ratio","setup_s"):
                            f[k]=float(f[k]); need(math.isfinite(f[k]) and f[k]>=0,k)
                        need(f["n"]==n and f["nnz_a"]>0 and 0<f["min_pivot"]<=f["max_pivot"],"factor range")
                        need(math.isclose(f["fill_ratio"],(f["nnz_l"]+f["nnz_u"]+n)/f["nnz_a"],rel_tol=1e-12),"fill mismatch")
                        factors.append(dict(matrix=a.name,mode=mode,step=step,**f))
                    pending=[]
            need(step==len(records) and not pending,"incomplete diagnostics")
            need(set(audit.iterdir())=={audit/f"x-{k:03}.mtx" for k in range(1,len(records)+1)},"solution set")
            for rec in records:
                x=audit/f"x-{rec['step']:03}.mtx"; h=contact.sha(x)
                residual=contact.certify(Path(rec["a"]),Path(rec["b"]),x,n)
                need(math.isfinite(residual) and residual<=1e-8 and h==contact.sha(x),"physical residual")
                cert.append({"step":rec["step"],"sha256":h,"residual":residual})
            comparison=(signature,[v["sha256"] for v in cert])
            if previous is None: previous=comparison
            else: need(previous==comparison,"numerical or vector mismatch")
            results.append(dict(mode=mode,variant=label,summary=summary,certifications=cert))
            print("DIAGNOSTIC_RUN",a.name,mode,label,"PASS",flush=True)
    need(all(contact.sha(Path(k))==v for k,v in hashes.items()),"input changed")
    with (root/"factors.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=list(factors[0])); w.writeheader(); w.writerows(factors)
    (root/"report.json").write_text(json.dumps({"status":"PASS","runs":results,"factors":factors},indent=2,allow_nan=False))
    print("FACTOR_DIAGNOSTIC_GATE=PASS",a.name,"factor_records",len(factors),flush=True)
if __name__=="__main__": main()
