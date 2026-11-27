#!/usr/bin/env python3
"""Alternate old/new libraries, verify deterministic public client signatures.
Timings are descriptive CI observations, not a release-performance SLA."""
import os,subprocess,sys,statistics,json
from pathlib import Path
exe,baseline,integrated,output=sys.argv[1:]
data={"baseline":{},"integrated":{}}
signatures={}
for rep in range(5):
    for name,lib in ([("baseline",baseline),("integrated",integrated)] if rep%2==0 else [("integrated",integrated),("baseline",baseline)]):
        env=dict(os.environ,LD_LIBRARY_PATH=lib)
        link=subprocess.check_output(["ldd",exe],env=env,text=True)
        assert lib in link,link
        stdout=subprocess.check_output([exe],env=env,text=True)
        assert "BASELINE CLIENT ABI INFO MATCH" in stdout
        seen=set()
        for row in stdout.splitlines():
            if not row.startswith("BENCH "):continue
            _,method,n,count,signature,seconds=row.split()
            key=method+":"+n;assert key not in seen;seen.add(key)
            assert float(seconds)>0
            if key in signatures:assert signatures[key]==signature,(name,key,signatures[key],signature)
            else:signatures[key]=signature
            data[name].setdefault(key,[]).append(float(seconds))
        assert len(seen)==14,seen
report=[]
for key in sorted(signatures):
    old=data["baseline"][key];new=data["integrated"][key]
    assert len(old)==len(new)==5
    ratio=statistics.median(new)/statistics.median(old)
    row=dict(case=key,baseline_seconds=old,integrated_seconds=new,median_ratio=ratio,signature=signatures[key])
    report.append(row);print("PERF",key,"baseline",statistics.median(old),"integrated",statistics.median(new),"ratio",ratio,"REVIEW" if ratio>1.20 else "OK")
Path(output).write_text(json.dumps(report,indent=2)+"\n")
print("BASELINE CLIENT NUMERICAL SIGNATURES PASS")
