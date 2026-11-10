#!/usr/bin/env bash
# Synthetic only; checks read-only 12/14/18 geometric overlap calculations.
set -euo pipefail
here="$(cd -- "$(dirname -- "$0")" && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT
python3 - "$tmp" <<'PY'
from pathlib import Path
import math, sys
p=Path(sys.argv[1]); n=24
cols=[[float(i==j) for i in range(n)] for j in range(20)]
def write(name, vectors):
    content='%%MatrixMarket matrix array real general\n% synthetic only\n'+f'{n} {len(vectors)}\n'
    content+=''.join('%.17g\n'%x for col in vectors for x in col)
    (p/name).write_text(content, encoding='ascii')
write('base.mtx', cols[:12])
write('offline.mtx', cols[:18])
write('contained.mtx', cols[:14])
write('disjoint.mtx', cols[:12]+cols[18:20])
write('mixed.mtx', cols[:12]+[[(a+b)/math.sqrt(2) for a,b in zip(cols[i],cols[i+6])] for i in (12,13)])
write('wrong_prefix.mtx', [cols[18]]+cols[1:14])
write('dependent.mtx', cols[:13]+[cols[12]])
PY
python3 - "$here" "$tmp" <<'PY'
import pathlib, subprocess, sys
here=pathlib.Path(sys.argv[1]); tmp=pathlib.Path(sys.argv[2])
for suffix, expected in [('contained',1.0),('disjoint',0.0),('mixed',0.5)]:
    run=subprocess.run([sys.executable, str(here/'mfamg_v8_auto14_overlap.py'),
       '--baseline', str(tmp/'base.mtx'), '--generated', str(tmp/(suffix+'.mtx')),
       '--offline-control', str(tmp/'offline.mtx')],capture_output=True,text=True)
    assert run.returncode == 0, run.stderr
    report=dict(line.split('=',1) for line in run.stdout.splitlines())
    assert report['V8_AUTO14_OVERLAP']=='OK'
    assert abs(float(report['V8_AUTO14_OVERLAP_MEAN_SQUARED_PROJECTION'])-expected)<1e-10, report
    assert abs(float(report['V8_AUTO14_OVERLAP_PRINCIPAL_COSINE_MIN']) - expected**0.5)<1e-10, report
    assert abs(float(report['V8_AUTO14_OVERLAP_PRINCIPAL_COSINE_MAX']) - expected**0.5)<1e-10, report
    assert str(tmp) not in run.stdout, 'exposed private path'
for suffix in ('wrong_prefix','dependent'):
    run=subprocess.run([sys.executable, str(here/'mfamg_v8_auto14_overlap.py'),
       '--baseline',str(tmp/'base.mtx'), '--generated',str(tmp/(suffix+'.mtx')),
       '--offline-control',str(tmp/'offline.mtx')],capture_output=True,text=True)
    assert run.returncode == 2, run.stdout
    assert 'V8_AUTO14_OVERLAP=BLOCKED' in run.stderr,run.stderr
    assert str(tmp) not in run.stderr, 'exposed private path'
# Fail closed if a privately supplied mode file changes while being inspected.
sys.path.insert(0, str(here))
import mfamg_v8_auto14_overlap as audit
original_reader = audit.read_modes
changed = [False]
def mutating_reader(path):
    result = original_reader(path)
    if str(path).endswith('contained.mtx') and not changed[0]:
        changed[0] = True
        with (tmp/'base.mtx').open('ab') as stream:
            stream.write(b'% concurrent modification\n')
    return result
audit.read_modes = mutating_reader
try:
    try:
        audit.analyze(tmp/'base.mtx', tmp/'contained.mtx', tmp/'offline.mtx')
    except ValueError as e:
        assert 'changed during post-hoc overlap audit' in str(e), str(e)
    else:
        raise AssertionError('a concurrently changed mode input was accepted')
finally:
    audit.read_modes = original_reader
assert changed[0]
PY
echo V8_AUTO14_OVERLAP_SMOKE_OK
