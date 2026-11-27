#!/usr/bin/env python3
"""Materialize an opt-in native RCI source package from the pinned generators.

Run only in a clean disposable checkout. The resulting tree contains C sources
and headers, not a build-time source mutation or a Python/Node dependency.
"""
from pathlib import Path
import difflib
import os
import re
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
Q=ROOT/"qualification/rci-restore"
P=ROOT/"src/rci_native"
I=ROOT/"include"
os.chdir(ROOT)
def run(*args,**kw):
    return subprocess.run(args,check=True,**kw)
def git(*args):
    return subprocess.check_output(["git",*args],text=True).strip()
assert not git("status","--porcelain","--untracked-files=no"), "clean checkout required"
assert not P.exists(), "already materialized"
original={ROOT/p:(ROOT/p).read_text() for p in git("ls-files","src").splitlines() if p.endswith(".c")}
def generate(complex_mode,full):
    env=dict(os.environ)
    env.pop("LIS_RCI_QUAL_COMPLEX",None)
    if complex_mode: env["LIS_RCI_QUAL_COMPLEX"]="1"
    run(sys.executable,str(Q/("qualify_full_native.py" if full else "qualify_all_native.py")),env=env)
    changed={p:p.read_text() for p,s in original.items() if p.read_text()!=s}
    inc={p.name:p.read_text() for p in Q.glob("*continuation.inc")}
    for p in changed:p.write_text(original[p])
    for p in Q.glob("*continuation.inc"):
        if not git("ls-files",str(p.relative_to(ROOT))):p.unlink()
    return changed,inc
real_src,real_inc=generate(False,True)
_,complex_inc=generate(True,False)
P.mkdir()
# Wrap only generator changes, preserving the default source branch verbatim.
# For full lifecycle only real-double is qualified.
for p,s in real_src.items():
    cond="defined(USE_NATIVE_RCI)"
    if p.name in ("lis_solver.c","lis_esolver.c"):
        cond+=" && !defined(_COMPLEX) && !defined(_LONG__DOUBLE)"
    s=s.replace("../../qualification/rci-restore/","../rci_native/")
    a=original[p].splitlines(keepends=True);b=s.splitlines(keepends=True)
    out=[]
    for op,i,j,k,l in difflib.SequenceMatcher(None,a,b,autojunk=False).get_opcodes():
        if op=="equal":out+=a[i:j]
        else:
            out+=["\n#if "+cond+"\n",*b[k:l]]
            if a[i:j]:out+=["\n#else\n",*a[i:j]]
            out+=["\n#endif /* native RCI */\n"]
    text="".join(out)
    if p.name=="lis_krylov.c":
        text='#ifdef HAVE_CONFIG_H\n#include "lis_config.h"\n#endif\n'+text
    p.write_text(text)
# Generated bodies include the native source attribution and license.
notice="/* Generated continuation; native LIS source copyright and license apply.\n"
notice+=" * Sources: "+", ".join(sorted(str(p.relative_to(ROOT)) for p in real_src))+"\n */\n"
for name,content in sorted(real_inc.items()):
    other=complex_inc.get(name,content)
    if other==content:
        (P/name).write_text(notice+content)
    else:
        stem=name.removesuffix(".inc")
        (P/(stem+"_real.inc")).write_text(notice+content)
        (P/(stem+"_complex.inc")).write_text(notice+other)
        (P/name).write_text('#ifdef _COMPLEX\n#include "'+stem+'_complex.inc"\n#else\n#include "'+stem+'_real.inc"\n#endif\n')
# Separate the shared status and operation constants from the legacy raw API.
legacy=I/"lis_rci.h"
s=legacy.read_text()
m=re.search(r"typedef enum \{ LIS_RCI_DONE=0,.*?\} lis_rci_status;",s,re.S)
assert m
(I/"lis_rci_status.h").write_text("#ifndef LIS_RCI_STATUS_H\n#define LIS_RCI_STATUS_H\n"+m[0]+"\n#endif\n")
legacy.write_text(s[:m.start()]+'#include "lis_rci_status.h"'+s[m.end():])
restore=Q/"lis_rci_restore.h"
s=restore.read_text()
m=re.search(r"enum lis_rci_opaque_operation \{.*?\};",s,re.S)
assert m
(I/"lis_rci_native_protocol.h").write_text('#ifndef LIS_RCI_NATIVE_PROTOCOL_H\n#define LIS_RCI_NATIVE_PROTOCOL_H\n#include "lis_rci_status.h"\n'+m[0]+"\n#endif\n")
restore.write_text(s[:m.start()]+'#include "lis_rci_native_protocol.h"'+s[m.end():])
private_names=("lis_rci_native_linear_restore","lis_rci_native_gcro_restore","lis_rci_native_fgcro_restore",
               "lis_rci_native_feast","lis_rci_krylov_restore","lis_rci_native_session","lis_rci_native_full",
               "lis_rci_native_full_feast","lis_rci_native_full_internal")
for stem in private_names:
    for suffix in (".c",".h"):
        f=Q/(stem+suffix)
        if not f.exists():continue
        s=f.read_text().replace('#include "lis_rci_restore.h"','#include "lis_rci_native_protocol.h"')
        if f.name in ("lis_rci_native_session.h","lis_rci_native_full.h"):
            (I/f.name).write_text(s.replace("not installed","optional installed"))
        else:(P/f.name).write_text(s)
# Feature is independently opt-in and leaves old RCI switches/ABI intact.
ac=ROOT/"configure.ac"
s=ac.read_text()
anchor='AS_CASE([$enable_rci], [yes|no], [], [AC_MSG_ERROR([--enable-rci accepts yes or no])])'
assert s.count(anchor)==1
s=s.replace(anchor,anchor+'''
AC_ARG_ENABLE([native-rci], AS_HELP_STRING([--enable-native-rci],[enable experimental shared native RCI cores [[default=no]]]), [], [enable_native_rci=no])
AS_CASE([$enable_native_rci], [yes|no], [], [AC_MSG_ERROR([--enable-native-rci accepts yes or no])])
''')
anchor='AM_CONDITIONAL([ENABLE_RCI], [test "x$enable_rci" = xyes])'
assert s.count(anchor)==1
s=s.replace(anchor,'''
if test "x$enable_native_rci" = xyes; then
  if test "x$enable_mpi" = xyes || test "x$enable_quad" = xyes; then
    AC_MSG_ERROR([native RCI requires non-MPI LIS without double-double solver mode])
  fi
  AC_DEFINE([USE_NATIVE_RCI], [1], [Define to 1 for the experimental shared native RCI cores.])
fi
AM_CONDITIONAL([ENABLE_NATIVE_RCI], [test "x$enable_native_rci" = xyes])
AM_CONDITIONAL([ENABLE_NATIVE_RCI_FULL], [test "x$enable_native_rci" = xyes && test "x$enable_complex" != xyes && test "x$enable_longdouble" != xyes])
'''+anchor)
anchor='if test "x$enable_rci" = xyes; then\n  AC_MSG_CHECKING([for nonfinite-aware'
assert s.count(anchor)==1
s=s.replace(anchor,'if test "x$enable_rci" = xyes || test "x$enable_native_rci" = xyes; then\n  AC_MSG_CHECKING([for nonfinite-aware')
ac.write_text(s)
im=I/"Makefile.am"
s=im.read_text()
s+='''
if ENABLE_NATIVE_RCI
nobase_include_HEADERS += lis_rci_status.h lis_rci_native_protocol.h lis_rci_native_session.h
else
if ENABLE_RCI
nobase_include_HEADERS += lis_rci_status.h
endif
endif
if ENABLE_NATIVE_RCI_FULL
nobase_include_HEADERS += lis_rci_native_full.h
endif
'''
im.write_text(s)
am=ROOT/"Makefile.am"
s=am.read_text()
s+="\nEXTRA_DIST += "+" \\\n    ".join(str(p.relative_to(ROOT)) for p in sorted(P.iterdir()))+" \\\n    doc/native-rci.md\n"
am.write_text(s)
(ROOT/"doc/native-rci.md").write_text("""# Experimental shared native RCI

Configure with \`--enable-native-rci\` (default off). The legacy
\`--enable-rci\` flag and legacy raw-vector API remain independent.
The optional native session API uses LIS vector handles, except Krylov
which uses its original real-double array ABI. See lis_rci_native_session.h
for ownership, request versions, split-complex shift conventions and status.
The full lifecycle header is installed in real-double builds only.
Its documented restrictions are checked before native setup.
MPI, double-double solver mode, and fast/finite-only math are unsupported.
Complex and extended scalar builds support the prepared native session.
Krylov itself stays real double. Do not use private restore/test symbols.
This is experimental, not a stable new ABI. Generated continuations are
materialized from pinned native LIS source by the qualification generators;
normal builds require neither Python nor Node.
""")
print("MATERIALIZED",len(real_src),"native source dispatch files,",len(list(P.iterdir())),"private files")
