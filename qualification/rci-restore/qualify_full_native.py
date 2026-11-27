#!/usr/bin/env python3
"""Lower the actual LIS solve setup/cleanup around the unified core.
Preserve classic routines unchanged as the independent public reference.
"""
import re,subprocess,sys
import generate_linear_continuation as g
subprocess.run([sys.executable,str(g.DIR/"qualify_unified_native.py")],check=True)
g.REAL_ONLY=True
p="src/solver/lis_solver.c";e="src/esolver/lis_esolver.c"
sources={p:g.get_source(p,"4feb18f8dbe399fee5cc3848b4edb77ce9a060d8"),e:g.get_source(e,"3ba869be3d881ce93ff46d418b3cdcbf06978243")}
specs=[("outer","lis_solve",p),("kernel","lis_solve_kernel",p),("eigen","lis_gesolve_once",e)]
frames={}
for name,func,path in specs:
    par,body=g.extract(g.strip_comments(sources[path]),func)
    body=g.no_quad(body)
    fields={}
    for arg in g.args(par):
        m=re.fullmatch(r"(.*?)\s*([A-Za-z_]\w*)",arg)
        assert m,arg
        fields[m[2]]=(m[1].strip(),"")
    # Keep nested synchronous locals on their original C scope. Only the
    # leading automatic declarations can survive the two audited yield sites.
    decl=re.compile(r"^[ \t]*(LIS_Comm|LIS_MATRIX|LIS_VECTOR|LIS_PRECON|LIS_SOLVER|LIS_ESOLVER|LIS_SCALAR|LIS_REAL|LIS_INT|double|char)\s+([^;]+);",re.M)
    head,sep,tail=body.partition("LIS_DEBUG_FUNC_IN;")
    assert sep
    def hoist(m):
        out=[]
        for arg in g.args(m[2]):
            n=re.fullmatch(r"(\**)\s*(\w+)(\[\d+\])?\s*(?:=\s*(.*))?",arg,re.S)
            assert n,(name,arg)
            assert n[2] not in fields,(name,n[2])
            fields[n[2]]=(m[1]+(" "+n[1] if n[1] else ""),n[3] or "")
            if n[4]:out.append(n[2]+" = "+n[4]+";")
        return "\n".join(out)
    body=decl.sub(hoist,head)+sep+tail
    # OpenMP requires a local loop induction variable, not a frame field.
    # The audited loops perform scaling only; no continuation lies inside.
    loop=re.compile(r"(#pragma omp parallel for\s*\n\s*#endif\s*\n\s*)for\(i=0;i<n;i\+\+\)\s*\{")
    pos=0;parts=[];count=0
    for m in loop.finditer(body):
        if m.start()<pos:continue
        end=g.closing(body,m.end()-1,"{","}")
        block=body[m.end():end]
        parts.append(body[pos:m.start()])
        parts.append(m[1]+"for(LIS_INT full_i=0;full_i<n;full_i++){"+re.sub(r"\bi\b","full_i",block)+"}")
        pos=end+1;count+=1
    body="".join(parts)+body[pos:]
    assert count==(6 if name=="kernel" else 0),(name,count)
    lex=re.compile(r'"(?:\\.|[^"\\])*"|'+r"'(?:\\.|[^'\\])*'|"+r"->|\.|[A-Za-z_]\w*|[^\s]")
    prev="";pos=0;out=[]
    for m in lex.finditer(body):
        out.append(body[pos:m.start()])
        tok=m[0];out.append("F->"+tok if tok in fields and prev not in (".","->") else tok)
        prev=tok;pos=m.end()
    body="".join(out)+body[pos:]
    if name=="outer":
        body,n=re.subn(r"F->err\s*=\s*lis_solve_kernel\s*\(\s*F->A\s*,\s*F->effective_b\s*,\s*F->x\s*,\s*F->solver\s*,\s*F->precon\s*\);",
          "C->kernel.A=F->A; C->kernel.b=F->effective_b; C->kernel.x=F->x; C->kernel.solver=F->solver; C->kernel.precon=F->precon;\nF->pc=1; S->stage=2; return FULL_AGAIN;\nresume_1: F->pc=0; F->err=S->completion;",body)
        assert n==1
        assert not re.search(r"\blis_solve_kernel\s*\(",body)
        sites=1
    elif name=="kernel":
        body,n=re.subn(r"F->err\s*=\s*lis_solver_execute\[F->nsolver\]\(F->solver\);",
          "F->pc=1; return full_native_begin_linear(S,F->solver,F->nsolver);\nresume_1: F->pc=0; F->err=S->completion;",body)
        assert n==1
        body,n=re.subn(r"F->residual_err\s*=\s*F->shell_scale\s*\?\s*lis_matvec\(F->A,F->x,F->t\)\s*:\s*lis_matvec\(F->A,F->xx,F->t\);",
          "F->pc=2; return full_native_cert(S,F->A,F->shell_scale?F->x:F->xx,F->t);\nresume_2: F->pc=0; F->residual_err=S->completion;",body)
        assert n==1
        body,n=re.subn(r"F->t\s*=\s*NULL;\s*F->xx->precision\s*=\s*LIS_PRECISION_DEFAULT;",
          "F->t=NULL; F->xx->precision=LIS_PRECISION_DEFAULT;\nif(S->forced){F->err=S->forced_code;F->solver->retcode=F->err;F->solver->resid=NAN;goto residual_done;}",body)
        assert n==1
        assert not re.search(r"\blis_matvec\s*\(|\blis_solver_execute\s*\[",body)
        sites=2
    else:
        body,n=re.subn(r"F->err\s*=\s*lis_esolver_execute\[F->nesolver\]\(F->esolver\);",
          "F->pc=1; return full_native_begin_feast(S,F->esolver);\nresume_1: F->pc=0; F->err=S->completion;",body)
        assert n==1
        assert not re.search(r"\blis_esolver_execute\s*\[",body)
        sites=1
    finish="full_native_finish(S," if name!="kernel" else "full_linear_parent(S,"
    body=re.sub(r"\breturn\s+([^;]+);",lambda m:m[0] if m[1]=="FULL_AGAIN" or m[1].startswith(("full_native_begin_","full_native_cert(")) else "return "+finish+m[1]+");",body)
    frames[name]=(fields,body,func,sites)
def emit(names,group):
    out=["/* Generated from pinned original LIS functions, qualification only. */"]
    for name in names:
        fields,_,_,_=frames[name]
        out.append("struct full_frame_"+name+" { int pc;")
        out.extend(" "+t+" "+n+a+";" for n,(t,a) in fields.items())
        out.append("};")
    out.append("struct full_frames_"+group+" {")
    out.extend(" struct full_frame_"+name+" "+name+";" for name in names)
    out.append("};")
    if group=="linear":
        out.append("static int full_linear_parent(lis_rci_native_full *S,LIS_INT rc){S->completion=rc;S->stage=1;return FULL_AGAIN;}")
    for name in names:
        _,body,func,sites=frames[name]
        out += ["#undef __FUNC__", '#define __FUNC__ "'+func+'"',
         "static int full_step_"+name+"(lis_rci_native_full *S){",
         "struct full_frames_"+group+" *C=S->frame; struct full_frame_"+name+" *F=&C->"+name+";",
         "switch(F->pc){case 0:break;"+"".join("case "+str(i)+":goto resume_"+str(i)+";" for i in range(1,sites+1))+"default:return full_native_finish(S,LIS_BREAKDOWN);}",body,"}"]
    (g.DIR/("full_"+group+"_continuation.inc")).write_text("\n".join(out)+"\n")
emit(["outer","kernel"],"linear")
emit(["eigen"],"feast")
(g.ROOT/p).write_text(sources[p]+'\n#include "../../qualification/rci-restore/lis_rci_native_full.c"\n')
(g.ROOT/e).write_text(sources[e]+'\n#include "../../qualification/rci-restore/lis_rci_native_full_feast.c"\n')
print("FULL NATIVE LIFECYCLE GENERATED: solve=1 kernel=2 eigen=1")
