# MFAMG v8: portable POSIX-awk log-to-CSV summary.
# Usage: awk -f test/mfamg_replay_v8_report.awk v8_case1.log ...
# Input may include status=4/MAXITER. A valid REPLAY_RESULT is required.
# Numeric log fields are emitted as recorded; NEVER infer convergence from rc.
BEGIN {
  OFS=","
  print "case_id,snapshot,n,modes,factor,method,recycle_k,restart,maxiter,process_exit,solve_rc,status,iter,reported_residual,true_relative,wall_seconds,recycle_dim,builds,updates,uses,rejects,intra_cycles,validation_rc,orth_error,relation_error,log_file"
}
function getval(key,line,    parts,n,i,eq) {
  n=split(line,parts,/[[:space:]]+/)
  for(i=1;i<=n;i++) {
    eq=index(parts[i],"=")
    if(eq && substr(parts[i],1,eq-1)==key)
      return substr(parts[i],eq+1)
  }
  return "NA"
}
function base(path,    s) {
  s=path
  sub(/^.*\//,"",s)
  return s
}
function reset(path) {
  current=path
  filebase=base(path)
  case_id=filebase
  snapshot=n=modes=factor=method=recycle_k=restart=maxiter="NA"
  process_exit=solve_rc=status=iter=reported_residual=true_relative=wall_seconds="NA"
  recycle_dim=builds=updates=uses=rejects="NA"
  validation_rc=orth_error=relation_error="NA"
  intra_cycles=0
  has_result=0
}
function emit() {
  if(current=="") return
  if(!has_result) {
    print "ERROR: missing REPLAY_RESULT in " current > "/dev/stderr"
    invalid=1
  }
  print case_id,snapshot,n,modes,factor,method,recycle_k,restart,maxiter, \
        process_exit,solve_rc,status,iter,reported_residual,true_relative, \
        wall_seconds,recycle_dim,builds,updates,uses,rejects,intra_cycles, \
        validation_rc,orth_error,relation_error,filebase
}
FNR==1 {
  if(NR>1) emit()
  reset(FILENAME)
}
/^REPLAY_SWEEP_META / {
  case_id=getval("case",$0)
  snapshot=getval("snapshot",$0)
  factor=getval("factor",$0)
}
/^REPLAY n=/ {
  n=getval("n",$0)
  modes=getval("modes",$0)
  restart=getval("restart",$0)
  maxiter=getval("maxiter",$0)
  method=getval("method",$0)
  recycle_k=getval("recycle_k",$0)
}
/^REPLAY_RESULT / {
  if(has_result) {
    print "ERROR: repeated REPLAY_RESULT in " FILENAME > "/dev/stderr"
    invalid=1
  }
  has_result=1
  solve_rc=getval("solve_rc",$0)
  status=getval("status",$0)
  iter=getval("iter",$0)
  reported_residual=getval("reported_residual",$0)
  true_relative=getval("true_relative",$0)
  wall_seconds=getval("wall_seconds",$0)
}
/^REPLAY_RECYCLE dim=/ { recycle_dim=getval("dim",$0) }
/^REPLAY_RECYCLE_STATS / {
  builds=getval("builds",$0)
  updates=getval("updates",$0)
  uses=getval("uses",$0)
  rejects=getval("rejects",$0)
}
/^REPLAY_RECYCLE_CHECK / {
  validation_rc=getval("rc",$0)
  orth_error=getval("orth_error",$0)
  relation_error=getval("relation_error",$0)
}
/^LIS_FGCRODR_INTRA_RECYCLE iter=/ { intra_cycles++ }
/^REPLAY_SWEEP_EXIT / { process_exit=getval("process_exit",$0) }
END {
  emit()
  if(invalid) exit 2
}
