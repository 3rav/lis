# Portable AWK assessment for real MFAMG v8 contact snapshot replay CSV.
# This only evaluates numerical evidence, NOT permission to merge/use in CCX.
# Usage: awk -f test/mfamg_replay_v8_assess.awk v8_summary.csv
BEGIN {
  FS=",";
  required="31_12m_f2_fgmres 31_12m_f2_fgcrodr_k0 31_12m_f2_fgcrodr_k8 116_12m_f2_fgmres 116_12m_f2_fgcrodr_k8";
  nreq=split(required,need," ");
}
NR==1 {
  for(i=1;i<=NF;i++) column[$i]=i;
  if(!("case_id" in column) || !("status" in column) || !("true_relative" in column) || !("intra_cycles" in column) || !("relation_error" in column)) bad=1;
  next;
}
NR>1 {
  id=$1;
  if(seen[id]++) bad=1;
  for(key in column) value[id,key]=$(column[key]);
}
function field(id,key) {return value[id,key];}
function isnum(v) {return v ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$/;}
function sane(id) {
  return isnum(field(id,"status")) && isnum(field(id,"solve_rc")) && isnum(field(id,"iter")) &&
         isnum(field(id,"true_relative")) && field(id,"true_relative")+0>=0 &&
         isnum(field(id,"reported_residual")) && isnum(field(id,"process_exit"));
}
function same_setup(a,b) {
  return field(a,"n")==field(b,"n") && field(a,"modes")==field(b,"modes") &&
         field(a,"factor")==field(b,"factor") && field(a,"restart")==field(b,"restart") &&
         field(a,"maxiter")==field(b,"maxiter");
}
function validated(id) {
  return field(id,"validation_rc")=="0" && isnum(field(id,"orth_error")) &&
         isnum(field(id,"relation_error")) && field(id,"orth_error")+0<=1e-8 &&
         field(id,"relation_error")+0<=1e-8 && isnum(field(id,"intra_cycles")) &&
         field(id,"intra_cycles")+0>=1;
}
function result(name,val) {print "ASSESS_" name "=" val;}
END {
  for(i=1;i<=nreq;i++) if(!seen[need[i]] || !sane(need[i])) bad=1;
  if(bad) {
    result("NUMERIC_GATE","INCOMPLETE_OR_INVALID_DATA");
    result("PRODUCTION_GATE","BLOCKED");
    exit 2;
  }
  b="31_12m_f2_fgmres";
  zero="31_12m_f2_fgcrodr_k0";
  rec="31_12m_f2_fgcrodr_k8";
  control="116_12m_f2_fgmres";
  controlrec="116_12m_f2_fgcrodr_k8";
  ref=1.832690841870e-3;
  brel=field(b,"true_relative")+0;
  base_ok=(field(b,"n")=="60265" && field(b,"status")=="4" &&
           field(b,"iter")+0>=999 && field(b,"iter")+0<=1001 &&
           brel/ref>=0.95 && brel/ref<=1.05);
  parity=(same_setup(b,zero) && field(b,"status")==field(zero,"status") &&
          field(b,"iter")==field(zero,"iter") &&
          field(b,"true_relative")==field(zero,"true_relative") &&
          field(b,"reported_residual")==field(zero,"reported_residual") &&
          field(b,"solve_rc")==field(zero,"solve_rc") &&
          field(b,"process_exit")==field(zero,"process_exit"));
  cbase_ok=(field(control,"status")=="0" && field(control,"process_exit")=="0" &&
            field(control,"solve_rc")=="0" && field(control,"true_relative")+0<=1e-5 &&
            field(control,"iter")+0>=530 && field(control,"iter")+0<=655);
  crec_ok=(same_setup(control,controlrec) && field(controlrec,"status")=="0" &&
           field(controlrec,"process_exit")=="0" && field(controlrec,"solve_rc")=="0" &&
           field(controlrec,"true_relative")+0<=1e-5 &&
           field(controlrec,"iter")+0<=1.05*(field(control,"iter")+0) &&
           validated(controlrec));
  rrel=field(rec,"true_relative")+0;
  recycle_numerical=((field(rec,"status")=="0" && rrel<=1e-5) ||
                     (field(rec,"status")=="4" && field(rec,"iter")==field(b,"iter") &&
                      rrel<=0.90*brel));
  recycle_ok=(same_setup(b,rec) && recycle_numerical && validated(rec));
  result("31_BASELINE",base_ok?"REPRODUCED":"REVIEW_MISMATCH");
  result("31_K0_PARITY",parity?"PASS":"FAIL");
  result("31_K8_NUMERICS",recycle_ok?"IMPROVED":"NOT_DEMONSTRATED");
  result("116_BASELINE",cbase_ok?"REPRODUCED":"REVIEW_MISMATCH");
  result("116_K8_CONTROL",crec_ok?"PASS":"FAIL");
  if(brel>0) result("31_K8_RESIDUAL_RATIO",sprintf("%.6g",rrel/brel));
  if(field(control,"iter")+0>0)
    result("116_K8_ITER_RATIO",sprintf("%.6g",field(controlrec,"iter")/(field(control,"iter")+0)));
  allok=base_ok && parity && recycle_ok && cbase_ok && crec_ok;
  result("NUMERIC_GATE",allok?"PASS_CONTACT_REPLAY_ONLY":"FAIL_OR_INCONCLUSIVE");
  result("PRODUCTION_GATE","BLOCKED_NEEDS_REPEATED_TIMING_AND_NEWTON_SEQUENCE");
}
