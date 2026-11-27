#!/usr/bin/env python3
"""Build a disposable library with all native cores plus versioned facade."""
import subprocess,sys
from pathlib import Path
here=Path(__file__).resolve().parent
subprocess.run([sys.executable,str(here/"qualify_all_native.py")],check=True)
p=here.parents[1]/"src/solver/lis_solver_cg.c"
p.write_text(p.read_text()+'\n#include "../../qualification/rci-restore/lis_rci_native_session.c"\n')
print("UNIFIED NATIVE SESSION PATCH READY")
