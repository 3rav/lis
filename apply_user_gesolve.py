#!/usr/bin/env python3
# Apply stage-1 LIS_MATRIX_USER generalized-eigensolver support.
#
# Target: current LIS 2.1.13 master containing the official LIS_MATRIX_USER API.
#
# Scope:
#   - make lis_gesolve() handle LIS_MATRIX_USER safely
#   - reject explicit storage conversion for user matrices
#   - reject matrix shifts for user matrices until a generic composed/shifted
#     operator exists
#   - print "user/shell" instead of indexing the explicit-storage name table
#   - add a generalized eigenproblem regression test for GPI and GII using
#     LIS_MATRIX_USER for both A and B

from pathlib import Path
import sys

ROOT = Path.cwd()

required = [
    "include/lis.h",
    "src/esolver/lis_esolver.c",
    "test/Makefile.am",
]
for rel in required:
    if not (ROOT / rel).is_file():
        sys.exit(f"ERROR: {rel} not found. Run this script from the LIS repository root.")

lis_h = (ROOT / "include/lis.h").read_text(encoding="utf-8")
esolver_path = ROOT / "src/esolver/lis_esolver.c"
makefile_path = ROOT / "test/Makefile.am"
test_path = ROOT / "test/getestuser.c"

if '#define LIS_VERSION\t"2.1.13"' not in lis_h and '#define LIS_VERSION "2.1.13"' not in lis_h:
    sys.exit("ERROR: expected LIS 2.1.13 base was not found in include/lis.h")
if "#define LIS_MATRIX_USER 22" not in lis_h:
    sys.exit("ERROR: official LIS_MATRIX_USER definition was not found")
if "lis_matrix_set_user" not in lis_h:
    sys.exit("ERROR: official lis_matrix_set_user() API was not found")

esolver = esolver_path.read_text(encoding="utf-8")
makefile = makefile_path.read_text(encoding="utf-8")

guard_marker = "matrix shifts are unavailable for LIS_MATRIX_USER"
print_marker = 'matrix storage format : user/shell'
make_marker = "getestuser_SOURCES = getestuser.c"

already = (
    guard_marker in esolver
    and print_marker in esolver
    and "getestuser" in makefile
    and test_path.is_file()
)
if already:
    print("Stage-1 user-matrix generalized-eigensolver changes are already present.")
    sys.exit(0)

# ---------------------------------------------------------------------------
# 1. lis_gesolve(): reject operations that require materialized storage.
# ---------------------------------------------------------------------------
option_anchor = '''\trval = esolver->options[LIS_EOPTIONS_RVAL];
\tif( nesolver < 1 || nesolver > LIS_ESOLVER_LEN )
'''

option_replacement = '''\trval = esolver->options[LIS_EOPTIONS_RVAL];

\t/*
\t * A user-defined matrix is an operator, not materialized LIS storage.
\t * Keep the first generalized-eigensolver step operator-only.  Explicit
\t * storage conversion and shifted/composed operators are separate features.
\t */
\tif( (A->matrix_type==LIS_MATRIX_USER ||
\t     (B!=NULL && B->matrix_type==LIS_MATRIX_USER)) && estorage>0 )
\t{
\t\tLIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
\t\t           "storage conversion is unavailable for LIS_MATRIX_USER\\n");
\t\treturn LIS_ERR_NOT_IMPLEMENTED;
\t}
\tif( (A->matrix_type==LIS_MATRIX_USER ||
\t     (B!=NULL && B->matrix_type==LIS_MATRIX_USER)) &&
\t    (oshift!=0.0 || ishift!=0.0) )
\t{
\t\tLIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
\t\t           "matrix shifts are unavailable for LIS_MATRIX_USER\\n");
\t\treturn LIS_ERR_NOT_IMPLEMENTED;
\t}

\tif( nesolver < 1 || nesolver > LIS_ESOLVER_LEN )
'''

if guard_marker not in esolver:
    if option_anchor not in esolver:
        sys.exit("ERROR: lis_gesolve option anchor not found; upstream code changed")
    esolver = esolver.replace(option_anchor, option_replacement, 1)

# ---------------------------------------------------------------------------
# 2. lis_gesolve(): avoid indexing lis_estoragename[] with LIS_MATRIX_USER=22.
# ---------------------------------------------------------------------------
print_anchor = '''\tif( A->matrix_type==LIS_MATRIX_BSR || A->matrix_type==LIS_MATRIX_BSC )
\t  {
\t    if ( output ) lis_printf(comm,"matrix storage format : %s(%D x %D)\\n", lis_estoragename[A->matrix_type-1],eblock,eblock);
\t  }
\telse
\t  {
\t    if ( output ) lis_printf(comm,"matrix storage format : %s\\n", lis_estoragename[A->matrix_type-1]);
\t  }
'''

print_replacement = '''\tif( A->matrix_type==LIS_MATRIX_BSR || A->matrix_type==LIS_MATRIX_BSC )
\t  {
\t    if ( output ) lis_printf(comm,"matrix storage format : %s(%D x %D)\\n", lis_estoragename[A->matrix_type-1],eblock,eblock);
\t  }
\telse if( A->matrix_type==LIS_MATRIX_USER )
\t  {
\t    if ( output ) lis_printf(comm,"matrix storage format : user/shell\\n");
\t  }
\telse
\t  {
\t    if ( output ) lis_printf(comm,"matrix storage format : %s\\n", lis_estoragename[A->matrix_type-1]);
\t  }
'''

if print_marker not in esolver:
    if print_anchor not in esolver:
        sys.exit("ERROR: lis_gesolve storage-print anchor not found; upstream code changed")
    esolver = esolver.replace(print_anchor, print_replacement, 1)

# ---------------------------------------------------------------------------
# 3. Add getestuser to the automake test programs.
# ---------------------------------------------------------------------------
program_anchor = "test_PROGRAMS = etest1 etest2 etest3 etest4 etest5 etest5b etest6 etest7 getest1 getest5 getest5b spmvtest1 spmvtest2 spmvtest2b spmvtest3 spmvtest3b spmvtest4 spmvtest5 test1 test2 test2b test3 test3b test3c test4 test5 test6 test7 testuser"
program_replacement = program_anchor + " getestuser"

if "getestuser" not in makefile:
    if program_anchor not in makefile:
        sys.exit("ERROR: expected test_PROGRAMS line not found; upstream code changed")
    makefile = makefile.replace(program_anchor, program_replacement, 1)

source_anchor = "testuser_SOURCES = testuser.c"
source_replacement = source_anchor + "\ngetestuser_SOURCES = getestuser.c"

if make_marker not in makefile:
    if source_anchor not in makefile:
        sys.exit("ERROR: testuser_SOURCES anchor not found; upstream code changed")
    makefile = makefile.replace(source_anchor, source_replacement, 1)

# ---------------------------------------------------------------------------
# 4. Generalized eigenproblem regression test.
#    A = diag(2,6,12,20), B = diag(1,2,3,4)
#    generalized eigenvalues = 2,3,4,5
#    GPI should converge to 5; GII should converge to 2.
# ---------------------------------------------------------------------------
test_source = r'''#include "lis.h"
#include <stdio.h>
#include <math.h>

typedef struct
{
	LIS_INT n;
	LIS_INT is;
	LIS_INT which;
} USER_DIAG;

static LIS_SCALAR user_diag_value(USER_DIAG *ctx, LIS_INT global_i)
{
	static const LIS_SCALAR adiag[4] = {2.0, 6.0, 12.0, 20.0};
	static const LIS_SCALAR bdiag[4] = {1.0, 2.0, 3.0, 4.0};

	return ctx->which==0 ? adiag[global_i] : bdiag[global_i];
}

static LIS_INT user_diag_matvec(void *user_data,
                                const LIS_SCALAR *x,
                                LIS_SCALAR *y)
{
	USER_DIAG *ctx = (USER_DIAG *)user_data;
	LIS_INT i;

	for(i=0;i<ctx->n;i++)
	{
		y[i] = user_diag_value(ctx,ctx->is+i) * x[i];
	}
	return LIS_SUCCESS;
}

/* The test operators are real diagonal matrices, so A^H == A and B^H == B. */
static LIS_INT user_diag_matvech(void *user_data,
                                 const LIS_SCALAR *x,
                                 LIS_SCALAR *y)
{
	return user_diag_matvec(user_data,x,y);
}

static LIS_INT create_user_diag(LIS_Comm comm,
                                LIS_INT global_n,
                                LIS_INT which,
                                USER_DIAG *ctx,
                                LIS_MATRIX *A)
{
	LIS_INT err, local_n, is, ie;

	err = lis_matrix_create(comm,A);
	if( err ) return err;

	err = lis_matrix_set_size(*A,0,global_n);
	if( err )
	{
		lis_matrix_destroy(*A);
		*A = NULL;
		return err;
	}

	err = lis_matrix_get_size(*A,&local_n,&global_n);
	if( err )
	{
		lis_matrix_destroy(*A);
		*A = NULL;
		return err;
	}
	err = lis_matrix_get_range(*A,&is,&ie);
	if( err )
	{
		lis_matrix_destroy(*A);
		*A = NULL;
		return err;
	}

	ctx->n = local_n;
	ctx->is = is;
	ctx->which = which;
	(void)ie;

	err = lis_matrix_set_user(*A,ctx,user_diag_matvec,user_diag_matvech);
	if( err )
	{
		lis_matrix_destroy(*A);
		*A = NULL;
		return err;
	}

	return LIS_SUCCESS;
}

static LIS_INT run_case(LIS_Comm comm,
                        LIS_MATRIX A,
                        LIS_MATRIX B,
                        LIS_VECTOR x,
                        const char *method,
                        LIS_SCALAR expected)
{
	LIS_ESOLVER esolver;
	LIS_SCALAR evalue;
	LIS_REAL error;
	LIS_INT err;
	char option[128];

	err = lis_esolver_create(&esolver);
	if( err ) return err;

	snprintf(option,sizeof(option),
	         "-e %s -etol 1.0e-10 -emaxiter 1000 -eprint mem",
	         method);

	err = lis_esolver_set_option(option,esolver);
	if( err )
	{
		lis_esolver_destroy(esolver);
		return err;
	}

	err = lis_esolver_set_optionC(esolver);
	if( err )
	{
		lis_esolver_destroy(esolver);
		return err;
	}

	err = lis_gesolve(A,B,x,&evalue,esolver);
	if( err )
	{
		fprintf(stderr,"%s: lis_gesolve failed: %d\n",method,(int)err);
		lis_esolver_destroy(esolver);
		return err;
	}

	error = (LIS_REAL)fabs(evalue-expected);
	lis_printf(comm,
	           "LIS_MATRIX_USER generalized %s eigenvalue = %.15e, error = %.3e\n",
	           method,(double)evalue,(double)error);

	lis_esolver_destroy(esolver);

	if( error > 1.0e-7 )
	{
		fprintf(stderr,"%s: eigenvalue verification failed\n",method);
		return LIS_FAILS;
	}

	return LIS_SUCCESS;
}

int main(int argc, char *argv[])
{
	LIS_Comm comm;
	LIS_MATRIX A = NULL, B = NULL;
	LIS_VECTOR x = NULL;
	USER_DIAG actx, bctx;
	LIS_INT err;
	const LIS_INT global_n = 4;

	err = lis_initialize(&argc,&argv);
	if( err ) return 1;

	comm = LIS_COMM_WORLD;

	err = create_user_diag(comm,global_n,0,&actx,&A);
	if( err ) goto fail;

	err = create_user_diag(comm,global_n,1,&bctx,&B);
	if( err ) goto fail;

	err = lis_vector_duplicate(A,&x);
	if( err ) goto fail;

	/* GPI converges to the largest generalized eigenvalue: 5. */
	err = run_case(comm,A,B,x,"gpi",5.0);
	if( err ) goto fail;

	/* GII converges to the smallest generalized eigenvalue: 2. */
	err = run_case(comm,A,B,x,"gii",2.0);
	if( err ) goto fail;

	lis_printf(comm,"LIS_MATRIX_USER generalized eigensolver test PASSED\n");

	lis_vector_destroy(x);
	lis_matrix_destroy(A);
	lis_matrix_destroy(B);
	lis_finalize();
	return 0;

fail:
	fprintf(stderr,"getestuser failed: %d\n",(int)err);
	if( x ) lis_vector_destroy(x);
	if( A ) lis_matrix_destroy(A);
	if( B ) lis_matrix_destroy(B);
	lis_finalize();
	return 1;
}
'''

if test_path.exists():
    existing = test_path.read_text(encoding="utf-8")
    if existing != test_source:
        sys.exit("ERROR: test/getestuser.c already exists with different content")
else:
    test_path.write_text(test_source, encoding="utf-8")

# ---------------------------------------------------------------------------
# Write modified tracked files only after every anchor has been validated.
# ---------------------------------------------------------------------------
esolver_path.write_text(esolver, encoding="utf-8")
makefile_path.write_text(makefile, encoding="utf-8")

# Sanity checks.
final_esolver = esolver_path.read_text(encoding="utf-8")
final_makefile = makefile_path.read_text(encoding="utf-8")

checks = [
    (guard_marker in final_esolver, "user-matrix shift guard"),
    ("storage conversion is unavailable for LIS_MATRIX_USER" in final_esolver,
     "user-matrix storage guard"),
    (print_marker in final_esolver, "user/shell storage output"),
    ("getestuser" in final_makefile, "getestuser Makefile entry"),
    (test_path.is_file(), "getestuser source"),
]
missing = [name for ok, name in checks if not ok]
if missing:
    sys.exit("ERROR: sanity check failed: " + ", ".join(missing))

print("Stage-1 LIS_MATRIX_USER generalized-eigensolver patch applied.")
print("Changed:")
print("  src/esolver/lis_esolver.c")
print("  test/Makefile.am")
print("  test/getestuser.c")
