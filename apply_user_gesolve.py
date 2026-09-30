#!/usr/bin/env python3
# Stage 1: enable safe LIS_MATRIX_USER use in generalized eigensolver path.
# Target base: current LIS 2.1.13 master with Akira's official LIS_MATRIX_USER.

from pathlib import Path
import sys

ROOT = Path.cwd()

for rel in ("include/lis.h", "src/esolver/lis_esolver.c", "test/Makefile.am"):
    if not (ROOT / rel).is_file():
        sys.exit(f"ERROR: {rel} not found; run from LIS repository root")

lis_h = (ROOT / "include/lis.h").read_bytes().decode("utf-8")
esolver_path = ROOT / "src/esolver/lis_esolver.c"
makefile_path = ROOT / "test/Makefile.am"
test_path = ROOT / "test/getestuser.c"

if '#define LIS_VERSION\t"2.1.13"' not in lis_h:
    sys.exit("ERROR: expected LIS 2.1.13 base")
if "#define LIS_MATRIX_USER 22" not in lis_h:
    sys.exit("ERROR: LIS_MATRIX_USER=22 not found")
if "lis_matrix_set_user" not in lis_h:
    sys.exit("ERROR: lis_matrix_set_user() API not found")

esolver = esolver_path.read_bytes().decode("utf-8")
makefile = makefile_path.read_bytes().decode("utf-8")
esolver_eol = "\r\n" if "\r\n" in esolver else "\n"
makefile_eol = "\r\n" if "\r\n" in makefile else "\n"

guard_marker = "matrix shifts are unavailable for LIS_MATRIX_USER"
print_marker = "matrix storage format : user/shell"

# ----------------------------------------------------------------------
# 1. Insert user-matrix guards after the single RVAL assignment.
# ----------------------------------------------------------------------
if guard_marker not in esolver:
    anchor = "\trval = esolver->options[LIS_EOPTIONS_RVAL];"
    if esolver.count(anchor) != 1:
        sys.exit(
            "ERROR: expected exactly one "
            "rval = esolver->options[LIS_EOPTIONS_RVAL] assignment"
        )

    insertion = r'''
	/*
	 * LIS_MATRIX_USER is an operator, not materialized LIS storage.
	 * Explicit storage conversion and matrix shifts are separate features.
	 */
	if( A->matrix_type==LIS_MATRIX_USER && estorage>0 )
	{
		LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
		           "storage conversion is unavailable for LIS_MATRIX_USER\n");
		return LIS_ERR_NOT_IMPLEMENTED;
	}
	if( (A->matrix_type==LIS_MATRIX_USER ||
	     (B!=NULL && B->matrix_type==LIS_MATRIX_USER)) &&
	    (oshift!=0.0 || ishift!=0.0) )
	{
		LIS_SETERR(LIS_ERR_NOT_IMPLEMENTED,
		           "matrix shifts are unavailable for LIS_MATRIX_USER\n");
		return LIS_ERR_NOT_IMPLEMENTED;
	}
'''
    insertion = insertion.rstrip("\n").replace("\n", esolver_eol)
    esolver = esolver.replace(anchor, anchor + esolver_eol + insertion, 1)

# ----------------------------------------------------------------------
# 2. Avoid lis_estoragename[A->matrix_type-1] for USER=22.
#    Replace only the generic one-line output inside lis_gesolve().
# ----------------------------------------------------------------------
if print_marker not in esolver:
    old = '\t    if ( output ) lis_printf(comm,"matrix storage format : %s\\n", lis_estoragename[A->matrix_type-1]);'
    new = (
        '\t    if( A->matrix_type==LIS_MATRIX_USER )\n'
        '\t      {\n'
        '\t        if ( output ) lis_printf(comm,"matrix storage format : user/shell\\n");\n'
        '\t      }\n'
        '\t    else\n'
        '\t      {\n'
        '\t        if ( output ) lis_printf(comm,"matrix storage format : %s\\n", lis_estoragename[A->matrix_type-1]);\n'
        '\t      }'
    )
    if esolver.count(old) != 1:
        sys.exit(
            "ERROR: expected exactly one generic matrix-storage output line "
            "in lis_gesolve()"
        )
    new = new.replace("\n", esolver_eol)
    esolver = esolver.replace(old, new, 1)

# ----------------------------------------------------------------------
# 3. Register new test in test/Makefile.am.
# ----------------------------------------------------------------------
program_line = (
    "test_PROGRAMS = etest1 etest2 etest3 etest4 etest5 etest5b "
    "etest6 etest7 getest1 getest5 getest5b spmvtest1 spmvtest2 "
    "spmvtest2b spmvtest3 spmvtest3b spmvtest4 spmvtest5 test1 "
    "test2 test2b test3 test3b test3c test4 test5 test6 test7 testuser"
)

if "getestuser_SOURCES = getestuser.c" not in makefile:
    if makefile.count(program_line) != 1:
        sys.exit("ERROR: expected test_PROGRAMS line not found uniquely")
    makefile = makefile.replace(program_line, program_line + " getestuser", 1)

    source_line = "testuser_SOURCES = testuser.c"
    if makefile.count(source_line) != 1:
        sys.exit("ERROR: testuser_SOURCES line not found uniquely")
    makefile = makefile.replace(
        source_line,
        source_line + makefile_eol + "getestuser_SOURCES = getestuser.c",
        1,
    )

# ----------------------------------------------------------------------
# 4. Add generalized user-matrix regression test.
#    A = diag(2,6,12,20), B = diag(1,2,3,4)
#    eigenvalues = 2,3,4,5
#    GPI -> largest 5, GII -> smallest 2.
# ----------------------------------------------------------------------
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
	LIS_INT err, local_n, is, ie, gn;

	err = lis_matrix_create(comm,A);
	if( err ) return err;

	err = lis_matrix_set_size(*A,0,global_n);
	if( err ) goto fail;

	err = lis_matrix_get_size(*A,&local_n,&gn);
	if( err ) goto fail;

	err = lis_matrix_get_range(*A,&is,&ie);
	if( err ) goto fail;

	ctx->n = local_n;
	ctx->is = is;
	ctx->which = which;
	(void)ie;
	(void)gn;

	err = lis_matrix_set_user(*A,ctx,user_diag_matvec,user_diag_matvech);
	if( err ) goto fail;

	return LIS_SUCCESS;

fail:
	lis_matrix_destroy(*A);
	*A = NULL;
	return err;
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
	if( err ) goto out;

	err = lis_gesolve(A,B,x,&evalue,esolver);
	if( err )
	{
		fprintf(stderr,"%s: lis_gesolve failed: %d\n",method,(int)err);
		goto out;
	}

	error = (LIS_REAL)fabs(evalue-expected);
	lis_printf(comm,
	           "LIS_MATRIX_USER generalized %s eigenvalue = %.15e, error = %.3e\n",
	           method,(double)evalue,(double)error);

	if( error > 1.0e-7 )
	{
		fprintf(stderr,"%s: eigenvalue verification failed\n",method);
		err = LIS_FAILS;
	}

out:
	lis_esolver_destroy(esolver);
	return err;
}

LIS_INT main(int argc, char *argv[])
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

	err = run_case(comm,A,B,x,"gpi",5.0);
	if( err ) goto fail;

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
    if test_path.read_bytes().decode("utf-8") != test_source:
        sys.exit("ERROR: test/getestuser.c already exists with different content")
else:
    test_path.write_bytes(test_source.encode("utf-8"))

# Make sure the snippets introduced by this generator contain no trailing
# whitespace. Do not normalize the whole upstream file.
for marker in (
    "storage conversion is unavailable for LIS_MATRIX_USER",
    "matrix shifts are unavailable for LIS_MATRIX_USER",
    "matrix storage format : user/shell",
):
    if marker not in esolver:
        sys.exit("ERROR: expected generated marker is missing: " + marker)

# Write tracked files while preserving the upstream newline convention.
esolver_path.write_bytes(esolver.encode("utf-8"))
makefile_path.write_bytes(makefile.encode("utf-8"))

# Final sanity checks.
final_esolver = esolver_path.read_bytes().decode("utf-8")
final_makefile = makefile_path.read_bytes().decode("utf-8")

checks = [
    guard_marker in final_esolver,
    print_marker in final_esolver,
    "getestuser" in final_makefile,
    test_path.is_file(),
]
if not all(checks):
    sys.exit("ERROR: final sanity checks failed")

print("Stage-1 generalized LIS_MATRIX_USER changes applied.")
print("Changed:")
print("  src/esolver/lis_esolver.c")
print("  test/Makefile.am")
print("  test/getestuser.c")
