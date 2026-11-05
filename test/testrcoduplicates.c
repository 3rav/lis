/*
 * RCO regression: updating existing entries must not realloc full rows.
 * Numeric verification uses the assembled CSR matrix and SpMV.
 */
#include <stdio.h>
#include "lis.h"

#define CHECK(expr) do { \
    LIS_INT rc=(expr); \
    if(rc!=LIS_SUCCESS) { \
        fprintf(stderr,"%s failed: %ld (line %d)\n",#expr,(long)rc,__LINE__); \
        failed=1; \
        goto cleanup; \
    } \
} while(0)

#define ASSERT_TRUE(expr) do { \
    if(!(expr)) { \
        fprintf(stderr,"Assertion failed: %s (line %d)\n",#expr,__LINE__); \
        failed=1; \
        goto cleanup; \
    } \
} while(0)

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR x=NULL,y=NULL;
    LIS_INT i,err;
    LIS_INT *initial_index;
    LIS_SCALAR val;
    int failed=0;

    err=lis_initialize(&argc,&argv);
    if(err) return 1;

    CHECK(lis_matrix_create(LIS_COMM_WORLD,&A));
    CHECK(lis_matrix_set_size(A,0,3));
    CHECK(lis_matrix_malloc(A,1,NULL));

    CHECK(lis_matrix_set_value(LIS_INS_VALUE,0,0,2.0,A));
    CHECK(lis_matrix_set_value(LIS_INS_VALUE,1,1,3.0,A));
    CHECK(lis_matrix_set_value(LIS_INS_VALUE,2,2,4.0,A));
    ASSERT_TRUE(A->w_row[0]==1 && A->w_nnz[0]==1);
    initial_index=A->w_index[0];

    for(i=0;i<200;i++) {
        CHECK(lis_matrix_set_value(LIS_ADD_VALUE,0,0,0.25,A));
        CHECK(lis_matrix_set_value(LIS_INS_VALUE,1,1,3.0,A));
    }
    ASSERT_TRUE(A->w_nnz[0]==1 && A->w_row[0]==1);
    ASSERT_TRUE(A->w_nnz[1]==1 && A->w_row[1]==1);
    ASSERT_TRUE(A->w_index[0]==initial_index);

    CHECK(lis_matrix_set_value(LIS_INS_VALUE,0,1,-1.0,A));
    ASSERT_TRUE(A->w_row[0]==2 && A->w_nnz[0]==2);
    for(i=0;i<100;i++) {
        CHECK(lis_matrix_set_value(LIS_ADD_VALUE,0,1,0.5,A));
    }
    ASSERT_TRUE(A->w_row[0]==2 && A->w_nnz[0]==2);

    CHECK(lis_matrix_assemble(A));
    CHECK(lis_vector_duplicate(A,&x));
    CHECK(lis_vector_duplicate(A,&y));
    CHECK(lis_vector_set_value(LIS_INS_VALUE,0,1.0,x));
    CHECK(lis_vector_set_value(LIS_INS_VALUE,1,2.0,x));
    CHECK(lis_vector_set_value(LIS_INS_VALUE,2,3.0,x));
    CHECK(lis_matvec(A,x,y));

    CHECK(lis_vector_get_value(y,0,&val));
    ASSERT_TRUE(val==(LIS_SCALAR)150.0); /* 52*1 + 49*2 */
    CHECK(lis_vector_get_value(y,1,&val));
    ASSERT_TRUE(val==(LIS_SCALAR)6.0);
    CHECK(lis_vector_get_value(y,2,&val));
    ASSERT_TRUE(val==(LIS_SCALAR)12.0);

    puts("RCO_DUPLICATE_UPDATE_PASSED");

cleanup:
    if(y) lis_vector_destroy(y);
    if(x) lis_vector_destroy(x);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return failed;
}
