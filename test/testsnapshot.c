#include "lis.h"
#include "lis_snapshot.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if defined(_COMPLEX)
int main(void) { return 0; }
#else

static int check_matrix(const char *prefix)
{
    char name[512],header[128];
    FILE *f;
    long long n,m,nnz,i,j;
    double value;
    double total=0.0;
    long long k;
    snprintf(name,sizeof(name),"%s_A.mtx",prefix);
    f=fopen(name,"r");
    if(!f) return 1;
    if(!fgets(header,sizeof(header),f) ||
       strcmp(header,"%%MatrixMarket matrix coordinate real general\n")!=0 ||
       fscanf(f,"%lld %lld %lld",&n,&m,&nnz)!=3 ||
       n!=2 || m!=2 || nnz!=4)
    {
        fclose(f);
        return 1;
    }
    for(k=0;k<nnz;k++)
    {
        if(fscanf(f,"%lld %lld %lf",&i,&j,&value)!=3)
        {
            fclose(f);
            return 1;
        }
        total+=value;
    }
    fclose(f);
    return fabs(total-5.0)>1e-12;
}

static int check_vector(const char *prefix,const char *suffix,int cols)
{
    char name[512],header[128];
    FILE *f;
    long long n,k;
    double value;
    int i;
    snprintf(name,sizeof(name),"%s_%s.mtx",prefix,suffix);
    f=fopen(name,"r");
    if(!f) return 1;
    if(!fgets(header,sizeof(header),f) ||
       strcmp(header,"%%MatrixMarket matrix array real general\n")!=0 ||
       fscanf(f,"%lld %lld",&n,&k)!=2 ||
       n!=2 || k!=cols)
    {
        fclose(f);
        return 1;
    }
    for(i=0;i<2*cols;i++)
        if(fscanf(f,"%lf",&value)!=1 || !isfinite(value))
        {
            fclose(f);
            return 1;
        }
    fclose(f);
    return 0;
}

static void cleanup(const char *prefix)
{
    const char *names[]={"A","b","x0","x","B"};
    char name[512];
    int i;
    for(i=0;i<5;i++)
    {
        snprintf(name,sizeof(name),"%s_%s.mtx",prefix,names[i]);
        remove(name);
    }
}

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x0=NULL,x=NULL,z0=NULL,z1=NULL;
    LIS_VECTOR basis[2];
    LIS_INT err=0;
    int failed=0,keep=(argc>1);
    const char *prefix=keep ? argv[1] : "lis-snapshot-c-test";

    err=lis_initialize(&argc,&argv);
    if(err) return 1;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto fail;
    err=lis_matrix_set_size(A,0,2);
    if(err) goto fail;
    err=lis_matrix_set_value(LIS_INS_VALUE,0,0,(LIS_SCALAR)4.0,A);
    if(err) goto fail;
    err=lis_matrix_set_value(LIS_INS_VALUE,0,1,(LIS_SCALAR)-1.0,A);
    if(err) goto fail;
    err=lis_matrix_set_value(LIS_INS_VALUE,1,0,(LIS_SCALAR)-1.0,A);
    if(err) goto fail;
    err=lis_matrix_set_value(LIS_INS_VALUE,1,1,(LIS_SCALAR)3.0,A);
    if(err) goto fail;
    err=lis_matrix_assemble(A);
    if(err) goto fail;
    err=lis_vector_duplicate(A,&b);
    if(err) goto fail;
    err=lis_vector_duplicate(A,&x0);
    if(err) goto fail;
    err=lis_vector_duplicate(A,&x);
    if(err) goto fail;
    err=lis_vector_duplicate(A,&z0);
    if(err) goto fail;
    err=lis_vector_duplicate(A,&z1);
    if(err) goto fail;
    lis_vector_set_value(LIS_INS_VALUE,0,(LIS_SCALAR)2.0,b);
    lis_vector_set_value(LIS_INS_VALUE,1,(LIS_SCALAR)5.0,b);
    lis_vector_set_all((LIS_SCALAR)0.0,x0);
    lis_vector_set_value(LIS_INS_VALUE,0,(LIS_SCALAR)1.0,x);
    lis_vector_set_value(LIS_INS_VALUE,1,(LIS_SCALAR)2.0,x);
    lis_vector_set_all((LIS_SCALAR)1.0,z0);
    lis_vector_set_value(LIS_INS_VALUE,0,(LIS_SCALAR)1.0,z1);
    lis_vector_set_value(LIS_INS_VALUE,1,(LIS_SCALAR)-1.0,z1);
    basis[0]=z0;
    basis[1]=z1;
    err=lis_snapshot_export_mm(A,b,x0,x,2,basis,prefix);
    if(err || check_matrix(prefix) || check_vector(prefix,"b",1) ||
       check_vector(prefix,"x0",1) || check_vector(prefix,"x",1) ||
       check_vector(prefix,"B",2))
        goto fail;
    /* Invalid basis must fail before writing anything to the new prefix. */
    err=lis_snapshot_export_mm(A,b,x0,x,-1,basis,prefix);
    if(err!=LIS_ERR_ILL_ARG) goto fail;
    printf("LIS_GENERIC_SNAPSHOT=PASSED\n");
    goto finish;

fail:
    fprintf(stderr,"LIS_GENERIC_SNAPSHOT=FAILED err=%d\n",(int)err);
    failed=1;

finish:
    if(!keep || failed) cleanup(prefix);
    if(z1) lis_vector_destroy(z1);
    if(z0) lis_vector_destroy(z0);
    if(x) lis_vector_destroy(x);
    if(x0) lis_vector_destroy(x0);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return failed;
}
#endif
