#include "lis.h"
#include "lis_snapshot.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if defined(_COMPLEX)
int main(void) { return 0; }
#else

static int header_matches(char *header,const char *expected)
{
    /* Allow both Windows CRLF and Unix LF when reading Matrix Market text. */
    header[strcspn(header,"\r\n")]='\0';
    return strcmp(header,expected)==0;
}

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
    header[0]='\0';
    if(!fgets(header,sizeof(header),f) ||
       !header_matches(header,"%%MatrixMarket matrix coordinate real general") ||
       fscanf(f,"%lld %lld %lld",&n,&m,&nnz)!=3 ||
       n!=2 || m!=2 || nnz!=4)
    {
        fprintf(stderr,"snapshot matrix header/dimension mismatch, header='%s'\n",header);
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
    header[0]='\0';
    if(!fgets(header,sizeof(header),f) ||
       !header_matches(header,"%%MatrixMarket matrix array real general") ||
       fscanf(f,"%lld %lld",&n,&k)!=2 ||
       n!=2 || k!=cols)
    {
        fprintf(stderr,"snapshot vector '%s' header/dimension mismatch, header='%s'\n",suffix,header);
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

/* A non-CalculiX USER matrix with a generic row callback. */
static LIS_INT user_matvec(void *userdata,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    (void)userdata;
    y[0]=(LIS_SCALAR)4.0*x[0]-x[1];
    y[1]=-x[0]+(LIS_SCALAR)3.0*x[1];
    return LIS_SUCCESS;
}
static LIS_INT user_row(void *userdata,LIS_INT row,LIS_INT *nnz,
                        LIS_INT *cols,LIS_SCALAR *vals)
{
    (void)userdata;
    if(row<0 || row>1 || nnz==NULL) return LIS_ERR_ILL_ARG;
    if(cols==NULL && vals==NULL)
    {
        *nnz=2;
        return LIS_SUCCESS;
    }
    if(cols==NULL || vals==NULL || *nnz<2) return LIS_ERR_ILL_ARG;
    cols[0]=0;
    cols[1]=1;
    vals[0]=(row==0)?(LIS_SCALAR)4.0:(LIS_SCALAR)-1.0;
    vals[1]=(row==0)?(LIS_SCALAR)-1.0:(LIS_SCALAR)3.0;
    *nnz=2;
    return LIS_SUCCESS;
}
static int check_user_snapshot(LIS_VECTOR b,LIS_VECTOR x0,LIS_VECTOR x,
                               LIS_VECTOR basis[],const char *prefix)
{
    LIS_MATRIX A=NULL;
    char name[512];
    int failed=1;
    LIS_INT err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto finish;
    err=lis_matrix_set_size(A,0,2);
    if(err) goto finish;
    err=lis_matrix_set_user(A,NULL,user_matvec,user_matvec);
    if(err) goto finish;
    err=lis_matrix_set_user_row(A,user_row);
    if(err) goto finish;
    err=lis_snapshot_export_mm(A,b,x0,x,2,basis,prefix);
    if(err) fprintf(stderr,"USER export failed rc=%d\n",(int)err);
    if(err || check_matrix(prefix) || check_vector(prefix,"B",2))
        goto finish;
    failed=0;
finish:
    /* USER and CSR must have identical standard Matrix Market content. */
    snprintf(name,sizeof(name),"%s_A.mtx",prefix);
    remove(name);
    snprintf(name,sizeof(name),"%s_b.mtx",prefix);
    remove(name);
    snprintf(name,sizeof(name),"%s_x0.mtx",prefix);
    remove(name);
    snprintf(name,sizeof(name),"%s_x.mtx",prefix);
    remove(name);
    snprintf(name,sizeof(name),"%s_B.mtx",prefix);
    remove(name);
    if(A) lis_matrix_destroy(A);
    return failed;
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
    if(err) fprintf(stderr,"CSR export failed rc=%d\n",(int)err);
    if(err || check_matrix(prefix) || check_vector(prefix,"b",1) ||
       check_vector(prefix,"x0",1) || check_vector(prefix,"x",1) ||
       check_vector(prefix,"B",2))
        goto fail;
    /* Invalid basis must fail before writing anything to the new prefix. */
    err=lis_snapshot_export_mm(A,b,x0,x,-1,basis,prefix);
    if(err!=LIS_ERR_ILL_ARG) goto fail;
    {
        char user_prefix[512];
        snprintf(user_prefix,sizeof(user_prefix),"%s-user",prefix);
        if(check_user_snapshot(b,x0,x,basis,user_prefix))
            goto fail;
    }
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
