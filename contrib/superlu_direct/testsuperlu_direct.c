#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lislib.h"
#include "lis_superlu_direct.h"

#define N 6
static const double M[N][N]={
    {5,1,0,0,0,0.7},
    {0.4,4,2,0,0,0},
    {0,0.7,6,-0.3,0,0},
    {0,0,0.5,3,0.3,0},
    {0,0,0,1.3,9,1.1},
    {0.6,0,0,0,0.4,2}
};
int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_VECTOR b=NULL,bt=NULL,x=NULL;
    LIS_INT e,i,j,t;
    const double scales[N]={1e-4,1e3,1e-2,1e2,1e-3,1e4};
    const double v[N]={0.6,-0.25,1.2,2.0,-1.0,0.3};
    double maxerr=0.0,max_relative_residual=0.0;

    e=lis_initialize(&argc,&argv);if(e)return 1;
    e=lis_superlu_direct_register();if(e)goto fail;
    e=lis_matrix_create(LIS_COMM_WORLD,&A);if(e)goto fail;
    e=lis_matrix_set_size(A,0,N);if(e)goto fail;
    for(i=0;i<N;i++)for(j=0;j<N;j++)if(M[i][j]!=0)
    {
        e=lis_matrix_set_value(LIS_INS_VALUE,i,j,M[i][j]*scales[i],A);
        if(e)goto fail;
    }
    e=lis_matrix_set_type(A,LIS_MATRIX_CSR);if(e)goto fail;
    e=lis_matrix_assemble(A);if(e)goto fail;
    e=lis_vector_duplicate(A,&b);if(e)goto fail;
    e=lis_vector_duplicate(A,&bt);if(e)goto fail;
    e=lis_vector_duplicate(A,&x);if(e)goto fail;
    for(i=0;i<N;i++)
    {
        b->value[i]=bt->value[i]=0;
        for(j=0;j<N;j++)
        {
            b->value[i]+=(LIS_SCALAR)(M[i][j]*scales[i]*v[j]);
            bt->value[i]+=(LIS_SCALAR)(M[j][i]*scales[j]*v[j]);
        }
    }
    e=lis_solver_create(&solver);if(e)goto fail;
    solver->A=A;
    e=lis_solver_set_option("-p superlu -i fgmres -scale none",solver);
    if(e)goto fail;
    e=lis_precon_create(solver,&pc);if(e){pc=NULL;goto fail;}
    solver->precon=pc;
    for(t=0;t<2;t++)
    {
        e=lis_psolve(solver,b,x);if(e)goto fail;
        for(i=0;i<N;i++)
        {
            double r=-b->value[i];
            for(j=0;j<N;j++)r+=M[i][j]*scales[i]*x->value[j];
            if(fabs(r)/(fabs(b->value[i])+1)>max_relative_residual)
                max_relative_residual=fabs(r)/(fabs(b->value[i])+1);
            if(fabs((double)x->value[i]-v[i])>maxerr)
                maxerr=fabs((double)x->value[i]-v[i]);
        }
        e=lis_psolveh(solver,bt,x);if(e)goto fail;
        for(i=0;i<N;i++)
        {
            double r=-bt->value[i];
            for(j=0;j<N;j++)r+=M[j][i]*scales[j]*x->value[j];
            if(fabs(r)/(fabs(bt->value[i])+1)>max_relative_residual)
                max_relative_residual=fabs(r)/(fabs(bt->value[i])+1);
            if(fabs((double)x->value[i]-v[i])>maxerr)
                maxerr=fabs((double)x->value[i]-v[i]);
        }
    }
    printf("SUPERLU_PLUGIN_TEST max_solution_error=%.12g max_scaled_residual=%.12g forward_and_transpose_reused=2\n",
           maxerr,max_relative_residual);
    e=(maxerr<1e-8 && max_relative_residual<1e-10)?
       LIS_SUCCESS:LIS_BREAKDOWN;
fail:
    if(e)fprintf(stderr,"SUPERLU_PLUGIN_TEST FAIL code=%d\n",(int)e);
    if(solver)solver->precon=NULL;
    if(pc)lis_precon_destroy(pc);
    if(solver)lis_solver_destroy(solver);
    if(x)lis_vector_destroy(x);
    if(bt)lis_vector_destroy(bt);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    lis_finalize();
    return e ? 1:0;
}
