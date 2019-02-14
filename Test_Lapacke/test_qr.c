#include <stdio.h>
#include <math.h>
#include "/opt/lapack-3.8.0/include/cblas.h"
#include "/opt/lapack-3.8.0/include/lapacke.h"

void NC(double *A,double *b,double *res)
{
    double A_normal[4],b_normal[2];
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                2,2,3,
                1.0,A,2,A,2,(double)0,A_normal,2);
    cblas_dgemv(CblasRowMajor,CblasTrans,
                3,2,
                1.0,A,2,b,1,(double)0,b_normal,1);
    // printf("A_normal:\n");
    // for(int i=0;i<4;i++)
    // {
    //     printf("%e\n",A_normal[i]);
    // }
    // printf("b_normal:\n");
    // for(int i=0;i<2;i++)
    // {
    //     printf("%e\n",b_normal[i]);
    // }
    lapack_int info;
    lapack_int ipiv[2];
    char uplo='L';
    info=LAPACKE_dpotrf(LAPACK_ROW_MAJOR,uplo,2,A_normal,2);
    info=LAPACKE_dpotrs(LAPACK_ROW_MAJOR,uplo,2,1,A_normal,2,b_normal,1);
    cblas_dcopy(2,b_normal,1,res,1);// move results 
}

void NG(double *A,double *b,double *res)
{
    double A_normal[4],b_normal[2];
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                2,2,3,
                1.0,A,2,A,2,(double)0,A_normal,2);
    cblas_dgemv(CblasRowMajor,CblasTrans,
                3,2,
                1.0,A,2,b,1,(double)0,b_normal,1);
    lapack_int info;
    lapack_int ipiv[2];
    char tran='N';
    info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR,2,2,A_normal,2,ipiv);
    info=LAPACKE_dgetrs(LAPACK_ROW_MAJOR,tran,2,1,A_normal,2,ipiv,b_normal,1);
    cblas_dcopy(2,b_normal,1,res,1);// move results 
}

void QR(double *A,double *b,double *res)
{
    lapack_int info;
    LAPACKE_dgels(LAPACK_ROW_MAJOR,'N',3,2,1,A,2,b,1);
    cblas_dcopy(2,b,1,res,1);// move results     
}

int main(int argc, char const *argv[])
{
    double k=1e9;
    double A[6]={k/sqrt(2)+0.5,-k/sqrt(6)+sqrt(3)/2.0,
                 -k/sqrt(2),   k/sqrt(6),
                 k/sqrt(2)-0.5,-k/sqrt(6)-sqrt(3)/2.0};
    double b[3]={3.0,2.0,-1.0};
    double x_exact[2]={1.0,sqrt(3)};
    // Cholesky method solving Normal Equations
    double resNC[2];
    NC(A,b,resNC);
    // Partially pivoted Guass method solving Normal Equation
    double resNG[3];
    NG(A,b,resNG);
    // QR method
    double resQR[3];
    QR(A,b,resQR);
    // Check results
    printf("Methods\t||x-x_exact||_2\n");
    cblas_daxpy(2,(double)-1,x_exact,1,resNC,1);   
    printf("NC:\t%e\n",cblas_dnrm2(2,resNC,1));
    cblas_daxpy(2,(double)-1,x_exact,1,resNG,1);   
    printf("NG:\t%e\n",cblas_dnrm2(2,resNG,1));
    cblas_daxpy(2,(double)-1,x_exact,1,resQR,1);   
    printf("QR:\t%e\n",cblas_dnrm2(2,resQR,1));
    return 0;
}
