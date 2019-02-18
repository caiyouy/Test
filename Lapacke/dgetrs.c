#include <stdio.h>
#include "/opt/lapack-3.8.0/include/lapacke.h"

int main (int argc, const char * argv[])
{
   lapack_int info;
   lapack_int m=2,n=2,lda=2;
   lapack_int ipiv[2];
   double a[2*2]={0.001,1.00,1.00,2.00};
   info=LAPACKE_dgetrf(LAPACK_COL_MAJOR, m, n, a, lda, ipiv);
   char tran='N'; // no transportation
   double b[2]={1.00,3.00};
   info=LAPACKE_dgetrs(LAPACK_COL_MAJOR,tran,n,1,a,m,ipiv,b,n);
   //output
   // printf("Matrix a:\n");
   // for(int i=0;i<4;i++)
   // {
   //    printf("%e\n",a[i]);
   // }
   // printf("ipiv: \n");
   // for(int i=0;i<2;i++)
   // {
   //    printf("%d\n",ipiv[i]);
   // }
   printf("b: \n");
   for(int i=0;i<2;i++)
   {
      printf("%e\n",b[i]);
   }
   return(info);
}

