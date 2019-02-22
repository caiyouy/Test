#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

double weno5_FV(double *u)
{
    double stencil_A[3]={-1.0/24, 1.0/12, 23.0/24};
    double stencil_B[3]={-1.0/24,13.0/12, -1.0/24};
    double stencil_C[3]={23.0/24, 1.0/12, -1.0/24};
    double gamma[3]={-9.0/80,49.0/40,-9.0/80};
    double epsilon=1e-6;
    // Approximation from three different stencils
    double res_A=stencil_A[0]*u[0]+stencil_A[1]*u[1]+stencil_A[2]*u[2];
    double res_B=stencil_B[0]*u[1]+stencil_B[1]*u[2]+stencil_B[2]*u[3];
    double res_C=stencil_C[0]*u[2]+stencil_C[1]*u[3]+stencil_C[2]*u[4];
    // Compute smoothness indicators
    double beta_A=13.0/12.0*pow(u[0]-2.0*u[1]+u[2],2)+1.0/4.0*pow(u[0]-4.0*u[1]+3*u[2],2);
    double beta_B=13.0/12.0*pow(u[1]-2.0*u[2]+u[3],2)+1.0/4.0*pow(u[1]-u[3],2);
    double beta_C=13.0/12.0*pow(u[2]-2.0*u[3]+u[4],2)+1.0/4.0*pow(3.0*u[2]-4.0*u[3]+u[4],2);
    // Compute weights
    double tmp_A=gamma[0]/pow(epsilon+beta_A,2);
    double tmp_B=gamma[1]/pow(epsilon+beta_B,2);
    double tmp_C=gamma[2]/pow(epsilon+beta_C,2);
    double tmp=tmp_A+tmp_B+tmp_C;
    double weight_A=tmp_A/tmp;
    double weight_B=tmp_B/tmp;
    double weight_C=tmp_C/tmp;
    // res
    return (weight_A*res_A+weight_B*res_B+weight_C*res_C);
}

//Fixed stencil (Five cells)
double tradition_FV(double *u)
{
    double stencil_A[3]={-1.0/24, 1.0/12, 23.0/24};
    double stencil_B[3]={-1.0/24,13.0/12, -1.0/24};
    double stencil_C[3]={23.0/24, 1.0/12, -1.0/24};
    double gamma[3]={-9.0/80,49.0/40,-9.0/80};
    // Approximation from three different stencils
    double res_A=stencil_A[0]*u[0]+stencil_A[1]*u[1]+stencil_A[2]*u[2];
    double res_B=stencil_B[0]*u[1]+stencil_B[1]*u[2]+stencil_B[2]*u[3];
    double res_C=stencil_C[0]*u[2]+stencil_C[1]*u[3]+stencil_C[2]*u[4];
    // res
    return (gamma[0]*res_A+gamma[1]*res_B+gamma[2]*res_C);
}

double initial(double lb,double ub)
{
    if(ub<=0) //interval in 2x
        return (ub-lb)*(ub+lb);
    if(lb>0)
        return -20*(ub-lb);
    return -lb*lb-20*ub;
}

int main()
{
    // Computation Region approximate [-0.5,0.5]
    double n=50; // Number of nodes
    double h=1.0/n;
    double left=(-24-0.4965-0.5)*h;
    double *u=(double *)malloc(n*sizeof(double));
    //Initialization
    for(int i=0;i<n;i++)
    {
        double lb=left+i*h;
        u[i]=initial(lb,lb+h)/h;
        // printf("%lf\n",u[i]);
    }
    //Reconstruction for u_{i}
    FILE *fp=fopen("res","w");
    for(int i=2;i<n-2;i++) // Simplify for boundary
    {
        double in[5];
        memcpy(in,u+i-2,5*sizeof(double));
        double res_WENO=weno5_FV(in);
        double res_tradition=tradition_FV(in);
        fprintf(fp,"%e %e %e\n",left+(i+0.5)*h,res_WENO,res_tradition);
    }
    fclose(fp);
    return 0;
}
