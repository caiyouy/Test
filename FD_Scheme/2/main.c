#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

double Func_initial(double x)
{
    if(x<=0)
        return 1.0;
    else
        return 0.0;
}

double weno5_FV(double *u)
{
    double stencil_A[3]={ 1.0/3.0,-7.0/6.0, 11.0/6.0};
    double stencil_B[3]={-1.0/6.0, 5.0/6.0,  1.0/3.0};
    double stencil_C[3]={ 1.0/3.0, 5.0/6.0, -1.0/6.0};
    double gamma[3]={1.0/10.0, 3.0/5.0, 3.0/10.0};
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
    double stencil_A[3]={ 1.0/3.0,-7.0/6.0, 11.0/6.0};
    double stencil_B[3]={-1.0/6.0, 5.0/6.0,  1.0/3.0};
    double stencil_C[3]={ 1.0/3.0, 5.0/6.0, -1.0/6.0};
    double gamma[3]={1.0/10.0, 3.0/5.0, 3.0/10.0};
    // Approximation from three different stencils
    double res_A=stencil_A[0]*u[0]+stencil_A[1]*u[1]+stencil_A[2]*u[2];
    double res_B=stencil_B[0]*u[1]+stencil_B[1]*u[2]+stencil_B[2]*u[3];
    double res_C=stencil_C[0]*u[2]+stencil_C[1]*u[3]+stencil_C[2]*u[4];
    // res
    return (gamma[0]*res_A+gamma[1]*res_B+gamma[2]*res_C);
}

double gauss_quadrature(double (*Func)(double),double lb,double ub)
{
    double points[3]={-sqrt(3.0/5.0),0,sqrt(3.0/5.0)};
    double weights[3]={5.0/9.0,8.0/9.0,5.0/9.0};
    // Linear transformation
    double a=(ub-lb)/2.0;
    double b=(ub+lb)/2.0;
    double res=0;
    for(int i=0;i<3;i++)
    {
        res+=Func(a*points[i]+b)*weights[i];
    }
    return res*a;
}

int main()
{
    // Computation Region [-0.5,0.5]
    double n=50; // Number of nodes
    double h=1.0/n;
    double *u=(double *)malloc(n*sizeof(double));
    //Initialization
    for(int i=0;i<n;i++)
    {
        double lb=-0.5+i*h;
        u[i]=gauss_quadrature(&Func_initial,lb,lb+h)/h;
        // printf("%lf\n",u[i]);
    }
    //Reconstruction for u_{i+0.5}
    FILE *fp=fopen("res","w");
    for(int i=3;i<n-1;i++) // Simplify for boundary
    {
        double in[5];
        memcpy(in,u+i-3,5*sizeof(double));
        double res_WENO=weno5_FV(in);
        double res_tradition=tradition_FV(in);
        fprintf(fp,"%e %e %e\n",-0.5+i*h,res_WENO,res_tradition);
    }
    fclose(fp);
    return 0;
}