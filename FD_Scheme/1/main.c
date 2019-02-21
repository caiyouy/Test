#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

double Func_initial(double x)
{
    if(fabs(x)<=1)
        return 1.0-fabs(x);
    else
        return 0;
}

void Lax_Friedrichs(double lb, double rb, int n, double CFL)
{
    //Space and Time step size
    double h=(rb-lb)/n;
    double k=CFL*h;
    double T=0.8;
    double *u=(double *)malloc((n+1)*sizeof(double));
    double *temp=(double *)malloc((n+1)*sizeof(double));    
    //Initialization
    for(int i=0;i<=n;i++)
    {
        u[i]=Func_initial(lb+h*i);
        // printf("%e\n",u[i]);
    }
    //Time Develop
    for(int i=0;i<T/k;i++)
    {
        //Boundary condition
        temp[0]=0; //lb: u=0
        for(int j=1;j<n;j++)
        {
            temp[j]=(u[j+1]+u[j-1])/2.0-CFL/2.0*(u[j+1]-u[j-1]); //LF
        }
        temp[n]=temp[n-1]; //rb: u[n+1]=u[n]
        memcpy(u,temp,(n+1)*sizeof(double));
    }
    //Output
    char name[20];
    sprintf(name, "LF_CFL%1.1lf.csv",CFL);
    FILE *fp=fopen(name,"w");
    for(int i=0;i<=n;i++)
    {
        fprintf(fp,"%e,%e\n",lb+i*h,u[i]);
    }
    fclose(fp);
    free(u); free(temp);
}

void Leapfrog(double lb, double rb, int n, double CFL)
{
    //Space and Time step size
    double h=(rb-lb)/n;
    double k=CFL*h;
    double T=0.8;
    double *u_pre=(double *)malloc((n+1)*sizeof(double));
    double *u_now=(double *)malloc((n+1)*sizeof(double));
    double *u_nex=(double *)malloc((n+1)*sizeof(double));
    //Initialization
    for(int i=0;i<=n;i++)
    {
        u_pre[i]=Func_initial(lb+h*i);
    }
    //Time Develop
    //First Time Develop;
    u_now[0]=0;
    for(int i=1;i<n;i++)
    {
        u_now[i]=u_pre[i]-CFL/2.0*(u_pre[i+1]-u_pre[i-1]);
    }
    u_now[n]=u_now[n-1];
    for(int i=1;i<T/k;i++)
    {
        //Boundary condition
        u_nex[0]=0; //lb: u=0
        for(int j=1;j<n;j++)
        {
            u_nex[j]=u_pre[j]-CFL*(u_now[j+1]-u_now[j-1]);
        }
        u_nex[n]=u_nex[n-1];
        memcpy(u_pre,u_now,(n+1)*sizeof(double));
        memcpy(u_now,u_nex,(n+1)*sizeof(double));
    }
    //Output
    char name[20];
    sprintf(name, "Leapfrog_CFL%1.1lf.csv",CFL);
    FILE *fp=fopen(name,"w");
    for(int i=0;i<=n;i++)
    {
        fprintf(fp,"%e,%e\n",lb+i*h,u_now[i]);
    }
    fclose(fp);
    free(u_pre);
    free(u_now);
    free(u_nex);
}

int main()
{
    // Computation Region [-2.0,3.0]
    double lb=-2.0;
    double rb= 3.0; 
    double n=50; // Number of nodes
    double CFL=0.8; // CFL number: k/h
    Lax_Friedrichs(lb,rb,n,CFL);
    CFL=1.6;
    Lax_Friedrichs(lb,rb,n,CFL);
    CFL=0.8;
    Leapfrog(lb,rb,n,CFL);
    return 0;
}