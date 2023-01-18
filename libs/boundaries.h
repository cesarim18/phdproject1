// It is not recommended to put function definitions  
// in a header file. Ideally there should be only 
// function declarations. Purpose of this code is 
// to only demonstrate working of header files. 
#ifndef _BOUNDARIESH_
#define _BOUNDARIESH_

#include "flux_extra.h"
// ext 61040
void leftboundaries(double *W1,double *W2,double *W3,double Gh1[],double Gh2[],
double Gh3[],double *R0,double *G0,double th[],double da[],int D_th,int D_s,
double beta,double gs,double gt,double u_gh,int cont[],double *Ao)
{
	int j,k;
	j = 2;

	double PAs,PA0[3],PA1[2];
	PAs = IVs(gs,1.0);
	PA0[0] = IVs(gs,1.0);
	PA0[1] = 2.0*IVs(gs,2.0);
	PA0[2] = IVs(gs,3.0);
	PA1[0] = IVsVs(gs,1.0);
	PA1[1] = IVsVs(gs,2.0);

	double PAt[2],PB0[2],PB1;
	PAt[0] = IVt(gt,1.0); 
	PAt[1] = IVt(gt,2.0);
	PB0[0] = IVt(gt,3.0);
	PB0[1] = IVt(gt,4.0);
	PB1 = IVsVt(gs,gt,3.0);

	int i;
	double Qr,Qs,Qt,maxL,minL,V[15],F[3];
	double R,Ga,VAs,VA0,VAt,VB0;

	for(k=0;k<D_th;k++){
		i = j + D_s*k;
		Qr = ga(W1[i],R0[i],sin(th[k]),da[j]);
		Qs = W2[i];
		Qt = W3[i];
		Variables(V,Qr,Qs,Qt,R0[i],Ao[i],G0[i],th[k],da[j],PAs,PAt,PA0,PB0,beta);
		Axial_Req(F,&maxL,&minL,V,PA1,PB1,PAs,PA0);
		
		if(maxL >= 0.0){
			// if (val == TRUE) Inflow
			Gh1[k] = 1.0;				// Dirichlet for R/R_o
			Gh2[k] = u_gh;				// Dirichlet for u
			Gh3[k] = 0.0;				// Dirichlet for u
			cont[k] = 1;
		}else{
			// if (val != TRUE) Outflow
			Gh1[k] = Qr;				// Neumann for R/R_o
			Gh2[k] = V[6];				// Neumann for u
			Gh3[k] = 0.0;				// Neumann for w
			cont[k] = 0;
		}

		R = Gh1[k]*R0[i-2];
		Ga = R*sin(th[k])*da[j-2];
		VAs = As(R,Ga,PAs);
		VA0 = A0(R,Ga,PA0);
		VAt = At(R,Ga,PAt);
		VB0 = B0(R,Ga,PB0);

		W1[i-2] = Q0(R,Ga);
		W2[i-2] = (W1[i-2]/VAs)*VA0*Gh2[k];
		W3[i-2] = (W1[i-2]/VAt)*VB0*Gh3[k];

		R = Gh1[k]*R0[i-1];
		Ga = R*sin(th[k])*da[j-1];
		VAs = As(R,Ga,PAs);
		VA0 = A0(R,Ga,PA0);
		VAt = At(R,Ga,PAt);
		VB0 = B0(R,Ga,PB0);

		W1[i-1] = Q0(R,Ga);
		W2[i-1] = (W1[i-1]/VAs)*VA0*Gh2[k];
		W3[i-1] = (W1[i-1]/VAt)*VB0*Gh3[k];
	}	
}

void rightboundaries(double *W1,double *W2,double *W3,double Gh1[],double Gh2[],
double Gh3[],double *R0,double *G0,double th[],double da[],int D_th,int D_s,
double beta,double gs,double gt,int cont[],double *Ao)
{
	int j,k;
	j = D_s-3;

	double PAs,PA0[3],PA1[2];
	PAs = IVs(gs,1.0);
	PA0[0] = IVs(gs,1.0);
	PA0[1] = 2.0*IVs(gs,2.0);
	PA0[2] = IVs(gs,3.0);
	PA1[0] = IVsVs(gs,1.0);
	PA1[1] = IVsVs(gs,2.0);

	double PAt[2],PB0[2],PB1;
	PAt[0] = IVt(gt,1.0); 
	PAt[1] = IVt(gt,2.0);
	PB0[0] = IVt(gt,3.0);
	PB0[1] = IVt(gt,4.0);
	PB1 = IVsVt(gs,gt,3.0);

	int i;
	double Qr,Qs,Qt,maxL,minL,V[15],F[3];
	double R,Ga,VAs,VA0,VAt,VB0;
	for(k=0;k<D_th;k++){
		i = j + D_s*k;
		Qr = ga(W1[i],R0[i],sin(th[k]),da[j]);
		Qs = W2[i];
		Qt = W3[i];
		Variables(V,Qr,Qs,Qt,R0[i],Ao[i],G0[i],th[k],da[j],PAs,PAt,PA0,PB0,beta);
		Axial_Req(F,&maxL,&minL,V,PA1,PB1,PAs,PA0);
		
		if(maxL >= 0.0){
			// if (val == TRUE) Outflow
			Gh1[k] = Qr;		// Neumann for R/R_o
			Gh2[k] = V[6];		// Neumann for u
			Gh3[k] = V[10];		// Neumann for w
			cont[k] = 0;
		}else{
			// if (val != TRUE) Inflow
			Gh1[k] = 1.0;		// Dirichlet for R/R_o
			Gh2[k] = 0.0;		// Dirichlet for u
			Gh3[k] = 0.0;		// Dirichlet for w
			cont[k] = 1;
		}

		R = Gh1[k]*R0[i+1];
		Ga = R*sin(th[k])*da[j+1];
		VAs = As(R,Ga,PAs);
		VA0 = A0(R,Ga,PA0);
		VAt = At(R,Ga,PAt);
		VB0 = B0(R,Ga,PB0);

		W1[i+1] = Q0(R,Ga);
		W2[i+1] = (W1[i+1]/VAs)*VA0*Gh2[k];
		W3[i+1] = (W1[i+1]/VAt)*VB0*Gh3[k];

		R = Gh1[k]*R0[i+2];
		Ga = R*sin(th[k])*da[j+2];
		VAs = As(R,Ga,PAs);
		VA0 = A0(R,Ga,PA0);
		VAt = At(R,Ga,PAt);
		VB0 = B0(R,Ga,PB0);

		W1[i+2] = Q0(R,Ga);
		W2[i+2] = (W1[i+2]/VAs)*VA0*Gh2[k];
		W3[i+2] = (W1[i+2]/VAt)*VB0*Gh3[k];
	}	
}

void periodicboundaries(double *W,int dim_s,int dim_th)
{
	int j;
	#pragma omp parallel for
	for(j=0;j<dim_s-4;j++){
		int i0,i1,i3m,i4m;
		int i2,i3,i1m,i2m;
		 i0 = (j+2) + 0*dim_s;
		 i1 = (j+2) + 1*dim_s;
		i3m = (j+2) + (dim_th-3)*dim_s;
		i4m = (j+2) + (dim_th-4)*dim_s;
		W[i0] = W[i4m];
		W[i1] = W[i3m];

		 i2 = (j+2) + 2*dim_s;
		 i3 = (j+2) + 3*dim_s;
		i1m = (j+2) + (dim_th-1)*dim_s;
		i2m = (j+2) + (dim_th-2)*dim_s;
		W[i1m] = W[i3];
		W[i2m] = W[i2];
	}
}

void free_eastboundaries(double *v,int dim_th,int dim_s)
{
	int i,k;
	for(k=0;k<dim_th;k++){
		i = 2 + k*dim_s;
		v[i-1] = v[i];
		v[i-2] = v[i];
	}
}

void free_westboundaries(double *v,int dim_th,int dim_s)
{
	int i,k;
	for(k=0;k<dim_th;k++){
		i = (dim_s-3) + k*dim_s;
		v[i+1] = v[i];
		v[i+2] = v[i];
	}
}

void free_northboundaries(double *v,int dim_th,int dim_s)
{
	int i,j;
	for(j=0;j<dim_s;j++){
		i = j + (dim_th-3)*dim_s;
		v[i+dim_s] = v[i];
		v[i+2*dim_s] = v[i];
	}
}

void free_southboundaries(double *v,int dim_th,int dim_s)
{
	int i,j;
	for(j=0;j<dim_s;j++){
		i = j + 2*dim_s;
		v[i-dim_s] = v[i];
		v[i-2*dim_s] = v[i];
	}
}

#endif
