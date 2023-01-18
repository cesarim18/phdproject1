// It is not recommended to put function definitions  
// in a header file. Ideally there should be only 
// function declarations. Purpose of this code is 
// to only demonstrate working of header files. 
#ifndef _PERTURBATIONH_
#define _PERTURBATIONH_

#include "profiles.h"
#include "flux_extra.h"

void perturbation(double x[],double y[],double z[],double a[],double th[],double da[],int dim_s,int dim_th,double *R,double *A0)
{
	int i_star,j_star,k_star;
	j_star = 0.5*(dim_s-4) + 2;
	k_star = 0.25*(dim_th-4) + 2;
	i_star = j_star + dim_s*k_star;
	printf("	(j,k) = (%d,%d)\n",j_star,k_star);
	double th_k,a_j;
	a_j = a[j_star];
	th_k = th[k_star];

	double xp,yp,zp,rp;
	rp = R[i_star];
	xp = x[j_star] - rp*sin(th_k)*sin(a_j);
	yp = y[j_star] + rp*cos(th_k);
	zp = z[j_star] + rp*sin(th_k)*cos(a_j);
	printf("	(xp,yp,zp,rp) = (%e,%e,%e,%e)\n",xp,yp,zp,rp);
	
	int i,j,k;
	double r_i,da_j,r_pert,Vx,Vy,Vz;
	for(k=0;k<dim_th;k++){
		for(j=0;j<dim_s;j++){
			i = j+dim_s*k;
			r_i = R[i];
			a_j = a[j];
			da_j = da[j];
			th_k = th[k];

			Vx = x[j] - r_i*sin(th_k)*sin(a_j);
			Vy = y[j] + r_i*cos(th_k);
			Vz = z[j] + r_i*sin(th_k)*cos(a_j);

			r_pert = sqrt(0.25*pow(Vx-xp,2.0) + pow(Vy-yp,2.0) + pow(Vz-zp,2.0))/rp;
			if(r_pert < 1.0){
				r_i = (1.0 + 0.5*sin(0.5*M_PI*(1.0-r_pert)))*r_i;
				printf("	perturbation at (j,k) = (%d,%d), (Vx,Vy,Vz) = (%e,%e,%e),r_pert = %.15e,th_k = %e\n",j,k,Vx,Vy,Vz,r_pert,th_k);
			}A0[i] = Q0(1.15*r_i,1.15*r_i*sin(th_k)*da_j);
		}	
	}
} 

void aneurysmG(double x[],double y[],double z[],double a[],double th[],int dim_s,int dim_th,double *R,double *G,double xp,double yp,double zp,double rp)
{
	double th_k,a_j;
	int i,j,k;
	double r_i,r_pert,Vx,Vy,Vz,G_i;

	for(k=0;k<dim_th;k++){
		for(j=0;j<dim_s;j++){
			i = j+dim_s*k;
			r_i = R[i];
			G_i = G[i];
			a_j = a[j];
			th_k = th[k];

			Vx = x[j] - r_i*sin(th_k)*sin(a_j);
			Vy = y[j] + r_i*cos(th_k);
			Vz = z[j] + r_i*sin(th_k)*cos(a_j);

			r_pert = sqrt(pow(Vx-xp,2.0) + pow(Vy-yp,2.0) + pow(Vz-zp,2.0))/rp;
			if(r_pert <= 1.0) G_i = (1.0 + 0.5*sin(0.5*M_PI*(1.0-r_pert)))*G_i;
			G[i] = G_i;
		}
	}
}

void aneurysmR(double x[],double y[],double z[],double a[],double th[],int dim_s,int dim_th,double *R,double xp,double yp,double zp,double rp)
{
	double th_k,a_j;
	int i,j,k;
	double r_i,r_pert,Vx,Vy,Vz;

	for(k=0;k<dim_th;k++){
		for(j=0;j<dim_s;j++){
			i = j+dim_s*k;
			r_i = R[i];
			a_j = a[j];
			th_k = th[k];

			Vx = x[j] - r_i*sin(th_k)*sin(a_j);
			Vy = y[j] + r_i*cos(th_k);
			Vz = z[j] + r_i*sin(th_k)*cos(a_j);

			r_pert = sqrt(pow(Vx-xp,2.0) + pow(Vy-yp,2.0) + pow(Vz-zp,2.0))/rp;
			if(r_pert <= 1.0) r_i = (1.0 - 0.75*sin(0.5*M_PI*(1.0-r_pert)))*r_i;
			R[i] = r_i;
		}	
	}
} 

void bulgingR(int dim_s,int dim_th,double *R,double s[],double Ls)
{
	int i,j,k;
	double s_j,r_i,pert,Lim1,Lim2,Lim3;
	Lim1 = 0.05*Ls;
	Lim2 = 0.2*Ls;
	Lim3 = 0.4*Ls;
	for(k=0;k<dim_th;k++){
		for(j=0;j<dim_s;j++){
			i = j+dim_s*k;
			s_j = s[j];
			r_i = R[i];

			if(Lim1 <= s_j && s_j <= Lim2){
				pert = 1.0 - (Lim2 - s_j)/(Lim2 - Lim1);
				r_i = (1.0 + 0.5*sin(0.5*M_PI*pert))*r_i;
			}if(Lim2 < s_j && s_j <= Lim3){
				pert = 1.0 - (s_j - Lim2)/(Lim3 - Lim2);
				r_i = (1.0 + 0.5*sin(0.5*M_PI*pert))*r_i;
			}R[i] = r_i;
		}	
	}
} 

void oclussionR(int dim_s,int dim_th,double *R,double s[],double Ls)
{
	int i,j,k;
	double s_j,r_i,pert,Lim1,Lim2,Lim3;
	Lim1 = 0.05*Ls;
	Lim2 = 0.2*Ls;
	Lim3 = 0.4*Ls;
	for(k=0;k<dim_th;k++){
		for(j=0;j<dim_s;j++){
			i = j+dim_s*k;
			s_j = s[j];
			r_i = R[i];

			if(Lim1 <= s_j && s_j <= Lim2){
				pert = 1.0 - (Lim2 - s_j)/(Lim2 - Lim1);
				r_i = (1.0 - 0.5*sin(0.5*M_PI*pert))*r_i;
			}if(Lim2 < s_j && s_j <= Lim3){
				pert = 1.0 - (s_j - Lim2)/(Lim3 - Lim2);
				r_i = (1.0 - 0.5*sin(0.5*M_PI*pert))*r_i;
			}R[i] = r_i;
		}	
	}
} 


void bulgingG(int dim_s,int dim_th,double *G,double s[],double Ls)
{
	int i,j,k;
	double s_j,G_i,pert,Lim1,Lim2,Lim3;
	Lim1 = 0.2*Ls;
	Lim2 = 0.35*Ls;
	Lim3 = 0.5*Ls;
	for(k=0;k<dim_th;k++){
		for(j=0;j<dim_s;j++){
			i = j+dim_s*k;
			s_j = s[j];
			G_i = G[i];

			if(Lim1 <= s_j && s_j <= Lim2){
				pert = 1.0 - (Lim2 - s_j)/(Lim2 - Lim1);
				G_i = (1.0 + 0.5*sin(0.5*M_PI*pert))*G_i;
			}if(Lim2 < s_j && s_j <= Lim3){
				pert = 1.0 - (s_j - Lim2)/(Lim3 - Lim2);
				G_i = (1.0 + 0.5*sin(0.5*M_PI*pert))*G_i;
			}G[i] = G_i;
		}	
	}
} 
#endif
