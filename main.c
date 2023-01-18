/*
	Casa Presidenta - No. 

	Blood Flow Simulation

	gcc main.c -fopenmp -lm -o salida
	export OMP_NUM_THREADS=8
	time ./salida

	New:
	clang main.c -fopenmp -lm -o salida
	
	MODIFICA LA VARIABLE FILE EN CADA EJECUCION QUE REALICES
	PARA QUE TE GUARDE POR SEPARADO LOS ARCHIVOS
*/

// LIBRARIES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <complex.h>
#include <sys/stat.h> 
#include <sys/types.h>
#include <omp.h>
#include "libs/print_files.h"
#include "libs/print_initial_condition.h"
#include "libs/print_final_condition.h"
#include "libs/flux.h"
#include "libs/flux_extra.h"
#include "libs/find_root.h"
#include "libs/profiles.h"
#include "libs/initial_R0G0_V2.h"
#include "libs/carpets.h"
#include "libs/boundaries.h"
#include "libs/post_process.h"
#include "libs/perturbation.h"
#include "libs/create_vectors.h"
#include "libs/cardiac_cycle.h"

void initial_cond(double [],double [],int ,int ,double *,double *);
void initial_cond2(double ,double [],int ,int ,double *,double *,double );
void heatequation(double [],int );
void printfrc(double [],double [],int ,int );
void printfarc(double [],double [],double [],int ,int );
void arc(double [],double [],double ,int );
void diffarc1(double [],double [],int ,double );
void diffarc2(double [],double [],int ,double );
void print_inflow(int [],int ,int );
void print_outflow(int [],int ,int );
void Imprimir_Boundaries(int [],int [],int ,int );
void sum_z(double [],double [],int, double );
void print_RoGo(double *,double *,int ,int ,int );

int main()
{
	int nthreads;
	nthreads = omp_get_max_threads();
	printf("Using %u OpenMP threads.\n",nthreads);

	int file;
	file = 1;
	create_folders(file);

	int N_s,N_th,movies;
	N_s = 100;
	N_th = 60;
	movies = 100; // In frames per Tfinal units

	double Tfinal;
	Tfinal = 2.0;// In seconds

	// Eccentricity for geometry
	double ex1,ex2,eps,ex0;
	ex0 = 0.5;
	ex1 = 0.0;
	ex2 = 0.0;
	eps = 1.0;//0.01;0.0001;

	int heat,geometry;
	geometry = 1; // 0 for aorta, 1 for horizontal, 2 for twisted aorta
	double Ls,barrier;
	switch(geometry)
	{
		case 0:{
			heat = 10;
			Ls = (7.0357+0.8+0.9+6.4737+15.2+1.8+0.7+0.7+4.3+4.3)/100.0; // In m
			printf("Ls = %.15e\n",Ls);
			barrier = 0.3*Ls;
		}break;
		case 1:{
			Ls = 50.0/100.0; // In m
			printf("Ls = %.15e\n",Ls);
		}break;
		case 2:{
			Ls = 100.0/100.0; // In m
		}break;
	}
	if(geometry == 0){
		heat = 10;
		Ls = (7.0357+0.8+0.9+6.4737+15.2+1.8+0.7+0.7+4.3+4.3)/100.0; // In m
		printf("Ls = %.15e\n",Ls);
		barrier = 0.3*Ls;
	}if(geometry == 1){
	heat = 10;

	Ls = (7.0357+0.8+0.9+6.4737+15.2+1.8+0.7+0.7+4.3+4.3)/100.0; // In m
	printf("Ls = %.15e\n",Ls);
	barrier = 0.3*Ls;

	int tapering,chipoteG,chipoteR,seccionR,seccionG;
	tapering = 1;
	chipoteG = 0;
	chipoteR = 0;
	seccionR = 0;
	seccionG = 0;

	double g_s,g_th;
	g_s = 9.0;
	g_th = 2.0;

	int pert;
	pert = 1;

	double parVin;
	parVin = 1.0;

	// Scales
	double Length_char,Radius_char,Angle_char;
	Length_char = 1.0;				// In m
	Radius_char = 0.01;				// In m
	Angle_char = 1.0;				// In rad

	double Axial_char,Radial_char,Angular_char;
	Axial_char = 1.0;				// In m/seg
	Radial_char = 0.01;				// In m/seg
	Angular_char = 1.0;				// In rad/seg

	double Time_char;
	Time_char = 1.0;				// In seg

	//--------------------------------------------------------------------------

	// index for loops
	int i,j,k;
	int js,ks,ic;
	int bt = 0;

	double max_a_0,max_b_0,min_R_0,max_a_1,max_b_1,min_R_1;
	double A1,Us1,cs1,Ga1;

	double u_gh,v_inflow;
	u_gh = 0.0;

	double max_a,max_b,min_R,wtime;
	max_a = 0.0;
	max_b = 0.0;
	min_R = 0.0;
	wtime = 0.0;

	// parameters
	double beta,rho,nu,g;
	beta = 2.0;
	rho = 1050.0;	 				// In Kg/m^3
	nu  = 1000.0*2.6*pow(10.0,-6.0)*rho;	// In m^2/seg  // MODIFICAR ESTE VALOR
	g = 9.81;

	// domain parameters
	double ds,dth;
	ds = (1.0/N_s)*Ls;
	dth = (1.0/N_th)*2.0*M_PI;

	int D_s,D_th;
	D_s = N_s+4;
	D_th = N_th+4;

	double W1_gh[D_th];
	double W2_gh[D_th];
	double W3_gh[D_th];

	// time variables
	double dt,dt_0,dt_1,T;
	dt = 0.0;
	dt_0 = 0.0;
	dt_1 = 0.0;
	T = 0.0;

	// Frames
	int savedata;
	savedata = 0;

	// numerical scheme parameters
	double CFL,v_mm,fac_visc;
	CFL = 0.45;
	v_mm = 1.5;
	fac_visc = 1.0;

	double par[13];

	// numerical scheme initiation
	par[0] = 0.0;
	par[1] = 0.0;
	par[2] = 0.0;
	par[3] = v_mm;
	par[4] = ds;
	par[5] = dth;
	par[6] = beta;
	par[7] = g;
	par[8] = fac_visc;
	par[9] = eps;
	par[10] = g_s;
	par[11] = g_th;
	par[12] = nu/rho;

	double max_wtime = 0.0;
	double min_wtime = 1.0;
	double sum_wtime = 0.0;

	double th0[D_th],thp[D_th],thm[D_th],thA[2*D_th+1];
	for(k=0;k<2*D_th+1;k++){
		thA[k] = -2.0*dth + k*0.5*dth;
	}for(k=0;k<D_th;k++){
		thm[k] = thA[2*k];
		th0[k] = thA[2*k+1];
		thp[k] = thA[2*(k+1)];
	}

	//--------------------------------------------------------------------------
	// Creating x and z variables
	//--------------------------------------------------------------------------
	double sA[2*D_s+1],alA[2*D_s+1],daA[2*D_s+1],ddA[2*D_s+1];
	for(j=0;j<2*D_s+1;j++){
		 sA[j] = (0.5*j-2.0)*ds;
		alA[j] = 0.0;
		daA[j] = 0.0;
		ddA[j] = 0.0;
	}

	if(geometry == 0){
		// aorta conditions
		int elec = 2*(3*N_s/10 + 2);
		printf("%d\n",elec);
		sA[elec] = 0.3*Ls;
		arc(alA,sA,barrier,2*D_s+1);
		diffarc1(daA,alA,2*D_s+1,0.5*ds);
		daA[elec] = daA[elec-1];
		diffarc2(ddA,alA,2*D_s+1,0.5*ds);
		printfarc(alA,daA,ddA,2*D_s+1,file);		
		for(j=0;j<heat;j++){
			heatequation(alA,2*D_s+1);
			diffarc1(daA,alA,2*D_s+1,0.5*ds);
			diffarc2(ddA,alA,2*D_s+1,0.5*ds);
			printfarc(alA,daA,ddA,2*D_s+1,file);
		}
		printf("after heat equation\n");
		for(j=0;j<2*D_s+1;j++){
			if(j<10) printf("      j = %d , ",j);
			if(j>9 && j<100) printf("     j = %d , ",j);
			if(j>99 &&  j<1000) printf("    j = %d , ",j);

			if(sA[j]<0.0) printf("%.15e, ",sA[j]/Ls);
			else printf(" %.15e, ",sA[j]/Ls);
			if(alA[j]<0.0) printf("%.15e, ",alA[j]);
			else printf(" %.15e, ",alA[j]);
			if(daA[j]<0.0) printf("%.15e, ",daA[j]);
			else printf(" %.15e, ",daA[j]);
			if(ddA[j]<0.0) printf("%.15e\n",ddA[j]);
			else printf(" %.15e\n",ddA[j]);
		}
	}
	//return 0;
	double s0[D_s],sp[D_s],sm[D_s];
	double al0[D_s],da0[D_s],dd0[D_s];
	double alp[D_s],dap[D_s],ddp[D_s];
	double alm[D_s],dam[D_s],ddm[D_s];
	for(j=0;j<D_s;j++){
		 sm[j] =  sA[2*j];
		alm[j] = alA[2*j];
		dam[j] = daA[2*j];
		ddm[j] = ddA[2*j];

		 s0[j] =  sA[2*j+1];
		al0[j] = alA[2*j+1];
		da0[j] = daA[2*j+1];
		dd0[j] = ddA[2*j+1];

		 sp[j] =  sA[2*(j+1)];
		alp[j] = alA[2*(j+1)];
		dap[j] = daA[2*(j+1)];
		ddp[j] = ddA[2*(j+1)];
	}

	double vec_r0[D_s],vec_c0[D_s];
	double vec_rm[D_s],vec_cm[D_s];
	double vec_rp[D_s],vec_cp[D_s];
	initial_radius_promvel(s0,D_s,tapering,vec_r0,vec_c0);
	initial_radius_promvel(sm,D_s,tapering,vec_rm,vec_cm);
	initial_radius_promvel(sp,D_s,tapering,vec_rp,vec_cp);

	for(j=0;j<D_s;j++){
		printf("%.15e, %.15e, %.15e, ",sm[j],vec_rm[j],vec_cm[j]);
		printf("%.15e, %.15e, %.15e, ",s0[j],vec_r0[j],vec_c0[j]);
		printf("%.15e, %.15e, %.15e\n",sp[j],vec_rp[j],vec_cp[j]);
	}//return 0;

	// for post processing
	double x_arc[2*D_s+1],y_arc[2*D_s+1],z_arc[2*D_s+1];
	x_arc[4] = sA[4];
	y_arc[4] = 0.0;
	z_arc[4] = 0.0;
	for(j=5;j<2*D_s+1;j++){
		x_arc[j] = x_arc[j-1] + 0.5*ds*cos(alA[j]);
		y_arc[j] = y_arc[j-1];
		z_arc[j] = z_arc[j-1] + 0.5*ds*sin(alA[j]);
	}for(j=3;j>-1;j--){
		x_arc[j] = x_arc[j+1] - 0.5*ds*cos(alA[j]);
		y_arc[j] = y_arc[j+1];
		z_arc[j] = z_arc[j+1] - 0.5*ds*sin(alA[j]);
	}

	double x_arc_m[D_s],y_arc_m[D_s],z_arc_m[D_s];
	double x_arc_0[D_s],y_arc_0[D_s],z_arc_0[D_s];
	double x_arc_p[D_s],y_arc_p[D_s],z_arc_p[D_s];
	for(j=0;j<D_s;j++){
		x_arc_m[j] = x_arc[2*j];
		y_arc_m[j] = y_arc[2*j];
		z_arc_m[j] = z_arc[2*j];

		x_arc_0[j] = x_arc[2*j+1];
		y_arc_0[j] = y_arc[2*j+1];
		z_arc_0[j] = z_arc[2*j+1];

		x_arc_p[j] = x_arc[2*(j+1)];
		y_arc_p[j] = y_arc[2*(j+1)];
		z_arc_p[j] = z_arc[2*(j+1)];
	}

	// -------------------------------------------------------------------------
	// wall elasticity properties
	// Young modulus reference per s and th
	int left_cont[D_th],right_cont[D_th];
	double *Ro,*RoW,*RoE,*RoS,*RoN;
	Crea1(&Ro,D_s*D_th);
	Crea1(&RoW,D_s*D_th);
	Crea1(&RoE,D_s*D_th);
	Crea1(&RoS,D_s*D_th);
	Crea1(&RoN,D_s*D_th);

	double *Go,*GoW,*GoE,*GoS,*GoN;
	Crea1(&Go,D_s*D_th);
	Crea1(&GoW,D_s*D_th);
	Crea1(&GoE,D_s*D_th);
	Crea1(&GoS,D_s*D_th);
	Crea1(&GoN,D_s*D_th);

	R0G0_elipse(RoE,GoE,vec_rp,vec_cp,th0,D_s,D_th,ex0);
	R0G0_elipse(RoW,GoW,vec_rm,vec_cm,th0,D_s,D_th,ex0);
	R0G0_elipse(RoN,GoN,vec_r0,vec_c0,thp,D_s,D_th,ex0);
	R0G0_elipse(RoS,GoS,vec_r0,vec_c0,thm,D_s,D_th,ex0);

	for(i=0;i<D_s*D_th;i++){
		//RoW[i] = 1.2*RoW[i];
		//RoE[i] = 1.2*RoE[i];
		//RoN[i] = 1.2*RoN[i];
		//RoS[i] = 1.2*RoS[i];
		
		Ro[i] = 0.0;
		Ro[i] = 0.25*(RoW[i]+RoE[i]+RoS[i]+RoN[i]);
		Go[i] = 0.0;
		Go[i] = 0.25*(GoW[i]+GoE[i]+GoS[i]+GoN[i]);
	}
	periodicboundaries(Go,D_s,D_th);
	periodicboundaries(Ro,D_s,D_th);

	double min_R_i;
	min_R_i = min_vec(Ro,D_s*D_th);
	printf("min Ro = %.15e\n",min_R_i);
	if(min_R_i == 0.0) return 0;
	print_RoGo(Go,Ro,D_s,D_th,file);
	//return 0;

	if(chipoteG == 1){
		double xp,yp,zp,rp;
		int s_p,th_p;
		s_p = N_s/8+2;
		th_p = N_th/2+5;
		rp = Ro[s_p + D_s*th_p];
		xp = x_arc_0[s_p] - rp*sin(th0[th_p])*sin(al0[s_p]);
		yp = y_arc_0[s_p] + rp*cos(th0[th_p]);
		zp = z_arc_0[s_p] + rp*sin(th0[th_p])*cos(al0[s_p]);

		aneurysmG(x_arc_p,y_arc_p,z_arc_p,alp,th0,D_s,D_th,RoW,GoW,xp,yp,zp,rp);
		aneurysmG(x_arc_m,y_arc_m,z_arc_m,alm,th0,D_s,D_th,RoE,GoE,xp,yp,zp,rp);
		aneurysmG(x_arc_0,y_arc_0,z_arc_0,al0,thp,D_s,D_th,RoN,GoN,xp,yp,zp,rp);
		aneurysmG(x_arc_0,y_arc_0,z_arc_0,al0,thm,D_s,D_th,RoS,GoS,xp,yp,zp,rp);

		periodicboundaries(GoW,D_s,D_th);
		periodicboundaries(GoE,D_s,D_th);
		periodicboundaries(GoS,D_s,D_th);
		periodicboundaries(GoN,D_s,D_th);

		for(i=0;i<D_s*D_th;i++){
			Go[i] = 0.0;
			Go[i] = 0.25*(GoW[i]+GoE[i]+GoS[i]+GoN[i]);
		}periodicboundaries(Go,D_s,D_th);
	}

	if(chipoteR == 1){
		int s_p,th_p;
		s_p = N_s/8+2;
		th_p = N_th/2+5;

		double xp,yp,zp,rp;
		rp = Ro[s_p + D_s*th_p];
		xp = x_arc_0[s_p] - rp*sin(th0[th_p])*sin(al0[s_p]);
		yp = y_arc_0[s_p] + rp*cos(th0[th_p]);
		zp = z_arc_0[s_p] + rp*sin(th0[th_p])*cos(al0[s_p]);

		aneurysmR(x_arc_p,y_arc_p,z_arc_p,alp,th0,D_s,D_th,RoE,xp,yp,zp,rp);
		aneurysmR(x_arc_m,y_arc_m,z_arc_m,alm,th0,D_s,D_th,RoW,xp,yp,zp,rp);
		aneurysmR(x_arc_0,y_arc_0,z_arc_0,al0,thp,D_s,D_th,RoN,xp,yp,zp,rp);
		aneurysmR(x_arc_0,y_arc_0,z_arc_0,al0,thm,D_s,D_th,RoS,xp,yp,zp,rp);

		periodicboundaries(RoE,D_s,D_th);
		periodicboundaries(RoW,D_s,D_th);
		periodicboundaries(RoN,D_s,D_th);
		periodicboundaries(RoS,D_s,D_th);

		for(i=0;i<D_s*D_th;i++){
			Ro[i] = 0.0;
			Ro[i] = 0.25*(RoW[i]+RoE[i]+RoS[i]+RoN[i]);
		}periodicboundaries(Ro,D_s,D_th);
	}

	if(seccionR == 1){
		oclussionR(D_s,D_th,RoE,sp,Ls);
		oclussionR(D_s,D_th,RoW,sm,Ls);
		oclussionR(D_s,D_th,RoS,s0,Ls);
		oclussionR(D_s,D_th,RoN,s0,Ls);

/*
		bulgingR(D_s,D_th,RoE,sp,Ls);
		bulgingR(D_s,D_th,RoW,sm,Ls);
		bulgingR(D_s,D_th,RoS,s0,Ls);
		bulgingR(D_s,D_th,RoN,s0,Ls);
*/
		periodicboundaries(RoE,D_s,D_th);
		periodicboundaries(RoW,D_s,D_th);
		periodicboundaries(RoN,D_s,D_th);
		periodicboundaries(RoS,D_s,D_th);

		for(i=0;i<D_s*D_th;i++){
			Ro[i] = 0.0;
			Ro[i] = 0.25*(RoW[i]+RoE[i]+RoS[i]+RoN[i]);
		}periodicboundaries(Ro,D_s,D_th);		
	}

	if(seccionG == 1){
		bulgingG(D_s,D_th,GoE,sp,Ls);
		bulgingG(D_s,D_th,GoW,sm,Ls);
		bulgingG(D_s,D_th,GoS,s0,Ls);
		bulgingG(D_s,D_th,GoN,s0,Ls);

		periodicboundaries(GoE,D_s,D_th);
		periodicboundaries(GoW,D_s,D_th);
		periodicboundaries(GoN,D_s,D_th);
		periodicboundaries(GoS,D_s,D_th);

		for(i=0;i<D_s*D_th;i++){
			Go[i] = 0.0;
			Go[i] = 0.25*(GoW[i]+GoE[i]+GoS[i]+GoN[i]);
		}periodicboundaries(Go,D_s,D_th);
	}


	double *Ao,*AoW,*AoE,*AoS,*AoN;
	Crea1(&Ao,D_s*D_th);
	Crea1(&AoW,D_s*D_th);
	Crea1(&AoE,D_s*D_th);
	Crea1(&AoS,D_s*D_th);
	Crea1(&AoN,D_s*D_th);

	double Gamma;
	initial_cond(th0,da0,D_s,D_th,Ro,Ao);
	initial_cond(th0,dap,D_s,D_th,RoE,AoE);
	initial_cond(th0,dam,D_s,D_th,RoW,AoW);
	initial_cond(thp,da0,D_s,D_th,RoN,AoN);
	initial_cond(thm,da0,D_s,D_th,RoS,AoS);
	print_RoGo(Go,Ro,D_s,D_th,file);
	//return 0;

	//--------------------------------------------------------------------------
	// creating the W variables RKFLUX W1 = A, W2 = A*U_s, W3 = A*U_th
	//--------------------------------------------------------------------------
	double *W01,*W02,*W03,*W11,*W12,*W13;
	Crea1(&W01,D_s*D_th);
	Crea1(&W02,D_s*D_th);
	Crea1(&W03,D_s*D_th);
	Crea1(&W11,D_s*D_th);
	Crea1(&W12,D_s*D_th);
	Crea1(&W13,D_s*D_th);
	
	for(i=0;i<D_s*D_th;i++){
		W01[i] = 0.0; // A
		W02[i] = 0.0; // psi_s_0*A*u
		W03[i] = 0.0; // A*L = psi_th_0*A*A*omega
		W11[i] = 0.0;
		W12[i] = 0.0;
		W13[i] = 0.0;
	}

	//--------------------------------------------------------------------------
	// Creating the post-processing variables
	//--------------------------------------------------------------------------
	double r_ic;
	int dim_post;
	dim_post = D_s*(N_th+1);

	double *R,*A,*u,*w,*L,*Ut,*P,*PT,*Qs,*Qt,*VT;
	Crea1(&R,dim_post);
	Crea1(&A,dim_post);
	Crea1(&u,dim_post);
	Crea1(&w,dim_post);
	Crea1(&L,dim_post);
	Crea1(&Ut,dim_post);
	Crea1(&P,dim_post);
	Crea1(&PT,dim_post);
	Crea1(&Qs,dim_post);
	Crea1(&Qt,dim_post);
	Crea1(&VT,dim_post);

	double *Vx,*Vy,*Vz;
	Crea1(&Vx,dim_post);
	Crea1(&Vy,dim_post);
	Crea1(&Vz,dim_post);

	double *Vex,*Vey,*Vez;
	Crea1(&Vex,dim_post);
	Crea1(&Vey,dim_post);
	Crea1(&Vez,dim_post);
	
	//--------------------------------------------------------------------------
	// Read data to create cardiac cycle
	//--------------------------------------------------------------------------
	FILE *entrada = NULL;
	entrada = fopen("cardiac_cycle_data/datos01.dat","r");
	if(entrada == NULL){
		puts("No se pudo abrir el archivo.");
		exit(0);
	}
	
	int N_cycle;
	fscanf(entrada,"%d",&N_cycle);
	double Time_CC[N_cycle];
	double Val_CC[N_cycle];

	for(k=0;k<N_cycle;k++){
		fscanf(entrada,"%lf",&Time_CC[k]);
		fscanf(entrada,"%lf",&Val_CC[k]);
	}fclose(entrada);

	double T_cycle;
	T_cycle = Time_CC[N_cycle-1];
	int Nf = 15;
	_Complex coeff[Nf+1];
	for(k=0;k<Nf+1;k++) coeff[k] = coefficients(N_cycle,T_cycle,k,Val_CC);
	
	//--------------------------------------------------------------------------
	// Initial condition for A
	//--------------------------------------------------------------------------
	double z_sum[D_s];
	sum_z(z_sum,al0,D_s,ds);
	if(pert == 0){
		for(i=0;i<D_s*D_th;i++) W01[i] = 1.44*Ao[i];
		// initial condition for Q1 and Q2
		double PAs;
		PAs = IVs(g_s,1.0);
		double PA0[3];
		PA0[0] = IVs(g_s,1.0);
		PA0[1] = 2.0*IVs(g_s,2.0);
		PA0[2] = IVs(g_s,3.0);
		double A_ic,R_ic,G_ic,p_ic;
		printf(" estoy por crear la condicion inicial\n");

		for(k=0;k<N_th+1;k++){
			for(j=0;j<D_s;j++){
				A_ic = W01[j+(k+2)*D_s];
				R_ic = findingroot(A_ic,sin(th0[k+2]),da0[j]);
				G_ic = R_ic*sin(th0[k+2])*da0[j];
				p_ic = A0(R_ic,G_ic,PA0)/As(R_ic,G_ic,PAs);
				W02[j+(k+2)*D_s] = p_ic*A_ic*0.0;
				W03[j+(k+2)*D_s] = A_ic*0.0;
			}printf("k = %d\n",k);
		}printf(" ya cree la condicion inicial\n");
	}else perturbation(x_arc_0,y_arc_0,z_arc_0,al0,th0,da0,D_s,D_th,Ro,W01);

	//--------------------------------------------------------------------------
	// Periodic condition to A
	//--------------------------------------------------------------------------
	periodicboundaries(W01,D_s,D_th);
	periodicboundaries(W02,D_s,D_th);
	periodicboundaries(W03,D_s,D_th);

	//--------------------------------------------------------------------------
	// creating the C of the RKFLUX
	//--------------------------------------------------------------------------
	double *CW01,*CW02,*CW03,*CW11,*CW12,*CW13;
	Crea1(&CW01,D_s*D_th);
	Crea1(&CW02,D_s*D_th);
	Crea1(&CW03,D_s*D_th);
	Crea1(&CW11,D_s*D_th);
	Crea1(&CW12,D_s*D_th);
	Crea1(&CW13,D_s*D_th);
	for(j=0;j<D_s*D_th;j++){
		CW01[j] = 0.0;
		CW02[j] = 0.0;
		CW03[j] = 0.0;
		CW11[j] = 0.0;
		CW12[j] = 0.0;
		CW13[j] = 0.0;
	}

	//--------------------------------------------------------------------------
	// printing the initial condition
	//--------------------------------------------------------------------------
	create_post_processing(D_s,N_th+1,W01,W02,W03,th0,da0,x_arc_0,y_arc_0,
	z_arc_0,al0,R,A,u,w,L,Ut,P,PT,Qs,Qt,VT,Vx,Vy,Vz,Vex,Vey,Vez,dth,Go,Ro,Ao,
	beta,g,g_s,g_th,z_sum);
	print_IC(R,A,u,w,L,Ut,P,PT,Qs,Qt,VT,Vx,Vy,Vz,Vex,Vey,Vez,D_s,N_th+1,file);
	printf("\n");		
	min_R = min_vec(R,D_s*(N_th+1));

	//--------------------------------------------------------------------------
	// First Post-Processing Data
	//--------------------------------------------------------------------------

	if(bt<10) printf("      bt = %d",bt);
	if(bt>9 && bt<100) printf("     bt = %d",bt);
	if(bt>99 &&  bt<1000) printf("    bt = %d",bt);
	if(bt>999 && bt<10000) printf("   bt = %d",bt);	
	if(bt>9999) printf("  bt = %d",bt);	
	printf("   dt = %e",dt);
	printf("    T = %e",T);
	printf("    Tb = %e",1.0*savedata*Tfinal/movies);
	printf("     SD = %d",savedata);
	print_vec_vel(Vex,Vey,Vez,D_s,N_th+1,file);
	print_vec_pos(Vx,Vy,Vz,D_s,N_th+1,file);
	print_3D(R,A,u,w,L,Ut,P,PT,Qs,Qt,VT,D_s,N_th+1,file);
	savedata = savedata + 1;
	printf("\n");

	double max_C1,max_C2,max_C3;
	int cont_cont;

	do{
		cont_cont = 0;
		wtime = omp_get_wtime();

		//----------------------------------------------------------------------
		// Ghost Velocity Routine
		//----------------------------------------------------------------------
		u_gh = ghost(coeff,Nf,T,T_cycle);
		v_inflow = parVin*u_gh;

		//----------------------------------------------------------------------
		// First Boundary Conditions
		//----------------------------------------------------------------------

		// Periodic Boundaries
		periodicboundaries(W01,D_s,D_th);
		periodicboundaries(W02,D_s,D_th);
		periodicboundaries(W03,D_s,D_th);

		// Left Boundaries
		leftboundaries(W01,W02,W03,W1_gh,W2_gh,W3_gh,Ro,Go,th0,da0,D_th,
		D_s,beta,g_s,g_th,v_inflow,left_cont,Ao);

		// Right Boundaries
		rightboundaries(W01,W02,W03,W1_gh,W2_gh,W3_gh,Ro,Go,th0,da0,D_th,
		D_s,beta,g_s,g_th,right_cont,Ao);
		
		Imprimir_Boundaries(left_cont,right_cont,D_th,file);
		//----------------------------------------------------------------------
		// First Step
		//----------------------------------------------------------------------
		par[0] = 0.0;
		par[1] = 0.0;
		par[2] = min_R;
		NumFluxG(CW01,CW02,CW03,W01,W02,W03,D_s,D_th,Ro,RoE,RoW,RoN,RoS,Go,GoE,
		GoW,GoN,GoS,par,thp,th0,thm,alp,al0,alm,dap,da0,dam,ddp,dd0,ddm,Ao,AoE,
		AoW,AoN,AoS);
		//return 0;
		max_a_0 = par[0];
		max_b_0 = par[1];
		min_R_0 = par[2];
		dt_0 = CFL/(max_a_0/ds + max_b_0/dth);
		#pragma omp parallel for
		for(ic=0;ic<D_s*D_th;ic++){
			W11[ic] = W01[ic] + dt_0*CW01[ic];
			W12[ic] = W02[ic] + dt_0*CW02[ic];
			W13[ic] = W03[ic] + dt_0*CW03[ic];
		}

		//----------------------------------------------------------------------
		// Second Boundary Conditions
		//----------------------------------------------------------------------
		// Periodic Boundaries
		periodicboundaries(W11,D_s,D_th);
		periodicboundaries(W12,D_s,D_th);
		periodicboundaries(W13,D_s,D_th);

		// Left Boundaries
		leftboundaries(W11,W12,W13,W1_gh,W2_gh,W3_gh,Ro,Go,th0,da0,D_th,D_s,
		beta,g_s,g_th,v_inflow,left_cont,Ao);

		// Right Boundaries
		rightboundaries(W11,W12,W13,W1_gh,W2_gh,W3_gh,Ro,Go,th0,da0,D_th,D_s,
		beta,g_s,g_th,right_cont,Ao);
		
		Imprimir_Boundaries(left_cont,right_cont,D_th,file);

		//----------------------------------------------------------------------
		// Second Step
		//----------------------------------------------------------------------
		par[0] = 0.0;
		par[1] = 0.0;
		par[2] = min_R_0;
		//return 0;
		NumFluxG(CW11,CW12,CW13,W11,W12,W13,D_s,D_th,Ro,RoE,RoW,RoN,RoS,Go,GoE,
		GoW,GoN,GoS,par,thp,th0,thm,alp,al0,alm,dap,da0,dam,ddp,dd0,ddm,Ao,AoE,
		AoW,AoN,AoS);
		max_a_1 = par[0];
		max_b_1 = par[1];
		min_R_1 = par[2];
		dt_1 = CFL/(max_a_1/ds + max_b_1/dth);
		#pragma omp parallel for
		for(ic=0;ic<D_s*D_th;ic++){
			W01[ic] = 0.5*(W01[ic] + W11[ic] + dt_1*CW11[ic]);
			W02[ic] = 0.5*(W02[ic] + W12[ic] + dt_1*CW12[ic]);
			W03[ic] = 0.5*(W03[ic] + W13[ic] + dt_1*CW13[ic]);
		}

		//----------------------------------------------------------------------
		// Final Variables
		//----------------------------------------------------------------------
		dt = 0.5*(dt_0 + dt_1);
		T = T + dt;
		wtime = omp_get_wtime()-wtime;
		bt = bt + 1;
		max_a = 0.5*(max_a_0 + max_a_1);
		max_b = 0.5*(max_b_0 + max_b_1);
		min_R = min_R_1;

		ImprimirValores(dt,max_a,max_b,min_R,u_gh,wtime,T,file);
		sum_wtime = sum_wtime + wtime;
		max_wtime = max(wtime,max_wtime);
		min_wtime = min(wtime,min_wtime);

		//----------------------------------------------------------------------
		// Post-Processing Data
		//----------------------------------------------------------------------

		if(T>=(1.0*savedata*Tfinal/movies)){
			if(bt<10) printf("      bt = %d",bt);
			if(bt>9 && bt<100) printf("     bt = %d",bt);
			if(bt>99 &&  bt<1000) printf("    bt = %d",bt);
			if(bt>999 && bt<10000) printf("   bt = %d",bt);	
			if(bt>9999) printf("  bt = %d",bt);	
			printf("   dt = %e",dt);
			printf("   R = %e",min_R/min_R_i);
			printf("   T = %e",T);
			printf("   Te = %e",wtime);

			printf("   Tb = %e",1.0*savedata*Tfinal/movies);
			printf("   SD = %d",savedata);
			create_post_processing(D_s,N_th+1,W01,W02,W03,th0,da0,x_arc_0,
			y_arc_0,z_arc_0,al0,R,A,u,w,L,Ut,P,PT,Qs,Qt,VT,Vx,Vy,Vz,Vex,Vey,Vez,
			dth,Go,Ro,Ao,beta,g,g_s,g_th,z_sum);
			print_vec_vel(Vex,Vey,Vez,D_s,N_th+1,file);
			print_vec_pos(Vx,Vy,Vz,D_s,N_th+1,file);
			print_3D(R,A,u,w,L,Ut,P,PT,Qs,Qt,VT,D_s,N_th+1,file);
			print_inflow(left_cont,D_th,file);
			print_outflow(right_cont,D_th,file);
			savedata = savedata + 1;
			printf("\n");
		}
	}while(T<Tfinal && dt>pow(10.0,-15.0) && min_R > 0.0);
	printf("\n");
	printf("The max loop-time was %.15e\n",max_wtime);
	printf("The min loop-time was %.15e\n",min_wtime);
	printf("The mean loop-time was %.15e\n",sum_wtime/bt);
	printf("\n");
	
	printf("T=%.15e, Tfinal=%.15e\n",T,Tfinal);
	printf("dt=%.15e, minR=%.15e\n",dt,min_R);
	print_IF(R,A,u,w,L,Ut,P,PT,Qs,Qt,VT,Vx,Vy,Vz,Vex,Vey,Vez,D_s,N_th+1,file);
	printf("\n");		
	
	//--------------------------------------------------------------------------
	// free memory
	//--------------------------------------------------------------------------
	free(R);
	free(A);
	free(u);
	free(w);
	free(L);
	free(P);
	free(PT);
	free(Qs);
	free(Qt);
	free(Ut);
	free(VT);
	
	free(Vx);
	free(Vy);
	free(Vz);

	free(Vex);
	free(Vey);
	free(Vez);

	free(Go);
	free(GoE);
	free(GoW);
	free(GoN);
	free(GoS);

	free(Ro);
	free(RoE);
	free(RoW);
	free(RoN);
	free(RoS);

	free(Ao);
	free(AoE);
	free(AoW);
	free(AoN);
	free(AoS);

	free(W01);
	free(W02);
	free(W03);

	free(CW01);
	free(CW02);
	free(CW03);

	free(W11);
	free(W12);
	free(W13);

	free(CW11);
	free(CW12);
	free(CW13);
}

void initial_cond(double th[],double da[],int D_s,int D_th,double *Ro,double *Ao)
{
	int i,j,k;
	for(k=0;k<D_th;k++){
		for(j=0;j<D_s;j++){
			i = j + D_s*k;
			Ao[i] = Q0(Ro[i],Ro[i]*sin(th[k])*da[j]);
		}
	}
}

void initial_cond2(double Po,double z_sum[],int D_s,int D_th,double *Go,double *Ao,double g)
{
	int i,j,k;
	double zo;
	for(j=0;j<D_s;j++){
		zo = z_sum[j];
		for(k=0;k<D_th;k++){
			i = j + D_s*k;
			Ao[i] = Ao[i]*(1.0 + (Po - g*zo)/Go[i]);
		}
	}
}

void sum_z(double z_sum[],double al[],int D_s,double ds)
{
	z_sum[0] = -0.5*(sin(al[0]) + 2.0*sin(al[1]) + sin(al[2]))*ds;
	z_sum[1] = -0.5*(sin(al[0]) + sin(al[1]))*ds;
	z_sum[2] = 0.0;
	int j;
	double sum = 0.0;
	for(j=3;j<D_s;j++){
		sum = sum + 0.5*(sin(al[j-1]) + sin(al[j]))*ds;
		z_sum[j] = sum;
	}
}

void arc(double al[],double s[],double barrier,int D_s)
{
	int j;
	double stop;
	stop = 1.0;
	for(j=0;j<D_s;j++){
		if(s[j]<stop*barrier) al[j] = 0.5*M_PI*(1.0 - 2.0*s[j]/barrier);
		else al[j] = 0.5*(1.0-2.0*stop)*M_PI;
		//if(s[j]<0.0) al[j] = 0.5*M_PI;
	}
}

void heatequation(double vec[],int D_s)
{
	int j;
	double vecn[D_s];
	for(j=0;j<5;j++) vecn[j] = vec[j];
	for(j=5;j<D_s-5;j++) vecn[j] = 0.25*(vec[j-1] + 2.0*vec[j] + vec[j+1]);
	for(j=D_s-5;j<D_s;j++) vecn[j] = vec[j];
	for(j=0;j<D_s;j++) vec[j] = vecn[j];
}

void diffarc1(double dvec[],double vec[],int D_s,double ds)
{
	int j;
	for(j=1;j<D_s-5;j++) dvec[j] = 0.5*(vec[j+1] - vec[j-1])/ds;
	dvec[0] = dvec[1];
}

void diffarc2(double dvec[],double vec[],int D_s,double ds)
{
	int j;
	for(j=1;j<D_s-5;j++) dvec[j] = (vec[j+1] - 2.0*vec[j] + vec[j-1])/(ds*ds);
	dvec[0] = dvec[1];
}

void printfarc(double vec[],double dvec[],double ddvec[],int dim,int carpet)
{
	int j;
	char eps_name[100];
	FILE *salida = NULL;
	sprintf(eps_name,"Example/Ex%d/Parameters/al.dat",carpet);
	salida = fopen(eps_name,"a");
	if(salida==NULL){
		puts("No se pudo abrir el archivo sobre dt.");
		exit(0);
	}for(j=0;j<dim;j++) fprintf(salida,"%.15e\n",vec[j]);
	fclose(salida);
	salida = NULL;

	FILE *salida1 = NULL;
	sprintf(eps_name,"Example/Ex%d/Parameters/da.dat",carpet);
	salida1 = fopen(eps_name,"a");
	if(salida1==NULL){
		puts("No se pudo abrir el archivo sobre dt.");
		exit(0);
	}for(j=0;j<dim;j++) fprintf(salida1,"%.15e\n",dvec[j]);
	fclose(salida1);
	salida1 = NULL;

	FILE *salida2 = NULL;
	sprintf(eps_name,"Example/Ex%d/Parameters/dd.dat",carpet);
	salida2 = fopen(eps_name,"a");
	if(salida2==NULL){
		puts("No se pudo abrir el archivo sobre dt.");
		exit(0);
	}for(j=0;j<dim;j++) fprintf(salida2,"%.15e\n",ddvec[j]);
	fclose(salida2);
	salida2 = NULL;
}

void printfrc(double vecr[],double vecc[],int dim,int carpet)
{
	int j;
	char eps_name[100];
	FILE *salida = NULL;
	sprintf(eps_name,"Example/Ex%d/Parameters/Radius.dat",carpet);
	salida = fopen(eps_name,"a");
	if(salida==NULL){
		puts("No se pudo abrir el archivo sobre Radius.");
		exit(0);
	}for(j=0;j<dim;j++) fprintf(salida,"%.15e\n",vecr[j]);
	fclose(salida);
	salida = NULL;

	FILE *salida1 = NULL;
	sprintf(eps_name,"Example/Ex%d/Parameters/Velocity.dat",carpet);
	salida1 = fopen(eps_name,"a");
	if(salida1==NULL){
		puts("No se pudo abrir el archivo sobre Velocity.");
		exit(0);
	}for(j=0;j<dim;j++) fprintf(salida1,"%.15e\n",vecc[j]);
	fclose(salida1);
	salida1 = NULL;
}


void print_inflow(int vec[],int D,int carpet)
{
	int j;
	char eps_name[100];
	FILE *salida = NULL;
	sprintf(eps_name,"Example/Ex%d/Parameters/Inflow.dat",carpet);
	salida = fopen(eps_name,"a");
	if(salida==NULL){
		puts("No se pudo abrir el archivo sobre inflow.");
		exit(0);
	}for(j=0;j<D;j++) fprintf(salida,"%d\n",vec[j]);
	fclose(salida);
	salida = NULL;
}

void Imprimir_Boundaries(int vL[],int vR[],int dim,int carpet)
{
	int j;
	char eps_name[100];
	FILE *salida = NULL;
	sprintf(eps_name,"Example/Ex%d/Parameters/Boundaries.dat",carpet);
	salida = fopen(eps_name,"a");
	if(salida==NULL){
		puts("No se pudo abrir el archivo sobre boundaries.");
		exit(0);
	}int sumL=0,sumR=0;
	for(j=0;j<dim;j++){
		sumL = sumL + vL[j];
		sumR = sumR + vR[j];
	}fprintf(salida,"(%d,%d)\n",sumL,sumR);
	fclose(salida);
	salida = NULL;
}

void print_outflow(int vec[],int D,int carpet)
{
	int j;
	char eps_name[100];
	FILE *salida = NULL;
	sprintf(eps_name,"Example/Ex%d/Parameters/Outflow.dat",carpet);
	salida = fopen(eps_name,"a");
	if(salida==NULL){
		puts("No se pudo abrir el archivo sobre outflow.");
		exit(0);
	}for(j=0;j<D;j++) fprintf(salida,"%d\n",vec[j]);
	fclose(salida);
	salida = NULL;
}

void print_RoGo(double *Go,double *Ro,int dim_s,int dim_th,int carpet)
{
	int j,k;
	char eps_name[100];

	FILE *Go_F = NULL;
	sprintf(eps_name,"Example/Ex%d/IC/IC_Go.dat",carpet);
	Go_F = fopen(eps_name,"a");
	if(Go_F==NULL){
		puts("No se puede abrir velocity_y!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=2;k<dim_th-1;k++){
		for(j=0;j<dim_s;j++){
			fprintf(Go_F,"%.15e\n",Go[j+dim_s*k]);
		}
	}//printf(" Save Vy");
	fclose(Go_F);
	Go_F = NULL;

	FILE *Ro_F = NULL;
	sprintf(eps_name,"Example/Ex%d/IC/IC_Ro.dat",carpet);
	Ro_F = fopen(eps_name,"a");
	if(Ro_F==NULL){
		puts("No se puede abrir velocity_y!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=2;k<dim_th-1;k++){
		for(j=0;j<dim_s;j++){
			fprintf(Ro_F,"%.15e\n",Ro[j+dim_s*k]);
		}
	}//printf(" Save Vy");
	fclose(Ro_F);
	Ro_F = NULL;
}
