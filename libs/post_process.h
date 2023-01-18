// It is not recommended to put function definitions  
// in a header file. Ideally there should be only 
// function declarations. Purpose of this code is 
// to only demonstrate working of header files. 
#ifndef _POST_PROCESSH_
#define _POST_PROCESSH_

#include "flux_extra.h"

void create_post_processing(int ds,int dth,double *W1,double *W2,double *W3,
double th[],double da[],double x[],double y[],double z[],double a[],double *R,
double *A,double *u,double *w,double *L,double *Ut,double *P,double *PT,
double *Qs,double *Qt,double *VT,double *Vx,double *Vy,double *Vz,double *Vex,
double *Vey,double *Vez,double delta_th,double *Go,double *Ro,double *Ao,
double beta,double g,double gs,double gt,double z_sum[])
{
	int j,k;
	int ic,is;
	double r_ic;
	double rc,rp,rm;
	double Rp,Rpp,Cr,Per;
	double xp,yp,zp;
	double Qr,Ga;
	double PAs,PAt[2],PA0[3],PB0[2];
	PAs = IVs(gs,1.0);
	PA0[0] = IVs(gs,1.0);
	PA0[1] = 2.0*IVs(gs,2.0);
	PA0[2] = IVs(gs,3.0);
	PAt[0] = IVt(gt,1.0); 
	PAt[1] = IVt(gt,2.0);
	PB0[0] = IVt(gt,3.0);
	PB0[1] = IVt(gt,4.0);
	double cs,VAs,ct,VAt;

	for(k=0;k<dth;k++){
		for(j=0;j<ds;j++){
			ic = j + k*ds;
			is = j + (k+2)*ds;
			 A[ic] = W1[is];
			Qs[ic] = W2[is];
			Qt[ic] = W3[is];
			 P[ic] = Go[is]*(pow(W1[is]/Ao[is],beta/2.0) - 1.0);

			rp = findingroot(W1[is],sin(th[k+3]),da[j]);
			rc = findingroot(W1[is],sin(th[k+2]),da[j]);
			rm = findingroot(W1[is],sin(th[k+1]),da[j]);

			 R[ic] = rc;
			PT[ic] = P[ic] + g*z_sum[j];

			Qr = rc/Ro[is];
			Ga = rc*sin(th[k+2])*da[j];
			cs = W1[is]/As(rc,Ga,PAs);
			ct = W1[is]/At(rc,Ga,PAt);
			VAs = A0(rc,Ga,PA0);
			VAt = B0(rc,Ga,PB0);

			u[ic] = W2[is]/(cs*VAs);
			w[ic] = W3[is]/(ct*VAt);
			L[ic] = W3[is]/W1[is];

			Ut[ic] = (ct/A[ic])*(PAt[1] - PB0[0]*Ga)*pow(rc,3.0)*w[ic]; 

			Vx[ic] = x[j] - rc*sin(th[k+2])*sin(a[j]);
			Vy[ic] = y[j] + rc*cos(th[k+2]);
			Vz[ic] = z[j] + rc*sin(th[k+2])*cos(a[j]);

			Rp = 0.5*(rp-rm)/delta_th;
			xp = (sin(th[k+2])*Rp/rc + cos(th[k+2]))*(-sin(a[j]));
			yp = -sin(th[k+2]) + cos(th[k+2])*Rp/rc;
			zp = (sin(th[k+2])*Rp/rc + cos(th[k+2]))*(cos(a[j]));

			Rpp = (rp - 2.0*rc + rm)/pow(delta_th,2.0);
			Per = sqrt(pow(rc,2.0)+pow(Rp,2.0));
			Cr = pow(Per,3.0)/(pow(rc,2.0)+2.0*pow(Rp,2.0)-rc*Rpp);

			Vex[ic] = (Cr*xp/Per)*Ut[ic] + cos(a[j])*u[ic];
			Vey[ic] = (Cr*yp/Per)*Ut[ic];
			Vez[ic] = (Cr*zp/Per)*Ut[ic] + sin(a[j])*u[ic];
			 VT[ic] = sqrt(pow(u[ic],2.0) + pow((Cr/rc)*Ut[ic],2.0));
// Cr es el radio de curvatura, R_c en el articulo.
// rc es el radio en la posicion que queremos calcular la velocidad.
		}	
	}
}

#endif
