// It is not recommended to put function definitions  
// in a header file. Ideally there should be only 
// function declarations. Purpose of this code is 
// to only demonstrate working of header files. 

#ifndef _FLUX_EXTRAH_
#define _FLUX_EXTRAH_

#include "minmod.h"
#include "find_root.h"
#include "profiles.h"

double ga(double A,double R0,double sint,double alfaps)
{
	return findingroot(A,sint,alfaps)/R0;
}

void Reconstruction(double diff,double V0,double V1,double V2,double V3,
double V4,double vmm,double *Vb1,double *Vb2,double *Vb3,double *Vb4)
{
	*Vb1 = V1 + 0.5*diff*minmod(V2,V1,V0,diff,vmm);
	*Vb2 = V2 - 0.5*diff*minmod(V3,V2,V1,diff,vmm);
	*Vb3 = V2 + 0.5*diff*minmod(V3,V2,V1,diff,vmm);
	*Vb4 = V3 - 0.5*diff*minmod(V4,V3,V2,diff,vmm);
}

double ApbarGo(double A,double Ao,double Go,double beta,double Aphat)
{
	
	return (A*Go*(pow(A/Ao,beta/2.0) - 1.0) - Aphat)/Go;
}

double ApbarGo2(double A,double Aphat,double p,double Go)
{
	
	return (A*p-Aphat)/Go;
}

double Aphat(double beta,double Go,double Ao,double A)
{
	double pressure;
	pressure = Go*(pow(A/Ao,beta/2.0)-1.0);
	return beta*pressure*A/(beta+2.0) - beta*Go*(Ao-A);
	//return beta*Go*Ao*(pow(A/Ao,beta/2.0 + 1.0) - 1.0)/(beta+2.0);
}

double Aphat2(double beta,double Go,double Ao,double A,double p)
{
	return beta*A*p/(beta+2.0) - beta*Go*(Ao-A)/(beta+2.0);
}

double Pressure(double beta,double Go,double Ao,double A)
{
	return Go*(pow(A/Ao,beta/2.0) - 1.0);
}

double dAphat(double beta,double Go,double A,double Ao)
{
	return (beta/2.0)*Go*pow(A/Ao,beta/2.0);
}

void Variables(double *V,double Qr,double Qs,double Qt,double Ro,double Ao,
double Go,double s,double da,double PAs,double PAt[],double PA0[],double PB0[],
double beta)
{
	// Radius
	V[0] = Qr*Ro;				// R
	V[1] = V[0]*sin(s)*da;		// Ga
	V[2] = Q0(V[0],V[1]);		// A

	// Axial
	V[3] = As(V[0],V[1],PAs);	// As(R,Ga,gs)
	V[4] = V[2]/V[3];			// c_s = A/As
	V[5] = A0(V[0],V[1],PA0);	// A0(R,Ga,gs)
	V[6] = Qs/(V[4]*V[5]);		// u

	// Angular
	V[7] = At(V[0],V[1],PAt);	// At(R,Ga,gt)
	V[8] = V[2]/V[7];			// c_t = A/At
	V[9] = B0(V[0],V[1],PB0);	// B0(R,Ga,gt)
	V[10] = Qt/(V[8]*V[9]);		// omega

	// Jacobian for derivatives
	V[11] = Jac(V[0],V[1]);

	// Pressure
	V[13] = dAphat(beta,Go,V[2],Ao);

	// Pressure
	V[14] = Pressure(beta,Go,Ao,V[2]);

	// Conservative Pressure
	V[12] = Aphat2(beta,Go,Ao,V[2],V[14]);
}

void lambda(double *lam,double *p,double pre,double h)
{
	lam[0] = p[0]*h;
	lam[1] = p[1]*h - sqrt(p[3]*pre + p[2]*pow(h,2.0));
	lam[2] = p[1]*h + sqrt(p[3]*pre + p[2]*pow(h,2.0));
}

double Disc(double a,double da,double b,double db,double c,double dc,double J,double A)
{
	return pow(A,2.0)*(a*J*(b*dc - 2.0*c*db) + 0.25*pow(b*da - a*db - 2.0*c*J,2.0))/pow(J*a*b,2.0);
}

double Extr(double a,double da,double b,double db,double c,double J,double A)
{
	return 0.5*A*(b*da - a*db + 2.0*c*J)/(J*a*b);
}

void fluxF(double *vec,double *A,double u,double o,double p,double cs,double ct)
{
	vec[0] = A[0]*u;
	vec[1] = cs*cs*A[1]*pow(u,2.0) + p;
	vec[2] = cs*ct*A[2]*u*o;
}

void create_lam_s(double p[],double A[],double D[],double V[])
{
	p[0] = V[4]*A[2]/V[9];
	// c_s*B1/B0
	p[1] = Extr(V[3],D[0],V[5],D[1],A[1],V[11],V[2]);
	p[2] = Disc(V[3],D[0],V[5],D[1],A[1],D[2],V[11],V[2]);
	p[3] = V[3]/V[5];
	// As/A0
}

void Axial_Req(double *F,double *maxL,double *minL,double *V,double PA1[],double PB1,double PAs,double PA0[])
{
	double A[3];
	A[0] = V[2];
	A[1] = A1(V[0],V[1],PA1);
	A[2] = B1(V[0],V[1],PB1);
	fluxF(F,A,V[6],V[10],V[12],V[4],V[8]);

	double D[3];
	D[0] = dAs(V[0],V[1],PAs);
	D[1] = dA0(V[0],V[1],PA0);
	D[2] = dA1(V[0],V[1],PA1);

	double p[4],lam[3];
	create_lam_s(p,A,D,V);
	lambda(lam,p,V[13],V[6]);
	*maxL = maximum(lam[0],lam[1],lam[2]);
	*minL = minimum(lam[0],lam[1],lam[2]);
}

void create_lam_th(double p[],double B[],double D[],double V[])
{
	p[0] = V[8]*B[1]/V[5];
	// c_t*A2/A0
	p[1] = Extr(V[7],D[0],V[9],D[1],B[2],V[11],V[2]);
	p[2] = Disc(V[7],D[0],V[9],D[1],B[2],D[2],V[11],V[2]);
	p[3] = V[7]/V[9];
	// At/B0	
}

void fluxG(double *vec,double *A,double u,double o,double p,double cs,double ct)
{
	vec[0] = A[0]*o;
	vec[1] = cs*ct*A[1]*u*o;
	vec[2] = ct*ct*A[2]*pow(o,2.0) + p;
}

void Angular_Req(double *F,double *maxL,double *minL,double *V,double PA2[],double PB2[],double PAt[],double PB0[])
{
	double B[3];
	B[0] = V[2];
	B[1] = A2(V[0],V[1],PA2);
	B[2] = B2(V[0],V[1],PB2);
	fluxG(F,B,V[6],V[10],V[12],V[4],V[8]);

	double D[3];
	D[0] = dAt(V[0],V[1],PAt);
	D[1] = dB0(V[0],V[1],PB0);
	D[2] = dB2(V[0],V[1],PB2);

	double p[4],lam[3];
	create_lam_th(p,B,D,V);
	lambda(lam,p,V[13],V[10]);
	*maxL = maximum(lam[0],lam[1],lam[2]);
	*minL = minimum(lam[0],lam[1],lam[2]);
}
 			
#endif
