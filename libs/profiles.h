// It is not recommended to put function definitions  
// in a header file. Ideally there should be only 
// function declarations. Purpose of this code is 
// to only demonstrate working of header files. 

#ifndef _PROFILESH_
#define _PROFILESH_

double IA(double n)
{
	return 1.0/n;
}

double IVs(double g,double n)
{
	return g*IA(n+1.0)*IA(g+n+1.0);
}

double IVt(double g,double n)
{
	return (2.0*g+n+2.0)*IA(g+1.0)*IA(g+n)*IA(g+n+1.0);
}

double IVsVs(double g,double n)
{
	return 2.0*g*g*IA(n+1.0)*IA(g+n+1.0)*IA(2.0*g+n+1.0);
}

double IVtVt(double g,double n)
{
	return (n*(n+1.0) + 2.0*g*(5.0*g+3.0*n+2.0))*pow(IA(g+1.0),2.0)*IA(2.0*g+n-1.0)*IA(2.0*g+n)+IA(2.0*g+n+1.0);
}

double IVsVt(double gs,double gt,double n)
{
	return gs*(pow(n+1.0,2.0) + gt*(3.0*gt+n+3.0) + gs*(2.0*gt+n+1.0))*IA(gt+1.0)*IA(gt+n)*IA(gt+n+1.0)*IA(gt+gs+n)*IA(gt+gs+n+1.0);
}

double Jac(double r,double ga)
{
	return r*(1.0-ga);
}

double Q0(double r,double ga)
{
	return (IA(2.0) - IA(3.0)*ga)*r*r;
}

double As(double r,double ga,double p)
{
	return p*r*r;
}

double dAs(double r,double ga,double p)
{
	return 2.0*p*r;
}

double At(double r,double ga,double p[])
{
	return (p[0] - p[1]*ga)*r*r;
}

double dAt(double r,double ga,double p[])
{
	return (2.0*p[0] - 3.0*p[1]*ga)*r;
}

double A0(double r,double ga,double p[])
{
	return (p[0] - p[1]*ga + p[2]*ga*ga)*r*r;
}

double dA0(double r,double ga,double p[])
{
	return (2.0*p[0] - 3.0*p[1]*ga + 4.0*p[2]*ga*ga)*r;
}

double A1(double r,double ga,double p[])
{
	return (p[0] - p[1]*ga)*r*r;
}

double dA1(double r,double ga,double p[])
{
	return (2.0*p[0] - 3.0*p[1]*ga)*r;
}

double A2(double r,double ga,double p[])
{
	return (p[0] - p[1]*ga + p[2]*ga*ga)*r*r;
}

double B0(double r,double ga,double p[])
{
	return (p[0] - p[1]*ga)*r*r*r*r;
}

double dB0(double r,double ga,double p[])
{
	return (4.0*p[0] - 5.0*p[1]*ga)*r*r*r;
}

double B1(double r,double ga,double p)
{
	return p*r*r*r*r;
}

double B2(double r,double ga,double p[])
{
	return (p[0] - p[1]*ga)*r*r*r*r;
}

double dB2(double r,double ga,double p[])
{
	return (4.0*p[0] - 5.0*p[1]*ga)*r*r*r;
}

double Visc_Us(double R,double cs,double u,double J)
{
	return cs*pow(J/R,2.0)*u;
}

double Visc_Ut(double R,double ct,double w,double J)
{
	return ct*R*J*w;
}

double Prod(double R,double A,double cs,double u,double gs)
{
	return IVsVs(gs,2.0)*pow(cs,2.0)*pow(R,3.0)*pow(u,2.0);
}
#endif
