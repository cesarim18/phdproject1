// It is not recommended to put function definitions  
// in a header file. Ideally there should be only 
// function declarations. Purpose of this code is 
// to only demonstrate working of header files. 
#ifndef _CARDIAC_CYCLEH_
#define _CARDIAC_CYCLEH_

#include "create_vectors.h"
#include <complex.h>
#include "maxmin.h"

_Complex coefficients(int Nc,double Tc,int k,double V[])
{
	_Complex cof;
	cof = 0.0 + 0.0*I;
	int j;
	double t_j;
	for(j=0;j<Nc;j++){
		t_j = 1.0*(j+1)*Tc/Nc;
		cof = cof + cexp(-2.0*M_PI*I*k*t_j/Tc)*V[j];
	}return (1.0/Tc)*(Tc/(1.0*Nc))*cof;
}

double ghost(_Complex coeff[],int Nf,double T,double Tc)
{
	double u_gh;
	u_gh = creal(coeff[0]);
	int k;
	
	for(k=1;k<Nf+1;k++){
		u_gh = u_gh + 2.0*creal(coeff[k]*cexp(2.0*M_PI*I*k*T/Tc));
	}return u_gh;
}

#endif


/*
cardiac_cycle = load('x point01 step.csv');
N_cycle = size(cardiac_cycle,1);
T_cycle = cardiac_cycle(N_cycle,1);
N_fourier = 15;
coeff = zeros(N_fourier+1,1);
for n = 0:N_fourier
    a_n = 0;
    for j=1:N_cycle
        t_j = T_cycle/N_cycle*(j-1);
        a_n = a_n + 1/T_cycle*(T_cycle/N_cycle)*exp(-sqrt(-1)*n*2*pi/T_cycle*t_j)*cardiac_cycle(j,2);
    end
    coeff(n+1,1) = a_n;
end
cardiac_cycle_fourier = zeros(N_cycle,2);
cardiac_cycle_fourier(:,1) = cardiac_cycle(:,1);

a_n = coeff(1,1);
for j=1:N_cycle
    cardiac_cycle_fourier(j,2) = a_n;
end
for n=1:N_fourier
    a_n = coeff(n+1,1);
    for j=1:N_cycle
        t_j = T_cycle/N_cycle*(j-1);
        cardiac_cycle_fourier(j,2) = cardiac_cycle_fourier(j,2)+2*real(a_n*exp(sqrt(-1)*n*2*pi/T_cycle*t_j));
    end
end

u_ghost_max = 0;
for j = 1:N_cycle-1
    u_ghost_max = max(u_ghost_max,(cardiac_cycle_fourier(j+1,2)-cardiac_cycle_fourier(j,2))/(T_cycle/N_cycle));
end
u_ghost_max

*/

