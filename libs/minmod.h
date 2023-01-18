#ifndef _MINMODH_
#define _MINMODH_

#include <math.h>
#include "maxmin.h"

int sgn(double x)
{
	if(x>0.0) return 1;
	if(x<0.0) return -1;
	return 0;
}

double minmod2(double a,double b,double c)
{
	double aux2;
	if(a==0.0 || b==0.0 || c==0.0 || sgn(a)+sgn(b)==0.0 || sgn(b)+sgn(c)==0.0 || sgn(c)+sgn(a)==0.0) aux2=0.0;
	else{
		aux2=fabs(a);
		if(fabs(b)<=aux2) aux2=fabs(b);
		if(fabs(c)<=aux2) aux2=fabs(c);
		aux2=sgn(a)*aux2;
	}return aux2;
}

double minmod3(double a,double b,double c)
{
	int aux;
	aux = sgn(a) + sgn(b) + sgn(c);
	if(fabs(aux) == 3){
		if(sgn(a) == -1) return maximum(a,b,c);
		if(sgn(a) == 1) return minimum(a,b,c); 
	}else return 0;
}

double minmod(double a,double b,double c,double dx,double theta)
{
	double aux1,aux2,aux3,mm;
	aux1 = theta*(b-c)/dx;
	aux2 =   0.5*(a-c)/dx;
	aux3 = theta*(a-b)/dx;
	mm = minmod3(aux1,aux2,aux3);
	return mm;
}

#endif
