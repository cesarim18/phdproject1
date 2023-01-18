#ifndef _MAXMINH_
#define _MAXMINH_

#include <math.h>

double min(double a,double b){
	double aux;
	aux=a;
	if(b<aux) aux=b;
	return aux;
}

double max(double a,double b){
	double aux;
	aux=a;
	if(b>aux) aux=b;
	return aux;
}

double maximum(double a,double b,double c)
{
	double aux1,aux2;
	aux1=max(a,b);
	aux2=max(aux1,c);
	return aux2;
}

double minimum(double a,double b,double c)
{
	double aux1,aux2;
	aux1=min(a,b);
	aux2=min(aux1,c);
	return aux2;
}

double min_vec(double *vec,int dim)
{
	int i;
	double min;
	min = vec[0];
	for(i=1;i<dim;i++){
		if(vec[i] < min) min = vec[i];
	}return min;
}

double max_vec(double *vec,int dim)
{
	int i;
	double max;
	max = vec[0];
	for(i=1;i<dim;i++){
		if(vec[i] > max) max = vec[i];
	}return max;
}
#endif
