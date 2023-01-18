// It is not recommended to put function definitions  
// in a header file. Ideally there should be only 
// function declarations. Purpose of this code is 
// to only demonstrate working of header files. 
#ifndef _FIND_ROOTH_
#define _FIND_ROOTH_

#include "profiles.h"

double fa1(double R,double sint,double alfaps)
{
	return Q0(R,R*sint*alfaps);
}

double fb1(double R,double sint,double alfaps)
{
	return Jac(R,R*sint*alfaps);
}

double newtonraphson(double area,double sint,double alfaps,double tol)
{
	int i=0,ite=5;
	double xo,xant,dx,fx,error;
	xo=sqrt(2.0*area);
	error = 0.0;
	do{
		xant = xo;
		xo = xant - (fa1(xant,sint,alfaps) - area)/fb1(xant,sint,alfaps);
		if(i>0){
			error = fabs((xo-xant)/xo);
			if(error<=tol) break;
		}i++;
	}while(i<=ite);
	return xo;	
}

double findingroot(double area,double sint,double alfaps)
{
	//return sqrt(2.0*area);
	return newtonraphson(area,sint,alfaps,pow(10.0,-15.0));
}

#endif
