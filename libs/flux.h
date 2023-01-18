// It is not recommended to put function definitions  
// in a header file. Ideally there should be only 
// function declarations. Purpose of this code is 
// to only demonstrate working of header files. 

#ifndef _FLUXH_
#define _FLUXH_

#include <math.h>
#include "maxmin.h"
#include "flux_extra.h"
#include "find_root.h"

void NumFluxG(double *C1,double *C2,double *C3,double *W1,double *W2,double *W3,
int D_s,int D_th,double *Ro,double *RoE,double *RoW,double *RoN,double *RoS,
double *Go,double *GoE,double *GoW,double *GoN,double *GoS,double par[],
double thp[],double th0[], double thm[],double alp[],double al0[],double alm[],
double dap[],double da0[],double dam[],double ddp[],double dd0[],double ddm[],
double *Ao,double *AoE,double *AoW,double *AoN,double *AoS)
{
	int i;
	double max_a,max_b,min_R,v_mm,delta_s,delta_th,beta,g,fac,eps;
	double gs,gt,nu;
	    max_a = par[0];
	    max_b = par[1];
	    min_R = par[2];
		 v_mm = par[3];
	  delta_s = par[4];
	 delta_th = par[5];
		 beta = par[6];
		    g = par[7];
		  fac = par[8];
		  eps = par[9];
		   gs = par[10];
		   gt = par[11];
		   nu = par[12];

	double PAs,PAo[3],PA1[2],PA2[3];
	PAs = IVs(gs,1.0);
	PAo[0] = IVs(gs,1.0);
	PAo[1] = 2.0*IVs(gs,2.0);
	PAo[2] = IVs(gs,3.0);
	PA1[0] = IVsVs(gs,1.0);
	PA1[1] = IVsVs(gs,2.0);
	PA2[0] = IVsVt(gs,gt,1.0);
	PA2[1] = 2.0*IVsVt(gs,gt,2.0);
	PA2[2] = IVsVt(gs,gt,3.0);

	double PAt[2],PB0[2],PB1,PB2[2];
	PAt[0] = IVt(gt,1.0); 
	PAt[1] = IVt(gt,2.0);
	PB0[0] = IVt(gt,3.0);
	PB0[1] = IVt(gt,4.0);
	PB1 = IVsVt(gs,gt,3.0);
	PB2[0] = IVtVt(gt,3.0);
	PB2[1] = IVtVt(gt,4.0);

	for(i=0;i<D_s*D_th;i++){
		C1[i]=0.0;
		C2[i]=0.0;
		C3[i]=0.0;
	}

	#pragma omp parallel for shared(max_a,max_b,min_R)
	for(i=0;i<(D_s-4)*(D_th-4);i++){

		int j,k,ic; //index

		// values for minmods
		double V0,V1,V2,V3,V4;
		double VV1,VV2,VV3,VV4;
	
		j = i%(D_s-4) + 2;
		k = i/(D_s-4) + 2;
		ic = j + D_s*k;

		// reconstruction of A in s-direction
		double QrWp,QrE0,QrW0,QrEm;
		V4 = ga(W1[ic+2],Ro[ic+2],sin(th0[k]),da0[j+2]);
		V3 = ga(W1[ic+1],Ro[ic+1],sin(th0[k]),da0[j+1]);
		V2 = ga(W1[ ic ],Ro[ ic ],sin(th0[k]),da0[ j ]);
		V1 = ga(W1[ic-1],Ro[ic-1],sin(th0[k]),da0[j-1]);
		V0 = ga(W1[ic-2],Ro[ic-2],sin(th0[k]),da0[j-2]);
		Reconstruction(delta_s,V0,V1,V2,V3,V4,v_mm,&QrEm,&QrW0,&QrE0,&QrWp);
	
		// reconstruction of A in theta-direction
		double QrSp,QrN0,QrS0,QrNm;
		V4 = ga(W1[ic+2*D_s],Ro[ic+2*D_s],sin(th0[k+2]),da0[j]);
		V3 = ga(W1[ ic+D_s ],Ro[ ic+D_s ],sin(th0[k+1]),da0[j]);
		V2 = ga(W1[   ic   ],Ro[   ic   ],sin(th0[ k ]),da0[j]);
		V1 = ga(W1[ ic-D_s ],Ro[ ic-D_s ],sin(th0[k-1]),da0[j]);
		V0 = ga(W1[ic-2*D_s],Ro[ic-2*D_s],sin(th0[k-2]),da0[j]);
		Reconstruction(delta_th,V0,V1,V2,V3,V4,v_mm,&QrNm,&QrS0,&QrN0,&QrSp);

		// reconstruction of Qs in s-direction
		double QsWp,QsE0,QsW0,QsEm;
		V4 = W2[ic+2];
		V3 = W2[ic+1];
		V2 = W2[ ic ];
		V1 = W2[ic-1];
		V0 = W2[ic-2];
		Reconstruction(delta_s,V0,V1,V2,V3,V4,v_mm,&QsEm,&QsW0,&QsE0,&QsWp);

		// reconstruction of Qs in theta-direction
		double QsSp,QsN0,QsS0,QsNm;
		V4 = W2[ic+2*D_s];
		V3 = W2[ ic+D_s ];
		V2 = W2[   ic   ];
		V1 = W2[ ic-D_s ];
		V0 = W2[ic-2*D_s];
		Reconstruction(delta_th,V0,V1,V2,V3,V4,v_mm,&QsNm,&QsS0,&QsN0,&QsSp);

		// reconstruction of Qt in s-direction
		double QtWp,QtE0,QtW0,QtEm;
		V4 = W3[ic+2];
		V3 = W3[ic+1];
		V2 = W3[ ic ];
		V1 = W3[ic-1];
		V0 = W3[ic-2];
		Reconstruction(delta_s,V0,V1,V2,V3,V4,v_mm,&QtEm,&QtW0,&QtE0,&QtWp);

		// reconstruction of Qth in theta-direction
		double QtSp,QtN0,QtS0,QtNm;
		V4 = W3[ic+2*D_s];
		V3 = W3[ ic+D_s ];	
		V2 = W3[   ic   ];
		V1 = W3[ ic-D_s ];
		V0 = W3[ic-2*D_s];
		Reconstruction(delta_th,V0,V1,V2,V3,V4,v_mm,&QtNm,&QtS0,&QtN0,&QtSp);

		// Calculi of Variables in Axial Direction
		double VWp[15],VE0[15],VW0[15],VEm[15];
		Variables(VWp,QrWp,QsWp,QtWp,RoW[ic+1],AoW[ic+1],GoW[ic+1],th0[k],dam[j+1],PAs,PAt,PAo,PB0,beta);
		Variables(VE0,QrE0,QsE0,QtE0,RoE[ ic ],AoE[ ic ],GoE[ ic ],th0[k],dap[ j ],PAs,PAt,PAo,PB0,beta);
		Variables(VW0,QrW0,QsW0,QtW0,RoW[ ic ],AoW[ ic ],GoW[ ic ],th0[k],dam[ j ],PAs,PAt,PAo,PB0,beta);
		Variables(VEm,QrEm,QsEm,QtEm,RoE[ic-1],AoE[ic-1],GoE[ic-1],th0[k],dap[j-1],PAs,PAt,PAo,PB0,beta);

		// Calculi of Variables in Angular Direction
		double VSp[15],VN0[15],VS0[15],VNm[15];
		Variables(VSp,QrSp,QsSp,QtSp,RoS[ic+D_s],AoS[ic+D_s],GoS[ic+D_s],thm[k+1],da0[j],PAs,PAt,PAo,PB0,beta);
		Variables(VN0,QrN0,QsN0,QtN0,RoN[  ic  ],AoN[  ic  ],GoN[  ic  ],thp[ k ],da0[j],PAs,PAt,PAo,PB0,beta);
		Variables(VS0,QrS0,QsS0,QtS0,RoS[  ic  ],AoS[  ic  ],GoS[  ic  ],thm[ k ],da0[j],PAs,PAt,PAo,PB0,beta);
		Variables(VNm,QrNm,QsNm,QtNm,RoN[ic-D_s],AoN[ic-D_s],GoN[ic-D_s],thp[k-1],da0[j],PAs,PAt,PAo,PB0,beta);

		// Calculi of Flux and Eigenvalues in Axial Direction
		double FWp[3],FE0[3],FW0[3],FEm[3];
		double maxWp,maxE0,maxW0,maxEm;
		double minWp,minE0,minW0,minEm;
		double aEmax,aEmin,aWmax,aWmin;

		Axial_Req(FWp,&maxWp,&minWp,VWp,PA1,PB1,PAs,PAo);
		Axial_Req(FE0,&maxE0,&minE0,VE0,PA1,PB1,PAs,PAo);
		aEmax = maximum(max(maxWp,maxE0),max(VWp[6],VE0[6]),0.0);
		aEmin = minimum(min(minWp,minE0),min(VWp[6],VE0[6]),0.0);

		Axial_Req(FW0,&maxW0,&minW0,VW0,PA1,PB1,PAs,PAo);
		Axial_Req(FEm,&maxEm,&minEm,VEm,PA1,PB1,PAs,PAo);
		aWmax = maximum(max(maxW0,maxEm),max(VW0[6],VEm[6]),0.0);
		aWmin = minimum(min(minW0,minEm),min(VW0[6],VEm[6]),0.0);

		// Calculi of Flux and Eigenvalues in Angular Direction
		double GSp[3],GN0[3],GS0[3],GNm[3];
		double maxSp,maxN0,maxS0,maxNm;
		double minSp,minN0,minS0,minNm;
		double bNmax,bNmin,bSmax,bSmin;

		Angular_Req(GSp,&maxSp,&minSp,VSp,PA2,PB2,PAt,PB0);
		Angular_Req(GN0,&maxN0,&minN0,VN0,PA2,PB2,PAt,PB0);

		bNmax = maximum(max(maxSp,maxN0),max(VSp[10],VN0[10]),0.0);
		bNmin = minimum(min(minSp,minN0),min(VSp[10],VN0[10]),0.0);

		Angular_Req(GS0,&maxS0,&minS0,VS0,PA2,PB2,PAt,PB0);
		Angular_Req(GNm,&maxNm,&minNm,VNm,PA2,PB2,PAt,PB0);

		bSmax = maximum(max(maxS0,maxNm),max(VS0[10],VNm[10]),0.0);
		bSmin = minimum(min(minS0,minNm),min(VS0[10],VNm[10]),0.0);

		// ---------------------------------------------------------------------
		//
		//							Numerical Fluxes
		//
		// ---------------------------------------------------------------------
		double HFW[3],HFE[3],HF[3];
		double HGN[3],HGS[3],HG[3];
		double H[3];

		HFE[0] = (aEmax*FE0[0] - aEmin*FWp[0])/(aEmax - aEmin) + (aEmax*aEmin/(aEmax - aEmin))*(VWp[2] - VE0[2]);
		HFW[0] = (aWmax*FEm[0] - aWmin*FW0[0])/(aWmax - aWmin) + (aWmax*aWmin/(aWmax - aWmin))*(VW0[2] - VEm[2]);
		 HF[0] = -(HFE[0] - HFW[0])/delta_s;
		HGN[0] = (bNmax*GN0[0] - bNmin*GSp[0])/(bNmax - bNmin) + (bNmax*bNmin/(bNmax - bNmin))*(VSp[2] - VN0[2]);
		HGS[0] = (bSmax*GNm[0] - bSmin*GS0[0])/(bSmax - bSmin) + (bSmax*bSmin/(bSmax - bSmin))*(VS0[2] - VNm[2]);
		 HG[0] = -(HGN[0] - HGS[0])/delta_th;
		  H[0] = HF[0] + HG[0];

		HFE[1] = (aEmax*FE0[1] - aEmin*FWp[1])/(aEmax - aEmin) + (aEmax*aEmin/(aEmax - aEmin))*(QsWp - QsE0);
		HFW[1] = (aWmax*FEm[1] - aWmin*FW0[1])/(aWmax - aWmin) + (aWmax*aWmin/(aWmax - aWmin))*(QsW0 - QsEm);
		 HF[1] = -(HFE[1] - HFW[1])/delta_s;
		HGN[1] = (bNmax*GN0[1] - bNmin*GSp[1])/(bNmax - bNmin) + (bNmax*bNmin/(bNmax - bNmin))*(QsSp - QsN0);
		HGS[1] = (bSmax*GNm[1] - bSmin*GS0[1])/(bSmax - bSmin) + (bSmax*bSmin/(bSmax - bSmin))*(QsS0 - QsNm);
		 HG[1] = -(HGN[1] - HGS[1])/delta_th;
		  H[1] = HF[1] + HG[1];

		HFE[2] = (aEmax*FE0[2] - aEmin*FWp[2])/(aEmax - aEmin) + (aEmax*aEmin/(aEmax - aEmin))*(QtWp - QtE0);
		HFW[2] = (aWmax*FEm[2] - aWmin*FW0[2])/(aWmax - aWmin) + (aWmax*aWmin/(aWmax - aWmin))*(QtW0 - QtEm);
		 HF[2] = -(HFE[2] - HFW[2])/delta_s;
		HGN[2] = (bNmax*GN0[2] - bNmin*GSp[2])/(bNmax - bNmin) + (bNmax*bNmin/(bNmax - bNmin))*(QtSp - QtN0);
		HGS[2] = (bSmax*GNm[2] - bSmin*GS0[2])/(bSmax - bSmin) + (bSmax*bSmin/(bSmax - bSmin))*(QtS0 - QtNm);
		 HG[2] = -(HGN[2] - HGS[2])/delta_th;
		  H[2] = HF[2] + HG[2];

		// ---------------------------------------------------------------------
		//
		//								Source terms
		//
		// ---------------------------------------------------------------------
		double ValN,ValS,ValE,ValW,ValC;

		double S2_Go,S2_Ao,S2_Pre,S2_Grav,S2_Prod,S2_Visc;
		double S3_Go,S3_Ao,S3_Pre,S3_Prod,S3_Visc;
		double S2_Ro,S2_ex,S3_Ro,S3_ex;

		// pressure in S2
		// Go derivative
		ValE = ApbarGo2(VE0[2],VE0[12],VE0[14],GoE[ic]);
		ValW = ApbarGo2(VW0[2],VW0[12],VW0[14],GoW[ic]);
		ValC = 0.5*(ValE + ValW);
		S2_Go = ValC*(GoE[ic] - GoW[ic])/delta_s;

		// Ao derivatives
		ValE = Jac(VE0[0]/QrE0,VE0[1]/QrE0)*VE0[12]/AoE[ic];
		ValW = Jac(VW0[0]/QrW0,VW0[1]/QrW0)*VW0[12]/AoW[ic];
		ValC = 0.5*(ValE + ValW);
		S2_Ro = ValC*(RoE[ic] - RoW[ic])/delta_s;

		ValE = pow(RoE[ic],3.0)*VE0[12]*ddp[j]/AoE[ic];
		ValW = pow(RoW[ic],3.0)*VW0[12]*ddm[j]/AoW[ic];
		S2_ex = 0.5*(ValE + ValW)*sin(th0[k])/3.0;

		S2_Ao = S2_Ro - S2_ex;

/*
		ValE = VE0[12]/AoE[ic];
		ValW = VW0[12]/AoW[ic];
		ValC = 0.5*(ValE + ValW);
		S2_Ao = ValC*(AoE[ic] - AoW[ic])/delta_s;
*/

		S2_Pre = -S2_Go + S2_Ao;

		// pressure in S3
		// Go derivative
		ValN = ApbarGo2(VN0[2],VN0[12],VN0[14],GoN[ic]);
		ValS = ApbarGo2(VS0[2],VS0[12],VS0[14],GoS[ic]);
		ValC = 0.5*(ValN + ValS);
		S3_Go = ValC*(GoN[ic] - GoS[ic])/delta_th;

		// Ao derivatives
		ValN = Jac(VN0[0]/QrN0,VN0[1]/QrN0)*VN0[12]/AoN[ic];
		ValS = Jac(VS0[0]/QrS0,VS0[1]/QrS0)*VS0[12]/AoS[ic];
		ValC = 0.5*(ValN + ValS);
		S3_Ro = ValC*(RoN[ic] - RoS[ic])/delta_th;

		ValN = pow(RoN[ic],3.0)*VN0[12]*cos(thp[k])/AoN[ic];
		ValS = pow(RoS[ic],3.0)*VS0[12]*cos(thm[k])/AoS[ic];
		S3_ex = 0.5*(ValN + ValS)*da0[j]/3.0;

		S3_Ao = S3_Ro - S3_ex;

/*
		ValN = VN0[12]/AoN[ic];
		ValS = VS0[12]/AoS[ic];
		ValC = 0.5*(ValN + ValS);
		S3_Ao = ValC*(AoN[ic] - AoS[ic])/delta_th;
*/
		S3_Pre = -S3_Go + S3_Ao;

		// gravity
		ValN = sin(al0[j])*VN0[2];
		ValS = sin(al0[j])*VS0[2];
		ValE = sin(alp[j])*VE0[2];
		ValW = sin(alm[j])*VW0[2];

		ValC = 0.25*(ValN + ValS + ValE + ValW);
		S2_Grav = g*ValC;

		// products
		ValN = Prod(VN0[0],VN0[2],VN0[4],VN0[6],gs)*sin(thp[k])*dd0[j];
		ValS = Prod(VS0[0],VS0[2],VS0[4],VS0[6],gs)*sin(thm[k])*dd0[j];
		ValE = Prod(VE0[0],VE0[2],VE0[4],VE0[6],gs)*sin(th0[k])*ddp[j];
		ValW = Prod(VW0[0],VW0[2],VW0[4],VW0[6],gs)*sin(th0[k])*ddm[j];
		ValC = 0.25*(ValN + ValS + ValE + ValW);		
		S2_Prod = ValC;

		ValN = Prod(VN0[0],VN0[2],VN0[4],VN0[6],gs)*cos(thp[k])*da0[j];
		ValS = Prod(VS0[0],VS0[2],VS0[4],VS0[6],gs)*cos(thm[k])*da0[j];
		ValE = Prod(VE0[0],VE0[2],VE0[4],VE0[6],gs)*cos(th0[k])*dap[j];
		ValW = Prod(VW0[0],VW0[2],VW0[4],VW0[6],gs)*cos(th0[k])*dam[j];
		ValC = 0.25*(ValN + ValS + ValE + ValW);		
		S3_Prod = ValC;

		// Viscosity
		ValN = Visc_Us(VN0[0],VN0[4],VN0[6],VN0[11]);
		ValS = Visc_Us(VS0[0],VS0[4],VS0[6],VS0[11]);
		ValE = Visc_Us(VE0[0],VE0[4],VE0[6],VE0[11]);
		ValW = Visc_Us(VW0[0],VW0[4],VW0[6],VW0[11]);
		ValC = 0.25*(ValN + ValS + ValE + ValW);
		S2_Visc = nu*gs*ValC;

		ValN = Visc_Ut(VN0[0],VN0[8],VN0[10],VN0[11]);
		ValS = Visc_Ut(VS0[0],VS0[8],VS0[10],VS0[11]);
		ValE = Visc_Ut(VE0[0],VE0[8],VE0[10],VE0[11]);
		ValW = Visc_Ut(VW0[0],VW0[8],VW0[10],VW0[11]);
		ValC = 0.25*(ValN + ValS + ValE + ValW);
		S3_Visc = nu*ValC/(gt+1.0);

		// finalice
		C1[ic] = H[0];
		C2[ic] = H[1] + S2_Pre - S2_Prod - S2_Visc - S2_Grav;
		C3[ic] = H[2] + S3_Pre - S3_Prod - S3_Visc;
		
		min_R = minimum(min_R,min(VSp[0],VN0[0]),min(VS0[0],VNm[0]));
		min_R = minimum(min_R,min(VWp[0],VE0[0]),min(VW0[0],VEm[0]));

		max_a = maximum(max_a,max(aEmax,-aEmin),max(aWmax,-aWmin));
		max_b = maximum(max_b,max(bNmax,-bNmin),max(bSmax,-bSmin));
	}par[0] = max_a;
	par[1] = max_b;
	par[2] = min_R;
}

#endif
