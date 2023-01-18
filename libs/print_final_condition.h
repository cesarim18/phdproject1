// It is not recommended to put function definitions  
// in a header file. Ideally there should be only 
// function declarations. Purpose of this code is 
// to only demonstrate working of header files. 
#ifndef _PRINT_FINAL_CONDITIONH_
#define _PRINT_FINAL_CONDITIONH_

void print_IF(double *R,double *A,double *u,double *w,double *L,double *Ut,
double *P,double *PT,double *Qs,double *Qt,double *VT,double *Vx,double *Vy,
double *Vz,double *Vex,double *Vey,double *Vez,int Ns,int Nth,int exa)
{
	int j,k;
	char eps_name[100];

	FILE *radius = NULL;
	sprintf(eps_name,"Example/Ex%d/IF/IF_R.dat",exa);
	radius = fopen(eps_name,"a");
	if(radius==NULL){
		puts("No se puede abrir IC_R!!! Ayuda Gerardo :(");
		exit(0);
	}
	for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(radius,"%.15e\n",R[j+Ns*k]);
		}
	}fclose(radius);
	radius = NULL;

/*
	FILE *area = NULL;
	sprintf(eps_name,"Example/Ex%d/IF/IF_A.dat",exa);
	area = fopen(eps_name,"a");
	if(area==NULL){
		puts("No se puede abrir IC_R!!! Ayuda Gerardo :(");
		exit(0);
	}
	for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(area,"%.15e\n",A[j+Ns*k]);
		}
	}fclose(area);
	area = NULL;
*/

	FILE *velu = NULL;
	sprintf(eps_name,"Example/Ex%d/IF/IF_u.dat",exa);
	velu = fopen(eps_name,"a");
	if(velu==NULL){
		puts("No se puede abrir el fileu1!!! Ayuda Gerardo :(");
		exit(0);
	}
	for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(velu,"%.15e\n",u[j+Ns*k]);
		}
	}fclose(velu);
	velu = NULL;

	FILE *velw = NULL;
	sprintf(eps_name,"Example/Ex%d/IF/IF_w.dat",exa);
	velw = fopen(eps_name,"a");
	if(velw==NULL){
		puts("No se puede abrir el fileu1!!! Ayuda Gerardo :(");
		exit(0);
	}
	for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(velw,"%.15e\n",w[j+Ns*k]);
		}
	}fclose(velw);
	velw = NULL;

	FILE *velUt = NULL;
	sprintf(eps_name,"Example/Ex%d/IF/IF_Ut.dat",exa);
	velUt = fopen(eps_name,"a");
	if(velUt==NULL){
		puts("No se puede abrir el fileu1!!! Ayuda Gerardo :(");
		exit(0);
	}
	for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(velUt,"%.15e\n",Ut[j+Ns*k]);
		}
	}fclose(velUt);
	velUt = NULL;

	FILE *momL = NULL;
	sprintf(eps_name,"Example/Ex%d/IF/IF_L.dat",exa);
	momL = fopen(eps_name,"a");
	if(momL==NULL){
		puts("No se puede abrir el fileu1!!! Ayuda Gerardo :(");
		exit(0);
	}
	for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(momL,"%.15e\n",L[j+Ns*k]);
		}
	}fclose(momL);
	momL = NULL;

	FILE *pressure = NULL;
	sprintf(eps_name,"Example/Ex%d/IF/IF_P.dat",exa);
	pressure = fopen(eps_name,"a");
	if(pressure==NULL){
		puts("No se puede abrir vessel_x!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(pressure,"%.15e\n",P[j+Ns*k]);
		}
	}fclose(pressure);
	pressure = NULL;

/*
	FILE *pressureT = NULL;
	sprintf(eps_name,"Example/Ex%d/IF/IF_PT.dat",exa);
	pressureT = fopen(eps_name,"a");
	if(pressureT==NULL){
		puts("No se puede abrir vessel_x!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(pressureT,"%.15e\n",PT[j+Ns*k]);
		}
	}fclose(pressureT);
	pressureT = NULL;
*/

/*
	FILE *cons_Qs = NULL;
	sprintf(eps_name,"Example/Ex%d/IF/IF_Qs.dat",exa);
	cons_Qs = fopen(eps_name,"a");
	if(cons_Qs==NULL){
		puts("No se puede abrir vessel_x!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(cons_Qs,"%.15e\n",Qs[j+Ns*k]);
		}
	}fclose(cons_Qs);
	cons_Qs = NULL;
*/

/*
	FILE *cons_Qt = NULL;
	sprintf(eps_name,"Example/Ex%d/IF/IF_Qt.dat",exa);
	cons_Qt = fopen(eps_name,"a");
	if(cons_Qt==NULL){
		puts("No se puede abrir vessel_x!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(cons_Qt,"%.15e\n",Qt[j+Ns*k]);
		}
	}fclose(cons_Qt);
	cons_Qt = NULL;
*/

/*
	FILE *velocity_total = NULL;
	sprintf(eps_name,"Example/Ex%d/IF/IF_VT.dat",exa);
	velocity_total = fopen(eps_name,"a");
	if(velocity_total==NULL){
		puts("No se puede abrir vessel_x!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(velocity_total,"%.15e\n",VT[j+Ns*k]);
		}
	}fclose(velocity_total);
	velocity_total = NULL;
*/

	FILE *vessel_x = NULL;
	sprintf(eps_name,"Example/Ex%d/IF/IF_Vx.dat",exa);
	vessel_x = fopen(eps_name,"a");
	if(vessel_x==NULL){
		puts("No se puede abrir vessel_x!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(vessel_x,"%.15e\n",Vx[j+Ns*k]);
		}
	}fclose(vessel_x);
	vessel_x = NULL;

	FILE *vessel_y = NULL;
	sprintf(eps_name,"Example/Ex%d/IF/IF_Vy.dat",exa);
	vessel_y = fopen(eps_name,"a");
	if(vessel_y==NULL){
		puts("No se puede abrir vessel_y!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(vessel_y,"%.15e\n",Vy[j+Ns*k]);
		}
	}fclose(vessel_y);
	vessel_y = NULL;

	FILE *vessel_z = NULL;
	sprintf(eps_name,"Example/Ex%d/IF/IF_Vz.dat",exa);
	vessel_z = fopen(eps_name,"a");
	if(vessel_z==NULL){
		puts("No se puede abrir vessel_z!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(vessel_z,"%.15e\n",Vz[j+Ns*k]);
		}
	}fclose(vessel_z);
	vessel_z = NULL;

	FILE *velocity_x = NULL;
	sprintf(eps_name,"Example/Ex%d/IF/IF_Velx.dat",exa);
	velocity_x = fopen(eps_name,"a");
	if(velocity_x==NULL){
		puts("No se puede abrir vessel_x!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(velocity_x,"%.15e\n",Vex[j+Ns*k]);
		}
	}fclose(velocity_x);
	velocity_x = NULL;

	FILE *velocity_y = NULL;
	sprintf(eps_name,"Example/Ex%d/IF/IF_Vely.dat",exa);
	velocity_y = fopen(eps_name,"a");
	if(velocity_y==NULL){
		puts("No se puede abrir vessel_x!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(velocity_y,"%.15e\n",Vey[j+Ns*k]);
		}
	}fclose(velocity_y);
	velocity_y = NULL;

	FILE *velocity_z = NULL;
	sprintf(eps_name,"Example/Ex%d/IF/IF_Velz.dat",exa);
	velocity_z = fopen(eps_name,"a");
	if(velocity_z==NULL){
		puts("No se puede abrir vessel_x!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(velocity_z,"%.15e\n",Vez[j+Ns*k]);
		}
	}fclose(velocity_z);
	velocity_z = NULL;

	printf("\n			INITIAL CONDITION HAS BEEN SAVED\n");
}

#endif
