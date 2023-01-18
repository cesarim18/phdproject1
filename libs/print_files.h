// It is not recommended to put function definitions  
// in a header file. Ideally there should be only 
// function declarations. Purpose of this code is 
// to only demonstrate working of header files. 
#ifndef _PRINT_FILESH_
#define _PRINT_FILESH_

void ImprimirValores(double dt,double max_a,double max_b,double min_R,double u_gh,double wtime,double time,int carpet)
{
	char eps_name[100];

	FILE *salida_dt = NULL;
	sprintf(eps_name,"Example/Ex%d/Parameters/dt.dat",carpet);
	salida_dt = fopen(eps_name,"a");
	if(salida_dt==NULL){
		puts("No se pudo abrir el archivo sobre dt.");
		exit(0);
	}fprintf(salida_dt,"%.15e\n",dt);
	fclose(salida_dt);
    salida_dt = NULL;

	FILE *salida_max_a = NULL;
	sprintf(eps_name,"Example/Ex%d/Parameters/max_a.dat",carpet);
	salida_max_a = fopen(eps_name,"a");
	if(salida_max_a==NULL){
		puts("No se pudo abrir el archivo sobre max_a.");
		exit(0);
	}fprintf(salida_max_a,"%.15e\n",max_a);
	fclose(salida_max_a);
    salida_max_a = NULL;

	FILE *salida_max_b = NULL;
	sprintf(eps_name,"Example/Ex%d/Parameters/max_b.dat",carpet);
	salida_max_b = fopen(eps_name,"a");
	if(salida_max_b==NULL){
		puts("No se pudo abrir el archivo sobre max_b.");
		exit(0);
	}fprintf(salida_max_b,"%.15e\n",max_b);
	fclose(salida_max_b);
    salida_max_b = NULL;

	FILE *salida_min_R = NULL;
	sprintf(eps_name,"Example/Ex%d/Parameters/min_R.dat",carpet);
	salida_min_R = fopen(eps_name,"a");
	if(salida_min_R==NULL){
		puts("No se pudo abrir el archivo sobre min_R.");
		exit(0);
	}fprintf(salida_min_R,"%.15e\n",min_R);
	fclose(salida_min_R);
    salida_min_R = NULL;

	FILE *salida_u_gh = NULL;
	sprintf(eps_name,"Example/Ex%d/Parameters/u_gh.dat",carpet);
	salida_u_gh = fopen(eps_name,"a");
	if(salida_u_gh==NULL){
		puts("No se pudo abrir el archivo sobre u_gh.");
		exit(0);
	}fprintf(salida_u_gh,"%.15e\n",u_gh);
	fclose(salida_u_gh);
    salida_u_gh = NULL;

	FILE *salida_wtime = NULL;
	sprintf(eps_name,"Example/Ex%d/Parameters/wtime.dat",carpet);
	salida_wtime = fopen(eps_name,"a");
	if(salida_wtime==NULL){
		puts("No se pudo abrir el archivo sobre wtime.");
		exit(0);
	}fprintf(salida_wtime,"%.15e\n",wtime);
	fclose(salida_wtime);
    salida_wtime = NULL;

	FILE *salida_time = NULL;
	sprintf(eps_name,"Example/Ex%d/Parameters/time.dat",carpet);
	salida_time = fopen(eps_name,"a");
	if(salida_time==NULL){
		puts("No se pudo abrir el archivo sobre time.");
		exit(0);
	}fprintf(salida_time,"%.15e\n",time);
	fclose(salida_time);
    salida_time = NULL;}

void print_3D(double *R,double *A,double *u,double *w,double *L,double *Ut,
double *P,double *PT,double *Qs,double *Qt,double *VT,int Ns,int Nth,int carpet)
{
//	print_3D(R,A,u,w,L,P,PT,Qs,Qt,VT,D_s,N_th+1,file);

	int j,k;
	char eps_name[100];

	FILE *radius = NULL;
	sprintf(eps_name,"Example/Ex%d/Files/R.dat",carpet);
	radius = fopen(eps_name,"a");
	if(radius==NULL){
		puts("No se puede abrir el fileu1!!! Ayuda Gerardo :(");
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
	sprintf(eps_name,"Example/Ex%d/Files/A.dat",carpet);
	area = fopen(eps_name,"a");
	if(area==NULL){
		puts("No se puede abrir el fileu1!!! Ayuda Gerardo :(");
		exit(0);
	}
	for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(area,"%.15e\n",A[j+Ns*k]);
		}
	}fclose(area);
	area = NULL;
*/

	FILE *velus = NULL;
	sprintf(eps_name,"Example/Ex%d/Files/u.dat",carpet);
	velus = fopen(eps_name,"a");
	if(velus==NULL){
		puts("No se puede abrir el fileu1!!! Ayuda Gerardo :(");
		exit(0);
	}
	for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(velus,"%.15e\n",u[j+Ns*k]);
		}
	}//printf(" Save Us");
	fclose(velus);
	velus = NULL;

	FILE *velw = NULL;
	sprintf(eps_name,"Example/Ex%d/Files/w.dat",carpet);
	velw = fopen(eps_name,"a");
	if(velw==NULL){
		puts("No se puede abrir el fileu1!!! Ayuda Gerardo :(");
		exit(0);
	}
	for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(velw,"%.15e\n",w[j+Ns*k]);
		}
	}//printf(" Save Ht");
	fclose(velw);
	velw = NULL;

	FILE *velUt = NULL;
	sprintf(eps_name,"Example/Ex%d/Files/Ut.dat",carpet);
	velUt = fopen(eps_name,"a");
	if(velUt==NULL){
		puts("No se puede abrir el fileu1!!! Ayuda Gerardo :(");
		exit(0);
	}
	for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(velUt,"%.15e\n",Ut[j+Ns*k]);
		}
	}//printf(" Save Ht");
	fclose(velUt);
	velUt = NULL;

	FILE *velHt = NULL;
	sprintf(eps_name,"Example/Ex%d/Files/L.dat",carpet);
	velHt = fopen(eps_name,"a");
	if(velHt==NULL){
		puts("No se puede abrir el fileu1!!! Ayuda Gerardo :(");
		exit(0);
	}
	for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(velHt,"%.15e\n",L[j+Ns*k]);
		}
	}//printf(" Save Ht");
	fclose(velHt);
	velHt = NULL;

/*
	FILE *cons_Qs = NULL;
	sprintf(eps_name,"Example/Ex%d/Files/Qs.dat",carpet);
	cons_Qs = fopen(eps_name,"a");
	if(cons_Qs==NULL){
		puts("No se puede abrir el fileu1!!! Ayuda Gerardo :(");
		exit(0);
	}
	for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(cons_Qs,"%.15e\n",Qs[j+Ns*k]);
		}
	}//printf(" Save Ro");
	fclose(cons_Qs);
	cons_Qs = NULL;
*/

/*
	FILE *cons_Qt = NULL;
	sprintf(eps_name,"Example/Ex%d/Files/Qt.dat",carpet);
	cons_Qt = fopen(eps_name,"a");
	if(cons_Qt==NULL){
		puts("No se puede abrir el fileu1!!! Ayuda Gerardo :(");
		exit(0);
	}
	for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(cons_Qt,"%.15e\n",Qt[j+Ns*k]);
		}
	}//printf(" Save Ro");
	fclose(cons_Qt);
	cons_Qt = NULL;
*/

/*
	FILE *veltot = NULL;
	sprintf(eps_name,"Example/Ex%d/Files/VT.dat",carpet);
	veltot = fopen(eps_name,"a");
	if(veltot==NULL){
		puts("No se puede abrir el fileu1!!! Ayuda Gerardo :(");
		exit(0);
	}
	for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(veltot,"%.15e\n",VT[j+Ns*k]);
		}
	}//printf(" Save Ro");
	fclose(veltot);
	veltot = NULL;
*/

	FILE *pressure = NULL;
	sprintf(eps_name,"Example/Ex%d/Files/P.dat",carpet);
	pressure = fopen(eps_name,"a");
	if(pressure==NULL){
		puts("No se puede abrir el fileu1!!! Ayuda Gerardo :(");
		exit(0);
	}
	for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(pressure,"%.15e\n",P[j+Ns*k]);
		}
	}//printf(" Save Ro");
	fclose(pressure);
	pressure = NULL;

	FILE *pressureT = NULL;
	sprintf(eps_name,"Example/Ex%d/Files/PT.dat",carpet);
	pressureT = fopen(eps_name,"a");
	if(pressureT==NULL){
		puts("No se puede abrir el fileu1!!! Ayuda Gerardo :(");
		exit(0);
	}
	for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(pressureT,"%.15e\n",PT[j+Ns*k]);
		}
	}//printf(" Save Ro");
	fclose(pressureT);
	pressureT = NULL;
}

void print_vec_pos(double *Vx,double *Vy,double *Vz,int Ns,int Nth,int carpet)
{
	int j,k;
	char eps_name[100];

	FILE *vessel_x = NULL;
	sprintf(eps_name,"Example/Ex%d/Files/Vx.dat",carpet);
	vessel_x = fopen(eps_name,"a");
	if(vessel_x==NULL){
		puts("No se puede abrir vessel_x!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(vessel_x,"%.15e\n",Vx[j+Ns*k]);
		}
	}//printf(" Save Vx");
	fclose(vessel_x);
	vessel_x = NULL;

	FILE *vessel_y = NULL;
	sprintf(eps_name,"Example/Ex%d/Files/Vy.dat",carpet);
	vessel_y = fopen(eps_name,"a");
	if(vessel_y==NULL){
		puts("No se puede abrir vessel_y!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(vessel_y,"%.15e\n",Vy[j+Ns*k]);
		}
	}//printf(" Save Vy");
	fclose(vessel_y);
	vessel_y = NULL;

	FILE *vessel_z = NULL;
	sprintf(eps_name,"Example/Ex%d/Files/Vz.dat",carpet);
	vessel_z = fopen(eps_name,"a");
	if(vessel_z==NULL){
		puts("No se puede abrir vessel_z!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(vessel_z,"%.15e\n",Vz[j+Ns*k]);
		}
	}//printf(" Save Vz");
	fclose(vessel_z);
	vessel_z = NULL;
}

void print_vec_vel(double *Vx,double *Vy,double *Vz,int Ns,int Nth,int carpet)
{
	int j,k;
	char eps_name[100];

	FILE *vessel_x = NULL;
	sprintf(eps_name,"Example/Ex%d/Files/Velx.dat",carpet);
	vessel_x = fopen(eps_name,"a");
	if(vessel_x==NULL){
		puts("No se puede abrir velocity_x!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(vessel_x,"%.15e\n",Vx[j+Ns*k]);
		}
	}//printf(" Save Vx");
	fclose(vessel_x);
	vessel_x = NULL;

	FILE *vessel_y = NULL;
	sprintf(eps_name,"Example/Ex%d/Files/Vely.dat",carpet);
	vessel_y = fopen(eps_name,"a");
	if(vessel_y==NULL){
		puts("No se puede abrir velocity_y!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(vessel_y,"%.15e\n",Vy[j+Ns*k]);
		}
	}//printf(" Save Vy");
	fclose(vessel_y);
	vessel_y = NULL;

	FILE *vessel_z = NULL;
	sprintf(eps_name,"Example/Ex%d/Files/Velz.dat",carpet);
	vessel_z = fopen(eps_name,"a");
	if(vessel_z==NULL){
		puts("No se puede abrir velocity_z!!! Ayuda Gerardo :(");
		exit(0);
	}for(k=0;k<Nth;k++){
		for(j=0;j<Ns;j++){
			fprintf(vessel_z,"%.15e\n",Vz[j+Ns*k]);
		}
	}//printf(" Save Vz");
	fclose(vessel_z);
	vessel_z = NULL;
}


#endif
