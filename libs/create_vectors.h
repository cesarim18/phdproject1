// It is not recommended to put function definitions  
// in a header file. Ideally there should be only 
// function declarations. Purpose of this code is 
// to only demonstrate working of header files. 
#ifndef _CREATE_VECTORSH_
#define _CREATE_VECTORSH_

bool Crea1(double **A,int n)
{
	(*A)=(double*)malloc(n*sizeof(double));
	if((*A)==NULL){
		puts("No hay espacio suficiente en la RAM.");
		system("pause");
		exit(0);
	}
}

bool Crea(double ***A,int m,int n)
{
    (*A)=(double**)malloc(m*sizeof(double*));
    if((*A)==NULL){
        puts("1No hay espacio suficiente en la RAM.");
        exit(0);
    }int i,j;
	for(i=0;i<m;i++){
        (*A)[i]=(double*)malloc((n)*sizeof(double));
        if((*A)[i]==NULL){
            printf("2No hay espacio suficiente en la RAM.");
            for(j=0;j<i;j++){
                free((*A)[j]);
            }free((*A));
            exit(0);
        }
    }
}
#endif
