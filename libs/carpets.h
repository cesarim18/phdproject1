// It is not recommended to put function definitions  
// in a header file. Ideally there should be only 
// function declarations. Purpose of this code is 
// to only demonstrate working of header files. 
#ifndef _CARPETSH_
#define _CARPETSH_

void create_folders(int file)
{
	char carpet[100];
	sprintf(carpet,"Example");
	mkdir(carpet,0777);
	sprintf(carpet,"Example/Ex%d",file);
	mkdir(carpet,0777);
	sprintf(carpet,"Example/Ex%d/Videos",file);
	mkdir(carpet,0777);
	sprintf(carpet,"Example/Ex%d/Files",file);
	mkdir(carpet,0777);
	sprintf(carpet,"Example/Ex%d/IC",file);
	mkdir(carpet,0777);
	sprintf(carpet,"Example/Ex%d/IF",file);
	mkdir(carpet,0777);
	sprintf(carpet,"Example/Ex%d/Parameters",file);
	mkdir(carpet,0777);
}

#endif
