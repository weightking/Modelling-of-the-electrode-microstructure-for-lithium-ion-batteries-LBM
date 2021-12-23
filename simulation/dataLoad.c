#include <stdio.h>
#include <stdlib.h>

void dataload(FILE *a, double ***b, int mInitial,int m, int nInitial, int n, int qInitial, int q){
	int i;
	int j;
	int k;

	if(a==NULL){
		printf("creat file failed\n");
		exit(1);                                    //exited abnormly
	}
	for (i=mInitial;i<m;i++){
		for (j=nInitial;j<n;j++){
			for (k=qInitial;k<q;k++){
				fscanf(a, "%lf",&b[i][j][k]);
			}
		}
	}
	fclose(a);
}