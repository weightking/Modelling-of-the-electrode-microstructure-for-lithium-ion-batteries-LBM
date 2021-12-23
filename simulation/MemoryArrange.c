#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

double ***memoryarrange(int m, int n, int q){
	int i;
	int j;
	double ***a=NULL;
	a=(double ***)malloc(m*sizeof(double **));
	if (a==NULL){
		printf("can not arrange space for a");
		exit(1);                                    //exited abnormly
	}
	for (i=0;i<m;i++){
		a[i]=(double **)malloc(n*sizeof(double*));
		if (a[i]==NULL){
			printf("can not arrange space for a[i]");
			exit(1);                                //exited abnormly
		}
		for (j=0;j<n;j++){
			a[i][j]=(double *)malloc(q*sizeof(double));
			if (a[i][j]==NULL){
				printf("can not arrange space for a[i][j]");
				exit(1);
			}
		}
	}
	return a;
}

double **memoryarrangeTwoDimension(int m, int n){
	int i;
	double **a=NULL;
	a=(double **)malloc(m*sizeof(double *));
	if (a==NULL){
		printf("can not arrange space for a");
		exit(1);                                    //exited abnormly
	}
	for (i=0;i<m;i++){
		a[i]=(double *)malloc(n*sizeof(double));
		if (a[i]==NULL){
			printf("can not arrange space for a[i]");
			exit(1);                                //exited abnormly
		}
	}
	return a;
}

double *memoryarrangeOneDimension(int m){
	double *a=NULL;
	a=(double *)malloc(m*sizeof(double));
	if (a==NULL){
		printf("can not arrange space for a");
		exit(1);
	}
	return a;
}