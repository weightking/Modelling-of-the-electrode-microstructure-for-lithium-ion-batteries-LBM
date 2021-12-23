#include <stdio.h>
#include <stdlib.h>

void freememory(double ***a,int m, int n, int q) {
	int i;
	int j;
	for(i = 0; i < m; i++){
		for (j=0; j<n; j++){
			free(a[i][j]);
			a[i][j]=NULL;
		}
	}
	for (i=0; i<m; i++){
		free(a[i]);
		a[i]=NULL;
	}
	free(a);
	a=NULL;
}

void freememoryTwoDimension(double **a, int m) {
	int i;
	for(i = 0; i < m; i++){
		free(a[i]);
		a[i]=NULL;
	}
	free(a);
	a=NULL;
}

void freememoryOneDimension(double *a){
	free(a);
	a=NULL;
}