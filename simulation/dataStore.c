#include <stdio.h>
#include <stdlib.h>
void fileprintout(FILE *a, double ***b, int m, int n, int q){
	int i;
	int j;
	int k;
	if(a==NULL){
		printf("creat file failed\n");
		exit(1);                                    //exited abnormly
	} 
	for (i=0;i<m;i++){
		for (j=0;j<n;j++){
			for (k=0;k<q;k++){
				fprintf(a,"%.10lf ",b[i][j][k]);
			}
		}
		fprintf(a,"\n");
	}
	fclose(a);
}