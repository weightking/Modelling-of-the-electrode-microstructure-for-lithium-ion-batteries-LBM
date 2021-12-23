#include <stdio.h>
#include <stdlib.h>
#include "MemoryArrange.h"

double ***gWall=NULL;

void gWallArrangeTwoPhase(int m, int n, int q){
	gWall=memoryarrange(m,n,q);
}

void gWallDefineTwoPhase(int mInitial, int m, int nInitial, int n, int qInitial, int q, double contactAngle){
	int i;
	int j;
	int k;
	for (i=mInitial;i<m;i++){
		for (j=nInitial;j<n;j++){
			for (k=qInitial;k<q;k++){
				gWall[i][j][k]=contactAngle;
			}
		}
	}
}