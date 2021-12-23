#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <malloc.h>
#include <math.h>
#include "MemoryArrange.h"
#include "matrixMove.h"
#include "memoryFree.h"

double ***domainPorous=NULL;
double ***domainTwoPhase = NULL;
double ***domain=NULL;
double ***Domain1=NULL;
double ***Domain2=NULL;
double ***Domain3=NULL;
double ***Domain4=NULL;
double ***Domain5=NULL;
double ***Domain6=NULL;
double ***Domain7=NULL;
double ***Domain8=NULL;
double ***Domain9=NULL;
double ***Domain10=NULL;
double ***Domain11=NULL;
double ***Domain12=NULL;
double ***Domain13=NULL;
double ***Domain14=NULL;
double ***Domain15=NULL;
double ***Domain16=NULL;
double ***Domain17=NULL;
double ***Domain18=NULL;
double ***Domain19 = NULL;
double ***Domain20 = NULL;
double ***Domain21 = NULL;
double ***Domain22 = NULL;
double ***Domain23 = NULL;
double ***Domain24 = NULL;
double ***Domain25 = NULL;
double ***Domain26 = NULL;

void domainArrange(int m, int n, int q){
	domainPorous=memoryarrange(m,n,q);
	domain=memoryarrange(m,n,q);
	Domain1=memoryarrange(m,n,q);
	Domain2=memoryarrange(m,n,q);
	Domain3=memoryarrange(m,n,q);
	Domain4=memoryarrange(m,n,q);
	Domain5=memoryarrange(m,n,q);
	Domain6=memoryarrange(m,n,q);
//	Domain7=memoryarrange(m,n,q);
//	Domain8=memoryarrange(m,n,q);
//	Domain9=memoryarrange(m,n,q);
//	Domain10=memoryarrange(m,n,q);
//	Domain11=memoryarrange(m,n,q);
//	Domain12=memoryarrange(m,n,q);
//	Domain13=memoryarrange(m,n,q);
//	Domain14=memoryarrange(m,n,q);
//	Domain15=memoryarrange(m,n,q);
//	Domain16=memoryarrange(m,n,q);
//	Domain17=memoryarrange(m,n,q);
//	Domain18=memoryarrange(m,n,q);
//	Domain19 = memoryarrange(m, n, q);
//	Domain20 = memoryarrange(m, n, q);
//	Domain21 = memoryarrange(m, n, q);
//	Domain22 = memoryarrange(m, n, q);
//	Domain23 = memoryarrange(m, n, q);
//	Domain24 = memoryarrange(m, n, q);
//	Domain25 = memoryarrange(m, n, q);
//	Domain26 = memoryarrange(m, n, q);
}

void domainShift(int m, int n, int q){
	shift1(domain,Domain1,m,n,q);
	shift2(domain,Domain2,m,n,q);
	shift3(domain,Domain3,m,n,q);
	shift4(domain,Domain4,m,n,q);
	shift5(domain,Domain5,m,n,q);
	shift6(domain,Domain6,m,n,q);
//	shift7(domain,Domain7,m,n,q);
//	shift8(domain,Domain8,m,n,q);
//	shift9(domain,Domain9,m,n,q);
//	shift10(domain,Domain10,m,n,q);
//	shift11(domain,Domain11,m,n,q);
//	shift12(domain,Domain12,m,n,q);
//	shift13(domain,Domain13,m,n,q);
//	shift14(domain,Domain14,m,n,q);
//	shift15(domain,Domain15,m,n,q);
//	shift16(domain,Domain16,m,n,q);
//	shift17(domain,Domain17,m,n,q);
//	shift18(domain,Domain18,m,n,q);
//	shift19(domain, Domain19, m, n, q);
//	shift20(domain, Domain20, m, n, q);
//	shift21(domain, Domain21, m, n, q);
//	shift22(domain, Domain22, m, n, q);
//	shift23(domain, Domain23, m, n, q);
//	shift24(domain, Domain24, m, n, q);
//	shift25(domain, Domain25, m, n, q);
//	shift26(domain, Domain26, m, n, q);
}

void domainDefine(int mInitial, int m, int nInitial, int n, int qInitial, int q, double phase){
	int i;
	int j;
	int k;
	for (i=mInitial;i<m;i++){
		for (j=nInitial;j<n;j++){
			for (k=qInitial;k<q;k++){
				domain[i][j][k]=phase;
			}
		}
	}
}

void domainDefineInletHoles(int qInitial, int q, int Length, int Width, int number, int radius, double phase){
	int i;
	int j;
	int k;
	int n;
	int coordX;
	int coordY;
	double distance1;
	double distance2;
	double distance3;
	double distance4;
	double distance5;
	double distance6;
	double distance7;
	double distance8;
	double distance9;
	srand((unsigned)time(NULL));
	for (n=0;n<number;n++){
		coordX=rand()%Length;
		coordY=rand()%Width;
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				distance1=((i-coordX)*(i-coordX)+(j-coordY)*(j-coordY));
				distance2=((i-Length-coordX)*(i-Length-coordX)+(j-Width-coordY)*(j-Width-coordY));
				distance3=((i-coordX)*(i-coordX)+(j-Width-coordY)*(j-Width-coordY));
				distance4=((i-Length-coordX)*(i-Length-coordX)+(j-coordY)*(j-coordY));
				distance5=((i+Length-coordX)*(i+Length-coordX)+(j-coordY)*(j-coordY));
				distance6=((i-coordX)*(i-coordX)+(j+Width-coordY)*(j+Width-coordY));
				distance7=((i-Length-coordX)*(i-Length-coordX)+(j+Width-coordY)*(j+Width-coordY));
				distance8=((i+Length-coordX)*(i+Length-coordX)+(j-Width-coordY)*(j-Width-coordY));
				distance9=((i+Length-coordX)*(i+Length-coordX)+(j+Width-coordY)*(j+Width-coordY));
				if (distance1<=(radius*radius)||distance2<=(radius*radius)||distance3<=(radius*radius)||distance4<=(radius*radius)||distance5<=(radius*radius)||distance6<=(radius*radius)||distance7<=(radius*radius)||distance8<=(radius*radius)||distance9<=(radius*radius)){
					for (k=qInitial;k<q;k++){
						domain[i][j][k]=phase;
					}
				}
			}
		}
	}
}

void domainFreeMemory(int m, int n, int q){
	freememory(domainPorous,m,n,q);
	freememory(domain,m,n,q);
	freememory(Domain1,m,n,q);
	freememory(Domain2,m,n,q);
	freememory(Domain3,m,n,q);
	freememory(Domain4,m,n,q);
	freememory(Domain5,m,n,q);
	freememory(Domain6,m,n,q);
//	freememory(Domain7,m,n,q);
//	freememory(Domain8,m,n,q);
//	freememory(Domain9,m,n,q);
//	freememory(Domain10,m,n,q);
//	freememory(Domain11,m,n,q);
//	freememory(Domain12,m,n,q);
//	freememory(Domain13,m,n,q);
//	freememory(Domain14,m,n,q);
//	freememory(Domain15,m,n,q);
//	freememory(Domain16,m,n,q);
//	freememory(Domain17,m,n,q);
//	freememory(Domain18,m,n,q);
//	freememory(Domain19, m, n, q);
//	freememory(Domain20, m, n, q);
//	freememory(Domain21, m, n, q);
//	freememory(Domain22, m, n, q);
//	freememory(Domain23, m, n, q);
//	freememory(Domain24, m, n, q);
//	freememory(Domain25, m, n, q);
//	freememory(Domain26, m, n, q);
}