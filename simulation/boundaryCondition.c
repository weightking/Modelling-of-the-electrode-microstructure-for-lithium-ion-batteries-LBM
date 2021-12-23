#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "MemoryArrange.h"
#include "memoryFree.h"
#include "twoPhaseField.h"
#include "domain.h"
#include "concentrationField.h"
#include "ghostFlowField.h"
#include "electrodePotential.h"
#include "electrolytePotential.h"
#include "OutMemoryArrange.h"
#include "equilibrium.h"
#include "matrixMove.h"

double **First0=NULL;
double **First1=NULL;
double **First2=NULL;
double **First3=NULL;
double **First4=NULL;
double **First5=NULL;
double **First6=NULL;
double **First7=NULL;
double **First8=NULL;
double **First9=NULL;
double **First10=NULL;
double **First11=NULL;
double **First12=NULL;
double **First13=NULL;
double **First14=NULL;
double **First15=NULL;
double **First16=NULL;
double **First17=NULL;
double **First18=NULL;
double **FirstFx=NULL;
double **FirstFy=NULL;
double **FirstFz=NULL;

double **Last0=NULL;
double **Last1=NULL;
double **Last2=NULL;
double **Last3=NULL;
double **Last4=NULL;
double **Last5=NULL;
double **Last6=NULL;
double **Last7=NULL;
double **Last8=NULL;
double **Last9=NULL;
double **Last10=NULL;
double **Last11=NULL;
double **Last12=NULL;
double **Last13=NULL;
double **Last14=NULL;
double **Last15=NULL;
double **Last16=NULL;
double **Last17=NULL;
double **Last18=NULL;
double **LastFx=NULL;
double **LastFy=NULL;
double **LastFz=NULL;
double **LastRho=NULL;
double **LastUx=NULL;
double **LastUy=NULL;
double **LastUz=NULL;
double **LastPressure2=NULL;
double **LastPressure1=NULL;

double **Last0Concentration=NULL;
double **Last1Concentration=NULL;
double **Last2Concentration=NULL;
double **Last3Concentration=NULL;
double **Last4Concentration=NULL;
double **Last5Concentration=NULL;
double **Last6Concentration=NULL;
double **Last7Concentration=NULL;
double **Last8Concentration=NULL;
double **Last9Concentration=NULL;
double **Last10Concentration=NULL;
double **Last11Concentration=NULL;
double **Last12Concentration=NULL;
double **Last13Concentration=NULL;
double **Last14Concentration=NULL;
double **Last15Concentration=NULL;
double **Last16Concentration=NULL;
double **Last17Concentration=NULL;
double **Last18Concentration=NULL;

double **First0Ghost=NULL;
double **First1Ghost=NULL;
double **First2Ghost=NULL;
double **First3Ghost=NULL;
double **First4Ghost=NULL;
double **First5Ghost=NULL;
double **First6Ghost=NULL;
double **First7Ghost=NULL;
double **First8Ghost=NULL;
double **First9Ghost=NULL;
double **First10Ghost=NULL;
double **First11Ghost=NULL;
double **First12Ghost=NULL;
double **First13Ghost=NULL;
double **First14Ghost=NULL;
double **First15Ghost=NULL;
double **First16Ghost=NULL;
double **First17Ghost=NULL;
double **First18Ghost=NULL;

double **Last0Ghost=NULL;
double **Last1Ghost=NULL;
double **Last2Ghost=NULL;
double **Last3Ghost=NULL;
double **Last4Ghost=NULL;
double **Last5Ghost=NULL;
double **Last6Ghost=NULL;
double **Last7Ghost=NULL;
double **Last8Ghost=NULL;
double **Last9Ghost=NULL;
double **Last10Ghost=NULL;
double **Last11Ghost=NULL;
double **Last12Ghost=NULL;
double **Last13Ghost=NULL;
double **Last14Ghost=NULL;
double **Last15Ghost=NULL;
double **Last16Ghost=NULL;
double **Last17Ghost=NULL;
double **Last18Ghost=NULL;

void FirstlineArrangeTwoPhase(int n, int q){
	First0=memoryarrangeTwoDimension(n,q);
	First1=memoryarrangeTwoDimension(n,q);
	First2=memoryarrangeTwoDimension(n,q);
	First3=memoryarrangeTwoDimension(n,q);
	First4=memoryarrangeTwoDimension(n,q);
	First5=memoryarrangeTwoDimension(n,q);
	First6=memoryarrangeTwoDimension(n,q);
	First7=memoryarrangeTwoDimension(n,q);
	First8=memoryarrangeTwoDimension(n,q);
	First9=memoryarrangeTwoDimension(n,q);
	First10=memoryarrangeTwoDimension(n,q);
	First11=memoryarrangeTwoDimension(n,q);
	First12=memoryarrangeTwoDimension(n,q);
	First13=memoryarrangeTwoDimension(n,q);
	First14=memoryarrangeTwoDimension(n,q);
	First15=memoryarrangeTwoDimension(n,q);
	First16=memoryarrangeTwoDimension(n,q);
	First17=memoryarrangeTwoDimension(n,q);
	First18=memoryarrangeTwoDimension(n,q);
	FirstFx=memoryarrangeTwoDimension(n,q);
	FirstFy=memoryarrangeTwoDimension(n,q);
	FirstFz=memoryarrangeTwoDimension(n,q);	
}

void LastlineArrangeTwoPhase(int m, int n){
	Last0=memoryarrangeTwoDimension(m,n);
	Last1=memoryarrangeTwoDimension(m,n);
	Last2=memoryarrangeTwoDimension(m,n);
	Last3=memoryarrangeTwoDimension(m,n);
	Last4=memoryarrangeTwoDimension(m,n);
	Last5=memoryarrangeTwoDimension(m,n);
	Last6=memoryarrangeTwoDimension(m,n);
	Last7=memoryarrangeTwoDimension(m,n);
	Last8=memoryarrangeTwoDimension(m,n);
	Last9=memoryarrangeTwoDimension(m,n);
	Last10=memoryarrangeTwoDimension(m,n);
	Last11=memoryarrangeTwoDimension(m,n);
	Last12=memoryarrangeTwoDimension(m,n);
	Last13=memoryarrangeTwoDimension(m,n);
	Last14=memoryarrangeTwoDimension(m,n);
	Last15=memoryarrangeTwoDimension(m,n);
	Last16=memoryarrangeTwoDimension(m,n);
	Last17=memoryarrangeTwoDimension(m,n);
	Last18=memoryarrangeTwoDimension(m,n);
	LastFx=memoryarrangeTwoDimension(m,n);
	LastFy=memoryarrangeTwoDimension(m,n);
	LastFz=memoryarrangeTwoDimension(m,n);
	LastRho=memoryarrangeTwoDimension(m,n);
	LastUx=memoryarrangeTwoDimension(m,n);
	LastUy=memoryarrangeTwoDimension(m,n);
	LastUz=memoryarrangeTwoDimension(m,n);
	LastPressure1=memoryarrangeTwoDimension(m,n);
	LastPressure2=memoryarrangeTwoDimension(m,n);
}

void LastlineArrangeConcentration(int n, int q){
	Last0Concentration=memoryarrangeTwoDimension(n,q);
	Last1Concentration=memoryarrangeTwoDimension(n,q);
	Last2Concentration=memoryarrangeTwoDimension(n,q);
	Last3Concentration=memoryarrangeTwoDimension(n,q);
	Last4Concentration=memoryarrangeTwoDimension(n,q);
	Last5Concentration=memoryarrangeTwoDimension(n,q);
	Last6Concentration=memoryarrangeTwoDimension(n,q);
	Last7Concentration=memoryarrangeTwoDimension(n,q);
	Last8Concentration=memoryarrangeTwoDimension(n,q);
	Last9Concentration=memoryarrangeTwoDimension(n,q);
	Last10Concentration=memoryarrangeTwoDimension(n,q);
	Last11Concentration=memoryarrangeTwoDimension(n,q);
	Last12Concentration=memoryarrangeTwoDimension(n,q);
	Last13Concentration=memoryarrangeTwoDimension(n,q);
	Last14Concentration=memoryarrangeTwoDimension(n,q);
	Last15Concentration=memoryarrangeTwoDimension(n,q);
	Last16Concentration=memoryarrangeTwoDimension(n,q);
	Last17Concentration=memoryarrangeTwoDimension(n,q);
	Last18Concentration=memoryarrangeTwoDimension(n,q);
}

void FirstlineArrangeVelocityGhost(int n, int q){
	First0Ghost=memoryarrangeTwoDimension(n,q);
	First1Ghost=memoryarrangeTwoDimension(n,q);
	First2Ghost=memoryarrangeTwoDimension(n,q);
	First3Ghost=memoryarrangeTwoDimension(n,q);
	First4Ghost=memoryarrangeTwoDimension(n,q);
	First5Ghost=memoryarrangeTwoDimension(n,q);
	First6Ghost=memoryarrangeTwoDimension(n,q);
	First7Ghost=memoryarrangeTwoDimension(n,q);
	First8Ghost=memoryarrangeTwoDimension(n,q);
	First9Ghost=memoryarrangeTwoDimension(n,q);
	First10Ghost=memoryarrangeTwoDimension(n,q);
	First11Ghost=memoryarrangeTwoDimension(n,q);
	First12Ghost=memoryarrangeTwoDimension(n,q);
	First13Ghost=memoryarrangeTwoDimension(n,q);
	First14Ghost=memoryarrangeTwoDimension(n,q);
	First15Ghost=memoryarrangeTwoDimension(n,q);
	First16Ghost=memoryarrangeTwoDimension(n,q);
	First17Ghost=memoryarrangeTwoDimension(n,q);
	First18Ghost=memoryarrangeTwoDimension(n,q);
}

void LastlineArrangeVelocityGhost(int n, int q){
	Last0Ghost=memoryarrangeTwoDimension(n,q);
	Last1Ghost=memoryarrangeTwoDimension(n,q);
	Last2Ghost=memoryarrangeTwoDimension(n,q);
	Last3Ghost=memoryarrangeTwoDimension(n,q);
	Last4Ghost=memoryarrangeTwoDimension(n,q);
	Last5Ghost=memoryarrangeTwoDimension(n,q);
	Last6Ghost=memoryarrangeTwoDimension(n,q);
	Last7Ghost=memoryarrangeTwoDimension(n,q);
	Last8Ghost=memoryarrangeTwoDimension(n,q);
	Last9Ghost=memoryarrangeTwoDimension(n,q);
	Last10Ghost=memoryarrangeTwoDimension(n,q);
	Last11Ghost=memoryarrangeTwoDimension(n,q);
	Last12Ghost=memoryarrangeTwoDimension(n,q);
	Last13Ghost=memoryarrangeTwoDimension(n,q);
	Last14Ghost=memoryarrangeTwoDimension(n,q);
	Last15Ghost=memoryarrangeTwoDimension(n,q);
	Last16Ghost=memoryarrangeTwoDimension(n,q);
	Last17Ghost=memoryarrangeTwoDimension(n,q);
	Last18Ghost=memoryarrangeTwoDimension(n,q);
}

double lengthSurfacelineCalculation(int layerNumber, int initial, int m, int finish, int n, double ***uz){
	int i;
	int j;
	int k;
	double av=0.0;
//#pragma omp parallel for reduction(+: av)  
	for (i=initial;i<m;i++){
		for (j=finish;j<n;j++){
			av=av+uz[i][j][layerNumber];
		}
	}
	av=av/(m-initial)/(n-finish);
	return av;
}

void SurfacelineSave(int layerNumber1,int m,int n,int q,double **Save0,double **Save1,double **Save2,double **Save3,double **Save4,double **Save5,double **Save6,
	double **Save7,double **Save8,double **Save9,double **Save10,double **Save11,double **Save12,double **Save13,double **Save14,double **Save15,double **Save16,double **Save17,double **Save18,double **SaveRho,double **SaveUx,double **SaveUy,double **SaveUz,double **SavePressure1,double **SavePressure2,
	double ***In0,double ***In1,double ***In2,double ***In3,double ***In4,double ***In5,double ***In6,double ***In7,double ***In8,double ***In9,
	double ***In10,double ***In11,double ***In12,double ***In13,double ***In14,double ***In15,double ***In16,double ***In17,double ***In18,double ***Rho,double ***ux,double ***uy,double ***uz,double ***pressure){
		int i;
		int j;
		int k;
#pragma omp parallel private(i,j)
		{
#pragma omp for schedule(dynamic)
			for (i=0;i<m;i++){
				for (j=0;j<n;j++){
					Save0[i][j]=In0[i][j][layerNumber1];
					Save1[i][j]=In1[i][j][layerNumber1];
					Save2[i][j]=In2[i][j][layerNumber1];
					Save3[i][j]=In3[i][j][layerNumber1];
					Save4[i][j]=In4[i][j][layerNumber1];
					Save5[i][j]=In5[i][j][layerNumber1];
					Save6[i][j]=In6[i][j][layerNumber1];
					Save7[i][j]=In7[i][j][layerNumber1];
					Save8[i][j]=In8[i][j][layerNumber1];
					Save9[i][j]=In9[i][j][layerNumber1];
					Save10[i][j]=In10[i][j][layerNumber1];
					Save11[i][j]=In11[i][j][layerNumber1];
					Save12[i][j]=In12[i][j][layerNumber1];
					Save13[i][j]=In13[i][j][layerNumber1];
					Save14[i][j]=In14[i][j][layerNumber1];
					Save15[i][j]=In15[i][j][layerNumber1];
					Save16[i][j]=In16[i][j][layerNumber1];
					Save17[i][j]=In17[i][j][layerNumber1];
					Save18[i][j]=In18[i][j][layerNumber1];
					SaveRho[i][j]=Rho[i][j][layerNumber1+1];
					SaveUx[i][j]=ux[i][j][layerNumber1];
					SaveUy[i][j]=uy[i][j][layerNumber1];
					SaveUz[i][j]=uz[i][j][layerNumber1];
					SavePressure2[i][j]=pressure[i][j][layerNumber1];
					SavePressure1[i][j]=pressure[i][j][layerNumber1+1];
				}
			}
		}
}

void LastsurfaceLineMemoryFreeTwoPhase(int m){
	freememoryTwoDimension(Last0,m);
	freememoryTwoDimension(Last1,m);
	freememoryTwoDimension(Last2,m);
	freememoryTwoDimension(Last3,m);
	freememoryTwoDimension(Last4,m);
	freememoryTwoDimension(Last5,m);
	freememoryTwoDimension(Last6,m);
	freememoryTwoDimension(Last7,m);
	freememoryTwoDimension(Last8,m);
	freememoryTwoDimension(Last9,m);
	freememoryTwoDimension(Last10,m);
	freememoryTwoDimension(Last11,m);
	freememoryTwoDimension(Last12,m);
	freememoryTwoDimension(Last13,m);
	freememoryTwoDimension(Last14,m);
	freememoryTwoDimension(Last15,m);
	freememoryTwoDimension(Last16,m);
	freememoryTwoDimension(Last17,m);
	freememoryTwoDimension(Last18,m);
	freememoryTwoDimension(LastFx,m);
	freememoryTwoDimension(LastFy,m);
	freememoryTwoDimension(LastFz,m);
	freememoryTwoDimension(LastRho,m);
	freememoryTwoDimension(LastUx,m);
	freememoryTwoDimension(LastUy,m);
	freememoryTwoDimension(LastUz,m);
	freememoryTwoDimension(LastPressure1,m);
	freememoryTwoDimension(LastPressure2,m);
}

void FirstsurfaceLineMemoryFreeTwoPhase(int n){
	freememoryTwoDimension(First0,n);
	freememoryTwoDimension(First1,n);
	freememoryTwoDimension(First2,n);
	freememoryTwoDimension(First3,n);
	freememoryTwoDimension(First4,n);
	freememoryTwoDimension(First5,n);
	freememoryTwoDimension(First6,n);
	freememoryTwoDimension(First7,n);
	freememoryTwoDimension(First8,n);
	freememoryTwoDimension(First9,n);
	freememoryTwoDimension(First10,n);
	freememoryTwoDimension(First11,n);
	freememoryTwoDimension(First12,n);
	freememoryTwoDimension(First13,n);
	freememoryTwoDimension(First14,n);
	freememoryTwoDimension(First15,n);
	freememoryTwoDimension(First16,n);
	freememoryTwoDimension(First17,n);
	freememoryTwoDimension(First18,n);
	freememoryTwoDimension(FirstFx,n);
	freememoryTwoDimension(FirstFy,n);
	freememoryTwoDimension(FirstFz,n);
}

void surfaceLineMemoryFreeConcentration(int n){
	freememoryTwoDimension(Last0Concentration,n);
	freememoryTwoDimension(Last1Concentration,n);
	freememoryTwoDimension(Last2Concentration,n);
	freememoryTwoDimension(Last3Concentration,n);
	freememoryTwoDimension(Last4Concentration,n);
	freememoryTwoDimension(Last5Concentration,n);
	freememoryTwoDimension(Last6Concentration,n);
	freememoryTwoDimension(Last7Concentration,n);
	freememoryTwoDimension(Last8Concentration,n);
	freememoryTwoDimension(Last9Concentration,n);
	freememoryTwoDimension(Last10Concentration,n);
	freememoryTwoDimension(Last11Concentration,n);
	freememoryTwoDimension(Last12Concentration,n);
	freememoryTwoDimension(Last13Concentration,n);
	freememoryTwoDimension(Last14Concentration,n);
	freememoryTwoDimension(Last15Concentration,n);
	freememoryTwoDimension(Last16Concentration,n);
	freememoryTwoDimension(Last17Concentration,n);
	freememoryTwoDimension(Last18Concentration,n);
}

void surfaceLineMemoryFreeVelocityGhost(int n){
	freememoryTwoDimension(Last0Ghost,n);
	freememoryTwoDimension(Last1Ghost,n);
	freememoryTwoDimension(Last2Ghost,n);
	freememoryTwoDimension(Last3Ghost,n);
	freememoryTwoDimension(Last4Ghost,n);
	freememoryTwoDimension(Last5Ghost,n);
	freememoryTwoDimension(Last6Ghost,n);
	freememoryTwoDimension(Last7Ghost,n);
	freememoryTwoDimension(Last8Ghost,n);
	freememoryTwoDimension(Last9Ghost,n);
	freememoryTwoDimension(Last10Ghost,n);
	freememoryTwoDimension(Last11Ghost,n);
	freememoryTwoDimension(Last12Ghost,n);
	freememoryTwoDimension(Last13Ghost,n);
	freememoryTwoDimension(Last14Ghost,n);
	freememoryTwoDimension(Last15Ghost,n);
	freememoryTwoDimension(Last16Ghost,n);
	freememoryTwoDimension(Last17Ghost,n);
	freememoryTwoDimension(Last18Ghost,n);

	freememoryTwoDimension(First0Ghost,n);
	freememoryTwoDimension(First1Ghost,n);
	freememoryTwoDimension(First2Ghost,n);
	freememoryTwoDimension(First3Ghost,n);
	freememoryTwoDimension(First4Ghost,n);
	freememoryTwoDimension(First5Ghost,n);
	freememoryTwoDimension(First6Ghost,n);
	freememoryTwoDimension(First7Ghost,n);
	freememoryTwoDimension(First8Ghost,n);
	freememoryTwoDimension(First9Ghost,n);
	freememoryTwoDimension(First10Ghost,n);
	freememoryTwoDimension(First11Ghost,n);
	freememoryTwoDimension(First12Ghost,n);
	freememoryTwoDimension(First13Ghost,n);
	freememoryTwoDimension(First14Ghost,n);
	freememoryTwoDimension(First15Ghost,n);
	freememoryTwoDimension(First16Ghost,n);
	freememoryTwoDimension(First17Ghost,n);
	freememoryTwoDimension(First18Ghost,n);
}

void lengthFirstConstantVelocityBoundaryTwoPhase(int nInitial,int n,int qInitial,int q,double InletUx,double InletUy,double InletUz,double solid,double b, double R, double T, double a){
	int i;
	int j;
	int k;
	double bRho4;
	double BRho4;
#pragma omp parallel private(j,k,bRho4,BRho4)
	{
#pragma omp for schedule(dynamic)
		for (j=nInitial;j<n;j++){
			for (k=qInitial;k<q;k++){
				if (domain[1][j][k]==0){
					rho[1][j][k]=((fIn0[1][j][k]+fIn2[1][j][k]+fIn3[1][j][k]+fIn5[1][j][k]+fIn6[1][j][k]+fIn9[1][j][k]+fIn12[1][j][k]+fIn15[1][j][k]+fIn18[1][j][k])
						+2.0*(fIn4[1][j][k]+fIn10[1][j][k]+fIn11[1][j][k]+fIn16[1][j][k]+fIn17[1][j][k]))/(1.0-InletUx);
					ux[1][j][k]=InletUx;
					uy[1][j][k]=InletUy;
					uz[1][j][k]=InletUz;
					bRho4=b*rho[1][j][k]/4.0;
					BRho4=1.0-bRho4;
					//pressure[1][j][k]=rho[1][j][k]*R*T/(1-b*rho[1][j][k])-a*alpha*rho[1][j][k]*rho[1][j][k]/(1+2*b*rho[1][j][k]-b*b*rho[1][j][k]*rho[1][j][k]);
					pressure[1][j][k]=rho[1][j][k]*R*T*(1+bRho4+bRho4*bRho4-bRho4*bRho4*bRho4)/(BRho4*BRho4*BRho4)-a*rho[1][j][k]*rho[1][j][k];
					fIn1[1][j][k]=fIn4[1][j][k]+1.0/3*rho[1][j][k]*InletUx;
					fIn7[1][j][k]=fIn10[1][j][k]+1.0/6*rho[1][j][k]*InletUx;
					fIn8[1][j][k]=fIn11[1][j][k]+1.0/6*rho[1][j][k]*InletUx;
					fIn13[1][j][k]=fIn16[1][j][k]+1.0/6*rho[1][j][k]*InletUx;
					fIn14[1][j][k]=fIn17[1][j][k]+1.0/6*rho[1][j][k]*InletUx;
				}
				if (domain[1][j][k]!=0){
					rho[1][j][k]=solid;
					ux[1][j][k]=0.0;
					uy[1][j][k]=0.0;
					uz[1][j][k]=0.0;
				}
				rho[0][j][k]=rho[1][j][k];
				ux[0][j][k]=ux[1][j][k];
				uy[0][j][k]=uy[1][j][k];
				uz[0][j][k]=uz[1][j][k];
				pressure[0][j][k]=pressure[1][j][k];
			}
		}
	}
}

void HeightFirstConstantVelocityBoundaryTwoPhase(int mInitial,int m,int nInitial,int n,double InletSideUz,double InletSideUx,double InletSideUy,double solid,double b, double R, double T, double a){
	int i;
	int j;
	int k;
	double bRho4;
	double BRho4;
#pragma omp parallel private(i,j,bRho4,BRho4)
	{
#pragma omp for schedule(dynamic)
		for (i=mInitial;i<m;i++){
			for (j=nInitial;j<n;j++){
				if (domain[i][j][1]==0){
					rho[i][j][1]=((fIn0[i][j][1]+fIn1[i][j][1]+fIn2[i][j][1]+fIn4[i][j][1]+fIn5[i][j][1]+fIn7[i][j][1]+fIn10[i][j][1]+fIn13[i][j][1]+fIn16[i][j][1])
						+2.0*(fIn6[i][j][1]+fIn11[i][j][1]+fIn12[i][j][1]+fIn14[i][j][1]+fIn15[i][j][1]))/(1.0-InletSideUz);
					ux[i][j][1]=InletSideUx;
					uy[i][j][1]=InletSideUy;
					uz[i][j][1]=InletSideUz;
					fIn3[i][j][1]=fIn6[i][j][1]+1.0/3*rho[i][j][1]*InletSideUz;
					fIn8[i][j][1]=fIn11[i][j][1]+1.0/6*rho[i][j][1]*(InletSideUx+InletSideUz)-0.5*(fIn1[i][j][1]+fIn7[i][j][1]+fIn13[i][j][1]-(fIn4[i][j][1]+fIn10[i][j][1]+fIn16[i][j][1]))+1.0/3*rho[i][j][1]*InletSideUx;
//					fIn8[i][j][1]=fIn11[i][j][1]+1.0/6*rho[i][j][1]*InletSideUz;
					fIn9[i][j][1]=fIn12[i][j][1]+1.0/6*rho[i][j][1]*(InletSideUy+InletSideUz)-0.5*(fIn2[i][j][1]+fIn7[i][j][1]+fIn16[i][j][1]-(fIn5[i][j][1]+fIn10[i][j][1]+fIn13[i][j][1]))+1.0/3*rho[i][j][1]*InletSideUy;
//					fIn9[i][j][1]=fIn12[i][j][1]+1.0/6*rho[i][j][1]*InletSideUz;
					fIn17[i][j][1]=fIn14[i][j][1]-1.0/6*rho[i][j][1]*(InletSideUx-InletSideUz)+(0.5*(fIn1[i][j][1]+fIn7[i][j][1]+fIn13[i][j][1]-(fIn4[i][j][1]+fIn10[i][j][1]+fIn16[i][j][1]))-1.0/3*rho[i][j][1]*InletSideUx);
//					fIn17[i][j][1]=fIn14[i][j][1]+1.0/6*rho[i][j][1]*InletSideUz;
					fIn15[i][j][1]=fIn18[i][j][1]-1.0/6*rho[i][j][1]*(InletSideUy-InletSideUz)+(0.5*(fIn2[i][j][1]+fIn7[i][j][1]+fIn16[i][j][1]-(fIn5[i][j][1]+fIn10[i][j][1]+fIn13[i][j][1]))-1.0/3*rho[i][j][1]*InletSideUy);
//					fIn15[i][j][1]=fIn18[i][j][1]+1.0/6*rho[i][j][1]*InletSideUz;
					bRho4=b*rho[i][j][1]/4.0;
					BRho4=1.0-bRho4;
					//pressure[i][j][1]=rho[i][j][1]*R*T/(1-b*rho[i][j][1])-a*alpha*rho[i][j][1]*rho[i][j][1]/(1+2*b*rho[i][j][1]-b*b*rho[i][j][1]*rho[i][j][1]);
					pressure[i][j][1]=rho[i][j][1]*R*T*(1+bRho4+bRho4*bRho4-bRho4*bRho4*bRho4)/(BRho4*BRho4*BRho4)-a*rho[i][j][1]*rho[i][j][1];
				}
				if (domain[i][j][1]!=0){
					rho[i][j][1]=solid;
					ux[i][j][1]=0.0;
					uy[i][j][1]=0.0;
					uz[i][j][1]=0.0;
				}
				rho[i][j][0]=rho[i][j][1];
				ux[i][j][0]=ux[i][j][1];
				uy[i][j][0]=uy[i][j][1];
				uz[i][j][0]=uz[i][j][1];
				pressure[i][j][0]=pressure[i][j][1];
			}
		}
	}
}

void HeightFirstConstantPressureBoundaryTwoPhase(int mInitial,int m,int nInitial,int n,double InletDensity,double InletSideUx,double InletSideUy,double solid,double b, double R, double T, double a){
	int i;
	int j;
	int k;
	double bRho4;
	double BRho4;
#pragma omp parallel private(i,j,bRho4,BRho4)
	{
#pragma omp for schedule(dynamic)
		for (i=mInitial;i<m;i++){
			for (j=nInitial;j<n;j++){
				if (domain[i][j][1]==0){
					rho[i][j][1]=InletDensity;
					ux[i][j][1]=InletSideUx;
					uy[i][j][1]=InletSideUy;
					uz[i][j][1]=1.0-((fIn0[i][j][1]+fIn1[i][j][1]+fIn2[i][j][1]+fIn4[i][j][1]+fIn5[i][j][1]+fIn7[i][j][1]+fIn10[i][j][1]+fIn13[i][j][1]+fIn16[i][j][1])
						+2.0*(fIn6[i][j][1]+fIn11[i][j][1]+fIn12[i][j][1]+fIn14[i][j][1]+fIn15[i][j][1]))/InletDensity;
					fIn3[i][j][1]=fIn6[i][j][1]+1.0/3*rho[i][j][1]*uz[i][j][1];
//					fIn8[i][j][1]=fIn11[i][j][1]+1.0/6*rho[i][j][1]*(InletSideUx+uz[i][j][1])-0.5*(fIn1[i][j][1]+fIn7[i][j][1]+fIn13[i][j][1]-(fIn4[i][j][1]+fIn10[i][j][1]+fIn16[i][j][1]))+1.0/3*rho[i][j][1]*InletSideUx;
					fIn8[i][j][1]=fIn11[i][j][1]+1.0/6*rho[i][j][1]*uz[i][j][1];
//					fIn9[i][j][1]=fIn12[i][j][1]+1.0/6*rho[i][j][1]*(InletSideUy+uz[i][j][1])-0.5*(fIn2[i][j][1]+fIn7[i][j][1]+fIn16[i][j][1]-(fIn5[i][j][1]+fIn10[i][j][1]+fIn13[i][j][1]))+1.0/3*rho[i][j][1]*InletSideUy;
					fIn9[i][j][1]=fIn12[i][j][1]+1.0/6*rho[i][j][1]*uz[i][j][1];
//					fIn17[i][j][1]=fIn14[i][j][1]-1.0/6*rho[i][j][1]*(InletSideUx-uz[i][j][1])+(0.5*(fIn1[i][j][1]+fIn7[i][j][1]+fIn13[i][j][1]-(fIn4[i][j][1]+fIn10[i][j][1]+fIn16[i][j][1]))-1.0/3*rho[i][j][1]*InletSideUx);
					fIn17[i][j][1]=fIn14[i][j][1]+1.0/6*rho[i][j][1]*uz[i][j][1];
//					fIn15[i][j][1]=fIn18[i][j][1]-1.0/6*rho[i][j][1]*(InletSideUy-uz[i][j][1])+(0.5*(fIn2[i][j][1]+fIn7[i][j][1]+fIn16[i][j][1]-(fIn5[i][j][1]+fIn10[i][j][1]+fIn13[i][j][1]))-1.0/3*rho[i][j][1]*InletSideUy);
					fIn15[i][j][1]=fIn18[i][j][1]+1.0/6*rho[i][j][1]*uz[i][j][1];
					bRho4=b*rho[i][j][1]/4.0;
					BRho4=1.0-bRho4;
					//pressure[i][j][1]=rho[i][j][1]*R*T/(1-b*rho[i][j][1])-a*alpha*rho[i][j][1]*rho[i][j][1]/(1+2*b*rho[i][j][1]-b*b*rho[i][j][1]*rho[i][j][1]);
					pressure[i][j][1]=rho[i][j][1]*R*T*(1+bRho4+bRho4*bRho4-bRho4*bRho4*bRho4)/(BRho4*BRho4*BRho4)-a*rho[i][j][1]*rho[i][j][1];
				}
				if (domain[i][j][1]!=0){
					rho[i][j][1]=solid;
					ux[i][j][1]=0.0;
					uy[i][j][1]=0.0;
					uz[i][j][1]=0.0;
				}
				rho[i][j][0]=rho[i][j][1];
				ux[i][j][0]=ux[i][j][1];
				uy[i][j][0]=uy[i][j][1];
				uz[i][j][0]=uz[i][j][1];
				pressure[i][j][0]=pressure[i][j][1];
			}
		}
	}
}

void HeightLastConstantVelocityBoundaryTwoPhase(int q,int mInitial,int m,int nInitial,int n,double InletSideUz,double InletSideUx,double InletSideUy,double solid,double b, double R, double T, double a){
	int i;
	int j;
	int k;
	double bRho4;
	double BRho4;
#pragma omp parallel private(i,j,bRho4,BRho4)
	{
#pragma omp for schedule(dynamic)
		for (i=mInitial;i<m;i++){
			for (j=nInitial;j<n;j++){
				if (domain[i][j][q-2]==0){
					rho[i][j][q-2]=((fIn0[i][j][q-2]+fIn1[i][j][q-2]+fIn2[i][j][q-2]+fIn4[i][j][q-2]+fIn5[i][j][q-2]+fIn7[i][j][q-2]+fIn10[i][j][q-2]+fIn13[i][j][q-2]+fIn16[i][j][q-2])
						+2.0*(fIn3[i][j][q-2]+fIn8[i][j][q-2]+fIn9[i][j][q-2]+fIn17[i][j][q-2]+fIn18[i][j][q-2]))/(1.0+InletSideUz);
					ux[i][j][q-2]=InletSideUx;
					uy[i][j][q-2]=InletSideUy;
					uz[i][j][q-2]=InletSideUz;
					fIn6[i][j][q-2]=fIn3[i][j][q-2]-1.0/3*rho[i][j][q-2]*InletSideUz;
					fIn11[i][j][q-2]=fIn8[i][j][q-2]-1.0/6*rho[i][j][q-2]*(InletSideUx+InletSideUz)+0.5*(fIn1[i][j][q-2]+fIn7[i][j][q-2]+fIn13[i][j][q-2]-(fIn4[i][j][q-2]+fIn10[i][j][q-2]+fIn16[i][j][q-2]))-1.0/3*rho[i][j][q-2]*InletSideUx;
					fIn12[i][j][q-2]=fIn9[i][j][q-2]-1.0/6*rho[i][j][q-2]*(InletSideUy+InletSideUz)+0.5*(fIn2[i][j][q-2]+fIn7[i][j][q-2]+fIn16[i][j][q-2]-(fIn5[i][j][q-2]+fIn10[i][j][q-2]+fIn13[i][j][q-2]))-1.0/3*rho[i][j][q-2]*InletSideUy;
					fIn14[i][j][q-2]=fIn17[i][j][q-2]+1.0/6*rho[i][j][q-2]*(InletSideUx-InletSideUz)-(0.5*(fIn1[i][j][q-2]+fIn7[i][j][q-2]+fIn13[i][j][q-2]-(fIn4[i][j][q-2]+fIn10[i][j][q-2]+fIn16[i][j][q-2]))-1.0/3*rho[i][j][q-2]*InletSideUx);
					fIn15[i][j][q-2]=fIn18[i][j][q-2]+1.0/6*rho[i][j][q-2]*(InletSideUy-InletSideUz)-(0.5*(fIn2[i][j][q-2]+fIn7[i][j][q-2]+fIn16[i][j][q-2]-(fIn5[i][j][q-2]+fIn10[i][j][q-2]+fIn13[i][j][q-2]))-1.0/3*rho[i][j][q-2]*InletSideUy);
					bRho4=b*rho[i][j][q-2]/4.0;
					BRho4=1.0-bRho4;
					//pressure[i][j][1]=rho[i][j][1]*R*T/(1-b*rho[i][j][1])-a*alpha*rho[i][j][1]*rho[i][j][1]/(1+2*b*rho[i][j][1]-b*b*rho[i][j][1]*rho[i][j][1]);
					pressure[i][j][q-2]=rho[i][j][q-2]*R*T*(1+bRho4+bRho4*bRho4-bRho4*bRho4*bRho4)/(BRho4*BRho4*BRho4)-a*rho[i][j][q-2]*rho[i][j][q-2];
				}
				if (domain[i][j][q-2]!=0){
					rho[i][j][q-2]=solid;
					ux[i][j][q-2]=0.0;
					uy[i][j][q-2]=0.0;
					uz[i][j][q-2]=0.0;
				}
				rho[i][j][q-1]=rho[i][j][q-2];
				ux[i][j][q-1]=ux[i][j][q-2];
				uy[i][j][q-1]=uy[i][j][q-2];
				uz[i][j][q-1]=uz[i][j][q-2];
				pressure[i][j][q-1]=pressure[i][j][q-2];
			}
		}
	}
}

void HeightLastConstantPressureBoundaryTwoPhase(int q,int mInitial,int m,int nInitial,int n,double outPutDensity,double InletSideUx,double InletSideUy,double solid,double b, double R, double T, double a){
	int i;
	int j;
	int k;
	double bRho4;
	double BRho4;
#pragma omp parallel private(i,j,bRho4,BRho4)
	{
#pragma omp for schedule(dynamic)
		for (i=mInitial;i<m;i++){
			for (j=nInitial;j<n;j++){
				if (domain[i][j][q-2]==0){
					rho[i][j][q-2]=outPutDensity;
					ux[i][j][q-2]=InletSideUx;
					uy[i][j][q-2]=InletSideUy;
					uz[i][j][q-2]=((fIn0[i][j][q-2]+fIn1[i][j][q-2]+fIn2[i][j][q-2]+fIn4[i][j][q-2]+fIn5[i][j][q-2]+fIn7[i][j][q-2]+fIn10[i][j][q-2]+fIn13[i][j][q-2]+fIn16[i][j][q-2])
						+2.0*(fIn3[i][j][q-2]+fIn8[i][j][q-2]+fIn9[i][j][q-2]+fIn17[i][j][q-2]+fIn18[i][j][q-2]))/outPutDensity-1.0;
//					fIn6[i][j][q-2]=fIn3[i][j][q-2]-1.0/3*rho[i][j][q-2]*uz[i][j][q-2];
					fIn6[i][j][q-2]=fIn3[i][j][q-2]-1.0/3*rho[i][j][q-2]*uz[i][j][q-2];
//					fIn11[i][j][q-2]=fIn8[i][j][q-2]-1.0/6*rho[i][j][q-2]*(InletSideUx+uz[i][j][q-2])+0.5*(fIn1[i][j][q-2]+fIn7[i][j][q-2]+fIn13[i][j][q-2]-(fIn4[i][j][q-2]+fIn10[i][j][q-2]+fIn16[i][j][q-2]))-1.0/3*rho[i][j][q-2]*InletSideUx;
					fIn11[i][j][q-2]=fIn8[i][j][q-2]-1.0/6*rho[i][j][q-2]*uz[i][j][q-2];
//					fIn12[i][j][q-2]=fIn9[i][j][q-2]-1.0/6*rho[i][j][q-2]*(InletSideUy+uz[i][j][q-2])+0.5*(fIn2[i][j][q-2]+fIn7[i][j][q-2]+fIn16[i][j][q-2]-(fIn5[i][j][q-2]+fIn10[i][j][q-2]+fIn13[i][j][q-2]))-1.0/3*rho[i][j][q-2]*InletSideUy;
					fIn12[i][j][q-2]=fIn9[i][j][q-2]-1.0/6*rho[i][j][q-2]*uz[i][j][q-2];
//					fIn14[i][j][q-2]=fIn17[i][j][q-2]+1.0/6*rho[i][j][q-2]*(InletSideUx-uz[i][j][q-2])-(0.5*(fIn1[i][j][q-2]+fIn7[i][j][q-2]+fIn13[i][j][q-2]-(fIn4[i][j][q-2]+fIn10[i][j][q-2]+fIn16[i][j][q-2]))-1.0/3*rho[i][j][q-2]*InletSideUx);
					fIn14[i][j][q-2]=fIn17[i][j][q-2]-1.0/6*rho[i][j][q-2]*uz[i][j][q-2];
//					fIn15[i][j][q-2]=fIn18[i][j][q-2]+1.0/6*rho[i][j][q-2]*(InletSideUy-uz[i][j][q-2])-(0.5*(fIn2[i][j][q-2]+fIn7[i][j][q-2]+fIn16[i][j][q-2]-(fIn5[i][j][q-2]+fIn10[i][j][q-2]+fIn13[i][j][q-2]))-1.0/3*rho[i][j][q-2]*InletSideUy);
					fIn15[i][j][q-2]=fIn18[i][j][q-2]-1.0/6*rho[i][j][q-2]*uz[i][j][q-2];
					bRho4=b*rho[i][j][q-2]/4.0;
					BRho4=1.0-bRho4;
					//pressure[i][j][1]=rho[i][j][1]*R*T/(1-b*rho[i][j][1])-a*alpha*rho[i][j][1]*rho[i][j][1]/(1+2*b*rho[i][j][1]-b*b*rho[i][j][1]*rho[i][j][1]);
					pressure[i][j][q-2]=rho[i][j][q-2]*R*T*(1+bRho4+bRho4*bRho4-bRho4*bRho4*bRho4)/(BRho4*BRho4*BRho4)-a*rho[i][j][q-2]*rho[i][j][q-2];
				}
				if (domain[i][j][q-2]!=0){
					rho[i][j][q-2]=solid;
					ux[i][j][q-2]=0.0;
					uy[i][j][q-2]=0.0;
					uz[i][j][q-2]=0.0;
				}
				rho[i][j][q-1]=rho[i][j][q-2];
				ux[i][j][q-1]=ux[i][j][q-2];
				uy[i][j][q-1]=uy[i][j][q-2];
				uz[i][j][q-1]=uz[i][j][q-2];
				pressure[i][j][q-1]=pressure[i][j][q-2];
			}
		}
	}
}

void HeightLastVelocityConvectiveTwoPhase(int q,int mInitial,int m,int nInitial,int n,double uzAV,double solid,double b, double R, double T, double a){
	int i;
	int j;
	double bRho4;
	double BRho4;
	double av2=0.0;
#pragma omp parallel private(i,j,bRho4,BRho4)
	{
#pragma omp for schedule(dynamic)
		for (i=mInitial;i<m;i++){
			for (j=nInitial;j<n;j++){
				if (domain[i][j][q-2]==0){
					fIn0[i][j][q-2]=(Last0[i][j]+uzAV*fIn0[i][j][q-3])/(1.0+uzAV);
					//					fIn0[m-2][j][k]=(Last0[j][k]+ux[m-3][j][k]*fIn0[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn1[i][j][q-2]=(Last1[i][j]+uzAV*fIn1[i][j][q-3])/(1.0+uzAV);
					//					fIn1[m-2][j][k]=(Last1[j][k]+ux[m-3][j][k]*fIn1[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn2[i][j][q-2]=(Last2[i][j]+uzAV*fIn2[i][j][q-3])/(1.0+uzAV);
					//					fIn2[m-2][j][k]=(Last2[j][k]+ux[m-3][j][k]*fIn2[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn3[i][j][q-2]=(Last3[i][j]+uzAV*fIn3[i][j][q-3])/(1.0+uzAV);
					//					fIn3[m-2][j][k]=(Last3[j][k]+ux[m-3][j][k]*fIn3[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn4[i][j][q-2]=(Last4[i][j]+uzAV*fIn4[i][j][q-3])/(1.0+uzAV);
					//					fIn4[m-2][j][k]=(Last4[j][k]+ux[m-3][j][k]*fIn4[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn5[i][j][q-2]=(Last5[i][j]+uzAV*fIn5[i][j][q-3])/(1.0+uzAV);
					//					fIn5[m-2][j][k]=(Last5[j][k]+ux[m-3][j][k]*fIn5[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn6[i][j][q-2]=(Last6[i][j]+uzAV*fIn6[i][j][q-3])/(1.0+uzAV);
					//					fIn6[m-2][j][k]=(Last6[j][k]+ux[m-3][j][k]*fIn6[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn7[i][j][q-2]=(Last7[i][j]+uzAV*fIn7[i][j][q-3])/(1.0+uzAV);
					//					fIn7[m-2][j][k]=(Last7[j][k]+ux[m-3][j][k]*fIn7[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn8[i][j][q-2]=(Last8[i][j]+uzAV*fIn8[i][j][q-3])/(1.0+uzAV);
					//					fIn8[m-2][j][k]=(Last8[j][k]+ux[m-3][j][k]*fIn8[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn9[i][j][q-2]=(Last9[i][j]+uzAV*fIn9[i][j][q-3])/(1.0+uzAV);
					//					fIn9[m-2][j][k]=(Last9[j][k]+ux[m-3][j][k]*fIn9[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn10[i][j][q-2]=(Last10[i][j]+uzAV*fIn10[i][j][q-3])/(1.0+uzAV);
					//					fIn10[m-2][j][k]=(Last10[j][k]+ux[m-3][j][k]*fIn10[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn11[i][j][q-2]=(Last11[i][j]+uzAV*fIn11[i][j][q-3])/(1.0+uzAV);
					//					fIn11[m-2][j][k]=(Last11[j][k]+ux[m-3][j][k]*fIn11[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn12[i][j][q-2]=(Last12[i][j]+uzAV*fIn12[i][j][q-3])/(1.0+uzAV);
					//					fIn12[m-2][j][k]=(Last12[j][k]+ux[m-3][j][k]*fIn12[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn13[i][j][q-2]=(Last13[i][j]+uzAV*fIn13[i][j][q-3])/(1.0+uzAV);
					//					fIn13[m-2][j][k]=(Last13[j][k]+ux[m-3][j][k]*fIn13[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn14[i][j][q-2]=(Last14[i][j]+uzAV*fIn14[i][j][q-3])/(1.0+uzAV);
					//					fIn14[m-2][j][k]=(Last14[j][k]+ux[m-3][j][k]*fIn14[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn15[i][j][q-2]=(Last15[i][j]+uzAV*fIn15[i][j][q-3])/(1.0+uzAV);
					//					fIn15[m-2][j][k]=(Last15[j][k]+ux[m-3][j][k]*fIn15[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn16[i][j][q-2]=(Last16[i][j]+uzAV*fIn16[i][j][q-3])/(1.0+uzAV);
					//					fIn16[m-2][j][k]=(Last16[j][k]+ux[m-3][j][k]*fIn16[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn17[i][j][q-2]=(Last17[i][j]+uzAV*fIn17[i][j][q-3])/(1.0+uzAV);
					//					fIn17[m-2][j][k]=(Last17[j][k]+ux[m-3][j][k]*fIn17[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn18[i][j][q-2]=(Last18[i][j]+uzAV*fIn18[i][j][q-3])/(1.0+uzAV);
					//					fIn18[m-2][j][k]=(Last18[j][k]+ux[m-3][j][k]*fIn18[m-3][j][k])/(1.0+ux[m-3][j][k]);
					rho[i][j][q-2]=fIn0[i][j][q-2]+fIn1[i][j][q-2]+fIn2[i][j][q-2]+fIn3[i][j][q-2]+fIn4[i][j][q-2]+fIn5[i][j][q-2]+fIn6[i][j][q-2]+fIn7[i][j][q-2]+fIn8[i][j][q-2]
					+fIn9[i][j][q-2]+fIn10[i][j][q-2]+fIn11[i][j][q-2]+fIn12[i][j][q-2]+fIn13[i][j][q-2]+fIn14[i][j][q-2]+fIn15[i][j][q-2]+fIn16[i][j][q-2]+fIn17[i][j][q-2]+fIn18[i][j][q-2];
//					ux[m-2][j][k]=(cxNS[0]*fIn0[m-2][j][k]+cxNS[1]*fIn1[m-2][j][k]+cxNS[2]*fIn2[m-2][j][k]+cxNS[3]*fIn3[m-2][j][k]+cxNS[4]*fIn4[m-2][j][k]+cxNS[5]*fIn5[m-2][j][k]+cxNS[6]*fIn6[m-2][j][k]+cxNS[7]*fIn7[m-2][j][k]+cxNS[8]*fIn8[m-2][j][k]\
//						+cxNS[9]*fIn9[m-2][j][k]+cxNS[10]*fIn10[m-2][j][k]+cxNS[11]*fIn11[m-2][j][k]+cxNS[12]*fIn12[m-2][j][k]+cxNS[13]*fIn13[m-2][j][k]+cxNS[14]*fIn14[m-2][j][k]+cxNS[15]*fIn15[m-2][j][k]+cxNS[16]*fIn16[m-2][j][k]+cxNS[17]*fIn17[m-2][j][k]+cxNS[18]*fIn18[m-2][j][k])/rho[m-2][j][k]+Fx[m-2][j][k]/2.0/rho[m-2][j][k];
					ux[i][j][q-2]=(LastUx[i][j]+uzAV*ux[i][j][q-3])/(1.0+uzAV);
//					ux[m-2][j][k]=(LastUx[j][k]+ux[m-3][j][k]*ux[m-3][j][k])/(1.0+ux[m-3][j][k]);
//					uy[m-2][j][k]=(cyNS[0]*fIn0[m-2][j][k]+cyNS[1]*fIn1[m-2][j][k]+cyNS[2]*fIn2[m-2][j][k]+cyNS[3]*fIn3[m-2][j][k]+cyNS[4]*fIn4[m-2][j][k]+cyNS[5]*fIn5[m-2][j][k]+cyNS[6]*fIn6[m-2][j][k]+cyNS[7]*fIn7[m-2][j][k]+cyNS[8]*fIn8[m-2][j][k]\
//						+cyNS[9]*fIn9[m-2][j][k]+cyNS[10]*fIn10[m-2][j][k]+cyNS[11]*fIn11[m-2][j][k]+cyNS[12]*fIn12[m-2][j][k]+cyNS[13]*fIn13[m-2][j][k]+cyNS[14]*fIn14[m-2][j][k]+cyNS[15]*fIn15[m-2][j][k]+cyNS[16]*fIn16[m-2][j][k]+cyNS[17]*fIn17[m-2][j][k]+cyNS[18]*fIn18[m-2][j][k])/rho[m-2][j][k]+Fy[m-2][j][k]/2.0/rho[m-2][j][k];
					uy[i][j][q-2]=(LastUy[i][j]+uzAV*uy[i][j][q-3])/(1.0+uzAV);
//					uy[m-2][j][k]=(LastUy[j][k]+ux[m-3][j][k]*uy[m-3][j][k])/(1.0+ux[m-3][j][k]);
//					uz[m-2][j][k]=(czNS[0]*fIn0[m-2][j][k]+czNS[1]*fIn1[m-2][j][k]+czNS[2]*fIn2[m-2][j][k]+czNS[3]*fIn3[m-2][j][k]+czNS[4]*fIn4[m-2][j][k]+czNS[5]*fIn5[m-2][j][k]+czNS[6]*fIn6[m-2][j][k]+czNS[7]*fIn7[m-2][j][k]+czNS[8]*fIn8[m-2][j][k]\
//						+czNS[9]*fIn9[m-2][j][k]+czNS[10]*fIn10[m-2][j][k]+czNS[11]*fIn11[m-2][j][k]+czNS[12]*fIn12[m-2][j][k]+czNS[13]*fIn13[m-2][j][k]+czNS[14]*fIn14[m-2][j][k]+czNS[15]*fIn15[m-2][j][k]+czNS[16]*fIn16[m-2][j][k]+czNS[17]*fIn17[m-2][j][k]+czNS[18]*fIn18[m-2][j][k])/rho[m-2][j][k]+Fz[m-2][j][k]/2.0/rho[m-2][j][k];
					uz[i][j][q-2]=(LastUz[i][j]+uzAV*uz[i][j][q-3])/(1.0+uzAV);
//					uz[m-2][j][k]=(LastUz[j][k]+ux[m-3][j][k]*uz[m-3][j][k])/(1.0+ux[m-3][j][k]);
					bRho4=b*rho[i][j][q-2]/4.0;
					BRho4=1.0-bRho4;
					//pressure[m-1][j][k]=rho[m-1][j][k]*R*T/(1-b*rho[m-1][j][k])-a*alpha*rho[m-1][j][k]*rho[m-1][j][k]/(1+2*b*rho[m-1][j][k]-b*b*rho[m-1][j][k]*rho[m-1][j][k]);
					//pressure[m-2][j][k]=rho[m-2][j][k]*R*T*(1+bRho4+bRho4*bRho4-bRho4*bRho4*bRho4)/(BRho4*BRho4*BRho4)-a*rho[m-2][j][k]*rho[m-2][j][k];
					pressure[i][j][q-2]=(LastPressure2[i][j]+uzAV*pressure[i][j][q-3])/(1.0+uzAV);
					//pressure[m-2][j][k]=(LastPressure2[j][k]+ux[m-3][j][k]*pressure[m-3][j][k])/(1.0+ux[m-3][j][k]);
				}
				if (domain[i][j][q-2]!=0){
					rho[i][j][q-2]=solid;
					ux[i][j][q-2]=0.0;
					uy[i][j][q-2]=0.0;
					uz[i][j][q-2]=0.0;
				}
			}
		}
	}
//#pragma omp parallel for reduction(+: av2)  
	for (i=mInitial;i<m;i++){
		for (j=nInitial;j<n;j++){
			av2=av2+uz[i][j][q-2];
		}
	}
	av2=av2/(m-mInitial)/(n-nInitial);
#pragma omp parallel private(i,j,bRho4,BRho4)
	{
#pragma omp for schedule(dynamic)
		for (i=mInitial;i<m;i++){
			for (j=nInitial;j<n;j++){
				if (domain[i][j][q-1]==0){
					rho[i][j][q-1]=(LastRho[i][j]+av2*rho[i][j][q-2])/(1.0+av2);
//					rho[m-1][j][k]=(LastRho[j][k]+ux[m-3][j][k]*rho[m-2][j][k])/(1.0+ux[m-3][j][k]);
					ux[i][j][q-1]=ux[i][j][q-2];
					uy[i][j][q-1]=uy[i][j][q-2];
					uz[i][j][q-1]=uz[i][j][q-2];
					bRho4=b*rho[i][j][q-1]/4.0;
					BRho4=1.0-bRho4;
//					pressure[m-1][j][k]=rho[m-1][j][k]*R*T*(1+bRho4+bRho4*bRho4-bRho4*bRho4*bRho4)/(BRho4*BRho4*BRho4)-a*rho[m-1][j][k]*rho[m-1][j][k];
					pressure[i][j][q-1]=(LastPressure1[i][j]+av2*pressure[i][j][q-2])/(1.0+av2);
//					pressure[m-1][j][k]=(LastPressure1[j][k]+ux[m-3][j][k]*pressure[m-2][j][k])/(1.0+ux[m-3][j][k]);
				}
				if (domain[i][j][q-1]!=0){
					rho[i][j][q-1]=solid;
					ux[i][j][q-1]=0.0;
					uy[i][j][q-1]=0.0;
					uz[i][j][q-1]=0.0;
				}
			}
		}
	}
}

void HeightLastNeumannVelocityTwoPhase(int q,int mInitial,int m,int nInitial,int n,double solid){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=mInitial;i<m;i++){
			for (j=nInitial;j<n;j++){
				if (domain[i][j][q-2]==0){
					fIn0[i][j][q-2]=fIn0[i][j][q-3];
					fIn1[i][j][q-2]=fIn1[i][j][q-3];
					fIn2[i][j][q-2]=fIn2[i][j][q-3];
					fIn3[i][j][q-2]=fIn3[i][j][q-3];
					fIn4[i][j][q-2]=fIn4[i][j][q-3];
					fIn5[i][j][q-2]=fIn5[i][j][q-3];
					fIn6[i][j][q-2]=fIn6[i][j][q-3];
					fIn7[i][j][q-2]=fIn7[i][j][q-3];
					fIn8[i][j][q-2]=fIn8[i][j][q-3];
					fIn9[i][j][q-2]=fIn9[i][j][q-3];
					fIn10[i][j][q-2]=fIn10[i][j][q-3];
					fIn11[i][j][q-2]=fIn11[i][j][q-3];
					fIn12[i][j][q-2]=fIn12[i][j][q-3];
					fIn13[i][j][q-2]=fIn13[i][j][q-3];
					fIn14[i][j][q-2]=fIn14[i][j][q-3];
					fIn15[i][j][q-2]=fIn15[i][j][q-3];
					fIn16[i][j][q-2]=fIn16[i][j][q-3];
					fIn17[i][j][q-2]=fIn17[i][j][q-3];
					fIn18[i][j][q-2]=fIn18[i][j][q-3];
					rho[i][j][q-2]=fIn0[i][j][q-2]+fIn1[i][j][q-2]+fIn2[i][j][q-2]+fIn3[i][j][q-2]+fIn4[i][j][q-2]+fIn5[i][j][q-2]+fIn6[i][j][q-2]+fIn7[i][j][q-2]+fIn8[i][j][q-2]
					+fIn9[i][j][q-2]+fIn10[i][j][q-2]+fIn11[i][j][q-2]+fIn12[i][j][q-2]+fIn13[i][j][q-2]+fIn14[i][j][q-2]+fIn15[i][j][q-2]+fIn16[i][j][q-2]+fIn17[i][j][q-2]+fIn18[i][j][q-2];
					ux[i][j][q-2]=(cxNS[0]*fIn0[i][j][q-2]+cxNS[1]*fIn1[i][j][q-2]+cxNS[2]*fIn2[i][j][q-2]+cxNS[3]*fIn3[i][j][q-2]+cxNS[4]*fIn4[i][j][q-2]+cxNS[5]*fIn5[i][j][q-2]+cxNS[6]*fIn6[i][j][q-2]+cxNS[7]*fIn7[i][j][q-2]+cxNS[8]*fIn8[i][j][q-2]\
						+cxNS[9]*fIn9[i][j][q-2]+cxNS[10]*fIn10[i][j][q-2]+cxNS[11]*fIn11[i][j][q-2]+cxNS[12]*fIn12[i][j][q-2]+cxNS[13]*fIn13[i][j][q-2]+cxNS[14]*fIn14[i][j][q-2]+cxNS[15]*fIn15[i][j][q-2]+cxNS[16]*fIn16[i][j][q-2]+cxNS[17]*fIn17[i][j][q-2]+cxNS[18]*fIn18[i][j][q-2])/rho[i][j][q-2]+Fx[i][j][q-2]/2.0/rho[i][j][q-2];
					uy[i][j][q-2]=(cyNS[0]*fIn0[i][j][q-2]+cyNS[1]*fIn1[i][j][q-2]+cyNS[2]*fIn2[i][j][q-2]+cyNS[3]*fIn3[i][j][q-2]+cyNS[4]*fIn4[i][j][q-2]+cyNS[5]*fIn5[i][j][q-2]+cyNS[6]*fIn6[i][j][q-2]+cyNS[7]*fIn7[i][j][q-2]+cyNS[8]*fIn8[i][j][q-2]\
						+cyNS[9]*fIn9[i][j][q-2]+cyNS[10]*fIn10[i][j][q-2]+cyNS[11]*fIn11[i][j][q-2]+cyNS[12]*fIn12[i][j][q-2]+cyNS[13]*fIn13[i][j][q-2]+cyNS[14]*fIn14[i][j][q-2]+cyNS[15]*fIn15[i][j][q-2]+cyNS[16]*fIn16[i][j][q-2]+cyNS[17]*fIn17[i][j][q-2]+cyNS[18]*fIn18[i][j][q-2])/rho[i][j][q-2]+Fy[i][j][q-2]/2.0/rho[i][j][q-2];
					uz[i][j][q-2]=(czNS[0]*fIn0[i][j][q-2]+czNS[1]*fIn1[i][j][q-2]+czNS[2]*fIn2[i][j][q-2]+czNS[3]*fIn3[i][j][q-2]+czNS[4]*fIn4[i][j][q-2]+czNS[5]*fIn5[i][j][q-2]+czNS[6]*fIn6[i][j][q-2]+czNS[7]*fIn7[i][j][q-2]+czNS[8]*fIn8[i][j][q-2]\
						+czNS[9]*fIn9[i][j][q-2]+czNS[10]*fIn10[i][j][q-2]+czNS[11]*fIn11[i][j][q-2]+czNS[12]*fIn12[i][j][q-2]+czNS[13]*fIn13[i][j][q-2]+czNS[14]*fIn14[i][j][q-2]+czNS[15]*fIn15[i][j][q-2]+czNS[16]*fIn16[i][j][q-2]+czNS[17]*fIn17[i][j][q-2]+czNS[18]*fIn18[i][j][q-2])/rho[i][j][q-2]+Fz[i][j][q-2]/2.0/rho[i][j][q-2];
					pressure[i][j][q-2]=pressure[i][j][q-3];
				}
				if (domain[i][j][q-2]!=0){
					rho[i][j][q-2]=solid;
					ux[i][j][q-2]=0.0;
					uy[i][j][q-2]=0.0;
					uz[i][j][q-2]=0.0;
				}
			}
		}
	}
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=mInitial;i<m;i++){
			for (j=nInitial;j<n;j++){
				if (domain[i][j][q-1]==0){
					rho[i][j][q-1]=rho[i][j][q-2];
					ux[i][j][q-1]=ux[i][j][q-2];
					uy[i][j][q-1]=uy[i][j][q-2];
					uz[i][j][q-1]=uz[i][j][q-2];
					pressure[i][j][q-1]=pressure[i][j][q-2];
				}
				if (domain[i][j][q-1]!=0){
					rho[i][j][q-1]=solid;
					ux[i][j][q-1]=0.0;
					uy[i][j][q-1]=0.0;
					uz[i][j][q-1]=0.0;
				}
			}
		}
	}
}

void lengthFirstForceVelocityConvectiveTwoPhase(int nInitial,int n,int qInitial,int q,double uxAV,double G){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=nInitial;j<n;j++){
			for (k=qInitial;k<q;k++){
				if (domain[0][j][k]==0){
					Fx[0][j][k]=(FirstFx[j][k]+uxAV*(Fx[1][j][k]-G))/(1.0+uxAV)+G;
					//Fx[0][j][k]=(FirstFx[j][k]+ux[1][j][k]*(Fx[1][j][k]-G))/(1.0+ux[1][j][k])+G;
					Fy[0][j][k]=(FirstFy[j][k]+uxAV*Fy[1][j][k])/(1.0+uxAV);
					//Fy[0][j][k]=(FirstFy[j][k]+ux[1][j][k]*Fy[1][j][k])/(1.0+ux[1][j][k]);
					Fz[0][j][k]=(FirstFz[j][k]+uxAV*Fz[1][j][k])/(1.0+uxAV);
					//Fz[0][j][k]=(FirstFz[j][k]+ux[1][j][k]*Fz[1][j][k])/(1.0+ux[1][j][k]);
				}
			}
		}
	}
}

void lengthFirstVelocityConvectiveTwoPhase(int nInitial,int n,int qInitial,int q,double uxAV,double solid,double b, double R, double T, double a){
	int j;
	int k;
	double bRho4;
	double BRho4;
#pragma omp parallel private(j,k,bRho4,BRho4)
	{
#pragma omp for schedule(dynamic)
		for (j=nInitial;j<n;j++){
			for (k=qInitial;k<q;k++){
				if (domain[0][j][k]==0){
					fIn0[0][j][k]=(First0[j][k]+uxAV*fIn0[1][j][k])/(1.0+uxAV);
					fIn1[0][j][k]=(First1[j][k]+uxAV*fIn1[1][j][k])/(1.0+uxAV);
					fIn2[0][j][k]=(First2[j][k]+uxAV*fIn2[1][j][k])/(1.0+uxAV);
					fIn3[0][j][k]=(First3[j][k]+uxAV*fIn3[1][j][k])/(1.0+uxAV);
					fIn4[0][j][k]=(First4[j][k]+uxAV*fIn4[1][j][k])/(1.0+uxAV);
					fIn5[0][j][k]=(First5[j][k]+uxAV*fIn5[1][j][k])/(1.0+uxAV);
					fIn6[0][j][k]=(First6[j][k]+uxAV*fIn6[1][j][k])/(1.0+uxAV);
					fIn7[0][j][k]=(First7[j][k]+uxAV*fIn7[1][j][k])/(1.0+uxAV);
					fIn8[0][j][k]=(First8[j][k]+uxAV*fIn8[1][j][k])/(1.0+uxAV);
					fIn9[0][j][k]=(First9[j][k]+uxAV*fIn9[1][j][k])/(1.0+uxAV);
					fIn10[0][j][k]=(First10[j][k]+uxAV*fIn10[1][j][k])/(1.0+uxAV);
					fIn11[0][j][k]=(First11[j][k]+uxAV*fIn11[1][j][k])/(1.0+uxAV);
					fIn12[0][j][k]=(First12[j][k]+uxAV*fIn12[1][j][k])/(1.0+uxAV);
					fIn13[0][j][k]=(First13[j][k]+uxAV*fIn13[1][j][k])/(1.0+uxAV);
					fIn14[0][j][k]=(First14[j][k]+uxAV*fIn14[1][j][k])/(1.0+uxAV);
					fIn15[0][j][k]=(First15[j][k]+uxAV*fIn15[1][j][k])/(1.0+uxAV);
					fIn16[0][j][k]=(First16[j][k]+uxAV*fIn16[1][j][k])/(1.0+uxAV);
					fIn17[0][j][k]=(First17[j][k]+uxAV*fIn17[1][j][k])/(1.0+uxAV);
					fIn18[0][j][k]=(First18[j][k]+uxAV*fIn18[1][j][k])/(1.0+uxAV);
					rho[0][j][k]=fIn0[0][j][k]+fIn1[0][j][k]+fIn2[0][j][k]+fIn3[0][j][k]+fIn4[0][j][k]+fIn5[0][j][k]+fIn6[0][j][k]+fIn7[0][j][k]+fIn8[0][j][k]
					+fIn9[0][j][k]+fIn10[0][j][k]+fIn11[0][j][k]+fIn12[0][j][k]+fIn13[0][j][k]+fIn14[0][j][k]+fIn15[0][j][k]+fIn16[0][j][k]+fIn17[0][j][k]+fIn18[0][j][k];
					ux[0][j][k]=(cxNS[0]*fIn0[0][j][k]+cxNS[1]*fIn1[0][j][k]+cxNS[2]*fIn2[0][j][k]+cxNS[3]*fIn3[0][j][k]+cxNS[4]*fIn4[0][j][k]+cxNS[5]*fIn5[0][j][k]+cxNS[6]*fIn6[0][j][k]+cxNS[7]*fIn7[0][j][k]+cxNS[8]*fIn8[0][j][k]\
						+cxNS[9]*fIn9[0][j][k]+cxNS[10]*fIn10[0][j][k]+cxNS[11]*fIn11[0][j][k]+cxNS[12]*fIn12[0][j][k]+cxNS[13]*fIn13[0][j][k]+cxNS[14]*fIn14[0][j][k]+cxNS[15]*fIn15[0][j][k]+cxNS[16]*fIn16[0][j][k]+cxNS[17]*fIn17[0][j][k]+cxNS[18]*fIn18[0][j][k])/rho[0][j][k]+Fx[0][j][k]/2.0/rho[0][j][k];
					uy[0][j][k]=(cyNS[0]*fIn0[0][j][k]+cyNS[1]*fIn1[0][j][k]+cyNS[2]*fIn2[0][j][k]+cyNS[3]*fIn3[0][j][k]+cyNS[4]*fIn4[0][j][k]+cyNS[5]*fIn5[0][j][k]+cyNS[6]*fIn6[0][j][k]+cyNS[7]*fIn7[0][j][k]+cyNS[8]*fIn8[0][j][k]\
						+cyNS[9]*fIn9[0][j][k]+cyNS[10]*fIn10[0][j][k]+cyNS[11]*fIn11[0][j][k]+cyNS[12]*fIn12[0][j][k]+cyNS[13]*fIn13[0][j][k]+cyNS[14]*fIn14[0][j][k]+cyNS[15]*fIn15[0][j][k]+cyNS[16]*fIn16[0][j][k]+cyNS[17]*fIn17[0][j][k]+cyNS[18]*fIn18[0][j][k])/rho[0][j][k]+Fy[0][j][k]/2.0/rho[0][j][k];
					uz[0][j][k]=(czNS[0]*fIn0[0][j][k]+czNS[0]*fIn1[0][j][k]+czNS[2]*fIn2[0][j][k]+czNS[3]*fIn3[0][j][k]+czNS[4]*fIn4[0][j][k]+czNS[5]*fIn5[0][j][k]+czNS[6]*fIn6[0][j][k]+czNS[7]*fIn7[0][j][k]+czNS[8]*fIn8[0][j][k]\
						+czNS[9]*fIn9[0][j][k]+czNS[10]*fIn10[0][j][k]+czNS[11]*fIn11[0][j][k]+czNS[12]*fIn12[0][j][k]+czNS[13]*fIn13[0][j][k]+czNS[14]*fIn14[0][j][k]+czNS[15]*fIn15[0][j][k]+czNS[16]*fIn16[0][j][k]+czNS[17]*fIn17[0][j][k]+czNS[18]*fIn18[0][j][k])/rho[0][j][k]+Fz[0][j][k]/2.0/rho[0][j][k];
					bRho4=b*rho[0][j][k]/4.0;
					BRho4=1.0-bRho4;
					//pressure[0][j][k]=rho[0][j][k]*R*T/(1-b*rho[0][j][k])-a*alpha*rho[0][j][k]*rho[0][j][k]/(1+2*b*rho[0][j][k]-b*b*rho[0][j][k]*rho[0][j][k]);
					pressure[0][j][k]=rho[0][j][k]*R*T*(1+bRho4+bRho4*bRho4-bRho4*bRho4*bRho4)/(BRho4*BRho4*BRho4)-a*rho[0][j][k]*rho[0][j][k];
				}
				if (domain[0][j][k]!=0){
					rho[0][j][k]=solid;
					ux[0][j][k]=0.0;
					uy[0][j][k]=0.0;
					uz[0][j][k]=0.0;
				}
			}
		}
	}
}

void lengthFirstVelocityConvectiveGhost(int nInitial,int n,int qInitial,int q,double uxAV,double solid,double gasCritical){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=nInitial;j<n;j++){
			for (k=qInitial;k<q;k++){
				if (domain[0][j][k]==0&&rho[0][j][k]<gasCritical){
					fIn0Ghost[0][j][k]=(First0Ghost[j][k]+uxAV*fIn0Ghost[1][j][k])/(1.0+uxAV);
					fIn1Ghost[0][j][k]=(First1Ghost[j][k]+uxAV*fIn1Ghost[1][j][k])/(1.0+uxAV);
					fIn2Ghost[0][j][k]=(First2Ghost[j][k]+uxAV*fIn2Ghost[1][j][k])/(1.0+uxAV);
					fIn3Ghost[0][j][k]=(First3Ghost[j][k]+uxAV*fIn3Ghost[1][j][k])/(1.0+uxAV);
					fIn4Ghost[0][j][k]=(First4Ghost[j][k]+uxAV*fIn4Ghost[1][j][k])/(1.0+uxAV);
					fIn5Ghost[0][j][k]=(First5Ghost[j][k]+uxAV*fIn5Ghost[1][j][k])/(1.0+uxAV);
					fIn6Ghost[0][j][k]=(First6Ghost[j][k]+uxAV*fIn6Ghost[1][j][k])/(1.0+uxAV);
					fIn7Ghost[0][j][k]=(First7Ghost[j][k]+uxAV*fIn7Ghost[1][j][k])/(1.0+uxAV);
					fIn8Ghost[0][j][k]=(First8Ghost[j][k]+uxAV*fIn8Ghost[1][j][k])/(1.0+uxAV);
					fIn9Ghost[0][j][k]=(First9Ghost[j][k]+uxAV*fIn9Ghost[1][j][k])/(1.0+uxAV);
					fIn10Ghost[0][j][k]=(First10Ghost[j][k]+uxAV*fIn10Ghost[1][j][k])/(1.0+uxAV);
					fIn11Ghost[0][j][k]=(First11Ghost[j][k]+uxAV*fIn11Ghost[1][j][k])/(1.0+uxAV);
					fIn12Ghost[0][j][k]=(First12Ghost[j][k]+uxAV*fIn12Ghost[1][j][k])/(1.0+uxAV);
					fIn13Ghost[0][j][k]=(First13Ghost[j][k]+uxAV*fIn13Ghost[1][j][k])/(1.0+uxAV);
					fIn14Ghost[0][j][k]=(First14Ghost[j][k]+uxAV*fIn14Ghost[1][j][k])/(1.0+uxAV);
					fIn15Ghost[0][j][k]=(First15Ghost[j][k]+uxAV*fIn15Ghost[1][j][k])/(1.0+uxAV);
					fIn16Ghost[0][j][k]=(First16Ghost[j][k]+uxAV*fIn16Ghost[1][j][k])/(1.0+uxAV);
					fIn17Ghost[0][j][k]=(First17Ghost[j][k]+uxAV*fIn17Ghost[1][j][k])/(1.0+uxAV);
					fIn18Ghost[0][j][k]=(First18Ghost[j][k]+uxAV*fIn18Ghost[1][j][k])/(1.0+uxAV);
					rhoGhost[0][j][k]=fIn0Ghost[0][j][k]+fIn1Ghost[0][j][k]+fIn2Ghost[0][j][k]+fIn3Ghost[0][j][k]+fIn4Ghost[0][j][k]+fIn5Ghost[0][j][k]+fIn6Ghost[0][j][k]+fIn7Ghost[0][j][k]+fIn8Ghost[0][j][k]
					+fIn9Ghost[0][j][k]+fIn10Ghost[0][j][k]+fIn11Ghost[0][j][k]+fIn12Ghost[0][j][k]+fIn13Ghost[0][j][k]+fIn14Ghost[0][j][k]+fIn15Ghost[0][j][k]+fIn16Ghost[0][j][k]+fIn17Ghost[0][j][k]+fIn18Ghost[0][j][k];
					uxGhost[0][j][k]=(cxNS[0]*fIn0Ghost[0][j][k]+cxNS[1]*fIn1Ghost[0][j][k]+cxNS[2]*fIn2Ghost[0][j][k]+cxNS[3]*fIn3Ghost[0][j][k]+cxNS[4]*fIn4Ghost[0][j][k]+cxNS[5]*fIn5Ghost[0][j][k]+cxNS[6]*fIn6Ghost[0][j][k]+cxNS[7]*fIn7Ghost[0][j][k]+cxNS[8]*fIn8Ghost[0][j][k]\
						+cxNS[9]*fIn9Ghost[0][j][k]+cxNS[10]*fIn10Ghost[0][j][k]+cxNS[11]*fIn11Ghost[0][j][k]+cxNS[12]*fIn12Ghost[0][j][k]+cxNS[13]*fIn13Ghost[0][j][k]+cxNS[14]*fIn14Ghost[0][j][k]+cxNS[15]*fIn15Ghost[0][j][k]+cxNS[16]*fIn16Ghost[0][j][k]+cxNS[17]*fIn17Ghost[0][j][k]+cxNS[18]*fIn18Ghost[0][j][k])/rhoGhost[0][j][k]+bodyForce[0][j][k]/2.0/rhoGhost[0][j][k];
					uyGhost[0][j][k]=(cyNS[0]*fIn0Ghost[0][j][k]+cyNS[1]*fIn1Ghost[0][j][k]+cyNS[2]*fIn2Ghost[0][j][k]+cyNS[3]*fIn3Ghost[0][j][k]+cyNS[4]*fIn4Ghost[0][j][k]+cyNS[5]*fIn5Ghost[0][j][k]+cyNS[6]*fIn6Ghost[0][j][k]+cyNS[7]*fIn7Ghost[0][j][k]+cyNS[8]*fIn8Ghost[0][j][k]\
						+cyNS[9]*fIn9Ghost[0][j][k]+cyNS[10]*fIn10Ghost[0][j][k]+cyNS[11]*fIn11Ghost[0][j][k]+cyNS[12]*fIn12Ghost[0][j][k]+cyNS[13]*fIn13Ghost[0][j][k]+cyNS[14]*fIn14Ghost[0][j][k]+cyNS[15]*fIn15Ghost[0][j][k]+cyNS[16]*fIn16Ghost[0][j][k]+cyNS[17]*fIn17Ghost[0][j][k]+cyNS[18]*fIn18Ghost[0][j][k])/rhoGhost[0][j][k];
					uzGhost[0][j][k]=(czNS[0]*fIn0Ghost[0][j][k]+czNS[0]*fIn1Ghost[0][j][k]+czNS[2]*fIn2Ghost[0][j][k]+czNS[3]*fIn3Ghost[0][j][k]+czNS[4]*fIn4Ghost[0][j][k]+czNS[5]*fIn5Ghost[0][j][k]+czNS[6]*fIn6Ghost[0][j][k]+czNS[7]*fIn7Ghost[0][j][k]+czNS[8]*fIn8Ghost[0][j][k]\
						+czNS[9]*fIn9Ghost[0][j][k]+czNS[10]*fIn10Ghost[0][j][k]+czNS[11]*fIn11Ghost[0][j][k]+czNS[12]*fIn12Ghost[0][j][k]+czNS[13]*fIn13Ghost[0][j][k]+czNS[14]*fIn14Ghost[0][j][k]+czNS[15]*fIn15Ghost[0][j][k]+czNS[16]*fIn16Ghost[0][j][k]+czNS[17]*fIn17Ghost[0][j][k]+czNS[18]*fIn18Ghost[0][j][k])/rhoGhost[0][j][k];
				}
				if (domain[0][j][k]!=0||rho[0][j][k]>=gasCritical){
					rho[0][j][k]=solid;
					ux[0][j][k]=0.0;
					uy[0][j][k]=0.0;
					uz[0][j][k]=0.0;
				}
			}
		}
	}
}

void lengthFirstNeumannInletTwoPhase(int nInitial,int n,int qInitial,int q,double solid){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=nInitial;j<n;j++){
			for (k=qInitial;k<q;k++){
				if (domain[1][j][k]==0){
					fIn0[1][j][k]=fIn0[2][j][k];
					fIn1[1][j][k]=fIn1[2][j][k];
					fIn2[1][j][k]=fIn2[2][j][k];
					fIn3[1][j][k]=fIn3[2][j][k];
					fIn4[1][j][k]=fIn4[2][j][k];
					fIn5[1][j][k]=fIn5[2][j][k];
					fIn6[1][j][k]=fIn6[2][j][k];
					fIn7[1][j][k]=fIn7[2][j][k];
					fIn8[1][j][k]=fIn8[2][j][k];
					fIn9[1][j][k]=fIn9[2][j][k];
					fIn10[1][j][k]=fIn10[2][j][k];
					fIn11[1][j][k]=fIn11[2][j][k];
					fIn12[1][j][k]=fIn12[2][j][k];
					fIn13[1][j][k]=fIn13[2][j][k];
					fIn14[1][j][k]=fIn14[2][j][k];
					fIn15[1][j][k]=fIn15[2][j][k];
					fIn16[1][j][k]=fIn16[2][j][k];
					fIn17[1][j][k]=fIn17[2][j][k];
					fIn18[1][j][k]=fIn18[2][j][k];
					rho[1][j][k]=rho[2][j][k];
					pressure[1][j][k]=pressure[2][j][k];
					ux[1][j][k]=(cxNS[0]*fIn0[1][j][k]+cxNS[1]*fIn1[1][j][k]+cxNS[2]*fIn2[1][j][k]+cxNS[3]*fIn3[1][j][k]+cxNS[4]*fIn4[1][j][k]+cxNS[5]*fIn5[1][j][k]+cxNS[6]*fIn6[1][j][k]+cxNS[7]*fIn7[1][j][k]+cxNS[8]*fIn8[1][j][k]\
						+cxNS[9]*fIn9[1][j][k]+cxNS[10]*fIn10[1][j][k]+cxNS[11]*fIn11[1][j][k]+cxNS[12]*fIn12[1][j][k]+cxNS[13]*fIn13[1][j][k]+cxNS[14]*fIn14[1][j][k]+cxNS[15]*fIn15[1][j][k]+cxNS[16]*fIn16[1][j][k]+cxNS[17]*fIn17[1][j][k]+cxNS[18]*fIn18[1][j][k])/rho[1][j][k]+Fx[1][j][k]/2.0/rho[1][j][k];
					uy[1][j][k]=(cyNS[0]*fIn0[1][j][k]+cyNS[1]*fIn1[1][j][k]+cyNS[2]*fIn2[1][j][k]+cyNS[3]*fIn3[1][j][k]+cyNS[4]*fIn4[1][j][k]+cyNS[5]*fIn5[1][j][k]+cyNS[6]*fIn6[1][j][k]+cyNS[7]*fIn7[1][j][k]+cyNS[8]*fIn8[1][j][k]\
						+cyNS[9]*fIn9[1][j][k]+cyNS[10]*fIn10[1][j][k]+cyNS[11]*fIn11[1][j][k]+cyNS[12]*fIn12[1][j][k]+cyNS[13]*fIn13[1][j][k]+cyNS[14]*fIn14[1][j][k]+cyNS[15]*fIn15[1][j][k]+cyNS[16]*fIn16[1][j][k]+cyNS[17]*fIn17[1][j][k]+cyNS[18]*fIn18[1][j][k])/rho[1][j][k]+Fy[1][j][k]/2.0/rho[1][j][k];
					uz[1][j][k]=(czNS[0]*fIn0[1][j][k]+czNS[1]*fIn1[1][j][k]+czNS[2]*fIn2[1][j][k]+czNS[3]*fIn3[1][j][k]+czNS[4]*fIn4[1][j][k]+czNS[5]*fIn5[1][j][k]+czNS[6]*fIn6[1][j][k]+czNS[7]*fIn7[1][j][k]+czNS[8]*fIn8[1][j][k]\
						+czNS[9]*fIn9[1][j][k]+czNS[10]*fIn10[1][j][k]+czNS[11]*fIn11[1][j][k]+czNS[12]*fIn12[1][j][k]+czNS[13]*fIn13[1][j][k]+czNS[14]*fIn14[1][j][k]+czNS[15]*fIn15[1][j][k]+czNS[16]*fIn16[1][j][k]+czNS[17]*fIn17[1][j][k]+czNS[18]*fIn18[1][j][k])/rho[1][j][k]+Fz[1][j][k]/2.0/rho[1][j][k];
				}
				if (domain[1][j][k]!=0){
					rho[1][j][k]=solid;
					ux[1][j][k]=0.0;
					uy[1][j][k]=0.0;
					uz[1][j][k]=0.0;
				}
			}
		}
	}
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=nInitial;j<n;j++){
			for (k=qInitial;k<q;k++){
				if (domain[0][j][k]==0){
					rho[0][j][k]=rho[1][j][k];
					pressure[0][j][k]=pressure[1][j][k];
					ux[0][j][k]=ux[1][j][k];
					uy[0][j][k]=uy[1][j][k];
					uz[0][j][k]=uz[1][j][k];
				}
				if (domain[0][j][k]!=0){
					rho[0][j][k]=solid;
					ux[0][j][k]=0.0;
					uy[0][j][k]=0.0;
					uz[0][j][k]=0.0;
				}
			}
		}
	}
}

void lengthFirstConstantPressureTwoPhase(int nInitial,int n,int qInitial,int q,double rhoConstant,double InletUy,double InletUz,double solid,double b, double R, double T, double a){
	int j;
	int k;
	double bRho4;
	double BRho4;
#pragma omp parallel private(j,k,bRho4,BRho4)
	{
#pragma omp for schedule(dynamic)
		for (j=nInitial;j<n;j++){
			for (k=qInitial;k<q;k++){
				if (domain[1][j][k]==0){
					rho[1][j][k]=rhoConstant;
					ux[1][j][k]=1-((fIn0[1][j][k]+fIn2[1][j][k]+fIn3[1][j][k]+fIn5[1][j][k]+fIn6[1][j][k]+fIn9[1][j][k]+fIn12[1][j][k]+fIn15[1][j][k]+fIn18[1][j][k])
						+2*(fIn4[1][j][k]+fIn10[1][j][k]+fIn11[1][j][k]+fIn16[1][j][k]+fIn17[1][j][k]))/rho[1][j][k];
					uy[1][j][k]=InletUy;
					uz[1][j][k]=InletUz;
					bRho4=b*rho[1][j][k]/4.0;
					BRho4=1.0-bRho4;
					//pressure[1][j][k]=rho[1][j][k]*R*T/(1-b*rho[1][j][k])-a*alpha*rho[1][j][k]*rho[1][j][k]/(1+2*b*rho[1][j][k]-b*b*rho[1][j][k]*rho[1][j][k]);
					pressure[1][j][k]=rho[1][j][k]*R*T*(1+bRho4+bRho4*bRho4-bRho4*bRho4*bRho4)/(BRho4*BRho4*BRho4)-a*rho[1][j][k]*rho[1][j][k];
					fIn1[1][j][k]=fIn4[1][j][k]+1/3*rho[1][j][k]*ux[1][j][k];
					fIn7[1][j][k]=fIn10[1][j][k]+1/6*rho[1][j][k]*ux[1][j][k];
					fIn8[1][j][k]=fIn11[1][j][k]+1/6*rho[1][j][k]*ux[1][j][k];
					fIn13[1][j][k]=fIn16[1][j][k]+1/6*rho[1][j][k]*ux[1][j][k];
					fIn14[1][j][k]=fIn17[1][j][k]+1/6*rho[1][j][k]*ux[1][j][k];
				}
				if (domain[1][j][k]!=0){
					rho[1][j][k]=solid;
					ux[1][j][k]=0.0;
					uy[1][j][k]=0.0;
					uz[1][j][k]=0.0;
				}
				rho[0][j][k]=rho[1][j][k];
				ux[0][j][k]=ux[1][j][k];
				uy[0][j][k]=uy[1][j][k];
				uz[0][j][k]=uz[1][j][k];
				pressure[0][j][k]=pressure[1][j][k];
			}
		}
	}
}

void lengthForceSave(int layerNumber1,int n,int q,double **SaveFx,double **SaveFy,double **SaveFz,double ***Fx,double ***Fy,double ***Fz,double G){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=0;j<n;j++){
			for (k=0;k<q;k++){
				SaveFx[j][k]=Fx[layerNumber1][j][k]-G;
				SaveFy[j][k]=Fy[layerNumber1][j][k];
				SaveFz[j][k]=Fz[layerNumber1][j][k];
			}
		}
	}
}

void lengthLastForceVelocityConvectiveTwoPhase(int m,int nInitial,int n,int qInitial,int q,double uxAV,double G){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=nInitial;j<n;j++){
			for (k=qInitial;k<q;k++){
				if (domain[m-1][j][k]==0){
					Fx[m-2][j][k]=(LastFx[j][k]+uxAV*(Fx[m-3][j][k]-G))/(1.0+uxAV)+G;
					//Fx[m-1][j][k]=(LastFx[j][k]+ux[m-2][j][k]*(Fx[m-2][j][k]-G))/(1.0+ux[m-2][j][k])+G;
					Fy[m-2][j][k]=(LastFy[j][k]+uxAV*Fy[m-3][j][k])/(1.0+uxAV);
					//Fy[m-1][j][k]=(LastFy[j][k]+ux[m-2][j][k]*Fy[m-2][j][k])/(1.0+ux[m-2][j][k]);
					Fz[m-2][j][k]=(LastFz[j][k]+uxAV*Fz[m-3][j][k])/(1.0+uxAV);
					//Fz[m-1][j][k]=(LastFz[j][k]+ux[m-2][j][k]*Fz[m-2][j][k])/(1.0+ux[m-2][j][k]);
				}
			}
		}
	}
}

void lengthLastVelocityConvectiveTwoPhase(int m,int nInitial,int n,int qInitial,int q,double uxAV,double solid,double b, double R, double T, double a){
	int j;
	int k;
	double bRho4;
	double BRho4;
	double av2=0.0;
#pragma omp parallel private(j,k,bRho4,BRho4)
	{
#pragma omp for schedule(dynamic)
		for (j=nInitial;j<n;j++){
			for (k=qInitial;k<q;k++){
				if (domain[m-2][j][k]==0){
					fIn0[m-2][j][k]=(Last0[j][k]+uxAV*fIn0[m-3][j][k])/(1.0+uxAV);
					//					fIn0[m-2][j][k]=(Last0[j][k]+ux[m-3][j][k]*fIn0[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn1[m-2][j][k]=(Last1[j][k]+uxAV*fIn1[m-3][j][k])/(1.0+uxAV);
					//					fIn1[m-2][j][k]=(Last1[j][k]+ux[m-3][j][k]*fIn1[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn2[m-2][j][k]=(Last2[j][k]+uxAV*fIn2[m-3][j][k])/(1.0+uxAV);
					//					fIn2[m-2][j][k]=(Last2[j][k]+ux[m-3][j][k]*fIn2[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn3[m-2][j][k]=(Last3[j][k]+uxAV*fIn3[m-3][j][k])/(1.0+uxAV);
					//					fIn3[m-2][j][k]=(Last3[j][k]+ux[m-3][j][k]*fIn3[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn4[m-2][j][k]=(Last4[j][k]+uxAV*fIn4[m-3][j][k])/(1.0+uxAV);
					//					fIn4[m-2][j][k]=(Last4[j][k]+ux[m-3][j][k]*fIn4[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn5[m-2][j][k]=(Last5[j][k]+uxAV*fIn5[m-3][j][k])/(1.0+uxAV);
					//					fIn5[m-2][j][k]=(Last5[j][k]+ux[m-3][j][k]*fIn5[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn6[m-2][j][k]=(Last6[j][k]+uxAV*fIn6[m-3][j][k])/(1.0+uxAV);
					//					fIn6[m-2][j][k]=(Last6[j][k]+ux[m-3][j][k]*fIn6[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn7[m-2][j][k]=(Last7[j][k]+uxAV*fIn7[m-3][j][k])/(1.0+uxAV);
					//					fIn7[m-2][j][k]=(Last7[j][k]+ux[m-3][j][k]*fIn7[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn8[m-2][j][k]=(Last8[j][k]+uxAV*fIn8[m-3][j][k])/(1.0+uxAV);
					//					fIn8[m-2][j][k]=(Last8[j][k]+ux[m-3][j][k]*fIn8[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn9[m-2][j][k]=(Last9[j][k]+uxAV*fIn9[m-3][j][k])/(1.0+uxAV);
					//					fIn9[m-2][j][k]=(Last9[j][k]+ux[m-3][j][k]*fIn9[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn10[m-2][j][k]=(Last10[j][k]+uxAV*fIn10[m-3][j][k])/(1.0+uxAV);
					//					fIn10[m-2][j][k]=(Last10[j][k]+ux[m-3][j][k]*fIn10[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn11[m-2][j][k]=(Last11[j][k]+uxAV*fIn11[m-3][j][k])/(1.0+uxAV);
					//					fIn11[m-2][j][k]=(Last11[j][k]+ux[m-3][j][k]*fIn11[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn12[m-2][j][k]=(Last12[j][k]+uxAV*fIn12[m-3][j][k])/(1.0+uxAV);
					//					fIn12[m-2][j][k]=(Last12[j][k]+ux[m-3][j][k]*fIn12[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn13[m-2][j][k]=(Last13[j][k]+uxAV*fIn13[m-3][j][k])/(1.0+uxAV);
					//					fIn13[m-2][j][k]=(Last13[j][k]+ux[m-3][j][k]*fIn13[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn14[m-2][j][k]=(Last14[j][k]+uxAV*fIn14[m-3][j][k])/(1.0+uxAV);
					//					fIn14[m-2][j][k]=(Last14[j][k]+ux[m-3][j][k]*fIn14[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn15[m-2][j][k]=(Last15[j][k]+uxAV*fIn15[m-3][j][k])/(1.0+uxAV);
					//					fIn15[m-2][j][k]=(Last15[j][k]+ux[m-3][j][k]*fIn15[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn16[m-2][j][k]=(Last16[j][k]+uxAV*fIn16[m-3][j][k])/(1.0+uxAV);
					//					fIn16[m-2][j][k]=(Last16[j][k]+ux[m-3][j][k]*fIn16[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn17[m-2][j][k]=(Last17[j][k]+uxAV*fIn17[m-3][j][k])/(1.0+uxAV);
					//					fIn17[m-2][j][k]=(Last17[j][k]+ux[m-3][j][k]*fIn17[m-3][j][k])/(1.0+ux[m-3][j][k]);
					fIn18[m-2][j][k]=(Last18[j][k]+uxAV*fIn18[m-3][j][k])/(1.0+uxAV);
					//					fIn18[m-2][j][k]=(Last18[j][k]+ux[m-3][j][k]*fIn18[m-3][j][k])/(1.0+ux[m-3][j][k]);
					rho[m-2][j][k]=fIn0[m-2][j][k]+fIn1[m-2][j][k]+fIn2[m-2][j][k]+fIn3[m-2][j][k]+fIn4[m-2][j][k]+fIn5[m-2][j][k]+fIn6[m-2][j][k]+fIn7[m-2][j][k]+fIn8[m-2][j][k]
					+fIn9[m-2][j][k]+fIn10[m-2][j][k]+fIn11[m-2][j][k]+fIn12[m-2][j][k]+fIn13[m-2][j][k]+fIn14[m-2][j][k]+fIn15[m-2][j][k]+fIn16[m-2][j][k]+fIn17[m-2][j][k]+fIn18[m-2][j][k];
//					ux[m-2][j][k]=(cxNS[0]*fIn0[m-2][j][k]+cxNS[1]*fIn1[m-2][j][k]+cxNS[2]*fIn2[m-2][j][k]+cxNS[3]*fIn3[m-2][j][k]+cxNS[4]*fIn4[m-2][j][k]+cxNS[5]*fIn5[m-2][j][k]+cxNS[6]*fIn6[m-2][j][k]+cxNS[7]*fIn7[m-2][j][k]+cxNS[8]*fIn8[m-2][j][k]\
//						+cxNS[9]*fIn9[m-2][j][k]+cxNS[10]*fIn10[m-2][j][k]+cxNS[11]*fIn11[m-2][j][k]+cxNS[12]*fIn12[m-2][j][k]+cxNS[13]*fIn13[m-2][j][k]+cxNS[14]*fIn14[m-2][j][k]+cxNS[15]*fIn15[m-2][j][k]+cxNS[16]*fIn16[m-2][j][k]+cxNS[17]*fIn17[m-2][j][k]+cxNS[18]*fIn18[m-2][j][k])/rho[m-2][j][k]+Fx[m-2][j][k]/2.0/rho[m-2][j][k];
					ux[m-2][j][k]=(LastUx[j][k]+uxAV*ux[m-3][j][k])/(1.0+uxAV);
//					ux[m-2][j][k]=(LastUx[j][k]+ux[m-3][j][k]*ux[m-3][j][k])/(1.0+ux[m-3][j][k]);
//					uy[m-2][j][k]=(cyNS[0]*fIn0[m-2][j][k]+cyNS[1]*fIn1[m-2][j][k]+cyNS[2]*fIn2[m-2][j][k]+cyNS[3]*fIn3[m-2][j][k]+cyNS[4]*fIn4[m-2][j][k]+cyNS[5]*fIn5[m-2][j][k]+cyNS[6]*fIn6[m-2][j][k]+cyNS[7]*fIn7[m-2][j][k]+cyNS[8]*fIn8[m-2][j][k]\
//						+cyNS[9]*fIn9[m-2][j][k]+cyNS[10]*fIn10[m-2][j][k]+cyNS[11]*fIn11[m-2][j][k]+cyNS[12]*fIn12[m-2][j][k]+cyNS[13]*fIn13[m-2][j][k]+cyNS[14]*fIn14[m-2][j][k]+cyNS[15]*fIn15[m-2][j][k]+cyNS[16]*fIn16[m-2][j][k]+cyNS[17]*fIn17[m-2][j][k]+cyNS[18]*fIn18[m-2][j][k])/rho[m-2][j][k]+Fy[m-2][j][k]/2.0/rho[m-2][j][k];
					uy[m-2][j][k]=(LastUy[j][k]+uxAV*uy[m-3][j][k])/(1.0+uxAV);
//					uy[m-2][j][k]=(LastUy[j][k]+ux[m-3][j][k]*uy[m-3][j][k])/(1.0+ux[m-3][j][k]);
//					uz[m-2][j][k]=(czNS[0]*fIn0[m-2][j][k]+czNS[1]*fIn1[m-2][j][k]+czNS[2]*fIn2[m-2][j][k]+czNS[3]*fIn3[m-2][j][k]+czNS[4]*fIn4[m-2][j][k]+czNS[5]*fIn5[m-2][j][k]+czNS[6]*fIn6[m-2][j][k]+czNS[7]*fIn7[m-2][j][k]+czNS[8]*fIn8[m-2][j][k]\
//						+czNS[9]*fIn9[m-2][j][k]+czNS[10]*fIn10[m-2][j][k]+czNS[11]*fIn11[m-2][j][k]+czNS[12]*fIn12[m-2][j][k]+czNS[13]*fIn13[m-2][j][k]+czNS[14]*fIn14[m-2][j][k]+czNS[15]*fIn15[m-2][j][k]+czNS[16]*fIn16[m-2][j][k]+czNS[17]*fIn17[m-2][j][k]+czNS[18]*fIn18[m-2][j][k])/rho[m-2][j][k]+Fz[m-2][j][k]/2.0/rho[m-2][j][k];
					uz[m-2][j][k]=(LastUz[j][k]+uxAV*uz[m-3][j][k])/(1.0+uxAV);
//					uz[m-2][j][k]=(LastUz[j][k]+ux[m-3][j][k]*uz[m-3][j][k])/(1.0+ux[m-3][j][k]);
					bRho4=b*rho[m-2][j][k]/4.0;
					BRho4=1.0-bRho4;
					//pressure[m-1][j][k]=rho[m-1][j][k]*R*T/(1-b*rho[m-1][j][k])-a*alpha*rho[m-1][j][k]*rho[m-1][j][k]/(1+2*b*rho[m-1][j][k]-b*b*rho[m-1][j][k]*rho[m-1][j][k]);
					//pressure[m-2][j][k]=rho[m-2][j][k]*R*T*(1+bRho4+bRho4*bRho4-bRho4*bRho4*bRho4)/(BRho4*BRho4*BRho4)-a*rho[m-2][j][k]*rho[m-2][j][k];
					pressure[m-2][j][k]=(LastPressure2[j][k]+uxAV*pressure[m-3][j][k])/(1.0+uxAV);
					//pressure[m-2][j][k]=(LastPressure2[j][k]+ux[m-3][j][k]*pressure[m-3][j][k])/(1.0+ux[m-3][j][k]);
				}
				if (domain[m-2][j][k]!=0){
					rho[m-2][j][k]=solid;
					ux[m-2][j][k]=0.0;
					uy[m-2][j][k]=0.0;
					uz[m-2][j][k]=0.0;
				}
			}
		}
	}
#pragma omp parallel for reduction(+: av2)  
	for (j=nInitial;j<n;j++){
		for (k=qInitial;k<q;k++){
			av2=av2+ux[m-2][j][k];
		}
	}
	av2=av2/(n-nInitial)/(q-qInitial);
#pragma omp parallel private(j,k,bRho4,BRho4)
	{
#pragma omp for schedule(dynamic)
		for (j=nInitial;j<n;j++){
			for (k=qInitial;k<q;k++){
				if (domain[m-1][j][k]==0){
					rho[m-1][j][k]=(LastRho[j][k]+av2*rho[m-2][j][k])/(1.0+av2);
//					rho[m-1][j][k]=(LastRho[j][k]+ux[m-3][j][k]*rho[m-2][j][k])/(1.0+ux[m-3][j][k]);
					ux[m-1][j][k]=ux[m-2][j][k];
					uy[m-1][j][k]=uy[m-2][j][k];
					uz[m-1][j][k]=uz[m-2][j][k];
					bRho4=b*rho[m-1][j][k]/4.0;
					BRho4=1.0-bRho4;
//					pressure[m-1][j][k]=rho[m-1][j][k]*R*T*(1+bRho4+bRho4*bRho4-bRho4*bRho4*bRho4)/(BRho4*BRho4*BRho4)-a*rho[m-1][j][k]*rho[m-1][j][k];
					pressure[m-1][j][k]=(LastPressure1[j][k]+av2*pressure[m-2][j][k])/(1.0+av2);
//					pressure[m-1][j][k]=(LastPressure1[j][k]+ux[m-3][j][k]*pressure[m-2][j][k])/(1.0+ux[m-3][j][k]);
				}
				if (domain[m-1][j][k]!=0){
					rho[m-1][j][k]=solid;
					ux[m-1][j][k]=0.0;
					uy[m-1][j][k]=0.0;
					uz[m-1][j][k]=0.0;
				}
			}
		}
	}
}

void lengthLastVelocityConvectiveGhost(int m,int nInitial,int n,int qInitial,int q,double uxAV,double solid,double gasCritical){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=nInitial;j<n;j++){
			for (k=qInitial;k<q;k++){
				if (domain[m-1][j][k]==0&&rho[m-1][j][k]<gasCritical){
					fIn0Ghost[m-1][j][k]=(Last0Ghost[j][k]+uxAV*fIn0Ghost[m-2][j][k])/(1.0+uxAV);
					fIn1Ghost[m-1][j][k]=(Last1Ghost[j][k]+uxAV*fIn1Ghost[m-2][j][k])/(1.0+uxAV);
					fIn2Ghost[m-1][j][k]=(Last2Ghost[j][k]+uxAV*fIn2Ghost[m-2][j][k])/(1.0+uxAV);
					fIn3Ghost[m-1][j][k]=(Last3Ghost[j][k]+uxAV*fIn3Ghost[m-2][j][k])/(1.0+uxAV);
					fIn4Ghost[m-1][j][k]=(Last4Ghost[j][k]+uxAV*fIn4Ghost[m-2][j][k])/(1.0+uxAV);
					fIn5Ghost[m-1][j][k]=(Last5Ghost[j][k]+uxAV*fIn5Ghost[m-2][j][k])/(1.0+uxAV);
					fIn6Ghost[m-1][j][k]=(Last6Ghost[j][k]+uxAV*fIn6Ghost[m-2][j][k])/(1.0+uxAV);
					fIn7Ghost[m-1][j][k]=(Last7Ghost[j][k]+uxAV*fIn7Ghost[m-2][j][k])/(1.0+uxAV);
					fIn8Ghost[m-1][j][k]=(Last8Ghost[j][k]+uxAV*fIn8Ghost[m-2][j][k])/(1.0+uxAV);
					fIn9Ghost[m-1][j][k]=(Last9Ghost[j][k]+uxAV*fIn9Ghost[m-2][j][k])/(1.0+uxAV);
					fIn10Ghost[m-1][j][k]=(Last10Ghost[j][k]+uxAV*fIn10Ghost[m-2][j][k])/(1.0+uxAV);
					fIn11Ghost[m-1][j][k]=(Last11Ghost[j][k]+uxAV*fIn11Ghost[m-2][j][k])/(1.0+uxAV);
					fIn12Ghost[m-1][j][k]=(Last12Ghost[j][k]+uxAV*fIn12Ghost[m-2][j][k])/(1.0+uxAV);
					fIn13Ghost[m-1][j][k]=(Last13Ghost[j][k]+uxAV*fIn13Ghost[m-2][j][k])/(1.0+uxAV);
					fIn14Ghost[m-1][j][k]=(Last14Ghost[j][k]+uxAV*fIn14Ghost[m-2][j][k])/(1.0+uxAV);
					fIn15Ghost[m-1][j][k]=(Last15Ghost[j][k]+uxAV*fIn15Ghost[m-2][j][k])/(1.0+uxAV);
					fIn16Ghost[m-1][j][k]=(Last16Ghost[j][k]+uxAV*fIn16Ghost[m-2][j][k])/(1.0+uxAV);
					fIn17Ghost[m-1][j][k]=(Last17Ghost[j][k]+uxAV*fIn17Ghost[m-2][j][k])/(1.0+uxAV);
					fIn18Ghost[m-1][j][k]=(Last18Ghost[j][k]+uxAV*fIn18Ghost[m-2][j][k])/(1.0+uxAV);
					rhoGhost[m-1][j][k]=fIn0Ghost[m-1][j][k]+fIn1Ghost[m-1][j][k]+fIn2Ghost[m-1][j][k]+fIn3Ghost[m-1][j][k]+fIn4Ghost[m-1][j][k]+fIn5Ghost[m-1][j][k]+fIn6Ghost[m-1][j][k]+fIn7Ghost[m-1][j][k]+fIn8Ghost[m-1][j][k]
					+fIn9Ghost[m-1][j][k]+fIn10Ghost[m-1][j][k]+fIn11Ghost[m-1][j][k]+fIn12Ghost[m-1][j][k]+fIn13Ghost[m-1][j][k]+fIn14Ghost[m-1][j][k]+fIn15Ghost[m-1][j][k]+fIn16Ghost[m-1][j][k]+fIn17Ghost[m-1][j][k]+fIn18Ghost[m-1][j][k];
					uxGhost[m-1][j][k]=(cxNS[0]*fIn0Ghost[m-1][j][k]+cxNS[1]*fIn1Ghost[m-1][j][k]+cxNS[2]*fIn2Ghost[m-1][j][k]+cxNS[3]*fIn3Ghost[m-1][j][k]+cxNS[4]*fIn4Ghost[m-1][j][k]+cxNS[5]*fIn5Ghost[m-1][j][k]+cxNS[6]*fIn6Ghost[m-1][j][k]+cxNS[7]*fIn7Ghost[m-1][j][k]+cxNS[8]*fIn8Ghost[m-1][j][k]\
						+cxNS[9]*fIn9Ghost[m-1][j][k]+cxNS[10]*fIn10Ghost[m-1][j][k]+cxNS[11]*fIn11Ghost[m-1][j][k]+cxNS[12]*fIn12Ghost[m-1][j][k]+cxNS[13]*fIn13Ghost[m-1][j][k]+cxNS[14]*fIn14Ghost[m-1][j][k]+cxNS[15]*fIn15Ghost[m-1][j][k]+cxNS[16]*fIn16Ghost[m-1][j][k]+cxNS[17]*fIn17Ghost[m-1][j][k]+cxNS[18]*fIn18Ghost[m-1][j][k])/rhoGhost[m-1][j][k]+bodyForce[m-1][j][k]/2.0/rhoGhost[m-1][j][k];
					uyGhost[m-1][j][k]=(cyNS[0]*fIn0Ghost[m-1][j][k]+cyNS[1]*fIn1Ghost[m-1][j][k]+cyNS[2]*fIn2Ghost[m-1][j][k]+cyNS[3]*fIn3Ghost[m-1][j][k]+cyNS[4]*fIn4Ghost[m-1][j][k]+cyNS[5]*fIn5Ghost[m-1][j][k]+cyNS[6]*fIn6Ghost[m-1][j][k]+cyNS[7]*fIn7Ghost[m-1][j][k]+cyNS[8]*fIn8Ghost[m-1][j][k]\
						+cyNS[9]*fIn9Ghost[m-1][j][k]+cyNS[10]*fIn10Ghost[m-1][j][k]+cyNS[11]*fIn11Ghost[m-1][j][k]+cyNS[12]*fIn12Ghost[m-1][j][k]+cyNS[13]*fIn13Ghost[m-1][j][k]+cyNS[14]*fIn14Ghost[m-1][j][k]+cyNS[15]*fIn15Ghost[m-1][j][k]+cyNS[16]*fIn16Ghost[m-1][j][k]+cyNS[17]*fIn17Ghost[m-1][j][k]+cyNS[18]*fIn18Ghost[m-1][j][k])/rhoGhost[m-1][j][k];
					uzGhost[m-1][j][k]=(czNS[0]*fIn0Ghost[m-1][j][k]+czNS[1]*fIn1Ghost[m-1][j][k]+czNS[2]*fIn2Ghost[m-1][j][k]+czNS[3]*fIn3Ghost[m-1][j][k]+czNS[4]*fIn4Ghost[m-1][j][k]+czNS[5]*fIn5Ghost[m-1][j][k]+czNS[6]*fIn6Ghost[m-1][j][k]+czNS[7]*fIn7Ghost[m-1][j][k]+czNS[8]*fIn8Ghost[m-1][j][k]\
						+czNS[9]*fIn9Ghost[m-1][j][k]+czNS[10]*fIn10Ghost[m-1][j][k]+czNS[11]*fIn11Ghost[m-1][j][k]+czNS[12]*fIn12Ghost[m-1][j][k]+czNS[13]*fIn13Ghost[m-1][j][k]+czNS[14]*fIn14Ghost[m-1][j][k]+czNS[15]*fIn15Ghost[m-1][j][k]+czNS[16]*fIn16Ghost[m-1][j][k]+czNS[17]*fIn17Ghost[m-1][j][k]+czNS[18]*fIn18Ghost[m-1][j][k])/rhoGhost[m-1][j][k];
				}
				if (domain[m-1][j][k]!=0||rho[m-1][j][k]>=gasCritical){
					rhoGhost[m-1][j][k]=solid;
					uxGhost[m-1][j][k]=0.0;
					uyGhost[m-1][j][k]=0.0;
					uzGhost[m-1][j][k]=0.0;
				}
			}
		}
	}
}

void lengthLastNeumannBoundaryTwoPhase(int m,int nInitial,int n,int qInitial,int q,double solid){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=nInitial;j<n;j++){
			for (k=qInitial;k<q;k++){
				if (domain[m-2][j][k]==0){
					fIn0[m-2][j][k]=fIn0[m-3][j][k];
					fIn1[m-2][j][k]=fIn1[m-3][j][k];
					fIn2[m-2][j][k]=fIn2[m-3][j][k];
					fIn3[m-2][j][k]=fIn3[m-3][j][k];
					fIn4[m-2][j][k]=fIn4[m-3][j][k];
					fIn5[m-2][j][k]=fIn5[m-3][j][k];
					fIn6[m-2][j][k]=fIn6[m-3][j][k];
					fIn7[m-2][j][k]=fIn7[m-3][j][k];
					fIn8[m-2][j][k]=fIn8[m-3][j][k];
					fIn9[m-2][j][k]=fIn9[m-3][j][k];
					fIn10[m-2][j][k]=fIn10[m-3][j][k];
					fIn11[m-2][j][k]=fIn11[m-3][j][k];
					fIn12[m-2][j][k]=fIn12[m-3][j][k];
					fIn13[m-2][j][k]=fIn13[m-3][j][k];
					fIn14[m-2][j][k]=fIn14[m-3][j][k];
					fIn15[m-2][j][k]=fIn15[m-3][j][k];
					fIn16[m-2][j][k]=fIn16[m-3][j][k];
					fIn17[m-2][j][k]=fIn17[m-3][j][k];
					fIn18[m-2][j][k]=fIn18[m-3][j][k];
					rho[m-2][j][k]=rho[m-3][j][k];
					pressure[m-2][j][k]=pressure[m-3][j][k];
					ux[m-2][j][k]=(cxNS[0]*fIn0[m-2][j][k]+cxNS[1]*fIn1[m-2][j][k]+cxNS[2]*fIn2[m-2][j][k]+cxNS[3]*fIn3[m-2][j][k]+cxNS[4]*fIn4[m-2][j][k]+cxNS[5]*fIn5[m-2][j][k]+cxNS[6]*fIn6[m-2][j][k]+cxNS[7]*fIn7[m-2][j][k]+cxNS[8]*fIn8[m-2][j][k]\
						+cxNS[9]*fIn9[m-2][j][k]+cxNS[10]*fIn10[m-2][j][k]+cxNS[11]*fIn11[m-2][j][k]+cxNS[12]*fIn12[m-2][j][k]+cxNS[13]*fIn13[m-2][j][k]+cxNS[14]*fIn14[m-2][j][k]+cxNS[15]*fIn15[m-2][j][k]+cxNS[16]*fIn16[m-2][j][k]+cxNS[17]*fIn17[m-2][j][k]+cxNS[18]*fIn18[m-2][j][k])/rho[m-2][j][k]+Fx[m-2][j][k]/2.0/rho[m-2][j][k];
					uy[m-2][j][k]=(cyNS[0]*fIn0[m-2][j][k]+cyNS[1]*fIn1[m-2][j][k]+cyNS[2]*fIn2[m-2][j][k]+cyNS[3]*fIn3[m-2][j][k]+cyNS[4]*fIn4[m-2][j][k]+cyNS[5]*fIn5[m-2][j][k]+cyNS[6]*fIn6[m-2][j][k]+cyNS[7]*fIn7[m-2][j][k]+cyNS[8]*fIn8[m-2][j][k]\
						+cyNS[9]*fIn9[m-2][j][k]+cyNS[10]*fIn10[m-2][j][k]+cyNS[11]*fIn11[m-2][j][k]+cyNS[12]*fIn12[m-2][j][k]+cyNS[13]*fIn13[m-2][j][k]+cyNS[14]*fIn14[m-2][j][k]+cyNS[15]*fIn15[m-2][j][k]+cyNS[16]*fIn16[m-2][j][k]+cyNS[17]*fIn17[m-2][j][k]+cyNS[18]*fIn18[m-2][j][k])/rho[m-2][j][k]+Fy[m-2][j][k]/2.0/rho[m-2][j][k];
					uz[m-2][j][k]=(czNS[0]*fIn0[m-2][j][k]+czNS[1]*fIn1[m-2][j][k]+czNS[2]*fIn2[m-2][j][k]+czNS[3]*fIn3[m-2][j][k]+czNS[4]*fIn4[m-2][j][k]+czNS[5]*fIn5[m-2][j][k]+czNS[6]*fIn6[m-2][j][k]+czNS[7]*fIn7[m-2][j][k]+czNS[8]*fIn8[m-2][j][k]\
						+czNS[9]*fIn9[m-2][j][k]+czNS[10]*fIn10[m-2][j][k]+czNS[11]*fIn11[m-2][j][k]+czNS[12]*fIn12[m-2][j][k]+czNS[13]*fIn13[m-2][j][k]+czNS[14]*fIn14[m-2][j][k]+czNS[15]*fIn15[m-2][j][k]+czNS[16]*fIn16[m-2][j][k]+czNS[17]*fIn17[m-2][j][k]+czNS[18]*fIn18[m-2][j][k])/rho[m-2][j][k]+Fz[m-2][j][k]/2.0/rho[m-2][j][k];
				}
				if (domain[m-2][j][k]!=0){
					rho[m-2][j][k]=solid;
					ux[m-2][j][k]=0.0;
					uy[m-2][j][k]=0.0;
					uz[m-2][j][k]=0.0;
				}
			}
		}
	}
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=nInitial;j<n;j++){
			for (k=qInitial;k<q;k++){
				if (domain[m-1][j][k]==0){
					rho[m-1][j][k]=rho[m-2][j][k];
					pressure[m-1][j][k]=pressure[m-2][j][k];
					ux[m-1][j][k]=ux[m-2][j][k];
					uy[m-1][j][k]=uy[m-2][j][k];
					uz[m-1][j][k]=uz[m-2][j][k];
				}
				if (domain[m-1][j][k]!=0){
					rho[m-1][j][k]=solid;
					ux[m-1][j][k]=0.0;
					uy[m-1][j][k]=0.0;
					uz[m-1][j][k]=0.0;
				}
			}
		}
	}
}

void lengthLastConstantVelocityBoundaryTwoPhase(int m,int nInitial,int n,int qInitial,int q,double InletUx,double InletUy,double InletUz,double solid,double b,double R,double T,double a){
	int j;
	int k;
	double bRho4;
	double BRho4;
#pragma omp parallel private(j,k,bRho4,BRho4)
	{
#pragma omp for schedule(dynamic)
		for (j=nInitial;j<n;j++){
			for (k=qInitial;k<q;k++){
				if (domain[m-2][j][k]==0){
					rho[m-2][j][k]=((fIn0[m-2][j][k]+fIn2[m-2][j][k]+fIn3[m-2][j][k]+fIn5[m-2][j][k]+fIn6[m-2][j][k]+fIn9[m-2][j][k]+fIn12[m-2][j][k]+fIn15[m-2][j][k]+fIn18[m-2][j][k])
						+2.0*(fIn1[m-2][j][k]+fIn7[m-2][j][k]+fIn8[m-2][j][k]+fIn13[m-2][j][k]+fIn14[m-2][j][k]))/(1+InletUx);
					ux[m-2][j][k]=InletUx;
					uy[m-2][j][k]=InletUy;
					uz[m-2][j][k]=InletUz;
					bRho4=b*rho[m-2][j][k]/4.0;
					BRho4=1.0-bRho4;
					//pressure[m-2][j][k]=rho[m-2][j][k]*R*T/(1-b*rho[m-2][j][k])-a*alpha*rho[m-2][j][k]*rho[m-2][j][k]/(1+2*b*rho[m-2][j][k]-b*b*rho[m-2][j][k]*rho[m-2][j][k]);
					pressure[m-2][j][k]=rho[m-2][j][k]*R*T*(1+bRho4+bRho4*bRho4-bRho4*bRho4*bRho4)/(BRho4*BRho4*BRho4)-a*rho[m-2][j][k]*rho[m-2][j][k];
					fIn4[m-2][j][k]=fIn1[m-2][j][k]-1.0/3*rho[m-2][j][k]*InletUx;
					fIn10[m-2][j][k]=fIn7[m-2][j][k]-1.0/6*rho[m-2][j][k]*InletUx;
					fIn11[m-2][j][k]=fIn8[m-2][j][k]-1.0/6*rho[m-2][j][k]*InletUx;
					fIn16[m-2][j][k]=fIn13[m-2][j][k]-1.0/6*rho[m-2][j][k]*InletUx;
					fIn17[m-2][j][k]=fIn14[m-2][j][k]-1.0/6*rho[m-2][j][k]*InletUx;
				}
				if (domain[m-2][j][k]!=0){
					rho[m-2][j][k]=solid;
					ux[m-2][j][k]=0.0;
					uy[m-2][j][k]=0.0;
					uz[m-2][j][k]=0.0;
				}
				rho[m-1][j][k]=rho[m-2][j][k];
				ux[m-1][j][k]=ux[m-2][j][k];
				uy[m-1][j][k]=uy[m-2][j][k];
				uz[m-1][j][k]=uz[m-2][j][k];
				pressure[m-1][j][k]=pressure[m-2][j][k];
			}
		}
	}
}

void lengthLastConstantPressureBoundaryTwoPhase(int m,int nInitial,int n,int qInitial,int q,double rhoConstant,double InletUy,double InletUz,double solid,double b,double R,double T,double a){
	int j;
	int k;
	double bRho4;
	double BRho4;
#pragma omp parallel private(j,k,bRho4,BRho4)
	{
#pragma omp for schedule(dynamic)
		for (j=nInitial;j<n;j++){
			for (k=qInitial;k<q;k++){
				if (domain[m-2][j][k]==0){
					rho[m-2][j][k]=rhoConstant;
					ux[m-2][j][k]=((fIn0[m-2][j][k]+fIn2[m-2][j][k]+fIn3[m-2][j][k]+fIn5[m-2][j][k]+fIn6[m-2][j][k]+fIn9[m-2][j][k]+fIn12[m-2][j][k]+fIn15[m-2][j][k]+fIn18[m-2][j][k])
						+2*(fIn1[m-2][j][k]+fIn7[m-2][j][k]+fIn8[m-2][j][k]+fIn13[m-2][j][k]+fIn14[m-2][j][k]))/rho[m-2][j][k]-1;
					uy[m-2][j][k]=InletUy;
					uz[m-2][j][k]=InletUz;
					bRho4=b*rho[m-2][j][k]/4.0;
					BRho4=1.0-bRho4;
					//pressure[m-2][j][k]=rho[m-2][j][k]*R*T/(1-b*rho[m-2][j][k])-a*alpha*rho[m-2][j][k]*rho[m-2][j][k]/(1+2*b*rho[m-2][j][k]-b*b*rho[m-2][j][k]*rho[m-2][j][k]);
					pressure[m-2][j][k]=rho[m-2][j][k]*R*T*(1+bRho4+bRho4*bRho4-bRho4*bRho4*bRho4)/(BRho4*BRho4*BRho4)-a*rho[m-2][j][k]*rho[m-2][j][k];
					fIn4[m-2][j][k]=fIn1[m-2][j][k]-1/3*rho[m-2][j][k]*ux[m-2][j][k];
					fIn10[m-2][j][k]=fIn7[m-2][j][k]-1/6*rho[m-2][j][k]*ux[m-2][j][k];
					fIn11[m-2][j][k]=fIn8[m-2][j][k]-1/6*rho[m-2][j][k]*ux[m-2][j][k];
					fIn16[m-2][j][k]=fIn13[m-2][j][k]-1/6*rho[m-2][j][k]*ux[m-2][j][k];
					fIn17[m-2][j][k]=fIn14[m-2][j][k]-1/6*rho[m-2][j][k]*ux[m-2][j][k];
				}
				if (domain[m-2][j][k]!=0){
					rho[m-2][j][k]=solid;
					ux[m-2][j][k]=0.0;
					uy[m-2][j][k]=0.0;
					uz[m-2][j][k]=0.0;
				}
				rho[m-1][j][k]=rho[m-2][j][k];
				ux[m-1][j][k]=ux[m-2][j][k];
				uy[m-1][j][k]=uy[m-2][j][k];
				uz[m-1][j][k]=uz[m-2][j][k];
				pressure[m-1][j][k]=pressure[m-2][j][k];
			}
		}
	}
}

void lengthFirstDirichletConcentration(int Width,int Height,double InletConcentration,double gasCritical){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
				if (domain[0][j][k]==0&&rho[0][j][k]<gasCritical){
					Concentration[0][j][k]=InletConcentration;
					cIn1[0][j][k]=-Out4[0][j][k]+2*tNS[4]*Concentration[0][j][k]*(1+4.5*((cxNS[4]*ux[0][j][k]+cyNS[4]*uy[0][j][k]+czNS[4]*uz[0][j][k])*(cxNS[4]*ux[0][j][k]+cyNS[4]*uy[0][j][k]+czNS[4]*uz[0][j][k]))-1.5*(ux[0][j][k]*ux[0][j][k]+uy[0][j][k]*uy[0][j][k]+uz[0][j][k]*uz[0][j][k]));
					cIn7[0][j][k]=-Out10[0][j][k]+2*tNS[10]*Concentration[0][j][k]*(1+4.5*((cxNS[10]*ux[0][j][k]+cyNS[10]*uy[0][j][k]+czNS[10]*uz[0][j][k])*(cxNS[10]*ux[0][j][k]+cyNS[10]*uy[0][j][k]+czNS[10]*uz[0][j][k]))-1.5*(ux[0][j][k]*ux[0][j][k]+uy[0][j][k]*uy[0][j][k]+uz[0][j][k]*uz[0][j][k]));
					cIn8[0][j][k]=-Out11[0][j][k]+2*tNS[11]*Concentration[0][j][k]*(1+4.5*((cxNS[11]*ux[0][j][k]+cyNS[11]*uy[0][j][k]+czNS[11]*uz[0][j][k])*(cxNS[11]*ux[0][j][k]+cyNS[11]*uy[0][j][k]+czNS[11]*uz[0][j][k]))-1.5*(ux[0][j][k]*ux[0][j][k]+uy[0][j][k]*uy[0][j][k]+uz[0][j][k]*uz[0][j][k]));
					cIn13[0][j][k]=-Out16[0][j][k]+2*tNS[16]*Concentration[0][j][k]*(1+4.5*((cxNS[16]*ux[0][j][k]+cyNS[16]*uy[0][j][k]+czNS[16]*uz[0][j][k])*(cxNS[16]*ux[0][j][k]+cyNS[16]*uy[0][j][k]+czNS[16]*uz[0][j][k]))-1.5*(ux[0][j][k]*ux[0][j][k]+uy[0][j][k]*uy[0][j][k]+uz[0][j][k]*uz[0][j][k]));
					cIn14[0][j][k]=-Out17[0][j][k]+2*tNS[17]*Concentration[0][j][k]*(1+4.5*((cxNS[17]*ux[0][j][k]+cyNS[17]*uy[0][j][k]+czNS[17]*uz[0][j][k])*(cxNS[17]*ux[0][j][k]+cyNS[17]*uy[0][j][k]+czNS[17]*uz[0][j][k]))-1.5*(ux[0][j][k]*ux[0][j][k]+uy[0][j][k]*uy[0][j][k]+uz[0][j][k]*uz[0][j][k]));
				}
				if (domain[0][j][k]!=0||rho[0][j][k]>=gasCritical){
					Concentration[0][j][k]=0.0;
				}
			}
		}
	}
}

void lengthLastConvectiveConcentration(int Length,int nInitial,int Width,int qInitial,int Height,double uxAV,double gasCritical){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=nInitial;j<Width;j++){
			for (k=qInitial;k<Height;k++){
				if (domain[Length-1][j][k]==0&&rho[Length-1][j][k]<gasCritical){
					cIn0[Length-1][j][k]=(Last0Concentration[j][k]+uxAV*cIn0[Length-2][j][k])/(1+uxAV);
					cIn1[Length-1][j][k]=(Last1Concentration[j][k]+uxAV*cIn1[Length-2][j][k])/(1+uxAV);
					cIn2[Length-1][j][k]=(Last2Concentration[j][k]+uxAV*cIn2[Length-2][j][k])/(1+uxAV);
					cIn3[Length-1][j][k]=(Last3Concentration[j][k]+uxAV*cIn3[Length-2][j][k])/(1+uxAV);
					cIn4[Length-1][j][k]=(Last4Concentration[j][k]+uxAV*cIn4[Length-2][j][k])/(1+uxAV);
					cIn5[Length-1][j][k]=(Last5Concentration[j][k]+uxAV*cIn5[Length-2][j][k])/(1+uxAV);
					cIn6[Length-1][j][k]=(Last6Concentration[j][k]+uxAV*cIn6[Length-2][j][k])/(1+uxAV);
					cIn7[Length-1][j][k]=(Last7Concentration[j][k]+uxAV*cIn7[Length-2][j][k])/(1+uxAV);
					cIn8[Length-1][j][k]=(Last8Concentration[j][k]+uxAV*cIn8[Length-2][j][k])/(1+uxAV);
					cIn9[Length-1][j][k]=(Last9Concentration[j][k]+uxAV*cIn9[Length-2][j][k])/(1+uxAV);
					cIn10[Length-1][j][k]=(Last10Concentration[j][k]+uxAV*cIn10[Length-2][j][k])/(1+uxAV);
					cIn11[Length-1][j][k]=(Last11Concentration[j][k]+uxAV*cIn11[Length-2][j][k])/(1+uxAV);
					cIn12[Length-1][j][k]=(Last12Concentration[j][k]+uxAV*cIn12[Length-2][j][k])/(1+uxAV);
					cIn13[Length-1][j][k]=(Last13Concentration[j][k]+uxAV*cIn13[Length-2][j][k])/(1+uxAV);
					cIn14[Length-1][j][k]=(Last14Concentration[j][k]+uxAV*cIn14[Length-2][j][k])/(1+uxAV);
					cIn15[Length-1][j][k]=(Last15Concentration[j][k]+uxAV*cIn15[Length-2][j][k])/(1+uxAV);
					cIn16[Length-1][j][k]=(Last16Concentration[j][k]+uxAV*cIn16[Length-2][j][k])/(1+uxAV);
					cIn17[Length-1][j][k]=(Last17Concentration[j][k]+uxAV*cIn17[Length-2][j][k])/(1+uxAV);
					cIn18[Length-1][j][k]=(Last18Concentration[j][k]+uxAV*cIn18[Length-2][j][k])/(1+uxAV);
					Concentration[Length-1][j][k]=cIn0[Length-1][j][k]+cIn1[Length-1][j][k]+cIn2[Length-1][j][k]+cIn3[Length-1][j][k]+cIn4[Length-1][j][k]+cIn5[Length-1][j][k]+cIn6[Length-1][j][k]+cIn7[Length-1][j][k]+cIn8[Length-1][j][k]
					+cIn9[Length-1][j][k]+cIn10[Length-1][j][k]+cIn11[Length-1][j][k]+cIn12[Length-1][j][k]+cIn13[Length-1][j][k]+cIn14[Length-1][j][k]+cIn15[Length-1][j][k]+cIn16[Length-1][j][k]+cIn17[Length-1][j][k]+cIn18[Length-1][j][k];
				}
				if (domain[Length-1][j][k]!=0||rho[Length-1][j][k]>=gasCritical){
					Concentration[Length-1][j][k]=0.0;
				}
			}
		}
	}
}

void lengthLastNeumannConcentration(int Length,int Width,int Height){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
					cIn0[Length-1][j][k]=cIn0[Length-2][j][k];
					cIn1[Length-1][j][k]=cIn1[Length-2][j][k];
					cIn2[Length-1][j][k]=cIn2[Length-2][j][k];
					cIn3[Length-1][j][k]=cIn3[Length-2][j][k];
					cIn4[Length-1][j][k]=cIn4[Length-2][j][k];
					cIn5[Length-1][j][k]=cIn5[Length-2][j][k];
					cIn6[Length-1][j][k]=cIn6[Length-2][j][k];
					cIn7[Length-1][j][k]=cIn7[Length-2][j][k];
					cIn8[Length-1][j][k]=cIn8[Length-2][j][k];
					cIn9[Length-1][j][k]=cIn9[Length-2][j][k];
					cIn10[Length-1][j][k]=cIn10[Length-2][j][k];
					cIn11[Length-1][j][k]=cIn11[Length-2][j][k];
					cIn12[Length-1][j][k]=cIn12[Length-2][j][k];
					cIn13[Length-1][j][k]=cIn13[Length-2][j][k];
					cIn14[Length-1][j][k]=cIn14[Length-2][j][k];
					cIn15[Length-1][j][k]=cIn15[Length-2][j][k];
					cIn16[Length-1][j][k]=cIn16[Length-2][j][k];
					cIn17[Length-1][j][k]=cIn17[Length-2][j][k];
					cIn18[Length-1][j][k]=cIn18[Length-2][j][k];
					Concentration[Length-1][j][k]=Concentration[Length-2][j][k];
			}
		}
	}
}

void gasLiquidBoundaryConcentration(int m, int n, int q, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical){
						if (rho1[i][j][k]<=gasCritical){
							cIn1[i][j][k]=-Out4[i][j][k]+2*tNS[1]*Concentration[i][j][k];
						}
						if (rho2[i][j][k]<=gasCritical){
							cIn2[i][j][k]=-Out5[i][j][k]+2*tNS[2]*Concentration[i][j][k];
						}
						if (rho3[i][j][k]<=gasCritical){
							cIn3[i][j][k]=-Out6[i][j][k]+2*tNS[3]*Concentration[i][j][k];
						}
						if (rho4[i][j][k]<=gasCritical){
							cIn4[i][j][k]=-Out1[i][j][k]+2*tNS[4]*Concentration[i][j][k];
						}
						if (rho5[i][j][k]<=gasCritical){
							cIn5[i][j][k]=-Out2[i][j][k]+2*tNS[5]*Concentration[i][j][k];
						}
						if (rho6[i][j][k]<=gasCritical){
							cIn6[i][j][k]=-Out3[i][j][k]+2*tNS[6]*Concentration[i][j][k];
						}
						if (rho7[i][j][k]<=gasCritical){
							cIn7[i][j][k]=-Out10[i][j][k]+2*tNS[7]*Concentration[i][j][k];
						}
						if (rho8[i][j][k]<=gasCritical){
							cIn8[i][j][k]=-Out11[i][j][k]+2*tNS[8]*Concentration[i][j][k];
						}
						if (rho9[i][j][k]<=gasCritical){
							cIn9[i][j][k]=-Out12[i][j][k]+2*tNS[9]*Concentration[i][j][k];
						}
						if (rho10[i][j][k]<=gasCritical){
							cIn10[i][j][k]=-Out7[i][j][k]+2*tNS[10]*Concentration[i][j][k];
						}
						if (rho11[i][j][k]<=gasCritical){
							cIn11[i][j][k]=-Out8[i][j][k]+2*tNS[11]*Concentration[i][j][k];
						}
						if (rho12[i][j][k]<=gasCritical){
							cIn12[i][j][k]=-Out9[i][j][k]+2*tNS[12]*Concentration[i][j][k];
						}
						if (rho13[i][j][k]<=gasCritical){
							cIn13[i][j][k]=-Out16[i][j][k]+2*tNS[13]*Concentration[i][j][k];
						}
						if (rho14[i][j][k]<=gasCritical){
							cIn14[i][j][k]=-Out17[i][j][k]+2*tNS[14]*Concentration[i][j][k];
						}
						if (rho15[i][j][k]<=gasCritical){
							cIn15[i][j][k]=-Out18[i][j][k]+2*tNS[15]*Concentration[i][j][k];
						}
						if (rho16[i][j][k]<=gasCritical){
							cIn16[i][j][k]=-Out13[i][j][k]+2*tNS[16]*Concentration[i][j][k];
						}
						if (rho17[i][j][k]<=gasCritical){
							cIn17[i][j][k]=-Out14[i][j][k]+2*tNS[17]*Concentration[i][j][k];
						}
						if (rho18[i][j][k]<=gasCritical){
							cIn18[i][j][k]=-Out15[i][j][k]+2*tNS[18]*Concentration[i][j][k];
						}
					}
				}
			}
		}
	}
}

void gasLiquidBoundaryConcentrationD3Q7(int m, int n, int q, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical){
						if (rho1[i][j][k]<=gasCritical){
//							cIn1[i][j][k]=-Out4[i][j][k]+2*tNSD3Q7[1]*Concentration[i][j][k];
							cIn1[i][j][k]=Out4[i][j][k];
						}
						if (rho2[i][j][k]<=gasCritical){
//							cIn2[i][j][k]=-Out5[i][j][k]+2*tNSD3Q7[2]*Concentration[i][j][k];
							cIn2[i][j][k]=Out5[i][j][k];
						}
						if (rho3[i][j][k]<=gasCritical){
//							cIn3[i][j][k]=-Out6[i][j][k]+2*tNSD3Q7[3]*Concentration[i][j][k];
							cIn3[i][j][k]=Out6[i][j][k];
						}
						if (rho4[i][j][k]<=gasCritical){
//							cIn4[i][j][k]=-Out1[i][j][k]+2*tNSD3Q7[4]*Concentration[i][j][k];
							cIn4[i][j][k]=Out1[i][j][k];
						}
						if (rho5[i][j][k]<=gasCritical){
//							cIn5[i][j][k]=-Out2[i][j][k]+2*tNSD3Q7[5]*Concentration[i][j][k];
							cIn5[i][j][k]=Out2[i][j][k];
						}
						if (rho6[i][j][k]<=gasCritical){
//							cIn6[i][j][k]=-Out3[i][j][k]+2*tNSD3Q7[6]*Concentration[i][j][k];
							cIn6[i][j][k]=Out3[i][j][k];
						}
					}
				}
			}
		}
	}
}

void gasLiquidBoundaryConcentrationD3Q15(int m, int n, int q, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<m; i++){
			for (j = 0; j<n; j++){
				for (k = 0; k<q; k++){
					if (domain[i][j][k] == 0 && rho[i][j][k]>gasCritical){
						if (rho1[i][j][k] <= gasCritical){
							cIn1[i][j][k] = Out2[i][j][k];
						}
						if (rho2[i][j][k] <= gasCritical){
							cIn3[i][j][k] = Out4[i][j][k];
						}
						if (rho3[i][j][k] <= gasCritical){
							cIn5[i][j][k] = Out6[i][j][k];
						}
						if (rho4[i][j][k] <= gasCritical){
							cIn2[i][j][k] = Out1[i][j][k];
						}
						if (rho5[i][j][k] <= gasCritical){
							cIn4[i][j][k] = Out3[i][j][k];
						}
						if (rho6[i][j][k] <= gasCritical){
							cIn6[i][j][k] = Out5[i][j][k];
						}
						if (rho19[i][j][k] <= gasCritical){
							cIn7[i][j][k] = Out14[i][j][k];
						}
						if (rho20[i][j][k] <= gasCritical){
							cIn8[i][j][k] = Out13[i][j][k];
						}
						if (rho21[i][j][k] <= gasCritical){
							cIn9[i][j][k] = Out12[i][j][k];
						}
						if (rho22[i][j][k] <= gasCritical){
							cIn10[i][j][k] = Out11[i][j][k];
						}
						if (rho23[i][j][k] <= gasCritical){
							cIn11[i][j][k] = Out10[i][j][k];
						}
						if (rho24[i][j][k] <= gasCritical){
							cIn12[i][j][k] = Out9[i][j][k];
						}
						if (rho25[i][j][k] <= gasCritical){
							cIn13[i][j][k] = Out8[i][j][k];
						}
						if (rho26[i][j][k] <= gasCritical){
							cIn14[i][j][k] = Out7[i][j][k];
						}
					}
				}
			}
		}
	}
}

void gasLiquidBoundaryConcentrationFirst(int m, int n, int q, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical){
						if (rho1[i][j][k]<=gasCritical){
							cIn1First[i][j][k]=-Out4[i][j][k]+2*tNS[1]*ConcentrationFirst[i][j][k];
						}
						if (rho2[i][j][k]<=gasCritical){
							cIn2First[i][j][k]=-Out5[i][j][k]+2*tNS[2]*ConcentrationFirst[i][j][k];
						}
						if (rho3[i][j][k]<=gasCritical){
							cIn3First[i][j][k]=-Out6[i][j][k]+2*tNS[3]*ConcentrationFirst[i][j][k];
						}
						if (rho4[i][j][k]<=gasCritical){
							cIn4First[i][j][k]=-Out1[i][j][k]+2*tNS[4]*ConcentrationFirst[i][j][k];
						}
						if (rho5[i][j][k]<=gasCritical){
							cIn5First[i][j][k]=-Out2[i][j][k]+2*tNS[5]*ConcentrationFirst[i][j][k];
						}
						if (rho6[i][j][k]<=gasCritical){
							cIn6First[i][j][k]=-Out3[i][j][k]+2*tNS[6]*ConcentrationFirst[i][j][k];
						}
						if (rho7[i][j][k]<=gasCritical){
							cIn7First[i][j][k]=-Out10[i][j][k]+2*tNS[7]*ConcentrationFirst[i][j][k];
						}
						if (rho8[i][j][k]<=gasCritical){
							cIn8First[i][j][k]=-Out11[i][j][k]+2*tNS[8]*ConcentrationFirst[i][j][k];
						}
						if (rho9[i][j][k]<=gasCritical){
							cIn9First[i][j][k]=-Out12[i][j][k]+2*tNS[9]*ConcentrationFirst[i][j][k];
						}
						if (rho10[i][j][k]<=gasCritical){
							cIn10First[i][j][k]=-Out7[i][j][k]+2*tNS[10]*ConcentrationFirst[i][j][k];
						}
						if (rho11[i][j][k]<=gasCritical){
							cIn11First[i][j][k]=-Out8[i][j][k]+2*tNS[11]*ConcentrationFirst[i][j][k];
						}
						if (rho12[i][j][k]<=gasCritical){
							cIn12First[i][j][k]=-Out9[i][j][k]+2*tNS[12]*ConcentrationFirst[i][j][k];
						}
						if (rho13[i][j][k]<=gasCritical){
							cIn13First[i][j][k]=-Out16[i][j][k]+2*tNS[13]*ConcentrationFirst[i][j][k];
						}
						if (rho14[i][j][k]<=gasCritical){
							cIn14First[i][j][k]=-Out17[i][j][k]+2*tNS[14]*ConcentrationFirst[i][j][k];
						}
						if (rho15[i][j][k]<=gasCritical){
							cIn15First[i][j][k]=-Out18[i][j][k]+2*tNS[15]*ConcentrationFirst[i][j][k];
						}
						if (rho16[i][j][k]<=gasCritical){
							cIn16First[i][j][k]=-Out13[i][j][k]+2*tNS[16]*ConcentrationFirst[i][j][k];
						}
						if (rho17[i][j][k]<=gasCritical){
							cIn17First[i][j][k]=-Out14[i][j][k]+2*tNS[17]*ConcentrationFirst[i][j][k];
						}
						if (rho18[i][j][k]<=gasCritical){
							cIn18First[i][j][k]=-Out15[i][j][k]+2*tNS[18]*ConcentrationFirst[i][j][k];
						}
					}
				}
			}
		}
	}
}

void gasLiquidBoundaryConcentrationFirstD3Q7(int m, int n, int q, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical){
						if (rho1[i][j][k]<=gasCritical){
//							cIn1First[i][j][k]=-Out4[i][j][k]+2*tNSD3Q7[1]*ConcentrationFirst[i][j][k];
							cIn1First[i][j][k]=Out4[i][j][k];
						}
						if (rho2[i][j][k]<=gasCritical){
//							cIn2First[i][j][k]=-Out5[i][j][k]+2*tNSD3Q7[2]*ConcentrationFirst[i][j][k];
							cIn2First[i][j][k]=Out5[i][j][k];
						}
						if (rho3[i][j][k]<=gasCritical){
//							cIn3First[i][j][k]=-Out6[i][j][k]+2*tNSD3Q7[3]*ConcentrationFirst[i][j][k];
							cIn3First[i][j][k]=Out6[i][j][k];
						}
						if (rho4[i][j][k]<=gasCritical){
//							cIn4First[i][j][k]=-Out1[i][j][k]+2*tNSD3Q7[4]*ConcentrationFirst[i][j][k];
							cIn4First[i][j][k]=Out1[i][j][k];
						}
						if (rho5[i][j][k]<=gasCritical){
//							cIn5First[i][j][k]=-Out2[i][j][k]+2*tNSD3Q7[5]*ConcentrationFirst[i][j][k];
							cIn5First[i][j][k]=Out2[i][j][k];
						}
						if (rho6[i][j][k]<=gasCritical){
//							cIn6First[i][j][k]=-Out3[i][j][k]+2*tNSD3Q7[6]*ConcentrationFirst[i][j][k];
							cIn6First[i][j][k]=Out3[i][j][k];
						}
					}
				}
			}
		}
	}
}

void gasLiquidBoundaryConcentrationSecond(int m, int n, int q, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical){
						if (rho1[i][j][k]<=gasCritical){
							cIn1Second[i][j][k]=-Out4[i][j][k]+2*tNS[1]*ConcentrationSecond[i][j][k];
						}
						if (rho2[i][j][k]<=gasCritical){
							cIn2Second[i][j][k]=-Out5[i][j][k]+2*tNS[2]*ConcentrationSecond[i][j][k];
						}
						if (rho3[i][j][k]<=gasCritical){
							cIn3Second[i][j][k]=-Out6[i][j][k]+2*tNS[3]*ConcentrationSecond[i][j][k];
						}
						if (rho4[i][j][k]<=gasCritical){
							cIn4Second[i][j][k]=-Out1[i][j][k]+2*tNS[4]*ConcentrationSecond[i][j][k];
						}
						if (rho5[i][j][k]<=gasCritical){
							cIn5Second[i][j][k]=-Out2[i][j][k]+2*tNS[5]*ConcentrationSecond[i][j][k];
						}
						if (rho6[i][j][k]<=gasCritical){
							cIn6Second[i][j][k]=-Out3[i][j][k]+2*tNS[6]*ConcentrationSecond[i][j][k];
						}
						if (rho7[i][j][k]<=gasCritical){
							cIn7Second[i][j][k]=-Out10[i][j][k]+2*tNS[7]*ConcentrationSecond[i][j][k];
						}
						if (rho8[i][j][k]<=gasCritical){
							cIn8Second[i][j][k]=-Out11[i][j][k]+2*tNS[8]*ConcentrationSecond[i][j][k];
						}
						if (rho9[i][j][k]<=gasCritical){
							cIn9Second[i][j][k]=-Out12[i][j][k]+2*tNS[9]*ConcentrationSecond[i][j][k];
						}
						if (rho10[i][j][k]<=gasCritical){
							cIn10Second[i][j][k]=-Out7[i][j][k]+2*tNS[10]*ConcentrationSecond[i][j][k];
						}
						if (rho11[i][j][k]<=gasCritical){
							cIn11Second[i][j][k]=-Out8[i][j][k]+2*tNS[11]*ConcentrationSecond[i][j][k];
						}
						if (rho12[i][j][k]<=gasCritical){
							cIn12Second[i][j][k]=-Out9[i][j][k]+2*tNS[12]*ConcentrationSecond[i][j][k];
						}
						if (rho13[i][j][k]<=gasCritical){
							cIn13Second[i][j][k]=-Out16[i][j][k]+2*tNS[13]*ConcentrationSecond[i][j][k];
						}
						if (rho14[i][j][k]<=gasCritical){
							cIn14Second[i][j][k]=-Out17[i][j][k]+2*tNS[14]*ConcentrationSecond[i][j][k];
						}
						if (rho15[i][j][k]<=gasCritical){
							cIn15Second[i][j][k]=-Out18[i][j][k]+2*tNS[15]*ConcentrationSecond[i][j][k];
						}
						if (rho16[i][j][k]<=gasCritical){
							cIn16Second[i][j][k]=-Out13[i][j][k]+2*tNS[16]*ConcentrationSecond[i][j][k];
						}
						if (rho17[i][j][k]<=gasCritical){
							cIn17Second[i][j][k]=-Out14[i][j][k]+2*tNS[17]*ConcentrationSecond[i][j][k];
						}
						if (rho18[i][j][k]<=gasCritical){
							cIn18Second[i][j][k]=-Out15[i][j][k]+2*tNS[18]*ConcentrationSecond[i][j][k];
						}
					}
				}
			}
		}
	}
}

void gasLiquidBoundaryConcentrationSecondD3Q7(int m, int n, int q, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical){
						if (rho1[i][j][k]<=gasCritical){
//							cIn1Second[i][j][k]=-Out4[i][j][k]+2*tNSD3Q7[1]*ConcentrationSecond[i][j][k];
							cIn1Second[i][j][k]=Out4[i][j][k];
						}
						if (rho2[i][j][k]<=gasCritical){
//							cIn2Second[i][j][k]=-Out5[i][j][k]+2*tNSD3Q7[2]*ConcentrationSecond[i][j][k];
							cIn2Second[i][j][k]=Out5[i][j][k];
						}
						if (rho3[i][j][k]<=gasCritical){
//							cIn3Second[i][j][k]=-Out6[i][j][k]+2*tNSD3Q7[3]*ConcentrationSecond[i][j][k];
							cIn3Second[i][j][k]=Out6[i][j][k];
						}
						if (rho4[i][j][k]<=gasCritical){
//							cIn4Second[i][j][k]=-Out1[i][j][k]+2*tNSD3Q7[4]*ConcentrationSecond[i][j][k];
							cIn4Second[i][j][k]=Out1[i][j][k];
						}
						if (rho5[i][j][k]<=gasCritical){
//							cIn5Second[i][j][k]=-Out2[i][j][k]+2*tNSD3Q7[5]*ConcentrationSecond[i][j][k];
							cIn5Second[i][j][k]=Out2[i][j][k];
						}
						if (rho6[i][j][k]<=gasCritical){
//							cIn6Second[i][j][k]=-Out3[i][j][k]+2*tNSD3Q7[6]*ConcentrationSecond[i][j][k];
							cIn6Second[i][j][k]=Out3[i][j][k];
						}
					}
				}
			}
		}
	}
}

void gasLiquidBoundaryVelocityGhost(int m, int n, int q, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical){
						if (rho1[i][j][k]<=gasCritical){
							fIn1Ghost[i][j][k]=Out4[i][j][k];
						}
						if (rho2[i][j][k]<=gasCritical){
							fIn2Ghost[i][j][k]=Out5[i][j][k];
						}
						if (rho3[i][j][k]<=gasCritical){
							fIn3Ghost[i][j][k]=Out6[i][j][k];
						}
						if (rho4[i][j][k]<=gasCritical){
							fIn4Ghost[i][j][k]=Out1[i][j][k];
						}
						if (rho5[i][j][k]<=gasCritical){
							fIn5Ghost[i][j][k]=Out2[i][j][k];
						}
						if (rho6[i][j][k]<=gasCritical){
							fIn6Ghost[i][j][k]=Out3[i][j][k];
						}
						if (rho7[i][j][k]<=gasCritical){
							fIn7Ghost[i][j][k]=Out10[i][j][k];
						}
						if (rho8[i][j][k]<=gasCritical){
							fIn8Ghost[i][j][k]=Out11[i][j][k];
						}
						if (rho9[i][j][k]<=gasCritical){
							fIn9Ghost[i][j][k]=Out12[i][j][k];
						}
						if (rho10[i][j][k]<=gasCritical){
							fIn10Ghost[i][j][k]=Out7[i][j][k];
						}
						if (rho11[i][j][k]<=gasCritical){
							fIn11Ghost[i][j][k]=Out8[i][j][k];
						}
						if (rho12[i][j][k]<=gasCritical){
							fIn12Ghost[i][j][k]=Out9[i][j][k];
						}
						if (rho13[i][j][k]<=gasCritical){
							fIn13Ghost[i][j][k]=Out16[i][j][k];
						}
						if (rho14[i][j][k]<=gasCritical){
							fIn14Ghost[i][j][k]=Out17[i][j][k];
						}
						if (rho15[i][j][k]<=gasCritical){
							fIn15Ghost[i][j][k]=Out18[i][j][k];
						}
						if (rho16[i][j][k]<=gasCritical){
							fIn16Ghost[i][j][k]=Out13[i][j][k];
						}
						if (rho17[i][j][k]<=gasCritical){
							fIn17Ghost[i][j][k]=Out14[i][j][k];
						}
						if (rho18[i][j][k]<=gasCritical){
							fIn18Ghost[i][j][k]=Out15[i][j][k];
						}
					}
				}
			}
		}
	}
}

void gasLiquidBoundaryElectrolytePotential(int m, int n, int q, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical){
						if (rho1[i][j][k]<=gasCritical){
							electrolyteIn1[i][j][k]=-Out4[i][j][k]+2*tNS[1]*electrolytePotential[i][j][k];
						}
						if (rho2[i][j][k]<=gasCritical){
							electrolyteIn2[i][j][k]=-Out5[i][j][k]+2*tNS[2]*electrolytePotential[i][j][k];
						}
						if (rho3[i][j][k]<=gasCritical){
							electrolyteIn3[i][j][k]=-Out6[i][j][k]+2*tNS[3]*electrolytePotential[i][j][k];
						}
						if (rho4[i][j][k]<=gasCritical){
							electrolyteIn4[i][j][k]=-Out1[i][j][k]+2*tNS[4]*electrolytePotential[i][j][k];
						}
						if (rho5[i][j][k]<=gasCritical){
							electrolyteIn5[i][j][k]=-Out2[i][j][k]+2*tNS[5]*electrolytePotential[i][j][k];
						}
						if (rho6[i][j][k]<=gasCritical){
							electrolyteIn6[i][j][k]=-Out3[i][j][k]+2*tNS[6]*electrolytePotential[i][j][k];
						}
						if (rho7[i][j][k]<=gasCritical){
							electrolyteIn7[i][j][k]=-Out10[i][j][k]+2*tNS[7]*electrolytePotential[i][j][k];
						}
						if (rho8[i][j][k]<=gasCritical){
							electrolyteIn8[i][j][k]=-Out11[i][j][k]+2*tNS[8]*electrolytePotential[i][j][k];
						}
						if (rho9[i][j][k]<=gasCritical){
							electrolyteIn9[i][j][k]=-Out12[i][j][k]+2*tNS[9]*electrolytePotential[i][j][k];
						}
						if (rho10[i][j][k]<=gasCritical){
							electrolyteIn10[i][j][k]=-Out7[i][j][k]+2*tNS[10]*electrolytePotential[i][j][k];
						}
						if (rho11[i][j][k]<=gasCritical){
							electrolyteIn11[i][j][k]=-Out8[i][j][k]+2*tNS[11]*electrolytePotential[i][j][k];
						}
						if (rho12[i][j][k]<=gasCritical){
							electrolyteIn12[i][j][k]=-Out9[i][j][k]+2*tNS[12]*electrolytePotential[i][j][k];
						}
						if (rho13[i][j][k]<=gasCritical){
							electrolyteIn13[i][j][k]=-Out16[i][j][k]+2*tNS[13]*electrolytePotential[i][j][k];
						}
						if (rho14[i][j][k]<=gasCritical){
							electrolyteIn14[i][j][k]=-Out17[i][j][k]+2*tNS[14]*electrolytePotential[i][j][k];
						}
						if (rho15[i][j][k]<=gasCritical){
							electrolyteIn15[i][j][k]=-Out18[i][j][k]+2*tNS[15]*electrolytePotential[i][j][k];
						}
						if (rho16[i][j][k]<=gasCritical){
							electrolyteIn16[i][j][k]=-Out13[i][j][k]+2*tNS[16]*electrolytePotential[i][j][k];
						}
						if (rho17[i][j][k]<=gasCritical){
							electrolyteIn17[i][j][k]=-Out14[i][j][k]+2*tNS[17]*electrolytePotential[i][j][k];
						}
						if (rho18[i][j][k]<=gasCritical){
							electrolyteIn18[i][j][k]=-Out15[i][j][k]+2*tNS[18]*electrolytePotential[i][j][k];
						}
					}
				}
			}
		}
	}
}

void gasLiquidBoundaryElectrolytePotentialD3Q7(int m, int n, int q, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical){
						if (rho1[i][j][k]<=gasCritical){
//							electrolyteIn1[i][j][k]=-Out4[i][j][k]+2*tNSD3Q7[1]*electrolytePotential[i][j][k];
							electrolyteIn1[i][j][k]=Out4[i][j][k];
						}
						if (rho2[i][j][k]<=gasCritical){
//							electrolyteIn2[i][j][k]=-Out5[i][j][k]+2*tNSD3Q7[2]*electrolytePotential[i][j][k];
							electrolyteIn2[i][j][k]=Out5[i][j][k];
						}
						if (rho3[i][j][k]<=gasCritical){
//							electrolyteIn3[i][j][k]=-Out6[i][j][k]+2*tNSD3Q7[3]*electrolytePotential[i][j][k];
							electrolyteIn3[i][j][k]=Out6[i][j][k];
						}
						if (rho4[i][j][k]<=gasCritical){
//							electrolyteIn4[i][j][k]=-Out1[i][j][k]+2*tNSD3Q7[4]*electrolytePotential[i][j][k];
							electrolyteIn4[i][j][k]=Out1[i][j][k];
						}
						if (rho5[i][j][k]<=gasCritical){
//							electrolyteIn5[i][j][k]=-Out2[i][j][k]+2*tNSD3Q7[5]*electrolytePotential[i][j][k];
							electrolyteIn5[i][j][k]=Out2[i][j][k];
						}
						if (rho6[i][j][k]<=gasCritical){
//							electrolyteIn6[i][j][k]=-Out3[i][j][k]+2*tNSD3Q7[6]*electrolytePotential[i][j][k];
							electrolyteIn6[i][j][k]=Out3[i][j][k];
						}
					}
				}
			}
		}
	}
}

void lengthLastNeumannVelocityGhost(int Length,int Width,int Height){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
					fIn0Ghost[Length-1][j][k]=fIn0Ghost[Length-2][j][k];
					fIn1Ghost[Length-1][j][k]=fIn1Ghost[Length-2][j][k];
					fIn2Ghost[Length-1][j][k]=fIn2Ghost[Length-2][j][k];
					fIn3Ghost[Length-1][j][k]=fIn3Ghost[Length-2][j][k];
					fIn4Ghost[Length-1][j][k]=fIn4Ghost[Length-2][j][k];
					fIn5Ghost[Length-1][j][k]=fIn5Ghost[Length-2][j][k];
					fIn6Ghost[Length-1][j][k]=fIn6Ghost[Length-2][j][k];
					fIn7Ghost[Length-1][j][k]=fIn7Ghost[Length-2][j][k];
					fIn8Ghost[Length-1][j][k]=fIn8Ghost[Length-2][j][k];
					fIn9Ghost[Length-1][j][k]=fIn9Ghost[Length-2][j][k];
					fIn10Ghost[Length-1][j][k]=fIn10Ghost[Length-2][j][k];
					fIn11Ghost[Length-1][j][k]=fIn11Ghost[Length-2][j][k];
					fIn12Ghost[Length-1][j][k]=fIn12Ghost[Length-2][j][k];
					fIn13Ghost[Length-1][j][k]=fIn13Ghost[Length-2][j][k];
					fIn14Ghost[Length-1][j][k]=fIn14Ghost[Length-2][j][k];
					fIn15Ghost[Length-1][j][k]=fIn15Ghost[Length-2][j][k];
					fIn16Ghost[Length-1][j][k]=fIn16Ghost[Length-2][j][k];
					fIn17Ghost[Length-1][j][k]=fIn17Ghost[Length-2][j][k];
					fIn18Ghost[Length-1][j][k]=fIn18Ghost[Length-2][j][k];
					uxGhost[Length-1][j][k]=uxGhost[Length-2][j][k];
					uyGhost[Length-1][j][k]=uyGhost[Length-2][j][k];
					uzGhost[Length-1][j][k]=uzGhost[Length-2][j][k];					
			}
		}
	}
}

void lengthFirstNeumannVelocityGhost(int Length,int Width,int Height){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
					fIn0Ghost[0][j][k]=fIn0Ghost[1][j][k];
					fIn1Ghost[0][j][k]=fIn1Ghost[1][j][k];
					fIn2Ghost[0][j][k]=fIn2Ghost[1][j][k];
					fIn3Ghost[0][j][k]=fIn3Ghost[1][j][k];
					fIn4Ghost[0][j][k]=fIn4Ghost[1][j][k];
					fIn5Ghost[0][j][k]=fIn5Ghost[1][j][k];
					fIn6Ghost[0][j][k]=fIn6Ghost[1][j][k];
					fIn7Ghost[0][j][k]=fIn7Ghost[1][j][k];
					fIn8Ghost[0][j][k]=fIn8Ghost[1][j][k];
					fIn9Ghost[0][j][k]=fIn9Ghost[1][j][k];
					fIn10Ghost[0][j][k]=fIn10Ghost[1][j][k];
					fIn11Ghost[0][j][k]=fIn11Ghost[1][j][k];
					fIn12Ghost[0][j][k]=fIn12Ghost[1][j][k];
					fIn13Ghost[0][j][k]=fIn13Ghost[1][j][k];
					fIn14Ghost[0][j][k]=fIn14Ghost[1][j][k];
					fIn15Ghost[0][j][k]=fIn15Ghost[1][j][k];
					fIn16Ghost[0][j][k]=fIn16Ghost[1][j][k];
					fIn17Ghost[0][j][k]=fIn17Ghost[1][j][k];
					fIn18Ghost[0][j][k]=fIn18Ghost[1][j][k];
					uxGhost[0][j][k]=uxGhost[1][j][k];
					uyGhost[0][j][k]=uyGhost[1][j][k];
					uzGhost[0][j][k]=uzGhost[1][j][k];					
			}
		}
	}
}

void widthFirstNeumannVelocityGhost(int Length,int Width,int Height){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (k=0;k<Height;k++){
					fIn0Ghost[i][0][k]=fIn0Ghost[i][1][k];
					fIn1Ghost[i][0][k]=fIn1Ghost[i][1][k];
					fIn2Ghost[i][0][k]=fIn2Ghost[i][1][k];
					fIn3Ghost[i][0][k]=fIn3Ghost[i][1][k];
					fIn4Ghost[i][0][k]=fIn4Ghost[i][1][k];
					fIn5Ghost[i][0][k]=fIn5Ghost[i][1][k];
					fIn6Ghost[i][0][k]=fIn6Ghost[i][1][k];
					fIn7Ghost[i][0][k]=fIn7Ghost[i][1][k];
					fIn8Ghost[i][0][k]=fIn8Ghost[i][1][k];
					fIn9Ghost[i][0][k]=fIn9Ghost[i][1][k];
					fIn10Ghost[i][0][k]=fIn10Ghost[i][1][k];
					fIn11Ghost[i][0][k]=fIn11Ghost[i][1][k];
					fIn12Ghost[i][0][k]=fIn12Ghost[i][1][k];
					fIn13Ghost[i][0][k]=fIn13Ghost[i][1][k];
					fIn14Ghost[i][0][k]=fIn14Ghost[i][1][k];
					fIn15Ghost[i][0][k]=fIn15Ghost[i][1][k];
					fIn16Ghost[i][0][k]=fIn16Ghost[i][1][k];
					fIn17Ghost[i][0][k]=fIn17Ghost[i][1][k];
					fIn18Ghost[i][0][k]=fIn18Ghost[i][1][k];
					uxGhost[i][0][k]=uxGhost[i][1][k];
					uyGhost[i][0][k]=uyGhost[i][1][k];
					uzGhost[i][0][k]=uzGhost[i][1][k];					
			}
		}
	}
}

void widthFirstNeumannConcentration(int Length,int Width,int Height){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (k=0;k<Height;k++){
					cIn0[i][0][k]=cIn0[i][1][k];
					cIn1[i][0][k]=cIn1[i][1][k];
					cIn2[i][0][k]=cIn2[i][1][k];
					cIn3[i][0][k]=cIn3[i][1][k];
					cIn4[i][0][k]=cIn4[i][1][k];
					cIn5[i][0][k]=cIn5[i][1][k];
					cIn6[i][0][k]=cIn6[i][1][k];
					cIn7[i][0][k]=cIn7[i][1][k];
					cIn8[i][0][k]=cIn8[i][1][k];
					cIn9[i][0][k]=cIn9[i][1][k];
					cIn10[i][0][k]=cIn10[i][1][k];
					cIn11[i][0][k]=cIn11[i][1][k];
					cIn12[i][0][k]=cIn12[i][1][k];
					cIn13[i][0][k]=cIn13[i][1][k];
					cIn14[i][0][k]=cIn14[i][1][k];
					cIn15[i][0][k]=cIn15[i][1][k];
					cIn16[i][0][k]=cIn16[i][1][k];
					cIn17[i][0][k]=cIn17[i][1][k];
					cIn18[i][0][k]=cIn18[i][1][k];
					Concentration[i][0][k]=Concentration[i][1][k];	
			}
		}
	}
}

void widthFirstNeumannConcentrationElectrode(int Length,int Width,int Height){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (k=0;k<Height;k++){
					cIn0First[i][0][k]=cIn0First[i][1][k];
					cIn1First[i][0][k]=cIn1First[i][1][k];
					cIn2First[i][0][k]=cIn2First[i][1][k];
					cIn3First[i][0][k]=cIn3First[i][1][k];
					cIn4First[i][0][k]=cIn4First[i][1][k];
					cIn5First[i][0][k]=cIn5First[i][1][k];
					cIn6First[i][0][k]=cIn6First[i][1][k];
					cIn7First[i][0][k]=cIn7First[i][1][k];
					cIn8First[i][0][k]=cIn8First[i][1][k];
					cIn9First[i][0][k]=cIn9First[i][1][k];
					cIn10First[i][0][k]=cIn10First[i][1][k];
					cIn11First[i][0][k]=cIn11First[i][1][k];
					cIn12First[i][0][k]=cIn12First[i][1][k];
					cIn13First[i][0][k]=cIn13First[i][1][k];
					cIn14First[i][0][k]=cIn14First[i][1][k];
					cIn15First[i][0][k]=cIn15First[i][1][k];
					cIn16First[i][0][k]=cIn16First[i][1][k];
					cIn17First[i][0][k]=cIn17First[i][1][k];
					cIn18First[i][0][k]=cIn18First[i][1][k];
					ConcentrationFirst[i][0][k]=ConcentrationFirst[i][1][k];	
			}
		}
	}
}

void widthFirstNeumannConcentrationSecond(int Length,int Width,int Height){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (k=0;k<Height;k++){
					cIn0Second[i][0][k]=cIn0Second[i][1][k];
					cIn1Second[i][0][k]=cIn1Second[i][1][k];
					cIn2Second[i][0][k]=cIn2Second[i][1][k];
					cIn3Second[i][0][k]=cIn3Second[i][1][k];
					cIn4Second[i][0][k]=cIn4Second[i][1][k];
					cIn5Second[i][0][k]=cIn5Second[i][1][k];
					cIn6Second[i][0][k]=cIn6Second[i][1][k];
					cIn7Second[i][0][k]=cIn7Second[i][1][k];
					cIn8Second[i][0][k]=cIn8Second[i][1][k];
					cIn9Second[i][0][k]=cIn9Second[i][1][k];
					cIn10Second[i][0][k]=cIn10Second[i][1][k];
					cIn11Second[i][0][k]=cIn11Second[i][1][k];
					cIn12Second[i][0][k]=cIn12Second[i][1][k];
					cIn13Second[i][0][k]=cIn13Second[i][1][k];
					cIn14Second[i][0][k]=cIn14Second[i][1][k];
					cIn15Second[i][0][k]=cIn15Second[i][1][k];
					cIn16Second[i][0][k]=cIn16Second[i][1][k];
					cIn17Second[i][0][k]=cIn17Second[i][1][k];
					cIn18Second[i][0][k]=cIn18Second[i][1][k];
					ConcentrationSecond[i][0][k]=ConcentrationSecond[i][1][k];	
			}
		}
	}
}

void HeightLastConstantVelocityBoundaryGhost(int q,int mInitial,int m,int nInitial,int n,double InletSideUz,double InletSideUx,double InletSideUy,double solid,double b, double R, double T, double a){
	int i;
	int j;
	int k;
	double bRho4;
	double BRho4;
#pragma omp parallel private(i,j,bRho4,BRho4)
	{
#pragma omp for schedule(dynamic)
		for (i=mInitial;i<m;i++){
			for (j=nInitial;j<n;j++){
				if (domain[i][j][q-2]==0){
					rhoGhost[i][j][q-2]=((fIn0Ghost[i][j][q-2]+fIn1Ghost[i][j][q-2]+fIn2Ghost[i][j][q-2]+fIn4Ghost[i][j][q-2]+fIn5Ghost[i][j][q-2]+fIn7Ghost[i][j][q-2]+fIn10Ghost[i][j][q-2]+fIn13Ghost[i][j][q-2]+fIn16Ghost[i][j][q-2])
						+2.0*(fIn3Ghost[i][j][q-2]+fIn8Ghost[i][j][q-2]+fIn9Ghost[i][j][q-2]+fIn17Ghost[i][j][q-2]+fIn18Ghost[i][j][q-2]))/(1+InletSideUz);
					uxGhost[i][j][q-2]=InletSideUx;
					uyGhost[i][j][q-2]=InletSideUy;
					uzGhost[i][j][q-2]=InletSideUz;
					fIn6Ghost[i][j][q-2]=fIn3Ghost[i][j][q-2]-1.0/3*rhoGhost[i][j][q-2]*InletSideUz;
					fIn11Ghost[i][j][q-2]=fIn8Ghost[i][j][q-2]-1.0/6*rhoGhost[i][j][q-2]*(InletSideUx+InletSideUz)+0.5*(fIn1Ghost[i][j][q-2]+fIn7Ghost[i][j][q-2]+fIn13Ghost[i][j][q-2]-(fIn4Ghost[i][j][q-2]+fIn10Ghost[i][j][q-2]+fIn16Ghost[i][j][q-2]))-1.0/3*rhoGhost[i][j][q-2]*InletSideUx;
					fIn12Ghost[i][j][q-2]=fIn9Ghost[i][j][q-2]-1.0/6*rhoGhost[i][j][q-2]*(InletSideUy+InletSideUz)+0.5*(fIn2Ghost[i][j][q-2]+fIn7Ghost[i][j][q-2]+fIn16Ghost[i][j][q-2]-(fIn5Ghost[i][j][q-2]+fIn10Ghost[i][j][q-2]+fIn13Ghost[i][j][q-2]))-1.0/3*rhoGhost[i][j][q-2]*InletSideUy;
					fIn14Ghost[i][j][q-2]=fIn17Ghost[i][j][q-2]+1.0/6*rhoGhost[i][j][q-2]*(InletSideUx-InletSideUz)-(0.5*(fIn1Ghost[i][j][q-2]+fIn7Ghost[i][j][q-2]+fIn13Ghost[i][j][q-2]-(fIn4Ghost[i][j][q-2]+fIn10Ghost[i][j][q-2]+fIn16Ghost[i][j][q-2]))-1.0/3*rhoGhost[i][j][q-2]*InletSideUx);
					fIn15Ghost[i][j][q-2]=fIn18Ghost[i][j][q-2]+1.0/6*rhoGhost[i][j][q-2]*(InletSideUy-InletSideUz)-(0.5*(fIn2Ghost[i][j][q-2]+fIn7Ghost[i][j][q-2]+fIn16Ghost[i][j][q-2]-(fIn5Ghost[i][j][q-2]+fIn10Ghost[i][j][q-2]+fIn13Ghost[i][j][q-2]))-1.0/3*rhoGhost[i][j][q-2]*InletSideUy);
				}
				if (domain[i][j][q-2]!=0){
					rhoGhost[i][j][q-2]=solid;
					uxGhost[i][j][q-2]=0.0;
					uyGhost[i][j][q-2]=0.0;
					uzGhost[i][j][q-2]=0.0;
				}
				rhoGhost[i][j][q-1]=rhoGhost[i][j][q-2];
				uxGhost[i][j][q-1]=uxGhost[i][j][q-2];
				uyGhost[i][j][q-1]=uyGhost[i][j][q-2];
				uzGhost[i][j][q-1]=uzGhost[i][j][q-2];
			}
		}
	}
}

void WidthLastConstantVelocityBoundaryGhost(int q, int mInitial, int m, int nInitial, int n, double InletSideUz, double InletSideUx, double InletSideUy, double solid){
	int i;
	int j;
	int k;
	double bRho4;
	double BRho4;
#pragma omp parallel private(i,k,bRho4,BRho4)
	{
#pragma omp for schedule(dynamic)
		for (i = mInitial; i<m; i++){
			for (k = nInitial; k<n; k++){
				if (domain[i][q-2][k] == 0){
					rhoGhost[i][q-2][k] = ((fIn0Ghost[i][q-2][k] + fIn1Ghost[i][q-2][k] + fIn3Ghost[i][q-2][k] + fIn4Ghost[i][q-2][k] + fIn6Ghost[i][q-2][k] + fIn8Ghost[i][q-2][k] + fIn11Ghost[i][q-2][k] + fIn14Ghost[i][q-2][k] + fIn17Ghost[i][q-2][k])
						+ 2.0*(fIn2Ghost[i][q-2][k] + fIn7Ghost[i][q-2][k] + fIn9Ghost[i][q-2][k] + fIn15Ghost[i][q-2][k] + fIn16Ghost[i][q-2][k])) / (1 + InletSideUy);
					uxGhost[i][q-2][k] = InletSideUx;
					uyGhost[i][q-2][k] = InletSideUy;
					uzGhost[i][q-2][k] = InletSideUz;
					fIn5Ghost[i][q-2][k] = fIn2Ghost[i][q-2][k] - 1.0 / 3 * rhoGhost[i][q-2][k] * InletSideUy;
					fIn10Ghost[i][q-2][k] = fIn7Ghost[i][q-2][k] - 1.0 / 6 * rhoGhost[i][q-2][k] * (InletSideUx + InletSideUy);
					fIn12Ghost[i][q-2][k] = fIn9Ghost[i][q-2][k] - 1.0 / 6 * rhoGhost[i][q-2][k] * (InletSideUy + InletSideUz);
					fIn13Ghost[i][q-2][k] = fIn16Ghost[i][q-2][k] + 1.0 / 6 * rhoGhost[i][q-2][k] * (InletSideUx - InletSideUy);
					fIn18Ghost[i][q-2][k] = fIn15Ghost[i][q-2][k] - 1.0 / 6 * rhoGhost[i][q-2][k] * (InletSideUy - InletSideUz);
				}
				if (domain[i][q-2][k] != 0){
					rhoGhost[i][q-2][k] = solid;
					uxGhost[i][q-2][k] = 0.0;
					uyGhost[i][q-2][k] = 0.0;
					uzGhost[i][q-2][k] = 0.0;
				}
				rhoGhost[i][q-1][k] = rhoGhost[i][q-2][k];
				uxGhost[i][q-1][k] = uxGhost[i][q-2][k];
				uyGhost[i][q-1][k] = uyGhost[i][q-2][k];
				uzGhost[i][q-1][k] = uzGhost[i][q-2][k];
			}
		}
	}
}

void WidthLastConstantPressureBoundaryGhost(int q, int mInitial, int m, int nInitial, int n, double outPutDensity, double InletSideUx, double InletSideUz, double solid, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i = mInitial; i<m; i++){
			for (k = nInitial; k<n; k++){
				if (domain[i][q-2][k] == 0 && rho[i][q-2][k]>gasCritical){
					rhoGhost[i][q - 2][k] = outPutDensity;
					uxGhost[i][q-2][k] = InletSideUx;
					uyGhost[i][q - 2][k] = ((fIn0Ghost[i][q - 2][k] + fIn1Ghost[i][q - 2][k] + fIn3Ghost[i][q - 2][k] + fIn4Ghost[i][q - 2][k] + fIn6Ghost[i][q - 2][k] + fIn8Ghost[i][q - 2][k] + fIn11Ghost[i][q - 2][k] + fIn14Ghost[i][q - 2][k] + fIn17Ghost[i][q - 2][k])
						+ 2.0*(fIn2Ghost[i][q - 2][k] + fIn7Ghost[i][q - 2][k] + fIn9Ghost[i][q - 2][k] + fIn15Ghost[i][q - 2][k] + fIn16Ghost[i][q - 2][k])) / outPutDensity - 1;
					uzGhost[i][q-2][k] = InletSideUz;
					fIn5Ghost[i][q - 2][k] = fIn2Ghost[i][q - 2][k] - 1.0 / 3 * rhoGhost[i][q - 2][k] * uyGhost[i][q - 2][k];
					fIn10Ghost[i][q - 2][k] = fIn7Ghost[i][q - 2][k] - 1.0 / 6 * rhoGhost[i][q - 2][k] * (InletSideUx + uyGhost[i][q - 2][k]);
					fIn12Ghost[i][q - 2][k] = fIn9Ghost[i][q - 2][k] - 1.0 / 6 * rhoGhost[i][q - 2][k] * (uyGhost[i][q - 2][k] + InletSideUz);
					fIn13Ghost[i][q - 2][k] = fIn16Ghost[i][q - 2][k] + 1.0 / 6 * rhoGhost[i][q - 2][k] * (InletSideUx - uyGhost[i][q - 2][k]);
					fIn18Ghost[i][q - 2][k] = fIn15Ghost[i][q - 2][k] - 1.0 / 6 * rhoGhost[i][q - 2][k] * (uyGhost[i][q - 2][k] - InletSideUz);
				}
				if (domain[i][q-2][k] != 0 || rho[i][q-2][k] <= gasCritical){
					rhoGhost[i][q-2][k] = solid;
					uxGhost[i][q-2][k] = 0.0;
					uyGhost[i][q-2][k] = 0.0;
					uzGhost[i][q-2][k] = 0.0;
				}
				rhoGhost[i][q-1][k] = rhoGhost[i][q-2][k];
				uxGhost[i][q-1][k] = uxGhost[i][q-2][k];
				uyGhost[i][q-1][k] = uyGhost[i][q-2][k];
				uzGhost[i][q-1][k] = uzGhost[i][q-2][k];
			}
		}
	}
}

void HeightFirstConstantVelocityBoundaryGhost(int mInitial, int m, int nInitial, int n, double InletSideUz, double InletSideUx, double InletSideUy, double solid, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i = mInitial; i<m; i++){
			for (j = nInitial; j<n; j++){
				if (domain[i][j][1] == 0 && rho[i][j][1]>gasCritical){
					rhoGhost[i][j][1] = ((fIn0Ghost[i][j][1] + fIn1Ghost[i][j][1] + fIn2Ghost[i][j][1] + fIn4Ghost[i][j][1] + fIn5Ghost[i][j][1] + fIn7Ghost[i][j][1] + fIn10Ghost[i][j][1] + fIn13Ghost[i][j][1] + fIn16Ghost[i][j][1])
						+ 2.0*(fIn6Ghost[i][j][1] + fIn11Ghost[i][j][1] + fIn12Ghost[i][j][1] + fIn14Ghost[i][j][1] + fIn15Ghost[i][j][1])) / (1 - InletSideUz);
					uxGhost[i][j][1] = InletSideUx;
					uyGhost[i][j][1] = InletSideUy;
					uzGhost[i][j][1] = InletSideUz;
					fIn3Ghost[i][j][1] = fIn6Ghost[i][j][1] + 1.0 / 3 * rhoGhost[i][j][1] * uzGhost[i][j][1];
//					fIn8Ghost[i][j][1]=fIn11Ghost[i][j][1]+1.0/6*rhoGhost[i][j][1]*(InletSideUx+uzGhost[i][j][1])-0.5*(fIn1Ghost[i][j][1]+fIn7Ghost[i][j][1]+fIn13Ghost[i][j][1]-(fIn4Ghost[i][j][1]+fIn10Ghost[i][j][1]+fIn16Ghost[i][j][1]))+1.0/3*rhoGhost[i][j][1]*InletSideUx;
					fIn8Ghost[i][j][1] = fIn11Ghost[i][j][1] + 1.0 / 6 * rhoGhost[i][j][1] * uzGhost[i][j][1];
//					fIn9Ghost[i][j][1]=fIn12Ghost[i][j][1]+1.0/6*rhoGhost[i][j][1]*(InletSideUy+uzGhost[i][j][1])-0.5*(fIn2Ghost[i][j][1]+fIn7Ghost[i][j][1]+fIn16Ghost[i][j][1]-(fIn5Ghost[i][j][1]+fIn10Ghost[i][j][1]+fIn13Ghost[i][j][1]))+1.0/3*rhoGhost[i][j][1]*InletSideUy;
					fIn9Ghost[i][j][1] = fIn12Ghost[i][j][1] + 1.0 / 6 * rhoGhost[i][j][1] * uzGhost[i][j][1];
//					fIn17Ghost[i][j][1]=fIn14Ghost[i][j][1]-1.0/6*rhoGhost[i][j][1]*(InletSideUx-uzGhost[i][j][1])+(0.5*(fIn1Ghost[i][j][1]+fIn7Ghost[i][j][1]+fIn13Ghost[i][j][1]-(fIn4Ghost[i][j][1]+fIn10Ghost[i][j][1]+fIn16Ghost[i][j][1]))-1.0/3*rhoGhost[i][j][1]*InletSideUx);
					fIn17Ghost[i][j][1] = fIn14Ghost[i][j][1] + 1.0 / 6 * rhoGhost[i][j][1] * uzGhost[i][j][1];
//					fIn18Ghost[i][j][1]=fIn15Ghost[i][j][1]-1.0/6*rhoGhost[i][j][1]*(InletSideUy-uzGhost[i][j][1])+(0.5*(fIn2Ghost[i][j][1]+fIn7Ghost[i][j][1]+fIn16Ghost[i][j][1]-(fIn5Ghost[i][j][1]+fIn10Ghost[i][j][1]+fIn13Ghost[i][j][1]))-1.0/3*rhoGhost[i][j][1]*InletSideUy);
					fIn18Ghost[i][j][1] = fIn15Ghost[i][j][1] + 1.0 / 6 * rhoGhost[i][j][1] * uzGhost[i][j][1];
				}
				if (domain[i][j][1] != 0 || rho[i][j][1] <= gasCritical){
					rhoGhost[i][j][1] = solid;
					uxGhost[i][j][1] = 0.0;
					uyGhost[i][j][1] = 0.0;
					uzGhost[i][j][1] = 0.0;
				}
				rhoGhost[i][j][0] = rhoGhost[i][j][1];
				uxGhost[i][j][0] = uxGhost[i][j][1];
				uyGhost[i][j][0] = uyGhost[i][j][1];
				uzGhost[i][j][0] = uzGhost[i][j][1];
			}
		}
	}
}

void HeightFirstConstantPressureBoundaryGhost(int mInitial, int m, int nInitial, int n, double InletDensity, double InletSideUx, double InletSideUy, double solid, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=mInitial;i<m;i++){
			for (j=nInitial;j<n;j++){
				if (domain[i][j][1]==0&&rho[i][j][1]>gasCritical){
					rhoGhost[i][j][1]=InletDensity;
					uxGhost[i][j][1]=InletSideUx;
					uyGhost[i][j][1]=InletSideUy;
					uzGhost[i][j][1]=1-((fIn0Ghost[i][j][1]+fIn1Ghost[i][j][1]+fIn2Ghost[i][j][1]+fIn4Ghost[i][j][1]+fIn5Ghost[i][j][1]+fIn7Ghost[i][j][1]+fIn10Ghost[i][j][1]+fIn13Ghost[i][j][1]+fIn16Ghost[i][j][1])
						+2.0*(fIn6Ghost[i][j][1]+fIn11Ghost[i][j][1]+fIn12Ghost[i][j][1]+fIn14Ghost[i][j][1]+fIn15Ghost[i][j][1]))/InletDensity;
					fIn3Ghost[i][j][1]=fIn6Ghost[i][j][1]+1.0/3*rhoGhost[i][j][1]*uzGhost[i][j][1];
//					fIn8Ghost[i][j][1]=fIn11Ghost[i][j][1]+1.0/6*rhoGhost[i][j][1]*(InletSideUx+uzGhost[i][j][1])-0.5*(fIn1Ghost[i][j][1]+fIn7Ghost[i][j][1]+fIn13Ghost[i][j][1]-(fIn4Ghost[i][j][1]+fIn10Ghost[i][j][1]+fIn16Ghost[i][j][1]))+1.0/3*rhoGhost[i][j][1]*InletSideUx;
					fIn8Ghost[i][j][1]=fIn11Ghost[i][j][1]+1.0/6*rhoGhost[i][j][1]*uzGhost[i][j][1];
//					fIn9Ghost[i][j][1]=fIn12Ghost[i][j][1]+1.0/6*rhoGhost[i][j][1]*(InletSideUy+uzGhost[i][j][1])-0.5*(fIn2Ghost[i][j][1]+fIn7Ghost[i][j][1]+fIn16Ghost[i][j][1]-(fIn5Ghost[i][j][1]+fIn10Ghost[i][j][1]+fIn13Ghost[i][j][1]))+1.0/3*rhoGhost[i][j][1]*InletSideUy;
					fIn9Ghost[i][j][1]=fIn12Ghost[i][j][1]+1.0/6*rhoGhost[i][j][1]*uzGhost[i][j][1];
//					fIn17Ghost[i][j][1]=fIn14Ghost[i][j][1]-1.0/6*rhoGhost[i][j][1]*(InletSideUx-uzGhost[i][j][1])+(0.5*(fIn1Ghost[i][j][1]+fIn7Ghost[i][j][1]+fIn13Ghost[i][j][1]-(fIn4Ghost[i][j][1]+fIn10Ghost[i][j][1]+fIn16Ghost[i][j][1]))-1.0/3*rhoGhost[i][j][1]*InletSideUx);
					fIn17Ghost[i][j][1]=fIn14Ghost[i][j][1]+1.0/6*rhoGhost[i][j][1]*uzGhost[i][j][1];
//					fIn18Ghost[i][j][1]=fIn15Ghost[i][j][1]-1.0/6*rhoGhost[i][j][1]*(InletSideUy-uzGhost[i][j][1])+(0.5*(fIn2Ghost[i][j][1]+fIn7Ghost[i][j][1]+fIn16Ghost[i][j][1]-(fIn5Ghost[i][j][1]+fIn10Ghost[i][j][1]+fIn13Ghost[i][j][1]))-1.0/3*rhoGhost[i][j][1]*InletSideUy);
					fIn18Ghost[i][j][1]=fIn15Ghost[i][j][1]+1.0/6*rhoGhost[i][j][1]*uzGhost[i][j][1];
				}
				if (domain[i][j][1]!=0||rho[i][j][1]<=gasCritical){
					rhoGhost[i][j][1]=solid;
					uxGhost[i][j][1]=0.0;
					uyGhost[i][j][1]=0.0;
					uzGhost[i][j][1]=0.0;
				}
				rhoGhost[i][j][0]=rhoGhost[i][j][1];
				uxGhost[i][j][0]=uxGhost[i][j][1];
				uyGhost[i][j][0]=uyGhost[i][j][1];
				uzGhost[i][j][0]=uzGhost[i][j][1];
			}
		}
	}
}

void HeightLastConstantPressureBoundaryGhost(int q, int mInitial, int m, int nInitial, int n, double outPutDensity, double InletSideUx, double InletSideUy, double solid, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=mInitial;i<m;i++){
			for (j=nInitial;j<n;j++){
				if (domain[i][j][q-2]==0&&rho[i][j][q-2]>gasCritical){
					uxGhost[i][j][q-2]=InletSideUx;
					uyGhost[i][j][q-2]=InletSideUy;
					uzGhost[i][j][q-2]=((fIn0Ghost[i][j][q-2]+fIn1Ghost[i][j][q-2]+fIn2Ghost[i][j][q-2]+fIn4Ghost[i][j][q-2]+fIn5Ghost[i][j][q-2]+fIn7Ghost[i][j][q-2]+fIn10Ghost[i][j][q-2]+fIn13Ghost[i][j][q-2]+fIn16Ghost[i][j][q-2])
						+2.0*(fIn3Ghost[i][j][q-2]+fIn8Ghost[i][j][q-2]+fIn9Ghost[i][j][q-2]+fIn17Ghost[i][j][q-2]+fIn18Ghost[i][j][q-2]))/outPutDensity-1;
					fIn6Ghost[i][j][q-2]=fIn3Ghost[i][j][q-2]-1.0/3*rhoGhost[i][j][q-2]*uzGhost[i][j][q-2];
					fIn11Ghost[i][j][q-2]=fIn8Ghost[i][j][q-2]-1.0/6*rhoGhost[i][j][q-2]*(InletSideUx+uzGhost[i][j][q-2])+0.5*(fIn1Ghost[i][j][q-2]+fIn7Ghost[i][j][q-2]+fIn13Ghost[i][j][q-2]-(fIn4Ghost[i][j][q-2]+fIn10Ghost[i][j][q-2]+fIn16Ghost[i][j][q-2]))-1.0/3*rhoGhost[i][j][q-2]*InletSideUx;
					fIn12Ghost[i][j][q-2]=fIn9Ghost[i][j][q-2]-1.0/6*rhoGhost[i][j][q-2]*(InletSideUy+uzGhost[i][j][q-2])+0.5*(fIn2Ghost[i][j][q-2]+fIn7Ghost[i][j][q-2]+fIn16Ghost[i][j][q-2]-(fIn5Ghost[i][j][q-2]+fIn10Ghost[i][j][q-2]+fIn13Ghost[i][j][q-2]))-1.0/3*rhoGhost[i][j][q-2]*InletSideUy;
					fIn14Ghost[i][j][q-2]=fIn17Ghost[i][j][q-2]+1.0/6*rhoGhost[i][j][q-2]*(InletSideUx-uzGhost[i][j][q-2])-(0.5*(fIn1Ghost[i][j][q-2]+fIn7Ghost[i][j][q-2]+fIn13Ghost[i][j][q-2]-(fIn4Ghost[i][j][q-2]+fIn10Ghost[i][j][q-2]+fIn16Ghost[i][j][q-2]))-1.0/3*rhoGhost[i][j][q-2]*InletSideUx);
					fIn15Ghost[i][j][q-2]=fIn18Ghost[i][j][q-2]+1.0/6*rhoGhost[i][j][q-2]*(InletSideUy-uzGhost[i][j][q-2])-(0.5*(fIn2Ghost[i][j][q-2]+fIn7Ghost[i][j][q-2]+fIn16Ghost[i][j][q-2]-(fIn5Ghost[i][j][q-2]+fIn10Ghost[i][j][q-2]+fIn13Ghost[i][j][q-2]))-1.0/3*rhoGhost[i][j][q-2]*InletSideUy);
				}
				if (domain[i][j][q-2]!=0||rho[i][j][q-2]<=gasCritical){
					rhoGhost[i][j][q-2]=solid;
					uxGhost[i][j][q-2]=0.0;
					uyGhost[i][j][q-2]=0.0;
					uzGhost[i][j][q-2]=0.0;
				}
				rhoGhost[i][j][q-1]=rhoGhost[i][j][q-2];
				uxGhost[i][j][q-1]=uxGhost[i][j][q-2];
				uyGhost[i][j][q-1]=uyGhost[i][j][q-2];
				uzGhost[i][j][q-1]=uzGhost[i][j][q-2];
			}
		}
	}
}



void HeightFirstDirichletConcentration(int Length,int Width,double InletConcentration,double gasCritical,double constant){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				if (domain[i][j][1]==0&&rho[i][j][1]>gasCritical){
					Concentration[i][j][1]=InletConcentration;
					cIn3[i][j][1]=-Out6[i][j][1]+2*tNS[3]*Concentration[i][j][1]*(1+4.5*((cxNS[3]*uxGhost[i][j][1]*constant+cyNS[3]*uyGhost[i][j][1]*constant+czNS[3]*uzGhost[i][j][1]*constant)*(cxNS[3]*uxGhost[i][j][1]*constant+cyNS[3]*uyGhost[i][j][1]*constant+czNS[3]*uzGhost[i][j][1]*constant))-1.5*(uxGhost[i][j][1]*uxGhost[i][j][1]*constant*constant+uyGhost[i][j][1]*uyGhost[i][j][1]*constant*constant+uzGhost[i][j][1]*uzGhost[i][j][1]*constant*constant));
					cIn8[i][j][1]=-Out11[i][j][1]+2*tNS[8]*Concentration[i][j][1]*(1+4.5*((cxNS[8]*uxGhost[i][j][1]*constant+cyNS[8]*uyGhost[i][j][1]*constant+czNS[8]*uzGhost[i][j][1]*constant)*(cxNS[8]*uxGhost[i][j][1]*constant+cyNS[8]*uyGhost[i][j][1]*constant+czNS[8]*uzGhost[i][j][1]*constant))-1.5*(uxGhost[i][j][1]*uxGhost[i][j][1]*constant*constant+uyGhost[i][j][1]*uyGhost[i][j][1]*constant*constant+uzGhost[i][j][1]*uzGhost[i][j][1]*constant*constant));
					cIn9[i][j][1]=-Out12[i][j][1]+2*tNS[9]*Concentration[i][j][1]*(1+4.5*((cxNS[9]*uxGhost[i][j][1]*constant+cyNS[9]*uyGhost[i][j][1]*constant+czNS[9]*uzGhost[i][j][1]*constant)*(cxNS[9]*uxGhost[i][j][1]*constant+cyNS[9]*uyGhost[i][j][1]*constant+czNS[9]*uzGhost[i][j][1]*constant))-1.5*(uxGhost[i][j][1]*uxGhost[i][j][1]*constant*constant+uyGhost[i][j][1]*uyGhost[i][j][1]*constant*constant+uzGhost[i][j][1]*uzGhost[i][j][1]*constant*constant));
					cIn17[i][j][1]=-Out14[i][j][1]+2*tNS[17]*Concentration[i][j][1]*(1+4.5*((cxNS[17]*uxGhost[i][j][1]*constant+cyNS[17]*uyGhost[i][j][1]*constant+czNS[17]*uzGhost[i][j][1]*constant)*(cxNS[17]*uxGhost[i][j][1]*constant+cyNS[17]*uyGhost[i][j][1]*constant+czNS[17]*uzGhost[i][j][1]*constant))-1.5*(uxGhost[i][j][1]*uxGhost[i][j][1]*constant*constant+uyGhost[i][j][1]*uyGhost[i][j][1]*constant*constant+uzGhost[i][j][1]*uzGhost[i][j][1]*constant*constant));
					cIn18[i][j][1]=-Out15[i][j][1]+2*tNS[18]*Concentration[i][j][1]*(1+4.5*((cxNS[18]*uxGhost[i][j][1]*constant+cyNS[18]*uyGhost[i][j][1]*constant+czNS[18]*uzGhost[i][j][1]*constant)*(cxNS[18]*uxGhost[i][j][1]*constant+cyNS[18]*uyGhost[i][j][1]*constant+czNS[18]*uzGhost[i][j][1]*constant))-1.5*(uxGhost[i][j][1]*uxGhost[i][j][1]*constant*constant+uyGhost[i][j][1]*uyGhost[i][j][1]*constant*constant+uzGhost[i][j][1]*uzGhost[i][j][1]*constant*constant));
				}
				if (domain[i][j][1]!=0||rho[i][j][1]<=gasCritical){
					Concentration[i][j][1]=0.0;
				}
				Concentration[i][j][0]=Concentration[i][j][1];
				cIn3[i][j][0]=cIn3[i][j][1];
				cIn8[i][j][0]=cIn8[i][j][1];
				cIn9[i][j][0]=cIn9[i][j][1];
				cIn17[i][j][0]=cIn17[i][j][1];
				cIn18[i][j][0]=cIn18[i][j][1];
			}
		}
	}
}

void HeightFirstDirichletConcentrationD3Q7(int Length,int Width,double InletConcentration,double gasCritical,double constant){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				if (domain[i][j][1]==0&&rho[i][j][1]>gasCritical){
					Concentration[i][j][1]=InletConcentration;
//					cIn3[i][j][1]=-Out6[i][j][1]+2*tNSD3Q7[3]*Concentration[i][j][1];
					cIn0[i][j][1] = Concentration[i][j][1] * tNSD3Q7[0];
					cIn1[i][j][1] = Concentration[i][j][1] * tNSD3Q7[1];
					cIn2[i][j][1] = Concentration[i][j][1] * tNSD3Q7[2];
					cIn3[i][j][1] = Concentration[i][j][1] * tNSD3Q7[3];
					cIn4[i][j][1] = Concentration[i][j][1] * tNSD3Q7[4];
					cIn5[i][j][1] = Concentration[i][j][1] * tNSD3Q7[5];
					cIn6[i][j][1] = Concentration[i][j][1] * tNSD3Q7[6];
				}
				if (domain[i][j][1]!=0||rho[i][j][1]<=gasCritical){
					Concentration[i][j][1]=0.0;
				}
				Concentration[i][j][0]=Concentration[i][j][1];
				cIn0[i][j][0] = cIn0[i][j][1];
				cIn1[i][j][0] = cIn1[i][j][1];
				cIn2[i][j][0] = cIn2[i][j][1];
				cIn3[i][j][0] = cIn3[i][j][1];
				cIn4[i][j][0] = cIn4[i][j][1];
				cIn5[i][j][0] = cIn5[i][j][1];
				cIn6[i][j][0] = cIn6[i][j][1];
			}
		}
	}
}

void HeightLastDirichletConcentrationD3Q7(int Length, int Width, double InletConcentration, double gasCritical, int constant){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (j = 0; j<Width; j++){
				if (domain[i][j][constant] == 0 && rho[i][j][constant]>gasCritical){
					Concentration[i][j][constant] = InletConcentration;
					cIn0[i][j][constant] = Concentration[i][j][constant] * tNSD3Q7[0];
					cIn1[i][j][constant] = Concentration[i][j][constant] * tNSD3Q7[1];
					cIn2[i][j][constant] = Concentration[i][j][constant] * tNSD3Q7[2];
					cIn3[i][j][constant] = Concentration[i][j][constant] * tNSD3Q7[3];
					cIn4[i][j][constant] = Concentration[i][j][constant] * tNSD3Q7[4];
					cIn5[i][j][constant] = Concentration[i][j][constant] * tNSD3Q7[5];
					cIn6[i][j][constant] = Concentration[i][j][constant] * tNSD3Q7[6];
				}
				if (domain[i][j][constant] != 0 || rho[i][j][constant] <= gasCritical){
					Concentration[i][j][constant] = 0.0;
				}
				Concentration[i][j][constant + 1] = Concentration[i][j][constant];
				cIn0[i][j][constant+1] = cIn0[i][j][constant];
				cIn1[i][j][constant+1] = cIn1[i][j][constant];
				cIn2[i][j][constant+1] = cIn2[i][j][constant];
				cIn3[i][j][constant+1] = cIn3[i][j][constant];
				cIn4[i][j][constant+1] = cIn4[i][j][constant];
				cIn5[i][j][constant+1] = cIn5[i][j][constant];
				cIn6[i][j][constant+1] = cIn6[i][j][constant];
			}
		}
	}
}

void WidthLastDirichletConcentrationD3Q7(int Length, int initialHeight, int Height, double InletConcentration, double gasCritical, int constant){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (k = initialHeight; k<Height; k++){
				if (domain[i][constant][k] == 0 && rho[i][constant][k]>gasCritical){
					Concentration[i][constant][k] = InletConcentration;
					cIn0[i][constant][k] = Concentration[i][constant][k] * tNSD3Q7[0];
					cIn1[i][constant][k] = Concentration[i][constant][k] * tNSD3Q7[1];
					cIn2[i][constant][k] = Concentration[i][constant][k] * tNSD3Q7[2];
					cIn3[i][constant][k] = Concentration[i][constant][k] * tNSD3Q7[3];
					cIn4[i][constant][k] = Concentration[i][constant][k] * tNSD3Q7[4];
					cIn5[i][constant][k] = Concentration[i][constant][k] * tNSD3Q7[5];
					cIn6[i][constant][k] = Concentration[i][constant][k] * tNSD3Q7[6];
				}
				if (domain[i][constant][k] != 0 || rho[i][constant][k] <= gasCritical){
					Concentration[i][constant][k] = 0.0;
				}
				Concentration[i][constant+1][k] = Concentration[i][constant][k];
				cIn0[i][constant+1][k] = cIn0[i][constant][k];
				cIn1[i][constant+1][k] = cIn1[i][constant][k];
				cIn2[i][constant+1][k] = cIn2[i][constant][k];
				cIn3[i][constant+1][k] = cIn3[i][constant][k];
				cIn4[i][constant+1][k] = cIn4[i][constant][k];
				cIn5[i][constant+1][k] = cIn5[i][constant][k];
				cIn6[i][constant+1][k] = cIn6[i][constant][k];
			}
		}
	}
}

void WidthLastDirichletConcentrationFirstD3Q7(int Length, int initialHeight, int Height, double InletConcentration, double gasCritical, int constant){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (k = initialHeight; k<Height; k++){
				if (domain[i][constant][k] == 0 && rho[i][constant][k]>gasCritical){
					ConcentrationFirst[i][constant][k] = InletConcentration;
					cIn0First[i][constant][k] = ConcentrationFirst[i][constant][k] * tNSD3Q7[0];
					cIn1First[i][constant][k] = ConcentrationFirst[i][constant][k] * tNSD3Q7[1];
					cIn2First[i][constant][k] = ConcentrationFirst[i][constant][k] * tNSD3Q7[2];
					cIn3First[i][constant][k] = ConcentrationFirst[i][constant][k] * tNSD3Q7[3];
					cIn4First[i][constant][k] = ConcentrationFirst[i][constant][k] * tNSD3Q7[4];
					cIn5First[i][constant][k] = ConcentrationFirst[i][constant][k] * tNSD3Q7[5];
					cIn6First[i][constant][k] = ConcentrationFirst[i][constant][k] * tNSD3Q7[6];
				}
				if (domain[i][constant][k] != 0 || rho[i][constant][k] <= gasCritical){
					ConcentrationFirst[i][constant][k] = 0.0;
				}
				ConcentrationFirst[i][constant + 1][k] = ConcentrationFirst[i][constant][k];
				cIn0First[i][constant + 1][k] = cIn0First[i][constant][k];
				cIn1First[i][constant + 1][k] = cIn1First[i][constant][k];
				cIn2First[i][constant + 1][k] = cIn2First[i][constant][k];
				cIn3First[i][constant + 1][k] = cIn3First[i][constant][k];
				cIn4First[i][constant + 1][k] = cIn4First[i][constant][k];
				cIn5First[i][constant + 1][k] = cIn5First[i][constant][k];
				cIn6First[i][constant + 1][k] = cIn6First[i][constant][k];
			}
		}
	}
}

void HeightFirstDirichletConcentrationFirst(int Length,int Width,double InletConcentration,double gasCritical,double constant){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				if (domain[i][j][1]==0&&rho[i][j][1]>gasCritical){
					ConcentrationFirst[i][j][1]=InletConcentration;
					cIn3First[i][j][1]=-Out6[i][j][1]+2*tNS[3]*ConcentrationFirst[i][j][1]*(1+4.5*((cxNS[3]*uxGhost[i][j][1]*constant+cyNS[3]*uyGhost[i][j][1]*constant+czNS[3]*uzGhost[i][j][1]*constant)*(cxNS[3]*uxGhost[i][j][1]*constant+cyNS[3]*uyGhost[i][j][1]*constant+czNS[3]*uzGhost[i][j][1]*constant))-1.5*(uxGhost[i][j][1]*uxGhost[i][j][1]*constant*constant+uyGhost[i][j][1]*uyGhost[i][j][1]*constant*constant+uzGhost[i][j][1]*uzGhost[i][j][1]*constant*constant));
					cIn8First[i][j][1]=-Out11[i][j][1]+2*tNS[8]*ConcentrationFirst[i][j][1]*(1+4.5*((cxNS[8]*uxGhost[i][j][1]*constant+cyNS[8]*uyGhost[i][j][1]*constant+czNS[8]*uzGhost[i][j][1]*constant)*(cxNS[8]*uxGhost[i][j][1]*constant+cyNS[8]*uyGhost[i][j][1]*constant+czNS[8]*uzGhost[i][j][1]*constant))-1.5*(uxGhost[i][j][1]*uxGhost[i][j][1]*constant*constant+uyGhost[i][j][1]*uyGhost[i][j][1]*constant*constant+uzGhost[i][j][1]*uzGhost[i][j][1]*constant*constant));
					cIn9First[i][j][1]=-Out12[i][j][1]+2*tNS[9]*ConcentrationFirst[i][j][1]*(1+4.5*((cxNS[9]*uxGhost[i][j][1]*constant+cyNS[9]*uyGhost[i][j][1]*constant+czNS[9]*uzGhost[i][j][1]*constant)*(cxNS[9]*uxGhost[i][j][1]*constant+cyNS[9]*uyGhost[i][j][1]*constant+czNS[9]*uzGhost[i][j][1]*constant))-1.5*(uxGhost[i][j][1]*uxGhost[i][j][1]*constant*constant+uyGhost[i][j][1]*uyGhost[i][j][1]*constant*constant+uzGhost[i][j][1]*uzGhost[i][j][1]*constant*constant));
					cIn17First[i][j][1]=-Out14[i][j][1]+2*tNS[17]*ConcentrationFirst[i][j][1]*(1+4.5*((cxNS[17]*uxGhost[i][j][1]*constant+cyNS[17]*uyGhost[i][j][1]*constant+czNS[17]*uzGhost[i][j][1]*constant)*(cxNS[17]*uxGhost[i][j][1]*constant+cyNS[17]*uyGhost[i][j][1]*constant+czNS[17]*uzGhost[i][j][1]*constant))-1.5*(uxGhost[i][j][1]*uxGhost[i][j][1]*constant*constant+uyGhost[i][j][1]*uyGhost[i][j][1]*constant*constant+uzGhost[i][j][1]*uzGhost[i][j][1]*constant*constant));
					cIn18First[i][j][1]=-Out15[i][j][1]+2*tNS[18]*ConcentrationFirst[i][j][1]*(1+4.5*((cxNS[18]*uxGhost[i][j][1]*constant+cyNS[18]*uyGhost[i][j][1]*constant+czNS[18]*uzGhost[i][j][1]*constant)*(cxNS[18]*uxGhost[i][j][1]*constant+cyNS[18]*uyGhost[i][j][1]*constant+czNS[18]*uzGhost[i][j][1]*constant))-1.5*(uxGhost[i][j][1]*uxGhost[i][j][1]*constant*constant+uyGhost[i][j][1]*uyGhost[i][j][1]*constant*constant+uzGhost[i][j][1]*uzGhost[i][j][1]*constant*constant));
				}
				if (domain[i][j][1]!=0||rho[i][j][1]<=gasCritical){
					ConcentrationFirst[i][j][1]=0.0;
				}
				ConcentrationFirst[i][j][0]=ConcentrationFirst[i][j][1];
				cIn3First[i][j][0]=cIn3First[i][j][1];
				cIn8First[i][j][0]=cIn8First[i][j][1];
				cIn9First[i][j][0]=cIn9First[i][j][1];
				cIn17First[i][j][0]=cIn17First[i][j][1];
				cIn18First[i][j][0]=cIn18First[i][j][1];
			}
		}
	}
}

void HeightFirstDirichletConcentrationFirstD3Q7(int Length,int Width,double InletConcentration,double gasCritical,double constant){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				if (domain[i][j][1]==0&&rho[i][j][1]>gasCritical){
					ConcentrationFirst[i][j][1]=InletConcentration;
//					cIn3First[i][j][1]=-Out6[i][j][1]+2*tNSD3Q7[3]*ConcentrationFirst[i][j][1];
					cIn0First[i][j][1] = ConcentrationFirst[i][j][1] * tNSD3Q7[0];
					cIn1First[i][j][1] = ConcentrationFirst[i][j][1] * tNSD3Q7[1];
					cIn2First[i][j][1] = ConcentrationFirst[i][j][1] * tNSD3Q7[2];
					cIn3First[i][j][1] = ConcentrationFirst[i][j][1] * tNSD3Q7[3];
					cIn4First[i][j][1] = ConcentrationFirst[i][j][1] * tNSD3Q7[4];
					cIn5First[i][j][1] = ConcentrationFirst[i][j][1] * tNSD3Q7[5];
					cIn6First[i][j][1] = ConcentrationFirst[i][j][1] * tNSD3Q7[6];
				}
				if (domain[i][j][1]!=0||rho[i][j][1]<=gasCritical){
					ConcentrationFirst[i][j][1]=0.0;
				}
				ConcentrationFirst[i][j][0]=ConcentrationFirst[i][j][1];
				cIn0First[i][j][0] = cIn0First[i][j][1];
				cIn1First[i][j][0] = cIn1First[i][j][1];
				cIn2First[i][j][0] = cIn2First[i][j][1];
				cIn3First[i][j][0] = cIn3First[i][j][1];
				cIn4First[i][j][0] = cIn4First[i][j][1];
				cIn5First[i][j][0] = cIn5First[i][j][1];
				cIn6First[i][j][0] = cIn6First[i][j][1];
			}
		}
	}
}

void HeightLastDirichletConcentrationFirstD3Q7(int Length, int Width, double InletConcentration, double gasCritical, int constant){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (j = 0; j<Width; j++){
				if (domain[i][j][constant] == 0 && rho[i][j][constant]>gasCritical){
					ConcentrationFirst[i][j][constant] = InletConcentration;
					cIn0First[i][j][constant] = ConcentrationFirst[i][j][constant] * tNSD3Q7[0];
					cIn1First[i][j][constant] = ConcentrationFirst[i][j][constant] * tNSD3Q7[1];
					cIn2First[i][j][constant] = ConcentrationFirst[i][j][constant] * tNSD3Q7[2];
					cIn3First[i][j][constant] = ConcentrationFirst[i][j][constant] * tNSD3Q7[3];
					cIn4First[i][j][constant] = ConcentrationFirst[i][j][constant] * tNSD3Q7[4];
					cIn5First[i][j][constant] = ConcentrationFirst[i][j][constant] * tNSD3Q7[5];
					cIn6First[i][j][constant] = ConcentrationFirst[i][j][constant] * tNSD3Q7[6];
				}
				if (domain[i][j][constant] != 0 || rho[i][j][constant] <= gasCritical){
					ConcentrationFirst[i][j][constant] = 0.0;
				}
				ConcentrationFirst[i][j][constant + 1] = ConcentrationFirst[i][j][constant];
				cIn0First[i][j][constant + 1] = cIn0First[i][j][constant];
				cIn1First[i][j][constant + 1] = cIn1First[i][j][constant];
				cIn2First[i][j][constant + 1] = cIn2First[i][j][constant];
				cIn3First[i][j][constant + 1] = cIn3First[i][j][constant];
				cIn4First[i][j][constant + 1] = cIn4First[i][j][constant];
				cIn5First[i][j][constant + 1] = cIn5First[i][j][constant];
				cIn6First[i][j][constant + 1] = cIn6First[i][j][constant];
			}
		}
	}
}

void HeightFirstDirichletConcentrationSecond(int Length,int Width,double InletConcentration,double gasCritical,double constant){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				if (domain[i][j][1]==0&&rho[i][j][1]>gasCritical){
					ConcentrationSecond[i][j][1]=InletConcentration;
					cIn3Second[i][j][1]=-Out6[i][j][1]+2*tNS[3]*ConcentrationSecond[i][j][1]*(1+4.5*((cxNS[3]*uxGhost[i][j][1]*constant+cyNS[3]*uyGhost[i][j][1]*constant+czNS[3]*uzGhost[i][j][1]*constant)*(cxNS[3]*uxGhost[i][j][1]*constant+cyNS[3]*uyGhost[i][j][1]*constant+czNS[3]*uzGhost[i][j][1]*constant))-1.5*(uxGhost[i][j][1]*uxGhost[i][j][1]*constant*constant+uyGhost[i][j][1]*uyGhost[i][j][1]*constant*constant+uzGhost[i][j][1]*uzGhost[i][j][1]*constant*constant));
					cIn8Second[i][j][1]=-Out11[i][j][1]+2*tNS[8]*ConcentrationSecond[i][j][1]*(1+4.5*((cxNS[8]*uxGhost[i][j][1]*constant+cyNS[8]*uyGhost[i][j][1]*constant+czNS[8]*uzGhost[i][j][1]*constant)*(cxNS[8]*uxGhost[i][j][1]*constant+cyNS[8]*uyGhost[i][j][1]*constant+czNS[8]*uzGhost[i][j][1]*constant))-1.5*(uxGhost[i][j][1]*uxGhost[i][j][1]*constant*constant+uyGhost[i][j][1]*uyGhost[i][j][1]*constant*constant+uzGhost[i][j][1]*uzGhost[i][j][1]*constant*constant));
					cIn9Second[i][j][1]=-Out12[i][j][1]+2*tNS[9]*ConcentrationSecond[i][j][1]*(1+4.5*((cxNS[9]*uxGhost[i][j][1]*constant+cyNS[9]*uyGhost[i][j][1]*constant+czNS[9]*uzGhost[i][j][1]*constant)*(cxNS[9]*uxGhost[i][j][1]*constant+cyNS[9]*uyGhost[i][j][1]*constant+czNS[9]*uzGhost[i][j][1]*constant))-1.5*(uxGhost[i][j][1]*uxGhost[i][j][1]*constant*constant+uyGhost[i][j][1]*uyGhost[i][j][1]*constant*constant+uzGhost[i][j][1]*uzGhost[i][j][1]*constant*constant));
					cIn17Second[i][j][1]=-Out14[i][j][1]+2*tNS[17]*ConcentrationSecond[i][j][1]*(1+4.5*((cxNS[17]*uxGhost[i][j][1]*constant+cyNS[17]*uyGhost[i][j][1]*constant+czNS[17]*uzGhost[i][j][1]*constant)*(cxNS[17]*uxGhost[i][j][1]*constant+cyNS[17]*uyGhost[i][j][1]*constant+czNS[17]*uzGhost[i][j][1]*constant))-1.5*(uxGhost[i][j][1]*uxGhost[i][j][1]*constant*constant+uyGhost[i][j][1]*uyGhost[i][j][1]*constant*constant+uzGhost[i][j][1]*uzGhost[i][j][1]*constant*constant));
					cIn18Second[i][j][1]=-Out15[i][j][1]+2*tNS[18]*ConcentrationSecond[i][j][1]*(1+4.5*((cxNS[18]*uxGhost[i][j][1]*constant+cyNS[18]*uyGhost[i][j][1]*constant+czNS[18]*uzGhost[i][j][1]*constant)*(cxNS[18]*uxGhost[i][j][1]*constant+cyNS[18]*uyGhost[i][j][1]*constant+czNS[18]*uzGhost[i][j][1]*constant))-1.5*(uxGhost[i][j][1]*uxGhost[i][j][1]*constant*constant+uyGhost[i][j][1]*uyGhost[i][j][1]*constant*constant+uzGhost[i][j][1]*uzGhost[i][j][1]*constant*constant));
				}
				if (domain[i][j][1]!=0||rho[i][j][1]<=gasCritical){
					ConcentrationSecond[i][j][1]=0.0;
				}
				ConcentrationSecond[i][j][0]=ConcentrationSecond[i][j][1];
				cIn3Second[i][j][0]=cIn3Second[i][j][1];
				cIn8Second[i][j][0]=cIn8Second[i][j][1];
				cIn9Second[i][j][0]=cIn9Second[i][j][1];
				cIn17Second[i][j][0]=cIn17Second[i][j][1];
				cIn18Second[i][j][0]=cIn18Second[i][j][1];
			}
		}
	}
}

void HeightFirstDirichletConcentrationSecondD3Q7(int Length,int Width,double InletConcentration,double gasCritical,double constant){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				if (domain[i][j][1]==0&&rho[i][j][1]>gasCritical){
					ConcentrationSecond[i][j][1]=InletConcentration;
//					cIn3Second[i][j][1]=-Out6[i][j][1]+2*tNSD3Q7[3]*ConcentrationSecond[i][j][1];
					cIn0Second[i][j][1] = ConcentrationSecond[i][j][1] * tNSD3Q7[0];
					cIn1Second[i][j][1] = ConcentrationSecond[i][j][1] * tNSD3Q7[1];
					cIn2Second[i][j][1] = ConcentrationSecond[i][j][1] * tNSD3Q7[2];
					cIn3Second[i][j][1] = ConcentrationSecond[i][j][1] * tNSD3Q7[3];
					cIn4Second[i][j][1] = ConcentrationSecond[i][j][1] * tNSD3Q7[4];
					cIn5Second[i][j][1] = ConcentrationSecond[i][j][1] * tNSD3Q7[5];
					cIn6Second[i][j][1] = ConcentrationSecond[i][j][1] * tNSD3Q7[6];
				}
				if (domain[i][j][1]!=0||rho[i][j][1]<=gasCritical){
					ConcentrationSecond[i][j][1]=0.0;
				}
				ConcentrationSecond[i][j][0]=ConcentrationSecond[i][j][1];
				cIn0Second[i][j][0] = cIn0Second[i][j][1];
				cIn1Second[i][j][0] = cIn1Second[i][j][1];
				cIn2Second[i][j][0] = cIn2Second[i][j][1];
				cIn3Second[i][j][0] = cIn3Second[i][j][1];
				cIn4Second[i][j][0] = cIn4Second[i][j][1];
				cIn5Second[i][j][0] = cIn5Second[i][j][1];
				cIn6Second[i][j][0] = cIn6Second[i][j][1];
			}
		}
	}
}

void HeightLastDirichletConcentrationSecondD3Q7(int Length, int Width, double InletConcentration, double gasCritical, int constant){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (j = 0; j<Width; j++){
				if (domain[i][j][constant] == 0 && rho[i][j][constant]>gasCritical){
					ConcentrationSecond[i][j][constant] = InletConcentration;
					cIn0Second[i][j][constant] = ConcentrationSecond[i][j][constant] * tNSD3Q7[0];
					cIn1Second[i][j][constant] = ConcentrationSecond[i][j][constant] * tNSD3Q7[1];
					cIn2Second[i][j][constant] = ConcentrationSecond[i][j][constant] * tNSD3Q7[2];
					cIn3Second[i][j][constant] = ConcentrationSecond[i][j][constant] * tNSD3Q7[3];
					cIn4Second[i][j][constant] = ConcentrationSecond[i][j][constant] * tNSD3Q7[4];
					cIn5Second[i][j][constant] = ConcentrationSecond[i][j][constant] * tNSD3Q7[5];
					cIn6Second[i][j][constant] = ConcentrationSecond[i][j][constant] * tNSD3Q7[6];
				}
				if (domain[i][j][constant] != 0 || rho[i][j][constant] <= gasCritical){
					ConcentrationSecond[i][j][constant] = 0.0;
				}
				ConcentrationSecond[i][j][constant + 1] = ConcentrationSecond[i][j][constant];
				cIn0Second[i][j][constant + 1] = cIn0Second[i][j][constant];
				cIn1Second[i][j][constant + 1] = cIn1Second[i][j][constant];
				cIn2Second[i][j][constant + 1] = cIn2Second[i][j][constant];
				cIn3Second[i][j][constant + 1] = cIn3Second[i][j][constant];
				cIn4Second[i][j][constant + 1] = cIn4Second[i][j][constant];
				cIn5Second[i][j][constant + 1] = cIn5Second[i][j][constant];
				cIn6Second[i][j][constant + 1] = cIn6Second[i][j][constant];
			}
		}
	}
}

void HeightLastDirichletConcentration(int q,int Length,int Width,double OutletConcentration,double gasCritical){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				if (domain[i][j][q-2]==0&&rho[i][j][q-2]<gasCritical){
					Concentration[i][j][q-2]=OutletConcentration;
					cIn6[i][j][q-2]=-Out3[i][j][q-2]+2*tNS[6]*Concentration[i][j][q-2]*(1+4.5*((cxNS[6]*uxGhost[i][j][q-2]+cyNS[6]*uyGhost[i][j][q-2]+czNS[6]*uzGhost[i][j][q-2])*(cxNS[6]*uxGhost[i][j][q-2]+cyNS[6]*uyGhost[i][j][q-2]+czNS[6]*uzGhost[i][j][q-2]))-1.5*(uxGhost[i][j][q-2]*ux[i][j][q-2]+uyGhost[i][j][q-2]*uyGhost[i][j][q-2]+uzGhost[i][j][q-2]*uzGhost[i][j][q-2]));
					cIn11[i][j][q-2]=-Out8[i][j][q-2]+2*tNS[11]*Concentration[i][j][q-2]*(1+4.5*((cxNS[11]*uxGhost[i][j][q-2]+cyNS[11]*uyGhost[i][j][q-2]+czNS[11]*uzGhost[i][j][q-2])*(cxNS[11]*uxGhost[i][j][q-2]+cyNS[11]*uyGhost[i][j][q-2]+czNS[11]*uzGhost[i][j][q-2]))-1.5*(uxGhost[i][j][q-2]*ux[i][j][q-2]+uyGhost[i][j][q-2]*uyGhost[i][j][q-2]+uzGhost[i][j][q-2]*uzGhost[i][j][q-2]));
					cIn12[i][j][q-2]=-Out9[i][j][q-2]+2*tNS[12]*Concentration[i][j][q-2]*(1+4.5*((cxNS[12]*uxGhost[i][j][q-2]+cyNS[12]*uyGhost[i][j][q-2]+czNS[12]*uzGhost[i][j][q-2])*(cxNS[12]*uxGhost[i][j][q-2]+cyNS[12]*uyGhost[i][j][q-2]+czNS[12]*uzGhost[i][j][q-2]))-1.5*(uxGhost[i][j][q-2]*ux[i][j][q-2]+uyGhost[i][j][q-2]*uyGhost[i][j][q-2]+uzGhost[i][j][q-2]*uzGhost[i][j][q-2]));
					cIn14[i][j][q-2]=-Out17[i][j][q-2]+2*tNS[14]*Concentration[i][j][q-2]*(1+4.5*((cxNS[14]*uxGhost[i][j][q-2]+cyNS[14]*uyGhost[i][j][q-2]+czNS[14]*uzGhost[i][j][q-2])*(cxNS[14]*uxGhost[i][j][q-2]+cyNS[14]*uyGhost[i][j][q-2]+czNS[14]*uzGhost[i][j][q-2]))-1.5*(uxGhost[i][j][q-2]*ux[i][j][q-2]+uyGhost[i][j][q-2]*uyGhost[i][j][q-2]+uzGhost[i][j][q-2]*uzGhost[i][j][q-2]));
					cIn15[i][j][q-2]=-Out18[i][j][q-2]+2*tNS[15]*Concentration[i][j][q-2]*(1+4.5*((cxNS[15]*uxGhost[i][j][q-2]+cyNS[15]*uyGhost[i][j][q-2]+czNS[15]*uzGhost[i][j][q-2])*(cxNS[15]*uxGhost[i][j][q-2]+cyNS[15]*uyGhost[i][j][q-2]+czNS[15]*uzGhost[i][j][q-2]))-1.5*(uxGhost[i][j][q-2]*ux[i][j][q-2]+uyGhost[i][j][q-2]*uyGhost[i][j][q-2]+uzGhost[i][j][q-2]*uzGhost[i][j][q-2]));
				}
				if (domain[i][j][q-2]!=0||rho[i][j][q-2]>=gasCritical){
					Concentration[i][j][q-2]=0.0;
				}
				Concentration[i][j][q-1]=Concentration[i][j][q-2];
				cIn6[i][j][q-1]=cIn6[i][j][q-2];
				cIn11[i][j][q-1]=cIn11[i][j][q-2];
				cIn12[i][j][q-1]=cIn12[i][j][q-2];
				cIn14[i][j][q-1]=cIn14[i][j][q-2];
				cIn15[i][j][q-1]=cIn15[i][j][q-2];
			}
		}
	}
}

void HeightLastNeumannConcentration(int Length,int Width,int Height){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
					cIn0[i][j][Height-1]=cIn0[i][j][Height-2];
					cIn1[i][j][Height-1]=cIn1[i][j][Height-2];
					cIn2[i][j][Height-1]=cIn2[i][j][Height-2];
					cIn3[i][j][Height-1]=cIn3[i][j][Height-2];
					cIn4[i][j][Height-1]=cIn4[i][j][Height-2];
					cIn5[i][j][Height-1]=cIn5[i][j][Height-2];
					cIn6[i][j][Height-1]=cIn6[i][j][Height-2];
					cIn7[i][j][Height-1]=cIn7[i][j][Height-2];
					cIn8[i][j][Height-1]=cIn8[i][j][Height-2];
					cIn9[i][j][Height-1]=cIn9[i][j][Height-2];
					cIn10[i][j][Height-1]=cIn10[i][j][Height-2];
					cIn11[i][j][Height-1]=cIn11[i][j][Height-2];
					cIn12[i][j][Height-1]=cIn12[i][j][Height-2];
					cIn13[i][j][Height-1]=cIn13[i][j][Height-2];
					cIn14[i][j][Height-1]=cIn14[i][j][Height-2];
					cIn15[i][j][Height-1]=cIn15[i][j][Height-2];
					cIn16[i][j][Height-1]=cIn16[i][j][Height-2];
					cIn17[i][j][Height-1]=cIn17[i][j][Height-2];
					cIn18[i][j][Height-1]=cIn18[i][j][Height-2];
					Concentration[i][j][Height-1]=Concentration[i][j][Height-2];	
			}
		}
	}
}

void HeightLastNeumannConcentrationD3Q7(int Length,int Width,int Height){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
					cIn0[i][j][Height-1]=cIn0[i][j][Height-2];
					cIn1[i][j][Height-1]=cIn1[i][j][Height-2];
					cIn2[i][j][Height-1]=cIn2[i][j][Height-2];
					cIn3[i][j][Height-1]=cIn3[i][j][Height-2];
					cIn4[i][j][Height-1]=cIn4[i][j][Height-2];
					cIn5[i][j][Height-1]=cIn5[i][j][Height-2];
					cIn6[i][j][Height-1]=cIn6[i][j][Height-2];
					Concentration[i][j][Height-1]=Concentration[i][j][Height-2];	
			}
		}
	}
}

void WidthLastNeumannConcentrationD3Q7(int Length, int Width, int initialHeight,int Height){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (k = initialHeight; k<Height; k++){
				cIn0[i][Width - 1][k] = cIn0[i][Width - 2][k];
				cIn1[i][Width - 1][k] = cIn1[i][Width - 2][k];
				cIn2[i][Width - 1][k] = cIn2[i][Width - 2][k];
				cIn3[i][Width - 1][k] = cIn3[i][Width - 2][k];
				cIn4[i][Width - 1][k] = cIn4[i][Width - 2][k];
				cIn5[i][Width - 1][k] = cIn5[i][Width - 2][k];
				cIn6[i][Width - 1][k] = cIn6[i][Width - 2][k];
				Concentration[i][Width - 1][k] = Concentration[i][Width - 2][k];
			}
		}
	}
}

void WidthLastNeumannConcentrationFirstD3Q7(int Length, int Width, int initialHeight, int Height){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (k = initialHeight; k<Height; k++){
				cIn0First[i][Width - 1][k] = cIn0First[i][Width - 2][k];
				cIn1First[i][Width - 1][k] = cIn1First[i][Width - 2][k];
				cIn2First[i][Width - 1][k] = cIn2First[i][Width - 2][k];
				cIn3First[i][Width - 1][k] = cIn3First[i][Width - 2][k];
				cIn4First[i][Width - 1][k] = cIn4First[i][Width - 2][k];
				cIn5First[i][Width - 1][k] = cIn5First[i][Width - 2][k];
				cIn6First[i][Width - 1][k] = cIn6First[i][Width - 2][k];
				ConcentrationFirst[i][Width - 1][k] = ConcentrationFirst[i][Width - 2][k];
			}
		}
	}
}

void HeightFirstNeumannConcentrationD3Q7(int Length, int Width, int Height){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (j = 0; j<Width; j++){
				cIn0[i][j][0] = cIn0[i][j][1];
				cIn1[i][j][0] = cIn1[i][j][1];
				cIn2[i][j][0] = cIn2[i][j][1];
				cIn3[i][j][0] = cIn3[i][j][1];
				cIn4[i][j][0] = cIn4[i][j][1];
				cIn5[i][j][0] = cIn5[i][j][1];
				cIn6[i][j][0] = cIn6[i][j][1];
				Concentration[i][j][0] = Concentration[i][j][1];
			}
		}
	}
}

void HeightLastNeumannConcentrationFirst(int Length,int Width,int Height){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
					cIn0First[i][j][Height-1]=cIn0First[i][j][Height-2];
					cIn1First[i][j][Height-1]=cIn1First[i][j][Height-2];
					cIn2First[i][j][Height-1]=cIn2First[i][j][Height-2];
					cIn3First[i][j][Height-1]=cIn3First[i][j][Height-2];
					cIn4First[i][j][Height-1]=cIn4First[i][j][Height-2];
					cIn5First[i][j][Height-1]=cIn5First[i][j][Height-2];
					cIn6First[i][j][Height-1]=cIn6First[i][j][Height-2];
					cIn7First[i][j][Height-1]=cIn7First[i][j][Height-2];
					cIn8First[i][j][Height-1]=cIn8First[i][j][Height-2];
					cIn9First[i][j][Height-1]=cIn9First[i][j][Height-2];
					cIn10First[i][j][Height-1]=cIn10First[i][j][Height-2];
					cIn11First[i][j][Height-1]=cIn11First[i][j][Height-2];
					cIn12First[i][j][Height-1]=cIn12First[i][j][Height-2];
					cIn13First[i][j][Height-1]=cIn13First[i][j][Height-2];
					cIn14First[i][j][Height-1]=cIn14First[i][j][Height-2];
					cIn15First[i][j][Height-1]=cIn15First[i][j][Height-2];
					cIn16First[i][j][Height-1]=cIn16First[i][j][Height-2];
					cIn17First[i][j][Height-1]=cIn17First[i][j][Height-2];
					cIn18First[i][j][Height-1]=cIn18First[i][j][Height-2];
					ConcentrationFirst[i][j][Height-1]=ConcentrationFirst[i][j][Height-2];	
			}
		}
	}
}

void HeightLastNeumannConcentrationFirstD3Q7(int Length,int Width,int Height){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
					cIn0First[i][j][Height-1]=cIn0First[i][j][Height-2];
					cIn1First[i][j][Height-1]=cIn1First[i][j][Height-2];
					cIn2First[i][j][Height-1]=cIn2First[i][j][Height-2];
					cIn3First[i][j][Height-1]=cIn3First[i][j][Height-2];
					cIn4First[i][j][Height-1]=cIn4First[i][j][Height-2];
					cIn5First[i][j][Height-1]=cIn5First[i][j][Height-2];
					cIn6First[i][j][Height-1]=cIn6First[i][j][Height-2];
					ConcentrationFirst[i][j][Height-1]=ConcentrationFirst[i][j][Height-2];	
			}
		}
	}
}

void HeightFirstNeumannConcentrationFirstD3Q7(int Length, int Width, int Height){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (j = 0; j<Width; j++){
				cIn0First[i][j][0] = cIn0First[i][j][1];
				cIn1First[i][j][0] = cIn1First[i][j][1];
				cIn2First[i][j][0] = cIn2First[i][j][1];
				cIn3First[i][j][0] = cIn3First[i][j][1];
				cIn4First[i][j][0] = cIn4First[i][j][1];
				cIn5First[i][j][0] = cIn5First[i][j][1];
				cIn6First[i][j][0] = cIn6First[i][j][1];
				ConcentrationFirst[i][j][0] = ConcentrationFirst[i][j][1];
			}
		}
	}
}

void HeightLastNeumannConcentrationSecond(int Length,int Width,int Height){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
					cIn0Second[i][j][Height-1]=cIn0Second[i][j][Height-2];
					cIn1Second[i][j][Height-1]=cIn1Second[i][j][Height-2];
					cIn2Second[i][j][Height-1]=cIn2Second[i][j][Height-2];
					cIn3Second[i][j][Height-1]=cIn3Second[i][j][Height-2];
					cIn4Second[i][j][Height-1]=cIn4Second[i][j][Height-2];
					cIn5Second[i][j][Height-1]=cIn5Second[i][j][Height-2];
					cIn6Second[i][j][Height-1]=cIn6Second[i][j][Height-2];
					cIn7Second[i][j][Height-1]=cIn7Second[i][j][Height-2];
					cIn8Second[i][j][Height-1]=cIn8Second[i][j][Height-2];
					cIn9Second[i][j][Height-1]=cIn9Second[i][j][Height-2];
					cIn10Second[i][j][Height-1]=cIn10Second[i][j][Height-2];
					cIn11Second[i][j][Height-1]=cIn11Second[i][j][Height-2];
					cIn12Second[i][j][Height-1]=cIn12Second[i][j][Height-2];
					cIn13Second[i][j][Height-1]=cIn13Second[i][j][Height-2];
					cIn14Second[i][j][Height-1]=cIn14Second[i][j][Height-2];
					cIn15Second[i][j][Height-1]=cIn15Second[i][j][Height-2];
					cIn16Second[i][j][Height-1]=cIn16Second[i][j][Height-2];
					cIn17Second[i][j][Height-1]=cIn17Second[i][j][Height-2];
					cIn18Second[i][j][Height-1]=cIn18Second[i][j][Height-2];
					ConcentrationSecond[i][j][Height-1]=ConcentrationSecond[i][j][Height-2];	
			}
		}
	}
}

void HeightLastNeumannConcentrationSecondD3Q7(int Length,int Width,int Height){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
					cIn0Second[i][j][Height-1]=cIn0Second[i][j][Height-2];
					cIn1Second[i][j][Height-1]=cIn1Second[i][j][Height-2];
					cIn2Second[i][j][Height-1]=cIn2Second[i][j][Height-2];
					cIn3Second[i][j][Height-1]=cIn3Second[i][j][Height-2];
					cIn4Second[i][j][Height-1]=cIn4Second[i][j][Height-2];
					cIn5Second[i][j][Height-1]=cIn5Second[i][j][Height-2];
					cIn6Second[i][j][Height-1]=cIn6Second[i][j][Height-2];
					ConcentrationSecond[i][j][Height-1]=ConcentrationSecond[i][j][Height-2];	
			}
		}
	}
}

void widthFirstNeumannElectrodePotential(int Length,int Width,int Height){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (k=0;k<Height;k++){
					electrodeIn0[i][0][k]=electrodeIn0[i][1][k];
					electrodeIn1[i][0][k]=electrodeIn1[i][1][k];
					electrodeIn2[i][0][k]=electrodeIn2[i][1][k];
					electrodeIn3[i][0][k]=electrodeIn3[i][1][k];
					electrodeIn4[i][0][k]=electrodeIn4[i][1][k];
					electrodeIn5[i][0][k]=electrodeIn5[i][1][k];
					electrodeIn6[i][0][k]=electrodeIn6[i][1][k];
					electrodeIn7[i][0][k]=electrodeIn7[i][1][k];
					electrodeIn8[i][0][k]=electrodeIn8[i][1][k];
					electrodeIn9[i][0][k]=electrodeIn9[i][1][k];
					electrodeIn10[i][0][k]=electrodeIn10[i][1][k];
					electrodeIn11[i][0][k]=electrodeIn11[i][1][k];
					electrodeIn12[i][0][k]=electrodeIn12[i][1][k];
					electrodeIn13[i][0][k]=electrodeIn13[i][1][k];
					electrodeIn14[i][0][k]=electrodeIn14[i][1][k];
					electrodeIn15[i][0][k]=electrodeIn15[i][1][k];
					electrodeIn16[i][0][k]=electrodeIn16[i][1][k];
					electrodeIn17[i][0][k]=electrodeIn17[i][1][k];
					electrodeIn18[i][0][k]=electrodeIn18[i][1][k];
					electrodePotential[i][0][k]=electrodePotential[i][1][k];	
			}
		}
	}
}

void widthFirstNeumannElectrodePotentialD3Q7(int Length,int Width,int Height){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (k=0;k<Height;k++){
					electrodeIn0[i][0][k]=electrodeIn0[i][1][k];
					electrodeIn1[i][0][k]=electrodeIn1[i][1][k];
					electrodeIn2[i][0][k]=electrodeIn2[i][1][k];
					electrodeIn3[i][0][k]=electrodeIn3[i][1][k];
					electrodeIn4[i][0][k]=electrodeIn4[i][1][k];
					electrodeIn5[i][0][k]=electrodeIn5[i][1][k];
					electrodeIn6[i][0][k]=electrodeIn6[i][1][k];
					electrodePotential[i][0][k]=electrodePotential[i][1][k];	
			}
		}
	}
}

void widthFirstDirichletElectrodePotentialD3Q7(int Length,int Width,int Height, double InletElectrodePotential){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (k=0;k<Height;k++){
				if (domain[i][1][k]==1){
					electrodePotential[i][1][k]=InletElectrodePotential;
					electrodeIn2[i][1][k]=-Out5[i][1][k]+2*tNSD3Q7[2]*electrodePotential[i][1][k];
				}
				if (domain[i][1][k]==0){
					electrodePotential[i][1][k]=0.0;
				}
				electrodePotential[i][0][k]=electrodePotential[i][1][k];
				electrodeIn2[i][0][k]=electrodeIn2[i][1][k];
			}
		}
	}
}

void lengthFirstNeumannElectrodePotential(int Length,int Width,int Height){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
					electrodeIn0[0][j][k]=electrodeIn0[1][j][k];
					electrodeIn1[0][j][k]=electrodeIn1[1][j][k];
					electrodeIn2[0][j][k]=electrodeIn2[1][j][k];
					electrodeIn3[0][j][k]=electrodeIn3[1][j][k];
					electrodeIn4[0][j][k]=electrodeIn4[1][j][k];
					electrodeIn5[0][j][k]=electrodeIn5[1][j][k];
					electrodeIn6[0][j][k]=electrodeIn6[1][j][k];
					electrodeIn7[0][j][k]=electrodeIn7[1][j][k];
					electrodeIn8[0][j][k]=electrodeIn8[1][j][k];
					electrodeIn9[0][j][k]=electrodeIn9[1][j][k];
					electrodeIn10[0][j][k]=electrodeIn10[1][j][k];
					electrodeIn11[0][j][k]=electrodeIn11[1][j][k];
					electrodeIn12[0][j][k]=electrodeIn12[1][j][k];
					electrodeIn13[0][j][k]=electrodeIn13[1][j][k];
					electrodeIn14[0][j][k]=electrodeIn14[1][j][k];
					electrodeIn15[0][j][k]=electrodeIn15[1][j][k];
					electrodeIn16[0][j][k]=electrodeIn16[1][j][k];
					electrodeIn17[0][j][k]=electrodeIn17[1][j][k];
					electrodeIn18[0][j][k]=electrodeIn18[1][j][k];
					electrodePotential[0][j][k]=electrodePotential[1][j][k];	
			}
		}
	}
}

void lengthFirstNeumannElectrodePotentialD3Q7(int Length,int Width,int Height){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
					electrodeIn0[0][j][k]=electrodeIn0[1][j][k];
					electrodeIn1[0][j][k]=electrodeIn1[1][j][k];
					electrodeIn2[0][j][k]=electrodeIn2[1][j][k];
					electrodeIn3[0][j][k]=electrodeIn3[1][j][k];
					electrodeIn4[0][j][k]=electrodeIn4[1][j][k];
					electrodeIn5[0][j][k]=electrodeIn5[1][j][k];
					electrodeIn6[0][j][k]=electrodeIn6[1][j][k];
					electrodePotential[0][j][k]=electrodePotential[1][j][k];	
			}
		}
	}
}

void lengthLastNeumannElectrodePotential(int Length,int Width,int Height){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
					electrodeIn0[Length-1][j][k]=electrodeIn0[Length-2][j][k];
					electrodeIn1[Length-1][j][k]=electrodeIn1[Length-2][j][k];
					electrodeIn2[Length-1][j][k]=electrodeIn2[Length-2][j][k];
					electrodeIn3[Length-1][j][k]=electrodeIn3[Length-2][j][k];
					electrodeIn4[Length-1][j][k]=electrodeIn4[Length-2][j][k];
					electrodeIn5[Length-1][j][k]=electrodeIn5[Length-2][j][k];
					electrodeIn6[Length-1][j][k]=electrodeIn6[Length-2][j][k];
					electrodeIn7[Length-1][j][k]=electrodeIn7[Length-2][j][k];
					electrodeIn8[Length-1][j][k]=electrodeIn8[Length-2][j][k];
					electrodeIn9[Length-1][j][k]=electrodeIn9[Length-2][j][k];
					electrodeIn10[Length-1][j][k]=electrodeIn10[Length-2][j][k];
					electrodeIn11[Length-1][j][k]=electrodeIn11[Length-2][j][k];
					electrodeIn12[Length-1][j][k]=electrodeIn12[Length-2][j][k];
					electrodeIn13[Length-1][j][k]=electrodeIn13[Length-2][j][k];
					electrodeIn14[Length-1][j][k]=electrodeIn14[Length-2][j][k];
					electrodeIn15[Length-1][j][k]=electrodeIn15[Length-2][j][k];
					electrodeIn16[Length-1][j][k]=electrodeIn16[Length-2][j][k];
					electrodeIn17[Length-1][j][k]=electrodeIn17[Length-2][j][k];
					electrodeIn18[Length-1][j][k]=electrodeIn18[Length-2][j][k];
					electrodePotential[Length-1][j][k]=electrodePotential[Length-2][j][k];	
			}
		}
	}
}

void lengthLastNeumannElectrodePotentialD3Q7(int Length,int Width,int Height){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
					electrodeIn0[Length-1][j][k]=electrodeIn0[Length-2][j][k];
					electrodeIn1[Length-1][j][k]=electrodeIn1[Length-2][j][k];
					electrodeIn2[Length-1][j][k]=electrodeIn2[Length-2][j][k];
					electrodeIn3[Length-1][j][k]=electrodeIn3[Length-2][j][k];
					electrodeIn4[Length-1][j][k]=electrodeIn4[Length-2][j][k];
					electrodeIn5[Length-1][j][k]=electrodeIn5[Length-2][j][k];
					electrodeIn6[Length-1][j][k]=electrodeIn6[Length-2][j][k];
					electrodePotential[Length-1][j][k]=electrodePotential[Length-2][j][k];	
			}
		}
	}
}

void HeightFirstNeumannElectrodePotential(int Length,int Width,int Height){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
					electrodeIn0[i][j][0]=electrodeIn0[i][j][1];
					electrodeIn1[i][j][0]=electrodeIn1[i][j][1];
					electrodeIn2[i][j][0]=electrodeIn2[i][j][1];
					electrodeIn3[i][j][0]=electrodeIn3[i][j][1];
					electrodeIn4[i][j][0]=electrodeIn4[i][j][1];
					electrodeIn5[i][j][0]=electrodeIn5[i][j][1];
					electrodeIn6[i][j][0]=electrodeIn6[i][j][1];
					electrodeIn7[i][j][0]=electrodeIn7[i][j][1];
					electrodeIn8[i][j][0]=electrodeIn8[i][j][1];
					electrodeIn9[i][j][0]=electrodeIn9[i][j][1];
					electrodeIn10[i][j][0]=electrodeIn10[i][j][1];
					electrodeIn11[i][j][0]=electrodeIn11[i][j][1];
					electrodeIn12[i][j][0]=electrodeIn12[i][j][1];
					electrodeIn13[i][j][0]=electrodeIn13[i][j][1];
					electrodeIn14[i][j][0]=electrodeIn14[i][j][1];
					electrodeIn15[i][j][0]=electrodeIn15[i][j][1];
					electrodeIn16[i][j][0]=electrodeIn16[i][j][1];
					electrodeIn17[i][j][0]=electrodeIn17[i][j][1];
					electrodeIn18[i][j][0]=electrodeIn18[i][j][1];
					electrodePotential[i][j][0]=electrodePotential[i][j][1];	
			}
		}
	}
}

void HeightFirstNeumannElectrodePotentialD3Q7(int Length,int Width,int Height){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				electrodeIn0[i][j][Height] = electrodeIn0[i][j][Height+1];
				electrodeIn1[i][j][Height] = electrodeIn1[i][j][Height+1];
				electrodeIn2[i][j][Height] = electrodeIn2[i][j][Height+1];
				electrodeIn3[i][j][Height] = electrodeIn3[i][j][Height+1];
				electrodeIn4[i][j][Height] = electrodeIn4[i][j][Height+1];
				electrodeIn5[i][j][Height] = electrodeIn5[i][j][Height+1];
				electrodeIn6[i][j][Height] = electrodeIn6[i][j][Height+1];
				electrodePotential[i][j][Height] = electrodePotential[i][j][Height+1];
			}
		}
	}
}

void HeightLastNeumannElectrodePotential(int Length,int Width,int Height){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
					electrodeIn0[i][j][Height-1]=electrodeIn0[i][j][Height-2];
					electrodeIn1[i][j][Height-1]=electrodeIn1[i][j][Height-2];
					electrodeIn2[i][j][Height-1]=electrodeIn2[i][j][Height-2];
					electrodeIn3[i][j][Height-1]=electrodeIn3[i][j][Height-2];
					electrodeIn4[i][j][Height-1]=electrodeIn4[i][j][Height-2];
					electrodeIn5[i][j][Height-1]=electrodeIn5[i][j][Height-2];
					electrodeIn6[i][j][Height-1]=electrodeIn6[i][j][Height-2];
					electrodeIn7[i][j][Height-1]=electrodeIn7[i][j][Height-2];
					electrodeIn8[i][j][Height-1]=electrodeIn8[i][j][Height-2];
					electrodeIn9[i][j][Height-1]=electrodeIn9[i][j][Height-2];
					electrodeIn10[i][j][Height-1]=electrodeIn10[i][j][Height-2];
					electrodeIn11[i][j][Height-1]=electrodeIn11[i][j][Height-2];
					electrodeIn12[i][j][Height-1]=electrodeIn12[i][j][Height-2];
					electrodeIn13[i][j][Height-1]=electrodeIn13[i][j][Height-2];
					electrodeIn14[i][j][Height-1]=electrodeIn14[i][j][Height-2];
					electrodeIn15[i][j][Height-1]=electrodeIn15[i][j][Height-2];
					electrodeIn16[i][j][Height-1]=electrodeIn16[i][j][Height-2];
					electrodeIn17[i][j][Height-1]=electrodeIn17[i][j][Height-2];
					electrodeIn18[i][j][Height-1]=electrodeIn18[i][j][Height-2];
					electrodePotential[i][j][Height-1]=electrodePotential[i][j][Height-2];	
			}
		}
	}
}

void HeightLastNeumannElectrodePotentialD3Q7(int Length,int Width,int Height){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
					electrodeIn0[i][j][Height-1]=electrodeIn0[i][j][Height-2];
					electrodeIn1[i][j][Height-1]=electrodeIn1[i][j][Height-2];
					electrodeIn2[i][j][Height-1]=electrodeIn2[i][j][Height-2];
					electrodeIn3[i][j][Height-1]=electrodeIn3[i][j][Height-2];
					electrodeIn4[i][j][Height-1]=electrodeIn4[i][j][Height-2];
					electrodeIn5[i][j][Height-1]=electrodeIn5[i][j][Height-2];
					electrodeIn6[i][j][Height-1]=electrodeIn6[i][j][Height-2];
					electrodePotential[i][j][Height-1]=electrodePotential[i][j][Height-2];	
			}
		}
	}
}

void widthLastDirichletelectrodePotential(int Length,int Width, int Height,double constant){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (k=0;k<Height;k++){
				if (domain[i][Width-1][k]==1){
					electrodePotential[i][Width-1][k]=electrodePotential[i][Width-2][k]+constant;
					electrodeIn5[i][Width-1][k]=-Out2[i][Width-1][k]+2*tNS[5]*electrodePotential[i][Width-1][k];
					electrodeIn10[i][Width-1][k]=-Out8[i][Width-1][k]+2*tNS[10]*electrodePotential[i][Width-1][k];
					electrodeIn12[i][Width-1][k]=-Out9[i][Width-1][k]+2*tNS[12]*electrodePotential[i][Width-1][k];
					electrodeIn13[i][Width-1][k]=-Out16[i][Width-1][k]+2*tNS[13]*electrodePotential[i][Width-1][k];
					electrodeIn18[i][Width-1][k]=-Out15[i][Width-1][k]+2*tNS[18]*electrodePotential[i][Width-1][k];
				}
			}
		}
	}
}

void widthLastDirichletelectrodePotentialD3Q7(int Length,int Width, int Initial, int Height,double constant){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (k=Initial;k<Height;k++){
				if (domain[i][Width-1][k]==2){
//					electrodePotential[i][Width-1][k]=InletElectrodePotential;
					electrodeIn0[i][Width - 1][k] = electrodeIn0[i][Width - 1][k] + tNSD3Q7[0] * constant;
					electrodeIn1[i][Width - 1][k] = electrodeIn1[i][Width - 1][k] + tNSD3Q7[1] * constant;
					electrodeIn2[i][Width - 1][k] = electrodeIn2[i][Width - 1][k] + tNSD3Q7[2] * constant;
					electrodeIn3[i][Width - 1][k] = electrodeIn3[i][Width - 1][k] + tNSD3Q7[3] * constant;
					electrodeIn4[i][Width - 1][k] = electrodeIn4[i][Width - 1][k] + tNSD3Q7[4] * constant;
					electrodeIn5[i][Width - 1][k] = Out2[i][Width - 1][k] + tNSD3Q7[5] * constant;
					electrodeIn6[i][Width - 1][k] = electrodeIn6[i][Width - 1][k] + tNSD3Q7[6] * constant;
					electrodePotential[i][Width - 1][k] = electrodeIn0[i][Width - 1][k] + electrodeIn1[i][Width - 1][k] + electrodeIn2[i][Width - 1][k] + electrodeIn3[i][Width - 1][k] + electrodeIn4[i][Width - 1][k] + electrodeIn5[i][Width - 1][k] + electrodeIn6[i][Width - 1][k];
				}
			}
		}
	}
}

void widthLastNeumannelectrodePotentialD3Q7(int Length, int Width, int Initial, int Height){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (k = Initial; k<Height; k++){
				if (domain[i][Width - 1][k] == 1){
					electrodeIn0[i][Width - 1][k] = electrodeIn0[i][Width - 2][k];
					electrodeIn1[i][Width - 1][k] = electrodeIn1[i][Width - 2][k];
					electrodeIn2[i][Width - 1][k] = electrodeIn2[i][Width - 2][k];
					electrodeIn3[i][Width - 1][k] = electrodeIn3[i][Width - 2][k];
					electrodeIn4[i][Width - 1][k] = electrodeIn4[i][Width - 2][k];
					electrodeIn5[i][Width - 1][k] = electrodeIn5[i][Width - 2][k];
					electrodeIn6[i][Width - 1][k] = electrodeIn6[i][Width - 2][k];
					electrodePotential[i][Width - 1][k] = electrodePotential[i][Width - 2][k];
				}
			}
		}
	}
}

void widthLastNeumannelectrolytePotentialD3Q7(int Length, int Width, int Initial, int Height){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (k = Initial; k<Height; k++){
				if (domain[i][Width - 1][k] == 0){
					electrolyteIn0[i][Width - 1][k] = electrolyteIn0[i][Width - 2][k];
					electrolyteIn1[i][Width - 1][k] = electrolyteIn1[i][Width - 2][k];
					electrolyteIn2[i][Width - 1][k] = electrolyteIn2[i][Width - 2][k];
					electrolyteIn3[i][Width - 1][k] = electrolyteIn3[i][Width - 2][k];
					electrolyteIn4[i][Width - 1][k] = electrolyteIn4[i][Width - 2][k];
					electrolyteIn5[i][Width - 1][k] = electrolyteIn5[i][Width - 2][k];
					electrolyteIn6[i][Width - 1][k] = electrolyteIn6[i][Width - 2][k];
					electrolytePotential[i][Width - 1][k] = electrolytePotential[i][Width - 2][k];
				}
			}
		}
	}
}

void HeightFirstNeumannElectrolytePotential(int Length,int Width,int Height){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
					electrolyteIn0[i][j][0]=electrolyteIn0[i][j][1];
					electrolyteIn1[i][j][0]=electrolyteIn1[i][j][1];
					electrolyteIn2[i][j][0]=electrolyteIn2[i][j][1];
					electrolyteIn3[i][j][0]=electrolyteIn3[i][j][1];
					electrolyteIn4[i][j][0]=electrolyteIn4[i][j][1];
					electrolyteIn5[i][j][0]=electrolyteIn5[i][j][1];
					electrolyteIn6[i][j][0]=electrolyteIn6[i][j][1];
					electrolyteIn7[i][j][0]=electrolyteIn7[i][j][1];
					electrolyteIn8[i][j][0]=electrolyteIn8[i][j][1];
					electrolyteIn9[i][j][0]=electrolyteIn9[i][j][1];
					electrolyteIn10[i][j][0]=electrolyteIn10[i][j][1];
					electrolyteIn11[i][j][0]=electrolyteIn11[i][j][1];
					electrolyteIn12[i][j][0]=electrolyteIn12[i][j][1];
					electrolyteIn13[i][j][0]=electrolyteIn13[i][j][1];
					electrolyteIn14[i][j][0]=electrolyteIn14[i][j][1];
					electrolyteIn15[i][j][0]=electrolyteIn15[i][j][1];
					electrolyteIn16[i][j][0]=electrolyteIn16[i][j][1];
					electrolyteIn17[i][j][0]=electrolyteIn17[i][j][1];
					electrolyteIn18[i][j][0]=electrolyteIn18[i][j][1];
					electrolytePotential[i][j][0]=electrolytePotential[i][j][1];	
			}
		}
	}
}

void HeightFirstNeumannElectrolytePotentialD3Q7(int Length,int Width,int Height){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				electrolyteIn0[i][j][Height] = electrolyteIn0[i][j][Height+1];
				electrolyteIn1[i][j][Height] = electrolyteIn1[i][j][Height+1];
				electrolyteIn2[i][j][Height] = electrolyteIn2[i][j][Height+1];
				electrolyteIn3[i][j][Height] = electrolyteIn3[i][j][Height+1];
				electrolyteIn4[i][j][Height] = electrolyteIn4[i][j][Height+1];
				electrolyteIn5[i][j][Height] = electrolyteIn5[i][j][Height+1];
				electrolyteIn6[i][j][Height] = electrolyteIn6[i][j][Height+1];
				electrolytePotential[i][j][Height] = electrolytePotential[i][j][Height+1];
			}
		}
	}
}

void HeightLastNeumannElectrolytePotential(int Length,int Width,int Height){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
					electrolyteIn0[i][j][Height-1]=electrolyteIn0[i][j][Height-2];
					electrolyteIn1[i][j][Height-1]=electrolyteIn1[i][j][Height-2];
					electrolyteIn2[i][j][Height-1]=electrolyteIn2[i][j][Height-2];
					electrolyteIn3[i][j][Height-1]=electrolyteIn3[i][j][Height-2];
					electrolyteIn4[i][j][Height-1]=electrolyteIn4[i][j][Height-2];
					electrolyteIn5[i][j][Height-1]=electrolyteIn5[i][j][Height-2];
					electrolyteIn6[i][j][Height-1]=electrolyteIn6[i][j][Height-2];
					electrolyteIn7[i][j][Height-1]=electrolyteIn7[i][j][Height-2];
					electrolyteIn8[i][j][Height-1]=electrolyteIn8[i][j][Height-2];
					electrolyteIn9[i][j][Height-1]=electrolyteIn9[i][j][Height-2];
					electrolyteIn10[i][j][Height-1]=electrolyteIn10[i][j][Height-2];
					electrolyteIn11[i][j][Height-1]=electrolyteIn11[i][j][Height-2];
					electrolyteIn12[i][j][Height-1]=electrolyteIn12[i][j][Height-2];
					electrolyteIn13[i][j][Height-1]=electrolyteIn13[i][j][Height-2];
					electrolyteIn14[i][j][Height-1]=electrolyteIn14[i][j][Height-2];
					electrolyteIn15[i][j][Height-1]=electrolyteIn15[i][j][Height-2];
					electrolyteIn16[i][j][Height-1]=electrolyteIn16[i][j][Height-2];
					electrolyteIn17[i][j][Height-1]=electrolyteIn17[i][j][Height-2];
					electrolyteIn18[i][j][Height-1]=electrolyteIn18[i][j][Height-2];
					electrolytePotential[i][j][Height-1]=electrolytePotential[i][j][Height-2];	
			}
		}
	}
}

void HeightLastNeumannElectrolytePotentialD3Q7(int Length,int Width,int Height){
	int i;
	int j;
#pragma omp parallel private(i,j)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
					electrolyteIn0[i][j][Height-1]=electrolyteIn0[i][j][Height-2];
					electrolyteIn1[i][j][Height-1]=electrolyteIn1[i][j][Height-2];
					electrolyteIn2[i][j][Height-1]=electrolyteIn2[i][j][Height-2];
					electrolyteIn3[i][j][Height-1]=electrolyteIn3[i][j][Height-2];
					electrolyteIn4[i][j][Height-1]=electrolyteIn4[i][j][Height-2];
					electrolyteIn5[i][j][Height-1]=electrolyteIn5[i][j][Height-2];
					electrolyteIn6[i][j][Height-1]=electrolyteIn6[i][j][Height-2];
					electrolytePotential[i][j][Height-1]=electrolytePotential[i][j][Height-2];	
			}
		}
	}
}

void widthFirstNeumannElectrolytePotential(int Length,int Width,int Height){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (k=0;k<Height;k++){
					electrolyteIn0[i][0][k]=electrolyteIn0[i][1][k];
					electrolyteIn1[i][0][k]=electrolyteIn1[i][1][k];
					electrolyteIn2[i][0][k]=electrolyteIn2[i][1][k];
					electrolyteIn3[i][0][k]=electrolyteIn3[i][1][k];
					electrolyteIn4[i][0][k]=electrolyteIn4[i][1][k];
					electrolyteIn5[i][0][k]=electrolyteIn5[i][1][k];
					electrolyteIn6[i][0][k]=electrolyteIn6[i][1][k];
					electrolyteIn7[i][0][k]=electrolyteIn7[i][1][k];
					electrolyteIn8[i][0][k]=electrolyteIn8[i][1][k];
					electrolyteIn9[i][0][k]=electrolyteIn9[i][1][k];
					electrolyteIn10[i][0][k]=electrolyteIn10[i][1][k];
					electrolyteIn11[i][0][k]=electrolyteIn11[i][1][k];
					electrolyteIn12[i][0][k]=electrolyteIn12[i][1][k];
					electrolyteIn13[i][0][k]=electrolyteIn13[i][1][k];
					electrolyteIn14[i][0][k]=electrolyteIn14[i][1][k];
					electrolyteIn15[i][0][k]=electrolyteIn15[i][1][k];
					electrolyteIn16[i][0][k]=electrolyteIn16[i][1][k];
					electrolyteIn17[i][0][k]=electrolyteIn17[i][1][k];
					electrolyteIn18[i][0][k]=electrolyteIn18[i][1][k];
					electrolytePotential[i][0][k]=electrolytePotential[i][1][k];	
			}
		}
	}
}

void widthFirstDirichletElectrolytePotentialD3Q7(int Length,int Width,int Height, double InletElectrolytePotential, double gasCritical){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (k=0;k<Height;k++){
				if (domain[i][1][k]==0&&rho[i][1][k]>gasCritical){
					electrolytePotential[i][1][k]=InletElectrolytePotential;
					electrolyteIn0[i][1][k] = electrolytePotential[i][1][k] * tNSD3Q7[0];
					electrolyteIn1[i][1][k] = electrolytePotential[i][1][k] * tNSD3Q7[1];
					electrolyteIn2[i][1][k] = electrolytePotential[i][1][k] * tNSD3Q7[2];
					electrolyteIn3[i][1][k] = electrolytePotential[i][1][k] * tNSD3Q7[3];
					electrolyteIn4[i][1][k] = electrolytePotential[i][1][k] * tNSD3Q7[4];
					electrolyteIn5[i][1][k] = electrolytePotential[i][1][k] * tNSD3Q7[5];
					electrolyteIn6[i][1][k] = electrolytePotential[i][1][k] * tNSD3Q7[6];
				}
				if (domain[i][1][k]==1||rho[i][1][k]<=gasCritical){
					electrolytePotential[i][1][k]=0.0;
				}
				electrolytePotential[i][0][k] = electrolytePotential[i][1][k];
				electrolyteIn0[i][0][k] = electrolyteIn0[i][1][k];
				electrolyteIn1[i][0][k] = electrolyteIn1[i][1][k];
				electrolyteIn2[i][0][k] = electrolyteIn2[i][1][k];
				electrolyteIn3[i][0][k] = electrolyteIn3[i][1][k];
				electrolyteIn4[i][0][k] = electrolyteIn4[i][1][k];
				electrolyteIn5[i][0][k] = electrolyteIn5[i][1][k];
				electrolyteIn6[i][0][k] = electrolyteIn6[i][1][k];
			}
		}
	}
}

void widthFirstConstantCurrentElectrolytePotentialD3Q7(int Length,int Width,int Initial, int Height, double constant, double gasCritical){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (k=Initial;k<Height;k++){
				if (domain[i][0][k]==0&&rho[i][0][k]>gasCritical){
					electrolyteIn0[i][0][k] = electrolyteIn0[i][0][k] + tNSD3Q7[0] * constant;
					electrolyteIn1[i][0][k] = electrolyteIn1[i][0][k] + tNSD3Q7[1] * constant;
					electrolyteIn2[i][0][k] = Out5[i][0][k]+tNSD3Q7[2]*constant;
					electrolyteIn3[i][0][k] = electrolyteIn3[i][0][k] + tNSD3Q7[3] * constant;
					electrolyteIn4[i][0][k] = electrolyteIn4[i][0][k] + tNSD3Q7[4] * constant;
					electrolyteIn5[i][0][k] = electrolyteIn5[i][0][k] + tNSD3Q7[5] * constant;
					electrolyteIn6[i][0][k] = electrolyteIn6[i][0][k] + tNSD3Q7[6] * constant;
					electrolytePotential[i][0][k] = electrolyteIn0[i][0][k] + electrolyteIn1[i][0][k] + electrolyteIn2[i][0][k] + electrolyteIn3[i][0][k] + electrolyteIn4[i][0][k] + electrolyteIn5[i][0][k] + electrolyteIn6[i][0][k];
				}
				if (domain[i][0][k] == 1 || domain[i][0][k] == 3 || rho[i][0][k] <= gasCritical){
					electrolytePotential[i][0][k]=0.0;
				}
			}
		}
	}
}

void widthFirstConstantCurrentConcentrationD3Q7(int Length, int Width, int Initial, int Height, double constant, double gasCritical, double transferenceNumber){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (k = Initial; k<Height; k++){
				if (domain[i][0][k] == 0 && rho[i][0][k]>gasCritical){
					cIn0[i][0][k] = cIn0[i][0][k] + tNSD3Q7[0] * constant*(1 - transferenceNumber);
					cIn1[i][0][k] = cIn1[i][0][k] + tNSD3Q7[1] * constant*(1 - transferenceNumber);
					cIn2[i][0][k] = Out5[i][0][k] + tNSD3Q7[2] * constant*(1 - transferenceNumber);
					cIn3[i][0][k] = cIn3[i][0][k] + tNSD3Q7[3] * constant*(1 - transferenceNumber);
					cIn4[i][0][k] = cIn4[i][0][k] + tNSD3Q7[4] * constant*(1 - transferenceNumber);
					cIn5[i][0][k] = cIn5[i][0][k] + tNSD3Q7[5] * constant*(1 - transferenceNumber);
					cIn6[i][0][k] = cIn6[i][0][k] + tNSD3Q7[6] * constant*(1 - transferenceNumber);
					Concentration[i][0][k] = cIn0[i][0][k] + cIn1[i][0][k] + cIn2[i][0][k] + cIn3[i][0][k] + cIn4[i][0][k] + cIn5[i][0][k] + cIn6[i][0][k];
				}
				if (domain[i][0][k] == 1 || domain[i][0][k] == 3 || rho[i][0][k] <= gasCritical){
					Concentration[i][0][k] = 0.0;
					cIn0[i][0][k] = 0.0;
					cIn1[i][0][k] = 0.0;
					cIn2[i][0][k] = 0.0;
					cIn3[i][0][k] = 0.0;
					cIn4[i][0][k] = 0.0;
					cIn5[i][0][k] = 0.0;
					cIn6[i][0][k] = 0.0;
				}
			}
		}
	}
}

void widthLastElectrodeConcentrationNeumannD3Q7(int Length, int Width, int Height){
	int i;
	int k;
#pragma omp parallel private(i,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (k = 0; k<Height; k++){
				if (domain[i][Width-6][k] == 1 && electrodeConnection[i][Width-7][k] == 1){
					cIn0First[i][Width - 6][k] = cIn0First[i][Width - 7][k];
					cIn1First[i][Width - 6][k] = cIn1First[i][Width - 7][k];
					cIn2First[i][Width - 6][k] = cIn2First[i][Width - 7][k];
					cIn3First[i][Width - 6][k] = cIn3First[i][Width - 7][k];
					cIn4First[i][Width - 6][k] = cIn4First[i][Width - 7][k];
					cIn5First[i][Width - 6][k] = cIn5First[i][Width - 7][k];
					cIn6First[i][Width - 6][k] = cIn6First[i][Width - 7][k];
					ConcentrationFirst[i][Width - 6][k] = ConcentrationFirst[i][Width - 7][k];
				}
			}
		}
	}
}

void lengthFirstNeumannElectrolytePotentialD3Q7(int Length,int Width,int Height){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
					electrolyteIn0[0][j][k]=electrolyteIn0[1][j][k];
					electrolyteIn1[0][j][k]=electrolyteIn1[1][j][k];
					electrolyteIn2[0][j][k]=electrolyteIn2[1][j][k];
					electrolyteIn3[0][j][k]=electrolyteIn3[1][j][k];
					electrolyteIn4[0][j][k]=electrolyteIn4[1][j][k];
					electrolyteIn5[0][j][k]=electrolyteIn5[1][j][k];
					electrolyteIn6[0][j][k]=electrolyteIn6[1][j][k];
					electrolytePotential[0][j][k]=electrolytePotential[1][j][k];	
			}
		}
	}
}

void lengthLastNeumannElectrolytePotentialD3Q7(int Length,int Width,int Height){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
					electrolyteIn0[Length-1][j][k]=electrolyteIn0[Length-2][j][k];
					electrolyteIn1[Length-1][j][k]=electrolyteIn1[Length-2][j][k];
					electrolyteIn2[Length-1][j][k]=electrolyteIn2[Length-2][j][k];
					electrolyteIn3[Length-1][j][k]=electrolyteIn3[Length-2][j][k];
					electrolyteIn4[Length-1][j][k]=electrolyteIn4[Length-2][j][k];
					electrolyteIn5[Length-1][j][k]=electrolyteIn5[Length-2][j][k];
					electrolyteIn6[Length-1][j][k]=electrolyteIn6[Length-2][j][k];
					electrolytePotential[Length-1][j][k]=electrolytePotential[Length-2][j][k];	
			}
		}
	}
}

void lengthFirstNeumannConcentrationD3Q7(int Length,int Width,int Height){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
					cIn0[0][j][k]=cIn0[1][j][k];
					cIn1[0][j][k]=cIn1[1][j][k];
					cIn2[0][j][k]=cIn2[1][j][k];
					cIn3[0][j][k]=cIn3[1][j][k];
					cIn4[0][j][k]=cIn4[1][j][k];
					cIn5[0][j][k]=cIn5[1][j][k];
					cIn6[0][j][k]=cIn6[1][j][k];
					Concentration[0][j][k]=Concentration[1][j][k];	
			}
		}
	}
}

void lengthLastNeumannConcentrationD3Q7(int Length,int Width,int Height){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
					cIn0[Length-1][j][k]=cIn0[Length-2][j][k];
					cIn1[Length-1][j][k]=cIn1[Length-2][j][k];
					cIn2[Length-1][j][k]=cIn2[Length-2][j][k];
					cIn3[Length-1][j][k]=cIn3[Length-2][j][k];
					cIn4[Length-1][j][k]=cIn4[Length-2][j][k];
					cIn5[Length-1][j][k]=cIn5[Length-2][j][k];
					cIn6[Length-1][j][k]=cIn6[Length-2][j][k];
					Concentration[Length-1][j][k]=Concentration[Length-2][j][k];	
			}
		}
	}
}

void lengthFirstNeumannConcentrationFirstD3Q7(int Length,int Width,int Height){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
					cIn0First[0][j][k]=cIn0First[1][j][k];
					cIn1First[0][j][k]=cIn1First[1][j][k];
					cIn2First[0][j][k]=cIn2First[1][j][k];
					cIn3First[0][j][k]=cIn3First[1][j][k];
					cIn4First[0][j][k]=cIn4First[1][j][k];
					cIn5First[0][j][k]=cIn5First[1][j][k];
					cIn6First[0][j][k]=cIn6First[1][j][k];
					ConcentrationFirst[0][j][k]=ConcentrationFirst[1][j][k];	
			}
		}
	}
}

void lengthLastNeumannConcentrationFirstD3Q7(int Length,int Width,int Height){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
					cIn0First[Length-1][j][k]=cIn0First[Length-2][j][k];
					cIn1First[Length-1][j][k]=cIn1First[Length-2][j][k];
					cIn2First[Length-1][j][k]=cIn2First[Length-2][j][k];
					cIn3First[Length-1][j][k]=cIn3First[Length-2][j][k];
					cIn4First[Length-1][j][k]=cIn4First[Length-2][j][k];
					cIn5First[Length-1][j][k]=cIn5First[Length-2][j][k];
					cIn6First[Length-1][j][k]=cIn6First[Length-2][j][k];
					ConcentrationFirst[Length-1][j][k]=ConcentrationFirst[Length-2][j][k];	
			}
		}
	}
}

void lengthFirstNeumannConcentrationSecondD3Q7(int Length,int Width,int Height){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
					cIn0Second[0][j][k]=cIn0Second[1][j][k];
					cIn1Second[0][j][k]=cIn1Second[1][j][k];
					cIn2Second[0][j][k]=cIn2Second[1][j][k];
					cIn3Second[0][j][k]=cIn3Second[1][j][k];
					cIn4Second[0][j][k]=cIn4Second[1][j][k];
					cIn5Second[0][j][k]=cIn5Second[1][j][k];
					cIn6Second[0][j][k]=cIn6Second[1][j][k];
					ConcentrationSecond[0][j][k]=ConcentrationSecond[1][j][k];	
			}
		}
	}
}

void lengthLastNeumannConcentrationSecondD3Q7(int Length,int Width,int Height){
	int j;
	int k;
#pragma omp parallel private(j,k)
	{
#pragma omp for schedule(dynamic)
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
					cIn0Second[Length-1][j][k]=cIn0Second[Length-2][j][k];
					cIn1Second[Length-1][j][k]=cIn1Second[Length-2][j][k];
					cIn2Second[Length-1][j][k]=cIn2Second[Length-2][j][k];
					cIn3Second[Length-1][j][k]=cIn3Second[Length-2][j][k];
					cIn4Second[Length-1][j][k]=cIn4Second[Length-2][j][k];
					cIn5Second[Length-1][j][k]=cIn5Second[Length-2][j][k];
					cIn6Second[Length-1][j][k]=cIn6Second[Length-2][j][k];
					ConcentrationSecond[Length-1][j][k]=ConcentrationSecond[Length-2][j][k];	
			}
		}
	}
}