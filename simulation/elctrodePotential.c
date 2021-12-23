#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "domain.h"
#include "MemoryArrange.h"
#include "OutMemoryArrange.h"
#include "matrixMove.h"
#include "memoryFree.h"
#include "twoPhaseField.h"
#include "concentrationField.h"
#include "electrolytePotential.h"

double ***electrodeIn0=NULL;
double ***electrodeIn1=NULL;
double ***electrodeIn2=NULL;
double ***electrodeIn3=NULL;
double ***electrodeIn4=NULL;
double ***electrodeIn5=NULL;
double ***electrodeIn6=NULL;
double ***electrodeIn7=NULL;
double ***electrodeIn8=NULL;
double ***electrodeIn9=NULL;
double ***electrodeIn10=NULL;
double ***electrodeIn11=NULL;
double ***electrodeIn12=NULL;
double ***electrodeIn13=NULL;
double ***electrodeIn14=NULL;
double ***electrodeIn15=NULL;
double ***electrodeIn16=NULL;
double ***electrodeIn17=NULL;
double ***electrodeIn18=NULL;
double ***electrodePotential=NULL;
double ***electrodeConnection=NULL;
double ***lastWallElectrode=NULL;
double ***localCurrent=NULL;
double ***reaction1 = NULL;
double ***reaction2 = NULL;
double ***reaction3 = NULL;
double ***reaction4 = NULL;
double ***reaction5 = NULL;
double ***reaction6 = NULL;

void FieldArrangeElectrodePotential(int Length, int Width, int Height){
	electrodeIn0=memoryarrange(Length,Width,Height);
	electrodeIn1=memoryarrange(Length,Width,Height);
	electrodeIn2=memoryarrange(Length,Width,Height);
	electrodeIn3=memoryarrange(Length,Width,Height);
	electrodeIn4=memoryarrange(Length,Width,Height);
	electrodeIn5=memoryarrange(Length,Width,Height);
	electrodeIn6=memoryarrange(Length,Width,Height);
//	electrodeIn7=memoryarrange(Length,Width,Height);
//	electrodeIn8=memoryarrange(Length,Width,Height);
//	electrodeIn9=memoryarrange(Length,Width,Height);
//	electrodeIn10=memoryarrange(Length,Width,Height);
//	electrodeIn11=memoryarrange(Length,Width,Height);
//	electrodeIn12=memoryarrange(Length,Width,Height);
//	electrodeIn13=memoryarrange(Length,Width,Height);
//	electrodeIn14=memoryarrange(Length,Width,Height);
//	electrodeIn15=memoryarrange(Length,Width,Height);
//	electrodeIn16=memoryarrange(Length,Width,Height);
//	electrodeIn17=memoryarrange(Length,Width,Height);
//	electrodeIn18=memoryarrange(Length,Width,Height);
	electrodePotential=memoryarrange(Length,Width,Height); 
	electrodeConnection=memoryarrange(Length,Width,Height);
	lastWallElectrode=memoryarrange(Length,Width,Height);
	localCurrent = memoryarrange(Length, Width, Height);
	reaction1 = memoryarrange(Length, Width, Height);
	reaction2 = memoryarrange(Length, Width, Height);
	reaction3 = memoryarrange(Length, Width, Height);
	reaction4 = memoryarrange(Length, Width, Height);
	reaction5 = memoryarrange(Length, Width, Height);
	reaction6 = memoryarrange(Length, Width, Height);
}

void electrodeConnectionRegion(int Length, int Width, int Height){
	int i;
	int j;
	int k;
	int k1;
	int k2;
	int Neightbour1;
	int Neightbour2;
	int Neightbour3;
	int Neightbour4;
	int Neightbour5;
	int Neightbour6;
	int Neightbour7;
	int Neightbour8;
	int Neightbour9;
	int Neightbour10;
	int Neightbour11;
	int Neightbour12;
	int Neightbour13;
	int Neightbour14;
	int Neightbour15;
	int Neightbour16;
	int Neightbour17;
	int Neightbour18;
	for (i=0;i<Length;i++){
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
				electrodeConnection[i][j][k]=0;
			}
		}
	}
	electrodeConnection[0][Width-1][0]=1;
	k2=0;
	k1=1;
	while (k1-k2>0){
		k2=k1;
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
					if (electrodeConnection[i][j][k]==0&&domain[i][j][k]==1){
						Neightbour1=neightbour1(i,j,k,Length,electrodeConnection);
						Neightbour2=neightbour2(i,j,k,Width,electrodeConnection);
						Neightbour3=neightbour3(i,j,k,Height,electrodeConnection);
						Neightbour4=neightbour4(i,j,k,Length,electrodeConnection);
						Neightbour5=neightbour5(i,j,k,Width,electrodeConnection);
						Neightbour6=neightbour6(i,j,k,Height,electrodeConnection);
						Neightbour7=neightbour7(i,j,k,Length,Width,electrodeConnection);
						Neightbour8=neightbour8(i,j,k,Length,Height,electrodeConnection);
						Neightbour9=neightbour9(i,j,k,Width,Height,electrodeConnection);
						Neightbour10=neightbour10(i,j,k,Length,Width,electrodeConnection);
						Neightbour11=neightbour11(i,j,k,Length,Height,electrodeConnection);
						Neightbour12=neightbour12(i,j,k,Width,Height,electrodeConnection);
						Neightbour13=neightbour13(i,j,k,Length,Width,electrodeConnection);
						Neightbour14=neightbour14(i,j,k,Length,Height,electrodeConnection);
						Neightbour15=neightbour15(i,j,k,Width,Height,electrodeConnection);
						Neightbour16=neightbour16(i,j,k,Length,Width,electrodeConnection);
						Neightbour17=neightbour17(i,j,k,Length,Height,electrodeConnection);
						Neightbour18=neightbour18(i,j,k,Width,Height,electrodeConnection);
						if (Neightbour1==1||Neightbour2==1||Neightbour3==1||Neightbour4==1||Neightbour5==1||Neightbour6==1||Neightbour7==1||Neightbour8==1||\
							Neightbour9==1||Neightbour10==1||Neightbour11==1||Neightbour12==1||Neightbour13==1||Neightbour14==1||Neightbour15==1||Neightbour16==1||Neightbour17==1||Neightbour18==1){
								electrodeConnection[i][j][k]=1;
								k1=k1+1;
						}
					}
				}
			}
		}
	}
}

void electrodeConnectionRegionD3Q7(int Length, int Width, int Height){
	int i;
	int j;
	int k;
	int k1;
	int k2;
	int Neightbour1;
	int Neightbour2;
	int Neightbour3;
	int Neightbour4;
	int Neightbour5;
	int Neightbour6;
	for (i=0;i<Length;i++){
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
				electrodeConnection[i][j][k]=0;
			}
		}
	}
	electrodeConnection[Length/2][Width-1][Height/2]=1;
	k2=0;
	k1=1;
	while (k1-k2>0){
		k2=k1;
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
					if (electrodeConnection[i][j][k] == 0 && (domain[i][j][k] == 1 || domain[i][j][k] == 2)){
						Neightbour1=neightbour1(i,j,k,Length,electrodeConnection);
						Neightbour2=neightbour2(i,j,k,Width,electrodeConnection);
						Neightbour3=neightbour3(i,j,k,Height,electrodeConnection);
						Neightbour4=neightbour4(i,j,k,Length,electrodeConnection);
						Neightbour5=neightbour5(i,j,k,Width,electrodeConnection);
						Neightbour6=neightbour6(i,j,k,Height,electrodeConnection);
						if (Neightbour1==1||Neightbour2==1||Neightbour3==1||Neightbour4==1||Neightbour5==1||Neightbour6==1){
								electrodeConnection[i][j][k]=1;
								k1=k1+1;
						}
					}
				}
			}
		}
	}
}

void FieldInitialelectrodePotential(int Length, int Width, int Height, double electrodePotentialInitial){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
					if (domain[i][j][k]==1&&electrodeConnection[i][j][k]==1){
						electrodePotential[i][j][k]=electrodePotentialInitial;
						electrodeIn0[i][j][k]=electrodePotential[i][j][k]*tNS[0];
						electrodeIn1[i][j][k]=electrodePotential[i][j][k]*tNS[1];
						electrodeIn2[i][j][k]=electrodePotential[i][j][k]*tNS[2];
						electrodeIn3[i][j][k]=electrodePotential[i][j][k]*tNS[3];
						electrodeIn4[i][j][k]=electrodePotential[i][j][k]*tNS[4];
						electrodeIn5[i][j][k]=electrodePotential[i][j][k]*tNS[5];
						electrodeIn6[i][j][k]=electrodePotential[i][j][k]*tNS[6];
						electrodeIn7[i][j][k]=electrodePotential[i][j][k]*tNS[7];
						electrodeIn8[i][j][k]=electrodePotential[i][j][k]*tNS[8];
						electrodeIn9[i][j][k]=electrodePotential[i][j][k]*tNS[9];
						electrodeIn10[i][j][k]=electrodePotential[i][j][k]*tNS[10];
						electrodeIn11[i][j][k]=electrodePotential[i][j][k]*tNS[11];
						electrodeIn12[i][j][k]=electrodePotential[i][j][k]*tNS[12];
						electrodeIn13[i][j][k]=electrodePotential[i][j][k]*tNS[13];
						electrodeIn14[i][j][k]=electrodePotential[i][j][k]*tNS[14];
						electrodeIn15[i][j][k]=electrodePotential[i][j][k]*tNS[15];
						electrodeIn16[i][j][k]=electrodePotential[i][j][k]*tNS[16];
						electrodeIn17[i][j][k]=electrodePotential[i][j][k]*tNS[17];
						electrodeIn18[i][j][k]=electrodePotential[i][j][k]*tNS[18];
					}
					if (domain[i][j][k]!=1||electrodeConnection[i][j][k]!=1){
						electrodePotential[i][j][k]=0.0;
						electrodeIn0[i][j][k]=0.0;
						electrodeIn1[i][j][k]=0.0;
						electrodeIn2[i][j][k]=0.0;
						electrodeIn3[i][j][k]=0.0;
						electrodeIn4[i][j][k]=0.0;
						electrodeIn5[i][j][k]=0.0;
						electrodeIn6[i][j][k]=0.0;
						electrodeIn7[i][j][k]=0.0;
						electrodeIn8[i][j][k]=0.0;
						electrodeIn9[i][j][k]=0.0;
						electrodeIn10[i][j][k]=0.0;
						electrodeIn11[i][j][k]=0.0;
						electrodeIn12[i][j][k]=0.0;
						electrodeIn13[i][j][k]=0.0;
						electrodeIn14[i][j][k]=0.0;
						electrodeIn15[i][j][k]=0.0;
						electrodeIn16[i][j][k]=0.0;
						electrodeIn17[i][j][k]=0.0;
						electrodeIn18[i][j][k]=0.0;
					}
				}
			}
		}
	}
}

void FieldInitialelectrodePotentialD3Q7(int Length, int Width, int Initial, int Height, double electrodePotentialInitial){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=Initial;k<Height;k++){
					if ((domain[i][j][k] == 1 || domain[i][j][k] == 2) && electrodeConnection[i][j][k] == 1){
						electrodePotential[i][j][k]=electrodePotentialInitial;
						electrodeIn0[i][j][k]=electrodePotential[i][j][k]*tNSD3Q7[0];
						electrodeIn1[i][j][k]=electrodePotential[i][j][k]*tNSD3Q7[1];
						electrodeIn2[i][j][k]=electrodePotential[i][j][k]*tNSD3Q7[2];
						electrodeIn3[i][j][k]=electrodePotential[i][j][k]*tNSD3Q7[3];
						electrodeIn4[i][j][k]=electrodePotential[i][j][k]*tNSD3Q7[4];
						electrodeIn5[i][j][k]=electrodePotential[i][j][k]*tNSD3Q7[5];
						electrodeIn6[i][j][k]=electrodePotential[i][j][k]*tNSD3Q7[6];
					}
					if (domain[i][j][k] == 0 || electrodeConnection[i][j][k] != 1 || domain[i][j][k] == 3){
						electrodePotential[i][j][k]=0.0;
						electrodeIn0[i][j][k]=0.0;
						electrodeIn1[i][j][k]=0.0;
						electrodeIn2[i][j][k]=0.0;
						electrodeIn3[i][j][k]=0.0;
						electrodeIn4[i][j][k]=0.0;
						electrodeIn5[i][j][k]=0.0;
						electrodeIn6[i][j][k]=0.0;
					}
				}
			}
		}
	}
}

void CollisionSRTElectrodePotential(int Length,int Width,int Height,double CrelaxationTime){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
					if (domain[i][j][k]==1&&electrodeConnection[i][j][k]==1){
						Out0[i][j][k]=electrodeIn0[i][j][k]-CrelaxationTime*(electrodeIn0[i][j][k]-tNS[0]*electrodePotential[i][j][k]);
						Out1[i][j][k]=electrodeIn1[i][j][k]-CrelaxationTime*(electrodeIn1[i][j][k]-tNS[1]*electrodePotential[i][j][k]);
						Out2[i][j][k]=electrodeIn2[i][j][k]-CrelaxationTime*(electrodeIn2[i][j][k]-tNS[2]*electrodePotential[i][j][k]);
						Out3[i][j][k]=electrodeIn3[i][j][k]-CrelaxationTime*(electrodeIn3[i][j][k]-tNS[3]*electrodePotential[i][j][k]);
						Out4[i][j][k]=electrodeIn4[i][j][k]-CrelaxationTime*(electrodeIn4[i][j][k]-tNS[4]*electrodePotential[i][j][k]);					
						Out5[i][j][k]=electrodeIn5[i][j][k]-CrelaxationTime*(electrodeIn5[i][j][k]-tNS[5]*electrodePotential[i][j][k]);
						Out6[i][j][k]=electrodeIn6[i][j][k]-CrelaxationTime*(electrodeIn6[i][j][k]-tNS[6]*electrodePotential[i][j][k]);
						Out7[i][j][k]=electrodeIn7[i][j][k]-CrelaxationTime*(electrodeIn7[i][j][k]-tNS[7]*electrodePotential[i][j][k]);
						Out8[i][j][k]=electrodeIn8[i][j][k]-CrelaxationTime*(electrodeIn8[i][j][k]-tNS[8]*electrodePotential[i][j][k]);
						Out9[i][j][k]=electrodeIn9[i][j][k]-CrelaxationTime*(electrodeIn9[i][j][k]-tNS[9]*electrodePotential[i][j][k]);
						Out10[i][j][k]=electrodeIn10[i][j][k]-CrelaxationTime*(electrodeIn10[i][j][k]-tNS[10]*electrodePotential[i][j][k]);
						Out11[i][j][k]=electrodeIn11[i][j][k]-CrelaxationTime*(electrodeIn11[i][j][k]-tNS[11]*electrodePotential[i][j][k]);
						Out12[i][j][k]=electrodeIn12[i][j][k]-CrelaxationTime*(electrodeIn12[i][j][k]-tNS[12]*electrodePotential[i][j][k]);
						Out13[i][j][k]=electrodeIn13[i][j][k]-CrelaxationTime*(electrodeIn13[i][j][k]-tNS[13]*electrodePotential[i][j][k]);
						Out14[i][j][k]=electrodeIn14[i][j][k]-CrelaxationTime*(electrodeIn14[i][j][k]-tNS[14]*electrodePotential[i][j][k]);
						Out15[i][j][k]=electrodeIn15[i][j][k]-CrelaxationTime*(electrodeIn15[i][j][k]-tNS[15]*electrodePotential[i][j][k]);
						Out16[i][j][k]=electrodeIn16[i][j][k]-CrelaxationTime*(electrodeIn16[i][j][k]-tNS[16]*electrodePotential[i][j][k]);
						Out17[i][j][k]=electrodeIn17[i][j][k]-CrelaxationTime*(electrodeIn17[i][j][k]-tNS[17]*electrodePotential[i][j][k]);
						Out18[i][j][k]=electrodeIn18[i][j][k]-CrelaxationTime*(electrodeIn18[i][j][k]-tNS[18]*electrodePotential[i][j][k]);
					}
				}
			}
		}
	}
}

void CollisionSRTElectrodePotentialD3Q7(int Length,int Width,int Height,double CrelaxationTime){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
					if ((domain[i][j][k] == 1||domain[i][j][k] == 2) && electrodeConnection[i][j][k] == 1){
						Out0[i][j][k]=electrodeIn0[i][j][k]-CrelaxationTime*(electrodeIn0[i][j][k]-tNSD3Q7[0]*electrodePotential[i][j][k]);
						Out1[i][j][k]=electrodeIn1[i][j][k]-CrelaxationTime*(electrodeIn1[i][j][k]-tNSD3Q7[1]*electrodePotential[i][j][k]);
						Out2[i][j][k]=electrodeIn2[i][j][k]-CrelaxationTime*(electrodeIn2[i][j][k]-tNSD3Q7[2]*electrodePotential[i][j][k]);
						Out3[i][j][k]=electrodeIn3[i][j][k]-CrelaxationTime*(electrodeIn3[i][j][k]-tNSD3Q7[3]*electrodePotential[i][j][k]);
						Out4[i][j][k]=electrodeIn4[i][j][k]-CrelaxationTime*(electrodeIn4[i][j][k]-tNSD3Q7[4]*electrodePotential[i][j][k]);					
						Out5[i][j][k]=electrodeIn5[i][j][k]-CrelaxationTime*(electrodeIn5[i][j][k]-tNSD3Q7[5]*electrodePotential[i][j][k]);
						Out6[i][j][k]=electrodeIn6[i][j][k]-CrelaxationTime*(electrodeIn6[i][j][k]-tNSD3Q7[6]*electrodePotential[i][j][k]);
					}
				}
			}
		}
	}
}

void MemoryFreeElectrodePotential(int m, int n, int q){
	freememory(electrodePotential,m,n,q);
	freememory(electrodeConnection,m,n,q);
	freememory(electrodeIn0,m,n,q);
	freememory(electrodeIn1,m,n,q);
	freememory(electrodeIn2,m,n,q);
	freememory(electrodeIn3,m,n,q);
	freememory(electrodeIn4,m,n,q);
	freememory(electrodeIn5,m,n,q);
	freememory(electrodeIn6,m,n,q);
//	freememory(electrodeIn7,m,n,q);
//	freememory(electrodeIn8,m,n,q);
//	freememory(electrodeIn9,m,n,q);
//	freememory(electrodeIn10,m,n,q);
//	freememory(electrodeIn11,m,n,q);
//	freememory(electrodeIn12,m,n,q);
//	freememory(electrodeIn13,m,n,q);
//	freememory(electrodeIn14,m,n,q);
//	freememory(electrodeIn15,m,n,q);
//	freememory(electrodeIn16,m,n,q);
//	freememory(electrodeIn17,m,n,q);
//	freememory(electrodeIn18,m,n,q);
	freememory(lastWallElectrode,m,n,q);
	freememory(localCurrent,m,n,q);
	freememory(reaction1, m, n, q);
	freememory(reaction2, m, n, q);
	freememory(reaction3, m, n, q);
	freememory(reaction4, m, n, q);
	freememory(reaction5, m, n, q);
	freememory(reaction6, m, n, q);
}

void StreamElectrodePotential(int m, int n, int q){
	shift0(Out0,electrodeIn0,m,n,q);
	shift1(Out1,electrodeIn1,m,n,q);
	shift2(Out2,electrodeIn2,m,n,q);
	shift3(Out3,electrodeIn3,m,n,q);
	shift4(Out4,electrodeIn4,m,n,q);
	shift5(Out5,electrodeIn5,m,n,q);
	shift6(Out6,electrodeIn6,m,n,q);
	shift7(Out7,electrodeIn7,m,n,q);
	shift8(Out8,electrodeIn8,m,n,q);
	shift9(Out9,electrodeIn9,m,n,q);
	shift10(Out10,electrodeIn10,m,n,q);
	shift11(Out11,electrodeIn11,m,n,q);
	shift12(Out12,electrodeIn12,m,n,q);
	shift13(Out13,electrodeIn13,m,n,q);
	shift14(Out14,electrodeIn14,m,n,q);
	shift15(Out15,electrodeIn15,m,n,q);
	shift16(Out16,electrodeIn16,m,n,q);
	shift17(Out17,electrodeIn17,m,n,q);
	shift18(Out18,electrodeIn18,m,n,q);
}

void StreamElectrodePotentialD3Q7(int m, int n, int q){
	shift0(Out0,electrodeIn0,m,n,q);
	shift1(Out1,electrodeIn1,m,n,q);
	shift2(Out2,electrodeIn2,m,n,q);
	shift3(Out3,electrodeIn3,m,n,q);
	shift4(Out4,electrodeIn4,m,n,q);
	shift5(Out5,electrodeIn5,m,n,q);
	shift6(Out6,electrodeIn6,m,n,q);
}

void FieldCalculationElectrodePotential(int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==1&&electrodeConnection[i][j][k]==1){
						electrodePotential[i][j][k]=electrodeIn0[i][j][k]+electrodeIn1[i][j][k]+electrodeIn2[i][j][k]+electrodeIn3[i][j][k]+electrodeIn4[i][j][k]+electrodeIn5[i][j][k]+electrodeIn6[i][j][k]+electrodeIn7[i][j][k]+electrodeIn8[i][j][k]\
							+electrodeIn9[i][j][k]+electrodeIn10[i][j][k]+electrodeIn11[i][j][k]+electrodeIn12[i][j][k]+electrodeIn13[i][j][k]+electrodeIn14[i][j][k]+electrodeIn15[i][j][k]+electrodeIn16[i][j][k]+electrodeIn17[i][j][k]+electrodeIn18[i][j][k];
					}
					if (domain[i][j][k]!=1){
						electrodePotential[i][j][k]=0.0;
						electrodeIn0[i][j][k]=0.0;
						electrodeIn1[i][j][k]=0.0;
						electrodeIn2[i][j][k]=0.0;
						electrodeIn3[i][j][k]=0.0;
						electrodeIn4[i][j][k]=0.0;
						electrodeIn5[i][j][k]=0.0;
						electrodeIn6[i][j][k]=0.0;
						electrodeIn7[i][j][k]=0.0;
						electrodeIn8[i][j][k]=0.0;
						electrodeIn9[i][j][k]=0.0;
						electrodeIn10[i][j][k]=0.0;
						electrodeIn11[i][j][k]=0.0;
						electrodeIn12[i][j][k]=0.0;
						electrodeIn13[i][j][k]=0.0;
						electrodeIn14[i][j][k]=0.0;
						electrodeIn15[i][j][k]=0.0;
						electrodeIn16[i][j][k]=0.0;
						electrodeIn17[i][j][k]=0.0;
						electrodeIn18[i][j][k]=0.0;
					}
					if (domain[i][j][k]==0&&electrodeConnection[i][j][k]!=1){
						electrodeIn0[i][j][k]=electrodePotential[i][j][k]*tNS[0];
						electrodeIn1[i][j][k]=electrodePotential[i][j][k]*tNS[1];
						electrodeIn2[i][j][k]=electrodePotential[i][j][k]*tNS[2];
						electrodeIn3[i][j][k]=electrodePotential[i][j][k]*tNS[3];
						electrodeIn4[i][j][k]=electrodePotential[i][j][k]*tNS[4];
						electrodeIn5[i][j][k]=electrodePotential[i][j][k]*tNS[5];
						electrodeIn6[i][j][k]=electrodePotential[i][j][k]*tNS[6];
						electrodeIn7[i][j][k]=electrodePotential[i][j][k]*tNS[7];
						electrodeIn8[i][j][k]=electrodePotential[i][j][k]*tNS[8];
						electrodeIn9[i][j][k]=electrodePotential[i][j][k]*tNS[9];
						electrodeIn10[i][j][k]=electrodePotential[i][j][k]*tNS[10];
						electrodeIn11[i][j][k]=electrodePotential[i][j][k]*tNS[11];
						electrodeIn12[i][j][k]=electrodePotential[i][j][k]*tNS[12];
						electrodeIn13[i][j][k]=electrodePotential[i][j][k]*tNS[13];
						electrodeIn14[i][j][k]=electrodePotential[i][j][k]*tNS[14];
						electrodeIn15[i][j][k]=electrodePotential[i][j][k]*tNS[15];
						electrodeIn16[i][j][k]=electrodePotential[i][j][k]*tNS[16];
						electrodeIn17[i][j][k]=electrodePotential[i][j][k]*tNS[17];
						electrodeIn18[i][j][k]=electrodePotential[i][j][k]*tNS[18];
					}					
				}
			}
		}
	}
}

void FieldCalculationElectrodePotentialD3Q7(int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if ((domain[i][j][k] == 1 || domain[i][j][k] == 2) && electrodeConnection[i][j][k] == 1){
						electrodePotential[i][j][k]=electrodeIn0[i][j][k]+electrodeIn1[i][j][k]+electrodeIn2[i][j][k]+electrodeIn3[i][j][k]+electrodeIn4[i][j][k]+electrodeIn5[i][j][k]+electrodeIn6[i][j][k];
					}
					if (domain[i][j][k] == 0 || domain[i][j][k] == 3 || electrodeConnection[i][j][k] != 1){
						electrodePotential[i][j][k]=0.0;
						electrodeIn0[i][j][k]=0.0;
						electrodeIn1[i][j][k]=0.0;
						electrodeIn2[i][j][k]=0.0;
						electrodeIn3[i][j][k]=0.0;
						electrodeIn4[i][j][k]=0.0;
						electrodeIn5[i][j][k]=0.0;
						electrodeIn6[i][j][k]=0.0;
					}
/*					if (domain[i][j][k]==0&&electrodeConnection[i][j][k]!=1){
						electrodeIn0[i][j][k]=electrodePotential[i][j][k]*tNSD3Q7[0];
						electrodeIn1[i][j][k]=electrodePotential[i][j][k]*tNSD3Q7[1];
						electrodeIn2[i][j][k]=electrodePotential[i][j][k]*tNSD3Q7[2];
						electrodeIn3[i][j][k]=electrodePotential[i][j][k]*tNSD3Q7[3];
						electrodeIn4[i][j][k]=electrodePotential[i][j][k]*tNSD3Q7[4];
						electrodeIn5[i][j][k]=electrodePotential[i][j][k]*tNSD3Q7[5];
						electrodeIn6[i][j][k]=electrodePotential[i][j][k]*tNSD3Q7[6];
					}*/					
				}
			}
		}
	}
}

void lastElectrode(int Length, int Width, int Height){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
						lastWallElectrode[i][j][k]=electrodePotential[i][j][k];
				}
			}
		}
	}
}

void reactionElectrodeD3Q7(int m, int n, int q, double gasCritical, double RateConstant, double constant, double K, double exchange, double MaxConcentration, double lastTotalCurrentResidual){
	int i;
	int j;
	int k;
	double valueConcentration;
	double rateConstant;
	double electrolytePotential;
	double overpotential;
	double domainConnectionNeightbour;
	double OCV;
#pragma omp parallel private(i,j,k,valueConcentration,rateConstant,electrolytePotential,overpotential,domainConnectionNeightbour,OCV)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					localCurrent[i][j][k] = 0;
					reaction1[i][j][k] = 0;
					reaction2[i][j][k] = 0;
					reaction3[i][j][k] = 0;
					reaction4[i][j][k] = 0;
					reaction5[i][j][k] = 0;
					reaction6[i][j][k] = 0;
					if (domain[i][j][k]==1&&electrodeConnection[i][j][k]==1){
						if (Domain1[i][j][k] == 0){
							domainConnectionNeightbour=neightbour4(i,j,k,m,domainConnection);
							if (domainConnectionNeightbour==1){
							OCV = 7.976 - 5.5419*(wallNearConcentrationFirst[i][j][k] / MaxConcentration) + 5.2824*pow((wallNearConcentrationFirst[i][j][k] / MaxConcentration), 1.0700)-1.0556*pow(10.0, -4.0)*exp(124.7407*(wallNearConcentrationFirst[i][j][k] / MaxConcentration) - 114.2593) - 4.0446*pow((wallNearConcentrationFirst[i][j][k] / MaxConcentration), 0.0766);
							valueConcentration=neightbour4(i, j, k, m, wallNearConcentration);
							electrolytePotential=neightbour4(i, j, k, m, lastWallElectrolyte);
							overpotential = electrodePotential[i][j][k] / exchange - electrolytePotential / exchange -OCV;
							rateConstant = -RateConstant*(exp(0.5*96485.0*overpotential / 8.314 / K) - exp(-0.5*96485.0*overpotential / 8.314 / K));
							reaction1[i][j][k] = 1.0*constant*rateConstant *pow(valueConcentration, 0.5)*pow(wallNearConcentrationFirst[i][j][k], 0.5)*pow((MaxConcentration - wallNearConcentrationFirst[i][j][k]), 0.5);
							localCurrent[i][j][k] = localCurrent[i][j][k] + reaction1[i][j][k];
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * reaction1[i][j][k] * (1 + lastTotalCurrentResidual);
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * reaction1[i][j][k] * (1 + lastTotalCurrentResidual);
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * reaction1[i][j][k] * (1 + lastTotalCurrentResidual);
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * reaction1[i][j][k] * (1 + lastTotalCurrentResidual);
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * reaction1[i][j][k] * (1 + lastTotalCurrentResidual);
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * reaction1[i][j][k] * (1 + lastTotalCurrentResidual);
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * reaction1[i][j][k] * (1 + lastTotalCurrentResidual);
							}
						}
						if (Domain2[i][j][k] == 0){
							domainConnectionNeightbour=neightbour5(i,j,k,n,domainConnection);
							if (domainConnectionNeightbour==1){
							OCV = 7.976 - 5.5419*(wallNearConcentrationFirst[i][j][k] / MaxConcentration) + 5.2824*pow((wallNearConcentrationFirst[i][j][k] / MaxConcentration), 1.0700)-1.0556*pow(10.0, -4.0)*exp(124.7407*(wallNearConcentrationFirst[i][j][k] / MaxConcentration) - 114.2593) - 4.0446*pow((wallNearConcentrationFirst[i][j][k] / MaxConcentration), 0.0766);
							valueConcentration=neightbour5(i,j,k,n,wallNearConcentration);
							electrolytePotential=neightbour5(i, j, k, n, lastWallElectrolyte);
							overpotential = electrodePotential[i][j][k] / exchange - electrolytePotential / exchange -OCV;
							rateConstant = -RateConstant*(exp(0.5*96485.0*overpotential / 8.314 / K) - exp(-0.5*96485.0*overpotential / 8.314 / K));
							reaction2[i][j][k] = 1.0*constant*rateConstant *pow(valueConcentration, 0.5)*pow(wallNearConcentrationFirst[i][j][k], 0.5)*pow((MaxConcentration - wallNearConcentrationFirst[i][j][k]), 0.5);
							localCurrent[i][j][k] = localCurrent[i][j][k] + reaction2[i][j][k];
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * reaction2[i][j][k] * (1 + lastTotalCurrentResidual);
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * reaction2[i][j][k] * (1 + lastTotalCurrentResidual);
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * reaction2[i][j][k] * (1 + lastTotalCurrentResidual);
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * reaction2[i][j][k] * (1 + lastTotalCurrentResidual);
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * reaction2[i][j][k] * (1 + lastTotalCurrentResidual);
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * reaction2[i][j][k] * (1 + lastTotalCurrentResidual);
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * reaction2[i][j][k] * (1 + lastTotalCurrentResidual);
							}
						}
						if (Domain3[i][j][k] == 0){
							domainConnectionNeightbour=neightbour6(i,j,k,q,domainConnection);
							if (domainConnectionNeightbour==1){
							OCV = 7.976 - 5.5419*(wallNearConcentrationFirst[i][j][k] / MaxConcentration) + 5.2824*pow((wallNearConcentrationFirst[i][j][k] / MaxConcentration), 1.0700)-1.0556*pow(10.0, -4.0)*exp(124.7407*(wallNearConcentrationFirst[i][j][k] / MaxConcentration) - 114.2593) - 4.0446*pow((wallNearConcentrationFirst[i][j][k] / MaxConcentration), 0.0766);
							valueConcentration=neightbour6(i, j, k, q, wallNearConcentration);
							electrolytePotential=neightbour6(i, j, k, q, lastWallElectrolyte);
							overpotential = electrodePotential[i][j][k] / exchange - electrolytePotential / exchange -OCV;
							rateConstant = -RateConstant*(exp(0.5*96485.0*overpotential / 8.314 / K) - exp(-0.5*96485.0*overpotential / 8.314 / K));
							reaction3[i][j][k] = 1.0*constant*rateConstant *pow(valueConcentration, 0.5)*pow(wallNearConcentrationFirst[i][j][k], 0.5)*pow((MaxConcentration - wallNearConcentrationFirst[i][j][k]), 0.5);
							localCurrent[i][j][k] = localCurrent[i][j][k] + reaction3[i][j][k];
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * reaction3[i][j][k] * (1 + lastTotalCurrentResidual);
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * reaction3[i][j][k] * (1 + lastTotalCurrentResidual);
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * reaction3[i][j][k] * (1 + lastTotalCurrentResidual);
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * reaction3[i][j][k] * (1 + lastTotalCurrentResidual);
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * reaction3[i][j][k] * (1 + lastTotalCurrentResidual);
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * reaction3[i][j][k] * (1 + lastTotalCurrentResidual);
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * reaction3[i][j][k] * (1 + lastTotalCurrentResidual);
							}
						}
						if (Domain4[i][j][k] == 0){
							domainConnectionNeightbour=neightbour1(i,j,k,m,domainConnection);
							if (domainConnectionNeightbour==1){
							OCV = 7.976 - 5.5419*(wallNearConcentrationFirst[i][j][k] / MaxConcentration) + 5.2824*pow((wallNearConcentrationFirst[i][j][k] / MaxConcentration), 1.0700)-1.0556*pow(10.0, -4.0)*exp(124.7407*(wallNearConcentrationFirst[i][j][k] / MaxConcentration) - 114.2593) - 4.0446*pow((wallNearConcentrationFirst[i][j][k] / MaxConcentration), 0.0766);
							valueConcentration=neightbour1(i, j, k, m, wallNearConcentration);
							electrolytePotential=neightbour1(i, j, k, m, lastWallElectrolyte);
							overpotential = electrodePotential[i][j][k] / exchange - electrolytePotential / exchange -OCV;
							rateConstant = -RateConstant*(exp(0.5*96485.0*overpotential / 8.314 / K) - exp(-0.5*96485.0*overpotential / 8.314 / K));
							reaction4[i][j][k] = 1.0*constant*rateConstant *pow(valueConcentration, 0.5)*pow(wallNearConcentrationFirst[i][j][k], 0.5)*pow((MaxConcentration - wallNearConcentrationFirst[i][j][k]), 0.5);
							localCurrent[i][j][k] = localCurrent[i][j][k] + reaction4[i][j][k];
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * reaction4[i][j][k] * (1 + lastTotalCurrentResidual);
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * reaction4[i][j][k] * (1 + lastTotalCurrentResidual);
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * reaction4[i][j][k] * (1 + lastTotalCurrentResidual);
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * reaction4[i][j][k] * (1 + lastTotalCurrentResidual);
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * reaction4[i][j][k] * (1 + lastTotalCurrentResidual);
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * reaction4[i][j][k] * (1 + lastTotalCurrentResidual);
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * reaction4[i][j][k] * (1 + lastTotalCurrentResidual);
							}
						}
						if (Domain5[i][j][k] == 0 && j != n - 1){
							domainConnectionNeightbour=neightbour2(i,j,k,n,domainConnection);
							if (domainConnectionNeightbour==1){
							OCV = 7.976 - 5.5419*(wallNearConcentrationFirst[i][j][k] / MaxConcentration) + 5.2824*pow((wallNearConcentrationFirst[i][j][k] / MaxConcentration), 1.0700)-1.0556*pow(10.0, -4.0)*exp(124.7407*(wallNearConcentrationFirst[i][j][k] / MaxConcentration) - 114.2593) - 4.0446*pow((wallNearConcentrationFirst[i][j][k] / MaxConcentration), 0.0766);
							valueConcentration=neightbour2(i, j, k, n, wallNearConcentration);
							electrolytePotential=neightbour2(i, j, k, n, lastWallElectrolyte);
							overpotential = electrodePotential[i][j][k] / exchange - electrolytePotential / exchange -OCV;
							rateConstant = -RateConstant*(exp(0.5*96485.0*overpotential / 8.314 / K) - exp(-0.5*96485.0*overpotential / 8.314 / K));
							reaction5[i][j][k] = 1.0*constant*rateConstant *pow(valueConcentration, 0.5)*pow(wallNearConcentrationFirst[i][j][k], 0.5)*pow((MaxConcentration - wallNearConcentrationFirst[i][j][k]), 0.5);
							localCurrent[i][j][k] = localCurrent[i][j][k] + reaction5[i][j][k];
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * reaction5[i][j][k] * (1 + lastTotalCurrentResidual);
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * reaction5[i][j][k] * (1 + lastTotalCurrentResidual);
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * reaction5[i][j][k] * (1 + lastTotalCurrentResidual);
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * reaction5[i][j][k] * (1 + lastTotalCurrentResidual);
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * reaction5[i][j][k] * (1 + lastTotalCurrentResidual);
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * reaction5[i][j][k] * (1 + lastTotalCurrentResidual);
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * reaction5[i][j][k] * (1 + lastTotalCurrentResidual);
							}
						}
						if (Domain6[i][j][k] == 0){
							domainConnectionNeightbour=neightbour3(i,j,k,q,domainConnection);
							if (domainConnectionNeightbour==1){
							OCV = 7.976 - 5.5419*(wallNearConcentrationFirst[i][j][k] / MaxConcentration) + 5.2824*pow((wallNearConcentrationFirst[i][j][k] / MaxConcentration), 1.0700)-1.0556*pow(10.0, -4.0)*exp(124.7407*(wallNearConcentrationFirst[i][j][k] / MaxConcentration) - 114.2593) - 4.0446*pow((wallNearConcentrationFirst[i][j][k] / MaxConcentration), 0.0766);
							valueConcentration=neightbour3(i, j, k, q, wallNearConcentration);
							electrolytePotential=neightbour3(i, j, k, q, lastWallElectrolyte);
							overpotential = electrodePotential[i][j][k] / exchange - electrolytePotential / exchange -OCV;
							rateConstant = -RateConstant*(exp(0.5*96485.0*overpotential / 8.314 / K) - exp(-0.5*96485.0*overpotential / 8.314 / K));
							reaction6[i][j][k] = 1.0*constant*rateConstant *pow(valueConcentration, 0.5)*pow(wallNearConcentrationFirst[i][j][k], 0.5)*pow((MaxConcentration - wallNearConcentrationFirst[i][j][k]), 0.5);
							localCurrent[i][j][k] = localCurrent[i][j][k] + reaction6[i][j][k];
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * reaction6[i][j][k] * (1 + lastTotalCurrentResidual);
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * reaction6[i][j][k] * (1 + lastTotalCurrentResidual);
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * reaction6[i][j][k] * (1 + lastTotalCurrentResidual);
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * reaction6[i][j][k] * (1 + lastTotalCurrentResidual);
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * reaction6[i][j][k] * (1 + lastTotalCurrentResidual);
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * reaction6[i][j][k] * (1 + lastTotalCurrentResidual);
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * reaction6[i][j][k] * (1 + lastTotalCurrentResidual);
							}
						}
					}
				}
			}
		}
	}
}