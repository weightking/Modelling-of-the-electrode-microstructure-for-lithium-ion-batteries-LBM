#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "MemoryArrange.h"
#include "memoryFree.h"
#include "domain.h"
#include "twoPhaseField.h"

extern double ***cuNS0=NULL;
extern double ***cuNS1=NULL;
extern double ***cuNS2=NULL;
extern double ***cuNS3=NULL;
extern double ***cuNS4=NULL;
extern double ***cuNS5=NULL;
extern double ***cuNS6=NULL;
extern double ***cuNS7=NULL;
extern double ***cuNS8=NULL;
extern double ***cuNS9=NULL;
extern double ***cuNS10=NULL;
extern double ***cuNS11=NULL;
extern double ***cuNS12=NULL;
extern double ***cuNS13=NULL;
extern double ***cuNS14=NULL;
extern double ***cuNS15=NULL;
extern double ***cuNS16=NULL;
extern double ***cuNS17=NULL;
extern double ***cuNS18=NULL;
extern double ***ux2uy2uz2=NULL;

extern double ***Eq0=NULL;
extern double ***Eq1=NULL;
extern double ***Eq2=NULL;
extern double ***Eq3=NULL;
extern double ***Eq4=NULL;
extern double ***Eq5=NULL;
extern double ***Eq6=NULL;
extern double ***Eq7=NULL;
extern double ***Eq8=NULL;
extern double ***Eq9=NULL;
extern double ***Eq10=NULL;
extern double ***Eq11=NULL;
extern double ***Eq12=NULL;
extern double ***Eq13=NULL;
extern double ***Eq14=NULL;
extern double ***Eq15=NULL;
extern double ***Eq16=NULL;
extern double ***Eq17=NULL;
extern double ***Eq18=NULL;

void equilibriumArrange(int Length, int Width, int Height){
	cuNS0=memoryarrange(Length,Width,Height);
	cuNS1=memoryarrange(Length,Width,Height);
	cuNS2=memoryarrange(Length,Width,Height);
	cuNS3=memoryarrange(Length,Width,Height);
	cuNS4=memoryarrange(Length,Width,Height);
	cuNS5=memoryarrange(Length,Width,Height);
	cuNS6=memoryarrange(Length,Width,Height);
	cuNS7=memoryarrange(Length,Width,Height);
	cuNS8=memoryarrange(Length,Width,Height);
	cuNS9=memoryarrange(Length,Width,Height);
	cuNS10=memoryarrange(Length,Width,Height);
	cuNS11=memoryarrange(Length,Width,Height);
	cuNS12=memoryarrange(Length,Width,Height);
	cuNS13=memoryarrange(Length,Width,Height);
	cuNS14=memoryarrange(Length,Width,Height);
	cuNS15=memoryarrange(Length,Width,Height);
	cuNS16=memoryarrange(Length,Width,Height);
	cuNS17=memoryarrange(Length,Width,Height);
	cuNS18=memoryarrange(Length,Width,Height);
	ux2uy2uz2=memoryarrange(Length,Width,Height);
	Eq0=memoryarrange(Length,Width,Height);
	Eq1=memoryarrange(Length,Width,Height);
	Eq2=memoryarrange(Length,Width,Height);
	Eq3=memoryarrange(Length,Width,Height);
	Eq4=memoryarrange(Length,Width,Height);
	Eq5=memoryarrange(Length,Width,Height);
	Eq6=memoryarrange(Length,Width,Height);
	Eq7=memoryarrange(Length,Width,Height);
	Eq8=memoryarrange(Length,Width,Height);
	Eq9=memoryarrange(Length,Width,Height);
	Eq10=memoryarrange(Length,Width,Height);
	Eq11=memoryarrange(Length,Width,Height);
	Eq12=memoryarrange(Length,Width,Height);
	Eq13=memoryarrange(Length,Width,Height);
	Eq14=memoryarrange(Length,Width,Height);
	Eq15=memoryarrange(Length,Width,Height);
	Eq16=memoryarrange(Length,Width,Height);
	Eq17=memoryarrange(Length,Width,Height);
	Eq18=memoryarrange(Length,Width,Height);
}

void equilibriumCalculation(int Length, int Width, int Height, double ***Ux, double ***Uy, double ***Uz){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
					if (domain[i][j][k]==0){
						cuNS0[i][j][k]=(3.0*(cxNS[0]*Ux[i][j][k]+cyNS[0]*Uy[i][j][k]+czNS[0]*Uz[i][j][k]));
						cuNS1[i][j][k]=(3.0*(cxNS[1]*Ux[i][j][k]+cyNS[1]*Uy[i][j][k]+czNS[1]*Uz[i][j][k]));
						cuNS2[i][j][k]=(3.0*(cxNS[2]*Ux[i][j][k]+cyNS[2]*Uy[i][j][k]+czNS[2]*Uz[i][j][k]));
						cuNS3[i][j][k]=(3.0*(cxNS[3]*Ux[i][j][k]+cyNS[3]*Uy[i][j][k]+czNS[3]*Uz[i][j][k]));
						cuNS4[i][j][k]=(3.0*(cxNS[4]*Ux[i][j][k]+cyNS[4]*Uy[i][j][k]+czNS[4]*Uz[i][j][k]));
						cuNS5[i][j][k]=(3.0*(cxNS[5]*Ux[i][j][k]+cyNS[5]*Uy[i][j][k]+czNS[5]*Uz[i][j][k]));
						cuNS6[i][j][k]=(3.0*(cxNS[6]*Ux[i][j][k]+cyNS[6]*Uy[i][j][k]+czNS[6]*Uz[i][j][k]));
						cuNS7[i][j][k]=(3.0*(cxNS[7]*Ux[i][j][k]+cyNS[7]*Uy[i][j][k]+czNS[7]*Uz[i][j][k]));
						cuNS8[i][j][k]=(3.0*(cxNS[8]*Ux[i][j][k]+cyNS[8]*Uy[i][j][k]+czNS[8]*Uz[i][j][k]));
						cuNS9[i][j][k]=(3.0*(cxNS[9]*Ux[i][j][k]+cyNS[9]*Uy[i][j][k]+czNS[9]*Uz[i][j][k]));
						cuNS10[i][j][k]=(3.0*(cxNS[10]*Ux[i][j][k]+cyNS[10]*Uy[i][j][k]+czNS[10]*Uz[i][j][k]));
						cuNS11[i][j][k]=(3.0*(cxNS[11]*Ux[i][j][k]+cyNS[11]*Uy[i][j][k]+czNS[11]*Uz[i][j][k]));
						cuNS12[i][j][k]=(3.0*(cxNS[12]*Ux[i][j][k]+cyNS[12]*Uy[i][j][k]+czNS[12]*Uz[i][j][k]));
						cuNS13[i][j][k]=(3.0*(cxNS[13]*Ux[i][j][k]+cyNS[13]*Uy[i][j][k]+czNS[13]*Uz[i][j][k]));
						cuNS14[i][j][k]=(3.0*(cxNS[14]*Ux[i][j][k]+cyNS[14]*Uy[i][j][k]+czNS[14]*Uz[i][j][k]));
						cuNS15[i][j][k]=(3.0*(cxNS[15]*Ux[i][j][k]+cyNS[15]*Uy[i][j][k]+czNS[15]*Uz[i][j][k]));
						cuNS16[i][j][k]=(3.0*(cxNS[16]*Ux[i][j][k]+cyNS[16]*Uy[i][j][k]+czNS[16]*Uz[i][j][k]));
						cuNS17[i][j][k]=(3.0*(cxNS[17]*Ux[i][j][k]+cyNS[17]*Uy[i][j][k]+czNS[17]*Uz[i][j][k]));
						cuNS18[i][j][k]=(3.0*(cxNS[18]*Ux[i][j][k]+cyNS[18]*Uy[i][j][k]+czNS[18]*Uz[i][j][k]));
						ux2uy2uz2[i][j][k]=Ux[i][j][k]*Ux[i][j][k]+Uy[i][j][k]*Uy[i][j][k]+Uz[i][j][k]*Uz[i][j][k];
						Eq0[i][j][k]=tNS[0]*(1+cuNS0[i][j][k]+0.5*(cuNS0[i][j][k]*cuNS0[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
						Eq1[i][j][k]=tNS[1]*(1+cuNS1[i][j][k]+0.5*(cuNS1[i][j][k]*cuNS1[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
						Eq2[i][j][k]=tNS[2]*(1+cuNS2[i][j][k]+0.5*(cuNS2[i][j][k]*cuNS2[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
						Eq3[i][j][k]=tNS[3]*(1+cuNS3[i][j][k]+0.5*(cuNS3[i][j][k]*cuNS3[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
						Eq4[i][j][k]=tNS[4]*(1+cuNS4[i][j][k]+0.5*(cuNS4[i][j][k]*cuNS4[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
						Eq5[i][j][k]=tNS[5]*(1+cuNS5[i][j][k]+0.5*(cuNS5[i][j][k]*cuNS5[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
						Eq6[i][j][k]=tNS[6]*(1+cuNS6[i][j][k]+0.5*(cuNS6[i][j][k]*cuNS6[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
						Eq7[i][j][k]=tNS[7]*(1+cuNS7[i][j][k]+0.5*(cuNS7[i][j][k]*cuNS7[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
						Eq8[i][j][k]=tNS[8]*(1+cuNS8[i][j][k]+0.5*(cuNS8[i][j][k]*cuNS8[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
						Eq9[i][j][k]=tNS[9]*(1+cuNS9[i][j][k]+0.5*(cuNS9[i][j][k]*cuNS9[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
						Eq10[i][j][k]=tNS[10]*(1+cuNS10[i][j][k]+0.5*(cuNS10[i][j][k]*cuNS10[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
						Eq11[i][j][k]=tNS[11]*(1+cuNS11[i][j][k]+0.5*(cuNS11[i][j][k]*cuNS11[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
						Eq12[i][j][k]=tNS[12]*(1+cuNS12[i][j][k]+0.5*(cuNS12[i][j][k]*cuNS12[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
						Eq13[i][j][k]=tNS[13]*(1+cuNS13[i][j][k]+0.5*(cuNS13[i][j][k]*cuNS13[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
						Eq14[i][j][k]=tNS[14]*(1+cuNS14[i][j][k]+0.5*(cuNS14[i][j][k]*cuNS14[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
						Eq15[i][j][k]=tNS[15]*(1+cuNS15[i][j][k]+0.5*(cuNS15[i][j][k]*cuNS15[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
						Eq16[i][j][k]=tNS[16]*(1+cuNS16[i][j][k]+0.5*(cuNS16[i][j][k]*cuNS16[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
						Eq17[i][j][k]=tNS[17]*(1+cuNS17[i][j][k]+0.5*(cuNS17[i][j][k]*cuNS17[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
						Eq18[i][j][k]=tNS[18]*(1+cuNS18[i][j][k]+0.5*(cuNS18[i][j][k]*cuNS18[i][j][k])-1.5*(ux2uy2uz2[i][j][k]));
					}
					if (domain[i][j][k]!=0){
						Eq0[i][j][k]=0.0;
						Eq1[i][j][k]=0.0;
						Eq2[i][j][k]=0.0;
						Eq3[i][j][k]=0.0;
						Eq4[i][j][k]=0.0;
						Eq5[i][j][k]=0.0;
						Eq6[i][j][k]=0.0;
						Eq7[i][j][k]=0.0;
						Eq8[i][j][k]=0.0;
						Eq9[i][j][k]=0.0;
						Eq10[i][j][k]=0.0;
						Eq11[i][j][k]=0.0;
						Eq12[i][j][k]=0.0;
						Eq13[i][j][k]=0.0;
						Eq14[i][j][k]=0.0;
						Eq15[i][j][k]=0.0;
						Eq16[i][j][k]=0.0;
						Eq17[i][j][k]=0.0;
						Eq18[i][j][k]=0.0;
					}
				}
			}
		}
	}
}

void equilibriumConcentrationCalculation(int Length, int Width, int Height, double ***Ux, double ***Uy, double ***Uz){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
					if (domain[i][j][k]==0){
						cuNS0[i][j][k]=4.0*(cxNS[0]*Ux[i][j][k]+cyNS[0]*Uy[i][j][k]+czNS[0]*Uz[i][j][k]);
						cuNS1[i][j][k]=4.0*(cxNS[1]*Ux[i][j][k]+cyNS[1]*Uy[i][j][k]+czNS[1]*Uz[i][j][k]);
						cuNS2[i][j][k]=4.0*(cxNS[2]*Ux[i][j][k]+cyNS[2]*Uy[i][j][k]+czNS[2]*Uz[i][j][k]);
						cuNS3[i][j][k]=4.0*(cxNS[3]*Ux[i][j][k]+cyNS[3]*Uy[i][j][k]+czNS[3]*Uz[i][j][k]);
						cuNS4[i][j][k]=4.0*(cxNS[4]*Ux[i][j][k]+cyNS[4]*Uy[i][j][k]+czNS[4]*Uz[i][j][k]);
						cuNS5[i][j][k]=4.0*(cxNS[5]*Ux[i][j][k]+cyNS[5]*Uy[i][j][k]+czNS[5]*Uz[i][j][k]);
						cuNS6[i][j][k]=4.0*(cxNS[6]*Ux[i][j][k]+cyNS[6]*Uy[i][j][k]+czNS[6]*Uz[i][j][k]);
						Eq0[i][j][k]=tNSD3Q7[0]*(1.0+cuNS0[i][j][k]);
						Eq1[i][j][k]=tNSD3Q7[1]*(1.0+cuNS1[i][j][k]);
						Eq2[i][j][k]=tNSD3Q7[2]*(1.0+cuNS2[i][j][k]);
						Eq3[i][j][k]=tNSD3Q7[3]*(1.0+cuNS3[i][j][k]);
						Eq4[i][j][k]=tNSD3Q7[4]*(1.0+cuNS4[i][j][k]);
						Eq5[i][j][k]=tNSD3Q7[5]*(1.0+cuNS5[i][j][k]);
						Eq6[i][j][k]=tNSD3Q7[6]*(1.0+cuNS6[i][j][k]);
					}
					if (domain[i][j][k]!=0){
						Eq0[i][j][k]=0.0;
						Eq1[i][j][k]=0.0;
						Eq2[i][j][k]=0.0;
						Eq3[i][j][k]=0.0;
						Eq4[i][j][k]=0.0;
						Eq5[i][j][k]=0.0;
						Eq6[i][j][k]=0.0;
					}
				}
			}
		}
	}
}

void MemoryFreeEquilibrium(int m, int n, int q){
	freememory(cuNS0,m,n,q);
	freememory(cuNS1,m,n,q);
	freememory(cuNS2,m,n,q);
	freememory(cuNS3,m,n,q);
	freememory(cuNS4,m,n,q);
	freememory(cuNS5,m,n,q);
	freememory(cuNS6,m,n,q);
	freememory(cuNS7,m,n,q);
	freememory(cuNS8,m,n,q);
	freememory(cuNS9,m,n,q);
	freememory(cuNS10,m,n,q);
	freememory(cuNS11,m,n,q);
	freememory(cuNS12,m,n,q);
	freememory(cuNS13,m,n,q);
	freememory(cuNS14,m,n,q);
	freememory(cuNS15,m,n,q);
	freememory(cuNS16,m,n,q);
	freememory(cuNS17,m,n,q);
	freememory(cuNS18,m,n,q);
	freememory(ux2uy2uz2,m,n,q);
	freememory(Eq0,m,n,q);
	freememory(Eq1,m,n,q);
	freememory(Eq2,m,n,q);
	freememory(Eq3,m,n,q);
	freememory(Eq4,m,n,q);
	freememory(Eq5,m,n,q);
	freememory(Eq6,m,n,q);
	freememory(Eq7,m,n,q);
	freememory(Eq8,m,n,q);
	freememory(Eq9,m,n,q);
	freememory(Eq10,m,n,q);
	freememory(Eq11,m,n,q);
	freememory(Eq12,m,n,q);
	freememory(Eq13,m,n,q);
	freememory(Eq14,m,n,q);
	freememory(Eq15,m,n,q);
	freememory(Eq16,m,n,q);
	freememory(Eq17,m,n,q);
	freememory(Eq18,m,n,q);
}