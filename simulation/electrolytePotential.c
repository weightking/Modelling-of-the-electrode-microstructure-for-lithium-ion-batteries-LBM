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
#include "electrodePotential.h"

double ***electrolyteIn0=NULL;
double ***electrolyteIn1=NULL;
double ***electrolyteIn2=NULL;
double ***electrolyteIn3=NULL;
double ***electrolyteIn4=NULL;
double ***electrolyteIn5=NULL;
double ***electrolyteIn6=NULL;
double ***electrolyteIn7=NULL;
double ***electrolyteIn8=NULL;
double ***electrolyteIn9=NULL;
double ***electrolyteIn10=NULL;
double ***electrolyteIn11=NULL;
double ***electrolyteIn12=NULL;
double ***electrolyteIn13=NULL;
double ***electrolyteIn14=NULL;
double ***electrolyteIn15=NULL;
double ***electrolyteIn16=NULL;
double ***electrolyteIn17=NULL;
double ***electrolyteIn18=NULL;
double ***electrolytePotential=NULL;
double ***electrolyteRelaxationTime=NULL;
double ***electrolyteAdditionalTerm=NULL;
double ***lastWallElectrolyte=NULL;

void FieldArrangeelectrolytePotential(int Length, int Width, int Height){
	electrolyteIn0=memoryarrange(Length,Width,Height);
	electrolyteIn1=memoryarrange(Length,Width,Height);
	electrolyteIn2=memoryarrange(Length,Width,Height);
	electrolyteIn3=memoryarrange(Length,Width,Height);
	electrolyteIn4=memoryarrange(Length,Width,Height);
	electrolyteIn5=memoryarrange(Length,Width,Height);
	electrolyteIn6=memoryarrange(Length,Width,Height);
//	electrolyteIn7=memoryarrange(Length,Width,Height);
//	electrolyteIn8=memoryarrange(Length,Width,Height);
//	electrolyteIn9=memoryarrange(Length,Width,Height);
//	electrolyteIn10=memoryarrange(Length,Width,Height);
//	electrolyteIn11=memoryarrange(Length,Width,Height);
//	electrolyteIn12=memoryarrange(Length,Width,Height);
//	electrolyteIn13=memoryarrange(Length,Width,Height);
//	electrolyteIn14=memoryarrange(Length,Width,Height);
//	electrolyteIn15=memoryarrange(Length,Width,Height);
//	electrolyteIn16=memoryarrange(Length,Width,Height);
//	electrolyteIn17=memoryarrange(Length,Width,Height);
//	electrolyteIn18=memoryarrange(Length,Width,Height);
	electrolytePotential=memoryarrange(Length,Width,Height); 
	electrolyteRelaxationTime=memoryarrange(Length,Width,Height);
	electrolyteAdditionalTerm=memoryarrange(Length,Width,Height);
	lastWallElectrolyte=memoryarrange(Length,Width,Height);
}

void FieldInitialelectrolytePotential(int Length, int Width, int Height, double electrolytePotentialInitial, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						electrolytePotential[i][j][k]=electrolytePotentialInitial;
						electrolyteIn0[i][j][k]=electrolytePotential[i][j][k]*tNS[0];
						electrolyteIn1[i][j][k]=electrolytePotential[i][j][k]*tNS[1];
						electrolyteIn2[i][j][k]=electrolytePotential[i][j][k]*tNS[2];
						electrolyteIn3[i][j][k]=electrolytePotential[i][j][k]*tNS[3];
						electrolyteIn4[i][j][k]=electrolytePotential[i][j][k]*tNS[4];
						electrolyteIn5[i][j][k]=electrolytePotential[i][j][k]*tNS[5];
						electrolyteIn6[i][j][k]=electrolytePotential[i][j][k]*tNS[6];
						electrolyteIn7[i][j][k]=electrolytePotential[i][j][k]*tNS[7];
						electrolyteIn8[i][j][k]=electrolytePotential[i][j][k]*tNS[8];
						electrolyteIn9[i][j][k]=electrolytePotential[i][j][k]*tNS[9];
						electrolyteIn10[i][j][k]=electrolytePotential[i][j][k]*tNS[10];
						electrolyteIn11[i][j][k]=electrolytePotential[i][j][k]*tNS[11];
						electrolyteIn12[i][j][k]=electrolytePotential[i][j][k]*tNS[12];
						electrolyteIn13[i][j][k]=electrolytePotential[i][j][k]*tNS[13];
						electrolyteIn14[i][j][k]=electrolytePotential[i][j][k]*tNS[14];
						electrolyteIn15[i][j][k]=electrolytePotential[i][j][k]*tNS[15];
						electrolyteIn16[i][j][k]=electrolytePotential[i][j][k]*tNS[16];
						electrolyteIn17[i][j][k]=electrolytePotential[i][j][k]*tNS[17];
						electrolyteIn18[i][j][k]=electrolytePotential[i][j][k]*tNS[18];
					}
					if (domain[i][j][k]!=0||rho[i][j][k]<=gasCritical||domainConnection[i][j][k]!=1){
						electrolytePotential[i][j][k]=0.0;
						electrolyteIn0[i][j][k]=0.0;
						electrolyteIn1[i][j][k]=0.0;
						electrolyteIn2[i][j][k]=0.0;
						electrolyteIn3[i][j][k]=0.0;
						electrolyteIn4[i][j][k]=0.0;
						electrolyteIn5[i][j][k]=0.0;
						electrolyteIn6[i][j][k]=0.0;
						electrolyteIn7[i][j][k]=0.0;
						electrolyteIn8[i][j][k]=0.0;
						electrolyteIn9[i][j][k]=0.0;
						electrolyteIn10[i][j][k]=0.0;
						electrolyteIn11[i][j][k]=0.0;
						electrolyteIn12[i][j][k]=0.0;
						electrolyteIn13[i][j][k]=0.0;
						electrolyteIn14[i][j][k]=0.0;
						electrolyteIn15[i][j][k]=0.0;
						electrolyteIn16[i][j][k]=0.0;
						electrolyteIn17[i][j][k]=0.0;
						electrolyteIn18[i][j][k]=0.0;
					}
				}
			}
		}
	}
}

void FieldInitialelectrolytePotentialD3Q7(int Length, int Width, int Initial,int Height, double electrolytePotentialInitial, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=Initial;k<Height;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						electrolytePotential[i][j][k]=electrolytePotentialInitial;
						electrolyteIn0[i][j][k]=electrolytePotential[i][j][k]*tNSD3Q7[0];
						electrolyteIn1[i][j][k]=electrolytePotential[i][j][k]*tNSD3Q7[1];
						electrolyteIn2[i][j][k]=electrolytePotential[i][j][k]*tNSD3Q7[2];
						electrolyteIn3[i][j][k]=electrolytePotential[i][j][k]*tNSD3Q7[3];
						electrolyteIn4[i][j][k]=electrolytePotential[i][j][k]*tNSD3Q7[4];
						electrolyteIn5[i][j][k]=electrolytePotential[i][j][k]*tNSD3Q7[5];
						electrolyteIn6[i][j][k]=electrolytePotential[i][j][k]*tNSD3Q7[6];
					}
					if (domain[i][j][k]!=0||rho[i][j][k]<=gasCritical||domainConnection[i][j][k]!=1){
						electrolytePotential[i][j][k]=0.0;
						electrolyteIn0[i][j][k]=0.0;
						electrolyteIn1[i][j][k]=0.0;
						electrolyteIn2[i][j][k]=0.0;
						electrolyteIn3[i][j][k]=0.0;
						electrolyteIn4[i][j][k]=0.0;
						electrolyteIn5[i][j][k]=0.0;
						electrolyteIn6[i][j][k]=0.0;
					}
				}
			}
		}
	}
}

void electrolyteRelaxationTimeCalculation(int Length, int Width, int Height, double constant, double constantMol, double constantLelectrolyteConductivity, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						electrolyteRelaxationTime[i][j][k] = 1.0 / ((1.12*(-8.2488 + 0.053248*constant - 0.000029871*constant*constant + 0.26235*Concentration[i][j][k] * constantMol - 0.0093063*Concentration[i][j][k] * constantMol * constant + 0.000008069*Concentration[i][j][k] * constantMol * constant*constant + 0.22002*Concentration[i][j][k] * constantMol * Concentration[i][j][k] * constantMol - 0.0001765*Concentration[i][j][k] * constantMol * Concentration[i][j][k] * constantMol * constant))* (1.12*(-8.2488 + 0.053248*constant - 0.000029871*constant*constant + 0.26235*Concentration[i][j][k] * constantMol - 0.0093063*Concentration[i][j][k] * constantMol * constant + 0.000008069*Concentration[i][j][k] * constantMol * constant*constant + 0.22002*Concentration[i][j][k] * constantMol * Concentration[i][j][k] * constantMol - 0.0001765*Concentration[i][j][k] * constantMol * Concentration[i][j][k] * constantMol * constant)) * Concentration[i][j][k] * constantMol*constantLelectrolyteConductivity * 4.0 + 0.5);
//						electrolyteRelaxationTime[i][j][k] = 1.0 / 1.43333;
					}
				}
			}
		}
	}
}

void additionalTermCalculation(int Length,int Width,int Height,double constant,double gasCritical,double concentrationDiffusivity,double concentrationFirstDiffusivity,double concentrationSecondDiffusivity,double concentrationThirdDiffusivity,int charge,int charge1,int charge2,int charge3){
	int i;
	int j;
	int k;
	double concentrationNeightbour1;
	double concentrationFirstNeightbour1;
	double concentrationSecondNeightbour1;
	double concentrationNeightbour2;
	double concentrationFirstNeightbour2;
	double concentrationSecondNeightbour2;
	double concentrationNeightbour3;
	double concentrationFirstNeightbour3;
	double concentrationSecondNeightbour3;
	double concentrationNeightbour4;
	double concentrationFirstNeightbour4;
	double concentrationSecondNeightbour4;
	double concentrationNeightbour5;
	double concentrationFirstNeightbour5;
	double concentrationSecondNeightbour5;
	double concentrationNeightbour6;
	double concentrationFirstNeightbour6;
	double concentrationSecondNeightbour6;
	double concentrationNeightbour7;
	double concentrationFirstNeightbour7;
	double concentrationSecondNeightbour7;
	double concentrationNeightbour8;
	double concentrationFirstNeightbour8;
	double concentrationSecondNeightbour8;
	double concentrationNeightbour9;
	double concentrationFirstNeightbour9;
	double concentrationSecondNeightbour9;
	double concentrationNeightbour10;
	double concentrationFirstNeightbour10;
	double concentrationSecondNeightbour10;
	double concentrationNeightbour11;
	double concentrationFirstNeightbour11;
	double concentrationSecondNeightbour11;
	double concentrationNeightbour12;
	double concentrationFirstNeightbour12;
	double concentrationSecondNeightbour12;
	double concentrationNeightbour13;
	double concentrationFirstNeightbour13;
	double concentrationSecondNeightbour13;
	double concentrationNeightbour14;
	double concentrationFirstNeightbour14;
	double concentrationSecondNeightbour14;
	double concentrationNeightbour15;
	double concentrationFirstNeightbour15;
	double concentrationSecondNeightbour15;
	double concentrationNeightbour16;
	double concentrationFirstNeightbour16;
	double concentrationSecondNeightbour16;
	double concentrationNeightbour17;
	double concentrationFirstNeightbour17;
	double concentrationSecondNeightbour17;
	double concentrationNeightbour18;
	double concentrationFirstNeightbour18;
	double concentrationSecondNeightbour18;
	double concentrationAdditionalTerm;
	double concentrationFirstAdditionalTerm;
	double concentrationSecondAdditionalTerm;
	double concentrationThirdAdditionalTerm;
#pragma omp parallel private(i,j,k,concentrationNeightbour1,concentrationFirstNeightbour1,concentrationSecondNeightbour1,concentrationNeightbour2,concentrationFirstNeightbour2,concentrationSecondNeightbour2,concentrationNeightbour3,concentrationFirstNeightbour3,concentrationSecondNeightbour3,concentrationNeightbour4,concentrationFirstNeightbour4,concentrationSecondNeightbour4,concentrationNeightbour5,concentrationFirstNeightbour5,concentrationSecondNeightbour5,concentrationNeightbour6,concentrationFirstNeightbour6,concentrationSecondNeightbour6,concentrationNeightbour7,concentrationFirstNeightbour7,concentrationSecondNeightbour7,concentrationNeightbour8,concentrationFirstNeightbour8,concentrationSecondNeightbour8,concentrationNeightbour9,concentrationFirstNeightbour9,concentrationSecondNeightbour9,concentrationNeightbour10,concentrationFirstNeightbour10,concentrationSecondNeightbour10,concentrationNeightbour11,concentrationFirstNeightbour11,concentrationSecondNeightbour11,concentrationNeightbour12,concentrationFirstNeightbour12,concentrationSecondNeightbour12,concentrationNeightbour13,concentrationFirstNeightbour13,concentrationSecondNeightbour13,concentrationNeightbour14,concentrationFirstNeightbour14,concentrationSecondNeightbour14,concentrationNeightbour15,concentrationFirstNeightbour15,concentrationSecondNeightbour15,concentrationNeightbour16,concentrationFirstNeightbour16,concentrationSecondNeightbour16,concentrationNeightbour17,concentrationFirstNeightbour17,concentrationSecondNeightbour17,concentrationNeightbour18,concentrationFirstNeightbour18,concentrationSecondNeightbour18,concentrationAdditionalTerm,concentrationFirstAdditionalTerm,concentrationSecondAdditionalTerm,concentrationThirdAdditionalTerm)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						concentrationNeightbour1=neightbour1(i,j,k,Length,Concentration);
						concentrationFirstNeightbour1=neightbour1(i,j,k,Length,ConcentrationFirst);
						concentrationSecondNeightbour1=neightbour1(i,j,k,Length,ConcentrationSecond);

						concentrationNeightbour2=neightbour2(i,j,k,Width,Concentration);
						concentrationFirstNeightbour2=neightbour2(i,j,k,Width,ConcentrationFirst);
						concentrationSecondNeightbour2=neightbour2(i,j,k,Width,ConcentrationSecond);

						concentrationNeightbour3=neightbour3(i,j,k,Height,Concentration);
						concentrationFirstNeightbour3=neightbour3(i,j,k,Height,ConcentrationFirst);
						concentrationSecondNeightbour3=neightbour3(i,j,k,Height,ConcentrationSecond);

						concentrationNeightbour4=neightbour4(i,j,k,Length,Concentration);
						concentrationFirstNeightbour4=neightbour4(i,j,k,Length,ConcentrationFirst);
						concentrationSecondNeightbour4=neightbour4(i,j,k,Length,ConcentrationSecond);

						concentrationNeightbour5=neightbour5(i,j,k,Width,Concentration);
						concentrationFirstNeightbour5=neightbour5(i,j,k,Width,ConcentrationFirst);
						concentrationSecondNeightbour5=neightbour5(i,j,k,Width,ConcentrationSecond);

						concentrationNeightbour6=neightbour6(i,j,k,Height,Concentration);
						concentrationFirstNeightbour6=neightbour6(i,j,k,Height,ConcentrationFirst);
						concentrationSecondNeightbour6=neightbour6(i,j,k,Height,ConcentrationSecond);

						concentrationNeightbour7=neightbour7(i,j,k,Length,Width,Concentration);
						concentrationFirstNeightbour7=neightbour7(i,j,k,Length,Width,ConcentrationFirst);
						concentrationSecondNeightbour7=neightbour7(i,j,k,Length,Width,ConcentrationSecond);

						concentrationNeightbour8=neightbour8(i,j,k,Length,Height,Concentration);
						concentrationFirstNeightbour8=neightbour8(i,j,k,Length,Height,ConcentrationFirst);
						concentrationSecondNeightbour8=neightbour8(i,j,k,Length,Height,ConcentrationSecond);

						concentrationNeightbour9=neightbour9(i,j,k,Width,Height,Concentration);
						concentrationFirstNeightbour9=neightbour9(i,j,k,Width,Height,ConcentrationFirst);
						concentrationSecondNeightbour9=neightbour9(i,j,k,Width,Height,ConcentrationSecond);

						concentrationNeightbour10=neightbour10(i,j,k,Length,Width,Concentration);
						concentrationFirstNeightbour10=neightbour10(i,j,k,Length,Width,ConcentrationFirst);
						concentrationSecondNeightbour10=neightbour10(i,j,k,Length,Width,ConcentrationSecond);

						concentrationNeightbour11=neightbour11(i,j,k,Length,Height,Concentration);
						concentrationFirstNeightbour11=neightbour11(i,j,k,Length,Height,ConcentrationFirst);
						concentrationSecondNeightbour11=neightbour11(i,j,k,Length,Height,ConcentrationSecond);

						concentrationNeightbour12=neightbour12(i,j,k,Width,Height,Concentration);
						concentrationFirstNeightbour12=neightbour12(i,j,k,Width,Height,ConcentrationFirst);
						concentrationSecondNeightbour12=neightbour12(i,j,k,Width,Height,ConcentrationSecond);

						concentrationNeightbour13=neightbour13(i,j,k,Length,Width,Concentration);
						concentrationFirstNeightbour13=neightbour13(i,j,k,Length,Width,ConcentrationFirst);
						concentrationSecondNeightbour13=neightbour13(i,j,k,Length,Width,ConcentrationSecond);

						concentrationNeightbour14=neightbour14(i,j,k,Length,Height,Concentration);
						concentrationFirstNeightbour14=neightbour14(i,j,k,Length,Height,ConcentrationFirst);
						concentrationSecondNeightbour14=neightbour14(i,j,k,Length,Height,ConcentrationSecond);

						concentrationNeightbour15=neightbour15(i,j,k,Width,Height,Concentration);
						concentrationFirstNeightbour15=neightbour15(i,j,k,Width,Height,ConcentrationFirst);
						concentrationSecondNeightbour15=neightbour15(i,j,k,Width,Height,ConcentrationSecond);

						concentrationNeightbour16=neightbour16(i,j,k,Length,Width,Concentration);
						concentrationFirstNeightbour16=neightbour16(i,j,k,Length,Width,ConcentrationFirst);
						concentrationSecondNeightbour16=neightbour16(i,j,k,Length,Width,ConcentrationSecond);

						concentrationNeightbour17=neightbour17(i,j,k,Length,Height,Concentration);
						concentrationFirstNeightbour17=neightbour17(i,j,k,Length,Height,ConcentrationFirst);
						concentrationSecondNeightbour17=neightbour17(i,j,k,Length,Height,ConcentrationSecond);

						concentrationNeightbour18=neightbour18(i,j,k,Width,Height,Concentration);
						concentrationFirstNeightbour18=neightbour18(i,j,k,Width,Height,ConcentrationFirst);
						concentrationSecondNeightbour18=neightbour18(i,j,k,Width,Height,ConcentrationSecond);

						if (Domain1[i][j][k]!=0||rho1[i][j][k]<=gasCritical){
							concentrationNeightbour4=concentrationNeightbour1;
							concentrationFirstNeightbour4=concentrationFirstNeightbour1;
							concentrationSecondNeightbour4=concentrationSecondNeightbour1;
						}
						if (Domain4[i][j][k]!=0||rho4[i][j][k]<=gasCritical){
							concentrationNeightbour1=concentrationNeightbour4;
							concentrationFirstNeightbour1=concentrationFirstNeightbour4;
							concentrationSecondNeightbour1=concentrationSecondNeightbour4;
						}
						if ((Domain1[i][j][k]!=0||rho1[i][j][k]<=gasCritical)&&(Domain4[i][j][k]!=0||rho4[i][j][k]<=gasCritical)){
							concentrationNeightbour4=Concentration[i][j][k];
							concentrationFirstNeightbour4=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour4=ConcentrationSecond[i][j][k];
							concentrationNeightbour1=Concentration[i][j][k];
							concentrationFirstNeightbour1=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour1=ConcentrationSecond[i][j][k];							
						}

						if (Domain2[i][j][k]!=0||rho2[i][j][k]<=gasCritical){
							concentrationNeightbour5=concentrationNeightbour2;
							concentrationFirstNeightbour5=concentrationFirstNeightbour2;
							concentrationSecondNeightbour5=concentrationSecondNeightbour2;
						}
						if (Domain5[i][j][k]!=0||rho5[i][j][k]<=gasCritical){
							concentrationNeightbour2=concentrationNeightbour5;
							concentrationFirstNeightbour2=concentrationFirstNeightbour5;
							concentrationSecondNeightbour2=concentrationSecondNeightbour5;
						}
						if ((Domain2[i][j][k]!=0||rho2[i][j][k]<=gasCritical)&&(Domain5[i][j][k]!=0||rho5[i][j][k]<=gasCritical)){
							concentrationNeightbour5=Concentration[i][j][k];
							concentrationFirstNeightbour5=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour5=ConcentrationSecond[i][j][k];
							concentrationNeightbour2=Concentration[i][j][k];
							concentrationFirstNeightbour2=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour2=ConcentrationSecond[i][j][k];							
						}

						if (Domain3[i][j][k]!=0||rho3[i][j][k]<=gasCritical){
							concentrationNeightbour6=concentrationNeightbour3;
							concentrationFirstNeightbour6=concentrationFirstNeightbour3;
							concentrationSecondNeightbour6=concentrationSecondNeightbour3;
						}
						if (Domain6[i][j][k]!=0||rho6[i][j][k]<=gasCritical){
							concentrationNeightbour3=concentrationNeightbour6;
							concentrationFirstNeightbour3=concentrationFirstNeightbour6;
							concentrationSecondNeightbour3=concentrationSecondNeightbour6;
						}
						if ((Domain3[i][j][k]!=0||rho3[i][j][k]<=gasCritical)&&(Domain6[i][j][k]!=0||rho6[i][j][k]<=gasCritical)){
							concentrationNeightbour6=Concentration[i][j][k];
							concentrationFirstNeightbour6=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour6=ConcentrationSecond[i][j][k];
							concentrationNeightbour3=Concentration[i][j][k];
							concentrationFirstNeightbour3=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour3=ConcentrationSecond[i][j][k];							
						}

						if (Domain7[i][j][k]!=0||rho7[i][j][k]<=gasCritical){
							concentrationNeightbour10=concentrationNeightbour7;
							concentrationFirstNeightbour10=concentrationFirstNeightbour7;
							concentrationSecondNeightbour10=concentrationSecondNeightbour7;
						}
						if (Domain10[i][j][k]!=0||rho10[i][j][k]<=gasCritical){
							concentrationNeightbour7=concentrationNeightbour10;
							concentrationFirstNeightbour7=concentrationFirstNeightbour10;
							concentrationSecondNeightbour7=concentrationSecondNeightbour10;
						}
						if ((Domain7[i][j][k]!=0||rho7[i][j][k]<=gasCritical)&&(Domain10[i][j][k]!=0||rho10[i][j][k]<=gasCritical)){
							concentrationNeightbour10=Concentration[i][j][k];
							concentrationFirstNeightbour10=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour10=ConcentrationSecond[i][j][k];
							concentrationNeightbour7=Concentration[i][j][k];
							concentrationFirstNeightbour7=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour7=ConcentrationSecond[i][j][k];							
						}

						if (Domain8[i][j][k]!=0||rho8[i][j][k]<=gasCritical){
							concentrationNeightbour11=concentrationNeightbour8;
							concentrationFirstNeightbour11=concentrationFirstNeightbour8;
							concentrationSecondNeightbour11=concentrationSecondNeightbour8;
						}
						if (Domain11[i][j][k]!=0||rho11[i][j][k]<=gasCritical){
							concentrationNeightbour8=concentrationNeightbour11;
							concentrationFirstNeightbour8=concentrationFirstNeightbour11;
							concentrationSecondNeightbour8=concentrationSecondNeightbour11;
						}
						if ((Domain8[i][j][k]!=0||rho8[i][j][k]<=gasCritical)&&(Domain11[i][j][k]!=0||rho11[i][j][k]<=gasCritical)){
							concentrationNeightbour11=Concentration[i][j][k];
							concentrationFirstNeightbour11=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour11=ConcentrationSecond[i][j][k];
							concentrationNeightbour8=Concentration[i][j][k];
							concentrationFirstNeightbour8=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour8=ConcentrationSecond[i][j][k];							
						}

						if (Domain9[i][j][k]!=0||rho9[i][j][k]<=gasCritical){
							concentrationNeightbour12=concentrationNeightbour9;
							concentrationFirstNeightbour12=concentrationFirstNeightbour9;
							concentrationSecondNeightbour12=concentrationSecondNeightbour9;
						}
						if (Domain12[i][j][k]!=0||rho12[i][j][k]<=gasCritical){
							concentrationNeightbour9=concentrationNeightbour12;
							concentrationFirstNeightbour9=concentrationFirstNeightbour12;
							concentrationSecondNeightbour9=concentrationSecondNeightbour12;
						}
						if ((Domain9[i][j][k]!=0||rho9[i][j][k]<=gasCritical)&&(Domain12[i][j][k]!=0||rho12[i][j][k]<=gasCritical)){
							concentrationNeightbour12=Concentration[i][j][k];
							concentrationFirstNeightbour12=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour12=ConcentrationSecond[i][j][k];
							concentrationNeightbour9=Concentration[i][j][k];
							concentrationFirstNeightbour9=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour9=ConcentrationSecond[i][j][k];							
						}

						if (Domain13[i][j][k]!=0||rho13[i][j][k]<=gasCritical){
							concentrationNeightbour16=concentrationNeightbour13;
							concentrationFirstNeightbour16=concentrationFirstNeightbour13;
							concentrationSecondNeightbour16=concentrationSecondNeightbour13;
						}
						if (Domain16[i][j][k]!=0||rho16[i][j][k]<=gasCritical){
							concentrationNeightbour13=concentrationNeightbour16;
							concentrationFirstNeightbour13=concentrationFirstNeightbour16;
							concentrationSecondNeightbour13=concentrationSecondNeightbour16;
						}
						if ((Domain13[i][j][k]!=0||rho13[i][j][k]<=gasCritical)&&(Domain16[i][j][k]!=0||rho16[i][j][k]<=gasCritical)){
							concentrationNeightbour16=Concentration[i][j][k];
							concentrationFirstNeightbour16=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour16=ConcentrationSecond[i][j][k];
							concentrationNeightbour13=Concentration[i][j][k];
							concentrationFirstNeightbour13=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour13=ConcentrationSecond[i][j][k];							
						}

						if (Domain14[i][j][k]!=0||rho14[i][j][k]<=gasCritical){
							concentrationNeightbour17=concentrationNeightbour14;
							concentrationFirstNeightbour17=concentrationFirstNeightbour14;
							concentrationSecondNeightbour17=concentrationSecondNeightbour14;
						}
						if (Domain17[i][j][k]!=0||rho17[i][j][k]<=gasCritical){
							concentrationNeightbour14=concentrationNeightbour17;
							concentrationFirstNeightbour14=concentrationFirstNeightbour17;
							concentrationSecondNeightbour14=concentrationSecondNeightbour17;
						}
						if ((Domain14[i][j][k]!=0||rho14[i][j][k]<=gasCritical)&&(Domain17[i][j][k]!=0||rho17[i][j][k]<=gasCritical)){
							concentrationNeightbour17=Concentration[i][j][k];
							concentrationFirstNeightbour17=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour17=ConcentrationSecond[i][j][k];
							concentrationNeightbour14=Concentration[i][j][k];
							concentrationFirstNeightbour14=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour14=ConcentrationSecond[i][j][k];							
						}

						if (Domain15[i][j][k]!=0||rho15[i][j][k]<=gasCritical){
							concentrationNeightbour18=concentrationNeightbour15;
							concentrationFirstNeightbour18=concentrationFirstNeightbour15;
							concentrationSecondNeightbour18=concentrationSecondNeightbour15;
						}
						if (Domain18[i][j][k]!=0||rho18[i][j][k]<=gasCritical){
							concentrationNeightbour15=concentrationNeightbour18;
							concentrationFirstNeightbour15=concentrationFirstNeightbour18;
							concentrationSecondNeightbour15=concentrationSecondNeightbour18;
						}
						if ((Domain15[i][j][k]!=0||rho15[i][j][k]<=gasCritical)&&(Domain18[i][j][k]!=0||rho18[i][j][k]<=gasCritical)){
							concentrationNeightbour18=Concentration[i][j][k];
							concentrationFirstNeightbour18=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour18=ConcentrationSecond[i][j][k];
							concentrationNeightbour15=Concentration[i][j][k];
							concentrationFirstNeightbour15=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour15=ConcentrationSecond[i][j][k];							
						}

						if (k==0){
							concentrationNeightbour6=Concentration[i][j][k];
							concentrationFirstNeightbour6=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour6=ConcentrationSecond[i][j][k];
							concentrationNeightbour11=concentrationNeightbour4;
							concentrationFirstNeightbour11=concentrationFirstNeightbour4;
							concentrationSecondNeightbour11=concentrationSecondNeightbour4;
							concentrationNeightbour12=concentrationNeightbour5;
							concentrationFirstNeightbour12=concentrationFirstNeightbour5;
							concentrationSecondNeightbour12=concentrationSecondNeightbour5;
							concentrationNeightbour14=concentrationNeightbour1;
							concentrationFirstNeightbour14=concentrationFirstNeightbour1;
							concentrationSecondNeightbour14=concentrationSecondNeightbour1;
							concentrationNeightbour15=concentrationNeightbour2;
							concentrationFirstNeightbour15=concentrationFirstNeightbour2;
							concentrationSecondNeightbour15=concentrationSecondNeightbour2;
						}
						if (k==Height-1){
							concentrationNeightbour3=Concentration[i][j][k];
							concentrationFirstNeightbour3=ConcentrationFirst[i][j][k];
							concentrationSecondNeightbour3=ConcentrationSecond[i][j][k];
							concentrationNeightbour8=concentrationNeightbour1;
							concentrationFirstNeightbour8=concentrationFirstNeightbour1;
							concentrationSecondNeightbour8=concentrationSecondNeightbour1;
							concentrationNeightbour9=concentrationNeightbour2;
							concentrationFirstNeightbour9=concentrationFirstNeightbour2;
							concentrationSecondNeightbour9=concentrationSecondNeightbour2;
							concentrationNeightbour17=concentrationNeightbour4;
							concentrationFirstNeightbour17=concentrationFirstNeightbour4;
							concentrationSecondNeightbour17=concentrationSecondNeightbour4;
							concentrationNeightbour18=concentrationNeightbour5;
							concentrationFirstNeightbour18=concentrationFirstNeightbour5;
							concentrationSecondNeightbour18=concentrationSecondNeightbour5;
						}

						concentrationAdditionalTerm=(concentrationNeightbour7+concentrationNeightbour10+concentrationNeightbour13+concentrationNeightbour16+concentrationNeightbour8+concentrationNeightbour11+concentrationNeightbour14+concentrationNeightbour17+concentrationNeightbour9+concentrationNeightbour12+concentrationNeightbour15+concentrationNeightbour18+2.0*concentrationNeightbour1+2.0*concentrationNeightbour4+2.0*concentrationNeightbour2+2.0*concentrationNeightbour5+2.0*concentrationNeightbour3+2.0*concentrationNeightbour6-24.0*Concentration[i][j][k])/6.0;
						concentrationFirstAdditionalTerm=(concentrationFirstNeightbour7+concentrationFirstNeightbour10+concentrationFirstNeightbour13+concentrationFirstNeightbour16+concentrationFirstNeightbour8+concentrationFirstNeightbour11+concentrationFirstNeightbour14+concentrationFirstNeightbour17+concentrationFirstNeightbour9+concentrationFirstNeightbour12+concentrationFirstNeightbour15+concentrationFirstNeightbour18+2.0*concentrationFirstNeightbour1+2.0*concentrationFirstNeightbour4+2.0*concentrationFirstNeightbour2+2.0*concentrationFirstNeightbour5+2.0*concentrationFirstNeightbour3+2.0*concentrationFirstNeightbour6-24.0*ConcentrationFirst[i][j][k])/6.0;
						concentrationSecondAdditionalTerm=(concentrationSecondNeightbour7+concentrationSecondNeightbour10+concentrationSecondNeightbour13+concentrationSecondNeightbour16+concentrationSecondNeightbour8+concentrationSecondNeightbour11+concentrationSecondNeightbour14+concentrationSecondNeightbour17+concentrationSecondNeightbour9+concentrationSecondNeightbour12+concentrationSecondNeightbour15+concentrationSecondNeightbour18+2.0*concentrationSecondNeightbour1+2.0*concentrationSecondNeightbour4+2.0*concentrationSecondNeightbour2+2.0*concentrationSecondNeightbour5+2.0*concentrationSecondNeightbour3+2.0*concentrationSecondNeightbour6-24.0*ConcentrationSecond[i][j][k])/6.0;
						concentrationThirdAdditionalTerm=(charge*concentrationAdditionalTerm+charge1*concentrationFirstAdditionalTerm+charge2*concentrationSecondAdditionalTerm)/(-charge3);
						electrolyteAdditionalTerm[i][j][k]=constant*(charge*concentrationDiffusivity*concentrationAdditionalTerm+charge1*concentrationFirstDiffusivity*concentrationFirstAdditionalTerm+charge2*concentrationSecondDiffusivity*concentrationSecondAdditionalTerm+charge3*concentrationThirdDiffusivity*concentrationThirdAdditionalTerm);
					}
				}
			}
		}
	}
}

void additionalTermCalculationD3Q7(int Length,int Width,int Height,double constant1,double constant2,double K,double gasCritical,double constantMol,double constantLelectrolyteConductivity,double exchange){
	int i;
	int j;
	int k;
	double concentrationNeightbour1;
	double concentrationNeightbour2;
	double concentrationNeightbour3;
	double concentrationNeightbour4;
	double concentrationNeightbour5;
	double concentrationNeightbour6;
	double concentrationAdditionalTermPartialFirstX;
	double concentrationAdditionalTermPartialFirstY;
	double concentrationAdditionalTermPartialFirstZ;
	double concentrationAdditionalTermPartialSecond;
	double otherAdditionalTermPart1Part;
	double otherAdditionalTermPart1;
	double otherAdditionalTermPart2;
	double otherAdditionalTerm;
	double Part1PartialFirstX;
	double Part1PartialFirstY;
	double Part1PartialFirstZ;
	double Part2PartialFirstX;
	double Part2PartialFirstY;
	double Part2PartialFirstZ;
	double otherAdditionalTermPartialFirstX;
	double otherAdditionalTermPartialFirstY;
	double otherAdditionalTermPartialFirstZ;
#pragma omp parallel private(i,j,k,concentrationNeightbour1,concentrationNeightbour2,concentrationNeightbour3,concentrationNeightbour4,concentrationNeightbour5,concentrationNeightbour6,concentrationAdditionalTermPartialFirstX,concentrationAdditionalTermPartialFirstY,concentrationAdditionalTermPartialFirstZ,concentrationAdditionalTermPartialSecond,otherAdditionalTermPart1Part,otherAdditionalTermPart1,otherAdditionalTermPart2,otherAdditionalTerm,Part1PartialFirstX,Part1PartialFirstY,Part1PartialFirstZ,Part2PartialFirstX,Part2PartialFirstY,Part2PartialFirstZ,otherAdditionalTermPartialFirstX,otherAdditionalTermPartialFirstY,otherAdditionalTermPartialFirstZ)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						concentrationNeightbour1=neightbour1(i,j,k,Length,Concentration)*constantMol;
						concentrationNeightbour2=neightbour2(i,j,k,Width,Concentration)*constantMol;
						concentrationNeightbour3=neightbour3(i,j,k,Height,Concentration)*constantMol;
						concentrationNeightbour4=neightbour4(i,j,k,Length,Concentration)*constantMol;
						concentrationNeightbour5=neightbour5(i,j,k,Width,Concentration)*constantMol;
						concentrationNeightbour6=neightbour6(i,j,k,Height,Concentration)*constantMol;

						if (Domain1[i][j][k]!=0||rho1[i][j][k]<=gasCritical){
							concentrationNeightbour4=Concentration[i][j][k]*constantMol;
						}
						if (Domain4[i][j][k]!=0||rho4[i][j][k]<=gasCritical){
							concentrationNeightbour1=Concentration[i][j][k]*constantMol;
						}
						if ((Domain1[i][j][k]!=0||rho1[i][j][k]<=gasCritical)&&(Domain4[i][j][k]!=0||rho4[i][j][k]<=gasCritical)){
							concentrationNeightbour4=Concentration[i][j][k]*constantMol;
							concentrationNeightbour1=Concentration[i][j][k]*constantMol;				
						}

						if (Domain2[i][j][k]!=0||rho2[i][j][k]<=gasCritical){
							concentrationNeightbour5=Concentration[i][j][k]*constantMol;
						}
						if (Domain5[i][j][k]!=0||rho5[i][j][k]<=gasCritical){
							concentrationNeightbour2=Concentration[i][j][k]*constantMol;
						}
						if ((Domain2[i][j][k]!=0||rho2[i][j][k]<=gasCritical)&&(Domain5[i][j][k]!=0||rho5[i][j][k]<=gasCritical)){
							concentrationNeightbour5=Concentration[i][j][k]*constantMol;
							concentrationNeightbour2=Concentration[i][j][k]*constantMol;							
						}

						if (Domain3[i][j][k]!=0||rho3[i][j][k]<=gasCritical){
							concentrationNeightbour6=Concentration[i][j][k]*constantMol;
						}
						if (Domain6[i][j][k]!=0||rho6[i][j][k]<=gasCritical){
							concentrationNeightbour3=Concentration[i][j][k]*constantMol;
						}
						if ((Domain3[i][j][k]!=0||rho3[i][j][k]<=gasCritical)&&(Domain6[i][j][k]!=0||rho6[i][j][k]<=gasCritical)){
							concentrationNeightbour6=Concentration[i][j][k]*constantMol;
							concentrationNeightbour3=Concentration[i][j][k]*constantMol;							
						}

						if (j==0){
							concentrationNeightbour5=Concentration[i][j][k]*constantMol;							
						}
						concentrationAdditionalTermPartialSecond = (2.0*concentrationNeightbour1 + 2.0*concentrationNeightbour4 + 2.0*concentrationNeightbour2 + 2.0*concentrationNeightbour5 + 2.0*concentrationNeightbour3 + 2.0*concentrationNeightbour6 - 12.0*Concentration[i][j][k]*constantMol) / 2.0;
						concentrationAdditionalTermPartialFirstX = (concentrationNeightbour1-concentrationNeightbour4)/2.0;
                        concentrationAdditionalTermPartialFirstY = (concentrationNeightbour2-concentrationNeightbour5)/2.0;
						concentrationAdditionalTermPartialFirstZ = (concentrationNeightbour3-concentrationNeightbour6)/2.0;
						otherAdditionalTermPart1Part = -8.2488+0.053248*K-0.000029871*K*K+0.26235*Concentration[i][j][k]*constantMol/1000.0-0.0093063*Concentration[i][j][k]*constantMol/1000.0*K+0.000008069*Concentration[i][j][k]*constantMol/1000.0*K*K+0.22002*Concentration[i][j][k]*constantMol/1000.0*Concentration[i][j][k]*constantMol/1000.0-0.0001765*Concentration[i][j][k]*constantMol/1000.0*Concentration[i][j][k]*constantMol/1000.0*K;
						otherAdditionalTermPart1 = constantLelectrolyteConductivity*exchange*1.2544*2*constant1*K/constant2*otherAdditionalTermPart1Part*otherAdditionalTermPart1Part/1000.0;
						otherAdditionalTermPart2 = 0.637/0.601*(0.601-0.24*sqrt(Concentration[i][j][k]*constantMol/1000.0)+0.982*(1-0.0052*(K-294))*Concentration[i][j][k]*constantMol/1000.0*sqrt(Concentration[i][j][k]*constantMol/1000.0));
						otherAdditionalTerm = otherAdditionalTermPart1*otherAdditionalTermPart2;
						Part1PartialFirstX = constantLelectrolyteConductivity*exchange*1.2544*2*constant1*K/constant2*2*otherAdditionalTermPart1Part*(0.26235*concentrationAdditionalTermPartialFirstX/1000.0-0.0093063*K*concentrationAdditionalTermPartialFirstX/1000.0+0.000008069*K*K*concentrationAdditionalTermPartialFirstX/1000.0+0.44004*Concentration[i][j][k]*constantMol/1000.0*concentrationAdditionalTermPartialFirstX/1000.0-0.0001765*2*K*Concentration[i][j][k]*constantMol/1000.0*concentrationAdditionalTermPartialFirstX/1000.0)/1000.0;
						Part1PartialFirstY = constantLelectrolyteConductivity*exchange*1.2544*2*constant1*K/constant2*2*otherAdditionalTermPart1Part*(0.26235*concentrationAdditionalTermPartialFirstY/1000.0-0.0093063*K*concentrationAdditionalTermPartialFirstY/1000.0+0.000008069*K*K*concentrationAdditionalTermPartialFirstY/1000.0+0.44004*Concentration[i][j][k]*constantMol/1000.0*concentrationAdditionalTermPartialFirstY/1000.0-0.0001765*2*K*Concentration[i][j][k]*constantMol/1000.0*concentrationAdditionalTermPartialFirstY/1000.0)/1000.0;
						Part1PartialFirstZ = constantLelectrolyteConductivity*exchange*1.2544*2*constant1*K/constant2*2*otherAdditionalTermPart1Part*(0.26235*concentrationAdditionalTermPartialFirstZ/1000.0-0.0093063*K*concentrationAdditionalTermPartialFirstZ/1000.0+0.000008069*K*K*concentrationAdditionalTermPartialFirstZ/1000.0+0.44004*Concentration[i][j][k]*constantMol/1000.0*concentrationAdditionalTermPartialFirstZ/1000.0-0.0001765*2*K*Concentration[i][j][k]*constantMol/1000.0*concentrationAdditionalTermPartialFirstZ/1000.0)/1000.0;
						Part2PartialFirstX = 0.637/0.601*(-0.24*0.5/sqrt(Concentration[i][j][k]*constantMol/1000.0)*concentrationAdditionalTermPartialFirstX/1000.0+0.982*1.5*(1-0.0052*(K-294))*sqrt(Concentration[i][j][k]*constantMol/1000.0)*concentrationAdditionalTermPartialFirstX/1000.0);
						Part2PartialFirstY = 0.637/0.601*(-0.24*0.5/sqrt(Concentration[i][j][k]*constantMol/1000.0)*concentrationAdditionalTermPartialFirstY/1000.0+0.982*1.5*(1-0.0052*(K-294))*sqrt(Concentration[i][j][k]*constantMol/1000.0)*concentrationAdditionalTermPartialFirstY/1000.0);
						Part2PartialFirstZ = 0.637/0.601*(-0.24*0.5/sqrt(Concentration[i][j][k]*constantMol/1000.0)*concentrationAdditionalTermPartialFirstZ/1000.0+0.982*1.5*(1-0.0052*(K-294))*sqrt(Concentration[i][j][k]*constantMol/1000.0)*concentrationAdditionalTermPartialFirstZ/1000.0);
						otherAdditionalTermPartialFirstX = Part1PartialFirstX*otherAdditionalTermPart2+otherAdditionalTermPart1*Part2PartialFirstX;
                        otherAdditionalTermPartialFirstY = Part1PartialFirstY*otherAdditionalTermPart2+otherAdditionalTermPart1*Part2PartialFirstY;
						otherAdditionalTermPartialFirstZ = Part1PartialFirstZ*otherAdditionalTermPart2+otherAdditionalTermPart1*Part2PartialFirstZ;
						electrolyteAdditionalTerm[i][j][k]= -(otherAdditionalTermPartialFirstX*concentrationAdditionalTermPartialFirstX+otherAdditionalTermPartialFirstY*concentrationAdditionalTermPartialFirstY+otherAdditionalTermPartialFirstZ*concentrationAdditionalTermPartialFirstZ)-otherAdditionalTerm*concentrationAdditionalTermPartialSecond;
					}
				}
			}
		}
	}
}

void CollisionSRTelectrolytePotential(int Length,int Width,int Height,double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						Out0[i][j][k]=electrolyteIn0[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn0[i][j][k]-tNS[0]*electrolytePotential[i][j][k]);//+tNS[0]*electrolyteAdditionalTerm[i][j][k];
						Out1[i][j][k]=electrolyteIn1[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn1[i][j][k]-tNS[1]*electrolytePotential[i][j][k]);//+tNS[1]*electrolyteAdditionalTerm[i][j][k];
						Out2[i][j][k]=electrolyteIn2[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn2[i][j][k]-tNS[2]*electrolytePotential[i][j][k]);//+tNS[2]*electrolyteAdditionalTerm[i][j][k];
						Out3[i][j][k]=electrolyteIn3[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn3[i][j][k]-tNS[3]*electrolytePotential[i][j][k]);//+tNS[3]*electrolyteAdditionalTerm[i][j][k];
						Out4[i][j][k]=electrolyteIn4[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn4[i][j][k]-tNS[4]*electrolytePotential[i][j][k]);//+tNS[4]*electrolyteAdditionalTerm[i][j][k];					
						Out5[i][j][k]=electrolyteIn5[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn5[i][j][k]-tNS[5]*electrolytePotential[i][j][k]);//+tNS[5]*electrolyteAdditionalTerm[i][j][k];
						Out6[i][j][k]=electrolyteIn6[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn6[i][j][k]-tNS[6]*electrolytePotential[i][j][k]);//+tNS[6]*electrolyteAdditionalTerm[i][j][k];
						Out7[i][j][k]=electrolyteIn7[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn7[i][j][k]-tNS[7]*electrolytePotential[i][j][k]);//+tNS[7]*electrolyteAdditionalTerm[i][j][k];
						Out8[i][j][k]=electrolyteIn8[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn8[i][j][k]-tNS[8]*electrolytePotential[i][j][k]);//+tNS[8]*electrolyteAdditionalTerm[i][j][k];
						Out9[i][j][k]=electrolyteIn9[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn9[i][j][k]-tNS[9]*electrolytePotential[i][j][k]);//+tNS[9]*electrolyteAdditionalTerm[i][j][k];
						Out10[i][j][k]=electrolyteIn10[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn10[i][j][k]-tNS[10]*electrolytePotential[i][j][k]);//+tNS[10]*electrolyteAdditionalTerm[i][j][k];
						Out11[i][j][k]=electrolyteIn11[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn11[i][j][k]-tNS[11]*electrolytePotential[i][j][k]);//+tNS[11]*electrolyteAdditionalTerm[i][j][k];
						Out12[i][j][k]=electrolyteIn12[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn12[i][j][k]-tNS[12]*electrolytePotential[i][j][k]);//+tNS[12]*electrolyteAdditionalTerm[i][j][k];
						Out13[i][j][k]=electrolyteIn13[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn13[i][j][k]-tNS[13]*electrolytePotential[i][j][k]);//+tNS[13]*electrolyteAdditionalTerm[i][j][k];
						Out14[i][j][k]=electrolyteIn14[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn14[i][j][k]-tNS[14]*electrolytePotential[i][j][k]);//+tNS[14]*electrolyteAdditionalTerm[i][j][k];
						Out15[i][j][k]=electrolyteIn15[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn15[i][j][k]-tNS[15]*electrolytePotential[i][j][k]);//+tNS[15]*electrolyteAdditionalTerm[i][j][k];
						Out16[i][j][k]=electrolyteIn16[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn16[i][j][k]-tNS[16]*electrolytePotential[i][j][k]);//+tNS[16]*electrolyteAdditionalTerm[i][j][k];
						Out17[i][j][k]=electrolyteIn17[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn17[i][j][k]-tNS[17]*electrolytePotential[i][j][k]);//+tNS[17]*electrolyteAdditionalTerm[i][j][k];
						Out18[i][j][k]=electrolyteIn18[i][j][k]-electrolyteRelaxationTime[i][j][k]*(electrolyteIn18[i][j][k]-tNS[18]*electrolytePotential[i][j][k]);//+tNS[18]*electrolyteAdditionalTerm[i][j][k];
					}
				}
			}
		}
	}
}

void CollisionSRTelectrolytePotentialD3Q7(int Length,int Width,int Height,double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						Out0[i][j][k] = electrolyteIn0[i][j][k] - electrolyteRelaxationTime[i][j][k] * (electrolyteIn0[i][j][k] - tNSD3Q7[0] * electrolytePotential[i][j][k]) +tNSD3Q7[0] * electrolyteAdditionalTerm[i][j][k];
						Out1[i][j][k] = electrolyteIn1[i][j][k] - electrolyteRelaxationTime[i][j][k] * (electrolyteIn1[i][j][k] - tNSD3Q7[1] * electrolytePotential[i][j][k]) +tNSD3Q7[1] * electrolyteAdditionalTerm[i][j][k];
						Out2[i][j][k] = electrolyteIn2[i][j][k] - electrolyteRelaxationTime[i][j][k] * (electrolyteIn2[i][j][k] - tNSD3Q7[2] * electrolytePotential[i][j][k]) +tNSD3Q7[2] * electrolyteAdditionalTerm[i][j][k];
						Out3[i][j][k] = electrolyteIn3[i][j][k] - electrolyteRelaxationTime[i][j][k] * (electrolyteIn3[i][j][k] - tNSD3Q7[3] * electrolytePotential[i][j][k]) +tNSD3Q7[3] * electrolyteAdditionalTerm[i][j][k];
						Out4[i][j][k] = electrolyteIn4[i][j][k] - electrolyteRelaxationTime[i][j][k] * (electrolyteIn4[i][j][k] - tNSD3Q7[4] * electrolytePotential[i][j][k]) +tNSD3Q7[4] * electrolyteAdditionalTerm[i][j][k];
						Out5[i][j][k] = electrolyteIn5[i][j][k] - electrolyteRelaxationTime[i][j][k] * (electrolyteIn5[i][j][k] - tNSD3Q7[5] * electrolytePotential[i][j][k]) +tNSD3Q7[5] * electrolyteAdditionalTerm[i][j][k];
						Out6[i][j][k] = electrolyteIn6[i][j][k] - electrolyteRelaxationTime[i][j][k] * (electrolyteIn6[i][j][k] - tNSD3Q7[6] * electrolytePotential[i][j][k]) +tNSD3Q7[6] * electrolyteAdditionalTerm[i][j][k];
					}
				}
			}
		}
	}
}

void MemoryFreeelectrolytePotential(int m, int n, int q){
	freememory(electrolytePotential,m,n,q);
	freememory(electrolyteIn0,m,n,q);
	freememory(electrolyteIn1,m,n,q);
	freememory(electrolyteIn2,m,n,q);
	freememory(electrolyteIn3,m,n,q);
	freememory(electrolyteIn4,m,n,q);
	freememory(electrolyteIn5,m,n,q);
	freememory(electrolyteIn6,m,n,q);
//	freememory(electrolyteIn7,m,n,q);
//	freememory(electrolyteIn8,m,n,q);
//	freememory(electrolyteIn9,m,n,q);
//	freememory(electrolyteIn10,m,n,q);
//	freememory(electrolyteIn11,m,n,q);
//	freememory(electrolyteIn12,m,n,q);
//	freememory(electrolyteIn13,m,n,q);
//	freememory(electrolyteIn14,m,n,q);
//	freememory(electrolyteIn15,m,n,q);
//	freememory(electrolyteIn16,m,n,q);
//	freememory(electrolyteIn17,m,n,q);
//	freememory(electrolyteIn18,m,n,q);
	freememory(electrolyteRelaxationTime,m,n,q);
	freememory(electrolyteAdditionalTerm,m,n,q);
	freememory(lastWallElectrolyte,m,n,q);
}

void StreamelectrolytePotential(int m, int n, int q){
	shift0(Out0,electrolyteIn0,m,n,q);
	shift1(Out1,electrolyteIn1,m,n,q);
	shift2(Out2,electrolyteIn2,m,n,q);
	shift3(Out3,electrolyteIn3,m,n,q);
	shift4(Out4,electrolyteIn4,m,n,q);
	shift5(Out5,electrolyteIn5,m,n,q);
	shift6(Out6,electrolyteIn6,m,n,q);
	shift7(Out7,electrolyteIn7,m,n,q);
	shift8(Out8,electrolyteIn8,m,n,q);
	shift9(Out9,electrolyteIn9,m,n,q);
	shift10(Out10,electrolyteIn10,m,n,q);
	shift11(Out11,electrolyteIn11,m,n,q);
	shift12(Out12,electrolyteIn12,m,n,q);
	shift13(Out13,electrolyteIn13,m,n,q);
	shift14(Out14,electrolyteIn14,m,n,q);
	shift15(Out15,electrolyteIn15,m,n,q);
	shift16(Out16,electrolyteIn16,m,n,q);
	shift17(Out17,electrolyteIn17,m,n,q);
	shift18(Out18,electrolyteIn18,m,n,q);
}

void StreamelectrolytePotentialD3Q7(int m, int n, int q){
	shift0(Out0,electrolyteIn0,m,n,q);
	shift1(Out1,electrolyteIn1,m,n,q);
	shift2(Out2,electrolyteIn2,m,n,q);
	shift3(Out3,electrolyteIn3,m,n,q);
	shift4(Out4,electrolyteIn4,m,n,q);
	shift5(Out5,electrolyteIn5,m,n,q);
	shift6(Out6,electrolyteIn6,m,n,q);
}

void FieldCalculationelectrolytePotential(int m, int n, int q, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						electrolytePotential[i][j][k]=electrolyteIn0[i][j][k]+electrolyteIn1[i][j][k]+electrolyteIn2[i][j][k]+electrolyteIn3[i][j][k]+electrolyteIn4[i][j][k]+electrolyteIn5[i][j][k]+electrolyteIn6[i][j][k]+electrolyteIn7[i][j][k]+electrolyteIn8[i][j][k]\
							+electrolyteIn9[i][j][k]+electrolyteIn10[i][j][k]+electrolyteIn11[i][j][k]+electrolyteIn12[i][j][k]+electrolyteIn13[i][j][k]+electrolyteIn14[i][j][k]+electrolyteIn15[i][j][k]+electrolyteIn16[i][j][k]+electrolyteIn17[i][j][k]+electrolyteIn18[i][j][k];
					}
					if (domain[i][j][k]!=0||rho[i][j][k]<=gasCritical){
						electrolytePotential[i][j][k]=0.0;
						electrolyteIn0[i][j][k]=0.0;
						electrolyteIn1[i][j][k]=0.0;
						electrolyteIn2[i][j][k]=0.0;
						electrolyteIn3[i][j][k]=0.0;
						electrolyteIn4[i][j][k]=0.0;
						electrolyteIn5[i][j][k]=0.0;
						electrolyteIn6[i][j][k]=0.0;
						electrolyteIn7[i][j][k]=0.0;
						electrolyteIn8[i][j][k]=0.0;
						electrolyteIn9[i][j][k]=0.0;
						electrolyteIn10[i][j][k]=0.0;
						electrolyteIn11[i][j][k]=0.0;
						electrolyteIn12[i][j][k]=0.0;
						electrolyteIn13[i][j][k]=0.0;
						electrolyteIn14[i][j][k]=0.0;
						electrolyteIn15[i][j][k]=0.0;
						electrolyteIn16[i][j][k]=0.0;
						electrolyteIn17[i][j][k]=0.0;
						electrolyteIn18[i][j][k]=0.0;
					}
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]!=1){
						electrolyteIn0[i][j][k]=electrolytePotential[i][j][k]*tNS[0];
						electrolyteIn1[i][j][k]=electrolytePotential[i][j][k]*tNS[1];
						electrolyteIn2[i][j][k]=electrolytePotential[i][j][k]*tNS[2];
						electrolyteIn3[i][j][k]=electrolytePotential[i][j][k]*tNS[3];
						electrolyteIn4[i][j][k]=electrolytePotential[i][j][k]*tNS[4];
						electrolyteIn5[i][j][k]=electrolytePotential[i][j][k]*tNS[5];
						electrolyteIn6[i][j][k]=electrolytePotential[i][j][k]*tNS[6];
						electrolyteIn7[i][j][k]=electrolytePotential[i][j][k]*tNS[7];
						electrolyteIn8[i][j][k]=electrolytePotential[i][j][k]*tNS[8];
						electrolyteIn9[i][j][k]=electrolytePotential[i][j][k]*tNS[9];
						electrolyteIn10[i][j][k]=electrolytePotential[i][j][k]*tNS[10];
						electrolyteIn11[i][j][k]=electrolytePotential[i][j][k]*tNS[11];
						electrolyteIn12[i][j][k]=electrolytePotential[i][j][k]*tNS[12];
						electrolyteIn13[i][j][k]=electrolytePotential[i][j][k]*tNS[13];
						electrolyteIn14[i][j][k]=electrolytePotential[i][j][k]*tNS[14];
						electrolyteIn15[i][j][k]=electrolytePotential[i][j][k]*tNS[15];
						electrolyteIn16[i][j][k]=electrolytePotential[i][j][k]*tNS[16];
						electrolyteIn17[i][j][k]=electrolytePotential[i][j][k]*tNS[17];
						electrolyteIn18[i][j][k]=electrolytePotential[i][j][k]*tNS[18];
					}					
				}
			}
		}
	}
}

void FieldCalculationelectrolytePotentialD3Q7(int m, int n, int q, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						electrolytePotential[i][j][k]=electrolyteIn0[i][j][k]+electrolyteIn1[i][j][k]+electrolyteIn2[i][j][k]+electrolyteIn3[i][j][k]+electrolyteIn4[i][j][k]+electrolyteIn5[i][j][k]+electrolyteIn6[i][j][k];
					}
					if (domain[i][j][k]!=0||rho[i][j][k]<=gasCritical||domainConnection[i][j][k]!=1){
						electrolytePotential[i][j][k]=0.0;
						electrolyteIn0[i][j][k]=0.0;
						electrolyteIn1[i][j][k]=0.0;
						electrolyteIn2[i][j][k]=0.0;
						electrolyteIn3[i][j][k]=0.0;
						electrolyteIn4[i][j][k]=0.0;
						electrolyteIn5[i][j][k]=0.0;
						electrolyteIn6[i][j][k]=0.0;
					}
/*					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]!=1){
						electrolyteIn0[i][j][k]=electrolytePotential[i][j][k]*tNSD3Q7[0];
						electrolyteIn1[i][j][k]=electrolytePotential[i][j][k]*tNSD3Q7[1];
						electrolyteIn2[i][j][k]=electrolytePotential[i][j][k]*tNSD3Q7[2];
						electrolyteIn3[i][j][k]=electrolytePotential[i][j][k]*tNSD3Q7[3];
						electrolyteIn4[i][j][k]=electrolytePotential[i][j][k]*tNSD3Q7[4];
						electrolyteIn5[i][j][k]=electrolytePotential[i][j][k]*tNSD3Q7[5];
						electrolyteIn6[i][j][k]=electrolytePotential[i][j][k]*tNSD3Q7[6];
					}*/					
				}
			}
		}
	}
}

void lastElectrolyte(int Length, int Width, int Height, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
						lastWallElectrolyte[i][j][k]=electrolytePotential[i][j][k];
				}
			}
		}
	}
}

void reactionElectrolyteD3Q7(int m, int n, int q, double lastTotalCurrentResidual){
	int i;
	int j;
	int k;
	double reactant;
#pragma omp parallel private(i,j,k,reactant)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&domainConnection[i][j][k]==1){
						if (Domain1[i][j][k]==1){
							reactant = -neightbour4(i, j, k, m, reaction4)*(1 + lastTotalCurrentResidual);
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * reactant;
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * reactant;
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * reactant;
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * reactant;
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * reactant;
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * reactant;
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * reactant;
							}
						if (Domain2[i][j][k] == 1 && j != 0){
							reactant = -neightbour5(i, j, k, n, reaction5)*(1 + lastTotalCurrentResidual);
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * reactant;
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * reactant;
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * reactant;
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * reactant;
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * reactant;
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * reactant;
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * reactant;
							}
						if (Domain3[i][j][k] == 1){
							reactant = -neightbour6(i, j, k, q, reaction6)*(1 + lastTotalCurrentResidual);
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * reactant;
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * reactant;
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * reactant;
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * reactant;
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * reactant;
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * reactant;
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * reactant;
							}
						if (Domain4[i][j][k] == 1){
							reactant = -neightbour1(i, j, k, m, reaction1)*(1 + lastTotalCurrentResidual);
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * reactant;
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * reactant;
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * reactant;
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * reactant;
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * reactant;
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * reactant;
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * reactant;
							}
						if (Domain5[i][j][k] == 1){
							reactant = -neightbour2(i, j, k, n, reaction2)*(1 + lastTotalCurrentResidual);
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * reactant;
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * reactant;
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * reactant;
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * reactant;
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * reactant;
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * reactant;
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * reactant;
							}
						if (Domain6[i][j][k] == 1){
							reactant = -neightbour3(i, j, k, q, reaction3)*(1 + lastTotalCurrentResidual);
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * reactant;
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * reactant;
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * reactant;
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * reactant;
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * reactant;
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * reactant;
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * reactant;
							}
					}
				}
			}
		}
	}
}