#include <stdio.h>
#include <stdlib.h>
#include "matrixMove.h"
#include "twoPhaseField.h"

#define jx (cxNS[0]*fIn0[i][j][k]+cxNS[1]*fIn1[i][j][k]+cxNS[2]*fIn2[i][j][k]+cxNS[3]*fIn3[i][j][k]+cxNS[4]*fIn4[i][j][k]+cxNS[5]*fIn5[i][j][k]+cxNS[6]*fIn6[i][j][k]+cxNS[7]*fIn7[i][j][k]+cxNS[8]*fIn8[i][j][k]\
	+cxNS[9]*fIn9[i][j][k]+cxNS[10]*fIn10[i][j][k]+cxNS[11]*fIn11[i][j][k]+cxNS[12]*fIn12[i][j][k]+cxNS[13]*fIn13[i][j][k]+cxNS[14]*fIn14[i][j][k]+cxNS[15]*fIn15[i][j][k]+cxNS[16]*fIn16[i][j][k]+cxNS[17]*fIn17[i][j][k]+cxNS[18]*fIn18[i][j][k])
#define jy (cyNS[0]*fIn0[i][j][k]+cyNS[1]*fIn1[i][j][k]+cyNS[2]*fIn2[i][j][k]+cyNS[3]*fIn3[i][j][k]+cyNS[4]*fIn4[i][j][k]+cyNS[5]*fIn5[i][j][k]+cyNS[6]*fIn6[i][j][k]+cyNS[7]*fIn7[i][j][k]+cyNS[8]*fIn8[i][j][k]\
	+cyNS[9]*fIn9[i][j][k]+cyNS[10]*fIn10[i][j][k]+cyNS[11]*fIn11[i][j][k]+cyNS[12]*fIn12[i][j][k]+cyNS[13]*fIn13[i][j][k]+cyNS[14]*fIn14[i][j][k]+cyNS[15]*fIn15[i][j][k]+cyNS[16]*fIn16[i][j][k]+cyNS[17]*fIn17[i][j][k]+cyNS[18]*fIn18[i][j][k])
#define jz (czNS[0]*fIn0[i][j][k]+czNS[1]*fIn1[i][j][k]+czNS[2]*fIn2[i][j][k]+czNS[3]*fIn3[i][j][k]+czNS[4]*fIn4[i][j][k]+czNS[5]*fIn5[i][j][k]+czNS[6]*fIn6[i][j][k]+czNS[7]*fIn7[i][j][k]+czNS[8]*fIn8[i][j][k]\
	+czNS[9]*fIn9[i][j][k]+czNS[10]*fIn10[i][j][k]+czNS[11]*fIn11[i][j][k]+czNS[12]*fIn12[i][j][k]+czNS[13]*fIn13[i][j][k]+czNS[14]*fIn14[i][j][k]+czNS[15]*fIn15[i][j][k]+czNS[16]*fIn16[i][j][k]+czNS[17]*fIn17[i][j][k]+czNS[18]*fIn18[i][j][k])
#define uTotX (jx/rho[i][j][k])
#define uTotY (jy/rho[i][j][k])
#define uTotZ (jz/rho[i][j][k])
#define Pressure rho[i][j][k]*R*T*(1+bRho4+bRho4*bRho4-bRho4*bRho4*bRho4)/(BRho4*BRho4*BRho4)-a*rho[i][j][k]*rho[i][j][k]
//#define Pressure rho[i][j][k]*R*T/(1-b*rho[i][j][k])-a*alpha*rho[i][j][k]*rho[i][j][k]/(1+2*b*rho[i][j][k]-b*b*rho[i][j][k]*rho[i][j][k])

void initialWallBoundaryTwoPhase(int mInitial,int m,int nInitial,int n,int qInitial,int q,double ***domain,
	double ***Domain1,double ***Domain2,double ***Domain3,double ***Domain4,double ***Domain5,double ***Domain6,double ***Domain7,double ***Domain8,double ***Domain9,double ***Domain10,
	double ***Domain11,double ***Domain12,double ***Domain13,double ***Domain14,double ***Domain15,double ***Domain16,double ***Domain17,double ***Domain18,double constant,double rhoLayer,
	double b, double R, double T, double a){
		int i;
		int j;
		int k;
		double bRho4;
		double BRho4;
#pragma omp parallel private(i,j,k,bRho4,BRho4)
		{
#pragma omp for schedule(dynamic)
			for (i=mInitial;i<m;i++){
				for (j=nInitial;j<n;j++){
					for (k=qInitial;k<q;k++){
						if (domain[i][j][k]==0&&(Domain1[i][j][k]==constant||Domain2[i][j][k]==constant||Domain3[i][j][k]==constant||Domain4[i][j][k]==constant||\
							Domain5[i][j][k]==constant||Domain6[i][j][k]==constant||Domain7[i][j][k]==constant||Domain8[i][j][k]==constant||Domain9[i][j][k]==constant||\
							Domain10[i][j][k]==constant||Domain11[i][j][k]==constant||Domain12[i][j][k]==constant||Domain13[i][j][k]==constant||Domain14[i][j][k]==constant||\
							Domain15[i][j][k]==constant||Domain16[i][j][k]==constant||Domain17[i][j][k]==constant||Domain18[i][j][k]==constant)){
								rho[i][j][k]=rhoLayer;
								bRho4=b*rho[i][j][k]/4.0;
								BRho4=1.0-bRho4;
								pressure[i][j][k]=Pressure;
								fIn0[i][j][k]=rho[i][j][k]*tNS[0];
								fIn1[i][j][k]=rho[i][j][k]*tNS[1];
								fIn2[i][j][k]=rho[i][j][k]*tNS[2];
								fIn3[i][j][k]=rho[i][j][k]*tNS[3];
								fIn4[i][j][k]=rho[i][j][k]*tNS[4];
								fIn5[i][j][k]=rho[i][j][k]*tNS[5];
								fIn6[i][j][k]=rho[i][j][k]*tNS[6];
								fIn7[i][j][k]=rho[i][j][k]*tNS[7];
								fIn8[i][j][k]=rho[i][j][k]*tNS[8];
								fIn9[i][j][k]=rho[i][j][k]*tNS[9];
								fIn10[i][j][k]=rho[i][j][k]*tNS[10];
								fIn11[i][j][k]=rho[i][j][k]*tNS[11];
								fIn12[i][j][k]=rho[i][j][k]*tNS[12];
								fIn13[i][j][k]=rho[i][j][k]*tNS[13];
								fIn14[i][j][k]=rho[i][j][k]*tNS[14];
								fIn15[i][j][k]=rho[i][j][k]*tNS[15];
								fIn16[i][j][k]=rho[i][j][k]*tNS[16];
								fIn17[i][j][k]=rho[i][j][k]*tNS[17];
								fIn18[i][j][k]=rho[i][j][k]*tNS[18];
								ux[i][j][k]=uTotX+Fx[i][j][k]/2.0/rho[i][j][k];
								uy[i][j][k]=uTotY+Fy[i][j][k]/2.0/rho[i][j][k];
								uz[i][j][k]=uTotZ+Fz[i][j][k]/2.0/rho[i][j][k];
						}
					}
				}
			}
		}
}

void initialWallNearestTwoPhase(int mInitial,int m,int nInitial,int n,int qInitial,int q,double ***domain,
	double ***Domain1,double ***Domain2,double ***Domain3,double ***Domain4,double ***Domain5,double ***Domain6,double ***Domain7,double ***Domain8,double ***Domain9,double ***Domain10,
	double ***Domain11,double ***Domain12,double ***Domain13,double ***Domain14,double ***Domain15,double ***Domain16,double ***Domain17,double ***Domain18,double constant,double rhoLayer,
	double b, double R, double T, double a, double GAS){
		int i;
		int j;
		int k;
		double bRho4;
		double BRho4;
#pragma omp parallel private(i,j,k,bRho4,BRho4)
		{
#pragma omp for schedule(dynamic)
			for (i=mInitial;i<m;i++){
				for (j=nInitial;j<n;j++){
					for (k=qInitial;k<q;k++){
						if (domain[i][j][k]==0&&rho[i][j][k]==GAS&&(Domain1[i][j][k]==constant||Domain2[i][j][k]==constant||Domain3[i][j][k]==constant||Domain4[i][j][k]==constant||\
							Domain5[i][j][k]==constant||Domain6[i][j][k]==constant||Domain7[i][j][k]==constant||Domain8[i][j][k]==constant||Domain9[i][j][k]==constant||\
							Domain10[i][j][k]==constant||Domain11[i][j][k]==constant||Domain12[i][j][k]==constant||Domain13[i][j][k]==constant||Domain14[i][j][k]==constant||\
							Domain15[i][j][k]==constant||Domain16[i][j][k]==constant||Domain17[i][j][k]==constant||Domain18[i][j][k]==constant)){
								rho[i][j][k]=rhoLayer;
								bRho4=b*rho[i][j][k]/4.0;
								BRho4=1.0-bRho4;
								pressure[i][j][k]=Pressure;
								fIn0[i][j][k]=rho[i][j][k]*tNS[0];
								fIn1[i][j][k]=rho[i][j][k]*tNS[1];
								fIn2[i][j][k]=rho[i][j][k]*tNS[2];
								fIn3[i][j][k]=rho[i][j][k]*tNS[3];
								fIn4[i][j][k]=rho[i][j][k]*tNS[4];
								fIn5[i][j][k]=rho[i][j][k]*tNS[5];
								fIn6[i][j][k]=rho[i][j][k]*tNS[6];
								fIn7[i][j][k]=rho[i][j][k]*tNS[7];
								fIn8[i][j][k]=rho[i][j][k]*tNS[8];
								fIn9[i][j][k]=rho[i][j][k]*tNS[9];
								fIn10[i][j][k]=rho[i][j][k]*tNS[10];
								fIn11[i][j][k]=rho[i][j][k]*tNS[11];
								fIn12[i][j][k]=rho[i][j][k]*tNS[12];
								fIn13[i][j][k]=rho[i][j][k]*tNS[13];
								fIn14[i][j][k]=rho[i][j][k]*tNS[14];
								fIn15[i][j][k]=rho[i][j][k]*tNS[15];
								fIn16[i][j][k]=rho[i][j][k]*tNS[16];
								fIn17[i][j][k]=rho[i][j][k]*tNS[17];
								fIn18[i][j][k]=rho[i][j][k]*tNS[18];
								ux[i][j][k]=uTotX+Fx[i][j][k]/2.0/rho[i][j][k];
								uy[i][j][k]=uTotY+Fy[i][j][k]/2.0/rho[i][j][k];
								uz[i][j][k]=uTotZ+Fz[i][j][k]/2.0/rho[i][j][k];
						}
					}
				}
			}
		}
}
