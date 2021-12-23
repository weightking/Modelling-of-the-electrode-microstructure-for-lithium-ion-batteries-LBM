#include <stdio.h>
#include <stdlib.h>
#include "twoPhaseField.h"
#include "domain.h"

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

double divergeNumberTwoPhase(int m, int n, int q, double FLUID){
	int i;
	int j;
	int k;
	double numberdivergence=0.0;
	for (i=0;i<m;i++){
		for (j=0;j<n;j++){
			for (k=0;k<q;k++){
				if (rho[i][j][k]<0||rho[i][j][k]>FLUID*3)
					numberdivergence=numberdivergence+1;
			}
		}
	}
	return numberdivergence;
}

void divergeFixTwoPhase(int m, int n, int q,double GAS, double FLUID, double b, double R, double T, double a){
	int i;
	int j;
	int k;
	double bRho4;
	double BRho4;
#pragma omp parallel private(i,j,k,bRho4,BRho4)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0){
						if (rho[i][j][k]<0){
							rho[i][j][k]=GAS;
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
							rho[i][j][k]=fIn0[i][j][k]+fIn1[i][j][k]+fIn2[i][j][k]+fIn3[i][j][k]+fIn4[i][j][k]+fIn5[i][j][k]+fIn6[i][j][k]+fIn7[i][j][k]+fIn8[i][j][k]\
								+fIn9[i][j][k]+fIn10[i][j][k]+fIn11[i][j][k]+fIn12[i][j][k]+fIn13[i][j][k]+fIn14[i][j][k]+fIn15[i][j][k]+fIn16[i][j][k]+fIn17[i][j][k]+fIn18[i][j][k];
							ux[i][j][k]=uTotX+Fx[i][j][k]/2.0/rho[i][j][k];
							uy[i][j][k]=uTotY+Fy[i][j][k]/2.0/rho[i][j][k];
							uz[i][j][k]=uTotZ+Fz[i][j][k]/2.0/rho[i][j][k];
							bRho4=b*rho[i][j][k]/4.0;
							BRho4=1.0-bRho4;
							pressure[i][j][k]=Pressure;
						}
						if (rho[i][j][k]>FLUID*3){
							rho[i][j][k]=FLUID;
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
							rho[i][j][k]=fIn0[i][j][k]+fIn1[i][j][k]+fIn2[i][j][k]+fIn3[i][j][k]+fIn4[i][j][k]+fIn5[i][j][k]+fIn6[i][j][k]+fIn7[i][j][k]+fIn8[i][j][k]\
								+fIn9[i][j][k]+fIn10[i][j][k]+fIn11[i][j][k]+fIn12[i][j][k]+fIn13[i][j][k]+fIn14[i][j][k]+fIn15[i][j][k]+fIn16[i][j][k]+fIn17[i][j][k]+fIn18[i][j][k];
							ux[i][j][k]=uTotX+Fx[i][j][k]/2.0/rho[i][j][k];
							uy[i][j][k]=uTotY+Fy[i][j][k]/2.0/rho[i][j][k];
							uz[i][j][k]=uTotZ+Fz[i][j][k]/2.0/rho[i][j][k];
							bRho4=b*rho[i][j][k]/4.0;
							BRho4=1.0-bRho4;
							pressure[i][j][k]=Pressure;
						}
					}
				}
			}
		}
	}
}