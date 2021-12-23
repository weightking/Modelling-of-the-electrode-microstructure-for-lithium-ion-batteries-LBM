#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "domain.h"
#include "OutMemoryArrange.h"
#include "twoPhaseField.h"
#include "concentrationField.h"
#include "ghostFlowField.h"
#include "electrodePotential.h"
#include "electrolytePotential.h"
#include "matrixMove.h"

void wallBoundaryTwoPhase(int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0){
						if (Domain1[i][j][k]!=0)
							fIn1[i][j][k]=Out4[i][j][k];
						if (Domain2[i][j][k]!=0)
							fIn2[i][j][k]=Out5[i][j][k];
						if (Domain3[i][j][k]!=0)
							fIn3[i][j][k]=Out6[i][j][k];
						if (Domain4[i][j][k]!=0)
							fIn4[i][j][k]=Out1[i][j][k];
						if (Domain5[i][j][k]!=0)
							fIn5[i][j][k]=Out2[i][j][k];
						if (Domain6[i][j][k]!=0)
							fIn6[i][j][k]=Out3[i][j][k];
						if (Domain7[i][j][k]!=0)
							fIn7[i][j][k]=Out10[i][j][k];
						if (Domain8[i][j][k]!=0)
							fIn8[i][j][k]=Out11[i][j][k];
						if (Domain9[i][j][k]!=0)
							fIn9[i][j][k]=Out12[i][j][k];
						if (Domain10[i][j][k]!=0)
							fIn10[i][j][k]=Out7[i][j][k];
						if (Domain11[i][j][k]!=0)
							fIn11[i][j][k]=Out8[i][j][k];
						if (Domain12[i][j][k]!=0)
							fIn12[i][j][k]=Out9[i][j][k];
						if (Domain13[i][j][k]!=0)
							fIn13[i][j][k]=Out16[i][j][k];
						if (Domain14[i][j][k]!=0)
							fIn14[i][j][k]=Out17[i][j][k];
						if (Domain15[i][j][k]!=0)
							fIn15[i][j][k]=Out18[i][j][k];
						if (Domain16[i][j][k]!=0)
							fIn16[i][j][k]=Out13[i][j][k];
						if (Domain17[i][j][k]!=0)
							fIn17[i][j][k]=Out14[i][j][k];
						if (Domain18[i][j][k]!=0)
							fIn18[i][j][k]=Out15[i][j][k];
					}
				}
			}
		}
	}
}

void wallBoundaryConcentration(int m, int n, int q, double gasCritical){
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
						if (Domain1[i][j][k]!=0)
							cIn1[i][j][k]=-Out4[i][j][k]+2*tNS[1]*Concentration[i][j][k];
//                            cIn1[i][j][k]=Out4[i][j][k];
						if (Domain2[i][j][k]!=0)
							cIn2[i][j][k]=-Out5[i][j][k]+2*tNS[2]*Concentration[i][j][k];
//                            cIn2[i][j][k]=Out5[i][j][k];
						if (Domain3[i][j][k]!=0)
							cIn3[i][j][k]=-Out6[i][j][k]+2*tNS[3]*Concentration[i][j][k];
//                            cIn3[i][j][k]=Out6[i][j][k];
						if (Domain4[i][j][k]!=0)
							cIn4[i][j][k]=-Out1[i][j][k]+2*tNS[4]*Concentration[i][j][k];
//                            cIn4[i][j][k]=Out1[i][j][k];
						if (Domain5[i][j][k]!=0)
							cIn5[i][j][k]=-Out2[i][j][k]+2*tNS[5]*Concentration[i][j][k];
//                            cIn5[i][j][k]=Out2[i][j][k];
						if (Domain6[i][j][k]!=0)
							cIn6[i][j][k]=-Out3[i][j][k]+2*tNS[6]*Concentration[i][j][k];
//                            cIn6[i][j][k]=Out3[i][j][k];
						if (Domain7[i][j][k]!=0)
							cIn7[i][j][k]=-Out10[i][j][k]+2*tNS[7]*Concentration[i][j][k];
//                            cIn7[i][j][k]=Out10[i][j][k];
						if (Domain8[i][j][k]!=0)
							cIn8[i][j][k]=-Out11[i][j][k]+2*tNS[8]*Concentration[i][j][k];
//                            cIn8[i][j][k]=Out11[i][j][k];
						if (Domain9[i][j][k]!=0)
							cIn9[i][j][k]=-Out12[i][j][k]+2*tNS[9]*Concentration[i][j][k];
//                            cIn9[i][j][k]=Out12[i][j][k];
						if (Domain10[i][j][k]!=0)
							cIn10[i][j][k]=-Out7[i][j][k]+2*tNS[10]*Concentration[i][j][k];
//                            cIn10[i][j][k]=Out7[i][j][k];
						if (Domain11[i][j][k]!=0)
							cIn11[i][j][k]=-Out8[i][j][k]+2*tNS[11]*Concentration[i][j][k];
//                            cIn11[i][j][k]=Out8[i][j][k];
						if (Domain12[i][j][k]!=0)
							cIn12[i][j][k]=-Out9[i][j][k]+2*tNS[12]*Concentration[i][j][k];
//                            cIn12[i][j][k]=Out9[i][j][k];
						if (Domain13[i][j][k]!=0)
							cIn13[i][j][k]=-Out16[i][j][k]+2*tNS[13]*Concentration[i][j][k];
//                            cIn13[i][j][k]=Out16[i][j][k];
						if (Domain14[i][j][k]!=0)
							cIn14[i][j][k]=-Out17[i][j][k]+2*tNS[14]*Concentration[i][j][k];
//                            cIn14[i][j][k]=Out17[i][j][k];
						if (Domain15[i][j][k]!=0)
							cIn15[i][j][k]=-Out18[i][j][k]+2*tNS[15]*Concentration[i][j][k];
//                            cIn15[i][j][k]=Out18[i][j][k];
						if (Domain16[i][j][k]!=0)
							cIn16[i][j][k]=-Out13[i][j][k]+2*tNS[16]*Concentration[i][j][k];
//                            cIn16[i][j][k]=Out13[i][j][k];
						if (Domain17[i][j][k]!=0)
							cIn17[i][j][k]=-Out14[i][j][k]+2*tNS[17]*Concentration[i][j][k];
//                            cIn17[i][j][k]=Out14[i][j][k];
						if (Domain18[i][j][k]!=0)
							cIn18[i][j][k]=-Out15[i][j][k]+2*tNS[18]*Concentration[i][j][k];
//                            cIn18[i][j][k]=Out15[i][j][k];
					}
				}
			}
		}
	}
}

void wallBoundaryConcentrationD3Q7(int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0){
						if (Domain1[i][j][k]!=0)
                            cIn1[i][j][k]=Out4[i][j][k];
						if (Domain2[i][j][k]!=0)
                            cIn2[i][j][k]=Out5[i][j][k];
						if (Domain3[i][j][k]!=0)
						    cIn3[i][j][k]=Out6[i][j][k];
						if (Domain4[i][j][k]!=0)
                            cIn4[i][j][k]=Out1[i][j][k];
						if (Domain5[i][j][k]!=0)
						    cIn5[i][j][k]=Out2[i][j][k];
						if (Domain6[i][j][k]!=0)
                            cIn6[i][j][k]=Out3[i][j][k];
					}
				}
			}
		}
	}
}

void wallBoundaryConcentrationD3Q15(int m, int n, int q, double gasCritical){
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
						if (Domain1[i][j][k] != 0)
							cIn1[i][j][k] = Out2[i][j][k];
						if (Domain2[i][j][k] != 0)
							cIn3[i][j][k] = Out4[i][j][k];
						if (Domain3[i][j][k] != 0)
							cIn5[i][j][k] = Out6[i][j][k];
						if (Domain4[i][j][k] != 0)
							cIn2[i][j][k] = Out1[i][j][k];
						if (Domain5[i][j][k] != 0)
							cIn4[i][j][k] = Out3[i][j][k];
						if (Domain6[i][j][k] != 0)
							cIn6[i][j][k] = Out5[i][j][k];
						if (Domain19[i][j][k] != 0)
							cIn7[i][j][k] = Out14[i][j][k];
						if (Domain20[i][j][k] != 0)
							cIn8[i][j][k] = Out13[i][j][k];
						if (Domain21[i][j][k] != 0)
							cIn9[i][j][k] = Out12[i][j][k];
						if (Domain22[i][j][k] != 0)
							cIn10[i][j][k] = Out11[i][j][k];
						if (Domain23[i][j][k] != 0)
							cIn11[i][j][k] = Out10[i][j][k];
						if (Domain24[i][j][k] != 0)
							cIn12[i][j][k] = Out9[i][j][k];
						if (Domain25[i][j][k] != 0)
							cIn13[i][j][k] = Out8[i][j][k];
						if (Domain26[i][j][k] != 0)
							cIn14[i][j][k] = Out7[i][j][k];
					}
				}
			}
		}
	}
}

void wallBoundaryConcentrationFirst(int m, int n, int q, double gasCritical){
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
						if (Domain1[i][j][k]!=0)
							cIn1First[i][j][k]=-Out4[i][j][k]+2*tNS[1]*ConcentrationFirst[i][j][k];
						if (Domain2[i][j][k]!=0)
							cIn2First[i][j][k]=-Out5[i][j][k]+2*tNS[2]*ConcentrationFirst[i][j][k];
						if (Domain3[i][j][k]!=0)
							cIn3First[i][j][k]=-Out6[i][j][k]+2*tNS[3]*ConcentrationFirst[i][j][k];
						if (Domain4[i][j][k]!=0)
							cIn4First[i][j][k]=-Out1[i][j][k]+2*tNS[4]*ConcentrationFirst[i][j][k];
						if (Domain5[i][j][k]!=0)
							cIn5First[i][j][k]=-Out2[i][j][k]+2*tNS[5]*ConcentrationFirst[i][j][k];
						if (Domain6[i][j][k]!=0)
							cIn6First[i][j][k]=-Out3[i][j][k]+2*tNS[6]*ConcentrationFirst[i][j][k];
						if (Domain7[i][j][k]!=0)
							cIn7First[i][j][k]=-Out10[i][j][k]+2*tNS[7]*ConcentrationFirst[i][j][k];
						if (Domain8[i][j][k]!=0)
							cIn8First[i][j][k]=-Out11[i][j][k]+2*tNS[8]*ConcentrationFirst[i][j][k];
						if (Domain9[i][j][k]!=0)
							cIn9First[i][j][k]=-Out12[i][j][k]+2*tNS[9]*ConcentrationFirst[i][j][k];
						if (Domain10[i][j][k]!=0)
							cIn10First[i][j][k]=-Out7[i][j][k]+2*tNS[10]*ConcentrationFirst[i][j][k];
						if (Domain11[i][j][k]!=0)
							cIn11First[i][j][k]=-Out8[i][j][k]+2*tNS[11]*ConcentrationFirst[i][j][k];
						if (Domain12[i][j][k]!=0)
							cIn12First[i][j][k]=-Out9[i][j][k]+2*tNS[12]*ConcentrationFirst[i][j][k];
						if (Domain13[i][j][k]!=0)
							cIn13First[i][j][k]=-Out16[i][j][k]+2*tNS[13]*ConcentrationFirst[i][j][k];
						if (Domain14[i][j][k]!=0)
							cIn14First[i][j][k]=-Out17[i][j][k]+2*tNS[14]*ConcentrationFirst[i][j][k];
						if (Domain15[i][j][k]!=0)
							cIn15First[i][j][k]=-Out18[i][j][k]+2*tNS[15]*ConcentrationFirst[i][j][k];
						if (Domain16[i][j][k]!=0)
							cIn16First[i][j][k]=-Out13[i][j][k]+2*tNS[16]*ConcentrationFirst[i][j][k];
						if (Domain17[i][j][k]!=0)
							cIn17First[i][j][k]=-Out14[i][j][k]+2*tNS[17]*ConcentrationFirst[i][j][k];
						if (Domain18[i][j][k]!=0)
							cIn18First[i][j][k]=-Out15[i][j][k]+2*tNS[18]*ConcentrationFirst[i][j][k];
					}
				}
			}
		}
	}
}

void wallBoundaryConcentrationElectrodeD3Q7(int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k] == 1 && electrodeConnection[i][j][k] == 1){
						if (Domain1[i][j][k]!=1)
//							cIn1First[i][j][k]=-Out4[i][j][k]+2*tNSD3Q7[1]*ConcentrationFirst[i][j][k];
                            cIn1First[i][j][k]=Out4[i][j][k];
						if (Domain2[i][j][k]!=1)
//							cIn2First[i][j][k]=-Out5[i][j][k]+2*tNSD3Q7[2]*ConcentrationFirst[i][j][k];
                            cIn2First[i][j][k]=Out5[i][j][k];
						if (Domain3[i][j][k]!=1)
//							cIn3First[i][j][k]=-Out6[i][j][k]+2*tNSD3Q7[3]*ConcentrationFirst[i][j][k];
                            cIn3First[i][j][k]=Out6[i][j][k];
						if (Domain4[i][j][k]!=1)
//							cIn4First[i][j][k]=-Out1[i][j][k]+2*tNSD3Q7[4]*ConcentrationFirst[i][j][k];
                            cIn4First[i][j][k]=Out1[i][j][k];
						if (Domain5[i][j][k]!=1)
//							cIn5First[i][j][k]=-Out2[i][j][k]+2*tNSD3Q7[5]*ConcentrationFirst[i][j][k];
                            cIn5First[i][j][k]=Out2[i][j][k];
						if (Domain6[i][j][k]!=1)
//							cIn6First[i][j][k]=-Out3[i][j][k]+2*tNSD3Q7[6]*ConcentrationFirst[i][j][k];
                            cIn6First[i][j][k]=Out3[i][j][k];
					}
				}
			}
		}
	}
}

void wallBoundaryConcentrationSecond(int m, int n, int q, double gasCritical){
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
						if (Domain1[i][j][k]!=0)
							cIn1Second[i][j][k]=-Out4[i][j][k]+2*tNS[1]*ConcentrationSecond[i][j][k];
						if (Domain2[i][j][k]!=0)
							cIn2Second[i][j][k]=-Out5[i][j][k]+2*tNS[2]*ConcentrationSecond[i][j][k];
						if (Domain3[i][j][k]!=0)
							cIn3Second[i][j][k]=-Out6[i][j][k]+2*tNS[3]*ConcentrationSecond[i][j][k];
						if (Domain4[i][j][k]!=0)
							cIn4Second[i][j][k]=-Out1[i][j][k]+2*tNS[4]*ConcentrationSecond[i][j][k];
						if (Domain5[i][j][k]!=0)
							cIn5Second[i][j][k]=-Out2[i][j][k]+2*tNS[5]*ConcentrationSecond[i][j][k];
						if (Domain6[i][j][k]!=0)
							cIn6Second[i][j][k]=-Out3[i][j][k]+2*tNS[6]*ConcentrationSecond[i][j][k];
						if (Domain7[i][j][k]!=0)
							cIn7Second[i][j][k]=-Out10[i][j][k]+2*tNS[7]*ConcentrationSecond[i][j][k];
						if (Domain8[i][j][k]!=0)
							cIn8Second[i][j][k]=-Out11[i][j][k]+2*tNS[8]*ConcentrationSecond[i][j][k];
						if (Domain9[i][j][k]!=0)
							cIn9Second[i][j][k]=-Out12[i][j][k]+2*tNS[9]*ConcentrationSecond[i][j][k];
						if (Domain10[i][j][k]!=0)
							cIn10Second[i][j][k]=-Out7[i][j][k]+2*tNS[10]*ConcentrationSecond[i][j][k];
						if (Domain11[i][j][k]!=0)
							cIn11Second[i][j][k]=-Out8[i][j][k]+2*tNS[11]*ConcentrationSecond[i][j][k];
						if (Domain12[i][j][k]!=0)
							cIn12Second[i][j][k]=-Out9[i][j][k]+2*tNS[12]*ConcentrationSecond[i][j][k];
						if (Domain13[i][j][k]!=0)
							cIn13Second[i][j][k]=-Out16[i][j][k]+2*tNS[13]*ConcentrationSecond[i][j][k];
						if (Domain14[i][j][k]!=0)
							cIn14Second[i][j][k]=-Out17[i][j][k]+2*tNS[14]*ConcentrationSecond[i][j][k];
						if (Domain15[i][j][k]!=0)
							cIn15Second[i][j][k]=-Out18[i][j][k]+2*tNS[15]*ConcentrationSecond[i][j][k];
						if (Domain16[i][j][k]!=0)
							cIn16Second[i][j][k]=-Out13[i][j][k]+2*tNS[16]*ConcentrationSecond[i][j][k];
						if (Domain17[i][j][k]!=0)
							cIn17Second[i][j][k]=-Out14[i][j][k]+2*tNS[17]*ConcentrationSecond[i][j][k];
						if (Domain18[i][j][k]!=0)
							cIn18Second[i][j][k]=-Out15[i][j][k]+2*tNS[18]*ConcentrationSecond[i][j][k];
					}
				}
			}
		}
	}
}

void wallBoundaryConcentrationSecondD3Q7(int m, int n, int q, double gasCritical){
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
						if (Domain1[i][j][k]!=0)
//							cIn1Second[i][j][k]=-Out4[i][j][k]+2*tNSD3Q7[1]*ConcentrationSecond[i][j][k];
                            cIn1Second[i][j][k]=Out4[i][j][k];
						if (Domain2[i][j][k]!=0)
//							cIn2Second[i][j][k]=-Out5[i][j][k]+2*tNSD3Q7[2]*ConcentrationSecond[i][j][k];
                            cIn2Second[i][j][k]=Out5[i][j][k];
						if (Domain3[i][j][k]!=0)
//							cIn3Second[i][j][k]=-Out6[i][j][k]+2*tNSD3Q7[3]*ConcentrationSecond[i][j][k];
                            cIn3Second[i][j][k]=Out6[i][j][k];
						if (Domain4[i][j][k]!=0)
//							cIn4Second[i][j][k]=-Out1[i][j][k]+2*tNSD3Q7[4]*ConcentrationSecond[i][j][k];
                            cIn4Second[i][j][k]=Out1[i][j][k];
						if (Domain5[i][j][k]!=0)
//							cIn5Second[i][j][k]=-Out2[i][j][k]+2*tNSD3Q7[5]*ConcentrationSecond[i][j][k];
                            cIn5Second[i][j][k]=Out2[i][j][k];
						if (Domain6[i][j][k]!=0)
//							cIn6Second[i][j][k]=-Out3[i][j][k]+2*tNSD3Q7[6]*ConcentrationSecond[i][j][k];
                            cIn6Second[i][j][k]=Out3[i][j][k];
					}
				}
			}
		}
	}
}

void wallBoundaryVelocityGhost(int m, int n, int q, double gasCritical){
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
						if (Domain1[i][j][k]!=0)
							fIn1Ghost[i][j][k]=Out4[i][j][k];
						if (Domain2[i][j][k]!=0)
							fIn2Ghost[i][j][k]=Out5[i][j][k];
						if (Domain3[i][j][k]!=0)
							fIn3Ghost[i][j][k]=Out6[i][j][k];
						if (Domain4[i][j][k]!=0)
							fIn4Ghost[i][j][k]=Out1[i][j][k];
						if (Domain5[i][j][k]!=0)
							fIn5Ghost[i][j][k]=Out2[i][j][k];
						if (Domain6[i][j][k]!=0)
							fIn6Ghost[i][j][k]=Out3[i][j][k];
						if (Domain7[i][j][k]!=0)
							fIn7Ghost[i][j][k]=Out10[i][j][k];
						if (Domain8[i][j][k]!=0)
							fIn8Ghost[i][j][k]=Out11[i][j][k];
						if (Domain9[i][j][k]!=0)
							fIn9Ghost[i][j][k]=Out12[i][j][k];
						if (Domain10[i][j][k]!=0)
							fIn10Ghost[i][j][k]=Out7[i][j][k];
						if (Domain11[i][j][k]!=0)
							fIn11Ghost[i][j][k]=Out8[i][j][k];
						if (Domain12[i][j][k]!=0)
							fIn12Ghost[i][j][k]=Out9[i][j][k];
						if (Domain13[i][j][k]!=0)
							fIn13Ghost[i][j][k]=Out16[i][j][k];
						if (Domain14[i][j][k]!=0)
							fIn14Ghost[i][j][k]=Out17[i][j][k];
						if (Domain15[i][j][k]!=0)
							fIn15Ghost[i][j][k]=Out18[i][j][k];
						if (Domain16[i][j][k]!=0)
							fIn16Ghost[i][j][k]=Out13[i][j][k];
						if (Domain17[i][j][k]!=0)
							fIn17Ghost[i][j][k]=Out14[i][j][k];
						if (Domain18[i][j][k]!=0)
							fIn18Ghost[i][j][k]=Out15[i][j][k];
					}
				}
			}
		}
	}
}

void wallSurfaceReactionBoundaryConcentration(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double RateConstant,double equilibriumPotential,double K,double exchange){
	int i;
	int j;
	int k;
	double rateConstant;
	double electrodePotential;
	double overpotential;
#pragma omp parallel private(i,j,k,rateConstant,electrodePotential,overpotential)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						if (Domain1[i][j][k]!=0){
							electrodePotential=neightbour4(i, j, k, m, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn1[i][j][k]=-Out4[i][j][k]+2*tNS[1]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain2[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour5(i, j, k, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn2[i][j][k]=-Out5[i][j][k]+2*tNS[2]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain3[i][j][k]!=0){
							electrodePotential=neightbour6(i, j, k, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn3[i][j][k]=-Out6[i][j][k]+2*tNS[3]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain4[i][j][k]!=0){
							electrodePotential=neightbour1(i, j, k, m, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn4[i][j][k]=-Out1[i][j][k]+2*tNS[4]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain5[i][j][k]!=0){
							electrodePotential=neightbour2(i, j, k, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn5[i][j][k]=-Out2[i][j][k]+2*tNS[5]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain6[i][j][k]!=0){
							electrodePotential=neightbour3(i, j, k, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn6[i][j][k]=-Out3[i][j][k]+2*tNS[6]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain7[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour10(i, j, k, m, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn7[i][j][k]=-Out10[i][j][k]+2*tNS[7]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain8[i][j][k]!=0){
							electrodePotential=neightbour11(i, j, k, m, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn8[i][j][k]=-Out11[i][j][k]+2*tNS[8]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain9[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour12(i, j, k, n, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn9[i][j][k]=-Out12[i][j][k]+2*tNS[9]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain10[i][j][k]!=0){
							electrodePotential=neightbour7(i, j, k, m, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn10[i][j][k]=-Out7[i][j][k]+2*tNS[10]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain11[i][j][k]!=0){
							electrodePotential=neightbour8(i, j, k, m, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn11[i][j][k]=-Out8[i][j][k]+2*tNS[11]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain12[i][j][k]!=0){
							electrodePotential=neightbour9(i, j, k, n, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn12[i][j][k]=-Out9[i][j][k]+2*tNS[12]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain13[i][j][k]!=0){
							electrodePotential=neightbour16(i, j, k, m, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn13[i][j][k]=-Out16[i][j][k]+2*tNS[13]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain14[i][j][k]!=0){
							electrodePotential=neightbour17(i, j, k, m, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn14[i][j][k]=-Out17[i][j][k]+2*tNS[14]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain15[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour18(i, j, k, n, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn15[i][j][k]=-Out18[i][j][k]+2*tNS[15]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain16[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour13(i, j, k, m, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn16[i][j][k]=-Out13[i][j][k]+2*tNS[16]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain17[i][j][k]!=0){
							electrodePotential=neightbour14(i, j, k, m, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn17[i][j][k]=-Out14[i][j][k]+2*tNS[17]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain18[i][j][k]!=0){
							electrodePotential=neightbour15(i, j, k, n, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn18[i][j][k]=-Out15[i][j][k]+2*tNS[18]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
					}
				}
			}
		}
	}
}

void wallSurfaceReactionBoundaryConcentrationD3Q7(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double RateConstant,double equilibriumPotential,double K,double exchange){
	int i;
	int j;
	int k;
	double rateConstant;
	double electrodePotential;
	double overpotential;
#pragma omp parallel private(i,j,k,rateConstant,electrodePotential,overpotential)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						if (Domain1[i][j][k]!=0){
							electrodePotential=neightbour4(i, j, k, m, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn1[i][j][k]=-Out4[i][j][k]+2*tNSD3Q7[1]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain2[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour5(i, j, k, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn2[i][j][k]=-Out5[i][j][k]+2*tNSD3Q7[2]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain3[i][j][k]!=0){
							electrodePotential=neightbour6(i, j, k, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn3[i][j][k]=-Out6[i][j][k]+2*tNSD3Q7[3]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain4[i][j][k]!=0){
							electrodePotential=neightbour1(i, j, k, m, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn4[i][j][k]=-Out1[i][j][k]+2*tNSD3Q7[4]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain5[i][j][k]!=0){
							electrodePotential=neightbour2(i, j, k, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn5[i][j][k]=-Out2[i][j][k]+2*tNSD3Q7[5]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain6[i][j][k]!=0){
							electrodePotential=neightbour3(i, j, k, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn6[i][j][k]=-Out3[i][j][k]+2*tNSD3Q7[6]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
					}
				}
			}
		}
	}
}

void wallSurfaceReactionBoundaryConcentrationFirst(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double RateConstant,double equilibriumPotential,double K,double exchange){
	int i;
	int j;
	int k;
	double rateConstant;
	double electrodePotential;
	double overpotential;
#pragma omp parallel private(i,j,k,rateConstant,electrodePotential,overpotential)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						if (Domain1[i][j][k]!=0){
							electrodePotential=neightbour4(i, j, k, m, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn1First[i][j][k]=-Out4[i][j][k]+2*tNS[1]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain2[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour5(i, j, k, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn2First[i][j][k]=-Out5[i][j][k]+2*tNS[2]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain3[i][j][k]!=0){
							electrodePotential=neightbour6(i, j, k, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn3First[i][j][k]=-Out6[i][j][k]+2*tNS[3]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain4[i][j][k]!=0){
							electrodePotential=neightbour1(i, j, k, m, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn4First[i][j][k]=-Out1[i][j][k]+2*tNS[4]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain5[i][j][k]!=0){
							electrodePotential=neightbour2(i, j, k, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn5First[i][j][k]=-Out2[i][j][k]+2*tNS[5]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain6[i][j][k]!=0){
							electrodePotential=neightbour3(i, j, k, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn6First[i][j][k]=-Out3[i][j][k]+2*tNS[6]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain7[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour10(i, j, k, m, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn7First[i][j][k]=-Out10[i][j][k]+2*tNS[7]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain8[i][j][k]!=0){
							electrodePotential=neightbour11(i, j, k, m, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn8First[i][j][k]=-Out11[i][j][k]+2*tNS[8]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain9[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour12(i, j, k, n, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn9First[i][j][k]=-Out12[i][j][k]+2*tNS[9]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain10[i][j][k]!=0){
							electrodePotential=neightbour7(i, j, k, m, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn10First[i][j][k]=-Out7[i][j][k]+2*tNS[10]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain11[i][j][k]!=0){
							electrodePotential=neightbour8(i, j, k, m, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn11First[i][j][k]=-Out8[i][j][k]+2*tNS[11]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain12[i][j][k]!=0){
							electrodePotential=neightbour9(i, j, k, n, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn12First[i][j][k]=-Out9[i][j][k]+2*tNS[12]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain13[i][j][k]!=0){
							electrodePotential=neightbour16(i, j, k, m, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn13First[i][j][k]=-Out16[i][j][k]+2*tNS[13]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain14[i][j][k]!=0){
							electrodePotential=neightbour17(i, j, k, m, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn14First[i][j][k]=-Out17[i][j][k]+2*tNS[14]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain15[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour18(i, j, k, n, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn15First[i][j][k]=-Out18[i][j][k]+2*tNS[15]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain16[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour13(i, j, k, m, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn16First[i][j][k]=-Out13[i][j][k]+2*tNS[16]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain17[i][j][k]!=0){
							electrodePotential=neightbour14(i, j, k, m, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn17First[i][j][k]=-Out14[i][j][k]+2*tNS[17]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain18[i][j][k]!=0){
							electrodePotential=neightbour15(i, j, k, n, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn18First[i][j][k]=-Out15[i][j][k]+2*tNS[18]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
					}
				}
			}
		}
	}
}

void wallSurfaceReactionBoundaryConcentrationFirstD3Q7(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double RateConstant,double equilibriumPotential,double K,double exchange){
	int i;
	int j;
	int k;
	double rateConstant;
	double electrodePotential;
	double overpotential;
#pragma omp parallel private(i,j,k,rateConstant,electrodePotential,overpotential)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						if (Domain1[i][j][k]!=0){
							electrodePotential=neightbour4(i, j, k, m, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn1First[i][j][k]=-Out4[i][j][k]+2*tNSD3Q7[1]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain2[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour5(i, j, k, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn2First[i][j][k]=-Out5[i][j][k]+2*tNSD3Q7[2]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain3[i][j][k]!=0){
							electrodePotential=neightbour6(i, j, k, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn3First[i][j][k]=-Out6[i][j][k]+2*tNSD3Q7[3]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain4[i][j][k]!=0){
							electrodePotential=neightbour1(i, j, k, m, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn4First[i][j][k]=-Out1[i][j][k]+2*tNSD3Q7[4]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain5[i][j][k]!=0){
							electrodePotential=neightbour2(i, j, k, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn5First[i][j][k]=-Out2[i][j][k]+2*tNSD3Q7[5]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
						if (Domain6[i][j][k]!=0){
							electrodePotential=neightbour3(i, j, k, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn6First[i][j][k]=-Out3[i][j][k]+2*tNSD3Q7[6]*0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant);
						}
					}
				}
			}
		}
	}
}

void wallSurfaceReactionBoundaryConcentrationSecond(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double diffuseCoefficiencySecond,double RateConstant,double equilibriumPotential,double K,double exchange){
	int i;
	int j;
	int k;
	double rateConstant;
	double electrodePotential;
	double overpotential;
#pragma omp parallel private(i,j,k,rateConstant,electrodePotential,overpotential)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						if (Domain1[i][j][k]!=0){
							electrodePotential=neightbour4(i, j, k, m, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn1Second[i][j][k]=-Out4[i][j][k]+2*tNS[1]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain2[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour5(i, j, k, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn2Second[i][j][k]=-Out5[i][j][k]+2*tNS[2]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain3[i][j][k]!=0){
							electrodePotential=neightbour6(i, j, k, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn3Second[i][j][k]=-Out6[i][j][k]+2*tNS[3]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain4[i][j][k]!=0){
							electrodePotential=neightbour1(i, j, k, m, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn4Second[i][j][k]=-Out1[i][j][k]+2*tNS[4]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain5[i][j][k]!=0){
							electrodePotential=neightbour2(i, j, k, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn5Second[i][j][k]=-Out2[i][j][k]+2*tNS[5]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain6[i][j][k]!=0){
							electrodePotential=neightbour3(i, j, k, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn6Second[i][j][k]=-Out3[i][j][k]+2*tNS[6]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain7[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour10(i, j, k, m, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn7Second[i][j][k]=-Out10[i][j][k]+2*tNS[7]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain8[i][j][k]!=0){
							electrodePotential=neightbour11(i, j, k, m, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn8Second[i][j][k]=-Out11[i][j][k]+2*tNS[8]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain9[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour12(i, j, k, n, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn9Second[i][j][k]=-Out12[i][j][k]+2*tNS[9]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain10[i][j][k]!=0){
							electrodePotential=neightbour7(i, j, k, m, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn10Second[i][j][k]=-Out7[i][j][k]+2*tNS[10]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain11[i][j][k]!=0){
							electrodePotential=neightbour8(i, j, k, m, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn11Second[i][j][k]=-Out8[i][j][k]+2*tNS[11]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain12[i][j][k]!=0){
							electrodePotential=neightbour9(i, j, k, n, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn12Second[i][j][k]=-Out9[i][j][k]+2*tNS[12]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain13[i][j][k]!=0){
							electrodePotential=neightbour16(i, j, k, m, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn13Second[i][j][k]=-Out16[i][j][k]+2*tNS[13]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain14[i][j][k]!=0){
							electrodePotential=neightbour17(i, j, k, m, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn14Second[i][j][k]=-Out17[i][j][k]+2*tNS[14]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain15[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour18(i, j, k, n, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn15Second[i][j][k]=-Out18[i][j][k]+2*tNS[15]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain16[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour13(i, j, k, m, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn16Second[i][j][k]=-Out13[i][j][k]+2*tNS[16]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain17[i][j][k]!=0){
							electrodePotential=neightbour14(i, j, k, m, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn17Second[i][j][k]=-Out14[i][j][k]+2*tNS[17]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain18[i][j][k]!=0){
							electrodePotential=neightbour15(i, j, k, n, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn18Second[i][j][k]=-Out15[i][j][k]+2*tNS[18]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
					}
				}
			}
		}
	}
}

void wallSurfaceReactionBoundaryConcentrationSecondD3Q7(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double diffuseCoefficiencySecond,double RateConstant,double equilibriumPotential,double K,double exchange){
	int i;
	int j;
	int k;
	double rateConstant;
	double electrodePotential;
	double overpotential;
#pragma omp parallel private(i,j,k,rateConstant,electrodePotential,overpotential)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						if (Domain1[i][j][k]!=0){
							electrodePotential=neightbour4(i, j, k, m, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn1Second[i][j][k]=-Out4[i][j][k]+2*tNSD3Q7[1]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain2[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour5(i, j, k, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn2Second[i][j][k]=-Out5[i][j][k]+2*tNSD3Q7[2]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain3[i][j][k]!=0){
							electrodePotential=neightbour6(i, j, k, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn3Second[i][j][k]=-Out6[i][j][k]+2*tNSD3Q7[3]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain4[i][j][k]!=0){
							electrodePotential=neightbour1(i, j, k, m, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn4Second[i][j][k]=-Out1[i][j][k]+2*tNSD3Q7[4]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain5[i][j][k]!=0){
							electrodePotential=neightbour2(i, j, k, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn5Second[i][j][k]=-Out2[i][j][k]+2*tNSD3Q7[5]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain6[i][j][k]!=0){
							electrodePotential=neightbour3(i, j, k, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							cIn6Second[i][j][k]=-Out3[i][j][k]+2*tNSD3Q7[6]*(ConcentrationSecond[i][j][k]-0.25*rateConstant/diffuseCoefficiencySecond*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
					}
				}
			}
		}
	}
}

void wallBoundaryElectrodePotential(int m, int n, int q){
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
						if (Domain1[i][j][k]!=1)
							electrodeIn1[i][j][k]=-Out4[i][j][k]+2*tNS[1]*electrodePotential[i][j][k];
						if (Domain2[i][j][k]!=1)
							electrodeIn2[i][j][k]=-Out5[i][j][k]+2*tNS[2]*electrodePotential[i][j][k];
						if (Domain3[i][j][k]!=1)
							electrodeIn3[i][j][k]=-Out6[i][j][k]+2*tNS[3]*electrodePotential[i][j][k];
						if (Domain4[i][j][k]!=1)
							electrodeIn4[i][j][k]=-Out1[i][j][k]+2*tNS[4]*electrodePotential[i][j][k];
						if (Domain5[i][j][k]!=1)
							electrodeIn5[i][j][k]=-Out2[i][j][k]+2*tNS[5]*electrodePotential[i][j][k];
						if (Domain6[i][j][k]!=1)
							electrodeIn6[i][j][k]=-Out3[i][j][k]+2*tNS[6]*electrodePotential[i][j][k];
						if (Domain7[i][j][k]!=1)
							electrodeIn7[i][j][k]=-Out10[i][j][k]+2*tNS[7]*electrodePotential[i][j][k];
						if (Domain8[i][j][k]!=1)
							electrodeIn8[i][j][k]=-Out11[i][j][k]+2*tNS[8]*electrodePotential[i][j][k];
						if (Domain9[i][j][k]!=1)
							electrodeIn9[i][j][k]=-Out12[i][j][k]+2*tNS[9]*electrodePotential[i][j][k];
						if (Domain10[i][j][k]!=1)
							electrodeIn10[i][j][k]=-Out7[i][j][k]+2*tNS[10]*electrodePotential[i][j][k];
						if (Domain11[i][j][k]!=1)
							electrodeIn11[i][j][k]=-Out8[i][j][k]+2*tNS[11]*electrodePotential[i][j][k];
						if (Domain12[i][j][k]!=1)
							electrodeIn12[i][j][k]=-Out9[i][j][k]+2*tNS[12]*electrodePotential[i][j][k];
						if (Domain13[i][j][k]!=1)
							electrodeIn13[i][j][k]=-Out16[i][j][k]+2*tNS[13]*electrodePotential[i][j][k];
						if (Domain14[i][j][k]!=1)
							electrodeIn14[i][j][k]=-Out17[i][j][k]+2*tNS[14]*electrodePotential[i][j][k];
						if (Domain15[i][j][k]!=1)
							electrodeIn15[i][j][k]=-Out18[i][j][k]+2*tNS[15]*electrodePotential[i][j][k];
						if (Domain16[i][j][k]!=1)
							electrodeIn16[i][j][k]=-Out13[i][j][k]+2*tNS[16]*electrodePotential[i][j][k];
						if (Domain17[i][j][k]!=1)
							electrodeIn17[i][j][k]=-Out14[i][j][k]+2*tNS[17]*electrodePotential[i][j][k];
						if (Domain18[i][j][k]!=1)
							electrodeIn18[i][j][k]=-Out15[i][j][k]+2*tNS[18]*electrodePotential[i][j][k];
					}
				}
			}
		}
	}
}

void wallBoundaryElectrodePotentialD3Q7(int m, int n, int q){
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
						if (Domain1[i][j][k] == 0 || Domain1[i][j][k] == 3)
                            electrodeIn1[i][j][k]=Out4[i][j][k];
						if (Domain2[i][j][k] == 0 || Domain2[i][j][k] == 3)
                            electrodeIn2[i][j][k]=Out5[i][j][k];
						if (Domain3[i][j][k] == 0 || Domain3[i][j][k] == 3)
                            electrodeIn3[i][j][k]=Out6[i][j][k];
						if (Domain4[i][j][k] == 0 || Domain4[i][j][k] == 3)
                            electrodeIn4[i][j][k]=Out1[i][j][k];
						if (Domain5[i][j][k] == 0 || Domain5[i][j][k] == 3)
                            electrodeIn5[i][j][k]=Out2[i][j][k];
						if (Domain6[i][j][k] == 0 || Domain6[i][j][k] == 3)
                            electrodeIn6[i][j][k]=Out3[i][j][k];
					}
				}
			}
		}
	}
}

void wallSurfaceReactionBoundaryElectrode(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double electrodeConductivity,double RateConstant,double constant,double equilibriumPotential,double K,double exchange){
	int i;
	int j;
	int k;
	double valueConcentration;
	double valueConcentrationFirst;
	double rateConstant;
	double electrolytePotential;
	double overpotential;
#pragma omp parallel private(i,j,k,valueConcentration,valueConcentrationFirst,rateConstant,electrolytePotential,overpotential)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==1&&electrodeConnection[i][j][k]==1){
						if (Domain1[i][j][k]!=1&&rho1[i][j][j]>gasCritical){
							valueConcentration=neightbour4(i, j, k, m, wallNearConcentration);
							valueConcentrationFirst=neightbour4(i,j,k,m,wallNearConcentrationFirst);
							electrolytePotential=neightbour4(i, j, k, m, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn1[i][j][k]=-Out4[i][j][k]+2*tNS[1]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain2[i][j][k]!=1&&rho2[i][j][j]>gasCritical&&j!=0){
							valueConcentration=neightbour5(i,j,k,n,wallNearConcentration);
							valueConcentrationFirst=neightbour5(i,j,k,n,wallNearConcentrationFirst);
							electrolytePotential=neightbour5(i, j, k, n, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn2[i][j][k]=-Out5[i][j][k]+2*tNS[2]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain3[i][j][k]!=1&&rho3[i][j][j]>gasCritical){
							valueConcentration=neightbour6(i, j, k, q, wallNearConcentration);
							valueConcentrationFirst=neightbour6(i, j, k, q, wallNearConcentrationFirst);
							electrolytePotential=neightbour6(i, j, k, q, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn3[i][j][k]=-Out6[i][j][k]+2*tNS[3]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain4[i][j][k]!=1&&rho4[i][j][j]>gasCritical){
							valueConcentration=neightbour1(i, j, k, m, wallNearConcentration);
                            valueConcentrationFirst=neightbour1(i, j, k, m, wallNearConcentrationFirst);
							electrolytePotential=neightbour1(i, j, k, m, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn4[i][j][k]=-Out1[i][j][k]+2*tNS[4]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain5[i][j][k]!=1&&rho5[i][j][j]>gasCritical){
							valueConcentration=neightbour2(i, j, k, n, wallNearConcentration);
                            valueConcentrationFirst=neightbour2(i, j, k, n, wallNearConcentrationFirst);
							electrolytePotential=neightbour2(i, j, k, n, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn5[i][j][k]=-Out2[i][j][k]+2*tNS[5]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain6[i][j][k]!=1&&rho6[i][j][j]>gasCritical){
							valueConcentration=neightbour3(i, j, k, q, wallNearConcentration);
                            valueConcentrationFirst=neightbour3(i, j, k, q, wallNearConcentrationFirst);
							electrolytePotential=neightbour3(i, j, k, q, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn6[i][j][k]=-Out3[i][j][k]+2*tNS[6]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain7[i][j][k]!=1&&rho7[i][j][j]>gasCritical&&j!=0){
							valueConcentration=neightbour10(i, j, k, m, n, wallNearConcentration);
                            valueConcentrationFirst=neightbour10(i, j, k, m, n, wallNearConcentrationFirst);
							electrolytePotential=neightbour10(i, j, k, m, n, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn7[i][j][k]=-Out10[i][j][k]+2*tNS[7]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain8[i][j][k]!=1&&rho8[i][j][j]>gasCritical){
							valueConcentration=neightbour11(i, j, k, m, q, wallNearConcentration);
                            valueConcentrationFirst=neightbour11(i, j, k, m, q, wallNearConcentrationFirst);
							electrolytePotential=neightbour11(i, j, k, m, q, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn8[i][j][k]=-Out11[i][j][k]+2*tNS[8]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain9[i][j][k]!=1&&rho9[i][j][j]>gasCritical&&j!=0){
							valueConcentration=neightbour12(i, j, k, n, q, wallNearConcentration);
                            valueConcentrationFirst=neightbour12(i, j, k, n, q, wallNearConcentrationFirst);
							electrolytePotential=neightbour12(i, j, k, n, q, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn9[i][j][k]=-Out12[i][j][k]+2*tNS[9]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain10[i][j][k]!=1&&rho10[i][j][j]>gasCritical){
							valueConcentration=neightbour7(i, j, k, m, n, wallNearConcentration);
                            valueConcentrationFirst=neightbour7(i, j, k, m, n, wallNearConcentrationFirst);
							electrolytePotential=neightbour7(i, j, k, m, n, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn10[i][j][k]=-Out7[i][j][k]+2*tNS[10]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain11[i][j][k]!=1&&rho11[i][j][j]>gasCritical){
							valueConcentration=neightbour8(i, j, k, m, q, wallNearConcentration);
                            valueConcentrationFirst=neightbour8(i, j, k, m, q, wallNearConcentrationFirst);
							electrolytePotential=neightbour8(i, j, k, m, q, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn11[i][j][k]=-Out8[i][j][k]+2*tNS[11]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain12[i][j][k]!=1&&rho12[i][j][j]>gasCritical){
							valueConcentration=neightbour9(i, j, k, n, q, wallNearConcentration);
                            valueConcentrationFirst=neightbour9(i, j, k, n, q, wallNearConcentrationFirst);
							electrolytePotential=neightbour9(i, j, k, n, q, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn12[i][j][k]=-Out9[i][j][k]+2*tNS[12]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain13[i][j][k]!=1&&rho13[i][j][j]>gasCritical){
							valueConcentration=neightbour16(i, j, k, m, n, wallNearConcentration);
                            valueConcentrationFirst=neightbour16(i, j, k, m, n, wallNearConcentrationFirst);
							electrolytePotential=neightbour16(i, j, k, m, n, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn13[i][j][k]=-Out16[i][j][k]+2*tNS[13]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain14[i][j][k]!=1&&rho14[i][j][j]>gasCritical){
							valueConcentration=neightbour17(i, j, k, m, q, wallNearConcentration);
                            valueConcentrationFirst=neightbour17(i, j, k, m, q, wallNearConcentrationFirst);
							electrolytePotential=neightbour17(i, j, k, m, q, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn14[i][j][k]=-Out17[i][j][k]+2*tNS[14]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain15[i][j][k]!=1&&rho15[i][j][j]>gasCritical&&j!=0){
							valueConcentration=neightbour18(i, j, k, n, q, wallNearConcentration);
                            valueConcentrationFirst=neightbour18(i, j, k, n, q, wallNearConcentrationFirst);
							electrolytePotential=neightbour18(i, j, k, n, q, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn15[i][j][k]=-Out18[i][j][k]+2*tNS[15]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain16[i][j][k]!=1&&rho16[i][j][j]>gasCritical&&j!=0){
							valueConcentration=neightbour13(i, j, k, m, n, wallNearConcentration);
                            valueConcentrationFirst=neightbour13(i, j, k, m, n, wallNearConcentrationFirst);
							electrolytePotential=neightbour13(i, j, k, m, n, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn16[i][j][k]=-Out13[i][j][k]+2*tNS[16]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain17[i][j][k]!=1&&rho17[i][j][j]>gasCritical){
							valueConcentration=neightbour14(i, j, k, m, q, wallNearConcentration);
                            valueConcentrationFirst=neightbour14(i, j, k, m, q, wallNearConcentrationFirst);
							electrolytePotential=neightbour14(i, j, k, m, q, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn17[i][j][k]=-Out14[i][j][k]+2*tNS[17]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain18[i][j][k]!=1&&rho18[i][j][j]>gasCritical){
							valueConcentration=neightbour15(i, j, k, n, q, wallNearConcentration);
                            valueConcentrationFirst=neightbour15(i, j, k, n, q, wallNearConcentrationFirst);
							electrolytePotential=neightbour15(i, j, k, n, q, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn18[i][j][k]=-Out15[i][j][k]+2*tNS[18]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
					}
				}
			}
		}
	}
}

void wallSurfaceReactionBoundaryElectrodeD3Q7(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double electrodeConductivity,double RateConstant,double constant,double equilibriumPotential,double K,double exchange){
	int i;
	int j;
	int k;
	double valueConcentration;
	double valueConcentrationFirst;
	double rateConstant;
	double electrolytePotential;
	double overpotential;
#pragma omp parallel private(i,j,k,valueConcentration,valueConcentrationFirst,rateConstant,electrolytePotential,overpotential)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==1&&electrodeConnection[i][j][k]==1){
						if (Domain1[i][j][k]!=1&&rho1[i][j][j]>gasCritical){
							valueConcentration=neightbour4(i, j, k, m, wallNearConcentration);
							valueConcentrationFirst=neightbour4(i,j,k,m,wallNearConcentrationFirst);
							electrolytePotential=neightbour4(i, j, k, m, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn1[i][j][k]=-Out4[i][j][k]+2*tNSD3Q7[1]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
/*						if (Domain2[i][j][k]!=1&&rho2[i][j][j]>gasCritical&&j!=0){
							valueConcentration=neightbour5(i,j,k,n,wallNearConcentration);
							valueConcentrationFirst=neightbour5(i,j,k,n,wallNearConcentrationFirst);
							electrolytePotential=neightbour5(i, j, k, n, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn2[i][j][k]=-Out5[i][j][k]+2*tNSD3Q7[2]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}*/
						if (Domain3[i][j][k]!=1&&rho3[i][j][j]>gasCritical){
							valueConcentration=neightbour6(i, j, k, q, wallNearConcentration);
							valueConcentrationFirst=neightbour6(i, j, k, q, wallNearConcentrationFirst);
							electrolytePotential=neightbour6(i, j, k, q, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn3[i][j][k]=-Out6[i][j][k]+2*tNSD3Q7[3]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain4[i][j][k]!=1&&rho4[i][j][j]>gasCritical){
							valueConcentration=neightbour1(i, j, k, m, wallNearConcentration);
                            valueConcentrationFirst=neightbour1(i, j, k, m, wallNearConcentrationFirst);
							electrolytePotential=neightbour1(i, j, k, m, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn4[i][j][k]=-Out1[i][j][k]+2*tNSD3Q7[4]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
/*						if (Domain5[i][j][k]!=1&&rho5[i][j][j]>gasCritical){
							valueConcentration=neightbour2(i, j, k, n, wallNearConcentration);
                            valueConcentrationFirst=neightbour2(i, j, k, n, wallNearConcentrationFirst);
							electrolytePotential=neightbour2(i, j, k, n, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn5[i][j][k]=-Out2[i][j][k]+2*tNSD3Q7[5]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}*/
						if (Domain6[i][j][k]!=1&&rho6[i][j][j]>gasCritical){
							valueConcentration=neightbour3(i, j, k, q, wallNearConcentration);
                            valueConcentrationFirst=neightbour3(i, j, k, q, wallNearConcentrationFirst);
							electrolytePotential=neightbour3(i, j, k, q, lastWallElectrolyte);
							overpotential=electrodePotential[i][j][k]/exchange-electrolytePotential/exchange-(equilibriumPotential+8.314*K/96485.0*log(valueConcentration/valueConcentrationFirst));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrodeIn6[i][j][k]=-Out3[i][j][k]+2*tNSD3Q7[6]*(electrodePotential[i][j][k]+0.5*constant*rateConstant/electrodeConductivity*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentration+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst-pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*valueConcentrationFirst+rateConstant*rateConstant*valueConcentration+rateConstant*rateConstant*valueConcentrationFirst+pow(rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentration+rateConstant*rateConstant*rateConstant*rateConstant*valueConcentrationFirst*valueConcentrationFirst+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst+2.0*rateConstant*rateConstant*rateConstant*rateConstant*valueConcentration*valueConcentrationFirst,0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
					}
				}
			}
		}
	}
}

void wallBoundaryElectrolytePotential(int m, int n, int q, double gasCritical){
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
						if (Domain1[i][j][k]!=0)
							electrolyteIn1[i][j][k]=-Out4[i][j][k]+2*tNS[1]*electrolytePotential[i][j][k];
						if (Domain2[i][j][k]!=0)
							electrolyteIn2[i][j][k]=-Out5[i][j][k]+2*tNS[2]*electrolytePotential[i][j][k];
						if (Domain3[i][j][k]!=0)
							electrolyteIn3[i][j][k]=-Out6[i][j][k]+2*tNS[3]*electrolytePotential[i][j][k];
						if (Domain4[i][j][k]!=0)
							electrolyteIn4[i][j][k]=-Out1[i][j][k]+2*tNS[4]*electrolytePotential[i][j][k];
						if (Domain5[i][j][k]!=0)
							electrolyteIn5[i][j][k]=-Out2[i][j][k]+2*tNS[5]*electrolytePotential[i][j][k];
						if (Domain6[i][j][k]!=0)
							electrolyteIn6[i][j][k]=-Out3[i][j][k]+2*tNS[6]*electrolytePotential[i][j][k];
						if (Domain7[i][j][k]!=0)
							electrolyteIn7[i][j][k]=-Out10[i][j][k]+2*tNS[7]*electrolytePotential[i][j][k];
						if (Domain8[i][j][k]!=0)
							electrolyteIn8[i][j][k]=-Out11[i][j][k]+2*tNS[8]*electrolytePotential[i][j][k];
						if (Domain9[i][j][k]!=0)
							electrolyteIn9[i][j][k]=-Out12[i][j][k]+2*tNS[9]*electrolytePotential[i][j][k];
						if (Domain10[i][j][k]!=0)
							electrolyteIn10[i][j][k]=-Out7[i][j][k]+2*tNS[10]*electrolytePotential[i][j][k];
						if (Domain11[i][j][k]!=0)
							electrolyteIn11[i][j][k]=-Out8[i][j][k]+2*tNS[11]*electrolytePotential[i][j][k];
						if (Domain12[i][j][k]!=0)
							electrolyteIn12[i][j][k]=-Out9[i][j][k]+2*tNS[12]*electrolytePotential[i][j][k];
						if (Domain13[i][j][k]!=0)
							electrolyteIn13[i][j][k]=-Out16[i][j][k]+2*tNS[13]*electrolytePotential[i][j][k];
						if (Domain14[i][j][k]!=0)
							electrolyteIn14[i][j][k]=-Out17[i][j][k]+2*tNS[14]*electrolytePotential[i][j][k];
						if (Domain15[i][j][k]!=0)
							electrolyteIn15[i][j][k]=-Out18[i][j][k]+2*tNS[15]*electrolytePotential[i][j][k];
						if (Domain16[i][j][k]!=0)
							electrolyteIn16[i][j][k]=-Out13[i][j][k]+2*tNS[16]*electrolytePotential[i][j][k];
						if (Domain17[i][j][k]!=0)
							electrolyteIn17[i][j][k]=-Out14[i][j][k]+2*tNS[17]*electrolytePotential[i][j][k];
						if (Domain18[i][j][k]!=0)
							electrolyteIn18[i][j][k]=-Out15[i][j][k]+2*tNS[18]*electrolytePotential[i][j][k];
					}
				}
			}
		}
	}
}

void wallBoundaryElectrolytePotentialD3Q7(int m, int n, int q, double gasCritical){
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
						if (Domain1[i][j][k]!=0)
                            electrolyteIn1[i][j][k]=Out4[i][j][k];
						if (Domain2[i][j][k]!=0)
                            electrolyteIn2[i][j][k]=Out5[i][j][k];
						if (Domain3[i][j][k]!=0)
                            electrolyteIn3[i][j][k]=Out6[i][j][k];
						if (Domain4[i][j][k]!=0)
                            electrolyteIn4[i][j][k]=Out1[i][j][k];
						if (Domain5[i][j][k]!=0)
                            electrolyteIn5[i][j][k]=Out2[i][j][k];
						if (Domain6[i][j][k]!=0)
                            electrolyteIn6[i][j][k]=Out3[i][j][k];
					}
				}
			}
		}
	}
}

void wallSurfaceReactionBoundaryElectrolyte(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double ***electrolyteRelaxationTime,double RateConstant,double constant,double equilibriumPotential,double K,double exchange){
	int i;
	int j;
	int k;
	double rateConstant;
	double electrodePotential;
	double overpotential;
#pragma omp parallel private(i,j,k,rateConstant,electrodePotential,overpotential)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						if (Domain1[i][j][k]!=0){
							electrodePotential=neightbour4(i, j, k, m, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn1[i][j][k]=-Out4[i][j][k]+2*tNS[1]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/3.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain2[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour5(i, j, k, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn2[i][j][k]=-Out5[i][j][k]+2*tNS[2]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/3.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain3[i][j][k]!=0){
							electrodePotential=neightbour6(i, j, k, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn3[i][j][k]=-Out6[i][j][k]+2*tNS[3]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/3.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain4[i][j][k]!=0){
							electrodePotential=neightbour1(i, j, k, m, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn4[i][j][k]=-Out1[i][j][k]+2*tNS[4]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/3.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain5[i][j][k]!=0){
							electrodePotential=neightbour2(i, j, k, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn5[i][j][k]=-Out2[i][j][k]+2*tNS[5]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/3.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain6[i][j][k]!=0){
							electrodePotential=neightbour3(i, j, k, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn6[i][j][k]=-Out3[i][j][k]+2*tNS[6]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/3.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain7[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour10(i, j, k, m, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn7[i][j][k]=-Out10[i][j][k]+2*tNS[7]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/3.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain8[i][j][k]!=0){
							electrodePotential=neightbour11(i, j, k, m, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn8[i][j][k]=-Out11[i][j][k]+2*tNS[8]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/3.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain9[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour12(i, j, k, n, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn9[i][j][k]=-Out12[i][j][k]+2*tNS[9]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/3.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain10[i][j][k]!=0){
							electrodePotential=neightbour7(i, j, k, m, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn10[i][j][k]=-Out7[i][j][k]+2*tNS[10]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/3.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain11[i][j][k]!=0){
							electrodePotential=neightbour8(i, j, k, m, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn11[i][j][k]=-Out8[i][j][k]+2*tNS[11]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/3.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain12[i][j][k]!=0){
							electrodePotential=neightbour9(i, j, k, n, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn12[i][j][k]=-Out9[i][j][k]+2*tNS[12]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/3.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain13[i][j][k]!=0){
							electrodePotential=neightbour16(i, j, k, m, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn13[i][j][k]=-Out16[i][j][k]+2*tNS[13]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/3.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain14[i][j][k]!=0){
							electrodePotential=neightbour17(i, j, k, m, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn14[i][j][k]=-Out17[i][j][k]+2*tNS[14]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/3.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain15[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour18(i, j, k, n, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn15[i][j][k]=-Out18[i][j][k]+2*tNS[15]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/3.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain16[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour13(i, j, k, m, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn16[i][j][k]=-Out13[i][j][k]+2*tNS[16]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/3.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain17[i][j][k]!=0){
							electrodePotential=neightbour14(i, j, k, m, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn17[i][j][k]=-Out14[i][j][k]+2*tNS[17]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/3.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain18[i][j][k]!=0){
							electrodePotential=neightbour15(i, j, k, n, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn18[i][j][k]=-Out15[i][j][k]+2*tNS[18]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/3.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
					}
				}
			}
		}
	}
}

void wallSurfaceReactionBoundaryElectrolyteD3Q7(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double ***electrolyteRelaxationTime,double RateConstant,double constant,double equilibriumPotential,double K,double exchange){
	int i;
	int j;
	int k;
	double rateConstant;
	double electrodePotential;
	double overpotential;
#pragma omp parallel private(i,j,k,rateConstant,electrodePotential,overpotential)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						if (Domain1[i][j][k]!=0){
							electrodePotential=neightbour4(i, j, k, m, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn1[i][j][k]=-Out4[i][j][k]+2*tNSD3Q7[1]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/4.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
/*						if (Domain2[i][j][k]!=0&&j!=0){
							electrodePotential=neightbour5(i, j, k, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn2[i][j][k]=-Out5[i][j][k]+2*tNSD3Q7[2]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/4.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}*/
						if (Domain3[i][j][k]!=0){
							electrodePotential=neightbour6(i, j, k, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn3[i][j][k]=-Out6[i][j][k]+2*tNSD3Q7[3]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/4.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
						if (Domain4[i][j][k]!=0){
							electrodePotential=neightbour1(i, j, k, m, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn4[i][j][k]=-Out1[i][j][k]+2*tNSD3Q7[4]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/4.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
/*						if (Domain5[i][j][k]!=0){
							electrodePotential=neightbour2(i, j, k, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn5[i][j][k]=-Out2[i][j][k]+2*tNSD3Q7[5]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/4.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}*/
						if (Domain6[i][j][k]!=0){
							electrodePotential=neightbour3(i, j, k, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=-RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							electrolyteIn6[i][j][k]=-Out3[i][j][k]+2*tNSD3Q7[6]*(electrolytePotential[i][j][k]-0.5*constant*rateConstant/((1.0/electrolyteRelaxationTime[i][j][k]-0.5)/4.0)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]-pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5)*pow(0.5*(8.0*diffuseCoefficiency*diffuseCoefficiency*wallNearConcentrationFirst[i][j][k]+rateConstant*rateConstant*wallNearConcentration[i][j][k]+rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]+pow(rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentration[i][j][k]+rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentrationFirst[i][j][k]*wallNearConcentrationFirst[i][j][k]+16.0*diffuseCoefficiency*diffuseCoefficiency*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k]+2.0*rateConstant*rateConstant*rateConstant*rateConstant*wallNearConcentration[i][j][k]*wallNearConcentrationFirst[i][j][k],0.5))/(4.0*diffuseCoefficiency*diffuseCoefficiency+rateConstant*rateConstant),0.5));
						}
					}
				}
			}
		}
	}
}