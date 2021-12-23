#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "domain.h"
#include "MemoryArrange.h"
#include "OutMemoryArrange.h"
#include "wallContactAngle.h"
#include "matrixMove.h"
#include "memoryFree.h"
#include "equilibrium.h"

#define viscosity 0.5/3
#define concentrationDiffusivity 0.00007
#define concentrationFirstDiffusivity 0.0013
#define concentrationSecondDiffusivity 0.031
#define concentrationThirdDiffusivity 0.00073333
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

#define fluidForceX1 Psi1*(-cxNS[1])
#define fluidForceY1 Psi1*(-cyNS[1])
#define fluidForceZ1 Psi1*(-czNS[1])
#define fluidForceX2 Psi2*(-cxNS[2])
#define fluidForceY2 Psi2*(-cyNS[2])
#define fluidForceZ2 Psi2*(-czNS[2])
#define fluidForceX3 Psi3*(-cxNS[3])
#define fluidForceY3 Psi3*(-cyNS[3])
#define fluidForceZ3 Psi3*(-czNS[3])
#define fluidForceX4 Psi4*(-cxNS[4])
#define fluidForceY4 Psi4*(-cyNS[4])
#define fluidForceZ4 Psi4*(-czNS[4])
#define fluidForceX5 Psi5*(-cxNS[5])
#define fluidForceY5 Psi5*(-cyNS[5])
#define fluidForceZ5 Psi5*(-czNS[5])
#define fluidForceX6 Psi6*(-cxNS[6])
#define fluidForceY6 Psi6*(-cyNS[6])
#define fluidForceZ6 Psi6*(-czNS[6])
#define fluidForceX7 Psi7*(-cxNS[7])
#define fluidForceY7 Psi7*(-cyNS[7])
#define fluidForceZ7 Psi7*(-czNS[7])
#define fluidForceX8 Psi8*(-cxNS[8])
#define fluidForceY8 Psi8*(-cyNS[8])
#define fluidForceZ8 Psi8*(-czNS[8])
#define fluidForceX9 Psi9*(-cxNS[9])
#define fluidForceY9 Psi9*(-cyNS[9])
#define fluidForceZ9 Psi9*(-czNS[9])
#define fluidForceX10 Psi10*(-cxNS[10])
#define fluidForceY10 Psi10*(-cyNS[10])
#define fluidForceZ10 Psi10*(-czNS[10])
#define fluidForceX11 Psi11*(-cxNS[11])
#define fluidForceY11 Psi11*(-cyNS[11])
#define fluidForceZ11 Psi11*(-czNS[11])
#define fluidForceX12 Psi12*(-cxNS[12])
#define fluidForceY12 Psi12*(-cyNS[12])
#define fluidForceZ12 Psi12*(-czNS[12])
#define fluidForceX13 Psi13*(-cxNS[13])
#define fluidForceY13 Psi13*(-cyNS[13])
#define fluidForceZ13 Psi13*(-czNS[13])
#define fluidForceX14 Psi14*(-cxNS[14])
#define fluidForceY14 Psi14*(-cyNS[14])
#define fluidForceZ14 Psi14*(-czNS[14])
#define fluidForceX15 Psi15*(-cxNS[15])
#define fluidForceY15 Psi15*(-cyNS[15])
#define fluidForceZ15 Psi15*(-czNS[15])
#define fluidForceX16 Psi16*(-cxNS[16])
#define fluidForceY16 Psi16*(-cyNS[16])
#define fluidForceZ16 Psi16*(-czNS[16])
#define fluidForceX17 Psi17*(-cxNS[17])
#define fluidForceY17 Psi17*(-cyNS[17])
#define fluidForceZ17 Psi17*(-czNS[17])
#define fluidForceX18 Psi18*(-cxNS[18])
#define fluidForceY18 Psi18*(-cyNS[18])
#define fluidForceZ18 Psi18*(-czNS[18])

#define Force0 ((1-0.5*A)*tNS[0]*(3.0*(cfNS0-uxFxuyFyuzFz)+3.0*cuNS0[i][j][k]*cfNS0))
#define Force1 ((1-0.5*A)*tNS[1]*(3.0*(cfNS1-uxFxuyFyuzFz)+3.0*cuNS1[i][j][k]*cfNS1))
#define Force2 ((1-0.5*A)*tNS[2]*(3.0*(cfNS2-uxFxuyFyuzFz)+3.0*cuNS2[i][j][k]*cfNS2))
#define Force3 ((1-0.5*A)*tNS[3]*(3.0*(cfNS3-uxFxuyFyuzFz)+3.0*cuNS3[i][j][k]*cfNS3))
#define Force4 ((1-0.5*A)*tNS[4]*(3.0*(cfNS4-uxFxuyFyuzFz)+3.0*cuNS4[i][j][k]*cfNS4))
#define Force5 ((1-0.5*A)*tNS[5]*(3.0*(cfNS5-uxFxuyFyuzFz)+3.0*cuNS5[i][j][k]*cfNS5))
#define Force6 ((1-0.5*A)*tNS[6]*(3.0*(cfNS6-uxFxuyFyuzFz)+3.0*cuNS6[i][j][k]*cfNS6))
#define Force7 ((1-0.5*A)*tNS[7]*(3.0*(cfNS7-uxFxuyFyuzFz)+3.0*cuNS7[i][j][k]*cfNS7))
#define Force8 ((1-0.5*A)*tNS[8]*(3.0*(cfNS8-uxFxuyFyuzFz)+3.0*cuNS8[i][j][k]*cfNS8))
#define Force9 ((1-0.5*A)*tNS[9]*(3.0*(cfNS9-uxFxuyFyuzFz)+3.0*cuNS9[i][j][k]*cfNS9))
#define Force10 ((1-0.5*A)*tNS[10]*(3.0*(cfNS10-uxFxuyFyuzFz)+3.0*cuNS10[i][j][k]*cfNS10))
#define Force11 ((1-0.5*A)*tNS[11]*(3.0*(cfNS11-uxFxuyFyuzFz)+3.0*cuNS11[i][j][k]*cfNS11))
#define Force12 ((1-0.5*A)*tNS[12]*(3.0*(cfNS12-uxFxuyFyuzFz)+3.0*cuNS12[i][j][k]*cfNS12))
#define Force13 ((1-0.5*A)*tNS[13]*(3.0*(cfNS13-uxFxuyFyuzFz)+3.0*cuNS13[i][j][k]*cfNS13))
#define Force14 ((1-0.5*A)*tNS[14]*(3.0*(cfNS14-uxFxuyFyuzFz)+3.0*cuNS14[i][j][k]*cfNS14))
#define Force15 ((1-0.5*A)*tNS[15]*(3.0*(cfNS15-uxFxuyFyuzFz)+3.0*cuNS15[i][j][k]*cfNS15))
#define Force16 ((1-0.5*A)*tNS[16]*(3.0*(cfNS16-uxFxuyFyuzFz)+3.0*cuNS16[i][j][k]*cfNS16))
#define Force17 ((1-0.5*A)*tNS[17]*(3.0*(cfNS17-uxFxuyFyuzFz)+3.0*cuNS17[i][j][k]*cfNS17))
#define Force18 ((1-0.5*A)*tNS[18]*(3.0*(cfNS18-uxFxuyFyuzFz)+3.0*cuNS18[i][j][k]*cfNS18))

#define vfIn0 (MInverse[0][0]*mfIn0+MInverse[0][1]*mfIn1+MInverse[0][2]*mfIn2+MInverse[0][3]*mfIn3+MInverse[0][4]*mfIn4+MInverse[0][5]*mfIn5+MInverse[0][6]*mfIn6+MInverse[0][7]*mfIn7+MInverse[0][8]*mfIn8\
	+MInverse[0][9]*mfIn9+MInverse[0][10]*mfIn10+MInverse[0][11]*mfIn11+MInverse[0][12]*mfIn12+MInverse[0][13]*mfIn13+MInverse[0][14]*mfIn14+MInverse[0][15]*mfIn15+MInverse[0][16]*mfIn16+MInverse[0][17]*mfIn17+MInverse[0][18]*mfIn18)
#define vfIn1 (MInverse[1][0]*mfIn0+MInverse[1][1]*mfIn1+MInverse[1][2]*mfIn2+MInverse[1][3]*mfIn3+MInverse[1][4]*mfIn4+MInverse[1][5]*mfIn5+MInverse[1][6]*mfIn6+MInverse[1][7]*mfIn7+MInverse[1][8]*mfIn8\
	+MInverse[1][9]*mfIn9+MInverse[1][10]*mfIn10+MInverse[1][11]*mfIn11+MInverse[1][12]*mfIn12+MInverse[1][13]*mfIn13+MInverse[1][14]*mfIn14+MInverse[1][15]*mfIn15+MInverse[1][16]*mfIn16+MInverse[1][17]*mfIn17+MInverse[1][18]*mfIn18)
#define vfIn2 (MInverse[2][0]*mfIn0+MInverse[2][1]*mfIn1+MInverse[2][2]*mfIn2+MInverse[2][3]*mfIn3+MInverse[2][4]*mfIn4+MInverse[2][5]*mfIn5+MInverse[2][6]*mfIn6+MInverse[2][7]*mfIn7+MInverse[2][8]*mfIn8\
	+MInverse[2][9]*mfIn9+MInverse[2][10]*mfIn10+MInverse[2][11]*mfIn11+MInverse[2][12]*mfIn12+MInverse[2][13]*mfIn13+MInverse[2][14]*mfIn14+MInverse[2][15]*mfIn15+MInverse[2][16]*mfIn16+MInverse[2][17]*mfIn17+MInverse[2][18]*mfIn18)
#define vfIn3 (MInverse[3][0]*mfIn0+MInverse[3][1]*mfIn1+MInverse[3][2]*mfIn2+MInverse[3][3]*mfIn3+MInverse[3][4]*mfIn4+MInverse[3][5]*mfIn5+MInverse[3][6]*mfIn6+MInverse[3][7]*mfIn7+MInverse[3][8]*mfIn8\
	+MInverse[3][9]*mfIn9+MInverse[3][10]*mfIn10+MInverse[3][11]*mfIn11+MInverse[3][12]*mfIn12+MInverse[3][13]*mfIn13+MInverse[3][14]*mfIn14+MInverse[3][15]*mfIn15+MInverse[3][16]*mfIn16+MInverse[3][17]*mfIn17+MInverse[3][18]*mfIn18)
#define vfIn4 (MInverse[4][0]*mfIn0+MInverse[4][1]*mfIn1+MInverse[4][2]*mfIn2+MInverse[4][3]*mfIn3+MInverse[4][4]*mfIn4+MInverse[4][5]*mfIn5+MInverse[4][6]*mfIn6+MInverse[4][7]*mfIn7+MInverse[4][8]*mfIn8\
	+MInverse[4][9]*mfIn9+MInverse[4][10]*mfIn10+MInverse[4][11]*mfIn11+MInverse[4][12]*mfIn12+MInverse[4][13]*mfIn13+MInverse[4][14]*mfIn14+MInverse[4][15]*mfIn15+MInverse[4][16]*mfIn16+MInverse[4][17]*mfIn17+MInverse[4][18]*mfIn18)
#define vfIn5 (MInverse[5][0]*mfIn0+MInverse[5][1]*mfIn1+MInverse[5][2]*mfIn2+MInverse[5][3]*mfIn3+MInverse[5][4]*mfIn4+MInverse[5][5]*mfIn5+MInverse[5][6]*mfIn6+MInverse[5][7]*mfIn7+MInverse[5][8]*mfIn8\
	+MInverse[5][9]*mfIn9+MInverse[5][10]*mfIn10+MInverse[5][11]*mfIn11+MInverse[5][12]*mfIn12+MInverse[5][13]*mfIn13+MInverse[5][14]*mfIn14+MInverse[5][15]*mfIn15+MInverse[5][16]*mfIn16+MInverse[5][17]*mfIn17+MInverse[5][18]*mfIn18)
#define vfIn6 (MInverse[6][0]*mfIn0+MInverse[6][1]*mfIn1+MInverse[6][2]*mfIn2+MInverse[6][3]*mfIn3+MInverse[6][4]*mfIn4+MInverse[6][5]*mfIn5+MInverse[6][6]*mfIn6+MInverse[6][7]*mfIn7+MInverse[6][8]*mfIn8\
	+MInverse[6][9]*mfIn9+MInverse[6][10]*mfIn10+MInverse[6][11]*mfIn11+MInverse[6][12]*mfIn12+MInverse[6][13]*mfIn13+MInverse[6][14]*mfIn14+MInverse[6][15]*mfIn15+MInverse[6][16]*mfIn16+MInverse[6][17]*mfIn17+MInverse[6][18]*mfIn18)
#define vfIn7 (MInverse[7][0]*mfIn0+MInverse[7][1]*mfIn1+MInverse[7][2]*mfIn2+MInverse[7][3]*mfIn3+MInverse[7][4]*mfIn4+MInverse[7][5]*mfIn5+MInverse[7][6]*mfIn6+MInverse[7][7]*mfIn7+MInverse[7][8]*mfIn8\
	+MInverse[7][9]*mfIn9+MInverse[7][10]*mfIn10+MInverse[7][11]*mfIn11+MInverse[7][12]*mfIn12+MInverse[7][13]*mfIn13+MInverse[7][14]*mfIn14+MInverse[7][15]*mfIn15+MInverse[7][16]*mfIn16+MInverse[7][17]*mfIn17+MInverse[7][18]*mfIn18)
#define vfIn8 (MInverse[8][0]*mfIn0+MInverse[8][1]*mfIn1+MInverse[8][2]*mfIn2+MInverse[8][3]*mfIn3+MInverse[8][4]*mfIn4+MInverse[8][5]*mfIn5+MInverse[8][6]*mfIn6+MInverse[8][7]*mfIn7+MInverse[8][8]*mfIn8\
	+MInverse[8][9]*mfIn9+MInverse[8][10]*mfIn10+MInverse[8][11]*mfIn11+MInverse[8][12]*mfIn12+MInverse[8][13]*mfIn13+MInverse[8][14]*mfIn14+MInverse[8][15]*mfIn15+MInverse[8][16]*mfIn16+MInverse[8][17]*mfIn17+MInverse[8][18]*mfIn18)
#define vfIn9 (MInverse[9][0]*mfIn0+MInverse[9][1]*mfIn1+MInverse[9][2]*mfIn2+MInverse[9][3]*mfIn3+MInverse[9][4]*mfIn4+MInverse[9][5]*mfIn5+MInverse[9][6]*mfIn6+MInverse[9][7]*mfIn7+MInverse[9][8]*mfIn8\
	+MInverse[9][9]*mfIn9+MInverse[9][10]*mfIn10+MInverse[9][11]*mfIn11+MInverse[9][12]*mfIn12+MInverse[9][13]*mfIn13+MInverse[9][14]*mfIn14+MInverse[9][15]*mfIn15+MInverse[9][16]*mfIn16+MInverse[9][17]*mfIn17+MInverse[9][18]*mfIn18)
#define vfIn10 (MInverse[10][0]*mfIn0+MInverse[10][1]*mfIn1+MInverse[10][2]*mfIn2+MInverse[10][3]*mfIn3+MInverse[10][4]*mfIn4+MInverse[10][5]*mfIn5+MInverse[10][6]*mfIn6+MInverse[10][7]*mfIn7+MInverse[10][8]*mfIn8\
	+MInverse[10][9]*mfIn9+MInverse[10][10]*mfIn10+MInverse[10][11]*mfIn11+MInverse[10][12]*mfIn12+MInverse[10][13]*mfIn13+MInverse[10][14]*mfIn14+MInverse[10][15]*mfIn15+MInverse[10][16]*mfIn16+MInverse[10][17]*mfIn17+MInverse[10][18]*mfIn18)
#define vfIn11 (MInverse[11][0]*mfIn0+MInverse[11][1]*mfIn1+MInverse[11][2]*mfIn2+MInverse[11][3]*mfIn3+MInverse[11][4]*mfIn4+MInverse[11][5]*mfIn5+MInverse[11][6]*mfIn6+MInverse[11][7]*mfIn7+MInverse[11][8]*mfIn8\
	+MInverse[11][9]*mfIn9+MInverse[11][10]*mfIn10+MInverse[11][11]*mfIn11+MInverse[11][12]*mfIn12+MInverse[11][13]*mfIn13+MInverse[11][14]*mfIn14+MInverse[11][15]*mfIn15+MInverse[11][16]*mfIn16+MInverse[11][17]*mfIn17+MInverse[11][18]*mfIn18)
#define vfIn12 (MInverse[12][0]*mfIn0+MInverse[12][1]*mfIn1+MInverse[12][2]*mfIn2+MInverse[12][3]*mfIn3+MInverse[12][4]*mfIn4+MInverse[12][5]*mfIn5+MInverse[12][6]*mfIn6+MInverse[12][7]*mfIn7+MInverse[12][8]*mfIn8\
	+MInverse[12][9]*mfIn9+MInverse[12][10]*mfIn10+MInverse[12][11]*mfIn11+MInverse[12][12]*mfIn12+MInverse[12][13]*mfIn13+MInverse[12][14]*mfIn14+MInverse[12][15]*mfIn15+MInverse[12][16]*mfIn16+MInverse[12][17]*mfIn17+MInverse[12][18]*mfIn18)
#define vfIn13 (MInverse[13][0]*mfIn0+MInverse[13][1]*mfIn1+MInverse[13][2]*mfIn2+MInverse[13][3]*mfIn3+MInverse[13][4]*mfIn4+MInverse[13][5]*mfIn5+MInverse[13][6]*mfIn6+MInverse[13][7]*mfIn7+MInverse[13][8]*mfIn8\
	+MInverse[13][9]*mfIn9+MInverse[13][10]*mfIn10+MInverse[13][11]*mfIn11+MInverse[13][12]*mfIn12+MInverse[13][13]*mfIn13+MInverse[13][14]*mfIn14+MInverse[13][15]*mfIn15+MInverse[13][16]*mfIn16+MInverse[13][17]*mfIn17+MInverse[13][18]*mfIn18)
#define vfIn14 (MInverse[14][0]*mfIn0+MInverse[14][1]*mfIn1+MInverse[14][2]*mfIn2+MInverse[14][3]*mfIn3+MInverse[14][4]*mfIn4+MInverse[14][5]*mfIn5+MInverse[14][6]*mfIn6+MInverse[14][7]*mfIn7+MInverse[14][8]*mfIn8\
	+MInverse[14][9]*mfIn9+MInverse[14][10]*mfIn10+MInverse[14][11]*mfIn11+MInverse[14][12]*mfIn12+MInverse[14][13]*mfIn13+MInverse[14][14]*mfIn14+MInverse[14][15]*mfIn15+MInverse[14][16]*mfIn16+MInverse[14][17]*mfIn17+MInverse[14][18]*mfIn18)
#define vfIn15 (MInverse[15][0]*mfIn0+MInverse[15][1]*mfIn1+MInverse[15][2]*mfIn2+MInverse[15][3]*mfIn3+MInverse[15][4]*mfIn4+MInverse[15][5]*mfIn5+MInverse[15][6]*mfIn6+MInverse[15][7]*mfIn7+MInverse[15][8]*mfIn8\
	+MInverse[15][9]*mfIn9+MInverse[15][10]*mfIn10+MInverse[15][11]*mfIn11+MInverse[15][12]*mfIn12+MInverse[15][13]*mfIn13+MInverse[15][14]*mfIn14+MInverse[15][15]*mfIn15+MInverse[15][16]*mfIn16+MInverse[15][17]*mfIn17+MInverse[15][18]*mfIn18)
#define vfIn16 (MInverse[16][0]*mfIn0+MInverse[16][1]*mfIn1+MInverse[16][2]*mfIn2+MInverse[16][3]*mfIn3+MInverse[16][4]*mfIn4+MInverse[16][5]*mfIn5+MInverse[16][6]*mfIn6+MInverse[16][7]*mfIn7+MInverse[16][8]*mfIn8\
	+MInverse[16][9]*mfIn9+MInverse[16][10]*mfIn10+MInverse[16][11]*mfIn11+MInverse[16][12]*mfIn12+MInverse[16][13]*mfIn13+MInverse[16][14]*mfIn14+MInverse[16][15]*mfIn15+MInverse[16][16]*mfIn16+MInverse[16][17]*mfIn17+MInverse[16][18]*mfIn18)
#define vfIn17 (MInverse[17][0]*mfIn0+MInverse[17][1]*mfIn1+MInverse[17][2]*mfIn2+MInverse[17][3]*mfIn3+MInverse[17][4]*mfIn4+MInverse[17][5]*mfIn5+MInverse[17][6]*mfIn6+MInverse[17][7]*mfIn7+MInverse[17][8]*mfIn8\
	+MInverse[17][9]*mfIn9+MInverse[17][10]*mfIn10+MInverse[17][11]*mfIn11+MInverse[17][12]*mfIn12+MInverse[17][13]*mfIn13+MInverse[17][14]*mfIn14+MInverse[17][15]*mfIn15+MInverse[17][16]*mfIn16+MInverse[17][17]*mfIn17+MInverse[17][18]*mfIn18)
#define vfIn18 (MInverse[18][0]*mfIn0+MInverse[18][1]*mfIn1+MInverse[18][2]*mfIn2+MInverse[18][3]*mfIn3+MInverse[18][4]*mfIn4+MInverse[18][5]*mfIn5+MInverse[18][6]*mfIn6+MInverse[18][7]*mfIn7+MInverse[18][8]*mfIn8\
	+MInverse[18][9]*mfIn9+MInverse[18][10]*mfIn10+MInverse[18][11]*mfIn11+MInverse[18][12]*mfIn12+MInverse[18][13]*mfIn13+MInverse[18][14]*mfIn14+MInverse[18][15]*mfIn15+MInverse[18][16]*mfIn16+MInverse[18][17]*mfIn17+MInverse[18][18]*mfIn18)

#define vfEq0 (MInverse[0][0]*mfEq0+MInverse[0][1]*mfEq1+MInverse[0][2]*mfEq2+MInverse[0][3]*mfEq3+MInverse[0][4]*mfEq4+MInverse[0][5]*mfEq5+MInverse[0][6]*mfEq6+MInverse[0][7]*mfEq7+MInverse[0][8]*mfEq8\
	+MInverse[0][9]*mfEq9+MInverse[0][10]*mfEq10+MInverse[0][11]*mfEq11+MInverse[0][12]*mfEq12+MInverse[0][13]*mfEq13+MInverse[0][14]*mfEq14+MInverse[0][15]*mfEq15+MInverse[0][16]*mfEq16+MInverse[0][17]*mfEq17+MInverse[0][18]*mfEq18)
#define vfEq1 (MInverse[1][0]*mfEq0+MInverse[1][1]*mfEq1+MInverse[1][2]*mfEq2+MInverse[1][3]*mfEq3+MInverse[1][4]*mfEq4+MInverse[1][5]*mfEq5+MInverse[1][6]*mfEq6+MInverse[1][7]*mfEq7+MInverse[1][8]*mfEq8\
	+MInverse[1][9]*mfEq9+MInverse[1][10]*mfEq10+MInverse[1][11]*mfEq11+MInverse[1][12]*mfEq12+MInverse[1][13]*mfEq13+MInverse[1][14]*mfEq14+MInverse[1][15]*mfEq15+MInverse[1][16]*mfEq16+MInverse[1][17]*mfEq17+MInverse[1][18]*mfEq18)
#define vfEq2 (MInverse[2][0]*mfEq0+MInverse[2][1]*mfEq1+MInverse[2][2]*mfEq2+MInverse[2][3]*mfEq3+MInverse[2][4]*mfEq4+MInverse[2][5]*mfEq5+MInverse[2][6]*mfEq6+MInverse[2][7]*mfEq7+MInverse[2][8]*mfEq8\
	+MInverse[2][9]*mfEq9+MInverse[2][10]*mfEq10+MInverse[2][11]*mfEq11+MInverse[2][12]*mfEq12+MInverse[2][13]*mfEq13+MInverse[2][14]*mfEq14+MInverse[2][15]*mfEq15+MInverse[2][16]*mfEq16+MInverse[2][17]*mfEq17+MInverse[2][18]*mfEq18)
#define vfEq3 (MInverse[3][0]*mfEq0+MInverse[3][1]*mfEq1+MInverse[3][2]*mfEq2+MInverse[3][3]*mfEq3+MInverse[3][4]*mfEq4+MInverse[3][5]*mfEq5+MInverse[3][6]*mfEq6+MInverse[3][7]*mfEq7+MInverse[3][8]*mfEq8\
	+MInverse[3][9]*mfEq9+MInverse[3][10]*mfEq10+MInverse[3][11]*mfEq11+MInverse[3][12]*mfEq12+MInverse[3][13]*mfEq13+MInverse[3][14]*mfEq14+MInverse[3][15]*mfEq15+MInverse[3][16]*mfEq16+MInverse[3][17]*mfEq17+MInverse[3][18]*mfEq18)
#define vfEq4 (MInverse[4][0]*mfEq0+MInverse[4][1]*mfEq1+MInverse[4][2]*mfEq2+MInverse[4][3]*mfEq3+MInverse[4][4]*mfEq4+MInverse[4][5]*mfEq5+MInverse[4][6]*mfEq6+MInverse[4][7]*mfEq7+MInverse[4][8]*mfEq8\
	+MInverse[4][9]*mfEq9+MInverse[4][10]*mfEq10+MInverse[4][11]*mfEq11+MInverse[4][12]*mfEq12+MInverse[4][13]*mfEq13+MInverse[4][14]*mfEq14+MInverse[4][15]*mfEq15+MInverse[4][16]*mfEq16+MInverse[4][17]*mfEq17+MInverse[4][18]*mfEq18)
#define vfEq5 (MInverse[5][0]*mfEq0+MInverse[5][1]*mfEq1+MInverse[5][2]*mfEq2+MInverse[5][3]*mfEq3+MInverse[5][4]*mfEq4+MInverse[5][5]*mfEq5+MInverse[5][6]*mfEq6+MInverse[5][7]*mfEq7+MInverse[5][8]*mfEq8\
	+MInverse[5][9]*mfEq9+MInverse[5][10]*mfEq10+MInverse[5][11]*mfEq11+MInverse[5][12]*mfEq12+MInverse[5][13]*mfEq13+MInverse[5][14]*mfEq14+MInverse[5][15]*mfEq15+MInverse[5][16]*mfEq16+MInverse[5][17]*mfEq17+MInverse[5][18]*mfEq18)
#define vfEq6 (MInverse[6][0]*mfEq0+MInverse[6][1]*mfEq1+MInverse[6][2]*mfEq2+MInverse[6][3]*mfEq3+MInverse[6][4]*mfEq4+MInverse[6][5]*mfEq5+MInverse[6][6]*mfEq6+MInverse[6][7]*mfEq7+MInverse[6][8]*mfEq8\
	+MInverse[6][9]*mfEq9+MInverse[6][10]*mfEq10+MInverse[6][11]*mfEq11+MInverse[6][12]*mfEq12+MInverse[6][13]*mfEq13+MInverse[6][14]*mfEq14+MInverse[6][15]*mfEq15+MInverse[6][16]*mfEq16+MInverse[6][17]*mfEq17+MInverse[6][18]*mfEq18)
#define vfEq7 (MInverse[7][0]*mfEq0+MInverse[7][1]*mfEq1+MInverse[7][2]*mfEq2+MInverse[7][3]*mfEq3+MInverse[7][4]*mfEq4+MInverse[7][5]*mfEq5+MInverse[7][6]*mfEq6+MInverse[7][7]*mfEq7+MInverse[7][8]*mfEq8\
	+MInverse[7][9]*mfEq9+MInverse[7][10]*mfEq10+MInverse[7][11]*mfEq11+MInverse[7][12]*mfEq12+MInverse[7][13]*mfEq13+MInverse[7][14]*mfEq14+MInverse[7][15]*mfEq15+MInverse[7][16]*mfEq16+MInverse[7][17]*mfEq17+MInverse[7][18]*mfEq18)
#define vfEq8 (MInverse[8][0]*mfEq0+MInverse[8][1]*mfEq1+MInverse[8][2]*mfEq2+MInverse[8][3]*mfEq3+MInverse[8][4]*mfEq4+MInverse[8][5]*mfEq5+MInverse[8][6]*mfEq6+MInverse[8][7]*mfEq7+MInverse[8][8]*mfEq8\
	+MInverse[8][9]*mfEq9+MInverse[8][10]*mfEq10+MInverse[8][11]*mfEq11+MInverse[8][12]*mfEq12+MInverse[8][13]*mfEq13+MInverse[8][14]*mfEq14+MInverse[8][15]*mfEq15+MInverse[8][16]*mfEq16+MInverse[8][17]*mfEq17+MInverse[8][18]*mfEq18)
#define vfEq9 (MInverse[9][0]*mfEq0+MInverse[9][1]*mfEq1+MInverse[9][2]*mfEq2+MInverse[9][3]*mfEq3+MInverse[9][4]*mfEq4+MInverse[9][5]*mfEq5+MInverse[9][6]*mfEq6+MInverse[9][7]*mfEq7+MInverse[9][8]*mfEq8\
	+MInverse[9][9]*mfEq9+MInverse[9][10]*mfEq10+MInverse[9][11]*mfEq11+MInverse[9][12]*mfEq12+MInverse[9][13]*mfEq13+MInverse[9][14]*mfEq14+MInverse[9][15]*mfEq15+MInverse[9][16]*mfEq16+MInverse[9][17]*mfEq17+MInverse[9][18]*mfEq18)
#define vfEq10 (MInverse[10][0]*mfEq0+MInverse[10][1]*mfEq1+MInverse[10][2]*mfEq2+MInverse[10][3]*mfEq3+MInverse[10][4]*mfEq4+MInverse[10][5]*mfEq5+MInverse[10][6]*mfEq6+MInverse[10][7]*mfEq7+MInverse[10][8]*mfEq8\
	+MInverse[10][9]*mfEq9+MInverse[10][10]*mfEq10+MInverse[10][11]*mfEq11+MInverse[10][12]*mfEq12+MInverse[10][13]*mfEq13+MInverse[10][14]*mfEq14+MInverse[10][15]*mfEq15+MInverse[10][16]*mfEq16+MInverse[10][17]*mfEq17+MInverse[10][18]*mfEq18)
#define vfEq11 (MInverse[11][0]*mfEq0+MInverse[11][1]*mfEq1+MInverse[11][2]*mfEq2+MInverse[11][3]*mfEq3+MInverse[11][4]*mfEq4+MInverse[11][5]*mfEq5+MInverse[11][6]*mfEq6+MInverse[11][7]*mfEq7+MInverse[11][8]*mfEq8\
	+MInverse[11][9]*mfEq9+MInverse[11][10]*mfEq10+MInverse[11][11]*mfEq11+MInverse[11][12]*mfEq12+MInverse[11][13]*mfEq13+MInverse[11][14]*mfEq14+MInverse[11][15]*mfEq15+MInverse[11][16]*mfEq16+MInverse[11][17]*mfEq17+MInverse[11][18]*mfEq18)
#define vfEq12 (MInverse[12][0]*mfEq0+MInverse[12][1]*mfEq1+MInverse[12][2]*mfEq2+MInverse[12][3]*mfEq3+MInverse[12][4]*mfEq4+MInverse[12][5]*mfEq5+MInverse[12][6]*mfEq6+MInverse[12][7]*mfEq7+MInverse[12][8]*mfEq8\
	+MInverse[12][9]*mfEq9+MInverse[12][10]*mfEq10+MInverse[12][11]*mfEq11+MInverse[12][12]*mfEq12+MInverse[12][13]*mfEq13+MInverse[12][14]*mfEq14+MInverse[12][15]*mfEq15+MInverse[12][16]*mfEq16+MInverse[12][17]*mfEq17+MInverse[12][18]*mfEq18)
#define vfEq13 (MInverse[13][0]*mfEq0+MInverse[13][1]*mfEq1+MInverse[13][2]*mfEq2+MInverse[13][3]*mfEq3+MInverse[13][4]*mfEq4+MInverse[13][5]*mfEq5+MInverse[13][6]*mfEq6+MInverse[13][7]*mfEq7+MInverse[13][8]*mfEq8\
	+MInverse[13][9]*mfEq9+MInverse[13][10]*mfEq10+MInverse[13][11]*mfEq11+MInverse[13][12]*mfEq12+MInverse[13][13]*mfEq13+MInverse[13][14]*mfEq14+MInverse[13][15]*mfEq15+MInverse[13][16]*mfEq16+MInverse[13][17]*mfEq17+MInverse[13][18]*mfEq18)
#define vfEq14 (MInverse[14][0]*mfEq0+MInverse[14][1]*mfEq1+MInverse[14][2]*mfEq2+MInverse[14][3]*mfEq3+MInverse[14][4]*mfEq4+MInverse[14][5]*mfEq5+MInverse[14][6]*mfEq6+MInverse[14][7]*mfEq7+MInverse[14][8]*mfEq8\
	+MInverse[14][9]*mfEq9+MInverse[14][10]*mfEq10+MInverse[14][11]*mfEq11+MInverse[14][12]*mfEq12+MInverse[14][13]*mfEq13+MInverse[14][14]*mfEq14+MInverse[14][15]*mfEq15+MInverse[14][16]*mfEq16+MInverse[14][17]*mfEq17+MInverse[14][18]*mfEq18)
#define vfEq15 (MInverse[15][0]*mfEq0+MInverse[15][1]*mfEq1+MInverse[15][2]*mfEq2+MInverse[15][3]*mfEq3+MInverse[15][4]*mfEq4+MInverse[15][5]*mfEq5+MInverse[15][6]*mfEq6+MInverse[15][7]*mfEq7+MInverse[15][8]*mfEq8\
	+MInverse[15][9]*mfEq9+MInverse[15][10]*mfEq10+MInverse[15][11]*mfEq11+MInverse[15][12]*mfEq12+MInverse[15][13]*mfEq13+MInverse[15][14]*mfEq14+MInverse[15][15]*mfEq15+MInverse[15][16]*mfEq16+MInverse[15][17]*mfEq17+MInverse[15][18]*mfEq18)
#define vfEq16 (MInverse[16][0]*mfEq0+MInverse[16][1]*mfEq1+MInverse[16][2]*mfEq2+MInverse[16][3]*mfEq3+MInverse[16][4]*mfEq4+MInverse[16][5]*mfEq5+MInverse[16][6]*mfEq6+MInverse[16][7]*mfEq7+MInverse[16][8]*mfEq8\
	+MInverse[16][9]*mfEq9+MInverse[16][10]*mfEq10+MInverse[16][11]*mfEq11+MInverse[16][12]*mfEq12+MInverse[16][13]*mfEq13+MInverse[16][14]*mfEq14+MInverse[16][15]*mfEq15+MInverse[16][16]*mfEq16+MInverse[16][17]*mfEq17+MInverse[16][18]*mfEq18)
#define vfEq17 (MInverse[17][0]*mfEq0+MInverse[17][1]*mfEq1+MInverse[17][2]*mfEq2+MInverse[17][3]*mfEq3+MInverse[17][4]*mfEq4+MInverse[17][5]*mfEq5+MInverse[17][6]*mfEq6+MInverse[17][7]*mfEq7+MInverse[17][8]*mfEq8\
	+MInverse[17][9]*mfEq9+MInverse[17][10]*mfEq10+MInverse[17][11]*mfEq11+MInverse[17][12]*mfEq12+MInverse[17][13]*mfEq13+MInverse[17][14]*mfEq14+MInverse[17][15]*mfEq15+MInverse[17][16]*mfEq16+MInverse[17][17]*mfEq17+MInverse[17][18]*mfEq18)
#define vfEq18 (MInverse[18][0]*mfEq0+MInverse[18][1]*mfEq1+MInverse[18][2]*mfEq2+MInverse[18][3]*mfEq3+MInverse[18][4]*mfEq4+MInverse[18][5]*mfEq5+MInverse[18][6]*mfEq6+MInverse[18][7]*mfEq7+MInverse[18][8]*mfEq8\
	+MInverse[18][9]*mfEq9+MInverse[18][10]*mfEq10+MInverse[18][11]*mfEq11+MInverse[18][12]*mfEq12+MInverse[18][13]*mfEq13+MInverse[18][14]*mfEq14+MInverse[18][15]*mfEq15+MInverse[18][16]*mfEq16+MInverse[18][17]*mfEq17+MInverse[18][18]*mfEq18)

#define vForce0 (MInverse[0][0]*mForce0+MInverse[0][1]*mForce1+MInverse[0][2]*mForce2+MInverse[0][3]*mForce3+MInverse[0][4]*mForce4+MInverse[0][5]*mForce5+MInverse[0][6]*mForce6+MInverse[0][7]*mForce7+MInverse[0][8]*mForce8\
	+MInverse[0][9]*mForce9+MInverse[0][10]*mForce10+MInverse[0][11]*mForce11+MInverse[0][12]*mForce12+MInverse[0][13]*mForce13+MInverse[0][14]*mForce14+MInverse[0][15]*mForce15+MInverse[0][16]*mForce16+MInverse[0][17]*mForce17+MInverse[0][18]*mForce18)
#define vForce1 (MInverse[1][0]*mForce0+MInverse[1][1]*mForce1+MInverse[1][2]*mForce2+MInverse[1][3]*mForce3+MInverse[1][4]*mForce4+MInverse[1][5]*mForce5+MInverse[1][6]*mForce6+MInverse[1][7]*mForce7+MInverse[1][8]*mForce8\
	+MInverse[1][9]*mForce9+MInverse[1][10]*mForce10+MInverse[1][11]*mForce11+MInverse[1][12]*mForce12+MInverse[1][13]*mForce13+MInverse[1][14]*mForce14+MInverse[1][15]*mForce15+MInverse[1][16]*mForce16+MInverse[1][17]*mForce17+MInverse[1][18]*mForce18)
#define vForce2 (MInverse[2][0]*mForce0+MInverse[2][1]*mForce1+MInverse[2][2]*mForce2+MInverse[2][3]*mForce3+MInverse[2][4]*mForce4+MInverse[2][5]*mForce5+MInverse[2][6]*mForce6+MInverse[2][7]*mForce7+MInverse[2][8]*mForce8\
	+MInverse[2][9]*mForce9+MInverse[2][10]*mForce10+MInverse[2][11]*mForce11+MInverse[2][12]*mForce12+MInverse[2][13]*mForce13+MInverse[2][14]*mForce14+MInverse[2][15]*mForce15+MInverse[2][16]*mForce16+MInverse[2][17]*mForce17+MInverse[2][18]*mForce18)
#define vForce3 (MInverse[3][0]*mForce0+MInverse[3][1]*mForce1+MInverse[3][2]*mForce2+MInverse[3][3]*mForce3+MInverse[3][4]*mForce4+MInverse[3][5]*mForce5+MInverse[3][6]*mForce6+MInverse[3][7]*mForce7+MInverse[3][8]*mForce8\
	+MInverse[3][9]*mForce9+MInverse[3][10]*mForce10+MInverse[3][11]*mForce11+MInverse[3][12]*mForce12+MInverse[3][13]*mForce13+MInverse[3][14]*mForce14+MInverse[3][15]*mForce15+MInverse[3][16]*mForce16+MInverse[3][17]*mForce17+MInverse[3][18]*mForce18)
#define vForce4 (MInverse[4][0]*mForce0+MInverse[4][1]*mForce1+MInverse[4][2]*mForce2+MInverse[4][3]*mForce3+MInverse[4][4]*mForce4+MInverse[4][5]*mForce5+MInverse[4][6]*mForce6+MInverse[4][7]*mForce7+MInverse[4][8]*mForce8\
	+MInverse[4][9]*mForce9+MInverse[4][10]*mForce10+MInverse[4][11]*mForce11+MInverse[4][12]*mForce12+MInverse[4][13]*mForce13+MInverse[4][14]*mForce14+MInverse[4][15]*mForce15+MInverse[4][16]*mForce16+MInverse[4][17]*mForce17+MInverse[4][18]*mForce18)
#define vForce5 (MInverse[5][0]*mForce0+MInverse[5][1]*mForce1+MInverse[5][2]*mForce2+MInverse[5][3]*mForce3+MInverse[5][4]*mForce4+MInverse[5][5]*mForce5+MInverse[5][6]*mForce6+MInverse[5][7]*mForce7+MInverse[5][8]*mForce8\
	+MInverse[5][9]*mForce9+MInverse[5][10]*mForce10+MInverse[5][11]*mForce11+MInverse[5][12]*mForce12+MInverse[5][13]*mForce13+MInverse[5][14]*mForce14+MInverse[5][15]*mForce15+MInverse[5][16]*mForce16+MInverse[5][17]*mForce17+MInverse[5][18]*mForce18)
#define vForce6 (MInverse[6][0]*mForce0+MInverse[6][1]*mForce1+MInverse[6][2]*mForce2+MInverse[6][3]*mForce3+MInverse[6][4]*mForce4+MInverse[6][5]*mForce5+MInverse[6][6]*mForce6+MInverse[6][7]*mForce7+MInverse[6][8]*mForce8\
	+MInverse[6][9]*mForce9+MInverse[6][10]*mForce10+MInverse[6][11]*mForce11+MInverse[6][12]*mForce12+MInverse[6][13]*mForce13+MInverse[6][14]*mForce14+MInverse[6][15]*mForce15+MInverse[6][16]*mForce16+MInverse[6][17]*mForce17+MInverse[6][18]*mForce18)
#define vForce7 (MInverse[7][0]*mForce0+MInverse[7][1]*mForce1+MInverse[7][2]*mForce2+MInverse[7][3]*mForce3+MInverse[7][4]*mForce4+MInverse[7][5]*mForce5+MInverse[7][6]*mForce6+MInverse[7][7]*mForce7+MInverse[7][8]*mForce8\
	+MInverse[7][9]*mForce9+MInverse[7][10]*mForce10+MInverse[7][11]*mForce11+MInverse[7][12]*mForce12+MInverse[7][13]*mForce13+MInverse[7][14]*mForce14+MInverse[7][15]*mForce15+MInverse[7][16]*mForce16+MInverse[7][17]*mForce17+MInverse[7][18]*mForce18)
#define vForce8 (MInverse[8][0]*mForce0+MInverse[8][1]*mForce1+MInverse[8][2]*mForce2+MInverse[8][3]*mForce3+MInverse[8][4]*mForce4+MInverse[8][5]*mForce5+MInverse[8][6]*mForce6+MInverse[8][7]*mForce7+MInverse[8][8]*mForce8\
	+MInverse[8][9]*mForce9+MInverse[8][10]*mForce10+MInverse[8][11]*mForce11+MInverse[8][12]*mForce12+MInverse[8][13]*mForce13+MInverse[8][14]*mForce14+MInverse[8][15]*mForce15+MInverse[8][16]*mForce16+MInverse[8][17]*mForce17+MInverse[8][18]*mForce18)
#define vForce9 (MInverse[9][0]*mForce0+MInverse[9][1]*mForce1+MInverse[9][2]*mForce2+MInverse[9][3]*mForce3+MInverse[9][4]*mForce4+MInverse[9][5]*mForce5+MInverse[9][6]*mForce6+MInverse[9][7]*mForce7+MInverse[9][8]*mForce8\
	+MInverse[9][9]*mForce9+MInverse[9][10]*mForce10+MInverse[9][11]*mForce11+MInverse[9][12]*mForce12+MInverse[9][13]*mForce13+MInverse[9][14]*mForce14+MInverse[9][15]*mForce15+MInverse[9][16]*mForce16+MInverse[9][17]*mForce17+MInverse[9][18]*mForce18)
#define vForce10 (MInverse[10][0]*mForce0+MInverse[10][1]*mForce1+MInverse[10][2]*mForce2+MInverse[10][3]*mForce3+MInverse[10][4]*mForce4+MInverse[10][5]*mForce5+MInverse[10][6]*mForce6+MInverse[10][7]*mForce7+MInverse[10][8]*mForce8\
	+MInverse[10][9]*mForce9+MInverse[10][10]*mForce10+MInverse[10][11]*mForce11+MInverse[10][12]*mForce12+MInverse[10][13]*mForce13+MInverse[10][14]*mForce14+MInverse[10][15]*mForce15+MInverse[10][16]*mForce16+MInverse[10][17]*mForce17+MInverse[10][18]*mForce18)
#define vForce11 (MInverse[11][0]*mForce0+MInverse[11][1]*mForce1+MInverse[11][2]*mForce2+MInverse[11][3]*mForce3+MInverse[11][4]*mForce4+MInverse[11][5]*mForce5+MInverse[11][6]*mForce6+MInverse[11][7]*mForce7+MInverse[11][8]*mForce8\
	+MInverse[11][9]*mForce9+MInverse[11][10]*mForce10+MInverse[11][11]*mForce11+MInverse[11][12]*mForce12+MInverse[11][13]*mForce13+MInverse[11][14]*mForce14+MInverse[11][15]*mForce15+MInverse[11][16]*mForce16+MInverse[11][17]*mForce17+MInverse[11][18]*mForce18)
#define vForce12 (MInverse[12][0]*mForce0+MInverse[12][1]*mForce1+MInverse[12][2]*mForce2+MInverse[12][3]*mForce3+MInverse[12][4]*mForce4+MInverse[12][5]*mForce5+MInverse[12][6]*mForce6+MInverse[12][7]*mForce7+MInverse[12][8]*mForce8\
	+MInverse[12][9]*mForce9+MInverse[12][10]*mForce10+MInverse[12][11]*mForce11+MInverse[12][12]*mForce12+MInverse[12][13]*mForce13+MInverse[12][14]*mForce14+MInverse[12][15]*mForce15+MInverse[12][16]*mForce16+MInverse[12][17]*mForce17+MInverse[12][18]*mForce18)
#define vForce13 (MInverse[13][0]*mForce0+MInverse[13][1]*mForce1+MInverse[13][2]*mForce2+MInverse[13][3]*mForce3+MInverse[13][4]*mForce4+MInverse[13][5]*mForce5+MInverse[13][6]*mForce6+MInverse[13][7]*mForce7+MInverse[13][8]*mForce8\
	+MInverse[13][9]*mForce9+MInverse[13][10]*mForce10+MInverse[13][11]*mForce11+MInverse[13][12]*mForce12+MInverse[13][13]*mForce13+MInverse[13][14]*mForce14+MInverse[13][15]*mForce15+MInverse[13][16]*mForce16+MInverse[13][17]*mForce17+MInverse[13][18]*mForce18)
#define vForce14 (MInverse[14][0]*mForce0+MInverse[14][1]*mForce1+MInverse[14][2]*mForce2+MInverse[14][3]*mForce3+MInverse[14][4]*mForce4+MInverse[14][5]*mForce5+MInverse[14][6]*mForce6+MInverse[14][7]*mForce7+MInverse[14][8]*mForce8\
	+MInverse[14][9]*mForce9+MInverse[14][10]*mForce10+MInverse[14][11]*mForce11+MInverse[14][12]*mForce12+MInverse[14][13]*mForce13+MInverse[14][14]*mForce14+MInverse[14][15]*mForce15+MInverse[14][16]*mForce16+MInverse[14][17]*mForce17+MInverse[14][18]*mForce18)
#define vForce15 (MInverse[15][0]*mForce0+MInverse[15][1]*mForce1+MInverse[15][2]*mForce2+MInverse[15][3]*mForce3+MInverse[15][4]*mForce4+MInverse[15][5]*mForce5+MInverse[15][6]*mForce6+MInverse[15][7]*mForce7+MInverse[15][8]*mForce8\
	+MInverse[15][9]*mForce9+MInverse[15][10]*mForce10+MInverse[15][11]*mForce11+MInverse[15][12]*mForce12+MInverse[15][13]*mForce13+MInverse[15][14]*mForce14+MInverse[15][15]*mForce15+MInverse[15][16]*mForce16+MInverse[15][17]*mForce17+MInverse[15][18]*mForce18)
#define vForce16 (MInverse[16][0]*mForce0+MInverse[16][1]*mForce1+MInverse[16][2]*mForce2+MInverse[16][3]*mForce3+MInverse[16][4]*mForce4+MInverse[16][5]*mForce5+MInverse[16][6]*mForce6+MInverse[16][7]*mForce7+MInverse[16][8]*mForce8\
	+MInverse[16][9]*mForce9+MInverse[16][10]*mForce10+MInverse[16][11]*mForce11+MInverse[16][12]*mForce12+MInverse[16][13]*mForce13+MInverse[16][14]*mForce14+MInverse[16][15]*mForce15+MInverse[16][16]*mForce16+MInverse[16][17]*mForce17+MInverse[16][18]*mForce18)
#define vForce17 (MInverse[17][0]*mForce0+MInverse[17][1]*mForce1+MInverse[17][2]*mForce2+MInverse[17][3]*mForce3+MInverse[17][4]*mForce4+MInverse[17][5]*mForce5+MInverse[17][6]*mForce6+MInverse[17][7]*mForce7+MInverse[17][8]*mForce8\
	+MInverse[17][9]*mForce9+MInverse[17][10]*mForce10+MInverse[17][11]*mForce11+MInverse[17][12]*mForce12+MInverse[17][13]*mForce13+MInverse[17][14]*mForce14+MInverse[17][15]*mForce15+MInverse[17][16]*mForce16+MInverse[17][17]*mForce17+MInverse[17][18]*mForce18)
#define vForce18 (MInverse[18][0]*mForce0+MInverse[18][1]*mForce1+MInverse[18][2]*mForce2+MInverse[18][3]*mForce3+MInverse[18][4]*mForce4+MInverse[18][5]*mForce5+MInverse[18][6]*mForce6+MInverse[18][7]*mForce7+MInverse[18][8]*mForce8\
	+MInverse[18][9]*mForce9+MInverse[18][10]*mForce10+MInverse[18][11]*mForce11+MInverse[18][12]*mForce12+MInverse[18][13]*mForce13+MInverse[18][14]*mForce14+MInverse[18][15]*mForce15+MInverse[18][16]*mForce16+MInverse[18][17]*mForce17+MInverse[18][18]*mForce18)

double ***fIn0=NULL;
double ***fIn1=NULL;
double ***fIn2=NULL;
double ***fIn3=NULL;
double ***fIn4=NULL;
double ***fIn5=NULL;
double ***fIn6=NULL;
double ***fIn7=NULL;
double ***fIn8=NULL;
double ***fIn9=NULL;
double ***fIn10=NULL;
double ***fIn11=NULL;
double ***fIn12=NULL;
double ***fIn13=NULL;
double ***fIn14=NULL;
double ***fIn15=NULL;
double ***fIn16=NULL;
double ***fIn17=NULL;
double ***fIn18=NULL;

double ***rho=NULL;
double ***rhoLastStep=NULL;
double ***rho1=NULL;
double ***rho2=NULL;
double ***rho3=NULL;
double ***rho4=NULL;
double ***rho5=NULL;
double ***rho6=NULL;
double ***rho7=NULL;
double ***rho8=NULL;
double ***rho9=NULL;
double ***rho10=NULL;
double ***rho11=NULL;
double ***rho12=NULL;
double ***rho13=NULL;
double ***rho14=NULL;
double ***rho15=NULL;
double ***rho16=NULL;
double ***rho17=NULL;
double ***rho18=NULL;
double ***rho19 = NULL;
double ***rho20 = NULL;
double ***rho21 = NULL;
double ***rho22 = NULL;
double ***rho23 = NULL;
double ***rho24 = NULL;
double ***rho25 = NULL;
double ***rho26 = NULL;

double ***ux=NULL;
double ***uy=NULL;
double ***uz=NULL;
double ***Fx=NULL;
double ***Fy=NULL;
double ***Fz=NULL;
double ***bodyForce=NULL;

double ***pressure=NULL;
double ***pressure1=NULL;
double ***pressure2=NULL;
double ***pressure3=NULL;
double ***pressure4=NULL;
double ***pressure5=NULL;
double ***pressure6=NULL;
double ***pressure7=NULL;
double ***pressure8=NULL;
double ***pressure9=NULL;
double ***pressure10=NULL;
double ***pressure11=NULL;
double ***pressure12=NULL;
double ***pressure13=NULL;
double ***pressure14=NULL;
double ***pressure15=NULL;
double ***pressure16=NULL;
double ***pressure17=NULL;
double ***pressure18=NULL;

void FieldArrangeTwoPhase(int m, int n, int q){
	fIn0=memoryarrange(m,n,q);
	fIn1=memoryarrange(m,n,q);
	fIn2=memoryarrange(m,n,q);
	fIn3=memoryarrange(m,n,q);
	fIn4=memoryarrange(m,n,q);
	fIn5=memoryarrange(m,n,q);
	fIn6=memoryarrange(m,n,q);
	fIn7=memoryarrange(m,n,q);
	fIn8=memoryarrange(m,n,q);
	fIn9=memoryarrange(m,n,q);
	fIn10=memoryarrange(m,n,q);
	fIn11=memoryarrange(m,n,q);
	fIn12=memoryarrange(m,n,q);
	fIn13=memoryarrange(m,n,q);
	fIn14=memoryarrange(m,n,q);
	fIn15=memoryarrange(m,n,q);
	fIn16=memoryarrange(m,n,q);
	fIn17=memoryarrange(m,n,q);
	fIn18=memoryarrange(m,n,q);
	rho=memoryarrange(m,n,q);
	rhoLastStep=memoryarrange(m,n,q);
	rho1=memoryarrange(m,n,q);
	rho2=memoryarrange(m,n,q);
	rho3=memoryarrange(m,n,q);
	rho4=memoryarrange(m,n,q);
	rho5=memoryarrange(m,n,q);
	rho6=memoryarrange(m,n,q);
	rho7=memoryarrange(m,n,q);
	rho8=memoryarrange(m,n,q);
	rho9=memoryarrange(m,n,q);
	rho10=memoryarrange(m,n,q);
	rho11=memoryarrange(m,n,q);
	rho12=memoryarrange(m,n,q);
	rho13=memoryarrange(m,n,q);
	rho14=memoryarrange(m,n,q);
	rho15=memoryarrange(m,n,q);
	rho16=memoryarrange(m,n,q);
	rho17=memoryarrange(m,n,q);
	rho18=memoryarrange(m,n,q);
//	rho19 = memoryarrange(m, n, q);
//	rho20 = memoryarrange(m, n, q);
//	rho21 = memoryarrange(m, n, q);
//	rho22 = memoryarrange(m, n, q);
//	rho23 = memoryarrange(m, n, q);
//	rho24 = memoryarrange(m, n, q);
//	rho25 = memoryarrange(m, n, q);
//	rho26 = memoryarrange(m, n, q);
	ux=memoryarrange(m,n,q);
	uy=memoryarrange(m,n,q);
	uz=memoryarrange(m,n,q);
	Fx=memoryarrange(m,n,q);
	Fy=memoryarrange(m,n,q);
	Fz=memoryarrange(m,n,q);
	bodyForce=memoryarrange(m,n,q);
	pressure=memoryarrange(m,n,q);
	pressure1=memoryarrange(m,n,q);
	pressure2=memoryarrange(m,n,q);
	pressure3=memoryarrange(m,n,q);
	pressure4=memoryarrange(m,n,q);
	pressure5=memoryarrange(m,n,q);
	pressure6=memoryarrange(m,n,q);
	pressure7=memoryarrange(m,n,q);
	pressure8=memoryarrange(m,n,q);
	pressure9=memoryarrange(m,n,q);
	pressure10=memoryarrange(m,n,q);
	pressure11=memoryarrange(m,n,q);
	pressure12=memoryarrange(m,n,q);
	pressure13=memoryarrange(m,n,q);
	pressure14=memoryarrange(m,n,q);
	pressure15=memoryarrange(m,n,q);
	pressure16=memoryarrange(m,n,q);
	pressure17=memoryarrange(m,n,q);
	pressure18=memoryarrange(m,n,q);
}

double tNS[]={1.0/3,1.0/18,1.0/18,1.0/18,1.0/18,1.0/18,1.0/18,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36};
double tNSD3Q7[]={1.0/4,1.0/8,1.0/8,1.0/8,1.0/8,1.0/8,1.0/8};
double cxNS[]={0.0,1.0,0.0,0.0,-1.0,0.0,0.0,1.0,1.0,0.0,-1.0,-1.0,0.0,1.0,1.0,0.0,-1.0,-1.0,0.0};
double cxNSD3Q7[]={0.0,1.0,0.0,0.0,-1.0,0.0,0.0};
double cyNS[]={0.0,0.0,1.0,0.0,0.0,-1.0,0.0,1.0,0.0,1.0,-1.0,0.0,-1.0,-1.0,0.0,1.0,1.0,0.0,-1.0};
double cyNSD3Q7[]={0.0,0.0,1.0,0.0,0.0,-1.0,0.0};
double czNS[]={0.0,0.0,0.0,1.0,0.0,0.0,-1.0,0.0,1.0,1.0,0.0,-1.0,-1.0,0.0,-1.0,-1.0,0.0,1.0,1.0};
double czNSD3Q7[]={0.0,0.0,0.0,1.0,0.0,0.0,-1.0};

double tNSD3Q15[] = {2.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/72,1.0/72,1.0/72,1.0/72,1.0/72,1.0/72,1.0/72,1.0/72};
double cxNSD3Q15[] = { 0.0, 1.0,-1.0,0.0, 0.0, 0.0, 0.0,1.0, 1.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0};
double cyNSD3Q15[] = { 0.0, 0.0, 0.0,1.0,-1.0, 0.0, 0.0,1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0};
double czNSD3Q15[] = { 0.0, 0.0, 0.0,0.0, 0.0, 1.0,-1.0,1.0,-1.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0};

double A_g[19][19]={{1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,1.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,1.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,1.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,1.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0/(viscosity*3+0.5),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0/(viscosity*3+0.5),0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.1,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0/(viscosity*3+0.5),0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0/(viscosity*3+0.5),0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0/(viscosity*3+0.5),0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.1,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.1,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.1}};

double A_gConcentration[7][7] = { 
{ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{ 0.0, 1.0 / (concentrationDiffusivity * 4 + 0.5), 0.0, 0.0, 0.0, 0.0, 0.0 },
{ 0.0, 0.0, 1.0 / (concentrationDiffusivity * 4 + 0.5), 0.0, 0.0, 0.0, 0.0 },
{ 0.0, 0.0, 0.0, 1.0 / (concentrationDiffusivity * 4 + 0.5), 0.0, 0.0, 0.0 },
{ 0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0},
{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0},
{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9}
};

double A_gConcentrationD3Q15[15][15] = {
	{1.0, 0.0, 0.0, 0.0,                                                  0.0,                                                         0.0,                                                 0.0,                                                         0.0,                                                 0.0,                                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 1.0, 0.0, 0.0,                                                  0.0,                                                         0.0,                                                 0.0,                                                         0.0,                                                 0.0,                                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 1.0, 0.0,                                                  0.0,                                                         0.0,                                                 0.0,                                                         0.0,                                                 0.0,                                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 1.0 / (concentrationDiffusivity * 3.0 + 0.5), 1.0 / (concentrationDiffusivity * 3.0 + 0.5)/2.0-1.0,0.0,                                                 0.0,                                                         0.0,                                                 0.0,                                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0,                                                  1.0,                                                         0.0,                                                 0.0,                                                         0.0,                                                 0.0,                                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0,                                                  0.0,                                                         1.0 / (concentrationDiffusivity * 3.0 + 0.5),1.0 / (concentrationDiffusivity * 3.0 + 0.5)/2.0-1.0,0.0,                                                 0.0,                                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0,                                                  0.0,                                                         0.0,                                                 1.0,                                                         0.0,                                                 0.0,                                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0,                                                  0.0,                                                         0.0,                                                 0.0,                                                         1.0 / (concentrationDiffusivity * 3.0 + 0.5),1.0 / (concentrationDiffusivity * 3.0 + 0.5)/2.0-1.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0,                                                  0.0,                                                         0.0,                                                 0.0,                                                         0.0,                                                 1.0,                                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0,                                                  0.0,                                                         0.0,                                                 0.0,                                                         0.0,                                                 0.0,                                                         1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0,                                                  0.0,                                                         0.0,                                                 0.0,                                                         0.0,                                                 0.0,                                                         0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0,                                                  0.0,                                                         0.0,                                                 0.0,                                                         0.0,                                                 0.0,                                                         0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0,                                                  0.0,                                                         0.0,                                                 0.0,                                                         0.0,                                                 0.0,                                                         0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0,                                                  0.0,                                                         0.0,                                                 0.0,                                                         0.0,                                                 0.0,                                                         0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
	{0.0, 0.0, 0.0, 0.0,                                                  0.0,                                                         0.0,                                                 0.0,                                                         0.0,                                                 0.0,                                                         0.0, 0.0, 0.0, 0.0, 0.0, 1.0}
};

double A=1.0/(viscosity*3+0.5);

double M[19][19]={{1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0},
{-30.0,-11.0,-11.0,-11.0,-11.0,-11.0,-11.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0},
{12.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0},
{0.0,1.0,0.0,0.0,-1.0,0.0,0.0,1.0,1.0,0.0,-1.0,-1.0,0.0,1.0,1.0,0.0,-1.0,-1.0,0.0},
{0.0,-4.0,0.0,0.0,4.0,0.0,0.0,1.0,1.0,0.0,-1.0,-1.0,0.0,1.0,1.0,0.0,-1.0,-1.0,0.0},
{0.0,0.0,1.0,0.0,0.0,-1.0,0.0,1.0,0.0,1.0,-1.0,0.0,-1.0,-1.0,0.0,1.0,1.0,0.0,-1.0},
{0.0,0.0,-4.0,0.0,0.0,4.0,0.0,1.0,0.0,1.0,-1.0,0.0,-1.0,-1.0,0.0,1.0,1.0,0.0,-1.0},
{0.0,0.0,0.0,1.0,0.0,0.0,-1.0,0.0,1.0,1.0,0.0,-1.0,-1.0,0.0,-1.0,-1.0,0.0,1.0,1.0},
{0.0,0.0,0.0,-4.0,0.0,0.0,4.0,0.0,1.0,1.0,0.0,-1.0,-1.0,0.0,-1.0,-1.0,0.0,1.0,1.0},
{0.0,2.0,-1.0,-1.0,2.0,-1.0,-1.0,1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0},
{0.0,-4.0,2.0,2.0,-4.0,2.0,2.0,1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0,1.0,1.0,-2.0},
{0.0,0.0,1.0,-1.0,0.0,1.0,-1.0,1.0,-1.0,0.0,1.0,-1.0,0.0,1.0,-1.0,0.0,1.0,-1.0,0.0},
{0.0,0.0,-2.0,2.0,0.0,-2.0,2.0,1.0,-1.0,0.0,1.0,-1.0,0.0,1.0,-1.0,0.0,1.0,-1.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,-1.0,0.0,0.0,-1.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,-1.0,0.0,0.0,-1.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,-1.0,0.0,0.0,-1.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,-1.0,0.0,-1.0,1.0,0.0,1.0,-1.0,0.0,-1.0,1.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,1.0,0.0,-1.0,1.0,0.0,1.0,-1.0,0.0,-1.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,-1.0,0.0,-1.0,1.0,0.0,-1.0,1.0,0.0,1.0,-1.0}};

double M_Concentration[7][7] = {
{1.0,1.0,1.0,1.0,1.0,1.0,1.0},
{0.0,1.0,0.0,0.0,-1.0,0.0,0.0},
{0.0,0.0,1.0,0.0,0.0,-1.0,0.0},
{0.0,0.0,0.0,1.0,0.0,0.0,-1.0},
{6.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0},
{0.0,2.0,-1.0,-1.0,2.0,-1.0,-1.0},
{0.0,0.0,1.0,-1.0,0.0,1.0,-1.0}
};

double M_ConcentrationD3Q15[15][15] = {
	{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
	{-2.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
	{16.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
	{ 0.0, 1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0},
	{ 0.0,-4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0},
	{ 0.0, 0.0, 0.0, 1.0,-1.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0},
	{ 0.0, 0.0, 0.0,-4.0, 4.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0},
	{ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0, 1.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0},
	{ 0.0, 0.0, 0.0, 0.0, 0.0,-4.0, 4.0, 1.0, 1.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0},
	{ 0.0, 2.0, 2.0,-1.0,-1.0,-1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{ 0.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0},
	{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0, 1.0, 1.0},
	{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0,-1.0, 1.0,-1.0, 1.0},
	{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0,-1.0, 1.0,-1.0, 1.0, 1.0,-1.0}
};

double MInverse_Concentration[7][7] = {
{0.14285714285, 0.0, 0.0, 0.0,0.14285714285, 0.0, 0.0 },
{0.14285714285, 0.5, 0.0, 0.0, -0.0238095238, 0.16666667, 0.0 },
{ 0.14285714285, 0.0, 0.5, 0.0, -0.0238095238, -0.083333333, 0.25 },
{ 0.14285714285, 0.0, 0.0, 0.5, -0.0238095238, -0.083333333, -0.25 },
{ 0.14285714285, -0.5, 0.0, 0.0, -0.0238095238, 0.16666667, 0.0 },
{ 0.14285714285, 0.0, -0.5, 0.0, -0.0238095238, -0.083333333, 0.25 },
{ 0.14285714285, 0.0, 0.0, -0.5, -0.0238095238, -0.08333333, -0.25 }
};

double MInverse_ConcentrationD3Q15[15][15] = {
	{ 0.0666667, -0.1111, 0.0444, 0.0, -3.4694e-18,-2.7756e-18, 2.7756e-18, 2.7756e-18, 6.9389e-19, 6.9389e-18, 7.4015e-18, 0.0,       -3.4694e-18, 0.0,        0.0},
	{ 0.0666667, -0.0556,-0.0111, 0.1, -0.1,       -3.3307e-18, 2.6368e-18, 3.3307e-18, 8.3267e-19, 0.1667,     7.4015e-18, 2.7756e-18,-4.1633e-18,-2.7756e-18, 2.7756e-18},
	{ 0.0666667, -0.0556,-0.0111,-0.1,  0.1,       -2.2204e-18, 2.9143e-18, 2.2204e-18, 5.5511e-19, 0.1667,     7.4015e-18,-2.7756e-18,-2.7756e-18, 2.7756e-18,-2.7756e-18},
	{ 0.0666667, -0.0556,-0.0111, 0.0,  0.0,        0.1,       -0.1,        5.5511e-18, 5.0885e-18,-0.0833,     0.25,       0.0,       -6.9389e-18, 0.0,        0.0},
	{ 0.0666667, -0.0556,-0.0111, 0.0,  0.0,       -0.1,        0.1,        0.0,        3.7007e-18,-0.0833,     0.25,       0.0,        0.0,        0.0,        0.0},
	{ 0.0666667, -0.0556,-0.0111, 0.0,  0.0,       -2.2204e-18,-5.5511e-19, 0.1,       -0.1,       -0.0833,    -0.25,       0.0,        0.0,        0.0,        0.0},
	{ 0.0666667, -0.0556,-0.0111, 0.0,  0.0,       -8.8818e-18,-2.2204e-18,-0.1,        0.1,       -0.0833,    -0.25,       0.0,        0.0,        0.0,        0.0},
	{ 0.0666667,  0.0556, 0.0028, 0.1,  0.025,      0.1,        0.025,      0.1,        0.025,     -6.9389e-18, 7.4015e-18, 0.125,      0.125,      0.125,      0.125},
	{ 0.0666667,  0.0556, 0.0028,-0.1, -0.025,      0.1,        0.025,      0.1,        0.025,      0.0,        7.4015e-18,-0.125,      0.125,     -0.125,     -0.125},
	{ 0.0666667,  0.0556, 0.0028, 0.1,  0.025,     -0.1,       -0.025,      0.1,        0.025,      0.0,        7.4015e-18,-0.125,     -0.125,      0.125,     -0.125},
	{ 0.0666667,  0.0556, 0.0028,-0.1, -0.025,     -0.1,       -0.025,      0.1,        0.025,     -6.9389e-18, 7.4015e-18, 0.125,     -0.125,     -0.125,      0.125},
	{ 0.0666667,  0.0556, 0.0028, 0.1,  0.025,      0.1,        0.025,     -0.1,       -0.025,      0.0,        7.4015e-18, 0.125,     -0.125,     -0.125,     -0.125},
	{ 0.0666667,  0.0556, 0.0028,-0.1, -0.025,      0.1,        0.025,     -0.1,       -0.025,      0.0,        7.4015e-18,-0.125,     -0.125,      0.125,      0.125},
	{ 0.0666667,  0.0556, 0.0028, 0.1,  0.025,     -0.1,       -0.025,     -0.1,       -0.025,     -6.9389e-18, 7.4015e-18,-0.125,      0.125,     -0.125,      0.125},
	{ 0.0666667,  0.0556, 0.0028,-0.1, -0.025,     -0.1,       -0.025,     -0.1,       -0.025,      0.0,        7.4015e-18, 0.125,      0.125,      0.125,     -0.125}
};

double MInverse[19][19]={{0.0526315789473684,	-0.0125313283208020,	0.0476190476190476,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0},
{0.0526315789473684,	-0.00459482038429407,	-0.0158730158730159,	0.100000000000000,	-0.100000000000000,	8.32667268468868e-18,	2.08166817117217e-18,	1.38777878078145e-17,	3.46944695195361e-18,	0.0555555555555555,	-0.0555555555555556,	-4.62592926927146e-19,	-2.31296463463573e-19,	0.0,	2.77555756156289e-17,	0.0,	1.38777878078145e-17,	-2.77555756156289e-17,	1.38777878078145e-17},
{0.0526315789473684,	-0.00459482038429407,	-0.0158730158730159,	3.33066907387547e-18,	8.32667268468867e-19,	0.100000000000000,	-0.100000000000000,	2.77555756156289e-18,	6.93889390390723e-19,	-0.0277777777777778, 0.0277777777777778,	0.0833333333333333,	-0.0833333333333334,	-5.55111512312578e-18,	-6.93889390390723e-18,	1.38777878078145e-17,	-9.71445146547012e-18,	6.93889390390724e-19,	1.04083408558608e-17},
{0.0526315789473684,	-0.00459482038429407,	-0.0158730158730159,	0.0,	0.0,	-1.38777878078145e-17,	-6.93889390390723e-18,	0.100000000000000,	-0.100000000000000,	-0.0277777777777778, 0.0277777777777778,	-0.0833333333333333,	0.0833333333333333,	0.0,	2.08166817117217e-17,	0.0,	0.0,	-3.46944695195361e-18,	1.04083408558608e-17},
{0.0526315789473684,	-0.00459482038429407,	-0.0158730158730159,	-0.100000000000000,	0.100000000000000,	0.0,	1.92592994438724e-34,	2.77555756156289e-18,	6.93889390390723e-19,	0.0555555555555556,	-0.0555555555555556,	5.08852219619864e-18,	2.54426109809932e-18,	6.93889390390723e-18,	0.0,	0.0,	3.46944695195361e-18,	3.46944695195361e-18,	0.0},
{0.0526315789473684,	-0.00459482038429407,	-0.0158730158730159,	2.22044604925031e-18,	5.55111512312578e-19,	-0.100000000000000,	0.100000000000000,	-6.24500451351651e-18,	-8.50014503228636e-18,	-0.0277777777777778,	0.0277777777777778,	0.0833333333333333,	-0.0833333333333333,	5.55111512312578e-18,	-8.67361737988404e-18,	0.0,	2.77555756156289e-18,	1.73472347597680e-19,	8.67361737988404e-19},
{0.0526315789473684,	-0.00459482038429407,	-0.0158730158730159,	0.0,	0.0,	-8.32667268468868e-18,	3.46944695195361e-18,	-0.100000000000000,	0.100000000000000,	-0.0277777777777778,	0.0277777777777778,	-0.0833333333333333,	0.0833333333333333,	0.0,	-6.93889390390723e-18,	0.0,	0.0,	3.46944695195361e-18,	-3.46944695195361e-18},
{0.0526315789473684,	0.00334168755221387,	0.00396825396825396,	0.100000000000000,	0.0250000000000000,	0.100000000000000,	0.0250000000000000,	-2.77555756156289e-17,	-6.93889390390723e-18,	0.0277777777777778,	0.0138888888888889,	0.0833333333333333,	0.0416666666666667,	0.250000000000000,	2.77555756156289e-17,	2.77555756156289e-17,	0.125000000000000,	-0.125000000000000,	0.0},
{0.0526315789473684,	0.00334168755221387,	0.00396825396825396,	0.100000000000000,	0.0250000000000000,	0.0,	0.0,	0.100000000000000,	0.0250000000000000,	0.0277777777777778,	0.0138888888888889,	-0.0833333333333333,	-0.0416666666666667,	0.0,	0.0,	0.250000000000000,	-0.125000000000000,	0.0,	0.125000000000000},
{0.0526315789473684,	0.00334168755221387,	0.00396825396825397,	0.0,	0.0,	0.100000000000000,	0.0250000000000000,	0.100000000000000,	0.0250000000000000,	-0.0555555555555555,	-0.0277777777777778,	-1.38777878078145e-17,	-6.93889390390723e-18,	0.0,	0.250000000000000,	0.0,	0.0,	0.125000000000000,	-0.125000000000000},
{0.0526315789473685,	0.00334168755221387,	0.00396825396825397,	-0.100000000000000,	-0.0250000000000000,	-0.100000000000000,	-0.0250000000000000,	2.77555756156289e-17,	6.93889390390723e-18,	0.0277777777777778,	0.0138888888888889,	0.0833333333333334,	0.0416666666666667,	0.250000000000000,	-2.77555756156289e-17,	-2.77555756156289e-17,	-0.125000000000000,	0.125000000000000,	0.0},
{0.0526315789473684,	0.00334168755221387,	0.00396825396825397,	-0.100000000000000,	-0.0250000000000000,	0.0,	0.0,	-0.100000000000000,	-0.0250000000000000,	0.0277777777777778,	0.0138888888888889,	-0.0833333333333333,	-0.0416666666666667,	0.0,	0.0,	0.250000000000000,	0.125000000000000,	0.0,	-0.125000000000000},
{0.0526315789473684,	0.00334168755221387,	0.00396825396825397,	0.0,	0.0,	-0.100000000000000,	-0.0250000000000000,	-0.100000000000000,	-0.0250000000000000,	-0.0555555555555556,	-0.0277777777777778,	1.38777878078145e-17,	6.93889390390723e-18,	0.0,	0.250000000000000,	0.0,	0.0,	-0.125000000000000,	0.125000000000000},
{0.0526315789473684,	0.00334168755221387,	0.00396825396825397,	0.100000000000000,	0.0250000000000000,	-0.100000000000000,	-0.0250000000000000,	0.0,	0.0,	0.0277777777777778,	0.0138888888888889,	0.0833333333333333,	0.0416666666666667,	-0.250000000000000,	0.0,	0.0,	0.125000000000000,	0.125000000000000,	0.0},
{0.0526315789473684,	0.00334168755221387,	0.00396825396825397,	0.100000000000000,	0.0250000000000000,	0.0,	0.0,	-0.100000000000000,	-0.0250000000000000,	0.0277777777777778,	0.0138888888888889,	-0.0833333333333333,	-0.0416666666666667,	0.0,	0.0,	-0.250000000000000,	-0.125000000000000,	0.0,	-0.125000000000000},
{0.0526315789473684,	0.00334168755221387,	0.00396825396825397,	0.0,	0.0,	0.100000000000000,	0.0250000000000000,	-0.100000000000000,	-0.0250000000000000,	-0.0555555555555556,	-0.0277777777777778,	-1.38777878078145e-17,	-6.93889390390723e-18,	0.0,	-0.250000000000000,	0.0,	0.0,	0.125000000000000,	0.125000000000000},
{0.0526315789473684,	0.00334168755221387,	0.00396825396825396,	-0.100000000000000,	-0.0250000000000000,	0.100000000000000,	0.0250000000000000,	0.0,	0.0,	0.0277777777777778,	0.0138888888888889,	0.0833333333333333,	0.0416666666666667,	-0.250000000000000,	0.0,	0.0,	-0.125000000000000,	-0.125000000000000,	0.0},
{0.0526315789473684,	0.00334168755221387,	0.00396825396825397,	-0.100000000000000,	-0.0250000000000000,	0.0,	-6.93889390390723e-18,	0.100000000000000,	0.0250000000000000,	0.0277777777777778,	0.0138888888888889,	-0.0833333333333333,	-0.0416666666666667,	0.0,	0.0,	-0.250000000000000,	0.125000000000000,	0.0,	0.125000000000000},
{0.0526315789473684,	0.00334168755221387,	0.00396825396825397,	0.0,	0.0,	-0.100000000000000,	-0.0250000000000000,	0.100000000000000,	0.0250000000000000,	-0.0555555555555555,	-0.0277777777777778,	1.38777878078145e-17,	6.93889390390723e-18,	0.0,	-0.250000000000000,	0.0,	0.0,	-0.125000000000000,	-0.125000000000000}};

double I[19][19]={{1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0}};

double mfIn0;
double mfIn1;
double mfIn2;
double mfIn3;
double mfIn4;
double mfIn5;
double mfIn6;
double mfIn7;
double mfIn8;
double mfIn9;
double mfIn10;
double mfIn11;
double mfIn12;
double mfIn13;
double mfIn14;
double mfIn15;
double mfIn16;
double mfIn17;
double mfIn18;
double mfEq0;
double mfEq1;
double mfEq2;
double mfEq3;
double mfEq4;
double mfEq5;
double mfEq6;
double mfEq7;
double mfEq8;
double mfEq9;
double mfEq10;
double mfEq11;
double mfEq12;
double mfEq13;
double mfEq14;
double mfEq15;
double mfEq16;
double mfEq17;
double mfEq18;
double mForce0;
double mForce1;
double mForce2;
double mForce3;
double mForce4;
double mForce5;
double mForce6;
double mForce7;
double mForce8;
double mForce9;
double mForce10;
double mForce11;
double mForce12;
double mForce13;
double mForce14;
double mForce15;
double mForce16;
double mForce17;
double mForce18;
double wallForceX1;
double wallForceY1;
double wallForceZ1;
double wallForceX2;
double wallForceY2;
double wallForceZ2;
double wallForceX3;
double wallForceY3;
double wallForceZ3;
double wallForceX4;
double wallForceY4;
double wallForceZ4;
double wallForceX5;
double wallForceY5;
double wallForceZ5;
double wallForceX6;
double wallForceY6;
double wallForceZ6;
double wallForceX7;
double wallForceY7;
double wallForceZ7;
double wallForceX8;
double wallForceY8;
double wallForceZ8;
double wallForceX9;
double wallForceY9;
double wallForceZ9;
double wallForceX10;
double wallForceY10;
double wallForceZ10;
double wallForceX11;
double wallForceY11;
double wallForceZ11;
double wallForceX12;
double wallForceY12;
double wallForceZ12;
double wallForceX13;
double wallForceY13;
double wallForceZ13;
double wallForceX14;
double wallForceY14;
double wallForceZ14;
double wallForceX15;
double wallForceY15;
double wallForceZ15;
double wallForceX16;
double wallForceY16;
double wallForceZ16;
double wallForceX17;
double wallForceY17;
double wallForceZ17;
double wallForceX18;
double wallForceY18;
double wallForceZ18;
double Factor0;
double Factor1;
double Factor2;
double Factor3;
double Factor4;
double Factor5;
double Factor6;
double Factor7;
double Factor8;
double Factor9;
double Factor10;
double Factor11;
double Factor12;
double Factor13;
double Factor14;
double Factor15;
double Factor16;
double Factor17;
double Factor18;
double Psi1;
double Psi2;
double Psi3;
double Psi4;
double Psi5;
double Psi6;
double Psi7;
double Psi8;
double Psi9;
double Psi10;
double Psi11;
double Psi12;
double Psi13;
double Psi14;
double Psi15;
double Psi16;
double Psi17;
double Psi18;
double sign1;
double sign2;
double sign3;
double sign4;
double sign5;
double sign6;
double sign7;
double sign8;
double sign9;
double sign10;
double sign11;
double sign12;
double sign13;
double sign14;
double sign15;
double sign16;
double sign17;
double sign18;
/******************************************************************************************************************************************/
void FieldInitialGasTwoPhase(int m, int n, int q, double GAS, double SOLID, double b, double R, double T, double a){
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
						rho[i][j][k]=GAS;
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
						Fx[i][j][k]=0.0;
						Fy[i][j][k]=0.0;
						Fz[i][j][k]=0.0;
						ux[i][j][k]=uTotX+Fx[i][j][k]/2.0/rho[i][j][k];
						uy[i][j][k]=uTotY+Fy[i][j][k]/2.0/rho[i][j][k];
						uz[i][j][k]=uTotZ+Fz[i][j][k]/2.0/rho[i][j][k];
					}
					if (domain[i][j][k]==1){
						rho[i][j][k]=SOLID;
						pressure[i][j][k]=0;
						ux[i][j][k]=0;
						uy[i][j][k]=0;
						uz[i][j][k]=0;
						fIn0[i][j][k]=0.0;
						fIn1[i][j][k]=0.0;
						fIn2[i][j][k]=0.0;
						fIn3[i][j][k]=0.0;
						fIn4[i][j][k]=0.0;
						fIn5[i][j][k]=0.0;
						fIn6[i][j][k]=0.0;
						fIn7[i][j][k]=0.0;
						fIn8[i][j][k]=0.0;
						fIn9[i][j][k]=0.0;
						fIn10[i][j][k]=0.0;
						fIn11[i][j][k]=0.0;
						fIn12[i][j][k]=0.0;
						fIn13[i][j][k]=0.0;
						fIn14[i][j][k]=0.0;
						fIn15[i][j][k]=0.0;
						fIn16[i][j][k]=0.0;
						fIn17[i][j][k]=0.0;
						fIn18[i][j][k]=0.0;
						Fx[i][j][k]=0.0;
						Fy[i][j][k]=0.0;
						Fz[i][j][k]=0.0;
					}
				}
			}
		}
	}
}
void FieldInitialLiquidTwoPhase(int m, int n, int q, double FLUID, double GAS, int WIDE, int radius, double b, double R, double T, double a){
	int i;
	int j;
	int k;
	double initialLiquid;
	double distance;
	double bRho4;
	double BRho4;
#pragma omp parallel private(i,j,k,bRho4,BRho4,initialLiquid,distance)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0){
//						distance=sqrt((i-m/2.0)*(i-m/2.0)+(j-n/2.0)*(j-n/2.0)+(k-q/2.0)*(k-q/2.0));
						initialLiquid=2.0*(abs(k-0.0)-WIDE)/5;
//						initialLiquid=2.0*(distance-radius)/5;
						rho[i][j][k]=(FLUID+GAS)/2-(FLUID-GAS)/2*tanh(initialLiquid);
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

void wallForceCalculationTwoPhase(double wallDensity){
	wallForceX1=wallDensity*tNS[1]*(-cxNS[1]);
	wallForceY1=wallDensity*tNS[1]*(-cyNS[1]);
	wallForceZ1=wallDensity*tNS[1]*(-czNS[1]);
	wallForceX2=wallDensity*tNS[2]*(-cxNS[2]);
	wallForceY2=wallDensity*tNS[2]*(-cyNS[2]);
	wallForceZ2=wallDensity*tNS[2]*(-czNS[2]);
	wallForceX3=wallDensity*tNS[3]*(-cxNS[3]);
	wallForceY3=wallDensity*tNS[3]*(-cyNS[3]);
	wallForceZ3=wallDensity*tNS[3]*(-czNS[3]);
	wallForceX4=wallDensity*tNS[4]*(-cxNS[4]);
	wallForceY4=wallDensity*tNS[4]*(-cyNS[4]);
	wallForceZ4=wallDensity*tNS[4]*(-czNS[4]);
	wallForceX5=wallDensity*tNS[5]*(-cxNS[5]);
	wallForceY5=wallDensity*tNS[5]*(-cyNS[5]);
	wallForceZ5=wallDensity*tNS[5]*(-czNS[5]);
	wallForceX6=wallDensity*tNS[6]*(-cxNS[6]);
	wallForceY6=wallDensity*tNS[6]*(-cyNS[6]);
	wallForceZ6=wallDensity*tNS[6]*(-czNS[6]);
	wallForceX7=wallDensity*tNS[7]*(-cxNS[7]);
	wallForceY7=wallDensity*tNS[7]*(-cyNS[7]);
	wallForceZ7=wallDensity*tNS[7]*(-czNS[7]);
	wallForceX8=wallDensity*tNS[8]*(-cxNS[8]);
	wallForceY8=wallDensity*tNS[8]*(-cyNS[8]);
	wallForceZ8=wallDensity*tNS[8]*(-czNS[8]);
	wallForceX9=wallDensity*tNS[9]*(-cxNS[9]);
	wallForceY9=wallDensity*tNS[9]*(-cyNS[9]);
	wallForceZ9=wallDensity*tNS[9]*(-czNS[9]);
	wallForceX10=wallDensity*tNS[10]*(-cxNS[10]);
	wallForceY10=wallDensity*tNS[10]*(-cyNS[10]);
	wallForceZ10=wallDensity*tNS[10]*(-czNS[10]);
	wallForceX11=wallDensity*tNS[11]*(-cxNS[11]);
	wallForceY11=wallDensity*tNS[11]*(-cyNS[11]);
	wallForceZ11=wallDensity*tNS[11]*(-czNS[11]);
	wallForceX12=wallDensity*tNS[12]*(-cxNS[12]);
	wallForceY12=wallDensity*tNS[12]*(-cyNS[12]);
	wallForceZ12=wallDensity*tNS[12]*(-czNS[12]);
	wallForceX13=wallDensity*tNS[13]*(-cxNS[13]);
	wallForceY13=wallDensity*tNS[13]*(-cyNS[13]);
	wallForceZ13=wallDensity*tNS[13]*(-czNS[13]);
	wallForceX14=wallDensity*tNS[14]*(-cxNS[14]);
	wallForceY14=wallDensity*tNS[14]*(-cyNS[14]);
	wallForceZ14=wallDensity*tNS[14]*(-czNS[14]);
	wallForceX15=wallDensity*tNS[15]*(-cxNS[15]);
	wallForceY15=wallDensity*tNS[15]*(-cyNS[15]);
	wallForceZ15=wallDensity*tNS[15]*(-czNS[15]);
	wallForceX16=wallDensity*tNS[16]*(-cxNS[16]);
	wallForceY16=wallDensity*tNS[16]*(-cyNS[16]);
	wallForceZ16=wallDensity*tNS[16]*(-czNS[16]);
	wallForceX17=wallDensity*tNS[17]*(-cxNS[17]);
	wallForceY17=wallDensity*tNS[17]*(-cyNS[17]);
	wallForceZ17=wallDensity*tNS[17]*(-czNS[17]);
	wallForceX18=wallDensity*tNS[18]*(-cxNS[18]);
	wallForceY18=wallDensity*tNS[18]*(-cyNS[18]);
	wallForceZ18=wallDensity*tNS[18]*(-czNS[18]);
}

void FactorCalculationTwoPhase(){
	Factor0=(I[0][0]-A_g[0][0]/2.0);
	Factor1=(I[1][1]-A_g[1][1]/2.0);
	Factor2=(I[2][2]-A_g[2][2]/2.0);
	Factor3=(I[3][3]-A_g[3][3]/2.0);
	Factor4=(I[4][4]-A_g[4][4]/2.0);
	Factor5=(I[5][5]-A_g[5][5]/2.0);
	Factor6=(I[6][6]-A_g[6][6]/2.0);
	Factor7=(I[7][7]-A_g[7][7]/2.0);
	Factor8=(I[8][8]-A_g[8][8]/2.0);
	Factor9=(I[9][9]-A_g[9][9]/2.0);
	Factor10=(I[10][10]-A_g[10][10]/2.0);
	Factor11=(I[11][11]-A_g[11][11]/2.0);
	Factor12=(I[12][12]-A_g[12][12]/2.0);
	Factor13=(I[13][13]-A_g[13][13]/2.0);
	Factor14=(I[14][14]-A_g[14][14]/2.0);
	Factor15=(I[15][15]-A_g[15][15]/2.0);
	Factor16=(I[16][16]-A_g[16][16]/2.0);
	Factor17=(I[17][17]-A_g[17][17]/2.0);
	Factor18=(I[18][18]-A_g[18][18]/2.0);
}

void bodyForceTwoPhase(int m, int n, int q, double G, double boundary){
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
						bodyForce[i][j][k]=G;
					}
				}
			}
		}
	}
}

void ForceCalculationTwoPhase(int m, int n, int q, double PHASE1){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k,sign1,sign2,sign3,sign4,sign5,sign6,sign7,sign8,sign9,sign10,sign11,sign12,sign13,sign14,sign15,sign16,sign17,sign18,\
	Psi1,Psi2,Psi3,Psi4,Psi5,Psi6,Psi7,Psi8,Psi9,Psi10,Psi11,Psi12,Psi13,Psi14,Psi15,Psi16,Psi17,Psi18)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0){
						if (pressure1[i][j][k]-1.0/3*rho1[i][j][k]>0)
							sign1=1.0;
						else sign1=-1.0;
						if (pressure2[i][j][k]-1.0/3*rho2[i][j][k]>0)
							sign2=1.0;
						else sign2=-1.0;
						if (pressure3[i][j][k]-1.0/3*rho3[i][j][k]>0)
							sign3=1.0;
						else sign3=-1.0;
						if (pressure4[i][j][k]-1.0/3*rho4[i][j][k]>0)
							sign4=1.0;
						else sign4=-1.0;
						if (pressure5[i][j][k]-1.0/3*rho5[i][j][k]>0)
							sign5=1.0;
						else sign5=-1.0;
						if (pressure6[i][j][k]-1.0/3*rho6[i][j][k]>0)
							sign6=1.0;
						else sign6=-1.0;
						if (pressure7[i][j][k]-1.0/3*rho7[i][j][k]>0)
							sign7=1.0;
						else sign7=-1.0;
						if (pressure8[i][j][k]-1.0/3*rho8[i][j][k]>0)
							sign8=1.0;
						else sign8=-1.0;
						if (pressure9[i][j][k]-1.0/3*rho9[i][j][k]>0)
							sign9=1.0;
						else sign9=-1.0;
						if (pressure10[i][j][k]-1.0/3*rho10[i][j][k]>0)
							sign10=1.0;
						else sign10=-1.0;
						if (pressure11[i][j][k]-1.0/3*rho11[i][j][k]>0)
							sign11=1.0;
						else sign11=-1.0;
						if (pressure12[i][j][k]-1.0/3*rho12[i][j][k]>0)
							sign12=1.0;
						else sign12=-1.0;
						if (pressure13[i][j][k]-1.0/3*rho13[i][j][k]>0)
							sign13=1.0;
						else sign13=-1.0;
						if (pressure14[i][j][k]-1.0/3*rho14[i][j][k]>0)
							sign14=1.0;
						else sign14=-1.0;
						if (pressure15[i][j][k]-1.0/3*rho15[i][j][k]>0)
							sign15=1.0;
						else sign15=-1.0;
						if (pressure16[i][j][k]-1.0/3*rho16[i][j][k]>0)
							sign16=1.0;
						else sign16=-1.0;
						if (pressure17[i][j][k]-1.0/3*rho17[i][j][k]>0)
							sign17=1.0;
						else sign17=-1.0;
						if (pressure18[i][j][k]-1.0/3*rho18[i][j][k]>0)
							sign18=1.0;
						else sign18=-1.0;
						Psi1=sign1*sqrt(2.0*fabs(pressure1[i][j][k]-1.0/3*rho1[i][j][k])/(6.0*1/18.0))*tNS[1];
						Psi2=sign2*sqrt(2.0*fabs(pressure2[i][j][k]-1.0/3*rho2[i][j][k])/(6.0*1/18.0))*tNS[2];
						Psi3=sign3*sqrt(2.0*fabs(pressure3[i][j][k]-1.0/3*rho3[i][j][k])/(6.0*1/18.0))*tNS[3];
						Psi4=sign4*sqrt(2.0*fabs(pressure4[i][j][k]-1.0/3*rho4[i][j][k])/(6.0*1/18.0))*tNS[4];
						Psi5=sign5*sqrt(2.0*fabs(pressure5[i][j][k]-1.0/3*rho5[i][j][k])/(6.0*1/18.0))*tNS[5];
						Psi6=sign6*sqrt(2.0*fabs(pressure6[i][j][k]-1.0/3*rho6[i][j][k])/(6.0*1/18.0))*tNS[6];
						Psi7=sign7*sqrt(2.0*fabs(pressure7[i][j][k]-1.0/3*rho7[i][j][k])/(6.0*1/18.0))*tNS[7];
						Psi8=sign8*sqrt(2.0*fabs(pressure8[i][j][k]-1.0/3*rho8[i][j][k])/(6.0*1/18.0))*tNS[8];
						Psi9=sign9*sqrt(2.0*fabs(pressure9[i][j][k]-1.0/3*rho9[i][j][k])/(6.0*1/18.0))*tNS[9];
						Psi10=sign10*sqrt(2.0*fabs(pressure10[i][j][k]-1.0/3*rho10[i][j][k])/(6.0*1/18.0))*tNS[10];
						Psi11=sign11*sqrt(2.0*fabs(pressure11[i][j][k]-1.0/3*rho11[i][j][k])/(6.0*1/18.0))*tNS[11];
						Psi12=sign12*sqrt(2.0*fabs(pressure12[i][j][k]-1.0/3*rho12[i][j][k])/(6.0*1/18.0))*tNS[12];
						Psi13=sign13*sqrt(2.0*fabs(pressure13[i][j][k]-1.0/3*rho13[i][j][k])/(6.0*1/18.0))*tNS[13];
						Psi14=sign14*sqrt(2.0*fabs(pressure14[i][j][k]-1.0/3*rho14[i][j][k])/(6.0*1/18.0))*tNS[14];
						Psi15=sign15*sqrt(2.0*fabs(pressure15[i][j][k]-1.0/3*rho15[i][j][k])/(6.0*1/18.0))*tNS[15];
						Psi16=sign16*sqrt(2.0*fabs(pressure16[i][j][k]-1.0/3*rho16[i][j][k])/(6.0*1/18.0))*tNS[16];
						Psi17=sign17*sqrt(2.0*fabs(pressure17[i][j][k]-1.0/3*rho17[i][j][k])/(6.0*1/18.0))*tNS[17];
						Psi18=sign18*sqrt(2.0*fabs(pressure18[i][j][k]-1.0/3*rho18[i][j][k])/(6.0*1/18.0))*tNS[18];
						if (Domain1[i][j][k]==PHASE1){           //judge the nearest point is wall or fluid space
							Fx[i][j][k]=wallForceX1*gWall[i][j][k];
							Fy[i][j][k]=wallForceY1*gWall[i][j][k];
							Fz[i][j][k]=wallForceZ1*gWall[i][j][k];
						}
						else {
							Fx[i][j][k]=fluidForceX1;
							Fy[i][j][k]=fluidForceY1;
							Fz[i][j][k]=fluidForceZ1;
						}
						if (Domain2[i][j][k]==PHASE1){            //judge the nearest point is wall or fluid space
							Fx[i][j][k]=Fx[i][j][k]+wallForceX2*gWall[i][j][k];
							Fy[i][j][k]=Fy[i][j][k]+wallForceY2*gWall[i][j][k];
							Fz[i][j][k]=Fz[i][j][k]+wallForceZ2*gWall[i][j][k];
						}
						else {
							Fx[i][j][k]=Fx[i][j][k]+fluidForceX2;
							Fy[i][j][k]=Fy[i][j][k]+fluidForceY2;
							Fz[i][j][k]=Fz[i][j][k]+fluidForceZ2;
						}
						if (Domain3[i][j][k]==PHASE1){            //judge the nearest point is wall or fluid space
							Fx[i][j][k]=Fx[i][j][k]+wallForceX3*gWall[i][j][k];
							Fy[i][j][k]=Fy[i][j][k]+wallForceY3*gWall[i][j][k];
							Fz[i][j][k]=Fz[i][j][k]+wallForceZ3*gWall[i][j][k];
						}
						else {
							Fx[i][j][k]=Fx[i][j][k]+fluidForceX3;
							Fy[i][j][k]=Fy[i][j][k]+fluidForceY3;
							Fz[i][j][k]=Fz[i][j][k]+fluidForceZ3;
						}
						if (Domain4[i][j][k]==PHASE1){            //judge the nearest point is wall or fluid space
							Fx[i][j][k]=Fx[i][j][k]+wallForceX4*gWall[i][j][k];
							Fy[i][j][k]=Fy[i][j][k]+wallForceY4*gWall[i][j][k];
							Fz[i][j][k]=Fz[i][j][k]+wallForceZ4*gWall[i][j][k];
						}
						else {
							Fx[i][j][k]=Fx[i][j][k]+fluidForceX4;
							Fy[i][j][k]=Fy[i][j][k]+fluidForceY4;
							Fz[i][j][k]=Fz[i][j][k]+fluidForceZ4;
						}
						if (Domain5[i][j][k]==PHASE1){            //judge the nearest point is wall or fluid space
							Fx[i][j][k]=Fx[i][j][k]+wallForceX5*gWall[i][j][k];
							Fy[i][j][k]=Fy[i][j][k]+wallForceY5*gWall[i][j][k];
							Fz[i][j][k]=Fz[i][j][k]+wallForceZ5*gWall[i][j][k];
						}
						else {
							Fx[i][j][k]=Fx[i][j][k]+fluidForceX5;
							Fy[i][j][k]=Fy[i][j][k]+fluidForceY5;
							Fz[i][j][k]=Fz[i][j][k]+fluidForceZ5;
						}
						if (Domain6[i][j][k]==PHASE1){            //judge the nearest point is wall or fluid space
							Fx[i][j][k]=Fx[i][j][k]+wallForceX6*gWall[i][j][k];
							Fy[i][j][k]=Fy[i][j][k]+wallForceY6*gWall[i][j][k];
							Fz[i][j][k]=Fz[i][j][k]+wallForceZ6*gWall[i][j][k];
						}
						else {
							Fx[i][j][k]=Fx[i][j][k]+fluidForceX6;
							Fy[i][j][k]=Fy[i][j][k]+fluidForceY6;
							Fz[i][j][k]=Fz[i][j][k]+fluidForceZ6;
						}
						if (Domain7[i][j][k]==PHASE1){            //judge the nearest point is wall or fluid space
							Fx[i][j][k]=Fx[i][j][k]+wallForceX7*gWall[i][j][k];
							Fy[i][j][k]=Fy[i][j][k]+wallForceY7*gWall[i][j][k];
							Fz[i][j][k]=Fz[i][j][k]+wallForceZ7*gWall[i][j][k];
						}
						else {
							Fx[i][j][k]=Fx[i][j][k]+fluidForceX7;
							Fy[i][j][k]=Fy[i][j][k]+fluidForceY7;
							Fz[i][j][k]=Fz[i][j][k]+fluidForceZ7;
						}
						if (Domain8[i][j][k]==PHASE1){            //judge the nearest point is wall or fluid space
							Fx[i][j][k]=Fx[i][j][k]+wallForceX8*gWall[i][j][k];
							Fy[i][j][k]=Fy[i][j][k]+wallForceY8*gWall[i][j][k];
							Fz[i][j][k]=Fz[i][j][k]+wallForceZ8*gWall[i][j][k];
						}
						else {
							Fx[i][j][k]=Fx[i][j][k]+fluidForceX8;
							Fy[i][j][k]=Fy[i][j][k]+fluidForceY8;
							Fz[i][j][k]=Fz[i][j][k]+fluidForceZ8;
						}
						if (Domain9[i][j][k]==PHASE1){           //judge the nearest point is wall or fluid space
							Fx[i][j][k]=Fx[i][j][k]+wallForceX9*gWall[i][j][k];
							Fy[i][j][k]=Fy[i][j][k]+wallForceY9*gWall[i][j][k];
							Fz[i][j][k]=Fz[i][j][k]+wallForceZ9*gWall[i][j][k];
						}
						else {
							Fx[i][j][k]=Fx[i][j][k]+fluidForceX9;
							Fy[i][j][k]=Fy[i][j][k]+fluidForceY9;
							Fz[i][j][k]=Fz[i][j][k]+fluidForceZ9;
						}
						if (Domain10[i][j][k]==PHASE1){           //judge the nearest point is wall or fluid space
							Fx[i][j][k]=Fx[i][j][k]+wallForceX10*gWall[i][j][k];
							Fy[i][j][k]=Fy[i][j][k]+wallForceY10*gWall[i][j][k];
							Fz[i][j][k]=Fz[i][j][k]+wallForceZ10*gWall[i][j][k];
						}
						else {
							Fx[i][j][k]=Fx[i][j][k]+fluidForceX10;
							Fy[i][j][k]=Fy[i][j][k]+fluidForceY10;
							Fz[i][j][k]=Fz[i][j][k]+fluidForceZ10;
						}
						if (Domain11[i][j][k]==PHASE1){           //judge the nearest point is wall or fluid space
							Fx[i][j][k]=Fx[i][j][k]+wallForceX11*gWall[i][j][k];
							Fy[i][j][k]=Fy[i][j][k]+wallForceY11*gWall[i][j][k];
							Fz[i][j][k]=Fz[i][j][k]+wallForceZ11*gWall[i][j][k];
						}
						else {
							Fx[i][j][k]=Fx[i][j][k]+fluidForceX11;
							Fy[i][j][k]=Fy[i][j][k]+fluidForceY11;
							Fz[i][j][k]=Fz[i][j][k]+fluidForceZ11;
						}
						if (Domain12[i][j][k]==PHASE1){           //judge the nearest point is wall or fluid space
							Fx[i][j][k]=Fx[i][j][k]+wallForceX12*gWall[i][j][k];
							Fy[i][j][k]=Fy[i][j][k]+wallForceY12*gWall[i][j][k];
							Fz[i][j][k]=Fz[i][j][k]+wallForceZ12*gWall[i][j][k];
						}
						else {
							Fx[i][j][k]=Fx[i][j][k]+fluidForceX12;
							Fy[i][j][k]=Fy[i][j][k]+fluidForceY12;
							Fz[i][j][k]=Fz[i][j][k]+fluidForceZ12;
						}
						if (Domain13[i][j][k]==PHASE1){           //judge the nearest point is wall or fluid space
							Fx[i][j][k]=Fx[i][j][k]+wallForceX13*gWall[i][j][k];
							Fy[i][j][k]=Fy[i][j][k]+wallForceY13*gWall[i][j][k];
							Fz[i][j][k]=Fz[i][j][k]+wallForceZ13*gWall[i][j][k];
						}
						else {
							Fx[i][j][k]=Fx[i][j][k]+fluidForceX13;
							Fy[i][j][k]=Fy[i][j][k]+fluidForceY13;
							Fz[i][j][k]=Fz[i][j][k]+fluidForceZ13;
						}
						if (Domain14[i][j][k]==PHASE1){           //judge the nearest point is wall or fluid space
							Fx[i][j][k]=Fx[i][j][k]+wallForceX14*gWall[i][j][k];
							Fy[i][j][k]=Fy[i][j][k]+wallForceY14*gWall[i][j][k];
							Fz[i][j][k]=Fz[i][j][k]+wallForceZ14*gWall[i][j][k];
						}
						else {
							Fx[i][j][k]=Fx[i][j][k]+fluidForceX14;
							Fy[i][j][k]=Fy[i][j][k]+fluidForceY14;
							Fz[i][j][k]=Fz[i][j][k]+fluidForceZ14;
						}
						if (Domain15[i][j][k]==PHASE1){           //judge the nearest point is wall or fluid space
							Fx[i][j][k]=Fx[i][j][k]+wallForceX15*gWall[i][j][k];
							Fy[i][j][k]=Fy[i][j][k]+wallForceY15*gWall[i][j][k];
							Fz[i][j][k]=Fz[i][j][k]+wallForceZ15*gWall[i][j][k];
						}
						else {
							Fx[i][j][k]=Fx[i][j][k]+fluidForceX15;
							Fy[i][j][k]=Fy[i][j][k]+fluidForceY15;
							Fz[i][j][k]=Fz[i][j][k]+fluidForceZ15;
						}
						if (Domain16[i][j][k]==PHASE1){           //judge the nearest point is wall or fluid space
							Fx[i][j][k]=Fx[i][j][k]+wallForceX16*gWall[i][j][k];
							Fy[i][j][k]=Fy[i][j][k]+wallForceY16*gWall[i][j][k];
							Fz[i][j][k]=Fz[i][j][k]+wallForceZ16*gWall[i][j][k];
						}
						else {
							Fx[i][j][k]=Fx[i][j][k]+fluidForceX16;
							Fy[i][j][k]=Fy[i][j][k]+fluidForceY16;
							Fz[i][j][k]=Fz[i][j][k]+fluidForceZ16;
						}
						if (Domain17[i][j][k]==PHASE1){           //judge the nearest point is wall or fluid space
							Fx[i][j][k]=Fx[i][j][k]+wallForceX17*gWall[i][j][k];
							Fy[i][j][k]=Fy[i][j][k]+wallForceY17*gWall[i][j][k];
							Fz[i][j][k]=Fz[i][j][k]+wallForceZ17*gWall[i][j][k];
						}
						else {
							Fx[i][j][k]=Fx[i][j][k]+fluidForceX17;
							Fy[i][j][k]=Fy[i][j][k]+fluidForceY17;
							Fz[i][j][k]=Fz[i][j][k]+fluidForceZ17;
						}
						if (Domain18[i][j][k]==PHASE1){           //judge the nearest point is wall or fluid space
							Fx[i][j][k]=Fx[i][j][k]+wallForceX18*gWall[i][j][k];
							Fy[i][j][k]=Fy[i][j][k]+wallForceY18*gWall[i][j][k];
							Fz[i][j][k]=Fz[i][j][k]+wallForceZ18*gWall[i][j][k];
						}
						else {
							Fx[i][j][k]=Fx[i][j][k]+fluidForceX18;
							Fy[i][j][k]=Fy[i][j][k]+fluidForceY18;
							Fz[i][j][k]=Fz[i][j][k]+fluidForceZ18;
						}
						Fx[i][j][k]=-sqrt(2.0*fabs(pressure[i][j][k]-1.0/3*rho[i][j][k])/(6.0*1/18.0))*Fx[i][j][k];
						Fy[i][j][k]=-sqrt(2.0*fabs(pressure[i][j][k]-1.0/3*rho[i][j][k])/(6.0*1/18.0))*Fy[i][j][k];
						Fz[i][j][k]=-sqrt(2.0*fabs(pressure[i][j][k]-1.0/3*rho[i][j][k])/(6.0*1/18.0))*Fz[i][j][k]+bodyForce[i][j][k];
					}
				}
			}
		}
	}
}

void CollisionSRTTwoPhase(int Length,int Width,int Height){
	int i;
	int j;
	int k;
	double cfNS0;
	double cfNS1;
	double cfNS2;
	double cfNS3;
	double cfNS4;
	double cfNS5;
	double cfNS6;
	double cfNS7;
	double cfNS8;
	double cfNS9;
	double cfNS10;
	double cfNS11;
	double cfNS12;
	double cfNS13;
	double cfNS14;
	double cfNS15;
	double cfNS16;
	double cfNS17;
	double cfNS18;
	double uxFxuyFyuzFz;
#pragma omp parallel private(i,j,k,cfNS0,cfNS1,cfNS2,cfNS3,cfNS4,cfNS5,cfNS6,cfNS7,cfNS8,cfNS9,cfNS10,cfNS11,cfNS12,cfNS13,cfNS14,cfNS15,cfNS16,cfNS17,cfNS18,uxFxuyFyuzFz)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
					if (domain[i][j][k]==0){
						cfNS0=cxNS[0]*Fx[i][j][k]+cyNS[0]*Fy[i][j][k]+czNS[0]*Fz[i][j][k];
						cfNS1=cxNS[1]*Fx[i][j][k]+cyNS[1]*Fy[i][j][k]+czNS[1]*Fz[i][j][k];
						cfNS2=cxNS[2]*Fx[i][j][k]+cyNS[2]*Fy[i][j][k]+czNS[2]*Fz[i][j][k];
						cfNS3=cxNS[3]*Fx[i][j][k]+cyNS[3]*Fy[i][j][k]+czNS[3]*Fz[i][j][k];
						cfNS4=cxNS[4]*Fx[i][j][k]+cyNS[4]*Fy[i][j][k]+czNS[4]*Fz[i][j][k];
						cfNS5=cxNS[5]*Fx[i][j][k]+cyNS[5]*Fy[i][j][k]+czNS[5]*Fz[i][j][k];
						cfNS6=cxNS[6]*Fx[i][j][k]+cyNS[6]*Fy[i][j][k]+czNS[6]*Fz[i][j][k];
						cfNS7=cxNS[7]*Fx[i][j][k]+cyNS[7]*Fy[i][j][k]+czNS[7]*Fz[i][j][k];
						cfNS8=cxNS[8]*Fx[i][j][k]+cyNS[8]*Fy[i][j][k]+czNS[8]*Fz[i][j][k];
						cfNS9=cxNS[9]*Fx[i][j][k]+cyNS[9]*Fy[i][j][k]+czNS[9]*Fz[i][j][k];
						cfNS10=cxNS[10]*Fx[i][j][k]+cyNS[10]*Fy[i][j][k]+czNS[10]*Fz[i][j][k];
						cfNS11=cxNS[11]*Fx[i][j][k]+cyNS[11]*Fy[i][j][k]+czNS[11]*Fz[i][j][k];
						cfNS12=cxNS[12]*Fx[i][j][k]+cyNS[12]*Fy[i][j][k]+czNS[12]*Fz[i][j][k];
						cfNS13=cxNS[13]*Fx[i][j][k]+cyNS[13]*Fy[i][j][k]+czNS[13]*Fz[i][j][k];
						cfNS14=cxNS[14]*Fx[i][j][k]+cyNS[14]*Fy[i][j][k]+czNS[14]*Fz[i][j][k];
						cfNS15=cxNS[15]*Fx[i][j][k]+cyNS[15]*Fy[i][j][k]+czNS[15]*Fz[i][j][k];
						cfNS16=cxNS[16]*Fx[i][j][k]+cyNS[16]*Fy[i][j][k]+czNS[16]*Fz[i][j][k];
						cfNS17=cxNS[17]*Fx[i][j][k]+cyNS[17]*Fy[i][j][k]+czNS[17]*Fz[i][j][k];
						cfNS18=cxNS[18]*Fx[i][j][k]+cyNS[18]*Fy[i][j][k]+czNS[18]*Fz[i][j][k];
						uxFxuyFyuzFz=ux[i][j][k]*Fx[i][j][k]+uy[i][j][k]*Fy[i][j][k]+uz[i][j][k]*Fz[i][j][k];
						Out0[i][j][k]=fIn0[i][j][k]-A*(fIn0[i][j][k]-Eq0[i][j][k]*rho[i][j][k])+Force0;
						Out1[i][j][k]=fIn1[i][j][k]-A*(fIn1[i][j][k]-Eq1[i][j][k]*rho[i][j][k])+Force1;
						Out2[i][j][k]=fIn2[i][j][k]-A*(fIn2[i][j][k]-Eq2[i][j][k]*rho[i][j][k])+Force2;
						Out3[i][j][k]=fIn3[i][j][k]-A*(fIn3[i][j][k]-Eq3[i][j][k]*rho[i][j][k])+Force3;
						Out4[i][j][k]=fIn4[i][j][k]-A*(fIn4[i][j][k]-Eq4[i][j][k]*rho[i][j][k])+Force4;
						Out5[i][j][k]=fIn5[i][j][k]-A*(fIn5[i][j][k]-Eq5[i][j][k]*rho[i][j][k])+Force5;
						Out6[i][j][k]=fIn6[i][j][k]-A*(fIn6[i][j][k]-Eq6[i][j][k]*rho[i][j][k])+Force6;
						Out7[i][j][k]=fIn7[i][j][k]-A*(fIn7[i][j][k]-Eq7[i][j][k]*rho[i][j][k])+Force7;
						Out8[i][j][k]=fIn8[i][j][k]-A*(fIn8[i][j][k]-Eq8[i][j][k]*rho[i][j][k])+Force8;
						Out9[i][j][k]=fIn9[i][j][k]-A*(fIn9[i][j][k]-Eq9[i][j][k]*rho[i][j][k])+Force9;
						Out10[i][j][k]=fIn10[i][j][k]-A*(fIn10[i][j][k]-Eq10[i][j][k]*rho[i][j][k])+Force10;
						Out11[i][j][k]=fIn11[i][j][k]-A*(fIn11[i][j][k]-Eq11[i][j][k]*rho[i][j][k])+Force11;
						Out12[i][j][k]=fIn12[i][j][k]-A*(fIn12[i][j][k]-Eq12[i][j][k]*rho[i][j][k])+Force12;
						Out13[i][j][k]=fIn13[i][j][k]-A*(fIn13[i][j][k]-Eq13[i][j][k]*rho[i][j][k])+Force13;
						Out14[i][j][k]=fIn14[i][j][k]-A*(fIn14[i][j][k]-Eq14[i][j][k]*rho[i][j][k])+Force14;
						Out15[i][j][k]=fIn15[i][j][k]-A*(fIn15[i][j][k]-Eq15[i][j][k]*rho[i][j][k])+Force15;
						Out16[i][j][k]=fIn16[i][j][k]-A*(fIn16[i][j][k]-Eq16[i][j][k]*rho[i][j][k])+Force16;
						Out17[i][j][k]=fIn17[i][j][k]-A*(fIn17[i][j][k]-Eq17[i][j][k]*rho[i][j][k])+Force17;
						Out18[i][j][k]=fIn18[i][j][k]-A*(fIn18[i][j][k]-Eq18[i][j][k]*rho[i][j][k])+Force18;
					}
					else{
						Out0[i][j][k]=0.0;
						Out1[i][j][k]=0.0;
						Out2[i][j][k]=0.0;
						Out3[i][j][k]=0.0;
						Out4[i][j][k]=0.0;
						Out5[i][j][k]=0.0;
						Out6[i][j][k]=0.0;
						Out7[i][j][k]=0.0;
						Out8[i][j][k]=0.0;
						Out9[i][j][k]=0.0;
						Out10[i][j][k]=0.0;
						Out11[i][j][k]=0.0;
						Out12[i][j][k]=0.0;
						Out13[i][j][k]=0.0;
						Out14[i][j][k]=0.0;
						Out15[i][j][k]=0.0;
						Out16[i][j][k]=0.0;
						Out17[i][j][k]=0.0;
						Out18[i][j][k]=0.0;
					}
				}
			}
		}
	}
}

void CollisionMRTTwoPhase(int m, int n, int q, double sigma){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k,mfIn0,mfIn1,mfIn2,mfIn3,mfIn4,mfIn5,mfIn6,mfIn7,mfIn8,mfIn9,\
	mfIn10,mfIn11,mfIn12,mfIn13,mfIn14,mfIn15,mfIn16,mfIn17,mfIn18,mfEq0,mfEq1,mfEq2,mfEq3,mfEq4,mfEq5,mfEq6,mfEq7,mfEq8,mfEq9,mfEq10,mfEq11,mfEq12,mfEq13,mfEq14,\
	mfEq15,mfEq16,mfEq17,mfEq18,mForce0,mForce1,mForce2,mForce3,mForce4,mForce5,mForce6,mForce7,mForce8,mForce9,mForce10,mForce11,mForce12,mForce13,mForce14,mForce15,\
	mForce16,mForce17,mForce18)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0){
						mfIn0=(A_g[0][0]*(M[0][0]*fIn0[i][j][k]+M[0][1]*fIn1[i][j][k]+M[0][2]*fIn2[i][j][k]+M[0][3]*fIn3[i][j][k]+M[0][4]*fIn4[i][j][k]+M[0][5]*fIn5[i][j][k]+M[0][6]*fIn6[i][j][k]+M[0][7]*fIn7[i][j][k]+M[0][8]*fIn8[i][j][k]\
							+M[0][9]*fIn9[i][j][k]+M[0][10]*fIn10[i][j][k]+M[0][11]*fIn11[i][j][k]+M[0][12]*fIn12[i][j][k]+M[0][13]*fIn13[i][j][k]+M[0][14]*fIn14[i][j][k]+M[0][15]*fIn15[i][j][k]+M[0][16]*fIn16[i][j][k]+M[0][17]*fIn17[i][j][k]+M[0][18]*fIn18[i][j][k]));
						mfIn1=(A_g[1][1]*(M[1][0]*fIn0[i][j][k]+M[1][1]*fIn1[i][j][k]+M[1][2]*fIn2[i][j][k]+M[1][3]*fIn3[i][j][k]+M[1][4]*fIn4[i][j][k]+M[1][5]*fIn5[i][j][k]+M[1][6]*fIn6[i][j][k]+M[1][7]*fIn7[i][j][k]+M[1][8]*fIn8[i][j][k]\
							+M[1][9]*fIn9[i][j][k]+M[1][10]*fIn10[i][j][k]+M[1][11]*fIn11[i][j][k]+M[1][12]*fIn12[i][j][k]+M[1][13]*fIn13[i][j][k]+M[1][14]*fIn14[i][j][k]+M[1][15]*fIn15[i][j][k]+M[1][16]*fIn16[i][j][k]+M[1][17]*fIn17[i][j][k]+M[1][18]*fIn18[i][j][k]));
						mfIn2=(A_g[2][2]*(M[2][0]*fIn0[i][j][k]+M[2][1]*fIn1[i][j][k]+M[2][2]*fIn2[i][j][k]+M[2][3]*fIn3[i][j][k]+M[2][4]*fIn4[i][j][k]+M[2][5]*fIn5[i][j][k]+M[2][6]*fIn6[i][j][k]+M[2][7]*fIn7[i][j][k]+M[2][8]*fIn8[i][j][k]\
							+M[2][9]*fIn9[i][j][k]+M[2][10]*fIn10[i][j][k]+M[2][11]*fIn11[i][j][k]+M[2][12]*fIn12[i][j][k]+M[2][13]*fIn13[i][j][k]+M[2][14]*fIn14[i][j][k]+M[2][15]*fIn15[i][j][k]+M[2][16]*fIn16[i][j][k]+M[2][17]*fIn17[i][j][k]+M[2][18]*fIn18[i][j][k]));
						mfIn3=(A_g[3][3]*(M[3][0]*fIn0[i][j][k]+M[3][1]*fIn1[i][j][k]+M[3][2]*fIn2[i][j][k]+M[3][3]*fIn3[i][j][k]+M[3][4]*fIn4[i][j][k]+M[3][5]*fIn5[i][j][k]+M[3][6]*fIn6[i][j][k]+M[3][7]*fIn7[i][j][k]+M[3][8]*fIn8[i][j][k]\
							+M[3][9]*fIn9[i][j][k]+M[3][10]*fIn10[i][j][k]+M[3][11]*fIn11[i][j][k]+M[3][12]*fIn12[i][j][k]+M[3][13]*fIn13[i][j][k]+M[3][14]*fIn14[i][j][k]+M[3][15]*fIn15[i][j][k]+M[3][16]*fIn16[i][j][k]+M[3][17]*fIn17[i][j][k]+M[3][18]*fIn18[i][j][k]));
						mfIn4=(A_g[4][4]*(M[4][0]*fIn0[i][j][k]+M[4][1]*fIn1[i][j][k]+M[4][2]*fIn2[i][j][k]+M[4][3]*fIn3[i][j][k]+M[4][4]*fIn4[i][j][k]+M[4][5]*fIn5[i][j][k]+M[4][6]*fIn6[i][j][k]+M[4][7]*fIn7[i][j][k]+M[4][8]*fIn8[i][j][k]\
							+M[4][9]*fIn9[i][j][k]+M[4][10]*fIn10[i][j][k]+M[4][11]*fIn11[i][j][k]+M[4][12]*fIn12[i][j][k]+M[4][13]*fIn13[i][j][k]+M[4][14]*fIn14[i][j][k]+M[4][15]*fIn15[i][j][k]+M[4][16]*fIn16[i][j][k]+M[4][17]*fIn17[i][j][k]+M[4][18]*fIn18[i][j][k]));
						mfIn5=(A_g[5][5]*(M[5][0]*fIn0[i][j][k]+M[5][1]*fIn1[i][j][k]+M[5][2]*fIn2[i][j][k]+M[5][3]*fIn3[i][j][k]+M[5][4]*fIn4[i][j][k]+M[5][5]*fIn5[i][j][k]+M[5][6]*fIn6[i][j][k]+M[5][7]*fIn7[i][j][k]+M[5][8]*fIn8[i][j][k]\
							+M[5][9]*fIn9[i][j][k]+M[5][10]*fIn10[i][j][k]+M[5][11]*fIn11[i][j][k]+M[5][12]*fIn12[i][j][k]+M[5][13]*fIn13[i][j][k]+M[5][14]*fIn14[i][j][k]+M[5][15]*fIn15[i][j][k]+M[5][16]*fIn16[i][j][k]+M[5][17]*fIn17[i][j][k]+M[5][18]*fIn18[i][j][k]));
						mfIn6=(A_g[6][6]*(M[6][0]*fIn0[i][j][k]+M[6][1]*fIn1[i][j][k]+M[6][2]*fIn2[i][j][k]+M[6][3]*fIn3[i][j][k]+M[6][4]*fIn4[i][j][k]+M[6][5]*fIn5[i][j][k]+M[6][6]*fIn6[i][j][k]+M[6][7]*fIn7[i][j][k]+M[6][8]*fIn8[i][j][k]\
							+M[6][9]*fIn9[i][j][k]+M[6][10]*fIn10[i][j][k]+M[6][11]*fIn11[i][j][k]+M[6][12]*fIn12[i][j][k]+M[6][13]*fIn13[i][j][k]+M[6][14]*fIn14[i][j][k]+M[6][15]*fIn15[i][j][k]+M[6][16]*fIn16[i][j][k]+M[6][17]*fIn17[i][j][k]+M[6][18]*fIn18[i][j][k]));
						mfIn7=(A_g[7][7]*(M[7][0]*fIn0[i][j][k]+M[7][1]*fIn1[i][j][k]+M[7][2]*fIn2[i][j][k]+M[7][3]*fIn3[i][j][k]+M[7][4]*fIn4[i][j][k]+M[7][5]*fIn5[i][j][k]+M[7][6]*fIn6[i][j][k]+M[7][7]*fIn7[i][j][k]+M[7][8]*fIn8[i][j][k]\
							+M[7][9]*fIn9[i][j][k]+M[7][10]*fIn10[i][j][k]+M[7][11]*fIn11[i][j][k]+M[7][12]*fIn12[i][j][k]+M[7][13]*fIn13[i][j][k]+M[7][14]*fIn14[i][j][k]+M[7][15]*fIn15[i][j][k]+M[7][16]*fIn16[i][j][k]+M[7][17]*fIn17[i][j][k]+M[7][18]*fIn18[i][j][k]));
						mfIn8=(A_g[8][8]*(M[8][0]*fIn0[i][j][k]+M[8][1]*fIn1[i][j][k]+M[8][2]*fIn2[i][j][k]+M[8][3]*fIn3[i][j][k]+M[8][4]*fIn4[i][j][k]+M[8][5]*fIn5[i][j][k]+M[8][6]*fIn6[i][j][k]+M[8][7]*fIn7[i][j][k]+M[8][8]*fIn8[i][j][k]\
							+M[8][9]*fIn9[i][j][k]+M[8][10]*fIn10[i][j][k]+M[8][11]*fIn11[i][j][k]+M[8][12]*fIn12[i][j][k]+M[8][13]*fIn13[i][j][k]+M[8][14]*fIn14[i][j][k]+M[8][15]*fIn15[i][j][k]+M[8][16]*fIn16[i][j][k]+M[8][17]*fIn17[i][j][k]+M[8][18]*fIn18[i][j][k]));
						mfIn9=(A_g[9][9]*(M[9][0]*fIn0[i][j][k]+M[9][1]*fIn1[i][j][k]+M[9][2]*fIn2[i][j][k]+M[9][3]*fIn3[i][j][k]+M[9][4]*fIn4[i][j][k]+M[9][5]*fIn5[i][j][k]+M[9][6]*fIn6[i][j][k]+M[9][7]*fIn7[i][j][k]+M[9][8]*fIn8[i][j][k]\
							+M[9][9]*fIn9[i][j][k]+M[9][10]*fIn10[i][j][k]+M[9][11]*fIn11[i][j][k]+M[9][12]*fIn12[i][j][k]+M[9][13]*fIn13[i][j][k]+M[9][14]*fIn14[i][j][k]+M[9][15]*fIn15[i][j][k]+M[9][16]*fIn16[i][j][k]+M[9][17]*fIn17[i][j][k]+M[9][18]*fIn18[i][j][k]));
						mfIn10=(A_g[10][10]*(M[10][0]*fIn0[i][j][k]+M[10][1]*fIn1[i][j][k]+M[10][2]*fIn2[i][j][k]+M[10][3]*fIn3[i][j][k]+M[10][4]*fIn4[i][j][k]+M[10][5]*fIn5[i][j][k]+M[10][6]*fIn6[i][j][k]+M[10][7]*fIn7[i][j][k]+M[10][8]*fIn8[i][j][k]\
							+M[10][9]*fIn9[i][j][k]+M[10][10]*fIn10[i][j][k]+M[10][11]*fIn11[i][j][k]+M[10][12]*fIn12[i][j][k]+M[10][13]*fIn13[i][j][k]+M[10][14]*fIn14[i][j][k]+M[10][15]*fIn15[i][j][k]+M[10][16]*fIn16[i][j][k]+M[10][17]*fIn17[i][j][k]+M[10][18]*fIn18[i][j][k]));
						mfIn11=(A_g[11][11]*(M[11][0]*fIn0[i][j][k]+M[11][1]*fIn1[i][j][k]+M[11][2]*fIn2[i][j][k]+M[11][3]*fIn3[i][j][k]+M[11][4]*fIn4[i][j][k]+M[11][5]*fIn5[i][j][k]+M[11][6]*fIn6[i][j][k]+M[11][7]*fIn7[i][j][k]+M[11][8]*fIn8[i][j][k]\
							+M[11][9]*fIn9[i][j][k]+M[11][10]*fIn10[i][j][k]+M[11][11]*fIn11[i][j][k]+M[11][12]*fIn12[i][j][k]+M[11][13]*fIn13[i][j][k]+M[11][14]*fIn14[i][j][k]+M[11][15]*fIn15[i][j][k]+M[11][16]*fIn16[i][j][k]+M[11][17]*fIn17[i][j][k]+M[11][18]*fIn18[i][j][k]));
						mfIn12=(A_g[12][12]*(M[12][0]*fIn0[i][j][k]+M[12][1]*fIn1[i][j][k]+M[12][2]*fIn2[i][j][k]+M[12][3]*fIn3[i][j][k]+M[12][4]*fIn4[i][j][k]+M[12][5]*fIn5[i][j][k]+M[12][6]*fIn6[i][j][k]+M[12][7]*fIn7[i][j][k]+M[12][8]*fIn8[i][j][k]\
							+M[12][9]*fIn9[i][j][k]+M[12][10]*fIn10[i][j][k]+M[12][11]*fIn11[i][j][k]+M[12][12]*fIn12[i][j][k]+M[12][13]*fIn13[i][j][k]+M[12][14]*fIn14[i][j][k]+M[12][15]*fIn15[i][j][k]+M[12][16]*fIn16[i][j][k]+M[12][17]*fIn17[i][j][k]+M[12][18]*fIn18[i][j][k]));
						mfIn13=(A_g[13][13]*(M[13][0]*fIn0[i][j][k]+M[13][1]*fIn1[i][j][k]+M[13][2]*fIn2[i][j][k]+M[13][3]*fIn3[i][j][k]+M[13][4]*fIn4[i][j][k]+M[13][5]*fIn5[i][j][k]+M[13][6]*fIn6[i][j][k]+M[13][7]*fIn7[i][j][k]+M[13][8]*fIn8[i][j][k]\
							+M[13][9]*fIn9[i][j][k]+M[13][10]*fIn10[i][j][k]+M[13][11]*fIn11[i][j][k]+M[13][12]*fIn12[i][j][k]+M[13][13]*fIn13[i][j][k]+M[13][14]*fIn14[i][j][k]+M[13][15]*fIn15[i][j][k]+M[13][16]*fIn16[i][j][k]+M[13][17]*fIn17[i][j][k]+M[13][18]*fIn18[i][j][k]));
						mfIn14=(A_g[14][14]*(M[14][0]*fIn0[i][j][k]+M[14][1]*fIn1[i][j][k]+M[14][2]*fIn2[i][j][k]+M[14][3]*fIn3[i][j][k]+M[14][4]*fIn4[i][j][k]+M[14][5]*fIn5[i][j][k]+M[14][6]*fIn6[i][j][k]+M[14][7]*fIn7[i][j][k]+M[14][8]*fIn8[i][j][k]\
							+M[14][9]*fIn9[i][j][k]+M[14][10]*fIn10[i][j][k]+M[14][11]*fIn11[i][j][k]+M[14][12]*fIn12[i][j][k]+M[14][13]*fIn13[i][j][k]+M[14][14]*fIn14[i][j][k]+M[14][15]*fIn15[i][j][k]+M[14][16]*fIn16[i][j][k]+M[14][17]*fIn17[i][j][k]+M[14][18]*fIn18[i][j][k]));
						mfIn15=(A_g[15][15]*(M[15][0]*fIn0[i][j][k]+M[15][1]*fIn1[i][j][k]+M[15][2]*fIn2[i][j][k]+M[15][3]*fIn3[i][j][k]+M[15][4]*fIn4[i][j][k]+M[15][5]*fIn5[i][j][k]+M[15][6]*fIn6[i][j][k]+M[15][7]*fIn7[i][j][k]+M[15][8]*fIn8[i][j][k]\
							+M[15][9]*fIn9[i][j][k]+M[15][10]*fIn10[i][j][k]+M[15][11]*fIn11[i][j][k]+M[15][12]*fIn12[i][j][k]+M[15][13]*fIn13[i][j][k]+M[15][14]*fIn14[i][j][k]+M[15][15]*fIn15[i][j][k]+M[15][16]*fIn16[i][j][k]+M[15][17]*fIn17[i][j][k]+M[15][18]*fIn18[i][j][k]));
						mfIn16=(A_g[16][16]*(M[16][0]*fIn0[i][j][k]+M[16][1]*fIn1[i][j][k]+M[16][2]*fIn2[i][j][k]+M[16][3]*fIn3[i][j][k]+M[16][4]*fIn4[i][j][k]+M[16][5]*fIn5[i][j][k]+M[16][6]*fIn6[i][j][k]+M[16][7]*fIn7[i][j][k]+M[16][8]*fIn8[i][j][k]\
							+M[16][9]*fIn9[i][j][k]+M[16][10]*fIn10[i][j][k]+M[16][11]*fIn11[i][j][k]+M[16][12]*fIn12[i][j][k]+M[16][13]*fIn13[i][j][k]+M[16][14]*fIn14[i][j][k]+M[16][15]*fIn15[i][j][k]+M[16][16]*fIn16[i][j][k]+M[16][17]*fIn17[i][j][k]+M[16][18]*fIn18[i][j][k]));
						mfIn17=(A_g[17][17]*(M[17][0]*fIn0[i][j][k]+M[17][1]*fIn1[i][j][k]+M[17][2]*fIn2[i][j][k]+M[17][3]*fIn3[i][j][k]+M[17][4]*fIn4[i][j][k]+M[17][5]*fIn5[i][j][k]+M[17][6]*fIn6[i][j][k]+M[17][7]*fIn7[i][j][k]+M[17][8]*fIn8[i][j][k]\
							+M[17][9]*fIn9[i][j][k]+M[17][10]*fIn10[i][j][k]+M[17][11]*fIn11[i][j][k]+M[17][12]*fIn12[i][j][k]+M[17][13]*fIn13[i][j][k]+M[17][14]*fIn14[i][j][k]+M[17][15]*fIn15[i][j][k]+M[17][16]*fIn16[i][j][k]+M[17][17]*fIn17[i][j][k]+M[17][18]*fIn18[i][j][k]));
						mfIn18=(A_g[18][18]*(M[18][0]*fIn0[i][j][k]+M[18][1]*fIn1[i][j][k]+M[18][2]*fIn2[i][j][k]+M[18][3]*fIn3[i][j][k]+M[18][4]*fIn4[i][j][k]+M[18][5]*fIn5[i][j][k]+M[18][6]*fIn6[i][j][k]+M[18][7]*fIn7[i][j][k]+M[18][8]*fIn8[i][j][k]\
							+M[18][9]*fIn9[i][j][k]+M[18][10]*fIn10[i][j][k]+M[18][11]*fIn11[i][j][k]+M[18][12]*fIn12[i][j][k]+M[18][13]*fIn13[i][j][k]+M[18][14]*fIn14[i][j][k]+M[18][15]*fIn15[i][j][k]+M[18][16]*fIn16[i][j][k]+M[18][17]*fIn17[i][j][k]+M[18][18]*fIn18[i][j][k]));
						mfEq0=A_g[0][0]*rho[i][j][k];
						mfEq1=A_g[1][1]*(-11.0*rho[i][j][k]+19.0*((rho[i][j][k]*ux[i][j][k])*(rho[i][j][k]*ux[i][j][k])+(rho[i][j][k]*uy[i][j][k])*(rho[i][j][k]*uy[i][j][k])+(rho[i][j][k]*uz[i][j][k])*(rho[i][j][k]*uz[i][j][k]))/rho[i][j][k]);
						mfEq2=A_g[2][2]*(3.0*rho[i][j][k]-11.0/2.0*((rho[i][j][k]*ux[i][j][k])*(rho[i][j][k]*ux[i][j][k])+(rho[i][j][k]*uy[i][j][k])*(rho[i][j][k]*uy[i][j][k])+(rho[i][j][k]*uz[i][j][k])*(rho[i][j][k]*uz[i][j][k]))/rho[i][j][k]);
						mfEq3=A_g[3][3]*(rho[i][j][k]*ux[i][j][k]);
						mfEq4=A_g[4][4]*(-2.0/3.0*rho[i][j][k]*ux[i][j][k]);
						mfEq5=A_g[5][5]*(rho[i][j][k]*uy[i][j][k]);
						mfEq6=A_g[6][6]*(-2.0/3.0*rho[i][j][k]*uy[i][j][k]);
						mfEq7=A_g[7][7]*(rho[i][j][k]*uz[i][j][k]);
						mfEq8=A_g[8][8]*(-2.0/3.0*rho[i][j][k]*uz[i][j][k]);
						mfEq9=A_g[9][9]*((2.0*(rho[i][j][k]*ux[i][j][k])*(rho[i][j][k]*ux[i][j][k])-((rho[i][j][k]*uy[i][j][k])*(rho[i][j][k]*uy[i][j][k])+(rho[i][j][k]*uz[i][j][k])*(rho[i][j][k]*uz[i][j][k])))/rho[i][j][k]);
						mfEq10=A_g[10][10]*(-1.0/2.0*((2.0*(rho[i][j][k]*ux[i][j][k])*(rho[i][j][k]*ux[i][j][k])-((rho[i][j][k]*uy[i][j][k])*(rho[i][j][k]*uy[i][j][k])+(rho[i][j][k]*uz[i][j][k])*(rho[i][j][k]*uz[i][j][k])))/rho[i][j][k]));
						mfEq11=A_g[11][11]*(((rho[i][j][k]*uy[i][j][k])*(rho[i][j][k]*uy[i][j][k])-(rho[i][j][k]*uz[i][j][k])*(rho[i][j][k]*uz[i][j][k]))/rho[i][j][k]);
						mfEq12=A_g[12][12]*(-1.0/2.0*((rho[i][j][k]*uy[i][j][k])*(rho[i][j][k]*uy[i][j][k])-(rho[i][j][k]*uz[i][j][k])*(rho[i][j][k]*uz[i][j][k]))/rho[i][j][k]);
						mfEq13=A_g[13][13]*(rho[i][j][k]*ux[i][j][k]*rho[i][j][k]*uy[i][j][k]/rho[i][j][k]);
						mfEq14=A_g[14][14]*(rho[i][j][k]*uy[i][j][k]*rho[i][j][k]*uz[i][j][k]/rho[i][j][k]);
						mfEq15=A_g[15][15]*(rho[i][j][k]*ux[i][j][k]*rho[i][j][k]*uz[i][j][k]/rho[i][j][k]);
						mfEq16=0.0;
						mfEq17=0.0;
						mfEq18=0.0;
						mForce0=(Factor0*0.0);
						mForce1=Factor1*((38.0*(ux[i][j][k]*Fx[i][j][k]+uy[i][j][k]*Fy[i][j][k]+uz[i][j][k]*Fz[i][j][k]))+114*sigma*(Fx[i][j][k]*Fx[i][j][k]+Fy[i][j][k]*Fy[i][j][k]+Fz[i][j][k]*Fz[i][j][k])/(sqrt(2*fabs(pressure[i][j][k]-1.0/3*rho[i][j][k])/(6.0*1/18))*sqrt(2*fabs(pressure[i][j][k]-1.0/3*rho[i][j][k])/(6.0*1/18)))/(1.0/A_g[1][1]-0.5));
						//mForce1=Factor1*(38.0*(ux[i][j][k]*Fx[i][j][k]+uy[i][j][k]*Fy[i][j][k]+uz[i][j][k]*Fz[i][j][k]));
						mForce2=Factor2*(-11.0*(ux[i][j][k]*Fx[i][j][k]+uy[i][j][k]*Fy[i][j][k]+uz[i][j][k]*Fz[i][j][k]));
						mForce3=Factor3*Fx[i][j][k];
						mForce4=Factor4*(-2.0/3.0*Fx[i][j][k]);
						mForce5=Factor5*Fy[i][j][k];
						mForce6=Factor6*(-2.0/3.0*Fy[i][j][k]);
						mForce7=Factor7*Fz[i][j][k];
						mForce8=Factor8*(-2.0/3.0*Fz[i][j][k]);
						mForce9=Factor9*(2.0*(2.0*ux[i][j][k]*Fx[i][j][k]-uy[i][j][k]*Fy[i][j][k]-uz[i][j][k]*Fz[i][j][k]));
						mForce10=Factor10*(-2.0*ux[i][j][k]*Fx[i][j][k]+uy[i][j][k]*Fy[i][j][k]+uz[i][j][k]*Fz[i][j][k]);
						mForce11=Factor11*(2.0*(uy[i][j][k]*Fy[i][j][k]-uz[i][j][k]*Fz[i][j][k]));
						mForce12=Factor12*(-uy[i][j][k]*Fy[i][j][k]+uz[i][j][k]*Fz[i][j][k]);
						mForce13=Factor13*(uy[i][j][k]*Fx[i][j][k]+ux[i][j][k]*Fy[i][j][k]);
						mForce14=Factor14*(uz[i][j][k]*Fy[i][j][k]+uy[i][j][k]*Fz[i][j][k]);
						mForce15=Factor15*(uz[i][j][k]*Fx[i][j][k]+ux[i][j][k]*Fz[i][j][k]);
						mForce16=Factor16*0.0;
						mForce17=Factor17*0.0;
						mForce18=Factor18*0.0;
						Out0[i][j][k]=fIn0[i][j][k]-(vfIn0-vfEq0)+vForce0;
						Out1[i][j][k]=fIn1[i][j][k]-(vfIn1-vfEq1)+vForce1;
						Out2[i][j][k]=fIn2[i][j][k]-(vfIn2-vfEq2)+vForce2;
						Out3[i][j][k]=fIn3[i][j][k]-(vfIn3-vfEq3)+vForce3;
						Out4[i][j][k]=fIn4[i][j][k]-(vfIn4-vfEq4)+vForce4;
						Out5[i][j][k]=fIn5[i][j][k]-(vfIn5-vfEq5)+vForce5;
						Out6[i][j][k]=fIn6[i][j][k]-(vfIn6-vfEq6)+vForce6;
						Out7[i][j][k]=fIn7[i][j][k]-(vfIn7-vfEq7)+vForce7;
						Out8[i][j][k]=fIn8[i][j][k]-(vfIn8-vfEq8)+vForce8;
						Out9[i][j][k]=fIn9[i][j][k]-(vfIn9-vfEq9)+vForce9;
						Out10[i][j][k]=fIn10[i][j][k]-(vfIn10-vfEq10)+vForce10;
						Out11[i][j][k]=fIn11[i][j][k]-(vfIn11-vfEq11)+vForce11;
						Out12[i][j][k]=fIn12[i][j][k]-(vfIn12-vfEq12)+vForce12;
						Out13[i][j][k]=fIn13[i][j][k]-(vfIn13-vfEq13)+vForce13;
						Out14[i][j][k]=fIn14[i][j][k]-(vfIn14-vfEq14)+vForce14;
						Out15[i][j][k]=fIn15[i][j][k]-(vfIn15-vfEq15)+vForce15;
						Out16[i][j][k]=fIn16[i][j][k]-(vfIn16-vfEq16)+vForce16;
						Out17[i][j][k]=fIn17[i][j][k]-(vfIn17-vfEq17)+vForce17;
						Out18[i][j][k]=fIn18[i][j][k]-(vfIn18-vfEq18)+vForce18;
					}
					else{
						Out0[i][j][k]=0.0;
						Out1[i][j][k]=0.0;
						Out2[i][j][k]=0.0;
						Out3[i][j][k]=0.0;
						Out4[i][j][k]=0.0;
						Out5[i][j][k]=0.0;
						Out6[i][j][k]=0.0;
						Out7[i][j][k]=0.0;
						Out8[i][j][k]=0.0;
						Out9[i][j][k]=0.0;
						Out10[i][j][k]=0.0;
						Out11[i][j][k]=0.0;
						Out12[i][j][k]=0.0;
						Out13[i][j][k]=0.0;
						Out14[i][j][k]=0.0;
						Out15[i][j][k]=0.0;
						Out16[i][j][k]=0.0;
						Out17[i][j][k]=0.0;
						Out18[i][j][k]=0.0;
					}
				}
			}
		}
	}
}

void FieldCalculationTwoPhase(int m, int n, int q, double SOLID, double b, double R, double T, double a){
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
						rho[i][j][k]=fIn0[i][j][k]+fIn1[i][j][k]+fIn2[i][j][k]+fIn3[i][j][k]+fIn4[i][j][k]+fIn5[i][j][k]+fIn6[i][j][k]+fIn7[i][j][k]+fIn8[i][j][k]\
							+fIn9[i][j][k]+fIn10[i][j][k]+fIn11[i][j][k]+fIn12[i][j][k]+fIn13[i][j][k]+fIn14[i][j][k]+fIn15[i][j][k]+fIn16[i][j][k]+fIn17[i][j][k]+fIn18[i][j][k];
						bRho4=b*rho[i][j][k]/4.0;
						BRho4=1.0-bRho4;
						pressure[i][j][k]=Pressure;
						ux[i][j][k]=uTotX+Fx[i][j][k]/2.0/rho[i][j][k];
						uy[i][j][k]=uTotY+Fy[i][j][k]/2.0/rho[i][j][k];
						uz[i][j][k]=uTotZ+Fz[i][j][k]/2.0/rho[i][j][k];
					}
					if (domain[i][j][k]==1){
						rho[i][j][k]=SOLID;
						pressure[i][j][k]=0.0;
						ux[i][j][k]=0.0;
						uy[i][j][k]=0.0;
						uz[i][j][k]=0.0;
						fIn0[i][j][k]=0.0;
						fIn1[i][j][k]=0.0;
						fIn2[i][j][k]=0.0;
						fIn3[i][j][k]=0.0;
						fIn4[i][j][k]=0.0;
						fIn5[i][j][k]=0.0;
						fIn6[i][j][k]=0.0;
						fIn7[i][j][k]=0.0;
						fIn8[i][j][k]=0.0;
						fIn9[i][j][k]=0.0;
						fIn10[i][j][k]=0.0;
						fIn11[i][j][k]=0.0;
						fIn12[i][j][k]=0.0;
						fIn13[i][j][k]=0.0;
						fIn14[i][j][k]=0.0;
						fIn15[i][j][k]=0.0;
						fIn16[i][j][k]=0.0;
						fIn17[i][j][k]=0.0;
						fIn18[i][j][k]=0.0;
					}
				}
			}
		}
	}
}

void rhoShift(int m, int n, int q){
	shift1(rho,rho1,m,n,q);
	shift2(rho,rho2,m,n,q);
	shift3(rho,rho3,m,n,q);
	shift4(rho,rho4,m,n,q);
	shift5(rho,rho5,m,n,q);
	shift6(rho,rho6,m,n,q);
	shift7(rho,rho7,m,n,q);
	shift8(rho,rho8,m,n,q);
	shift9(rho,rho9,m,n,q);
	shift10(rho,rho10,m,n,q);
	shift11(rho,rho11,m,n,q);
	shift12(rho,rho12,m,n,q);
	shift13(rho,rho13,m,n,q);
	shift14(rho,rho14,m,n,q);
	shift15(rho,rho15,m,n,q);
	shift16(rho,rho16,m,n,q);
	shift17(rho,rho17,m,n,q);
	shift18(rho,rho18,m,n,q);
//	shift19(rho, rho19, m, n, q);
//	shift20(rho, rho20, m, n, q);
//	shift21(rho, rho21, m, n, q);
//	shift22(rho, rho22, m, n, q);
//	shift23(rho, rho23, m, n, q);
//	shift24(rho, rho24, m, n, q);
//	shift25(rho, rho25, m, n, q);
//	shift26(rho, rho26, m, n, q);
}

void pressureShift(int m, int n, int q){
	shift1(pressure,pressure1,m,n,q);
	shift2(pressure,pressure2,m,n,q);
	shift3(pressure,pressure3,m,n,q);
	shift4(pressure,pressure4,m,n,q);
	shift5(pressure,pressure5,m,n,q);
	shift6(pressure,pressure6,m,n,q);
	shift7(pressure,pressure7,m,n,q);
	shift8(pressure,pressure8,m,n,q);
	shift9(pressure,pressure9,m,n,q);
	shift10(pressure,pressure10,m,n,q);
	shift11(pressure,pressure11,m,n,q);
	shift12(pressure,pressure12,m,n,q);
	shift13(pressure,pressure13,m,n,q);
	shift14(pressure,pressure14,m,n,q);
	shift15(pressure,pressure15,m,n,q);
	shift16(pressure,pressure16,m,n,q);
	shift17(pressure,pressure17,m,n,q);
	shift18(pressure,pressure18,m,n,q);
}

void StreamTwoPhase(int m, int n, int q){
	shift0(Out0,fIn0,m,n,q);
	shift1(Out1,fIn1,m,n,q);
	shift2(Out2,fIn2,m,n,q);
	shift3(Out3,fIn3,m,n,q);
	shift4(Out4,fIn4,m,n,q);
	shift5(Out5,fIn5,m,n,q);
	shift6(Out6,fIn6,m,n,q);
	shift7(Out7,fIn7,m,n,q);
	shift8(Out8,fIn8,m,n,q);
	shift9(Out9,fIn9,m,n,q);
	shift10(Out10,fIn10,m,n,q);
	shift11(Out11,fIn11,m,n,q);
	shift12(Out12,fIn12,m,n,q);
	shift13(Out13,fIn13,m,n,q);
	shift14(Out14,fIn14,m,n,q);
	shift15(Out15,fIn15,m,n,q);
	shift16(Out16,fIn16,m,n,q);
	shift17(Out17,fIn17,m,n,q);
	shift18(Out18,fIn18,m,n,q);
}

void MemoryFreeTwoPhase(int m, int n, int q){
	freememory(rho,m,n,q);
	freememory(rhoLastStep,m,n,q);
	freememory(rho1,m,n,q);
	freememory(rho2,m,n,q);
	freememory(rho3,m,n,q);
	freememory(rho4,m,n,q);
	freememory(rho5,m,n,q);
	freememory(rho6,m,n,q);
	freememory(rho7,m,n,q);
	freememory(rho8,m,n,q);
	freememory(rho9,m,n,q);
	freememory(rho10,m,n,q);
	freememory(rho11,m,n,q);
	freememory(rho12,m,n,q);
	freememory(rho13,m,n,q);
	freememory(rho14,m,n,q);
	freememory(rho15,m,n,q);
	freememory(rho16,m,n,q);
	freememory(rho17,m,n,q);
	freememory(rho18,m,n,q);
//	freememory(rho19, m, n, q);
//	freememory(rho20, m, n, q);
//	freememory(rho21, m, n, q);
//	freememory(rho22, m, n, q);
//	freememory(rho23, m, n, q);
//	freememory(rho24, m, n, q);
//	freememory(rho25, m, n, q);
//	freememory(rho26, m, n, q);
	freememory(fIn0,m,n,q);
	freememory(fIn1,m,n,q);
	freememory(fIn2,m,n,q);
	freememory(fIn3,m,n,q);
	freememory(fIn4,m,n,q);
	freememory(fIn5,m,n,q);
	freememory(fIn6,m,n,q);
	freememory(fIn7,m,n,q);
	freememory(fIn8,m,n,q);
	freememory(fIn9,m,n,q);
	freememory(fIn10,m,n,q);
	freememory(fIn11,m,n,q);
	freememory(fIn12,m,n,q);
	freememory(fIn13,m,n,q);
	freememory(fIn14,m,n,q);
	freememory(fIn15,m,n,q);
	freememory(fIn16,m,n,q);
	freememory(fIn17,m,n,q);
	freememory(fIn18,m,n,q);
	freememory(ux,m,n,q);
	freememory(uy,m,n,q);
	freememory(uz,m,n,q);
	freememory(Fx,m,n,q);
	freememory(Fy,m,n,q);
	freememory(Fz,m,n,q);
	freememory(bodyForce,m,n,q);
	freememory(pressure,m,n,q);
	freememory(pressure1,m,n,q);
	freememory(pressure2,m,n,q);
	freememory(pressure3,m,n,q);
	freememory(pressure4,m,n,q);
	freememory(pressure5,m,n,q);
	freememory(pressure6,m,n,q);
	freememory(pressure7,m,n,q);
	freememory(pressure8,m,n,q);
	freememory(pressure9,m,n,q);
	freememory(pressure10,m,n,q);
	freememory(pressure11,m,n,q);
	freememory(pressure12,m,n,q);
	freememory(pressure13,m,n,q);
	freememory(pressure14,m,n,q);
	freememory(pressure15,m,n,q);
	freememory(pressure16,m,n,q);
	freememory(pressure17,m,n,q);
	freememory(pressure18,m,n,q);
}

void rhoFieldStoreTwoPhase(int m, int n, int q){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					rhoLastStep[i][j][k]=rho[i][j][k];
				}
			}
		}
	}
}