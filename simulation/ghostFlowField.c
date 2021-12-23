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
#include "equilibrium.h"

#define viscosityGhost 0.5/3
#define jxGhost (cxNS[0]*fIn0Ghost[i][j][k]+cxNS[1]*fIn1Ghost[i][j][k]+cxNS[2]*fIn2Ghost[i][j][k]+cxNS[3]*fIn3Ghost[i][j][k]+cxNS[4]*fIn4Ghost[i][j][k]+cxNS[5]*fIn5Ghost[i][j][k]+cxNS[6]*fIn6Ghost[i][j][k]+cxNS[7]*fIn7Ghost[i][j][k]+cxNS[8]*fIn8Ghost[i][j][k]\
	+cxNS[9]*fIn9Ghost[i][j][k]+cxNS[10]*fIn10Ghost[i][j][k]+cxNS[11]*fIn11Ghost[i][j][k]+cxNS[12]*fIn12Ghost[i][j][k]+cxNS[13]*fIn13Ghost[i][j][k]+cxNS[14]*fIn14Ghost[i][j][k]+cxNS[15]*fIn15Ghost[i][j][k]+cxNS[16]*fIn16Ghost[i][j][k]+cxNS[17]*fIn17Ghost[i][j][k]+cxNS[18]*fIn18Ghost[i][j][k])
#define jyGhost (cyNS[0]*fIn0Ghost[i][j][k]+cyNS[1]*fIn1Ghost[i][j][k]+cyNS[2]*fIn2Ghost[i][j][k]+cyNS[3]*fIn3Ghost[i][j][k]+cyNS[4]*fIn4Ghost[i][j][k]+cyNS[5]*fIn5Ghost[i][j][k]+cyNS[6]*fIn6Ghost[i][j][k]+cyNS[7]*fIn7Ghost[i][j][k]+cyNS[8]*fIn8Ghost[i][j][k]\
	+cyNS[9]*fIn9Ghost[i][j][k]+cyNS[10]*fIn10Ghost[i][j][k]+cyNS[11]*fIn11Ghost[i][j][k]+cyNS[12]*fIn12Ghost[i][j][k]+cyNS[13]*fIn13Ghost[i][j][k]+cyNS[14]*fIn14Ghost[i][j][k]+cyNS[15]*fIn15Ghost[i][j][k]+cyNS[16]*fIn16Ghost[i][j][k]+cyNS[17]*fIn17Ghost[i][j][k]+cyNS[18]*fIn18Ghost[i][j][k])
#define jzGhost (czNS[0]*fIn0Ghost[i][j][k]+czNS[1]*fIn1Ghost[i][j][k]+czNS[2]*fIn2Ghost[i][j][k]+czNS[3]*fIn3Ghost[i][j][k]+czNS[4]*fIn4Ghost[i][j][k]+czNS[5]*fIn5Ghost[i][j][k]+czNS[6]*fIn6Ghost[i][j][k]+czNS[7]*fIn7Ghost[i][j][k]+czNS[8]*fIn8Ghost[i][j][k]\
	+czNS[9]*fIn9Ghost[i][j][k]+czNS[10]*fIn10Ghost[i][j][k]+czNS[11]*fIn11Ghost[i][j][k]+czNS[12]*fIn12Ghost[i][j][k]+czNS[13]*fIn13Ghost[i][j][k]+czNS[14]*fIn14Ghost[i][j][k]+czNS[15]*fIn15Ghost[i][j][k]+czNS[16]*fIn16Ghost[i][j][k]+czNS[17]*fIn17Ghost[i][j][k]+czNS[18]*fIn18Ghost[i][j][k])
#define uTotXGhost (jxGhost/rhoGhost[i][j][k])
#define uTotYGhost (jyGhost/rhoGhost[i][j][k])
#define uTotZGhost (jzGhost/rhoGhost[i][j][k])

#define Force0Ghost ((1-0.5*AGhost)*tNS[0]*(3.0*(cfNS0-uxFxuyFyuzFz)+3.0*cuNS0[i][j][k]*cfNS0))
#define Force1Ghost ((1-0.5*AGhost)*tNS[1]*(3.0*(cfNS1-uxFxuyFyuzFz)+3.0*cuNS1[i][j][k]*cfNS1))
#define Force2Ghost ((1-0.5*AGhost)*tNS[2]*(3.0*(cfNS2-uxFxuyFyuzFz)+3.0*cuNS2[i][j][k]*cfNS2))
#define Force3Ghost ((1-0.5*AGhost)*tNS[3]*(3.0*(cfNS3-uxFxuyFyuzFz)+3.0*cuNS3[i][j][k]*cfNS3))
#define Force4Ghost ((1-0.5*AGhost)*tNS[4]*(3.0*(cfNS4-uxFxuyFyuzFz)+3.0*cuNS4[i][j][k]*cfNS4))
#define Force5Ghost ((1-0.5*AGhost)*tNS[5]*(3.0*(cfNS5-uxFxuyFyuzFz)+3.0*cuNS5[i][j][k]*cfNS5))
#define Force6Ghost ((1-0.5*AGhost)*tNS[6]*(3.0*(cfNS6-uxFxuyFyuzFz)+3.0*cuNS6[i][j][k]*cfNS6))
#define Force7Ghost ((1-0.5*AGhost)*tNS[7]*(3.0*(cfNS7-uxFxuyFyuzFz)+3.0*cuNS7[i][j][k]*cfNS7))
#define Force8Ghost ((1-0.5*AGhost)*tNS[8]*(3.0*(cfNS8-uxFxuyFyuzFz)+3.0*cuNS8[i][j][k]*cfNS8))
#define Force9Ghost ((1-0.5*AGhost)*tNS[9]*(3.0*(cfNS9-uxFxuyFyuzFz)+3.0*cuNS9[i][j][k]*cfNS9))
#define Force10Ghost ((1-0.5*AGhost)*tNS[10]*(3.0*(cfNS10-uxFxuyFyuzFz)+3.0*cuNS10[i][j][k]*cfNS10))
#define Force11Ghost ((1-0.5*AGhost)*tNS[11]*(3.0*(cfNS11-uxFxuyFyuzFz)+3.0*cuNS11[i][j][k]*cfNS11))
#define Force12Ghost ((1-0.5*AGhost)*tNS[12]*(3.0*(cfNS12-uxFxuyFyuzFz)+3.0*cuNS12[i][j][k]*cfNS12))
#define Force13Ghost ((1-0.5*AGhost)*tNS[13]*(3.0*(cfNS13-uxFxuyFyuzFz)+3.0*cuNS13[i][j][k]*cfNS13))
#define Force14Ghost ((1-0.5*AGhost)*tNS[14]*(3.0*(cfNS14-uxFxuyFyuzFz)+3.0*cuNS14[i][j][k]*cfNS14))
#define Force15Ghost ((1-0.5*AGhost)*tNS[15]*(3.0*(cfNS15-uxFxuyFyuzFz)+3.0*cuNS15[i][j][k]*cfNS15))
#define Force16Ghost ((1-0.5*AGhost)*tNS[16]*(3.0*(cfNS16-uxFxuyFyuzFz)+3.0*cuNS16[i][j][k]*cfNS16))
#define Force17Ghost ((1-0.5*AGhost)*tNS[17]*(3.0*(cfNS17-uxFxuyFyuzFz)+3.0*cuNS17[i][j][k]*cfNS17))
#define Force18Ghost ((1-0.5*AGhost)*tNS[18]*(3.0*(cfNS18-uxFxuyFyuzFz)+3.0*cuNS18[i][j][k]*cfNS18))

#define vfGhostIn0 (MInverse[0][0]*mfIn0+MInverse[0][1]*mfIn1+MInverse[0][2]*mfIn2+MInverse[0][3]*mfIn3+MInverse[0][4]*mfIn4+MInverse[0][5]*mfIn5+MInverse[0][6]*mfIn6+MInverse[0][7]*mfIn7+MInverse[0][8]*mfIn8\
	+ MInverse[0][9] * mfIn9 + MInverse[0][10] * mfIn10 + MInverse[0][11] * mfIn11 + MInverse[0][12] * mfIn12 + MInverse[0][13] * mfIn13 + MInverse[0][14] * mfIn14 + MInverse[0][15] * mfIn15 + MInverse[0][16] * mfIn16 + MInverse[0][17] * mfIn17 + MInverse[0][18] * mfIn18)
#define vfGhostIn1 (MInverse[1][0]*mfIn0+MInverse[1][1]*mfIn1+MInverse[1][2]*mfIn2+MInverse[1][3]*mfIn3+MInverse[1][4]*mfIn4+MInverse[1][5]*mfIn5+MInverse[1][6]*mfIn6+MInverse[1][7]*mfIn7+MInverse[1][8]*mfIn8\
	+ MInverse[1][9] * mfIn9 + MInverse[1][10] * mfIn10 + MInverse[1][11] * mfIn11 + MInverse[1][12] * mfIn12 + MInverse[1][13] * mfIn13 + MInverse[1][14] * mfIn14 + MInverse[1][15] * mfIn15 + MInverse[1][16] * mfIn16 + MInverse[1][17] * mfIn17 + MInverse[1][18] * mfIn18)
#define vfGhostIn2 (MInverse[2][0]*mfIn0+MInverse[2][1]*mfIn1+MInverse[2][2]*mfIn2+MInverse[2][3]*mfIn3+MInverse[2][4]*mfIn4+MInverse[2][5]*mfIn5+MInverse[2][6]*mfIn6+MInverse[2][7]*mfIn7+MInverse[2][8]*mfIn8\
	+ MInverse[2][9] * mfIn9 + MInverse[2][10] * mfIn10 + MInverse[2][11] * mfIn11 + MInverse[2][12] * mfIn12 + MInverse[2][13] * mfIn13 + MInverse[2][14] * mfIn14 + MInverse[2][15] * mfIn15 + MInverse[2][16] * mfIn16 + MInverse[2][17] * mfIn17 + MInverse[2][18] * mfIn18)
#define vfGhostIn3 (MInverse[3][0]*mfIn0+MInverse[3][1]*mfIn1+MInverse[3][2]*mfIn2+MInverse[3][3]*mfIn3+MInverse[3][4]*mfIn4+MInverse[3][5]*mfIn5+MInverse[3][6]*mfIn6+MInverse[3][7]*mfIn7+MInverse[3][8]*mfIn8\
	+ MInverse[3][9] * mfIn9 + MInverse[3][10] * mfIn10 + MInverse[3][11] * mfIn11 + MInverse[3][12] * mfIn12 + MInverse[3][13] * mfIn13 + MInverse[3][14] * mfIn14 + MInverse[3][15] * mfIn15 + MInverse[3][16] * mfIn16 + MInverse[3][17] * mfIn17 + MInverse[3][18] * mfIn18)
#define vfGhostIn4 (MInverse[4][0]*mfIn0+MInverse[4][1]*mfIn1+MInverse[4][2]*mfIn2+MInverse[4][3]*mfIn3+MInverse[4][4]*mfIn4+MInverse[4][5]*mfIn5+MInverse[4][6]*mfIn6+MInverse[4][7]*mfIn7+MInverse[4][8]*mfIn8\
	+ MInverse[4][9] * mfIn9 + MInverse[4][10] * mfIn10 + MInverse[4][11] * mfIn11 + MInverse[4][12] * mfIn12 + MInverse[4][13] * mfIn13 + MInverse[4][14] * mfIn14 + MInverse[4][15] * mfIn15 + MInverse[4][16] * mfIn16 + MInverse[4][17] * mfIn17 + MInverse[4][18] * mfIn18)
#define vfGhostIn5 (MInverse[5][0]*mfIn0+MInverse[5][1]*mfIn1+MInverse[5][2]*mfIn2+MInverse[5][3]*mfIn3+MInverse[5][4]*mfIn4+MInverse[5][5]*mfIn5+MInverse[5][6]*mfIn6+MInverse[5][7]*mfIn7+MInverse[5][8]*mfIn8\
	+ MInverse[5][9] * mfIn9 + MInverse[5][10] * mfIn10 + MInverse[5][11] * mfIn11 + MInverse[5][12] * mfIn12 + MInverse[5][13] * mfIn13 + MInverse[5][14] * mfIn14 + MInverse[5][15] * mfIn15 + MInverse[5][16] * mfIn16 + MInverse[5][17] * mfIn17 + MInverse[5][18] * mfIn18)
#define vfGhostIn6 (MInverse[6][0]*mfIn0+MInverse[6][1]*mfIn1+MInverse[6][2]*mfIn2+MInverse[6][3]*mfIn3+MInverse[6][4]*mfIn4+MInverse[6][5]*mfIn5+MInverse[6][6]*mfIn6+MInverse[6][7]*mfIn7+MInverse[6][8]*mfIn8\
	+ MInverse[6][9] * mfIn9 + MInverse[6][10] * mfIn10 + MInverse[6][11] * mfIn11 + MInverse[6][12] * mfIn12 + MInverse[6][13] * mfIn13 + MInverse[6][14] * mfIn14 + MInverse[6][15] * mfIn15 + MInverse[6][16] * mfIn16 + MInverse[6][17] * mfIn17 + MInverse[6][18] * mfIn18)
#define vfGhostIn7 (MInverse[7][0]*mfIn0+MInverse[7][1]*mfIn1+MInverse[7][2]*mfIn2+MInverse[7][3]*mfIn3+MInverse[7][4]*mfIn4+MInverse[7][5]*mfIn5+MInverse[7][6]*mfIn6+MInverse[7][7]*mfIn7+MInverse[7][8]*mfIn8\
	+ MInverse[7][9] * mfIn9 + MInverse[7][10] * mfIn10 + MInverse[7][11] * mfIn11 + MInverse[7][12] * mfIn12 + MInverse[7][13] * mfIn13 + MInverse[7][14] * mfIn14 + MInverse[7][15] * mfIn15 + MInverse[7][16] * mfIn16 + MInverse[7][17] * mfIn17 + MInverse[7][18] * mfIn18)
#define vfGhostIn8 (MInverse[8][0]*mfIn0+MInverse[8][1]*mfIn1+MInverse[8][2]*mfIn2+MInverse[8][3]*mfIn3+MInverse[8][4]*mfIn4+MInverse[8][5]*mfIn5+MInverse[8][6]*mfIn6+MInverse[8][7]*mfIn7+MInverse[8][8]*mfIn8\
	+ MInverse[8][9] * mfIn9 + MInverse[8][10] * mfIn10 + MInverse[8][11] * mfIn11 + MInverse[8][12] * mfIn12 + MInverse[8][13] * mfIn13 + MInverse[8][14] * mfIn14 + MInverse[8][15] * mfIn15 + MInverse[8][16] * mfIn16 + MInverse[8][17] * mfIn17 + MInverse[8][18] * mfIn18)
#define vfGhostIn9 (MInverse[9][0]*mfIn0+MInverse[9][1]*mfIn1+MInverse[9][2]*mfIn2+MInverse[9][3]*mfIn3+MInverse[9][4]*mfIn4+MInverse[9][5]*mfIn5+MInverse[9][6]*mfIn6+MInverse[9][7]*mfIn7+MInverse[9][8]*mfIn8\
	+ MInverse[9][9] * mfIn9 + MInverse[9][10] * mfIn10 + MInverse[9][11] * mfIn11 + MInverse[9][12] * mfIn12 + MInverse[9][13] * mfIn13 + MInverse[9][14] * mfIn14 + MInverse[9][15] * mfIn15 + MInverse[9][16] * mfIn16 + MInverse[9][17] * mfIn17 + MInverse[9][18] * mfIn18)
#define vfGhostIn10 (MInverse[10][0]*mfIn0+MInverse[10][1]*mfIn1+MInverse[10][2]*mfIn2+MInverse[10][3]*mfIn3+MInverse[10][4]*mfIn4+MInverse[10][5]*mfIn5+MInverse[10][6]*mfIn6+MInverse[10][7]*mfIn7+MInverse[10][8]*mfIn8\
	+ MInverse[10][9] * mfIn9 + MInverse[10][10] * mfIn10 + MInverse[10][11] * mfIn11 + MInverse[10][12] * mfIn12 + MInverse[10][13] * mfIn13 + MInverse[10][14] * mfIn14 + MInverse[10][15] * mfIn15 + MInverse[10][16] * mfIn16 + MInverse[10][17] * mfIn17 + MInverse[10][18] * mfIn18)
#define vfGhostIn11 (MInverse[11][0]*mfIn0+MInverse[11][1]*mfIn1+MInverse[11][2]*mfIn2+MInverse[11][3]*mfIn3+MInverse[11][4]*mfIn4+MInverse[11][5]*mfIn5+MInverse[11][6]*mfIn6+MInverse[11][7]*mfIn7+MInverse[11][8]*mfIn8\
	+ MInverse[11][9] * mfIn9 + MInverse[11][10] * mfIn10 + MInverse[11][11] * mfIn11 + MInverse[11][12] * mfIn12 + MInverse[11][13] * mfIn13 + MInverse[11][14] * mfIn14 + MInverse[11][15] * mfIn15 + MInverse[11][16] * mfIn16 + MInverse[11][17] * mfIn17 + MInverse[11][18] * mfIn18)
#define vfGhostIn12 (MInverse[12][0]*mfIn0+MInverse[12][1]*mfIn1+MInverse[12][2]*mfIn2+MInverse[12][3]*mfIn3+MInverse[12][4]*mfIn4+MInverse[12][5]*mfIn5+MInverse[12][6]*mfIn6+MInverse[12][7]*mfIn7+MInverse[12][8]*mfIn8\
	+ MInverse[12][9] * mfIn9 + MInverse[12][10] * mfIn10 + MInverse[12][11] * mfIn11 + MInverse[12][12] * mfIn12 + MInverse[12][13] * mfIn13 + MInverse[12][14] * mfIn14 + MInverse[12][15] * mfIn15 + MInverse[12][16] * mfIn16 + MInverse[12][17] * mfIn17 + MInverse[12][18] * mfIn18)
#define vfGhostIn13 (MInverse[13][0]*mfIn0+MInverse[13][1]*mfIn1+MInverse[13][2]*mfIn2+MInverse[13][3]*mfIn3+MInverse[13][4]*mfIn4+MInverse[13][5]*mfIn5+MInverse[13][6]*mfIn6+MInverse[13][7]*mfIn7+MInverse[13][8]*mfIn8\
	+ MInverse[13][9] * mfIn9 + MInverse[13][10] * mfIn10 + MInverse[13][11] * mfIn11 + MInverse[13][12] * mfIn12 + MInverse[13][13] * mfIn13 + MInverse[13][14] * mfIn14 + MInverse[13][15] * mfIn15 + MInverse[13][16] * mfIn16 + MInverse[13][17] * mfIn17 + MInverse[13][18] * mfIn18)
#define vfGhostIn14 (MInverse[14][0]*mfIn0+MInverse[14][1]*mfIn1+MInverse[14][2]*mfIn2+MInverse[14][3]*mfIn3+MInverse[14][4]*mfIn4+MInverse[14][5]*mfIn5+MInverse[14][6]*mfIn6+MInverse[14][7]*mfIn7+MInverse[14][8]*mfIn8\
	+ MInverse[14][9] * mfIn9 + MInverse[14][10] * mfIn10 + MInverse[14][11] * mfIn11 + MInverse[14][12] * mfIn12 + MInverse[14][13] * mfIn13 + MInverse[14][14] * mfIn14 + MInverse[14][15] * mfIn15 + MInverse[14][16] * mfIn16 + MInverse[14][17] * mfIn17 + MInverse[14][18] * mfIn18)
#define vfGhostIn15 (MInverse[15][0]*mfIn0+MInverse[15][1]*mfIn1+MInverse[15][2]*mfIn2+MInverse[15][3]*mfIn3+MInverse[15][4]*mfIn4+MInverse[15][5]*mfIn5+MInverse[15][6]*mfIn6+MInverse[15][7]*mfIn7+MInverse[15][8]*mfIn8\
	+ MInverse[15][9] * mfIn9 + MInverse[15][10] * mfIn10 + MInverse[15][11] * mfIn11 + MInverse[15][12] * mfIn12 + MInverse[15][13] * mfIn13 + MInverse[15][14] * mfIn14 + MInverse[15][15] * mfIn15 + MInverse[15][16] * mfIn16 + MInverse[15][17] * mfIn17 + MInverse[15][18] * mfIn18)
#define vfGhostIn16 (MInverse[16][0]*mfIn0+MInverse[16][1]*mfIn1+MInverse[16][2]*mfIn2+MInverse[16][3]*mfIn3+MInverse[16][4]*mfIn4+MInverse[16][5]*mfIn5+MInverse[16][6]*mfIn6+MInverse[16][7]*mfIn7+MInverse[16][8]*mfIn8\
	+ MInverse[16][9] * mfIn9 + MInverse[16][10] * mfIn10 + MInverse[16][11] * mfIn11 + MInverse[16][12] * mfIn12 + MInverse[16][13] * mfIn13 + MInverse[16][14] * mfIn14 + MInverse[16][15] * mfIn15 + MInverse[16][16] * mfIn16 + MInverse[16][17] * mfIn17 + MInverse[16][18] * mfIn18)
#define vfGhostIn17 (MInverse[17][0]*mfIn0+MInverse[17][1]*mfIn1+MInverse[17][2]*mfIn2+MInverse[17][3]*mfIn3+MInverse[17][4]*mfIn4+MInverse[17][5]*mfIn5+MInverse[17][6]*mfIn6+MInverse[17][7]*mfIn7+MInverse[17][8]*mfIn8\
	+ MInverse[17][9] * mfIn9 + MInverse[17][10] * mfIn10 + MInverse[17][11] * mfIn11 + MInverse[17][12] * mfIn12 + MInverse[17][13] * mfIn13 + MInverse[17][14] * mfIn14 + MInverse[17][15] * mfIn15 + MInverse[17][16] * mfIn16 + MInverse[17][17] * mfIn17 + MInverse[17][18] * mfIn18)
#define vfGhostIn18 (MInverse[18][0]*mfIn0+MInverse[18][1]*mfIn1+MInverse[18][2]*mfIn2+MInverse[18][3]*mfIn3+MInverse[18][4]*mfIn4+MInverse[18][5]*mfIn5+MInverse[18][6]*mfIn6+MInverse[18][7]*mfIn7+MInverse[18][8]*mfIn8\
	+ MInverse[18][9] * mfIn9 + MInverse[18][10] * mfIn10 + MInverse[18][11] * mfIn11 + MInverse[18][12] * mfIn12 + MInverse[18][13] * mfIn13 + MInverse[18][14] * mfIn14 + MInverse[18][15] * mfIn15 + MInverse[18][16] * mfIn16 + MInverse[18][17] * mfIn17 + MInverse[18][18] * mfIn18)

#define vfGhostEq0 (MInverse[0][0]*mfEq0+MInverse[0][1]*mfEq1+MInverse[0][2]*mfEq2+MInverse[0][3]*mfEq3+MInverse[0][4]*mfEq4+MInverse[0][5]*mfEq5+MInverse[0][6]*mfEq6+MInverse[0][7]*mfEq7+MInverse[0][8]*mfEq8\
	+ MInverse[0][9] * mfEq9 + MInverse[0][10] * mfEq10 + MInverse[0][11] * mfEq11 + MInverse[0][12] * mfEq12 + MInverse[0][13] * mfEq13 + MInverse[0][14] * mfEq14 + MInverse[0][15] * mfEq15 + MInverse[0][16] * mfEq16 + MInverse[0][17] * mfEq17 + MInverse[0][18] * mfEq18)
#define vfGhostEq1 (MInverse[1][0]*mfEq0+MInverse[1][1]*mfEq1+MInverse[1][2]*mfEq2+MInverse[1][3]*mfEq3+MInverse[1][4]*mfEq4+MInverse[1][5]*mfEq5+MInverse[1][6]*mfEq6+MInverse[1][7]*mfEq7+MInverse[1][8]*mfEq8\
	+ MInverse[1][9] * mfEq9 + MInverse[1][10] * mfEq10 + MInverse[1][11] * mfEq11 + MInverse[1][12] * mfEq12 + MInverse[1][13] * mfEq13 + MInverse[1][14] * mfEq14 + MInverse[1][15] * mfEq15 + MInverse[1][16] * mfEq16 + MInverse[1][17] * mfEq17 + MInverse[1][18] * mfEq18)
#define vfGhostEq2 (MInverse[2][0]*mfEq0+MInverse[2][1]*mfEq1+MInverse[2][2]*mfEq2+MInverse[2][3]*mfEq3+MInverse[2][4]*mfEq4+MInverse[2][5]*mfEq5+MInverse[2][6]*mfEq6+MInverse[2][7]*mfEq7+MInverse[2][8]*mfEq8\
	+ MInverse[2][9] * mfEq9 + MInverse[2][10] * mfEq10 + MInverse[2][11] * mfEq11 + MInverse[2][12] * mfEq12 + MInverse[2][13] * mfEq13 + MInverse[2][14] * mfEq14 + MInverse[2][15] * mfEq15 + MInverse[2][16] * mfEq16 + MInverse[2][17] * mfEq17 + MInverse[2][18] * mfEq18)
#define vfGhostEq3 (MInverse[3][0]*mfEq0+MInverse[3][1]*mfEq1+MInverse[3][2]*mfEq2+MInverse[3][3]*mfEq3+MInverse[3][4]*mfEq4+MInverse[3][5]*mfEq5+MInverse[3][6]*mfEq6+MInverse[3][7]*mfEq7+MInverse[3][8]*mfEq8\
	+ MInverse[3][9] * mfEq9 + MInverse[3][10] * mfEq10 + MInverse[3][11] * mfEq11 + MInverse[3][12] * mfEq12 + MInverse[3][13] * mfEq13 + MInverse[3][14] * mfEq14 + MInverse[3][15] * mfEq15 + MInverse[3][16] * mfEq16 + MInverse[3][17] * mfEq17 + MInverse[3][18] * mfEq18)
#define vfGhostEq4 (MInverse[4][0]*mfEq0+MInverse[4][1]*mfEq1+MInverse[4][2]*mfEq2+MInverse[4][3]*mfEq3+MInverse[4][4]*mfEq4+MInverse[4][5]*mfEq5+MInverse[4][6]*mfEq6+MInverse[4][7]*mfEq7+MInverse[4][8]*mfEq8\
	+ MInverse[4][9] * mfEq9 + MInverse[4][10] * mfEq10 + MInverse[4][11] * mfEq11 + MInverse[4][12] * mfEq12 + MInverse[4][13] * mfEq13 + MInverse[4][14] * mfEq14 + MInverse[4][15] * mfEq15 + MInverse[4][16] * mfEq16 + MInverse[4][17] * mfEq17 + MInverse[4][18] * mfEq18)
#define vfGhostEq5 (MInverse[5][0]*mfEq0+MInverse[5][1]*mfEq1+MInverse[5][2]*mfEq2+MInverse[5][3]*mfEq3+MInverse[5][4]*mfEq4+MInverse[5][5]*mfEq5+MInverse[5][6]*mfEq6+MInverse[5][7]*mfEq7+MInverse[5][8]*mfEq8\
	+ MInverse[5][9] * mfEq9 + MInverse[5][10] * mfEq10 + MInverse[5][11] * mfEq11 + MInverse[5][12] * mfEq12 + MInverse[5][13] * mfEq13 + MInverse[5][14] * mfEq14 + MInverse[5][15] * mfEq15 + MInverse[5][16] * mfEq16 + MInverse[5][17] * mfEq17 + MInverse[5][18] * mfEq18)
#define vfGhostEq6 (MInverse[6][0]*mfEq0+MInverse[6][1]*mfEq1+MInverse[6][2]*mfEq2+MInverse[6][3]*mfEq3+MInverse[6][4]*mfEq4+MInverse[6][5]*mfEq5+MInverse[6][6]*mfEq6+MInverse[6][7]*mfEq7+MInverse[6][8]*mfEq8\
	+ MInverse[6][9] * mfEq9 + MInverse[6][10] * mfEq10 + MInverse[6][11] * mfEq11 + MInverse[6][12] * mfEq12 + MInverse[6][13] * mfEq13 + MInverse[6][14] * mfEq14 + MInverse[6][15] * mfEq15 + MInverse[6][16] * mfEq16 + MInverse[6][17] * mfEq17 + MInverse[6][18] * mfEq18)
#define vfGhostEq7 (MInverse[7][0]*mfEq0+MInverse[7][1]*mfEq1+MInverse[7][2]*mfEq2+MInverse[7][3]*mfEq3+MInverse[7][4]*mfEq4+MInverse[7][5]*mfEq5+MInverse[7][6]*mfEq6+MInverse[7][7]*mfEq7+MInverse[7][8]*mfEq8\
	+ MInverse[7][9] * mfEq9 + MInverse[7][10] * mfEq10 + MInverse[7][11] * mfEq11 + MInverse[7][12] * mfEq12 + MInverse[7][13] * mfEq13 + MInverse[7][14] * mfEq14 + MInverse[7][15] * mfEq15 + MInverse[7][16] * mfEq16 + MInverse[7][17] * mfEq17 + MInverse[7][18] * mfEq18)
#define vfGhostEq8 (MInverse[8][0]*mfEq0+MInverse[8][1]*mfEq1+MInverse[8][2]*mfEq2+MInverse[8][3]*mfEq3+MInverse[8][4]*mfEq4+MInverse[8][5]*mfEq5+MInverse[8][6]*mfEq6+MInverse[8][7]*mfEq7+MInverse[8][8]*mfEq8\
	+ MInverse[8][9] * mfEq9 + MInverse[8][10] * mfEq10 + MInverse[8][11] * mfEq11 + MInverse[8][12] * mfEq12 + MInverse[8][13] * mfEq13 + MInverse[8][14] * mfEq14 + MInverse[8][15] * mfEq15 + MInverse[8][16] * mfEq16 + MInverse[8][17] * mfEq17 + MInverse[8][18] * mfEq18)
#define vfGhostEq9 (MInverse[9][0]*mfEq0+MInverse[9][1]*mfEq1+MInverse[9][2]*mfEq2+MInverse[9][3]*mfEq3+MInverse[9][4]*mfEq4+MInverse[9][5]*mfEq5+MInverse[9][6]*mfEq6+MInverse[9][7]*mfEq7+MInverse[9][8]*mfEq8\
	+ MInverse[9][9] * mfEq9 + MInverse[9][10] * mfEq10 + MInverse[9][11] * mfEq11 + MInverse[9][12] * mfEq12 + MInverse[9][13] * mfEq13 + MInverse[9][14] * mfEq14 + MInverse[9][15] * mfEq15 + MInverse[9][16] * mfEq16 + MInverse[9][17] * mfEq17 + MInverse[9][18] * mfEq18)
#define vfGhostEq10 (MInverse[10][0]*mfEq0+MInverse[10][1]*mfEq1+MInverse[10][2]*mfEq2+MInverse[10][3]*mfEq3+MInverse[10][4]*mfEq4+MInverse[10][5]*mfEq5+MInverse[10][6]*mfEq6+MInverse[10][7]*mfEq7+MInverse[10][8]*mfEq8\
	+ MInverse[10][9] * mfEq9 + MInverse[10][10] * mfEq10 + MInverse[10][11] * mfEq11 + MInverse[10][12] * mfEq12 + MInverse[10][13] * mfEq13 + MInverse[10][14] * mfEq14 + MInverse[10][15] * mfEq15 + MInverse[10][16] * mfEq16 + MInverse[10][17] * mfEq17 + MInverse[10][18] * mfEq18)
#define vfGhostEq11 (MInverse[11][0]*mfEq0+MInverse[11][1]*mfEq1+MInverse[11][2]*mfEq2+MInverse[11][3]*mfEq3+MInverse[11][4]*mfEq4+MInverse[11][5]*mfEq5+MInverse[11][6]*mfEq6+MInverse[11][7]*mfEq7+MInverse[11][8]*mfEq8\
	+ MInverse[11][9] * mfEq9 + MInverse[11][10] * mfEq10 + MInverse[11][11] * mfEq11 + MInverse[11][12] * mfEq12 + MInverse[11][13] * mfEq13 + MInverse[11][14] * mfEq14 + MInverse[11][15] * mfEq15 + MInverse[11][16] * mfEq16 + MInverse[11][17] * mfEq17 + MInverse[11][18] * mfEq18)
#define vfGhostEq12 (MInverse[12][0]*mfEq0+MInverse[12][1]*mfEq1+MInverse[12][2]*mfEq2+MInverse[12][3]*mfEq3+MInverse[12][4]*mfEq4+MInverse[12][5]*mfEq5+MInverse[12][6]*mfEq6+MInverse[12][7]*mfEq7+MInverse[12][8]*mfEq8\
	+ MInverse[12][9] * mfEq9 + MInverse[12][10] * mfEq10 + MInverse[12][11] * mfEq11 + MInverse[12][12] * mfEq12 + MInverse[12][13] * mfEq13 + MInverse[12][14] * mfEq14 + MInverse[12][15] * mfEq15 + MInverse[12][16] * mfEq16 + MInverse[12][17] * mfEq17 + MInverse[12][18] * mfEq18)
#define vfGhostEq13 (MInverse[13][0]*mfEq0+MInverse[13][1]*mfEq1+MInverse[13][2]*mfEq2+MInverse[13][3]*mfEq3+MInverse[13][4]*mfEq4+MInverse[13][5]*mfEq5+MInverse[13][6]*mfEq6+MInverse[13][7]*mfEq7+MInverse[13][8]*mfEq8\
	+ MInverse[13][9] * mfEq9 + MInverse[13][10] * mfEq10 + MInverse[13][11] * mfEq11 + MInverse[13][12] * mfEq12 + MInverse[13][13] * mfEq13 + MInverse[13][14] * mfEq14 + MInverse[13][15] * mfEq15 + MInverse[13][16] * mfEq16 + MInverse[13][17] * mfEq17 + MInverse[13][18] * mfEq18)
#define vfGhostEq14 (MInverse[14][0]*mfEq0+MInverse[14][1]*mfEq1+MInverse[14][2]*mfEq2+MInverse[14][3]*mfEq3+MInverse[14][4]*mfEq4+MInverse[14][5]*mfEq5+MInverse[14][6]*mfEq6+MInverse[14][7]*mfEq7+MInverse[14][8]*mfEq8\
	+ MInverse[14][9] * mfEq9 + MInverse[14][10] * mfEq10 + MInverse[14][11] * mfEq11 + MInverse[14][12] * mfEq12 + MInverse[14][13] * mfEq13 + MInverse[14][14] * mfEq14 + MInverse[14][15] * mfEq15 + MInverse[14][16] * mfEq16 + MInverse[14][17] * mfEq17 + MInverse[14][18] * mfEq18)
#define vfGhostEq15 (MInverse[15][0]*mfEq0+MInverse[15][1]*mfEq1+MInverse[15][2]*mfEq2+MInverse[15][3]*mfEq3+MInverse[15][4]*mfEq4+MInverse[15][5]*mfEq5+MInverse[15][6]*mfEq6+MInverse[15][7]*mfEq7+MInverse[15][8]*mfEq8\
	+ MInverse[15][9] * mfEq9 + MInverse[15][10] * mfEq10 + MInverse[15][11] * mfEq11 + MInverse[15][12] * mfEq12 + MInverse[15][13] * mfEq13 + MInverse[15][14] * mfEq14 + MInverse[15][15] * mfEq15 + MInverse[15][16] * mfEq16 + MInverse[15][17] * mfEq17 + MInverse[15][18] * mfEq18)
#define vfGhostEq16 (MInverse[16][0]*mfEq0+MInverse[16][1]*mfEq1+MInverse[16][2]*mfEq2+MInverse[16][3]*mfEq3+MInverse[16][4]*mfEq4+MInverse[16][5]*mfEq5+MInverse[16][6]*mfEq6+MInverse[16][7]*mfEq7+MInverse[16][8]*mfEq8\
	+ MInverse[16][9] * mfEq9 + MInverse[16][10] * mfEq10 + MInverse[16][11] * mfEq11 + MInverse[16][12] * mfEq12 + MInverse[16][13] * mfEq13 + MInverse[16][14] * mfEq14 + MInverse[16][15] * mfEq15 + MInverse[16][16] * mfEq16 + MInverse[16][17] * mfEq17 + MInverse[16][18] * mfEq18)
#define vfGhostEq17 (MInverse[17][0]*mfEq0+MInverse[17][1]*mfEq1+MInverse[17][2]*mfEq2+MInverse[17][3]*mfEq3+MInverse[17][4]*mfEq4+MInverse[17][5]*mfEq5+MInverse[17][6]*mfEq6+MInverse[17][7]*mfEq7+MInverse[17][8]*mfEq8\
	+ MInverse[17][9] * mfEq9 + MInverse[17][10] * mfEq10 + MInverse[17][11] * mfEq11 + MInverse[17][12] * mfEq12 + MInverse[17][13] * mfEq13 + MInverse[17][14] * mfEq14 + MInverse[17][15] * mfEq15 + MInverse[17][16] * mfEq16 + MInverse[17][17] * mfEq17 + MInverse[17][18] * mfEq18)
#define vfGhostEq18 (MInverse[18][0]*mfEq0+MInverse[18][1]*mfEq1+MInverse[18][2]*mfEq2+MInverse[18][3]*mfEq3+MInverse[18][4]*mfEq4+MInverse[18][5]*mfEq5+MInverse[18][6]*mfEq6+MInverse[18][7]*mfEq7+MInverse[18][8]*mfEq8\
	+ MInverse[18][9] * mfEq9 + MInverse[18][10] * mfEq10 + MInverse[18][11] * mfEq11 + MInverse[18][12] * mfEq12 + MInverse[18][13] * mfEq13 + MInverse[18][14] * mfEq14 + MInverse[18][15] * mfEq15 + MInverse[18][16] * mfEq16 + MInverse[18][17] * mfEq17 + MInverse[18][18] * mfEq18)

#define vGhostForce0 (MInverse[0][0]*mForce0+MInverse[0][1]*mForce1+MInverse[0][2]*mForce2+MInverse[0][3]*mForce3+MInverse[0][4]*mForce4+MInverse[0][5]*mForce5+MInverse[0][6]*mForce6+MInverse[0][7]*mForce7+MInverse[0][8]*mForce8\
	+ MInverse[0][9] * mForce9 + MInverse[0][10] * mForce10 + MInverse[0][11] * mForce11 + MInverse[0][12] * mForce12 + MInverse[0][13] * mForce13 + MInverse[0][14] * mForce14 + MInverse[0][15] * mForce15 + MInverse[0][16] * mForce16 + MInverse[0][17] * mForce17 + MInverse[0][18] * mForce18)
#define vGhostForce1 (MInverse[1][0]*mForce0+MInverse[1][1]*mForce1+MInverse[1][2]*mForce2+MInverse[1][3]*mForce3+MInverse[1][4]*mForce4+MInverse[1][5]*mForce5+MInverse[1][6]*mForce6+MInverse[1][7]*mForce7+MInverse[1][8]*mForce8\
	+ MInverse[1][9] * mForce9 + MInverse[1][10] * mForce10 + MInverse[1][11] * mForce11 + MInverse[1][12] * mForce12 + MInverse[1][13] * mForce13 + MInverse[1][14] * mForce14 + MInverse[1][15] * mForce15 + MInverse[1][16] * mForce16 + MInverse[1][17] * mForce17 + MInverse[1][18] * mForce18)
#define vGhostForce2 (MInverse[2][0]*mForce0+MInverse[2][1]*mForce1+MInverse[2][2]*mForce2+MInverse[2][3]*mForce3+MInverse[2][4]*mForce4+MInverse[2][5]*mForce5+MInverse[2][6]*mForce6+MInverse[2][7]*mForce7+MInverse[2][8]*mForce8\
	+ MInverse[2][9] * mForce9 + MInverse[2][10] * mForce10 + MInverse[2][11] * mForce11 + MInverse[2][12] * mForce12 + MInverse[2][13] * mForce13 + MInverse[2][14] * mForce14 + MInverse[2][15] * mForce15 + MInverse[2][16] * mForce16 + MInverse[2][17] * mForce17 + MInverse[2][18] * mForce18)
#define vGhostForce3 (MInverse[3][0]*mForce0+MInverse[3][1]*mForce1+MInverse[3][2]*mForce2+MInverse[3][3]*mForce3+MInverse[3][4]*mForce4+MInverse[3][5]*mForce5+MInverse[3][6]*mForce6+MInverse[3][7]*mForce7+MInverse[3][8]*mForce8\
	+ MInverse[3][9] * mForce9 + MInverse[3][10] * mForce10 + MInverse[3][11] * mForce11 + MInverse[3][12] * mForce12 + MInverse[3][13] * mForce13 + MInverse[3][14] * mForce14 + MInverse[3][15] * mForce15 + MInverse[3][16] * mForce16 + MInverse[3][17] * mForce17 + MInverse[3][18] * mForce18)
#define vGhostForce4 (MInverse[4][0]*mForce0+MInverse[4][1]*mForce1+MInverse[4][2]*mForce2+MInverse[4][3]*mForce3+MInverse[4][4]*mForce4+MInverse[4][5]*mForce5+MInverse[4][6]*mForce6+MInverse[4][7]*mForce7+MInverse[4][8]*mForce8\
	+ MInverse[4][9] * mForce9 + MInverse[4][10] * mForce10 + MInverse[4][11] * mForce11 + MInverse[4][12] * mForce12 + MInverse[4][13] * mForce13 + MInverse[4][14] * mForce14 + MInverse[4][15] * mForce15 + MInverse[4][16] * mForce16 + MInverse[4][17] * mForce17 + MInverse[4][18] * mForce18)
#define vGhostForce5 (MInverse[5][0]*mForce0+MInverse[5][1]*mForce1+MInverse[5][2]*mForce2+MInverse[5][3]*mForce3+MInverse[5][4]*mForce4+MInverse[5][5]*mForce5+MInverse[5][6]*mForce6+MInverse[5][7]*mForce7+MInverse[5][8]*mForce8\
	+ MInverse[5][9] * mForce9 + MInverse[5][10] * mForce10 + MInverse[5][11] * mForce11 + MInverse[5][12] * mForce12 + MInverse[5][13] * mForce13 + MInverse[5][14] * mForce14 + MInverse[5][15] * mForce15 + MInverse[5][16] * mForce16 + MInverse[5][17] * mForce17 + MInverse[5][18] * mForce18)
#define vGhostForce6 (MInverse[6][0]*mForce0+MInverse[6][1]*mForce1+MInverse[6][2]*mForce2+MInverse[6][3]*mForce3+MInverse[6][4]*mForce4+MInverse[6][5]*mForce5+MInverse[6][6]*mForce6+MInverse[6][7]*mForce7+MInverse[6][8]*mForce8\
	+ MInverse[6][9] * mForce9 + MInverse[6][10] * mForce10 + MInverse[6][11] * mForce11 + MInverse[6][12] * mForce12 + MInverse[6][13] * mForce13 + MInverse[6][14] * mForce14 + MInverse[6][15] * mForce15 + MInverse[6][16] * mForce16 + MInverse[6][17] * mForce17 + MInverse[6][18] * mForce18)
#define vGhostForce7 (MInverse[7][0]*mForce0+MInverse[7][1]*mForce1+MInverse[7][2]*mForce2+MInverse[7][3]*mForce3+MInverse[7][4]*mForce4+MInverse[7][5]*mForce5+MInverse[7][6]*mForce6+MInverse[7][7]*mForce7+MInverse[7][8]*mForce8\
	+ MInverse[7][9] * mForce9 + MInverse[7][10] * mForce10 + MInverse[7][11] * mForce11 + MInverse[7][12] * mForce12 + MInverse[7][13] * mForce13 + MInverse[7][14] * mForce14 + MInverse[7][15] * mForce15 + MInverse[7][16] * mForce16 + MInverse[7][17] * mForce17 + MInverse[7][18] * mForce18)
#define vGhostForce8 (MInverse[8][0]*mForce0+MInverse[8][1]*mForce1+MInverse[8][2]*mForce2+MInverse[8][3]*mForce3+MInverse[8][4]*mForce4+MInverse[8][5]*mForce5+MInverse[8][6]*mForce6+MInverse[8][7]*mForce7+MInverse[8][8]*mForce8\
	+ MInverse[8][9] * mForce9 + MInverse[8][10] * mForce10 + MInverse[8][11] * mForce11 + MInverse[8][12] * mForce12 + MInverse[8][13] * mForce13 + MInverse[8][14] * mForce14 + MInverse[8][15] * mForce15 + MInverse[8][16] * mForce16 + MInverse[8][17] * mForce17 + MInverse[8][18] * mForce18)
#define vGhostForce9 (MInverse[9][0]*mForce0+MInverse[9][1]*mForce1+MInverse[9][2]*mForce2+MInverse[9][3]*mForce3+MInverse[9][4]*mForce4+MInverse[9][5]*mForce5+MInverse[9][6]*mForce6+MInverse[9][7]*mForce7+MInverse[9][8]*mForce8\
	+ MInverse[9][9] * mForce9 + MInverse[9][10] * mForce10 + MInverse[9][11] * mForce11 + MInverse[9][12] * mForce12 + MInverse[9][13] * mForce13 + MInverse[9][14] * mForce14 + MInverse[9][15] * mForce15 + MInverse[9][16] * mForce16 + MInverse[9][17] * mForce17 + MInverse[9][18] * mForce18)
#define vGhostForce10 (MInverse[10][0]*mForce0+MInverse[10][1]*mForce1+MInverse[10][2]*mForce2+MInverse[10][3]*mForce3+MInverse[10][4]*mForce4+MInverse[10][5]*mForce5+MInverse[10][6]*mForce6+MInverse[10][7]*mForce7+MInverse[10][8]*mForce8\
	+ MInverse[10][9] * mForce9 + MInverse[10][10] * mForce10 + MInverse[10][11] * mForce11 + MInverse[10][12] * mForce12 + MInverse[10][13] * mForce13 + MInverse[10][14] * mForce14 + MInverse[10][15] * mForce15 + MInverse[10][16] * mForce16 + MInverse[10][17] * mForce17 + MInverse[10][18] * mForce18)
#define vGhostForce11 (MInverse[11][0]*mForce0+MInverse[11][1]*mForce1+MInverse[11][2]*mForce2+MInverse[11][3]*mForce3+MInverse[11][4]*mForce4+MInverse[11][5]*mForce5+MInverse[11][6]*mForce6+MInverse[11][7]*mForce7+MInverse[11][8]*mForce8\
	+ MInverse[11][9] * mForce9 + MInverse[11][10] * mForce10 + MInverse[11][11] * mForce11 + MInverse[11][12] * mForce12 + MInverse[11][13] * mForce13 + MInverse[11][14] * mForce14 + MInverse[11][15] * mForce15 + MInverse[11][16] * mForce16 + MInverse[11][17] * mForce17 + MInverse[11][18] * mForce18)
#define vGhostForce12 (MInverse[12][0]*mForce0+MInverse[12][1]*mForce1+MInverse[12][2]*mForce2+MInverse[12][3]*mForce3+MInverse[12][4]*mForce4+MInverse[12][5]*mForce5+MInverse[12][6]*mForce6+MInverse[12][7]*mForce7+MInverse[12][8]*mForce8\
	+ MInverse[12][9] * mForce9 + MInverse[12][10] * mForce10 + MInverse[12][11] * mForce11 + MInverse[12][12] * mForce12 + MInverse[12][13] * mForce13 + MInverse[12][14] * mForce14 + MInverse[12][15] * mForce15 + MInverse[12][16] * mForce16 + MInverse[12][17] * mForce17 + MInverse[12][18] * mForce18)
#define vGhostForce13 (MInverse[13][0]*mForce0+MInverse[13][1]*mForce1+MInverse[13][2]*mForce2+MInverse[13][3]*mForce3+MInverse[13][4]*mForce4+MInverse[13][5]*mForce5+MInverse[13][6]*mForce6+MInverse[13][7]*mForce7+MInverse[13][8]*mForce8\
	+ MInverse[13][9] * mForce9 + MInverse[13][10] * mForce10 + MInverse[13][11] * mForce11 + MInverse[13][12] * mForce12 + MInverse[13][13] * mForce13 + MInverse[13][14] * mForce14 + MInverse[13][15] * mForce15 + MInverse[13][16] * mForce16 + MInverse[13][17] * mForce17 + MInverse[13][18] * mForce18)
#define vGhostForce14 (MInverse[14][0]*mForce0+MInverse[14][1]*mForce1+MInverse[14][2]*mForce2+MInverse[14][3]*mForce3+MInverse[14][4]*mForce4+MInverse[14][5]*mForce5+MInverse[14][6]*mForce6+MInverse[14][7]*mForce7+MInverse[14][8]*mForce8\
	+ MInverse[14][9] * mForce9 + MInverse[14][10] * mForce10 + MInverse[14][11] * mForce11 + MInverse[14][12] * mForce12 + MInverse[14][13] * mForce13 + MInverse[14][14] * mForce14 + MInverse[14][15] * mForce15 + MInverse[14][16] * mForce16 + MInverse[14][17] * mForce17 + MInverse[14][18] * mForce18)
#define vGhostForce15 (MInverse[15][0]*mForce0+MInverse[15][1]*mForce1+MInverse[15][2]*mForce2+MInverse[15][3]*mForce3+MInverse[15][4]*mForce4+MInverse[15][5]*mForce5+MInverse[15][6]*mForce6+MInverse[15][7]*mForce7+MInverse[15][8]*mForce8\
	+ MInverse[15][9] * mForce9 + MInverse[15][10] * mForce10 + MInverse[15][11] * mForce11 + MInverse[15][12] * mForce12 + MInverse[15][13] * mForce13 + MInverse[15][14] * mForce14 + MInverse[15][15] * mForce15 + MInverse[15][16] * mForce16 + MInverse[15][17] * mForce17 + MInverse[15][18] * mForce18)
#define vGhostForce16 (MInverse[16][0]*mForce0+MInverse[16][1]*mForce1+MInverse[16][2]*mForce2+MInverse[16][3]*mForce3+MInverse[16][4]*mForce4+MInverse[16][5]*mForce5+MInverse[16][6]*mForce6+MInverse[16][7]*mForce7+MInverse[16][8]*mForce8\
	+ MInverse[16][9] * mForce9 + MInverse[16][10] * mForce10 + MInverse[16][11] * mForce11 + MInverse[16][12] * mForce12 + MInverse[16][13] * mForce13 + MInverse[16][14] * mForce14 + MInverse[16][15] * mForce15 + MInverse[16][16] * mForce16 + MInverse[16][17] * mForce17 + MInverse[16][18] * mForce18)
#define vGhostForce17 (MInverse[17][0]*mForce0+MInverse[17][1]*mForce1+MInverse[17][2]*mForce2+MInverse[17][3]*mForce3+MInverse[17][4]*mForce4+MInverse[17][5]*mForce5+MInverse[17][6]*mForce6+MInverse[17][7]*mForce7+MInverse[17][8]*mForce8\
	+ MInverse[17][9] * mForce9 + MInverse[17][10] * mForce10 + MInverse[17][11] * mForce11 + MInverse[17][12] * mForce12 + MInverse[17][13] * mForce13 + MInverse[17][14] * mForce14 + MInverse[17][15] * mForce15 + MInverse[17][16] * mForce16 + MInverse[17][17] * mForce17 + MInverse[17][18] * mForce18)
#define vGhostForce18 (MInverse[18][0]*mForce0+MInverse[18][1]*mForce1+MInverse[18][2]*mForce2+MInverse[18][3]*mForce3+MInverse[18][4]*mForce4+MInverse[18][5]*mForce5+MInverse[18][6]*mForce6+MInverse[18][7]*mForce7+MInverse[18][8]*mForce8\
	+ MInverse[18][9] * mForce9 + MInverse[18][10] * mForce10 + MInverse[18][11] * mForce11 + MInverse[18][12] * mForce12 + MInverse[18][13] * mForce13 + MInverse[18][14] * mForce14 + MInverse[18][15] * mForce15 + MInverse[18][16] * mForce16 + MInverse[18][17] * mForce17 + MInverse[18][18] * mForce18)

double ***fIn0Ghost=NULL;
double ***fIn1Ghost=NULL;
double ***fIn2Ghost=NULL;
double ***fIn3Ghost=NULL;
double ***fIn4Ghost=NULL;
double ***fIn5Ghost=NULL;
double ***fIn6Ghost=NULL;
double ***fIn7Ghost=NULL;
double ***fIn8Ghost=NULL;
double ***fIn9Ghost=NULL;
double ***fIn10Ghost=NULL;
double ***fIn11Ghost=NULL;
double ***fIn12Ghost=NULL;
double ***fIn13Ghost=NULL;
double ***fIn14Ghost=NULL;
double ***fIn15Ghost=NULL;
double ***fIn16Ghost=NULL;
double ***fIn17Ghost=NULL;
double ***fIn18Ghost=NULL;
double ***rhoGhost=NULL;
double ***uxGhost=NULL;
double ***uyGhost=NULL;
double ***uzGhost=NULL;

double AGhost=1.0/(viscosityGhost*3+0.5);

void FieldArrangeGhostVelocity(int Length, int Width, int Height){
	fIn0Ghost=memoryarrange(Length,Width,Height);
	fIn1Ghost=memoryarrange(Length,Width,Height);
	fIn2Ghost=memoryarrange(Length,Width,Height);
	fIn3Ghost=memoryarrange(Length,Width,Height);
	fIn4Ghost=memoryarrange(Length,Width,Height);
	fIn5Ghost=memoryarrange(Length,Width,Height);
	fIn6Ghost=memoryarrange(Length,Width,Height);
	fIn7Ghost=memoryarrange(Length,Width,Height);
	fIn8Ghost=memoryarrange(Length,Width,Height);
	fIn9Ghost=memoryarrange(Length,Width,Height);
	fIn10Ghost=memoryarrange(Length,Width,Height);
	fIn11Ghost=memoryarrange(Length,Width,Height);
	fIn12Ghost=memoryarrange(Length,Width,Height);
	fIn13Ghost=memoryarrange(Length,Width,Height);
	fIn14Ghost=memoryarrange(Length,Width,Height);
	fIn15Ghost=memoryarrange(Length,Width,Height);
	fIn16Ghost=memoryarrange(Length,Width,Height);
	fIn17Ghost=memoryarrange(Length,Width,Height);
	fIn18Ghost=memoryarrange(Length,Width,Height);
	rhoGhost=memoryarrange(Length,Width,Height); 
	uxGhost=memoryarrange(Length,Width,Height); 
	uyGhost=memoryarrange(Length,Width,Height); 
	uzGhost=memoryarrange(Length,Width,Height); 
}

void FieldInitialGhostVelocity(int Length, int Width, int Height, double GAS, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical){
						rhoGhost[i][j][k]=GAS;
						fIn0Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[0];
						fIn1Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[1];
						fIn2Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[2];
						fIn3Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[3];
						fIn4Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[4];
						fIn5Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[5];
						fIn6Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[6];
						fIn7Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[7];
						fIn8Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[8];
						fIn9Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[9];
						fIn10Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[10];
						fIn11Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[11];
						fIn12Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[12];
						fIn13Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[13];
						fIn14Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[14];
						fIn15Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[15];
						fIn16Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[16];
						fIn17Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[17];
						fIn18Ghost[i][j][k]=rhoGhost[i][j][k]*tNS[18];
						uxGhost[i][j][k]=uTotXGhost;
						uyGhost[i][j][k]=uTotYGhost;
						uzGhost[i][j][k]=uTotZGhost;
					}
					if (domain[i][j][k]!=0||rho[i][j][k]<=gasCritical){
						rhoGhost[i][j][k]=0.0;
						uxGhost[i][j][k]=0.0;
						uyGhost[i][j][k]=0.0;
						uzGhost[i][j][k]=0.0;
						fIn0Ghost[i][j][k]=0.0;
						fIn1Ghost[i][j][k]=0.0;
						fIn2Ghost[i][j][k]=0.0;
						fIn3Ghost[i][j][k]=0.0;
						fIn4Ghost[i][j][k]=0.0;
						fIn5Ghost[i][j][k]=0.0;
						fIn6Ghost[i][j][k]=0.0;
						fIn7Ghost[i][j][k]=0.0;
						fIn8Ghost[i][j][k]=0.0;
						fIn9Ghost[i][j][k]=0.0;
						fIn10Ghost[i][j][k]=0.0;
						fIn11Ghost[i][j][k]=0.0;
						fIn12Ghost[i][j][k]=0.0;
						fIn13Ghost[i][j][k]=0.0;
						fIn14Ghost[i][j][k]=0.0;
						fIn15Ghost[i][j][k]=0.0;
						fIn16Ghost[i][j][k]=0.0;
						fIn17Ghost[i][j][k]=0.0;
						fIn18Ghost[i][j][k]=0.0;
					}
				}
			}
		}
	}
}

void CollisionSRTGhostVelocity(int Length,int Width,int Height,double gasCritical){
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
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical){
						cfNS0=czNS[0]*bodyForce[i][j][k];
						cfNS1=czNS[1]*bodyForce[i][j][k];
						cfNS2=czNS[2]*bodyForce[i][j][k];
						cfNS3=czNS[3]*bodyForce[i][j][k];
						cfNS4=czNS[4]*bodyForce[i][j][k];
						cfNS5=czNS[5]*bodyForce[i][j][k];
						cfNS6=czNS[6]*bodyForce[i][j][k];
						cfNS7=czNS[7]*bodyForce[i][j][k];
						cfNS8=czNS[8]*bodyForce[i][j][k];
						cfNS9=czNS[9]*bodyForce[i][j][k];
						cfNS10=czNS[10]*bodyForce[i][j][k];
						cfNS11=czNS[11]*bodyForce[i][j][k];
						cfNS12=czNS[12]*bodyForce[i][j][k];
						cfNS13=czNS[13]*bodyForce[i][j][k];
						cfNS14=czNS[14]*bodyForce[i][j][k];
						cfNS15=czNS[15]*bodyForce[i][j][k];
						cfNS16=czNS[16]*bodyForce[i][j][k];
						cfNS17=czNS[17]*bodyForce[i][j][k];
						cfNS18=czNS[18]*bodyForce[i][j][k];
						uxFxuyFyuzFz=uzGhost[i][j][k]*bodyForce[i][j][k];
						Out0[i][j][k]=fIn0Ghost[i][j][k]-AGhost*(fIn0Ghost[i][j][k]-Eq0[i][j][k]*rhoGhost[i][j][k])+Force0Ghost;
						Out1[i][j][k]=fIn1Ghost[i][j][k]-AGhost*(fIn1Ghost[i][j][k]-Eq1[i][j][k]*rhoGhost[i][j][k])+Force1Ghost;
						Out2[i][j][k]=fIn2Ghost[i][j][k]-AGhost*(fIn2Ghost[i][j][k]-Eq2[i][j][k]*rhoGhost[i][j][k])+Force2Ghost;
						Out3[i][j][k]=fIn3Ghost[i][j][k]-AGhost*(fIn3Ghost[i][j][k]-Eq3[i][j][k]*rhoGhost[i][j][k])+Force3Ghost;
						Out4[i][j][k]=fIn4Ghost[i][j][k]-AGhost*(fIn4Ghost[i][j][k]-Eq4[i][j][k]*rhoGhost[i][j][k])+Force4Ghost;
						Out5[i][j][k]=fIn5Ghost[i][j][k]-AGhost*(fIn5Ghost[i][j][k]-Eq5[i][j][k]*rhoGhost[i][j][k])+Force5Ghost;
						Out6[i][j][k]=fIn6Ghost[i][j][k]-AGhost*(fIn6Ghost[i][j][k]-Eq6[i][j][k]*rhoGhost[i][j][k])+Force6Ghost;
						Out7[i][j][k]=fIn7Ghost[i][j][k]-AGhost*(fIn7Ghost[i][j][k]-Eq7[i][j][k]*rhoGhost[i][j][k])+Force7Ghost;
						Out8[i][j][k]=fIn8Ghost[i][j][k]-AGhost*(fIn8Ghost[i][j][k]-Eq8[i][j][k]*rhoGhost[i][j][k])+Force8Ghost;
						Out9[i][j][k]=fIn9Ghost[i][j][k]-AGhost*(fIn9Ghost[i][j][k]-Eq9[i][j][k]*rhoGhost[i][j][k])+Force9Ghost;
						Out10[i][j][k]=fIn10Ghost[i][j][k]-AGhost*(fIn10Ghost[i][j][k]-Eq10[i][j][k]*rhoGhost[i][j][k])+Force10Ghost;
						Out11[i][j][k]=fIn11Ghost[i][j][k]-AGhost*(fIn11Ghost[i][j][k]-Eq11[i][j][k]*rhoGhost[i][j][k])+Force11Ghost;
						Out12[i][j][k]=fIn12Ghost[i][j][k]-AGhost*(fIn12Ghost[i][j][k]-Eq12[i][j][k]*rhoGhost[i][j][k])+Force12Ghost;
						Out13[i][j][k]=fIn13Ghost[i][j][k]-AGhost*(fIn13Ghost[i][j][k]-Eq13[i][j][k]*rhoGhost[i][j][k])+Force13Ghost;
						Out14[i][j][k]=fIn14Ghost[i][j][k]-AGhost*(fIn14Ghost[i][j][k]-Eq14[i][j][k]*rhoGhost[i][j][k])+Force14Ghost;
						Out15[i][j][k]=fIn15Ghost[i][j][k]-AGhost*(fIn15Ghost[i][j][k]-Eq15[i][j][k]*rhoGhost[i][j][k])+Force15Ghost;
						Out16[i][j][k]=fIn16Ghost[i][j][k]-AGhost*(fIn16Ghost[i][j][k]-Eq16[i][j][k]*rhoGhost[i][j][k])+Force16Ghost;
						Out17[i][j][k]=fIn17Ghost[i][j][k]-AGhost*(fIn17Ghost[i][j][k]-Eq17[i][j][k]*rhoGhost[i][j][k])+Force17Ghost;
						Out18[i][j][k]=fIn18Ghost[i][j][k]-AGhost*(fIn18Ghost[i][j][k]-Eq18[i][j][k]*rhoGhost[i][j][k])+Force18Ghost;
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

void CollisionMRTGhostVelocity(int m, int n, int q, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k,mfIn0,mfIn1,mfIn2,mfIn3,mfIn4,mfIn5,mfIn6,mfIn7,mfIn8,mfIn9,\
	mfIn10, mfIn11, mfIn12, mfIn13, mfIn14, mfIn15, mfIn16, mfIn17, mfIn18, mfEq0, mfEq1, mfEq2, mfEq3, mfEq4, mfEq5, mfEq6, mfEq7, mfEq8, mfEq9, mfEq10, mfEq11, mfEq12, mfEq13, mfEq14, \
	mfEq15, mfEq16, mfEq17, mfEq18, mForce0, mForce1, mForce2, mForce3, mForce4, mForce5, mForce6, mForce7, mForce8, mForce9, mForce10, mForce11, mForce12, mForce13, mForce14, mForce15, \
	mForce16, mForce17, mForce18)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<m; i++){
			for (j = 0; j<n; j++){
				for (k = 0; k<q; k++){
					if (domain[i][j][k] == 0 && rho[i][j][k]>gasCritical){
						mfIn0 = (A_g[0][0] * (M[0][0] * fIn0Ghost[i][j][k] + M[0][1] * fIn1Ghost[i][j][k] + M[0][2] * fIn2Ghost[i][j][k] + M[0][3] * fIn3Ghost[i][j][k] + M[0][4] * fIn4Ghost[i][j][k] + M[0][5] * fIn5Ghost[i][j][k] + M[0][6] * fIn6Ghost[i][j][k] + M[0][7] * fIn7Ghost[i][j][k] + M[0][8] * fIn8Ghost[i][j][k]\
							+ M[0][9] * fIn9Ghost[i][j][k] + M[0][10] * fIn10Ghost[i][j][k] + M[0][11] * fIn11Ghost[i][j][k] + M[0][12] * fIn12Ghost[i][j][k] + M[0][13] * fIn13Ghost[i][j][k] + M[0][14] * fIn14Ghost[i][j][k] + M[0][15] * fIn15Ghost[i][j][k] + M[0][16] * fIn16Ghost[i][j][k] + M[0][17] * fIn17Ghost[i][j][k] + M[0][18] * fIn18Ghost[i][j][k]));
						mfIn1 = (A_g[1][1] * (M[1][0] * fIn0Ghost[i][j][k] + M[1][1] * fIn1Ghost[i][j][k] + M[1][2] * fIn2Ghost[i][j][k] + M[1][3] * fIn3Ghost[i][j][k] + M[1][4] * fIn4Ghost[i][j][k] + M[1][5] * fIn5Ghost[i][j][k] + M[1][6] * fIn6Ghost[i][j][k] + M[1][7] * fIn7Ghost[i][j][k] + M[1][8] * fIn8Ghost[i][j][k]\
							+ M[1][9] * fIn9Ghost[i][j][k] + M[1][10] * fIn10Ghost[i][j][k] + M[1][11] * fIn11Ghost[i][j][k] + M[1][12] * fIn12Ghost[i][j][k] + M[1][13] * fIn13Ghost[i][j][k] + M[1][14] * fIn14Ghost[i][j][k] + M[1][15] * fIn15Ghost[i][j][k] + M[1][16] * fIn16Ghost[i][j][k] + M[1][17] * fIn17Ghost[i][j][k] + M[1][18] * fIn18Ghost[i][j][k]));
						mfIn2 = (A_g[2][2] * (M[2][0] * fIn0Ghost[i][j][k] + M[2][1] * fIn1Ghost[i][j][k] + M[2][2] * fIn2Ghost[i][j][k] + M[2][3] * fIn3Ghost[i][j][k] + M[2][4] * fIn4Ghost[i][j][k] + M[2][5] * fIn5Ghost[i][j][k] + M[2][6] * fIn6Ghost[i][j][k] + M[2][7] * fIn7Ghost[i][j][k] + M[2][8] * fIn8Ghost[i][j][k]\
							+ M[2][9] * fIn9Ghost[i][j][k] + M[2][10] * fIn10Ghost[i][j][k] + M[2][11] * fIn11Ghost[i][j][k] + M[2][12] * fIn12Ghost[i][j][k] + M[2][13] * fIn13Ghost[i][j][k] + M[2][14] * fIn14Ghost[i][j][k] + M[2][15] * fIn15Ghost[i][j][k] + M[2][16] * fIn16Ghost[i][j][k] + M[2][17] * fIn17Ghost[i][j][k] + M[2][18] * fIn18Ghost[i][j][k]));
						mfIn3 = (A_g[3][3] * (M[3][0] * fIn0Ghost[i][j][k] + M[3][1] * fIn1Ghost[i][j][k] + M[3][2] * fIn2Ghost[i][j][k] + M[3][3] * fIn3Ghost[i][j][k] + M[3][4] * fIn4Ghost[i][j][k] + M[3][5] * fIn5Ghost[i][j][k] + M[3][6] * fIn6Ghost[i][j][k] + M[3][7] * fIn7Ghost[i][j][k] + M[3][8] * fIn8Ghost[i][j][k]\
							+ M[3][9] * fIn9Ghost[i][j][k] + M[3][10] * fIn10Ghost[i][j][k] + M[3][11] * fIn11Ghost[i][j][k] + M[3][12] * fIn12Ghost[i][j][k] + M[3][13] * fIn13Ghost[i][j][k] + M[3][14] * fIn14Ghost[i][j][k] + M[3][15] * fIn15Ghost[i][j][k] + M[3][16] * fIn16Ghost[i][j][k] + M[3][17] * fIn17Ghost[i][j][k] + M[3][18] * fIn18Ghost[i][j][k]));
						mfIn4 = (A_g[4][4] * (M[4][0] * fIn0Ghost[i][j][k] + M[4][1] * fIn1Ghost[i][j][k] + M[4][2] * fIn2Ghost[i][j][k] + M[4][3] * fIn3Ghost[i][j][k] + M[4][4] * fIn4Ghost[i][j][k] + M[4][5] * fIn5Ghost[i][j][k] + M[4][6] * fIn6Ghost[i][j][k] + M[4][7] * fIn7Ghost[i][j][k] + M[4][8] * fIn8Ghost[i][j][k]\
							+ M[4][9] * fIn9Ghost[i][j][k] + M[4][10] * fIn10Ghost[i][j][k] + M[4][11] * fIn11Ghost[i][j][k] + M[4][12] * fIn12Ghost[i][j][k] + M[4][13] * fIn13Ghost[i][j][k] + M[4][14] * fIn14Ghost[i][j][k] + M[4][15] * fIn15Ghost[i][j][k] + M[4][16] * fIn16Ghost[i][j][k] + M[4][17] * fIn17Ghost[i][j][k] + M[4][18] * fIn18Ghost[i][j][k]));
						mfIn5 = (A_g[5][5] * (M[5][0] * fIn0Ghost[i][j][k] + M[5][1] * fIn1Ghost[i][j][k] + M[5][2] * fIn2Ghost[i][j][k] + M[5][3] * fIn3Ghost[i][j][k] + M[5][4] * fIn4Ghost[i][j][k] + M[5][5] * fIn5Ghost[i][j][k] + M[5][6] * fIn6Ghost[i][j][k] + M[5][7] * fIn7Ghost[i][j][k] + M[5][8] * fIn8Ghost[i][j][k]\
							+ M[5][9] * fIn9Ghost[i][j][k] + M[5][10] * fIn10Ghost[i][j][k] + M[5][11] * fIn11Ghost[i][j][k] + M[5][12] * fIn12Ghost[i][j][k] + M[5][13] * fIn13Ghost[i][j][k] + M[5][14] * fIn14Ghost[i][j][k] + M[5][15] * fIn15Ghost[i][j][k] + M[5][16] * fIn16Ghost[i][j][k] + M[5][17] * fIn17Ghost[i][j][k] + M[5][18] * fIn18Ghost[i][j][k]));
						mfIn6 = (A_g[6][6] * (M[6][0] * fIn0Ghost[i][j][k] + M[6][1] * fIn1Ghost[i][j][k] + M[6][2] * fIn2Ghost[i][j][k] + M[6][3] * fIn3Ghost[i][j][k] + M[6][4] * fIn4Ghost[i][j][k] + M[6][5] * fIn5Ghost[i][j][k] + M[6][6] * fIn6Ghost[i][j][k] + M[6][7] * fIn7Ghost[i][j][k] + M[6][8] * fIn8Ghost[i][j][k]\
							+ M[6][9] * fIn9Ghost[i][j][k] + M[6][10] * fIn10Ghost[i][j][k] + M[6][11] * fIn11Ghost[i][j][k] + M[6][12] * fIn12Ghost[i][j][k] + M[6][13] * fIn13Ghost[i][j][k] + M[6][14] * fIn14Ghost[i][j][k] + M[6][15] * fIn15Ghost[i][j][k] + M[6][16] * fIn16Ghost[i][j][k] + M[6][17] * fIn17Ghost[i][j][k] + M[6][18] * fIn18Ghost[i][j][k]));
						mfIn7 = (A_g[7][7] * (M[7][0] * fIn0Ghost[i][j][k] + M[7][1] * fIn1Ghost[i][j][k] + M[7][2] * fIn2Ghost[i][j][k] + M[7][3] * fIn3Ghost[i][j][k] + M[7][4] * fIn4Ghost[i][j][k] + M[7][5] * fIn5Ghost[i][j][k] + M[7][6] * fIn6Ghost[i][j][k] + M[7][7] * fIn7Ghost[i][j][k] + M[7][8] * fIn8Ghost[i][j][k]\
							+ M[7][9] * fIn9Ghost[i][j][k] + M[7][10] * fIn10Ghost[i][j][k] + M[7][11] * fIn11Ghost[i][j][k] + M[7][12] * fIn12Ghost[i][j][k] + M[7][13] * fIn13Ghost[i][j][k] + M[7][14] * fIn14Ghost[i][j][k] + M[7][15] * fIn15Ghost[i][j][k] + M[7][16] * fIn16Ghost[i][j][k] + M[7][17] * fIn17Ghost[i][j][k] + M[7][18] * fIn18Ghost[i][j][k]));
						mfIn8 = (A_g[8][8] * (M[8][0] * fIn0Ghost[i][j][k] + M[8][1] * fIn1Ghost[i][j][k] + M[8][2] * fIn2Ghost[i][j][k] + M[8][3] * fIn3Ghost[i][j][k] + M[8][4] * fIn4Ghost[i][j][k] + M[8][5] * fIn5Ghost[i][j][k] + M[8][6] * fIn6Ghost[i][j][k] + M[8][7] * fIn7Ghost[i][j][k] + M[8][8] * fIn8Ghost[i][j][k]\
							+ M[8][9] * fIn9Ghost[i][j][k] + M[8][10] * fIn10Ghost[i][j][k] + M[8][11] * fIn11Ghost[i][j][k] + M[8][12] * fIn12Ghost[i][j][k] + M[8][13] * fIn13Ghost[i][j][k] + M[8][14] * fIn14Ghost[i][j][k] + M[8][15] * fIn15Ghost[i][j][k] + M[8][16] * fIn16Ghost[i][j][k] + M[8][17] * fIn17Ghost[i][j][k] + M[8][18] * fIn18Ghost[i][j][k]));
						mfIn9 = (A_g[9][9] * (M[9][0] * fIn0Ghost[i][j][k] + M[9][1] * fIn1Ghost[i][j][k] + M[9][2] * fIn2Ghost[i][j][k] + M[9][3] * fIn3Ghost[i][j][k] + M[9][4] * fIn4Ghost[i][j][k] + M[9][5] * fIn5Ghost[i][j][k] + M[9][6] * fIn6Ghost[i][j][k] + M[9][7] * fIn7Ghost[i][j][k] + M[9][8] * fIn8Ghost[i][j][k]\
							+ M[9][9] * fIn9Ghost[i][j][k] + M[9][10] * fIn10Ghost[i][j][k] + M[9][11] * fIn11Ghost[i][j][k] + M[9][12] * fIn12Ghost[i][j][k] + M[9][13] * fIn13Ghost[i][j][k] + M[9][14] * fIn14Ghost[i][j][k] + M[9][15] * fIn15Ghost[i][j][k] + M[9][16] * fIn16Ghost[i][j][k] + M[9][17] * fIn17Ghost[i][j][k] + M[9][18] * fIn18Ghost[i][j][k]));
						mfIn10 = (A_g[10][10] * (M[10][0] * fIn0Ghost[i][j][k] + M[10][1] * fIn1Ghost[i][j][k] + M[10][2] * fIn2Ghost[i][j][k] + M[10][3] * fIn3Ghost[i][j][k] + M[10][4] * fIn4Ghost[i][j][k] + M[10][5] * fIn5Ghost[i][j][k] + M[10][6] * fIn6Ghost[i][j][k] + M[10][7] * fIn7Ghost[i][j][k] + M[10][8] * fIn8Ghost[i][j][k]\
							+ M[10][9] * fIn9Ghost[i][j][k] + M[10][10] * fIn10Ghost[i][j][k] + M[10][11] * fIn11Ghost[i][j][k] + M[10][12] * fIn12Ghost[i][j][k] + M[10][13] * fIn13Ghost[i][j][k] + M[10][14] * fIn14Ghost[i][j][k] + M[10][15] * fIn15Ghost[i][j][k] + M[10][16] * fIn16Ghost[i][j][k] + M[10][17] * fIn17Ghost[i][j][k] + M[10][18] * fIn18Ghost[i][j][k]));
						mfIn11 = (A_g[11][11] * (M[11][0] * fIn0Ghost[i][j][k] + M[11][1] * fIn1Ghost[i][j][k] + M[11][2] * fIn2Ghost[i][j][k] + M[11][3] * fIn3Ghost[i][j][k] + M[11][4] * fIn4Ghost[i][j][k] + M[11][5] * fIn5Ghost[i][j][k] + M[11][6] * fIn6Ghost[i][j][k] + M[11][7] * fIn7Ghost[i][j][k] + M[11][8] * fIn8Ghost[i][j][k]\
							+ M[11][9] * fIn9Ghost[i][j][k] + M[11][10] * fIn10Ghost[i][j][k] + M[11][11] * fIn11Ghost[i][j][k] + M[11][12] * fIn12Ghost[i][j][k] + M[11][13] * fIn13Ghost[i][j][k] + M[11][14] * fIn14Ghost[i][j][k] + M[11][15] * fIn15Ghost[i][j][k] + M[11][16] * fIn16Ghost[i][j][k] + M[11][17] * fIn17Ghost[i][j][k] + M[11][18] * fIn18Ghost[i][j][k]));
						mfIn12 = (A_g[12][12] * (M[12][0] * fIn0Ghost[i][j][k] + M[12][1] * fIn1Ghost[i][j][k] + M[12][2] * fIn2Ghost[i][j][k] + M[12][3] * fIn3Ghost[i][j][k] + M[12][4] * fIn4Ghost[i][j][k] + M[12][5] * fIn5Ghost[i][j][k] + M[12][6] * fIn6Ghost[i][j][k] + M[12][7] * fIn7Ghost[i][j][k] + M[12][8] * fIn8Ghost[i][j][k]\
							+ M[12][9] * fIn9Ghost[i][j][k] + M[12][10] * fIn10Ghost[i][j][k] + M[12][11] * fIn11Ghost[i][j][k] + M[12][12] * fIn12Ghost[i][j][k] + M[12][13] * fIn13Ghost[i][j][k] + M[12][14] * fIn14Ghost[i][j][k] + M[12][15] * fIn15Ghost[i][j][k] + M[12][16] * fIn16Ghost[i][j][k] + M[12][17] * fIn17Ghost[i][j][k] + M[12][18] * fIn18Ghost[i][j][k]));
						mfIn13 = (A_g[13][13] * (M[13][0] * fIn0Ghost[i][j][k] + M[13][1] * fIn1Ghost[i][j][k] + M[13][2] * fIn2Ghost[i][j][k] + M[13][3] * fIn3Ghost[i][j][k] + M[13][4] * fIn4Ghost[i][j][k] + M[13][5] * fIn5Ghost[i][j][k] + M[13][6] * fIn6Ghost[i][j][k] + M[13][7] * fIn7Ghost[i][j][k] + M[13][8] * fIn8Ghost[i][j][k]\
							+ M[13][9] * fIn9Ghost[i][j][k] + M[13][10] * fIn10Ghost[i][j][k] + M[13][11] * fIn11Ghost[i][j][k] + M[13][12] * fIn12Ghost[i][j][k] + M[13][13] * fIn13Ghost[i][j][k] + M[13][14] * fIn14Ghost[i][j][k] + M[13][15] * fIn15Ghost[i][j][k] + M[13][16] * fIn16Ghost[i][j][k] + M[13][17] * fIn17Ghost[i][j][k] + M[13][18] * fIn18Ghost[i][j][k]));
						mfIn14 = (A_g[14][14] * (M[14][0] * fIn0Ghost[i][j][k] + M[14][1] * fIn1Ghost[i][j][k] + M[14][2] * fIn2Ghost[i][j][k] + M[14][3] * fIn3Ghost[i][j][k] + M[14][4] * fIn4Ghost[i][j][k] + M[14][5] * fIn5Ghost[i][j][k] + M[14][6] * fIn6Ghost[i][j][k] + M[14][7] * fIn7Ghost[i][j][k] + M[14][8] * fIn8Ghost[i][j][k]\
							+ M[14][9] * fIn9Ghost[i][j][k] + M[14][10] * fIn10Ghost[i][j][k] + M[14][11] * fIn11Ghost[i][j][k] + M[14][12] * fIn12Ghost[i][j][k] + M[14][13] * fIn13Ghost[i][j][k] + M[14][14] * fIn14Ghost[i][j][k] + M[14][15] * fIn15Ghost[i][j][k] + M[14][16] * fIn16Ghost[i][j][k] + M[14][17] * fIn17Ghost[i][j][k] + M[14][18] * fIn18Ghost[i][j][k]));
						mfIn15 = (A_g[15][15] * (M[15][0] * fIn0Ghost[i][j][k] + M[15][1] * fIn1Ghost[i][j][k] + M[15][2] * fIn2Ghost[i][j][k] + M[15][3] * fIn3Ghost[i][j][k] + M[15][4] * fIn4Ghost[i][j][k] + M[15][5] * fIn5Ghost[i][j][k] + M[15][6] * fIn6Ghost[i][j][k] + M[15][7] * fIn7Ghost[i][j][k] + M[15][8] * fIn8Ghost[i][j][k]\
							+ M[15][9] * fIn9Ghost[i][j][k] + M[15][10] * fIn10Ghost[i][j][k] + M[15][11] * fIn11Ghost[i][j][k] + M[15][12] * fIn12Ghost[i][j][k] + M[15][13] * fIn13Ghost[i][j][k] + M[15][14] * fIn14Ghost[i][j][k] + M[15][15] * fIn15Ghost[i][j][k] + M[15][16] * fIn16Ghost[i][j][k] + M[15][17] * fIn17Ghost[i][j][k] + M[15][18] * fIn18Ghost[i][j][k]));
						mfIn16 = (A_g[16][16] * (M[16][0] * fIn0Ghost[i][j][k] + M[16][1] * fIn1Ghost[i][j][k] + M[16][2] * fIn2Ghost[i][j][k] + M[16][3] * fIn3Ghost[i][j][k] + M[16][4] * fIn4Ghost[i][j][k] + M[16][5] * fIn5Ghost[i][j][k] + M[16][6] * fIn6Ghost[i][j][k] + M[16][7] * fIn7Ghost[i][j][k] + M[16][8] * fIn8Ghost[i][j][k]\
							+ M[16][9] * fIn9Ghost[i][j][k] + M[16][10] * fIn10Ghost[i][j][k] + M[16][11] * fIn11Ghost[i][j][k] + M[16][12] * fIn12Ghost[i][j][k] + M[16][13] * fIn13Ghost[i][j][k] + M[16][14] * fIn14Ghost[i][j][k] + M[16][15] * fIn15Ghost[i][j][k] + M[16][16] * fIn16Ghost[i][j][k] + M[16][17] * fIn17Ghost[i][j][k] + M[16][18] * fIn18Ghost[i][j][k]));
						mfIn17 = (A_g[17][17] * (M[17][0] * fIn0Ghost[i][j][k] + M[17][1] * fIn1Ghost[i][j][k] + M[17][2] * fIn2Ghost[i][j][k] + M[17][3] * fIn3Ghost[i][j][k] + M[17][4] * fIn4Ghost[i][j][k] + M[17][5] * fIn5Ghost[i][j][k] + M[17][6] * fIn6Ghost[i][j][k] + M[17][7] * fIn7Ghost[i][j][k] + M[17][8] * fIn8Ghost[i][j][k]\
							+ M[17][9] * fIn9Ghost[i][j][k] + M[17][10] * fIn10Ghost[i][j][k] + M[17][11] * fIn11Ghost[i][j][k] + M[17][12] * fIn12Ghost[i][j][k] + M[17][13] * fIn13Ghost[i][j][k] + M[17][14] * fIn14Ghost[i][j][k] + M[17][15] * fIn15Ghost[i][j][k] + M[17][16] * fIn16Ghost[i][j][k] + M[17][17] * fIn17Ghost[i][j][k] + M[17][18] * fIn18Ghost[i][j][k]));
						mfIn18 = (A_g[18][18] * (M[18][0] * fIn0Ghost[i][j][k] + M[18][1] * fIn1Ghost[i][j][k] + M[18][2] * fIn2Ghost[i][j][k] + M[18][3] * fIn3Ghost[i][j][k] + M[18][4] * fIn4Ghost[i][j][k] + M[18][5] * fIn5Ghost[i][j][k] + M[18][6] * fIn6Ghost[i][j][k] + M[18][7] * fIn7Ghost[i][j][k] + M[18][8] * fIn8Ghost[i][j][k]\
							+ M[18][9] * fIn9Ghost[i][j][k] + M[18][10] * fIn10Ghost[i][j][k] + M[18][11] * fIn11Ghost[i][j][k] + M[18][12] * fIn12Ghost[i][j][k] + M[18][13] * fIn13Ghost[i][j][k] + M[18][14] * fIn14Ghost[i][j][k] + M[18][15] * fIn15Ghost[i][j][k] + M[18][16] * fIn16Ghost[i][j][k] + M[18][17] * fIn17Ghost[i][j][k] + M[18][18] * fIn18Ghost[i][j][k]));
						mfEq0 = A_g[0][0] * rhoGhost[i][j][k];
						mfEq1 = A_g[1][1] * (-11.0*rhoGhost[i][j][k] + 19.0*((rhoGhost[i][j][k] * uxGhost[i][j][k])*(rhoGhost[i][j][k] * uxGhost[i][j][k]) + (rhoGhost[i][j][k] * uyGhost[i][j][k])*(rhoGhost[i][j][k] * uyGhost[i][j][k]) + (rhoGhost[i][j][k] * uzGhost[i][j][k])*(rhoGhost[i][j][k] * uzGhost[i][j][k])) / rhoGhost[i][j][k]);
						mfEq2 = A_g[2][2] * (3.0*rhoGhost[i][j][k] - 11.0 / 2.0*((rhoGhost[i][j][k] * uxGhost[i][j][k])*(rhoGhost[i][j][k] * uxGhost[i][j][k]) + (rhoGhost[i][j][k] * uyGhost[i][j][k])*(rhoGhost[i][j][k] * uyGhost[i][j][k]) + (rhoGhost[i][j][k] * uzGhost[i][j][k])*(rhoGhost[i][j][k] * uzGhost[i][j][k])) / rhoGhost[i][j][k]);
						mfEq3 = A_g[3][3] * (rhoGhost[i][j][k] * uxGhost[i][j][k]);
						mfEq4 = A_g[4][4] * (-2.0 / 3.0*rhoGhost[i][j][k] * uxGhost[i][j][k]);
						mfEq5 = A_g[5][5] * (rhoGhost[i][j][k] * uyGhost[i][j][k]);
						mfEq6 = A_g[6][6] * (-2.0 / 3.0*rhoGhost[i][j][k] * uyGhost[i][j][k]);
						mfEq7 = A_g[7][7] * (rhoGhost[i][j][k] * uzGhost[i][j][k]);
						mfEq8 = A_g[8][8] * (-2.0 / 3.0*rhoGhost[i][j][k] * uzGhost[i][j][k]);
						mfEq9 = A_g[9][9] * ((2.0*(rhoGhost[i][j][k] * uxGhost[i][j][k])*(rhoGhost[i][j][k] * uxGhost[i][j][k]) - ((rhoGhost[i][j][k] * uyGhost[i][j][k])*(rhoGhost[i][j][k] * uyGhost[i][j][k]) + (rhoGhost[i][j][k] * uzGhost[i][j][k])*(rhoGhost[i][j][k] * uzGhost[i][j][k]))) / rhoGhost[i][j][k]);
						mfEq10 = A_g[10][10] * (-1.0 / 2.0*((2.0*(rhoGhost[i][j][k] * uxGhost[i][j][k])*(rhoGhost[i][j][k] * uxGhost[i][j][k]) - ((rhoGhost[i][j][k] * uyGhost[i][j][k])*(rhoGhost[i][j][k] * uyGhost[i][j][k]) + (rhoGhost[i][j][k] * uzGhost[i][j][k])*(rhoGhost[i][j][k] * uzGhost[i][j][k]))) / rhoGhost[i][j][k]));
						mfEq11 = A_g[11][11] * (((rhoGhost[i][j][k] * uyGhost[i][j][k])*(rhoGhost[i][j][k] * uyGhost[i][j][k]) - (rhoGhost[i][j][k] * uzGhost[i][j][k])*(rhoGhost[i][j][k] * uzGhost[i][j][k])) / rhoGhost[i][j][k]);
						mfEq12 = A_g[12][12] * (-1.0 / 2.0*((rhoGhost[i][j][k] * uyGhost[i][j][k])*(rhoGhost[i][j][k] * uyGhost[i][j][k]) - (rhoGhost[i][j][k] * uzGhost[i][j][k])*(rhoGhost[i][j][k] * uzGhost[i][j][k])) / rhoGhost[i][j][k]);
						mfEq13 = A_g[13][13] * (rhoGhost[i][j][k] * uxGhost[i][j][k] * rhoGhost[i][j][k] * uyGhost[i][j][k] / rhoGhost[i][j][k]);
						mfEq14 = A_g[14][14] * (rhoGhost[i][j][k] * uyGhost[i][j][k] * rhoGhost[i][j][k] * uzGhost[i][j][k] / rhoGhost[i][j][k]);
						mfEq15 = A_g[15][15] * (rhoGhost[i][j][k] * uxGhost[i][j][k] * rhoGhost[i][j][k] * uzGhost[i][j][k] / rhoGhost[i][j][k]);
						mfEq16 = 0.0;
						mfEq17 = 0.0;
						mfEq18 = 0.0;
						mForce0 = (Factor0*0.0);
						mForce1=Factor1*(38.0*(uzGhost[i][j][k]*bodyForce[i][j][k]));
						mForce2 = Factor2*(-11.0*(uzGhost[i][j][k] * bodyForce[i][j][k]));
						mForce3 = Factor3*0;
						mForce4 = Factor4*(-2.0 / 3.0*0);
						mForce5 = Factor5*0;
						mForce6 = Factor6*(-2.0 / 3.0*0);
						mForce7 = Factor7*bodyForce[i][j][k];
						mForce8 = Factor8*(-2.0 / 3.0*bodyForce[i][j][k]);
						mForce9 = Factor9*(2.0*(2.0*uxGhost[i][j][k] * 0 - uyGhost[i][j][k] * 0 - uzGhost[i][j][k] * bodyForce[i][j][k]));
						mForce10 = Factor10*(-2.0*uxGhost[i][j][k] * 0 + uyGhost[i][j][k] * 0 + uzGhost[i][j][k] * bodyForce[i][j][k]);
						mForce11 = Factor11*(2.0*(uyGhost[i][j][k] * 0 - uzGhost[i][j][k] * bodyForce[i][j][k]));
						mForce12 = Factor12*(-uyGhost[i][j][k] * 0 + uzGhost[i][j][k] * bodyForce[i][j][k]);
						mForce13 = Factor13*(uyGhost[i][j][k] * 0 + uxGhost[i][j][k] * 0);
						mForce14 = Factor14*(uzGhost[i][j][k] * 0 + uyGhost[i][j][k] * bodyForce[i][j][k]);
						mForce15 = Factor15*(uzGhost[i][j][k] * 0 + uxGhost[i][j][k] * bodyForce[i][j][k]);
						mForce16 = Factor16*0.0;
						mForce17 = Factor17*0.0;
						mForce18 = Factor18*0.0;
						Out0[i][j][k] = fIn0Ghost[i][j][k] - (vfGhostIn0 - vfGhostEq0) + vGhostForce0;
						Out1[i][j][k] = fIn1Ghost[i][j][k] - (vfGhostIn1 - vfGhostEq1) + vGhostForce1;
						Out2[i][j][k] = fIn2Ghost[i][j][k] - (vfGhostIn2 - vfGhostEq2) + vGhostForce2;
						Out3[i][j][k] = fIn3Ghost[i][j][k] - (vfGhostIn3 - vfGhostEq3) + vGhostForce3;
						Out4[i][j][k] = fIn4Ghost[i][j][k] - (vfGhostIn4 - vfGhostEq4) + vGhostForce4;
						Out5[i][j][k] = fIn5Ghost[i][j][k] - (vfGhostIn5 - vfGhostEq5) + vGhostForce5;
						Out6[i][j][k] = fIn6Ghost[i][j][k] - (vfGhostIn6 - vfGhostEq6) + vGhostForce6;
						Out7[i][j][k] = fIn7Ghost[i][j][k] - (vfGhostIn7 - vfGhostEq7) + vGhostForce7;
						Out8[i][j][k] = fIn8Ghost[i][j][k] - (vfGhostIn8 - vfGhostEq8) + vGhostForce8;
						Out9[i][j][k] = fIn9Ghost[i][j][k] - (vfGhostIn9 - vfGhostEq9) + vGhostForce9;
						Out10[i][j][k] = fIn10Ghost[i][j][k] - (vfGhostIn10 - vfGhostEq10) + vGhostForce10;
						Out11[i][j][k] = fIn11Ghost[i][j][k] - (vfGhostIn11 - vfGhostEq11) + vGhostForce11;
						Out12[i][j][k] = fIn12Ghost[i][j][k] - (vfGhostIn12 - vfGhostEq12) + vGhostForce12;
						Out13[i][j][k] = fIn13Ghost[i][j][k] - (vfGhostIn13 - vfGhostEq13) + vGhostForce13;
						Out14[i][j][k] = fIn14Ghost[i][j][k] - (vfGhostIn14 - vfGhostEq14) + vGhostForce14;
						Out15[i][j][k] = fIn15Ghost[i][j][k] - (vfGhostIn15 - vfGhostEq15) + vGhostForce15;
						Out16[i][j][k] = fIn16Ghost[i][j][k] - (vfGhostIn16 - vfGhostEq16) + vGhostForce16;
						Out17[i][j][k] = fIn17Ghost[i][j][k] - (vfGhostIn17 - vfGhostEq17) + vGhostForce17;
						Out18[i][j][k] = fIn18Ghost[i][j][k] - (vfGhostIn18 - vfGhostEq18) + vGhostForce18;
					}
					else{
						Out0[i][j][k] = 0.0;
						Out1[i][j][k] = 0.0;
						Out2[i][j][k] = 0.0;
						Out3[i][j][k] = 0.0;
						Out4[i][j][k] = 0.0;
						Out5[i][j][k] = 0.0;
						Out6[i][j][k] = 0.0;
						Out7[i][j][k] = 0.0;
						Out8[i][j][k] = 0.0;
						Out9[i][j][k] = 0.0;
						Out10[i][j][k] = 0.0;
						Out11[i][j][k] = 0.0;
						Out12[i][j][k] = 0.0;
						Out13[i][j][k] = 0.0;
						Out14[i][j][k] = 0.0;
						Out15[i][j][k] = 0.0;
						Out16[i][j][k] = 0.0;
						Out17[i][j][k] = 0.0;
						Out18[i][j][k] = 0.0;
					}
				}
			}
		}
	}
}

void MemoryFreeGhostVelocity(int m, int n, int q){
	freememory(rhoGhost,m,n,q);
	freememory(uxGhost,m,n,q);
	freememory(uyGhost,m,n,q);
	freememory(uzGhost,m,n,q);
	freememory(fIn0Ghost,m,n,q);
	freememory(fIn1Ghost,m,n,q);
	freememory(fIn2Ghost,m,n,q);
	freememory(fIn3Ghost,m,n,q);
	freememory(fIn4Ghost,m,n,q);
	freememory(fIn5Ghost,m,n,q);
	freememory(fIn6Ghost,m,n,q);
	freememory(fIn7Ghost,m,n,q);
	freememory(fIn8Ghost,m,n,q);
	freememory(fIn9Ghost,m,n,q);
	freememory(fIn10Ghost,m,n,q);
	freememory(fIn11Ghost,m,n,q);
	freememory(fIn12Ghost,m,n,q);
	freememory(fIn13Ghost,m,n,q);
	freememory(fIn14Ghost,m,n,q);
	freememory(fIn15Ghost,m,n,q);
	freememory(fIn16Ghost,m,n,q);
	freememory(fIn17Ghost,m,n,q);
	freememory(fIn18Ghost,m,n,q);
}

void StreamGhostVelocity(int m, int n, int q){
	shift0(Out0,fIn0Ghost,m,n,q);
	shift1(Out1,fIn1Ghost,m,n,q);
	shift2(Out2,fIn2Ghost,m,n,q);
	shift3(Out3,fIn3Ghost,m,n,q);
	shift4(Out4,fIn4Ghost,m,n,q);
	shift5(Out5,fIn5Ghost,m,n,q);
	shift6(Out6,fIn6Ghost,m,n,q);
	shift7(Out7,fIn7Ghost,m,n,q);
	shift8(Out8,fIn8Ghost,m,n,q);
	shift9(Out9,fIn9Ghost,m,n,q);
	shift10(Out10,fIn10Ghost,m,n,q);
	shift11(Out11,fIn11Ghost,m,n,q);
	shift12(Out12,fIn12Ghost,m,n,q);
	shift13(Out13,fIn13Ghost,m,n,q);
	shift14(Out14,fIn14Ghost,m,n,q);
	shift15(Out15,fIn15Ghost,m,n,q);
	shift16(Out16,fIn16Ghost,m,n,q);
	shift17(Out17,fIn17Ghost,m,n,q);
	shift18(Out18,fIn18Ghost,m,n,q);
}

void FieldCalculationGhostVelocity(int m, int n, int q, double gasCritical){
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
						rhoGhost[i][j][k]=fIn0Ghost[i][j][k]+fIn1Ghost[i][j][k]+fIn2Ghost[i][j][k]+fIn3Ghost[i][j][k]+fIn4Ghost[i][j][k]+fIn5Ghost[i][j][k]+fIn6Ghost[i][j][k]+fIn7Ghost[i][j][k]+fIn8Ghost[i][j][k]\
							+fIn9Ghost[i][j][k]+fIn10Ghost[i][j][k]+fIn11Ghost[i][j][k]+fIn12Ghost[i][j][k]+fIn13Ghost[i][j][k]+fIn14Ghost[i][j][k]+fIn15Ghost[i][j][k]+fIn16Ghost[i][j][k]+fIn17Ghost[i][j][k]+fIn18Ghost[i][j][k];
						uxGhost[i][j][k]=uTotXGhost;
						uyGhost[i][j][k]=uTotYGhost;
						uzGhost[i][j][k]=uTotZGhost+bodyForce[i][j][k]/2.0/rhoGhost[i][j][k];
					}
					if (domain[i][j][k]!=0||rho[i][j][k]<=gasCritical){
						rhoGhost[i][j][k]=0.0;
						uxGhost[i][j][k]=0.0;
						uyGhost[i][j][k]=0.0;
						uzGhost[i][j][k]=0.0;
						fIn0Ghost[i][j][k]=0.0;
						fIn1Ghost[i][j][k]=0.0;
						fIn2Ghost[i][j][k]=0.0;
						fIn3Ghost[i][j][k]=0.0;
						fIn4Ghost[i][j][k]=0.0;
						fIn5Ghost[i][j][k]=0.0;
						fIn6Ghost[i][j][k]=0.0;
						fIn7Ghost[i][j][k]=0.0;
						fIn8Ghost[i][j][k]=0.0;
						fIn9Ghost[i][j][k]=0.0;
						fIn10Ghost[i][j][k]=0.0;
						fIn11Ghost[i][j][k]=0.0;
						fIn12Ghost[i][j][k]=0.0;
						fIn13Ghost[i][j][k]=0.0;
						fIn14Ghost[i][j][k]=0.0;
						fIn15Ghost[i][j][k]=0.0;
						fIn16Ghost[i][j][k]=0.0;
						fIn17Ghost[i][j][k]=0.0;
						fIn18Ghost[i][j][k]=0.0;
					}
				}
			}
		}
	}
}

void refillGhostVelocity(int m, int n, int q, double gasCritical,double ***rho1, double ***rho2, double gasDensity){
	int i;
	int j;
	int k;
	int number;
	double valueRho;
	double valueRho1;
#pragma omp parallel private(i,j,k,number,valueRho)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					number=0;
					if (domain[i][j][k]==0&&rho1[i][j][k]>=gasCritical&&rho2[i][j][k]<gasCritical){	
							fIn0Ghost[i][j][k]=0.0;
							fIn1Ghost[i][j][k]=0.0;
							fIn2Ghost[i][j][k]=0.0;
							fIn3Ghost[i][j][k]=0.0;
							fIn4Ghost[i][j][k]=0.0;
							fIn5Ghost[i][j][k]=0.0;
							fIn6Ghost[i][j][k]=0.0;
							fIn7Ghost[i][j][k]=0.0;
							fIn8Ghost[i][j][k]=0.0;
							fIn9Ghost[i][j][k]=0.0;
							fIn10Ghost[i][j][k]=0.0;
							fIn11Ghost[i][j][k]=0.0;
							fIn12Ghost[i][j][k]=0.0;
							fIn13Ghost[i][j][k]=0.0;
							fIn14Ghost[i][j][k]=0.0;
							fIn15Ghost[i][j][k]=0.0;
							fIn16Ghost[i][j][k]=0.0;
							fIn17Ghost[i][j][k]=0.0;
							fIn18Ghost[i][j][k]=0.0;
						valueRho=neightbour4(i, j, k, m, rho1);
						valueRho1=neightbour4(i, j, k, m, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]+neightbour4(i,j,k,m,fIn0Ghost);
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]+neightbour4(i,j,k,m,fIn1Ghost);
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]+neightbour4(i,j,k,m,fIn2Ghost);
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]+neightbour4(i,j,k,m,fIn3Ghost);
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]+neightbour4(i,j,k,m,fIn4Ghost);
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]+neightbour4(i,j,k,m,fIn5Ghost);
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]+neightbour4(i,j,k,m,fIn6Ghost);
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]+neightbour4(i,j,k,m,fIn7Ghost);
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]+neightbour4(i,j,k,m,fIn8Ghost);
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]+neightbour4(i,j,k,m,fIn9Ghost);
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]+neightbour4(i,j,k,m,fIn10Ghost);
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]+neightbour4(i,j,k,m,fIn11Ghost);
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]+neightbour4(i,j,k,m,fIn12Ghost);
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]+neightbour4(i,j,k,m,fIn13Ghost);
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]+neightbour4(i,j,k,m,fIn14Ghost);
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]+neightbour4(i,j,k,m,fIn15Ghost);
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]+neightbour4(i,j,k,m,fIn16Ghost);
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]+neightbour4(i,j,k,m,fIn17Ghost);
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]+neightbour4(i,j,k,m,fIn18Ghost);
						}
						valueRho=neightbour5(i, j, k, n, rho1);
						valueRho1=neightbour5(i, j, k, n, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){					
							number=number+1;
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]+neightbour5(i,j,k,n,fIn0Ghost);
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]+neightbour5(i,j,k,n,fIn1Ghost);
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]+neightbour5(i,j,k,n,fIn2Ghost);
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]+neightbour5(i,j,k,n,fIn3Ghost);
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]+neightbour5(i,j,k,n,fIn4Ghost);
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]+neightbour5(i,j,k,n,fIn5Ghost);
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]+neightbour5(i,j,k,n,fIn6Ghost);
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]+neightbour5(i,j,k,n,fIn7Ghost);
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]+neightbour5(i,j,k,n,fIn8Ghost);
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]+neightbour5(i,j,k,n,fIn9Ghost);
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]+neightbour5(i,j,k,n,fIn10Ghost);
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]+neightbour5(i,j,k,n,fIn11Ghost);
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]+neightbour5(i,j,k,n,fIn12Ghost);
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]+neightbour5(i,j,k,n,fIn13Ghost);
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]+neightbour5(i,j,k,n,fIn14Ghost);
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]+neightbour5(i,j,k,n,fIn15Ghost);
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]+neightbour5(i,j,k,n,fIn16Ghost);
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]+neightbour5(i,j,k,n,fIn17Ghost);
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]+neightbour5(i,j,k,n,fIn18Ghost);
						}
						valueRho=neightbour6(i, j, k, q, rho1);
						valueRho1=neightbour6(i, j, k, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]+neightbour6(i,j,k,q,fIn0Ghost);
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]+neightbour6(i,j,k,q,fIn1Ghost);
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]+neightbour6(i,j,k,q,fIn2Ghost);
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]+neightbour6(i,j,k,q,fIn3Ghost);
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]+neightbour6(i,j,k,q,fIn4Ghost);
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]+neightbour6(i,j,k,q,fIn5Ghost);
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]+neightbour6(i,j,k,q,fIn6Ghost);
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]+neightbour6(i,j,k,q,fIn7Ghost);
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]+neightbour6(i,j,k,q,fIn8Ghost);
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]+neightbour6(i,j,k,q,fIn9Ghost);
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]+neightbour6(i,j,k,q,fIn10Ghost);
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]+neightbour6(i,j,k,q,fIn11Ghost);
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]+neightbour6(i,j,k,q,fIn12Ghost);
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]+neightbour6(i,j,k,q,fIn13Ghost);
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]+neightbour6(i,j,k,q,fIn14Ghost);
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]+neightbour6(i,j,k,q,fIn15Ghost);
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]+neightbour6(i,j,k,q,fIn16Ghost);
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]+neightbour6(i,j,k,q,fIn17Ghost);
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]+neightbour6(i,j,k,q,fIn18Ghost);							
						}							
						valueRho=neightbour1(i, j, k, m, rho1);
						valueRho1=neightbour1(i, j, k, m, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]+neightbour1(i,j,k,m,fIn0Ghost);
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]+neightbour1(i,j,k,m,fIn1Ghost);
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]+neightbour1(i,j,k,m,fIn2Ghost);
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]+neightbour1(i,j,k,m,fIn3Ghost);
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]+neightbour1(i,j,k,m,fIn4Ghost);
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]+neightbour1(i,j,k,m,fIn5Ghost);
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]+neightbour1(i,j,k,m,fIn6Ghost);
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]+neightbour1(i,j,k,m,fIn7Ghost);
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]+neightbour1(i,j,k,m,fIn8Ghost);
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]+neightbour1(i,j,k,m,fIn9Ghost);
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]+neightbour1(i,j,k,m,fIn10Ghost);
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]+neightbour1(i,j,k,m,fIn11Ghost);
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]+neightbour1(i,j,k,m,fIn12Ghost);
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]+neightbour1(i,j,k,m,fIn13Ghost);
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]+neightbour1(i,j,k,m,fIn14Ghost);
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]+neightbour1(i,j,k,m,fIn15Ghost);
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]+neightbour1(i,j,k,m,fIn16Ghost);
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]+neightbour1(i,j,k,m,fIn17Ghost);
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]+neightbour1(i,j,k,m,fIn18Ghost);
						}
						valueRho=neightbour2(i, j, k, n, rho1);
						valueRho1=neightbour2(i, j, k, n, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]+neightbour2(i,j,k,n,fIn0Ghost);
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]+neightbour2(i,j,k,n,fIn1Ghost);
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]+neightbour2(i,j,k,n,fIn2Ghost);
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]+neightbour2(i,j,k,n,fIn3Ghost);
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]+neightbour2(i,j,k,n,fIn4Ghost);
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]+neightbour2(i,j,k,n,fIn5Ghost);
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]+neightbour2(i,j,k,n,fIn6Ghost);
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]+neightbour2(i,j,k,n,fIn7Ghost);
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]+neightbour2(i,j,k,n,fIn8Ghost);
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]+neightbour2(i,j,k,n,fIn9Ghost);
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]+neightbour2(i,j,k,n,fIn10Ghost);
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]+neightbour2(i,j,k,n,fIn11Ghost);
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]+neightbour2(i,j,k,n,fIn12Ghost);
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]+neightbour2(i,j,k,n,fIn13Ghost);
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]+neightbour2(i,j,k,n,fIn14Ghost);
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]+neightbour2(i,j,k,n,fIn15Ghost);
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]+neightbour2(i,j,k,n,fIn16Ghost);
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]+neightbour2(i,j,k,n,fIn17Ghost);
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]+neightbour2(i,j,k,n,fIn18Ghost);
						}
						valueRho=neightbour3(i, j, k, q, rho1);
						valueRho1=neightbour3(i, j, k, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]+neightbour3(i,j,k,q,fIn0Ghost);
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]+neightbour3(i,j,k,q,fIn1Ghost);
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]+neightbour3(i,j,k,q,fIn2Ghost);
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]+neightbour3(i,j,k,q,fIn3Ghost);
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]+neightbour3(i,j,k,q,fIn4Ghost);
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]+neightbour3(i,j,k,q,fIn5Ghost);
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]+neightbour3(i,j,k,q,fIn6Ghost);
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]+neightbour3(i,j,k,q,fIn7Ghost);
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]+neightbour3(i,j,k,q,fIn8Ghost);
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]+neightbour3(i,j,k,q,fIn9Ghost);
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]+neightbour3(i,j,k,q,fIn10Ghost);
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]+neightbour3(i,j,k,q,fIn11Ghost);
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]+neightbour3(i,j,k,q,fIn12Ghost);
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]+neightbour3(i,j,k,q,fIn13Ghost);
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]+neightbour3(i,j,k,q,fIn14Ghost);
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]+neightbour3(i,j,k,q,fIn15Ghost);
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]+neightbour3(i,j,k,q,fIn16Ghost);
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]+neightbour3(i,j,k,q,fIn17Ghost);
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]+neightbour3(i,j,k,q,fIn18Ghost);
						}
						valueRho=neightbour10(i, j, k, m, n, rho1);
						valueRho1=neightbour10(i, j, k, m, n, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn0Ghost);
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn1Ghost);
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn2Ghost);
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn3Ghost);
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn4Ghost);
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn5Ghost);
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn6Ghost);
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn7Ghost);
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn8Ghost);
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn9Ghost);
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn10Ghost);
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn11Ghost);
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn12Ghost);
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn13Ghost);
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn14Ghost);
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn15Ghost);
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn16Ghost);
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn17Ghost);
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]+neightbour10(i,j,k,m,n,fIn18Ghost);
						}
						valueRho=neightbour11(i, j, k, m, q, rho1);
						valueRho1=neightbour11(i, j, k, m, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn0Ghost);
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn1Ghost);
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn2Ghost);
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn3Ghost);
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn4Ghost);
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn5Ghost);
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn6Ghost);
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn7Ghost);
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn8Ghost);
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn9Ghost);
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn10Ghost);
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn11Ghost);
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn12Ghost);
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn13Ghost);
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn14Ghost);
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn15Ghost);
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn16Ghost);
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn17Ghost);
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]+neightbour11(i,j,k,m,q,fIn18Ghost);
						}
						valueRho=neightbour12(i, j, k, n, q, rho1);
						valueRho1=neightbour12(i, j, k, n, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn0Ghost);
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn1Ghost);
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn2Ghost);
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn3Ghost);
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn4Ghost);
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn5Ghost);
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn6Ghost);
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn7Ghost);
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn8Ghost);
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn9Ghost);
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn10Ghost);
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn11Ghost);
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn12Ghost);
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn13Ghost);
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn14Ghost);
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn15Ghost);
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn16Ghost);
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn17Ghost);
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]+neightbour12(i,j,k,n,q,fIn18Ghost);
						}
						valueRho=neightbour7(i, j, k, m, n, rho1);
						valueRho1=neightbour7(i, j, k, m, n, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn0Ghost);
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn1Ghost);
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn2Ghost);
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn3Ghost);
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn4Ghost);
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn5Ghost);
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn6Ghost);
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn7Ghost);
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn8Ghost);
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn9Ghost);
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn10Ghost);
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn11Ghost);
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn12Ghost);
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn13Ghost);
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn14Ghost);
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn15Ghost);
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn16Ghost);
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn17Ghost);
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]+neightbour7(i,j,k,m,n,fIn18Ghost);
						}
						valueRho=neightbour8(i, j, k, m, q, rho1);
						valueRho1=neightbour8(i, j, k, m, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn0Ghost);
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn1Ghost);
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn2Ghost);
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn3Ghost);
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn4Ghost);
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn5Ghost);
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn6Ghost);
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn7Ghost);
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn8Ghost);
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn9Ghost);
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn10Ghost);
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn11Ghost);
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn12Ghost);
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn13Ghost);
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn14Ghost);
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn15Ghost);
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn16Ghost);
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn17Ghost);
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]+neightbour8(i,j,k,m,q,fIn18Ghost);
						}
						valueRho=neightbour9(i, j, k, n, q, rho1);
						valueRho1=neightbour9(i, j, k, n, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn0Ghost);
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn1Ghost);
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn2Ghost);
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn3Ghost);
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn4Ghost);
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn5Ghost);
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn6Ghost);
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn7Ghost);
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn8Ghost);
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn9Ghost);
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn10Ghost);
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn11Ghost);
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn12Ghost);
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn13Ghost);
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn14Ghost);
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn15Ghost);
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn16Ghost);
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn17Ghost);
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]+neightbour9(i,j,k,n,q,fIn18Ghost);
						}
						valueRho=neightbour16(i, j, k, m, n, rho1);
						valueRho1=neightbour16(i, j, k, m, n, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn0Ghost);
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn1Ghost);
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn2Ghost);
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn3Ghost);
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn4Ghost);
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn5Ghost);
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn6Ghost);
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn7Ghost);
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn8Ghost);
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn9Ghost);
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn10Ghost);
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn11Ghost);
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn12Ghost);
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn13Ghost);
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn14Ghost);
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn15Ghost);
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn16Ghost);
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn17Ghost);
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]+neightbour16(i,j,k,m,n,fIn18Ghost);
						}
						valueRho=neightbour17(i, j, k, m, q, rho1);
						valueRho1=neightbour17(i, j, k, m, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn0Ghost);
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn1Ghost);
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn2Ghost);
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn3Ghost);
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn4Ghost);
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn5Ghost);
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn6Ghost);
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn7Ghost);
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn8Ghost);
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn9Ghost);
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn10Ghost);
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn11Ghost);
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn12Ghost);
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn13Ghost);
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn14Ghost);
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn15Ghost);
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn16Ghost);
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn17Ghost);
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]+neightbour17(i,j,k,m,q,fIn18Ghost);
						}
						valueRho=neightbour18(i, j, k, n, q, rho1);
						valueRho1=neightbour18(i, j, k, n, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn0Ghost);
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn1Ghost);
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn2Ghost);
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn3Ghost);
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn4Ghost);
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn5Ghost);
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn6Ghost);
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn7Ghost);
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn8Ghost);
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn9Ghost);
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn10Ghost);
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn11Ghost);
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn12Ghost);
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn13Ghost);
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn14Ghost);
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn15Ghost);
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn16Ghost);
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn17Ghost);
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]+neightbour18(i,j,k,n,q,fIn18Ghost);
						}
						valueRho=neightbour13(i, j, k, m, n, rho1);
						valueRho1=neightbour13(i, j, k, m, n, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn0Ghost);
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn1Ghost);
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn2Ghost);
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn3Ghost);
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn4Ghost);
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn5Ghost);
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn6Ghost);
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn7Ghost);
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn8Ghost);
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn9Ghost);
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn10Ghost);
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn11Ghost);
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn12Ghost);
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn13Ghost);
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn14Ghost);
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn15Ghost);
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn16Ghost);
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn17Ghost);
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]+neightbour13(i,j,k,m,n,fIn18Ghost);
						}
						valueRho=neightbour14(i, j, k, m, q, rho1);
						valueRho1=neightbour14(i, j, k, m, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn0Ghost);
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn1Ghost);
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn2Ghost);
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn3Ghost);
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn4Ghost);
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn5Ghost);
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn6Ghost);
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn7Ghost);
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn8Ghost);
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn9Ghost);
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn10Ghost);
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn11Ghost);
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn12Ghost);
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn13Ghost);
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn14Ghost);
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn15Ghost);
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn16Ghost);
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn17Ghost);
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]+neightbour14(i,j,k,m,q,fIn18Ghost);
						}
						valueRho=neightbour15(i, j, k, n, q, rho1);
						valueRho1=neightbour15(i, j, k, n, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn0Ghost);
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn1Ghost);
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn2Ghost);
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn3Ghost);
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn4Ghost);
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn5Ghost);
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn6Ghost);
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn7Ghost);
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn8Ghost);
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn9Ghost);
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn10Ghost);
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn11Ghost);
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn12Ghost);
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn13Ghost);
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn14Ghost);
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn15Ghost);
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn16Ghost);
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn17Ghost);
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]+neightbour15(i,j,k,n,q,fIn18Ghost);
						}						
						if (number==0){
							fIn0Ghost[i][j][k]=gasDensity*tNS[0];
							fIn1Ghost[i][j][k]=gasDensity*tNS[1];
							fIn2Ghost[i][j][k]=gasDensity*tNS[2];
							fIn3Ghost[i][j][k]=gasDensity*tNS[3];
							fIn4Ghost[i][j][k]=gasDensity*tNS[4];
							fIn5Ghost[i][j][k]=gasDensity*tNS[5];
							fIn6Ghost[i][j][k]=gasDensity*tNS[6];
							fIn7Ghost[i][j][k]=gasDensity*tNS[7];
							fIn8Ghost[i][j][k]=gasDensity*tNS[8];
							fIn9Ghost[i][j][k]=gasDensity*tNS[9];
							fIn10Ghost[i][j][k]=gasDensity*tNS[10];
							fIn11Ghost[i][j][k]=gasDensity*tNS[11];
							fIn12Ghost[i][j][k]=gasDensity*tNS[12];
							fIn13Ghost[i][j][k]=gasDensity*tNS[13];
							fIn14Ghost[i][j][k]=gasDensity*tNS[14];
							fIn15Ghost[i][j][k]=gasDensity*tNS[15];
							fIn16Ghost[i][j][k]=gasDensity*tNS[16];
							fIn17Ghost[i][j][k]=gasDensity*tNS[17];
							fIn18Ghost[i][j][k]=gasDensity*tNS[18];
							rhoGhost[i][j][k]=gasDensity;
							uxGhost[i][j][k]=0.0;
							uyGhost[i][j][k]=0.0;
							uzGhost[i][j][k]=0.0;
						}
						else{
							fIn0Ghost[i][j][k]=fIn0Ghost[i][j][k]/number;
							fIn1Ghost[i][j][k]=fIn1Ghost[i][j][k]/number;
							fIn2Ghost[i][j][k]=fIn2Ghost[i][j][k]/number;
							fIn3Ghost[i][j][k]=fIn3Ghost[i][j][k]/number;
							fIn4Ghost[i][j][k]=fIn4Ghost[i][j][k]/number;
							fIn5Ghost[i][j][k]=fIn5Ghost[i][j][k]/number;
							fIn6Ghost[i][j][k]=fIn6Ghost[i][j][k]/number;
							fIn7Ghost[i][j][k]=fIn7Ghost[i][j][k]/number;
							fIn8Ghost[i][j][k]=fIn8Ghost[i][j][k]/number;
							fIn9Ghost[i][j][k]=fIn9Ghost[i][j][k]/number;
							fIn10Ghost[i][j][k]=fIn10Ghost[i][j][k]/number;
							fIn11Ghost[i][j][k]=fIn11Ghost[i][j][k]/number;
							fIn12Ghost[i][j][k]=fIn12Ghost[i][j][k]/number;
							fIn13Ghost[i][j][k]=fIn13Ghost[i][j][k]/number;
							fIn14Ghost[i][j][k]=fIn14Ghost[i][j][k]/number;
							fIn15Ghost[i][j][k]=fIn15Ghost[i][j][k]/number;
							fIn16Ghost[i][j][k]=fIn16Ghost[i][j][k]/number;
							fIn17Ghost[i][j][k]=fIn17Ghost[i][j][k]/number;
							fIn18Ghost[i][j][k]=fIn18Ghost[i][j][k]/number;
							rhoGhost[i][j][k]=fIn0Ghost[i][j][k]+fIn1Ghost[i][j][k]+fIn2Ghost[i][j][k]+fIn3Ghost[i][j][k]+fIn4Ghost[i][j][k]+fIn5Ghost[i][j][k]+fIn6Ghost[i][j][k]+fIn7Ghost[i][j][k]+fIn8Ghost[i][j][k]\
								+fIn9Ghost[i][j][k]+fIn10Ghost[i][j][k]+fIn11Ghost[i][j][k]+fIn12Ghost[i][j][k]+fIn13Ghost[i][j][k]+fIn14Ghost[i][j][k]+fIn15Ghost[i][j][k]+fIn16Ghost[i][j][k]+fIn17Ghost[i][j][k]+fIn18Ghost[i][j][k];
							uxGhost[i][j][k]=uTotXGhost;
							uyGhost[i][j][k]=uTotYGhost;
							uzGhost[i][j][k]=uTotZGhost+bodyForce[i][j][k]/2.0/rhoGhost[i][j][k];
						}
					}
				}
			}
		}
	}
}