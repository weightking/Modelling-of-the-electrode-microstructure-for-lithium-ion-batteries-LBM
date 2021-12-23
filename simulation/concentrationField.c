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
#include "matrixMove.h"
#include "electrodePotential.h"
#include "electrolytePotential.h"

#define vfConcentrationIn0 (MInverse_Concentration[0][0]*mfIn0+MInverse_Concentration[0][1]*mfIn1+MInverse_Concentration[0][2]*mfIn2+MInverse_Concentration[0][3]*mfIn3+MInverse_Concentration[0][4]*mfIn4+MInverse_Concentration[0][5]*mfIn5+MInverse_Concentration[0][6]*mfIn6)
#define vfConcentrationIn1 (MInverse_Concentration[1][0]*mfIn0+MInverse_Concentration[1][1]*mfIn1+MInverse_Concentration[1][2]*mfIn2+MInverse_Concentration[1][3]*mfIn3+MInverse_Concentration[1][4]*mfIn4+MInverse_Concentration[1][5]*mfIn5+MInverse_Concentration[1][6]*mfIn6)
#define vfConcentrationIn2 (MInverse_Concentration[2][0]*mfIn0+MInverse_Concentration[2][1]*mfIn1+MInverse_Concentration[2][2]*mfIn2+MInverse_Concentration[2][3]*mfIn3+MInverse_Concentration[2][4]*mfIn4+MInverse_Concentration[2][5]*mfIn5+MInverse_Concentration[2][6]*mfIn6)
#define vfConcentrationIn3 (MInverse_Concentration[3][0]*mfIn0+MInverse_Concentration[3][1]*mfIn1+MInverse_Concentration[3][2]*mfIn2+MInverse_Concentration[3][3]*mfIn3+MInverse_Concentration[3][4]*mfIn4+MInverse_Concentration[3][5]*mfIn5+MInverse_Concentration[3][6]*mfIn6)
#define vfConcentrationIn4 (MInverse_Concentration[4][0]*mfIn0+MInverse_Concentration[4][1]*mfIn1+MInverse_Concentration[4][2]*mfIn2+MInverse_Concentration[4][3]*mfIn3+MInverse_Concentration[4][4]*mfIn4+MInverse_Concentration[4][5]*mfIn5+MInverse_Concentration[4][6]*mfIn6)
#define vfConcentrationIn5 (MInverse_Concentration[5][0]*mfIn0+MInverse_Concentration[5][1]*mfIn1+MInverse_Concentration[5][2]*mfIn2+MInverse_Concentration[5][3]*mfIn3+MInverse_Concentration[5][4]*mfIn4+MInverse_Concentration[5][5]*mfIn5+MInverse_Concentration[5][6]*mfIn6)
#define vfConcentrationIn6 (MInverse_Concentration[6][0]*mfIn0+MInverse_Concentration[6][1]*mfIn1+MInverse_Concentration[6][2]*mfIn2+MInverse_Concentration[6][3]*mfIn3+MInverse_Concentration[6][4]*mfIn4+MInverse_Concentration[6][5]*mfIn5+MInverse_Concentration[6][6]*mfIn6)

#define vfConcentrationEq0 (MInverse_Concentration[0][0]*mfEq0+MInverse_Concentration[0][1]*mfEq1+MInverse_Concentration[0][2]*mfEq2+MInverse_Concentration[0][3]*mfEq3+MInverse_Concentration[0][4]*mfEq4+MInverse_Concentration[0][5]*mfEq5+MInverse_Concentration[0][6]*mfEq6)
#define vfConcentrationEq1 (MInverse_Concentration[1][0]*mfEq0+MInverse_Concentration[1][1]*mfEq1+MInverse_Concentration[1][2]*mfEq2+MInverse_Concentration[1][3]*mfEq3+MInverse_Concentration[1][4]*mfEq4+MInverse_Concentration[1][5]*mfEq5+MInverse_Concentration[1][6]*mfEq6)
#define vfConcentrationEq2 (MInverse_Concentration[2][0]*mfEq0+MInverse_Concentration[2][1]*mfEq1+MInverse_Concentration[2][2]*mfEq2+MInverse_Concentration[2][3]*mfEq3+MInverse_Concentration[2][4]*mfEq4+MInverse_Concentration[2][5]*mfEq5+MInverse_Concentration[2][6]*mfEq6)
#define vfConcentrationEq3 (MInverse_Concentration[3][0]*mfEq0+MInverse_Concentration[3][1]*mfEq1+MInverse_Concentration[3][2]*mfEq2+MInverse_Concentration[3][3]*mfEq3+MInverse_Concentration[3][4]*mfEq4+MInverse_Concentration[3][5]*mfEq5+MInverse_Concentration[3][6]*mfEq6)
#define vfConcentrationEq4 (MInverse_Concentration[4][0]*mfEq0+MInverse_Concentration[4][1]*mfEq1+MInverse_Concentration[4][2]*mfEq2+MInverse_Concentration[4][3]*mfEq3+MInverse_Concentration[4][4]*mfEq4+MInverse_Concentration[4][5]*mfEq5+MInverse_Concentration[4][6]*mfEq6)
#define vfConcentrationEq5 (MInverse_Concentration[5][0]*mfEq0+MInverse_Concentration[5][1]*mfEq1+MInverse_Concentration[5][2]*mfEq2+MInverse_Concentration[5][3]*mfEq3+MInverse_Concentration[5][4]*mfEq4+MInverse_Concentration[5][5]*mfEq5+MInverse_Concentration[5][6]*mfEq6)
#define vfConcentrationEq6 (MInverse_Concentration[6][0]*mfEq0+MInverse_Concentration[6][1]*mfEq1+MInverse_Concentration[6][2]*mfEq2+MInverse_Concentration[6][3]*mfEq3+MInverse_Concentration[6][4]*mfEq4+MInverse_Concentration[6][5]*mfEq5+MInverse_Concentration[6][6]*mfEq6)

#define vfConcentrationD3Q15In0 (MInverse_ConcentrationD3Q15[0][0]*mfIn0+MInverse_ConcentrationD3Q15[0][1]*mfIn1+MInverse_ConcentrationD3Q15[0][2]*mfIn2+MInverse_ConcentrationD3Q15[0][3]*mfIn3+MInverse_ConcentrationD3Q15[0][4]*mfIn4+MInverse_ConcentrationD3Q15[0][5]*mfIn5+MInverse_ConcentrationD3Q15[0][6]*mfIn6+MInverse_ConcentrationD3Q15[0][7]*mfIn7+MInverse_ConcentrationD3Q15[0][8]*mfIn8+MInverse_ConcentrationD3Q15[0][9]*mfIn9+MInverse_ConcentrationD3Q15[0][10]*mfIn10+MInverse_ConcentrationD3Q15[0][11]*mfIn11+MInverse_ConcentrationD3Q15[0][12]*mfIn12+MInverse_ConcentrationD3Q15[0][13]*mfIn13+MInverse_ConcentrationD3Q15[0][14]*mfIn14)
#define vfConcentrationD3Q15In1 (MInverse_ConcentrationD3Q15[1][0]*mfIn0+MInverse_ConcentrationD3Q15[1][1]*mfIn1+MInverse_ConcentrationD3Q15[1][2]*mfIn2+MInverse_ConcentrationD3Q15[1][3]*mfIn3+MInverse_ConcentrationD3Q15[1][4]*mfIn4+MInverse_ConcentrationD3Q15[1][5]*mfIn5+MInverse_ConcentrationD3Q15[1][6]*mfIn6+MInverse_ConcentrationD3Q15[1][7]*mfIn7+MInverse_ConcentrationD3Q15[1][8]*mfIn8+MInverse_ConcentrationD3Q15[1][9]*mfIn9+MInverse_ConcentrationD3Q15[1][10]*mfIn10+MInverse_ConcentrationD3Q15[1][11]*mfIn11+MInverse_ConcentrationD3Q15[1][12]*mfIn12+MInverse_ConcentrationD3Q15[1][13]*mfIn13+MInverse_ConcentrationD3Q15[1][14]*mfIn14)
#define vfConcentrationD3Q15In2 (MInverse_ConcentrationD3Q15[2][0]*mfIn0+MInverse_ConcentrationD3Q15[2][1]*mfIn1+MInverse_ConcentrationD3Q15[2][2]*mfIn2+MInverse_ConcentrationD3Q15[2][3]*mfIn3+MInverse_ConcentrationD3Q15[2][4]*mfIn4+MInverse_ConcentrationD3Q15[2][5]*mfIn5+MInverse_ConcentrationD3Q15[2][6]*mfIn6+MInverse_ConcentrationD3Q15[2][7]*mfIn7+MInverse_ConcentrationD3Q15[2][8]*mfIn8+MInverse_ConcentrationD3Q15[2][9]*mfIn9+MInverse_ConcentrationD3Q15[2][10]*mfIn10+MInverse_ConcentrationD3Q15[2][11]*mfIn11+MInverse_ConcentrationD3Q15[2][12]*mfIn12+MInverse_ConcentrationD3Q15[2][13]*mfIn13+MInverse_ConcentrationD3Q15[2][14]*mfIn14)
#define vfConcentrationD3Q15In3 (MInverse_ConcentrationD3Q15[3][0]*mfIn0+MInverse_ConcentrationD3Q15[3][1]*mfIn1+MInverse_ConcentrationD3Q15[3][2]*mfIn2+MInverse_ConcentrationD3Q15[3][3]*mfIn3+MInverse_ConcentrationD3Q15[3][4]*mfIn4+MInverse_ConcentrationD3Q15[3][5]*mfIn5+MInverse_ConcentrationD3Q15[3][6]*mfIn6+MInverse_ConcentrationD3Q15[3][7]*mfIn7+MInverse_ConcentrationD3Q15[3][8]*mfIn8+MInverse_ConcentrationD3Q15[3][9]*mfIn9+MInverse_ConcentrationD3Q15[3][10]*mfIn10+MInverse_ConcentrationD3Q15[3][11]*mfIn11+MInverse_ConcentrationD3Q15[3][12]*mfIn12+MInverse_ConcentrationD3Q15[3][13]*mfIn13+MInverse_ConcentrationD3Q15[3][14]*mfIn14)
#define vfConcentrationD3Q15In4 (MInverse_ConcentrationD3Q15[4][0]*mfIn0+MInverse_ConcentrationD3Q15[4][1]*mfIn1+MInverse_ConcentrationD3Q15[4][2]*mfIn2+MInverse_ConcentrationD3Q15[4][3]*mfIn3+MInverse_ConcentrationD3Q15[4][4]*mfIn4+MInverse_ConcentrationD3Q15[4][5]*mfIn5+MInverse_ConcentrationD3Q15[4][6]*mfIn6+MInverse_ConcentrationD3Q15[4][7]*mfIn7+MInverse_ConcentrationD3Q15[4][8]*mfIn8+MInverse_ConcentrationD3Q15[4][9]*mfIn9+MInverse_ConcentrationD3Q15[4][10]*mfIn10+MInverse_ConcentrationD3Q15[4][11]*mfIn11+MInverse_ConcentrationD3Q15[4][12]*mfIn12+MInverse_ConcentrationD3Q15[4][13]*mfIn13+MInverse_ConcentrationD3Q15[4][14]*mfIn14)
#define vfConcentrationD3Q15In5 (MInverse_ConcentrationD3Q15[5][0]*mfIn0+MInverse_ConcentrationD3Q15[5][1]*mfIn1+MInverse_ConcentrationD3Q15[5][2]*mfIn2+MInverse_ConcentrationD3Q15[5][3]*mfIn3+MInverse_ConcentrationD3Q15[5][4]*mfIn4+MInverse_ConcentrationD3Q15[5][5]*mfIn5+MInverse_ConcentrationD3Q15[5][6]*mfIn6+MInverse_ConcentrationD3Q15[5][7]*mfIn7+MInverse_ConcentrationD3Q15[5][8]*mfIn8+MInverse_ConcentrationD3Q15[5][9]*mfIn9+MInverse_ConcentrationD3Q15[5][10]*mfIn10+MInverse_ConcentrationD3Q15[5][11]*mfIn11+MInverse_ConcentrationD3Q15[5][12]*mfIn12+MInverse_ConcentrationD3Q15[5][13]*mfIn13+MInverse_ConcentrationD3Q15[5][14]*mfIn14)
#define vfConcentrationD3Q15In6 (MInverse_ConcentrationD3Q15[6][0]*mfIn0+MInverse_ConcentrationD3Q15[6][1]*mfIn1+MInverse_ConcentrationD3Q15[6][2]*mfIn2+MInverse_ConcentrationD3Q15[6][3]*mfIn3+MInverse_ConcentrationD3Q15[6][4]*mfIn4+MInverse_ConcentrationD3Q15[6][5]*mfIn5+MInverse_ConcentrationD3Q15[6][6]*mfIn6+MInverse_ConcentrationD3Q15[6][7]*mfIn7+MInverse_ConcentrationD3Q15[6][8]*mfIn8+MInverse_ConcentrationD3Q15[6][9]*mfIn9+MInverse_ConcentrationD3Q15[6][10]*mfIn10+MInverse_ConcentrationD3Q15[6][11]*mfIn11+MInverse_ConcentrationD3Q15[6][12]*mfIn12+MInverse_ConcentrationD3Q15[6][13]*mfIn13+MInverse_ConcentrationD3Q15[6][14]*mfIn14)
#define vfConcentrationD3Q15In7 (MInverse_ConcentrationD3Q15[7][0]*mfIn0+MInverse_ConcentrationD3Q15[7][1]*mfIn1+MInverse_ConcentrationD3Q15[7][2]*mfIn2+MInverse_ConcentrationD3Q15[7][3]*mfIn3+MInverse_ConcentrationD3Q15[7][4]*mfIn4+MInverse_ConcentrationD3Q15[7][5]*mfIn5+MInverse_ConcentrationD3Q15[7][6]*mfIn6+MInverse_ConcentrationD3Q15[7][7]*mfIn7+MInverse_ConcentrationD3Q15[7][8]*mfIn8+MInverse_ConcentrationD3Q15[7][9]*mfIn9+MInverse_ConcentrationD3Q15[7][10]*mfIn10+MInverse_ConcentrationD3Q15[7][11]*mfIn11+MInverse_ConcentrationD3Q15[7][12]*mfIn12+MInverse_ConcentrationD3Q15[7][13]*mfIn13+MInverse_ConcentrationD3Q15[7][14]*mfIn14)
#define vfConcentrationD3Q15In8 (MInverse_ConcentrationD3Q15[8][0]*mfIn0+MInverse_ConcentrationD3Q15[8][1]*mfIn1+MInverse_ConcentrationD3Q15[8][2]*mfIn2+MInverse_ConcentrationD3Q15[8][3]*mfIn3+MInverse_ConcentrationD3Q15[8][4]*mfIn4+MInverse_ConcentrationD3Q15[8][5]*mfIn5+MInverse_ConcentrationD3Q15[8][6]*mfIn6+MInverse_ConcentrationD3Q15[8][7]*mfIn7+MInverse_ConcentrationD3Q15[8][8]*mfIn8+MInverse_ConcentrationD3Q15[8][9]*mfIn9+MInverse_ConcentrationD3Q15[8][10]*mfIn10+MInverse_ConcentrationD3Q15[8][11]*mfIn11+MInverse_ConcentrationD3Q15[8][12]*mfIn12+MInverse_ConcentrationD3Q15[8][13]*mfIn13+MInverse_ConcentrationD3Q15[8][14]*mfIn14)
#define vfConcentrationD3Q15In9 (MInverse_ConcentrationD3Q15[9][0]*mfIn0+MInverse_ConcentrationD3Q15[9][1]*mfIn1+MInverse_ConcentrationD3Q15[9][2]*mfIn2+MInverse_ConcentrationD3Q15[9][3]*mfIn3+MInverse_ConcentrationD3Q15[9][4]*mfIn4+MInverse_ConcentrationD3Q15[9][5]*mfIn5+MInverse_ConcentrationD3Q15[9][6]*mfIn6+MInverse_ConcentrationD3Q15[9][7]*mfIn7+MInverse_ConcentrationD3Q15[9][8]*mfIn8+MInverse_ConcentrationD3Q15[9][9]*mfIn9+MInverse_ConcentrationD3Q15[9][10]*mfIn10+MInverse_ConcentrationD3Q15[9][11]*mfIn11+MInverse_ConcentrationD3Q15[9][12]*mfIn12+MInverse_ConcentrationD3Q15[9][13]*mfIn13+MInverse_ConcentrationD3Q15[9][14]*mfIn14)
#define vfConcentrationD3Q15In10 (MInverse_ConcentrationD3Q15[10][0]*mfIn0+MInverse_ConcentrationD3Q15[10][1]*mfIn1+MInverse_ConcentrationD3Q15[10][2]*mfIn2+MInverse_ConcentrationD3Q15[10][3]*mfIn3+MInverse_ConcentrationD3Q15[10][4]*mfIn4+MInverse_ConcentrationD3Q15[10][5]*mfIn5+MInverse_ConcentrationD3Q15[10][6]*mfIn6+MInverse_ConcentrationD3Q15[10][7]*mfIn7+MInverse_ConcentrationD3Q15[10][8]*mfIn8+MInverse_ConcentrationD3Q15[10][9]*mfIn9+MInverse_ConcentrationD3Q15[10][10]*mfIn10+MInverse_ConcentrationD3Q15[10][11]*mfIn11+MInverse_ConcentrationD3Q15[10][12]*mfIn12+MInverse_ConcentrationD3Q15[10][13]*mfIn13+MInverse_ConcentrationD3Q15[10][14]*mfIn14)
#define vfConcentrationD3Q15In11 (MInverse_ConcentrationD3Q15[11][0]*mfIn0+MInverse_ConcentrationD3Q15[11][1]*mfIn1+MInverse_ConcentrationD3Q15[11][2]*mfIn2+MInverse_ConcentrationD3Q15[11][3]*mfIn3+MInverse_ConcentrationD3Q15[11][4]*mfIn4+MInverse_ConcentrationD3Q15[11][5]*mfIn5+MInverse_ConcentrationD3Q15[11][6]*mfIn6+MInverse_ConcentrationD3Q15[11][7]*mfIn7+MInverse_ConcentrationD3Q15[11][8]*mfIn8+MInverse_ConcentrationD3Q15[11][9]*mfIn9+MInverse_ConcentrationD3Q15[11][10]*mfIn10+MInverse_ConcentrationD3Q15[11][11]*mfIn11+MInverse_ConcentrationD3Q15[11][12]*mfIn12+MInverse_ConcentrationD3Q15[11][13]*mfIn13+MInverse_ConcentrationD3Q15[11][14]*mfIn14)
#define vfConcentrationD3Q15In12 (MInverse_ConcentrationD3Q15[12][0]*mfIn0+MInverse_ConcentrationD3Q15[12][1]*mfIn1+MInverse_ConcentrationD3Q15[12][2]*mfIn2+MInverse_ConcentrationD3Q15[12][3]*mfIn3+MInverse_ConcentrationD3Q15[12][4]*mfIn4+MInverse_ConcentrationD3Q15[12][5]*mfIn5+MInverse_ConcentrationD3Q15[12][6]*mfIn6+MInverse_ConcentrationD3Q15[12][7]*mfIn7+MInverse_ConcentrationD3Q15[12][8]*mfIn8+MInverse_ConcentrationD3Q15[12][9]*mfIn9+MInverse_ConcentrationD3Q15[12][10]*mfIn10+MInverse_ConcentrationD3Q15[12][11]*mfIn11+MInverse_ConcentrationD3Q15[12][12]*mfIn12+MInverse_ConcentrationD3Q15[12][13]*mfIn13+MInverse_ConcentrationD3Q15[12][14]*mfIn14)
#define vfConcentrationD3Q15In13 (MInverse_ConcentrationD3Q15[13][0]*mfIn0+MInverse_ConcentrationD3Q15[13][1]*mfIn1+MInverse_ConcentrationD3Q15[13][2]*mfIn2+MInverse_ConcentrationD3Q15[13][3]*mfIn3+MInverse_ConcentrationD3Q15[13][4]*mfIn4+MInverse_ConcentrationD3Q15[13][5]*mfIn5+MInverse_ConcentrationD3Q15[13][6]*mfIn6+MInverse_ConcentrationD3Q15[13][7]*mfIn7+MInverse_ConcentrationD3Q15[13][8]*mfIn8+MInverse_ConcentrationD3Q15[13][9]*mfIn9+MInverse_ConcentrationD3Q15[13][10]*mfIn10+MInverse_ConcentrationD3Q15[13][11]*mfIn11+MInverse_ConcentrationD3Q15[13][12]*mfIn12+MInverse_ConcentrationD3Q15[13][13]*mfIn13+MInverse_ConcentrationD3Q15[13][14]*mfIn14)
#define vfConcentrationD3Q15In14 (MInverse_ConcentrationD3Q15[14][0]*mfIn0+MInverse_ConcentrationD3Q15[14][1]*mfIn1+MInverse_ConcentrationD3Q15[14][2]*mfIn2+MInverse_ConcentrationD3Q15[14][3]*mfIn3+MInverse_ConcentrationD3Q15[14][4]*mfIn4+MInverse_ConcentrationD3Q15[14][5]*mfIn5+MInverse_ConcentrationD3Q15[14][6]*mfIn6+MInverse_ConcentrationD3Q15[14][7]*mfIn7+MInverse_ConcentrationD3Q15[14][8]*mfIn8+MInverse_ConcentrationD3Q15[14][9]*mfIn9+MInverse_ConcentrationD3Q15[14][10]*mfIn10+MInverse_ConcentrationD3Q15[14][11]*mfIn11+MInverse_ConcentrationD3Q15[14][12]*mfIn12+MInverse_ConcentrationD3Q15[14][13]*mfIn13+MInverse_ConcentrationD3Q15[14][14]*mfIn14)

#define vfConcentrationD3Q15Eq0 (MInverse_ConcentrationD3Q15[0][0]*mfEq0+MInverse_ConcentrationD3Q15[0][1]*mfEq1+MInverse_ConcentrationD3Q15[0][2]*mfEq2+MInverse_ConcentrationD3Q15[0][3]*mfEq3+MInverse_ConcentrationD3Q15[0][4]*mfEq4+MInverse_ConcentrationD3Q15[0][5]*mfEq5+MInverse_ConcentrationD3Q15[0][6]*mfEq6+MInverse_ConcentrationD3Q15[0][7]*mfEq7+MInverse_ConcentrationD3Q15[0][8]*mfEq8+MInverse_ConcentrationD3Q15[0][9]*mfEq9+MInverse_ConcentrationD3Q15[0][10]*mfEq10+MInverse_ConcentrationD3Q15[0][11]*mfEq11+MInverse_ConcentrationD3Q15[0][12]*mfEq12+MInverse_ConcentrationD3Q15[0][13]*mfEq13+MInverse_ConcentrationD3Q15[0][14]*mfEq14)
#define vfConcentrationD3Q15Eq1 (MInverse_ConcentrationD3Q15[1][0]*mfEq0+MInverse_ConcentrationD3Q15[1][1]*mfEq1+MInverse_ConcentrationD3Q15[1][2]*mfEq2+MInverse_ConcentrationD3Q15[1][3]*mfEq3+MInverse_ConcentrationD3Q15[1][4]*mfEq4+MInverse_ConcentrationD3Q15[1][5]*mfEq5+MInverse_ConcentrationD3Q15[1][6]*mfEq6+MInverse_ConcentrationD3Q15[1][7]*mfEq7+MInverse_ConcentrationD3Q15[1][8]*mfEq8+MInverse_ConcentrationD3Q15[1][9]*mfEq9+MInverse_ConcentrationD3Q15[1][10]*mfEq10+MInverse_ConcentrationD3Q15[1][11]*mfEq11+MInverse_ConcentrationD3Q15[1][12]*mfEq12+MInverse_ConcentrationD3Q15[1][13]*mfEq13+MInverse_ConcentrationD3Q15[1][14]*mfEq14)
#define vfConcentrationD3Q15Eq2 (MInverse_ConcentrationD3Q15[2][0]*mfEq0+MInverse_ConcentrationD3Q15[2][1]*mfEq1+MInverse_ConcentrationD3Q15[2][2]*mfEq2+MInverse_ConcentrationD3Q15[2][3]*mfEq3+MInverse_ConcentrationD3Q15[2][4]*mfEq4+MInverse_ConcentrationD3Q15[2][5]*mfEq5+MInverse_ConcentrationD3Q15[2][6]*mfEq6+MInverse_ConcentrationD3Q15[2][7]*mfEq7+MInverse_ConcentrationD3Q15[2][8]*mfEq8+MInverse_ConcentrationD3Q15[2][9]*mfEq9+MInverse_ConcentrationD3Q15[2][10]*mfEq10+MInverse_ConcentrationD3Q15[2][11]*mfEq11+MInverse_ConcentrationD3Q15[2][12]*mfEq12+MInverse_ConcentrationD3Q15[2][13]*mfEq13+MInverse_ConcentrationD3Q15[2][14]*mfEq14)
#define vfConcentrationD3Q15Eq3 (MInverse_ConcentrationD3Q15[3][0]*mfEq0+MInverse_ConcentrationD3Q15[3][1]*mfEq1+MInverse_ConcentrationD3Q15[3][2]*mfEq2+MInverse_ConcentrationD3Q15[3][3]*mfEq3+MInverse_ConcentrationD3Q15[3][4]*mfEq4+MInverse_ConcentrationD3Q15[3][5]*mfEq5+MInverse_ConcentrationD3Q15[3][6]*mfEq6+MInverse_ConcentrationD3Q15[3][7]*mfEq7+MInverse_ConcentrationD3Q15[3][8]*mfEq8+MInverse_ConcentrationD3Q15[3][9]*mfEq9+MInverse_ConcentrationD3Q15[3][10]*mfEq10+MInverse_ConcentrationD3Q15[3][11]*mfEq11+MInverse_ConcentrationD3Q15[3][12]*mfEq12+MInverse_ConcentrationD3Q15[3][13]*mfEq13+MInverse_ConcentrationD3Q15[3][14]*mfEq14)
#define vfConcentrationD3Q15Eq4 (MInverse_ConcentrationD3Q15[4][0]*mfEq0+MInverse_ConcentrationD3Q15[4][1]*mfEq1+MInverse_ConcentrationD3Q15[4][2]*mfEq2+MInverse_ConcentrationD3Q15[4][3]*mfEq3+MInverse_ConcentrationD3Q15[4][4]*mfEq4+MInverse_ConcentrationD3Q15[4][5]*mfEq5+MInverse_ConcentrationD3Q15[4][6]*mfEq6+MInverse_ConcentrationD3Q15[4][7]*mfEq7+MInverse_ConcentrationD3Q15[4][8]*mfEq8+MInverse_ConcentrationD3Q15[4][9]*mfEq9+MInverse_ConcentrationD3Q15[4][10]*mfEq10+MInverse_ConcentrationD3Q15[4][11]*mfEq11+MInverse_ConcentrationD3Q15[4][12]*mfEq12+MInverse_ConcentrationD3Q15[4][13]*mfEq13+MInverse_ConcentrationD3Q15[4][14]*mfEq14)
#define vfConcentrationD3Q15Eq5 (MInverse_ConcentrationD3Q15[5][0]*mfEq0+MInverse_ConcentrationD3Q15[5][1]*mfEq1+MInverse_ConcentrationD3Q15[5][2]*mfEq2+MInverse_ConcentrationD3Q15[5][3]*mfEq3+MInverse_ConcentrationD3Q15[5][4]*mfEq4+MInverse_ConcentrationD3Q15[5][5]*mfEq5+MInverse_ConcentrationD3Q15[5][6]*mfEq6+MInverse_ConcentrationD3Q15[5][7]*mfEq7+MInverse_ConcentrationD3Q15[5][8]*mfEq8+MInverse_ConcentrationD3Q15[5][9]*mfEq9+MInverse_ConcentrationD3Q15[5][10]*mfEq10+MInverse_ConcentrationD3Q15[5][11]*mfEq11+MInverse_ConcentrationD3Q15[5][12]*mfEq12+MInverse_ConcentrationD3Q15[5][13]*mfEq13+MInverse_ConcentrationD3Q15[5][14]*mfEq14)
#define vfConcentrationD3Q15Eq6 (MInverse_ConcentrationD3Q15[6][0]*mfEq0+MInverse_ConcentrationD3Q15[6][1]*mfEq1+MInverse_ConcentrationD3Q15[6][2]*mfEq2+MInverse_ConcentrationD3Q15[6][3]*mfEq3+MInverse_ConcentrationD3Q15[6][4]*mfEq4+MInverse_ConcentrationD3Q15[6][5]*mfEq5+MInverse_ConcentrationD3Q15[6][6]*mfEq6+MInverse_ConcentrationD3Q15[6][7]*mfEq7+MInverse_ConcentrationD3Q15[6][8]*mfEq8+MInverse_ConcentrationD3Q15[6][9]*mfEq9+MInverse_ConcentrationD3Q15[6][10]*mfEq10+MInverse_ConcentrationD3Q15[6][11]*mfEq11+MInverse_ConcentrationD3Q15[6][12]*mfEq12+MInverse_ConcentrationD3Q15[6][13]*mfEq13+MInverse_ConcentrationD3Q15[6][14]*mfEq14)
#define vfConcentrationD3Q15Eq7 (MInverse_ConcentrationD3Q15[7][0]*mfEq0+MInverse_ConcentrationD3Q15[7][1]*mfEq1+MInverse_ConcentrationD3Q15[7][2]*mfEq2+MInverse_ConcentrationD3Q15[7][3]*mfEq3+MInverse_ConcentrationD3Q15[7][4]*mfEq4+MInverse_ConcentrationD3Q15[7][5]*mfEq5+MInverse_ConcentrationD3Q15[7][6]*mfEq6+MInverse_ConcentrationD3Q15[7][7]*mfEq7+MInverse_ConcentrationD3Q15[7][8]*mfEq8+MInverse_ConcentrationD3Q15[7][9]*mfEq9+MInverse_ConcentrationD3Q15[7][10]*mfEq10+MInverse_ConcentrationD3Q15[7][11]*mfEq11+MInverse_ConcentrationD3Q15[7][12]*mfEq12+MInverse_ConcentrationD3Q15[7][13]*mfEq13+MInverse_ConcentrationD3Q15[7][14]*mfEq14)
#define vfConcentrationD3Q15Eq8 (MInverse_ConcentrationD3Q15[8][0]*mfEq0+MInverse_ConcentrationD3Q15[8][1]*mfEq1+MInverse_ConcentrationD3Q15[8][2]*mfEq2+MInverse_ConcentrationD3Q15[8][3]*mfEq3+MInverse_ConcentrationD3Q15[8][4]*mfEq4+MInverse_ConcentrationD3Q15[8][5]*mfEq5+MInverse_ConcentrationD3Q15[8][6]*mfEq6+MInverse_ConcentrationD3Q15[8][7]*mfEq7+MInverse_ConcentrationD3Q15[8][8]*mfEq8+MInverse_ConcentrationD3Q15[8][9]*mfEq9+MInverse_ConcentrationD3Q15[8][10]*mfEq10+MInverse_ConcentrationD3Q15[8][11]*mfEq11+MInverse_ConcentrationD3Q15[8][12]*mfEq12+MInverse_ConcentrationD3Q15[8][13]*mfEq13+MInverse_ConcentrationD3Q15[8][14]*mfEq14)
#define vfConcentrationD3Q15Eq9 (MInverse_ConcentrationD3Q15[9][0]*mfEq0+MInverse_ConcentrationD3Q15[9][1]*mfEq1+MInverse_ConcentrationD3Q15[9][2]*mfEq2+MInverse_ConcentrationD3Q15[9][3]*mfEq3+MInverse_ConcentrationD3Q15[9][4]*mfEq4+MInverse_ConcentrationD3Q15[9][5]*mfEq5+MInverse_ConcentrationD3Q15[9][6]*mfEq6+MInverse_ConcentrationD3Q15[9][7]*mfEq7+MInverse_ConcentrationD3Q15[9][8]*mfEq8+MInverse_ConcentrationD3Q15[9][9]*mfEq9+MInverse_ConcentrationD3Q15[9][10]*mfEq10+MInverse_ConcentrationD3Q15[9][11]*mfEq11+MInverse_ConcentrationD3Q15[9][12]*mfEq12+MInverse_ConcentrationD3Q15[9][13]*mfEq13+MInverse_ConcentrationD3Q15[9][14]*mfEq14)
#define vfConcentrationD3Q15Eq10 (MInverse_ConcentrationD3Q15[10][0]*mfEq0+MInverse_ConcentrationD3Q15[10][1]*mfEq1+MInverse_ConcentrationD3Q15[10][2]*mfEq2+MInverse_ConcentrationD3Q15[10][3]*mfEq3+MInverse_ConcentrationD3Q15[10][4]*mfEq4+MInverse_ConcentrationD3Q15[10][5]*mfEq5+MInverse_ConcentrationD3Q15[10][6]*mfEq6+MInverse_ConcentrationD3Q15[10][7]*mfEq7+MInverse_ConcentrationD3Q15[10][8]*mfEq8+MInverse_ConcentrationD3Q15[10][9]*mfEq9+MInverse_ConcentrationD3Q15[10][10]*mfEq10+MInverse_ConcentrationD3Q15[10][11]*mfEq11+MInverse_ConcentrationD3Q15[10][12]*mfEq12+MInverse_ConcentrationD3Q15[10][13]*mfEq13+MInverse_ConcentrationD3Q15[10][14]*mfEq14)
#define vfConcentrationD3Q15Eq11 (MInverse_ConcentrationD3Q15[11][0]*mfEq0+MInverse_ConcentrationD3Q15[11][1]*mfEq1+MInverse_ConcentrationD3Q15[11][2]*mfEq2+MInverse_ConcentrationD3Q15[11][3]*mfEq3+MInverse_ConcentrationD3Q15[11][4]*mfEq4+MInverse_ConcentrationD3Q15[11][5]*mfEq5+MInverse_ConcentrationD3Q15[11][6]*mfEq6+MInverse_ConcentrationD3Q15[11][7]*mfEq7+MInverse_ConcentrationD3Q15[11][8]*mfEq8+MInverse_ConcentrationD3Q15[11][9]*mfEq9+MInverse_ConcentrationD3Q15[11][10]*mfEq10+MInverse_ConcentrationD3Q15[11][11]*mfEq11+MInverse_ConcentrationD3Q15[11][12]*mfEq12+MInverse_ConcentrationD3Q15[11][13]*mfEq13+MInverse_ConcentrationD3Q15[11][14]*mfEq14)
#define vfConcentrationD3Q15Eq12 (MInverse_ConcentrationD3Q15[12][0]*mfEq0+MInverse_ConcentrationD3Q15[12][1]*mfEq1+MInverse_ConcentrationD3Q15[12][2]*mfEq2+MInverse_ConcentrationD3Q15[12][3]*mfEq3+MInverse_ConcentrationD3Q15[12][4]*mfEq4+MInverse_ConcentrationD3Q15[12][5]*mfEq5+MInverse_ConcentrationD3Q15[12][6]*mfEq6+MInverse_ConcentrationD3Q15[12][7]*mfEq7+MInverse_ConcentrationD3Q15[12][8]*mfEq8+MInverse_ConcentrationD3Q15[12][9]*mfEq9+MInverse_ConcentrationD3Q15[12][10]*mfEq10+MInverse_ConcentrationD3Q15[12][11]*mfEq11+MInverse_ConcentrationD3Q15[12][12]*mfEq12+MInverse_ConcentrationD3Q15[12][13]*mfEq13+MInverse_ConcentrationD3Q15[12][14]*mfEq14)
#define vfConcentrationD3Q15Eq13 (MInverse_ConcentrationD3Q15[13][0]*mfEq0+MInverse_ConcentrationD3Q15[13][1]*mfEq1+MInverse_ConcentrationD3Q15[13][2]*mfEq2+MInverse_ConcentrationD3Q15[13][3]*mfEq3+MInverse_ConcentrationD3Q15[13][4]*mfEq4+MInverse_ConcentrationD3Q15[13][5]*mfEq5+MInverse_ConcentrationD3Q15[13][6]*mfEq6+MInverse_ConcentrationD3Q15[13][7]*mfEq7+MInverse_ConcentrationD3Q15[13][8]*mfEq8+MInverse_ConcentrationD3Q15[13][9]*mfEq9+MInverse_ConcentrationD3Q15[13][10]*mfEq10+MInverse_ConcentrationD3Q15[13][11]*mfEq11+MInverse_ConcentrationD3Q15[13][12]*mfEq12+MInverse_ConcentrationD3Q15[13][13]*mfEq13+MInverse_ConcentrationD3Q15[13][14]*mfEq14)
#define vfConcentrationD3Q15Eq14 (MInverse_ConcentrationD3Q15[14][0]*mfEq0+MInverse_ConcentrationD3Q15[14][1]*mfEq1+MInverse_ConcentrationD3Q15[14][2]*mfEq2+MInverse_ConcentrationD3Q15[14][3]*mfEq3+MInverse_ConcentrationD3Q15[14][4]*mfEq4+MInverse_ConcentrationD3Q15[14][5]*mfEq5+MInverse_ConcentrationD3Q15[14][6]*mfEq6+MInverse_ConcentrationD3Q15[14][7]*mfEq7+MInverse_ConcentrationD3Q15[14][8]*mfEq8+MInverse_ConcentrationD3Q15[14][9]*mfEq9+MInverse_ConcentrationD3Q15[14][10]*mfEq10+MInverse_ConcentrationD3Q15[14][11]*mfEq11+MInverse_ConcentrationD3Q15[14][12]*mfEq12+MInverse_ConcentrationD3Q15[14][13]*mfEq13+MInverse_ConcentrationD3Q15[14][14]*mfEq14)

double ***cIn0=NULL;
double ***cIn1=NULL;
double ***cIn2=NULL;
double ***cIn3=NULL;
double ***cIn4=NULL;
double ***cIn5=NULL;
double ***cIn6=NULL;
double ***cIn7=NULL;
double ***cIn8=NULL;
double ***cIn9=NULL;
double ***cIn10=NULL;
double ***cIn11=NULL;
double ***cIn12=NULL;
double ***cIn13=NULL;
double ***cIn14=NULL;
double ***cIn15=NULL;
double ***cIn16=NULL;
double ***cIn17=NULL;
double ***cIn18=NULL;
double ***Concentration=NULL;
double ***domainConnection=NULL;
double ***wallNearConcentration=NULL;
double ***wallNearConcentrationFirst=NULL;
double ***convectTermConcentration = NULL;
double ***diffusionTermConcentration = NULL;
double ***reactantSurface1 = NULL;
double ***reactantSurface2 = NULL;
double ***reactantSurface3 = NULL;
double ***reactantSurface4 = NULL;
double ***reactantSurface5 = NULL;
double ***reactantSurface6 = NULL;
double ***massTransCoeMatrix = NULL;
double ***electrolyteConcentrationRelaxationTime = NULL;

double ***cIn0First=NULL;
double ***cIn1First=NULL;
double ***cIn2First=NULL;
double ***cIn3First=NULL;
double ***cIn4First=NULL;
double ***cIn5First=NULL;
double ***cIn6First=NULL;
double ***cIn7First=NULL;
double ***cIn8First=NULL;
double ***cIn9First=NULL;
double ***cIn10First=NULL;
double ***cIn11First=NULL;
double ***cIn12First=NULL;
double ***cIn13First=NULL;
double ***cIn14First=NULL;
double ***cIn15First=NULL;
double ***cIn16First=NULL;
double ***cIn17First=NULL;
double ***cIn18First=NULL;
double ***ConcentrationFirst=NULL;
double ***convectTermConcentrationFirst = NULL;
double ***diffusionTermConcentrationFirst = NULL;

double ***cIn0Second=NULL;
double ***cIn1Second=NULL;
double ***cIn2Second=NULL;
double ***cIn3Second=NULL;
double ***cIn4Second=NULL;
double ***cIn5Second=NULL;
double ***cIn6Second=NULL;
double ***cIn7Second=NULL;
double ***cIn8Second=NULL;
double ***cIn9Second=NULL;
double ***cIn10Second=NULL;
double ***cIn11Second=NULL;
double ***cIn12Second=NULL;
double ***cIn13Second=NULL;
double ***cIn14Second=NULL;
double ***cIn15Second=NULL;
double ***cIn16Second=NULL;
double ***cIn17Second=NULL;
double ***cIn18Second=NULL;
double ***ConcentrationSecond=NULL;
double ***convectTermConcentrationSecond = NULL;

void FieldArrangeConcentration(int Length, int Width, int Height){
	cIn0=memoryarrange(Length,Width,Height);
	cIn1=memoryarrange(Length,Width,Height);
	cIn2=memoryarrange(Length,Width,Height);
	cIn3=memoryarrange(Length,Width,Height);
	cIn4=memoryarrange(Length,Width,Height);
	cIn5=memoryarrange(Length,Width,Height);
	cIn6=memoryarrange(Length,Width,Height);
	cIn7=memoryarrange(Length,Width,Height);
	cIn8=memoryarrange(Length,Width,Height);
	cIn9=memoryarrange(Length,Width,Height);
	cIn10=memoryarrange(Length,Width,Height);
	cIn11=memoryarrange(Length,Width,Height);
	cIn12=memoryarrange(Length,Width,Height);
	cIn13=memoryarrange(Length,Width,Height);
	cIn14=memoryarrange(Length,Width,Height);
	cIn15=memoryarrange(Length,Width,Height);
	cIn16=memoryarrange(Length,Width,Height);
	cIn17=memoryarrange(Length,Width,Height);
	cIn18=memoryarrange(Length,Width,Height);
	Concentration=memoryarrange(Length,Width,Height); 
	domainConnection=memoryarrange(Length,Width,Height);
	wallNearConcentration=memoryarrange(Length,Width,Height);
	wallNearConcentrationFirst=memoryarrange(Length,Width,Height);
	convectTermConcentration = memoryarrange(Length, Width, Height);
	diffusionTermConcentration = memoryarrange(Length, Width, Height);
	reactantSurface1 = memoryarrange(Length, Width, Height);
	reactantSurface2 = memoryarrange(Length, Width, Height);
	reactantSurface3 = memoryarrange(Length, Width, Height);
	reactantSurface4 = memoryarrange(Length, Width, Height);
	reactantSurface5 = memoryarrange(Length, Width, Height);
	reactantSurface6 = memoryarrange(Length, Width, Height);
	massTransCoeMatrix = memoryarrange(Length, Width, Height);
	electrolyteConcentrationRelaxationTime = memoryarrange(Length, Width, Height);
}

void FieldArrangeConcentrationFirst(int Length, int Width, int Height){
	cIn0First=memoryarrange(Length,Width,Height);
	cIn1First=memoryarrange(Length,Width,Height);
	cIn2First=memoryarrange(Length,Width,Height);
	cIn3First=memoryarrange(Length,Width,Height);
	cIn4First=memoryarrange(Length,Width,Height);
	cIn5First=memoryarrange(Length,Width,Height);
	cIn6First=memoryarrange(Length,Width,Height);
	cIn7First=memoryarrange(Length,Width,Height);
	cIn8First=memoryarrange(Length,Width,Height);
	cIn9First=memoryarrange(Length,Width,Height);
	cIn10First=memoryarrange(Length,Width,Height);
	cIn11First=memoryarrange(Length,Width,Height);
	cIn12First=memoryarrange(Length,Width,Height);
	cIn13First=memoryarrange(Length,Width,Height);
	cIn14First=memoryarrange(Length,Width,Height);
	cIn15First=memoryarrange(Length,Width,Height);
	cIn16First=memoryarrange(Length,Width,Height);
	cIn17First=memoryarrange(Length,Width,Height);
	cIn18First=memoryarrange(Length,Width,Height);
	ConcentrationFirst=memoryarrange(Length,Width,Height); 
	convectTermConcentrationFirst = memoryarrange(Length, Width, Height);
	diffusionTermConcentrationFirst = memoryarrange(Length, Width, Height);
}

void FieldArrangeConcentrationSecond(int Length, int Width, int Height){
	cIn0Second=memoryarrange(Length,Width,Height);
	cIn1Second=memoryarrange(Length,Width,Height);
	cIn2Second=memoryarrange(Length,Width,Height);
	cIn3Second=memoryarrange(Length,Width,Height);
	cIn4Second=memoryarrange(Length,Width,Height);
	cIn5Second=memoryarrange(Length,Width,Height);
	cIn6Second=memoryarrange(Length,Width,Height);
	cIn7Second=memoryarrange(Length,Width,Height);
	cIn8Second=memoryarrange(Length,Width,Height);
	cIn9Second=memoryarrange(Length,Width,Height);
	cIn10Second=memoryarrange(Length,Width,Height);
	cIn11Second=memoryarrange(Length,Width,Height);
	cIn12Second=memoryarrange(Length,Width,Height);
	cIn13Second=memoryarrange(Length,Width,Height);
	cIn14Second=memoryarrange(Length,Width,Height);
	cIn15Second=memoryarrange(Length,Width,Height);
	cIn16Second=memoryarrange(Length,Width,Height);
	cIn17Second=memoryarrange(Length,Width,Height);
	cIn18Second=memoryarrange(Length,Width,Height);
	ConcentrationSecond=memoryarrange(Length,Width,Height); 
	convectTermConcentrationSecond = memoryarrange(Length, Width, Height);
}

void domainConnectionRegion(int Length, int Width, int Height, double gasCritical){
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
				domainConnection[i][j][k]=0;
			}
		}
	}
	domainConnection[0][0][0]=1;
	k2=0;
	k1=1;
	while (k1-k2>0){
		k2=k1;
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
					if (domainConnection[i][j][k]==0&&domain[i][j][k]==0&&rho[i][j][k]>gasCritical){
						Neightbour1=neightbour1(i,j,k,Length,domainConnection);
						Neightbour2=neightbour2(i,j,k,Width,domainConnection);
						Neightbour3=neightbour3(i,j,k,Height,domainConnection);
						Neightbour4=neightbour4(i,j,k,Length,domainConnection);
						Neightbour5=neightbour5(i,j,k,Width,domainConnection);
						Neightbour6=neightbour6(i,j,k,Height,domainConnection);
						Neightbour7=neightbour7(i,j,k,Length,Width,domainConnection);
						Neightbour8=neightbour8(i,j,k,Length,Height,domainConnection);
						Neightbour9=neightbour9(i,j,k,Width,Height,domainConnection);
						Neightbour10=neightbour10(i,j,k,Length,Width,domainConnection);
						Neightbour11=neightbour11(i,j,k,Length,Height,domainConnection);
						Neightbour12=neightbour12(i,j,k,Width,Height,domainConnection);
						Neightbour13=neightbour13(i,j,k,Length,Width,domainConnection);
						Neightbour14=neightbour14(i,j,k,Length,Height,domainConnection);
						Neightbour15=neightbour15(i,j,k,Width,Height,domainConnection);
						Neightbour16=neightbour16(i,j,k,Length,Width,domainConnection);
						Neightbour17=neightbour17(i,j,k,Length,Height,domainConnection);
						Neightbour18=neightbour18(i,j,k,Width,Height,domainConnection);
						if (Neightbour1==1||Neightbour2==1||Neightbour3==1||Neightbour4==1||Neightbour5==1||Neightbour6==1||Neightbour7==1||Neightbour8==1||\
							Neightbour9==1||Neightbour10==1||Neightbour11==1||Neightbour12==1||Neightbour13==1||Neightbour14==1||Neightbour15==1||Neightbour16==1||Neightbour17==1||Neightbour18==1){
								domainConnection[i][j][k]=1;
								k1=k1+1;
						}
					}
				}
			}
		}
	}
}

void domainConnectionRegionD3Q7(int Length, int Width, int Height, double gasCritical){
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
				domainConnection[i][j][k]=0;
			}
		}
	}
	domainConnection[0][0][0]=1;
	k2=0;
	k1=1;
	while (k1-k2>0){
		k2=k1;
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
					if (domainConnection[i][j][k]==0&&domain[i][j][k]==0&&rho[i][j][k]>gasCritical){
						Neightbour1=neightbour1(i,j,k,Length,domainConnection);
						Neightbour2=neightbour2(i,j,k,Width,domainConnection);
						Neightbour3=neightbour3(i,j,k,Height,domainConnection);
						Neightbour4=neightbour4(i,j,k,Length,domainConnection);
						Neightbour5=neightbour5(i,j,k,Width,domainConnection);
						Neightbour6=neightbour6(i,j,k,Height,domainConnection);
						if (Neightbour1==1||Neightbour2==1||Neightbour3==1||Neightbour4==1||Neightbour5==1||Neightbour6==1){
								domainConnection[i][j][k]=1;
								k1=k1+1;
						}
					}
				}
			}
		}
	}
}

void FieldInitialConcentration(int Length, int Width, int Height, double ConcentrationInitial, double gasCritical, double ***ux, double ***uy, double ***uz){
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
						Concentration[i][j][k]=ConcentrationInitial;
						cIn0[i][j][k]=Concentration[i][j][k]*tNS[0];//*(1+3*(cxNS[0]*ux[i][j][k]+cyNS[0]*uy[i][j][k]+czNS[0]*uz[i][j][k])+4.5*((cxNS[0]*ux[i][j][k]+cyNS[0]*uy[i][j][k]+czNS[0]*uz[i][j][k])*(cxNS[0]*ux[i][j][k]+cyNS[0]*uy[i][j][k]+czNS[0]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn1[i][j][k]=Concentration[i][j][k]*tNS[1];//*(1+3*(cxNS[1]*ux[i][j][k]+cyNS[1]*uy[i][j][k]+czNS[1]*uz[i][j][k])+4.5*((cxNS[1]*ux[i][j][k]+cyNS[1]*uy[i][j][k]+czNS[1]*uz[i][j][k])*(cxNS[1]*ux[i][j][k]+cyNS[1]*uy[i][j][k]+czNS[1]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn2[i][j][k]=Concentration[i][j][k]*tNS[2];//*(1+3*(cxNS[2]*ux[i][j][k]+cyNS[2]*uy[i][j][k]+czNS[2]*uz[i][j][k])+4.5*((cxNS[2]*ux[i][j][k]+cyNS[2]*uy[i][j][k]+czNS[2]*uz[i][j][k])*(cxNS[2]*ux[i][j][k]+cyNS[2]*uy[i][j][k]+czNS[2]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn3[i][j][k]=Concentration[i][j][k]*tNS[3];//*(1+3*(cxNS[3]*ux[i][j][k]+cyNS[3]*uy[i][j][k]+czNS[3]*uz[i][j][k])+4.5*((cxNS[3]*ux[i][j][k]+cyNS[3]*uy[i][j][k]+czNS[3]*uz[i][j][k])*(cxNS[3]*ux[i][j][k]+cyNS[3]*uy[i][j][k]+czNS[3]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn4[i][j][k]=Concentration[i][j][k]*tNS[4];//*(1+3*(cxNS[4]*ux[i][j][k]+cyNS[4]*uy[i][j][k]+czNS[4]*uz[i][j][k])+4.5*((cxNS[4]*ux[i][j][k]+cyNS[4]*uy[i][j][k]+czNS[4]*uz[i][j][k])*(cxNS[4]*ux[i][j][k]+cyNS[4]*uy[i][j][k]+czNS[4]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn5[i][j][k]=Concentration[i][j][k]*tNS[5];//*(1+3*(cxNS[5]*ux[i][j][k]+cyNS[5]*uy[i][j][k]+czNS[5]*uz[i][j][k])+4.5*((cxNS[5]*ux[i][j][k]+cyNS[5]*uy[i][j][k]+czNS[5]*uz[i][j][k])*(cxNS[5]*ux[i][j][k]+cyNS[5]*uy[i][j][k]+czNS[5]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn6[i][j][k]=Concentration[i][j][k]*tNS[6];//*(1+3*(cxNS[6]*ux[i][j][k]+cyNS[6]*uy[i][j][k]+czNS[6]*uz[i][j][k])+4.5*((cxNS[6]*ux[i][j][k]+cyNS[6]*uy[i][j][k]+czNS[6]*uz[i][j][k])*(cxNS[6]*ux[i][j][k]+cyNS[6]*uy[i][j][k]+czNS[6]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn7[i][j][k]=Concentration[i][j][k]*tNS[7];//*(1+3*(cxNS[7]*ux[i][j][k]+cyNS[7]*uy[i][j][k]+czNS[7]*uz[i][j][k])+4.5*((cxNS[7]*ux[i][j][k]+cyNS[7]*uy[i][j][k]+czNS[7]*uz[i][j][k])*(cxNS[7]*ux[i][j][k]+cyNS[7]*uy[i][j][k]+czNS[7]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn8[i][j][k]=Concentration[i][j][k]*tNS[8];//*(1+3*(cxNS[8]*ux[i][j][k]+cyNS[8]*uy[i][j][k]+czNS[8]*uz[i][j][k])+4.5*((cxNS[8]*ux[i][j][k]+cyNS[8]*uy[i][j][k]+czNS[8]*uz[i][j][k])*(cxNS[8]*ux[i][j][k]+cyNS[8]*uy[i][j][k]+czNS[8]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn9[i][j][k]=Concentration[i][j][k]*tNS[9];//*(1+3*(cxNS[9]*ux[i][j][k]+cyNS[9]*uy[i][j][k]+czNS[9]*uz[i][j][k])+4.5*((cxNS[9]*ux[i][j][k]+cyNS[9]*uy[i][j][k]+czNS[9]*uz[i][j][k])*(cxNS[9]*ux[i][j][k]+cyNS[9]*uy[i][j][k]+czNS[9]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn10[i][j][k]=Concentration[i][j][k]*tNS[10];//*(1+3*(cxNS[10]*ux[i][j][k]+cyNS[10]*uy[i][j][k]+czNS[10]*uz[i][j][k])+4.5*((cxNS[10]*ux[i][j][k]+cyNS[10]*uy[i][j][k]+czNS[10]*uz[i][j][k])*(cxNS[10]*ux[i][j][k]+cyNS[10]*uy[i][j][k]+czNS[10]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn11[i][j][k]=Concentration[i][j][k]*tNS[11];//*(1+3*(cxNS[11]*ux[i][j][k]+cyNS[11]*uy[i][j][k]+czNS[11]*uz[i][j][k])+4.5*((cxNS[11]*ux[i][j][k]+cyNS[11]*uy[i][j][k]+czNS[11]*uz[i][j][k])*(cxNS[11]*ux[i][j][k]+cyNS[11]*uy[i][j][k]+czNS[11]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn12[i][j][k]=Concentration[i][j][k]*tNS[12];//*(1+3*(cxNS[12]*ux[i][j][k]+cyNS[12]*uy[i][j][k]+czNS[12]*uz[i][j][k])+4.5*((cxNS[12]*ux[i][j][k]+cyNS[12]*uy[i][j][k]+czNS[12]*uz[i][j][k])*(cxNS[12]*ux[i][j][k]+cyNS[12]*uy[i][j][k]+czNS[12]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn13[i][j][k]=Concentration[i][j][k]*tNS[13];//*(1+3*(cxNS[13]*ux[i][j][k]+cyNS[13]*uy[i][j][k]+czNS[13]*uz[i][j][k])+4.5*((cxNS[13]*ux[i][j][k]+cyNS[13]*uy[i][j][k]+czNS[13]*uz[i][j][k])*(cxNS[13]*ux[i][j][k]+cyNS[13]*uy[i][j][k]+czNS[13]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn14[i][j][k]=Concentration[i][j][k]*tNS[14];//*(1+3*(cxNS[14]*ux[i][j][k]+cyNS[14]*uy[i][j][k]+czNS[14]*uz[i][j][k])+4.5*((cxNS[14]*ux[i][j][k]+cyNS[14]*uy[i][j][k]+czNS[14]*uz[i][j][k])*(cxNS[14]*ux[i][j][k]+cyNS[14]*uy[i][j][k]+czNS[14]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn15[i][j][k]=Concentration[i][j][k]*tNS[15];//*(1+3*(cxNS[15]*ux[i][j][k]+cyNS[15]*uy[i][j][k]+czNS[15]*uz[i][j][k])+4.5*((cxNS[15]*ux[i][j][k]+cyNS[15]*uy[i][j][k]+czNS[15]*uz[i][j][k])*(cxNS[15]*ux[i][j][k]+cyNS[15]*uy[i][j][k]+czNS[15]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn16[i][j][k]=Concentration[i][j][k]*tNS[16];//*(1+3*(cxNS[16]*ux[i][j][k]+cyNS[16]*uy[i][j][k]+czNS[16]*uz[i][j][k])+4.5*((cxNS[16]*ux[i][j][k]+cyNS[16]*uy[i][j][k]+czNS[16]*uz[i][j][k])*(cxNS[16]*ux[i][j][k]+cyNS[16]*uy[i][j][k]+czNS[16]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn17[i][j][k]=Concentration[i][j][k]*tNS[17];//*(1+3*(cxNS[17]*ux[i][j][k]+cyNS[17]*uy[i][j][k]+czNS[17]*uz[i][j][k])+4.5*((cxNS[17]*ux[i][j][k]+cyNS[17]*uy[i][j][k]+czNS[17]*uz[i][j][k])*(cxNS[17]*ux[i][j][k]+cyNS[17]*uy[i][j][k]+czNS[17]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn18[i][j][k]=Concentration[i][j][k]*tNS[18];//*(1+3*(cxNS[18]*ux[i][j][k]+cyNS[18]*uy[i][j][k]+czNS[18]*uz[i][j][k])+4.5*((cxNS[18]*ux[i][j][k]+cyNS[18]*uy[i][j][k]+czNS[18]*uz[i][j][k])*(cxNS[18]*ux[i][j][k]+cyNS[18]*uy[i][j][k]+czNS[18]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
					}
					if (domain[i][j][k]!=0||rho[i][j][k]<=gasCritical||domainConnection[i][j][k]!=1){
						Concentration[i][j][k]=0.0;
						cIn0[i][j][k]=0.0;
						cIn1[i][j][k]=0.0;
						cIn2[i][j][k]=0.0;
						cIn3[i][j][k]=0.0;
						cIn4[i][j][k]=0.0;
						cIn5[i][j][k]=0.0;
						cIn6[i][j][k]=0.0;
						cIn7[i][j][k]=0.0;
						cIn8[i][j][k]=0.0;
						cIn9[i][j][k]=0.0;
						cIn10[i][j][k]=0.0;
						cIn11[i][j][k]=0.0;
						cIn12[i][j][k]=0.0;
						cIn13[i][j][k]=0.0;
						cIn14[i][j][k]=0.0;
						cIn15[i][j][k]=0.0;
						cIn16[i][j][k]=0.0;
						cIn17[i][j][k]=0.0;
						cIn18[i][j][k]=0.0;
					}
				}
			}
		}
	}
}

void FieldInitialConcentrationElectrolyteD3Q7(int Length, int Width, int Height, double ConcentrationInitial, double gasCritical){
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
						Concentration[i][j][k]=ConcentrationInitial;
						cIn0[i][j][k]=Concentration[i][j][k]*tNSD3Q7[0];
						cIn1[i][j][k]=Concentration[i][j][k]*tNSD3Q7[1];
						cIn2[i][j][k]=Concentration[i][j][k]*tNSD3Q7[2];
						cIn3[i][j][k]=Concentration[i][j][k]*tNSD3Q7[3];
						cIn4[i][j][k]=Concentration[i][j][k]*tNSD3Q7[4];
						cIn5[i][j][k]=Concentration[i][j][k]*tNSD3Q7[5];
						cIn6[i][j][k]=Concentration[i][j][k]*tNSD3Q7[6];
					}
					if (domain[i][j][k]!=0||rho[i][j][k]<=gasCritical||domainConnection[i][j][k]!=1){
						Concentration[i][j][k]=0.0;
						cIn0[i][j][k]=0.0;
						cIn1[i][j][k]=0.0;
						cIn2[i][j][k]=0.0;
						cIn3[i][j][k]=0.0;
						cIn4[i][j][k]=0.0;
						cIn5[i][j][k]=0.0;
						cIn6[i][j][k]=0.0;
					}
				}
			}
		}
	}
}

void FieldInitialConcentrationD3Q15(int Length, int Width, int Height, double ConcentrationInitial, double gasCritical, double ***ux, double ***uy, double ***uz){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (j = 0; j<Width; j++){
				for (k = 0; k<Height; k++){
					if (domain[i][j][k] == 0 && rho[i][j][k]>gasCritical&&domainConnection[i][j][k] == 1){
						Concentration[i][j][k] = ConcentrationInitial;
						cIn0[i][j][k] = Concentration[i][j][k] * tNSD3Q15[0];
						cIn1[i][j][k] = Concentration[i][j][k] * tNSD3Q15[1];
						cIn2[i][j][k] = Concentration[i][j][k] * tNSD3Q15[2];
						cIn3[i][j][k] = Concentration[i][j][k] * tNSD3Q15[3];
						cIn4[i][j][k] = Concentration[i][j][k] * tNSD3Q15[4];
						cIn5[i][j][k] = Concentration[i][j][k] * tNSD3Q15[5];
						cIn6[i][j][k] = Concentration[i][j][k] * tNSD3Q15[6];
						cIn7[i][j][k] = Concentration[i][j][k] * tNSD3Q15[7];
						cIn8[i][j][k] = Concentration[i][j][k] * tNSD3Q15[8];
						cIn9[i][j][k] = Concentration[i][j][k] * tNSD3Q15[9];
						cIn10[i][j][k] = Concentration[i][j][k] * tNSD3Q15[10];
						cIn11[i][j][k] = Concentration[i][j][k] * tNSD3Q15[11];
						cIn12[i][j][k] = Concentration[i][j][k] * tNSD3Q15[12];
						cIn13[i][j][k] = Concentration[i][j][k] * tNSD3Q15[13];
						cIn14[i][j][k] = Concentration[i][j][k] * tNSD3Q15[14];
					}
					if (domain[i][j][k] != 0 || rho[i][j][k] <= gasCritical || domainConnection[i][j][k] != 1){
						Concentration[i][j][k] = 0.0;
						cIn0[i][j][k] = 0.0;
						cIn1[i][j][k] = 0.0;
						cIn2[i][j][k] = 0.0;
						cIn3[i][j][k] = 0.0;
						cIn4[i][j][k] = 0.0;
						cIn5[i][j][k] = 0.0;
						cIn6[i][j][k] = 0.0;
						cIn7[i][j][k] = 0.0;
						cIn8[i][j][k] = 0.0;
						cIn9[i][j][k] = 0.0;
						cIn10[i][j][k] = 0.0;
						cIn11[i][j][k] = 0.0;
						cIn12[i][j][k] = 0.0;
						cIn13[i][j][k] = 0.0;
						cIn14[i][j][k] = 0.0;
					}
				}
			}
		}
	}
}

void FieldInitialConcentrationFirst(int Length, int Width, int Height, double ConcentrationInitial, double gasCritical, double ***ux, double ***uy, double ***uz){
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
						ConcentrationFirst[i][j][k]=ConcentrationInitial;
						cIn0First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[0];//*(1+3*(cxNS[0]*ux[i][j][k]+cyNS[0]*uy[i][j][k]+czNS[0]*uz[i][j][k])+4.5*((cxNS[0]*ux[i][j][k]+cyNS[0]*uy[i][j][k]+czNS[0]*uz[i][j][k])*(cxNS[0]*ux[i][j][k]+cyNS[0]*uy[i][j][k]+czNS[0]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn1First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[1];//*(1+3*(cxNS[1]*ux[i][j][k]+cyNS[1]*uy[i][j][k]+czNS[1]*uz[i][j][k])+4.5*((cxNS[1]*ux[i][j][k]+cyNS[1]*uy[i][j][k]+czNS[1]*uz[i][j][k])*(cxNS[1]*ux[i][j][k]+cyNS[1]*uy[i][j][k]+czNS[1]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn2First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[2];//*(1+3*(cxNS[2]*ux[i][j][k]+cyNS[2]*uy[i][j][k]+czNS[2]*uz[i][j][k])+4.5*((cxNS[2]*ux[i][j][k]+cyNS[2]*uy[i][j][k]+czNS[2]*uz[i][j][k])*(cxNS[2]*ux[i][j][k]+cyNS[2]*uy[i][j][k]+czNS[2]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn3First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[3];//*(1+3*(cxNS[3]*ux[i][j][k]+cyNS[3]*uy[i][j][k]+czNS[3]*uz[i][j][k])+4.5*((cxNS[3]*ux[i][j][k]+cyNS[3]*uy[i][j][k]+czNS[3]*uz[i][j][k])*(cxNS[3]*ux[i][j][k]+cyNS[3]*uy[i][j][k]+czNS[3]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn4First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[4];//*(1+3*(cxNS[4]*ux[i][j][k]+cyNS[4]*uy[i][j][k]+czNS[4]*uz[i][j][k])+4.5*((cxNS[4]*ux[i][j][k]+cyNS[4]*uy[i][j][k]+czNS[4]*uz[i][j][k])*(cxNS[4]*ux[i][j][k]+cyNS[4]*uy[i][j][k]+czNS[4]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn5First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[5];//*(1+3*(cxNS[5]*ux[i][j][k]+cyNS[5]*uy[i][j][k]+czNS[5]*uz[i][j][k])+4.5*((cxNS[5]*ux[i][j][k]+cyNS[5]*uy[i][j][k]+czNS[5]*uz[i][j][k])*(cxNS[5]*ux[i][j][k]+cyNS[5]*uy[i][j][k]+czNS[5]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn6First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[6];//*(1+3*(cxNS[6]*ux[i][j][k]+cyNS[6]*uy[i][j][k]+czNS[6]*uz[i][j][k])+4.5*((cxNS[6]*ux[i][j][k]+cyNS[6]*uy[i][j][k]+czNS[6]*uz[i][j][k])*(cxNS[6]*ux[i][j][k]+cyNS[6]*uy[i][j][k]+czNS[6]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn7First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[7];//*(1+3*(cxNS[7]*ux[i][j][k]+cyNS[7]*uy[i][j][k]+czNS[7]*uz[i][j][k])+4.5*((cxNS[7]*ux[i][j][k]+cyNS[7]*uy[i][j][k]+czNS[7]*uz[i][j][k])*(cxNS[7]*ux[i][j][k]+cyNS[7]*uy[i][j][k]+czNS[7]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn8First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[8];//*(1+3*(cxNS[8]*ux[i][j][k]+cyNS[8]*uy[i][j][k]+czNS[8]*uz[i][j][k])+4.5*((cxNS[8]*ux[i][j][k]+cyNS[8]*uy[i][j][k]+czNS[8]*uz[i][j][k])*(cxNS[8]*ux[i][j][k]+cyNS[8]*uy[i][j][k]+czNS[8]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn9First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[9];//*(1+3*(cxNS[9]*ux[i][j][k]+cyNS[9]*uy[i][j][k]+czNS[9]*uz[i][j][k])+4.5*((cxNS[9]*ux[i][j][k]+cyNS[9]*uy[i][j][k]+czNS[9]*uz[i][j][k])*(cxNS[9]*ux[i][j][k]+cyNS[9]*uy[i][j][k]+czNS[9]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn10First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[10];//*(1+3*(cxNS[10]*ux[i][j][k]+cyNS[10]*uy[i][j][k]+czNS[10]*uz[i][j][k])+4.5*((cxNS[10]*ux[i][j][k]+cyNS[10]*uy[i][j][k]+czNS[10]*uz[i][j][k])*(cxNS[10]*ux[i][j][k]+cyNS[10]*uy[i][j][k]+czNS[10]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn11First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[11];//*(1+3*(cxNS[11]*ux[i][j][k]+cyNS[11]*uy[i][j][k]+czNS[11]*uz[i][j][k])+4.5*((cxNS[11]*ux[i][j][k]+cyNS[11]*uy[i][j][k]+czNS[11]*uz[i][j][k])*(cxNS[11]*ux[i][j][k]+cyNS[11]*uy[i][j][k]+czNS[11]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn12First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[12];//*(1+3*(cxNS[12]*ux[i][j][k]+cyNS[12]*uy[i][j][k]+czNS[12]*uz[i][j][k])+4.5*((cxNS[12]*ux[i][j][k]+cyNS[12]*uy[i][j][k]+czNS[12]*uz[i][j][k])*(cxNS[12]*ux[i][j][k]+cyNS[12]*uy[i][j][k]+czNS[12]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn13First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[13];//*(1+3*(cxNS[13]*ux[i][j][k]+cyNS[13]*uy[i][j][k]+czNS[13]*uz[i][j][k])+4.5*((cxNS[13]*ux[i][j][k]+cyNS[13]*uy[i][j][k]+czNS[13]*uz[i][j][k])*(cxNS[13]*ux[i][j][k]+cyNS[13]*uy[i][j][k]+czNS[13]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn14First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[14];//*(1+3*(cxNS[14]*ux[i][j][k]+cyNS[14]*uy[i][j][k]+czNS[14]*uz[i][j][k])+4.5*((cxNS[14]*ux[i][j][k]+cyNS[14]*uy[i][j][k]+czNS[14]*uz[i][j][k])*(cxNS[14]*ux[i][j][k]+cyNS[14]*uy[i][j][k]+czNS[14]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn15First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[15];//*(1+3*(cxNS[15]*ux[i][j][k]+cyNS[15]*uy[i][j][k]+czNS[15]*uz[i][j][k])+4.5*((cxNS[15]*ux[i][j][k]+cyNS[15]*uy[i][j][k]+czNS[15]*uz[i][j][k])*(cxNS[15]*ux[i][j][k]+cyNS[15]*uy[i][j][k]+czNS[15]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn16First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[16];//*(1+3*(cxNS[16]*ux[i][j][k]+cyNS[16]*uy[i][j][k]+czNS[16]*uz[i][j][k])+4.5*((cxNS[16]*ux[i][j][k]+cyNS[16]*uy[i][j][k]+czNS[16]*uz[i][j][k])*(cxNS[16]*ux[i][j][k]+cyNS[16]*uy[i][j][k]+czNS[16]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn17First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[17];//*(1+3*(cxNS[17]*ux[i][j][k]+cyNS[17]*uy[i][j][k]+czNS[17]*uz[i][j][k])+4.5*((cxNS[17]*ux[i][j][k]+cyNS[17]*uy[i][j][k]+czNS[17]*uz[i][j][k])*(cxNS[17]*ux[i][j][k]+cyNS[17]*uy[i][j][k]+czNS[17]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn18First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[18];//*(1+3*(cxNS[18]*ux[i][j][k]+cyNS[18]*uy[i][j][k]+czNS[18]*uz[i][j][k])+4.5*((cxNS[18]*ux[i][j][k]+cyNS[18]*uy[i][j][k]+czNS[18]*uz[i][j][k])*(cxNS[18]*ux[i][j][k]+cyNS[18]*uy[i][j][k]+czNS[18]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
					}
					if (domain[i][j][k]!=0||rho[i][j][k]<=gasCritical||domainConnection[i][j][k]!=1){
						ConcentrationFirst[i][j][k]=0.0;
						cIn0First[i][j][k]=0.0;
						cIn1First[i][j][k]=0.0;
						cIn2First[i][j][k]=0.0;
						cIn3First[i][j][k]=0.0;
						cIn4First[i][j][k]=0.0;
						cIn5First[i][j][k]=0.0;
						cIn6First[i][j][k]=0.0;
						cIn7First[i][j][k]=0.0;
						cIn8First[i][j][k]=0.0;
						cIn9First[i][j][k]=0.0;
						cIn10First[i][j][k]=0.0;
						cIn11First[i][j][k]=0.0;
						cIn12First[i][j][k]=0.0;
						cIn13First[i][j][k]=0.0;
						cIn14First[i][j][k]=0.0;
						cIn15First[i][j][k]=0.0;
						cIn16First[i][j][k]=0.0;
						cIn17First[i][j][k]=0.0;
						cIn18First[i][j][k]=0.0;
					}
				}
			}
		}
	}
}

void FieldInitialConcentrationElectrodeD3Q7(int Length, int Width, int Height, double ConcentrationInitial, double gasCritical){
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
						ConcentrationFirst[i][j][k]=ConcentrationInitial;
						cIn0First[i][j][k]=ConcentrationFirst[i][j][k]*tNSD3Q7[0];//*(1+3*(cxNS[0]*ux[i][j][k]+cyNS[0]*uy[i][j][k]+czNS[0]*uz[i][j][k])+4.5*((cxNS[0]*ux[i][j][k]+cyNS[0]*uy[i][j][k]+czNS[0]*uz[i][j][k])*(cxNS[0]*ux[i][j][k]+cyNS[0]*uy[i][j][k]+czNS[0]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn1First[i][j][k]=ConcentrationFirst[i][j][k]*tNSD3Q7[1];//*(1+3*(cxNS[1]*ux[i][j][k]+cyNS[1]*uy[i][j][k]+czNS[1]*uz[i][j][k])+4.5*((cxNS[1]*ux[i][j][k]+cyNS[1]*uy[i][j][k]+czNS[1]*uz[i][j][k])*(cxNS[1]*ux[i][j][k]+cyNS[1]*uy[i][j][k]+czNS[1]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn2First[i][j][k]=ConcentrationFirst[i][j][k]*tNSD3Q7[2];//*(1+3*(cxNS[2]*ux[i][j][k]+cyNS[2]*uy[i][j][k]+czNS[2]*uz[i][j][k])+4.5*((cxNS[2]*ux[i][j][k]+cyNS[2]*uy[i][j][k]+czNS[2]*uz[i][j][k])*(cxNS[2]*ux[i][j][k]+cyNS[2]*uy[i][j][k]+czNS[2]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn3First[i][j][k]=ConcentrationFirst[i][j][k]*tNSD3Q7[3];//*(1+3*(cxNS[3]*ux[i][j][k]+cyNS[3]*uy[i][j][k]+czNS[3]*uz[i][j][k])+4.5*((cxNS[3]*ux[i][j][k]+cyNS[3]*uy[i][j][k]+czNS[3]*uz[i][j][k])*(cxNS[3]*ux[i][j][k]+cyNS[3]*uy[i][j][k]+czNS[3]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn4First[i][j][k]=ConcentrationFirst[i][j][k]*tNSD3Q7[4];//*(1+3*(cxNS[4]*ux[i][j][k]+cyNS[4]*uy[i][j][k]+czNS[4]*uz[i][j][k])+4.5*((cxNS[4]*ux[i][j][k]+cyNS[4]*uy[i][j][k]+czNS[4]*uz[i][j][k])*(cxNS[4]*ux[i][j][k]+cyNS[4]*uy[i][j][k]+czNS[4]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn5First[i][j][k]=ConcentrationFirst[i][j][k]*tNSD3Q7[5];//*(1+3*(cxNS[5]*ux[i][j][k]+cyNS[5]*uy[i][j][k]+czNS[5]*uz[i][j][k])+4.5*((cxNS[5]*ux[i][j][k]+cyNS[5]*uy[i][j][k]+czNS[5]*uz[i][j][k])*(cxNS[5]*ux[i][j][k]+cyNS[5]*uy[i][j][k]+czNS[5]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn6First[i][j][k]=ConcentrationFirst[i][j][k]*tNSD3Q7[6];//*(1+3*(cxNS[6]*ux[i][j][k]+cyNS[6]*uy[i][j][k]+czNS[6]*uz[i][j][k])+4.5*((cxNS[6]*ux[i][j][k]+cyNS[6]*uy[i][j][k]+czNS[6]*uz[i][j][k])*(cxNS[6]*ux[i][j][k]+cyNS[6]*uy[i][j][k]+czNS[6]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
					}
					if (domain[i][j][k]!=1||electrodeConnection[i][j][k]!=1){
						ConcentrationFirst[i][j][k] = 0.0;
						cIn0First[i][j][k] = 0.0;
						cIn1First[i][j][k] = 0.0;
						cIn2First[i][j][k] = 0.0;
						cIn3First[i][j][k] = 0.0;
						cIn4First[i][j][k] = 0.0;
						cIn5First[i][j][k] = 0.0;
						cIn6First[i][j][k] = 0.0;
					}
				}
			}
		}
	}
}

void FieldInitialConcentrationSecond(int Length, int Width, int Height, double ConcentrationInitial, double gasCritical, double ***ux, double ***uy, double ***uz){
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
						ConcentrationSecond[i][j][k]=ConcentrationInitial;
						cIn0Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[0];//*(1+3*(cxNS[0]*ux[i][j][k]+cyNS[0]*uy[i][j][k]+czNS[0]*uz[i][j][k])+4.5*((cxNS[0]*ux[i][j][k]+cyNS[0]*uy[i][j][k]+czNS[0]*uz[i][j][k])*(cxNS[0]*ux[i][j][k]+cyNS[0]*uy[i][j][k]+czNS[0]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn1Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[1];//*(1+3*(cxNS[1]*ux[i][j][k]+cyNS[1]*uy[i][j][k]+czNS[1]*uz[i][j][k])+4.5*((cxNS[1]*ux[i][j][k]+cyNS[1]*uy[i][j][k]+czNS[1]*uz[i][j][k])*(cxNS[1]*ux[i][j][k]+cyNS[1]*uy[i][j][k]+czNS[1]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn2Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[2];//*(1+3*(cxNS[2]*ux[i][j][k]+cyNS[2]*uy[i][j][k]+czNS[2]*uz[i][j][k])+4.5*((cxNS[2]*ux[i][j][k]+cyNS[2]*uy[i][j][k]+czNS[2]*uz[i][j][k])*(cxNS[2]*ux[i][j][k]+cyNS[2]*uy[i][j][k]+czNS[2]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn3Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[3];//*(1+3*(cxNS[3]*ux[i][j][k]+cyNS[3]*uy[i][j][k]+czNS[3]*uz[i][j][k])+4.5*((cxNS[3]*ux[i][j][k]+cyNS[3]*uy[i][j][k]+czNS[3]*uz[i][j][k])*(cxNS[3]*ux[i][j][k]+cyNS[3]*uy[i][j][k]+czNS[3]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn4Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[4];//*(1+3*(cxNS[4]*ux[i][j][k]+cyNS[4]*uy[i][j][k]+czNS[4]*uz[i][j][k])+4.5*((cxNS[4]*ux[i][j][k]+cyNS[4]*uy[i][j][k]+czNS[4]*uz[i][j][k])*(cxNS[4]*ux[i][j][k]+cyNS[4]*uy[i][j][k]+czNS[4]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn5Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[5];//*(1+3*(cxNS[5]*ux[i][j][k]+cyNS[5]*uy[i][j][k]+czNS[5]*uz[i][j][k])+4.5*((cxNS[5]*ux[i][j][k]+cyNS[5]*uy[i][j][k]+czNS[5]*uz[i][j][k])*(cxNS[5]*ux[i][j][k]+cyNS[5]*uy[i][j][k]+czNS[5]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn6Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[6];//*(1+3*(cxNS[6]*ux[i][j][k]+cyNS[6]*uy[i][j][k]+czNS[6]*uz[i][j][k])+4.5*((cxNS[6]*ux[i][j][k]+cyNS[6]*uy[i][j][k]+czNS[6]*uz[i][j][k])*(cxNS[6]*ux[i][j][k]+cyNS[6]*uy[i][j][k]+czNS[6]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn7Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[7];//*(1+3*(cxNS[7]*ux[i][j][k]+cyNS[7]*uy[i][j][k]+czNS[7]*uz[i][j][k])+4.5*((cxNS[7]*ux[i][j][k]+cyNS[7]*uy[i][j][k]+czNS[7]*uz[i][j][k])*(cxNS[7]*ux[i][j][k]+cyNS[7]*uy[i][j][k]+czNS[7]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn8Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[8];//*(1+3*(cxNS[8]*ux[i][j][k]+cyNS[8]*uy[i][j][k]+czNS[8]*uz[i][j][k])+4.5*((cxNS[8]*ux[i][j][k]+cyNS[8]*uy[i][j][k]+czNS[8]*uz[i][j][k])*(cxNS[8]*ux[i][j][k]+cyNS[8]*uy[i][j][k]+czNS[8]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn9Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[9];//*(1+3*(cxNS[9]*ux[i][j][k]+cyNS[9]*uy[i][j][k]+czNS[9]*uz[i][j][k])+4.5*((cxNS[9]*ux[i][j][k]+cyNS[9]*uy[i][j][k]+czNS[9]*uz[i][j][k])*(cxNS[9]*ux[i][j][k]+cyNS[9]*uy[i][j][k]+czNS[9]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn10Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[10];//*(1+3*(cxNS[10]*ux[i][j][k]+cyNS[10]*uy[i][j][k]+czNS[10]*uz[i][j][k])+4.5*((cxNS[10]*ux[i][j][k]+cyNS[10]*uy[i][j][k]+czNS[10]*uz[i][j][k])*(cxNS[10]*ux[i][j][k]+cyNS[10]*uy[i][j][k]+czNS[10]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn11Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[11];//*(1+3*(cxNS[11]*ux[i][j][k]+cyNS[11]*uy[i][j][k]+czNS[11]*uz[i][j][k])+4.5*((cxNS[11]*ux[i][j][k]+cyNS[11]*uy[i][j][k]+czNS[11]*uz[i][j][k])*(cxNS[11]*ux[i][j][k]+cyNS[11]*uy[i][j][k]+czNS[11]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn12Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[12];//*(1+3*(cxNS[12]*ux[i][j][k]+cyNS[12]*uy[i][j][k]+czNS[12]*uz[i][j][k])+4.5*((cxNS[12]*ux[i][j][k]+cyNS[12]*uy[i][j][k]+czNS[12]*uz[i][j][k])*(cxNS[12]*ux[i][j][k]+cyNS[12]*uy[i][j][k]+czNS[12]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn13Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[13];//*(1+3*(cxNS[13]*ux[i][j][k]+cyNS[13]*uy[i][j][k]+czNS[13]*uz[i][j][k])+4.5*((cxNS[13]*ux[i][j][k]+cyNS[13]*uy[i][j][k]+czNS[13]*uz[i][j][k])*(cxNS[13]*ux[i][j][k]+cyNS[13]*uy[i][j][k]+czNS[13]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn14Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[14];//*(1+3*(cxNS[14]*ux[i][j][k]+cyNS[14]*uy[i][j][k]+czNS[14]*uz[i][j][k])+4.5*((cxNS[14]*ux[i][j][k]+cyNS[14]*uy[i][j][k]+czNS[14]*uz[i][j][k])*(cxNS[14]*ux[i][j][k]+cyNS[14]*uy[i][j][k]+czNS[14]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn15Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[15];//*(1+3*(cxNS[15]*ux[i][j][k]+cyNS[15]*uy[i][j][k]+czNS[15]*uz[i][j][k])+4.5*((cxNS[15]*ux[i][j][k]+cyNS[15]*uy[i][j][k]+czNS[15]*uz[i][j][k])*(cxNS[15]*ux[i][j][k]+cyNS[15]*uy[i][j][k]+czNS[15]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn16Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[16];//*(1+3*(cxNS[16]*ux[i][j][k]+cyNS[16]*uy[i][j][k]+czNS[16]*uz[i][j][k])+4.5*((cxNS[16]*ux[i][j][k]+cyNS[16]*uy[i][j][k]+czNS[16]*uz[i][j][k])*(cxNS[16]*ux[i][j][k]+cyNS[16]*uy[i][j][k]+czNS[16]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn17Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[17];//*(1+3*(cxNS[17]*ux[i][j][k]+cyNS[17]*uy[i][j][k]+czNS[17]*uz[i][j][k])+4.5*((cxNS[17]*ux[i][j][k]+cyNS[17]*uy[i][j][k]+czNS[17]*uz[i][j][k])*(cxNS[17]*ux[i][j][k]+cyNS[17]*uy[i][j][k]+czNS[17]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn18Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[18];//*(1+3*(cxNS[18]*ux[i][j][k]+cyNS[18]*uy[i][j][k]+czNS[18]*uz[i][j][k])+4.5*((cxNS[18]*ux[i][j][k]+cyNS[18]*uy[i][j][k]+czNS[18]*uz[i][j][k])*(cxNS[18]*ux[i][j][k]+cyNS[18]*uy[i][j][k]+czNS[18]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
					}
					if (domain[i][j][k]!=0||rho[i][j][k]<=gasCritical||domainConnection[i][j][k]!=1){
						ConcentrationSecond[i][j][k]=0.0;
						cIn0Second[i][j][k]=0.0;
						cIn1Second[i][j][k]=0.0;
						cIn2Second[i][j][k]=0.0;
						cIn3Second[i][j][k]=0.0;
						cIn4Second[i][j][k]=0.0;
						cIn5Second[i][j][k]=0.0;
						cIn6Second[i][j][k]=0.0;
						cIn7Second[i][j][k]=0.0;
						cIn8Second[i][j][k]=0.0;
						cIn9Second[i][j][k]=0.0;
						cIn10Second[i][j][k]=0.0;
						cIn11Second[i][j][k]=0.0;
						cIn12Second[i][j][k]=0.0;
						cIn13Second[i][j][k]=0.0;
						cIn14Second[i][j][k]=0.0;
						cIn15Second[i][j][k]=0.0;
						cIn16Second[i][j][k]=0.0;
						cIn17Second[i][j][k]=0.0;
						cIn18Second[i][j][k]=0.0;
					}
				}
			}
		}
	}
}

void FieldInitialConcentrationSecondD3Q7(int Length, int Width, int Height, double ConcentrationInitial, double gasCritical){
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
						ConcentrationSecond[i][j][k]=ConcentrationInitial;
						cIn0Second[i][j][k]=ConcentrationSecond[i][j][k]*tNSD3Q7[0];//*(1+3*(cxNS[0]*ux[i][j][k]+cyNS[0]*uy[i][j][k]+czNS[0]*uz[i][j][k])+4.5*((cxNS[0]*ux[i][j][k]+cyNS[0]*uy[i][j][k]+czNS[0]*uz[i][j][k])*(cxNS[0]*ux[i][j][k]+cyNS[0]*uy[i][j][k]+czNS[0]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn1Second[i][j][k]=ConcentrationSecond[i][j][k]*tNSD3Q7[1];//*(1+3*(cxNS[1]*ux[i][j][k]+cyNS[1]*uy[i][j][k]+czNS[1]*uz[i][j][k])+4.5*((cxNS[1]*ux[i][j][k]+cyNS[1]*uy[i][j][k]+czNS[1]*uz[i][j][k])*(cxNS[1]*ux[i][j][k]+cyNS[1]*uy[i][j][k]+czNS[1]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn2Second[i][j][k]=ConcentrationSecond[i][j][k]*tNSD3Q7[2];//*(1+3*(cxNS[2]*ux[i][j][k]+cyNS[2]*uy[i][j][k]+czNS[2]*uz[i][j][k])+4.5*((cxNS[2]*ux[i][j][k]+cyNS[2]*uy[i][j][k]+czNS[2]*uz[i][j][k])*(cxNS[2]*ux[i][j][k]+cyNS[2]*uy[i][j][k]+czNS[2]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn3Second[i][j][k]=ConcentrationSecond[i][j][k]*tNSD3Q7[3];//*(1+3*(cxNS[3]*ux[i][j][k]+cyNS[3]*uy[i][j][k]+czNS[3]*uz[i][j][k])+4.5*((cxNS[3]*ux[i][j][k]+cyNS[3]*uy[i][j][k]+czNS[3]*uz[i][j][k])*(cxNS[3]*ux[i][j][k]+cyNS[3]*uy[i][j][k]+czNS[3]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn4Second[i][j][k]=ConcentrationSecond[i][j][k]*tNSD3Q7[4];//*(1+3*(cxNS[4]*ux[i][j][k]+cyNS[4]*uy[i][j][k]+czNS[4]*uz[i][j][k])+4.5*((cxNS[4]*ux[i][j][k]+cyNS[4]*uy[i][j][k]+czNS[4]*uz[i][j][k])*(cxNS[4]*ux[i][j][k]+cyNS[4]*uy[i][j][k]+czNS[4]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn5Second[i][j][k]=ConcentrationSecond[i][j][k]*tNSD3Q7[5];//*(1+3*(cxNS[5]*ux[i][j][k]+cyNS[5]*uy[i][j][k]+czNS[5]*uz[i][j][k])+4.5*((cxNS[5]*ux[i][j][k]+cyNS[5]*uy[i][j][k]+czNS[5]*uz[i][j][k])*(cxNS[5]*ux[i][j][k]+cyNS[5]*uy[i][j][k]+czNS[5]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
						cIn6Second[i][j][k]=ConcentrationSecond[i][j][k]*tNSD3Q7[6];//*(1+3*(cxNS[6]*ux[i][j][k]+cyNS[6]*uy[i][j][k]+czNS[6]*uz[i][j][k])+4.5*((cxNS[6]*ux[i][j][k]+cyNS[6]*uy[i][j][k]+czNS[6]*uz[i][j][k])*(cxNS[6]*ux[i][j][k]+cyNS[6]*uy[i][j][k]+czNS[6]*uz[i][j][k]))-1.5*(ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]));
					}
					if (domain[i][j][k]!=0||rho[i][j][k]<=gasCritical||domainConnection[i][j][k]!=1){
						ConcentrationSecond[i][j][k]=0.0;
						cIn0Second[i][j][k]=0.0;
						cIn1Second[i][j][k]=0.0;
						cIn2Second[i][j][k]=0.0;
						cIn3Second[i][j][k]=0.0;
						cIn4Second[i][j][k]=0.0;
						cIn5Second[i][j][k]=0.0;
						cIn6Second[i][j][k]=0.0;
					}
				}
			}
		}
	}
}

void FieldInitialReactantSurface(int Length, int Width, int Height){
	int i;
	int j;
	int k;
	for (i = 0; i < Length; i++){
		for (j = 0; j < Width; j++){
			for (k = 0; k < Height; k++){
				reactantSurface1[i][j][k] = 0;
				reactantSurface2[i][j][k] = 0;
				reactantSurface3[i][j][k] = 0;
				reactantSurface4[i][j][k] = 0;
				reactantSurface5[i][j][k] = 0;
				reactantSurface6[i][j][k] = 0;
			}
		}
	}
}

void CollisionSRTConcentration(int Length,int Width,int Height,double CrelaxationTime,double gasCritical){
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
						Out0[i][j][k]=cIn0[i][j][k]-CrelaxationTime*(cIn0[i][j][k]-Eq0[i][j][k]*Concentration[i][j][k]);
						Out1[i][j][k]=cIn1[i][j][k]-CrelaxationTime*(cIn1[i][j][k]-Eq1[i][j][k]*Concentration[i][j][k]);
						Out2[i][j][k]=cIn2[i][j][k]-CrelaxationTime*(cIn2[i][j][k]-Eq2[i][j][k]*Concentration[i][j][k]);
						Out3[i][j][k]=cIn3[i][j][k]-CrelaxationTime*(cIn3[i][j][k]-Eq3[i][j][k]*Concentration[i][j][k]);
						Out4[i][j][k]=cIn4[i][j][k]-CrelaxationTime*(cIn4[i][j][k]-Eq4[i][j][k]*Concentration[i][j][k]);					
						Out5[i][j][k]=cIn5[i][j][k]-CrelaxationTime*(cIn5[i][j][k]-Eq5[i][j][k]*Concentration[i][j][k]);
						Out6[i][j][k]=cIn6[i][j][k]-CrelaxationTime*(cIn6[i][j][k]-Eq6[i][j][k]*Concentration[i][j][k]);
						Out7[i][j][k]=cIn7[i][j][k]-CrelaxationTime*(cIn7[i][j][k]-Eq7[i][j][k]*Concentration[i][j][k]);
						Out8[i][j][k]=cIn8[i][j][k]-CrelaxationTime*(cIn8[i][j][k]-Eq8[i][j][k]*Concentration[i][j][k]);
						Out9[i][j][k]=cIn9[i][j][k]-CrelaxationTime*(cIn9[i][j][k]-Eq9[i][j][k]*Concentration[i][j][k]);
						Out10[i][j][k]=cIn10[i][j][k]-CrelaxationTime*(cIn10[i][j][k]-Eq10[i][j][k]*Concentration[i][j][k]);
						Out11[i][j][k]=cIn11[i][j][k]-CrelaxationTime*(cIn11[i][j][k]-Eq11[i][j][k]*Concentration[i][j][k]);
						Out12[i][j][k]=cIn12[i][j][k]-CrelaxationTime*(cIn12[i][j][k]-Eq12[i][j][k]*Concentration[i][j][k]);
						Out13[i][j][k]=cIn13[i][j][k]-CrelaxationTime*(cIn13[i][j][k]-Eq13[i][j][k]*Concentration[i][j][k]);
						Out14[i][j][k]=cIn14[i][j][k]-CrelaxationTime*(cIn14[i][j][k]-Eq14[i][j][k]*Concentration[i][j][k]);
						Out15[i][j][k]=cIn15[i][j][k]-CrelaxationTime*(cIn15[i][j][k]-Eq15[i][j][k]*Concentration[i][j][k]);
						Out16[i][j][k]=cIn16[i][j][k]-CrelaxationTime*(cIn16[i][j][k]-Eq16[i][j][k]*Concentration[i][j][k]);
						Out17[i][j][k]=cIn17[i][j][k]-CrelaxationTime*(cIn17[i][j][k]-Eq17[i][j][k]*Concentration[i][j][k]);
						Out18[i][j][k]=cIn18[i][j][k]-CrelaxationTime*(cIn18[i][j][k]-Eq18[i][j][k]*Concentration[i][j][k]);
					}
				}
			}
		}
	}
}

void CollisionMRTConcentrationD3Q7(int Length, int Width, int Height, double gasCritical, double ***Ux, double ***Uy, double ***Uz, double a){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k,mfIn0,mfIn1,mfIn2,mfIn3,mfIn4,mfIn5,mfIn6,mfEq0, mfEq1, mfEq2, mfEq3, mfEq4, mfEq5, mfEq6)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (j = 0; j<Width; j++){
				for (k = 0; k<Height; k++){
					if (domain[i][j][k] == 0 && rho[i][j][k]>gasCritical&&domainConnection[i][j][k] == 1){
						mfIn0 = (A_gConcentration[0][0] * (M_Concentration[0][0] * cIn0[i][j][k] + M_Concentration[0][1] * cIn1[i][j][k] + M_Concentration[0][2] * cIn2[i][j][k] + M_Concentration[0][3] * cIn3[i][j][k] + M_Concentration[0][4] * cIn4[i][j][k] + M_Concentration[0][5] * cIn5[i][j][k] + M_Concentration[0][6] * cIn6[i][j][k]));
						mfIn1 = (A_gConcentration[1][1] * (M_Concentration[1][0] * cIn0[i][j][k] + M_Concentration[1][1] * cIn1[i][j][k] + M_Concentration[1][2] * cIn2[i][j][k] + M_Concentration[1][3] * cIn3[i][j][k] + M_Concentration[1][4] * cIn4[i][j][k] + M_Concentration[1][5] * cIn5[i][j][k] + M_Concentration[1][6] * cIn6[i][j][k]));
						mfIn2 = (A_gConcentration[2][2] * (M_Concentration[2][0] * cIn0[i][j][k] + M_Concentration[2][1] * cIn1[i][j][k] + M_Concentration[2][2] * cIn2[i][j][k] + M_Concentration[2][3] * cIn3[i][j][k] + M_Concentration[2][4] * cIn4[i][j][k] + M_Concentration[2][5] * cIn5[i][j][k] + M_Concentration[2][6] * cIn6[i][j][k]));
						mfIn3 = (A_gConcentration[3][3] * (M_Concentration[3][0] * cIn0[i][j][k] + M_Concentration[3][1] * cIn1[i][j][k] + M_Concentration[3][2] * cIn2[i][j][k] + M_Concentration[3][3] * cIn3[i][j][k] + M_Concentration[3][4] * cIn4[i][j][k] + M_Concentration[3][5] * cIn5[i][j][k] + M_Concentration[3][6] * cIn6[i][j][k]));
						mfIn4 = (A_gConcentration[4][4] * (M_Concentration[4][0] * cIn0[i][j][k] + M_Concentration[4][1] * cIn1[i][j][k] + M_Concentration[4][2] * cIn2[i][j][k] + M_Concentration[4][3] * cIn3[i][j][k] + M_Concentration[4][4] * cIn4[i][j][k] + M_Concentration[4][5] * cIn5[i][j][k] + M_Concentration[4][6] * cIn6[i][j][k]));
						mfIn5 = (A_gConcentration[5][5] * (M_Concentration[5][0] * cIn0[i][j][k] + M_Concentration[5][1] * cIn1[i][j][k] + M_Concentration[5][2] * cIn2[i][j][k] + M_Concentration[5][3] * cIn3[i][j][k] + M_Concentration[5][4] * cIn4[i][j][k] + M_Concentration[5][5] * cIn5[i][j][k] + M_Concentration[5][6] * cIn6[i][j][k]));
						mfIn6 = (A_gConcentration[6][6] * (M_Concentration[6][0] * cIn0[i][j][k] + M_Concentration[6][1] * cIn1[i][j][k] + M_Concentration[6][2] * cIn2[i][j][k] + M_Concentration[6][3] * cIn3[i][j][k] + M_Concentration[6][4] * cIn4[i][j][k] + M_Concentration[6][5] * cIn5[i][j][k] + M_Concentration[6][6] * cIn6[i][j][k]));

						mfEq0 = A_gConcentration[0][0] * Concentration[i][j][k];
						mfEq1 = A_gConcentration[1][1] * Concentration[i][j][k] * (M_Concentration[1][0] * Eq0[i][j][k] + M_Concentration[1][1] * Eq1[i][j][k] + M_Concentration[1][2] * Eq2[i][j][k] + M_Concentration[1][3] * Eq3[i][j][k] + M_Concentration[1][4] * Eq4[i][j][k] + M_Concentration[1][5] * Eq5[i][j][k] + M_Concentration[1][6] * Eq6[i][j][k]);
						mfEq2 = A_gConcentration[2][2] * Concentration[i][j][k] * (M_Concentration[2][0] * Eq0[i][j][k] + M_Concentration[2][1] * Eq1[i][j][k] + M_Concentration[2][2] * Eq2[i][j][k] + M_Concentration[2][3] * Eq3[i][j][k] + M_Concentration[2][4] * Eq4[i][j][k] + M_Concentration[2][5] * Eq5[i][j][k] + M_Concentration[2][6] * Eq6[i][j][k]);
						mfEq3 = A_gConcentration[3][3] * Concentration[i][j][k] * (M_Concentration[3][0] * Eq0[i][j][k] + M_Concentration[3][1] * Eq1[i][j][k] + M_Concentration[3][2] * Eq2[i][j][k] + M_Concentration[3][3] * Eq3[i][j][k] + M_Concentration[3][4] * Eq4[i][j][k] + M_Concentration[3][5] * Eq5[i][j][k] + M_Concentration[3][6] * Eq6[i][j][k]);
						mfEq4 = A_gConcentration[4][4] * Concentration[i][j][k] * (M_Concentration[4][0] * Eq0[i][j][k] + M_Concentration[4][1] * Eq1[i][j][k] + M_Concentration[4][2] * Eq2[i][j][k] + M_Concentration[4][3] * Eq3[i][j][k] + M_Concentration[4][4] * Eq4[i][j][k] + M_Concentration[4][5] * Eq5[i][j][k] + M_Concentration[4][6] * Eq6[i][j][k]);
						mfEq5 = A_gConcentration[5][5] * Concentration[i][j][k] * (M_Concentration[5][0] * Eq0[i][j][k] + M_Concentration[5][1] * Eq1[i][j][k] + M_Concentration[5][2] * Eq2[i][j][k] + M_Concentration[5][3] * Eq3[i][j][k] + M_Concentration[5][4] * Eq4[i][j][k] + M_Concentration[5][5] * Eq5[i][j][k] + M_Concentration[5][6] * Eq6[i][j][k]);
						mfEq6 = A_gConcentration[6][6] * Concentration[i][j][k] * (M_Concentration[6][0] * Eq0[i][j][k] + M_Concentration[6][1] * Eq1[i][j][k] + M_Concentration[6][2] * Eq2[i][j][k] + M_Concentration[6][3] * Eq3[i][j][k] + M_Concentration[6][4] * Eq4[i][j][k] + M_Concentration[6][5] * Eq5[i][j][k] + M_Concentration[6][6] * Eq6[i][j][k]);

						Out0[i][j][k] = cIn0[i][j][k] - (vfConcentrationIn0 - vfConcentrationEq0);
						Out1[i][j][k] = cIn1[i][j][k] - (vfConcentrationIn1 - vfConcentrationEq1);
						Out2[i][j][k] = cIn2[i][j][k] - (vfConcentrationIn2 - vfConcentrationEq2);
						Out3[i][j][k] = cIn3[i][j][k] - (vfConcentrationIn3 - vfConcentrationEq3);
						Out4[i][j][k] = cIn4[i][j][k] - (vfConcentrationIn4 - vfConcentrationEq4);
						Out5[i][j][k] = cIn5[i][j][k] - (vfConcentrationIn5 - vfConcentrationEq5);
						Out6[i][j][k] = cIn6[i][j][k] - (vfConcentrationIn6 - vfConcentrationEq6);
					}
					else{
						Out0[i][j][k] = 0.0;
						Out1[i][j][k] = 0.0;
						Out2[i][j][k] = 0.0;
						Out3[i][j][k] = 0.0;
						Out4[i][j][k] = 0.0;
						Out5[i][j][k] = 0.0;
						Out6[i][j][k] = 0.0;
					}
				}
			}
		}
	}
}

void CollisionMRTConcentrationD3Q15(int Length, int Width, int Height, double gasCritical, double ***Ux, double ***Uy, double ***Uz, double alpha1, double alpha2){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k,mfIn0,mfIn1,mfIn2,mfIn3,mfIn4,mfIn5,mfIn6,mfIn7,mfIn8,mfIn9,mfIn10,mfIn11,mfIn12,mfIn13,mfIn14,mfEq0, mfEq1, mfEq2, mfEq3, mfEq4, mfEq5, mfEq6,mfEq7,mfEq8,mfEq9,mfEq10,mfEq11,mfEq12,mfEq13,mfEq14)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (j = 0; j<Width; j++){
				for (k = 0; k<Height; k++){
					if (domain[i][j][k] == 0 && rho[i][j][k]>gasCritical&&domainConnection[i][j][k] == 1){
						mfIn0 = (A_gConcentrationD3Q15[0][0] * (M_ConcentrationD3Q15[0][0] * cIn0[i][j][k] + M_ConcentrationD3Q15[0][1] * cIn1[i][j][k] + M_ConcentrationD3Q15[0][2] * cIn2[i][j][k] + M_ConcentrationD3Q15[0][3] * cIn3[i][j][k] + M_ConcentrationD3Q15[0][4] * cIn4[i][j][k] + M_ConcentrationD3Q15[0][5] * cIn5[i][j][k] + M_ConcentrationD3Q15[0][6] * cIn6[i][j][k] + M_ConcentrationD3Q15[0][7] * cIn7[i][j][k] + M_ConcentrationD3Q15[0][8] * cIn8[i][j][k]\
							+ M_ConcentrationD3Q15[0][9] * cIn9[i][j][k] + M_ConcentrationD3Q15[0][10] * cIn10[i][j][k] + M_ConcentrationD3Q15[0][11] * cIn11[i][j][k] + M_ConcentrationD3Q15[0][12] * cIn12[i][j][k] + M_ConcentrationD3Q15[0][13] * cIn13[i][j][k] + M_ConcentrationD3Q15[0][14] * cIn14[i][j][k]));
						mfIn1 = (A_gConcentrationD3Q15[1][1] * (M_ConcentrationD3Q15[1][0] * cIn0[i][j][k] + M_ConcentrationD3Q15[1][1] * cIn1[i][j][k] + M_ConcentrationD3Q15[1][2] * cIn2[i][j][k] + M_ConcentrationD3Q15[1][3] * cIn3[i][j][k] + M_ConcentrationD3Q15[1][4] * cIn4[i][j][k] + M_ConcentrationD3Q15[1][5] * cIn5[i][j][k] + M_ConcentrationD3Q15[1][6] * cIn6[i][j][k] + M_ConcentrationD3Q15[1][7] * cIn7[i][j][k] + M_ConcentrationD3Q15[1][8] * cIn8[i][j][k]\
							+ M_ConcentrationD3Q15[1][9] * cIn9[i][j][k] + M_ConcentrationD3Q15[1][10] * cIn10[i][j][k] + M_ConcentrationD3Q15[1][11] * cIn11[i][j][k] + M_ConcentrationD3Q15[1][12] * cIn12[i][j][k] + M_ConcentrationD3Q15[1][13] * cIn13[i][j][k] + M_ConcentrationD3Q15[1][14] * cIn14[i][j][k]));
						mfIn2 = (A_gConcentrationD3Q15[2][2] * (M_ConcentrationD3Q15[2][0] * cIn0[i][j][k] + M_ConcentrationD3Q15[2][1] * cIn1[i][j][k] + M_ConcentrationD3Q15[2][2] * cIn2[i][j][k] + M_ConcentrationD3Q15[2][3] * cIn3[i][j][k] + M_ConcentrationD3Q15[2][4] * cIn4[i][j][k] + M_ConcentrationD3Q15[2][5] * cIn5[i][j][k] + M_ConcentrationD3Q15[2][6] * cIn6[i][j][k] + M_ConcentrationD3Q15[2][7] * cIn7[i][j][k] + M_ConcentrationD3Q15[2][8] * cIn8[i][j][k]\
							+ M_ConcentrationD3Q15[2][9] * cIn9[i][j][k] + M_ConcentrationD3Q15[2][10] * cIn10[i][j][k] + M_ConcentrationD3Q15[2][11] * cIn11[i][j][k] + M_ConcentrationD3Q15[2][12] * cIn12[i][j][k] + M_ConcentrationD3Q15[2][13] * cIn13[i][j][k] + M_ConcentrationD3Q15[2][14] * cIn14[i][j][k]));
						mfIn3 = (A_gConcentrationD3Q15[3][3] * (M_ConcentrationD3Q15[3][0] * cIn0[i][j][k] + M_ConcentrationD3Q15[3][1] * cIn1[i][j][k] + M_ConcentrationD3Q15[3][2] * cIn2[i][j][k] + M_ConcentrationD3Q15[3][3] * cIn3[i][j][k] + M_ConcentrationD3Q15[3][4] * cIn4[i][j][k] + M_ConcentrationD3Q15[3][5] * cIn5[i][j][k] + M_ConcentrationD3Q15[3][6] * cIn6[i][j][k] + M_ConcentrationD3Q15[3][7] * cIn7[i][j][k] + M_ConcentrationD3Q15[3][8] * cIn8[i][j][k]\
							+ M_ConcentrationD3Q15[3][9] * cIn9[i][j][k] + M_ConcentrationD3Q15[3][10] * cIn10[i][j][k] + M_ConcentrationD3Q15[3][11] * cIn11[i][j][k] + M_ConcentrationD3Q15[3][12] * cIn12[i][j][k] + M_ConcentrationD3Q15[3][13] * cIn13[i][j][k] + M_ConcentrationD3Q15[3][14] * cIn14[i][j][k]) +A_gConcentrationD3Q15[3][4] * (M_ConcentrationD3Q15[4][0] * cIn0[i][j][k] + M_ConcentrationD3Q15[4][1] * cIn1[i][j][k] + M_ConcentrationD3Q15[4][2] * cIn2[i][j][k] + M_ConcentrationD3Q15[4][3] * cIn3[i][j][k] + M_ConcentrationD3Q15[4][4] * cIn4[i][j][k] + M_ConcentrationD3Q15[4][5] * cIn5[i][j][k] + M_ConcentrationD3Q15[4][6] * cIn6[i][j][k] + M_ConcentrationD3Q15[4][7] * cIn7[i][j][k] + M_ConcentrationD3Q15[4][8] * cIn8[i][j][k]\
							+ M_ConcentrationD3Q15[4][9] * cIn9[i][j][k] + M_ConcentrationD3Q15[4][10] * cIn10[i][j][k] + M_ConcentrationD3Q15[4][11] * cIn11[i][j][k] + M_ConcentrationD3Q15[4][12] * cIn12[i][j][k] + M_ConcentrationD3Q15[4][13] * cIn13[i][j][k] + M_ConcentrationD3Q15[4][14] * cIn14[i][j][k]));
						mfIn4 = (A_gConcentrationD3Q15[4][4] * (M_ConcentrationD3Q15[4][0] * cIn0[i][j][k] + M_ConcentrationD3Q15[4][1] * cIn1[i][j][k] + M_ConcentrationD3Q15[4][2] * cIn2[i][j][k] + M_ConcentrationD3Q15[4][3] * cIn3[i][j][k] + M_ConcentrationD3Q15[4][4] * cIn4[i][j][k] + M_ConcentrationD3Q15[4][5] * cIn5[i][j][k] + M_ConcentrationD3Q15[4][6] * cIn6[i][j][k] + M_ConcentrationD3Q15[4][7] * cIn7[i][j][k] + M_ConcentrationD3Q15[4][8] * cIn8[i][j][k]\
							+ M_ConcentrationD3Q15[4][9] * cIn9[i][j][k] + M_ConcentrationD3Q15[4][10] * cIn10[i][j][k] + M_ConcentrationD3Q15[4][11] * cIn11[i][j][k] + M_ConcentrationD3Q15[4][12] * cIn12[i][j][k] + M_ConcentrationD3Q15[4][13] * cIn13[i][j][k] + M_ConcentrationD3Q15[4][14] * cIn14[i][j][k]));
						mfIn5 = (A_gConcentrationD3Q15[5][5] * (M_ConcentrationD3Q15[5][0] * cIn0[i][j][k] + M_ConcentrationD3Q15[5][1] * cIn1[i][j][k] + M_ConcentrationD3Q15[5][2] * cIn2[i][j][k] + M_ConcentrationD3Q15[5][3] * cIn3[i][j][k] + M_ConcentrationD3Q15[5][4] * cIn4[i][j][k] + M_ConcentrationD3Q15[5][5] * cIn5[i][j][k] + M_ConcentrationD3Q15[5][6] * cIn6[i][j][k] + M_ConcentrationD3Q15[5][7] * cIn7[i][j][k] + M_ConcentrationD3Q15[5][8] * cIn8[i][j][k]\
							+ M_ConcentrationD3Q15[5][9] * cIn9[i][j][k] + M_ConcentrationD3Q15[5][10] * cIn10[i][j][k] + M_ConcentrationD3Q15[5][11] * cIn11[i][j][k] + M_ConcentrationD3Q15[5][12] * cIn12[i][j][k] + M_ConcentrationD3Q15[5][13] * cIn13[i][j][k] + M_ConcentrationD3Q15[5][14] * cIn14[i][j][k]) +A_gConcentrationD3Q15[5][6] * (M_ConcentrationD3Q15[6][0] * cIn0[i][j][k] + M_ConcentrationD3Q15[6][1] * cIn1[i][j][k] + M_ConcentrationD3Q15[6][2] * cIn2[i][j][k] + M_ConcentrationD3Q15[6][3] * cIn3[i][j][k] + M_ConcentrationD3Q15[6][4] * cIn4[i][j][k] + M_ConcentrationD3Q15[6][5] * cIn5[i][j][k] + M_ConcentrationD3Q15[6][6] * cIn6[i][j][k] + M_ConcentrationD3Q15[6][7] * cIn7[i][j][k] + M_ConcentrationD3Q15[6][8] * cIn8[i][j][k]\
							+ M_ConcentrationD3Q15[6][9] * cIn9[i][j][k] + M_ConcentrationD3Q15[6][10] * cIn10[i][j][k] + M_ConcentrationD3Q15[6][11] * cIn11[i][j][k] + M_ConcentrationD3Q15[6][12] * cIn12[i][j][k] + M_ConcentrationD3Q15[6][13] * cIn13[i][j][k] + M_ConcentrationD3Q15[6][14] * cIn14[i][j][k]));
						mfIn6 = (A_gConcentrationD3Q15[6][6] * (M_ConcentrationD3Q15[6][0] * cIn0[i][j][k] + M_ConcentrationD3Q15[6][1] * cIn1[i][j][k] + M_ConcentrationD3Q15[6][2] * cIn2[i][j][k] + M_ConcentrationD3Q15[6][3] * cIn3[i][j][k] + M_ConcentrationD3Q15[6][4] * cIn4[i][j][k] + M_ConcentrationD3Q15[6][5] * cIn5[i][j][k] + M_ConcentrationD3Q15[6][6] * cIn6[i][j][k] + M_ConcentrationD3Q15[6][7] * cIn7[i][j][k] + M_ConcentrationD3Q15[6][8] * cIn8[i][j][k]\
							+ M_ConcentrationD3Q15[6][9] * cIn9[i][j][k] + M_ConcentrationD3Q15[6][10] * cIn10[i][j][k] + M_ConcentrationD3Q15[6][11] * cIn11[i][j][k] + M_ConcentrationD3Q15[6][12] * cIn12[i][j][k] + M_ConcentrationD3Q15[6][13] * cIn13[i][j][k] + M_ConcentrationD3Q15[6][14] * cIn14[i][j][k]));
						mfIn7 = (A_gConcentrationD3Q15[7][7] * (M_ConcentrationD3Q15[7][0] * cIn0[i][j][k] + M_ConcentrationD3Q15[7][1] * cIn1[i][j][k] + M_ConcentrationD3Q15[7][2] * cIn2[i][j][k] + M_ConcentrationD3Q15[7][3] * cIn3[i][j][k] + M_ConcentrationD3Q15[7][4] * cIn4[i][j][k] + M_ConcentrationD3Q15[7][5] * cIn5[i][j][k] + M_ConcentrationD3Q15[7][6] * cIn6[i][j][k] + M_ConcentrationD3Q15[7][7] * cIn7[i][j][k] + M_ConcentrationD3Q15[7][8] * cIn8[i][j][k]\
							+ M_ConcentrationD3Q15[7][9] * cIn9[i][j][k] + M_ConcentrationD3Q15[7][10] * cIn10[i][j][k] + M_ConcentrationD3Q15[7][11] * cIn11[i][j][k] + M_ConcentrationD3Q15[7][12] * cIn12[i][j][k] + M_ConcentrationD3Q15[7][13] * cIn13[i][j][k] + M_ConcentrationD3Q15[7][14] * cIn14[i][j][k]) +A_gConcentrationD3Q15[7][8] * (M_ConcentrationD3Q15[8][0] * cIn0[i][j][k] + M_ConcentrationD3Q15[8][1] * cIn1[i][j][k] + M_ConcentrationD3Q15[8][2] * cIn2[i][j][k] + M_ConcentrationD3Q15[8][3] * cIn3[i][j][k] + M_ConcentrationD3Q15[8][4] * cIn4[i][j][k] + M_ConcentrationD3Q15[8][5] * cIn5[i][j][k] + M_ConcentrationD3Q15[8][6] * cIn6[i][j][k] + M_ConcentrationD3Q15[8][7] * cIn7[i][j][k] + M_ConcentrationD3Q15[8][8] * cIn8[i][j][k]\
							+ M_ConcentrationD3Q15[8][9] * cIn9[i][j][k] + M_ConcentrationD3Q15[8][10] * cIn10[i][j][k] + M_ConcentrationD3Q15[8][11] * cIn11[i][j][k] + M_ConcentrationD3Q15[8][12] * cIn12[i][j][k] + M_ConcentrationD3Q15[8][13] * cIn13[i][j][k] + M_ConcentrationD3Q15[8][14] * cIn14[i][j][k]));
						mfIn8 = (A_gConcentrationD3Q15[8][8] * (M_ConcentrationD3Q15[8][0] * cIn0[i][j][k] + M_ConcentrationD3Q15[8][1] * cIn1[i][j][k] + M_ConcentrationD3Q15[8][2] * cIn2[i][j][k] + M_ConcentrationD3Q15[8][3] * cIn3[i][j][k] + M_ConcentrationD3Q15[8][4] * cIn4[i][j][k] + M_ConcentrationD3Q15[8][5] * cIn5[i][j][k] + M_ConcentrationD3Q15[8][6] * cIn6[i][j][k] + M_ConcentrationD3Q15[8][7] * cIn7[i][j][k] + M_ConcentrationD3Q15[8][8] * cIn8[i][j][k]\
							+ M_ConcentrationD3Q15[8][9] * cIn9[i][j][k] + M_ConcentrationD3Q15[8][10] * cIn10[i][j][k] + M_ConcentrationD3Q15[8][11] * cIn11[i][j][k] + M_ConcentrationD3Q15[8][12] * cIn12[i][j][k] + M_ConcentrationD3Q15[8][13] * cIn13[i][j][k] + M_ConcentrationD3Q15[8][14] * cIn14[i][j][k]));
						mfIn9 = (A_gConcentrationD3Q15[9][9] * (M_ConcentrationD3Q15[9][0] * cIn0[i][j][k] + M_ConcentrationD3Q15[9][1] * cIn1[i][j][k] + M_ConcentrationD3Q15[9][2] * cIn2[i][j][k] + M_ConcentrationD3Q15[9][3] * cIn3[i][j][k] + M_ConcentrationD3Q15[9][4] * cIn4[i][j][k] + M_ConcentrationD3Q15[9][5] * cIn5[i][j][k] + M_ConcentrationD3Q15[9][6] * cIn6[i][j][k] + M_ConcentrationD3Q15[9][7] * cIn7[i][j][k] + M_ConcentrationD3Q15[9][8] * cIn8[i][j][k]\
							+ M_ConcentrationD3Q15[9][9] * cIn9[i][j][k] + M_ConcentrationD3Q15[9][10] * cIn10[i][j][k] + M_ConcentrationD3Q15[9][11] * cIn11[i][j][k] + M_ConcentrationD3Q15[9][12] * cIn12[i][j][k] + M_ConcentrationD3Q15[9][13] * cIn13[i][j][k] + M_ConcentrationD3Q15[9][14] * cIn14[i][j][k]));
						mfIn10 = (A_gConcentrationD3Q15[10][10] * (M_ConcentrationD3Q15[10][0] * cIn0[i][j][k] + M_ConcentrationD3Q15[10][1] * cIn1[i][j][k] + M_ConcentrationD3Q15[10][2] * cIn2[i][j][k] + M_ConcentrationD3Q15[10][3] * cIn3[i][j][k] + M_ConcentrationD3Q15[10][4] * cIn4[i][j][k] + M_ConcentrationD3Q15[10][5] * cIn5[i][j][k] + M_ConcentrationD3Q15[10][6] * cIn6[i][j][k] + M_ConcentrationD3Q15[10][7] * cIn7[i][j][k] + M_ConcentrationD3Q15[10][8] * cIn8[i][j][k]\
							+ M_ConcentrationD3Q15[10][9] * cIn9[i][j][k] + M_ConcentrationD3Q15[10][10] * cIn10[i][j][k] + M_ConcentrationD3Q15[10][11] * cIn11[i][j][k] + M_ConcentrationD3Q15[10][12] * cIn12[i][j][k] + M_ConcentrationD3Q15[10][13] * cIn13[i][j][k] + M_ConcentrationD3Q15[10][14] * cIn14[i][j][k]));
						mfIn11 = (A_gConcentrationD3Q15[11][11] * (M_ConcentrationD3Q15[11][0] * cIn0[i][j][k] + M_ConcentrationD3Q15[11][1] * cIn1[i][j][k] + M_ConcentrationD3Q15[11][2] * cIn2[i][j][k] + M_ConcentrationD3Q15[11][3] * cIn3[i][j][k] + M_ConcentrationD3Q15[11][4] * cIn4[i][j][k] + M_ConcentrationD3Q15[11][5] * cIn5[i][j][k] + M_ConcentrationD3Q15[11][6] * cIn6[i][j][k] + M_ConcentrationD3Q15[11][7] * cIn7[i][j][k] + M_ConcentrationD3Q15[11][8] * cIn8[i][j][k]\
							+ M_ConcentrationD3Q15[11][9] * cIn9[i][j][k] + M_ConcentrationD3Q15[11][10] * cIn10[i][j][k] + M_ConcentrationD3Q15[11][11] * cIn11[i][j][k] + M_ConcentrationD3Q15[11][12] * cIn12[i][j][k] + M_ConcentrationD3Q15[11][13] * cIn13[i][j][k] + M_ConcentrationD3Q15[11][14] * cIn14[i][j][k]));
						mfIn12 = (A_gConcentrationD3Q15[12][12] * (M_ConcentrationD3Q15[12][0] * cIn0[i][j][k] + M_ConcentrationD3Q15[12][1] * cIn1[i][j][k] + M_ConcentrationD3Q15[12][2] * cIn2[i][j][k] + M_ConcentrationD3Q15[12][3] * cIn3[i][j][k] + M_ConcentrationD3Q15[12][4] * cIn4[i][j][k] + M_ConcentrationD3Q15[12][5] * cIn5[i][j][k] + M_ConcentrationD3Q15[12][6] * cIn6[i][j][k] + M_ConcentrationD3Q15[12][7] * cIn7[i][j][k] + M_ConcentrationD3Q15[12][8] * cIn8[i][j][k]\
							+ M_ConcentrationD3Q15[12][9] * cIn9[i][j][k] + M_ConcentrationD3Q15[12][10] * cIn10[i][j][k] + M_ConcentrationD3Q15[12][11] * cIn11[i][j][k] + M_ConcentrationD3Q15[12][12] * cIn12[i][j][k] + M_ConcentrationD3Q15[12][13] * cIn13[i][j][k] + M_ConcentrationD3Q15[12][14] * cIn14[i][j][k]));
						mfIn13 = (A_gConcentrationD3Q15[13][13] * (M_ConcentrationD3Q15[13][0] * cIn0[i][j][k] + M_ConcentrationD3Q15[13][1] * cIn1[i][j][k] + M_ConcentrationD3Q15[13][2] * cIn2[i][j][k] + M_ConcentrationD3Q15[13][3] * cIn3[i][j][k] + M_ConcentrationD3Q15[13][4] * cIn4[i][j][k] + M_ConcentrationD3Q15[13][5] * cIn5[i][j][k] + M_ConcentrationD3Q15[13][6] * cIn6[i][j][k] + M_ConcentrationD3Q15[13][7] * cIn7[i][j][k] + M_ConcentrationD3Q15[13][8] * cIn8[i][j][k]\
							+ M_ConcentrationD3Q15[13][9] * cIn9[i][j][k] + M_ConcentrationD3Q15[13][10] * cIn10[i][j][k] + M_ConcentrationD3Q15[13][11] * cIn11[i][j][k] + M_ConcentrationD3Q15[13][12] * cIn12[i][j][k] + M_ConcentrationD3Q15[13][13] * cIn13[i][j][k] + M_ConcentrationD3Q15[13][14] * cIn14[i][j][k]));
						mfIn14 = (A_gConcentrationD3Q15[14][14] * (M_ConcentrationD3Q15[14][0] * cIn0[i][j][k] + M_ConcentrationD3Q15[14][1] * cIn1[i][j][k] + M_ConcentrationD3Q15[14][2] * cIn2[i][j][k] + M_ConcentrationD3Q15[14][3] * cIn3[i][j][k] + M_ConcentrationD3Q15[14][4] * cIn4[i][j][k] + M_ConcentrationD3Q15[14][5] * cIn5[i][j][k] + M_ConcentrationD3Q15[14][6] * cIn6[i][j][k] + M_ConcentrationD3Q15[14][7] * cIn7[i][j][k] + M_ConcentrationD3Q15[14][8] * cIn8[i][j][k]\
							+ M_ConcentrationD3Q15[14][9] * cIn9[i][j][k] + M_ConcentrationD3Q15[14][10] * cIn10[i][j][k] + M_ConcentrationD3Q15[14][11] * cIn11[i][j][k] + M_ConcentrationD3Q15[14][12] * cIn12[i][j][k] + M_ConcentrationD3Q15[14][13] * cIn13[i][j][k] + M_ConcentrationD3Q15[14][14] * cIn14[i][j][k]));

						mfEq0 = A_gConcentrationD3Q15[0][0] * Concentration[i][j][k];
						mfEq1 = A_gConcentrationD3Q15[1][1] * 5.0/6.0*alpha1*Concentration[i][j][k];
						mfEq2 = A_gConcentrationD3Q15[2][2] * 1.0/24.0*alpha2*Concentration[i][j][k];
						mfEq3 = A_gConcentrationD3Q15[3][3] * 1.0*Concentration[i][j][k] * Ux[i][j][k]+A_gConcentrationD3Q15[3][4] * (-1.0*Concentration[i][j][k] * Ux[i][j][k]);
						mfEq4 = A_gConcentrationD3Q15[4][4] * (-1.0*Concentration[i][j][k]*Ux[i][j][k]);
						mfEq5 = A_gConcentrationD3Q15[5][5] * 1.0*Concentration[i][j][k] * Uy[i][j][k]+A_gConcentrationD3Q15[5][6] * (-1.0*Concentration[i][j][k]*Uy[i][j][k]);
						mfEq6 = A_gConcentrationD3Q15[6][6] * (-1.0*Concentration[i][j][k]*Uy[i][j][k]);
						mfEq7 = A_gConcentrationD3Q15[7][7] * 1.0*Concentration[i][j][k] * Uz[i][j][k]+A_gConcentrationD3Q15[7][8] * (-1.0*Concentration[i][j][k]*Uz[i][j][k]);
						mfEq8 = A_gConcentrationD3Q15[8][8] * (-1.0*Concentration[i][j][k]*Uz[i][j][k]);
						mfEq9 = 0.0;
						mfEq10 = 0.0;
						mfEq11 = 0.0;
						mfEq12 = 0.0;
						mfEq13 = 0.0;
						mfEq14 = 0.0;

						Out0[i][j][k] = cIn0[i][j][k] -(vfConcentrationD3Q15In0 - vfConcentrationD3Q15Eq0);
						Out1[i][j][k] = cIn1[i][j][k] -(vfConcentrationD3Q15In1 - vfConcentrationD3Q15Eq1);
						Out2[i][j][k] = cIn2[i][j][k] -(vfConcentrationD3Q15In2 - vfConcentrationD3Q15Eq2);
						Out3[i][j][k] = cIn3[i][j][k] -(vfConcentrationD3Q15In3 - vfConcentrationD3Q15Eq3);
						Out4[i][j][k] = cIn4[i][j][k] -(vfConcentrationD3Q15In4 - vfConcentrationD3Q15Eq4);
						Out5[i][j][k] = cIn5[i][j][k] -(vfConcentrationD3Q15In5 - vfConcentrationD3Q15Eq5);
						Out6[i][j][k] = cIn6[i][j][k] -(vfConcentrationD3Q15In6 - vfConcentrationD3Q15Eq6);
						Out7[i][j][k] = cIn7[i][j][k] -(vfConcentrationD3Q15In7 - vfConcentrationD3Q15Eq7);
						Out8[i][j][k] = cIn8[i][j][k] -(vfConcentrationD3Q15In8 - vfConcentrationD3Q15Eq8);
						Out9[i][j][k] = cIn9[i][j][k] -(vfConcentrationD3Q15In9 - vfConcentrationD3Q15Eq9);
						Out10[i][j][k] = cIn10[i][j][k] -(vfConcentrationD3Q15In10 - vfConcentrationD3Q15Eq10);
						Out11[i][j][k] = cIn11[i][j][k] -(vfConcentrationD3Q15In11 - vfConcentrationD3Q15Eq11);
						Out12[i][j][k] = cIn12[i][j][k] -(vfConcentrationD3Q15In12 - vfConcentrationD3Q15Eq12);
						Out13[i][j][k] = cIn13[i][j][k] -(vfConcentrationD3Q15In13 - vfConcentrationD3Q15Eq13);
						Out14[i][j][k] = cIn14[i][j][k] -(vfConcentrationD3Q15In14 - vfConcentrationD3Q15Eq14);
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
					}
				}
			}
		}
	}
}

void electrolyteConcentrationRelaxationTimeCalculation(int Length, int Width, int Height, double constant, double gasCritical,double constantMol,double constant1){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (j = 0; j<Width; j++){
				for (k = 0; k<Height; k++){
					if (domain[i][j][k] == 0 && rho[i][j][k]>gasCritical&&domainConnection[i][j][k] == 1){
//						electrolyteConcentrationRelaxationTime[i][j][k] = 1.0 / (pow(10.0, (-4.43 - (54.0 / (constant - (229.0 + 5.0*Concentration[i][j][k]*constantMol))) - 0.22*Concentration[i][j][k]*constantMol))*constant1 * 4.0 + 0.5);
						electrolyteConcentrationRelaxationTime[i][j][k] = 1.0/(0.5105*4+0.5);
					}
				}
			}
		}
	}
}

void CollisionSRTConcentrationD3Q7(int Length, int Width, int Height, double gasCritical){
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
						Out0[i][j][k] = cIn0[i][j][k] - electrolyteConcentrationRelaxationTime[i][j][k] * (cIn0[i][j][k] - tNSD3Q7[0] * Concentration[i][j][k]);// +diffusionTermConcentration[i][j][k]);
						Out1[i][j][k] = cIn1[i][j][k] - electrolyteConcentrationRelaxationTime[i][j][k] * (cIn1[i][j][k] - tNSD3Q7[1] * Concentration[i][j][k]);// +diffusionTermConcentration[i][j][k]);
						Out2[i][j][k] = cIn2[i][j][k] - electrolyteConcentrationRelaxationTime[i][j][k] * (cIn2[i][j][k] - tNSD3Q7[2] * Concentration[i][j][k]);// +diffusionTermConcentration[i][j][k]);
						Out3[i][j][k] = cIn3[i][j][k] - electrolyteConcentrationRelaxationTime[i][j][k] * (cIn3[i][j][k] - tNSD3Q7[3] * Concentration[i][j][k]);// +diffusionTermConcentration[i][j][k]);
						Out4[i][j][k] = cIn4[i][j][k] - electrolyteConcentrationRelaxationTime[i][j][k] * (cIn4[i][j][k] - tNSD3Q7[4] * Concentration[i][j][k]);// +diffusionTermConcentration[i][j][k]);
						Out5[i][j][k] = cIn5[i][j][k] - electrolyteConcentrationRelaxationTime[i][j][k] * (cIn5[i][j][k] - tNSD3Q7[5] * Concentration[i][j][k]);// +diffusionTermConcentration[i][j][k]);
						Out6[i][j][k] = cIn6[i][j][k] - electrolyteConcentrationRelaxationTime[i][j][k] * (cIn6[i][j][k] - tNSD3Q7[6] * Concentration[i][j][k]);// +diffusionTermConcentration[i][j][k]);
					}
				}
			}
		}
	}
}

void CollisionSRTConcentrationFirst(int Length,int Width,int Height,double CrelaxationTime,double gasCritical){
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
						Out0[i][j][k]=cIn0First[i][j][k]-CrelaxationTime*(cIn0First[i][j][k]-Eq0[i][j][k]*ConcentrationFirst[i][j][k]);
						Out1[i][j][k]=cIn1First[i][j][k]-CrelaxationTime*(cIn1First[i][j][k]-Eq1[i][j][k]*ConcentrationFirst[i][j][k]);
						Out2[i][j][k]=cIn2First[i][j][k]-CrelaxationTime*(cIn2First[i][j][k]-Eq2[i][j][k]*ConcentrationFirst[i][j][k]);
						Out3[i][j][k]=cIn3First[i][j][k]-CrelaxationTime*(cIn3First[i][j][k]-Eq3[i][j][k]*ConcentrationFirst[i][j][k]);
						Out4[i][j][k]=cIn4First[i][j][k]-CrelaxationTime*(cIn4First[i][j][k]-Eq4[i][j][k]*ConcentrationFirst[i][j][k]);
						Out5[i][j][k]=cIn5First[i][j][k]-CrelaxationTime*(cIn5First[i][j][k]-Eq5[i][j][k]*ConcentrationFirst[i][j][k]);
						Out6[i][j][k]=cIn6First[i][j][k]-CrelaxationTime*(cIn6First[i][j][k]-Eq6[i][j][k]*ConcentrationFirst[i][j][k]);
						Out7[i][j][k]=cIn7First[i][j][k]-CrelaxationTime*(cIn7First[i][j][k]-Eq7[i][j][k]*ConcentrationFirst[i][j][k]);
						Out8[i][j][k]=cIn8First[i][j][k]-CrelaxationTime*(cIn8First[i][j][k]-Eq8[i][j][k]*ConcentrationFirst[i][j][k]);
						Out9[i][j][k]=cIn9First[i][j][k]-CrelaxationTime*(cIn9First[i][j][k]-Eq9[i][j][k]*ConcentrationFirst[i][j][k]);
						Out10[i][j][k]=cIn10First[i][j][k]-CrelaxationTime*(cIn10First[i][j][k]-Eq10[i][j][k]*ConcentrationFirst[i][j][k]);
						Out11[i][j][k]=cIn11First[i][j][k]-CrelaxationTime*(cIn11First[i][j][k]-Eq11[i][j][k]*ConcentrationFirst[i][j][k]);
						Out12[i][j][k]=cIn12First[i][j][k]-CrelaxationTime*(cIn12First[i][j][k]-Eq12[i][j][k]*ConcentrationFirst[i][j][k]);
						Out13[i][j][k]=cIn13First[i][j][k]-CrelaxationTime*(cIn13First[i][j][k]-Eq13[i][j][k]*ConcentrationFirst[i][j][k]);
						Out14[i][j][k]=cIn14First[i][j][k]-CrelaxationTime*(cIn14First[i][j][k]-Eq14[i][j][k]*ConcentrationFirst[i][j][k]);
						Out15[i][j][k]=cIn15First[i][j][k]-CrelaxationTime*(cIn15First[i][j][k]-Eq15[i][j][k]*ConcentrationFirst[i][j][k]);
						Out16[i][j][k]=cIn16First[i][j][k]-CrelaxationTime*(cIn16First[i][j][k]-Eq16[i][j][k]*ConcentrationFirst[i][j][k]);
						Out17[i][j][k]=cIn17First[i][j][k]-CrelaxationTime*(cIn17First[i][j][k]-Eq17[i][j][k]*ConcentrationFirst[i][j][k]);
						Out18[i][j][k]=cIn18First[i][j][k]-CrelaxationTime*(cIn18First[i][j][k]-Eq18[i][j][k]*ConcentrationFirst[i][j][k]);
					}
				}
			}
		}
	}
}

void CollisionSRTConcentrationElectrodeD3Q7(int Length, int Width, int Height, double CrelaxationTime, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
					if (domain[i][j][k] == 1 && electrodeConnection[i][j][k] == 1){
						Out0[i][j][k] = cIn0First[i][j][k] - CrelaxationTime*(cIn0First[i][j][k] - tNSD3Q7[0] * ConcentrationFirst[i][j][k]);
						Out1[i][j][k] = cIn1First[i][j][k] - CrelaxationTime*(cIn1First[i][j][k] - tNSD3Q7[1] * ConcentrationFirst[i][j][k]);
						Out2[i][j][k] = cIn2First[i][j][k] - CrelaxationTime*(cIn2First[i][j][k] - tNSD3Q7[2] * ConcentrationFirst[i][j][k]);
						Out3[i][j][k] = cIn3First[i][j][k] - CrelaxationTime*(cIn3First[i][j][k] - tNSD3Q7[3] * ConcentrationFirst[i][j][k]);
						Out4[i][j][k] = cIn4First[i][j][k] - CrelaxationTime*(cIn4First[i][j][k] - tNSD3Q7[4] * ConcentrationFirst[i][j][k]);
						Out5[i][j][k] = cIn5First[i][j][k] - CrelaxationTime*(cIn5First[i][j][k] - tNSD3Q7[5] * ConcentrationFirst[i][j][k]);
						Out6[i][j][k] = cIn6First[i][j][k] - CrelaxationTime*(cIn6First[i][j][k] - tNSD3Q7[6] * ConcentrationFirst[i][j][k]);
					}
				}
			}
		}
	}
}

void CollisionMRTConcentrationElectrodeD3Q7(int Length, int Width, int Height){
	int i;
	int j;
	int k;
	double mfIn0;
	double mfIn1;
	double mfIn2;
	double mfIn3;
	double mfIn4;
	double mfIn5;
	double mfIn6;
	double mfEq0;
	double mfEq1;
	double mfEq2;
	double mfEq3;
	double mfEq4;
	double mfEq5;
	double mfEq6;
#pragma omp parallel private(i,j,k,mfIn0,mfIn1,mfIn2,mfIn3,mfIn4,mfIn5,mfIn6,mfEq0, mfEq1, mfEq2, mfEq3, mfEq4, mfEq5, mfEq6)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (j = 0; j<Width; j++){
				for (k = 0; k<Height; k++){
					if (domain[i][j][k] == 1 && electrodeConnection[i][j][k] == 1){
						mfIn0 = (A_gConcentration[0][0] * (M_Concentration[0][0] * cIn0First[i][j][k] + M_Concentration[0][1] * cIn1First[i][j][k] + M_Concentration[0][2] * cIn2First[i][j][k] + M_Concentration[0][3] * cIn3First[i][j][k] + M_Concentration[0][4] * cIn4First[i][j][k] + M_Concentration[0][5] * cIn5First[i][j][k] + M_Concentration[0][6] * cIn6First[i][j][k]));
						mfIn1 = (A_gConcentration[1][1] * (M_Concentration[1][0] * cIn0First[i][j][k] + M_Concentration[1][1] * cIn1First[i][j][k] + M_Concentration[1][2] * cIn2First[i][j][k] + M_Concentration[1][3] * cIn3First[i][j][k] + M_Concentration[1][4] * cIn4First[i][j][k] + M_Concentration[1][5] * cIn5First[i][j][k] + M_Concentration[1][6] * cIn6First[i][j][k]));
						mfIn2 = (A_gConcentration[2][2] * (M_Concentration[2][0] * cIn0First[i][j][k] + M_Concentration[2][1] * cIn1First[i][j][k] + M_Concentration[2][2] * cIn2First[i][j][k] + M_Concentration[2][3] * cIn3First[i][j][k] + M_Concentration[2][4] * cIn4First[i][j][k] + M_Concentration[2][5] * cIn5First[i][j][k] + M_Concentration[2][6] * cIn6First[i][j][k]));
						mfIn3 = (A_gConcentration[3][3] * (M_Concentration[3][0] * cIn0First[i][j][k] + M_Concentration[3][1] * cIn1First[i][j][k] + M_Concentration[3][2] * cIn2First[i][j][k] + M_Concentration[3][3] * cIn3First[i][j][k] + M_Concentration[3][4] * cIn4First[i][j][k] + M_Concentration[3][5] * cIn5First[i][j][k] + M_Concentration[3][6] * cIn6First[i][j][k]));
						mfIn4 = (A_gConcentration[4][4] * (M_Concentration[4][0] * cIn0First[i][j][k] + M_Concentration[4][1] * cIn1First[i][j][k] + M_Concentration[4][2] * cIn2First[i][j][k] + M_Concentration[4][3] * cIn3First[i][j][k] + M_Concentration[4][4] * cIn4First[i][j][k] + M_Concentration[4][5] * cIn5First[i][j][k] + M_Concentration[4][6] * cIn6First[i][j][k]));
						mfIn5 = (A_gConcentration[5][5] * (M_Concentration[5][0] * cIn0First[i][j][k] + M_Concentration[5][1] * cIn1First[i][j][k] + M_Concentration[5][2] * cIn2First[i][j][k] + M_Concentration[5][3] * cIn3First[i][j][k] + M_Concentration[5][4] * cIn4First[i][j][k] + M_Concentration[5][5] * cIn5First[i][j][k] + M_Concentration[5][6] * cIn6First[i][j][k]));
						mfIn6 = (A_gConcentration[6][6] * (M_Concentration[6][0] * cIn0First[i][j][k] + M_Concentration[6][1] * cIn1First[i][j][k] + M_Concentration[6][2] * cIn2First[i][j][k] + M_Concentration[6][3] * cIn3First[i][j][k] + M_Concentration[6][4] * cIn4First[i][j][k] + M_Concentration[6][5] * cIn5First[i][j][k] + M_Concentration[6][6] * cIn6First[i][j][k]));

						mfEq0 = A_gConcentration[0][0] * ConcentrationFirst[i][j][k];
						mfEq1 = A_gConcentration[1][1] * ConcentrationFirst[i][j][k] * (M_Concentration[1][0] * tNSD3Q7[0] + M_Concentration[1][1] * tNSD3Q7[1] + M_Concentration[1][2] * tNSD3Q7[2] + M_Concentration[1][3] * tNSD3Q7[3] + M_Concentration[1][4] * tNSD3Q7[4] + M_Concentration[1][5] * tNSD3Q7[5] + M_Concentration[1][6] * tNSD3Q7[6]);
						mfEq2 = A_gConcentration[2][2] * ConcentrationFirst[i][j][k] * (M_Concentration[2][0] * tNSD3Q7[0] + M_Concentration[2][1] * tNSD3Q7[1] + M_Concentration[2][2] * tNSD3Q7[2] + M_Concentration[2][3] * tNSD3Q7[3] + M_Concentration[2][4] * tNSD3Q7[4] + M_Concentration[2][5] * tNSD3Q7[5] + M_Concentration[2][6] * tNSD3Q7[6]);
						mfEq3 = A_gConcentration[3][3] * ConcentrationFirst[i][j][k] * (M_Concentration[3][0] * tNSD3Q7[0] + M_Concentration[3][1] * tNSD3Q7[1] + M_Concentration[3][2] * tNSD3Q7[2] + M_Concentration[3][3] * tNSD3Q7[3] + M_Concentration[3][4] * tNSD3Q7[4] + M_Concentration[3][5] * tNSD3Q7[5] + M_Concentration[3][6] * tNSD3Q7[6]);
						mfEq4 = A_gConcentration[4][4] * ConcentrationFirst[i][j][k] * (M_Concentration[4][0] * tNSD3Q7[0] + M_Concentration[4][1] * tNSD3Q7[1] + M_Concentration[4][2] * tNSD3Q7[2] + M_Concentration[4][3] * tNSD3Q7[3] + M_Concentration[4][4] * tNSD3Q7[4] + M_Concentration[4][5] * tNSD3Q7[5] + M_Concentration[4][6] * tNSD3Q7[6]);
						mfEq5 = A_gConcentration[5][5] * ConcentrationFirst[i][j][k] * (M_Concentration[5][0] * tNSD3Q7[0] + M_Concentration[5][1] * tNSD3Q7[1] + M_Concentration[5][2] * tNSD3Q7[2] + M_Concentration[5][3] * tNSD3Q7[3] + M_Concentration[5][4] * tNSD3Q7[4] + M_Concentration[5][5] * tNSD3Q7[5] + M_Concentration[5][6] * tNSD3Q7[6]);
						mfEq6 = A_gConcentration[6][6] * ConcentrationFirst[i][j][k] * (M_Concentration[6][0] * tNSD3Q7[0] + M_Concentration[6][1] * tNSD3Q7[1] + M_Concentration[6][2] * tNSD3Q7[2] + M_Concentration[6][3] * tNSD3Q7[3] + M_Concentration[6][4] * tNSD3Q7[4] + M_Concentration[6][5] * tNSD3Q7[5] + M_Concentration[6][6] * tNSD3Q7[6]);

						Out0[i][j][k] = cIn0First[i][j][k] - (vfConcentrationIn0 - vfConcentrationEq0);
						Out1[i][j][k] = cIn1First[i][j][k] - (vfConcentrationIn1 - vfConcentrationEq1);
						Out2[i][j][k] = cIn2First[i][j][k] - (vfConcentrationIn2 - vfConcentrationEq2);
						Out3[i][j][k] = cIn3First[i][j][k] - (vfConcentrationIn3 - vfConcentrationEq3);
						Out4[i][j][k] = cIn4First[i][j][k] - (vfConcentrationIn4 - vfConcentrationEq4);
						Out5[i][j][k] = cIn5First[i][j][k] - (vfConcentrationIn5 - vfConcentrationEq5);
						Out6[i][j][k] = cIn6First[i][j][k] - (vfConcentrationIn6 - vfConcentrationEq6);
					}
				}
			}
		}
	}
}

void CollisionSRTConcentrationSecond(int Length,int Width,int Height,double CrelaxationTime,double gasCritical){
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
						Out0[i][j][k]=cIn0Second[i][j][k]-CrelaxationTime*(cIn0Second[i][j][k]-Eq0[i][j][k]*ConcentrationSecond[i][j][k]);
						Out1[i][j][k]=cIn1Second[i][j][k]-CrelaxationTime*(cIn1Second[i][j][k]-Eq1[i][j][k]*ConcentrationSecond[i][j][k]);
						Out2[i][j][k]=cIn2Second[i][j][k]-CrelaxationTime*(cIn2Second[i][j][k]-Eq2[i][j][k]*ConcentrationSecond[i][j][k]);
						Out3[i][j][k]=cIn3Second[i][j][k]-CrelaxationTime*(cIn3Second[i][j][k]-Eq3[i][j][k]*ConcentrationSecond[i][j][k]);
						Out4[i][j][k]=cIn4Second[i][j][k]-CrelaxationTime*(cIn4Second[i][j][k]-Eq4[i][j][k]*ConcentrationSecond[i][j][k]);
						Out5[i][j][k]=cIn5Second[i][j][k]-CrelaxationTime*(cIn5Second[i][j][k]-Eq5[i][j][k]*ConcentrationSecond[i][j][k]);
						Out6[i][j][k]=cIn6Second[i][j][k]-CrelaxationTime*(cIn6Second[i][j][k]-Eq6[i][j][k]*ConcentrationSecond[i][j][k]);
						Out7[i][j][k]=cIn7Second[i][j][k]-CrelaxationTime*(cIn7Second[i][j][k]-Eq7[i][j][k]*ConcentrationSecond[i][j][k]);
						Out8[i][j][k]=cIn8Second[i][j][k]-CrelaxationTime*(cIn8Second[i][j][k]-Eq8[i][j][k]*ConcentrationSecond[i][j][k]);
						Out9[i][j][k]=cIn9Second[i][j][k]-CrelaxationTime*(cIn9Second[i][j][k]-Eq9[i][j][k]*ConcentrationSecond[i][j][k]);
						Out10[i][j][k]=cIn10Second[i][j][k]-CrelaxationTime*(cIn10Second[i][j][k]-Eq10[i][j][k]*ConcentrationSecond[i][j][k]);
						Out11[i][j][k]=cIn11Second[i][j][k]-CrelaxationTime*(cIn11Second[i][j][k]-Eq11[i][j][k]*ConcentrationSecond[i][j][k]);
						Out12[i][j][k]=cIn12Second[i][j][k]-CrelaxationTime*(cIn12Second[i][j][k]-Eq12[i][j][k]*ConcentrationSecond[i][j][k]);
						Out13[i][j][k]=cIn13Second[i][j][k]-CrelaxationTime*(cIn13Second[i][j][k]-Eq13[i][j][k]*ConcentrationSecond[i][j][k]);
						Out14[i][j][k]=cIn14Second[i][j][k]-CrelaxationTime*(cIn14Second[i][j][k]-Eq14[i][j][k]*ConcentrationSecond[i][j][k]);
						Out15[i][j][k]=cIn15Second[i][j][k]-CrelaxationTime*(cIn15Second[i][j][k]-Eq15[i][j][k]*ConcentrationSecond[i][j][k]);
						Out16[i][j][k]=cIn16Second[i][j][k]-CrelaxationTime*(cIn16Second[i][j][k]-Eq16[i][j][k]*ConcentrationSecond[i][j][k]);
						Out17[i][j][k]=cIn17Second[i][j][k]-CrelaxationTime*(cIn17Second[i][j][k]-Eq17[i][j][k]*ConcentrationSecond[i][j][k]);
						Out18[i][j][k]=cIn18Second[i][j][k]-CrelaxationTime*(cIn18Second[i][j][k]-Eq18[i][j][k]*ConcentrationSecond[i][j][k]);
					}
				}
			}
		}
	}
}

void CollisionSRTConcentrationSecondD3Q7(int Length, int Width, int Height, double CrelaxationTime, double gasCritical){
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
						Out0[i][j][k] = cIn0Second[i][j][k] - CrelaxationTime*(cIn0Second[i][j][k] - tNSD3Q7[0] * ConcentrationSecond[i][j][k]);
						Out1[i][j][k] = cIn1Second[i][j][k] - CrelaxationTime*(cIn1Second[i][j][k] - tNSD3Q7[1] * ConcentrationSecond[i][j][k]);
						Out2[i][j][k] = cIn2Second[i][j][k] - CrelaxationTime*(cIn2Second[i][j][k] - tNSD3Q7[2] * ConcentrationSecond[i][j][k]);
						Out3[i][j][k] = cIn3Second[i][j][k] - CrelaxationTime*(cIn3Second[i][j][k] - tNSD3Q7[3] * ConcentrationSecond[i][j][k]);
						Out4[i][j][k] = cIn4Second[i][j][k] - CrelaxationTime*(cIn4Second[i][j][k] - tNSD3Q7[4] * ConcentrationSecond[i][j][k]);
						Out5[i][j][k] = cIn5Second[i][j][k] - CrelaxationTime*(cIn5Second[i][j][k] - tNSD3Q7[5] * ConcentrationSecond[i][j][k]);
						Out6[i][j][k] = cIn6Second[i][j][k] - CrelaxationTime*(cIn6Second[i][j][k] - tNSD3Q7[6] * ConcentrationSecond[i][j][k]);
					}
				}
			}
		}
	}
}

void MemoryFreeConcentration(int m, int n, int q){
	freememory(Concentration,m,n,q);
	freememory(domainConnection,m,n,q);
	freememory(wallNearConcentration,m,n,q);
	freememory(wallNearConcentrationFirst,m,n,q);
	freememory(convectTermConcentration, m, n, q);
	freememory(diffusionTermConcentration, m, n, q);
	freememory(reactantSurface1, m, n, q);
	freememory(reactantSurface2, m, n, q);
	freememory(reactantSurface3, m, n, q);
	freememory(reactantSurface4, m, n, q);
	freememory(reactantSurface5, m, n, q);
	freememory(reactantSurface6, m, n, q);
	freememory(massTransCoeMatrix, m, n, q);
	freememory(electrolyteConcentrationRelaxationTime,m,n,q);
	freememory(cIn0,m,n,q);
	freememory(cIn1,m,n,q);
	freememory(cIn2,m,n,q);
	freememory(cIn3,m,n,q);
	freememory(cIn4,m,n,q);
	freememory(cIn5,m,n,q);
	freememory(cIn6,m,n,q);
	freememory(cIn7,m,n,q);
	freememory(cIn8,m,n,q);
	freememory(cIn9,m,n,q);
	freememory(cIn10,m,n,q);
	freememory(cIn11,m,n,q);
	freememory(cIn12,m,n,q);
	freememory(cIn13,m,n,q);
	freememory(cIn14,m,n,q);
	freememory(cIn15,m,n,q);
	freememory(cIn16,m,n,q);
	freememory(cIn17,m,n,q);
	freememory(cIn18,m,n,q);
}

void MemoryFreeConcentrationFirst(int m, int n, int q){
	freememory(ConcentrationFirst,m,n,q);
	freememory(convectTermConcentrationFirst, m, n, q);
	freememory(diffusionTermConcentrationFirst,m,n,q);
	freememory(cIn0First,m,n,q);
	freememory(cIn1First,m,n,q);
	freememory(cIn2First,m,n,q);
	freememory(cIn3First,m,n,q);
	freememory(cIn4First,m,n,q);
	freememory(cIn5First,m,n,q);
	freememory(cIn6First,m,n,q);
	freememory(cIn7First,m,n,q);
	freememory(cIn8First,m,n,q);
	freememory(cIn9First,m,n,q);
	freememory(cIn10First,m,n,q);
	freememory(cIn11First,m,n,q);
	freememory(cIn12First,m,n,q);
	freememory(cIn13First,m,n,q);
	freememory(cIn14First,m,n,q);
	freememory(cIn15First,m,n,q);
	freememory(cIn16First,m,n,q);
	freememory(cIn17First,m,n,q);
	freememory(cIn18First,m,n,q);
}

void MemoryFreeConcentrationSecond(int m, int n, int q){
	freememory(ConcentrationSecond,m,n,q);
	freememory(convectTermConcentrationSecond, m, n, q);
	freememory(cIn0Second,m,n,q);
	freememory(cIn1Second,m,n,q);
	freememory(cIn2Second,m,n,q);
	freememory(cIn3Second,m,n,q);
	freememory(cIn4Second,m,n,q);
	freememory(cIn5Second,m,n,q);
	freememory(cIn6Second,m,n,q);
	freememory(cIn7Second,m,n,q);
	freememory(cIn8Second,m,n,q);
	freememory(cIn9Second,m,n,q);
	freememory(cIn10Second,m,n,q);
	freememory(cIn11Second,m,n,q);
	freememory(cIn12Second,m,n,q);
	freememory(cIn13Second,m,n,q);
	freememory(cIn14Second,m,n,q);
	freememory(cIn15Second,m,n,q);
	freememory(cIn16Second,m,n,q);
	freememory(cIn17Second,m,n,q);
	freememory(cIn18Second,m,n,q);
}

void StreamConcentration(int m, int n, int q){
	shift0(Out0,cIn0,m,n,q);
	shift1(Out1,cIn1,m,n,q);
	shift2(Out2,cIn2,m,n,q);
	shift3(Out3,cIn3,m,n,q);
	shift4(Out4,cIn4,m,n,q);
	shift5(Out5,cIn5,m,n,q);
	shift6(Out6,cIn6,m,n,q);
	shift7(Out7,cIn7,m,n,q);
	shift8(Out8,cIn8,m,n,q);
	shift9(Out9,cIn9,m,n,q);
	shift10(Out10,cIn10,m,n,q);
	shift11(Out11,cIn11,m,n,q);
	shift12(Out12,cIn12,m,n,q);
	shift13(Out13,cIn13,m,n,q);
	shift14(Out14,cIn14,m,n,q);
	shift15(Out15,cIn15,m,n,q);
	shift16(Out16,cIn16,m,n,q);
	shift17(Out17,cIn17,m,n,q);
	shift18(Out18,cIn18,m,n,q);
}

void StreamConcentrationD3Q7(int m, int n, int q){
	shift0(Out0,cIn0,m,n,q);
	shift1(Out1,cIn1,m,n,q);
	shift2(Out2,cIn2,m,n,q);
	shift3(Out3,cIn3,m,n,q);
	shift4(Out4,cIn4,m,n,q);
	shift5(Out5,cIn5,m,n,q);
	shift6(Out6,cIn6,m,n,q);
}

void StreamConcentrationD3Q15(int m, int n, int q){
	shift0(Out0, cIn0, m, n, q);
	shift1(Out1, cIn1, m, n, q);
	shift4(Out2, cIn2, m, n, q);
	shift2(Out3, cIn3, m, n, q);
	shift5(Out4, cIn4, m, n, q);
	shift3(Out5, cIn5, m, n, q);
	shift6(Out6, cIn6, m, n, q);
	shift19(Out7, cIn7, m, n, q);
	shift20(Out8, cIn8, m, n, q);
	shift21(Out9, cIn9, m, n, q);
	shift22(Out10, cIn10, m, n, q);
	shift23(Out11, cIn11, m, n, q);
	shift24(Out12, cIn12, m, n, q);
	shift25(Out13, cIn13, m, n, q);
	shift26(Out14, cIn14, m, n, q);
}

void StreamConcentrationFirst(int m, int n, int q){
	shift0(Out0,cIn0First,m,n,q);
	shift1(Out1,cIn1First,m,n,q);
	shift2(Out2,cIn2First,m,n,q);
	shift3(Out3,cIn3First,m,n,q);
	shift4(Out4,cIn4First,m,n,q);
	shift5(Out5,cIn5First,m,n,q);
	shift6(Out6,cIn6First,m,n,q);
	shift7(Out7,cIn7First,m,n,q);
	shift8(Out8,cIn8First,m,n,q);
	shift9(Out9,cIn9First,m,n,q);
	shift10(Out10,cIn10First,m,n,q);
	shift11(Out11,cIn11First,m,n,q);
	shift12(Out12,cIn12First,m,n,q);
	shift13(Out13,cIn13First,m,n,q);
	shift14(Out14,cIn14First,m,n,q);
	shift15(Out15,cIn15First,m,n,q);
	shift16(Out16,cIn16First,m,n,q);
	shift17(Out17,cIn17First,m,n,q);
	shift18(Out18,cIn18First,m,n,q);
}

void StreamConcentrationElectrodeD3Q7(int m, int n, int q){
	shift0(Out0,cIn0First,m,n,q);
	shift1(Out1,cIn1First,m,n,q);
	shift2(Out2,cIn2First,m,n,q);
	shift3(Out3,cIn3First,m,n,q);
	shift4(Out4,cIn4First,m,n,q);
	shift5(Out5,cIn5First,m,n,q);
	shift6(Out6,cIn6First,m,n,q);
}

void StreamConcentrationSecond(int m, int n, int q){
	shift0(Out0,cIn0Second,m,n,q);
	shift1(Out1,cIn1Second,m,n,q);
	shift2(Out2,cIn2Second,m,n,q);
	shift3(Out3,cIn3Second,m,n,q);
	shift4(Out4,cIn4Second,m,n,q);
	shift5(Out5,cIn5Second,m,n,q);
	shift6(Out6,cIn6Second,m,n,q);
	shift7(Out7,cIn7Second,m,n,q);
	shift8(Out8,cIn8Second,m,n,q);
	shift9(Out9,cIn9Second,m,n,q);
	shift10(Out10,cIn10Second,m,n,q);
	shift11(Out11,cIn11Second,m,n,q);
	shift12(Out12,cIn12Second,m,n,q);
	shift13(Out13,cIn13Second,m,n,q);
	shift14(Out14,cIn14Second,m,n,q);
	shift15(Out15,cIn15Second,m,n,q);
	shift16(Out16,cIn16Second,m,n,q);
	shift17(Out17,cIn17Second,m,n,q);
	shift18(Out18,cIn18Second,m,n,q);
}

void StreamConcentrationSecondD3Q7(int m, int n, int q){
	shift0(Out0,cIn0Second,m,n,q);
	shift1(Out1,cIn1Second,m,n,q);
	shift2(Out2,cIn2Second,m,n,q);
	shift3(Out3,cIn3Second,m,n,q);
	shift4(Out4,cIn4Second,m,n,q);
	shift5(Out5,cIn5Second,m,n,q);
	shift6(Out6,cIn6Second,m,n,q);
}

void FieldCalculationConcentration(int m, int n, int q, double gasCritical){
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
						Concentration[i][j][k]=cIn0[i][j][k]+cIn1[i][j][k]+cIn2[i][j][k]+cIn3[i][j][k]+cIn4[i][j][k]+cIn5[i][j][k]+cIn6[i][j][k]+cIn7[i][j][k]+cIn8[i][j][k]\
							+cIn9[i][j][k]+cIn10[i][j][k]+cIn11[i][j][k]+cIn12[i][j][k]+cIn13[i][j][k]+cIn14[i][j][k]+cIn15[i][j][k]+cIn16[i][j][k]+cIn17[i][j][k]+cIn18[i][j][k];
					}
					if (domain[i][j][k]!=0||rho[i][j][k]<=gasCritical){
						Concentration[i][j][k]=0.0;
						cIn0[i][j][k]=0.0;
						cIn1[i][j][k]=0.0;
						cIn2[i][j][k]=0.0;
						cIn3[i][j][k]=0.0;
						cIn4[i][j][k]=0.0;
						cIn5[i][j][k]=0.0;
						cIn6[i][j][k]=0.0;
						cIn7[i][j][k]=0.0;
						cIn8[i][j][k]=0.0;
						cIn9[i][j][k]=0.0;
						cIn10[i][j][k]=0.0;
						cIn11[i][j][k]=0.0;
						cIn12[i][j][k]=0.0;
						cIn13[i][j][k]=0.0;
						cIn14[i][j][k]=0.0;
						cIn15[i][j][k]=0.0;
						cIn16[i][j][k]=0.0;
						cIn17[i][j][k]=0.0;
						cIn18[i][j][k]=0.0;
					}
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]!=1){
						cIn0[i][j][k]=Concentration[i][j][k]*tNS[0];
						cIn1[i][j][k]=Concentration[i][j][k]*tNS[1];
						cIn2[i][j][k]=Concentration[i][j][k]*tNS[2];
						cIn3[i][j][k]=Concentration[i][j][k]*tNS[3];
						cIn4[i][j][k]=Concentration[i][j][k]*tNS[4];
						cIn5[i][j][k]=Concentration[i][j][k]*tNS[5];
						cIn6[i][j][k]=Concentration[i][j][k]*tNS[6];
						cIn7[i][j][k]=Concentration[i][j][k]*tNS[7];
						cIn8[i][j][k]=Concentration[i][j][k]*tNS[8];
						cIn9[i][j][k]=Concentration[i][j][k]*tNS[9];
						cIn10[i][j][k]=Concentration[i][j][k]*tNS[10];
						cIn11[i][j][k]=Concentration[i][j][k]*tNS[11];
						cIn12[i][j][k]=Concentration[i][j][k]*tNS[12];
						cIn13[i][j][k]=Concentration[i][j][k]*tNS[13];
						cIn14[i][j][k]=Concentration[i][j][k]*tNS[14];
						cIn15[i][j][k]=Concentration[i][j][k]*tNS[15];
						cIn16[i][j][k]=Concentration[i][j][k]*tNS[16];
						cIn17[i][j][k]=Concentration[i][j][k]*tNS[17];
						cIn18[i][j][k]=Concentration[i][j][k]*tNS[18];
					}					
				}
			}
		}
	}
}

void FieldCalculationConcentrationD3Q7(int m, int n, int q, double gasCritical){
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
						Concentration[i][j][k]=cIn0[i][j][k]+cIn1[i][j][k]+cIn2[i][j][k]+cIn3[i][j][k]+cIn4[i][j][k]+cIn5[i][j][k]+cIn6[i][j][k];
					}
					if (domain[i][j][k]!=0||rho[i][j][k]<=gasCritical||domainConnection[i][j][k]!=1){
						Concentration[i][j][k]=0.0;
						cIn0[i][j][k]=0.0;
						cIn1[i][j][k]=0.0;
						cIn2[i][j][k]=0.0;
						cIn3[i][j][k]=0.0;
						cIn4[i][j][k]=0.0;
						cIn5[i][j][k]=0.0;
						cIn6[i][j][k]=0.0;
					}
/*					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]!=1){
						cIn0[i][j][k]=Concentration[i][j][k]*tNSD3Q7[0];
						cIn1[i][j][k]=Concentration[i][j][k]*tNSD3Q7[1];
						cIn2[i][j][k]=Concentration[i][j][k]*tNSD3Q7[2];
						cIn3[i][j][k]=Concentration[i][j][k]*tNSD3Q7[3];
						cIn4[i][j][k]=Concentration[i][j][k]*tNSD3Q7[4];
						cIn5[i][j][k]=Concentration[i][j][k]*tNSD3Q7[5];
						cIn6[i][j][k]=Concentration[i][j][k]*tNSD3Q7[6];
					}*/					
				}
			}
		}
	}
}

void FieldCalculationConcentrationD3Q15(int m, int n, int q, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<m; i++){
			for (j = 0; j<n; j++){
				for (k = 0; k<q; k++){
					if (domain[i][j][k] == 0 && rho[i][j][k]>gasCritical&&domainConnection[i][j][k] == 1){
						Concentration[i][j][k] = cIn0[i][j][k] + cIn1[i][j][k] + cIn2[i][j][k] + cIn3[i][j][k] + cIn4[i][j][k] + cIn5[i][j][k] + cIn6[i][j][k] + cIn7[i][j][k] + cIn8[i][j][k] + cIn9[i][j][k] + cIn10[i][j][k] + cIn11[i][j][k] + cIn12[i][j][k] + cIn13[i][j][k] + cIn14[i][j][k];
					}
					if (domain[i][j][k] != 0 || rho[i][j][k] <= gasCritical || domainConnection[i][j][k] != 1){
						Concentration[i][j][k] = 0.0;
						cIn0[i][j][k] = 0.0;
						cIn1[i][j][k] = 0.0;
						cIn2[i][j][k] = 0.0;
						cIn3[i][j][k] = 0.0;
						cIn4[i][j][k] = 0.0;
						cIn5[i][j][k] = 0.0;
						cIn6[i][j][k] = 0.0;
						cIn7[i][j][k] = 0.0;
						cIn8[i][j][k] = 0.0;
						cIn9[i][j][k] = 0.0;
						cIn10[i][j][k] = 0.0;
						cIn11[i][j][k] = 0.0;
						cIn12[i][j][k] = 0.0;
						cIn13[i][j][k] = 0.0;
						cIn14[i][j][k] = 0.0;
					}
				}
			}
		}
	}
}

void FieldCalculationConcentrationFirst(int m, int n, int q, double gasCritical){
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
						ConcentrationFirst[i][j][k]=cIn0First[i][j][k]+cIn1First[i][j][k]+cIn2First[i][j][k]+cIn3First[i][j][k]+cIn4First[i][j][k]+cIn5First[i][j][k]+cIn6First[i][j][k]+cIn7First[i][j][k]+cIn8First[i][j][k]\
							+cIn9First[i][j][k]+cIn10First[i][j][k]+cIn11First[i][j][k]+cIn12First[i][j][k]+cIn13First[i][j][k]+cIn14First[i][j][k]+cIn15First[i][j][k]+cIn16First[i][j][k]+cIn17First[i][j][k]+cIn18First[i][j][k];
					}
					if (domain[i][j][k]!=0||rho[i][j][k]<=gasCritical){
						ConcentrationFirst[i][j][k]=0.0;
						cIn0First[i][j][k]=0.0;
						cIn1First[i][j][k]=0.0;
						cIn2First[i][j][k]=0.0;
						cIn3First[i][j][k]=0.0;
						cIn4First[i][j][k]=0.0;
						cIn5First[i][j][k]=0.0;
						cIn6First[i][j][k]=0.0;
						cIn7First[i][j][k]=0.0;
						cIn8First[i][j][k]=0.0;
						cIn9First[i][j][k]=0.0;
						cIn10First[i][j][k]=0.0;
						cIn11First[i][j][k]=0.0;
						cIn12First[i][j][k]=0.0;
						cIn13First[i][j][k]=0.0;
						cIn14First[i][j][k]=0.0;
						cIn15First[i][j][k]=0.0;
						cIn16First[i][j][k]=0.0;
						cIn17First[i][j][k]=0.0;
						cIn18First[i][j][k]=0.0;
					}
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]!=1){
						cIn0First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[0];
						cIn1First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[1];
						cIn2First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[2];
						cIn3First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[3];
						cIn4First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[4];
						cIn5First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[5];
						cIn6First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[6];
						cIn7First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[7];
						cIn8First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[8];
						cIn9First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[9];
						cIn10First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[10];
						cIn11First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[11];
						cIn12First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[12];
						cIn13First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[13];
						cIn14First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[14];
						cIn15First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[15];
						cIn16First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[16];
						cIn17First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[17];
						cIn18First[i][j][k]=ConcentrationFirst[i][j][k]*tNS[18];
					}					
				}
			}
		}
	}
}

void FieldCalculationConcentrationElectrodeD3Q7(int m, int n, int q){
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
						ConcentrationFirst[i][j][k]=cIn0First[i][j][k]+cIn1First[i][j][k]+cIn2First[i][j][k]+cIn3First[i][j][k]+cIn4First[i][j][k]+cIn5First[i][j][k]+cIn6First[i][j][k];
					}
					if (domain[i][j][k]!=1){
						ConcentrationFirst[i][j][k]=0.0;
						cIn0First[i][j][k]=0.0;
						cIn1First[i][j][k]=0.0;
						cIn2First[i][j][k]=0.0;
						cIn3First[i][j][k]=0.0;
						cIn4First[i][j][k]=0.0;
						cIn5First[i][j][k]=0.0;
						cIn6First[i][j][k]=0.0;
					}
				}
			}
		}
	}
}

void FieldCalculationConcentrationSecond(int m, int n, int q, double gasCritical){
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
						ConcentrationSecond[i][j][k]=cIn0Second[i][j][k]+cIn1Second[i][j][k]+cIn2Second[i][j][k]+cIn3Second[i][j][k]+cIn4Second[i][j][k]+cIn5Second[i][j][k]+cIn6Second[i][j][k]+cIn7Second[i][j][k]+cIn8Second[i][j][k]\
							+cIn9Second[i][j][k]+cIn10Second[i][j][k]+cIn11Second[i][j][k]+cIn12Second[i][j][k]+cIn13Second[i][j][k]+cIn14Second[i][j][k]+cIn15Second[i][j][k]+cIn16Second[i][j][k]+cIn17Second[i][j][k]+cIn18Second[i][j][k];
					}
					if (domain[i][j][k]!=0||rho[i][j][k]<=gasCritical){
						ConcentrationSecond[i][j][k]=0.0;
						cIn0Second[i][j][k]=0.0;
						cIn1Second[i][j][k]=0.0;
						cIn2Second[i][j][k]=0.0;
						cIn3Second[i][j][k]=0.0;
						cIn4Second[i][j][k]=0.0;
						cIn5Second[i][j][k]=0.0;
						cIn6Second[i][j][k]=0.0;
						cIn7Second[i][j][k]=0.0;
						cIn8Second[i][j][k]=0.0;
						cIn9Second[i][j][k]=0.0;
						cIn10Second[i][j][k]=0.0;
						cIn11Second[i][j][k]=0.0;
						cIn12Second[i][j][k]=0.0;
						cIn13Second[i][j][k]=0.0;
						cIn14Second[i][j][k]=0.0;
						cIn15Second[i][j][k]=0.0;
						cIn16Second[i][j][k]=0.0;
						cIn17Second[i][j][k]=0.0;
						cIn18Second[i][j][k]=0.0;
					}
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]!=1){
						cIn0Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[0];
						cIn1Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[1];
						cIn2Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[2];
						cIn3Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[3];
						cIn4Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[4];
						cIn5Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[5];
						cIn6Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[6];
						cIn7Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[7];
						cIn8Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[8];
						cIn9Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[9];
						cIn10Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[10];
						cIn11Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[11];
						cIn12Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[12];
						cIn13Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[13];
						cIn14Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[14];
						cIn15Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[15];
						cIn16Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[16];
						cIn17Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[17];
						cIn18Second[i][j][k]=ConcentrationSecond[i][j][k]*tNS[18];
					}					
				}
			}
		}
	}
}

void FieldCalculationConcentrationSecondD3Q7(int m, int n, int q, double gasCritical){
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
						ConcentrationSecond[i][j][k]=cIn0Second[i][j][k]+cIn1Second[i][j][k]+cIn2Second[i][j][k]+cIn3Second[i][j][k]+cIn4Second[i][j][k]+cIn5Second[i][j][k]+cIn6Second[i][j][k];
					}
					if (domain[i][j][k]!=0||rho[i][j][k]<=gasCritical||domainConnection[i][j][k]!=1){
						ConcentrationSecond[i][j][k]=0.0;
						cIn0Second[i][j][k]=0.0;
						cIn1Second[i][j][k]=0.0;
						cIn2Second[i][j][k]=0.0;
						cIn3Second[i][j][k]=0.0;
						cIn4Second[i][j][k]=0.0;
						cIn5Second[i][j][k]=0.0;
						cIn6Second[i][j][k]=0.0;
					}
/*					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]!=1){
						cIn0Second[i][j][k]=ConcentrationSecond[i][j][k]*tNSD3Q7[0];
						cIn1Second[i][j][k]=ConcentrationSecond[i][j][k]*tNSD3Q7[1];
						cIn2Second[i][j][k]=ConcentrationSecond[i][j][k]*tNSD3Q7[2];
						cIn3Second[i][j][k]=ConcentrationSecond[i][j][k]*tNSD3Q7[3];
						cIn4Second[i][j][k]=ConcentrationSecond[i][j][k]*tNSD3Q7[4];
						cIn5Second[i][j][k]=ConcentrationSecond[i][j][k]*tNSD3Q7[5];
						cIn6Second[i][j][k]=ConcentrationSecond[i][j][k]*tNSD3Q7[6];
					}*/					
				}
			}
		}
	}
}

void reactionConcentrationD3Q7(int m, int n, int q, double gasCritical, double latticeF, double transferenceNumber, double totalCurrentResidual){
	int i;
	int j;
	int k;
	double reactant;
	double residual;
#pragma omp parallel private(i,j,k,reactant,residual)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=1;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						residual = 1.0 + totalCurrentResidual;
						if (Domain1[i][j][k] == 1){
							reactant = -neightbour4(i, j, k, m, reaction4) / latticeF;
							reactant = reactant / residual;
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * (reactant*(1 - transferenceNumber));
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * (reactant*(1 - transferenceNumber));
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * (reactant*(1 - transferenceNumber));
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * (reactant*(1 - transferenceNumber));
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * (reactant*(1 - transferenceNumber));
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * (reactant*(1 - transferenceNumber));
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * (reactant*(1 - transferenceNumber));
							}
						if (Domain2[i][j][k] == 1 && j != 0){
							reactant = -neightbour5(i, j, k, n, reaction5) / latticeF;
							reactant = reactant / residual;
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * (reactant*(1 - transferenceNumber));
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * (reactant*(1 - transferenceNumber));
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * (reactant*(1 - transferenceNumber));
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * (reactant*(1 - transferenceNumber));
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * (reactant*(1 - transferenceNumber));
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * (reactant*(1 - transferenceNumber));
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * (reactant*(1 - transferenceNumber));
							}
						if (Domain3[i][j][k] == 1){
							reactant = -neightbour6(i, j, k, q, reaction6) / latticeF;
							reactant = reactant / residual;
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * (reactant*(1 - transferenceNumber));
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * (reactant*(1 - transferenceNumber));
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * (reactant*(1 - transferenceNumber));
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * (reactant*(1 - transferenceNumber));
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * (reactant*(1 - transferenceNumber));
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * (reactant*(1 - transferenceNumber));
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * (reactant*(1 - transferenceNumber));
							}
						if (Domain4[i][j][k] == 1){
							reactant = -neightbour1(i, j, k, m, reaction1) / latticeF;
							reactant = reactant / residual;
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * (reactant*(1 - transferenceNumber));
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * (reactant*(1 - transferenceNumber));
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * (reactant*(1 - transferenceNumber));
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * (reactant*(1 - transferenceNumber));
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * (reactant*(1 - transferenceNumber));
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * (reactant*(1 - transferenceNumber));
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * (reactant*(1 - transferenceNumber));
							}
						if (Domain5[i][j][k] == 1){
							reactant = -neightbour2(i, j, k, n, reaction2) / latticeF;
							reactant = reactant / residual;
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * (reactant*(1 - transferenceNumber));
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * (reactant*(1 - transferenceNumber));
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * (reactant*(1 - transferenceNumber));
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * (reactant*(1 - transferenceNumber));
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * (reactant*(1 - transferenceNumber));
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * (reactant*(1 - transferenceNumber));
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * (reactant*(1 - transferenceNumber));
							}
						if (Domain6[i][j][k] == 1){
							reactant = -neightbour3(i, j, k, q, reaction3) / latticeF;
							reactant = reactant / residual;
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * (reactant*(1 - transferenceNumber));
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * (reactant*(1 - transferenceNumber));
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * (reactant*(1 - transferenceNumber));
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * (reactant*(1 - transferenceNumber));
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * (reactant*(1 - transferenceNumber));
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * (reactant*(1 - transferenceNumber));
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * (reactant*(1 - transferenceNumber));
							}
					}
				}
			}
		}
	}
}

void reactionConcentrationElectrodeD3Q7(int m, int n, int q, double gasCritical, double latticeF, double totalCurrentResidual){
	int i;
	int j;
	int k;
	double residual;
#pragma omp parallel private(i,j,k,residual)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < m; i++){
			for (j = 1; j < n; j++){
				for (k = 0; k<q; k++){
					if (domain[i][j][k] == 1 && electrodeConnection[i][j][k] == 1){
						residual = totalCurrentResidual+1.0;
						if (Domain1[i][j][k] == 0){
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * (reaction1[i][j][k] / latticeF) / residual;
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * (reaction1[i][j][k] / latticeF) / residual;
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * (reaction1[i][j][k] / latticeF) / residual;
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * (reaction1[i][j][k] / latticeF) / residual;
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * (reaction1[i][j][k] / latticeF) / residual;
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * (reaction1[i][j][k] / latticeF) / residual;
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * (reaction1[i][j][k] / latticeF) / residual;
						}
						if (Domain2[i][j][k] == 0){
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * (reaction2[i][j][k] / latticeF) / residual;
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * (reaction2[i][j][k] / latticeF) / residual;
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * (reaction2[i][j][k] / latticeF) / residual;
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * (reaction2[i][j][k] / latticeF) / residual;
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * (reaction2[i][j][k] / latticeF) / residual;
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * (reaction2[i][j][k] / latticeF) / residual;
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * (reaction2[i][j][k] / latticeF) / residual;
						}
						if (Domain3[i][j][k] == 0){
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * (reaction3[i][j][k] / latticeF) / residual;
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * (reaction3[i][j][k] / latticeF) / residual;
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * (reaction3[i][j][k] / latticeF) / residual;
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * (reaction3[i][j][k] / latticeF) / residual;
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * (reaction3[i][j][k] / latticeF) / residual;
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * (reaction3[i][j][k] / latticeF) / residual;
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * (reaction3[i][j][k] / latticeF) / residual;
						}
						if (Domain4[i][j][k] == 0){
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * (reaction4[i][j][k] / latticeF) / residual;
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * (reaction4[i][j][k] / latticeF) / residual;
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * (reaction4[i][j][k] / latticeF) / residual;
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * (reaction4[i][j][k] / latticeF) / residual;
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * (reaction4[i][j][k] / latticeF) / residual;
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * (reaction4[i][j][k] / latticeF) / residual;
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * (reaction4[i][j][k] / latticeF) / residual;
						}
						if (Domain5[i][j][k] == 0 && j != n - 1){
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * (reaction5[i][j][k] / latticeF) / residual;
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * (reaction5[i][j][k] / latticeF) / residual;
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * (reaction5[i][j][k] / latticeF) / residual;
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * (reaction5[i][j][k] / latticeF) / residual;
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * (reaction5[i][j][k] / latticeF) / residual;
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * (reaction5[i][j][k] / latticeF) / residual;
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * (reaction5[i][j][k] / latticeF) / residual;
						}
						if (Domain6[i][j][k] == 0){
							Out0[i][j][k] = Out0[i][j][k] + tNSD3Q7[0] * (reaction6[i][j][k] / latticeF) / residual;
							Out1[i][j][k] = Out1[i][j][k] + tNSD3Q7[1] * (reaction6[i][j][k] / latticeF) / residual;
							Out2[i][j][k] = Out2[i][j][k] + tNSD3Q7[2] * (reaction6[i][j][k] / latticeF) / residual;
							Out3[i][j][k] = Out3[i][j][k] + tNSD3Q7[3] * (reaction6[i][j][k] / latticeF) / residual;
							Out4[i][j][k] = Out4[i][j][k] + tNSD3Q7[4] * (reaction6[i][j][k] / latticeF) / residual;
							Out5[i][j][k] = Out5[i][j][k] + tNSD3Q7[5] * (reaction6[i][j][k] / latticeF) / residual;
							Out6[i][j][k] = Out6[i][j][k] + tNSD3Q7[6] * (reaction6[i][j][k] / latticeF) / residual;
						}
					}
				}
			}
		}
	}
}

void reactionConcentrationSecondD3Q7(int m, int n, int q, double gasCritical,double RateConstant,double equilibriumPotential,double K,double exchange){
	int i;
	int j;
	int k;
	double rateConstant;
	double electrodePotential;
	double overpotential;
	double reactant;
	double electrodePotentialConnection;
#pragma omp parallel private(i,j,k,rateConstant,electrodePotential,overpotential,reactant,electrodePotentialConnection)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					if (domain[i][j][k]==0&&rho[i][j][k]>gasCritical&&domainConnection[i][j][k]==1){
						if (wallNearConcentration[i][j][k] == 0.5){
							wallNearConcentration[i][j][k] = 0.499999;
						}
						if (wallNearConcentrationFirst[i][j][k] == 0.0){
							wallNearConcentrationFirst[i][j][k] = 0.000001;
						}
						if (Domain1[i][j][k]!=0){
							electrodePotentialConnection=neightbour4(i,j,k,m,electrodeConnection);
							if (electrodePotentialConnection==1){
							electrodePotential=neightbour4(i, j, k, m, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							reactant = 2.0*rateConstant*pow(wallNearConcentration[i][j][k], 0.5)*pow(wallNearConcentrationFirst[i][j][k], 0.5);
							Out0[i][j][k]=Out0[i][j][k]+tNSD3Q7[0]*reactant;
							Out1[i][j][k]=Out1[i][j][k]+tNSD3Q7[1]*reactant;
							Out2[i][j][k]=Out2[i][j][k]+tNSD3Q7[2]*reactant;
							Out3[i][j][k]=Out3[i][j][k]+tNSD3Q7[3]*reactant;
							Out4[i][j][k]=Out4[i][j][k]+tNSD3Q7[4]*reactant;
							Out5[i][j][k]=Out5[i][j][k]+tNSD3Q7[5]*reactant;
							Out6[i][j][k]=Out6[i][j][k]+tNSD3Q7[6]*reactant;
							}
						}
						if (Domain2[i][j][k]!=0&&j!=0){
							electrodePotentialConnection=neightbour5(i,j,k,n,electrodeConnection);
							if (electrodePotentialConnection==1){
							electrodePotential=neightbour5(i, j, k, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							reactant = 2.0*rateConstant*pow(wallNearConcentration[i][j][k], 0.5)*pow(wallNearConcentrationFirst[i][j][k], 0.5);
							Out0[i][j][k]=Out0[i][j][k]+tNSD3Q7[0]*reactant;
							Out1[i][j][k]=Out1[i][j][k]+tNSD3Q7[1]*reactant;
							Out2[i][j][k]=Out2[i][j][k]+tNSD3Q7[2]*reactant;
							Out3[i][j][k]=Out3[i][j][k]+tNSD3Q7[3]*reactant;
							Out4[i][j][k]=Out4[i][j][k]+tNSD3Q7[4]*reactant;
							Out5[i][j][k]=Out5[i][j][k]+tNSD3Q7[5]*reactant;
							Out6[i][j][k]=Out6[i][j][k]+tNSD3Q7[6]*reactant;
							}
						}
						if (Domain3[i][j][k]!=0){
							electrodePotentialConnection=neightbour6(i,j,k,q,electrodeConnection);
							if (electrodePotentialConnection==1){
							electrodePotential=neightbour6(i, j, k, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							reactant = 2.0*rateConstant*pow(wallNearConcentration[i][j][k], 0.5)*pow(wallNearConcentrationFirst[i][j][k], 0.5);
							Out0[i][j][k]=Out0[i][j][k]+tNSD3Q7[0]*reactant;
							Out1[i][j][k]=Out1[i][j][k]+tNSD3Q7[1]*reactant;
							Out2[i][j][k]=Out2[i][j][k]+tNSD3Q7[2]*reactant;
							Out3[i][j][k]=Out3[i][j][k]+tNSD3Q7[3]*reactant;
							Out4[i][j][k]=Out4[i][j][k]+tNSD3Q7[4]*reactant;
							Out5[i][j][k]=Out5[i][j][k]+tNSD3Q7[5]*reactant;
							Out6[i][j][k]=Out6[i][j][k]+tNSD3Q7[6]*reactant;
							}
						}
						if (Domain4[i][j][k]!=0){
							electrodePotentialConnection=neightbour1(i,j,k,m,electrodeConnection);
							if (electrodePotentialConnection==1){
							electrodePotential=neightbour1(i, j, k, m, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							reactant = 2.0*rateConstant*pow(wallNearConcentration[i][j][k], 0.5)*pow(wallNearConcentrationFirst[i][j][k], 0.5);
							Out0[i][j][k]=Out0[i][j][k]+tNSD3Q7[0]*reactant;
							Out1[i][j][k]=Out1[i][j][k]+tNSD3Q7[1]*reactant;
							Out2[i][j][k]=Out2[i][j][k]+tNSD3Q7[2]*reactant;
							Out3[i][j][k]=Out3[i][j][k]+tNSD3Q7[3]*reactant;
							Out4[i][j][k]=Out4[i][j][k]+tNSD3Q7[4]*reactant;
							Out5[i][j][k]=Out5[i][j][k]+tNSD3Q7[5]*reactant;
							Out6[i][j][k]=Out6[i][j][k]+tNSD3Q7[6]*reactant;
							}
						}
						if (Domain5[i][j][k]!=0){
							electrodePotentialConnection=neightbour2(i,j,k,n,electrodeConnection);
							if (electrodePotentialConnection==1){
							electrodePotential=neightbour2(i, j, k, n, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							reactant = 2.0*rateConstant*pow(wallNearConcentration[i][j][k], 0.5)*pow(wallNearConcentrationFirst[i][j][k], 0.5);
							Out0[i][j][k]=Out0[i][j][k]+tNSD3Q7[0]*reactant;
							Out1[i][j][k]=Out1[i][j][k]+tNSD3Q7[1]*reactant;
							Out2[i][j][k]=Out2[i][j][k]+tNSD3Q7[2]*reactant;
							Out3[i][j][k]=Out3[i][j][k]+tNSD3Q7[3]*reactant;
							Out4[i][j][k]=Out4[i][j][k]+tNSD3Q7[4]*reactant;
							Out5[i][j][k]=Out5[i][j][k]+tNSD3Q7[5]*reactant;
							Out6[i][j][k]=Out6[i][j][k]+tNSD3Q7[6]*reactant;
							}
						}
						if (Domain6[i][j][k]!=0){
							electrodePotentialConnection=neightbour3(i,j,k,q,electrodeConnection);
							if (electrodePotentialConnection==1){
							electrodePotential=neightbour3(i, j, k, q, lastWallElectrode);
							overpotential=electrodePotential/exchange-electrolytePotential[i][j][k]/exchange-(equilibriumPotential+8.314*K/96485.0*log(wallNearConcentration[i][j][k]/wallNearConcentrationFirst[i][j][k]));
							rateConstant=RateConstant*(exp(0.5*96485.0*overpotential/8.314/K)-exp(-0.5*96485.0*overpotential/8.314/K));
							reactant = 2.0*rateConstant*pow(wallNearConcentration[i][j][k], 0.5)*pow(wallNearConcentrationFirst[i][j][k], 0.5);
							Out0[i][j][k]=Out0[i][j][k]+tNSD3Q7[0]*reactant;
							Out1[i][j][k]=Out1[i][j][k]+tNSD3Q7[1]*reactant;
							Out2[i][j][k]=Out2[i][j][k]+tNSD3Q7[2]*reactant;
							Out3[i][j][k]=Out3[i][j][k]+tNSD3Q7[3]*reactant;
							Out4[i][j][k]=Out4[i][j][k]+tNSD3Q7[4]*reactant;
							Out5[i][j][k]=Out5[i][j][k]+tNSD3Q7[5]*reactant;
							Out6[i][j][k]=Out6[i][j][k]+tNSD3Q7[6]*reactant;
							}
						}
					}
				}
			}
		}
	}
}

void refillConcentration(int m, int n, int q, double gasCritical,double ***rho1, double ***rho2, double ***ux, double ***uy, double ***uz){
	int i;
	int j;
	int k;
	int number;
	double valueRho;
	double valueRho1;
#pragma omp parallel private(i,j,k,number)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<q;k++){
					number=0;
					if (domain[i][j][k]==0&&rho1[i][j][k]>=gasCritical&&rho2[i][j][k]<gasCritical){		
						cIn0[i][j][k]=0.0;
						cIn1[i][j][k]=0.0;
						cIn2[i][j][k]=0.0;
						cIn3[i][j][k]=0.0;
						cIn4[i][j][k]=0.0;
						cIn5[i][j][k]=0.0;
						cIn6[i][j][k]=0.0;
						cIn7[i][j][k]=0.0;
						cIn8[i][j][k]=0.0;
						cIn9[i][j][k]=0.0;
						cIn10[i][j][k]=0.0;
						cIn11[i][j][k]=0.0;
						cIn12[i][j][k]=0.0;
						cIn13[i][j][k]=0.0;
						cIn14[i][j][k]=0.0;
						cIn15[i][j][k]=0.0;
						cIn16[i][j][k]=0.0;
						cIn17[i][j][k]=0.0;
						cIn18[i][j][k]=0.0;
						Concentration[i][j][k]=0.0;
						valueRho=neightbour4(i, j, k, m, rho1);
						valueRho1=neightbour4(i, j, k, m, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							cIn0[i][j][k]=cIn0[i][j][k]+neightbour4(i, j, k, m, cIn0);
							cIn1[i][j][k]=cIn1[i][j][k]+neightbour4(i, j, k, m, cIn1);
							cIn2[i][j][k]=cIn2[i][j][k]+neightbour4(i, j, k, m, cIn2);
							cIn3[i][j][k]=cIn3[i][j][k]+neightbour4(i, j, k, m, cIn3);
							cIn4[i][j][k]=cIn4[i][j][k]+neightbour4(i, j, k, m, cIn4);
							cIn5[i][j][k]=cIn5[i][j][k]+neightbour4(i, j, k, m, cIn5);
							cIn6[i][j][k]=cIn6[i][j][k]+neightbour4(i, j, k, m, cIn6);
							cIn7[i][j][k]=cIn7[i][j][k]+neightbour4(i, j, k, m, cIn7);
							cIn8[i][j][k]=cIn8[i][j][k]+neightbour4(i, j, k, m, cIn8);
							cIn9[i][j][k]=cIn9[i][j][k]+neightbour4(i, j, k, m, cIn9);
							cIn10[i][j][k]=cIn10[i][j][k]+neightbour4(i, j, k, m, cIn10);
							cIn11[i][j][k]=cIn11[i][j][k]+neightbour4(i, j, k, m, cIn11);
							cIn12[i][j][k]=cIn12[i][j][k]+neightbour4(i, j, k, m, cIn12);
							cIn13[i][j][k]=cIn13[i][j][k]+neightbour4(i, j, k, m, cIn13);
							cIn14[i][j][k]=cIn14[i][j][k]+neightbour4(i, j, k, m, cIn14);
							cIn15[i][j][k]=cIn15[i][j][k]+neightbour4(i, j, k, m, cIn15);
							cIn16[i][j][k]=cIn16[i][j][k]+neightbour4(i, j, k, m, cIn16);
							cIn17[i][j][k]=cIn17[i][j][k]+neightbour4(i, j, k, m, cIn17);
							cIn18[i][j][k]=cIn18[i][j][k]+neightbour4(i, j, k, m, cIn18);
							Concentration[i][j][k]=Concentration[i][j][k]+neightbour4(i, j, k, m, Concentration);
						}
						valueRho=neightbour5(i, j, k, n, rho1);
						valueRho1=neightbour5(i, j, k, n, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							cIn0[i][j][k]=cIn0[i][j][k]+neightbour5(i, j, k, n, cIn0);
							cIn1[i][j][k]=cIn1[i][j][k]+neightbour5(i, j, k, n, cIn1);
							cIn2[i][j][k]=cIn2[i][j][k]+neightbour5(i, j, k, n, cIn2);
							cIn3[i][j][k]=cIn3[i][j][k]+neightbour5(i, j, k, n, cIn3);
							cIn4[i][j][k]=cIn4[i][j][k]+neightbour5(i, j, k, n, cIn4);
							cIn5[i][j][k]=cIn5[i][j][k]+neightbour5(i, j, k, n, cIn5);
							cIn6[i][j][k]=cIn6[i][j][k]+neightbour5(i, j, k, n, cIn6);
							cIn7[i][j][k]=cIn7[i][j][k]+neightbour5(i, j, k, n, cIn7);
							cIn8[i][j][k]=cIn8[i][j][k]+neightbour5(i, j, k, n, cIn8);
							cIn9[i][j][k]=cIn9[i][j][k]+neightbour5(i, j, k, n, cIn9);
							cIn10[i][j][k]=cIn10[i][j][k]+neightbour5(i, j, k, n, cIn10);
							cIn11[i][j][k]=cIn11[i][j][k]+neightbour5(i, j, k, n, cIn11);
							cIn12[i][j][k]=cIn12[i][j][k]+neightbour5(i, j, k, n, cIn12);
							cIn13[i][j][k]=cIn13[i][j][k]+neightbour5(i, j, k, n, cIn13);
							cIn14[i][j][k]=cIn14[i][j][k]+neightbour5(i, j, k, n, cIn14);
							cIn15[i][j][k]=cIn15[i][j][k]+neightbour5(i, j, k, n, cIn15);
							cIn16[i][j][k]=cIn16[i][j][k]+neightbour5(i, j, k, n, cIn16);
							cIn17[i][j][k]=cIn17[i][j][k]+neightbour5(i, j, k, n, cIn17);
							cIn18[i][j][k]=cIn18[i][j][k]+neightbour5(i, j, k, n, cIn18);
							Concentration[i][j][k]=Concentration[i][j][k]+neightbour5(i, j, k, n, Concentration);
						}
						valueRho=neightbour6(i, j, k, q, rho1);
						valueRho1=neightbour6(i, j, k, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							cIn0[i][j][k]=cIn0[i][j][k]+neightbour6(i, j, k, q, cIn0);
							cIn1[i][j][k]=cIn1[i][j][k]+neightbour6(i, j, k, q, cIn1);
							cIn2[i][j][k]=cIn2[i][j][k]+neightbour6(i, j, k, q, cIn2);
							cIn3[i][j][k]=cIn3[i][j][k]+neightbour6(i, j, k, q, cIn3);
							cIn4[i][j][k]=cIn4[i][j][k]+neightbour6(i, j, k, q, cIn4);
							cIn5[i][j][k]=cIn5[i][j][k]+neightbour6(i, j, k, q, cIn5);
							cIn6[i][j][k]=cIn6[i][j][k]+neightbour6(i, j, k, q, cIn6);
							cIn7[i][j][k]=cIn7[i][j][k]+neightbour6(i, j, k, q, cIn7);
							cIn8[i][j][k]=cIn8[i][j][k]+neightbour6(i, j, k, q, cIn8);
							cIn9[i][j][k]=cIn9[i][j][k]+neightbour6(i, j, k, q, cIn9);
							cIn10[i][j][k]=cIn10[i][j][k]+neightbour6(i, j, k, q, cIn10);
							cIn11[i][j][k]=cIn11[i][j][k]+neightbour6(i, j, k, q, cIn11);
							cIn12[i][j][k]=cIn12[i][j][k]+neightbour6(i, j, k, q, cIn12);
							cIn13[i][j][k]=cIn13[i][j][k]+neightbour6(i, j, k, q, cIn13);
							cIn14[i][j][k]=cIn14[i][j][k]+neightbour6(i, j, k, q, cIn14);
							cIn15[i][j][k]=cIn15[i][j][k]+neightbour6(i, j, k, q, cIn15);
							cIn16[i][j][k]=cIn16[i][j][k]+neightbour6(i, j, k, q, cIn16);
							cIn17[i][j][k]=cIn17[i][j][k]+neightbour6(i, j, k, q, cIn17);
							cIn18[i][j][k]=cIn18[i][j][k]+neightbour6(i, j, k, q, cIn18);
							Concentration[i][j][k]=Concentration[i][j][k]+neightbour6(i, j, k, q, Concentration);
						}							
						valueRho=neightbour1(i, j, k, m, rho1);
						valueRho1=neightbour1(i, j, k, m, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							cIn0[i][j][k]=cIn0[i][j][k]+neightbour1(i, j, k, m, cIn0);
							cIn1[i][j][k]=cIn1[i][j][k]+neightbour1(i, j, k, m, cIn1);
							cIn2[i][j][k]=cIn2[i][j][k]+neightbour1(i, j, k, m, cIn2);
							cIn3[i][j][k]=cIn3[i][j][k]+neightbour1(i, j, k, m, cIn3);
							cIn4[i][j][k]=cIn4[i][j][k]+neightbour1(i, j, k, m, cIn4);
							cIn5[i][j][k]=cIn5[i][j][k]+neightbour1(i, j, k, m, cIn5);
							cIn6[i][j][k]=cIn6[i][j][k]+neightbour1(i, j, k, m, cIn6);
							cIn7[i][j][k]=cIn7[i][j][k]+neightbour1(i, j, k, m, cIn7);
							cIn8[i][j][k]=cIn8[i][j][k]+neightbour1(i, j, k, m, cIn8);
							cIn9[i][j][k]=cIn9[i][j][k]+neightbour1(i, j, k, m, cIn9);
							cIn10[i][j][k]=cIn10[i][j][k]+neightbour1(i, j, k, m, cIn10);
							cIn11[i][j][k]=cIn11[i][j][k]+neightbour1(i, j, k, m, cIn11);
							cIn12[i][j][k]=cIn12[i][j][k]+neightbour1(i, j, k, m, cIn12);
							cIn13[i][j][k]=cIn13[i][j][k]+neightbour1(i, j, k, m, cIn13);
							cIn14[i][j][k]=cIn14[i][j][k]+neightbour1(i, j, k, m, cIn14);
							cIn15[i][j][k]=cIn15[i][j][k]+neightbour1(i, j, k, m, cIn15);
							cIn16[i][j][k]=cIn16[i][j][k]+neightbour1(i, j, k, m, cIn16);
							cIn17[i][j][k]=cIn17[i][j][k]+neightbour1(i, j, k, m, cIn17);
							cIn18[i][j][k]=cIn18[i][j][k]+neightbour1(i, j, k, m, cIn18);
							Concentration[i][j][k]=Concentration[i][j][k]+neightbour1(i, j, k, m, Concentration);
						}
						valueRho=neightbour2(i, j, k, n, rho1);
						valueRho1=neightbour2(i, j, k, n, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							cIn0[i][j][k]=cIn0[i][j][k]+neightbour2(i, j, k, n, cIn0);
							cIn1[i][j][k]=cIn1[i][j][k]+neightbour2(i, j, k, n, cIn1);
							cIn2[i][j][k]=cIn2[i][j][k]+neightbour2(i, j, k, n, cIn2);
							cIn3[i][j][k]=cIn3[i][j][k]+neightbour2(i, j, k, n, cIn3);
							cIn4[i][j][k]=cIn4[i][j][k]+neightbour2(i, j, k, n, cIn4);
							cIn5[i][j][k]=cIn5[i][j][k]+neightbour2(i, j, k, n, cIn5);
							cIn6[i][j][k]=cIn6[i][j][k]+neightbour2(i, j, k, n, cIn6);
							cIn7[i][j][k]=cIn7[i][j][k]+neightbour2(i, j, k, n, cIn7);
							cIn8[i][j][k]=cIn8[i][j][k]+neightbour2(i, j, k, n, cIn8);
							cIn9[i][j][k]=cIn9[i][j][k]+neightbour2(i, j, k, n, cIn9);
							cIn10[i][j][k]=cIn10[i][j][k]+neightbour2(i, j, k, n, cIn10);
							cIn11[i][j][k]=cIn11[i][j][k]+neightbour2(i, j, k, n, cIn11);
							cIn12[i][j][k]=cIn12[i][j][k]+neightbour2(i, j, k, n, cIn12);
							cIn13[i][j][k]=cIn13[i][j][k]+neightbour2(i, j, k, n, cIn13);
							cIn14[i][j][k]=cIn14[i][j][k]+neightbour2(i, j, k, n, cIn14);
							cIn15[i][j][k]=cIn15[i][j][k]+neightbour2(i, j, k, n, cIn15);
							cIn16[i][j][k]=cIn16[i][j][k]+neightbour2(i, j, k, n, cIn16);
							cIn17[i][j][k]=cIn17[i][j][k]+neightbour2(i, j, k, n, cIn17);
							cIn18[i][j][k]=cIn18[i][j][k]+neightbour2(i, j, k, n, cIn18);
							Concentration[i][j][k]=Concentration[i][j][k]+neightbour2(i, j, k, n, Concentration);
						}
						valueRho=neightbour3(i, j, k, q, rho1);
						valueRho1=neightbour3(i, j, k, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							cIn0[i][j][k]=cIn0[i][j][k]+neightbour3(i, j, k, q, cIn0);
							cIn1[i][j][k]=cIn1[i][j][k]+neightbour3(i, j, k, q, cIn1);
							cIn2[i][j][k]=cIn2[i][j][k]+neightbour3(i, j, k, q, cIn2);
							cIn3[i][j][k]=cIn3[i][j][k]+neightbour3(i, j, k, q, cIn3);
							cIn4[i][j][k]=cIn4[i][j][k]+neightbour3(i, j, k, q, cIn4);
							cIn5[i][j][k]=cIn5[i][j][k]+neightbour3(i, j, k, q, cIn5);
							cIn6[i][j][k]=cIn6[i][j][k]+neightbour3(i, j, k, q, cIn6);
							cIn7[i][j][k]=cIn7[i][j][k]+neightbour3(i, j, k, q, cIn7);
							cIn8[i][j][k]=cIn8[i][j][k]+neightbour3(i, j, k, q, cIn8);
							cIn9[i][j][k]=cIn9[i][j][k]+neightbour3(i, j, k, q, cIn9);
							cIn10[i][j][k]=cIn10[i][j][k]+neightbour3(i, j, k, q, cIn10);
							cIn11[i][j][k]=cIn11[i][j][k]+neightbour3(i, j, k, q, cIn11);
							cIn12[i][j][k]=cIn12[i][j][k]+neightbour3(i, j, k, q, cIn12);
							cIn13[i][j][k]=cIn13[i][j][k]+neightbour3(i, j, k, q, cIn13);
							cIn14[i][j][k]=cIn14[i][j][k]+neightbour3(i, j, k, q, cIn14);
							cIn15[i][j][k]=cIn15[i][j][k]+neightbour3(i, j, k, q, cIn15);
							cIn16[i][j][k]=cIn16[i][j][k]+neightbour3(i, j, k, q, cIn16);
							cIn17[i][j][k]=cIn17[i][j][k]+neightbour3(i, j, k, q, cIn17);
							cIn18[i][j][k]=cIn18[i][j][k]+neightbour3(i, j, k, q, cIn18);
							Concentration[i][j][k]=Concentration[i][j][k]+neightbour3(i, j, k, q, Concentration);
						}
						valueRho=neightbour10(i, j, k, m, n, rho1);
						valueRho1=neightbour10(i, j, k, m, n, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							cIn0[i][j][k]=cIn0[i][j][k]+neightbour10(i, j, k, m, n, cIn0);
							cIn1[i][j][k]=cIn1[i][j][k]+neightbour10(i, j, k, m, n, cIn1);
							cIn2[i][j][k]=cIn2[i][j][k]+neightbour10(i, j, k, m, n, cIn2);
							cIn3[i][j][k]=cIn3[i][j][k]+neightbour10(i, j, k, m, n, cIn3);
							cIn4[i][j][k]=cIn4[i][j][k]+neightbour10(i, j, k, m, n, cIn4);
							cIn5[i][j][k]=cIn5[i][j][k]+neightbour10(i, j, k, m, n, cIn5);
							cIn6[i][j][k]=cIn6[i][j][k]+neightbour10(i, j, k, m, n, cIn6);
							cIn7[i][j][k]=cIn7[i][j][k]+neightbour10(i, j, k, m, n, cIn7);
							cIn8[i][j][k]=cIn8[i][j][k]+neightbour10(i, j, k, m, n, cIn8);
							cIn9[i][j][k]=cIn9[i][j][k]+neightbour10(i, j, k, m, n, cIn9);
							cIn10[i][j][k]=cIn10[i][j][k]+neightbour10(i, j, k, m, n, cIn10);
							cIn11[i][j][k]=cIn11[i][j][k]+neightbour10(i, j, k, m, n, cIn11);
							cIn12[i][j][k]=cIn12[i][j][k]+neightbour10(i, j, k, m, n, cIn12);
							cIn13[i][j][k]=cIn13[i][j][k]+neightbour10(i, j, k, m, n, cIn13);
							cIn14[i][j][k]=cIn14[i][j][k]+neightbour10(i, j, k, m, n, cIn14);
							cIn15[i][j][k]=cIn15[i][j][k]+neightbour10(i, j, k, m, n, cIn15);
							cIn16[i][j][k]=cIn16[i][j][k]+neightbour10(i, j, k, m, n, cIn16);
							cIn17[i][j][k]=cIn17[i][j][k]+neightbour10(i, j, k, m, n, cIn17);
							cIn18[i][j][k]=cIn18[i][j][k]+neightbour10(i, j, k, m, n, cIn18);
							Concentration[i][j][k]=Concentration[i][j][k]+neightbour10(i, j, k, m, n, Concentration);
						}
						valueRho=neightbour11(i, j, k, m, q, rho1);
						valueRho1=neightbour11(i, j, k, m, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							cIn0[i][j][k]=cIn0[i][j][k]+neightbour11(i, j, k, m, q, cIn0);
							cIn1[i][j][k]=cIn1[i][j][k]+neightbour11(i, j, k, m, q, cIn1);
							cIn2[i][j][k]=cIn2[i][j][k]+neightbour11(i, j, k, m, q, cIn2);
							cIn3[i][j][k]=cIn3[i][j][k]+neightbour11(i, j, k, m, q, cIn3);
							cIn4[i][j][k]=cIn4[i][j][k]+neightbour11(i, j, k, m, q, cIn4);
							cIn5[i][j][k]=cIn5[i][j][k]+neightbour11(i, j, k, m, q, cIn5);
							cIn6[i][j][k]=cIn6[i][j][k]+neightbour11(i, j, k, m, q, cIn6);
							cIn7[i][j][k]=cIn7[i][j][k]+neightbour11(i, j, k, m, q, cIn7);
							cIn8[i][j][k]=cIn8[i][j][k]+neightbour11(i, j, k, m, q, cIn8);
							cIn9[i][j][k]=cIn9[i][j][k]+neightbour11(i, j, k, m, q, cIn9);
							cIn10[i][j][k]=cIn10[i][j][k]+neightbour11(i, j, k, m, q, cIn10);
							cIn11[i][j][k]=cIn11[i][j][k]+neightbour11(i, j, k, m, q, cIn11);
							cIn12[i][j][k]=cIn12[i][j][k]+neightbour11(i, j, k, m, q, cIn12);
							cIn13[i][j][k]=cIn13[i][j][k]+neightbour11(i, j, k, m, q, cIn13);
							cIn14[i][j][k]=cIn14[i][j][k]+neightbour11(i, j, k, m, q, cIn14);
							cIn15[i][j][k]=cIn15[i][j][k]+neightbour11(i, j, k, m, q, cIn15);
							cIn16[i][j][k]=cIn16[i][j][k]+neightbour11(i, j, k, m, q, cIn16);
							cIn17[i][j][k]=cIn17[i][j][k]+neightbour11(i, j, k, m, q, cIn17);
							cIn18[i][j][k]=cIn18[i][j][k]+neightbour11(i, j, k, m, q, cIn18);
							Concentration[i][j][k]=Concentration[i][j][k]+neightbour11(i, j, k, m, q, Concentration);
						}
						valueRho=neightbour12(i, j, k, n, q, rho1);
						valueRho1=neightbour12(i, j, k, n, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							cIn0[i][j][k]=cIn0[i][j][k]+neightbour12(i, j, k, n, q, cIn0);
							cIn1[i][j][k]=cIn1[i][j][k]+neightbour12(i, j, k, n, q, cIn1);
							cIn2[i][j][k]=cIn2[i][j][k]+neightbour12(i, j, k, n, q, cIn2);
							cIn3[i][j][k]=cIn3[i][j][k]+neightbour12(i, j, k, n, q, cIn3);
							cIn4[i][j][k]=cIn4[i][j][k]+neightbour12(i, j, k, n, q, cIn4);
							cIn5[i][j][k]=cIn5[i][j][k]+neightbour12(i, j, k, n, q, cIn5);
							cIn6[i][j][k]=cIn6[i][j][k]+neightbour12(i, j, k, n, q, cIn6);
							cIn7[i][j][k]=cIn7[i][j][k]+neightbour12(i, j, k, n, q, cIn7);
							cIn8[i][j][k]=cIn8[i][j][k]+neightbour12(i, j, k, n, q, cIn8);
							cIn9[i][j][k]=cIn9[i][j][k]+neightbour12(i, j, k, n, q, cIn9);
							cIn10[i][j][k]=cIn10[i][j][k]+neightbour12(i, j, k, n, q, cIn10);
							cIn11[i][j][k]=cIn11[i][j][k]+neightbour12(i, j, k, n, q, cIn11);
							cIn12[i][j][k]=cIn12[i][j][k]+neightbour12(i, j, k, n, q, cIn12);
							cIn13[i][j][k]=cIn13[i][j][k]+neightbour12(i, j, k, n, q, cIn13);
							cIn14[i][j][k]=cIn14[i][j][k]+neightbour12(i, j, k, n, q, cIn14);
							cIn15[i][j][k]=cIn15[i][j][k]+neightbour12(i, j, k, n, q, cIn15);
							cIn16[i][j][k]=cIn16[i][j][k]+neightbour12(i, j, k, n, q, cIn16);
							cIn17[i][j][k]=cIn17[i][j][k]+neightbour12(i, j, k, n, q, cIn17);
							cIn18[i][j][k]=cIn18[i][j][k]+neightbour12(i, j, k, n, q, cIn18);
							Concentration[i][j][k]=Concentration[i][j][k]+neightbour12(i, j, k, n, q, Concentration);
						}
						valueRho=neightbour7(i, j, k, m, n, rho1);
						valueRho1=neightbour7(i, j, k, m, n, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							cIn0[i][j][k]=cIn0[i][j][k]+neightbour7(i, j, k, m, n, cIn0);
							cIn1[i][j][k]=cIn1[i][j][k]+neightbour7(i, j, k, m, n, cIn1);
							cIn2[i][j][k]=cIn2[i][j][k]+neightbour7(i, j, k, m, n, cIn2);
							cIn3[i][j][k]=cIn3[i][j][k]+neightbour7(i, j, k, m, n, cIn3);
							cIn4[i][j][k]=cIn4[i][j][k]+neightbour7(i, j, k, m, n, cIn4);
							cIn5[i][j][k]=cIn5[i][j][k]+neightbour7(i, j, k, m, n, cIn5);
							cIn6[i][j][k]=cIn6[i][j][k]+neightbour7(i, j, k, m, n, cIn6);
							cIn7[i][j][k]=cIn7[i][j][k]+neightbour7(i, j, k, m, n, cIn7);
							cIn8[i][j][k]=cIn8[i][j][k]+neightbour7(i, j, k, m, n, cIn8);
							cIn9[i][j][k]=cIn9[i][j][k]+neightbour7(i, j, k, m, n, cIn9);
							cIn10[i][j][k]=cIn10[i][j][k]+neightbour7(i, j, k, m, n, cIn10);
							cIn11[i][j][k]=cIn11[i][j][k]+neightbour7(i, j, k, m, n, cIn11);
							cIn12[i][j][k]=cIn12[i][j][k]+neightbour7(i, j, k, m, n, cIn12);
							cIn13[i][j][k]=cIn13[i][j][k]+neightbour7(i, j, k, m, n, cIn13);
							cIn14[i][j][k]=cIn14[i][j][k]+neightbour7(i, j, k, m, n, cIn14);
							cIn15[i][j][k]=cIn15[i][j][k]+neightbour7(i, j, k, m, n, cIn15);
							cIn16[i][j][k]=cIn16[i][j][k]+neightbour7(i, j, k, m, n, cIn16);
							cIn17[i][j][k]=cIn17[i][j][k]+neightbour7(i, j, k, m, n, cIn17);
							cIn18[i][j][k]=cIn18[i][j][k]+neightbour7(i, j, k, m, n, cIn18);
							Concentration[i][j][k]=Concentration[i][j][k]+neightbour7(i, j, k, m, n, Concentration);
						}
						valueRho=neightbour8(i, j, k, m, q, rho1);
						valueRho1=neightbour8(i, j, k, m, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							cIn0[i][j][k]=cIn0[i][j][k]+neightbour8(i, j, k, m, q, cIn0);
							cIn1[i][j][k]=cIn1[i][j][k]+neightbour8(i, j, k, m, q, cIn1);
							cIn2[i][j][k]=cIn2[i][j][k]+neightbour8(i, j, k, m, q, cIn2);
							cIn3[i][j][k]=cIn3[i][j][k]+neightbour8(i, j, k, m, q, cIn3);
							cIn4[i][j][k]=cIn4[i][j][k]+neightbour8(i, j, k, m, q, cIn4);
							cIn5[i][j][k]=cIn5[i][j][k]+neightbour8(i, j, k, m, q, cIn5);
							cIn6[i][j][k]=cIn6[i][j][k]+neightbour8(i, j, k, m, q, cIn6);
							cIn7[i][j][k]=cIn7[i][j][k]+neightbour8(i, j, k, m, q, cIn7);
							cIn8[i][j][k]=cIn8[i][j][k]+neightbour8(i, j, k, m, q, cIn8);
							cIn9[i][j][k]=cIn9[i][j][k]+neightbour8(i, j, k, m, q, cIn9);
							cIn10[i][j][k]=cIn10[i][j][k]+neightbour8(i, j, k, m, q, cIn10);
							cIn11[i][j][k]=cIn11[i][j][k]+neightbour8(i, j, k, m, q, cIn11);
							cIn12[i][j][k]=cIn12[i][j][k]+neightbour8(i, j, k, m, q, cIn12);
							cIn13[i][j][k]=cIn13[i][j][k]+neightbour8(i, j, k, m, q, cIn13);
							cIn14[i][j][k]=cIn14[i][j][k]+neightbour8(i, j, k, m, q, cIn14);
							cIn15[i][j][k]=cIn15[i][j][k]+neightbour8(i, j, k, m, q, cIn15);
							cIn16[i][j][k]=cIn16[i][j][k]+neightbour8(i, j, k, m, q, cIn16);
							cIn17[i][j][k]=cIn17[i][j][k]+neightbour8(i, j, k, m, q, cIn17);
							cIn18[i][j][k]=cIn18[i][j][k]+neightbour8(i, j, k, m, q, cIn18);
							Concentration[i][j][k]=Concentration[i][j][k]+neightbour8(i, j, k, m, q, Concentration);
						}
						valueRho=neightbour9(i, j, k, n, q, rho1);
						valueRho1=neightbour9(i, j, k, n, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							cIn0[i][j][k]=cIn0[i][j][k]+neightbour9(i, j, k, n, q, cIn0);
							cIn1[i][j][k]=cIn1[i][j][k]+neightbour9(i, j, k, n, q, cIn1);
							cIn2[i][j][k]=cIn2[i][j][k]+neightbour9(i, j, k, n, q, cIn2);
							cIn3[i][j][k]=cIn3[i][j][k]+neightbour9(i, j, k, n, q, cIn3);
							cIn4[i][j][k]=cIn4[i][j][k]+neightbour9(i, j, k, n, q, cIn4);
							cIn5[i][j][k]=cIn5[i][j][k]+neightbour9(i, j, k, n, q, cIn5);
							cIn6[i][j][k]=cIn6[i][j][k]+neightbour9(i, j, k, n, q, cIn6);
							cIn7[i][j][k]=cIn7[i][j][k]+neightbour9(i, j, k, n, q, cIn7);
							cIn8[i][j][k]=cIn8[i][j][k]+neightbour9(i, j, k, n, q, cIn8);
							cIn9[i][j][k]=cIn9[i][j][k]+neightbour9(i, j, k, n, q, cIn9);
							cIn10[i][j][k]=cIn10[i][j][k]+neightbour9(i, j, k, n, q, cIn10);
							cIn11[i][j][k]=cIn11[i][j][k]+neightbour9(i, j, k, n, q, cIn11);
							cIn12[i][j][k]=cIn12[i][j][k]+neightbour9(i, j, k, n, q, cIn12);
							cIn13[i][j][k]=cIn13[i][j][k]+neightbour9(i, j, k, n, q, cIn13);
							cIn14[i][j][k]=cIn14[i][j][k]+neightbour9(i, j, k, n, q, cIn14);
							cIn15[i][j][k]=cIn15[i][j][k]+neightbour9(i, j, k, n, q, cIn15);
							cIn16[i][j][k]=cIn16[i][j][k]+neightbour9(i, j, k, n, q, cIn16);
							cIn17[i][j][k]=cIn17[i][j][k]+neightbour9(i, j, k, n, q, cIn17);
							cIn18[i][j][k]=cIn18[i][j][k]+neightbour9(i, j, k, n, q, cIn18);
							Concentration[i][j][k]=Concentration[i][j][k]+neightbour9(i, j, k, n, q, Concentration);
						}
						valueRho=neightbour16(i, j, k, m, n, rho1);
						valueRho1=neightbour16(i, j, k, m, n, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							cIn0[i][j][k]=cIn0[i][j][k]+neightbour16(i, j, k, m, n, cIn0);
							cIn1[i][j][k]=cIn1[i][j][k]+neightbour16(i, j, k, m, n, cIn1);
							cIn2[i][j][k]=cIn2[i][j][k]+neightbour16(i, j, k, m, n, cIn2);
							cIn3[i][j][k]=cIn3[i][j][k]+neightbour16(i, j, k, m, n, cIn3);
							cIn4[i][j][k]=cIn4[i][j][k]+neightbour16(i, j, k, m, n, cIn4);
							cIn5[i][j][k]=cIn5[i][j][k]+neightbour16(i, j, k, m, n, cIn5);
							cIn6[i][j][k]=cIn6[i][j][k]+neightbour16(i, j, k, m, n, cIn6);
							cIn7[i][j][k]=cIn7[i][j][k]+neightbour16(i, j, k, m, n, cIn7);
							cIn8[i][j][k]=cIn8[i][j][k]+neightbour16(i, j, k, m, n, cIn8);
							cIn9[i][j][k]=cIn9[i][j][k]+neightbour16(i, j, k, m, n, cIn9);
							cIn10[i][j][k]=cIn10[i][j][k]+neightbour16(i, j, k, m, n, cIn10);
							cIn11[i][j][k]=cIn11[i][j][k]+neightbour16(i, j, k, m, n, cIn11);
							cIn12[i][j][k]=cIn12[i][j][k]+neightbour16(i, j, k, m, n, cIn12);
							cIn13[i][j][k]=cIn13[i][j][k]+neightbour16(i, j, k, m, n, cIn13);
							cIn14[i][j][k]=cIn14[i][j][k]+neightbour16(i, j, k, m, n, cIn14);
							cIn15[i][j][k]=cIn15[i][j][k]+neightbour16(i, j, k, m, n, cIn15);
							cIn16[i][j][k]=cIn16[i][j][k]+neightbour16(i, j, k, m, n, cIn16);
							cIn17[i][j][k]=cIn17[i][j][k]+neightbour16(i, j, k, m, n, cIn17);
							cIn18[i][j][k]=cIn18[i][j][k]+neightbour16(i, j, k, m, n, cIn18);
							Concentration[i][j][k]=Concentration[i][j][k]+neightbour16(i, j, k, m, n, Concentration);
						}
						valueRho=neightbour17(i, j, k, m, q, rho1);
						valueRho1=neightbour17(i, j, k, m, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							cIn0[i][j][k]=cIn0[i][j][k]+neightbour17(i, j, k, m, q, cIn0);
							cIn1[i][j][k]=cIn1[i][j][k]+neightbour17(i, j, k, m, q, cIn1);
							cIn2[i][j][k]=cIn2[i][j][k]+neightbour17(i, j, k, m, q, cIn2);
							cIn3[i][j][k]=cIn3[i][j][k]+neightbour17(i, j, k, m, q, cIn3);
							cIn4[i][j][k]=cIn4[i][j][k]+neightbour17(i, j, k, m, q, cIn4);
							cIn5[i][j][k]=cIn5[i][j][k]+neightbour17(i, j, k, m, q, cIn5);
							cIn6[i][j][k]=cIn6[i][j][k]+neightbour17(i, j, k, m, q, cIn6);
							cIn7[i][j][k]=cIn7[i][j][k]+neightbour17(i, j, k, m, q, cIn7);
							cIn8[i][j][k]=cIn8[i][j][k]+neightbour17(i, j, k, m, q, cIn8);
							cIn9[i][j][k]=cIn9[i][j][k]+neightbour17(i, j, k, m, q, cIn9);
							cIn10[i][j][k]=cIn10[i][j][k]+neightbour17(i, j, k, m, q, cIn10);
							cIn11[i][j][k]=cIn11[i][j][k]+neightbour17(i, j, k, m, q, cIn11);
							cIn12[i][j][k]=cIn12[i][j][k]+neightbour17(i, j, k, m, q, cIn12);
							cIn13[i][j][k]=cIn13[i][j][k]+neightbour17(i, j, k, m, q, cIn13);
							cIn14[i][j][k]=cIn14[i][j][k]+neightbour17(i, j, k, m, q, cIn14);
							cIn15[i][j][k]=cIn15[i][j][k]+neightbour17(i, j, k, m, q, cIn15);
							cIn16[i][j][k]=cIn16[i][j][k]+neightbour17(i, j, k, m, q, cIn16);
							cIn17[i][j][k]=cIn17[i][j][k]+neightbour17(i, j, k, m, q, cIn17);
							cIn18[i][j][k]=cIn18[i][j][k]+neightbour17(i, j, k, m, q, cIn18);
							Concentration[i][j][k]=Concentration[i][j][k]+neightbour17(i, j, k, m, q, Concentration);
						}
						valueRho=neightbour18(i, j, k, n, q, rho1);
						valueRho1=neightbour18(i, j, k, n, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							cIn0[i][j][k]=cIn0[i][j][k]+neightbour18(i, j, k, n, q, cIn0);
							cIn1[i][j][k]=cIn1[i][j][k]+neightbour18(i, j, k, n, q, cIn1);
							cIn2[i][j][k]=cIn2[i][j][k]+neightbour18(i, j, k, n, q, cIn2);
							cIn3[i][j][k]=cIn3[i][j][k]+neightbour18(i, j, k, n, q, cIn3);
							cIn4[i][j][k]=cIn4[i][j][k]+neightbour18(i, j, k, n, q, cIn4);
							cIn5[i][j][k]=cIn5[i][j][k]+neightbour18(i, j, k, n, q, cIn5);
							cIn6[i][j][k]=cIn6[i][j][k]+neightbour18(i, j, k, n, q, cIn6);
							cIn7[i][j][k]=cIn7[i][j][k]+neightbour18(i, j, k, n, q, cIn7);
							cIn8[i][j][k]=cIn8[i][j][k]+neightbour18(i, j, k, n, q, cIn8);
							cIn9[i][j][k]=cIn9[i][j][k]+neightbour18(i, j, k, n, q, cIn9);
							cIn10[i][j][k]=cIn10[i][j][k]+neightbour18(i, j, k, n, q, cIn10);
							cIn11[i][j][k]=cIn11[i][j][k]+neightbour18(i, j, k, n, q, cIn11);
							cIn12[i][j][k]=cIn12[i][j][k]+neightbour18(i, j, k, n, q, cIn12);
							cIn13[i][j][k]=cIn13[i][j][k]+neightbour18(i, j, k, n, q, cIn13);
							cIn14[i][j][k]=cIn14[i][j][k]+neightbour18(i, j, k, n, q, cIn14);
							cIn15[i][j][k]=cIn15[i][j][k]+neightbour18(i, j, k, n, q, cIn15);
							cIn16[i][j][k]=cIn16[i][j][k]+neightbour18(i, j, k, n, q, cIn16);
							cIn17[i][j][k]=cIn17[i][j][k]+neightbour18(i, j, k, n, q, cIn17);
							cIn18[i][j][k]=cIn18[i][j][k]+neightbour18(i, j, k, n, q, cIn18);
							Concentration[i][j][k]=Concentration[i][j][k]+neightbour18(i, j, k, n, q, Concentration);
						}
						valueRho=neightbour13(i, j, k, m, n, rho1);
						valueRho1=neightbour13(i, j, k, m, n, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							cIn0[i][j][k]=cIn0[i][j][k]+neightbour13(i, j, k, m, n, cIn0);
							cIn1[i][j][k]=cIn1[i][j][k]+neightbour13(i, j, k, m, n, cIn1);
							cIn2[i][j][k]=cIn2[i][j][k]+neightbour13(i, j, k, m, n, cIn2);
							cIn3[i][j][k]=cIn3[i][j][k]+neightbour13(i, j, k, m, n, cIn3);
							cIn4[i][j][k]=cIn4[i][j][k]+neightbour13(i, j, k, m, n, cIn4);
							cIn5[i][j][k]=cIn5[i][j][k]+neightbour13(i, j, k, m, n, cIn5);
							cIn6[i][j][k]=cIn6[i][j][k]+neightbour13(i, j, k, m, n, cIn6);
							cIn7[i][j][k]=cIn7[i][j][k]+neightbour13(i, j, k, m, n, cIn7);
							cIn8[i][j][k]=cIn8[i][j][k]+neightbour13(i, j, k, m, n, cIn8);
							cIn9[i][j][k]=cIn9[i][j][k]+neightbour13(i, j, k, m, n, cIn9);
							cIn10[i][j][k]=cIn10[i][j][k]+neightbour13(i, j, k, m, n, cIn10);
							cIn11[i][j][k]=cIn11[i][j][k]+neightbour13(i, j, k, m, n, cIn11);
							cIn12[i][j][k]=cIn12[i][j][k]+neightbour13(i, j, k, m, n, cIn12);
							cIn13[i][j][k]=cIn13[i][j][k]+neightbour13(i, j, k, m, n, cIn13);
							cIn14[i][j][k]=cIn14[i][j][k]+neightbour13(i, j, k, m, n, cIn14);
							cIn15[i][j][k]=cIn15[i][j][k]+neightbour13(i, j, k, m, n, cIn15);
							cIn16[i][j][k]=cIn16[i][j][k]+neightbour13(i, j, k, m, n, cIn16);
							cIn17[i][j][k]=cIn17[i][j][k]+neightbour13(i, j, k, m, n, cIn17);
							cIn18[i][j][k]=cIn18[i][j][k]+neightbour13(i, j, k, m, n, cIn18);
							Concentration[i][j][k]=Concentration[i][j][k]+neightbour13(i, j, k, m, n, Concentration);
						}
						valueRho=neightbour14(i, j, k, m, q, rho1);
						valueRho1=neightbour14(i, j, k, m, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							cIn0[i][j][k]=cIn0[i][j][k]+neightbour14(i, j, k, m, q, cIn0);
							cIn1[i][j][k]=cIn1[i][j][k]+neightbour14(i, j, k, m, q, cIn1);
							cIn2[i][j][k]=cIn2[i][j][k]+neightbour14(i, j, k, m, q, cIn2);
							cIn3[i][j][k]=cIn3[i][j][k]+neightbour14(i, j, k, m, q, cIn3);
							cIn4[i][j][k]=cIn4[i][j][k]+neightbour14(i, j, k, m, q, cIn4);
							cIn5[i][j][k]=cIn5[i][j][k]+neightbour14(i, j, k, m, q, cIn5);
							cIn6[i][j][k]=cIn6[i][j][k]+neightbour14(i, j, k, m, q, cIn6);
							cIn7[i][j][k]=cIn7[i][j][k]+neightbour14(i, j, k, m, q, cIn7);
							cIn8[i][j][k]=cIn8[i][j][k]+neightbour14(i, j, k, m, q, cIn8);
							cIn9[i][j][k]=cIn9[i][j][k]+neightbour14(i, j, k, m, q, cIn9);
							cIn10[i][j][k]=cIn10[i][j][k]+neightbour14(i, j, k, m, q, cIn10);
							cIn11[i][j][k]=cIn11[i][j][k]+neightbour14(i, j, k, m, q, cIn11);
							cIn12[i][j][k]=cIn12[i][j][k]+neightbour14(i, j, k, m, q, cIn12);
							cIn13[i][j][k]=cIn13[i][j][k]+neightbour14(i, j, k, m, q, cIn13);
							cIn14[i][j][k]=cIn14[i][j][k]+neightbour14(i, j, k, m, q, cIn14);
							cIn15[i][j][k]=cIn15[i][j][k]+neightbour14(i, j, k, m, q, cIn15);
							cIn16[i][j][k]=cIn16[i][j][k]+neightbour14(i, j, k, m, q, cIn16);
							cIn17[i][j][k]=cIn17[i][j][k]+neightbour14(i, j, k, m, q, cIn17);
							cIn18[i][j][k]=cIn18[i][j][k]+neightbour14(i, j, k, m, q, cIn18);
							Concentration[i][j][k]=Concentration[i][j][k]+neightbour14(i, j, k, m, q, Concentration);
						}
						valueRho=neightbour15(i, j, k, n, q, rho1);
						valueRho1=neightbour15(i, j, k, n, q, rho2);
						if ((valueRho<gasCritical&&valueRho>0)&&(valueRho1<gasCritical&&valueRho1>0)){
							number=number+1;
							cIn0[i][j][k]=cIn0[i][j][k]+neightbour15(i, j, k, n, q, cIn0);
							cIn1[i][j][k]=cIn1[i][j][k]+neightbour15(i, j, k, n, q, cIn1);
							cIn2[i][j][k]=cIn2[i][j][k]+neightbour15(i, j, k, n, q, cIn2);
							cIn3[i][j][k]=cIn3[i][j][k]+neightbour15(i, j, k, n, q, cIn3);
							cIn4[i][j][k]=cIn4[i][j][k]+neightbour15(i, j, k, n, q, cIn4);
							cIn5[i][j][k]=cIn5[i][j][k]+neightbour15(i, j, k, n, q, cIn5);
							cIn6[i][j][k]=cIn6[i][j][k]+neightbour15(i, j, k, n, q, cIn6);
							cIn7[i][j][k]=cIn7[i][j][k]+neightbour15(i, j, k, n, q, cIn7);
							cIn8[i][j][k]=cIn8[i][j][k]+neightbour15(i, j, k, n, q, cIn8);
							cIn9[i][j][k]=cIn9[i][j][k]+neightbour15(i, j, k, n, q, cIn9);
							cIn10[i][j][k]=cIn10[i][j][k]+neightbour15(i, j, k, n, q, cIn10);
							cIn11[i][j][k]=cIn11[i][j][k]+neightbour15(i, j, k, n, q, cIn11);
							cIn12[i][j][k]=cIn12[i][j][k]+neightbour15(i, j, k, n, q, cIn12);
							cIn13[i][j][k]=cIn13[i][j][k]+neightbour15(i, j, k, n, q, cIn13);
							cIn14[i][j][k]=cIn14[i][j][k]+neightbour15(i, j, k, n, q, cIn14);
							cIn15[i][j][k]=cIn15[i][j][k]+neightbour15(i, j, k, n, q, cIn15);
							cIn16[i][j][k]=cIn16[i][j][k]+neightbour15(i, j, k, n, q, cIn16);
							cIn17[i][j][k]=cIn17[i][j][k]+neightbour15(i, j, k, n, q, cIn17);
							cIn18[i][j][k]=cIn18[i][j][k]+neightbour15(i, j, k, n, q, cIn18);
							Concentration[i][j][k]=Concentration[i][j][k]+neightbour15(i, j, k, n, q, Concentration);
						}						
						if (number==0){
							Concentration[i][j][k]=0.0;
						}
						else{
							cIn0[i][j][k]=cIn0[i][j][k]/number;
							cIn1[i][j][k]=cIn1[i][j][k]/number;
							cIn2[i][j][k]=cIn2[i][j][k]/number;
							cIn3[i][j][k]=cIn3[i][j][k]/number;
							cIn4[i][j][k]=cIn4[i][j][k]/number;
							cIn5[i][j][k]=cIn5[i][j][k]/number;
							cIn6[i][j][k]=cIn6[i][j][k]/number;
							cIn7[i][j][k]=cIn7[i][j][k]/number;
							cIn8[i][j][k]=cIn8[i][j][k]/number;
							cIn9[i][j][k]=cIn9[i][j][k]/number;
							cIn10[i][j][k]=cIn10[i][j][k]/number;
							cIn11[i][j][k]=cIn11[i][j][k]/number;
							cIn12[i][j][k]=cIn12[i][j][k]/number;
							cIn13[i][j][k]=cIn13[i][j][k]/number;
							cIn14[i][j][k]=cIn14[i][j][k]/number;
							cIn15[i][j][k]=cIn15[i][j][k]/number;
							cIn16[i][j][k]=cIn16[i][j][k]/number;
							cIn17[i][j][k]=cIn17[i][j][k]/number;
							cIn18[i][j][k]=cIn18[i][j][k]/number;
							Concentration[i][j][k]=Concentration[i][j][k]/number;
						}
					}
				}
			}
		}
	}
}

void lastConcentration(int Length, int Width, int Height, double gasCritical){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i=0;i<Length;i++){
			for (j=0;j<Width;j++){
				for (k=0;k<Height;k++){
						wallNearConcentration[i][j][k]=Concentration[i][j][k];						
						wallNearConcentrationFirst[i][j][k] = ConcentrationFirst[i][j][k];
				}
			}
		}
	}
}

void additionalTermCalculationConcentrationD3Q7(int Length, int Width, int Height, double gasCritical, double ***Ux, double ***Uy, double ***Uz,double constant){
	int i;
	int j;
	int k;
	double concentrationNeightbour1;
	double concentrationNeightbour2;
	double concentrationNeightbour3;
	double concentrationNeightbour4;
	double concentrationNeightbour5;
	double concentrationNeightbour6;
#pragma omp parallel private(i,j,k,concentrationNeightbour1,concentrationNeightbour2,concentrationNeightbour3,concentrationNeightbour4,concentrationNeightbour5,concentrationNeightbour6)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (j = 0; j<Width; j++){
				for (k = 0; k<Height; k++){
					if (domain[i][j][k] == 0 && rho[i][j][k]>gasCritical&&domainConnection[i][j][k] == 1){
						concentrationNeightbour1 = neightbour1(i, j, k, Length, Concentration);
						concentrationNeightbour2 = neightbour2(i, j, k, Width, Concentration);
						concentrationNeightbour3 = neightbour3(i, j, k, Height, Concentration);
						concentrationNeightbour4 = neightbour4(i, j, k, Length, Concentration);
						concentrationNeightbour5 = neightbour5(i, j, k, Width, Concentration);
						concentrationNeightbour6 = neightbour6(i, j, k, Height, Concentration);

						if (Domain1[i][j][k] != 0 || rho1[i][j][k] <= gasCritical){
							concentrationNeightbour4 = Concentration[i][j][k];
						}
						if (Domain4[i][j][k] != 0 || rho4[i][j][k] <= gasCritical){
							concentrationNeightbour1 = Concentration[i][j][k];
						}
						if ((Domain1[i][j][k] != 0 || rho1[i][j][k] <= gasCritical) && (Domain4[i][j][k] != 0 || rho4[i][j][k] <= gasCritical)){
							concentrationNeightbour4 = Concentration[i][j][k];
							concentrationNeightbour1 = Concentration[i][j][k];
						}

						if (Domain2[i][j][k] != 0 || rho2[i][j][k] <= gasCritical){
							concentrationNeightbour5 = Concentration[i][j][k];
						}
						if (Domain5[i][j][k] != 0 || rho5[i][j][k] <= gasCritical){
							concentrationNeightbour2 = Concentration[i][j][k];
						}
						if ((Domain2[i][j][k] != 0 || rho2[i][j][k] <= gasCritical) && (Domain5[i][j][k] != 0 || rho5[i][j][k] <= gasCritical)){
							concentrationNeightbour5 = Concentration[i][j][k];
							concentrationNeightbour2 = Concentration[i][j][k];
						}

						if (Domain3[i][j][k] != 0 || rho3[i][j][k] <= gasCritical){
							concentrationNeightbour6 = Concentration[i][j][k];
						}
						if (Domain6[i][j][k] != 0 || rho6[i][j][k] <= gasCritical){
							concentrationNeightbour3 = Concentration[i][j][k];
						}
						if ((Domain3[i][j][k] != 0 || rho3[i][j][k] <= gasCritical) && (Domain6[i][j][k] != 0 || rho6[i][j][k] <= gasCritical)){
							concentrationNeightbour6 = Concentration[i][j][k];
							concentrationNeightbour3 = Concentration[i][j][k];
						}

						if (k == 0){
							concentrationNeightbour6 = Concentration[i][j][k];
						}
						if (k == Height - 1){
							concentrationNeightbour3 = Concentration[i][j][k];
						}
						convectTermConcentration[i][j][k] = (-Ux[i][j][k] * 0.5*(concentrationNeightbour1 - concentrationNeightbour4) - Uy[i][j][k] * 0.5*(concentrationNeightbour2 - concentrationNeightbour5) - Uz[i][j][k] * 0.5*(concentrationNeightbour3 - concentrationNeightbour6))*constant;
					}
				}
			}
		}
	}
}

void additionalTermCalculationConcentrationFirstD3Q7(int Length, int Width, int Height, double gasCritical, double ***Ux, double ***Uy, double ***Uz, double constant){
	int i;
	int j;
	int k;
	double concentrationNeightbour1;
	double concentrationNeightbour2;
	double concentrationNeightbour3;
	double concentrationNeightbour4;
	double concentrationNeightbour5;
	double concentrationNeightbour6;
#pragma omp parallel private(i,j,k,concentrationNeightbour1,concentrationNeightbour2,concentrationNeightbour3,concentrationNeightbour4,concentrationNeightbour5,concentrationNeightbour6)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (j = 0; j<Width; j++){
				for (k = 0; k<Height; k++){
					if (domain[i][j][k] == 0 && rho[i][j][k]>gasCritical&&domainConnection[i][j][k] == 1){
						concentrationNeightbour1 = neightbour1(i, j, k, Length, ConcentrationFirst);
						concentrationNeightbour2 = neightbour2(i, j, k, Width, ConcentrationFirst);
						concentrationNeightbour3 = neightbour3(i, j, k, Height, ConcentrationFirst);
						concentrationNeightbour4 = neightbour4(i, j, k, Length, ConcentrationFirst);
						concentrationNeightbour5 = neightbour5(i, j, k, Width, ConcentrationFirst);
						concentrationNeightbour6 = neightbour6(i, j, k, Height, ConcentrationFirst);

						if (Domain1[i][j][k] != 0 || rho1[i][j][k] <= gasCritical){
							concentrationNeightbour4 = ConcentrationFirst[i][j][k];
						}
						if (Domain4[i][j][k] != 0 || rho4[i][j][k] <= gasCritical){
							concentrationNeightbour1 = ConcentrationFirst[i][j][k];
						}
						if ((Domain1[i][j][k] != 0 || rho1[i][j][k] <= gasCritical) && (Domain4[i][j][k] != 0 || rho4[i][j][k] <= gasCritical)){
							concentrationNeightbour4 = ConcentrationFirst[i][j][k];
							concentrationNeightbour1 = ConcentrationFirst[i][j][k];
						}

						if (Domain2[i][j][k] != 0 || rho2[i][j][k] <= gasCritical){
							concentrationNeightbour5 = ConcentrationFirst[i][j][k];
						}
						if (Domain5[i][j][k] != 0 || rho5[i][j][k] <= gasCritical){
							concentrationNeightbour2 = ConcentrationFirst[i][j][k];
						}
						if ((Domain2[i][j][k] != 0 || rho2[i][j][k] <= gasCritical) && (Domain5[i][j][k] != 0 || rho5[i][j][k] <= gasCritical)){
							concentrationNeightbour5 = ConcentrationFirst[i][j][k];
							concentrationNeightbour2 = ConcentrationFirst[i][j][k];
						}

						if (Domain3[i][j][k] != 0 || rho3[i][j][k] <= gasCritical){
							concentrationNeightbour6 = ConcentrationFirst[i][j][k];
						}
						if (Domain6[i][j][k] != 0 || rho6[i][j][k] <= gasCritical){
							concentrationNeightbour3 = ConcentrationFirst[i][j][k];
						}
						if ((Domain3[i][j][k] != 0 || rho3[i][j][k] <= gasCritical) && (Domain6[i][j][k] != 0 || rho6[i][j][k] <= gasCritical)){
							concentrationNeightbour6 = ConcentrationFirst[i][j][k];
							concentrationNeightbour3 = ConcentrationFirst[i][j][k];
						}

						if (k == 0){
							concentrationNeightbour6 = ConcentrationFirst[i][j][k];
						}
						if (k == Height - 1){
							concentrationNeightbour3 = ConcentrationFirst[i][j][k];
						}
						convectTermConcentrationFirst[i][j][k] = (-Ux[i][j][k] * 0.5*(concentrationNeightbour1 - concentrationNeightbour4) - Uy[i][j][k] * 0.5*(concentrationNeightbour2 - concentrationNeightbour5) - Uz[i][j][k] * 0.5*(concentrationNeightbour3 - concentrationNeightbour6))*constant;
					}
				}
			}
		}
	}
}

void additionalTermCalculationConcentrationSecondD3Q7(int Length, int Width, int Height, double gasCritical, double ***Ux, double ***Uy, double ***Uz,double constant){
	int i;
	int j;
	int k;
	double concentrationNeightbour1;
	double concentrationNeightbour2;
	double concentrationNeightbour3;
	double concentrationNeightbour4;
	double concentrationNeightbour5;
	double concentrationNeightbour6;
#pragma omp parallel private(i,j,k,concentrationNeightbour1,concentrationNeightbour2,concentrationNeightbour3,concentrationNeightbour4,concentrationNeightbour5,concentrationNeightbour6)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (j = 0; j<Width; j++){
				for (k = 0; k<Height; k++){
					if (domain[i][j][k] == 0 && rho[i][j][k]>gasCritical&&domainConnection[i][j][k] == 1){
						concentrationNeightbour1 = neightbour1(i, j, k, Length, ConcentrationSecond);
						concentrationNeightbour2 = neightbour2(i, j, k, Width, ConcentrationSecond);
						concentrationNeightbour3 = neightbour3(i, j, k, Height, ConcentrationSecond);
						concentrationNeightbour4 = neightbour4(i, j, k, Length, ConcentrationSecond);
						concentrationNeightbour5 = neightbour5(i, j, k, Width, ConcentrationSecond);
						concentrationNeightbour6 = neightbour6(i, j, k, Height, ConcentrationSecond);

						if (Domain1[i][j][k] != 0 || rho1[i][j][k] <= gasCritical){
							concentrationNeightbour4 = ConcentrationSecond[i][j][k];
						}
						if (Domain4[i][j][k] != 0 || rho4[i][j][k] <= gasCritical){
							concentrationNeightbour1 = ConcentrationSecond[i][j][k];
						}
						if ((Domain1[i][j][k] != 0 || rho1[i][j][k] <= gasCritical) && (Domain4[i][j][k] != 0 || rho4[i][j][k] <= gasCritical)){
							concentrationNeightbour4 = ConcentrationSecond[i][j][k];
							concentrationNeightbour1 = ConcentrationSecond[i][j][k];
						}

						if (Domain2[i][j][k] != 0 || rho2[i][j][k] <= gasCritical){
							concentrationNeightbour5 = ConcentrationSecond[i][j][k];
						}
						if (Domain5[i][j][k] != 0 || rho5[i][j][k] <= gasCritical){
							concentrationNeightbour2 = ConcentrationSecond[i][j][k];
						}
						if ((Domain2[i][j][k] != 0 || rho2[i][j][k] <= gasCritical) && (Domain5[i][j][k] != 0 || rho5[i][j][k] <= gasCritical)){
							concentrationNeightbour5 = ConcentrationSecond[i][j][k];
							concentrationNeightbour2 = ConcentrationSecond[i][j][k];
						}

						if (Domain3[i][j][k] != 0 || rho3[i][j][k] <= gasCritical){
							concentrationNeightbour6 = ConcentrationSecond[i][j][k];
						}
						if (Domain6[i][j][k] != 0 || rho6[i][j][k] <= gasCritical){
							concentrationNeightbour3 = ConcentrationSecond[i][j][k];
						}
						if ((Domain3[i][j][k] != 0 || rho3[i][j][k] <= gasCritical) && (Domain6[i][j][k] != 0 || rho6[i][j][k] <= gasCritical)){
							concentrationNeightbour6 = ConcentrationSecond[i][j][k];
							concentrationNeightbour3 = ConcentrationSecond[i][j][k];
						}

						if (k == 0){
							concentrationNeightbour6 = ConcentrationSecond[i][j][k];
						}
						if (k == Height - 1){
							concentrationNeightbour3 = ConcentrationSecond[i][j][k];
						}
						convectTermConcentrationSecond[i][j][k] = (-Ux[i][j][k] * 0.5*(concentrationNeightbour1 - concentrationNeightbour4) - Uy[i][j][k] * 0.5*(concentrationNeightbour2 - concentrationNeightbour5) - Uz[i][j][k] * 0.5*(concentrationNeightbour3 - concentrationNeightbour6))*constant;
					}
				}
			}
		}
	}
}

void additionalDiffusionTermCalculationD3Q7(int Length, int Width, int Height, double constant1, double constant2, double gasCritical){
	int i;
	int j;
	int k;
	double concentrationNeightbour1;
	double concentrationFirstNeightbour1;
	double concentrationNeightbour2;
	double concentrationFirstNeightbour2;
	double concentrationNeightbour3;
	double concentrationFirstNeightbour3;
	double concentrationNeightbour4;
	double concentrationFirstNeightbour4;
	double concentrationNeightbour5;
	double concentrationFirstNeightbour5;
	double concentrationNeightbour6;
	double concentrationFirstNeightbour6;
	double concentrationAdditionalTerm;
	double concentrationFirstAdditionalTerm;
#pragma omp parallel private(i,j,k,concentrationNeightbour1,concentrationFirstNeightbour1,concentrationNeightbour2,concentrationFirstNeightbour2,concentrationNeightbour3,concentrationFirstNeightbour3,concentrationNeightbour4,concentrationFirstNeightbour4,concentrationNeightbour5,concentrationFirstNeightbour5,concentrationNeightbour6,concentrationFirstNeightbour6,concentrationAdditionalTerm,concentrationFirstAdditionalTerm)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (j = 0; j<Width; j++){
				for (k = 0; k<Height; k++){
					if (domain[i][j][k] == 0 && rho[i][j][k]>gasCritical&&domainConnection[i][j][k] == 1){
						concentrationNeightbour1 = neightbour1(i, j, k, Length, Concentration);
						concentrationFirstNeightbour1 = neightbour1(i, j, k, Length, ConcentrationFirst);

						concentrationNeightbour2 = neightbour2(i, j, k, Width, Concentration);
						concentrationFirstNeightbour2 = neightbour2(i, j, k, Width, ConcentrationFirst);

						concentrationNeightbour3 = neightbour3(i, j, k, Height, Concentration);
						concentrationFirstNeightbour3 = neightbour3(i, j, k, Height, ConcentrationFirst);

						concentrationNeightbour4 = neightbour4(i, j, k, Length, Concentration);
						concentrationFirstNeightbour4 = neightbour4(i, j, k, Length, ConcentrationFirst);

						concentrationNeightbour5 = neightbour5(i, j, k, Width, Concentration);
						concentrationFirstNeightbour5 = neightbour5(i, j, k, Width, ConcentrationFirst);

						concentrationNeightbour6 = neightbour6(i, j, k, Height, Concentration);
						concentrationFirstNeightbour6 = neightbour6(i, j, k, Height, ConcentrationFirst);

						if (Domain1[i][j][k] != 0 || rho1[i][j][k] <= gasCritical){
							concentrationNeightbour4 = Concentration[i][j][k];
							concentrationFirstNeightbour4 = ConcentrationFirst[i][j][k];
						}
						if (Domain4[i][j][k] != 0 || rho4[i][j][k] <= gasCritical){
							concentrationNeightbour1 = Concentration[i][j][k];
							concentrationFirstNeightbour1 = ConcentrationFirst[i][j][k];
						}
						if ((Domain1[i][j][k] != 0 || rho1[i][j][k] <= gasCritical) && (Domain4[i][j][k] != 0 || rho4[i][j][k] <= gasCritical)){
							concentrationNeightbour4 = Concentration[i][j][k];
							concentrationFirstNeightbour4 = ConcentrationFirst[i][j][k];
							concentrationNeightbour1 = Concentration[i][j][k];
							concentrationFirstNeightbour1 = ConcentrationFirst[i][j][k];
						}

						if (Domain2[i][j][k] != 0 || rho2[i][j][k] <= gasCritical){
							concentrationNeightbour5 = Concentration[i][j][k];
							concentrationFirstNeightbour5 = ConcentrationFirst[i][j][k];
						}
						if (Domain5[i][j][k] != 0 || rho5[i][j][k] <= gasCritical){
							concentrationNeightbour2 = Concentration[i][j][k];
							concentrationFirstNeightbour2 = ConcentrationFirst[i][j][k];
						}
						if ((Domain2[i][j][k] != 0 || rho2[i][j][k] <= gasCritical) && (Domain5[i][j][k] != 0 || rho5[i][j][k] <= gasCritical)){
							concentrationNeightbour5 = Concentration[i][j][k];
							concentrationFirstNeightbour5 = ConcentrationFirst[i][j][k];
							concentrationNeightbour2 = Concentration[i][j][k];
							concentrationFirstNeightbour2 = ConcentrationFirst[i][j][k];
						}

						if (Domain3[i][j][k] != 0 || rho3[i][j][k] <= gasCritical){
							concentrationNeightbour6 = Concentration[i][j][k];
							concentrationFirstNeightbour6 = ConcentrationFirst[i][j][k];
						}
						if (Domain6[i][j][k] != 0 || rho6[i][j][k] <= gasCritical){
							concentrationNeightbour3 = Concentration[i][j][k];
							concentrationFirstNeightbour3 = ConcentrationFirst[i][j][k];
						}
						if ((Domain3[i][j][k] != 0 || rho3[i][j][k] <= gasCritical) && (Domain6[i][j][k] != 0 || rho6[i][j][k] <= gasCritical)){
							concentrationNeightbour6 = Concentration[i][j][k];
							concentrationFirstNeightbour6 = ConcentrationFirst[i][j][k];
							concentrationNeightbour3 = Concentration[i][j][k];
							concentrationFirstNeightbour3 = ConcentrationFirst[i][j][k];
						}

						if (k == 0){
							concentrationNeightbour6 = Concentration[i][j][k];
							concentrationFirstNeightbour6 = ConcentrationFirst[i][j][k];
						}
						if (k == Height - 1){
							concentrationNeightbour3 = Concentration[i][j][k];
							concentrationFirstNeightbour3 = ConcentrationFirst[i][j][k];
						}
						if (j == 0){
							concentrationNeightbour5 = Concentration[i][j][k];
							concentrationFirstNeightbour5 = ConcentrationFirst[i][j][k];
						}
						concentrationAdditionalTerm = (2.0*concentrationNeightbour1 + 2.0*concentrationNeightbour4 + 2.0*concentrationNeightbour2 + 2.0*concentrationNeightbour5 + 2.0*concentrationNeightbour3 + 2.0*concentrationNeightbour6 - 12.0*Concentration[i][j][k]) / 2.0;
						concentrationFirstAdditionalTerm = (2.0*concentrationFirstNeightbour1 + 2.0*concentrationFirstNeightbour4 + 2.0*concentrationFirstNeightbour2 + 2.0*concentrationFirstNeightbour5 + 2.0*concentrationFirstNeightbour3 + 2.0*concentrationFirstNeightbour6 - 12.0*ConcentrationFirst[i][j][k]) / 2.0;
						diffusionTermConcentration[i][j][k] = -constant1*concentrationAdditionalTerm;
						diffusionTermConcentrationFirst[i][j][k] = -constant2*concentrationFirstAdditionalTerm;
					}
				}
			}
		}
	}
}

void massTransMatrix(int Length, int Width, int Height, double a, double b){
	int i;
	int j;
	int k;
#pragma omp parallel private(i,j,k)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i<Length; i++){
			for (j = 0; j<Width; j++){
				for (k = 0; k<Height; k++){
					if (domain[i][j][k] == 0){
//						massTransCoeMatrix[i][j][k] = 1.0 / (1.0 / (1.0 + exp(-267000 * (pow((Ux[i][j][k] * Ux[i][j][k] + Uy[i][j][k] * Uy[i][j][k] + Uz[i][j][k] * Uz[i][j][k]), 0.5)) + 42.302)) + 0.235) / 0.00008511*a / b;
						massTransCoeMatrix[i][j][k] = 29447;
					}
					if (domain[i][j][k] != 0){
						massTransCoeMatrix[i][j][k] = 0.0;
					}
				}
			}
		}
	}
}
