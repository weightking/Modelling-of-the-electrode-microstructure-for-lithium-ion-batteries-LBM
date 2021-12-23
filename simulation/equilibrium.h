#ifndef _EQUILIBRIUM_H
#define _EQUILIBRIUM_H

extern double ***cuNS0;
extern double ***cuNS1;
extern double ***cuNS2;
extern double ***cuNS3;
extern double ***cuNS4;
extern double ***cuNS5;
extern double ***cuNS6;
extern double ***cuNS7;
extern double ***cuNS8;
extern double ***cuNS9;
extern double ***cuNS10;
extern double ***cuNS11;
extern double ***cuNS12;
extern double ***cuNS13;
extern double ***cuNS14;
extern double ***cuNS15;
extern double ***cuNS16;
extern double ***cuNS17;
extern double ***cuNS18;
extern double ***ux2uy2uz2;
extern double ***Eq0;
extern double ***Eq1;
extern double ***Eq2;
extern double ***Eq3;
extern double ***Eq4;
extern double ***Eq5;
extern double ***Eq6;
extern double ***Eq7;
extern double ***Eq8;
extern double ***Eq9;
extern double ***Eq10;
extern double ***Eq11;
extern double ***Eq12;
extern double ***Eq13;
extern double ***Eq14;
extern double ***Eq15;
extern double ***Eq16;
extern double ***Eq17;
extern double ***Eq18;

extern void equilibriumArrange(int Length, int Width, int Height);
extern void equilibriumCalculation(int Length, int Width, int Height, double ***Ux, double ***Uy, double ***Uz);
extern void equilibriumConcentrationCalculation(int Length, int Width, int Height, double ***Ux, double ***Uy, double ***Uz);
extern void MemoryFreeEquilibrium(int m, int n, int q);

#endif