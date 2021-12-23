#ifndef _DIVERGENCE_H
#define _DEVERGENCE_H

extern double divergeNumberTwoPhase(int m, int n, int q, double FLUID);
extern void divergeFixTwoPhase(int m, int n, int q,double GAS, double FLUID, double b, double R, double T, double a);

#endif