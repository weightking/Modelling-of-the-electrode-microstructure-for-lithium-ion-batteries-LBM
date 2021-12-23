#ifndef _DOMAIN_H
#define _DOMAIN_H

extern double ***domainPorous;
extern double ***domainTwoPhase;
extern double ***domain;
extern double ***Domain1;
extern double ***Domain2;
extern double ***Domain3;
extern double ***Domain4;
extern double ***Domain5;
extern double ***Domain6;
extern double ***Domain7;
extern double ***Domain8;
extern double ***Domain9;
extern double ***Domain10;
extern double ***Domain11;
extern double ***Domain12;
extern double ***Domain13;
extern double ***Domain14;
extern double ***Domain15;
extern double ***Domain16;
extern double ***Domain17;
extern double ***Domain18;
extern double ***Domain19;
extern double ***Domain20;
extern double ***Domain21;
extern double ***Domain22;
extern double ***Domain23;
extern double ***Domain24;
extern double ***Domain25;
extern double ***Domain26;

extern void domainArrange(int m, int n, int q);

extern void domainShift(int m, int n, int q);

extern void domainDefine(int mInitial, int m, int nInitial, int n, int qInitial, int q, double phase);

extern void domainDefineInletHoles(int qInitial, int q, int Length, int Width, int number, int radius, double phase);

extern void domainFreeMemory(int m, int n, int q);
#endif