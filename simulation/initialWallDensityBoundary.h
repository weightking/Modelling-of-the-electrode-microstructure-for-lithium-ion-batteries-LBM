#ifndef _INITIALWALLDENSITYBOUNDARY_H
#define _INITIALWALLDENSITYBOUNDARY_H

extern void initialWallBoundaryTwoPhase(int mInitial,int m,int nInitial,int n,int qInitial,int q,double ***domain,
	double ***Domain1,double ***Domain2,double ***Domain3,double ***Domain4,double ***Domain5,double ***Domain6,double ***Domain7,double ***Domain8,double ***Domain9,double ***Domain10,
	double ***Domain11,double ***Domain12,double ***Domain13,double ***Domain14,double ***Domain15,double ***Domain16,double ***Domain17,double ***Domain18,double constant,double rhoLayer,
	double b, double R, double T, double a);

extern void initialWallNearestTwoPhase(int mInitial,int m,int nInitial,int n,int qInitial,int q,double ***domain,
	double ***Domain1,double ***Domain2,double ***Domain3,double ***Domain4,double ***Domain5,double ***Domain6,double ***Domain7,double ***Domain8,double ***Domain9,double ***Domain10,
	double ***Domain11,double ***Domain12,double ***Domain13,double ***Domain14,double ***Domain15,double ***Domain16,double ***Domain17,double ***Domain18,double constant,double rhoLayer,
	double b, double R, double T, double a, double GAS);

#endif