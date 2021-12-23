#ifndef _WALLCONTACTANGLE_H
#define _WALLCONTACTANGLE_H

extern double ***gWall;

extern void gWallArrangeTwoPhase(int m, int n, int q);
extern void gWallDefineTwoPhase(int mInitial, int m, int nInitial, int n, int qInitial, int q, double contactAngle);

#endif