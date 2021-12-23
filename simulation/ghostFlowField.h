#ifndef _GHOSTFLOWFIELD_H
#define _GHOSTFLOWFIELD_H

extern double ***fIn0Ghost;
extern double ***fIn1Ghost;
extern double ***fIn2Ghost;
extern double ***fIn3Ghost;
extern double ***fIn4Ghost;
extern double ***fIn5Ghost;
extern double ***fIn6Ghost;
extern double ***fIn7Ghost;
extern double ***fIn8Ghost;
extern double ***fIn9Ghost;
extern double ***fIn10Ghost;
extern double ***fIn11Ghost;
extern double ***fIn12Ghost;
extern double ***fIn13Ghost;
extern double ***fIn14Ghost;
extern double ***fIn15Ghost;
extern double ***fIn16Ghost;
extern double ***fIn17Ghost;
extern double ***fIn18Ghost;
extern double ***rhoGhost;
extern double ***uxGhost;
extern double ***uyGhost;
extern double ***uzGhost;

double AGhost;

extern void FieldArrangeGhostVelocity(int Length, int Width, int Height);
extern void FieldInitialGhostVelocity(int Length, int Width, int Height, double GAS, double gasCritical);
extern void CollisionSRTGhostVelocity(int Length,int Width,int Height,double gasCritical);
extern void CollisionMRTGhostVelocity(int m, int n, int q, double gasCritical);
extern void MemoryFreeGhostVelocity(int m, int n, int q);
extern void StreamGhostVelocity(int m, int n, int q);
extern void FieldCalculationGhostVelocity(int m, int n, int q, double gasCritical);
extern void refillGhostVelocity(int m, int n, int q, double gasCritical,double ***rho1, double ***rho2, double gasDensity);

#endif