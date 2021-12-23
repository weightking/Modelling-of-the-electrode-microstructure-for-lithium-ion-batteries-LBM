#ifndef _ELECTRODEPOTENTIAL_H
#define _ELECTRODEPOTENTIAL_H

extern double ***electrodeIn0;
extern double ***electrodeIn1;
extern double ***electrodeIn2;
extern double ***electrodeIn3;
extern double ***electrodeIn4;
extern double ***electrodeIn5;
extern double ***electrodeIn6;
extern double ***electrodeIn7;
extern double ***electrodeIn8;
extern double ***electrodeIn9;
extern double ***electrodeIn10;
extern double ***electrodeIn11;
extern double ***electrodeIn12;
extern double ***electrodeIn13;
extern double ***electrodeIn14;
extern double ***electrodeIn15;
extern double ***electrodeIn16;
extern double ***electrodeIn17;
extern double ***electrodeIn18;
extern double ***electrodePotential;
extern double ***electrodeConnection;
extern double ***lastWallElectrode;
extern double ***localCurrent;
extern double ***reaction1;
extern double ***reaction2;
extern double ***reaction3;
extern double ***reaction4;
extern double ***reaction5;
extern double ***reaction6;

extern void FieldArrangeElectrodePotential(int Length, int Width, int Height);
extern void electrodeConnectionRegion(int Length, int Width, int Height);
extern void electrodeConnectionRegionD3Q7(int Length, int Width, int Height);
extern void FieldInitialelectrodePotential(int Length, int Width, int Height, double electrodePotentialInitial);
extern void FieldInitialelectrodePotentialD3Q7(int Length, int Width, int Initial, int Height, double electrodePotentialInitial);
extern void CollisionSRTElectrodePotential(int Length,int Width,int Height,double CrelaxationTime);
extern void CollisionSRTElectrodePotentialD3Q7(int Length,int Width,int Height,double CrelaxationTime);
extern void MemoryFreeElectrodePotential(int m, int n, int q);
extern void StreamElectrodePotential(int m, int n, int q);
extern void StreamElectrodePotentialD3Q7(int m, int n, int q);
extern void FieldCalculationElectrodePotential(int m, int n, int q);
extern void FieldCalculationElectrodePotentialD3Q7(int m, int n, int q);
extern void lastElectrode(int Length, int Width, int Height);
extern void reactionElectrodeD3Q7(int m, int n, int q, double gasCritical, double RateConstant, double constant, double K, double exchange, double MaxConcentrationdouble, double lastTotalCurrentResidual);

#endif