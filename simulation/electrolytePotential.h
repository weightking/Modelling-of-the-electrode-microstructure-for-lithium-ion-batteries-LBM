#ifndef _ELECTROLYTEPOTENTIAL_H
#define _ELECTROLYTEPOTENTIAL_H

extern double ***electrolyteIn0;
extern double ***electrolyteIn1;
extern double ***electrolyteIn2;
extern double ***electrolyteIn3;
extern double ***electrolyteIn4;
extern double ***electrolyteIn5;
extern double ***electrolyteIn6;
extern double ***electrolyteIn7;
extern double ***electrolyteIn8;
extern double ***electrolyteIn9;
extern double ***electrolyteIn10;
extern double ***electrolyteIn11;
extern double ***electrolyteIn12;
extern double ***electrolyteIn13;
extern double ***electrolyteIn14;
extern double ***electrolyteIn15;
extern double ***electrolyteIn16;
extern double ***electrolyteIn17;
extern double ***electrolyteIn18;
extern double ***electrolytePotential;
extern double ***electrolyteRelaxationTime;
extern double ***electrolyteAdditionalTerm;
extern double ***lastWallElectrolyte;

extern void FieldArrangeelectrolytePotential(int Length, int Width, int Height);
extern void FieldInitialelectrolytePotential(int Length, int Width, int Height, double electrolytePotentialInitial, double gasCritical);
extern void FieldInitialelectrolytePotentialD3Q7(int Length, int Width, int Initial, int Height, double electrolytePotentialInitial, double gasCritical);
extern void electrolyteRelaxationTimeCalculation(int Length, int Width, int Height, double constant, double constantMol, double constantLelectrolyteConductivity, double gasCritical);
extern void additionalTermCalculation(int Length,int Width,int Height,double constant,double gasCritical,double concentrationDiffusivity,double concentrationFirstDiffusivity,double concentrationSecondDiffusivity,double concentrationThirdDiffusivity,int charge,int charge1,int charge2,int charge3);
extern void additionalTermCalculationD3Q7(int Length,int Width,int Height,double constant1,double constant2,double K,double gasCritical,double constantMol,double constantLelectrolyteConductivity,double exchange);
extern void CollisionSRTelectrolytePotential(int Length,int Width,int Height,double gasCritical);
extern void CollisionSRTelectrolytePotentialD3Q7(int Length,int Width,int Height,double gasCritical);
extern void MemoryFreeelectrolytePotential(int m, int n, int q);
extern void StreamelectrolytePotential(int m, int n, int q);
extern void StreamelectrolytePotentialD3Q7(int m, int n, int q);
extern void FieldCalculationelectrolytePotential(int m, int n, int q, double gasCritical);
extern void FieldCalculationelectrolytePotentialD3Q7(int m, int n, int q, double gasCritical);
extern void lastElectrolyte(int Length, int Width, int Height, double gasCritical);
extern void reactionElectrolyteD3Q7(int m, int n, int q, double lastTotalCurrentResidual);

#endif