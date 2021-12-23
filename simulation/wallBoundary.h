#ifndef _WALLBOUNDARY_H
#define _WALLBOUNDARY_H

extern void wallBoundaryTwoPhase(int m, int n, int q);
extern void wallBoundaryConcentration(int m, int n, int q, double gasCritical);
extern void wallBoundaryConcentrationD3Q7(int m, int n, int q);
extern void wallBoundaryConcentrationD3Q15(int m, int n, int q, double gasCritical);
extern void wallBoundaryConcentrationFirst(int m, int n, int q, double gasCritical);
extern void wallBoundaryConcentrationElectrodeD3Q7(int m, int n, int q);
extern void wallBoundaryConcentrationSecond(int m, int n, int q, double gasCritical);
extern void wallBoundaryConcentrationSecondD3Q7(int m, int n, int q, double gasCritical);
extern void wallBoundaryVelocityGhost(int m, int n, int q, double gasCritical);
extern void wallSurfaceReactionBoundaryConcentration(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double RateConstant,double equilibriumPotential,double K,double exchange);
extern void wallSurfaceReactionBoundaryConcentrationD3Q7(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double RateConstant,double equilibriumPotential,double K,double exchange);
extern void wallSurfaceReactionBoundaryConcentrationFirst(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double RateConstant,double equilibriumPotential,double K,double exchange);
extern void wallSurfaceReactionBoundaryConcentrationFirstD3Q7(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double RateConstant,double equilibriumPotential,double K,double exchange);
extern void wallSurfaceReactionBoundaryConcentrationSecond(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double diffuseCoefficiencySecond,double RateConstant,double equilibriumPotential,double K,double exchange);
extern void wallSurfaceReactionBoundaryConcentrationSecondD3Q7(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double diffuseCoefficiencySecond,double RateConstant,double equilibriumPotential,double K,double exchange);
extern void wallBoundaryElectrodePotential(int m, int n, int q);
extern void wallBoundaryElectrodePotentialD3Q7(int m, int n, int q);
extern void wallSurfaceReactionBoundaryElectrode(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double electrodeConductivity,double RateConstant,double constant,double equilibriumPotential,double K,double exchange);
extern void wallSurfaceReactionBoundaryElectrodeD3Q7(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double electrodeConductivity,double RateConstant,double constant,double equilibriumPotential,double K,double exchange);
extern void wallBoundaryElectrolytePotential(int m, int n, int q, double gasCritical);
extern void wallBoundaryElectrolytePotentialD3Q7(int m, int n, int q, double gasCritical);
extern void wallSurfaceReactionBoundaryElectrolyte(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double ***electrolyteRelaxationTime,double RateConstant,double constant,double equilibriumPotential,double K,double exchange);
extern void wallSurfaceReactionBoundaryElectrolyteD3Q7(int m, int n, int q, double gasCritical,double diffuseCoefficiency,double ***electrolyteRelaxationTime,double RateConstant,double constant,double equilibriumPotential,double K,double exchange);

#endif