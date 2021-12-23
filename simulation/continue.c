#include <stdio.h>
#include <stdlib.h>
#include "dataLoad.h"
#include "twoPhaseField.h"
#include "concentrationField.h"
#include "ghostFlowField.h"
#include "electrodePotential.h"
#include "electrolytePotential.h"
#include "domain.h"

FILE *InitialDomain;

FILE *InitialRho;
FILE *InitialUx;
FILE *InitialUy;
FILE *InitialUz;
FILE *InitialFx;
FILE *InitialFy;
FILE *InitialFz;
FILE *InitialPressure;
FILE *InitialFIn0;
FILE *InitialFIn1;
FILE *InitialFIn2;
FILE *InitialFIn3;
FILE *InitialFIn4;
FILE *InitialFIn5;
FILE *InitialFIn6;
FILE *InitialFIn7;
FILE *InitialFIn8;
FILE *InitialFIn9;
FILE *InitialFIn10;
FILE *InitialFIn11;
FILE *InitialFIn12;
FILE *InitialFIn13;
FILE *InitialFIn14;
FILE *InitialFIn15;
FILE *InitialFIn16;
FILE *InitialFIn17;
FILE *InitialFIn18;

FILE *InitialReactantSurface1;
FILE *InitialReactantSurface2;
FILE *InitialReactantSurface3;
FILE *InitialReactantSurface4;
FILE *InitialReactantSurface5;
FILE *InitialReactantSurface6;

FILE *InitialConcentration;
FILE *InitialCIn0;
FILE *InitialCIn1;
FILE *InitialCIn2;
FILE *InitialCIn3;
FILE *InitialCIn4;
FILE *InitialCIn5;
FILE *InitialCIn6;
FILE *InitialCIn7;
FILE *InitialCIn8;
FILE *InitialCIn9;
FILE *InitialCIn10;
FILE *InitialCIn11;
FILE *InitialCIn12;
FILE *InitialCIn13;
FILE *InitialCIn14;
FILE *InitialCIn15;
FILE *InitialCIn16;
FILE *InitialCIn17;
FILE *InitialCIn18;

FILE *InitialConcentrationFirst;
FILE *InitialCIn0First;
FILE *InitialCIn1First;
FILE *InitialCIn2First;
FILE *InitialCIn3First;
FILE *InitialCIn4First;
FILE *InitialCIn5First;
FILE *InitialCIn6First;
FILE *InitialCIn7First;
FILE *InitialCIn8First;
FILE *InitialCIn9First;
FILE *InitialCIn10First;
FILE *InitialCIn11First;
FILE *InitialCIn12First;
FILE *InitialCIn13First;
FILE *InitialCIn14First;
FILE *InitialCIn15First;
FILE *InitialCIn16First;
FILE *InitialCIn17First;
FILE *InitialCIn18First;

FILE *InitialConcentrationSecond;
FILE *InitialCIn0Second;
FILE *InitialCIn1Second;
FILE *InitialCIn2Second;
FILE *InitialCIn3Second;
FILE *InitialCIn4Second;
FILE *InitialCIn5Second;
FILE *InitialCIn6Second;
FILE *InitialCIn7Second;
FILE *InitialCIn8Second;
FILE *InitialCIn9Second;
FILE *InitialCIn10Second;
FILE *InitialCIn11Second;
FILE *InitialCIn12Second;
FILE *InitialCIn13Second;
FILE *InitialCIn14Second;
FILE *InitialCIn15Second;
FILE *InitialCIn16Second;
FILE *InitialCIn17Second;
FILE *InitialCIn18Second;

FILE *InitialRhoGhost;
FILE *InitialUxGhost;
FILE *InitialUyGhost;
FILE *InitialUzGhost;
FILE *InitialFIn0Ghost;
FILE *InitialFIn1Ghost;
FILE *InitialFIn2Ghost;
FILE *InitialFIn3Ghost;
FILE *InitialFIn4Ghost;
FILE *InitialFIn5Ghost;
FILE *InitialFIn6Ghost;
FILE *InitialFIn7Ghost;
FILE *InitialFIn8Ghost;
FILE *InitialFIn9Ghost;
FILE *InitialFIn10Ghost;
FILE *InitialFIn11Ghost;
FILE *InitialFIn12Ghost;
FILE *InitialFIn13Ghost;
FILE *InitialFIn14Ghost;
FILE *InitialFIn15Ghost;
FILE *InitialFIn16Ghost;
FILE *InitialFIn17Ghost;
FILE *InitialFIn18Ghost;

FILE *InitialElectrodePotential;
FILE *InitialElectrodeIn0;
FILE *InitialElectrodeIn1;
FILE *InitialElectrodeIn2;
FILE *InitialElectrodeIn3;
FILE *InitialElectrodeIn4;
FILE *InitialElectrodeIn5;
FILE *InitialElectrodeIn6;
FILE *InitialElectrodeIn7;
FILE *InitialElectrodeIn8;
FILE *InitialElectrodeIn9;
FILE *InitialElectrodeIn10;
FILE *InitialElectrodeIn11;
FILE *InitialElectrodeIn12;
FILE *InitialElectrodeIn13;
FILE *InitialElectrodeIn14;
FILE *InitialElectrodeIn15;
FILE *InitialElectrodeIn16;
FILE *InitialElectrodeIn17;
FILE *InitialElectrodeIn18;

FILE *InitialElectrolytePotential;
FILE *InitialElectrolyteIn0;
FILE *InitialElectrolyteIn1;
FILE *InitialElectrolyteIn2;
FILE *InitialElectrolyteIn3;
FILE *InitialElectrolyteIn4;
FILE *InitialElectrolyteIn5;
FILE *InitialElectrolyteIn6;
FILE *InitialElectrolyteIn7;
FILE *InitialElectrolyteIn8;
FILE *InitialElectrolyteIn9;
FILE *InitialElectrolyteIn10;
FILE *InitialElectrolyteIn11;
FILE *InitialElectrolyteIn12;
FILE *InitialElectrolyteIn13;
FILE *InitialElectrolyteIn14;
FILE *InitialElectrolyteIn15;
FILE *InitialElectrolyteIn16;
FILE *InitialElectrolyteIn17;
FILE *InitialElectrolyteIn18;

void fileAdressDomain(){
	InitialDomain=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\domain.txt","r");
}

void fileAddressTwoPhase(){
	InitialRho=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\rhoFlow15000.txt","r");
	InitialPressure=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\pressureFlow15000.txt","r");
	InitialUx=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\uxFlow15000.txt","r");
	InitialUy=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\uyFlow15000.txt","r");
	InitialUz=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\uzFlow15000.txt","r");
	InitialFx=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\FxFlow15000.txt","r");
	InitialFy=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\FyFlow15000.txt","r");
	InitialFz=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\FzFlow15000.txt","r");
	InitialFIn0=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn0Flow)15000.txt","r");
	InitialFIn1=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn1Flow)15000.txt","r");
	InitialFIn2=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn2Flow)15000.txt","r");
	InitialFIn3=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn3Flow)15000.txt","r");
	InitialFIn4=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn4Flow)15000.txt","r");
	InitialFIn5=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn5Flow)15000.txt","r");
	InitialFIn6=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn6Flow)15000.txt","r");
	InitialFIn7=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn7Flow)15000.txt","r");
	InitialFIn8=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn8Flow)15000.txt","r");
	InitialFIn9=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn9Flow)15000.txt","r");
	InitialFIn10=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn10Flow)15000.txt","r");
	InitialFIn11=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn11Flow)15000.txt","r");
	InitialFIn12=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn12Flow)15000.txt","r");
	InitialFIn13=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn13Flow)15000.txt","r");
	InitialFIn14=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn14Flow)15000.txt","r");
	InitialFIn15=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn15Flow)15000.txt","r");
	InitialFIn16=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn16Flow)15000.txt","r");
	InitialFIn17=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn17Flow)15000.txt","r");
	InitialFIn18=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn18Flow)15000.txt","r");
}

void fileAddressReactantSurface(){
	InitialReactantSurface1 = fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(reactant1Surface)15000.txt", "r");
	InitialReactantSurface2 = fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(reactant2Surface)15000.txt", "r");
	InitialReactantSurface3 = fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(reactant3Surface)15000.txt", "r");
	InitialReactantSurface4 = fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(reactant4Surface)15000.txt", "r");
	InitialReactantSurface5 = fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(reactant5Surface)15000.txt", "r");
	InitialReactantSurface6 = fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(reactant6Surface)15000.txt", "r");
}

void fileAddressConcentration(){
	InitialConcentration=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\ConcentrationFlow15000.txt","r");
	InitialCIn0=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn0Flow)15000.txt","r");
	InitialCIn1=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn1Flow)15000.txt","r");
	InitialCIn2=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn2Flow)15000.txt","r");
	InitialCIn3=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn3Flow)15000.txt","r");
	InitialCIn4=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn4Flow)15000.txt","r");
	InitialCIn5=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn5Flow)15000.txt","r");
	InitialCIn6=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn6Flow)15000.txt","r");
}

void fileAddressConcentrationFirst(){
	InitialConcentrationFirst=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\ConcentrationFirstFlow15000.txt","r");
	InitialCIn0First=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn0FirstFlow)15000.txt","r");
	InitialCIn1First=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn1FirstFlow)15000.txt","r");
	InitialCIn2First=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn2FirstFlow)15000.txt","r");
	InitialCIn3First=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn3FirstFlow)15000.txt","r");
	InitialCIn4First=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn4FirstFlow)15000.txt","r");
	InitialCIn5First=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn5FirstFlow)15000.txt","r");
	InitialCIn6First=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn6FirstFlow)15000.txt","r");
}

void fileAddressConcentrationSecond(){
	InitialConcentrationSecond=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\ConcentrationSecondFlow15000.txt","r");
	InitialCIn0Second=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(cIn0SecondFlow)15000.txt","r");
	InitialCIn1Second=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(cIn1SecondFlow)15000.txt","r");
	InitialCIn2Second=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(cIn2SecondFlow)15000.txt","r");
	InitialCIn3Second=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(cIn3SecondFlow)15000.txt","r");
	InitialCIn4Second=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(cIn4SecondFlow)15000.txt","r");
	InitialCIn5Second=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(cIn5SecondFlow)15000.txt","r");
	InitialCIn6Second=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(cIn6SecondFlow)15000.txt","r");
}

void fileAddressVelocityGhost(){
	InitialRhoGhost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\rhoGhostFlow30.txt","r");
	InitialUxGhost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\uxGhostFlow30.txt","r");
	InitialUyGhost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\uyGhostFlow30.txt","r");
	InitialUzGhost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\uzGhostFlow30.txt","r");
	InitialFIn0Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn0GhostFlow)30.txt","r");
	InitialFIn1Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn1GhostFlow)30.txt","r");
	InitialFIn2Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn2GhostFlow)30.txt","r");
	InitialFIn3Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn3GhostFlow)30.txt","r");
	InitialFIn4Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn4GhostFlow)30.txt","r");
	InitialFIn5Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn5GhostFlow)30.txt","r");
	InitialFIn6Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn6GhostFlow)30.txt","r");
	InitialFIn7Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn7GhostFlow)30.txt","r");
	InitialFIn8Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn8GhostFlow)30.txt","r");
	InitialFIn9Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn9GhostFlow)30.txt","r");
	InitialFIn10Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn10GhostFlow)30.txt","r");
	InitialFIn11Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn11GhostFlow)30.txt","r");
	InitialFIn12Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn12GhostFlow)30.txt","r");
	InitialFIn13Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn13GhostFlow)30.txt","r");
	InitialFIn14Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn14GhostFlow)30.txt","r");
	InitialFIn15Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn15GhostFlow)30.txt","r");
	InitialFIn16Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn16GhostFlow)30.txt","r");
	InitialFIn17Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn17GhostFlow)30.txt","r");
	InitialFIn18Ghost=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\(fIn18GhostFlow)30.txt","r");
}

void fileAddressElectrodePotential(){
	InitialElectrodePotential=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\ElectrodePotentialFlow15000.txt","r");
	InitialElectrodeIn0=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(ElectrodeIn0Flow)15000.txt","r");
	InitialElectrodeIn1=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(ElectrodeIn1Flow)15000.txt","r");
	InitialElectrodeIn2=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(ElectrodeIn2Flow)15000.txt","r");
	InitialElectrodeIn3=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(ElectrodeIn3Flow)15000.txt","r");
	InitialElectrodeIn4=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(ElectrodeIn4Flow)15000.txt","r");
	InitialElectrodeIn5=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(ElectrodeIn5Flow)15000.txt","r");
	InitialElectrodeIn6=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(ElectrodeIn6Flow)15000.txt","r");
}

void fileAddressElectrolytePotential(){
	InitialElectrolytePotential=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\ElectrolytePotentialFlow15000.txt","r");
	InitialElectrolyteIn0=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(ElectrolyteIn0Flow)15000.txt","r");
	InitialElectrolyteIn1=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(ElectrolyteIn1Flow)15000.txt","r");
	InitialElectrolyteIn2=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(ElectrolyteIn2Flow)15000.txt","r");
	InitialElectrolyteIn3=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(ElectrolyteIn3Flow)15000.txt","r");
	InitialElectrolyteIn4=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(ElectrolyteIn4Flow)15000.txt","r");
	InitialElectrolyteIn5=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(ElectrolyteIn5Flow)15000.txt","r");
	InitialElectrolyteIn6=fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(ElectrolyteIn6Flow)15000.txt","r");
}

void ContinueLoadDomain(int m, int n, int q){
	dataload(InitialDomain,domain,0,m,0,n,0,q);
}

void ContinueLoadTwoPhase(int m, int n, int q){
	dataload(InitialRho, rho, 0,m, 0,n, 0,q);
	dataload(InitialPressure, pressure, 0,m, 0,n, 0,q);
	dataload(InitialUx, ux, 0,m, 0,n, 0,q);
	dataload(InitialUy, uy, 0,m, 0,n, 0,q);
	dataload(InitialUz, uz, 0,m, 0,n, 0,q);
	dataload(InitialFx, Fx, 0,m, 0,n, 0,q);
	dataload(InitialFy, Fy, 0,m, 0,n, 0,q);
	dataload(InitialFz, Fz, 0,m, 0,n, 0,q);
	dataload(InitialFIn0, fIn0, 0,m, 0,n, 0,q);
	dataload(InitialFIn1, fIn1, 0,m, 0,n, 0,q);
	dataload(InitialFIn2, fIn2, 0,m, 0,n, 0,q);
	dataload(InitialFIn3, fIn3, 0,m, 0,n, 0,q);
	dataload(InitialFIn4, fIn4, 0,m, 0,n, 0,q);
	dataload(InitialFIn5, fIn5, 0,m, 0,n, 0,q);
	dataload(InitialFIn6, fIn6, 0,m, 0,n, 0,q);
	dataload(InitialFIn7, fIn7, 0,m, 0,n, 0,q);
	dataload(InitialFIn8, fIn8, 0,m, 0,n, 0,q);
	dataload(InitialFIn9, fIn9, 0,m, 0,n, 0,q);
	dataload(InitialFIn10, fIn10, 0,m, 0,n, 0,q);
	dataload(InitialFIn11, fIn11, 0,m, 0,n, 0,q);
	dataload(InitialFIn12, fIn12, 0,m, 0,n, 0,q);
	dataload(InitialFIn13, fIn13, 0,m, 0,n, 0,q);
	dataload(InitialFIn14, fIn14, 0,m, 0,n, 0,q);
	dataload(InitialFIn15, fIn15, 0,m, 0,n, 0,q);
	dataload(InitialFIn16, fIn16, 0,m, 0,n, 0,q);
	dataload(InitialFIn17, fIn17, 0,m, 0,n, 0,q);
	dataload(InitialFIn18, fIn18, 0,m, 0,n, 0,q);
}

void ContinueLoadReactantSurface(int m, int n, int q){
	dataload(InitialReactantSurface1, reactantSurface1, 0, m, 0, n, 0, q);
	dataload(InitialReactantSurface2, reactantSurface2, 0, m, 0, n, 0, q);
	dataload(InitialReactantSurface3, reactantSurface3, 0, m, 0, n, 0, q);
	dataload(InitialReactantSurface4, reactantSurface4, 0, m, 0, n, 0, q);
	dataload(InitialReactantSurface5, reactantSurface5, 0, m, 0, n, 0, q);
	dataload(InitialReactantSurface6, reactantSurface6, 0, m, 0, n, 0, q);
}

void ContinueLoadConcentration(int m, int n, int q){
	dataload(InitialConcentration, Concentration, 0,m, 0,n, 0,q);
	dataload(InitialCIn0, cIn0, 0,m, 0,n, 0,q);
	dataload(InitialCIn1, cIn1, 0,m, 0,n, 0,q);
	dataload(InitialCIn2, cIn2, 0,m, 0,n, 0,q);
	dataload(InitialCIn3, cIn3, 0,m, 0,n, 0,q);
	dataload(InitialCIn4, cIn4, 0,m, 0,n, 0,q);
	dataload(InitialCIn5, cIn5, 0,m, 0,n, 0,q);
	dataload(InitialCIn6, cIn6, 0,m, 0,n, 0,q);
}

void ContinueLoadConcentrationFirst(int m, int n, int q){
	dataload(InitialConcentrationFirst, ConcentrationFirst, 0,m, 0,n, 0,q);
	dataload(InitialCIn0First, cIn0First, 0,m, 0,n, 0,q);
	dataload(InitialCIn1First, cIn1First, 0,m, 0,n, 0,q);
	dataload(InitialCIn2First, cIn2First, 0,m, 0,n, 0,q);
	dataload(InitialCIn3First, cIn3First, 0,m, 0,n, 0,q);
	dataload(InitialCIn4First, cIn4First, 0,m, 0,n, 0,q);
	dataload(InitialCIn5First, cIn5First, 0,m, 0,n, 0,q);
	dataload(InitialCIn6First, cIn6First, 0,m, 0,n, 0,q);
}

void ContinueLoadConcentrationSecond(int m, int n, int q){
	dataload(InitialConcentrationSecond, ConcentrationSecond, 0,m, 0,n, 0,q);
	dataload(InitialCIn0Second, cIn0Second, 0,m, 0,n, 0,q);
	dataload(InitialCIn1Second, cIn1Second, 0,m, 0,n, 0,q);
	dataload(InitialCIn2Second, cIn2Second, 0,m, 0,n, 0,q);
	dataload(InitialCIn3Second, cIn3Second, 0,m, 0,n, 0,q);
	dataload(InitialCIn4Second, cIn4Second, 0,m, 0,n, 0,q);
	dataload(InitialCIn5Second, cIn5Second, 0,m, 0,n, 0,q);
	dataload(InitialCIn6Second, cIn6Second, 0,m, 0,n, 0,q);
}

void ContinueLoadVelocityGhost(int m, int n, int q){
	dataload(InitialRhoGhost, rhoGhost, 0,m, 0,n, 0,q);
	dataload(InitialUxGhost, uxGhost, 0,m, 0,n, 0,q);
	dataload(InitialUyGhost, uyGhost, 0,m, 0,n, 0,q);
	dataload(InitialUzGhost, uzGhost, 0,m, 0,n, 0,q);
	dataload(InitialFIn0Ghost, fIn0Ghost, 0,m, 0,n, 0,q);
	dataload(InitialFIn1Ghost, fIn1Ghost, 0,m, 0,n, 0,q);
	dataload(InitialFIn2Ghost, fIn2Ghost, 0,m, 0,n, 0,q);
	dataload(InitialFIn3Ghost, fIn3Ghost, 0,m, 0,n, 0,q);
	dataload(InitialFIn4Ghost, fIn4Ghost, 0,m, 0,n, 0,q);
	dataload(InitialFIn5Ghost, fIn5Ghost, 0,m, 0,n, 0,q);
	dataload(InitialFIn6Ghost, fIn6Ghost, 0,m, 0,n, 0,q);
	dataload(InitialFIn7Ghost, fIn7Ghost, 0,m, 0,n, 0,q);
	dataload(InitialFIn8Ghost, fIn8Ghost, 0,m, 0,n, 0,q);
	dataload(InitialFIn9Ghost, fIn9Ghost, 0,m, 0,n, 0,q);
	dataload(InitialFIn10Ghost, fIn10Ghost, 0,m, 0,n, 0,q);
	dataload(InitialFIn11Ghost, fIn11Ghost, 0,m, 0,n, 0,q);
	dataload(InitialFIn12Ghost, fIn12Ghost, 0,m, 0,n, 0,q);
	dataload(InitialFIn13Ghost, fIn13Ghost, 0,m, 0,n, 0,q);
	dataload(InitialFIn14Ghost, fIn14Ghost, 0,m, 0,n, 0,q);
	dataload(InitialFIn15Ghost, fIn15Ghost, 0,m, 0,n, 0,q);
	dataload(InitialFIn16Ghost, fIn16Ghost, 0,m, 0,n, 0,q);
	dataload(InitialFIn17Ghost, fIn17Ghost, 0,m, 0,n, 0,q);
	dataload(InitialFIn18Ghost, fIn18Ghost, 0,m, 0,n, 0,q);
}

void ContinueLoadElectrodePotential(int m, int n, int q){
	dataload(InitialElectrodePotential, electrodePotential, 0,m, 0,n, 0,q);
	dataload(InitialElectrodeIn0, electrodeIn0, 0,m, 0,n, 0,q);
	dataload(InitialElectrodeIn1, electrodeIn1, 0,m, 0,n, 0,q);
	dataload(InitialElectrodeIn2, electrodeIn2, 0,m, 0,n, 0,q);
	dataload(InitialElectrodeIn3, electrodeIn3, 0,m, 0,n, 0,q);
	dataload(InitialElectrodeIn4, electrodeIn4, 0,m, 0,n, 0,q);
	dataload(InitialElectrodeIn5, electrodeIn5, 0,m, 0,n, 0,q);
	dataload(InitialElectrodeIn6, electrodeIn6, 0,m, 0,n, 0,q);
}

void ContinueLoadElectrolytePotential(int m, int n, int q){
	dataload(InitialElectrolytePotential, electrolytePotential, 0,m, 0,n, 0,q);
	dataload(InitialElectrolyteIn0, electrolyteIn0, 0,m, 0,n, 0,q);
	dataload(InitialElectrolyteIn1, electrolyteIn1, 0,m, 0,n, 0,q);
	dataload(InitialElectrolyteIn2, electrolyteIn2, 0,m, 0,n, 0,q);
	dataload(InitialElectrolyteIn3, electrolyteIn3, 0,m, 0,n, 0,q);
	dataload(InitialElectrolyteIn4, electrolyteIn4, 0,m, 0,n, 0,q);
	dataload(InitialElectrolyteIn5, electrolyteIn5, 0,m, 0,n, 0,q);
	dataload(InitialElectrolyteIn6, electrolyteIn6, 0,m, 0,n, 0,q);
}



