#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "twoPhaseField.h"
#include "concentrationField.h"
#include "ghostFlowField.h"
#include "electrodePotential.h"
#include "electrolytePotential.h"
#include "dataStore.h"
#include "OutMemoryArrange.h"
#include "equilibrium.h"
#include "domain.h"

FILE *fpRho;
FILE *fpUx;
FILE *fpUy;
FILE *fpUz;
FILE *fpFx;
FILE *fpFy;
FILE *fpFz;
FILE *fpPressure;
FILE *fpFIn0;
FILE *fpFIn1;
FILE *fpFIn2;
FILE *fpFIn3;
FILE *fpFIn4;
FILE *fpFIn5;
FILE *fpFIn6;
FILE *fpFIn7;
FILE *fpFIn8;
FILE *fpFIn9;
FILE *fpFIn10;
FILE *fpFIn11;
FILE *fpFIn12;
FILE *fpFIn13;
FILE *fpFIn14;
FILE *fpFIn15;
FILE *fpFIn16;
FILE *fpFIn17;
FILE *fpFIn18;
FILE *fpDomain;

FILE *fpConcentration;
FILE *fpcIn0;
FILE *fpcIn1;
FILE *fpcIn2;
FILE *fpcIn3;
FILE *fpcIn4;
FILE *fpcIn5;
FILE *fpcIn6;
FILE *fpcIn7;
FILE *fpcIn8;
FILE *fpcIn9;
FILE *fpcIn10;
FILE *fpcIn11;
FILE *fpcIn12;
FILE *fpcIn13;
FILE *fpcIn14;
FILE *fpcIn15;
FILE *fpcIn16;
FILE *fpcIn17;
FILE *fpcIn18;

FILE *fpConcentrationFirst;
FILE *fpcIn0First;
FILE *fpcIn1First;
FILE *fpcIn2First;
FILE *fpcIn3First;
FILE *fpcIn4First;
FILE *fpcIn5First;
FILE *fpcIn6First;
FILE *fpcIn7First;
FILE *fpcIn8First;
FILE *fpcIn9First;
FILE *fpcIn10First;
FILE *fpcIn11First;
FILE *fpcIn12First;
FILE *fpcIn13First;
FILE *fpcIn14First;
FILE *fpcIn15First;
FILE *fpcIn16First;
FILE *fpcIn17First;
FILE *fpcIn18First;

FILE *fpConcentrationSecond;
FILE *fpcIn0Second;
FILE *fpcIn1Second;
FILE *fpcIn2Second;
FILE *fpcIn3Second;
FILE *fpcIn4Second;
FILE *fpcIn5Second;
FILE *fpcIn6Second;
FILE *fpcIn7Second;
FILE *fpcIn8Second;
FILE *fpcIn9Second;
FILE *fpcIn10Second;
FILE *fpcIn11Second;
FILE *fpcIn12Second;
FILE *fpcIn13Second;
FILE *fpcIn14Second;
FILE *fpcIn15Second;
FILE *fpcIn16Second;
FILE *fpcIn17Second;
FILE *fpcIn18Second;

FILE *fpRhoGhost;
FILE *fpUxGhost;
FILE *fpUyGhost;
FILE *fpUzGhost;
FILE *fpFIn0Ghost;
FILE *fpFIn1Ghost;
FILE *fpFIn2Ghost;
FILE *fpFIn3Ghost;
FILE *fpFIn4Ghost;
FILE *fpFIn5Ghost;
FILE *fpFIn6Ghost;
FILE *fpFIn7Ghost;
FILE *fpFIn8Ghost;
FILE *fpFIn9Ghost;
FILE *fpFIn10Ghost;
FILE *fpFIn11Ghost;
FILE *fpFIn12Ghost;
FILE *fpFIn13Ghost;
FILE *fpFIn14Ghost;
FILE *fpFIn15Ghost;
FILE *fpFIn16Ghost;
FILE *fpFIn17Ghost;
FILE *fpFIn18Ghost;

FILE *fpElectrodePotential;
FILE *fpelectrodeIn0;
FILE *fpelectrodeIn1;
FILE *fpelectrodeIn2;
FILE *fpelectrodeIn3;
FILE *fpelectrodeIn4;
FILE *fpelectrodeIn5;
FILE *fpelectrodeIn6;
FILE *fpelectrodeIn7;
FILE *fpelectrodeIn8;
FILE *fpelectrodeIn9;
FILE *fpelectrodeIn10;
FILE *fpelectrodeIn11;
FILE *fpelectrodeIn12;
FILE *fpelectrodeIn13;
FILE *fpelectrodeIn14;
FILE *fpelectrodeIn15;
FILE *fpelectrodeIn16;
FILE *fpelectrodeIn17;
FILE *fpelectrodeIn18;

FILE *fpElectrolytePotential;
FILE *fpelectrolyteIn0;
FILE *fpelectrolyteIn1;
FILE *fpelectrolyteIn2;
FILE *fpelectrolyteIn3;
FILE *fpelectrolyteIn4;
FILE *fpelectrolyteIn5;
FILE *fpelectrolyteIn6;
FILE *fpelectrolyteIn7;
FILE *fpelectrolyteIn8;
FILE *fpelectrolyteIn9;
FILE *fpelectrolyteIn10;
FILE *fpelectrolyteIn11;
FILE *fpelectrolyteIn12;
FILE *fpelectrolyteIn13;
FILE *fpelectrolyteIn14;
FILE *fpelectrolyteIn15;
FILE *fpelectrolyteIn16;
FILE *fpelectrolyteIn17;
FILE *fpelectrolyteIn18;

FILE *fpreactantSurface1;
FILE *fpreactantSurface2;
FILE *fpreactantSurface3;
FILE *fpreactantSurface4;
FILE *fpreactantSurface5;
FILE *fpreactantSurface6;

char numstep[25];
char filenameRho[200];
char filenameUx[200];
char filenameUy[200];
char filenameUz[200];
char filenameFx[200];
char filenameFy[200];
char filenameFz[200];
char filenamePressure[200];
char filenamefIn0[200];
char filenamefIn1[200];
char filenamefIn2[200];
char filenamefIn3[200];
char filenamefIn4[200];
char filenamefIn5[200];
char filenamefIn6[200];
char filenamefIn7[200];
char filenamefIn8[200];
char filenamefIn9[200];
char filenamefIn10[200];
char filenamefIn11[200];
char filenamefIn12[200];
char filenamefIn13[200];
char filenamefIn14[200];
char filenamefIn15[200];
char filenamefIn16[200];
char filenamefIn17[200];
char filenamefIn18[200];
char filenameDomain[200];
char p[]=".txt";

char filenameConcentration[200];
char filenamecIn0[200];
char filenamecIn1[200];
char filenamecIn2[200];
char filenamecIn3[200];
char filenamecIn4[200];
char filenamecIn5[200];
char filenamecIn6[200];
char filenamecIn7[200];
char filenamecIn8[200];
char filenamecIn9[200];
char filenamecIn10[200];
char filenamecIn11[200];
char filenamecIn12[200];
char filenamecIn13[200];
char filenamecIn14[200];
char filenamecIn15[200];
char filenamecIn16[200];
char filenamecIn17[200];
char filenamecIn18[200];

char filenameConcentrationFirst[200];
char filenamecIn0First[200];
char filenamecIn1First[200];
char filenamecIn2First[200];
char filenamecIn3First[200];
char filenamecIn4First[200];
char filenamecIn5First[200];
char filenamecIn6First[200];
char filenamecIn7First[200];
char filenamecIn8First[200];
char filenamecIn9First[200];
char filenamecIn10First[200];
char filenamecIn11First[200];
char filenamecIn12First[200];
char filenamecIn13First[200];
char filenamecIn14First[200];
char filenamecIn15First[200];
char filenamecIn16First[200];
char filenamecIn17First[200];
char filenamecIn18First[200];

char filenameConcentrationSecond[200];
char filenamecIn0Second[200];
char filenamecIn1Second[200];
char filenamecIn2Second[200];
char filenamecIn3Second[200];
char filenamecIn4Second[200];
char filenamecIn5Second[200];
char filenamecIn6Second[200];
char filenamecIn7Second[200];
char filenamecIn8Second[200];
char filenamecIn9Second[200];
char filenamecIn10Second[200];
char filenamecIn11Second[200];
char filenamecIn12Second[200];
char filenamecIn13Second[200];
char filenamecIn14Second[200];
char filenamecIn15Second[200];
char filenamecIn16Second[200];
char filenamecIn17Second[200];
char filenamecIn18Second[200];

char filenameRhoGhost[200];
char filenameUxGhost[200];
char filenameUyGhost[200];
char filenameUzGhost[200];
char filenamefIn0Ghost[200];
char filenamefIn1Ghost[200];
char filenamefIn2Ghost[200];
char filenamefIn3Ghost[200];
char filenamefIn4Ghost[200];
char filenamefIn5Ghost[200];
char filenamefIn6Ghost[200];
char filenamefIn7Ghost[200];
char filenamefIn8Ghost[200];
char filenamefIn9Ghost[200];
char filenamefIn10Ghost[200];
char filenamefIn11Ghost[200];
char filenamefIn12Ghost[200];
char filenamefIn13Ghost[200];
char filenamefIn14Ghost[200];
char filenamefIn15Ghost[200];
char filenamefIn16Ghost[200];
char filenamefIn17Ghost[200];
char filenamefIn18Ghost[200];

char filenameElectrodePotential[200];
char filenameelectrodeIn0[200];
char filenameelectrodeIn1[200];
char filenameelectrodeIn2[200];
char filenameelectrodeIn3[200];
char filenameelectrodeIn4[200];
char filenameelectrodeIn5[200];
char filenameelectrodeIn6[200];
char filenameelectrodeIn7[200];
char filenameelectrodeIn8[200];
char filenameelectrodeIn9[200];
char filenameelectrodeIn10[200];
char filenameelectrodeIn11[200];
char filenameelectrodeIn12[200];
char filenameelectrodeIn13[200];
char filenameelectrodeIn14[200];
char filenameelectrodeIn15[200];
char filenameelectrodeIn16[200];
char filenameelectrodeIn17[200];
char filenameelectrodeIn18[200];

char filenameElectrolytePotential[200];
char filenameelectrolyteIn0[200];
char filenameelectrolyteIn1[200];
char filenameelectrolyteIn2[200];
char filenameelectrolyteIn3[200];
char filenameelectrolyteIn4[200];
char filenameelectrolyteIn5[200];
char filenameelectrolyteIn6[200];
char filenameelectrolyteIn7[200];
char filenameelectrolyteIn8[200];
char filenameelectrolyteIn9[200];
char filenameelectrolyteIn10[200];
char filenameelectrolyteIn11[200];
char filenameelectrolyteIn12[200];
char filenameelectrolyteIn13[200];
char filenameelectrolyteIn14[200];
char filenameelectrolyteIn15[200];
char filenameelectrolyteIn16[200];
char filenameelectrolyteIn17[200];
char filenameelectrolyteIn18[200];

char filenamereactantSurface1[200];
char filenamereactantSurface2[200];
char filenamereactantSurface3[200];
char filenamereactantSurface4[200];
char filenamereactantSurface5[200];
char filenamereactantSurface6[200];
void resultsWriteDomain(int m, int n, int q){
		strcpy(filenameDomain,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\electrodesimulation.txt");
		strcat(filenameDomain,p);
		fpDomain=fopen(filenameDomain,"w");
		fileprintout(fpDomain,domain,m,n,q);
}

void resultsWriteReactantSurface(int m, int n, int q, int step, int INTERVAL){
	if (step%INTERVAL == 0){
		_itoa(step, numstep, 10);
		strcpy(filenamereactantSurface1, "C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(reactant1Surface)");
		strcat(filenamereactantSurface1, numstep);
		strcat(filenamereactantSurface1, p);
		fpreactantSurface1 = fopen(filenamereactantSurface1, "w");
		fileprintout(fpreactantSurface1, reaction1, m, n, q);

		strcpy(filenamereactantSurface2, "C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(reactant2Surface)");
		strcat(filenamereactantSurface2, numstep);
		strcat(filenamereactantSurface2, p);
		fpreactantSurface2 = fopen(filenamereactantSurface2, "w");
		fileprintout(fpreactantSurface2, reaction2, m, n, q);

		strcpy(filenamereactantSurface3, "C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(reactant3Surface)");
		strcat(filenamereactantSurface3, numstep);
		strcat(filenamereactantSurface3, p);
		fpreactantSurface3 = fopen(filenamereactantSurface3, "w");
		fileprintout(fpreactantSurface3, reaction3, m, n, q);

		strcpy(filenamereactantSurface4, "C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(reactant4Surface)");
		strcat(filenamereactantSurface4, numstep);
		strcat(filenamereactantSurface4, p);
		fpreactantSurface4 = fopen(filenamereactantSurface4, "w");
		fileprintout(fpreactantSurface4, reaction4, m, n, q);

		strcpy(filenamereactantSurface5, "C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(reactant5Surface)");
		strcat(filenamereactantSurface5, numstep);
		strcat(filenamereactantSurface5, p);
		fpreactantSurface5 = fopen(filenamereactantSurface5, "w");
		fileprintout(fpreactantSurface5, reaction5, m, n, q);

		strcpy(filenamereactantSurface6, "C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(reactant6Surface)");
		strcat(filenamereactantSurface6, numstep);
		strcat(filenamereactantSurface6, p);
		fpreactantSurface6 = fopen(filenamereactantSurface6, "w");
		fileprintout(fpreactantSurface6, reaction6, m, n, q);
	}
}

void resultsWriteTwoPhase(int m,int n,int q,int step,int INTERVAL){
	if (step%INTERVAL==0){
		_itoa(step,numstep,10);
		strcpy(filenameRho,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\rhoFlow");
		strcat(filenameRho,numstep);
		strcat(filenameRho,p);
		fpRho=fopen(filenameRho,"w");
		fileprintout(fpRho,rho,m,n,q);

		strcpy(filenameUx,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\uxFlow");
		strcat(filenameUx,numstep);
		strcat(filenameUx,p);
		fpUx=fopen(filenameUx,"w");
		fileprintout(fpUx,ux,m,n,q);

		strcpy(filenameUy,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\uyFlow");
		strcat(filenameUy,numstep);
		strcat(filenameUy,p);
		fpUy=fopen(filenameUy,"w");
		fileprintout(fpUy,uy,m,n,q);

		strcpy(filenameUz,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\uzFlow");
		strcat(filenameUz,numstep);
		strcat(filenameUz,p);
		fpUz=fopen(filenameUz,"w");
		fileprintout(fpUz,uz,m,n,q);

		strcpy(filenameFx,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\FxFlow");
		strcat(filenameFx,numstep);
		strcat(filenameFx,p);
		fpFx=fopen(filenameFx,"w");
		fileprintout(fpFx,Fx,m,n,q);

		strcpy(filenameFy,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\FyFlow");
		strcat(filenameFy,numstep);
		strcat(filenameFy,p);
		fpFy=fopen(filenameFy,"w");
		fileprintout(fpFy,Fy,m,n,q);

		strcpy(filenameFz,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\FzFlow");
		strcat(filenameFz,numstep);
		strcat(filenameFz,p);
		fpFz=fopen(filenameFz,"w");
		fileprintout(fpFz,Fz,m,n,q);

		strcpy(filenamePressure,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\pressureFlow");
		strcat(filenamePressure,numstep);
		strcat(filenamePressure,p);
		fpPressure=fopen(filenamePressure,"w");
		fileprintout(fpPressure,pressure,m,n,q);

		strcpy(filenamefIn0,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn0Flow)");
		strcat(filenamefIn0,numstep);
		strcat(filenamefIn0,p);
		fpFIn0=fopen(filenamefIn0,"w");
		fileprintout(fpFIn0,fIn0,m,n,q);

		strcpy(filenamefIn1,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn1Flow)");
		strcat(filenamefIn1,numstep);
		strcat(filenamefIn1,p);
		fpFIn1=fopen(filenamefIn1,"w");
		fileprintout(fpFIn1,fIn1,m,n,q);

		strcpy(filenamefIn2,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn2Flow)");
		strcat(filenamefIn2,numstep);
		strcat(filenamefIn2,p);
		fpFIn2=fopen(filenamefIn2,"w");
		fileprintout(fpFIn2,fIn2,m,n,q);

		strcpy(filenamefIn3,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn3Flow)");
		strcat(filenamefIn3,numstep);
		strcat(filenamefIn3,p);
		fpFIn3=fopen(filenamefIn3,"w");
		fileprintout(fpFIn3,fIn3,m,n,q);

		strcpy(filenamefIn4,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn4Flow)");
		strcat(filenamefIn4,numstep);
		strcat(filenamefIn4,p);
		fpFIn4=fopen(filenamefIn4,"w");
		fileprintout(fpFIn4,fIn4,m,n,q);

		strcpy(filenamefIn5,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn5Flow)");
		strcat(filenamefIn5,numstep);
		strcat(filenamefIn5,p);
		fpFIn5=fopen(filenamefIn5,"w");
		fileprintout(fpFIn5,fIn5,m,n,q);

		strcpy(filenamefIn6,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn6Flow)");
		strcat(filenamefIn6,numstep);
		strcat(filenamefIn6,p);
		fpFIn6=fopen(filenamefIn6,"w");
		fileprintout(fpFIn6,fIn6,m,n,q);

		strcpy(filenamefIn7,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn7Flow)");
		strcat(filenamefIn7,numstep);
		strcat(filenamefIn7,p);
		fpFIn7=fopen(filenamefIn7,"w");
		fileprintout(fpFIn7,fIn7,m,n,q);

		strcpy(filenamefIn8,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn8Flow)");
		strcat(filenamefIn8,numstep);
		strcat(filenamefIn8,p);
		fpFIn8=fopen(filenamefIn8,"w");
		fileprintout(fpFIn8,fIn8,m,n,q);

		strcpy(filenamefIn9,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn9Flow)");
		strcat(filenamefIn9,numstep);
		strcat(filenamefIn9,p);
		fpFIn9=fopen(filenamefIn9,"w");
		fileprintout(fpFIn9,fIn9,m,n,q);

		strcpy(filenamefIn10,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn10Flow)");
		strcat(filenamefIn10,numstep);
		strcat(filenamefIn10,p);
		fpFIn10=fopen(filenamefIn10,"w");
		fileprintout(fpFIn10,fIn10,m,n,q);

		strcpy(filenamefIn11,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn11Flow)");
		strcat(filenamefIn11,numstep);
		strcat(filenamefIn11,p);
		fpFIn11=fopen(filenamefIn11,"w");
		fileprintout(fpFIn11,fIn11,m,n,q);

		strcpy(filenamefIn12,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn12Flow)");
		strcat(filenamefIn12,numstep);
		strcat(filenamefIn12,p);
		fpFIn12=fopen(filenamefIn12,"w");
		fileprintout(fpFIn12,fIn12,m,n,q);

		strcpy(filenamefIn13,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn13Flow)");
		strcat(filenamefIn13,numstep);
		strcat(filenamefIn13,p);
		fpFIn13=fopen(filenamefIn13,"w");
		fileprintout(fpFIn13,fIn13,m,n,q);

		strcpy(filenamefIn14,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn14Flow)");
		strcat(filenamefIn14,numstep);
		strcat(filenamefIn14,p);
		fpFIn14=fopen(filenamefIn14,"w");
		fileprintout(fpFIn14,fIn14,m,n,q);

		strcpy(filenamefIn15,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn15Flow)");
		strcat(filenamefIn15,numstep);
		strcat(filenamefIn15,p);
		fpFIn15=fopen(filenamefIn15,"w");
		fileprintout(fpFIn15,fIn15,m,n,q);

		strcpy(filenamefIn16,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn16Flow)");
		strcat(filenamefIn16,numstep);
		strcat(filenamefIn16,p);
		fpFIn16=fopen(filenamefIn16,"w");
		fileprintout(fpFIn16,fIn16,m,n,q);

		strcpy(filenamefIn17,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn17Flow)");
		strcat(filenamefIn17,numstep);
		strcat(filenamefIn17,p);
		fpFIn17=fopen(filenamefIn17,"w");
		fileprintout(fpFIn17,fIn17,m,n,q);

		strcpy(filenamefIn18,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn18Flow)");
		strcat(filenamefIn18,numstep);
		strcat(filenamefIn18,p);
		fpFIn18=fopen(filenamefIn18,"w");
		fileprintout(fpFIn18,fIn18,m,n,q);
	}
}

void resultsWriteConcentration(int m,int n,int q,int step,int INTERVAL){
	if (step%INTERVAL==0){
		_itoa(step,numstep,10);
		strcpy(filenameConcentration,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\concentrationFlow");
		strcat(filenameConcentration,numstep);
		strcat(filenameConcentration,p);
		fpConcentration=fopen(filenameConcentration,"w");
		fileprintout(fpConcentration, Concentration, m, n, q);

		strcpy(filenamecIn0,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn0Flow)");
		strcat(filenamecIn0,numstep);
		strcat(filenamecIn0,p);
		fpcIn0=fopen(filenamecIn0,"w");
		fileprintout(fpcIn0,cIn0,m,n,q);

		strcpy(filenamecIn1,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn1Flow)");
		strcat(filenamecIn1,numstep);
		strcat(filenamecIn1,p);
		fpcIn1=fopen(filenamecIn1,"w");
		fileprintout(fpcIn1,cIn1,m,n,q);

		strcpy(filenamecIn2,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn2Flow)");
		strcat(filenamecIn2,numstep);
		strcat(filenamecIn2,p);
		fpcIn2=fopen(filenamecIn2,"w");
		fileprintout(fpcIn2,cIn2,m,n,q);

		strcpy(filenamecIn3,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn3Flow)");
		strcat(filenamecIn3,numstep);
		strcat(filenamecIn3,p);
		fpcIn3=fopen(filenamecIn3,"w");
		fileprintout(fpcIn3,cIn3,m,n,q);

		strcpy(filenamecIn4,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn4Flow)");
		strcat(filenamecIn4,numstep);
		strcat(filenamecIn4,p);
		fpcIn4=fopen(filenamecIn4,"w");
		fileprintout(fpcIn4,cIn4,m,n,q);

		strcpy(filenamecIn5,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn5Flow)");
		strcat(filenamecIn5,numstep);
		strcat(filenamecIn5,p);
		fpcIn5=fopen(filenamecIn5,"w");
		fileprintout(fpcIn5,cIn5,m,n,q);

		strcpy(filenamecIn6,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn6Flow)");
		strcat(filenamecIn6,numstep);
		strcat(filenamecIn6,p);
		fpcIn6=fopen(filenamecIn6,"w");
		fileprintout(fpcIn6,cIn6,m,n,q);
	}
}

void resultsWriteConcentrationFirst(int m,int n,int q,int step,int INTERVAL){
	if (step%INTERVAL==0){
		_itoa(step,numstep,10);
		strcpy(filenameConcentrationFirst,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\concentrationFirstFlow");
		strcat(filenameConcentrationFirst,numstep);
		strcat(filenameConcentrationFirst,p);
		fpConcentrationFirst=fopen(filenameConcentrationFirst,"w");
		fileprintout(fpConcentrationFirst,ConcentrationFirst,m,n,q);

		strcpy(filenamecIn0First,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn0FirstFlow)");
		strcat(filenamecIn0First,numstep);
		strcat(filenamecIn0First,p);
		fpcIn0First=fopen(filenamecIn0First,"w");
		fileprintout(fpcIn0First,cIn0First,m,n,q);

		strcpy(filenamecIn1First,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn1FirstFlow)");
		strcat(filenamecIn1First,numstep);
		strcat(filenamecIn1First,p);
		fpcIn1First=fopen(filenamecIn1First,"w");
		fileprintout(fpcIn1First,cIn1First,m,n,q);

		strcpy(filenamecIn2First,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn2FirstFlow)");
		strcat(filenamecIn2First,numstep);
		strcat(filenamecIn2First,p);
		fpcIn2First=fopen(filenamecIn2First,"w");
		fileprintout(fpcIn2First,cIn2First,m,n,q);

		strcpy(filenamecIn3First,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn3FirstFlow)");
		strcat(filenamecIn3First,numstep);
		strcat(filenamecIn3First,p);
		fpcIn3First=fopen(filenamecIn3First,"w");
		fileprintout(fpcIn3First,cIn3First,m,n,q);

		strcpy(filenamecIn4First,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn4FirstFlow)");
		strcat(filenamecIn4First,numstep);
		strcat(filenamecIn4First,p);
		fpcIn4First=fopen(filenamecIn4First,"w");
		fileprintout(fpcIn4First,cIn4First,m,n,q);

		strcpy(filenamecIn5First,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn5FirstFlow)");
		strcat(filenamecIn5First,numstep);
		strcat(filenamecIn5First,p);
		fpcIn5First=fopen(filenamecIn5First,"w");
		fileprintout(fpcIn5First,cIn5First,m,n,q);

		strcpy(filenamecIn6First,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn6FirstFlow)");
		strcat(filenamecIn6First,numstep);
		strcat(filenamecIn6First,p);
		fpcIn6First=fopen(filenamecIn6First,"w");
		fileprintout(fpcIn6First,cIn6First,m,n,q);
	}
}

void resultsWriteConcentrationSecond(int m,int n,int q,int step,int INTERVAL){
	if (step%INTERVAL==0){
		_itoa(step,numstep,10);
		strcpy(filenameConcentrationSecond,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\concentrationSecondFlow");
		strcat(filenameConcentrationSecond,numstep);
		strcat(filenameConcentrationSecond,p);
		fpConcentrationSecond=fopen(filenameConcentrationSecond,"w");
		fileprintout(fpConcentrationSecond, ConcentrationSecond, m, n, q);

		strcpy(filenamecIn0Second,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn0SecondFlow)");
		strcat(filenamecIn0Second,numstep);
		strcat(filenamecIn0Second,p);
		fpcIn0Second=fopen(filenamecIn0Second,"w");
		fileprintout(fpcIn0Second,cIn0Second,m,n,q);

		strcpy(filenamecIn1Second,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn1SecondFlow)");
		strcat(filenamecIn1Second,numstep);
		strcat(filenamecIn1Second,p);
		fpcIn1Second=fopen(filenamecIn1Second,"w");
		fileprintout(fpcIn1Second,cIn1Second,m,n,q);

		strcpy(filenamecIn2Second,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn2SecondFlow)");
		strcat(filenamecIn2Second,numstep);
		strcat(filenamecIn2Second,p);
		fpcIn2Second=fopen(filenamecIn2Second,"w");
		fileprintout(fpcIn2Second,cIn2Second,m,n,q);

		strcpy(filenamecIn3Second,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn3SecondFlow)");
		strcat(filenamecIn3Second,numstep);
		strcat(filenamecIn3Second,p);
		fpcIn3Second=fopen(filenamecIn3Second,"w");
		fileprintout(fpcIn3Second,cIn3Second,m,n,q);

		strcpy(filenamecIn4Second,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn4SecondFlow)");
		strcat(filenamecIn4Second,numstep);
		strcat(filenamecIn4Second,p);
		fpcIn4Second=fopen(filenamecIn4Second,"w");
		fileprintout(fpcIn4Second,cIn4Second,m,n,q);

		strcpy(filenamecIn5Second,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn5SecondFlow)");
		strcat(filenamecIn5Second,numstep);
		strcat(filenamecIn5Second,p);
		fpcIn5Second=fopen(filenamecIn5Second,"w");
		fileprintout(fpcIn5Second,cIn5Second,m,n,q);

		strcpy(filenamecIn6Second,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(cIn6SecondFlow)");
		strcat(filenamecIn6Second,numstep);
		strcat(filenamecIn6Second,p);
		fpcIn6Second=fopen(filenamecIn6Second,"w");
		fileprintout(fpcIn6Second,cIn6Second,m,n,q);
	}
}

void resultsWriteVelocityGhost(int m,int n,int q,int step,int INTERVAL){
	if (step%INTERVAL==0){
		_itoa(step,numstep,10);
		strcpy(filenameRhoGhost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\rhoGhostFlow");
		strcat(filenameRhoGhost,numstep);
		strcat(filenameRhoGhost,p);
		fpRhoGhost=fopen(filenameRhoGhost,"w");
		fileprintout(fpRhoGhost,rhoGhost,m,n,q);

		strcpy(filenameUxGhost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\uxGhostFlow");
		strcat(filenameUxGhost,numstep);
		strcat(filenameUxGhost,p);
		fpUxGhost=fopen(filenameUxGhost,"w");
		fileprintout(fpUxGhost,uxGhost,m,n,q);

		strcpy(filenameUyGhost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\uyGhostFlow");
		strcat(filenameUyGhost,numstep);
		strcat(filenameUyGhost,p);
		fpUyGhost=fopen(filenameUyGhost,"w");
		fileprintout(fpUyGhost,uyGhost,m,n,q);

		strcpy(filenameUzGhost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\uzGhostFlow");
		strcat(filenameUzGhost,numstep);
		strcat(filenameUzGhost,p);
		fpUzGhost=fopen(filenameUzGhost,"w");
		fileprintout(fpUzGhost,uzGhost,m,n,q);

		strcpy(filenamefIn0Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn0GhostFlow)");
		strcat(filenamefIn0Ghost,numstep);
		strcat(filenamefIn0Ghost,p);
		fpFIn0Ghost=fopen(filenamefIn0Ghost,"w");
		fileprintout(fpFIn0Ghost,fIn0Ghost,m,n,q);

		strcpy(filenamefIn1Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn1GhostFlow)");
		strcat(filenamefIn1Ghost,numstep);
		strcat(filenamefIn1Ghost,p);
		fpFIn1Ghost=fopen(filenamefIn1Ghost,"w");
		fileprintout(fpFIn1Ghost,fIn1Ghost,m,n,q);

		strcpy(filenamefIn2Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn2GhostFlow)");
		strcat(filenamefIn2Ghost,numstep);
		strcat(filenamefIn2Ghost,p);
		fpFIn2Ghost=fopen(filenamefIn2Ghost,"w");
		fileprintout(fpFIn2Ghost,fIn2Ghost,m,n,q);

		strcpy(filenamefIn3Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn3GhostFlow)");
		strcat(filenamefIn3Ghost,numstep);
		strcat(filenamefIn3Ghost,p);
		fpFIn3Ghost=fopen(filenamefIn3Ghost,"w");
		fileprintout(fpFIn3Ghost,fIn3Ghost,m,n,q);

		strcpy(filenamefIn4Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn4GhostFlow)");
		strcat(filenamefIn4Ghost,numstep);
		strcat(filenamefIn4Ghost,p);
		fpFIn4Ghost=fopen(filenamefIn4Ghost,"w");
		fileprintout(fpFIn4Ghost,fIn4Ghost,m,n,q);

		strcpy(filenamefIn5Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn5GhostFlow)");
		strcat(filenamefIn5Ghost,numstep);
		strcat(filenamefIn5Ghost,p);
		fpFIn5Ghost=fopen(filenamefIn5Ghost,"w");
		fileprintout(fpFIn5Ghost,fIn5Ghost,m,n,q);

		strcpy(filenamefIn6Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn6GhostFlow)");
		strcat(filenamefIn6Ghost,numstep);
		strcat(filenamefIn6Ghost,p);
		fpFIn6Ghost=fopen(filenamefIn6Ghost,"w");
		fileprintout(fpFIn6Ghost,fIn6Ghost,m,n,q);

		strcpy(filenamefIn7Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn7GhostFlow)");
		strcat(filenamefIn7Ghost,numstep);
		strcat(filenamefIn7Ghost,p);
		fpFIn7Ghost=fopen(filenamefIn7Ghost,"w");
		fileprintout(fpFIn7Ghost,fIn7Ghost,m,n,q);

		strcpy(filenamefIn8Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn8GhostFlow)");
		strcat(filenamefIn8Ghost,numstep);
		strcat(filenamefIn8Ghost,p);
		fpFIn8Ghost=fopen(filenamefIn8Ghost,"w");
		fileprintout(fpFIn8Ghost,fIn8Ghost,m,n,q);

		strcpy(filenamefIn9Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn9GhostFlow)");
		strcat(filenamefIn9Ghost,numstep);
		strcat(filenamefIn9Ghost,p);
		fpFIn9Ghost=fopen(filenamefIn9Ghost,"w");
		fileprintout(fpFIn9Ghost,fIn9Ghost,m,n,q);

		strcpy(filenamefIn10Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn10GhostFlow)");
		strcat(filenamefIn10Ghost,numstep);
		strcat(filenamefIn10Ghost,p);
		fpFIn10Ghost=fopen(filenamefIn10Ghost,"w");
		fileprintout(fpFIn10Ghost,fIn10Ghost,m,n,q);

		strcpy(filenamefIn11Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn11GhostFlow)");
		strcat(filenamefIn11Ghost,numstep);
		strcat(filenamefIn11Ghost,p);
		fpFIn11Ghost=fopen(filenamefIn11Ghost,"w");
		fileprintout(fpFIn11Ghost,fIn11Ghost,m,n,q);

		strcpy(filenamefIn12Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn12GhostFlow)");
		strcat(filenamefIn12Ghost,numstep);
		strcat(filenamefIn12Ghost,p);
		fpFIn12Ghost=fopen(filenamefIn12Ghost,"w");
		fileprintout(fpFIn12Ghost,fIn12Ghost,m,n,q);

		strcpy(filenamefIn13Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn13GhostFlow)");
		strcat(filenamefIn13Ghost,numstep);
		strcat(filenamefIn13Ghost,p);
		fpFIn13Ghost=fopen(filenamefIn13Ghost,"w");
		fileprintout(fpFIn13Ghost,fIn13Ghost,m,n,q);

		strcpy(filenamefIn14Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn14GhostFlow)");
		strcat(filenamefIn14Ghost,numstep);
		strcat(filenamefIn14Ghost,p);
		fpFIn14Ghost=fopen(filenamefIn14Ghost,"w");
		fileprintout(fpFIn14Ghost,fIn14Ghost,m,n,q);

		strcpy(filenamefIn15Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn15GhostFlow)");
		strcat(filenamefIn15Ghost,numstep);
		strcat(filenamefIn15Ghost,p);
		fpFIn15Ghost=fopen(filenamefIn15Ghost,"w");
		fileprintout(fpFIn15Ghost,fIn15Ghost,m,n,q);

		strcpy(filenamefIn16Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn16GhostFlow)");
		strcat(filenamefIn16Ghost,numstep);
		strcat(filenamefIn16Ghost,p);
		fpFIn16Ghost=fopen(filenamefIn16Ghost,"w");
		fileprintout(fpFIn16Ghost,fIn16Ghost,m,n,q);

		strcpy(filenamefIn17Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn17GhostFlow)");
		strcat(filenamefIn17Ghost,numstep);
		strcat(filenamefIn17Ghost,p);
		fpFIn17Ghost=fopen(filenamefIn17Ghost,"w");
		fileprintout(fpFIn17Ghost,fIn17Ghost,m,n,q);

		strcpy(filenamefIn18Ghost,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(fIn18GhostFlow)");
		strcat(filenamefIn18Ghost,numstep);
		strcat(filenamefIn18Ghost,p);
		fpFIn18Ghost=fopen(filenamefIn18Ghost,"w");
		fileprintout(fpFIn18Ghost,fIn18Ghost,m,n,q);
	}
}

void resultsWriteElectrodePotential(int m,int n,int q,int step,int INTERVAL){
	if (step%INTERVAL==0){
		_itoa(step,numstep,10);
		strcpy(filenameElectrodePotential,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\ElectrodePotentialFlow");
		strcat(filenameElectrodePotential,numstep);
		strcat(filenameElectrodePotential,p);
		fpElectrodePotential=fopen(filenameElectrodePotential,"w");
		fileprintout(fpElectrodePotential,electrodePotential,m,n,q);

		strcpy(filenameelectrodeIn0,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(electrodeIn0Flow)");
		strcat(filenameelectrodeIn0,numstep);
		strcat(filenameelectrodeIn0,p);
		fpelectrodeIn0=fopen(filenameelectrodeIn0,"w");
		fileprintout(fpelectrodeIn0,electrodeIn0,m,n,q);

		strcpy(filenameelectrodeIn1,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(electrodeIn1Flow)");
		strcat(filenameelectrodeIn1,numstep);
		strcat(filenameelectrodeIn1,p);
		fpelectrodeIn1=fopen(filenameelectrodeIn1,"w");
		fileprintout(fpelectrodeIn1,electrodeIn1,m,n,q);

		strcpy(filenameelectrodeIn2,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(electrodeIn2Flow)");
		strcat(filenameelectrodeIn2,numstep);
		strcat(filenameelectrodeIn2,p);
		fpelectrodeIn2=fopen(filenameelectrodeIn2,"w");
		fileprintout(fpelectrodeIn2,electrodeIn2,m,n,q);

		strcpy(filenameelectrodeIn3,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(electrodeIn3Flow)");
		strcat(filenameelectrodeIn3,numstep);
		strcat(filenameelectrodeIn3,p);
		fpelectrodeIn3=fopen(filenameelectrodeIn3,"w");
		fileprintout(fpelectrodeIn3,electrodeIn3,m,n,q);

		strcpy(filenameelectrodeIn4,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(electrodeIn4Flow)");
		strcat(filenameelectrodeIn4,numstep);
		strcat(filenameelectrodeIn4,p);
		fpelectrodeIn4=fopen(filenameelectrodeIn4,"w");
		fileprintout(fpelectrodeIn4,electrodeIn4,m,n,q);

		strcpy(filenameelectrodeIn5,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(electrodeIn5Flow)");
		strcat(filenameelectrodeIn5,numstep);
		strcat(filenameelectrodeIn5,p);
		fpelectrodeIn5=fopen(filenameelectrodeIn5,"w");
		fileprintout(fpelectrodeIn5,electrodeIn5,m,n,q);

		strcpy(filenameelectrodeIn6,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(electrodeIn6Flow)");
		strcat(filenameelectrodeIn6,numstep);
		strcat(filenameelectrodeIn6,p);
		fpelectrodeIn6=fopen(filenameelectrodeIn6,"w");
		fileprintout(fpelectrodeIn6,electrodeIn6,m,n,q);
	}
}

void resultsWriteElectrolytePotential(int m,int n,int q,int step,int INTERVAL){
	if (step%INTERVAL==0){
		_itoa(step,numstep,10);
		strcpy(filenameElectrolytePotential,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\ElectrolytePotentialFlow");
		strcat(filenameElectrolytePotential,numstep);
		strcat(filenameElectrolytePotential,p);
		fpElectrolytePotential=fopen(filenameElectrolytePotential,"w");
		fileprintout(fpElectrolytePotential,electrolytePotential,m,n,q);

		strcpy(filenameelectrolyteIn0,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(electrolyteIn0Flow)");
		strcat(filenameelectrolyteIn0,numstep);
		strcat(filenameelectrolyteIn0,p);
		fpelectrolyteIn0=fopen(filenameelectrolyteIn0,"w");
		fileprintout(fpelectrolyteIn0,electrolyteIn0,m,n,q);

		strcpy(filenameelectrolyteIn1,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(electrolyteIn1Flow)");
		strcat(filenameelectrolyteIn1,numstep);
		strcat(filenameelectrolyteIn1,p);
		fpelectrolyteIn1=fopen(filenameelectrolyteIn1,"w");
		fileprintout(fpelectrolyteIn1,electrolyteIn1,m,n,q);

		strcpy(filenameelectrolyteIn2,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(electrolyteIn2Flow)");
		strcat(filenameelectrolyteIn2,numstep);
		strcat(filenameelectrolyteIn2,p);
		fpelectrolyteIn2=fopen(filenameelectrolyteIn2,"w");
		fileprintout(fpelectrolyteIn2,electrolyteIn2,m,n,q);

		strcpy(filenameelectrolyteIn3,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(electrolyteIn3Flow)");
		strcat(filenameelectrolyteIn3,numstep);
		strcat(filenameelectrolyteIn3,p);
		fpelectrolyteIn3=fopen(filenameelectrolyteIn3,"w");
		fileprintout(fpelectrolyteIn3,electrolyteIn3,m,n,q);

		strcpy(filenameelectrolyteIn4,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(electrolyteIn4Flow)");
		strcat(filenameelectrolyteIn4,numstep);
		strcat(filenameelectrolyteIn4,p);
		fpelectrolyteIn4=fopen(filenameelectrolyteIn4,"w");
		fileprintout(fpelectrolyteIn4,electrolyteIn4,m,n,q);

		strcpy(filenameelectrolyteIn5,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(electrolyteIn5Flow)");
		strcat(filenameelectrolyteIn5,numstep);
		strcat(filenameelectrolyteIn5,p);
		fpelectrolyteIn5=fopen(filenameelectrolyteIn5,"w");
		fileprintout(fpelectrolyteIn5,electrolyteIn5,m,n,q);

		strcpy(filenameelectrolyteIn6,"C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\simulation\\results2C\\(electrolyteIn6Flow)");
		strcat(filenameelectrolyteIn6,numstep);
		strcat(filenameelectrolyteIn6,p);
		fpelectrolyteIn6=fopen(filenameelectrolyteIn6,"w");
		fileprintout(fpelectrolyteIn6,electrolyteIn6,m,n,q);
	}
}