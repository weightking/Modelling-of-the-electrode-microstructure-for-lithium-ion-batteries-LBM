#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <omp.h>

#include "matrixMove.h"
#include "MemoryArrange.h"
#include "memoryFree.h"
#include "dataStore.h"
#include "dataLoad.h"

#include "domain.h"
#include "wallContactAngle.h"
#include "OutMemoryArrange.h"
#include "twoPhaseField.h"
#include "equilibrium.h"
#include "concentrationField.h"
#include "ghostFlowField.h"
#include "electrodePotential.h"
#include "electrolytePotential.h"
#include "continue.h"
#include "initialWallDensityBoundary.h"
#include "wallBoundary.h"
#include "divergence.h"
#include "resultsStore.h"
#include "boundaryCondition.h"

#define threadsNumber 20
#define Length 64
#define Width 55
#define Height 64
#define pi 3.1415926
#define gasCritical 0.0
#define ConcentrationElectrolyteInitial 1.0
#define ConcentrationElectrodeInitial 7.97616256
#define electrodePotentialExchange 12.50
#define electrolytePotentialInitial 0.0
#define electrodePotentialInitial 52.4
//#define concentrationDiffusivity 0.0000066667
#define concentrationDiffusivity 0.000257
#define electrodeConductivity 0.5/3
#define solid 0.0
#define phase0 0.0
#define stepInitial 15000
#define interval 100
#define intervalPotential 1
#define intervalPotential2 1000
#define intervalElectrolyte 18
#define totalstep 16000

int main()
{
	clock_t start=clock();
	clock_t t3;
	clock_t t4;
	FILE *Initialfp;
	FILE *Initialfp1;
	/*************************************value for the operation of file***************************/
	int i;
	int j;
	int k;
	/********************************************value for loop*************************************/
	double avTime=0.0;
	int step;
	int stepElectrolyte;
	double totalCurrentResidual = 1000.0;
	double lastTotalCurrentResidual = 0.0;
	double totalCurrent = 0.0;
	double lastTotalCurrent = 0.0;
	double totalCurrentError = 1000.0;
	double lastTotalCurrentError = 0.0;
	double CrelaxationTime=1.0/(concentrationDiffusivity*4.0+0.5);
	double electrodeRelaxationTime=1.0/(electrodeConductivity*4.0+0.5);
	double constantMol = 1000.0;
	double constantLelectrolyteConductivity = 0.0135;
	double constantLelectrolyteDiffusivity = 675680000;
	/**********************************************************************************************************************/
	double FaradayConstant=96485.0;
	double universalGasConstant=8.314;
	double temperatureReaction=298.0;
	double reactionRateConstant = 0.0001;
	double K=298.0;
	double deltaElectrodePotential = -0.0000065875*5*4;
	double electrolyteCurrentDensityBoundary = 0.0000065875*5*4;
	double latticeF = 0.0048;
	double MaxConcentration=36.224;
	double transferenceNumber = 0.4;
	/**************************************************Bottom surface reaction***********************************************/
	omp_set_num_threads(threadsNumber);
	/****************************************arrange the number of threads for parallel running**************************************************/
	printf("reaction rate is:%0.8lf\n",reactionRateConstant);
	domainArrange(Length, Width, Height);
	OutMemoryArrange(Length, Width, Height);
	FieldArrangeTwoPhase(Length, Width, Height);
	FieldArrangeConcentration(Length, Width, Height);
	FieldArrangeConcentrationFirst(Length, Width, Height);
	FieldArrangeElectrodePotential(Length, Width, Height);
	FieldArrangeelectrolytePotential(Length, Width, Height);
	/********************************************************************************************************************************************/
	domainDefine(0, Length, 0, Width, 0, Height, phase0);
	Initialfp = fopen("C:\\duo\\lithium-Ion battery\\simulationForPaper\\ExperimentStructure\\structure\\electrodeStructure.txt", "r");
	dataload(Initialfp, domainPorous, 0, Length, 0, Width, 0, Height);
	for (i=0;i<Length;i++){
		for (j=0;j<Width;j++){
			for (k=0;k<Height;k++){
				domain[i][j][k]=domainPorous[i][j][k];
			}
		}
	}
	domainShift(Length, Width, Height);                                     //shift domain for justify the wall boundary layer
	/*************************************************load domain*******************************************************************/
	for (i = 0; i<Length; i++){
		for (j = 0; j<Width; j++){
			for (k = 0; k<Height; k++){
				rho[i][j][k] = 0.3;
			}
		}
	}
	/***********************************************load rho domain*****************************************************************/
	fileAddressConcentration();
	fileAddressConcentrationFirst();
	fileAddressElectrodePotential();
	fileAddressElectrolytePotential();
	ContinueLoadConcentration(Length,Width,Height);
	ContinueLoadConcentrationFirst(Length,Width,Height);
	ContinueLoadElectrodePotential(Length,Width,Height);
	ContinueLoadElectrolytePotential(Length, Width, Height);
	/************************************************continue Initial****************************************************************/
	domainConnectionRegionD3Q7(Length, Width, Height,gasCritical);
	electrodeConnectionRegionD3Q7(Length, Width, Height);
	/**********************************************connection of electrode and pore space********************************************/
/*	FieldInitialConcentrationElectrolyteD3Q7(Length, Width, Height, ConcentrationElectrolyteInitial, gasCritical);
	FieldInitialConcentrationElectrodeD3Q7(Length, Width, Height, ConcentrationElectrodeInitial, gasCritical);
	FieldInitialelectrodePotentialD3Q7(Length, Width, 0, Height, electrodePotentialInitial);
	FieldInitialelectrolytePotentialD3Q7(Length, Width, 0, Height, electrolytePotentialInitial, gasCritical);
	/******************************************************************initialization***********************************************/
	rhoShift(Length, Width, Height);                                                     //get rho of neighbour

	step=stepInitial;
	while (step<totalstep){
		t3=clock();
		electrolyteRelaxationTimeCalculation(Length, Width, Height, K, constantMol / 1000.0, constantLelectrolyteConductivity, gasCritical);
		additionalTermCalculationD3Q7(Length, Width, Height, universalGasConstant, FaradayConstant, K, gasCritical, constantMol, constantLelectrolyteConductivity, electrodePotentialExchange);
		electrolyteConcentrationRelaxationTimeCalculation(Length, Width, Height, K, gasCritical, constantMol / 1000.0, constantLelectrolyteDiffusivity);
		/***********************************************************************************************************************************/
		lastConcentration(Length, Width, Height, gasCritical);
		/*******************************************************save concentration from last step***********************************************************************/
		while ((fabs(totalCurrentError)>0.02 || fabs(totalCurrentResidual)>0.001) && (step%intervalPotential == 0)){
			lastElectrode(Length, Width, Height);
			lastElectrolyte(Length, Width, Height, gasCritical);
			totalCurrent = 0.0;
			if (step%intervalPotential2 == 0){
				lastTotalCurrentError = lastTotalCurrentError / 50.0;
			}
			/***************************************************************************************************************************************************************/
			CollisionSRTElectrodePotentialD3Q7(Length, Width, Height, electrodeRelaxationTime);
			reactionElectrodeD3Q7(Length, Width, Height, gasCritical, reactionRateConstant, latticeF, K, electrodePotentialExchange, MaxConcentration, lastTotalCurrentError*50.0);
			StreamElectrodePotentialD3Q7(Length, Width, Height);
			wallBoundaryElectrodePotentialD3Q7(Length, Width, Height);
			FieldCalculationElectrodePotentialD3Q7(Length, Width, Height);
			widthLastDirichletelectrodePotentialD3Q7(Length, Width, 0, Height, deltaElectrodePotential);
			/*************************************************electrodePotential Field collision and stream***********************************************************************/
			CollisionSRTelectrolytePotentialD3Q7(Length, Width, Height, gasCritical);
			reactionElectrolyteD3Q7(Length, Width, Height, lastTotalCurrentError*50.0);
			StreamelectrolytePotentialD3Q7(Length, Width, Height);
			wallBoundaryElectrolytePotentialD3Q7(Length, Width, Height, gasCritical);
			gasLiquidBoundaryElectrolytePotentialD3Q7(Length, Width, Height, gasCritical);
			FieldCalculationelectrolytePotentialD3Q7(Length, Width, Height, gasCritical);
			widthFirstConstantCurrentElectrolytePotentialD3Q7(Length, Width, 0, Height, electrolyteCurrentDensityBoundary, gasCritical);
			/*************************************************electrolytePotential Field collision and stream***********************************************************************/
			for (i = 0; i < Length; i++){
				for (j = 0; j < Width; j++){
					for (k = 0; k < Height; k++){
						totalCurrent = totalCurrent + localCurrent[i][j][k];
					}
				}
			}
			totalCurrentResidual = fabs(totalCurrent -lastTotalCurrent) / totalCurrent;
			totalCurrentError = (totalCurrent + deltaElectrodePotential*Length*Height) / -deltaElectrodePotential / Length / Height;
			printf("Residual of total Current: %f\n", totalCurrentResidual);
			printf("Error of total Current: %f\n", totalCurrentError);
			printf("progress: %f%%\n", 1.0*step / totalstep * 100);                                                                     //record the progrmme running progress
			printf("\n");
			lastTotalCurrentResidual = totalCurrentResidual;
			lastTotalCurrent = totalCurrent;
			lastTotalCurrentError = totalCurrentError;
			/*************************************************************************************************************************************************************************/
		}
		stepElectrolyte = 0;
		while (stepElectrolyte < intervalElectrolyte){
			CollisionSRTConcentrationD3Q7(Length, Width, Height, gasCritical);
			reactionConcentrationD3Q7(Length, Width, Height, gasCritical, latticeF*intervalElectrolyte, transferenceNumber, lastTotalCurrentError);
			StreamConcentrationD3Q7(Length, Width, Height);
			wallBoundaryConcentrationD3Q7(Length, Width, Height);
			FieldCalculationConcentrationD3Q7(Length, Width, Height, gasCritical);
			widthFirstConstantCurrentConcentrationD3Q7(Length, Width, 0, Height, electrolyteCurrentDensityBoundary / latticeF / intervalElectrolyte, gasCritical, transferenceNumber);
			stepElectrolyte = stepElectrolyte + 1;
		}
		/*************************************************concentration Field collision and stream***********************************************************************/
//		CollisionSRTConcentrationElectrodeD3Q7(Length, Width, Height, CrelaxationTime, gasCritical);
		CollisionMRTConcentrationElectrodeD3Q7(Length, Width, Height);
		reactionConcentrationElectrodeD3Q7(Length, Width, Height, gasCritical, latticeF, lastTotalCurrentError);
		StreamConcentrationElectrodeD3Q7(Length, Width, Height);
		wallBoundaryConcentrationElectrodeD3Q7(Length, Width, Height);
		FieldCalculationConcentrationElectrodeD3Q7(Length, Width, Height, gasCritical);
//		widthLastElectrodeConcentrationNeumannD3Q7(Length, Width, Height);
		/*************************************************concentrationFirst Field collision and stream***********************************************************************/		
		totalCurrentResidual = 1000.0;
		totalCurrentError = 1000.0;
		/*****************************************************************************************************************************************************************/
		t4=clock();
		avTime=avTime+(double)(t4-t3)/CLOCKS_PER_SEC;
		step = step + 1;
		printf("progress: %f%%\n", 1.0*step / totalstep * 100);
		resultsWriteConcentration(Length,Width,Height,step,interval);
		resultsWriteConcentrationFirst(Length,Width,Height,step,interval);
		resultsWriteElectrodePotential(Length,Width,Height,step,interval);
		resultsWriteElectrolytePotential(Length,Width,Height,step,interval);
		resultsWriteReactantSurface(Length, Width, Height, step, interval);
		/**************************************************write out********************************************************************************************/
	}
	resultsWriteDomain(Length, Width, Height);
	avTime=avTime/totalstep;
	printf("Time elapsed: %f\n",((double)clock()-start)/CLOCKS_PER_SEC);                                                       // record the progrmme running time 
	printf("Time elapsed: %f\n",avTime);                                                                                       // record the running time for every step

	domainFreeMemory(Length, Width, Height);
	MemoryFreeTwoPhase(Length, Width, Height);
	MemoryFreeConcentration(Length, Width, Height);
	MemoryFreeConcentrationFirst(Length, Width, Height);
	MemoryFreeElectrodePotential(Length, Width, Height);
	MemoryFreeelectrolytePotential(Length, Width, Height);
	OutMemoryFree(Length,Width,Height);
	/********************************************freeout the shift memory****************************************************************************************/
}