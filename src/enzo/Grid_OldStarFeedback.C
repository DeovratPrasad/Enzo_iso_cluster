/***********************************************************************
/
/  GRID: ADD Stellar Feedback to Elliptical Galaxies (Hernquist profile)
/
/  written by: Yuan Li
/  date:       Aug, 2015
/  modified1: 
/
/  PURPOSE: Stellar wind adds mass, and SN Ia adds heat
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, double *MassUnits, FLOAT Time);

int grid::OldStarFeedback()
{
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
                                             Vel3Num, TENum) == FAIL)   ///this or thisgrid
     ENZO_FAIL("Error in IdentifyPhysicalQuantities.");


//ClusterSMBHBCG is M_* here
//EllipticalGalaxyRe is Re in kpc
//OldStarFeedbackAlpha  alpha -19


  int i, j, k;
  float a=0;   

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1,
    TimeUnits = 1.0, VelocityUnits = 1.0;
  double MassUnits=1.0 ;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  float EnergyUnits;
  EnergyUnits = MassUnits* POW(VelocityUnits, 2.0);

  a=(EllipticalGalaxyRe*kpc_cm/1.8153)/LengthUnits;  // in code unit


  FLOAT r, x, y = 0, z = 0, rho_star=0;

   //double a_0 = 1.6*3.086e21, rho_0 = 3.3*1.67e-22; //NGC4472
   //double a_0 = 20.0*3.086e21, rho_0 = 4.0*1.67e-24; //A2029
   double a_0 = 20.0*3.086e21, rho_0 = 1.0*1.67e-24; //Phoenix
   //double a_0 = 1.2*3.086e21 , rho_0 = 5.0*1.67e-22; // NGC5044
   float cell_volume, cell_mass, ge_cell, ke_cell;
   float VelocityDispersionInKmPerSec= 300.0;

  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++) {

        /* Compute position */

        x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
        if (GridRank > 1)
          y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
        if (GridRank > 2)
          z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

          /* Find distance from center (0.5, 0.5, 0.5). */

        r = sqrt(POW(fabs(x-0.5), 2) +
                   POW(fabs(y-0.5), 2) +
                   POW(fabs(z-0.5), 2) );
        r = max(r, 0.1*CellWidth[0][0]);

	// Added by Deovrat Prasad.
         rho_star= rho_0/((r*LengthUnits/a_0)*pow((1.0 +(r*LengthUnits/a_0)),3)); // Hernquist density profile
         cell_volume = CellWidth[0][i]*CellWidth[1][j]*CellWidth[2][k];
         cell_mass = BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)]*cell_volume;
 
         ge_cell = (0.000475*rho_star*(dtFixed*TimeUnits)*SNIaFeedbackEnergy  //SNIa
                    +rho_star*(dtFixed*TimeUnits)*OldStarFeedbackAlpha*1.0e-19*POW(VelocityDispersionInKmPerSec*1.0e5,2)); //Stellar                                            
         ke_cell = 0.5*cell_mass*(pow(BaryonField[Vel1Num][GRIDINDEX_NOGHOST(i,j,k)],2)
                                +pow(BaryonField[Vel2Num][GRIDINDEX_NOGHOST(i,j,k)],2)
                                +pow(BaryonField[Vel3Num][GRIDINDEX_NOGHOST(i,j,k)],2));
 
	 // For PPM and DEF --- Deovrat 
         if (HydroMethod == PPM_DirectEuler){
                 if (DualEnergyFormalism)
			 BaryonField[GENum][GRIDINDEX_NOGHOST(i,j,k)] += ge_cell/(EnergyUnits*BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)]);
 
                 BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)] += (ge_cell+ke_cell)/(EnergyUnits*BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)]);
         //For ZEUS_Hydro -- Yuan
         }else if(HydroMethod == Zeus_Hydro){
                 BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)] += ge_cell/(EnergyUnits*BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)]);
         } else{
                 ENZO_FAIL("HydroMethod not defined for Old Star Feedback.");
         }
        BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)] += rho_star*(dtFixed*TimeUnits)*OldStarFeedbackAlpha*1.0e-19/DensityUnits;
  }

  return SUCCESS;

}
 
