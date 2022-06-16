/***********************************************************************
/
/ Gets the accretion rate for a jet powered by Bondi accretion. To remove
/ mass from the accretion zone, the accretion rate is divided by the total mass in
/ the zone, and the density in each cell is multiplied by this amount. Basically,
/ the mass is being removed in a density weighted manner.
/ Returns the accretion rate in code units
/
/  written by: Greg Meece
/  date:       April 2015
/
************************************************************************/

#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "units.h"
#include "Fluxes.h"
#include "GridList.h"
#include "phys_constants.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "ActiveParticle.h"
#include "ActiveParticle_AGNParticle.h"

#define NO_DEBUG_AP

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

float grid::AGNParticleGetBondiMassRate(ActiveParticleType* ThisParticle) {

   /* Return if this doesn't involve us */
   if (MyProcessorNumber != ProcessorNumber)
      return SUCCESS;

   // Cast this particle to an AGN particle so that we can access AGN specific
   // properties.
   ActiveParticleType_AGNParticle* tp = static_cast <ActiveParticleType_AGNParticle*>(ThisParticle);

   float max_radius = max(tp -> CoolingRadius, tp -> FeedbackRadius);

   float xsink = tp -> pos[0];
   float ysink = tp -> pos[1];
   float zsink = tp -> pos[2];

   if ((GridLeftEdge[0]    > xsink + max_radius) ||
         (GridLeftEdge[1]  > ysink + max_radius) ||
         (GridLeftEdge[2]  > zsink + max_radius) ||
         (GridRightEdge[0] < xsink - max_radius) ||
         (GridRightEdge[1] < ysink - max_radius) ||
         (GridRightEdge[2] < zsink - max_radius))
      return SUCCESS;

   // Figure out how to access grid quantities
   int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;

   if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                        Vel3Num, TENum) == FAIL) {
      ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
      }

   // Get the units
   float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
         VelocityUnits = 1, TimeUnits = 1, aUnits = 1;
   float MassUnits;

   GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
            &TimeUnits, &VelocityUnits, Time);
   MassUnits = DensityUnits * pow(LengthUnits, 3.0);

   float cell_volume = pow(CellWidth[0][0], 3.0); //Code
   float radius; //Code

   int index;

   // Some cgs constants
   float cm_per_kpc = 3.08567758e21;
   float seconds_per_year = 365.25 * 24.0 * 3600.0;
   float c_cgs = 2.99792458e+10;
   float g_cgs = 6.67259e-8;

   float g_code = g_cgs * DensityUnits * TimeUnits * TimeUnits;

   float cell_mass;
   float total_mass; // Total mass in the accretion zone
   float sound_speed; // Mass weighted average sound speed
   float velocity; // Mass weighted average velocity magnitude
   float temp, cs;
   float mean_density;
   float vx, vy, vz;
   float mdot;

   // Calculate and sum Mdot over all cells
   total_mass = 0.0;
   velocity = 0.0;
   sound_speed = 0.0;
   mdot = 0.0;

   // How much mass is removed during this time step?
   float delta_mass;
   delta_mass = mdot * dtFixed;

   // Calculate the temperature
   float* temperature = new float[GridDimension[0] * GridDimension[1] * GridDimension[2]];
   ComputeTemperatureField(temperature);


   // Add up the total mass, sound speed, and velocity within the sphere
   for (int k = GridStartIndex[2]; k < GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j < GridEndIndex[1]; j++) {
         for (int i = GridStartIndex[0]; i < GridEndIndex[0]; i++) {
            index = k * GridDimension[1] * GridDimension[0]
                  + j * GridDimension[0] + i;

            radius = pow(CellLeftEdge[0][i] + CellWidth[0][i] * 0.5 - xsink, 2.0)
                   + pow(CellLeftEdge[1][j] + CellWidth[1][j] * 0.5 - ysink, 2.0)
                   + pow(CellLeftEdge[2][k] + CellWidth[2][k] * 0.5 - zsink, 2.0);
            radius = sqrt(radius);

            //if (radius < tp -> CoolingRadius) {
	    if (radius < tp -> BondiRadius) {
               // Cell mass
               cell_mass = BaryonField[DensNum][index] * cell_volume; //code
               total_mass += cell_mass;

               // Velocity magnitude
               vx = BaryonField[Vel1Num][index];
               vy = BaryonField[Vel2Num][index];
               vz = BaryonField[Vel3Num][index];

               velocity += sqrt(vx*vx + vy*vy + vz*vz) * cell_mass;

               // Sound speed
               if (DualEnergyFormalism) {
                  temp = BaryonField[GENum][index] * (TemperatureUnits * (Gamma - 1.0) * Mu); // Temp in K
                  }

               else {
                  temp = BaryonField[TENum][index] - 0.5 * (vx*vx + vy*vy + vz*vz);
                  temp *= (TemperatureUnits * (Gamma - 1.0) * Mu);
                  }

               temp = temperature[index];

               if (temp < 0.0) {
                  printf("error:: temp, TE, vx, vy, vz: %"GSYM", %"GSYM", %"GSYM", %"GSYM", %"GSYM"\n", temp, BaryonField[TENum][index], vx, vy, vz);
                  exit(1);
                  }

               cs = sqrt(Gamma * kboltz * temp / (Mu * mh)); // Sounds speed this cell, cgs
               cs /= (LengthUnits / TimeUnits); // Code units
               sound_speed += cs * cell_mass;

               } // End total mass calculation

            } // End loop over i
         } // End loop over j
      } // End loop over k

   delete [] temperature;


   // Finish taking mass weighted averages
   sound_speed /= total_mass;
   velocity /= total_mass;
   mean_density = (3.0 * total_mass) / (4.0 * M_PI * pow(tp -> BondiRadius, 3.0));

   // Calculate the bondi accretion rate
   float mbh = 5.0e+10 * 1.99e+33;
   mbh /= MassUnits;

   mdot = (4.0 * M_PI * g_code * g_code * mbh * mbh * mean_density)
          / pow(sound_speed, 3.0);

   // Case 1: Constant boost factor
   if (AGNParticleFeedbackType == 3) {
      mdot *= AGNParticleBondiBoostFactor;
      }

   // Case 2: Booth and Shay: Density dependent boost factor
   else if (AGNParticleFeedbackType == 6) {
      float boost_factor, rho_0;

      rho_0 = AGNParticleBSBondiDensityCGS / DensityUnits;

      if (mean_density <= rho_0)
         boost_factor = 1.0;
      else
         boost_factor = pow(mean_density/rho_0, AGNParticleBSBondiBeta);

      mdot *= boost_factor;

      printf("rho, rho_0, boost_factor: %"FSYM" %"FSYM" %"FSYM"\n", mean_density, rho_0, boost_factor);
      
   } else if ( AGNParticleBondiAccretion == 1){ //added by DP to do Bondi in addition to cold mode.
	  mdot *= 1.0; 

   // Case 3: Bad feedback type selected.
   } else {
      printf("Error! Feedback type %"ISYM" is not Bondi accretion.\n", AGNParticleFeedbackType);
      exit(1);
      }

   //printf("mdot: %"GSYM" (cgs: %"GSYM")\n", mdot, mdot*(MassUnits/TimeUnits));
   //printf("total_mass: %"GSYM"\n", total_mass);
   //printf("velocity, sound_speed, mbh, rho: %"GSYM", %"GSYM", %"GSYM", %"GSYM"\n", velocity * (LengthUnits/TimeUnits), sound_speed * (LengthUnits/TimeUnits), mbh * MassUnits, mean_density * DensityUnits);

   // Calculate the fraction of the gas that condenses
   delta_mass = mdot * dtFixed;
   float condensation_fraction = delta_mass / total_mass;

   // Remove mass from the accretion zone
   for (int k = GridStartIndex[2]; k < GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j < GridEndIndex[1]; j++) {
         for (int i = GridStartIndex[0]; i < GridEndIndex[0]; i++) {
            index = k * GridDimension[1] * GridDimension[0]
                  + j * GridDimension[0] + i;

            radius = pow(CellLeftEdge[0][i] + CellWidth[0][i] * 0.5 - xsink, 2.0)
                   + pow(CellLeftEdge[1][j] + CellWidth[1][j] * 0.5 - ysink, 2.0)
                   + pow(CellLeftEdge[2][k] + CellWidth[2][k] * 0.5 - zsink, 2.0);
            radius = sqrt(radius);

            if (radius < tp -> BondiRadius) {
               // Condensed mass drops out
               BaryonField[DensNum][index] *= (1.0 - condensation_fraction);
               } // End total mass calculation

            } // End loop over i
         } // End loop over j
      } // End loop over k

   return mdot;
}

