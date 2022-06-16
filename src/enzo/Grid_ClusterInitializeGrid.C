/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A Cool Core Cluster
/
/  written by: Yuan Li and Greg Bryan
/  date:       Dec 2011 
/  modified1: Deovrat Prasad
/  date: Dec 2021
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
float gasdev();

static int ClusterParticleCount = 0;

static float CosmologySimulationInitialFractionHII   = 1.2e-5;
static float CosmologySimulationInitialFractionHeII  = 1.0e-14;
static float CosmologySimulationInitialFractionHeIII = 1.0e-17;
static float CosmologySimulationInitialFractionHM    = 2.0e-9;
static float CosmologySimulationInitialFractionH2I   = 2.0e-20;
static float CosmologySimulationInitialFractionH2II  = 3.0e-14;

int grid::ClusterInitializeGrid(int NumberOfSpheres,
			     FLOAT SphereRadius[MAX_SPHERES],
			     FLOAT SphereCoreRadius[MAX_SPHERES],
			     float SphereDensity[MAX_SPHERES],
			     float SphereTemperature[MAX_SPHERES],
			     FLOAT SpherePosition[MAX_SPHERES][MAX_DIMENSION],
			     float SphereVelocity[MAX_SPHERES][MAX_DIMENSION],
			     int   SphereType[MAX_SPHERES],
			     int   SphereUseParticles,
			     float UniformVelocity[MAX_DIMENSION],
			     int   SphereUseColour,
			     float InitialTemperature, 
			     float ClusterInitialSpinParameter,
			     int level)
{
  /* declarations */

  int dim, i, j, k, m, sphere;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

  /* create fields */

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  int ivel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1) 
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;
  if (MultiSpecies) {
    FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
    if (MultiSpecies > 1) {
      FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
      FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
      FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
    }
    if (MultiSpecies > 2) {
      FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
    }
  }
  int ColourNum = NumberOfBaryonFields;
  if (SphereUseColour)
    FieldType[NumberOfBaryonFields++] = Metallicity; /* fake it with metals */

  /* Set various units. */

  const double keV=1.1604e7;
  //double c200 = 5.0; //A2029
  double c200 = 4.7; // Phoenix
  //double c200 = 9.5; //N5044
  //double c200 = 7.5; //N4472
  float DensityUnits, LengthUnits, TemperatureUnits = 1, TimeUnits, 
    VelocityUnits, CriticalDensity = 1, BoxLength = 1, mu = 0.6;
  double MassUnits = 1;
  FLOAT a, dadt, ExpansionFactor = 1;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, &MassUnits, Time);
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
    ExpansionFactor = a/(1.0+InitialRedshift);
    CriticalDensity = 2.78e11*POW(HubbleConstantNow, 2); // in Msolar/Mpc^3
    BoxLength = ComovingBoxSize*ExpansionFactor/HubbleConstantNow;  // in Mpc
  } else {
    CriticalDensity = 2.78e11*pow(0.74,2); // in Msolar/Mpc^3 for h=0.74
    BoxLength = LengthUnits / Mpc_cm;
    HubbleConstantNow = 1.0;
    OmegaMatterNow = 1.0;
  }


  /* Set densities */

  float BaryonMeanDensity = SphereUseParticles ? 0.1 : 1.0;
  if (SphereUseParticles == 2) BaryonMeanDensity = 0.9;
  float ParticleMeanDensity = 1.0 - BaryonMeanDensity, ParticleCount = 0;

  /* Set the point source gravity parameters for the NFW profile. */
  if (PointSourceGravity == 2) {
    PointSourceGravityCoreRadius = SphereCoreRadius[0]*LengthUnits; // in CGS
//printf("begin calculating PointSourceGravityConstant");
    PointSourceGravityConstant = 4.0*pi*SphereDensity[0]*DensityUnits *
                pow(SphereCoreRadius[0]*LengthUnits, 3) *
               (log(1.0+c200) - c200/(1.0+c200))/SolarMass;// in Msun. Not Mvir, but Ms.
    BaryonMeanDensity = 0.15; // 15% baryon fraction
//printf("PointSourceGravityConstant= %"GSYM"\n", PointSourceGravityConstant);
  }
  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    NumberOfParticles = (SphereUseParticles > 0) ? 1 : 0;
    for (dim = 0; dim < GridRank; dim++)
      NumberOfParticles *= (GridEndIndex[dim] - GridStartIndex[dim] + 1);
    return SUCCESS;
  }

  /* Compute NFW profile-info. The parameters are SphereCoreRadius (which is
     the "knee radius") and SphereDensity (the overdensity ar the knee
     radius).  The pressure is computed by integrated the hydro-static
     equation and from this comes the temperature and the dm velocity
     dispersion. */

#define NFW_POINTS 2000
  float Allg[NFW_POINTS], NFWg[NFW_POINTS], NFWPressure[NFW_POINTS], NFWTemp[NFW_POINTS], x1,
        NFWDensity[NFW_POINTS], GasDensity[NFW_POINTS],NFWMass[NFW_POINTS], NFWSigma[NFW_POINTS];
  float NFWRadius[NFW_POINTS], Totp[NFW_POINTS], Kent[NFW_POINTS];
  double dpdr = 0, dpdr_old, rkpc, pg, fn, fnp, err; // Newton-Raphson variables - Deovrat
  //double a_0 = 20.0*3.086e21, rho_0 = 4.0*1.67e-24; //A2029
  double a_0 = 20.0*3.086e21, rho_0 = 1.0*1.67e-24; //Phoenix
  //double a_0 = 1.2*3.086e21 , rho_0 = 5.0*1.67e-22; // NGC5044
  //double a_0 = 1.6*3.086e21, rho_0 = 3.3*1.67e-22; // NGC4472

  double Konst, gmma=1.667, invgmma=0.6, sl=1.0e-10;
  sphere = 0;

  FILE *fptr = fopen("NFWProfile.out", "w");
  Konst = 0.617*mh*POW((1.17*mh),gmma-1.0)/(kboltz*1.16e7); // constant for Newton-Raphson
  printf("Before loop, c200=%"GSYM",a_0=%"GSYM",rho_0=%"GSYM",ClusterSMBHMass=%"GSYM"\n",c200, a_0, rho_0, ClusterSMBHMass);
  for (i = 0; i < NFW_POINTS; i++) {
    NFWRadius[i] = SphereRadius[sphere]*pow(10, -5*(float(i)/NFW_POINTS));
    x1 = NFWRadius[i]/SphereCoreRadius[sphere];
    NFWDensity[i] = SphereDensity[sphere]/(x1*(1.0+x1)*(1.0+x1));    // DM Density
    if (SphereType[sphere]>=6 && SphereType[sphere] <= 8) {  //aka 6, 7, 8: Perseus Cluster
       //rkpc=NFWRadius[i]*LengthUnits/(1.0e-3*Mpc);
       Kent[i] = 10.0 + 165.0*pow(NFWRadius[i]*LengthUnits/(0.1*Mpc_cm), 0.95); //A2029
    /* Set Gravity. NFW Dark Matter, BCG+BH */
       //Allg[i]=GravConst*PointSourceGravityConstant*SolarMass*((log(1.0+x1)-x1/(1.0+x1))
       //      /(log(1.0+c200)-c200/(1.0+c200)))/POW(NFWRadius[i]*LengthUnits, 2.0)
       //       + 2.0*3.14159265*GravConst*rho_0*pow(a_0,3)/pow((a_0+NFWRadius[i])*LengthUnits, 2) // using Hernquist profile for BCG
       //       + GravConst*SolarMass* ClusterSMBHMass/POW(NFWRadius[i]*LengthUnits, 2);
//
       //Totp[i] = -(GravConst*PointSourceGravityConstant*SolarMass*log(1.0+x1)/(log(1.0+c200) - c200/(1.0+c200))/(NFWRadius[i]*LengthUnits + sl*Mpc_cm) //NFW component
       //	 + 2.0*3.14159265*GravConst*rho_0*pow(a_0,3)/(a_0+NFWRadius[i]*LengthUnits) //Hernquist potential - BCG
       //           + GravConst*ClusterSMBHMass*SolarMass/(NFWRadius[i]*LengthUnits - 2.0*GravConst*ClusterSMBHMass*SolarMass/POW(clight,2))); //SMBH potential
  }else if (SphereType[sphere]==1) {  //galaxy
     Kent[i] = 5.0 + 85.0*pow((NFWRadius[i]*LengthUnits/(0.1*Mpc_cm)),1.1); //Galaxy 2.0*10^12 Msun
  }else if (SphereType[sphere]==2) {  //NGC4472//NGC5044
     //Kent[i] = 1.3 + 150.0*pow((NFWRadius[i]*LengthUnits/(0.1*Mpc_cm)),1.05); //NGC5044
     Kent[i] = 1.5 + 400.0*pow((NFWRadius[i]*LengthUnits/(0.1*Mpc_cm)),1.05); //NGC4472
  }else {
           ENZO_FAIL("Entropy for this SphereType not defined.");
  }
        Allg[i]=GravConst*PointSourceGravityConstant*SolarMass*
                ((log(1.0+x1)-x1/(1.0+x1)) /(log(1.0+c200)-c200/(1.0+c200)))/pow(NFWRadius[i]*LengthUnits, 2.0) 
                + 2.0*3.14159265*GravConst*rho_0*pow(a_0,3)/pow((a_0+NFWRadius[i])*LengthUnits, 2)  //using Hernquist profile
                +GravConst*SolarMass* ClusterSMBHMass / POW(NFWRadius[i]*LengthUnits, 2)  ;
 
	Totp[i] =  -(GravConst*PointSourceGravityConstant*SolarMass*log(1.0+x1)/(log(1.0+c200) - c200/(1.0+c200))/(NFWRadius[i]*LengthUnits + sl*Mpc_cm) //NFW
       		 + 2.0*3.14159265*GravConst*rho_0*pow(a_0,3)/(a_0+NFWRadius[i]*LengthUnits) //Hernquist profile - stellar profile
                 + GravConst*ClusterSMBHMass*SolarMass/(NFWRadius[i]*LengthUnits - 2.0*GravConst*ClusterSMBHMass*SolarMass/POW(clight,2)))    ; // SMBH  profile
	//printf("c200=%"GSYM",a_0=%"GSYM",rho_0=%"GSYM",ClusterSMBHMass=%"GSYM"\n",c200, a_0, rho_0, ClusterSMBHMass);

    if (i==0){
	// This section has been completely modified by Deovrat Prasad (Nov, 2018)
        GasDensity[i] = 1.0e-29; // Gas density estimate at the outer boundary - free parameter
        //NFWPressure[i]=kboltz*NFWTemp[i]*GasDensity[i]/ (mu * mh);  // in cgs
        NFWPressure[i]=Kent[i]*pow(GasDensity[i],gmma)/Konst;
        NFWTemp[i] = NFWPressure[i]*0.62*mh/(GasDensity[i]*kboltz);
    }
    // Using the value of gas pressure at the outer boundary we solve the
    // hydrostatic equilibrium equation using Newton-Raphson to construct
     // the pressure (and consequently density) profile for the Galaxy Cluster.
    if (i>0) {
	 err = 1.0e20;
         pg = NFWPressure[i-1];
         while(err >=1.0e-10){
		 //fn = pg - NFWPressure[i-1] + (0.5*mu*mh/kboltz)*(Totp[i] - Totp[i-1])*(pg/NFWTemp[i] + NFWPressure[i-1]/NFWTemp[i-1]);
                 //fnp = 1.0 + (0.5*mu*mh/kboltz)*(Totp[i] - Totp[i-1])/NFWTemp[i];
                 fn = pg - NFWPressure[i-1] + 0.5*(Totp[i] - Totp[i-1])*pow(Konst,invgmma)*(pow(pg/Kent[i], invgmma)
                      + pow(NFWPressure[i-1]/Kent[i-1], invgmma));
                 fnp = 1.0 + 0.5*(Totp[i]- Totp[i-1])*invgmma*pow(pg,(invgmma-1.0))*pow(Konst/Kent[i], invgmma);
		 pg = pg - fn/fnp;
                 err = fabs(fn/fnp)/pg;
         }
         NFWPressure[i] = pg;
    }
     GasDensity[i] = pow(NFWPressure[i]*Konst/Kent[i], invgmma);
     NFWTemp[i] = NFWPressure[i]*0.62*mh/(GasDensity[i]*kboltz);
    
  fprintf(fptr, "%"ISYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" \n", i, NFWRadius[i],
         NFWDensity[i], Allg[i], NFWPressure[i], NFWTemp[i], GasDensity[i]);
}  //end for
  fprintf(fptr, "CriticalDensity = %"GSYM" , DensityUnits = %"GSYM", TimeUnits=%"GSYM", LengthUnits= %"GSYM"\n", CriticalDensity, DensityUnits, TimeUnits, LengthUnits);
  fclose(fptr);

  /* Loop over the set-up twice, once to count the particles, the second
     time to initialize them. */

  int SetupLoopCount, npart = 0;
  for (SetupLoopCount = 0; SetupLoopCount < 1+min(SphereUseParticles, 1);
       SetupLoopCount++) {


  /* Set up the baryon field. */


  /* allocate fields */

  if (SetupLoopCount == 0)
      this->AllocateGrids();

  /* Loop over the mesh. */

  float density, dens1, Velocity[MAX_DIMENSION],
    temperature, temp1, sigma, sigma1, colour, gas_density, gas_dens1;
  FLOAT r, x, y = 0, z = 0;
  int n = 0;

  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++, n++) {

        /* Compute position */

        x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
        if (GridRank > 1)
          y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
        if (GridRank > 2)
          z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

        /* Loop over spheres. */

        density = 1.0;                              //Set background DM density
        dens1   = 0.0;                              //for DM
        gas_density = density * BaryonMeanDensity;  //Set background gas density
        gas_dens1 = gas_density;                     //for gas
        temperature = temp1 = InitialTemperature;
        sigma = sigma1 = 0;
        colour = 1.0e-10;
        for (dim = 0; dim < MAX_DIMENSION; dim++)
          Velocity[dim] = 0;
        for (sphere = 0; sphere < NumberOfSpheres; sphere++) {

          /* Find distance from center. */

          r = sqrt(POW(fabs(x-SpherePosition[sphere][0]), 2) +
                   POW(fabs(y-SpherePosition[sphere][1]), 2) +
                   POW(fabs(z-SpherePosition[sphere][2]), 2) );
          r = max(r, 0.1*CellWidth[0][0]);

          if (r < SphereRadius[sphere]) {

              FLOAT xpos, ypos, zpos, vc, rz;
              x1 = r/SphereCoreRadius[sphere];
              dens1 = SphereDensity[sphere]/(x1*(1.0+x1)*(1.0+x1));
            /* 6 Perseus with self-gravity */
            /* 3, 7, 8 Perseus no gas self-gravity, force HSE */
            for (m = 1; m < NFW_POINTS; m++)
                if (r > NFWRadius[m]) {
                  temp1 = NFWTemp[m] + (NFWTemp[m-1] - NFWTemp[m])*
                    (r - NFWRadius[m])/(NFWRadius[m-1] - NFWRadius[m]);               // in K
                    if (SphereType[sphere] == 6) {
                     gas_dens1 =mh*(0.0192/(1.0+POW(r*LengthUnits/(18.0e-3*Mpc_cm),3.0))+0.046/pow((1.0+pow(r*LengthUnits/(57.0e-3*Mpc_cm), 2.0)), 1.8)+
                     0.0048/POW((1.0+pow(r*LengthUnits/(200.0e-3*Mpc_cm), 2.0)), 1.1))/DensityUnits/0.88;
                    }
                     else{
                     gas_dens1 = (GasDensity[m] + (GasDensity[m-1] - GasDensity[m])*
                      (r - NFWRadius[m])/(NFWRadius[m-1] - NFWRadius[m]))/DensityUnits;  // in code unit
                     }
                  break;  // break when NFWRadius just drops below r
            }

              /* Loop over dims if using Zeus (since vel's face-centered). */

              for (dim = 0; dim < 1+(HydroMethod == Zeus_Hydro ? GridRank : 0);
                   dim++) {

                /* Compute position. */

              xpos = x-SpherePosition[sphere][0] - (dim == 1 ? 0.5*CellWidth[0][0] : 0.0);
              ypos = y-SpherePosition[sphere][1] - (dim == 2 ? 0.5*CellWidth[1][0] : 0.0);
              zpos = z-SpherePosition[sphere][2] - (dim == 3 ? 0.5*CellWidth[2][0] : 0.0);

              vc = ClusterInitialSpinParameter*sqrt(GravConst*PointSourceGravityConstant*SolarMass/(PointSourceGravityCoreRadius)); /*in GCS unit*/

              rz = sqrt(POW(fabs(xpos), 2) + pow(fabs(ypos), 2));
              rz = max(rz, 0.1*CellWidth[0][0]);
              
              if (r > 6.25e-4) {  //10kpc
                if (dim == 0 || dim == 1)
                  Velocity[0] = (-ypos*vc/rz+gasdev()*1.0e5*SphereVelocity[0][dim])/VelocityUnits; 
                if (dim == 0 || dim == 2)
                  Velocity[1] = (xpos*vc/rz+gasdev()*1.0e5*SphereVelocity[0][dim])/VelocityUnits;
                if (dim == 0 || dim == 3)
                  Velocity[2] = gasdev()*1.0e5*SphereVelocity[0][dim]/VelocityUnits;
              } 
              else {
                if (dim == 0 || dim == 1)
                  Velocity[0] = (-ypos*vc/rz)/VelocityUnits; 
                if (dim == 0 || dim == 2)
                  Velocity[1] = (xpos*vc/rz)/VelocityUnits;
                if (dim == 0 || dim == 3)
                  Velocity[2] = 0;                    
              } 
              } // end: loop over dims

            /* If the density is larger than the background (or the previous

               sphere), then set the velocity. */


            if (dens1 > density) {
              gas_density = gas_dens1;
              if (temp1 == InitialTemperature)
                temp1 = SphereTemperature[sphere];
              temperature = temp1;
            }

          } // end: if (r < SphereRadius)
        } // end: loop over spheres



        /* Set density. */
    BaryonField[0][n] = gas_density;


	/* Set Velocities. */

	for (dim = 0; dim < GridRank; dim++)
	  BaryonField[ivel+dim][n] = Velocity[dim] + UniformVelocity[dim];

	/* Set energy (thermal and then total if necessary). */

	BaryonField[1][n] = temperature/TemperatureUnits/
                            ((Gamma-1.0)*mu);

	if (DualEnergyFormalism)
	  BaryonField[2][n] = BaryonField[1][n];

	if (HydroMethod != Zeus_Hydro)
	  for (dim = 0; dim < GridRank; dim++)
	    BaryonField[1][n] += 0.5*POW(BaryonField[ivel+dim][n], 2);

      } // end loop over grid

  } // end loop SetupLoopCount

  return SUCCESS;
}

