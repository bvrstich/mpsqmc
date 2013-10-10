#ifndef MPSQMC_H
#define MPSQMC_H

#include "MPSstate.h"
#include "MPO.h"
#include "DMRG.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on September 4, 2013 */

class MPSQMC{

   public:
   
      //Constructor
      MPSQMC(MPO * theMPO, Random * RN, const int Dtrunc, const int Nwalkers, const double dtau);
      
      //Destructor
      ~MPSQMC();
      
      //Setup the trial wfn
      void SetupTrial();
      
      //Setup the walkers
      void SetupWalkers();
      
      //The energy based on projecting the set of walkers onto the trial wfn
      double EnergyFunction();
      
      //Let the walkers propagate for steps steps
      void Walk(const int steps);
      
   private:
   
      /*************************
      *** Global information ***
      *************************/
   
      //The problem MPO
      MPO * theMPO;
      
      //The number of terms in the MPO
      int nMPOterms;
      
      //The random number generator
      Random * RN;
      
      //The MPS truncation dimension
      int Dtrunc;
      
      //The total desired number of walkers
      int totalNDesiredWalkers;
      
      //The imaginary time step size (positive)
      double dtau;
      
      /*********************************************************************
      *** Some functions for propagation, population control and writing ***
      *********************************************************************/
      
      //Propagate my population of walkers for 1 timestep. Return the sum of the coeff of my walkers.
      double PropagateSeparately();
      
      //Add the current walker energy to walkerEnergyHistorySum[walker] & return the fluctuation metric for the current step
      double updateEnergyHistory(const int step);
      
      //Control the population of walkers based on scaling * weight
      void SeparatePopulationControl(const double scaling);
      
      //Write the projected energy, target energy, and fluctuation metric at MC time "step"
      void write(const int step, const double projectedEnergy, const double targetEnergy, const double fluctMetric);
      
      /****************************
      *** Trial wfn information ***
      ****************************/
      
      //Whether or not the trial wfn was set up
      bool bSetupTrial;
      
      //Trial wfn (one per thread)
      MPSstate ** Psi0;
      
      //MPO times trial wfn (one per thread)
      MPSstate ** HPsi0;
      
      //Specific MPO terms times trial wfn (one array per thread)
      MPSstate *** MPOtermsPsi0;
      
      /****************************************************************************************************
      *** Arrays to store the walkers, their coefficients, their overlaps, and their energy history sum ***
      ****************************************************************************************************/
      
      //Whether or not the walkers were set up
      bool bSetupWalkers;
      
      //The walkers
      MPSstate ** theWalkers;
      
      //The walkers: copy array
      MPSstate ** theWalkersCopyArray;
      
      //The coefficients
      double * walkerCoeff;
      
      //The coefficients: copy array
      double * walkerCoeffCopy;
      
      //The walker overlaps with the trial wfn
      double * walkerOverlap;
      
      //The walker overlaps with the trial wfn: copy array
      double * walkerOverlapCopy;
      
      //The sum of the energies at each previous time step per walker
      double * walkerEnergyHistorySum;
      
      //The sum of the energies at each previous time step per walker: copy array
      double * walkerEnergyHistorySumCopy;
      
      /***************************************************
      *** For MPI: some helper variables and functions ***
      ***************************************************/
      
      //My MPI rank
      int MPIrank;
      
      //The MPI size
      int MPIsize;
      
      //The current number of walkers in total
      int totalNCurrentWalkers;
      
      //The desired number of walkers on this process
      int myNDesiredWalkers;
      
      //The current number of walkers on this process
      int myNCurrentWalkers;
      
      //The max. allowed number of walkers on this process
      int myMaxNWalkers;
      
      //The max. number of openMP threads on this process
      int myNOpenMPthreads;
      
      //Workspace to temporarily store the overlaps <Psi_T | B(x) | Psi_i> in (one array per thread)
      double ** thePDF;
      
      //To keep track of the sum of the walker coeff
      double * sumWalkerCoeffPerThread;
      
      //The MPSstate is allocated and filled at MPIrank 0, and needs to be allocated and copied at the other ranks
      MPSstate * BroadcastCopyConstruct(MPSstate * pointer);
      
      //The desired walkers per rank
      int * NDesiredWalkersPerRank;
      
      //To control the populations on each rank
      void PopulationBalancing();
      
      //BubbleSort algorithm. Modifies order so that for all index: values[ order[index] ] >= values[ order[index+1] ].
      void BubbleSort(double * values, int * order, const int length);
      
};

#endif
