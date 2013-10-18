#ifndef MPSQMC2_H
#define MPSQMC2_H

#include "MPSstate.h"
#include "HeisenbergMPO.h"
#include "DMRG.h"
#include "TrotterHeisenberg.h"
#include "GridGenerator.h"
#include "Walker.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 15, 2013 */

class MPSQMC2{

   public:
   
      //Constructor
      MPSQMC2(HeisenbergMPO * theMPO, GridGenerator * theGrid, Random * RN, const int Dtrunc, const int Nwalkers, const double dtau);
      
      //Destructor
      ~MPSQMC2();
      
      //Let the walkers propagate for steps steps
      void Walk(const int steps);
      
   private:
   
      /******************************
      *** Constructor information ***
      ******************************/
   
      //The problem MPO
      HeisenbergMPO * theMPO;
      
      //The Trotter decomposition of the MPO
      TrotterHeisenberg * theTrotter;
      
      //The grid generator
      GridGenerator * theGrid;
      
      //The random number generator
      Random * RN;
      
      //The MPS truncation dimension
      int Dtrunc;
      
      //The total desired number of walkers
      int totalNDesiredWalkers;
      
      //The current number of walkers in total
      int totalNCurrentWalkers;
      
      //The imaginary time step size (positive)
      double dtau;
      
      /*********************************************************************
      *** Some functions for propagation, population control and writing ***
      *********************************************************************/
      
      //Propagate my population of walkers for 1 timestep. Return the sum of the coeff of my walkers.
      double PropagateSeparately();
      
      //Control the population of walkers based on scaling * weight
      void SeparatePopulationControl(const double scaling);
      
      //The MPSstate is allocated and filled at MPIrank 0, and needs to be allocated and copied to the other ranks
      MPSstate * BroadcastCopyConstruct(MPSstate * pointer);
      
      //To control the populations on each rank
      void PopulationBalancing();
      
      //Calculate the single walker projected energies, update the energy history, calculate the fluctuation metric, and the total projected energy
      double EnergyFunctionAndHistory(const int step, double * projectedEnergy, bool doFluctuationMetric=true);
      
      //BubbleSort algorithm. Modifies order so that for all index: values[ order[index] ] >= values[ order[index+1] ].
      void BubbleSort(double * values, int * order, const int length);
      
      //Write the projected energy, target energy, and fluctuation metric at MC time "step"
      void write(const int step, const double projectedEnergy, const double targetEnergy, const double fluctMetric);
      
      /****************************
      *** Trial wfn information ***
      ****************************/
      
      //Trial wfn (one per thread)
      MPSstate ** Psi0;
      
      //MPO times trial wfn (one per thread)
      MPSstate ** HPsi0;
      
      //Specific Trotter terms (hermitian conjugate!!) times trial wfn: TrotterTermsTimesPsi0[thread][coupling][ firstOp + d*d* secondOp ]
      MPSstate **** TrotterTermsTimesPsi0;
      
      //Number of non-zero couplings
      int nCouplings;
      
      //First index of coupling
      int * firstIndexCoupling;
      
      //Second index of coupling
      int * secondIndexCoupling;
      
      //Size of the two-site operator SVD decomposition (d*d)
      int trotterSVDsize;
      
      //Setup the trial wfn
      void SetupTrial();
      
      /******************
      *** The walkers ***
      ******************/
      
      //The walkers
      Walker ** theWalkers;
      
      //The walkers: copy array
      Walker ** theWalkersCopyArray;
      
      //Setup the walkers
      void SetupWalkers();
      
      /***************************************************
      *** For MPI: some helper variables and functions ***
      ***************************************************/
      
      //My MPI rank
      int MPIrank;
      
      //The MPI size
      int MPIsize;
      
      //The max. number of openMP threads on this process
      int * NThreadsPerRank;
      
      //The desired walkers per rank
      int * NDesiredWalkersPerRank;
      
      //The desired walkers per rank
      int * NCurrentWalkersPerRank;
      
      //The max. allowed number of walkers on this process
      int myMaxNWalkers;
      
      //Set the OpenMP and MPI load distribution
      void SetupOMPandMPILoadDistribution();
      
      /************************************************************************************************************
      *** Some variables which are needed at each MC step, so that they don't need to be reallocated every time ***
      ************************************************************************************************************/
      
      //Workspace to temporarily store the PDF in for different possible moves per Trotter term
      double ** thePDF;
      
      //Workspace to temporarily store the different left and right operator combinations per Trotter term in
      double ** theOperatorCombos;
      
      //To keep track of the sum of the walker coeff
      double * sumWalkerWeightPerThread;
      
};

#endif
