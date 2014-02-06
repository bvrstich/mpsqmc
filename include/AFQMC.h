#ifndef AFQMC_H
#define AFQMC_H

#include "MPSstate.h"
#include "TrotterJ1J2.h"
#include "Walker.h"
#include "J1J2MPO.h"
#include <vector>

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 15, 2013 */

class AFQMC{

   public:
   
      //constructor with input trialwavefunction
      AFQMC(J1J2MPO *theMPO,TrotterJ1J2 *theTrotter, Random * RN, MPSstate *Psi0_in,const int DW, const int Nwalkers, const double dtau);
      
      //Destructor
      ~AFQMC();
      
      //Let the walkers propagate for steps steps
      void Walk(const int steps);

      static void init(int,int,int,Random *);

      static void clear();

   private:
      
      //!backup storage for the 'error control'
      static MPSstate **stor;
   
      /******************************
      *** Constructor information ***
      ******************************/
   
      //The Trotter decomposition of the MPO
      TrotterJ1J2 *theTrotter;

      //The J1J2 MPO
      J1J2MPO *theMPO;

      //number of trotter terms
      int n_trot;
      
      //The random number generator
      Random * RN;

      //physical dimension of the sites
      int phys_d;
      
      //The MPS truncation dimension of the trial wavefunction
      int DT;

      //The MPS truncation dimension of the walkers
      int DW;
      
      //The total desired number of walkers
      int totalNDesiredWalkers;
      
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
      complex<double> gEP();
      
      //BubbleSort algorithm. Modifies order so that for all index: values[ order[index] ] >= values[ order[index+1] ].
      void BubbleSort(double * values, int * order, const int length);
      
      //Write the projected energy, target energy, and fluctuation metric at MC time "step"
      void write(const int step,const int nwalkers, const double projectedEnergy, const double targetEnergy);
      
      /****************************
      *** Trial wfn information ***
      ****************************/
      
      //Trial wfn (one per thread)
      MPSstate * Psi0;

      //Setup the trial wfn
      void SetupTrial();
      
      /******************
      *** The walkers ***
      ******************/
      
      //The walkers
      std::vector<Walker*> theWalkers;
      
      //Setup the walkers
      void SetupWalkers(bool);
      
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
      
};

#endif
