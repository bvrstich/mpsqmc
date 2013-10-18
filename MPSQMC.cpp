#include "MPItag.h" //Here USE_MPI_IN_MPSQMC can be defined to switch to non-MPI compilers.

#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include <omp.h>

#include "MPSQMC.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on September 5, 2013 */

using namespace std;

MPSQMC::MPSQMC(MPO * theMPO, Random * RN, const int Dtrunc, const int Nwalkers, const double dtau){

   this->theMPO = theMPO;
   this->RN = RN;
   this->Dtrunc = Dtrunc;
   this->dtau = dtau;
   this->nMPOterms = theMPO->gRN_nTerms();
   
   this->totalNDesiredWalkers = Nwalkers;
   this->totalNCurrentWalkers = Nwalkers;
   
   #ifdef _OPENMP
      myNOpenMPthreads = omp_get_max_threads();
   #else
      myNOpenMPthreads = 1;
   #endif

   #ifdef USE_MPI_IN_MPSQMC
      //MPI info
      MPIsize = MPI::COMM_WORLD.Get_size();
      MPIrank = MPI::COMM_WORLD.Get_rank();
      
      //Finding out how many threads are running for each process, as well as the total amount
      int * threadsPerRank = new int[MPIsize];
      threadsPerRank[MPIrank] = myNOpenMPthreads;
      for (int count=0; count<MPIsize; count++){
         MPI::COMM_WORLD.Bcast(threadsPerRank + count, 1, MPI::INT, count);
      }
      int totalNOpenMPthreads = 0;
      for (int count=0; count<MPIsize; count++){ totalNOpenMPthreads += threadsPerRank[count]; }
      
      //Distribute the workload according to the threads
      myNCurrentWalkers = Nwalkers / totalNOpenMPthreads;
      myNCurrentWalkers *= myNOpenMPthreads;
      
      //Find the cumulative threadcount: t[0] = t_0 ; t[1] = t_0 + t_1 ; t[2] = t_0 + t_1 + t_2 ; ... ; t[n] = sum(i=0..n) t_i
      for (int count=1; count<MPIsize; count++){ threadsPerRank[count] += threadsPerRank[count-1]; }
      
      //Distribute the remainder of the walkers over the first "remainder" threads:
      int remainder = Nwalkers % totalNOpenMPthreads;
      if (MPIrank>0){ remainder -= threadsPerRank[MPIrank-1]; }
      if (remainder > 0){ myNCurrentWalkers += min(remainder, myNOpenMPthreads); }
      
      //Find out the numbers for the other ranks
      NDesiredWalkersPerRank = new int[MPIsize];
      NDesiredWalkersPerRank[MPIrank] = myNCurrentWalkers;
      for (int count=0; count<MPIsize; count++){
         MPI::COMM_WORLD.Bcast(NDesiredWalkersPerRank + count, 1, MPI::INT, count);
      }
      
      if (MPIrank==0){
         cout << "There are " << MPIsize << " MPI processes." << endl;
         for (int cnt=0; cnt<MPIsize; cnt++){
            int nThreads = (cnt==0) ? threadsPerRank[cnt] : (threadsPerRank[cnt] - threadsPerRank[cnt-1]);
            cout << "   MPI process " << cnt << " has " << nThreads << " OMP threads and " << NDesiredWalkersPerRank[cnt] << " desired walkers." << endl;
         }
      }
      delete [] threadsPerRank;
   #else
      MPIsize = 1;
      MPIrank = 0;
      myNCurrentWalkers = Nwalkers;
   #endif
   
   myNDesiredWalkers = myNCurrentWalkers;
   myMaxNWalkers = max(1000,3*myNDesiredWalkers);
   bSetupTrial = false;
   bSetupWalkers = false;
   SetupTrial();
   SetupWalkers();

}

MPSQMC::~MPSQMC(){

   #ifdef USE_MPI_IN_MPSQMC
      delete [] NDesiredWalkersPerRank;
   #endif

   if (bSetupTrial){
      for (int cnt=0; cnt<myNOpenMPthreads; cnt++){
         for (int cnt2=1; cnt2<nMPOterms+1; cnt2++){ delete MPOtermsPsi0[cnt][cnt2]; }
         delete Psi0[cnt];
         delete HPsi0[cnt];
         delete [] MPOtermsPsi0[cnt];
      }
      delete [] Psi0;
      delete [] HPsi0;
      delete [] MPOtermsPsi0;
   }
   if (bSetupWalkers){
      for (int cnt=0; cnt<myNCurrentWalkers; cnt++){ delete theWalkers[cnt]; }
      delete [] theWalkers;
      delete [] theWalkersCopyArray;
      delete [] walkerCoeff;
      delete [] walkerCoeffCopy;
      delete [] walkerOverlap;
      delete [] walkerOverlapCopy;
      delete [] walkerEnergyHistorySum;
      delete [] walkerEnergyHistorySumCopy;
      
      for (int cnt=0; cnt<myNOpenMPthreads; cnt++){ delete [] thePDF[cnt]; }
      delete [] thePDF;
      
      delete [] sumWalkerCoeffPerThread;
   }
   
}

void MPSQMC::SetupTrial(){

   if (!bSetupTrial){
   
      bSetupTrial = true;
      
      //The lowest state for given D
      Psi0 = new MPSstate*[myNOpenMPthreads];
      if (MPIrank==0){
         Psi0[0] = new MPSstate(theMPO->gLength(),Dtrunc,theMPO->gPhys_d(),RN);
         DMRG * solver = new DMRG(Psi0[0],theMPO);
         solver->Solve();
         delete solver;
         Psi0[0]->LeftNormalize();
      }
      #ifdef USE_MPI_IN_MPSQMC
         Psi0[0] = BroadcastCopyConstruct(Psi0[0]);
      #endif
      for (int cnt=1; cnt<myNOpenMPthreads; cnt++){ Psi0[cnt] = new MPSstate(Psi0[0]); }
      
      HPsi0 = new MPSstate*[myNOpenMPthreads];
      if (MPIrank==0){
         HPsi0[0] = new MPSstate(theMPO->gLength(),Dtrunc,theMPO->gPhys_d(),RN);
         HPsi0[0]->ApplyMPO(theMPO, Psi0[0]);
         HPsi0[0]->CompressState(); //Compression only throws away Schmidt values which are numerically zero...
      }
      #ifdef USE_MPI_IN_MPSQMC
         HPsi0[0] = BroadcastCopyConstruct(HPsi0[0]);
      #endif
      for (int cnt=1; cnt<myNOpenMPthreads; cnt++){ HPsi0[cnt] = new MPSstate(HPsi0[0]); }

      MPOtermsPsi0 = new MPSstate**[myNOpenMPthreads];
      for (int cnt=0; cnt<myNOpenMPthreads; cnt++){
         MPOtermsPsi0[cnt] = new MPSstate*[nMPOterms+1];
         MPOtermsPsi0[cnt][0] = Psi0[cnt];
         for (int cnt2=1; cnt2<nMPOterms+1; cnt2++){
            if (cnt==0){
               MPOtermsPsi0[cnt][cnt2] = new MPSstate(Psi0[cnt]);
               MPOtermsPsi0[cnt][cnt2]->ApplyMPOtermHC(theMPO,cnt2-1);
            } else {
               MPOtermsPsi0[cnt][cnt2] = new MPSstate(MPOtermsPsi0[0][cnt2]);
            }
         }
      }
      
   }

}

MPSstate * MPSQMC::BroadcastCopyConstruct(MPSstate * pointer){

   #ifdef USE_MPI_IN_MPSQMC
      //Make the Psi storage exactly the same as rank 0 Psi storage
      int dim = theMPO->gLength()+1;
      int * VirtualDims = new int[dim];
      int truncDim = 0;
      if (MPIrank==0){
         for (int cnt=0; cnt<dim; cnt++){ VirtualDims[cnt] = pointer->gDimAtBound(cnt); }
         truncDim = pointer->gDtrunc();
      }
      MPI::COMM_WORLD.Bcast(VirtualDims, dim, MPI::INT, 0);
      MPI::COMM_WORLD.Bcast(&truncDim, 1, MPI::INT, 0);
      if (MPIrank>0){
         pointer = new MPSstate(theMPO->gLength(), truncDim, theMPO->gPhys_d(), VirtualDims, RN);
      }
      delete [] VirtualDims;
      
      //Broadcast the Psi storage from rank 0
      for (int site=0; site<theMPO->gLength(); site++){
         int dim = pointer->gDimAtBound(site) * pointer->gDimAtBound(site+1) * pointer->gPhys_d();
         double * storage = pointer->gMPStensor(site)->gStorage();
         MPI::COMM_WORLD.Bcast(storage, dim, MPI::DOUBLE, 0);
      }
   #endif
   
   return pointer;

}

void MPSQMC::SetupWalkers(){

   const bool copyTrial = true;

   if (!bSetupWalkers){
   
      bSetupWalkers = true;
      
      theWalkers = new MPSstate * [myMaxNWalkers];
      theWalkersCopyArray = new MPSstate * [myMaxNWalkers];
      walkerCoeff = new double[myMaxNWalkers];
      walkerCoeffCopy = new double[myMaxNWalkers];
      walkerOverlap = new double[myMaxNWalkers];
      walkerOverlapCopy = new double[myMaxNWalkers];
      walkerEnergyHistorySum = new double[myMaxNWalkers];
      walkerEnergyHistorySumCopy = new double[myMaxNWalkers];
      
      thePDF = new double*[myNOpenMPthreads];
      for (int cnt=0; cnt<myNOpenMPthreads; cnt++){ thePDF[cnt] = new double[nMPOterms+1]; }
      
      sumWalkerCoeffPerThread = new double[myNOpenMPthreads];
      
      for (int cnt=0; cnt<myNCurrentWalkers; cnt++){
         if (copyTrial){
            theWalkers[cnt] = new MPSstate(Psi0[0]);
         } else {
            theWalkers[cnt] = new MPSstate(theMPO->gLength(), Dtrunc, theMPO->gPhys_d(), RN);
         }
         theWalkersCopyArray[cnt] = NULL;
         walkerCoeff[cnt] = 1.0;
         walkerCoeffCopy[cnt] = 0.0;
         walkerOverlap[cnt] = 1.0;
         walkerOverlapCopy[cnt] = 0.0;
         walkerEnergyHistorySum[cnt] = 0.0;
         walkerEnergyHistorySumCopy[cnt] = 0.0;
      }
      for (int cnt=myNCurrentWalkers; cnt<myMaxNWalkers; cnt++){
         theWalkers[cnt] = NULL;
         theWalkersCopyArray[cnt] = NULL;
         walkerCoeff[cnt] = 0.0;
         walkerCoeffCopy[cnt] = 0.0;
         walkerOverlap[cnt] = 0.0;
         walkerOverlapCopy[cnt] = 0.0;
         walkerEnergyHistorySum[cnt] = 0.0;
         walkerEnergyHistorySumCopy[cnt] = 0.0;
      }
   
   }

}

void MPSQMC::Walk(const int steps){

   double projectedEnergy = EnergyFunction();
   if (MPIrank==0){
      ofstream output("energies.txt",ios::trunc);
      output << "#Step\t\tE_P\t\tE_T\t\tFluctMetric" << endl;
      output.close();
      
      cout << "Energy at start = " << projectedEnergy << endl;
      cout << "---------------------------------------------------------" << endl;
   }

   for (int step=1; step<=steps; step++){

      //Propagate the walkers of each rank separately --> no MPI in that function
      double mySumOfWalkerCoeff = PropagateSeparately();
      
      //Form the total sum of the walker coefficients
      #ifdef USE_MPI_IN_MPSQMC
         double totalSum = 0.0;
         MPI::COMM_WORLD.Allreduce(&mySumOfWalkerCoeff, &totalSum, 1, MPI::DOUBLE, MPI::SUM);
         double avgCoeff = totalSum / totalNCurrentWalkers;
      #else
         double avgCoeff = mySumOfWalkerCoeff / totalNCurrentWalkers;
      #endif
      
      //Update the energy history for each walker and return the fluctuation metric
      double fluctMetric = updateEnergyHistory(step);
      
      //Calculate the energy --> uses MPI
      projectedEnergy = EnergyFunction();
      
      //Calculate the scaling for population control
      double scaling = totalNDesiredWalkers / (totalNCurrentWalkers * avgCoeff);
      double targetEnergy = log(scaling)/dtau;
      
      if (MPIrank==0){
         cout << "        Step = " << step << endl;
         cout << "   # walkers = " << totalNCurrentWalkers << endl;
         cout << " avg(weight) = " << avgCoeff << endl;
         cout << "         E_P = " << projectedEnergy << endl;
         cout << "         E_T = " << targetEnergy << endl;
         cout << "       Omega = " << fluctMetric << endl;
         write(step, projectedEnergy, targetEnergy, fluctMetric);
         cout << "---------------------------------------------------------" << endl;
      }
      
      //Based on scaling, control the population on each rank separately --> no MPI in that function
      SeparatePopulationControl(scaling);
      
      //Population balancing --> uses MPI of course
      #ifdef USE_MPI_IN_MPSQMC
         PopulationBalancing();
      #endif
      
      //Share all global variables that may have changed
      #ifdef USE_MPI_IN_MPSQMC
         totalNCurrentWalkers = 0;
         MPI::COMM_WORLD.Allreduce(&myNCurrentWalkers, &totalNCurrentWalkers, 1, MPI::INT, MPI::SUM);
      #else
         totalNCurrentWalkers = myNCurrentWalkers;
      #endif

   }

}

double MPSQMC::PropagateSeparately(){

   for (int cnt=0; cnt<myNOpenMPthreads; cnt++){ sumWalkerCoeffPerThread[cnt] = 0.0; }
   
   #pragma omp parallel for schedule(static) default(none)
   for (int walker=0; walker<myNCurrentWalkers; walker++){
      
      #ifdef _OPENMP
         int myID = omp_get_thread_num();
      #else
         int myID = 0;
      #endif
      
      //Get PDF for IS (importance sampling)
      thePDF[myID][0] = 1.0;
      double NormOfPDF = thePDF[myID][0];
      for (int cnt=1; cnt<nMPOterms+1; cnt++){
         thePDF[myID][cnt] = max( -dtau * MPOtermsPsi0[myID][cnt]->InnerProduct(theWalkers[walker]) / walkerOverlap[walker] , 0.0 );
         NormOfPDF += thePDF[myID][cnt];
      }
      for (int cnt=0; cnt<nMPOterms+1; cnt++){ thePDF[myID][cnt] /= NormOfPDF; }
      
      //Get what you should do
      double randomNumber = RN->rand();
      int count = 0;
      double CumPDF = thePDF[myID][0];
      while ((CumPDF <= randomNumber) && (count<nMPOterms)){
         count += 1;
         CumPDF += thePDF[myID][count];
      }
      
      //Do what you should do
      if (count==0){
         //Walker norms only contribute to the overlap, not to the weights --> See PRB 55, 7464 (1997): end of paragraph IV-C.
         //Hence do nothing, walker and walker norm do not evolve.
      } else {
         theWalkers[walker]->ApplyMPOterm(theMPO,count-1);
         double val = theWalkers[walker]->LeftNormalize();
         val = MPOtermsPsi0[myID][0]->InnerProduct(theWalkers[walker]);
         if (val<0.0){ //Due to IS, only overlaps>0 are selected. Because norms are arbitrary (see comment in if statement), the factor -dtau can be left out too.
            theWalkers[walker]->ChangePhase();
            val *= -1.0;
         }
         walkerOverlap[walker] = val;
      }
      walkerCoeff[walker] *= NormOfPDF; //Weights are automatically positive
      sumWalkerCoeffPerThread[myID] += walkerCoeff[walker];
   }
   
   int sumOfWalkerCoeff = 0;
   for (int cnt=0; cnt<myNOpenMPthreads; cnt++){ sumOfWalkerCoeff += sumWalkerCoeffPerThread[cnt]; }
   return sumOfWalkerCoeff;

}

double MPSQMC::updateEnergyHistory(const int step){

   double sumOfIndividualAvgEnergiesSquared = 0.0; //Per MPI rank
   double sumOfIndividialAvgEnergies = 0.0; //Per MPI rank
   
   for (int walker=0; walker<myNCurrentWalkers; walker++){
      double myCurrentEnergy = HPsi0[0]->InnerProduct(theWalkers[walker]) / walkerOverlap[walker];
      walkerEnergyHistorySum[walker] += myCurrentEnergy; //EnergyHistory update
      double myAvgEnergy = walkerEnergyHistorySum[walker] / step; // e_j(t) : individual time average energy
      sumOfIndividialAvgEnergies += myAvgEnergy;
      sumOfIndividualAvgEnergiesSquared += myAvgEnergy * myAvgEnergy;
   }
   
   #ifdef USE_MPI_IN_MPSQMC
      double allE = 0.0;
      double allEsquared = 0.0;
      MPI::COMM_WORLD.Allreduce(&sumOfIndividialAvgEnergies,        &allE,        1, MPI::DOUBLE, MPI::SUM);
      MPI::COMM_WORLD.Allreduce(&sumOfIndividualAvgEnergiesSquared, &allEsquared, 1, MPI::DOUBLE, MPI::SUM);
   #else
      double allE = sumOfIndividialAvgEnergies;
      double allEsquared = sumOfIndividualAvgEnergiesSquared;
   #endif
   
   double fluctMetric = ( allEsquared - allE * allE / totalNCurrentWalkers ) / totalNCurrentWalkers;
   return fluctMetric;

}

void MPSQMC::SeparatePopulationControl(const double scaling){
   
   int newNumberOfWalkers = 0;
   for (int walker=0; walker<myNCurrentWalkers; walker++){
   
      int nCopies = (int) (walkerCoeff[walker] * scaling + RN->rand());
      double newCoeff = 1.0;
      
      if (newNumberOfWalkers + nCopies > myMaxNWalkers){
         cout << "MPSQMC::Walk (MPI rank " << MPIrank << ") -> The number of desired walkers exceeds the max. allowed number at MPI." << endl;
      }
      
      if (nCopies==0){ delete theWalkers[walker]; }
      if (nCopies>=1){
         theWalkersCopyArray[newNumberOfWalkers] = theWalkers[walker];
         walkerOverlapCopy[newNumberOfWalkers] = walkerOverlap[walker];
         walkerCoeffCopy[newNumberOfWalkers] = newCoeff;
         walkerEnergyHistorySumCopy[newNumberOfWalkers] = walkerEnergyHistorySum[walker];
         for (int cnt=1; cnt<nCopies; cnt++){
            theWalkersCopyArray[newNumberOfWalkers+cnt] = new MPSstate(theWalkers[walker]);
            walkerOverlapCopy[newNumberOfWalkers+cnt] = walkerOverlap[walker];
            walkerCoeffCopy[newNumberOfWalkers+cnt] = newCoeff;
            walkerEnergyHistorySumCopy[newNumberOfWalkers+cnt] = walkerEnergyHistorySum[walker];
         }
         newNumberOfWalkers += nCopies;
      }
   }
   myNCurrentWalkers = newNumberOfWalkers;
   //Swap the walker arrays
   MPSstate ** tempMPSstate = theWalkers;
   theWalkers = theWalkersCopyArray;
   theWalkersCopyArray = tempMPSstate;
   //Swap the overlap arrays
   double * tempOverlap = walkerOverlap;
   walkerOverlap = walkerOverlapCopy;
   walkerOverlapCopy = tempOverlap;
   //Swap the coeff arrays
   double * tempCoeff = walkerCoeff;
   walkerCoeff = walkerCoeffCopy;
   walkerCoeffCopy = tempCoeff;
   //Swap the energy history arrays
   double * tempEnergyHistory = walkerEnergyHistorySum;
   walkerEnergyHistorySum = walkerEnergyHistorySumCopy;
   walkerEnergyHistorySumCopy = tempEnergyHistory;

}

double MPSQMC::EnergyFunction(){

   double myNumerator = 0.0;
   double myDenominator = 0.0;
   for (int walker=0; walker<myNCurrentWalkers; walker++){
      myNumerator += walkerCoeff[walker] * HPsi0[0]->InnerProduct(theWalkers[walker]) / walkerOverlap[walker];
      myDenominator += walkerCoeff[walker];
   }
   #ifdef USE_MPI_IN_MPSQMC
      double totalNumerator = 0.0;
      double totalDenominator = 0.0;
      MPI::COMM_WORLD.Allreduce(&myNumerator,   &totalNumerator,   1, MPI::DOUBLE, MPI::SUM);
      MPI::COMM_WORLD.Allreduce(&myDenominator, &totalDenominator, 1, MPI::DOUBLE, MPI::SUM);
      double Energy = totalNumerator / totalDenominator;
   #else
      double Energy = myNumerator / myDenominator;
   #endif

   return Energy;

}

void MPSQMC::write(const int step, const double projectedEnergy, const double targetEnergy, const double fluctMetric){

   ofstream output("energies.txt",ios::app);
   output.precision(10);
   output << step << "\t\t" << projectedEnergy << "\t\t" << targetEnergy << "\t\t" << fluctMetric << endl;
   output.close();

}

void MPSQMC::BubbleSort(double * values, int * order, const int length){

   for (int cnt=0; cnt<length; cnt++){ order[cnt] = cnt; }
   bool allOK = false;
   while (!allOK){
      allOK = true;
      for (int cnt=0; cnt<length-1; cnt++){
         if ( values[ order[cnt] ] < values[ order[cnt+1] ] ){
            allOK = false;
            int temp = order[cnt];
            order[cnt] = order[cnt+1];
            order[cnt+1] = temp;
         }
      }
   } // Now for all index: values[ order[index] ] >= values[ order[index+1] ]

}

void MPSQMC::PopulationBalancing(){

   #ifdef USE_MPI_IN_MPSQMC
      const bool debuginfoprint = true;
      const double threshold_start = 0.1;
      const double threshold_stop  = 0.03;
      const bool clusterhasinfiniband = true;

      int * Ncurr = new int[MPIsize];
      double * Noffset = new double[MPIsize];
      bool oneFracDeviating = false;
      
      Ncurr[MPIrank] = myNCurrentWalkers;
      int Ntotal = 0;
      for (int count=0; count<MPIsize; count++){
         MPI::COMM_WORLD.Bcast(&Ncurr[count], 1, MPI::INT, count); //Now we know everybody's current population
         Ntotal += Ncurr[count]; //and the total population
      }
      
      const double fractionScaling = ((double) Ntotal) / totalNDesiredWalkers; //We proportionally want to distribute the total load
      for (int count=0; count<MPIsize; count++){
         Noffset[count] = Ncurr[count] - NDesiredWalkersPerRank[count] * fractionScaling; //So we calculate the offset from "equilibrium"
         double frac = fabs( Noffset[count] ) / ( NDesiredWalkersPerRank[count] * fractionScaling ); //and the percentage of offset
         if (frac > threshold_start){ oneFracDeviating = true; } //if one or more offsets are too large: redistribute       
      }
      
      if ((debuginfoprint) && (MPIrank==0)){
         for (int cnt=0; cnt<MPIsize; cnt++){
            cout << "MPSQMC::PopulationBalancing -> Nwalkers(" << cnt << ") = " << Ncurr[cnt] << " and Noffset(" << cnt << ") = " << Noffset[cnt] << endl;
         }
      }
      
      if (oneFracDeviating){
         
         if (debuginfoprint){
            double projectedEnergy = EnergyFunction();
            if (MPIrank==0){ cout << "MPSQMC::PopulationBalancing -> As a check: projected E before = " << projectedEnergy << endl; }
         }
         
         int * work = new int[MPIsize];
         int communication_round = 0;
      
         while (oneFracDeviating){
         
            communication_round++;
         
            //Do a bubble sort from large to small (positive to negative). Not optimal algo, optimal = quicksort (dlasrt_) --> but multi-array sort in lapack?
            BubbleSort(Noffset, work, MPIsize); // Now for all index: Noffset [ work[index] ] >= Noffset[ work[index+1] ]
            
            //Communicate walkers between rank work[comm] and rank work[MPIsize - 1 - comm] until "drained"
            for (int comm=0; comm< ((clusterhasinfiniband) ? MPIsize/2 : 1); comm++){
            
               const int sender = work[comm];
               const int receiver = work[MPIsize - 1 - comm];
               
               if ((Noffset[sender] > 0.0) && (Noffset[receiver] < 0.0)){
               
                  int amount = (int) min(Noffset[sender], -Noffset[receiver]); //This explains the "until drained" statement.
                  if (amount>0){
                  
                     if ((debuginfoprint) && (MPIrank==0)){
                        cout << "MPSQMC::PopulationBalancing -> Moving " << amount << " walkers from rank " << sender << " to rank " << receiver << " during communication round " << communication_round << "." << endl;
                     }
                  
                     if (MPIrank == sender){
                        for (int counter=0; counter<amount; counter++){
                           double Coeff = walkerCoeff[myNCurrentWalkers-1];
                           double Overlap = walkerOverlap[myNCurrentWalkers-1];
                           MPI::COMM_WORLD.Send(&Coeff, 1, MPI::DOUBLE, receiver, 0);
                           MPI::COMM_WORLD.Send(&Overlap, 1, MPI::DOUBLE, receiver, 0);
                           for (int site=0; site<theMPO->gLength(); site++){
                              int dim = Psi0[0]->gPhys_d() * Psi0[0]->gDimAtBound(site) * Psi0[0]->gDimAtBound(site+1);
                              double * SendStorage = theWalkers[myNCurrentWalkers-1]->gMPStensor(site)->gStorage();
                              MPI::COMM_WORLD.Send(SendStorage, dim, MPI::DOUBLE, receiver, 0);
                           }
                           delete theWalkers[myNCurrentWalkers-1];
                           myNCurrentWalkers--;
                        }
                     }
                     
                     if (MPIrank == receiver){
                        for (int counter=0; counter<amount; counter++){
                           double Coeff = 0.0;
                           double Overlap = 0.0;
                           MPI::Status status;
                           MPI::COMM_WORLD.Recv(&Coeff, 1, MPI::DOUBLE, sender, 0, status);
                           MPI::COMM_WORLD.Recv(&Overlap, 1, MPI::DOUBLE, sender, 0, status);
                           walkerCoeff[myNCurrentWalkers] = Coeff;
                           walkerOverlap[myNCurrentWalkers] = Overlap;
                           theWalkers[myNCurrentWalkers] = new MPSstate(Psi0[0]);
                           for (int site=0; site<theMPO->gLength(); site++){
                              int dim = Psi0[0]->gPhys_d() * Psi0[0]->gDimAtBound(site) * Psi0[0]->gDimAtBound(site+1);
                              double * RecvStorage = theWalkers[myNCurrentWalkers]->gMPStensor(site)->gStorage();
                              MPI::COMM_WORLD.Recv(RecvStorage, dim, MPI::DOUBLE, sender, 0, status);
                           }
                           myNCurrentWalkers++;
                        }
                     }
                     
                     Ncurr[sender] -= amount;
                     Ncurr[receiver] += amount;
                  }
                  
               }
            
            } //Everything got communicated in parallel
            
            MPI::COMM_WORLD.Barrier();
            
            //Determine the offsets again : now with threshold_stop!!!!
            oneFracDeviating = false;
            for (int count=0; count<MPIsize; count++){
               Noffset[count] = Ncurr[count] - NDesiredWalkersPerRank[count] * fractionScaling;
               double frac = fabs( Noffset[count] ) / ( NDesiredWalkersPerRank[count] * fractionScaling );
               if (frac > threshold_stop){ oneFracDeviating = true; }
            }
         
         }
         
         delete [] work;
         
         if (debuginfoprint){
            double projectedEnergy = EnergyFunction();
            if (MPIrank==0){ cout << "MPSQMC::PopulationBalancing -> As a check: projected E after  = " << projectedEnergy << endl; }
         }
      
      }

      delete [] Noffset;
      delete [] Ncurr;
   #endif

}

/*void MPSQMC::PopulationBalancing(){

   #ifdef USE_MPI_IN_MPSQMC
      const bool debuginfoprint = true;
      const double threshold = 0.1;
      
      if (debuginfoprint){
         double projectedEnergy = EnergyFunction();
         if (MPIrank==0){ cout << "MPSQMC::PopulationBalancing -> As a check: projected E before = " << projectedEnergy << endl; }
      }

      int * Ncurr = new int[MPIsize];
      double * Noffset = new double[MPIsize];
      bool oneFracDeviating = false;
      
      Ncurr[MPIrank] = myNCurrentWalkers;
      int Ntotal = 0;
      for (int count=0; count<MPIsize; count++){
         MPI::COMM_WORLD.Bcast(&Ncurr[count], 1, MPI::INT, count); //Now we know everybody's current population
         Ntotal += Ncurr[count]; //and the total population
      }
      
      const double fractionScaling = ((double) Ntotal) / totalNDesiredWalkers; //We proportionally want to distribute the total load
      for (int count=0; count<MPIsize; count++){
         Noffset[count] = Ncurr[count] - NDesiredWalkersPerRank[count] * fractionScaling; //So we calculate the offset from "equilibrium"
         double frac = fabs( Noffset[count] ) / ( NDesiredWalkersPerRank[count] * fractionScaling ); //and the percentage of offset
         if (frac > threshold){ oneFracDeviating = true; } //if one or more offsets are too large: redistribute       
      }
      
      if ((debuginfoprint) && (MPIrank==0)){
         for (int cnt=0; cnt<MPIsize; cnt++){
            cout << "MPSQMC::PopulationBalancing -> Nwalkers(" << cnt << ") = " << Ncurr[cnt] << " and Noffset(" << cnt << ") = " << Noffset[cnt] << endl;
         }
      }
      
      while (oneFracDeviating){
      
         //Do a bubble sort from large to small (positive to negative). Not optimal algo, optimal = quicksort (dlasrt_) --> but multi-array sort in lapack?
         int * work = new int[MPIsize];
         for (int cnt=0; cnt<MPIsize; cnt++){ work[cnt] = cnt; }
         bool allOK = false;
         while (!allOK){
            allOK = true;
            for (int cnt=0; cnt<MPIsize-1; cnt++){
               if ( Noffset[ work[cnt] ] < Noffset[ work[cnt+1] ] ){
                  allOK = false;
                  int temp = work[cnt];
                  work[cnt] = work[cnt+1];
                  work[cnt+1] = temp;
               }
            }
         } // Now for all index: Noffset [ work[index] ] >= Noffset[ work[index+1] ]
         
         //Communicate walkers between rank work[comm] and rank work[MPIsize - 1 - comm] until "drained"
         for (int comm=0; comm< MPIsize/2; comm++){
         
            const int sender = work[comm];
            const int receiver = work[MPIsize - 1 - comm];
            
            if ((Noffset[sender] > 0.0) && (Noffset[receiver] < 0.0)){
            
               int amount = (int) min(Noffset[sender], -Noffset[receiver]); //This explains the "until drained" statement.
               if (amount>0){
               
                  if ((debuginfoprint) && (MPIrank==0)){
                     cout << "MPSQMC::PopulationBalancing -> Moving " << amount << " walkers from rank " << sender << " to rank " << receiver << "." << endl;
                  }
               
                  if (MPIrank == sender){
                     for (int counter=0; counter<amount; counter++){
                        double Coeff = walkerCoeff[myNCurrentWalkers-1];
                        double Overlap = walkerOverlap[myNCurrentWalkers-1];
                        MPI::COMM_WORLD.Send(&Coeff, 1, MPI::DOUBLE, receiver, 0);
                        MPI::COMM_WORLD.Send(&Overlap, 1, MPI::DOUBLE, receiver, 0);
                        for (int site=0; site<theMPO->gLength(); site++){
                           int dim = Psi0[0]->gPhys_d() * Psi0[0]->gDimAtBound(site) * Psi0[0]->gDimAtBound(site+1);
                           double * SendStorage = theWalkers[myNCurrentWalkers-1]->gMPStensor(site)->gStorage();
                           MPI::COMM_WORLD.Send(SendStorage, dim, MPI::DOUBLE, receiver, 0);
                        }
                        delete theWalkers[myNCurrentWalkers-1];
                        myNCurrentWalkers--;
                     }
                  }
                  
                  if (MPIrank == receiver){
                     for (int counter=0; counter<amount; counter++){
                        double Coeff = 0.0;
                        double Overlap = 0.0;
                        MPI::Status status;
                        MPI::COMM_WORLD.Recv(&Coeff, 1, MPI::DOUBLE, sender, 0, status);
                        MPI::COMM_WORLD.Recv(&Overlap, 1, MPI::DOUBLE, sender, 0, status);
                        walkerCoeff[myNCurrentWalkers] = Coeff;
                        walkerOverlap[myNCurrentWalkers] = Overlap;
                        theWalkers[myNCurrentWalkers] = new MPSstate(Psi0[0]);
                        for (int site=0; site<theMPO->gLength(); site++){
                           int dim = Psi0[0]->gPhys_d() * Psi0[0]->gDimAtBound(site) * Psi0[0]->gDimAtBound(site+1);
                           double * RecvStorage = theWalkers[myNCurrentWalkers]->gMPStensor(site)->gStorage();
                           MPI::COMM_WORLD.Recv(RecvStorage, dim, MPI::DOUBLE, sender, 0, status);
                        }
                        myNCurrentWalkers++;
                     }
                  }
                  
                  Ncurr[sender] -= amount;
                  Ncurr[receiver] += amount;
               }
               
            }
         
         } //Everything got communicated in parallel
         
         //Determine the offsets again
         oneFracDeviating = false;
         for (int count=0; count<MPIsize; count++){
            Noffset[count] = Ncurr[count] - NDesiredWalkersPerRank[count] * fractionScaling;
            double frac = fabs( Noffset[count] ) / ( NDesiredWalkersPerRank[count] * fractionScaling );
            if (frac > threshold){ oneFracDeviating = true; }
         }
         
         delete [] work;
      
      }
      
      if (debuginfoprint){
         double projectedEnergy = EnergyFunction();
         if (MPIrank==0){ cout << "MPSQMC::PopulationBalancing -> As a check: projected E after  = " << projectedEnergy << endl; }
      }

      delete [] Noffset;
      delete [] Ncurr;
   #endif

}*/

/*void MPSQMC::PopulationBalancing(){

   //Very very very rudimentary population balancing. Better to set up a schedule what needs to be sent where, and then do it all simultaneously.

   #ifdef USE_MPI_IN_MPSQMC
      const bool debuginfoprint = false;
      
      if (debuginfoprint){
         double projectedEnergy = EnergyFunction();
         if (MPIrank==0){
            cout << "MPSQMC::PopulationBalancing -> As a check: projected E before = " << projectedEnergy << endl;
         }
      }
   
      const double thresholdUp = 1.1;
      const double thresholdLow = 0.9;

      int * Ncurr = new int[MPIsize];
      double * frac = new double[MPIsize];
      bool oneFracLarger = false;
      bool oneFracSmaller = false;
      
      Ncurr[MPIrank] = myNCurrentWalkers;
      for (int count=0; count<MPIsize; count++){
         MPI::COMM_WORLD.Bcast(&Ncurr[count], 1, MPI::INT, count);
         frac[count] = ((double) (Ncurr[count])) / NDesiredWalkersPerRank[count];
         if (frac[count] > thresholdUp){ oneFracLarger = true; }
         if (frac[count] < thresholdLow){ oneFracSmaller = true; }
      }
      
      if (debuginfoprint){
         if (MPIrank==0){
            for (int cnt=0; cnt<MPIsize; cnt++){
               cout << "MPSQMC::PopulationBalancing -> Nwalkers(" << cnt << ") = " << Ncurr[cnt] << endl;
            }
         }
      }
      
      while ((oneFracLarger) || (oneFracSmaller)){
      
         //Find out the smallest and largest fractions
         int largest = 0;
         int smallest = 0;
         for (int count=1; count<MPIsize; count++){
            if (frac[count] > frac[largest]){ largest = count; }
            if (frac[count] < frac[smallest]){ smallest = count; }
         }
         
         if (debuginfoprint){
            if (MPIrank==0){
               cout << "MPSQMC::PopulationBalancing -> Moving a walker from rank " << largest << " to rank " << smallest << "." << endl;
            }
         }
         
         //Send over coeff, overlap, and walker --> using the fact that the virtual dimensions of Psi0 never change during walker propagation
         if (MPIrank == largest){
            double Coeff = walkerCoeff[myNCurrentWalkers-1];
            double Overlap = walkerOverlap[myNCurrentWalkers-1];
            MPI::COMM_WORLD.Send(&Coeff, 1, MPI::DOUBLE, smallest, 0);
            MPI::COMM_WORLD.Send(&Overlap, 1, MPI::DOUBLE, smallest, 0);
            for (int site=0; site<theMPO->gLength(); site++){
               int dim = Psi0[0]->gPhys_d() * Psi0[0]->gDimAtBound(site) * Psi0[0]->gDimAtBound(site+1);
               double * SendStorage = theWalkers[myNCurrentWalkers-1]->gMPStensor(site)->gStorage();
               MPI::COMM_WORLD.Send(SendStorage, dim, MPI::DOUBLE, smallest, 0);
            }
            delete theWalkers[myNCurrentWalkers-1];
            myNCurrentWalkers--;
         }
         if (MPIrank == smallest){
            double Coeff = 0.0;
            double Overlap = 0.0;
            MPI::Status status;
            MPI::COMM_WORLD.Recv(&Coeff, 1, MPI::DOUBLE, largest, 0, status);
            MPI::COMM_WORLD.Recv(&Overlap, 1, MPI::DOUBLE, largest, 0, status);
            walkerCoeff[myNCurrentWalkers] = Coeff;
            walkerOverlap[myNCurrentWalkers] = Overlap;
            theWalkers[myNCurrentWalkers] = new MPSstate(Psi0[0]);
            for (int site=0; site<theMPO->gLength(); site++){
               int dim = Psi0[0]->gPhys_d() * Psi0[0]->gDimAtBound(site) * Psi0[0]->gDimAtBound(site+1);
               double * RecvStorage = theWalkers[myNCurrentWalkers]->gMPStensor(site)->gStorage();
               MPI::COMM_WORLD.Recv(RecvStorage, dim, MPI::DOUBLE, largest, 0, status);
            }
            myNCurrentWalkers++;
         }
         
         MPI::COMM_WORLD.Barrier();
         
         //Recalculate the fractions
         Ncurr[largest] -= 1;
         Ncurr[smallest] += 1;
         oneFracLarger = false;
         oneFracSmaller = false;
         for (int count=0; count<MPIsize; count++){
            frac[count] = ((double) (Ncurr[count])) / NDesiredWalkersPerRank[count];
            if (frac[count] > thresholdUp){ oneFracLarger = true; }
            if (frac[count] < thresholdLow){ oneFracSmaller = true; }
         }
         
      }
      
      if (debuginfoprint){
         double projectedEnergy = EnergyFunction();
         if (MPIrank==0){
            cout << "MPSQMC::PopulationBalancing -> As a check: projected E after  = " << projectedEnergy << endl;
         }
      }

      delete [] frac;
      delete [] Ncurr;
   #endif

}*/

