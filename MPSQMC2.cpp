#include "MPItag.h" //Here USE_MPI_IN_MPSQMC can be defined to switch to non-MPI compilers.

#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include <omp.h>

#include "MPSQMC2.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 15, 2013 */

using namespace std;

/**
 * constructor of the MPSQMC2 object, takes input parameters that define the QMC walk.
 * @param theMPO MPO containing relevant matrix elements, defines system being studied
 * @param theGrid Generates a grid for the random picking of terms from the two-site trotter MPO's
 * @param Dtrunc dimension of the trialstate
 * @param Dtrunc dimension of the walkers
 * @param Nwalkers number of Walker states
 * @param dtau time step of each evolution
 */
MPSQMC2::MPSQMC2(HeisenbergMPO * theMPO, GridGenerator * theGrid, Random * RN, const int DT,const int DW, const int Nwalkers, const double dtau){

   this->theMPO = theMPO;
   this->theGrid = theGrid;
   this->RN = RN;
   this->DT = DT;
   this->DW = DW;
   this->dtau = dtau;

   this->theTrotter = new TrotterHeisenberg(theMPO,dtau);

   this->totalNDesiredWalkers = Nwalkers;
   this->totalNCurrentWalkers = Nwalkers;

   SetupOMPandMPILoadDistribution();

   myMaxNWalkers = max(1000,3*NDesiredWalkersPerRank[MPIrank]);

   MPI::COMM_WORLD.Barrier();

   SetupTrial();
   SetupWalkers();

}

void MPSQMC2::SetupOMPandMPILoadDistribution(){

   //Find the maximum number of threads for this process
#ifdef _OPENMP
   const int myNOMPthreads = omp_get_max_threads();
#else
   const int myNOMPthreads = 1;
#endif

   //Find out about the COMM WORLD, and set the walker distribution over the processes
#ifdef USE_MPI_IN_MPSQMC
   MPIsize = MPI::COMM_WORLD.Get_size();
   MPIrank = MPI::COMM_WORLD.Get_rank();

   //Finding out how many threads are running for each process, as well as the total amount
   NThreadsPerRank          = new int[MPIsize];
   NThreadsPerRank[MPIrank] = myNOMPthreads;

   int totalNOMPthreads = 0;

   for (int count=0; count<MPIsize; count++){

      MPI::COMM_WORLD.Bcast(NThreadsPerRank + count, 1, MPI::INT, count);
      totalNOMPthreads += NThreadsPerRank[count];

   }

   //Distribute the workload according to the threads
   NDesiredWalkersPerRank = new int[MPIsize];

   int remainder = totalNDesiredWalkers;

   for(int count=0; count<MPIsize; count++){

      NDesiredWalkersPerRank[count] = (int)(((double) totalNDesiredWalkers * NThreadsPerRank[count]) / totalNOMPthreads);
      remainder -= NDesiredWalkersPerRank[count];

   }

   for (int count=0; count<remainder; count++)
      NDesiredWalkersPerRank[count] += 1;

   NCurrentWalkersPerRank = new int[MPIsize];

   for (int count=0; count<MPIsize; count++)
      NCurrentWalkersPerRank[count] = NDesiredWalkersPerRank[count];

   if (MPIrank==0){

      cout << "There are " << MPIsize << " MPI processes." << endl;

      for (int cnt=0; cnt<MPIsize; cnt++)
         cout << "   MPI rank " << cnt << " has " << NThreadsPerRank[cnt] << " threads and carries " << NDesiredWalkersPerRank[cnt] << " walkers." << endl;
      
   }

#else
   MPIsize = 1;
   MPIrank = 0;
   NThreadsPerRank                 = new int[MPIsize];
   NDesiredWalkersPerRank          = new int[MPIsize];
   NCurrentWalkersPerRank          = new int[MPIsize];
   NThreadsPerRank[MPIrank]        = myNOMPthreads;
   NDesiredWalkersPerRank[MPIrank] = totalNDesiredWalkers;
   NCurrentWalkersPerRank[MPIrank] = NDesiredWalkersPerRank[MPIrank];
#endif

}

MPSQMC2::~MPSQMC2(){

   //MPSQMC2::MPSQMC2
   delete theTrotter;

   //MPSQMC2::SetupTrial
   for (int cnt=0; cnt<NThreadsPerRank[MPIrank]; cnt++){
      for (int cnt2=0; cnt2<nCouplings; cnt2++){
         for (int cnt3=0; cnt3<trotterSVDsize*trotterSVDsize; cnt3++){ delete TrotterTermsTimesPsi0[cnt][cnt2][cnt3]; }
         delete [] TrotterTermsTimesPsi0[cnt][cnt2];
      }
      delete Psi0[cnt];
      delete HPsi0[cnt];
      delete [] TrotterTermsTimesPsi0[cnt];
   }
   delete [] Psi0;
   delete [] HPsi0;
   delete [] TrotterTermsTimesPsi0;
   delete [] firstIndexCoupling;
   delete [] secondIndexCoupling;

   //MPSQMC2::SetupWalkers
   for (int cnt=0; cnt<NCurrentWalkersPerRank[MPIrank]; cnt++){ delete theWalkers[cnt]; }
   delete [] theWalkers;
   delete [] theWalkersCopyArray;
   for (int cnt=0; cnt<NThreadsPerRank[MPIrank]; cnt++){
      delete [] thePDF[cnt];
      delete [] theOperatorCombos[cnt];
   }
   delete [] thePDF;
   delete [] theOperatorCombos;
   delete [] sumWalkerWeightPerThread;

   //MPSQMC2::SetupOMPandMPILoadDistribution
   delete [] NThreadsPerRank;
   delete [] NDesiredWalkersPerRank;
   delete [] NCurrentWalkersPerRank;

}

/**
 * construct the trial wavefunction by performing a DMRG calculation
 */
void MPSQMC2::SetupTrial(){

   //Rank 0 does the DMRG optimization, and the solution gets copied so every thread on every rank has 1 copy.
   Psi0 = new MPSstate * [NThreadsPerRank[MPIrank]];

   if (MPIrank==0){

      //create initial MPS guess
      Psi0[0] = new MPSstate(theMPO->gLength(),DT,theMPO->gPhys_d(),RN);

      //find optimal MPS for specific MPO using DMRG algorithm
      DMRG * solver = new DMRG(Psi0[0],theMPO);

      solver->Solve();

      delete solver;

      Psi0[0]->LeftNormalize();

   }

#ifdef USE_MPI_IN_MPSQMC
   Psi0[0] = BroadcastCopyConstruct(Psi0[0]);
#endif

   for(int cnt = 1;cnt < NThreadsPerRank[MPIrank];cnt++)
      Psi0[cnt] = new MPSstate(Psi0[0]);

   //Rank 0 calculates MPO times trial, and the result gets copied so every thread on every rank has 1 copy.
   HPsi0 = new MPSstate * [NThreadsPerRank[MPIrank]];

   if(MPIrank==0){

      HPsi0[0] = new MPSstate(theMPO->gLength(),DT,theMPO->gPhys_d(),RN);
      HPsi0[0]->ApplyMPO(theMPO, Psi0[0]);
      HPsi0[0]->CompressState(); //Compression only throws away Schmidt values which are numerically zero...

   }

#ifdef USE_MPI_IN_MPSQMC
   HPsi0[0] = BroadcastCopyConstruct(HPsi0[0]);
#endif

   for(int cnt = 1;cnt < NThreadsPerRank[MPIrank];cnt++)
      HPsi0[cnt] = new MPSstate(HPsi0[0]);

   //Find the number and the indices of the non-zero couplings, as well as the Trotter SVD size
   nCouplings = 0;

   for (int first=0; first<theMPO->gLength()-1; first++)
      for (int second=first+1; second<theMPO->gLength(); second++){

         double coupling = theTrotter->gCoupling(first,second);

         if (coupling != 0.0)
            nCouplings++;

      }

   firstIndexCoupling  = new int [nCouplings];
   secondIndexCoupling = new int [nCouplings];

   nCouplings = 0;

   for(int first = 0;first < theMPO->gLength() - 1;first++)
      for(int second = first + 1;second < theMPO->gLength();second++){

         double coupling = theTrotter->gCoupling(first,second);

         if(coupling != 0.0){

            firstIndexCoupling[nCouplings]  = first;
            secondIndexCoupling[nCouplings] = second;
            nCouplings++;

         }

      }

   trotterSVDsize = theMPO->gPhys_d() * theMPO->gPhys_d();

   //Multiply for all possible combinations of the SVD decomposition of the two-site Trotter terms, its hermitian conjugate into the trial
   TrotterTermsTimesPsi0 = new MPSstate***[ NThreadsPerRank[MPIrank] ];

   for(int cnt = 0;cnt < NThreadsPerRank[MPIrank];cnt++){

      TrotterTermsTimesPsi0[cnt] = new MPSstate ** [nCouplings];

      for(int cnt2 = 0;cnt2 < nCouplings;cnt2++){//loop over the different couplings <ij>

         TrotterTermsTimesPsi0[cnt][cnt2] = new MPSstate*[trotterSVDsize*trotterSVDsize];//d^2 x d^2 different MPS's

         for (int cnt3 = 0;cnt3 < trotterSVDsize;cnt3++)
            for (int cnt4 = 0;cnt4 < trotterSVDsize;cnt4++){

               if (cnt==0){//apply

                  TrotterTermsTimesPsi0[cnt][cnt2][cnt3 + trotterSVDsize * cnt4] = new MPSstate(Psi0[cnt]);
                  TrotterTermsTimesPsi0[cnt][cnt2][cnt3 + trotterSVDsize * cnt4]->ApplyTwoSiteTrotterTerm(theTrotter, firstIndexCoupling[cnt2], secondIndexCoupling[cnt2], cnt3, cnt4, true); //true for HC!!

               }
               else//copy
                  TrotterTermsTimesPsi0[cnt][cnt2][cnt3 + trotterSVDsize * cnt4] = new MPSstate(TrotterTermsTimesPsi0[0][cnt2][cnt3 + trotterSVDsize * cnt4]);

            }

      }

   }

}

/** 
 * different function
 */
MPSstate * MPSQMC2::BroadcastCopyConstruct(MPSstate * pointer){

#ifdef USE_MPI_IN_MPSQMC
   //Make the Psi storage exactly the same as rank 0 Psi storage
   int dim = theMPO->gLength()+1;
   int * VirtualDims = new int[dim];
   int truncDim = 0;
   
   if (MPIrank==0){

      for(int cnt = 0;cnt < dim;cnt++)
         VirtualDims[cnt] = pointer->gDimAtBound(cnt);

      truncDim = pointer->gDtrunc();

   }

   MPI::COMM_WORLD.Bcast(VirtualDims, dim, MPI::INT, 0);
   MPI::COMM_WORLD.Bcast(&truncDim, 1, MPI::INT, 0);

   if(MPIrank>0)
      pointer = new MPSstate(theMPO->gLength(), truncDim, theMPO->gPhys_d(), VirtualDims, RN);
   
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

/**
 * initialize the walkers
 */
void MPSQMC2::SetupWalkers(){

   const bool copyTrial = false;

   theWalkers          = new Walker * [myMaxNWalkers];
   theWalkersCopyArray = new Walker * [myMaxNWalkers];

   thePDF            = new double * [NThreadsPerRank[MPIrank]];
   theOperatorCombos = new double * [NThreadsPerRank[MPIrank]];

   for(int cnt = 0;cnt < NThreadsPerRank[MPIrank];cnt++){

      thePDF[cnt]            = new double [theGrid->gNPoints()];
      theOperatorCombos[cnt] = new double [trotterSVDsize * trotterSVDsize];

   }

   sumWalkerWeightPerThread = new double [NThreadsPerRank[MPIrank]];

   for (int cnt = 0;cnt < NCurrentWalkersPerRank[MPIrank];cnt++){

      if (copyTrial)
         theWalkers[cnt] = new Walker(Psi0[0], 1.0, 1.0, 0.0);
      else{

         theWalkers[cnt] = new Walker(theMPO->gLength(), DW, theMPO->gPhys_d(), RN);
         theWalkers[cnt]->setOverlap(Psi0[0]);

      }

      theWalkersCopyArray[cnt] = NULL;

   }

   for(int cnt=NCurrentWalkersPerRank[MPIrank]; cnt<myMaxNWalkers;cnt++){

      theWalkers[cnt] = NULL;
      theWalkersCopyArray[cnt] = NULL;

   }

}

void MPSQMC2::Walk(const int steps){

   double projectedEnergy = 0.0;

   EnergyFunctionAndHistory(1, &projectedEnergy, false);

   if (MPIrank==0){

      ofstream output("energies.txt",ios::trunc);
      output << "#Step\t\tE_P\t\tE_T\t\tFluctMetric" << endl;
      output.close();

      cout << "Energy at start = " << projectedEnergy << endl;
      cout << "---------------------------------------------------------" << endl;

   }

   if(trotterSVDsize != theGrid->gDim())
      cerr << "trotterSVDsize = " << trotterSVDsize << " and is not equal to theGrid->gDim() = " << theGrid->gDim() << endl;

   for (int step=1; step<=steps; step++){

      //Propagate the walkers of each rank separately --> no MPI in that function
      double mySumOfWalkerWeights = PropagateSeparately();

      //Form the total sum of the walker weights and calculate the scaling for population control
#ifdef USE_MPI_IN_MPSQMC

      double totalSumOfWalkerWeights = 0.0;
      MPI::COMM_WORLD.Allreduce(&mySumOfWalkerWeights, &totalSumOfWalkerWeights, 1, MPI::DOUBLE, MPI::SUM);
      double avgWeight = totalSumOfWalkerWeights / totalNCurrentWalkers;

#else

      double avgWeight = mySumOfWalkerWeights / totalNCurrentWalkers;

#endif
      double scaling = totalNDesiredWalkers / (totalNCurrentWalkers * avgWeight);
      double targetEnergy = log(scaling)/dtau;

      //Update the energy history for each walker, return the fluctuation metric and total projected energy --> uses MPI
      double fluctMetric = EnergyFunctionAndHistory(step, &projectedEnergy, true);

      if (MPIrank==0){
         cout << "        Step = " << step << endl;
         cout << "   # walkers = " << totalNCurrentWalkers << endl;
         cout << " avg(weight) = " << avgWeight << endl;
         cout << "         E_P = " << projectedEnergy << endl;
         cout << "         E_T = " << targetEnergy << endl;
         cout << "       Omega = " << fluctMetric << endl;
         write(step,totalNCurrentWalkers,projectedEnergy, targetEnergy, fluctMetric);
         cout << "---------------------------------------------------------" << endl;
      }

      //Based on scaling, first control the population on each rank separately, and then balance the population over the ranks (uses MPI)
      SeparatePopulationControl(scaling);

#ifdef USE_MPI_IN_MPSQMC
      PopulationBalancing();
#endif

   }

}

/**
 * Here the trotter terms, propagator terms are applied to every walker individually.
 */
double MPSQMC2::PropagateSeparately(){

   for(int cnt = 0;cnt < NThreadsPerRank[MPIrank];cnt++)
      sumWalkerWeightPerThread[cnt] = 0.0;

#pragma omp parallel for schedule(static) default(none)
   for(int walker=0; walker < NCurrentWalkersPerRank[MPIrank]; walker++){

#ifdef _OPENMP
      int myID = omp_get_thread_num();
#else
      int myID = 0;
#endif

      //Propagate with single site terms --> exp (dtau * h * SiZ * 0.5)  for Hterm = -h SiZ
      theWalkers[walker]->gState()->ApplyOneSiteTrotterTermEverywhere(theTrotter);
      theWalkers[walker]->setOverlap(Psi0[myID]); //Normalizes and correctly sets the overlap again

      //Loop over the non-zero couplings
      for(int trotterTerm = 0;trotterTerm < nCouplings;trotterTerm++){

         //Calculate the possible operator combinations
         for (int firstOp = 0;firstOp < trotterSVDsize;firstOp++)
            for (int secondOp = 0;secondOp < trotterSVDsize;secondOp++){

               int combinedIndex = firstOp + trotterSVDsize * secondOp;

               double expectation = TrotterTermsTimesPsi0[myID][trotterTerm][combinedIndex]->InnerProduct(theWalkers[walker]->gState());

               theOperatorCombos[myID][combinedIndex] = expectation;

            }

         //Calculate the PDF for possible moves --> theGrid->gDIM() should be equal to trotterSVDsize!!
         double PDFnorm = 0.0;

         for(int cnt = 0;cnt < theGrid->gNPoints();cnt++){

            thePDF[myID][cnt] = 0.0;

            for(int left = 0;left < trotterSVDsize;left++)
               for(int right = 0;right < trotterSVDsize;right++)
                  thePDF[myID][cnt] += theOperatorCombos[myID][left+trotterSVDsize*right] * theGrid->gCoOfPoint(cnt,left) * theGrid->gCoOfPoint(cnt,right);
            
            thePDF[myID][cnt] = max(0.0, thePDF[myID][cnt]);
            thePDF[myID][cnt] *= trotterSVDsize / ( theWalkers[walker]->gOverlap() * theGrid->gNPoints());

            PDFnorm += thePDF[myID][cnt];

         }

         double oneOverPDFnorm = 1.0 / PDFnorm;

         for(int cnt = 0;cnt < theGrid->gNPoints();cnt++)
            thePDF[myID][cnt] *= oneOverPDFnorm;

         theWalkers[walker]->multiplyWeight(PDFnorm);

         //Get what you should do
         double randomNumber = RN->rand();

         int count = 0;

         double CumulativePDF = thePDF[myID][0];

         while ( (CumulativePDF <= randomNumber) && (count < theGrid->gNPoints()-1) ){

            count += 1;
            CumulativePDF += thePDF[myID][count];

         }

         //Apply the specific term (count) and normalize again
         theWalkers[walker]->gState()->ApplyTwoSiteTrotterTerm(theTrotter, firstIndexCoupling[trotterTerm], secondIndexCoupling[trotterTerm], theGrid, count);
         theWalkers[walker]->setOverlap(Psi0[myID]);

      }

      //Propagate with single site terms --> exp (dtau * h * SiZ * 0.5)  for Hterm = -h SiZ
      theWalkers[walker]->gState()->ApplyOneSiteTrotterTermEverywhere(theTrotter);
      theWalkers[walker]->setOverlap(Psi0[myID]); //Normalizes and correctly sets the overlap again

      sumWalkerWeightPerThread[myID] += theWalkers[walker]->gWeight();

   }

   double sumOfWalkerWeights = 0;

   for (int cnt = 0;cnt < NThreadsPerRank[MPIrank];cnt++)
      sumOfWalkerWeights += sumWalkerWeightPerThread[cnt]; 

   return sumOfWalkerWeights;

}

void MPSQMC2::SeparatePopulationControl(const double scaling){

   const bool debugPrint = true;

   int newNumberOfWalkers = 0;
   double minWeight = 1.0;

   for(int walker=0;walker < NCurrentWalkersPerRank[MPIrank];walker++){

      double scaledWeight = theWalkers[walker]->gWeight() * scaling;

      if(scaledWeight < minWeight)
         minWeight = scaledWeight;

      int nCopies = 1;
      double newWeight = scaledWeight;

      if (newWeight < 0.25){ //Energy doesn't change statistically

         nCopies = (int) ( scaledWeight + RN->rand());
         newWeight = 1.0;

         if (debugPrint)
            cout << "Walker with weight " << scaledWeight << " will be " << nCopies << " times copied with weight 1 ." << endl;

      }

      if (newWeight > 1.5){ //Energy doesn't change

         nCopies = (int) ( scaledWeight + RN->rand());
         newWeight = scaledWeight / nCopies;

         if (debugPrint)
            cout << "Walker with weight " << scaledWeight << " will be " << nCopies << " times copied with weight " << newWeight << " ." << endl;

      }

      if (newNumberOfWalkers + nCopies > myMaxNWalkers)
         cerr << "MPSQMC::Walk (MPI rank " << MPIrank << ") -> The number of desired walkers exceeds the max. allowed number at MPI." << endl;

      if (nCopies==0)
         delete theWalkers[walker];

      if (nCopies>=1){

         theWalkers[walker]->setWeight(newWeight);
         theWalkersCopyArray[newNumberOfWalkers] = theWalkers[walker];

         for (int cnt=1; cnt<nCopies; cnt++)
            theWalkersCopyArray[newNumberOfWalkers+cnt] = new Walker(theWalkers[walker]);

         newNumberOfWalkers += nCopies;

      }
   }

   if (debugPrint)
      cout << "The min. encountered weight on rank " << MPIrank << " is " << minWeight << " ." << endl;

   NCurrentWalkersPerRank[MPIrank] = newNumberOfWalkers;

   //Swap the walker arrays
   Walker ** tempWalkers = theWalkers;
   theWalkers = theWalkersCopyArray;
   theWalkersCopyArray = tempWalkers;

#ifdef USE_MPI_IN_MPSQMC
   totalNCurrentWalkers = 0;
   MPI::COMM_WORLD.Allreduce(NCurrentWalkersPerRank + MPIrank, &totalNCurrentWalkers, 1, MPI::INT, MPI::SUM);
#else
   totalNCurrentWalkers = NCurrentWalkersPerRank[MPIrank];
#endif

}

double MPSQMC2::EnergyFunctionAndHistory(const int step, double * projectedEnergy, bool doFluctuationMetric){

   double fluctMetric_sumOfESq = 0.0; //Per MPI rank
   double fluctMetric_sumOfE   = 0.0; //Per MPI rank

   double projE_numerator   = 0.0; //Per MPI rank
   double projE_denominator = 0.0; //Per MPI rank

   for(int walker = 0;walker < NCurrentWalkersPerRank[MPIrank];walker++){

      const double walkerEnergy = theWalkers[walker]->calcIndividualProjectedEnergy(HPsi0[0]); // <Psi_T | H | walk > / <Psi_T | walk >

      //For the projected energy
      projE_numerator   += theWalkers[walker]->gWeight() * walkerEnergy;
      projE_denominator += theWalkers[walker]->gWeight();

      //For the energy history and the fluctuation metric
      if(doFluctuationMetric){

         theWalkers[walker]->addEnergyToHistory(walkerEnergy);
         const double walkerAvgEnergy = theWalkers[walker]->gEnergyHistory() / step; // e_j(t) : individual time average energy
         fluctMetric_sumOfE    += walkerAvgEnergy;
         fluctMetric_sumOfESq  += walkerAvgEnergy * walkerAvgEnergy;

      }

   }

#ifdef USE_MPI_IN_MPSQMC
   //For the fluctuation metric
   double fluctMetric = 0.0;
   if (doFluctuationMetric){
      double fluctMetric_allE   = 0.0;
      double fluctMetric_allEsq = 0.0;
      MPI::COMM_WORLD.Allreduce(&fluctMetric_sumOfE,   &fluctMetric_allE,   1, MPI::DOUBLE, MPI::SUM);
      MPI::COMM_WORLD.Allreduce(&fluctMetric_sumOfESq, &fluctMetric_allEsq, 1, MPI::DOUBLE, MPI::SUM);
      fluctMetric        = ( fluctMetric_allEsq - fluctMetric_allE * fluctMetric_allE / totalNCurrentWalkers ) / totalNCurrentWalkers;
   }

   //For the projected energy
   double totalNum = 0.0;
   double totalDen = 0.0;
   MPI::COMM_WORLD.Allreduce(&projE_numerator,   &totalNum, 1, MPI::DOUBLE, MPI::SUM);
   MPI::COMM_WORLD.Allreduce(&projE_denominator, &totalDen, 1, MPI::DOUBLE, MPI::SUM);
   projectedEnergy[0] = totalNum / totalDen;
#else
   
   //For the fluctuation metric
   double fluctMetric = 0.0;

   if (doFluctuationMetric)
      fluctMetric = ( fluctMetric_sumOfESq - fluctMetric_sumOfE * fluctMetric_sumOfE / totalNCurrentWalkers ) / totalNCurrentWalkers;

   //For the projected energy
   projectedEnergy[0] = projE_numerator / projE_denominator;
#endif

   return fluctMetric;

}

void MPSQMC2::write(const int step,const int nwalkers,const double projectedEnergy, const double targetEnergy, const double fluctMetric){

   ofstream output("energies.txt",ios::app);
   output.precision(10);
   output << step << "\t\t" << nwalkers << "\t" << projectedEnergy << "\t\t" << targetEnergy << "\t\t" << fluctMetric << endl;
   output.close();

}

void MPSQMC2::BubbleSort(double * values, int * order, const int length){

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

void MPSQMC2::PopulationBalancing(){

#ifdef USE_MPI_IN_MPSQMC
   const bool debuginfoprint = true;
   const double threshold_start = 0.1;
   const double threshold_stop  = 0.03;
   const bool clusterhasinfiniband = true;

   double * Noffset = new double[MPIsize];
   bool oneFracDeviating = false;
   const double fractionScaling = ((double) totalNCurrentWalkers) / totalNDesiredWalkers; //We proportionally want to distribute the total load

   for (int count = 0;count < MPIsize;count++){

      Noffset[count] = NCurrentWalkersPerRank[count] - NDesiredWalkersPerRank[count] * fractionScaling; //So we calculate the offset from "equilibrium"
      double frac = fabs( Noffset[count] ) / ( NDesiredWalkersPerRank[count] * fractionScaling ); //and the percentage of offset

      cout << count << "\t" << Noffset[count] << endl;

      if(frac > threshold_start) 
         oneFracDeviating = true; //if one or more offsets are too large: redistribute       

   }

   
   if ((debuginfoprint) && (MPIrank==0)){

      for (int cnt=0; cnt<MPIsize; cnt++)
         cout << "MPSQMC::PopulationBalancing -> Nwalkers(" << cnt << ") = " << NCurrentWalkersPerRank[cnt] << " and Noffset(" << cnt << ") = " << Noffset[cnt] << endl;

   }

   if (oneFracDeviating){

      if (debuginfoprint){
         double projectedEnergy = 0.0;
         EnergyFunctionAndHistory(1, &projectedEnergy, false);
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
                        double Weight = theWalkers[NCurrentWalkersPerRank[MPIrank]-1]->gWeight();
                        double Overlap = theWalkers[NCurrentWalkersPerRank[MPIrank]-1]->gOverlap();
                        double Ehistory = theWalkers[NCurrentWalkersPerRank[MPIrank]-1]->gEnergyHistory();
                        MPI::COMM_WORLD.Send(&Weight, 1, MPI::DOUBLE, receiver, 0);
                        MPI::COMM_WORLD.Send(&Overlap, 1, MPI::DOUBLE, receiver, 0);
                        MPI::COMM_WORLD.Send(&Ehistory, 1, MPI::DOUBLE, receiver, 0);
                        for (int site=0; site<theMPO->gLength(); site++){
                           int dim = Psi0[0]->gPhys_d() * Psi0[0]->gDimAtBound(site) * Psi0[0]->gDimAtBound(site+1);
                           double * SendStorage = theWalkers[NCurrentWalkersPerRank[MPIrank]-1]->gState()->gMPStensor(site)->gStorage();
                           MPI::COMM_WORLD.Send(SendStorage, dim, MPI::DOUBLE, receiver, 0);
                        }
                        delete theWalkers[NCurrentWalkersPerRank[MPIrank]-1];
                        NCurrentWalkersPerRank[MPIrank]  -= 1;
                        NCurrentWalkersPerRank[receiver] += 1;
                     }
                  }

                  if (MPIrank == receiver){
                     for (int counter=0; counter<amount; counter++){
                        double Weight = 0.0;
                        double Overlap = 0.0;
                        double Ehistory = 0.0;
                        MPI::Status status;
                        MPI::COMM_WORLD.Recv(&Weight, 1, MPI::DOUBLE, sender, 0, status);
                        MPI::COMM_WORLD.Recv(&Overlap, 1, MPI::DOUBLE, sender, 0, status);
                        MPI::COMM_WORLD.Recv(&Ehistory, 1, MPI::DOUBLE, sender, 0, status);
                        theWalkers[NCurrentWalkersPerRank[MPIrank]] = new Walker(Psi0[0], Overlap, Weight, Ehistory);
                        for (int site=0; site<theMPO->gLength(); site++){
                           int dim = Psi0[0]->gPhys_d() * Psi0[0]->gDimAtBound(site) * Psi0[0]->gDimAtBound(site+1);
                           double * RecvStorage = theWalkers[NCurrentWalkersPerRank[MPIrank]]->gState()->gMPStensor(site)->gStorage();
                           MPI::COMM_WORLD.Recv(RecvStorage, dim, MPI::DOUBLE, sender, 0, status);
                        }
                        NCurrentWalkersPerRank[MPIrank] += 1;
                        NCurrentWalkersPerRank[sender]  -= 1;
                     }
                  }

                  if ((MPIrank != sender) && (MPIrank != receiver)){
                     NCurrentWalkersPerRank[sender]   -= amount;
                     NCurrentWalkersPerRank[receiver] += amount;
                  }

               }

            }

         } //Everything got communicated in parallel

         MPI::COMM_WORLD.Barrier();

         //Determine the offsets again : now with threshold_stop!!!!
         oneFracDeviating = false;
         for (int count=0; count<MPIsize; count++){
            Noffset[count] = NCurrentWalkersPerRank[count] - NDesiredWalkersPerRank[count] * fractionScaling;
            double frac = fabs( Noffset[count] ) / ( NDesiredWalkersPerRank[count] * fractionScaling );
            if (frac > threshold_stop){ oneFracDeviating = true; }
         }

      }

      delete [] work;

      if (debuginfoprint){
         double projectedEnergy = 0.0;
         EnergyFunctionAndHistory(1, &projectedEnergy, false);
         if (MPIrank==0){ cout << "MPSQMC::PopulationBalancing -> As a check: projected E after  = " << projectedEnergy << endl; }
      }

   }

   delete [] Noffset;
#endif

}

