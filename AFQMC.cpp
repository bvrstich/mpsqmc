#include "MPItag.h" //Here USE_MPI_IN_MPSQMC can be defined to switch to non-MPI compilers.

#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include <omp.h>
#include <vector>

#include "AFQMC.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 15, 2013 */
using namespace std;

/**
 * constructor of the AFQMC object, takes input parameters that define the QMC walk.
 * @param theMPO MPO containing relevant matrix elements, defines system being studied
 * @param Psi0 input trialwavefunction
 * @param Dtrunc dimension of the walkers
 * @param Nwalkers number of Walker states
 * @param dtau time step of each evolution
 */
AFQMC::AFQMC(HeisenbergMPO * theMPO, Random * RN, MPSstate *Psi0_in,const int DW, const int Nwalkers, const double dtau){

   this->theMPO = theMPO;
   this->RN = RN;
   this->DT = Psi0_in->gDtrunc();
   this->DW = DW;
   this->dtau = dtau;
   this->phys_d = theMPO->gPhys_d();

   this->theTrotter = new TrotterHeisenberg(theMPO,dtau);
   n_trot = theTrotter->gn_trot();

   this->totalNDesiredWalkers = Nwalkers;

   SetupOMPandMPILoadDistribution();

   myMaxNWalkers = max(1000,3*NDesiredWalkersPerRank[MPIrank]);

#ifdef USE_MPI_IN_MPSQMC
   MPI::COMM_WORLD.Barrier();
#endif

   Psi0 = new MPSstate(Psi0_in);

   SetupTrial();

#ifdef USE_MPI_IN_MPSQMC
   MPI::COMM_WORLD.Barrier();
#endif

   SetupWalkers(true);

}

void AFQMC::SetupOMPandMPILoadDistribution(){

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
   NThreadsPerRank          = new int [MPIsize];
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

   if(MPIrank==0){

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

AFQMC::~AFQMC(){

   //AFQMC::AFQMC
   delete theTrotter;

   //AFQMC::SetupTrial
   delete Psi0;
   delete HPsi0;

   for(int k = 0;k < 3*n_trot;++k){

      delete VPsi0[k];
      delete V2Psi0[k];

   }

   delete [] VPsi0;
   delete [] V2Psi0;

   //AFQMC::SetupWalkers
   for (int cnt = 0;cnt < theWalkers.size(); cnt++)
      delete theWalkers[cnt];

   //AFQMC::SetupOMPandMPILoadDistribution
   delete [] NThreadsPerRank;
   delete [] NDesiredWalkersPerRank;
   delete [] NCurrentWalkersPerRank;

}

/**
 * construct the trial wavefunction
 */
void AFQMC::SetupTrial(){

   if(MPIrank==0)
      Psi0->LeftNormalize();

#ifdef USE_MPI_IN_MPSQMC
   Psi0 = BroadcastCopyConstruct(Psi0);
#endif

   //Rank 0 calculates MPO times trial, and the result gets copied so every thread on every rank has 1 copy.
   if(MPIrank==0){

      HPsi0 = new MPSstate(theMPO->gLength(),DT,theMPO->gPhys_d(),RN);
      HPsi0->ApplyMPO(false,theMPO, Psi0);
      HPsi0->CompressState(); //Compression only throws away Schmidt values which are numerically zero...

   }

#ifdef USE_MPI_IN_MPSQMC
   HPsi0 = BroadcastCopyConstruct(HPsi0);
#endif

   //now Apply the hermitian conjugate of the V's times trialstate
   if(MPIrank==0){

      VPsi0 = new MPSstate * [3*n_trot];
      V2Psi0 = new MPSstate * [3*n_trot];

      for(int r = 0;r < 3;++r)
         for(int k = 0;k < n_trot;++k){

         VPsi0[r*n_trot + k] = new MPSstate(theMPO->gLength(),DT,theMPO->gPhys_d(),RN);

         VPsi0[r*n_trot + k]->ApplyMPO(true,theTrotter->gV_Op(k,r) , Psi0);
         VPsi0[r*n_trot + k]->CompressState(); //Compression only throws away Schmidt values which are numerically zero...

         V2Psi0[r*n_trot + k] = new MPSstate(theMPO->gLength(),DT,theMPO->gPhys_d(),RN);

         V2Psi0[r*n_trot + k]->ApplyMPO(true,theTrotter->gV_Op(k,r) , VPsi0[r*n_trot + k]);
         V2Psi0[r*n_trot + k]->CompressState(); //Compression only throws away Schmidt values which are numerically zero...

      }

   }

#ifdef USE_MPI_IN_MPSQMC
   for(int r = 0;r < 3;++r)
      for(int k = 0;k < length;++k)
         VPsi0[r*n_trot + k] = BroadcastCopyConstruct(VPsi0[r*n_trot + k]);
#endif

}

/** 
 * different function
 */
MPSstate * AFQMC::BroadcastCopyConstruct(MPSstate * pointer){

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
void AFQMC::SetupWalkers(bool copyTrial){

   theWalkers.resize(NCurrentWalkersPerRank[MPIrank]);

   if(copyTrial){

      theWalkers[0] = new Walker(Psi0,1.0,n_trot);

      if(DW < DT)
         theWalkers[0]->gState()->CompressState(DW);//compress the state to Walker D

      theWalkers[0]->sOverlap(Psi0);
      theWalkers[0]->sEL(HPsi0);
      theWalkers[0]->sVL(VPsi0);

   }
   else{

      theWalkers[0] = new Walker(theMPO->gLength(), DW, theMPO->gPhys_d(), n_trot,RN);

      theWalkers[0]->sOverlap(Psi0);
      theWalkers[0]->sEL(HPsi0);
      theWalkers[0]->sVL(VPsi0);

   }

   for(int cnt = 1;cnt < theWalkers.size();cnt++){

      if(copyTrial)
         theWalkers[cnt] = new Walker(theWalkers[0]);
      else{

         theWalkers[cnt] = new Walker(theMPO->gLength(), DW, theMPO->gPhys_d(), n_trot,RN);
         
         theWalkers[cnt]->sOverlap(Psi0);
         theWalkers[cnt]->sEL(HPsi0);
         theWalkers[cnt]->sVL(VPsi0);

      }

   }

}

void AFQMC::Walk(const int steps){

   complex<double> EP = gEP();

   if (MPIrank==0){

      char filename[100];
      sprintf(filename,"output/Heisenberg1D/ener_L%dDT%dDW%d.txt",theMPO->gLength(),DT,DW);

      ofstream output(filename,ios::trunc);

      output << "#Step\t\tE_P\t\tE_T\t" << endl;
      output.close();

      cout << "Energy at start = " << EP << endl;
      cout << "---------------------------------------------------------" << endl;

   }

   for (int step=1; step<=steps; step++){

      //Propagate the walkers of each rank separately --> no MPI in that function
      double wsum = PropagateSeparately();

      //Form the total sum of the walker weights and calculate the scaling for population control
#ifdef USE_MPI_IN_MPSQMC
      MPI::COMM_WORLD.Barrier();

      double twsum = 0.0;
      MPI::COMM_WORLD.Allreduce(&wsum, &twsum, 1, MPI::DOUBLE, MPI::SUM);
      double avgw = twsum /(double) theWalkers.size();
#else
      double avgw = wsum / (double)theWalkers.size();
#endif
      double scaling = totalNDesiredWalkers / wsum;

      double ET = log(scaling)/dtau;

      //Update the energy history for each walker, return the fluctuation metric and total projected energy --> uses MPI
      EP = gEP();

      if (MPIrank==0){

         cout << "        Step = " << step << endl;
         cout << "   # walkers = " << theWalkers.size() << endl;
         cout << " avg(weight) = " << avgw << endl;
         cout << "         E_P = " << EP << endl;
         cout << "         E_T = " << ET << endl;
         cout << "---------------------------------------------------------" << endl;

         write(step,theWalkers.size(),std::real(EP), ET);

      }

      //Based on scaling, first control the population on each rank separately, and then balance the population over the ranks (uses MPI)
      SeparatePopulationControl(scaling);

      double min_en = 0.0;
      double min_ov = 1.0;

      for(int i = 0;i < theWalkers.size();++i){

         if(min_en > std::real(theWalkers[i]->gEL()))
            min_en = std::real(theWalkers[i]->gEL());

         if(min_ov > std::abs(theWalkers[i]->gOverlap()))
            min_ov = std::abs(theWalkers[i]->gOverlap());

      }

      cout << "Minimal Energy:\t" << min_en << endl;
      cout << "Minimal Overlap:\t" << min_ov << endl;

#ifdef USE_MPI_IN_MPSQMC
      MPI::COMM_WORLD.Barrier();
#endif

#ifdef USE_MPI_IN_MPSQMC
      PopulationBalancing();
#endif

   }

}

/**
 * Here the trotter terms, propagator terms are applied to every walker individually.
 */
double AFQMC::PropagateSeparately(){

   double sum = 0.0;

#pragma omp parallel for reduction(+: sum)
   for(int walker=0; walker < theWalkers.size(); walker++){

      //Propagate with single site terms --> exp (dtau * h * SiZ * 0.5)  for Hterm = -h SiZ, if a field is present
      if(theTrotter->gIsMagneticField())
         theWalkers[walker]->gState()->ApplyH1(theTrotter);

      //now loop over the auxiliary fields:
      for(int k = 0;k < n_trot;++k)
         for(int r = 0;r < 3;++r){

         double x = RN->normal();

         complex<double> shift = theWalkers[walker]->gVL(k,r);
         theWalkers[walker]->gState()->ApplyAF(k,r,(complex<double>)x + shift,theTrotter);

      }

      //Propagate with single site terms --> exp (dtau * h * SiZ * 0.5)  for Hterm = -h SiZ
      if(theTrotter->gIsMagneticField())
         theWalkers[walker]->gState()->ApplyH1(theTrotter);

      theWalkers[walker]->sOverlap(Psi0);

      complex<double> prev_EL = theWalkers[walker]->gEL();

      theWalkers[walker]->sEL(HPsi0);

      complex<double> EL = theWalkers[walker]->gEL();

      double scale = exp(-0.5 * dtau * std::real(EL + prev_EL));

      theWalkers[walker]->multWeight(scale);

      theWalkers[walker]->sVL(VPsi0);

      sum += theWalkers[walker]->gWeight();

   }

   return sum;

}

void AFQMC::SeparatePopulationControl(const double scaling){

   double minw = 1.0;
   double maxw = 1.0;

   double sum = 0.0;

   for(int walker = 0;walker < theWalkers.size();walker++){

      theWalkers[walker]->multWeight(scaling);

      double weight = theWalkers[walker]->gWeight();

      if(weight < minw)
         minw = weight;

      if(weight > maxw)
         maxw = weight;

      if (weight < 0.25){ //Energy doesn't change statistically

         int nCopies = (int) ( weight + RN->rand());

         if(nCopies == 0){

            cout << "Walker with weight " << weight << " will be deleted." << endl;

            delete theWalkers[walker];

            theWalkers.erase(theWalkers.begin() + walker);

         }
         else
            theWalkers[walker]->sWeight(1.0);

      }

      if(weight > 1.5){ //statically energy doesn't change

         int nCopies =(int) ( weight + RN->rand());
         double new_weight = weight / (double) nCopies;

         theWalkers[walker]->sWeight(new_weight);

         cout << "Walker with weight " << weight << " will be " << nCopies << " times copied with weight " << new_weight << "." << endl;

         for(int i = 1;i < nCopies;++i){

            Walker *nw = new Walker(theWalkers[walker]);

            theWalkers.push_back(nw);

         }

      }

      sum += weight;

   }

   cout << endl;
   cout << "total weight:\t" << sum << endl;
   cout << endl;

   cout << "The min. encountered weight on rank " << MPIrank << " is " << minw << " ." << endl;
   cout << "The max. encountered weight on rank " << MPIrank << " is " << maxw << " ." << endl;

#ifdef USE_MPI_IN_MPSQMC
   totalNCurrentWalkers = 0;
   MPI::COMM_WORLD.Allreduce(NCurrentWalkersPerRank + MPIrank, &totalNCurrentWalkers, 1, MPI::INT, MPI::SUM);
#endif

}

complex<double> AFQMC::gEP(){

   complex<double> projE_num = 0.0;
   complex<double> projE_den = 0.0;

   for(int walker = 0;walker < theWalkers.size();walker++){

      const complex<double> walkerEnergy = theWalkers[walker]->gEL(); // <Psi_T | H | walk > / <Psi_T | walk >

      //For the projected energy
      projE_num   += theWalkers[walker]->gWeight() * walkerEnergy;//* theWalkers[walker]->gOverlap();
      projE_den += theWalkers[walker]->gWeight();// * theWalkers[walker]->gOverlap();

   }

   complex<double> EP = 0.0;

#ifdef USE_MPI_IN_MPSQMC
   //For the projected energy
   double totalNum = 0.0;
   double totalDen = 0.0;

   MPI::COMM_WORLD.Allreduce(&projE_num, &totalNum, 1, MPI::DOUBLE, MPI::SUM);
   MPI::COMM_WORLD.Allreduce(&projE_den, &totalDen, 1, MPI::DOUBLE, MPI::SUM);

   EP = totalE_num / totalE_den;
#else
   EP = projE_num / projE_den;
#endif

   return EP;

}


void AFQMC::write(const int step,const int nwalkers,const double EP, const double ET){

   char filename[100];
   sprintf(filename,"output/Heisenberg1D/ener_L%dDT%dDW%d.txt",theMPO->gLength(),DT,DW);
   ofstream output(filename,ios::app);
   output.precision(10);
   output << step << "\t\t" << nwalkers << "\t" << EP << "\t\t" << ET << endl;
   output.close();

}

/*
   void AFQMC::BubbleSort(double * values, int * order, const int length){

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

   void AFQMC::PopulationBalancing(){

#ifdef USE_MPI_IN_MPSQMC
const bool debuginfoprint = true;
const double threshold_start = 0.1;
const double threshold_stop  = 0.03;
const bool clusterhasinfiniband = true;

double * Noffset = new double[MPIsize];
bool oneFracDeviating = false;
const double fractionScaling = ((double) totalNCurrentWalkers) / totalNDesiredWalkers; //We proportionally want to distribute the total load

//send and recieve to syncronize NCurrentWalkersPerRank!
if(MPIrank > 0)
MPI_Send(&NCurrentWalkersPerRank[MPIrank],1,MPI_INT,0,MPIrank,MPI_COMM_WORLD);

MPI::COMM_WORLD.Barrier();

//receive
if(MPIrank == 0)
for(int rank = 1;rank < MPIsize;++rank)
MPI_Recv(&NCurrentWalkersPerRank[rank],1,MPI_INT,rank,rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

//broadcast to other ranks:
MPI_Bcast(NCurrentWalkersPerRank,MPIsize,MPI_INT,0,MPI_COMM_WORLD);

for (int count = 0;count < MPIsize;count++){

Noffset[count] = NCurrentWalkersPerRank[count] - NDesiredWalkersPerRank[count] * fractionScaling; //So we calculate the offset from "equilibrium"
double frac = fabs( Noffset[count] ) / ( NDesiredWalkersPerRank[count] * fractionScaling ); //and the percentage of offset

if(frac > threshold_start) 
oneFracDeviating = true; //if one or more offsets are too large: redistribute       

}

if ((debuginfoprint) && (MPIrank==0)){

for (int cnt=0; cnt<MPIsize; cnt++)
cout << "MPSQMC::PopulationBalancing -> Nwalkers(" << cnt << ") = " << NCurrentWalkersPerRank[cnt] << " and Noffset(" << cnt << ") = " << Noffset[cnt] << endl;

}

if(oneFracDeviating){

if(debuginfoprint){

double projectedEnergy = 0.0;
EnergyFunctionAndHistory(1, &projectedEnergy, false);

if (MPIrank==0)
cout << "MPSQMC::PopulationBalancing -> As a check: projected E before = " << projectedEnergy << endl;

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
*/

void AFQMC::testProp(){

   double sum = 0.0;

#pragma omp parallel for reduction(+: sum)
   for(int walker = 0;walker < theWalkers.size();++walker){

      for(int k = 0;k < n_trot;++k)
         for(int r = 0;r < 3;++r){

            double x = RN->normal();
            theWalkers[walker]->gState()->ApplyAF(k,r,(complex<double>)x,theTrotter);

         }

      sum += std::real(theWalkers[walker]->gState()->InnerProduct(Psi0));

   }

   cout << sum/theWalkers.size() << endl;
   //cout << log((double)theWalkers.size()/sum)/dtau << endl;

   complex<double> tmp(0.0,0.0);

   for(int k = 0;k < n_trot;++k)
      for(int r = 0;r < 3;++r)
         tmp += Psi0->InnerProduct(V2Psi0[r*n_trot + k]);

   cout << "Energy sum aux fields:\t" << -0.5 * tmp/dtau << endl;

}
