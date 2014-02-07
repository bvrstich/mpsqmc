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
 * @param theTrotter trotter terms for the interaction
 * @param Psi0 input trialwavefunction
 * @param Dtrunc dimension of the walkers
 * @param Nwalkers number of Walker states
 * @param dtau time step of each evolution
 */
AFQMC::AFQMC(J1J2MPO *theMPO,TrotterJ1J2 *theTrotter,Random *RN, MPSstate *Psi0_in,const int DW, const int Nwalkers, const double dtau){

   this->RN = RN;
   this->DT = Psi0_in->gDtrunc();
   this->DW = DW;
   this->dtau = dtau;
   this->phys_d = theTrotter->gPhys_d();
   this->theTrotter = theTrotter;
   this->theMPO = theMPO;

   n_trot = theTrotter->gn_trot();

   this->totalNDesiredWalkers = Nwalkers;

   SetupOMPandMPILoadDistribution();

   myMaxNWalkers = max(1000,3*NDesiredWalkersPerRank[MPIrank]);

#ifdef USE_MPI_IN_MPSQMC
   MPI::COMM_WORLD.Barrier();
#endif

   Psi0 = new MPSstate(Psi0_in);

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

   //Distribute the workload according to the threads
   NDesiredWalkersPerRank = new int[MPIsize];

   int remainder = totalNDesiredWalkers;

   for(int count=0; count<MPIsize; count++){

      NDesiredWalkersPerRank[count] = (int)(((double) totalNDesiredWalkers ) / (double) MPIsize);
      remainder -= NDesiredWalkersPerRank[count];

   }

   for (int count=0; count<remainder; count++)
      NDesiredWalkersPerRank[count] += 1;

   if(MPIrank==0){

      cout << "There are " << MPIsize << " MPI processes." << endl;

      for (int cnt=0; cnt<MPIsize; cnt++)
         cout << "   MPI rank " << cnt << " carries " << NDesiredWalkersPerRank[cnt] << " walkers." << endl;
      
   }

#else
   MPIsize = 1;
   MPIrank = 0;
   NDesiredWalkersPerRank          = new int[MPIsize];
   NCurrentWalkersPerRank          = new int[MPIsize];
   NDesiredWalkersPerRank[MPIrank] = totalNDesiredWalkers;
#endif

}

AFQMC::~AFQMC(){

   //AFQMC::SetupTrial
   delete Psi0;

   //AFQMC::SetupWalkers
   for (int cnt = 0;cnt < theWalkers.size(); cnt++)
      delete theWalkers[cnt];

   //AFQMC::SetupOMPandMPILoadDistribution
   delete [] NDesiredWalkersPerRank;

}

/**
 * initialize the walkers
 */
void AFQMC::SetupWalkers(bool copyTrial){

   theWalkers.resize(NDesiredWalkersPerRank[MPIrank]);

   if(copyTrial){

      int L = sqrt(theTrotter->glength());

      double J2 = theTrotter->gJ2();
      int j2 = 10*J2;

      char filename[100];

      if(J2 == 10)
         sprintf(filename,"/home/bright/bestanden/programmas/dmrg/J1J2/%dx%d/J2=1.0/PsiW/DT=%d.mps",L,L,DT);
      else
         sprintf(filename,"/home/bright/bestanden/programmas/dmrg/J1J2/%dx%d/J2=0.%d/PsiW/DT=%d.mps",L,L,j2,DT);

      MPSstate input(filename,RN);
      stor = new MPSstate(filename,RN);

      theWalkers[0] = new Walker(&input,1.0,n_trot);

      theWalkers[0]->sOverlap(Psi0);
      theWalkers[0]->sEL(theMPO,Psi0);
      theWalkers[0]->sVL(theTrotter,Psi0);

   }
   else{

      theWalkers[0] = new Walker(theTrotter->glength(), DW, theTrotter->gPhys_d(), n_trot,RN);
      theWalkers[0]->gState()->normalize();

      theWalkers[0]->sOverlap(Psi0);
      theWalkers[0]->sEL(theMPO,Psi0);
      theWalkers[0]->sVL(theTrotter,Psi0);

   }

   for(int cnt = 1;cnt < theWalkers.size();cnt++){

      if(copyTrial)
         theWalkers[cnt] = new Walker(theWalkers[0]);
      else{

         theWalkers[cnt] = new Walker(theTrotter->glength(), DW, theTrotter->gPhys_d(), n_trot,RN);
         theWalkers[cnt]->gState()->normalize();

         theWalkers[cnt]->sOverlap(Psi0);
         theWalkers[cnt]->sEL(theMPO,Psi0);
         theWalkers[cnt]->sVL(theTrotter,Psi0);

      }

   }

}

void AFQMC::Walk(const int steps){

   complex<double> EP = gEP();

   if (MPIrank==0){

      int L = sqrt(theTrotter->glength());
      double J2 = theTrotter->gJ2();

      int j2 = 10*J2;

      char filename[100];

      if(j2 == 10)
         sprintf(filename,"output/J1J2/%dx%d/J2=1.0/DT%dDW%d.txt",L,L,DT,DW);
      else
         sprintf(filename,"output/J1J2/%dx%d/J2=0.%d/DT%dDW%d.txt",L,L,j2,DT,DW);

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

      MPI::COMM_WORLD.Barrier();

      totalNWalkers = 0;

      int nwalk = theWalkers.size();

      MPI::COMM_WORLD.Allreduce(&nwalk, &totalNWalkers, 1, MPI::INT, MPI::SUM);

      double avgw = twsum /(double) totalNWalkers;
#else
      double avgw = wsum / (double)theWalkers.size();
#endif
      double scaling = totalNDesiredWalkers / (avgw*totalNWalkers);

      double ET = log(scaling)/dtau;

      //Update the energy history for each walker, return the fluctuation metric and total projected energy --> uses MPI
      EP = gEP();

      if (MPIrank==0){

         cout << "        Step = " << step << endl;
         cout << "   # walkers = " << totalNWalkers << endl;
         cout << " avg(weight) = " << avgw << endl;
         cout << "         E_P = " << EP << endl;
         cout << "         E_T = " << ET << endl;
         cout << "---------------------------------------------------------" << endl;

         write(step,theWalkers.size(),std::real(EP), ET);

      }

      //Based on scaling, first control the population on each rank separately, and then balance the population over the ranks (uses MPI)
      SeparatePopulationControl(scaling);

#ifdef USE_MPI_IN_MPSQMC
      MPI::COMM_WORLD.Barrier();

      PopulationBalancing();
#endif

   }

}

/**
 * Here the trotter terms, propagator terms are applied to every walker individually.
 */
double AFQMC::PropagateSeparately(){

   double sum = 0.0;
   double width = sqrt(2.0/dtau);

   int num_rej = 0;

   for(int walker=0; walker < theWalkers.size(); walker++){

      //backup the state
      stor->copy(theWalkers[walker]->gState());

      //now loop over the auxiliary fields:
      for(int k = 0;k < n_trot;++k)
         for(int r = 0;r < 3;++r){

            double x = RN->normal();

            complex<double> shift = theWalkers[walker]->gVL(k,r);
            theWalkers[walker]->gState()->ApplyAF(k,r,(complex<double>)x + shift,theTrotter);

         }

      //get the local energy
      complex<double> tmpOverlap = theWalkers[walker]->gState()->InnerProduct(Psi0);
      complex<double> tmpEL = theWalkers[walker]->gState()->expectation(theMPO,Psi0)/tmpOverlap;

      if( (std::real(tmpEL) < std::real(theWalkers[walker]->gEL()) - width) || 

            (std::real(tmpEL) > std::real(theWalkers[walker]->gEL()) + width) ){//very rare event, will cause numerical unstability

         num_rej++;

         //copy the state back!
         theWalkers[walker]->gState()->copy(stor);

      }
      else{//go on

         theWalkers[walker]->sOverlap(Psi0);

         complex<double> prev_EL = theWalkers[walker]->gEL();

         theWalkers[walker]->sEL(tmpEL);

         complex<double> EL = theWalkers[walker]->gEL();

         double scale = exp(-0.5 * dtau * std::real(EL + prev_EL));

         theWalkers[walker]->multWeight(scale);

         theWalkers[walker]->sVL(theTrotter,Psi0);

         sum += theWalkers[walker]->gWeight();

      }

   }

   MPI::COMM_WORLD.Barrier();

   int tot_num_rej = 0;

   MPI::COMM_WORLD.Allreduce(&num_rej, &tot_num_rej, 1, MPI::INT, MPI::SUM);

   if(MPIrank == 0){

      cout << endl;
      cout << "Number of rejected moves:\t" << tot_num_rej << endl;
      cout << endl;

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

            //cout << "Walker with weight " << weight << " will be deleted." << endl;

            delete theWalkers[walker];

            theWalkers.erase(theWalkers.begin() + walker);

         }
         else
            theWalkers[walker]->sWeight(1.0);

      }
      //else{
      if(weight > 1.5){ //statically energy doesn't change

         int nCopies =(int) ( weight + RN->rand());
         double new_weight = weight / (double) nCopies;

         theWalkers[walker]->sWeight(new_weight);

         //cout << "Walker with weight " << weight << " will be " << nCopies << " times copied with weight " << new_weight << "." << endl;

         for(int i = 1;i < nCopies;++i){

            Walker *nw = new Walker(theWalkers[walker]);

            theWalkers.push_back(nw);

         }

      }

      sum += weight;

   }

}

complex<double> AFQMC::gEP(){

   complex<double> projE_num = 0.0;
   complex<double> projE_den = 0.0;

   for(int walker = 0;walker < theWalkers.size();walker++){

      const complex<double> walkerEnergy = theWalkers[walker]->gEL(); // <Psi_T | H | walk > / <Psi_T | walk >

      //For the projected energy
      projE_num   += theWalkers[walker]->gWeight() * walkerEnergy;
      projE_den += theWalkers[walker]->gWeight();

   }

   complex<double> EP = 0.0;

#ifdef USE_MPI_IN_MPSQMC
   //For the projected energy
   complex<double> totalNum = 0.0;
   complex<double> totalDen = 0.0;

   MPI::COMM_WORLD.Allreduce(&projE_num, &totalNum, 1, MPI::DOUBLE_COMPLEX, MPI::SUM);
   MPI::COMM_WORLD.Allreduce(&projE_den, &totalDen, 1, MPI::DOUBLE_COMPLEX, MPI::SUM);

   EP = totalNum / totalDen;
#else
   EP = projE_num / projE_den;
#endif

   return EP;

}

void AFQMC::write(const int step,const int nwalkers,const double EP, const double ET){

   char filename[100];

   int L = sqrt(theTrotter->glength());
   double J2 = theTrotter->gJ2();

   int j2 = 10*J2;

   if(j2 == 10)
      sprintf(filename,"output/J1J2/%dx%d/J2=1.0/DT%dDW%d.txt",L,L,DT,DW);
   else
      sprintf(filename,"output/J1J2/%dx%d/J2=0.%d/DT%dDW%d.txt",L,L,j2,DT,DW);

   ofstream output(filename,ios::app);
   output.precision(10);
   output << step << "\t\t" << nwalkers << "\t" << EP << "\t\t" << ET << endl;
   output.close();

}

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
   const double threshold_start = 0.1;
   const double threshold_stop  = 0.03;
   const bool clusterhasinfiniband = true;

   double * Noffset = new double[MPIsize];
   bool oneFracDeviating = false;

   totalNWalkers = 0;
   int nwalk = theWalkers.size();

   MPI::COMM_WORLD.Allreduce(&nwalk, &totalNWalkers, 1, MPI::INT, MPI::SUM);

   const double fractionScaling = (double) totalNWalkers / (double)totalNDesiredWalkers; //We proportionally want to distribute the total load

   MPI::COMM_WORLD.Barrier();

   Noffset[MPIrank] = theWalkers.size() - NDesiredWalkersPerRank[MPIrank] * fractionScaling; //So we calculate the offset from "equilibrium"

   MPI::COMM_WORLD.Barrier();

   //send and receive the offsets
   if(MPIrank > 0)
      MPI_Send(&Noffset[MPIrank],1,MPI_DOUBLE,0,MPIrank,MPI_COMM_WORLD);

   MPI::COMM_WORLD.Barrier();

   //receive
   if(MPIrank == 0)
      for(int rank = 1;rank < MPIsize;++rank)
         MPI_Recv(&Noffset[rank],1,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

   if(MPIrank == 0){
      cout << endl;
      for(int count = 0;count < MPIsize;++count){

         double frac = fabs( Noffset[count] ) / ( NDesiredWalkersPerRank[count] * fractionScaling ); //and the percentage of offset

         cout << "Deviation of average occupation on rank\t" << count << " = " << Noffset[count] << "\tFractional offset:\t" << frac << endl;

         if(frac > threshold_start) 
            oneFracDeviating = true; //if one or more offsets are too large: redistribute       

      }
      cout << endl;

   }

   //broadcast!
   MPI_Bcast(Noffset,MPIsize,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&oneFracDeviating,1,MPI::BOOL,0,MPI_COMM_WORLD);

   if(oneFracDeviating){

      int * work = new int[MPIsize];
      int communication_round = 0;

      while (oneFracDeviating ){

         communication_round++;

         //Do a bubble sort from large to small (positive to negative). Not optimal algo, optimal = quicksort (dlasrt_) --> but multi-array sort in lapack?
         BubbleSort(Noffset, work, MPIsize); // Now for all index: Noffset [ work[index] ] >= Noffset[ work[index+1] ]

         //Communicate walkers between rank work[comm] and rank work[MPIsize - 1 - comm] until "drained"
         for (int comm=0; comm < ((clusterhasinfiniband) ? MPIsize/2 : 1); comm++){

            const int sender = work[comm];
            const int receiver = work[MPIsize - 1 - comm];

            if ((Noffset[sender] > 0.0) && (Noffset[receiver] < 0.0)){

               int amount = (int) min(Noffset[sender], -Noffset[receiver]); //This explains the "until drained" statement.

               if (amount > 0){

                  if (MPIrank == sender){

                     for (int counter=0; counter<amount; counter++){

                        double weight = theWalkers[theWalkers.size() - 1]->gWeight();
                        complex<double> overlap = theWalkers[theWalkers.size()-1]->gOverlap();
                        complex<double> EL = theWalkers[theWalkers.size()-1]->gEL();
                        complex<double> *VL = theWalkers[theWalkers.size()-1]->gVL();

                        MPI::COMM_WORLD.Send(&weight, 1, MPI::DOUBLE, receiver, 0);
                        MPI::COMM_WORLD.Send(&overlap, 1, MPI::DOUBLE_COMPLEX, receiver, 0);
                        MPI::COMM_WORLD.Send(&EL, 1, MPI::DOUBLE_COMPLEX, receiver, 0);

                        int dim = 3*n_trot;
                        MPI::COMM_WORLD.Send(VL,dim , MPI::DOUBLE_COMPLEX, receiver, 0);

                        for (int site=0; site < Psi0->gLength(); site++){

                           int dim = Psi0->gPhys_d() * theWalkers[theWalkers.size() - 1]->gState()->gDimAtBound(site) * theWalkers[theWalkers.size() - 1]->gState()->gDimAtBound(site+1);

                           complex<double> * SendStorage = theWalkers[theWalkers.size() - 1]->gState()->gMPStensor(site)->gStorage();

                           MPI::COMM_WORLD.Send(SendStorage, dim, MPI::DOUBLE_COMPLEX, receiver, 0);

                        }

                        //remove the last element
                        delete theWalkers[theWalkers.size() - 1];
                        theWalkers.pop_back();

                     }

                  }

                  if (MPIrank == receiver){

                     for (int counter=0; counter<amount; counter++){

                        double weight = 0.0;
                        complex<double> overlap(0.0,0.0);
                        complex<double> EL(0.0,0.0);

                        MPI::Status status;

                        MPI::COMM_WORLD.Recv(&weight, 1, MPI::DOUBLE, sender, 0, status);
                        MPI::COMM_WORLD.Recv(&overlap, 1, MPI::DOUBLE_COMPLEX, sender, 0, status);
                        MPI::COMM_WORLD.Recv(&EL, 1, MPI::DOUBLE_COMPLEX, sender, 0, status);

                        int dim = 3*n_trot;
                        complex<double> *VL = new complex<double> [dim];

                        MPI::COMM_WORLD.Recv(VL,dim , MPI::DOUBLE_COMPLEX, sender, 0,status);

                        Walker *walk;
                        theWalkers.push_back(walk);

                        theWalkers[theWalkers.size() - 1] = new Walker(theWalkers[0]->gState(),weight, overlap, EL,VL,n_trot);

                        for(int site = 0;site < Psi0->gLength();site++){

                           int dim = Psi0->gPhys_d() * theWalkers[theWalkers.size() - 1]->gState()->gDimAtBound(site) * theWalkers[theWalkers.size() - 1]->gState()->gDimAtBound(site+1);

                           complex<double> * RecvStorage = theWalkers[theWalkers.size() - 1]->gState()->gMPStensor(site)->gStorage();

                           MPI::COMM_WORLD.Recv(RecvStorage, dim, MPI::DOUBLE_COMPLEX, sender, 0, status);

                        }

                        delete [] VL;

                     }

                  }

               }

            }

         } //Everything got communicated in parallel

         MPI::COMM_WORLD.Barrier();

         //Determine the offsets again : now with threshold_stop!!!!
         oneFracDeviating = false;

         Noffset[MPIrank] = theWalkers.size() - NDesiredWalkersPerRank[MPIrank] * fractionScaling; //So we calculate the offset from "equilibrium"

         MPI::COMM_WORLD.Barrier();

         //send and receive the offsets
         if(MPIrank > 0)
            MPI_Send(&Noffset[MPIrank],1,MPI_DOUBLE,0,MPIrank,MPI_COMM_WORLD);

         MPI::COMM_WORLD.Barrier();

         //receive
         if(MPIrank == 0)
            for(int rank = 1;rank < MPIsize;++rank)
               MPI_Recv(&Noffset[rank],1,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

         if(MPIrank == 0){
            for(int count = 0;count < MPIsize;++count){

               double frac = fabs( Noffset[count] ) / ( NDesiredWalkersPerRank[count] * fractionScaling ); //and the percentage of offset

               cout << "Deviation of average occupation on rank\t" << count << " = " << Noffset[count] << "\tFractional offset:\t" << frac << endl;

               if(frac > threshold_stop) 
                  oneFracDeviating = true; //if one or more offsets are too large: redistribute       
            }

            cout << endl;
         }

         //broadcast!
         MPI_Bcast(Noffset,MPIsize,MPI_DOUBLE,0,MPI_COMM_WORLD);
         MPI_Bcast(&oneFracDeviating,1,MPI::BOOL,0,MPI_COMM_WORLD);

      }

      delete [] work;

   }

   delete [] Noffset;

#endif

}
