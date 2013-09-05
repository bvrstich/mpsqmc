#include <iostream>

#include "MPItag.h" //Here USE_MPI_IN_MPSQMC can be defined to switch to non-MPI compilers.
#include "HeisenbergMPO.h"
#include "MPSstate.h"
#include "DMRG.h"
#include "MPSQMC.h"
#include "Random.h"

using namespace std;

void MPOcheck();
void DMRGcheck();
void RNcheck();
void MPSQMCcheck();

int main(void){

   #ifdef USE_MPI_IN_MPSQMC
      MPI::Init();
   #endif
   
   cout.precision(15);
   
   //MPOcheck();
   //DMRGcheck();
   //RNcheck();
   MPSQMCcheck();
   
   #ifdef USE_MPI_IN_MPSQMC
      MPI::Finalize();
   #endif
   
   return 0;

}

void MPSQMCcheck(){

   //Should work both with and without MPI
   
   /* The MPO */
   int length = 10;
   int d = 3;
   HeisenbergMPO theMPO(length,d);
   for (int cnt=0; cnt<length-1; cnt++){ theMPO.sCoupling(cnt,cnt+1,1.0); }
   theMPO.sField(0.0);
   theMPO.findNonZeroContributions();
   
   /* The MPSQMC */
   int Dtrunc = 2;
   int Nwalkers = 100;
   double dtau = 0.01;
   int nSteps = 10;
   Random RN;
   MPSQMC thePopulation(&theMPO, &RN, Dtrunc, Nwalkers, dtau);
   thePopulation.Walk(nSteps);

}

void RNcheck(){

   //Should work both with and without MPI
   Random RN;
   RN.test();

}

void DMRGcheck(){

   #ifdef USE_MPI_IN_MPSQMC
      int rank = MPI::COMM_WORLD.Get_rank();
   #else
      int rank = 0;
   #endif

   if (rank==0){
      int length = 10;
      int d = 3;
      HeisenbergMPO theMPO(length,d);
      for (int cnt=0; cnt<length-1; cnt++){ theMPO.sCoupling(cnt,cnt+1,1.0); }
      theMPO.sField(0.0);
      theMPO.findNonZeroContributions();
      
      int D = 10;
      Random RN;
      MPSstate Psi0(length, D, d, &RN);
      DMRG theSolver(&Psi0, &theMPO);
      double Energy = theSolver.Solve();
      cout << "The energy from DMRG = " << Energy << endl; //n.n. J=1 ; h=0 ; L=10 ; d=3 ; E(FCI)=-12.8945601322109
      
      double Norm = Psi0.RightNormalize();
      MPSstate Psi1(&Psi0);
      Psi1.ApplyMPO(&theMPO,&Psi0);
      Psi1.CompressState();
      double Energy2 = Psi1.InnerProduct(&Psi0);
      cout << "The energy from <Psi0|MPO|Psi0> = " << Energy2 << endl;
   }

}

void MPOcheck(){

   #ifdef USE_MPI_IN_MPSQMC
      int rank = MPI::COMM_WORLD.Get_rank();
   #else
      int rank = 0;
   #endif

   if (rank==0){
      int length = 6;
      int d = 6;
      HeisenbergMPO test(length,d);
      double value = 1.0;
      for (int cnt=0; cnt<length-1; cnt++){
         for (int cnt2=cnt+1; cnt2<length; cnt2++){
            value += 1.0;
            test.sCoupling(cnt,cnt2,value);
         }
      }
      test.sField(1.0);
      test.findNonZeroContributions();
      cout << test;
      
      /*****/
      
      int length2 = 7;
      int d2 = 3;
      HeisenbergMPO test2(length2,d2);
      double value2 = 1.0;
      for (int cnt=0; cnt<length2-1; cnt++){
         for (int cnt2=cnt+1; cnt2<length2; cnt2++){
            value2 += 1.0;
            test2.sCoupling(cnt,cnt2,value2);
         }
      }
      test2.sField(-0.0);
      test2.findNonZeroContributions();
      cout << test2;
   }

}


