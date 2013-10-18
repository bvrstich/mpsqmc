#include "MPItag.h" //Here USE_MPI_IN_MPSQMC can be defined to switch to non-MPI compilers.

#include <iostream>

#include "HeisenbergMPO.h"
#include "MPSstate.h"
#include "DMRG.h"
#include "MPSQMC2.h"
#include "Random.h"
#include "TrotterHeisenberg.h"
#include "GridGenerator.h"

using namespace std;

void MPOcheck();
void DMRGcheck();
void RNcheck();
void MPSQMCcheck();
void HeisenbergSquareLattice();

int main(void){

   #ifdef USE_MPI_IN_MPSQMC
      MPI::Init();
   #endif
   
   cout.precision(15);
   
   //MPOcheck();
   //DMRGcheck();
   //RNcheck();
   //MPSQMCcheck();
   HeisenbergSquareLattice();
   
   #ifdef USE_MPI_IN_MPSQMC
      MPI::Finalize();
   #endif
   
   return 0;

}

void HeisenbergSquareLattice(){

   #ifdef USE_MPI_IN_MPSQMC
      int rank = MPI::COMM_WORLD.Get_rank();
   #else
      int rank = 0;
   #endif

   /* The MPO */
   int base = 4;
   int length = base*base;
   int d = 2;
   bool useLadder = false;
   HeisenbergMPO theMPO(length,d,useLadder);
   for (int row=0; row<base; row++){
      for (int col=0; col<base; col++){
         int number = row + base * col;
         int neighbour1 = (row + 1       )%base + base * col;
         int neighbour2 = (row - 1 + base)%base + base * col;
         int neighbour3 = row                   + base * ((col - 1 + base)%base);
         int neighbour4 = row                   + base * ((col + 1       )%base);
         theMPO.sCoupling(number, neighbour1, 1.0);
         theMPO.sCoupling(number, neighbour2, 1.0);
         theMPO.sCoupling(number, neighbour3, 1.0);
         theMPO.sCoupling(number, neighbour4, 1.0);
      }
   }
   theMPO.sField(0.0);
   theMPO.findNonZeroContributions();
   
   Random RN;
   
   /* DMRG calculations */
   const bool doDMRG = true;
   if ((rank==0) && (doDMRG)){
      int D = 10;
      MPSstate Psi0(length, D, d, &RN);
      DMRG theSolver(&Psi0, &theMPO);
      double Energy = theSolver.Solve();
      cout << "The energy from DMRG = " << Energy << endl; //J=1 square 4x4, h=0, d=2 E("FCI") = -11.2284832084289
   }
   
   /*GridGenerator theGrid(4);
   theGrid.FillMarsaglia(4);*/
   
   /* The MPSQMC */
   /*int Dtrunc = 2;
   int Nwalkers = 1000;
   double dtau = 0.01;
   int nSteps = 10000;
   MPSQMC2 thePopulation(&theMPO, &theGrid, &RN, Dtrunc, Nwalkers, dtau);
   thePopulation.Walk(nSteps);*/

}

void MPSQMCcheck(){

   int length = 10;
   int d = 2;
   bool useLadder = false;
   HeisenbergMPO theMPO(length,d,useLadder);
   for (int cnt=0; cnt<length-1; cnt++){ theMPO.sCoupling(cnt,cnt+1,1.0); }
   theMPO.sField(0.0);
   theMPO.findNonZeroContributions();
   
   Random RN;
   
   if (false){
      int D = 2*2*2*2*2;
      MPSstate Psi0(length, D, d, &RN);
      DMRG theSolver(&Psi0, &theMPO);
      double Energy = theSolver.Solve();
      cout << "The energy from DMRG = " << Energy << endl; //E("FCI") = -4.25803520728288
   }
   
   GridGenerator theGrid(4);
   theGrid.FillSimple(4);

   int Dtrunc = 2;
   int Nwalkers = 1000;
   double dtau = 0.01;
   int nSteps = 10000;
   MPSQMC2 thePopulation(&theMPO, &theGrid, &RN, Dtrunc, Nwalkers, dtau);
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
      bool useLadder = false;
      HeisenbergMPO theMPO(length,d,useLadder);
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
   
      for (int geval=0; geval<2; geval++){
         int length = 6;
         int d = 6;
         bool useLadder = (geval==0) ? true : false;
         HeisenbergMPO test(length,d,useLadder);
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
      }
      
      /*****/
      
      for (int geval=0; geval<2; geval++){
         int length = 7;
         int d = 3;
         bool useLadder = (geval==0) ? true : false;
         HeisenbergMPO test(length,d,useLadder);
         double value = 1.0;
         for (int cnt=0; cnt<length-1; cnt++){
            for (int cnt2=cnt+1; cnt2<length; cnt2++){
               value += 1.0;
               test.sCoupling(cnt,cnt2,value);
            }
         }
         test.sField(-0.0);
         test.findNonZeroContributions();
         cout << test;
      }
      
   }

}


