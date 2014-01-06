#include "MPItag.h" //Here USE_MPI_IN_MPSQMC can be defined to switch to non-MPI compilers.

#include <iostream>

#include "Random.h"
#include "MPStensor.h"
#include "HeisenbergMPO.h"
#include "TwoSiteObject.h"
#include "MPSstate.h"
#include "AFMPO.h"
#include "TrotterHeisenberg.h"
#include "Walker.h"
#include "AFQMC.h"

using namespace std;

void MPOcheck();
void DMRGcheck();
void RNcheck();
void MPSQMCcheck();
void J1J2SquareLattice(const int d, const int L, const double J1, const double J2);

int main(int argc,char *argv[]){

#ifdef USE_MPI_IN_MPSQMC
   MPI::Init();
#endif

   cout.precision(15);

   int L = 16;
   int DT = 4;
   int DW = 4;
   int d = 2;

   HeisenbergMPO theMPO(L,d);

   double J1 = 1.0;
   double J2 = 0.0;

   int base = 4;

   /* The MPO */
   int length = base*base;

   for (int row=0; row<base; row++){
      for (int col=0; col<base; col++){
         int number = row + base * col;
         int neighbour1 = (row + 1       )%base + base * col;
         int neighbour2 = (row - 1 + base)%base + base * col;
         int neighbour3 = row                   + base * ((col - 1 + base)%base);
         int neighbour4 = row                   + base * ((col + 1       )%base);
         theMPO.sCoupling(number, neighbour1, J1);
         theMPO.sCoupling(number, neighbour2, J1);
         theMPO.sCoupling(number, neighbour3, J1);
         theMPO.sCoupling(number, neighbour4, J1);
         int neighbour5 = (row + 1       )%base + base * ((col + 1       )%base);
         int neighbour6 = (row + 1       )%base + base * ((col - 1 + base)%base);
         int neighbour7 = (row - 1 + base)%base + base * ((col - 1 + base)%base);
         int neighbour8 = (row - 1 + base)%base + base * ((col + 1       )%base);
         theMPO.sCoupling(number, neighbour5, J2);
         theMPO.sCoupling(number, neighbour6, J2);
         theMPO.sCoupling(number, neighbour7, J2);
         theMPO.sCoupling(number, neighbour8, J2);
      }
   }

   theMPO.sField(0.0);

   Random RN;
   MPSstate Psi0("debug.mps",&RN);

   int Nwalkers = 1000;
   double dtau = 0.01;
   int nSteps = 1;

   AFQMC thePopulation(&theMPO, &RN,&Psi0,DW, Nwalkers, dtau);
   thePopulation.Walk(nSteps);

#ifdef USE_MPI_IN_MPSQMC
   MPI::Finalize();
#endif

   return 0;

}
