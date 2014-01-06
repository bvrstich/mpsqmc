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

   int L = 20;
   int DT = 4;
   int DW = 4;
   int d = 2;

   HeisenbergMPO theMPO(L,d);

   for (int cnt = 0;cnt < L-1;cnt++)
      theMPO.sCoupling(cnt,cnt+1,1.0);

   theMPO.sField(0.0);

   Random RN;
   MPSstate Psi0("debug.mps",&RN);

   int Nwalkers = 1000;
   double dtau = 0.01;
   int nSteps = 1000;

   AFQMC thePopulation(&theMPO, &RN,&Psi0,DW, Nwalkers, dtau);
   thePopulation.Walk(nSteps);

#ifdef USE_MPI_IN_MPSQMC
   MPI::Finalize();
#endif

   return 0;

}
