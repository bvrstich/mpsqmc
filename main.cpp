#include "MPItag.h" //Here USE_MPI_IN_MPSQMC can be defined to switch to non-MPI compilers.

#include <iostream>

#include "Random.h"
#include "MPStensor.h"
#include "TwoSiteObject.h"
#include "MPSstate.h"
#include "AFMPO.h"
#include "TrotterJ1J2.h"
#include "Lapack.h"
#include "Walker.h"
#include "AFQMC.h"
#include "WorkSpace.h"

using namespace std;

int main(int argc,char *argv[]){

#ifdef USE_MPI_IN_MPSQMC
   MPI::Init();
#endif

   cout.precision(15);

   int L = atoi(argv[1]); 
   int J2 = atoi(argv[2]);
   int DT = atoi(argv[3]);
   int DW = atoi(argv[4]);
   int d = 2;

   Random RN;

   char filename[100];

   if(J2 == 10)
      sprintf(filename,"input/J1J2/%dx%d/J2=1.0/Psi0/DT=%d.mps",L,L,DT);
   else
      sprintf(filename,"input/J1J2/%dx%d/J2=0.%d/Psi0/DT=%d.mps",L,L,J2,DT);

   MPSstate Psi0(filename,&RN);

   int Nwalkers = 1000;
   double dtau = 0.01;
   int nSteps = 100000;

   TrotterJ1J2 theTrotter(L,d,(double)0.1*J2,dtau);

   AFQMC::init(L*L,DT,d,&RN);

   AFQMC thePopulation(&theTrotter,&RN,&Psi0,DW, Nwalkers, dtau);
   thePopulation.Walk(nSteps);

   AFQMC::clear();

#ifdef USE_MPI_IN_MPSQMC
   MPI::Finalize();
#endif

   return 0;

}
