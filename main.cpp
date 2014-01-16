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

void set2DHeis(int L,double J1,double J2,HeisenbergMPO &theMPO);
void set1DHeis(int L,double J,HeisenbergMPO &theMPO);

int main(int argc,char *argv[]){

#ifdef USE_MPI_IN_MPSQMC
   MPI::Init();
#endif

   cout.precision(15);

   int L = 16;
   int DT = 4;
   int DW = 2;
   int d = 2;

   HeisenbergMPO theMPO(L,d);

   set2DHeis(sqrt(L),1.0,0.0,theMPO);
   //set1DHeis(L,1.0,theMPO);

   Random RN;
   //MPSstate Psi0("1D_Heis_L=10_DT=4.mps",&RN);
   MPSstate Psi0("debug_2D.mps",&RN);
  
   int Nwalkers = 10000;
   double dtau = 0.001;
   int nSteps = 100000;

   AFQMC thePopulation(&theMPO, &RN,&Psi0,DW, Nwalkers, dtau);
   //thePopulation.testProp();
   thePopulation.Walk(nSteps);
/*
   MPSstate HPsi0(L,DT,d,&RN);
   MPSstate H2Psi0(L,DT,d,&RN);

   HPsi0.ApplyMPO(false,&theMPO,&Psi0);
   H2Psi0.ApplyMPO(false,&theMPO,&HPsi0);

   cout << HPsi0.InnerProduct(&Psi0) << endl;

   cout << 1.0 - dtau * HPsi0.InnerProduct(&Psi0) + 0.5 * dtau*dtau * H2Psi0.InnerProduct(&Psi0) << endl;
*/
#ifdef USE_MPI_IN_MPSQMC
   MPI::Finalize();
#endif

   return 0;

}

void set2DHeis(int L,double J1,double J2,HeisenbergMPO &theMPO){

   for (int row=0; row<L; row++){
      for (int col=0; col<L; col++){
         int number = row + L * col;
         int neighbour1 = (row + 1       )%L + L * col;
         int neighbour2 = (row - 1 + L)%L + L * col;
         int neighbour3 = row                   + L * ((col - 1 + L)%L);
         int neighbour4 = row                   + L * ((col + 1       )%L);
         theMPO.sCoupling(number, neighbour1, J1);
         theMPO.sCoupling(number, neighbour2, J1);
         theMPO.sCoupling(number, neighbour3, J1);
         theMPO.sCoupling(number, neighbour4, J1);
         int neighbour5 = (row + 1       )%L + L * ((col + 1       )%L);
         int neighbour6 = (row + 1       )%L + L * ((col - 1 + L)%L);
         int neighbour7 = (row - 1 + L)%L + L * ((col - 1 + L)%L);
         int neighbour8 = (row - 1 + L)%L + L * ((col + 1       )%L);
         theMPO.sCoupling(number, neighbour5, J2);
         theMPO.sCoupling(number, neighbour6, J2);
         theMPO.sCoupling(number, neighbour7, J2);
         theMPO.sCoupling(number, neighbour8, J2);
      }
   }

   theMPO.sField(0.0);

}

void set1DHeis(int L,double J,HeisenbergMPO &theMPO){

   for (int cnt = 0;cnt < L-1;cnt++)
      theMPO.sCoupling(cnt,cnt+1,J);

   theMPO.sField(0.0);

}
