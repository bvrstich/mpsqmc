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

   int L = 4;
   int DT = 4;
   int DW = 4;
   int d = 2;

   HeisenbergMPO theMPO(L*L,d);

   set2DHeis(L,1.0,0.0,theMPO);

   TrotterHeisenberg theTrotter(&theMPO,0.01);

   Random RN;

   MPSstate A(L*L,DT,d,&RN);
   MPSstate B(L*L,DT,d,&RN);

   MPSstate tmp(L*L,DT,d,&RN);
   MPSstate tmp2(L*L,DT,d,&RN);

   tmp.ApplyMPO(false,&theMPO,&A);

   cout << tmp.InnerProduct(&B) << endl;

   complex<double> val(0.0,0.0);

   for(int r = 0;r < 3;++r)
      for(int k = 0;k < theTrotter.gn_trot();++k){

         tmp.ApplyMPO(false,theTrotter.gV_Op(k,r),&A);
         tmp2.ApplyMPO(false,theTrotter.gV_Op(k,r),&tmp);

         val += tmp2.InnerProduct(&B);

      }

   cout << -0.5 * val / 0.01 << endl;

/*
   MPSstate Psi0("debug_2D.mps",&RN);

   int Nwalkers = 1000;
   double dtau = 0.01;
   int nSteps = 100000;

   AFQMC thePopulation(&theMPO, &RN,&Psi0,DW, Nwalkers, dtau);
   thePopulation.Walk(1);
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
