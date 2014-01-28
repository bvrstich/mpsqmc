#include "MPItag.h" //Here USE_MPI_IN_MPSQMC can be defined to switch to non-MPI compilers.

#include <iostream>

#include "Random.h"
#include "MPStensor.h"
#include "HeisenbergMPO.h"
#include "TwoSiteObject.h"
#include "MPSstate.h"
#include "AFMPO.h"
#include "TrotterHeisenberg.h"
#include "Lapack.h"
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
/*
   int L = 20;
   int DT = 32;
   int DW = 2;
   int d = 2;

   HeisenbergMPO theMPO(L,d);

   //set2DHeis(sqrt(L),1.0,0.0,theMPO);
   //set1DHeis(L,1.0,theMPO);

   ifstream input("J.in");

   int I,J;
   double value;

   while(input >> I >> J >> value){

      theMPO.sCoupling(I,J,value);
      theMPO.sCoupling(J,I,value);

   }

   theMPO.sField(0.0);

   Random RN;

   //MPSstate Psi0("input/Heisenberg2D/L4DT4.mps",&RN);
   //MPSstate Psi0("input/Heisenberg1D/L50D4.mps",&RN);
   MPSstate Psi0("debug_D32.mps",&RN);

   int Nwalkers = 1000;
   double dtau = 0.01;
   int nSteps = 100000;

   AFQMC thePopulation(&theMPO, &RN,&Psi0,DW, Nwalkers, dtau);
   thePopulation.Walk(nSteps);
*/
#ifdef USE_MPI_IN_MPSQMC
   MPI::Finalize();
#endif

   Random RN;

   int L = 20;
   int DTA = 64;
   int DTB = 8;
   int d = 2;

   MPSstate A(L,DTA,d,&RN);
   MPSstate B(L,DTB,d,&RN);

   HeisenbergMPO theMPO(L,d);

   ifstream input("J.in");

   int I,J;
   double value;

   while(input >> I >> J >> value){

      theMPO.sCoupling(I,J,value);
      theMPO.sCoupling(J,I,value);

   }

   theMPO.sField(0.0);

   int dimRA = A.gDimAtBound(1);
   int dimRB = B.gDimAtBound(1);
   int dimRmpo = theMPO.dimR(0);

   int DO = theMPO.gDtrunc();

   complex<double> *tmp2 = new complex<double> [d*DO*DTA*DTB];
   complex<double> *tmp1 = new complex<double> [d*DO*DTB*DTA];
   complex<double> *tmp = new complex<double> [DO*DTB*DTA];

   for(int i = 0;i < dimRA;++i)
      for(int o = 0;o < dimRmpo;++o)
         for(int s = 0;s < d;++s){

            tmp2[s*dimRmpo*dimRA + o*dimRA + i] = complex<double>(0.0,0.0);

            for(int s_ = 0;s_ < d;++s_)
               tmp2[s*dimRmpo*dimRA + o*dimRA + i] += theMPO(0,0,s_,s,o) * A[0](s_,0,i);

         }

   for(int o = 0;o < dimRmpo;++o)
      for(int i = 0;i < dimRA;++i)
         for(int j = 0;j < dimRB;++j){

            tmp[o*dimRA*dimRB + j*dimRA + i] = complex<double>(0.0,0.0);

            for(int s = 0;s < d;++s)
               tmp[o*dimRA*dimRB + j*dimRA + i] += tmp2[s*dimRmpo*dimRA + o*dimRA + i] * std::conj(B[0](s,0,j));

         }

   for(int site = 1;site < L;++site){

      int dimRA = A.gDimAtBound(site + 1);
      int dimLA = A.gDimAtBound(site);

      int dimRB = B.gDimAtBound(site + 1);
      int dimLB = B.gDimAtBound(site);

      int dimLmpo = theMPO.dimL(site);
      int dimRmpo = theMPO.dimR(site);

      //top
      for(int i = 0;i < dimRA;++i)
         for(int j = 0;j < dimLB;++j)
            for(int o = 0;o < dimLmpo;++o)
               for(int s = 0;s < d;++s){

                  tmp1[s*dimLmpo*dimRA*dimLB + o*dimRA*dimLB + j*dimRA + i] = complex<double>(0.0,0.0);

                  for(int k = 0;k < dimLA;++k)
                     tmp1[s*dimLmpo*dimRA*dimLB + o*dimRA*dimLB + j*dimRA + i] += tmp[o*dimLA*dimLB + j*dimLA + k] * A[site](s,k,i);

               }

      //operator
      for(int i = 0;i < dimRA;++i)
         for(int j = 0;j < dimLB;++j)
            for(int o = 0;o < dimRmpo;++o)
               for(int s = 0;s < d;++s){

                  tmp2[s*dimRmpo*dimRA*dimLB + o*dimRA*dimLB + j*dimRA + i] = complex<double>(0.0,0.0);

                  for(int p = 0;p < dimLmpo;++p)
                     for(int s_ = 0;s_ < d;++s_)
                        tmp2[s*dimRmpo*dimRA*dimLB + o*dimRA*dimLB + j*dimRA + i] += tmp1[s_*dimLmpo*dimRA*dimLB + p*dimRA*dimLB + j*dimRA + i] * theMPO(site,p,s_,s,o);

               }

      //bottom
      for(int i = 0;i < dimRA;++i)
         for(int j = 0;j < dimRB;++j)
            for(int o = 0;o < dimRmpo;++o){

               tmp[o*dimRA*dimRB + j*dimRA + i] = complex<double>(0.0,0.0);

               for(int s = 0;s < d;++s)
                  for(int k = 0;k < dimLB;++k)
                     tmp[o*dimRA*dimRB + j*dimRA + i] += tmp2[s*dimRmpo*dimRA*dimLB + o*dimRA*dimLB + k*dimRA + i] * std::conj(B[site](s,k,j));

            }

   }

   cout << tmp[0] << endl;

   MPSstate OA(L,DTA,d,&RN);

   OA.ApplyMPO(false,&theMPO,&A);
   cout << OA.InnerProduct(&B) << endl;

   cout << A.expectation(&theMPO,&B) << endl;

   delete [] tmp;
   delete [] tmp1;
   delete [] tmp2;

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

   for (int cnt = 0;cnt < L-1;cnt++){

      theMPO.sCoupling(cnt,cnt+1,J);
      theMPO.sCoupling(cnt+1,cnt,J);

   }

   theMPO.sField(0.0);

}
