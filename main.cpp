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

   Random RN;

   //MPSstate Psi0("debug_2D.mps",&RN);
   
   MPSstate Psi0(L*L,DT,d,&RN);
   MPSstate HPsi0(L*L,DT,d,&RN);

   HPsi0.ApplyMPO(false,&theMPO,&Psi0);

   int DO = theMPO.gDtrunc();

   complex<double> *tmp = new complex<double> [DT*DT*DO];
   complex<double> *tmp1 = new complex<double> [DT*DT*DO*d];
   complex<double> *tmp2 = new complex<double> [DT*DT*DO*d];

   int dimR = Psi0.gDimAtBound(1);
   int dimRmpo = theMPO.dimR(0);

   for(int i = 0;i < dimR;++i)
      for(int o = 0;o < dimRmpo;++o)
         for(int s = 0;s < d;++s){

            tmp2[i*DO*d + o*d + s] = complex<double>(0.0,0.0);

            for(int s_ = 0;s_ < d;++s_)
               tmp2[i*DO*d + o*d + s] += Psi0[0](s_,0,i) * theMPO(0,0,s_,s,o);

   }

   for(int s = 0;s < d;++s)
      for(int o = 0;o < dimRmpo;++o)
         for(int i = 0;i < dimR;++i)
            cout << s << "\t" << o << "\t" << i << "\t" << tmp2[i*DO*d + o*d + s] << endl;

   for(int o = 0;o < dimRmpo;++o)
      for(int i = 0;i < dimR;++i)
         for(int j = 0;j < dimR;++j){

            tmp[o*DT*DT + j*DT + i] = complex<double>(0.0,0.0);

            for(int s = 0;s < d;++s)
               tmp[o*DT*DT + j*DT + i] += tmp2[i*DO*d + o*d + s] * std::conj(Psi0[0](s,0,j));

         }

   for(int site = 1;site < L*L;++site){

      int dimR = Psi0.gDimAtBound(site + 1);
      int dimL = Psi0.gDimAtBound(site);

      int dimLmpo = theMPO.dimL(site);
      int dimRmpo = theMPO.dimR(site);

      //top
      for(int i = 0;i < dimR;++i)
         for(int j = 0;j < dimL;++j)
            for(int o = 0;o < dimLmpo;++o)
               for(int s = 0;s < d;++s){

                  tmp1[s*DO*DT*DT + o*DT*DT + j*DT + i] = complex<double>(0.0,0.0);

                  for(int k = 0;k < dimL;++k)
                     tmp1[s*DO*DT*DT + o*DT*DT + j*DT + i] += tmp[o*DT*DT + j*DT + k] * Psi0[site](s,k,i);

               }

      //operator
      for(int i = 0;i < dimR;++i)
         for(int j = 0;j < dimL;++j)
            for(int o = 0;o < dimRmpo;++o)
               for(int s = 0;s < d;++s){

                  tmp2[s*DO*DT*DT + o*DT*DT + j*DT + i] = complex<double>(0.0,0.0);

                  for(int p = 0;p < dimLmpo;++p)
                     for(int s_ = 0;s_ < d;++s_)
                        tmp2[s*DO*DT*DT + o*DT*DT + j*DT + i] += tmp1[s_*DO*DT*DT + p*DT*DT + j*DT + i] * theMPO(site,p,s_,s,o);

               }

      //bottom
      for(int i = 0;i < dimR;++i)
         for(int j = 0;j < dimR;++j)
            for(int o = 0;o < dimRmpo;++o){

               tmp[o*DT*DT + j*DT + i] = complex<double>(0.0,0.0);

               for(int s = 0;s < d;++s)
                  for(int k = 0;k < dimL;++k)
                     tmp[o*DT*DT + j*DT + i] += tmp2[s*DO*DT*DT + o*DT*DT + k*DT + i] * std::conj(Psi0[site](s,k,j));

               }

   }

   cout << tmp[0] << endl;

   delete [] tmp;
   delete [] tmp1;
   delete [] tmp2;

   cout << HPsi0.InnerProduct(&Psi0) << endl;
   cout << Psi0.expectation(&theMPO,&Psi0) << endl;
   /*
      int Nwalkers = 1000;
      double dtau = 0.01;
      int nSteps = 100000;

      AFQMC thePopulation(&theMPO, &RN,&Psi0,DW, Nwalkers, dtau);
      thePopulation.Walk(nSteps);
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
