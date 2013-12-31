#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "TrotterHeisenberg.h"
#include "OpSz.h"
#include "OpSx.h"
#include "OpSy.h"
#include "Lapack.h"

using namespace std;
const bool TrotterHeisenberg::debugPrint;

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 14, 2013 */

/**
 * constructor of the TrotterHeisenberg object. This contains the information about the trotter decomposition of one- and two-body terms.
 * @param theMPO MPO object containing the information about the system
 * @param dtau timestep
 */
TrotterHeisenberg::TrotterHeisenberg(HeisenbergMPO *theMPO, const double dtau){

   this->length = theMPO->gLength();
   this->phys_d = theMPO->gPhys_d();
   this->dtau = dtau;

   Sx = new OpSx(phys_d);
   Sy = new OpSy(phys_d);
   Sz = new OpSz(phys_d);

   //The couplings, and the different couplings
   Jv = new double [length * length];

   for(int i = 0;i < length;i++){

      Jv[ i*length + i ] = 0.0;

      for(int j = i + 1;j < length;j++){

         double theCoupling = theMPO->gCoupling(i,j);

         //couplingmatrix will contain the coupling J_[ij]
         Jv[ j*length + i ] = theCoupling;
         Jv[ i*length + j ] = theCoupling;

      }

   }

   //diagonalize the couplingmatrix
   char jobz = 'V';
   char uplo = 'U';

   Jeig = new double [length];

   int lwork = 3*length - 1;
   double * work = new double [lwork];
   int info;

   dsyev_(&jobz,&uplo,&length,Jv,&length,Jeig,work,&lwork,&info);

   delete [] work;

   //The magnetic field and single-site propagator
   this->theField = theMPO->gField();

   this->isMagneticField = (fabs(theField) < 1.0e-15) ? false : true;

   //make the single site propagators
   SzProp = new complex<double> [phys_d*phys_d];

   //e^Sz is diagonal in the physical site-basis
   for(int i = 0;i < phys_d;++i)
      for(int j = 0;j < phys_d;++j)
         SzProp[j*phys_d + i] = exp((*Sz)(i,j));

   //for e^Sx and e^Sy we need the eigenvalues and eigenvectors of Sx and Sy to construct the propagator
   complex<double> *mem = new complex<double> [phys_d*phys_d];
   double *mem_eig = new double [phys_d];

   int clwork = 3*phys_d - 1;

   complex<double> *cwork = new complex<double> [clwork];

   int rlwork = 3*phys_d - 2;

   double *rwork = new double [rlwork];

   for(int i = 0;i < phys_d;++i)
      for(int j = 0;j < phys_d;++j)
         mem[j*phys_d + i] = (*Sx)(i,j);

   zheev_(&jobz,&uplo,&phys_d,mem,&phys_d,mem_eig,cwork,&clwork,rwork,&info);

   //now construct the exponential:
   SxProp = new complex<double> [phys_d*phys_d];

   for(int i = 0;i < phys_d;++i)
      for(int j = 0;j < phys_d;++j){

         SxProp[j*phys_d + i] = complex<double>(0.0,0.0);

         for(int k = 0;k < phys_d;++k)
            SxProp[j*phys_d + i] += exp(mem_eig[k]) * mem[k*phys_d + i] * std::conj(mem[k*phys_d + j]);

      }

   //Finally Sy
   for(int i = 0;i < phys_d;++i)
      for(int j = 0;j < phys_d;++j)
         mem[j*phys_d + i] = (*Sy)(i,j);

   zheev_(&jobz,&uplo,&phys_d,mem,&phys_d,mem_eig,cwork,&clwork,rwork,&info);

   //construct exponential
   SyProp = new complex<double> [phys_d*phys_d];

   for(int i = 0;i < phys_d;++i)
      for(int j = 0;j < phys_d;++j){

         SyProp[j*phys_d + i] = complex<double>(0.0,0.0);

         for(int k = 0;k < phys_d;++k)
            SyProp[j*phys_d + i] += exp(mem_eig[k]) * mem[k*phys_d + i] * std::conj(mem[k*phys_d + j]);

      }

   delete [] mem;
   delete [] mem_eig;
   delete [] cwork;
   delete [] rwork;

}

TrotterHeisenberg::~TrotterHeisenberg(){

   delete Sx;
   delete Sy;
   delete Sz;

   delete [] Jv;
   delete [] Jeig;

   delete [] SxProp;
   delete [] SyProp;
   delete [] SzProp;

}

bool TrotterHeisenberg::gIsMagneticField() const {

   return isMagneticField; 

}

double TrotterHeisenberg::gField() const {

   return theField; 

}

/**
 * @return the j'th index of the k'th eigenvector of the couplingMatrix
 */
double TrotterHeisenberg::gJv(const int k, const int i) const {

   return Jv[ k*length + i ]; 

}

/**
 * @return the i'th eigenvalue of the couplingmatrix
 */
double TrotterHeisenberg::gJeig( const int i) const {

   return Jeig[i]; 

}


/**
 * @return element (i,j) of the single site propagotor e^{tau h S_z)
 */
complex<double> TrotterHeisenberg::gSzProp(const int i, const int j) const{

   return SzProp[j*phys_d + i];

}
