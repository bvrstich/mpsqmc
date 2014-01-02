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
   J = new complex<double> [length * length];

   for(int i = 0;i < length;i++){

      J[ i*length + i ] = 0.0;

      for(int j = i + 1;j < length;j++){

         complex<double> theCoupling = theMPO->gCoupling(i,j);

         //couplingmatrix will contain the coupling J_[ij]
         J[ j*length + i ] = theCoupling;
         J[ i*length + j ] = conj(theCoupling);

      }

   }

   //diagonalize the couplingmatrix
   char jobz = 'V';
   char uplo = 'U';

   Jeig = new double [length];

   int lwork = 3*length - 1;
   complex<double> * work = new complex<double> [lwork];

   int rlwork = 3*length - 2;
   double *rwork = new double [rlwork];
  
   int info;

   zheev_(&jobz,&uplo,&length,J,&length,Jeig,work,&lwork,rwork,&info);

   delete [] work;
   delete [] rwork;

   V = new complex<double> [length*length];

   //now transform the J elements with dtau and the eigenvalues for the propagator to form the transformation V
   for(int k = 0;k < length;++k){

      complex<double> tmp =  std::sqrt( (complex<double>)-2.0 * Jeig[k] * dtau);

      for(int i = 0;i < length;++i)
         V[k*length + i] = tmp * J[k*length + i];

   }

   //The magnetic field and single-site propagator
   this->theField = theMPO->gField();

   this->isMagneticField = (abs(theField) < 1.0e-15) ? false : true;

   //make the single site propagators
   H1Prop = new complex<double> [phys_d*phys_d];

   //e^Sz is diagonal in the physical site-basis
   for(int i = 0;i < phys_d;++i)
      for(int j = 0;j < phys_d;++j)
         H1Prop[j*phys_d + i] = exp(-0.5 * dtau * theField * (*Sz)(i,j));

   Sx_vec = new complex<double> [phys_d*phys_d];
   Sy_vec = new complex<double> [phys_d*phys_d];

   //for e^Sx and e^Sy we need the eigenvalues and eigenvectors of Sx and Sy to construct the propagator
   eig = new double [phys_d];

   lwork = 3*phys_d - 1;
   work = new complex<double> [lwork];

   rlwork = 3*phys_d - 2;
   rwork = new double [rlwork];

   for(int i = 0;i < phys_d;++i)
      for(int j = 0;j < phys_d;++j)
         Sx_vec[j*phys_d + i] = (*Sx)(i,j);

   zheev_(&jobz,&uplo,&phys_d,Sx_vec,&phys_d,eig,work,&lwork,rwork,&info);

   //Finally Sy
   for(int i = 0;i < phys_d;++i)
      for(int j = 0;j < phys_d;++j)
         Sy_vec[j*phys_d + i] = (*Sy)(i,j);

   zheev_(&jobz,&uplo,&phys_d,Sy_vec,&phys_d,eig,work,&lwork,rwork,&info);

   delete [] work;
   delete [] rwork;

   AFProp = new complex<double> [length * phys_d * phys_d];

   //finally construct the auxiliary field operators as MPO's:
   V_Op = new AFMPO * [3*length];

   for(int r = 0;r < 3;++r)
      for(int i = 0;i < length;++i)
         V_Op[r*length + i] = new AFMPO(length,phys_d,r,V + i*length);

}

TrotterHeisenberg::~TrotterHeisenberg(){

   delete Sx;
   delete Sy;
   delete Sz;

   delete [] J;
   delete [] Jeig;

   delete [] V;

   delete [] H1Prop;

   delete [] eig;
   delete [] Sx_vec;
   delete [] Sy_vec;

   delete [] AFProp;

   for(int r = 0;r < 3;++r)
      for(int i = 0;i < length;++i)
         delete V_Op[r*length + i];

   delete [] V_Op;

}

bool TrotterHeisenberg::gIsMagneticField() const {

   return isMagneticField; 

}

complex<double> TrotterHeisenberg::gField() const {

   return theField; 

}

/**
 * @return the j'th index of the k'th eigenvector of the couplingMatrix
 */
complex<double> TrotterHeisenberg::gJ(const int k, const int i) const {

   return J[ k*length + i ]; 

}

/**
 * @return the j'th index of the k'th eigenvector of the transformation matrix
 */
complex<double> TrotterHeisenberg::gV(const int k, const int i) const {

   return V[ k*length + i ]; 

}

/**
 * @return the i'th eigenvalue of the couplingmatrix
 */
double TrotterHeisenberg::gJeig( const int i) const {

   return Jeig[i]; 

}

/**
 * @return the timestep
 */
double TrotterHeisenberg::gtau() const {

   return dtau;

}

/**
 * @return the propagator matrix for the H1
 */
complex<double> TrotterHeisenberg::gH1Prop(int i,int j) const {

   return H1Prop[j*phys_d + i];

}

/**
 * @return the propagator matrix for the H1
 */
complex<double> TrotterHeisenberg::gAFProp(int site,int i,int j) const {

   return AFProp[site*phys_d*phys_d + j*phys_d + i];

}

/**
 * fill the AFProp array with the correct propagator
 * @param k index of eigenvectors of J
 * @param r type of operator, x,y or z
 * @param x auxiliary field variable
 */
void TrotterHeisenberg::fillAFProp(int k,int r,double x){

   if(r == 0){//x

      for(int site = 0;site < length;++site){

         for(int i = 0;i < phys_d;++i)
            for(int j = 0;j < phys_d;++j){

               AFProp[site*phys_d*phys_d + j*phys_d + i] = complex<double>(0.0,0.0);

               for(int l = 0;l < phys_d;++l)//loop over eigenvector of Sx--> l
                  AFProp[site*phys_d*phys_d + j*phys_d + i] += exp(x * V[k*length + site] * eig[l]) * Sx_vec[l*phys_d + i] * std::conj(Sx_vec[l*phys_d + j]);

            }

      }

   }
   else if(r == 1){//y

      for(int site = 0;site < length;++site){

         for(int i = 0;i < phys_d;++i)
            for(int j = 0;j < phys_d;++j){

               AFProp[site*phys_d*phys_d + j*phys_d + i] = complex<double>(0.0,0.0);

               for(int l = 0;l < phys_d;++l)//loop over eigenvector of Sx--> l
                  AFProp[site*phys_d*phys_d + j*phys_d + i] += exp(x * V[k*length + site] * eig[l]) * Sy_vec[l*phys_d + i] * std::conj(Sy_vec[l*phys_d + j]);

            }

      }

   }
   else{//z

      for(int site = 0;site < length;++site)
         for(int i = 0;i < phys_d;++i)
            for(int j = 0;j < phys_d;++j)
               AFProp[site*phys_d*phys_d + j*phys_d + i] = exp(x * V[k*length + site] * (*Sz)(i,j) );

   }

}

/**
 * @return the auxiliary field operator V associated with the k'th eigenvector of J, and operator Sr
 */
AFMPO *TrotterHeisenberg::gV_Op(int k,int r){

   return V_Op[r*length + k];

}
