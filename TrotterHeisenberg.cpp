#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "omp.h"

#include "TrotterHeisenberg.h"
#include "Random.h"
#include "OpSz.h"
#include "OpSx.h"
#include "OpSy.h"
#include "Lapack.h"

using namespace std;

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

   n_trot = 0;

   for(int i = 0;i < length;++i){

      if(fabs(Jeig[i]) > 1.0e-14)
         n_trot++;

   }

   V = new complex<double> [n_trot*length];

   //now transform the J elements with dtau and the eigenvalues for the propagator to form the transformation V
   int k = 0;

   for(int i = 0;i < length;++i){

      if(fabs(Jeig[i]) > 1.0e-14){

         complex<double> tmp =  std::sqrt( (complex<double>)-Jeig[i] * dtau);

         for(int j = 0;j < length;++j)
            V[k*length + j] = tmp * J[i*length + j];

         ++k;

      }

   }

   //The magnetic field and single-site propagator
   this->theField = theMPO->gField();

   this->isMagneticField = (abs(theField) < 1.0e-15) ? false : true;

   //make the single site propagators
   H1Prop = new complex<double> [phys_d*phys_d];

   //e^Sz is diagonal in the physical site-basis
   for(int i = 0;i < phys_d;++i){

      H1Prop[i*phys_d + i] = exp(-0.5 * dtau * theField * (*Sz)(i,i));

      for(int j = i + 1;j < phys_d;++j)
         H1Prop[j*phys_d + i] = H1Prop[i*phys_d + j] = complex<double>(0.0,0.0);

   }

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

   //Find the maximum number of threads for this process
#ifdef _OPENMP
   const int myNOMPthreads = omp_get_max_threads();
#else
   const int myNOMPthreads = 1;
#endif

   AFProp = new complex<double> * [myNOMPthreads];

   for(int thr = 0;thr < myNOMPthreads;++thr)
      AFProp[thr] = new complex<double> [length * phys_d * phys_d];

   //finally construct the auxiliary field operators as MPO's:
   V_Op = new AFMPO * [3*n_trot];

   for(int r = 0;r < 3;++r)
      for(int k = 0;k < n_trot;++k)
         V_Op[r*n_trot + k] = new AFMPO(length,phys_d,r,V + k*length);

   Sv = new complex<double> * [myNOMPthreads];

   for(int thr = 0;thr < myNOMPthreads;++thr)
      Sv[thr] = new complex<double> [phys_d * phys_d];

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

#ifdef _OPENMP
   const int myNOMPthreads = omp_get_max_threads();
#else
   const int myNOMPthreads = 1;
#endif

   for(int thr = 0;thr < myNOMPthreads;++thr){

      delete [] AFProp[thr];
      delete [] Sv[thr];

   }

   delete [] AFProp;
   delete [] Sv;

   for(int r = 0;r < 3;++r)
      for(int i = 0;i < n_trot;++i)
         delete V_Op[r*n_trot + i];

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
complex<double> TrotterHeisenberg::gAFProp(int myID,int site,int i,int j) const {

   return AFProp[myID][site*phys_d*phys_d + j*phys_d + i];

}

/**
 * fill the AFProp array with the correct propagator
 * @param k index of eigenvectors of J
 * @param r type of operator, x,y or z
 * @param x shifted auxiliary field variable
 */
void TrotterHeisenberg::fillAFProp(int myID,int k,int r,complex<double> x){

   if(r == 0){//x

      for(int site = 0;site < length;++site){

         for(int i = 0;i < phys_d;++i)
            for(int j = 0;j < phys_d;++j){

               AFProp[myID][site*phys_d*phys_d + j*phys_d + i] = complex<double>(0.0,0.0);

               for(int l = 0;l < phys_d;++l)//loop over eigenvector of Sx--> l
                  AFProp[myID][site*phys_d*phys_d + j*phys_d + i] += exp(  x * V[k*length + site] * eig[l] ) * Sx_vec[l*phys_d + i] * std::conj(Sx_vec[l*phys_d + j]);

            }

      }

   }
   else if(r == 1){//y

      for(int site = 0;site < length;++site){

         for(int i = 0;i < phys_d;++i)
            for(int j = 0;j < phys_d;++j){

               AFProp[myID][site*phys_d*phys_d + j*phys_d + i] = complex<double>(0.0,0.0);

               for(int l = 0;l < phys_d;++l)//loop over eigenvector of Sy--> l
                  AFProp[myID][site*phys_d*phys_d + j*phys_d + i] += exp( x * V[k*length + site] * eig[l]) * Sy_vec[l*phys_d + i] * std::conj(Sy_vec[l*phys_d + j]);

            }

      }

   }
   else{//z

      //Sz is diagonal in the basis, so this is easy
      for(int site = 0;site < length;++site)
         for(int i = 0;i < phys_d;++i){

            AFProp[myID][site*phys_d*phys_d + i*phys_d + i] = exp(x * V[k*length + site] * (*Sz)(i,i) );

            for(int j = i + 1;j < phys_d;++j)
               AFProp[myID][site*phys_d*phys_d + j*phys_d + i] = AFProp[myID][site*phys_d*phys_d + i*phys_d + j] = complex<double>(0.0,0.0);

         }

   }

}

/**
 * fill the AFProp array with the correct propagator: pick a random direction
 * @param k index of eigenvectors of J
 * @param x shifted auxiliary field variable
 */
void TrotterHeisenberg::fillAFProp(int myID,int k,complex<double> x,Random *RN){

   const double PI = 3.14159265358979323;

   //draw random (theta,phi) 
   double theta = RN->rand() * PI;
   double phi = RN->rand() * 2.0 * PI;

   for(int i = 0;i < phys_d;++i)
      for(int j = 0;j < phys_d;++j)
         Sv[myID][phys_d*j + i] = cos(phi)*sin(theta) * (*Sx)(i,j) + sin(phi)*sin(theta) * (*Sy)(i,j) + cos(theta) * (*Sz)(i,j);

   //diagonalize the couplingmatrix
   char jobz = 'V';
   char uplo = 'U';

   int lwork = 3*phys_d - 1;
   complex<double> * work = new complex<double> [lwork];

   int rlwork = 3*phys_d - 2;
   double *rwork = new double [rlwork];

   int info;

   zheev_(&jobz,&uplo,&phys_d,Sv[myID],&phys_d,eig,work,&lwork,rwork,&info);

   for(int site = 0;site < length;++site){

      for(int i = 0;i < phys_d;++i)
         for(int j = 0;j < phys_d;++j){

            AFProp[myID][site*phys_d*phys_d + j*phys_d + i] = complex<double>(0.0,0.0);

            for(int l = 0;l < phys_d;++l)//loop over eigenvector of Sx--> l
               AFProp[myID][site*phys_d*phys_d + j*phys_d + i] += exp(  x * V[k*n_trot + site] * eig[l] ) * Sv[myID][l*phys_d + i] * std::conj(Sv[myID][l*phys_d + j]);

         }

   }

   delete [] work;
   delete [] rwork;

}

/**
 * @return the auxiliary field operator V associated with the k'th eigenvector of J, and operator Sr
 */
AFMPO *TrotterHeisenberg::gV_Op(int k,int r){

   return V_Op[r*n_trot + k];

}

/**
 * @return the number of trotter terms
 */
int TrotterHeisenberg::gn_trot() const {

   return n_trot;

}
