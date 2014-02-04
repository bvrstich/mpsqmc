#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <omp.h>
#include "TrotterJ1J2.h"
#include "Random.h"
#include "OpSz.h"
#include "OpSx.h"
#include "OpSy.h"
#include "Lapack.h"

using namespace std;

/**
 * constructor of the TrotterJ1J2 object. This contains the information about the trotter decomposition of one- and two-body terms.
 * @param L dimension of the lattice
 * @param d physical dimension
 * @param J2 coupling strength
 * @param dtau timestep
 */
TrotterJ1J2::TrotterJ1J2(bool pbc,int L,int d,double J2,const double dtau){

   this->length = L*L;
   this->phys_d = d;
   this->dtau = dtau;
   this->J2 = J2;

   Sx = new OpSx(phys_d);
   Sy = new OpSy(phys_d);
   Sz = new OpSz(phys_d);

   //The couplings, and the different couplings
   J = new complex<double> [length * length];

   for(int i = 0;i < length;i++)
      for(int j = 0;j < length;j++)
         J[j*length + i] = 0.0;

   if(pbc){

      for (int row=0; row< L; row++)
         for (int col=0; col< L; col++){

            int number = row + L * col;

            int neighbour1 = (row + 1       )%L + L * col;
            int neighbour2 = (row - 1 + L)%L + L * col;
            int neighbour3 = row                   + L * ((col - 1 + L)%L);
            int neighbour4 = row                   + L * ((col + 1       )%L);

            J[number*length + neighbour1] = 1.0;
            J[number*length + neighbour2] = 1.0;
            J[number*length + neighbour3] = 1.0;
            J[number*length + neighbour4] = 1.0;

            int neighbour5 = (row + 1       )%L + L * ((col + 1       )%L);
            int neighbour6 = (row + 1       )%L + L * ((col - 1 + L)%L);
            int neighbour7 = (row - 1 + L)%L + L * ((col - 1 + L)%L);
            int neighbour8 = (row - 1 + L)%L + L * ((col + 1       )%L);

            J[number*length + neighbour5] = J2;
            J[number*length + neighbour6] = J2;
            J[number*length + neighbour7] = J2;
            J[number*length + neighbour8] = J2;

         }

   }
   else{

      for (int row=1; row< L - 1; row++)
         for (int col=1; col < L -1; col++){

            int number = row + L * col;

            int neighbour1 = row + 1 + L * col;
            int neighbour2 = row - 1 + L * col;
            int neighbour3 = row                   + L * ( col - 1 );
            int neighbour4 = row                   + L * ( col + 1 );

            J[number*length+neighbour1] = 1.0;
            J[number*length+neighbour2] = 1.0;
            J[number*length+neighbour3] = 1.0;
            J[number*length+neighbour4] = 1.0;

            J[neighbour1*length+number] = 1.0;
            J[neighbour2*length+number] = 1.0;
            J[neighbour3*length+number] = 1.0;
            J[neighbour4*length+number] = 1.0;

            int neighbour5 = row + 1 + L * (col + 1 );
            int neighbour6 = row + 1 + L * (col - 1 );
            int neighbour7 = row - 1 + L * (col - 1 );
            int neighbour8 = row - 1 + L * (col + 1 );

            J[number*length+neighbour5] = J2;
            J[number*length+neighbour6] = J2;
            J[number*length+neighbour7] = J2;
            J[number*length+neighbour8] = J2;

            J[neighbour5*length+number] = J2;
            J[neighbour6*length+number] = J2;
            J[neighbour7*length+number] = J2;
            J[neighbour8*length+number] = J2;

         }

      //edges
      for (int col=0; col < L -1; col++){

         //row = 0
         int number = L * col;
         int neighbour = L * ( col + 1 );

         J[number*length+neighbour] = 1.0;
         J[neighbour*length+number] = 1.0;

         //row = L - 1
         number = L - 1 + L * col;
         neighbour = L - 1 + L * ( col + 1 );

         J[number*length+neighbour] = 1.0;
         J[neighbour*length+number] = 1.0;

      }

      for (int row=0; row < L -1; row++){

         //col = 0
         int number = row;
         int neighbour = row + 1;

         J[number*length+neighbour] = 1.0;
         J[neighbour*length+number] = 1.0;

         //col = L - 1
         number = row + L * (L - 1);
         neighbour = row + 1 + L * (L - 1);

         J[number*length+neighbour] = 1.0;
         J[neighbour*length+number] = 1.0;

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

TrotterJ1J2::~TrotterJ1J2(){

   delete Sx;
   delete Sy;
   delete Sz;

   delete [] J;
   delete [] Jeig;

   delete [] V;

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

/**
 * @return the j'th index of the k'th eigenvector of the couplingMatrix
 */
complex<double> TrotterJ1J2::gJ(const int k, const int i) const {

   return J[ k*length + i ]; 

}

/**
 * @return the j'th index of the k'th eigenvector of the transformation matrix
 */
complex<double> TrotterJ1J2::gV(const int k, const int i) const {

   return V[ k*length + i ]; 

}

/**
 * @return the i'th eigenvalue of the couplingmatrix
 */
double TrotterJ1J2::gJeig( const int i) const {

   return Jeig[i]; 

}

/**
 * @return the timestep
 */
double TrotterJ1J2::gtau() const {

   return dtau;

}

/**
 * @return the propagator matrix for the AF
 */
complex<double> TrotterJ1J2::gAFProp(int myID,int site,int i,int j) const {

   return AFProp[myID][site*phys_d*phys_d + j*phys_d + i];

}

/**
 * fill the AFProp array with the correct propagator
 * @param k index of eigenvectors of J
 * @param r type of operator, x,y or z
 * @param x shifted auxiliary field variable
 */
void TrotterJ1J2::fillAFProp(int myID,int k,int r,complex<double> x){

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
void TrotterJ1J2::fillAFProp(int myID,int k,complex<double> x,Random *RN){

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
AFMPO *TrotterJ1J2::gV_Op(int k,int r){

   return V_Op[r*n_trot + k];

}

/**
 * @return the number of trotter terms
 */
int TrotterJ1J2::gn_trot() const {

   return n_trot;

}

/**
 * @return the length of the chain
 */
int TrotterJ1J2::glength() const {

   return length;

}

/**
 * @return the J2 coupling parameter
 */
double TrotterJ1J2::gJ2() const {

   return J2;

}

/**
 * @return the J2 coupling parameter
 */
int TrotterJ1J2::gPhys_d() const {

   return phys_d;

}
