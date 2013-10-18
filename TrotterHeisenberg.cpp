#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "TrotterHeisenberg.h"
#include "OpSz.h"
#include "OpSx.h"
#include "OpISy.h"
#include "Lapack.h"

using namespace std;
const bool TrotterHeisenberg::debugPrint;

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 14, 2013 */

TrotterHeisenberg::TrotterHeisenberg(HeisenbergMPO * theMPO, const double dtau){

   this->length = theMPO->gLength();
   this->phys_d = theMPO->gPhys_d();
   this->dtau = dtau;
   
   theSz = new OpSz(phys_d);
   theSx = new OpSx(phys_d);
   theISy = new OpISy(phys_d);
   
   //The couplings, and the different couplings
   couplingMx = new double[(length * (length - 1))/2];
   nDifferentCouplings = 0;
   fDifferentCouplings = new double[(length * (length - 1))/2];
   for (int first=0; first<length; first++){
      for (int second=first+1; second<length; second++){
         double theCoupling = theMPO->gCoupling(first,second);
         couplingMx[ (second*(second-1))/2 + first ] = theCoupling;
         if (theCoupling != 0.0){
            bool hasNoMatch = true;
            for (int cnt=0; cnt<nDifferentCouplings; cnt++){
               if (theCoupling == fDifferentCouplings[cnt]){ hasNoMatch = false; }
            }
            if (hasNoMatch){
               fDifferentCouplings[nDifferentCouplings] = theCoupling;
               nDifferentCouplings++;
            }
         }
      }
   }
   cout << "TrotterHeisenberg::TrotterHeisenberg --> The different spin couplings J_{ij}: ";
   for (int cnt=0; cnt<nDifferentCouplings; cnt++){ cout << fDifferentCouplings[cnt] << " \t "; }
   cout << endl;
   
   //The magnetic field and single-site propagator
   this->theField = theMPO->gField();
   this->isMagneticField = (theField == 0.0) ? false : true;
   if (isMagneticField){
      SingleSitePropagator = new double[phys_d];
      for (int cnt=0; cnt<phys_d; cnt++){
         SingleSitePropagator[cnt] = exp( 0.5 * theField * dtau * (*theSz)(cnt,cnt) ); //Mind the factor 0.5!!!!
      }
   }
   cout << "TrotterHeisenberg::TrotterHeisenberg --> The magnetic field h = " << theField << " and is it activated? : " << isMagneticField << endl;
   
   //Set the two-site propagators
   SetTheTwoSitePropagators();

}

void TrotterHeisenberg::SetTheTwoSitePropagators(){

   int n = phys_d * phys_d;

   int lwork = 7*n*n + 4*n;
   double * work2 = new double[n*n];
   double * work  = new double[lwork];
   for (int row1=0; row1<phys_d; row1++){
      for (int col1=0; col1<phys_d; col1++){
         for (int row2=0; row2<phys_d; row2++){
            for (int col2=0; col2<phys_d; col2++){
               work2[ row1 + phys_d * ( row2 + phys_d * ( col1 + phys_d * col2 ) ) ] //This is hermitian
                  = (*theSx)(row1,col1) * (*theSx)(row2,col2) - (*theISy)(row1,col1) * (*theISy)(row2,col2) + (*theSz)(row1,col1) * (*theSz)(row2,col2);
            }
         }
      }
   }
   
   double * EigenVecsTwoSite = new double[n*n];
   double * EigenValsTwoSite = new double[n];
   char jobz = 'V';
   char uplo = 'U';
   for (int cnt=0; cnt<n*n; cnt++){ EigenVecsTwoSite[cnt] = work2[cnt]; }
   int info;
   dsyev_(&jobz, &uplo, &n, EigenVecsTwoSite, &n, EigenValsTwoSite, work, &lwork, &info); // TwoSiteOp = EigenVecs2site * EigenVals2site * EigenVecs2site^T
   if (info!=0){ cerr << "TrotterHeisenberg::SetTheTwoSitePropagators --> dsyev_ output info = " << info << endl; }
   
   if (debugPrint){
      cout << "TrotterHeisenberg::SetTheTwoSitePropagators --> Eigenvalues of S1xS2x + S1yS2y + S1zS2z are ";
      for (int bla=0; bla<n; bla++){ cout << EigenValsTwoSite[bla] << " \t ";}
      cout << endl;
      cout << "                                                Remember that this can be checked with [ S(S+1) - S1(S1+1) - S2(S2+1) ]/2" << endl;
      double RMS = 0.0;
      for (int row=0; row<n; row++){
         for (int col=0; col<n; col++){
            work[row + n*col] = 0.0;
            for (int bla=0; bla<n; bla++){
               work[row + n*col] += EigenVecsTwoSite[row + n*bla] * EigenValsTwoSite[bla] * EigenVecsTwoSite[col + n*bla];
            }
            RMS += (work2[row + n*col] - work[row + n*col]) * (work2[row + n*col] - work[row + n*col]);
         }
      }
      cout << "                                                RMS deviation of the eigenvalue decomposition of S1xS2x + S1yS2y + S1zS2z = " << RMS << endl;
   }
   
   PropagatorPerCoupling = new double * [nDifferentCouplings];
   TwoSitePropU          = new double * [nDifferentCouplings];
   TwoSitePropVT         = new double * [nDifferentCouplings];
   TwoSitePropS          = new double * [nDifferentCouplings];
   int * iwork           = new int[8*n];
   
   for (int cnt=0; cnt<nDifferentCouplings; cnt++){
   
      PropagatorPerCoupling[cnt] = new double[n*n];
      TwoSitePropU[cnt]          = new double[n*n];
      TwoSitePropVT[cnt]         = new double[n*n];
      TwoSitePropS[cnt]          = new double[n];
      
      for (int col=0; col<n; col++){
         double preFactor = exp( - 0.5 * dtau * fDifferentCouplings[cnt] * EigenValsTwoSite[col] );
         for (int row=0; row<n; row++){
            work[row + n*col] = preFactor * EigenVecsTwoSite[row + n*col];
         }
      }
      
      char notrans = 'N';
      char trans = 'T';
      double alpha = 1.0;
      double beta = 0.0; //SET
      //PropagatorPerCoupling[cnt] = EigenVecs * exp( -dtau * fDifferentCouplings[cnt] * EigenVals) * EigenVecs^T
      dgemm_(&notrans,&trans,&n,&n,&n,&alpha,work,&n,work,&n,&beta,PropagatorPerCoupling[cnt],&n);
      
      //Propagator[row1 + d*row2, col1 + d*col2] --> Propagator[row1 + d*col1, row2 + d*col2]
      for (int row1=0; row1<phys_d; row1++){
         for (int col2=0; col2<phys_d; col2++){
            for (int row2=0; row2<phys_d; row2++){
               for (int col1=row2+1; col1<phys_d; col1++){
                  double temp = PropagatorPerCoupling[cnt][ row1 + phys_d * ( row2 + phys_d * ( col1 + phys_d * col2 ) ) ];
                  PropagatorPerCoupling[cnt][      row1 + phys_d * ( row2 + phys_d * ( col1 + phys_d * col2 ) ) ]
                     = PropagatorPerCoupling[cnt][ row1 + phys_d * ( col1 + phys_d * ( row2 + phys_d * col2 ) ) ];
                  PropagatorPerCoupling[cnt][ row1 + phys_d * ( col1 + phys_d * ( row2 + phys_d * col2 ) ) ] = temp;
               }
            }
         }
      }
      
      for (int bla=0; bla<n*n; bla++){ work2[bla] = PropagatorPerCoupling[cnt][bla]; }
      jobz = 'S';
      dgesdd_(&jobz, &n, &n, work2, &n, TwoSitePropS[cnt], TwoSitePropU[cnt], &n, TwoSitePropVT[cnt], &n, work, &lwork, iwork, &info);
      if (info!=0){ cerr << "TrotterHeisenberg::SetTheTwoSitePropagators --> dgesdd_ output info = " << info << endl; }
      //Now sum(j=1..n) U[row1 + d*col1 + d*d*j] S[j] VT[j + d*d*(row2 + d*col2)] = PropagatorPerCoupling[row1 + d*col1, row2 + d*col2]
      
      if (debugPrint){
         cout << "TrotterHeisenberg::SetTheTwoSitePropagators --> Singular values for J = " << fDifferentCouplings[cnt] << " are ";
         for (int bla=0; bla<n; bla++){ cout << TwoSitePropS[cnt][bla] << " \t ";}
         cout << endl;
         double RMS = 0.0;
         for (int row=0; row<n; row++){
            for (int col=0; col<n; col++){
               work2[row + n*col] = 0.0;
               for (int bla=0; bla<n; bla++){
                  work2[row + n*col] += TwoSitePropU[cnt][row + n*bla] * TwoSitePropS[cnt][bla] * TwoSitePropVT[cnt][bla + n*col];
               }
               RMS += (work2[row + n*col] - PropagatorPerCoupling[cnt][row + n*col]) * (work2[row + n*col] - PropagatorPerCoupling[cnt][row + n*col]);
            }
         }
         cout << "                                                SVD RMS for this J-value = " << RMS << endl;
         cout << "                                                The operator matrices :" << endl;
         
         for (int bla=0; bla<n; bla++){
            
            cout << "                                                ######   S = " << TwoSitePropS[cnt][bla] << endl;
            cout << "                                                ###      Left = " << endl;
            for (int row=0; row<phys_d; row++){
               cout << "                                                            ";
               for (int col=0; col<phys_d; col++){
                  cout << TwoSitePropU[cnt][row + phys_d * ( col + phys_d * bla)] << "\t";
               }
               cout << endl;
            }
            cout << "                                                ###      Right = " << endl;
            for (int row=0; row<phys_d; row++){
               cout << "                                                            ";
               for (int col=0; col<phys_d; col++){
                  cout << TwoSitePropVT[cnt][bla + phys_d*phys_d * (row + phys_d * col)] << "\t";
               }
               cout << endl;
            }
         }
         
      }
      
   }
   
   delete [] EigenVecsTwoSite;
   delete [] EigenValsTwoSite;
   delete [] iwork;
   delete [] work;
   delete [] work2;

}

TrotterHeisenberg::~TrotterHeisenberg(){

   delete theSz;
   delete theSx;
   delete theISy;
   
   delete [] couplingMx;
   delete [] fDifferentCouplings;

   if (isMagneticField){ delete [] SingleSitePropagator; }
   for (int cnt=0; cnt<nDifferentCouplings; cnt++){
      delete [] PropagatorPerCoupling[cnt];
      delete [] TwoSitePropU[cnt];
      delete [] TwoSitePropVT[cnt];
      delete [] TwoSitePropS[cnt];
   }
   delete [] PropagatorPerCoupling;
   delete [] TwoSitePropU;
   delete [] TwoSitePropVT;
   delete [] TwoSitePropS;

}

bool TrotterHeisenberg::gIsMagneticField() const{ return isMagneticField; }

double TrotterHeisenberg::gField() const{ return theField; }

double TrotterHeisenberg::gCoupling(const int i, const int j) const{

   if ((i==j) || (i<0) || (i>=length) || (j<0) || (j>=length)) {
      cerr << "TrotterHeisenberg::gCoupling  :  The combination of indices (" << i << "," << j << ") is not OK." << endl;
      return NAN;
   }
   
   if (i<j){ return couplingMx[ (j*(j-1))/2 + i ]; }
   return couplingMx[ (i*(i-1))/2 + j ];

}

double TrotterHeisenberg::gTwoSitePropSVD_Sing(const double J, const int k) const{

   if ((k<0) || (k>=phys_d*phys_d)){
      cerr << "TrotterHeisenberg:gSVD_Sing --> variable k out of bound; k = " << k << endl;
      return NAN;
   }
   
   for (int cnt=0; cnt<nDifferentCouplings; cnt++){
      if (J==fDifferentCouplings[cnt]){ return TwoSitePropS[cnt][k]; }
   }
   
   cerr << "TrotterHeisenberg:gSVD_Sing --> J was not found; J = " << J << endl;
   return NAN;

}
      
double TrotterHeisenberg::gTwoSitePropSVD_Left(const double J, const int k, const int i, const int j) const{

   if ((k<0) || (k>=phys_d*phys_d)){
      cerr << "TrotterHeisenberg:gSVD_Left --> variable k out of bound; k = " << k << endl;
      return NAN;
   }
   
   if ((i<0) || (j<0) || (i>=phys_d) || (j>=phys_d)){
      cerr << "TrotterHeisenberg:gSVD_Left --> variable i and/or j out of bound; i = " << i << " and j = " << j << endl;
      return NAN;
   }
   
   for (int cnt=0; cnt<nDifferentCouplings; cnt++){
      if (J==fDifferentCouplings[cnt]){ return TwoSitePropU[cnt][i + phys_d * ( j + phys_d * k ) ]; }
   }
   
   cerr << "TrotterHeisenberg:gSVD_Left --> J was not found; J = " << J << endl;
   return NAN;

}
      
double TrotterHeisenberg::gTwoSitePropSVD_Right(const double J, const int k, const int i, const int j) const{

   if ((k<0) || (k>=phys_d*phys_d)){
      cerr << "TrotterHeisenberg:gSVD_Right --> variable k out of bound; k = " << k << endl;
      return NAN;
   }
   
   if ((i<0) || (j<0) || (i>=phys_d) || (j>=phys_d)){
      cerr << "TrotterHeisenberg:gSVD_Right --> variable i and/or j out of bound; i = " << i << " and j = " << j << endl;
      return NAN;
   }
   
   for (int cnt=0; cnt<nDifferentCouplings; cnt++){
      if (J==fDifferentCouplings[cnt]){ return TwoSitePropVT[cnt][k + phys_d*phys_d * ( i + phys_d * j ) ]; }
   }
   
   cerr << "TrotterHeisenberg:gSVD_Right --> J was not found; J = " << J << endl;
   return NAN;

}

double TrotterHeisenberg::gSingleSiteProp(const int i, const int j) const{

   if ((i<0) || (i>=phys_d) || (j<0) || (j>=phys_d)){
      cerr << "TrotterHeisenberg:gSingleSiteProp --> variable i and/or j out of bound; i = " << i << " and j = " << j << endl;
      return NAN;
   }
   
   if (i!=j){ return 0.0; }
   
   if (isMagneticField){ return SingleSitePropagator[i]; }
   
   return 1.0; //When no magnetic field exp( h * other stuff) = identity

}


