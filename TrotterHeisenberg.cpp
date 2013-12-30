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

   theSz = new OpSz(phys_d);
   theSx = new OpSx(phys_d);
   theSy = new OpSy(phys_d);

   //The couplings, and the different couplings
   couplingMx = new double [(length * (length - 1))/2];

   for(int first = 0;first < length;first++)
      for(int second = first + 1;second < length;second++){

         double theCoupling = theMPO->gCoupling(first,second);

         //couplingmatrix will contain the coupling J_[ij]
         couplingMx[ (second*(second-1))/2 + first ] = theCoupling;

      }

   //The magnetic field and single-site propagator
   this->theField = theMPO->gField();

   this->isMagneticField = (fabs(theField) < 1.0e-15) ? false : true;

   if(isMagneticField){

      SingleSitePropagator = new complex<double> [phys_d];

      //e^Sz is diagonal in the site-basis, so just a d-dimensional vector:
      for (int cnt = 0;cnt < phys_d;cnt++)
         SingleSitePropagator[cnt] = exp( 0.5 * theField * dtau * (*theSz)(cnt,cnt) ); //Mind the factor 0.5!!!!

   }

}

TrotterHeisenberg::~TrotterHeisenberg(){

   delete theSz;
   delete theSx;
   delete theSy;

   delete [] couplingMx;

   if (isMagneticField)
      delete [] SingleSitePropagator;

}

bool TrotterHeisenberg::gIsMagneticField() const {

   return isMagneticField; 

}

double TrotterHeisenberg::gField() const {

   return theField; 
   
}

double TrotterHeisenberg::gCoupling(const int i, const int j) const {

   if((i==j) || (i<0) || (i>=length) || (j<0) || (j>=length)) {

      cerr << "TrotterHeisenberg::gCoupling  :  The combination of indices (" << i << "," << j << ") is not OK." << endl;

      return NAN;

   }

   if(i < j)
      return couplingMx[ (j*(j-1))/2 + i ]; 

   return couplingMx[ (i*(i-1))/2 + j ];

}

/**
 * @return element (i,j) of the single site propagotor e^{tau h S_z)
 */
complex<double> TrotterHeisenberg::gSingleSiteProp(const int i, const int j) const{

   if ((i<0) || (i>=phys_d) || (j<0) || (j>=phys_d)){

      cerr << "TrotterHeisenberg:gSingleSiteProp --> variable i and/or j out of bound; i = " << i << " and j = " << j << endl;

      return NAN;

   }

   //always diagonal
   if(i != j)
      return 0.0;

   if(isMagneticField)
      return SingleSitePropagator[i];

   return 1.0; //When no magnetic field exp( h * other stuff) = identity

}
