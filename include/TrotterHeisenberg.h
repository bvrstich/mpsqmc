#ifndef TROTTERHEISENBERG_H
#define TROTTERHEISENBERG_H

#include "HeisenbergMPO.h"
#include "OpSz.h"
#include "OpSx.h"
#include "OpISy.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 14, 2013 */

class TrotterHeisenberg{

   public:
   
      //Constructor
      TrotterHeisenberg(HeisenbergMPO *theMPO, double dtau);
      
      //Destructor
      ~TrotterHeisenberg();
      
      //Return whether the magnetic field differs from zero
      bool gIsMagneticField() const;
      
      //Return the magnetic field
      double gField() const;
      
      //Get a particular coupling element J_ij
      double gCoupling(const int i, const int j) const;
      
      //Get for the coupling J, the k-th singular value of the two-site propagator
      double gTwoSitePropSVD_Sing(const double J, const int k) const;
      
      //Get for the coupling J, the element Op1(i,j) corresponding to the k-th singular value of the two-site propagator
      double gTwoSitePropSVD_Left(const double J, const int k, const int i, const int j) const;
      
      //Get for the coupling J, the element Op2(i,j) corresponding to the k-th singular value of the two-site propagator
      double gTwoSitePropSVD_Right(const double J, const int k, const int i, const int j) const;
      
      //Get the single-site propagator element Op(i,j)
      double gSingleSiteProp(const int i, const int j) const;
      
   private:
   
      //Debug print
      static const bool debugPrint = true;
      
      //The chain length
      int length;
      
      //The physical index size
      int phys_d;
      
      //The magnetic field
      double theField;
      
      //Whether there is a magnetic field or not
      bool isMagneticField;
      
      //The coupling matrix
      double *couplingMx;
      
      //The number of different two-site couplings
      int nDifferentCouplings;

      //The different two-site couplings
      double *fDifferentCouplings;
      
      //The time step
      double dtau;
      
      //Sz
      OpSz *theSz;
      
      //Sx
      OpSx *theSx;
      
      //iSy
      OpISy *theISy;
      
      //If there is a magnetic field: construct single site propagator exp^{theField * dtau * S^z / 2}
      double *SingleSitePropagator;
      
      //The propagator per different coupling J = fDifferentCoupling[cnt]: exp^{ -J * dtau * vec{S}_1 * vec{S}_2 }
      //with the convention PropagatorPerCoupling[cnt][row1 + d*col1 + d^2*row2 * d^3*col2]
      double ** PropagatorPerCoupling;
      
      //SVD of PropagatorPerCoupling[row1 + d*col1 + d^2*row2 + d^3*col2] =  sum(j=1..d^2) U[row1 + d*col1 + d*d*j] S[j] VT[j + d*d*(row2 + d*col2)]
      double ** TwoSitePropU;
      double ** TwoSitePropVT;
      double ** TwoSitePropS;
      
      //Set PropagatorPerCoupling, TwoSitePropU, TwoSitePropVT, TwoSitePropS
      void SetTheTwoSitePropagators();
      
};

#endif
