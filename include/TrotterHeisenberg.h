#ifndef TROTTERHEISENBERG_H
#define TROTTERHEISENBERG_H

#include "HeisenbergMPO.h"
#include "OpSz.h"
#include "OpSx.h"
#include "OpSy.h"

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
      
      //Get the single-site propagator element Op(i,j)
      complex<double> gSingleSiteProp(const int i, const int j) const;
      
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
      
      //The time step
      double dtau;
      
      //Sz
      OpSz *theSz;
      
      //Sx
      OpSx *theSx;
      
      //iSy
      OpSy *theSy;
      
      //If there is a magnetic field: construct single site propagator exp^{theField * dtau * S^z / 2}
      complex<double> *SingleSitePropagator;
      
};

#endif
