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
      double gJv(const int i, const int j) const;

      double gJeig(const int i) const;
      
      //Get the single-site propagator element Op(i,j)
      complex<double> gSxProp(const int i, const int j) const;
      complex<double> gSyProp(const int i, const int j) const;
      complex<double> gSzProp(const int i, const int j) const;
      
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
      double *Jv;
      
      //The time step
      double dtau;

      //!single site operators
      OpSx *Sx;
      OpSy *Sy;
      OpSz *Sz;

      //!the single site propagators
      complex<double> *SxProp;
      complex<double> *SyProp;
      complex<double> *SzProp;
      
      double *Jeig;
      
};

#endif
