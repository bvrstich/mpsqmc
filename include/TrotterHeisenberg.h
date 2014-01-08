#ifndef TROTTERHEISENBERG_H
#define TROTTERHEISENBERG_H

#include "HeisenbergMPO.h"
#include "OpSz.h"
#include "OpSx.h"
#include "OpSy.h"
#include "AFMPO.h"

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
      complex<double> gField() const;

      double gtau() const;
      
      //Get a particular coupling element J_ij
      complex<double> gJ(const int i, const int j) const;

      complex<double> gV(const int i, const int j) const;

      double gJeig(const int i) const;
      
      //Get the single-site propagator element Op(i,j)
      complex<double> gH1Prop(const int i, const int j) const;

      complex<double> gAFProp(int myID,int site,const int i, const int j) const;

      void fillAFProp(int myID,int k,int r,complex<double> x);

      AFMPO *gV_Op(int k,int r);

      int gn_trot() const;
      
   private:
   
      //The chain length
      int length;
      
      //The physical index size
      int phys_d;
      
      //The magnetic field
      complex<double> theField;
      
      //Whether there is a magnetic field or not
      bool isMagneticField;
      
      //The coupling matrix
      complex<double> *J;
      
      //!the transformation matrix
      complex<double> *V;

      //!the number of non-zero eigenvalues of the coupling matrix: equals the number of trotter product terms
      int n_trot;

      AFMPO **V_Op;
      
      //The time step
      double dtau;

      //!single site operators
      OpSx *Sx;
      OpSy *Sy;
      OpSz *Sz;

      complex<double> *Sx_vec;
      complex<double> *Sy_vec;

      double *eig;

      //!the single site propagators
      complex<double> *H1Prop;

      //!the auxiliary-field propagators
      complex<double> **AFProp;
      
      double *Jeig;
      
};

#endif
