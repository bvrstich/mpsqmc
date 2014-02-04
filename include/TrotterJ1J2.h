#ifndef TROTTERJ1J2_H
#define TROTTERJ1J2_H

#include "OpSz.h"
#include "OpSx.h"
#include "OpSy.h"
#include "AFMPO.h"
#include "Random.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 14, 2013 */

class TrotterJ1J2{

   public:
   
      //Constructor
      TrotterJ1J2(bool,int,int,double J2, double dtau);
      
      //Destructor
      ~TrotterJ1J2();
      
      double gtau() const;
      
      //Get a particular coupling element J_ij
      complex<double> gJ(const int i, const int j) const;

      complex<double> gV(const int i, const int j) const;

      double gJeig(const int i) const;
      
      complex<double> gAFProp(int myID,int site,const int i, const int j) const;

      void fillAFProp(int myID,int k,int r,complex<double> x);

      void fillAFProp(int myID,int k,complex<double> x,Random *);

      AFMPO *gV_Op(int k,int r);

      int gn_trot() const;

      int glength() const;

      double gJ2() const;

      int gPhys_d() const;
      
   private:
   
      //The chain length
      int length;
      
      //The physical index size
      int phys_d;

      double J2;
      
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

      //!the auxiliary-field propagators
      complex<double> **AFProp;

      //!random direction operator
      complex<double> **Sv;
      
      double *Jeig;
      
};

#endif
