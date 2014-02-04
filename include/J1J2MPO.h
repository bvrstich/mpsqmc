#ifndef J1J2MPO_H
#define J1J2MPO_H

#include <stdlib.h>
#include <iostream>

#include "MPO.h"
#include "Op0.h"
#include "OpI.h"
#include "OpSx.h"
#include "OpSy.h"
#include "OpSz.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

using namespace std;

class J1J2MPO : public MPO{

   public:
   
      //Constructor
      J1J2MPO(const int length, const int phys_d,double J2);
      
      //Destructor
      virtual ~J1J2MPO();
      
      //Fill MPOdimensions
      void fillMPOdimensions();
      
      //Fill MPOprefactors and MPOoperators
      void fillMPOprefactorsMPOoperators();

      //Get the prefactor for a specific MPO term
      complex<double> gPrefactor(const int site, const int row, const int col) const;
      
      //Get the operator for a specific MPO term
      Operator * gOperator(const int site, const int row, const int col) const;
      
      //Get the coupling between spin i and j
      complex<double> gCoupling(const int i, const int j) const;

      complex<double> *gCouplingPointer(const int i, const int j);
      
      //Set the coupling between spin i and j
      void sCoupling(const int i, const int j, const complex<double> value);
      
      //Print its contents
      friend ostream& operator<<(ostream& os, const J1J2MPO& theMPO);
      
   private:

      complex<double> fOne;
      complex<double> fZero;

      complex<double> MinusH;

      //Operator 0
      Op0 *opZero;
      
      //Operator I
      OpI *opOne;
      
      //Operator S^x
      OpSx *opSx;
      
      //Operator S^y
      OpSy *opSy;
      
      //Operator S^z
      OpSz *opSz;
      
      //The coupling matrix --> Jij of [\sum_ij J_ij \vec{S}_i . \vec{S}_j]
      complex<double> *J;
      
};

#endif
