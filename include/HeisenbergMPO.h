#ifndef HEISENBERGMPO_H
#define HEISENBERGMPO_H

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

class HeisenbergMPO : public MPO{

   public:
   
      //Constructor
      HeisenbergMPO(const int length, const int phys_d);
      
      //Destructor
      ~HeisenbergMPO();
      
      //Fill MPOdimensions
      void fillMPOdimensions();
      
      //Fill MPOprefactors and MPOoperators
      void fillMPOprefactorsMPOoperators();

      void fillMPOprefactorsMPOoperators_debug();

      //Get the prefactor for a specific MPO term
      complex<double> gPrefactor(const int site, const int row, const int col) const;
      
      //Get the operator for a specific MPO term
      Operator * gOperator(const int site, const int row, const int col) const;
      
      //Get the coupling between spin i and j
      complex<double> gCoupling(const int i, const int j) const;
      
      //Get the magnetic field strength
      complex<double> gField() const;
      
      //Set the coupling between spin i and j
      void sCoupling(const int i, const int j, const complex<double> value);
      
      //Set the magnetic field strength
      void sField(const complex<double> value);
      
      //Print the Heisenberg operators
      void printOperators(ostream& os) const;
      
      //Print the operator tag
      void printOperatorTag(ostream& os, Operator * theOp) const;
      
      //Print its contents
      friend ostream& operator<<(ostream& os, const HeisenbergMPO& theMPO);

      double gJ2();
      
      void sJ2(double);
      
   private:

      double J2;

      complex<double> fOne;
      complex<double> fZero;

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
      
      //Minus the magnetic field --> -h of [- h \sum S^z_i]
      complex<double> MinusH;
      
      //The coupling matrix --> Jij of [\sum_ij J_ij \vec{S}_i . \vec{S}_j]
      complex<double> *couplingMatrix;
      
      //Same as before, but now the pointer to where the variable is stored --> allows for later changes to have influence on the MPO
      complex<double> *gCouplingPointer(const int i, const int j);
      
};

#endif
