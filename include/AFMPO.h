#ifndef AFMPO_H
#define AFMPO_H

#include <stdlib.h>
#include <iostream>

#include "MPO.h"
#include "Op0.h"
#include "OpI.h"
#include "OpSx.h"
#include "OpSy.h"
#include "OpSz.h"

/**
 * @author Brecht Verstichel
 * class written for the MPO's which represent the auxiliary field operators 'v'. These are 'one-body' operators which are 
 * represented as 2 dimensional MPO's.
 */

using namespace std;

class AFMPO : public MPO{

   public:
   
      //Constructor
      AFMPO(const int length, const int phys_d);
      
      //Destructor
      ~AFMPO();
      
      //Fill MPOdimensions
      void fillMPOdimensions();
      
      //Fill MPOprefactors and MPOoperators
      void fillMPOprefactorsMPOoperators();
      
      //Get the prefactor for a specific MPO term
      double gPrefactor(const int site, const int row, const int col) const;
      
      //Get the operator for a specific MPO term
      Operator * gOperator(const int site, const int row, const int col) const;
      
      //Get the coupling between spin i and j
      double gCoupling(const int i, const int j) const;
      
      //Get the magnetic field strength
      double gField() const;
      
      //Set the coupling between spin i and j
      void sCoupling(const int i, const int j, const double value);
      
      //Set the magnetic field strength
      void sField(const double value);
      
      //Print the Heisenberg operators
      void printOperators(ostream& os) const;
      
      //Print the operator tag
      void printOperatorTag(ostream& os, Operator * theOp) const;
      
      //Print its contents
      friend ostream& operator<<(ostream& os, const AFMPO& theMPO);
      
   private:
   
      //Value of zero
      double fZero;
      
      //Value of one
      double fOne;
      
      //Operator 0
      Op0 * opZero;
      
      //Operator I
      OpI * opOne;
      
      //Operator S^x
      OpSx * opSx;
      
      //Operator i*S^y
      OpSy * opSy;
      
      //Operator S^z
      OpSz * opSz;
      
      //Minus the magnetic field --> -h of [- h \sum S^z_i]
      double MinusH;
      
      //The coupling matrix --> Jij of [\sum_ij J_ij \vec{S}_i . \vec{S}_j]
      double * couplingMatrix;
      
      //Same as before, but now the pointer to where the variable is stored --> allows for later changes to have influence on the MPO
      double * gCouplingPointer(const int i, const int j);
      
};

#endif
