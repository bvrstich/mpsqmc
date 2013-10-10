#ifndef HEISENBERGMPO_H
#define HEISENBERGMPO_H

#include <stdlib.h>
#include <iostream>

#include "MPO.h"
#include "Op0.h"
#include "OpI.h"
#include "OpSx.h"
#include "OpISy.h"
#include "OpSz.h"
#include "OpSup.h"
#include "OpSdown.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

using namespace std;

class HeisenbergMPO : public MPO{

   public:
   
      //Constructor
      HeisenbergMPO(const int length, const int phys_d, const bool useLadder);
      
      //Destructor
      ~HeisenbergMPO();
      
      //Fill MPOdimensions
      void fillMPOdimensions();
      
      //Fill MPOprefactors and MPOoperators
      void fillMPOprefactorsMPOoperators();
      
      //Get the prefactor for a specific MPO term
      double gPrefactor(const int site, const int row, const int col) const;
      
      //Get the operator for a specific MPO term
      Operator * gOperator(const int site, const int row, const int col) const;
      
      //Get the number of non-zero pairs
      int gNnonzeroPairs() const;
      
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
      friend ostream& operator<<(ostream& os, const HeisenbergMPO& theMPO);
      
      //For the current filled in values of MinusH and couplingMatrix, find the non-zero contributions
      void findNonZeroContributions();
      
   private:
   
      //Whether to use ladder operators (S^+ S^- S^z) or the cartesian components (S^x i*S^y S^z)
      bool useLadder;
   
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
      OpISy * opISy;
      
      //Operator S^z
      OpSz * opSz;
      
      //Operator S^+/sqrt(2)
      OpSup * opSplus;
      
      //Operator S^-/sqrt(2)
      OpSdown * opSminus;
      
      //Minus the magnetic field --> -h of [- h \sum S^z_i]
      double MinusH;
      
      //The coupling matrix --> Jij of [\sum_ij J_ij \vec{S}_i . \vec{S}_j]
      double * couplingMatrix;
      double * minusCouplingMatrix;
      
      //Same as before, but now the pointer to where the variable is stored --> allows for later changes to have influence on the MPO
      double * gCouplingPointer(const int i, const int j);
      double * gMinusCouplingPointer(const int i, const int j);
      
};

#endif
