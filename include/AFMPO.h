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
      AFMPO(const int length, const int phys_d,int,complex<double> *);
      
      //Destructor
      ~AFMPO();
      
      //Fill MPOdimensions
      void fillMPOdimensions();
      
      //Fill MPOprefactors and MPOoperators
      void fillMPOprefactorsMPOoperators(int,complex<double> *v);
      
      //Get the prefactor for a specific MPO term
      complex<double> gPrefactor(const int site, const int row, const int col) const;
      
      //Get the operator for a specific MPO term
      Operator * gOperator(const int site, const int row, const int col) const;
      
      //Print its contents
      friend ostream& operator<<(ostream& os, const AFMPO& theMPO);
      
   private:
   
      //Value of zero
      complex<double> fZero;
      
      //Value of one
      complex<double> fOne;
      
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
      
};

#endif
