#ifndef MPO_H
#define MPO_H

#include <stdlib.h>
#include <iostream>

#include "Operator.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

using namespace std;

class MPO{

   public:
      
      //Get the size of the local Hilbert space
      int gPhys_d() const;
      
      //Get the chain length
      int gLength() const;
      
      //Get the MPO dimension left of the site
      int dimL(const int site) const;
      
      //Get the MPO dimension right of the site
      int dimR(const int site) const;

      int gDtrunc() const;

      complex<double> operator()(int,int,int,int,int);
      
      //Get a certain prefactor
      virtual complex<double> gPrefactor(const int site, const int row, const int col) const = 0;
      
      //Get the pointer to a certain Operator
      virtual Operator * gOperator(const int site, const int row, const int col) const = 0;
      
   protected:

      //The chain length
      int length;
      
      //The physical dimension (local Hilbert space size)
      int phys_d;
      
      //The MPO dimensions (array of length lenght+1, starting with 1 and ending with 1)
      int *MPOdimensions;

      //!the maximal dimension of the MPO
      int Dtrunc;
      
      //The MPO operator pointers : MPOoperators[site][row][col] is the pointer to that Operator
      Operator **** MPOoperators;
      
      //The MPO operator prefactor pointers : MPOprefactors[site][row][col] gives the pointer to the prefactor of the Operator at MPOoperators[site][row][col]
      complex<double> **** MPOprefactors;
      
};

#endif
