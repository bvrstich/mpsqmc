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
      
      //Get a certain prefactor
      virtual double gPrefactor(const int site, const int row, const int col) const = 0;
      
      //Get the pointer to a certain Operator
      virtual Operator * gOperator(const int site, const int row, const int col) const = 0;
      
      //Get the number of contributing MPO strings
      int gRN_nTerms() const;

      //Get the prefactor of string number count
      double gRN_prefactor(const int count) const;
      
      //Get the pointer to the operator of string number count, at site
      Operator * gRN_operator(const int count, const int site) const;
      
      //Get the hermitian conjugate index of a specific MPO term
      int gRN_HCfriend(const int index) const;
      
   protected:
   
      //The chain length
      int length;
      
      //The physical dimension (local Hilbert space size)
      int phys_d;
      
      //The MPO dimensions (array of length lenght+1, starting with 1 and ending with 1)
      int * MPOdimensions;
      
      //The MPO operator pointers : MPOoperators[site][row][col] is the pointer to that Operator
      Operator **** MPOoperators;
      
      //The MPO operator prefactor pointers : MPOprefactors[site][row][col] gives the pointer to the prefactor of the Operator at MPOoperators[site][row][col]
      double **** MPOprefactors;
      
      //For the random application of an MPO term: number of contributing strings
      int RN_nTerms;
      
      //For the random application of an MPO term: prefactors of the contributing strings, length RN_nTerms
      double * RN_prefactors;
      
      //For the random application of an MPO term: the strings of operators, length RN_nTerms x length
      Operator *** RN_operators;
      
      //For the random application of an MPO term: array of length RN_nTerms, with the indices (0<=index<RN_nTerms) of the HC MPO terms
      int * RN_HCfriends;
      
};

#endif
