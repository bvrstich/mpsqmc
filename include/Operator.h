#ifndef OPERATOR_H
#define OPERATOR_H

#include <stdlib.h>
#include <iostream>

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

using namespace std;

class Operator{

   public:
      
      //Get the size of the local Hilbert space
      int gPhys_d() const;
      
      //Get an element of of the operator; j is the column [C_i = sum_j OP(i,j) C_j]
      virtual double operator()(const int i, const int j) const = 0;
      
      //Print its contents
      friend ostream& operator<<(ostream& os, const Operator& theOp);
      
      //Am I the zero operator
      virtual bool AmIOp0() const = 0;
      
      //Am I the unity matrix
      virtual bool AmIOpI() const = 0;
      
   protected:
   
      //The physical dimension (local Hilbert space size)
      int phys_d;
      
      //Storage for what you want to store (can be used sparse though)
      double * storage;
      
};

#endif
