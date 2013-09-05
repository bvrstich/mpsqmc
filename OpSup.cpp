#include <math.h>
#include "OpSup.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

OpSup::OpSup(const int phys_d) : Operator(){

   this->phys_d = phys_d; //S = (phys_d - 1)/2
   storage = new double[phys_d - 1];
   //Order of local basis: smallest M to largest M. Important to know: this operator is actually S^+ / sqrt(2) !!!
   for (int count=0; count<phys_d-1; count++){ storage[count] = 0.5 * sqrt( 0.5 * ( phys_d*phys_d - 1.0 - (1 - phys_d + 2*count)*(3 - phys_d + 2*count) )); }

}

OpSup::~OpSup(){ delete [] storage; }

double OpSup::operator()(const int i, const int j) const{

   if ((i!=j+1) || (i>phys_d) || (j<0)){ return 0.0; }
   return storage[j];

}

