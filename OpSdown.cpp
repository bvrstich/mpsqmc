#include <math.h>
#include "OpSdown.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

OpSdown::OpSdown(const int phys_d) : Operator(){

   this->phys_d = phys_d; //S = (phys_d - 1)/2
   storage = new double[phys_d - 1];
   //Order of local basis: smallest M to largest M. Important to know: this operator is actually S^- / sqrt(2) !!!
   for (int count=0; count<phys_d-1; count++){ storage[count] = 0.5 * sqrt( 0.5 * ( phys_d*phys_d - 1.0 - (1 - phys_d + 2*count)*(3 - phys_d + 2*count) )); }

}

OpSdown::~OpSdown(){ delete [] storage; }

double OpSdown::operator()(const int i, const int j) const{

   if ((i+1!=j) || (j>phys_d) || (i<0)){ return 0.0; }
   return storage[i];

}

