#include <math.h>
#include "OpSx.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 10, 2013 */

OpSx::OpSx(const int phys_d) : Operator(){

   this->phys_d = phys_d; //S = (phys_d - 1)/2
   storage = new complex<double> [phys_d - 1]; //Upper and lower diagonal are equal --> only store once

   //Order of local basis: smallest M to largest M.
   for(int count = 0;count < phys_d-1;count++)
      storage[count] = complex<double>(0.25 * sqrt( phys_d*phys_d - 1.0 - (1 - phys_d + 2*count)*(3 - phys_d + 2*count) ),0.0);

}

OpSx::~OpSx(){ 
   
   delete [] storage; 
   
}

complex<double> OpSx::operator()(const int i, const int j) const{

   if( ((i!=j+1) && (i+1!=j)) || (j>phys_d) || (i>phys_d) || (i<0) || (j<0) )
      return complex<double>(0.0,0.0);

   if(i==j+1)
      return storage[j];

   return storage[i];

}
