#include "OpSz.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

OpSz::OpSz(const int phys_d) : Operator(){

   this->phys_d = phys_d; //S = (phys_d - 1)/2
   storage = new complex<double>[phys_d];

   for (int count=0; count<phys_d; count++)
      storage[count] = complex<double>(0.5 * (1 - phys_d) + count,0.0); //Order of local basis: smallest M to largest M

}

OpSz::~OpSz(){ 
   
   delete [] storage; 
   
}

complex<double> OpSz::operator()(const int i, const int j) const{

   if ((i!=j) || (i>phys_d) || (i<0))
      return 0.0;

   return storage[i];

}
