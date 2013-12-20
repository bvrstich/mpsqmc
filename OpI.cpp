#include "OpI.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

OpI::OpI(const int phys_d) : Operator(){

   this->phys_d = phys_d;

}

OpI::~OpI(){ }

complex<double> OpI::operator()(const int i, const int j) const{

   if ((i!=j) || (i>phys_d) || (i<0))
      return complex<double>(0.0,0.0); 

   return complex<double>(1.0,0.0);

}
