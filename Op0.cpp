#include "Op0.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

Op0::Op0(const int phys_d) : Operator(){

   this->phys_d = phys_d;

}

Op0::~Op0(){ }

double Op0::operator()(const int i, const int j) const{

   return 0.0;

}

