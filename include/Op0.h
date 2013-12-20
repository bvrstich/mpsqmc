#ifndef OP0_H
#define OP0_H

#include "Operator.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

class Op0 : public Operator{

   public:
   
      //Constructor
      Op0(const int phys_d);
      
      //Destructor
      ~Op0();
      
      //Get an element
      complex<double> operator()(const int i, const int j) const;
      
      //Am I the zero operator
      bool AmIOp0() const{ return true; }
      
      //Am I the unity matrix
      bool AmIOpI() const{ return false; }
      
};

#endif
