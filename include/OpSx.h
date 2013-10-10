#ifndef OPSX_H
#define OPSX_H

#include "Operator.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 10, 2013 */

class OpSx : public Operator{

   public:
   
      //Constructor
      OpSx(const int phys_d);
      
      //Destructor
      ~OpSx();
      
      //Get an element
      double operator()(const int i, const int j) const;
      
      //Am I the zero operator
      bool AmIOp0() const{ return false; }
      
      //Am I the unity matrix
      bool AmIOpI() const{ return false; }
      
};

#endif
