#ifndef OPISY_H
#define OPISY_H

#include "Operator.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 10, 2013 */

class OpSy : public Operator{

   public:
   
      //Constructor
      OpSy(const int phys_d);
      
      //Destructor
      ~OpSy();
      
      //Get an element
      complex<double> operator()(const int i, const int j) const;
      
      //Am I the zero operator
      bool AmIOp0() const{ return false; }
      
      //Am I the unity matrix
      bool AmIOpI() const{ return false; }
      
};

#endif
