#ifndef OPSZ_H
#define OPSZ_H

#include "Operator.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

class OpSz : public Operator{

   public:
   
      //Constructor
      OpSz(const int phys_d);
      
      //Destructor
      ~OpSz();
      
      //Get an element
      double operator()(const int i, const int j) const;
      
      //Am I the zero operator
      bool AmIOp0() const{ return false; }
      
      //Am I the unity matrix
      bool AmIOpI() const{ return false; }
      
};

#endif
