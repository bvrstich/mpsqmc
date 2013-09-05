#ifndef OPI_H
#define OPI_H

#include "Operator.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

class OpI : public Operator{

   public:
   
      //Constructor
      OpI(const int phys_d);
      
      //Destructor
      ~OpI();
      
      //Get an element
      double operator()(const int i, const int j) const;
      
      //Am I the zero operator
      bool AmIOp0() const{ return false; }
      
      //Am I the unity matrix
      bool AmIOpI() const{ return true; }
      
};

#endif
