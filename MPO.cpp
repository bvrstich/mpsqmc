#include "MPO.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

int MPO::gLength() const{
   
   return length; 
   
}

int MPO::gPhys_d() const {
   
   return phys_d; 
   
}

int MPO::dimL(const int site) const{

   if ((site<0) || (site>=length))
      return 0; 

   return MPOdimensions[site];
 
}

int MPO::dimR(const int site) const{

   if ((site<0) || (site>=length)) 
      return 0;

   return MPOdimensions[site+1];

}
