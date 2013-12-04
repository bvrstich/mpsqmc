#include "MPO.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

int MPO::gLength() const{ return length; }

int MPO::gPhys_d() const{ return phys_d; }

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

int MPO::gRN_nTerms() const{ return RN_nTerms; }

double MPO::gRN_prefactor(const int count) const{

   if ((count<0) || (count>=RN_nTerms)){ return 0.0; }
   return RN_prefactors[count];

}

Operator * MPO::gRN_operator(const int count, const int site) const{

   if ((count<0) || (count>=RN_nTerms) || (site<0) || (site>=length)){ return NULL; }
   return RN_operators[count][site];

}

int MPO::gRN_HCfriend(const int index) const{

   if ((index<0) || (index>=RN_nTerms)){ return -1; }
   return RN_HCfriends[index];

}


