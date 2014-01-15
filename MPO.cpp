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

int MPO::gDtrunc() const{

   return Dtrunc;

}

complex<double> MPO::operator()(int site,int vL,int su,int sd,int vR){

   return MPOprefactors[site][vL][vR][0] * (*MPOoperators[site][vL][vR])(su,sd);

}
