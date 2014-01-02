#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "AFMPO.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

AFMPO::AFMPO(const int length, const int phys_d,int r,complex<double> *v) : MPO(){

   this->length = length;
   this->phys_d = phys_d;

   this->opZero   = new Op0(phys_d);
   this->opOne    = new OpI(phys_d);

   this->opSx  = new OpSx(phys_d);
   this->opSy = new OpSy(phys_d);
   this->opSz     = new OpSz(phys_d);
   
   fZero = complex<double>(0.0,0.0);
   fOne = complex<double>(1.0,0.0);
   
   MPOdimensions = new int [length+1];
   fillMPOdimensions();
   
   MPOprefactors = new complex<double> *** [length];
   MPOoperators  = new Operator *** [length];

   for(int site = 0;site < length;site++){

      MPOprefactors[site] = new complex<double> ** [dimL(site)];
      MPOoperators[site]  = new Operator ** [dimL(site)];

      for(int row=0; row<dimL(site); row++){

         MPOprefactors[site][row] = new complex<double> * [dimR(site)];
         MPOoperators[site][row]  = new Operator * [dimR(site)];

      }

   }

   fillMPOprefactorsMPOoperators(r,v);

}

AFMPO::~AFMPO(){

   delete opZero;
   delete opOne;
   delete opSz;
   delete opSx;
   delete opSy;
   
   for (int site=0; site<length; site++){

      for (int row=0; row<dimL(site); row++){

         delete [] MPOprefactors[site][row];
         delete [] MPOoperators[site][row];

      }

      delete [] MPOprefactors[site];
      delete [] MPOoperators[site];

   }

   delete [] MPOprefactors;
   delete [] MPOoperators;
   
   delete [] MPOdimensions;
   
}

void AFMPO::fillMPOdimensions(){

   MPOdimensions[0] = MPOdimensions[length] = 1;

   for (int site = 0;site < length;site++)
      MPOdimensions[site+1] = 2;

}

complex<double> AFMPO::gPrefactor(const int site, const int row, const int col) const{

   if ((site<0) || (site>=length) || (row<0) || (row>=dimL(site)) || (col<0) || (col>=dimR(site)))
      return 0.0;

   return *MPOprefactors[site][row][col];

}

Operator * AFMPO::gOperator(const int site, const int row, const int col) const{

   if ((site<0) || (site>=length) || (row<0) || (row>=dimL(site)) || (col<0) || (col>=dimR(site)))
      return NULL; 

   if(abs(gPrefactor(site,row,col)) < 1.0e-15)
      return opZero;

   return MPOoperators[site][row][col];

}

void AFMPO::fillMPOprefactorsMPOoperators(int r,complex<double> *v){

   //first site
   MPOprefactors[0][0][0] = &fOne;
   MPOoperators[0][0][0] = opOne;

   MPOprefactors[0][0][1] = v;

   if(r == 0)
      MPOoperators[0][0][1] = opSx;
   else if(r == 1)
      MPOoperators[0][0][1] = opSy;
   else if(r == 2)
      MPOoperators[0][0][1] = opSz;

   for(int site = 1;site < length - 1;++site){
      
      MPOprefactors[site][0][0] = &fOne;
      MPOoperators[site][0][0] = opOne;

      MPOprefactors[site][0][1] = v + site;

      if(r == 0)
         MPOoperators[site][0][1] = opSx;
      else if(r == 1)
         MPOoperators[site][0][1] = opSy;
      else if(r == 2)
         MPOoperators[site][0][1] = opSz;
 
      MPOprefactors[site][1][0] = &fZero;
      MPOoperators[site][1][0] = opZero;

      MPOprefactors[site][1][1] = &fOne;
      MPOoperators[site][1][1] = opOne;

   }

   //last site
   MPOprefactors[length - 1][0][0] = v + length - 1;

   if(r == 0)
      MPOoperators[length-1][0][0] = opSx;
   else if(r == 1)
      MPOoperators[length-1][0][0] = opSy;
   else if(r == 2)
      MPOoperators[length-1][0][0] = opSz;

   MPOprefactors[length-1][1][0] = &fOne;
   MPOoperators[length-1][1][0] = opOne;

}
