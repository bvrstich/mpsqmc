#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "HeisenbergMPO.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

HeisenbergMPO::HeisenbergMPO(const int length, const int phys_d) : MPO(){

   if (length < 3) 
      cerr << "HeisenbergMPO::HeisenbergMPO  :  The length must be at least three for the program to work." << endl;

   this->length = length;
   this->phys_d = phys_d;
   this->opZero   = new Op0(phys_d);
   this->opOne    = new OpI(phys_d);

   this->opSx  = new OpSx(phys_d);
   this->opSy = new OpSy(phys_d);
   this->opSz     = new OpSz(phys_d);
   
   MinusH = -0.0;
   couplingMatrix      = new double[(length * (length - 1))/2];

   for (int count=0; count<(length * (length - 1))/2; count++)
      couplingMatrix[count] = 0.0; 
   
   fZero = 0.0;
   fOne = 1.0;
   
   MPOdimensions = new int [length+1];
   fillMPOdimensions();
   
   MPOprefactors = new double *** [length];
   MPOoperators  = new Operator *** [length];

   for(int site = 0;site < length;site++){

      MPOprefactors[site] = new double**[dimL(site)];
      MPOoperators[site]  = new Operator**[dimL(site)];

      for(int row=0; row<dimL(site); row++){

         MPOprefactors[site][row] = new double*[dimR(site)];
         MPOoperators[site][row]  = new Operator*[dimR(site)];

      }

   }

   fillMPOprefactorsMPOoperators();

}

HeisenbergMPO::~HeisenbergMPO(){

   delete opZero;
   delete opOne;
   delete opSz;
   delete opSx;
   delete opSy;
   
   delete [] couplingMatrix;
   
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

void HeisenbergMPO::fillMPOdimensions(){

   int locSwitch = length/2;

   MPOdimensions[0] = MPOdimensions[length] = 1;

   for (int site=0; site< locSwitch; site++)
      MPOdimensions[site+1] = 2 + 3*(site+1);

   for (int site=length-1; site>locSwitch; site--)
      MPOdimensions[site]   = 2 + 3*(length-site);

}

void HeisenbergMPO::fillMPOprefactorsMPOoperators(){

   int locSwitch = length/2;
   
   //First and last site are pretty easy (length of min. 3 is assumed)
   MPOprefactors[0][0][0] = &MinusH;
   MPOoperators[0][0][0] = opSz;
   MPOprefactors[0][0][1] = &fOne;
   MPOoperators[0][0][1] = opSx;
   MPOprefactors[0][0][2] = &fOne;
   MPOoperators[0][0][2] = opSy;
   MPOprefactors[0][0][3] = &fOne;
   MPOoperators[0][0][3] = opSz;
   MPOprefactors[0][0][4] = &fOne;
   MPOoperators[0][0][4] = opOne;
   
   MPOprefactors[length-1][0][0] = &fOne;
   MPOoperators[length-1][0][0] = opOne;
   MPOprefactors[length-1][1][0] = &fOne;
   MPOoperators[length-1][1][0] = opSx;
   MPOprefactors[length-1][2][0] = &fOne;
   MPOoperators[length-1][2][0] = opSy;
   MPOprefactors[length-1][3][0] = &fOne;
   MPOoperators[length-1][3][0] = opSz;
   MPOprefactors[length-1][4][0] = &MinusH;
   MPOoperators[length-1][4][0] = opSz;
   
   //What is to be done before the switch site, is also well known
   for (int site=1; site<locSwitch; site++){
   
      for (int row=0; row<dimL(site); row++){
         for (int col=0; col<dimR(site); col++){
            MPOprefactors[site][row][col] = &fZero;
            MPOoperators[site][row][col] = opZero;
         }
      }
      for (int row=0; row<dimL(site)-1; row++){
         MPOprefactors[site][row][row] = &fOne;
         MPOoperators[site][row][row] = opOne;
      }
      MPOprefactors[site][dimL(site)-1][0] = &MinusH;
      MPOoperators[site][dimL(site)-1][0] = opSz;
      MPOprefactors[site][dimL(site)-1][dimL(site)-1] = &fOne;
      MPOoperators[site][dimL(site)-1][dimL(site)-1] = opSx;
      MPOprefactors[site][dimL(site)-1][dimL(site)] = &fOne;
      MPOoperators[site][dimL(site)-1][dimL(site)] = opSy;
      MPOprefactors[site][dimL(site)-1][dimL(site)+1] = &fOne;
      MPOoperators[site][dimL(site)-1][dimL(site)+1] = opSz;
      MPOprefactors[site][dimL(site)-1][dimL(site)+2] = &fOne;
      MPOoperators[site][dimL(site)-1][dimL(site)+2] = opOne;

      for (int row=0; row<site; row++){

         double * value = gCouplingPointer(row,site);

         MPOprefactors[site][1 + 3*row + 0][0] = value;
         MPOoperators[site][ 1 + 3*row + 0][0] = opSx;
         MPOprefactors[site][1 + 3*row + 1][0] = value;
         MPOoperators[site][ 1 + 3*row + 1][0] = opSy;
         MPOprefactors[site][1 + 3*row + 2][0] = value;
         MPOoperators[site][ 1 + 3*row + 2][0] = opSz;

      }
   
   }
   
   //On the switch site, we change stored operators to complementary operators
   {
   
      for (int row=0; row<dimL(locSwitch); row++)
         for (int col=0; col<dimR(locSwitch); col++){

            MPOprefactors[locSwitch][row][col] = &fZero;
            MPOoperators[ locSwitch][row][col] = opZero;

         }

      MPOprefactors[locSwitch][0][0] = &fOne;
      MPOoperators[locSwitch][0][0] = opOne;
      MPOprefactors[locSwitch][dimL(locSwitch)-1][0] = &MinusH;
      MPOoperators[ locSwitch][dimL(locSwitch)-1][0] = opSz;

      for (int row=0; row<locSwitch; row++){

         double * value = gCouplingPointer(row,locSwitch);

         MPOprefactors[locSwitch][1 + 3*row + 0][0] = value;
         MPOoperators[ locSwitch][1 + 3*row + 0][0] = opSx;
         MPOprefactors[locSwitch][1 + 3*row + 1][0] = value;
         MPOoperators[ locSwitch][1 + 3*row + 1][0] = opSy;
         MPOprefactors[locSwitch][1 + 3*row + 2][0] = value;
         MPOoperators[ locSwitch][1 + 3*row + 2][0] = opSz;

         for (int col=0; col<length-1-locSwitch; col++){

            value = gCouplingPointer(row,locSwitch+1+col);

            for (int cnt=0; cnt<3; cnt++){

               MPOprefactors[locSwitch][1 + 3*row + cnt][1 + 3*col + cnt] = value;
               MPOoperators[ locSwitch][1 + 3*row + cnt][1 + 3*col + cnt] = opOne;

            }

         }

      }

      for (int col=0; col<length-1-locSwitch; col++){

         double * value = gCouplingPointer(locSwitch,locSwitch+1+col);

         MPOprefactors[locSwitch][dimL(locSwitch)-1][1 + 3*col + 0] = value;
         MPOoperators[ locSwitch][dimL(locSwitch)-1][1 + 3*col + 0] = opSx;
         MPOprefactors[locSwitch][dimL(locSwitch)-1][1 + 3*col + 1] = value;
         MPOoperators[ locSwitch][dimL(locSwitch)-1][1 + 3*col + 1] = opSy;
         MPOprefactors[locSwitch][dimL(locSwitch)-1][1 + 3*col + 2] = value;
         MPOoperators[ locSwitch][dimL(locSwitch)-1][1 + 3*col + 2] = opSz;

      }

      MPOprefactors[locSwitch][dimL(locSwitch)-1][1 + 3*(length-1-locSwitch) + 0] = &fOne;
      MPOoperators[ locSwitch][dimL(locSwitch)-1][1 + 3*(length-1-locSwitch) + 0] = opOne;
   
   }
   
   //After the switch site, we need to store complementary operators
   for (int site=locSwitch+1; site<length-1; site++){
   
      for (int row=0; row<dimL(site); row++)
         for (int col=0; col<dimR(site); col++){

            MPOprefactors[site][row][col] = &fZero;
            MPOoperators[site][row][col] = opZero;

         }
      
      MPOprefactors[site][0][0] = &fOne;
      MPOoperators[site][0][0] = opOne;
      MPOprefactors[site][dimL(site)-1][0] = &MinusH;
      MPOoperators[site][dimL(site)-1][0] = opSz;
      MPOprefactors[site][1][0] = &fOne;
      MPOoperators[ site][1][0] = opSx;
      MPOprefactors[site][2][0] = &fOne;
      MPOoperators[ site][2][0] = opSy;
      MPOprefactors[site][3][0] = &fOne;
      MPOoperators[ site][3][0] = opSz;

      for (int col=1; col<dimR(site); col++){

         MPOprefactors[site][3 + col][col] = &fOne;
         MPOoperators[ site][3 + col][col] = opOne;

      }

      for (int col=0; col<length-1-site; col++){

         double * value = gCouplingPointer(site,site+1+col);

         MPOprefactors[site][dimL(site)-1][1 + 3*col + 0] = value;
         MPOoperators[ site][dimL(site)-1][1 + 3*col + 0] = opSx;
         MPOprefactors[site][dimL(site)-1][1 + 3*col + 1] = value;
         MPOoperators[ site][dimL(site)-1][1 + 3*col + 1] = opSy;
         MPOprefactors[site][dimL(site)-1][1 + 3*col + 2] = value;
         MPOoperators[ site][dimL(site)-1][1 + 3*col + 2] = opSz;

      }
   
   }

}

void HeisenbergMPO::sField(const double value){
   
   MinusH = -value; 
   
}

double HeisenbergMPO::gField() const{
   
   return -MinusH; 
   
}

void HeisenbergMPO::sCoupling(const int i, const int j, const double value){

   if ((i==j) || (i<0) || (i>=length) || (j<0) || (j>=length)) {
      cerr << "HeisenbergMPO::sCoupling  :  The combination of indices (" << i << "," << j << ") is not OK." << endl;
      return;
   }
   
   if (i<j)
      couplingMatrix[ (j*(j-1))/2 + i ]      = value;
   else 
      couplingMatrix[ (i*(i-1))/2 + j ]      = value;

}

double HeisenbergMPO::gCoupling(const int i, const int j) const{

   if ((i==j) || (i<0) || (i>=length) || (j<0) || (j>=length)) {
      cerr << "HeisenbergMPO::gCoupling  :  The combination of indices (" << i << "," << j << ") is not OK." << endl;
      return NAN;
   }
   
   if (i<j){ return couplingMatrix[ (j*(j-1))/2 + i ]; }
   return couplingMatrix[ (i*(i-1))/2 + j ];

}

double * HeisenbergMPO::gCouplingPointer(const int i, const int j){

   if ((i==j) || (i<0) || (i>=length) || (j<0) || (j>=length)) {
      cerr << "HeisenbergMPO::gCouplingPointer  :  The combination of indices (" << i << "," << j << ") is not OK." << endl;
      return NULL;
   }
   
   if (i<j){ return couplingMatrix + (j*(j-1))/2 + i; }
   return couplingMatrix + (i*(i-1))/2 + j;

}

double HeisenbergMPO::gPrefactor(const int site, const int row, const int col) const{

   if ((site<0) || (site>=length) || (row<0) || (row>=dimL(site)) || (col<0) || (col>=dimR(site))){ return 0.0; }
   return MPOprefactors[site][row][col][0];

}

Operator * HeisenbergMPO::gOperator(const int site, const int row, const int col) const{

   if ((site<0) || (site>=length) || (row<0) || (row>=dimL(site)) || (col<0) || (col>=dimR(site))){ return NULL; }
   if (gPrefactor(site,row,col)==0.0){ return opZero; }
   return MPOoperators[site][row][col];

}

void HeisenbergMPO::printOperators(ostream& os) const{

   os << *opZero;
   os << *opOne;
   os << *opSz;
   os << *opSx;
   os << *opSy;

}

void HeisenbergMPO::printOperatorTag(ostream& os, Operator * theOp) const{

   if (theOp == opZero){   os << "---"; }
   if (theOp == opOne){    os << "I  "; }
   if (theOp == opSz){     os << "Sz ";  }
   if (theOp == opSx){  os << "Sx "; }
   if (theOp == opSy){ os << "Sy"; }
   

}

ostream& operator<<(ostream& os, const HeisenbergMPO& theMPO){

   theMPO.printOperators(os);
   
   os << "---------------------------------------------------------------------------------------" << endl;

   os << "HeisenbergMPO  with  L = " << theMPO.gLength() << "  and  d = " << theMPO.gPhys_d() << endl;
   os << "Field = " << theMPO.gField() << endl;
   os << "Non-zero coupling between sites (i,j) = value :" << endl;
   for (int i=0; i<theMPO.gLength()-1; i++){
      for (int j=i+1; j<theMPO.gLength(); j++){
         double value = theMPO.gCoupling(i,j);
         if (value != 0.0){ os << "\t\t(" << i << "," << j << ") = " << value << endl; }
      }
   }
   os << "MPO dimensions = \n\t\t1\t";
   for (int site=0; site<theMPO.gLength(); site++){ os << theMPO.dimR(site) << "\t"; }
   os << endl;
   
   os << "---------------------------------------------------------------------------------------" << endl;
   
   for (int site=0; site<theMPO.gLength(); site++){
      os << "MPO at site " << site << " is :" << endl;
      for (int row=0; row<theMPO.dimL(site); row++){
         os << "\t";
         for (int col=0; col<theMPO.dimR(site); col++){
            double value = theMPO.gPrefactor(site,row,col);
            if (value>=0.0){ os << "( " << value << ";"; }
            else { os << "(" << value << ";"; }
            theMPO.printOperatorTag(os, theMPO.gOperator(site,row,col));
            os << ")\t";
         }
         os << endl;
      }
      os << "---------------------------------------------------------------------------------------" << endl;
   }
   
   return os;
   
}
