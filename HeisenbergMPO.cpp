#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "HeisenbergMPO.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

HeisenbergMPO::HeisenbergMPO(const int length, const int phys_d, const bool useLadder) : MPO(){

   if (length < 3){ cerr << "HeisenbergMPO::HeisenbergMPO  :  The length must be at least three for the program to work." << endl; }

   this->length = length;
   this->phys_d = phys_d;
   this->useLadder = useLadder;
   this->opZero   = new Op0(phys_d);
   this->opOne    = new OpI(phys_d);
   this->opSz     = new OpSz(phys_d);
   if (useLadder){
      this->opSplus  = new OpSup(phys_d);
      this->opSminus = new OpSdown(phys_d);
   } else {
      this->opSx  = new OpSx(phys_d);
      this->opISy = new OpISy(phys_d);
   }
   
   MinusH = -0.0;
   couplingMatrix      = new double[(length * (length - 1))/2];
   minusCouplingMatrix = new double[(length * (length - 1))/2];
   for (int count=0; count<(length * (length - 1))/2; count++){ minusCouplingMatrix[count] = couplingMatrix[count] = 0.0; }
   
   fZero = 0.0;
   fOne = 1.0;
   
   MPOdimensions = new int[length+1];
   fillMPOdimensions();
   
   MPOprefactors = new double***[length];
   MPOoperators  = new Operator***[length];
   for (int site=0; site<length; site++){
      MPOprefactors[site] = new double**[dimL(site)];
      MPOoperators[site]  = new Operator**[dimL(site)];
      for (int row=0; row<dimL(site); row++){
         MPOprefactors[site][row] = new double*[dimR(site)];
         MPOoperators[site][row]  = new Operator*[dimR(site)];
      }
   }
   fillMPOprefactorsMPOoperators();
   
   RN_nTerms = 0; //Means not allocated

}

HeisenbergMPO::~HeisenbergMPO(){

   delete opZero;
   delete opOne;
   delete opSz;
   if (useLadder){
      delete opSplus;
      delete opSminus;
   } else {
      delete opSx;
      delete opISy;
   }
   delete [] couplingMatrix;
   delete [] minusCouplingMatrix;
   
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
   
   if (RN_nTerms > 0){
      delete [] RN_prefactors;
      for (int count=0; count<RN_nTerms; count++){ delete [] RN_operators[count]; }
      delete [] RN_operators;
      delete [] RN_HCfriends;
   }

}

void HeisenbergMPO::fillMPOdimensions(){

   int locSwitch = length/2;
   MPOdimensions[0] = MPOdimensions[length] = 1;
   for (int site=0; site<locSwitch; site++){        MPOdimensions[site+1] = 2 + 3*(site+1);      }
   for (int site=length-1; site>locSwitch; site--){ MPOdimensions[site]   = 2 + 3*(length-site); }

}

void HeisenbergMPO::fillMPOprefactorsMPOoperators(){

   int locSwitch = length/2;
   
   //First and last site are pretty easy (length of min. 3 is assumed)
   MPOprefactors[0][0][0] = &MinusH;
   MPOoperators[0][0][0] = opSz;
   MPOprefactors[0][0][1] = &fOne;
   MPOoperators[0][0][1] = (useLadder) ? (Operator *) opSplus : (Operator *) opSx;
   MPOprefactors[0][0][2] = &fOne;
   MPOoperators[0][0][2] = (useLadder) ? (Operator *) opSminus : (Operator *) opISy;
   MPOprefactors[0][0][3] = &fOne;
   MPOoperators[0][0][3] = opSz;
   MPOprefactors[0][0][4] = &fOne;
   MPOoperators[0][0][4] = opOne;
   
   MPOprefactors[length-1][0][0] = &fOne;
   MPOoperators[length-1][0][0] = opOne;
   MPOprefactors[length-1][1][0] = &fOne;
   MPOoperators[length-1][1][0] = (useLadder) ? (Operator *) opSminus : (Operator *) opSx;
   MPOprefactors[length-1][2][0] = &fOne;
   MPOoperators[length-1][2][0] = (useLadder) ? (Operator *) opSplus : (Operator *) opISy;
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
      MPOoperators[site][dimL(site)-1][dimL(site)-1] = (useLadder) ? (Operator *) opSplus : (Operator *) opSx;
      MPOprefactors[site][dimL(site)-1][dimL(site)] = &fOne;
      MPOoperators[site][dimL(site)-1][dimL(site)] = (useLadder) ? (Operator *) opSminus : (Operator *) opISy;
      MPOprefactors[site][dimL(site)-1][dimL(site)+1] = &fOne;
      MPOoperators[site][dimL(site)-1][dimL(site)+1] = opSz;
      MPOprefactors[site][dimL(site)-1][dimL(site)+2] = &fOne;
      MPOoperators[site][dimL(site)-1][dimL(site)+2] = opOne;
      for (int row=0; row<site; row++){
         double * value = gCouplingPointer(row,site);
         MPOprefactors[site][1 + 3*row + 0][0] = value;
         MPOoperators[site][ 1 + 3*row + 0][0] = (useLadder) ? (Operator *) opSminus : (Operator *) opSx;
         MPOprefactors[site][1 + 3*row + 1][0] = (useLadder) ? value : gMinusCouplingPointer(row,site);
         MPOoperators[site][ 1 + 3*row + 1][0] = (useLadder) ? (Operator *) opSplus : (Operator *) opISy;
         MPOprefactors[site][1 + 3*row + 2][0] = value;
         MPOoperators[site][ 1 + 3*row + 2][0] = opSz;
      }
   
   }
   
   //On the switch site, we change stored operators to complementary operators
   {
   
      for (int row=0; row<dimL(locSwitch); row++){
         for (int col=0; col<dimR(locSwitch); col++){
            MPOprefactors[locSwitch][row][col] = &fZero;
            MPOoperators[ locSwitch][row][col] = opZero;
         }
      }
      MPOprefactors[locSwitch][0][0] = &fOne;
      MPOoperators[locSwitch][0][0] = opOne;
      MPOprefactors[locSwitch][dimL(locSwitch)-1][0] = &MinusH;
      MPOoperators[ locSwitch][dimL(locSwitch)-1][0] = opSz;
      for (int row=0; row<locSwitch; row++){
         double * value = gCouplingPointer(row,locSwitch);
         MPOprefactors[locSwitch][1 + 3*row + 0][0] = value;
         MPOoperators[ locSwitch][1 + 3*row + 0][0] = (useLadder) ? (Operator *) opSminus : (Operator *) opSx;
         MPOprefactors[locSwitch][1 + 3*row + 1][0] = (useLadder) ? value : gMinusCouplingPointer(row,locSwitch);
         MPOoperators[ locSwitch][1 + 3*row + 1][0] = (useLadder) ? (Operator *) opSplus : (Operator *) opISy;
         MPOprefactors[locSwitch][1 + 3*row + 2][0] = value;
         MPOoperators[ locSwitch][1 + 3*row + 2][0] = opSz;
         for (int col=0; col<length-1-locSwitch; col++){
            value = gCouplingPointer(row,locSwitch+1+col);
            for (int cnt=0; cnt<3; cnt++){
               MPOprefactors[locSwitch][1 + 3*row + cnt][1 + 3*col + cnt] = (useLadder) ? value : ((cnt==1) ? gMinusCouplingPointer(row,locSwitch+1+col) : value);
               MPOoperators[ locSwitch][1 + 3*row + cnt][1 + 3*col + cnt] = opOne;
            }
         }
      }
      for (int col=0; col<length-1-locSwitch; col++){
         double * value = gCouplingPointer(locSwitch,locSwitch+1+col);
         MPOprefactors[locSwitch][dimL(locSwitch)-1][1 + 3*col + 0] = value;
         MPOoperators[ locSwitch][dimL(locSwitch)-1][1 + 3*col + 0] = (useLadder) ? (Operator *) opSplus : (Operator *) opSx;
         MPOprefactors[locSwitch][dimL(locSwitch)-1][1 + 3*col + 1] = (useLadder) ? value : gMinusCouplingPointer(locSwitch,locSwitch+1+col);
         MPOoperators[ locSwitch][dimL(locSwitch)-1][1 + 3*col + 1] = (useLadder) ? (Operator *) opSminus : (Operator *) opISy;
         MPOprefactors[locSwitch][dimL(locSwitch)-1][1 + 3*col + 2] = value;
         MPOoperators[ locSwitch][dimL(locSwitch)-1][1 + 3*col + 2] = opSz;
      }
      MPOprefactors[locSwitch][dimL(locSwitch)-1][1 + 3*(length-1-locSwitch) + 0] = &fOne;
      MPOoperators[ locSwitch][dimL(locSwitch)-1][1 + 3*(length-1-locSwitch) + 0] = opOne;
   
   }
   
   //After the switch site, we need to store complementary operators
   for (int site=locSwitch+1; site<length-1; site++){
   
      for (int row=0; row<dimL(site); row++){
         for (int col=0; col<dimR(site); col++){
            MPOprefactors[site][row][col] = &fZero;
            MPOoperators[site][row][col] = opZero;
         }
      }
      MPOprefactors[site][0][0] = &fOne;
      MPOoperators[site][0][0] = opOne;
      MPOprefactors[site][dimL(site)-1][0] = &MinusH;
      MPOoperators[site][dimL(site)-1][0] = opSz;
      MPOprefactors[site][1][0] = &fOne;
      MPOoperators[ site][1][0] = (useLadder) ? (Operator *) opSminus : (Operator *) opSx;
      MPOprefactors[site][2][0] = &fOne;
      MPOoperators[ site][2][0] = (useLadder) ? (Operator *) opSplus : (Operator *) opISy;
      MPOprefactors[site][3][0] = &fOne;
      MPOoperators[ site][3][0] = opSz;
      for (int col=1; col<dimR(site); col++){
         MPOprefactors[site][3 + col][col] = &fOne;
         MPOoperators[ site][3 + col][col] = opOne;
      }
      for (int col=0; col<length-1-site; col++){
         double * value = gCouplingPointer(site,site+1+col);
         MPOprefactors[site][dimL(site)-1][1 + 3*col + 0] = value;
         MPOoperators[ site][dimL(site)-1][1 + 3*col + 0] = (useLadder) ? (Operator *) opSplus : (Operator *) opSx;
         MPOprefactors[site][dimL(site)-1][1 + 3*col + 1] = (useLadder) ? value : gMinusCouplingPointer(site,site+1+col);
         MPOoperators[ site][dimL(site)-1][1 + 3*col + 1] = (useLadder) ? (Operator *) opSminus : (Operator *) opISy;
         MPOprefactors[site][dimL(site)-1][1 + 3*col + 2] = value;
         MPOoperators[ site][dimL(site)-1][1 + 3*col + 2] = opSz;
      }
   
   }

}

void HeisenbergMPO::sField(const double value){ MinusH = -value; }

double HeisenbergMPO::gField() const{ return -MinusH; }

void HeisenbergMPO::sCoupling(const int i, const int j, const double value){

   if ((i==j) || (i<0) || (i>=length) || (j<0) || (j>=length)) {
      cerr << "HeisenbergMPO::sCoupling  :  The combination of indices (" << i << "," << j << ") is not OK." << endl;
      return;
   }
   
   if (i<j){
      couplingMatrix[ (j*(j-1))/2 + i ]      = value;
      minusCouplingMatrix[ (j*(j-1))/2 + i ] = -value;
   } else {
      couplingMatrix[ (i*(i-1))/2 + j ]      = value;
      minusCouplingMatrix[ (i*(i-1))/2 + j ] = -value;
   }

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

double * HeisenbergMPO::gMinusCouplingPointer(const int i, const int j){

   if ((i==j) || (i<0) || (i>=length) || (j<0) || (j>=length)) {
      cerr << "HeisenbergMPO::gMinusCouplingPointer  :  The combination of indices (" << i << "," << j << ") is not OK." << endl;
      return NULL;
   }
   
   if (i<j){ return minusCouplingMatrix + (j*(j-1))/2 + i; }
   return minusCouplingMatrix + (i*(i-1))/2 + j;

}

int HeisenbergMPO::gNnonzeroPairs() const{

   int number = 0;
   for (int count=0; count<(length*(length-1))/2; count++){
      if (couplingMatrix[count] != 0.0){ number++; }
   }
   return number;

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
   if (useLadder){
      os << *opSplus;
      os << *opSminus;
   } else {
      os << *opSx;
      os << *opISy;
   }

}

void HeisenbergMPO::printOperatorTag(ostream& os, Operator * theOp) const{

   if (theOp == opZero){   os << "---"; }
   if (theOp == opOne){    os << "I  "; }
   if (theOp == opSz){     os << "Sz ";  }
   if (useLadder){
      if (theOp == opSplus){  os << "S+ "; }
      if (theOp == opSminus){ os << "S- "; }
   } else {
      if (theOp == opSx){  os << "Sx "; }
      if (theOp == opISy){ os << "iSy"; }
   }

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
   
   if (theMPO.gRN_nTerms()>0){
      os << "The MPO non-zero contributions :" << endl;
      for (int count=0; count<theMPO.gRN_nTerms(); count++){
         double value = theMPO.gRN_prefactor(count);
         if (value>=0.0){ os << "\t\tTerm " << count << "\t f =  " << value << "\t\t"; }
         else { os << "\t\tTerm " << count << "\t f = " << value << "\t\t"; }
         for (int site=0; site<theMPO.gLength(); site++){
            theMPO.printOperatorTag(os, theMPO.gRN_operator(count,site));
            os << " ";
         }
         os << "\t\tMy HC friend = " << theMPO.gRN_HCfriend(count) << endl;
      }
      os << "---------------------------------------------------------------------------------------" << endl;
   }

   return os;
   
}

void HeisenbergMPO::findNonZeroContributions(){

   //Erase a possible previous attempt
   if (RN_nTerms!=0){
      delete [] RN_prefactors;
      for (int count=0; count<RN_nTerms; count++){ delete [] RN_operators[count]; }
      delete [] RN_operators;
      delete [] RN_HCfriends;
      RN_nTerms = 0;
   }
   
   const bool debug_MPO = false;
   
   if (debug_MPO){
   
      cout << "This of course takes exponentially long..." << endl;
      //First find the number of non-zero contributions
      double * MPObonds = new double[length+1];
      for (int bond=0; bond<=length; bond++){ MPObonds[bond] = 0; }
      bool stop = false;
      while (!stop){
   
         bool allNonzero = true;
         for (int site=0; site<length; site++){
            if (gOperator(site, MPObonds[site], MPObonds[site+1])->AmIOp0()){ allNonzero = false; }
         }
         if (allNonzero){ RN_nTerms++; }
         int bondToAugment = length-1;
         bool wasAbleToAugment = false;
         while ((!wasAbleToAugment) && (bondToAugment>0)){
            if (MPObonds[bondToAugment] + 1 < dimL(bondToAugment)){
               wasAbleToAugment = true;
               MPObonds[bondToAugment] += 1;
            } else {
               MPObonds[bondToAugment] = 0;
               bondToAugment -= 1;
            }
         }
         if (bondToAugment == 0){ stop = true; }
         
      }
      cout << "Number of nonzero terms = " << RN_nTerms << endl;
      delete [] MPObonds;
      
   } else {
   
      for (int row=0; row<length-1; row++){
         for (int col=row+1; col<length; col++){
            if ( gCoupling(row,col) != 0.0 ){ RN_nTerms += 3; }
         }
      }
      
      if (gField() != 0.0){ RN_nTerms += length; }
   
   }
   
   //If necessary, allocate the necessary arrays
   if (RN_nTerms > 0){
      RN_prefactors = new double[RN_nTerms];
      RN_operators = new Operator**[RN_nTerms];
      for (int count=0; count<RN_nTerms; count++){
         RN_operators[count] = new Operator*[length];
      }
      RN_HCfriends = new int[RN_nTerms];
      for (int cnt=0; cnt<RN_nTerms; cnt++){ RN_HCfriends[cnt] = -1; }
   }
   
   //Then fill the operators and the prefactors
   if (RN_nTerms>0){
      if (debug_MPO){
   
         double * MPObonds = new double[length+1];
         int theCurrentTerm = 0;
         for (int bond=0; bond<=length; bond++){ MPObonds[bond] = 0; }
         bool stop = false;
         while (!stop){
            bool allNonzero = true;
            double factor = 1.0;
            for (int site=0; site<length; site++){
               if (gOperator(site, MPObonds[site], MPObonds[site+1])->AmIOp0() == true){ allNonzero = false; }
               factor *= gPrefactor(site, MPObonds[site], MPObonds[site+1]);
            }
            if (allNonzero){
               RN_prefactors[theCurrentTerm] = factor;
               for (int site=0; site<length; site++){
                  RN_operators[theCurrentTerm][site] = gOperator(site, MPObonds[site], MPObonds[site+1]);
               }
               theCurrentTerm++;
            }
            int bondToAugment = length-1;
            bool wasAbleToAugment = false;
            while ((!wasAbleToAugment) && (bondToAugment>0)){
               if (MPObonds[bondToAugment] + 1 < dimL(bondToAugment)){
                  wasAbleToAugment = true;
                  MPObonds[bondToAugment] += 1;
               } else {
                  MPObonds[bondToAugment] = 0;
                  bondToAugment -= 1;
               }
            }
            if (bondToAugment == 0){ stop = true; }
         }
         if (theCurrentTerm != RN_nTerms){ cerr << "HeisenbergMPO::findNonZeroContributions()  :  mismatch of the number of non-zero contributions" << endl; }
         delete [] MPObonds;
      
      } else {
      
         for (int count=0; count<RN_nTerms; count++){ for (int site=0; site<length; site++){ RN_operators[count][site] = opOne; } }
         int theCurrentTerm = 0;
         for (int row=0; row<length-1; row++){
            for (int col=row+1; col<length; col++){
               double value = gCoupling(row,col);
               if (value != 0.0){
                  if (useLadder){
                     RN_operators[theCurrentTerm  ][row] = opSplus;
                     RN_operators[theCurrentTerm  ][col] = opSminus;
                     RN_operators[theCurrentTerm+1][row] = opSminus;
                     RN_operators[theCurrentTerm+1][col] = opSplus;
                  } else {
                     RN_operators[theCurrentTerm  ][row] = opSx;
                     RN_operators[theCurrentTerm  ][col] = opSx;
                     RN_operators[theCurrentTerm+1][row] = opISy;
                     RN_operators[theCurrentTerm+1][col] = opISy;
                  }
                  RN_operators[theCurrentTerm+2][row] = opSz;
                  RN_operators[theCurrentTerm+2][col] = opSz;
                  RN_prefactors[theCurrentTerm  ] = value;
                  RN_prefactors[theCurrentTerm+1] = ((useLadder) ? 1 : -1) * value;
                  RN_prefactors[theCurrentTerm+2] = value;
                  theCurrentTerm += 3;
               }
            }
         }
         if (gField() != 0.0){
            for (int site=0; site<length; site++){
               RN_operators[theCurrentTerm][site] = opSz;
               RN_prefactors[theCurrentTerm] = MinusH;
               theCurrentTerm += 1;
            }
         }
         if (theCurrentTerm != RN_nTerms){ cerr << "HeisenbergMPO::findNonZeroContributions()  :  mismatch of the number of non-zero contributions" << endl; }
         
      }
   }
   
   //Then find the HC friends
   if (useLadder){
      Operator ** whatItShouldBe = new Operator*[length];
      for (int index=0; index<RN_nTerms; index++){
         if (RN_HCfriends[index] == -1){
      
            //Find what the HC looks like
            bool isDifferent = false;
            for (int site=0; site<length; site++){
               whatItShouldBe[site] = RN_operators[index][site];
               if (whatItShouldBe[site] == opSplus){
                  whatItShouldBe[site] = opSminus;
                  isDifferent = true;
               } else{
                  if (whatItShouldBe[site] == opSminus){
                     whatItShouldBe[site] = opSplus;
                     isDifferent = true;
                  }
               }
            }
            
            //Find the HC
            if (isDifferent==false){
               RN_HCfriends[index] = index;
            } else {
               const double desiredPrefactor = RN_prefactors[index];
               for (int index2=index+1; index2<RN_nTerms; index2++){
                  if (desiredPrefactor == RN_prefactors[index2]){
                     bool allTheSame = true;
                     for (int site=0; site<length; site++){
                        if (whatItShouldBe[site] != RN_operators[index2][site]){
                           allTheSame = false;
                           site = length;
                        }
                     }
                     if (allTheSame){
                        RN_HCfriends[index] = index2;
                        RN_HCfriends[index2] = index;
                        index2 = RN_nTerms;
                     }
                  }
               }
            }
            
         }
      }
      delete [] whatItShouldBe;
   } else {
      for (int index=0; index<RN_nTerms; index++){ RN_HCfriends[index] = index; }
   }

}

