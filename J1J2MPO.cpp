#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "J1J2MPO.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 29, 2013 */

J1J2MPO::J1J2MPO(int L, int phys_d,double J2) : MPO(){

   int length = L*L;

   if (length < 3) 
      cerr << "J1J2MPO::J1J2MPO  :  The length must be at least three for the program to work." << endl;

   this->length = length;
   this->phys_d = phys_d;
   this->opZero   = new Op0(phys_d);
   this->opOne    = new OpI(phys_d);

   this->opSx  = new OpSx(phys_d);
   this->opSy = new OpSy(phys_d);
   this->opSz     = new OpSz(phys_d);

   fOne = complex<double>(1.0,0.0);
   fZero = complex<double>(0.0,0.0);

   J = new complex<double> [length * length];

   for(int i = 0;i < length;i++)
      for(int j = 0;j < length;j++)
         J[j*length + i] = 0.0;

   for (int row=0; row< L; row++)
      for (int col=0; col< L; col++){

         int number = row + L * col;

         int neighbour1 = (row + 1       )%L + L * col;
         int neighbour2 = (row - 1 + L)%L + L * col;
         int neighbour3 = row                   + L * ((col - 1 + L)%L);
         int neighbour4 = row                   + L * ((col + 1       )%L);

         J[number*length + neighbour1] = 1.0;
         J[number*length + neighbour2] = 1.0;
         J[number*length + neighbour3] = 1.0;
         J[number*length + neighbour4] = 1.0;

         int neighbour5 = (row + 1       )%L + L * ((col + 1       )%L);
         int neighbour6 = (row + 1       )%L + L * ((col - 1 + L)%L);
         int neighbour7 = (row - 1 + L)%L + L * ((col - 1 + L)%L);
         int neighbour8 = (row - 1 + L)%L + L * ((col + 1       )%L);

         J[number*length + neighbour5] = J2;
         J[number*length + neighbour6] = J2;
         J[number*length + neighbour7] = J2;
         J[number*length + neighbour8] = J2;

      }
 
   MPOdimensions = new int [length+1];
   fillMPOdimensions();
   
   MPOprefactors = new complex<double> *** [length];
   MPOoperators  = new Operator *** [length];

   for(int site = 0;site < length;site++){

      MPOprefactors[site] = new complex<double> ** [dimL(site)];
      MPOoperators[site]  = new Operator**[dimL(site)];

      for(int row=0; row<dimL(site); row++){

         MPOprefactors[site][row] = new complex<double> * [dimR(site)];
         MPOoperators[site][row]  = new Operator * [dimR(site)];

      }

   }

   fillMPOprefactorsMPOoperators();

}

J1J2MPO::~J1J2MPO(){

   delete opZero;
   delete opOne;
   delete opSz;
   delete opSx;
   delete opSy;
   
   delete [] J;
   
   for (int site = 0;site < length;site++){

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

void J1J2MPO::fillMPOdimensions(){

   int locSwitch = length/2;

   MPOdimensions[0] = MPOdimensions[length] = 1;

   for (int site=0; site< locSwitch; site++)
      MPOdimensions[site+1] = 2 + 3*(site+1);

   for (int site=length-1; site>locSwitch; site--)
      MPOdimensions[site]   = 2 + 3*(length-site);

   Dtrunc = 1;

   for(int site = 1;site < length;++site)
      if(Dtrunc < MPOdimensions[site])
         Dtrunc = MPOdimensions[site];

}

void J1J2MPO::fillMPOprefactorsMPOoperators(){

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
   
      for (int row=0; row<dimL(site); row++)
         for (int col=0; col<dimR(site); col++){

            MPOprefactors[site][row][col] = &fZero;
            MPOoperators[site][row][col] = opZero;
         
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

      for(int row=0; row<site; row++){

         complex<double> * value = gCouplingPointer(row,site);

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

         complex<double> * value = gCouplingPointer(row,locSwitch);

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

         complex<double> * value = gCouplingPointer(locSwitch,locSwitch+1+col);

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

         complex<double> * value = gCouplingPointer(site,site+1+col);

         MPOprefactors[site][dimL(site)-1][1 + 3*col + 0] = value;
         MPOoperators[ site][dimL(site)-1][1 + 3*col + 0] = opSx;

         MPOprefactors[site][dimL(site)-1][1 + 3*col + 1] = value;
         MPOoperators[ site][dimL(site)-1][1 + 3*col + 1] = opSy;

         MPOprefactors[site][dimL(site)-1][1 + 3*col + 2] = value;
         MPOoperators[ site][dimL(site)-1][1 + 3*col + 2] = opSz;

      }
   
   }

}

complex<double> J1J2MPO::gCoupling(const int i, const int j) const{

   return J[j*length + i ];

}

complex<double> *J1J2MPO::gCouplingPointer(const int i, const int j){

   return J + j*length + i;

}

complex<double> J1J2MPO::gPrefactor(const int site, const int row, const int col) const{

   if ((site<0) || (site>=length) || (row<0) || (row>=dimL(site)) || (col<0) || (col>=dimR(site)))
      return 0.0; 

   return *MPOprefactors[site][row][col];

}

Operator * J1J2MPO::gOperator(const int site, const int row, const int col) const{

   if ((site<0) || (site>=length) || (row<0) || (row>=dimL(site)) || (col<0) || (col>=dimR(site)))
      return NULL;

   if(abs(gPrefactor(site,row,col)) < 1.0e-15) 
      return opZero;

   return MPOoperators[site][row][col];

}
