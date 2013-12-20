#ifndef TWOSITEOBJECT_H
#define TWOSITEOBJECT_H

#include <iostream>
#include <fstream>
#include <complex>

using std::ostream;
using std::ifstream;
using std::complex;

#include "MPStensor.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 12, 2013 */

class TwoSiteObject{

   public:
   
      //Constructor
      TwoSiteObject(const int DimL, const int DimR, const int phys_d);
      
      //Destructor
      ~TwoSiteObject();
      
      //Get the size of the local Hilbert space
      int gPhys_d() const;
      
      //Get the left virtual dimension
      int gDimL() const;
      
      //Get the right virtual dimension
      int gDimR() const;
      
      //Fill the two-site object. In case the storage is not large enough, it is recreated.
      void Compose(MPStensor * MPSleft, MPStensor * MPSright);
      
      //Get the size of the storage
      int gStorageSize() const;
      
      //Get a block of the two-site object
      complex<double> * gStorage(const int d_left, const int d_right);
      
      //Get the full storage
      complex<double> * gStorage();
      
      //Decompose the two-site object
      int Decompose(MPStensor * MPSleft, MPStensor * MPSright, const int Dtrunc, bool movingright, bool possiblyCompress);
      
   private:
      
      //The physical dimension (local Hilbert space size)
      int phys_d;
      
      //The left virtual dimension
      int DimL;
      
      //The right virtual dimension
      int DimR;
      
      //The storage size
      int storageSize;
      
      //The storage space
      complex<double> * storage;
      
      //Work spaces as large as the storage space
      complex<double> * work_large;
      complex<double> * work_large2;
      complex<double> * work_large3;
      
      //Other SVD work spaces and their sizes
      int SVD_dimMin;
      int SVD_lwork;
      int * SVD_iwork; //8*SVD_dimMin
      complex<double> * SVD_work; //SVD_lwork
      double *SVD_Svalues; //SVD_dimMin
      
      //Schmidt values smaller than this number are always thrown out
      static const double SchmidtTrunc;
      
};

#endif
