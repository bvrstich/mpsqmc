#ifndef MPSTENSOR_H
#define MPSTENSOR_H

#include <iostream>
#include <fstream>
#include <complex>

using std::ostream;
using std::ifstream;
using std::complex;

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 9, 2013 */

#include "Random.h"

class MPStensor{

   friend ostream &operator<<(ostream &output,const MPStensor &tensor);

   public:
   
      //Constructor
      MPStensor(const int dimL, const int dimR, const int phys_d, Random * RN);
      
      //Copy constructor
      MPStensor(MPStensor * toCopy);

      MPStensor(const char *filename,Random *RN);
      
      //Destructor
      ~MPStensor();
      
      //Get the left virtual dimension
      int gDimL() const;
      
      //Get the right virtual dimension
      int gDimR() const;
      
      //Get the size of the local Hilbert space
      int gPhys_d() const;
      
      //Get storage size
      int gStorageSize() const;
      
      //Get the pointer to the total storage
      complex<double> * gStorage();
      
      //Get the pointer to the block corresponding to a specific index of the local Hilbert space
      complex<double> * gStorage(const int d_val);

      complex<double> &operator()(int s,int i,int j) const;
      
      //Get the random number generator
      Random * gRN();
      
      //Left-normalization
      void QR(complex<double> * Rmx, complex<double> * mem, complex<double> * tau, complex<double> * work);
      
      //Right-normalization
      void LQ(complex<double> * Lmx, complex<double> * tau, complex<double> * work);
      
      //Multiply at the left with Lmx
      void LeftMultiply(complex<double> * Lmx, complex<double> * work);
      
      //Multiply at the right with Rmx
      void RightMultiply(complex<double> * Rmx, complex<double> * work);
      
      //Reset the virtual dimensions; the storage is only reset if the required storage size is larger than the current storage size
      void Reset(const int dimL, const int dimR);
     
      void copy(MPStensor *toCopy);
 
   private:
   
      //The random number generator
      Random * RN;
      
      //The left bond dimension
      int dimL;
      
      //The right bond dimension
      int dimR;
      
      //The size of the local Hilbert space
      int phys_d;
      
      //The storage size (so you don't always have to recreate this)
      int storageSize;
   
      //Storage for the MPS variables
      complex<double> * storage;
      
      //Fill storage with random complex numbers 0 <= val < 1
      void random();
      
};

#endif
