#ifndef MPSTENSOR_H
#define MPSTENSOR_H

#include <iostream>
#include <fstream>

using std::ostream;
using std::ifstream;

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
      double * gStorage();
      
      //Get the pointer to the block corresponding to a specific index of the local Hilbert space
      double * gStorage(const int d_val);
      
      //Get the random number generator
      Random * gRN();
      
      //Left-normalization
      void QR(double * Rmx, double * mem, double * tau, double * work);
      
      //Right-normalization
      void LQ(double * Lmx, double * tau, double * work);
      
      //Multiply at the left with Lmx
      void LeftMultiply(double * Lmx, double * work);
      
      //Multiply at the right with Rmx
      void RightMultiply(double * Rmx, double * work);
      
      //Reset the virtual dimensions; the storage is only reset if the required storage size is larger than the current storage size
      void Reset(const int dimL, const int dimR);
      
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
      double * storage;
      
      //Fill storage with random numbers 0 <= val < 1
      void random();
      
};

#endif
