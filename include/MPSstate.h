#ifndef MPSSTATE_H
#define MPSSTATE_H

#include <iostream>
#include <fstream>
#include <complex>

using std::ostream;
using std::ifstream;
using std::complex;

#include "MPStensor.h"
#include "MPO.h"
#include "TwoSiteObject.h"
#include "TrotterJ1J2.h"

#include "WorkSpace.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 9, 2013 */

class MPSstate{

   friend ostream &operator<<(ostream &output,MPSstate &mps);

   public:
   
      //Constructor
      MPSstate(const int length, const int Dtrunc, const int phys_d, Random * RN);
      
      //Second constructor: with given virtual dimension array
      MPSstate(const int length, const int Dtrunc, const int phys_d, int * VirtualDims, Random * RN);
      
      //Copy contructor
      MPSstate(MPSstate * toCopy);

      MPSstate(const char *filename,Random *RN);
 
      //Destructor
      ~MPSstate();
      
      //Get the virtual dimension truncation
      int gDtrunc() const;
      
      //Get the virtual dimension at boundary i: MPS chain [0 1 2 3] has boundaries [0 1 2 3 4]. Hence tensor i -> left bound i, right bound i+1.
      int gDimAtBound(const int bound) const;

      void printVdim() const;
      
      //Get the size of the local Hilbert space
      int gPhys_d() const;
      
      //Get the chain length
      int gLength() const;
      
      //Get the pointer to the MPS tensor at site
      MPStensor * gMPStensor(const int site);

      MPStensor &operator[](int site);
      
      //Get the random number generator
      Random * gRN();
      
      //Left-normalize the whole chain, return the norm
      complex<double> LeftNormalize();
      
      //Right-normalize the whole chain, return the norm
      complex<double> RightNormalize();

      //normalize the state, no canonicalization. return the norm
      complex<double> normalize();
      
      //Calculate the overlap
      complex<double> InnerProduct(MPSstate * OtherState);

      //Calculate <phi|O|psi>
      complex<double> expectation(MPO *,MPSstate * OtherState);

      void ChangePhase();
      
      //Multiply with a scalar
      void ScalarMultiplication(const complex<double> factor);
      
      //Compress the state so that max. truncD virtual states remain.
      void CompressState(const int truncD);
      
      //Compress lossless
      void CompressState(){

         CompressState(Dtrunc);

      }
      
      //Do H * Psi0 and store in this object
      void ApplyMPO(bool conj,MPO * theMPO, MPSstate * Psi0);
      
      //apply the auxiliary field
      void ApplyAF(int k,int r,complex<double> x,TrotterJ1J2 * theTrotter);

      void ApplyAF(int k,complex<double> x,TrotterJ1J2 * theTrotter);
      
      static void InitWork(int D,int DO,int d);
      static void ClearWork();
      
   private:
   
      //The random number generator
      Random * RN;
   
      //The chain length
      int length;
      
      //The bond dimension truncation
      int Dtrunc;
      
      //The physical dimension (local Hilbert space size)
      int phys_d;
      
      //The MPStensor array
      MPStensor ** theTensors;
   
      //The virtual dimensions in the chain (array with length "this->length"+1)
      int * VirtualD;
      
      //For compression: is the two-site object allocated?
      bool TwoSiteObjectAllocated;
      
      //For compression: the two-site object
      TwoSiteObject * the2siteObject;
      
      static WorkSpace *ws;
      
};

#endif
