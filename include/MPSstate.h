#ifndef MPSSTATE_H
#define MPSSTATE_H

#include "MPStensor.h"
#include "MPO.h"
#include "TwoSiteObject.h"
#include "TrotterHeisenberg.h"
#include "GridGenerator.h"

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
      
      //Get the size of the local Hilbert space
      int gPhys_d() const;
      
      //Get the chain length
      int gLength() const;
      
      //Get the pointer to the MPS tensor at site
      MPStensor * gMPStensor(const int site);
      
      //Get the random number generator
      Random * gRN();
      
      //Left-normalize the whole chain, return the norm
      double LeftNormalize();
      
      //Right-normalize the whole chain, return the norm
      double RightNormalize();
      
      //Calculate the overlap
      double InnerProduct(MPSstate * OtherState);
      
      //Multiply with a scalar
      void ScalarMultiplication(const double factor);
      
      //Together with CompressState and ScalarMultiplication, this provides the functionality for imaginary time evolution
      void ResetContentsAndStoreSumOf(MPSstate * state1, MPSstate * state2);
      
      //Change the phase
      void ChangePhase();
      
      //Compress the state so that max. truncD virtual states remain.
      void CompressState(const int truncD);
      
      //Compress lossless
      void CompressState(){ CompressState(Dtrunc); }
      
      //Do H * Psi0 and store in this object
      void ApplyMPO(MPO * theMPO, MPSstate * Psi0);
      
      //Apply a specific non-zero term of the MPO to this wfn
      void ApplyMPOterm(MPO * theMPO, const int SelectedTerm);
      
      //Apply the hermitian conjugate of a specific non-zero term of the MPO to this wfn
      void ApplyMPOtermHC(MPO * theMPO, const int SelectedTerm){ ApplyMPOterm(theMPO, theMPO->gRN_HCfriend(SelectedTerm)); }
      
      //Apply a particular two-site Trotter term
      void ApplyTwoSiteTrotterTerm(TrotterHeisenberg * theTrotter, const int firstSite, const int secondSite, const int leftSVDindex, const int rightSVDindex, const bool doHC = false);
      
      //Apply a particular two-site Trotter term
      void ApplyTwoSiteTrotterTerm(TrotterHeisenberg * theTrotter, const int firstSite, const int secondSite, GridGenerator * theGrid, const int gridPoint);
      
      //Apply the single-site Trotter term on each site
      void ApplyOneSiteTrotterTermEverywhere(TrotterHeisenberg * theTrotter);
      
      //Check whether the work arrays are allocated with at least size size
      void checkWork1(const int size);
      void checkWork2(const int size);
      void checkWork3(const int size);
      
      //Get pointer to the work arrays
      double * gWork1(){ return work1; }
      double * gWork2(){ return work2; }
      double * gWork3(){ return work3; }
      
      
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
      
      //For QR etc: workspaces
      double * work1;
      double * work2;
      double * work3;
      int sizeWork1;
      int sizeWork2;
      int sizeWork3;
      bool work1Allocated;
      bool work2Allocated;
      bool work3Allocated;
      
};

#endif

