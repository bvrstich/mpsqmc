#ifndef DMRG_H
#define DMRG_H

#include "MPStensor.h"
#include "MPSstate.h"
#include "MPO.h"
#include "TwoSiteObject.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 16, 2013 */

class DMRG{

   public:
   
      //Constructor
      DMRG(MPSstate * MPStrial, MPO * theMPO);
      
      //Destructor
      ~DMRG();
      
      //Solver: puts the lowest energy MPS with same bond dimensions as MPStrial in MPStrial; returns the energy
      double Solve();
      
   private:
   
      //The length of the MPS chain
      int length;
   
      //Partially contracted <Psi | MPO | Psi>
      double ** boundaryTerms;
      
      //Pointer to the given MPStrial
      MPSstate * MPSsolution;
      
      //Pointer to the given MPO
      MPO * theMPO;
      
      //The MPO * MPS --> can use the workspaces of this object to construct the boundaryTerms;
      MPSstate * MPOMPS;
      
      //The two-site object : for two-site effective Hamiltonian eigenvalue problem
      TwoSiteObject * the2siteObject;
      
      //Construct the left boundary operator at bound from the previous one
      void ConstructBoundaryTermMovingRight(const int bound);
      
      //Construct the right boundary operator at bound from the previous one
      void ConstructBoundaryTermMovingLeft(const int bound);
      
      //Do a left DMRG sweep
      double sweepleft();
      
      //Do a right DMRG sweep
      double sweepright();
      
      //Solve the local Hermitian two-site standard eigenvalue problem
      /*double SolveHeff(const int index);*/
      
      //Solve with Davidson's algorithm
      double SolveDAVIDSON(const int index);
      
      //Do the matrix-vector multiplication
      void matvec(double * vec, double * resultvec, const int index);
      
      //Fill the diagonal preconditioner for Davidson's algo
      void fillDiag(double * diag, const int index);
      
      //SolveHeff: are the arrays allocated
      bool SolveHeffAllocated;
      
      //SolveHeff: vectorsize
      int SolveHeffVectorSize;
      
      //SolveHeff: which eigenvalues
      char * which;
      
      //SolveHeff: residual
      double * resid;
      
      //SolveHeff: Arnoldi vectors
      double * v;
      
      //SolveHeff: ipntr
      int * ipntr;
      
      //SolveHeff: iparam
      int * iparam;
      
      //SolveHeff: workd
      double * workd;
      
      //SolveHeff: workl
      double * workl;
      
};

#endif
