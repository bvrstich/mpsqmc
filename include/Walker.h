#ifndef WALKER_H
#define WALKER_H

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 16, 2013 */

#include "MPSstate.h"

class Walker{

   public:
   
      //Constructor copying an MPSstate
      Walker(MPSstate * theState, const double weight,int n_trot);
      
      //Constructor copying an entire Walker
      Walker(Walker * theWalker);
      
      //Constructor creating a random MPSstate, with weight 1, overlap NAN
      Walker(const int length, const int Dtrunc, const int phys_d, int n_trot,Random * RN);
      
      //Destructor
      ~Walker();
      
      //Return the walker weight
      double gWeight() const;
      
      //Return the walker overlap
      complex<double> gOverlap() const;

      //Return the pointer to the MPSstate
      MPSstate * gState();

      //Set the walker weight
      void sWeight(const double weight);

      //Multiply factor into the walker weight
      void multWeight(const double factor);
      
      //Set the overlap with the trial wfn
      void sOverlap(MPSstate *Psi0);

      void sEL(MPO *theMPO,MPSstate *Psi0);

      void sEL(MPSstate *HPsi0);

      void sEL(complex<double> );

      void sVL(MPSstate **VPsi0);

      void update_weight(MPSstate *Psi0,complex<double> x,complex<double> shift);

      int gn_trot() const;

      complex<double> gEL() const;

      complex<double> gVL(int,int) const;
      
   private:
   
      //The MPS state
      MPSstate *theState;

      //nr of trotter terms
      int n_trot;
      
      //The walker overlap with the trial wfn
      complex<double> overlap;

      //!local energy
      complex<double> EL;

      //!local auxiliary operators: <PsiT|v|phi>/<PsiT|phi>
      complex<double> *VL;
      
      //The walker weight
      double weight;
      
};

#endif
