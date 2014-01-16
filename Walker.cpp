#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <algorithm>

#include "Walker.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 16, 2013 */

Walker::Walker(MPSstate * theState, const double weight,int n_trot){

   this->theState      = new MPSstate(theState);
   this->weight        = weight;
   this->n_trot        = n_trot;

   VL = new complex<double> [3*n_trot];

}

Walker::Walker(Walker * theWalker){

   this->theState      = new MPSstate(theWalker->gState());
   this->weight        = theWalker->gWeight();
   this->overlap       = theWalker->gOverlap();
   this->n_trot        = theWalker->gn_trot();

   this->EL = theWalker->gEL();

   VL = new complex<double> [3*n_trot];

   for(int r = 0;r < 3;++r)
      for(int k = 0;k < n_trot;++k)
         VL[r*n_trot + k] = theWalker->gVL(k,r);

}

Walker::Walker(const int length, const int Dtrunc, const int phys_d, int n_trot,Random * RN){

   this->theState      = new MPSstate(length, Dtrunc, phys_d, RN);
   this->weight        = 1.0;
   this->n_trot        = n_trot;

   VL = new complex<double> [3*n_trot];

}

Walker::~Walker(){

   delete theState;

   delete [] VL;

}

/** 
 * @return the number of trotter terms
 */
int Walker::gn_trot() const {

   return n_trot;

}

double Walker::gWeight() const{
   
   return weight; 
   
}

complex<double> Walker::gOverlap() const{
   
   return overlap; 
   
}

complex<double> Walker::gEL() const{
   
   return EL; 
   
}

complex<double> Walker::gVL(int k,int r) const{
   
   return VL[r*n_trot + k]; 
   
}

MPSstate * Walker::gState(){
   
   return theState; 
   
}

void Walker::sWeight(const double weight){
   
   this->weight = weight; 
   
}

void Walker::multWeight(const double factor){
   
   this->weight *= factor; 
   
}

/**
 * update the weight of the walker and set the new overlap and local energy.
 * @param dtau timestep
 */
void Walker::update_weight(MPSstate *Psi0,complex<double> x,complex<double> shift){

   complex<double> prev_overlap = overlap;

   theState->normalize();
   complex<double> val = theState->InnerProduct(Psi0);

   if (std::real(val) < 0.0){

      theState->ChangePhase();
      val *= -1;

   }

   this->overlap = val;

   complex<double> scal = overlap/prev_overlap * exp(x*shift - 0.5 * shift*shift);

   weight *= std::abs(scal);

}

/**
 * calculate the overlap with the trial, Psi0
 */
void Walker::sOverlap(MPSstate * Psi0){

   complex<double> norm = theState->normalize();

   weight *= std::real(norm);

   complex<double> val = theState->InnerProduct(Psi0);
/*
   if (std::real(val) < 0.0){

      theState->ChangePhase();
      val *= -1;

   }
*/
   this->overlap = val;

}

/** 
 * set the Local Energy: overlap has to be set first!
 */
void Walker::sEL(MPSstate * HPsi0){

   EL = theState->InnerProduct(HPsi0)/overlap;

}

/** 
 * set the v operator Energy: overlap has to be set first!
 */
void Walker::sVL(MPSstate ** VPsi0){

   for(int r = 0;r < 3;++r)
      for(int k = 0;k < n_trot;++k)
         VL[r*n_trot + k] = theState->InnerProduct(VPsi0[r*n_trot + k])/overlap;

}
