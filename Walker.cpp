#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <algorithm>

#include "Walker.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 16, 2013 */

Walker::Walker(MPSstate * theState, const double weight){

   this->theState      = new MPSstate(theState);
   this->weight        = weight;

   VL = new complex<double> [3*theState->gLength()];

}

Walker::Walker(Walker * theWalker){

   this->theState      = new MPSstate(theWalker->gState());
   this->weight        = theWalker->gWeight();
   this->overlap       = theWalker->gOverlap();

   this->EL = theWalker->gEL();

   VL = new complex<double> [3*theState->gLength()];

   for(int r = 0;r < 3;++r)
      for(int k = 0;k < theState->gLength();++k)
         VL[r*theState->gLength() + k] = theWalker->gVL(k,r);

}

Walker::Walker(const int length, const int Dtrunc, const int phys_d, Random * RN){

   this->theState      = new MPSstate(length, Dtrunc, phys_d, RN);
   this->weight        = 1.0;

   VL = new complex<double> [3*length];

}

Walker::~Walker(){

   delete theState;

   delete [] VL;

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
   
   return VL[r*theState->gLength() + k]; 
   
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
void Walker::update_weight(double dtau,MPSstate *Psi0,MPSstate *HPsi0){

   complex<double> prev_over = overlap;

   theState->normalize();
   overlap = theState->InnerProduct(Psi0);

   //energy reset and scaling
   complex<double> prev_EL = EL;

   EL = theState->InnerProduct(HPsi0)/overlap;

   prev_EL = 0.5 * ( prev_EL + EL );

   double exponent = -dtau * std::real(prev_EL);

   weight *= exp(exponent);

   //complex weight rescaling
   double theta = std::arg(overlap/prev_over);

   weight *= std::max(0.0,cos(theta));

}

/**
 * calculate the overlap with the trial, Psi0
 */
void Walker::sOverlap(MPSstate * Psi0){

   theState->normalize();
   overlap = theState->InnerProduct(Psi0);

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
      for(int k = 0;k < theState->gLength();++k)
         VL[r*theState->gLength() + k] = theState->InnerProduct(VPsi0[r*theState->gLength() + k])/overlap;

}
