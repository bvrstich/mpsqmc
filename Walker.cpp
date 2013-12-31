#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "Walker.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 16, 2013 */

Walker::Walker(MPSstate * theState, const complex<double> overlap, const double weight){

   this->theState      = new MPSstate(theState);
   this->weight        = weight;
   this->overlap       = overlap;

}

Walker::Walker(Walker * theWalker){

   this->theState      = new MPSstate(theWalker->gState());
   this->weight        = theWalker->gWeight();
   this->overlap       = theWalker->gOverlap();

}

Walker::Walker(const int length, const int Dtrunc, const int phys_d, Random * RN){

   this->theState      = new MPSstate(length, Dtrunc, phys_d, RN);
   this->weight        = 1.0;
   this->overlap       = NAN;

}

Walker::~Walker(){

   delete theState;

}

double Walker::gWeight() const{
   
   return weight; 
   
}

complex<double> Walker::gOverlap() const{
   
   return overlap; 
   
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
 * calculate the overlap with the trial, Psi0
 */
void Walker::sOverlap(MPSstate * Psi0){

   overlap = theState->InnerProduct(Psi0);

}

/** 
 * set the Local Energy: overlap has to be set first!
 */
void Walker::sEL(MPSstate * HPsi0){

   EL = theState->InnerProduct(HPsi0)/overlap;

}
