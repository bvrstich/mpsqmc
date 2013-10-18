#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "Walker.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 16, 2013 */

Walker::Walker(MPSstate * theState, const double overlap, const double weight, const double energyHistory){

   this->theState      = new MPSstate(theState);
   this->weight        = weight;
   this->overlap       = overlap;
   this->energyHistory = energyHistory;

}

Walker::Walker(Walker * theWalker){

   this->theState      = new MPSstate(theWalker->gState());
   this->weight        = theWalker->gWeight();
   this->overlap       = theWalker->gOverlap();
   this->energyHistory = theWalker->gEnergyHistory();

}

Walker::Walker(const int length, const int Dtrunc, const int phys_d, Random * RN){

   this->theState      = new MPSstate(length, Dtrunc, phys_d, RN);
   this->weight        = 1.0;
   this->overlap       = NAN;
   this->energyHistory = 0.0;

}

Walker::~Walker(){

   delete theState;

}

double Walker::gWeight() const{ return weight; }

double Walker::gOverlap() const{ return overlap; }

double Walker::gEnergyHistory() const{ return energyHistory; }

MPSstate * Walker::gState(){ return theState; }

void Walker::setWeight(const double weight){ this->weight = weight; }

void Walker::multiplyWeight(const double factor){ this->weight *= factor; }

void Walker::setOverlap(MPSstate * theTrial){

   double val = theState->LeftNormalize();
   val        = theTrial->InnerProduct(theState);
   
   if (val<0.0){ //Due to IS, only overlaps>0 are selected.
      theState->ChangePhase();
      val *= -1;
   }
   
   this->overlap = val;

}

void Walker::setEnergyHistory(const double energyHistory){ this->energyHistory = energyHistory; }

void Walker::addEnergyToHistory(const double term){ this->energyHistory += term; }

double Walker::calcIndividualProjectedEnergy(MPSstate * HamTimesTrial){

   //Set overlap must be called first !!!!
   return HamTimesTrial->InnerProduct(theState) / overlap;

}


