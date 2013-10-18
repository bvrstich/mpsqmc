#ifndef WALKER_H
#define WALKER_H

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 16, 2013 */

#include "MPSstate.h"

class Walker{

   public:
   
      //Constructor copying an MPSstate
      Walker(MPSstate * theState, const double overlap, const double weight, const double energyHistory);
      
      //Constructor copying an entire Walker
      Walker(Walker * theWalker);
      
      //Constructor creating a random MPSstate, with weight 1, overlap NAN, and energy history 0
      Walker(const int length, const int Dtrunc, const int phys_d, Random * RN);
      
      //Destructor
      ~Walker();
      
      //Return the walker weight
      double gWeight() const;
      
      //Return the walker overlap
      double gOverlap() const;

      //Return the energy history
      double gEnergyHistory() const;

      //Return the pointer to the MPSstate
      MPSstate * gState();

      //Set the walker weight
      void setWeight(const double weight);

      //Multiply factor into the walker weight
      void multiplyWeight(const double factor);
      
      //Set the overlap with the trial wfn
      void setOverlap(MPSstate * theTrial);

      //Set the energy history
      void setEnergyHistory(const double energyHistory);

      //Add a term to the energyHistory
      void addEnergyToHistory(const double term);
      
      //Calculate the walker's individual projected energy. Note that setOverlap should have been adjusted correctly!!!!
      double calcIndividualProjectedEnergy(MPSstate * HamTimesTrial);
      
   private:
   
      //The MPS state
      MPSstate * theState;
      
      //The walker overlap with the trial wfn
      double overlap;
      
      //The walker weight
      double weight;
      
      //The walker energy history: for all previous steps, contains the sum of the walker energies
      double energyHistory;

};

#endif
