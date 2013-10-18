#include <stdlib.h>
#include <math.h>
#include <iostream>

#include <boost/math/special_functions/erf.hpp>

#include "GridGenerator.h"

using namespace std;

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 15, 2013 */

GridGenerator::GridGenerator(const int dim){

   this->dim = dim;
   this->gridFilled = false;
   this->nStates = 0;

}

GridGenerator::~GridGenerator(){

   if (gridFilled){ delete [] grid; }

}

int GridGenerator::gDim() const{ return dim; }

int GridGenerator::gNPoints() const{ return nStates; }

double GridGenerator::gCoOfPoint(const int point, const int index) const{ return grid[ index + dim * point ]; }

ostream& operator<<(ostream& os, const GridGenerator& theGrid){

   const int dimension = theGrid.gDim();
   const int nStates   = theGrid.gNPoints();

   os << "GridGenerator print for Dim = " << dimension << endl;
   
   //Print the points, and calculate in the meanwhile the sum of the outer product of the vectors
   double * sumTest = new double[dimension*dimension];
   for (int cnt=0; cnt<dimension*dimension; cnt++){ sumTest[cnt] = 0.0; }
   
   os << "The points:" << endl;
   for (int cnt=0; cnt<nStates; cnt++){
      os << "      Point " << cnt+1 << "\t [";
      double sumSq = 0.0;
      for (int cnt2=0; cnt2<dimension; cnt2++){
         double val = theGrid.gCoOfPoint(cnt,cnt2);
         os << val << "\t";
         double prod = val*val;
         sumSq += prod;
         sumTest[ cnt2 * (1 + dimension) ] += prod;
         for (int cnt3=cnt2+1; cnt3<dimension; cnt3++){
            prod = val * theGrid.gCoOfPoint(cnt,cnt3);
            sumTest[ cnt2 + dimension * cnt3 ] += prod;
            sumTest[ cnt3 + dimension * cnt2 ] += prod;
         }
      }
      os << "] \t\t with 2-norm " << sqrt(sumSq) << endl;
   }
   
   //Compare the outer product of the vectors with nStates/dimension
   os << "#vectors / dim = " << ((double) nStates)/dimension << " vs. sum over the outer product of the vectors:" << endl;
   for (int row=0; row<dimension; row++){
      os << "      [ ";
      for (int col=0; col<dimension; col++){
         os << sumTest[ row + dimension * col ] << "\t";
      }
      os << "]" << endl;
   }
   
   delete [] sumTest;
   
   //Gather all different distances, how many times they occur, order them, and print them
   os << "Distances:" << endl;
   
   double * dists = new double[nStates * (nStates - 1)/2];
   int * counts  = new int[nStates * (nStates - 1)/2];
   int nPoss = 0;
   
   for (int cnt=0; cnt<nStates; cnt++){
      for (int cnt2=cnt+1; cnt2<nStates; cnt2++){
         double dist = 0.0;
         for (int co=0; co<dimension; co++){
            double diff = theGrid.gCoOfPoint(cnt,co) - theGrid.gCoOfPoint(cnt2,co);
            dist += diff * diff;
         }
         dist = sqrt(dist);
         bool isAlreadyIn = false;
         for (int cnt3=0; cnt3<nPoss; cnt3++){
            if (fabs(dists[cnt3] - dist) < 1e-8){
               isAlreadyIn = true;
               counts[cnt3] += 1;
               cnt3 = nPoss;
            }
         }
         if (!isAlreadyIn){
            dists[nPoss] = dist;
            counts[nPoss] = 1;
            nPoss++;
         }
      }
   }

   bool allOK = false;
   while (!allOK){
      allOK = true;
      for (int cnt=0; cnt<nPoss-1; cnt++){
         if ( dists[ cnt ] > dists[ cnt+1 ] ){
            allOK = false;
            double temp = dists[cnt];
            dists[cnt] = dists[cnt+1];
            dists[cnt+1] = temp;
            int temp2 = counts[cnt];
            counts[cnt] = counts[cnt+1];
            counts[cnt+1] = temp2;
         }
      }
   } // Now for all index: poss[ cnt ] <= poss [cnt+1 ]

   int totalCounts = 0;
   for (int cnt=0; cnt<nPoss; cnt++){
      os << "      Distance " << dists[cnt] << " has " << counts[cnt] << " counts." << endl;
      totalCounts += counts[cnt];
   }
   os << "Please note that 2 should be the max. distance, which occurs #vectors/2 = " << ((double)nStates)/2 << " number of times." << endl;
   os << "The total number of counts in the previous list is " << totalCounts << " and #vectors*(#vectors-1)/2 = " << nStates*(nStates-1)/2 << endl;
   delete [] dists;
   delete [] counts;

   return os;

}

void GridGenerator::FillMarsaglia(const int pointsPerCo){

   if (gridFilled){
      delete [] grid;
      gridFilled = false;
   }

   double * oneDimGrid = new double[pointsPerCo];
   for (int cnt=0; cnt<pointsPerCo; cnt++){ oneDimGrid[cnt] = 2*(cnt+0.5)/pointsPerCo - 1.0; }
   for (int cnt=0; cnt<pointsPerCo; cnt++){
      cout << "Point " << oneDimGrid[cnt] << " is transformed to ";
      oneDimGrid[cnt] = boost::math::erf_inv(oneDimGrid[cnt]);
      cout << oneDimGrid[cnt] << endl;
   }
   
   int statesEstimate = 1;
   for (int cnt=0; cnt<dim; cnt++){ statesEstimate *= pointsPerCo; } //Too much as [x x x x] generates the same point as [2x 2x 2x 2x]
   grid = new double[dim * statesEstimate];
   gridFilled = true;
   nStates = 0;
   
   double * attempt = new double[dim];
   for (int cnt=0; cnt<statesEstimate; cnt++){
      int countcopy = cnt;
      double norm = 0.0;
      for (int cnt2=0; cnt2<dim; cnt2++){
         int remainder = countcopy % pointsPerCo;
         attempt[cnt2] = oneDimGrid[remainder];
         norm += attempt[cnt2] * attempt[cnt2];
         countcopy = (countcopy - remainder)/pointsPerCo;
      }
      if (norm>1e-16){
         norm = 1.0/sqrt(norm);
         for (int cnt2=0; cnt2<dim; cnt2++){ attempt[cnt2] *= norm; }
         bool doCopy = true;
         for (int cntState=0; cntState<nStates; cntState++){
            double distSq = 0.0;
            for (int co=0; co<dim; co++){
               double diff = attempt[co] - grid[co + dim*cntState];
               distSq += diff*diff;
            }
            if (distSq<1e-16){
               doCopy = false;
               cntState = nStates;
            }
         }
         if (doCopy){
            for (int co=0; co<dim; co++){
               grid[co + dim*nStates] = attempt[co];
            }
            nStates++;
         }
      }
   }
   
   delete [] attempt;
   delete [] oneDimGrid;

}

void GridGenerator::FillSimple(const int denominator){

   if (gridFilled){
      delete [] grid;
      gridFilled = false;
   }

   int * currentState = new int[dim];
   
   //Calculate the number of states
   currentState[0] = denominator;
   for (int cnt=1; cnt<dim; cnt++){ currentState[cnt] = 0; }
   nStates = 2; //[+denom,0,0...] and [-denom,0,0...]
   bool couldGenerate = true;
   while (couldGenerate){
      couldGenerate = GenerateNextSimple(currentState, denominator);
      if (couldGenerate){
         int factor = 1;
         for (int cnt=0; cnt<dim; cnt++){
            if (currentState[cnt]>0){ factor *= 2; }
         }
         nStates += factor;
      }
   }
   
   //Allocate the grid
   grid = new double[dim * nStates];
   gridFilled = true;
   
   //Fill the grid
   int * nonZeros = new int[dim];
   currentState[0] = denominator;
   grid[ 0 + dim * 0 ] =  1.0;
   grid[ 0 + dim * 1 ] = -1.0;
   for (int cnt=1; cnt<dim; cnt++){
      currentState[cnt] = 0;
      grid[ cnt + dim * 0 ] = grid[ cnt + dim * 1 ] = 0.0;
   }
   nStates = 2;
   couldGenerate = true;
   while (couldGenerate){
      couldGenerate = GenerateNextSimple(currentState, denominator);
      if (couldGenerate){
         int factor = 1;
         int numNonZeros = 0;
         for (int cnt=0; cnt<dim; cnt++){
            if (currentState[cnt]>0){
               factor *= 2;
               nonZeros[numNonZeros] = cnt;
               numNonZeros += 1;
            }
         }
         for (int cnt=0; cnt<factor; cnt++){
            for (int cnt2=0; cnt2<dim; cnt2++){ grid[ cnt2 + dim * ( nStates + cnt ) ] = 0.0; }
            int countcopy = cnt;
            for (int cnt2=0; cnt2<numNonZeros; cnt2++){
               int remainder = countcopy % 2;
               grid[ nonZeros[cnt2] + dim * (nStates + cnt) ] = ((remainder == 0) ? 1 : -1) * sqrt( ( (double) currentState[nonZeros[cnt2]] ) / denominator );
               countcopy = (countcopy - remainder)/2;
            }
         }
         nStates += factor;
      }
   }
   
   delete [] nonZeros;
   delete [] currentState;

}

bool GridGenerator::GenerateNextSimple(int * state, const int denominator){

   if (state[dim-1] == denominator){ return false; }
   
   if (state[dim-1] == 0){
   
      int index = dim-2;
      while (state[index]==0){ index--; }
      state[index]   -= 1;
      state[index+1] += 1;
   
   } else {
   
      int index = dim-2;
      while (state[index]==0){ index--; }
      int lastsite    = state[dim-1];
      state[dim-1]   -= lastsite;
      state[index+1] += 1 + lastsite;
      state[index]   -= 1;
   
   }
   
   return true;

}

