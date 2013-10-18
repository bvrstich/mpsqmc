#ifndef GRIDGENERATOR_H
#define GRIDGENERATOR_H

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on October 15, 2013 */

using namespace std;

class GridGenerator{

   public:
   
      //Constructor
      GridGenerator(const int dim);
      
      //Destructor
      ~GridGenerator();
      
      //Get the Euclidean space dimension
      int gDim() const;
      
      //Get the number of grid points
      int gNPoints() const;

      //Get a specific coordinate of a specific point
      double gCoOfPoint(const int point, const int index) const;
      
      //Print its contents
      friend ostream& operator<<(ostream& os, const GridGenerator& theGrid);
      
      //Fill the grid with allowed points of the type +- sqrt( integer / denominator )
      void FillSimple(const int denominator);
      
      //Fill the grid according to the Marsaglia idea of uniform sampling on a sphere: pointsPerCo are selected uniformly on [0,1] and are transformed with the inverse cumulative distribution function of the basic gaussian to a "uniformly sampled" gaussian, after which they're used to generate points on the sphere
      void FillMarsaglia(const int pointsPerCo);
      
      
   private:
   
      //The dimension of the Euclidean space in which the hypersphere resides
      int dim;
      
      //Whether grid was already allocated or not
      bool gridFilled;
      
      //The amount of states in the grid
      int nStates;
      
      //The grid points: grid[ i + dim * k ] constains coordinate i of gridpoint k
      double * grid;
      
      //Given state, generate the new one for "FillSimple()". Returns false if no new state can be generated & generates only the ones with all + signs.
      bool GenerateNextSimple(int * state, const int denominator);
      
      
};

#endif
