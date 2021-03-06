#ifndef RANDOM_H
#define RANDOM_H

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on September 4, 2013 */

class Random{

   public:
   
      //Constructor
      Random();
      
      //Destructor
      ~Random();
      
      //Draw OpenMP and MPI thread safe random numbers from the uniform distribution (0,1)
      double rand();
      
      //Tester of the Random number generator
      void test();
      
   private:
   
      //Number of OpenMP threads
      int num_omp_threads;
      
      //My MPI rank
      int rank;
   
      //The Mersenne Twister RN generator
      boost::random::mt19937 *gen;
      
      //The uniform real distribution
      boost::random::uniform_real_distribution<double> **dists;
      
};

#endif
