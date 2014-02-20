#include "MPItag.h" //Here USE_MPI_IN_MPSQMC can be defined to switch to non-MPI compilers.

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <omp.h>
//#include </opt/intel/composer_xe_2013.0.079/compiler/include/omp.h>
#include <time.h>

#include "Random.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on September 4, 2013 */

using namespace std;

Random::Random(){
   
   #ifdef USE_MPI_IN_MPSQMC
      rank = MPI::COMM_WORLD.Get_rank();
   #else
      rank = 0;
   #endif
   
   #ifdef _OPENMP
      num_omp_threads = omp_get_max_threads();
   #else
      num_omp_threads = 1;
   #endif
   
   gen = new boost::random::mt19937 [num_omp_threads];

   for (int cnt=0; cnt<num_omp_threads; cnt++)
      gen[cnt].seed( time(0) + cnt*cnt*23 + 3571*rank );
   
   dists = new boost::random::uniform_real_distribution<double> * [num_omp_threads];

   for (int cnt=0; cnt<num_omp_threads; cnt++)
      dists[cnt] = new boost::random::uniform_real_distribution<double> (0,1); 

   gauss = new boost::random::normal_distribution<double> * [num_omp_threads];

   for (int cnt=0; cnt<num_omp_threads; cnt++)
      gauss[cnt] = new boost::random::normal_distribution<double> (0.0,1.0); 

}

Random::~Random(){

   for (int cnt=0; cnt<num_omp_threads; cnt++)
      delete dists[cnt];

   delete [] dists;

   for (int cnt=0; cnt<num_omp_threads; cnt++)
      delete gauss[cnt];

   delete [] gauss;

   delete [] gen;
   
}

double Random::rand(){

   #ifdef _OPENMP
      int tid = omp_get_thread_num();
   #else
      int tid = 0;
   #endif

   return (*dists[tid])(gen[tid]);

}

double Random::normal(){

   #ifdef _OPENMP
      int tid = omp_get_thread_num();
   #else
      int tid = 0;
   #endif

   return (*gauss[tid])(gen[tid]);

}


void Random::test(){

   sleep(3*rank);

   int factor = 10;
   double * RNs = new double[factor * num_omp_threads];
   int * tIDs = new int[factor * num_omp_threads];

   #pragma omp parallel for schedule(static) default(none) shared(RNs, tIDs, factor)
   for (int cnt=0; cnt<factor*num_omp_threads; cnt++){
   
      #ifdef _OPENMP
         int tid = omp_get_thread_num();
      #else
         int tid = 0;
      #endif
      
      RNs[cnt] = rand();
      tIDs[cnt] = tid;
   
   }
   
   for (int cnt=0; cnt<factor*num_omp_threads; cnt++){
      cout << "tID = " << tIDs[cnt] << " and RN = " << RNs[cnt] << endl;
   }
   
   delete [] RNs;
   delete [] tIDs;

}
