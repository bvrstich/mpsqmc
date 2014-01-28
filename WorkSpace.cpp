#include <omp.h>

#include "WorkSpace.h"

/**
 * allocate as many workspaces as needed by OMP
 * @param D truncation dimension of largest bond MPS
 * @param DO truncation dimension of MPO
 * @param d physical dimension
 */
WorkSpace::WorkSpace(int D,int DO,int d){
   
#ifdef _OPENMP
   num_omp_threads = omp_get_max_threads();
#else
   num_omp_threads = 1;
#endif

   work1 = new complex<double> * [num_omp_threads];
   work2 = new complex<double> * [num_omp_threads];
   work3 = new complex<double> * [num_omp_threads];

   dim1 = D*D*DO;
   dim2 = D*D*DO*d;
   dim3 = D*D*DO*d;

   for(int i = 0;i < num_omp_threads;++i){

      work1[i] = new complex<double> [dim1];
      work2[i] = new complex<double> [dim2];
      work3[i] = new complex<double> [dim3];

   }

}

WorkSpace::~WorkSpace(){

   for(int i = 0;i < num_omp_threads;++i){

      delete [] work1[i];
      delete [] work2[i];
      delete [] work3[i];

   }

   delete [] work1;
   delete [] work2;
   delete [] work3;

}

/**
 * set the workspaces on rank 'rank' to zero
 */
void WorkSpace::clean(int rank){

   for(int i = 0;i < dim1;++i)
      work1[rank][i] = 0.0;

   for(int i = 0;i < dim2;++i)
      work2[rank][i] = 0.0;

   for(int i = 0;i < dim3;++i)
      work3[rank][i] = 0.0;

}
