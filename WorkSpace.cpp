#include <omp.h>
#include <iostream>
#include "WorkSpace.h"

/**
 * allocate as many workspaces as needed by OMP
 * @param D truncation dimension of largest bond MPS
 * @param DO truncation dimension of MPO
 * @param d physical dimension
 */
WorkSpace::WorkSpace(int D,int DO,int d){
   
   dim1 = D*D*DO*d;
   dim2 = D*D*DO*d;
   dim3 = D*D*DO*d;

   work1 = new complex<double> [dim1];
   work2 = new complex<double> [dim2];
   work3 = new complex<double> [dim3];

}

WorkSpace::~WorkSpace(){

   delete [] work1;
   delete [] work2;
   delete [] work3;

}

/**
 * set the workspaces on rank 'rank' to zero
 */
void WorkSpace::clean(){

   for(int i = 0;i < dim1;++i)
      work1[i] = 0.0;

   for(int i = 0;i < dim2;++i)
      work2[i] = 0.0;

   for(int i = 0;i < dim3;++i)
      work3[i] = 0.0;

}
