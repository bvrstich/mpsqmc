#ifndef WORKSPACE_H
#define WORKSPACE_H

#include <complex>

using std::complex;

/**
 * class to organize the workspaces needed by the lapack functions in MPSstate
 */
class WorkSpace{

   public:
   
      //Constructor
      WorkSpace(int D,int DO,int d);
      
      //Destructor
      virtual ~WorkSpace();

      void clean(int rank);

      complex<double> **work1;
      complex<double> **work2;
      complex<double> **work3;

   private:
   
      //Number of OpenMP threads
      int num_omp_threads;
 
      int dim1;
      int dim2;
      int dim3;

     
};

#endif
