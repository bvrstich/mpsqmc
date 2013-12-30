#include <stdlib.h>
#include "MPStensor.h"
#include "Lapack.h"
#include "Random.h"
#include <iostream>

using std::cout;
using std::endl;
using std::ostream;
using std::ifstream;


/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 9, 2013 */

MPStensor::MPStensor(const int dimL, const int dimR, const int phys_d, Random * RN){

   this->dimL = dimL;
   this->dimR = dimR;
   this->phys_d = phys_d;
   this->RN = RN;
   this->storageSize = dimL * dimR * phys_d;
   storage = new complex<double> [storageSize];
   random();

}

MPStensor::MPStensor(MPStensor * toCopy){
   
   this->dimL = toCopy->gDimL();
   this->dimR = toCopy->gDimR();
   this->phys_d = toCopy->gPhys_d();
   this->RN = toCopy->gRN();
   this->storageSize = toCopy->gStorageSize();
   
   storage = new complex<double>[storageSize];
   int inc = 1;
   zcopy_(&storageSize,toCopy->gStorage(),&inc,storage,&inc);
   
}

//construct from file: read in REAL numbers!
MPStensor::MPStensor(const char *filename,Random *RN){

   ifstream in(filename);
   in >> dimL >> phys_d >> dimR >> storageSize;

   this->RN = RN;

   storage = new complex<double> [storageSize];

   for(int i = 0;i < storageSize;++i){

      double value;

      in >> i >> value;
      
      storage[i] = complex<double>(value,0.0);

   }

}

void MPStensor::Reset(const int dimL, const int dimR){

   this->dimL = dimL;
   this->dimR = dimR;

   if(dimL*dimR*phys_d > storageSize){

      storageSize = dimL*dimR*phys_d;
      delete [] storage;
      storage = new complex<double>[storageSize];
   }

   random();

}

MPStensor::~MPStensor(){
   
   delete [] storage;
   
}

void MPStensor::random(){

   for(int cnt = 0;cnt < storageSize;cnt++)
      storage[cnt] = complex<double>(RN->rand() - 0.5,RN->rand() - 0.5);

}

int MPStensor::gDimL() const {

   return dimL; 

}

int MPStensor::gDimR() const {
   
   return dimR; 
   
}

int MPStensor::gPhys_d() const {
   
   return phys_d; 
   
}

int MPStensor::gStorageSize() const {
   
   return storageSize; 
   
}

complex<double> * MPStensor::gStorage() {
   
   return storage; 
   
}


complex<double> * MPStensor::gStorage(const int d_val){

   if((d_val<0) || (d_val>=phys_d))
      return NULL; 

   return storage + d_val*dimL*dimR;

}

Random * MPStensor::gRN(){

   return RN; 

}

void MPStensor::QR(complex<double> *Rmx, complex<double> *mem, complex<double> * tau, complex<double> *work){
      
   int m = dimL*phys_d;
   int n = dimR;
   
   for(int d_val=0; d_val<phys_d; d_val++)
      for(int irow=0; irow<dimL; irow++)
         for(int icol=0; icol<dimR; icol++)
            mem[irow + dimL*(d_val + phys_d * icol)] = storage[irow + dimL * (icol + dimR * d_val)];
         
   int info;   
   zgeqrf_(&m,&n,mem,&m,tau,work,&n,&info);
      
   for(int irow=0; irow<dimR; irow++){

      for (int icol=irow; icol<dimR; icol++)
         Rmx[irow + dimR*icol] = mem[irow + m*icol];

     for (int icol=0; icol<irow; icol++)
         Rmx[irow + dimR*icol] = 0.0;
      
   }
      
   zungqr_(&m,&n,&n,mem,&m,tau,work,&n,&info);
   
   for(int d_val=0; d_val<phys_d; d_val++)
      for(int irow=0; irow<dimL; irow++)
         for(int icol=0; icol<dimR; icol++)
            storage[irow + dimL * (icol + dimR * d_val)] = mem[irow + dimL * (d_val + phys_d * icol)];

}

/**
 * access to the individual numbers in an easy way
 * @param s physical dim index
 * @param i virtual row index
 * @param j virtual column index
 */
complex<double> &MPStensor::operator()(int s,int i,int j) const {

   return storage[i + dimL * (j + dimR * s)];

}

void MPStensor::LQ(complex<double> * Lmx, complex<double> * tau, complex<double> * work){
      
   int m = dimL;
   int n = dimR*phys_d;

   int info;
      
   zgelqf_(&m,&n,storage,&m,tau,work,&m,&info);
      
   for (int irow=0; irow<dimL; irow++){

      for (int icol=0; icol<=irow; icol++)
         Lmx[irow + dimL*icol] = storage[irow + m*icol];
      
      for (int icol=irow+1; icol<dimL; icol++)
         Lmx[irow + dimL*icol] = 0.0;
      
   }
      
   zunglq_(&m,&n,&m,storage,&m,tau,work,&m,&info);

}

void MPStensor::LeftMultiply(complex<double> * Lmx, complex<double> * work){

   int m = dimL;
   int n = dimR;
   int k = dimL;

   char notrans = 'N';
   complex<double> alpha(1.0,0.0);
   complex<double> beta(0.0,0.0);
   int dim = m*n;
   int inc = 1;

   for (int d_val=0; d_val<phys_d; d_val++){

      zgemm_(&notrans, &notrans, &m, &n, &k, &alpha, Lmx, &m, gStorage(d_val), &k, &beta, work, &m);
      zcopy_(&dim, work, &inc, gStorage(d_val), &inc);

   }
   
}

void MPStensor::RightMultiply(complex<double> * Rmx, complex<double> * work){

   int m = dimL;
   int n = dimR;
   int k = dimR;

   char notrans = 'N';
   complex<double> alpha(1.0,0.0);
   complex<double> beta(0.0,0.0);
   int dim = m*n;
   int inc = 1;

   for (int d_val=0; d_val<phys_d; d_val++){

      zgemm_(&notrans, &notrans, &m, &n, &k, &alpha, gStorage(d_val), &m, Rmx, &k, &beta, work, &m);
      zcopy_(&dim, work, &inc, gStorage(d_val), &inc);

   }
   
}

ostream &operator<<(ostream &output,const MPStensor &tensor){

   output << tensor.dimL << "\t" << tensor.phys_d << "\t" << tensor.dimR << "\t" << tensor.storageSize << endl;

   for(int i = 0;i < tensor.storageSize;++i)
      output << i << "\t" << tensor.storage[i] << endl;

   return output;

}
