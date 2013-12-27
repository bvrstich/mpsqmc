#include <stdlib.h>
#include <algorithm>

#include "TwoSiteObject.h"
#include "MPStensor.h"
#include "Lapack.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 12, 2013 */

using std::min;
using std::max;

#include <iostream>
using namespace std;

const double TwoSiteObject::SchmidtTrunc = 1e-13;

TwoSiteObject::TwoSiteObject(const int DimL, const int DimR, const int phys_d){

   this->DimL = DimL;
   this->DimR = DimR;
   this->phys_d = phys_d;
   this->storageSize = DimL*DimR*phys_d*phys_d;
   
   storage = new complex<double> [storageSize];
   work_large = new complex<double> [storageSize];
   work_large2 = new complex<double> [storageSize];
   work_large3 = new complex<double>[storageSize];
   
   int dimLeft  = DimL*phys_d;
   int dimRight = DimR*phys_d;
   SVD_dimMin   = min(dimLeft,dimRight);
   SVD_lwork    = 3*SVD_dimMin*SVD_dimMin + max(max(dimLeft,dimRight),4*SVD_dimMin*SVD_dimMin+4*SVD_dimMin);
   SVD_lrwork   = 5*SVD_dimMin*SVD_dimMin + 7*SVD_dimMin;
   
   SVD_Svalues = new double[SVD_dimMin];
   SVD_work    = new complex<double> [SVD_lwork];
   SVD_rwork    = new double [SVD_lrwork];
   SVD_iwork   = new int[SVD_dimMin*8];

}

TwoSiteObject::~TwoSiteObject(){
   
   delete [] storage;
   delete [] work_large;
   delete [] work_large2;
   delete [] work_large3;
   
   delete [] SVD_Svalues;
   delete [] SVD_work;
   delete [] SVD_iwork;
   delete [] SVD_rwork;
   
}

int TwoSiteObject::gDimL() const{
   
   return DimL; 
   
}

int TwoSiteObject::gDimR() const{
   
   return DimR; 
   
}

int TwoSiteObject::gPhys_d() const{
   
   return phys_d; 
   
}

int TwoSiteObject::gStorageSize() const{
   
   return storageSize; 
   
}

complex<double> * TwoSiteObject::gStorage(const int d_left, const int d_right){
   
   if ((d_left<0) || (d_left>=phys_d) || (d_right<0) || (d_right>=phys_d))
      return NULL; 

   return storage + DimL*DimR*(d_left + phys_d*d_right);

}

complex<double> * TwoSiteObject::gStorage(){
   
   return storage; 
   
}

void TwoSiteObject::Compose(MPStensor * MPSleft, MPStensor * MPSright){

   if ((MPSleft->gDimR() != MPSright->gDimL()) || (MPSleft->gPhys_d() != phys_d) || (MPSright->gPhys_d() != phys_d)) 
      return; 

   //Check if the current arrays are long enough
   int requiredSpace = MPSleft->gDimL() * MPSright->gDimR() * phys_d * phys_d;
   int requiredDimMin = phys_d * min( MPSleft->gDimL() , MPSright->gDimR() );
   int requiredLWork = 3*requiredDimMin*requiredDimMin + max(phys_d*max(MPSleft->gDimL(),MPSright->gDimR()),4*requiredDimMin*requiredDimMin+4*requiredDimMin);
   int requiredLRWork = 5*requiredDimMin*requiredDimMin + 7*requiredDimMin;

   if(requiredSpace > storageSize){

      storageSize = requiredSpace;

      delete [] storage;
      delete [] work_large;
      delete [] work_large2;
      delete [] work_large3;

      storage = new complex<double>[storageSize];
      work_large = new complex<double>[storageSize];
      work_large2 = new complex<double>[storageSize];
      work_large3 = new complex<double>[storageSize];
      
   }

   if(requiredDimMin > SVD_dimMin){

      SVD_dimMin = requiredDimMin;

      delete [] SVD_iwork;
      delete [] SVD_Svalues;

      SVD_iwork = new int[8*SVD_dimMin];
      SVD_Svalues = new double[SVD_dimMin];

   }

   if(requiredLWork > SVD_lwork){

      SVD_lwork = requiredLWork;

      delete [] SVD_work;

      SVD_work = new complex<double> [SVD_lwork];

   }

   if(requiredLRWork > SVD_lrwork){

      SVD_lrwork = requiredLRWork;

      delete [] SVD_rwork;

      SVD_rwork = new double [SVD_lrwork];

   }

   this->DimL = MPSleft->gDimL();
   this->DimR = MPSright->gDimR();

   //Do the actual work
   char notrans = 'N';
   int DimM = MPSleft->gDimR();
   complex<double> alpha(1.0,0.0);
   complex<double> beta(0.0,0.0);

   for(int d_left=0; d_left<phys_d; d_left++)
      for(int d_right=0; d_right<phys_d; d_right++)
         zgemm_(&notrans,&notrans,&DimL,&DimR,&DimM,&alpha,MPSleft->gStorage(d_left),&DimL,MPSright->gStorage(d_right),&DimM,&beta,gStorage(d_left,d_right),&DimL);

}

int TwoSiteObject::Decompose(MPStensor * MPSleft, MPStensor * MPSright, const int Dtrunc, bool movingright, bool possiblyCompress){

   complex<double> * mem = work_large;

   for (int alpha=0; alpha<DimL; alpha++)
      for (int i_left=0; i_left<phys_d; i_left++)
         for (int beta=0; beta<DimR; beta++)
            for (int i_right=0; i_right<phys_d; i_right++)
               mem[alpha + DimL*( i_left + phys_d*( beta + DimR*i_right ))] = gStorage(i_left,i_right)[alpha + DimL*beta];

   int dimLeft = DimL*phys_d;
   int dimRight = DimR*phys_d;
   int dimMin = min(dimLeft,dimRight);

   char jobz = 'S';

   double * Svalues = SVD_Svalues;
   complex<double> * Uvalues = work_large2;
   complex<double> * VTvalues = work_large3;

   int info;

   //dgesdd is not thread-safe in every implementation (intel MKL is safe, Atlas is not safe)
   //#pragma omp critical
   zgesdd_(&jobz, &dimLeft, &dimRight, mem, &dimLeft, Svalues, Uvalues, &dimLeft, VTvalues, &dimMin, SVD_work, &SVD_lwork,SVD_rwork, SVD_iwork, &info);

   int dimTrunc = Dtrunc;

   if(possiblyCompress){

      int index_trunc = dimMin;

      for (int cnt=0; cnt<dimMin; cnt++)
         if (Svalues[cnt]<SchmidtTrunc){

            index_trunc = cnt;
            cnt = dimMin;

         }

      if (index_trunc < dimTrunc)
         dimTrunc = index_trunc;

   }

   if (dimTrunc != MPSleft->gDimR()){

      MPSleft->Reset(DimL,dimTrunc);
      MPSright->Reset(dimTrunc,DimR);

   }

   for(int i_phys=0; i_phys<phys_d; i_phys++){

      complex<double> * temp = MPSleft->gStorage(i_phys);

      for (int alpha=0; alpha<DimL; alpha++)
         for (int beta=0; beta<dimTrunc; beta++)
            temp[alpha + DimL*beta] = Uvalues[alpha + DimL*(i_phys + phys_d*beta)] * ((movingright)?1.0:Svalues[beta]);

   }

   for(int i_phys = 0;i_phys<phys_d; i_phys++){

      complex<double> * temp = MPSright->gStorage(i_phys);

      for (int alpha=0; alpha<dimTrunc; alpha++)
         for (int beta=0; beta<DimR; beta++)
            temp[alpha + dimTrunc*beta] = VTvalues[alpha + dimMin*(beta + DimR*i_phys)] * ((movingright)?Svalues[alpha]:1.0);

   }

   return dimTrunc;
}

ostream &operator<<(ostream &output,const TwoSiteObject &tso){

   for(int i = 0;i < tso.storageSize;++i)
      output << i << "\t" << tso.storage[i] << endl;

   return output;

}
