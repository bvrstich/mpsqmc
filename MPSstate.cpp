#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include "MPSstate.h"
#include "Lapack.h"
#include "Operator.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 9, 2013 */

using std::min;

MPSstate::MPSstate(const int length, const int Dtrunc, const int phys_d, Random * RN){

   this->length = length;
   this->Dtrunc = Dtrunc;
   this->phys_d = phys_d;
   this->RN = RN;
   
   VirtualD = new int[length+1];
   VirtualD[0] = 1;

   for (int cnt=1; cnt<length; cnt++)
      VirtualD[cnt] = min(VirtualD[cnt-1] * phys_d, Dtrunc);
   
   VirtualD[length] = 1;

   for (int cnt=length-1; cnt>0; cnt--)
      VirtualD[cnt] = min(min(VirtualD[cnt+1] * phys_d, Dtrunc), VirtualD[cnt]);
   
   theTensors = new MPStensor * [length];

   for (int cnt=0; cnt<length; cnt++)
      theTensors[cnt] = new MPStensor(VirtualD[cnt], VirtualD[cnt+1], phys_d, RN);
   
   TwoSiteObjectAllocated = false;
   work1Allocated = false;
   work2Allocated = false;
   work3Allocated = false;

}

MPSstate::MPSstate(const int length, const int Dtrunc, const int phys_d, int * VirtualDims, Random * RN){

   this->length = length;
   this->Dtrunc = Dtrunc;
   this->phys_d = phys_d;
   this->RN = RN;
   
   VirtualD = new int[length+1];

   for (int cnt=0; cnt<=length; cnt++)
      VirtualD[cnt] = VirtualDims[cnt];
   
   
   theTensors = new MPStensor * [length];

   for (int cnt=0; cnt<length; cnt++)
      theTensors[cnt] = new MPStensor(VirtualD[cnt], VirtualD[cnt+1], phys_d, RN);
   
   TwoSiteObjectAllocated = false;
   work1Allocated = false;
   work2Allocated = false;
   work3Allocated = false;

}

MPSstate::MPSstate(MPSstate * toCopy){

   this->length = toCopy->gLength();
   this->phys_d = toCopy->gPhys_d();
   this->Dtrunc = toCopy->gDtrunc();
   this->RN = toCopy->gRN();
   
   VirtualD = new int[length+1];

   for (int cnt=0; cnt<=length; cnt++)
      VirtualD[cnt] = toCopy->gDimAtBound(cnt);
   

   theTensors = new MPStensor * [length];

   for(int cnt=0; cnt<length; cnt++)
      theTensors[cnt] = new MPStensor(toCopy->gMPStensor(cnt));
 
   TwoSiteObjectAllocated = false;
   work1Allocated = false;
   work2Allocated = false;
   work3Allocated = false;

}

//construct from file: read in REAL numbers!
MPSstate::MPSstate(const char *filename,Random *RN){

   ifstream in(filename);
   in >> length >> Dtrunc >> phys_d;

   this->RN = RN;

   //get the dimensions
   VirtualD = new int [length + 1];

   for(int i = 0;i < length + 1;++i)
      in >> i >> VirtualD[i];

   theTensors = new MPStensor * [length];

   for(int cnt = 0;cnt < length;cnt++)
      theTensors[cnt] = new MPStensor(VirtualD[cnt],VirtualD[cnt+1],phys_d,RN);

   //now fill the rest from file
   for(int i = 0;i < length;++i){

      int storsize;

      in >> i >> storsize;

      for(int j = 0;j < storsize;++j){

         double value;

         in >> i >> j >> value;

         theTensors[i]->gStorage()[j] = complex<double>(value,0.0);

      }

   }

   //allocate the tensors
   TwoSiteObjectAllocated = false;
   work1Allocated = false;
   work2Allocated = false;
   work3Allocated = false;

}

MPSstate::~MPSstate(){
   
   for (int cnt=0; cnt<length; cnt++)
      delete theTensors[cnt];
   
   delete [] theTensors;

   delete [] VirtualD;
   
   if (TwoSiteObjectAllocated)
      delete the2siteObject; 

   if (work1Allocated)
      delete [] work1; 

   if(work2Allocated)
      delete [] work2; 

   if(work3Allocated)
      delete [] work3; 
   
}

int MPSstate::gLength() const{ 
   
   return length; 
   
}

int MPSstate::gPhys_d() const{
   
   return phys_d; 
   
}

int MPSstate::gDtrunc() const{
   
   return Dtrunc; 
   
}

int MPSstate::gDimAtBound(const int bound) const{

   if ((bound<0) || (bound>length))
      return 0;

   return VirtualD[bound];

}

MPStensor * MPSstate::gMPStensor(const int site){

   if ((site<0) || (site>=length)) 
      return NULL;

   return theTensors[site];

}

Random * MPSstate::gRN(){
   
   return RN; 
   
}

void MPSstate::checkWork1(const int size){

   if (!work1Allocated){

      sizeWork1 = size;
      work1Allocated = true;
      work1 = new complex<double>[sizeWork1];

   } 
   else {

      if (size > sizeWork1){

         sizeWork1 = size;
         delete [] work1;
         work1 = new complex<double>[sizeWork1];

      }
   }

}

void MPSstate::checkWork2(const int size){

   if (!work2Allocated){

      sizeWork2 = size;
      work2Allocated = true;
      work2 = new complex<double>[sizeWork2];

   } 
   else {

      if (size > sizeWork2){

         sizeWork2 = size;
         delete [] work2;
         work2 = new complex<double>[sizeWork2];

      }
   }

}

void MPSstate::checkWork3(const int size){

   if (!work3Allocated){

      sizeWork3 = size;
      work3Allocated = true;
      work3 = new complex<double>[sizeWork3];

   } 
   else {

      if (size > sizeWork3){

         sizeWork3 = size;
         delete [] work3;
         work3 = new complex<double>[sizeWork3];

      }
   }

}

complex<double> MPSstate::LeftNormalize(){

   checkWork1(Dtrunc*Dtrunc);
   checkWork2(Dtrunc*Dtrunc*phys_d);
   checkWork3(2*Dtrunc);

   for (int cnt=0; cnt< length; cnt++){

      theTensors[cnt]->QR(work1,work2,work3,work3+Dtrunc);

      if (cnt!=length-1)
         theTensors[cnt+1]->LeftMultiply(work1,work2);

   }
   
   return work1[0];

}

complex<double> MPSstate::RightNormalize(){

   checkWork1(Dtrunc*Dtrunc);
   checkWork2(Dtrunc*Dtrunc);
   checkWork3(2*Dtrunc);

   for (int cnt=length-1; cnt>=0; cnt--){
   
      theTensors[cnt]->LQ(work1,work3,work3+Dtrunc);

      if (cnt!=0)
         theTensors[cnt-1]->RightMultiply(work1,work2);
   
   }
   
   return work1[0];

}

complex<double> MPSstate::InnerProduct(MPSstate * OtherState){

   if((phys_d != OtherState->gPhys_d()) || (length != OtherState->gLength()))
      return 0.0;
   
   int theworksize = Dtrunc * OtherState->gDtrunc();
   checkWork1(theworksize);
   checkWork2(theworksize);
   checkWork3(theworksize);
   
   work1[0] = 1;

   char notrans = 'N';
   char herm = 'C';

   complex<double> one(1.0,0.0);
   complex<double> set(0.0,0.0);

   for(int cnt = 0;cnt < length;cnt++){
   
      int dimLthis = VirtualD[cnt];
      int dimRthis = VirtualD[cnt+1];

      int dimLother = OtherState->gDimAtBound(cnt);
      int dimRother = OtherState->gDimAtBound(cnt+1);
      
      for (int cnt2=0; cnt2<dimRthis*dimRother; cnt2++)
         work2[cnt2] = complex<double>(0.0,0.0);
      
      for (int d_val=0;d_val < phys_d;d_val++){

         zgemm_(&notrans, &notrans, &dimLother, &dimRthis, &dimLthis, &one, work1, &dimLother, theTensors[cnt]->gStorage(d_val), &dimLthis, &set, work3, &dimLother);
         zgemm_(&herm, &notrans, &dimRother, &dimRthis, &dimLother, &one, OtherState->gMPStensor(cnt)->gStorage(d_val), &dimLother, work3, &dimLother, &one, work2, &dimRother);

      }
      
      complex<double> * swap = work1;
      work1 = work2;
      work2 = swap;
      
   }
   
   return work1[0];

}


void MPSstate::ScalarMultiplication(const complex<double> factor){

   complex<double> site_fact = pow(factor,complex<double>(1.0/length,0.0));
   
   for (int site=0; site<length; site++){

      int dim = VirtualD[site] * VirtualD[site+1] * phys_d;

      int inc = 1;
      zscal_(&dim,&site_fact,theTensors[site]->gStorage(),&inc);

   }

}

void MPSstate::CompressState(const int truncD){

   if (!TwoSiteObjectAllocated){
      the2siteObject = new TwoSiteObject(Dtrunc, Dtrunc, phys_d);
      TwoSiteObjectAllocated = true;
   }

   for (int cnt=0; cnt<length-1; cnt++){
      the2siteObject->Compose(theTensors[cnt],theTensors[cnt+1]);
      VirtualD[cnt+1] = the2siteObject->Decompose(theTensors[cnt], theTensors[cnt+1], truncD, true, true);
   }
   
   for (int cnt=length-2; cnt>=0; cnt--){
      the2siteObject->Compose(theTensors[cnt],theTensors[cnt+1]);
      VirtualD[cnt+1] = the2siteObject->Decompose(theTensors[cnt], theTensors[cnt+1], truncD, false, true);
   }

   Dtrunc = truncD;

}

void MPSstate::printVdim() const {

   for(int i = 0;i < length + 1;++i)
      cout << i << "\t" << VirtualD[i] << endl;

}

void MPSstate::ApplyMPO(MPO * theMPO, MPSstate * Psi0){
   
   //Readjust the dimensions
   for(int cnt=1; cnt<length; cnt++) 
      VirtualD[cnt] = theMPO->dimL(cnt) * Psi0->gDimAtBound(cnt);

   for(int cnt=0; cnt<length; cnt++){

      theTensors[cnt]->Reset(VirtualD[cnt], VirtualD[cnt+1]); //Only reallocates storage if it's required.
      const int dim = VirtualD[cnt] * VirtualD[cnt+1] * phys_d;
      complex<double> * temp = theTensors[cnt]->gStorage();

      for (int cnt2=0; cnt2<dim; cnt2++) 
         temp[cnt2] = 0.0;  //Storage is put to zero. If "AmIOp0()==true", then nothing should be done.
   }

   Dtrunc = 1;

   for (int cnt=0; cnt<=length; cnt++)
      if (VirtualD[cnt] > Dtrunc)
         Dtrunc = VirtualD[cnt]; 
  
   for (int site=0; site<length; site++){

      const int dimLsmall = Psi0->gDimAtBound(site);
      const int dimRsmall = Psi0->gDimAtBound(site+1);

      const int dimLmpo = theMPO->dimL(site);
      const int dimRmpo = theMPO->dimR(site);
      
      for (int MPOleft = 0; MPOleft < dimLmpo; MPOleft++)
         for (int MPOright = 0; MPOright < dimRmpo; MPOright++){

            Operator * theOp = theMPO->gOperator(site,MPOleft,MPOright);

            if (!(theOp->AmIOp0())){

               const double factor = theMPO->gPrefactor(site,MPOleft,MPOright);

               for(int phys_up = 0;phys_up < phys_d;phys_up++){//physical index of the resultant MPStensor: i.e. upper index of the MPO

                  complex<double> * target = theTensors[site]->gStorage(phys_up);

                  if (theOp->AmIOpI()){

                     complex<double> * source = Psi0->gMPStensor(site)->gStorage(phys_up);

                     //copy the old tensor into the new larger one:
                     for (int row=0; row<dimLsmall; row++)
                        for (int col=0; col<dimRsmall; col++)
                           target[ (MPOleft * dimLsmall + row) + (dimLsmall * dimLmpo) * (MPOright * dimRsmall + col) ] = factor * source[ row + dimLsmall * col ];
                        
                     
                  }
                  else{

                     for(int phys_down = 0;phys_down < phys_d;phys_down++){//loop over the lower physical index of the MPO: this will be contracted with the phys index of the original MPStensor

                        const complex<double> factorbis = factor * (*theOp)(phys_up,phys_down);

                        if(std::abs(factorbis) > 1.0e-15){

                           complex<double> * source = Psi0->gMPStensor(site)->gStorage(phys_down);

                           //do the contraction of phys_down
                           for (int row=0; row<dimLsmall; row++)
                              for (int col=0; col<dimRsmall; col++)
                                 target[ (MPOleft * dimLsmall + row) + (dimLsmall * dimLmpo) * (MPOright * dimRsmall + col) ] += factorbis * source[ row + dimLsmall * col ];


                        }
                     }

                  }

               }//loop over phys_up close

            }//if not zero

         }//mpo left right loop

   }//loop over sites

}

void MPSstate::ApplyOneSiteTrotterTermEverywhere(TrotterHeisenberg * theTrotter){

   if (theTrotter->gIsMagneticField()){

      for (int site=0; site<length; site++){

         int sizeBlock = VirtualD[site] * VirtualD[site+1];
         int size = sizeBlock * phys_d;
         checkWork1(size);

         for (int cnt=0; cnt<size; cnt++)
            work1[cnt] = 0.0;

         for (int phys_up=0; phys_up<phys_d; phys_up++)
            for (int phys_down=0; phys_down<phys_d; phys_down++){

               complex<double> OperatorValue = theTrotter->gSingleSiteProp( phys_up, phys_down);

               if (std::abs(OperatorValue) > 1.0e-15){

                  int inc = 1;
                  zaxpy_(&sizeBlock, &OperatorValue, theTensors[site]->gStorage(phys_down), &inc, work1 + phys_up*sizeBlock, &inc);

               }

            }

         int inc = 1;
         zcopy_(&size, work1, &inc, theTensors[site]->gStorage(), &inc);

      }
   }

}

ostream &operator<<(ostream &output,MPSstate &mps){

   //first print the essentials
   output << mps.length << "\t" << mps.Dtrunc << "\t" << mps.phys_d << endl;

   //then the bond dimensions
   for(int i = 0;i < mps.length+1;++i)
      output << i << "\t" << mps.VirtualD[i] << endl;

   //finally the tensors themselves
   for(int i = 0;i < mps.length;++i){

      int storsize = mps.VirtualD[i] * mps.VirtualD[i + 1] * mps.phys_d;

      output << i << "\t" << storsize << endl;

      complex<double> *storage = mps.gMPStensor(i)->gStorage();

      for(int j = 0;j < storsize;++j)
         output << i << "\t" << j << "\t" << storage[j] << endl;

   }

   return output;

}
