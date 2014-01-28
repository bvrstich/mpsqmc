#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <omp.h>

#include "MPSstate.h"
#include "Lapack.h"
#include "Operator.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 9, 2013 */

using std::min;

WorkSpace *MPSstate::ws;

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

}

//construct from file: read in REAL numbers!
MPSstate::MPSstate(const char *filename,Random *RN){

   ifstream in(filename);
   
   in >> length >> Dtrunc >> phys_d;

   std::cout << length << "\t" << Dtrunc << "\t" << phys_d << endl;

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

}

MPSstate::~MPSstate(){
   
   for (int cnt=0; cnt<length; cnt++)
      delete theTensors[cnt];
   
   delete [] theTensors;

   delete [] VirtualD;
   
   if (TwoSiteObjectAllocated)
      delete the2siteObject; 

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

   return VirtualD[bound];

}

MPStensor * MPSstate::gMPStensor(const int site){

   return theTensors[site];

}

MPStensor &MPSstate::operator[](int site){

   return *theTensors[site];
}


Random * MPSstate::gRN(){
   
   return RN; 
   
}

complex<double> MPSstate::normalize(){

   complex<double> scal = std::sqrt(this->InnerProduct(this));

   this->ScalarMultiplication(1.0/scal);

   return scal;

}

complex<double> MPSstate::LeftNormalize(){

#ifdef _OPENMP
   int myID = omp_get_thread_num();
#else
   int myID = 0;
#endif
   
   ws->clean(myID);

   for(int cnt=0; cnt< length; cnt++){

      theTensors[cnt]->QR(ws->work1[myID],ws->work2[myID],ws->work3[myID],ws->work3[myID]+Dtrunc);

      if (cnt!=length-1)
         theTensors[cnt+1]->LeftMultiply(ws->work1[myID],ws->work2[myID]);

   }

   return ws->work1[myID][0];

}

complex<double> MPSstate::RightNormalize(){

#ifdef _OPENMP
   int myID = omp_get_thread_num();
#else
   int myID = 0;
#endif
   
   ws->clean(myID);

   for(int cnt=length-1; cnt>=0; cnt--){
   
      theTensors[cnt]->LQ(ws->work1[myID],ws->work3[myID],ws->work3[myID]+Dtrunc);

      if(cnt!=0)
         theTensors[cnt-1]->RightMultiply(ws->work1[myID],ws->work2[myID]);
   
   }
   
   return ws->work1[myID][0];

}

complex<double> MPSstate::expectation(MPO *theMPO,MPSstate * OtherState){

   int DO = theMPO->gDtrunc();
   int Dtrunc_other = OtherState->gDtrunc();
  
#ifdef _OPENMP
   int myID = omp_get_thread_num();
#else
   int myID = 0;
#endif
   
   //reset the workspaces to zero
   ws->clean(myID);

   char notrans = 'N';
   char trans = 'T';
   char herm = 'C';

   int inc = 1;

   complex<double> one(1.0,0.0);
   complex<double> zero(0.0,0.0);

   int dimRthis = this->gDimAtBound(1);
   int dimRother = OtherState->gDimAtBound(1);
   int dimRmpo = theMPO->dimR(0);

   int sizeBlock = dimRthis;
   int m,n,k;

   for(int o = 0;o < dimRmpo;++o)
      for(int s = 0;s < phys_d;++s)
         for(int s_ = 0;s_ < phys_d;++s_){

            complex<double> factor =  (*theMPO)(0,0,s_,s,o);

            if(std::abs(factor) > 1.0e-15)
               zaxpy_(&sizeBlock,&factor, theTensors[0]->gStorage(s_), &inc, ws->work2[myID] + s*dimRmpo*dimRthis + o*dimRthis, &inc);

         }

   m = dimRthis * dimRmpo;
   n = dimRother;
   k = phys_d;

   zgemm_(&notrans, &herm, &m, &n, &k, &one, ws->work2[myID] , &m, OtherState->gMPStensor(0)->gStorage() , &n, &zero, ws->work1[myID], &m);

   for(int site = 1;site < theMPO->gLength();++site){

      int dimRthis = this->gDimAtBound(site + 1);
      int dimLthis = this->gDimAtBound(site);

      int dimRother = OtherState->gDimAtBound(site + 1);
      int dimLother = OtherState->gDimAtBound(site);

      int dimLmpo = theMPO->dimL(site);
      int dimRmpo = theMPO->dimR(site);

      m = dimRthis * phys_d;
      n = dimLother * dimLmpo;
      k = dimLthis;

      zgemm_(&trans,&notrans, &m, &n, &k, &one, theTensors[site]->gStorage() , &k, ws->work1[myID] , &k, &zero, ws->work3[myID], &m);

      for(int i = 0;i < dimRthis;++i)
         for(int j = 0;j < dimLother;++j)
            for(int o = 0;o < dimRmpo;++o)
               for(int s = 0;s < phys_d;++s)
                  ws->work2[myID][s*dimLother*dimRthis*dimRmpo + j*dimRthis*dimRmpo + o*dimRthis + i] = 0.0;

      sizeBlock = dimRthis;

      for(int o = 0;o < dimRmpo;++o)
         for(int s = 0;s < phys_d;++s)
            for(int p = 0;p < dimLmpo;++p)
               for(int s_ = 0;s_ < phys_d;++s_){

                  complex<double> factor =  (*theMPO)(site,p,s_,s,o);

                  if(std::abs(factor) > 1.0e-15)
                     for(int j = 0;j < dimLother;++j)
                        zaxpy_(&sizeBlock,&factor, ws->work3[myID] + j*dimLmpo*dimRthis*phys_d + p*phys_d*dimRthis + s_*dimRthis, &inc, ws->work2[myID] + s*dimLother*dimRthis*dimRmpo + j*dimRthis*dimRmpo + o*dimRthis, &inc);

               }

      //permute index
      for(int s = 0;s < phys_d;++s)
         for(int i = 0;i < dimLother;++i)
            for(int j = 0;j < dimRother;++j)
               ws->work3[myID][s*dimLother*dimRother + i*dimRother + j] = OtherState->gMPStensor(site)->gStorage()[s*dimLother*dimRother + j*dimLother + i];
       
      m = dimRthis * dimRmpo;
      n = dimRother;
      k = dimLother * phys_d;

      zgemm_(&notrans,&herm, &m, &n, &k, &one, ws->work2[myID] , &m, ws->work3[myID] , &n, &zero, ws->work1[myID], &m);

   }

   return ws->work1[myID][0];

}

complex<double> MPSstate::InnerProduct(MPSstate * OtherState){
  
#ifdef _OPENMP
   int myID = omp_get_thread_num();
#else
   int myID = 0;
#endif
   
   //reset the workspaces to zero
   ws->clean(myID);

   ws->work1[myID][0] = 1;

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
         ws->work2[myID][cnt2] = complex<double>(0.0,0.0);

      for (int d_val=0;d_val < phys_d;d_val++){

         zgemm_(&notrans, &notrans, &dimLother, &dimRthis, &dimLthis, &one, ws->work1[myID], &dimLother, theTensors[cnt]->gStorage(d_val), &dimLthis, &set, ws->work3[myID], &dimLother);
         zgemm_(&herm, &notrans, &dimRother, &dimRthis, &dimLother, &one, OtherState->gMPStensor(cnt)->gStorage(d_val), &dimLother, ws->work3[myID], &dimLother, &one, ws->work2[myID], &dimRother);

      }

      complex<double> * swap = ws->work1[myID];
      ws->work1[myID] = ws->work2[myID];
      ws->work2[myID] = swap;

   }

   return ws->work1[myID][0];

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

   Dtrunc = 1;

   //new D
   for(int i = 0;i < length;++i)
      if(VirtualD[i] > Dtrunc)
         Dtrunc = VirtualD[i];

}

void MPSstate::printVdim() const {

   for(int i = 0;i < length + 1;++i)
      cout << i << "\t" << VirtualD[i] << endl;

}

/**
 * Apply an MPO to an MPS and get another MPS. H\Psi> = \phi>
 * @param conj if true use the hermitian conjugate of the MPO
 */
void MPSstate::ApplyMPO(bool conj,MPO * theMPO, MPSstate * Psi0){

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

               complex<double> factor = theMPO->gPrefactor(site,MPOleft,MPOright);

               //take the hermitian conjugate of the operator if requested
               if(conj)
                  factor = std::conj(factor);

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

                        complex<double> factorbis;

                        if(conj)
                           factorbis = factor * std::conj((*theOp)(phys_down,phys_up));
                        else
                           factorbis = factor * (*theOp)(phys_up,phys_down);

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

/**
 * Apply the single site part of the Hamiltonian e^{-dtH1} to the MPS. actually, there is a factor 1/2 here because of: e^-dtH = e^-dtH1/2 e^-dtH2 e^-dtH1/2
 */
void MPSstate::ApplyH1(TrotterHeisenberg * theTrotter){
  
#ifdef _OPENMP
   int myID = omp_get_thread_num();
#else
   int myID = 0;
#endif
   
   for (int site=0; site<length; site++){

      int sizeBlock = VirtualD[site] * VirtualD[site+1];
      int size = sizeBlock * phys_d;

      for (int cnt=0; cnt<size; cnt++)
         ws->work1[myID][cnt] = 0.0;

      for (int phys_up=0; phys_up<phys_d; phys_up++)
         for (int phys_down=0; phys_down<phys_d; phys_down++){

            complex<double> OperatorValue = theTrotter->gH1Prop( phys_up, phys_down);

            if (std::abs(OperatorValue) > 1.0e-15){

               int inc = 1;
               zaxpy_(&sizeBlock, &OperatorValue, theTensors[site]->gStorage(phys_down), &inc, ws->work1[myID] + phys_up*sizeBlock, &inc);

            }

         }

      int inc = 1;
      zcopy_(&size, ws->work1[myID], &inc, theTensors[site]->gStorage(), &inc);

   }

}

/**
 * Apply a specific auxiliary field operator to the MPS
 * @param k the index of the eigenvector of J_[ij], part of the index of the auxiliary field
 * @param r the 'type' of the single-site operator: 0=e^Sx, 1=e^Sy, 2=e^Sz
 * @param x stochastic variable drawn from a normal distribution, the auxiliary field
 * @param theTrotter the object containing the info about the propagators
 */
void MPSstate::ApplyAF(int k,int r,complex<double> x,TrotterHeisenberg * theTrotter){

#ifdef _OPENMP
   int myID = omp_get_thread_num();
#else
   int myID = 0;
#endif

   //fill the allocated memory with the correct values
   theTrotter->fillAFProp(myID,k,r,x);

   for(int site = 0;site < length;site++){

      int sizeBlock = VirtualD[site] * VirtualD[site+1];
      int size = sizeBlock * phys_d;

      for (int cnt=0; cnt<size; cnt++)
         ws->work1[myID][cnt] = 0.0;

      for (int phys_up=0; phys_up<phys_d; phys_up++)
         for (int phys_down=0; phys_down<phys_d; phys_down++){

            complex<double> OperatorValue = theTrotter->gAFProp(myID,site,phys_up,phys_down);

            if(std::abs(OperatorValue) > 1.0e-15){

               int inc = 1;
               zaxpy_(&sizeBlock, &OperatorValue, theTensors[site]->gStorage(phys_down), &inc, ws->work1[myID] + phys_up*sizeBlock, &inc);

            }

         }

      int inc = 1;
      zcopy_(&size, ws->work1[myID], &inc, theTensors[site]->gStorage(), &inc);

   }
}

/**
 * Apply a specific auxiliary field operator to the MPS
 * @param k the index of the eigenvector of J_[ij], part of the index of the auxiliary field
 * @param r the 'type' of the single-site operator: 0=e^Sx, 1=e^Sy, 2=e^Sz
 * @param x stochastic variable drawn from a normal distribution, the auxiliary field
 * @param theTrotter the object containing the info about the propagators
 */
void MPSstate::ApplyAF(int k,complex<double> x,TrotterHeisenberg * theTrotter){

#ifdef _OPENMP
   int myID = omp_get_thread_num();
#else
   int myID = 0;
#endif

   //fill the allocated memory with the correct values
   theTrotter->fillAFProp(myID,k,x,RN);

   for(int site = 0;site < length;site++){

      int sizeBlock = VirtualD[site] * VirtualD[site+1];
      int size = sizeBlock * phys_d;

      for (int cnt=0; cnt<size; cnt++)
         ws->work1[myID][cnt] = 0.0;

      for (int phys_up=0; phys_up<phys_d; phys_up++)
         for (int phys_down=0; phys_down<phys_d; phys_down++){

            complex<double> OperatorValue = theTrotter->gAFProp(myID,site,phys_up,phys_down);

            if(std::abs(OperatorValue) > 1.0e-15){

               int inc = 1;
               zaxpy_(&sizeBlock, &OperatorValue, theTensors[site]->gStorage(phys_down), &inc, ws->work1[myID] + phys_up*sizeBlock, &inc);

            }

         }

      int inc = 1;
      zcopy_(&size, ws->work1[myID], &inc, theTensors[site]->gStorage(), &inc);

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

/**
 * change the phase of the wavefunction, multiply tensor on site 0 with (-1)
 * needed for important sampling with constrained phase
 */
void MPSstate::ChangePhase(){

   const int site = 0;
   int dim = VirtualD[site] * VirtualD[site+1] * phys_d;
   complex<double> alpha(-1.0,0.0);
   int inc = 1;
   zscal_(&dim,&alpha,theTensors[site]->gStorage(),&inc);

}

/**
 * initialize the static workspace variables
 */
void MPSstate::InitWork(int D,int DO,int d){

   ws = new WorkSpace(D,DO,d);

}

void MPSstate::ClearWork(){

   delete ws;

}
