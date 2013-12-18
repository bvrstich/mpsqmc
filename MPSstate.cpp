#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include "Lapack.h"
#include "MPSstate.h"
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
   for (int cnt=1; cnt<length; cnt++){
      VirtualD[cnt] = min(VirtualD[cnt-1] * phys_d, Dtrunc);
   }
   VirtualD[length] = 1;
   for (int cnt=length-1; cnt>0; cnt--){
      VirtualD[cnt] = min(min(VirtualD[cnt+1] * phys_d, Dtrunc), VirtualD[cnt]);
   }
   
   theTensors = new MPStensor * [length];
   for (int cnt=0; cnt<length; cnt++){
      theTensors[cnt] = new MPStensor(VirtualD[cnt], VirtualD[cnt+1], phys_d, RN);
   }
   
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
   for (int cnt=0; cnt<=length; cnt++){
      VirtualD[cnt] = VirtualDims[cnt];
   }
   
   theTensors = new MPStensor * [length];
   for (int cnt=0; cnt<length; cnt++){
      theTensors[cnt] = new MPStensor(VirtualD[cnt], VirtualD[cnt+1], phys_d, RN);
   }
   
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
   for (int cnt=0; cnt<=length; cnt++){
      VirtualD[cnt] = toCopy->gDimAtBound(cnt);
   }
   
   theTensors = new MPStensor * [length];
   for (int cnt=0; cnt<length; cnt++){
      theTensors[cnt] = new MPStensor(toCopy->gMPStensor(cnt));
   }
   
   TwoSiteObjectAllocated = false;
   work1Allocated = false;
   work2Allocated = false;
   work3Allocated = false;

}

//construct from file
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

      for(int j = 0;j < storsize;++j)
         in >> i >> j >> theTensors[i]->gStorage()[j];

   }


   //allocate the tensors
   TwoSiteObjectAllocated = false;
   work1Allocated = false;
   work2Allocated = false;
   work3Allocated = false;

}

MPSstate::~MPSstate(){
   
   for (int cnt=0; cnt<length; cnt++){
      delete theTensors[cnt];
   }
   delete [] theTensors;
   
   delete [] VirtualD;
   
   if (TwoSiteObjectAllocated){ delete the2siteObject; }
   if (work1Allocated){ delete [] work1; }
   if (work2Allocated){ delete [] work2; }
   if (work3Allocated){ delete [] work3; }
   
}

int MPSstate::gLength() const{ return length; }

int MPSstate::gPhys_d() const{ return phys_d; }

int MPSstate::gDtrunc() const{ return Dtrunc; }

int MPSstate::gDimAtBound(const int bound) const{

   if ((bound<0) || (bound>length)){ return 0; }
   return VirtualD[bound];

}

MPStensor * MPSstate::gMPStensor(const int site){

   if ((site<0) || (site>=length)){ return NULL; }
   return theTensors[site];

}

Random * MPSstate::gRN(){ return RN; }

void MPSstate::checkWork1(const int size){

   if (!work1Allocated){
      sizeWork1 = size;
      work1Allocated = true;
      work1 = new double[sizeWork1];
   } else {
      if (size > sizeWork1){
         sizeWork1 = size;
         delete [] work1;
         work1 = new double[sizeWork1];
      }
   }

}

void MPSstate::checkWork2(const int size){

   if (!work2Allocated){
      sizeWork2 = size;
      work2Allocated = true;
      work2 = new double[sizeWork2];
   } else {
      if (size > sizeWork2){
         sizeWork2 = size;
         delete [] work2;
         work2 = new double[sizeWork2];
      }
   }

}

void MPSstate::checkWork3(const int size){

   if (!work3Allocated){
      sizeWork3 = size;
      work3Allocated = true;
      work3 = new double[sizeWork3];
   } else {
      if (size > sizeWork3){
         sizeWork3 = size;
         delete [] work3;
         work3 = new double[sizeWork3];
      }
   }

}

double MPSstate::LeftNormalize(){

   checkWork1(Dtrunc*Dtrunc);
   checkWork2(Dtrunc*Dtrunc*phys_d);
   checkWork3(2*Dtrunc);

   for (int cnt=0; cnt<length; cnt++){
      theTensors[cnt]->QR(work1,work2,work3,work3+Dtrunc);
      if (cnt!=length-1){ theTensors[cnt+1]->LeftMultiply(work1,work2); }
   }
   
   return work1[0];

}

double MPSstate::RightNormalize(){

   checkWork1(Dtrunc*Dtrunc);
   checkWork2(Dtrunc*Dtrunc);
   checkWork3(2*Dtrunc);

   for (int cnt=length-1; cnt>=0; cnt--){
   
      theTensors[cnt]->LQ(work1,work3,work3+Dtrunc);
      if (cnt!=0){ theTensors[cnt-1]->RightMultiply(work1,work2); }
   
   }
   
   return work1[0];

}

double MPSstate::InnerProduct(MPSstate * OtherState){

   if ((phys_d != OtherState->gPhys_d()) || (length != OtherState->gLength())){ return 0.0; }
   
   int theworksize = Dtrunc * OtherState->gDtrunc();
   checkWork1(theworksize);
   checkWork2(theworksize);
   checkWork3(theworksize);
   
   work1[0] = 1;
   char notrans = 'N';
   char trans = 'T';
   double one = 1.0;
   double set = 0.0;
   for (int cnt=0; cnt<length; cnt++){
   
      int dimLthis = VirtualD[cnt];
      int dimRthis = VirtualD[cnt+1];
      int dimLother = OtherState->gDimAtBound(cnt);
      int dimRother = OtherState->gDimAtBound(cnt+1);
      
      for (int cnt2=0; cnt2<dimRthis*dimRother; cnt2++){ work2[cnt2] = 0.0; }
      
      for (int d_val=0; d_val<phys_d; d_val++){
         dgemm_(&notrans, &notrans, &dimLother, &dimRthis, &dimLthis, &one, work1, &dimLother, theTensors[cnt]->gStorage(d_val), &dimLthis, &set, work3, &dimLother);
         dgemm_(&trans, &notrans, &dimRother, &dimRthis, &dimLother, &one, OtherState->gMPStensor(cnt)->gStorage(d_val), &dimLother, work3, &dimLother, &one, work2, &dimRother);
      }
      
      double * swap = work1;
      work1 = work2;
      work2 = swap;
      
   }
   
   return work1[0];

}

void MPSstate::ScalarMultiplication(const double factor){

   double fact = pow(fabs(factor),1.0/length);
   
   int sign = 0;
   if (factor>0.0){ sign = 1; }
   if (factor<0.0){ sign = -1; }
   
   for (int site=0; site<length; site++){
      int dim = VirtualD[site] * VirtualD[site+1] * phys_d;
      double alpha = ((site==0)?sign:1) * fact;
      int inc = 1;
      dscal_(&dim,&alpha,theTensors[site]->gStorage(),&inc);
   }

}

/**
 * change the phase of the wavefunction, multiply tensor on site 0 with (-1)
 * needed for important sampling with constrained path.
 */
void MPSstate::ChangePhase(){

   const int site = 0;
   int dim = VirtualD[site] * VirtualD[site+1] * phys_d;
   double alpha = -1.0;
   int inc = 1;
   dscal_(&dim,&alpha,theTensors[site]->gStorage(),&inc);

}

void MPSstate::ResetContentsAndStoreSumOf(MPSstate * state1, MPSstate * state2){

   if ((state1->gPhys_d() != state2->gPhys_d()) || (state1->gLength()!=state2->gLength())){ return; }

   //Set the parameters and (if necessary) reallocate the theTensors
   length = state1->gLength();
   Dtrunc = state1->gDtrunc() + state2->gDtrunc();
   phys_d = state1->gPhys_d();
   for (int cnt=1; cnt<length; cnt++){ VirtualD[cnt] = state1->gDimAtBound(cnt) + state2->gDimAtBound(cnt); }
   for (int cnt=0; cnt<length; cnt++){ theTensors[cnt]->Reset(VirtualD[cnt], VirtualD[cnt+1]); }
   
   //Fill the MPS tensors
   int index = 0; //Left end
   for (int d_val=0; d_val<phys_d; d_val++){
   
      int dim1 = state1->gDimAtBound(index+1);
      int dim2 = state2->gDimAtBound(index+1);
      double * temp = theTensors[index]->gStorage(d_val); // 1 x (dim1 + dim2)
      double * temp1 = state1->gMPStensor(index)->gStorage(d_val);
      double * temp2 = state2->gMPStensor(index)->gStorage(d_val);
      for (int cnt=0; cnt<dim1; cnt++){ temp[cnt]      = temp1[cnt]; }
      for (int cnt=0; cnt<dim2; cnt++){ temp[dim1+cnt] = temp2[cnt]; }
      
   }
   
   index = length-1; //Right end
   for (int d_val=0; d_val<phys_d; d_val++){
   
      int dim1 = state1->gDimAtBound(index);
      int dim2 = state2->gDimAtBound(index);
      double * temp = theTensors[index]->gStorage(d_val); // (dim1 + dim2) x 1
      double * temp1 = state1->gMPStensor(index)->gStorage(d_val);
      double * temp2 = state2->gMPStensor(index)->gStorage(d_val);
      for (int cnt=0; cnt<dim1; cnt++){ temp[cnt]      = temp1[cnt]; }
      for (int cnt=0; cnt<dim2; cnt++){ temp[dim1+cnt] = temp2[cnt]; }
      
   }
   
   for (int site=1; site<length-1; site++){ //Other sites
      for (int d_val=0; d_val<phys_d; d_val++){
         int dim1_left  = state1->gDimAtBound(site);
         int dim2_left  = state2->gDimAtBound(site);
         int dim1_right = state1->gDimAtBound(site+1);
         int dim2_right = state2->gDimAtBound(site+1);
         
         double * temp = theTensors[site]->gStorage(d_val); // (dim1_left + dim2_left) x (dim1_right + dim2_right)
         double * temp1 = state1->gMPStensor(site)->gStorage(d_val);
         double * temp2 = state2->gMPStensor(site)->gStorage(d_val);
      
         for (int cnt=0; cnt<dim1_left; cnt++){
            for (int cnt2=0;          cnt2<dim1_right;            cnt2++){ temp[cnt + (dim1_left+dim2_left)*cnt2] = temp1[cnt+dim1_left*cnt2]; }
            for (int cnt2=dim1_right; cnt2<dim1_right+dim2_right; cnt2++){ temp[cnt + (dim1_left+dim2_left)*cnt2] = 0.0;                       }
         }
         for (int cnt=0; cnt<dim2_left; cnt++){
            for (int cnt2=0; cnt2<dim1_right; cnt2++){ temp[(cnt + dim1_left) + (dim1_left+dim2_left)*cnt2             ] = 0.0;                         }
            for (int cnt2=0; cnt2<dim2_right; cnt2++){ temp[(cnt + dim1_left) + (dim1_left+dim2_left)*(dim1_right+cnt2)] = temp2[cnt + dim2_left*cnt2]; }
         }
      }
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

void MPSstate::ApplyMPO(MPO * theMPO, MPSstate * Psi0){
   
   //Readjust the dimensions
   for (int cnt=1; cnt<length; cnt++){ VirtualD[cnt] = theMPO->dimL(cnt) * Psi0->gDimAtBound(cnt); }
   for (int cnt=0; cnt<length; cnt++){
      theTensors[cnt]->Reset(VirtualD[cnt], VirtualD[cnt+1]); //Only reallocates storage if it's required.
      const int dim = VirtualD[cnt] * VirtualD[cnt+1] * phys_d;
      double * temp = theTensors[cnt]->gStorage();
      for (int cnt2=0; cnt2<dim; cnt2++){ temp[cnt2] = 0.0; } //Storage is put to zero. If "AmIOp0()==true", then nothing should be done.
   }
   Dtrunc = 1;
   for (int cnt=0; cnt<=length; cnt++){
      if (VirtualD[cnt]>Dtrunc){ Dtrunc = VirtualD[cnt]; }
   }
   
   //Apply the MPO: OpenMP parallellizable
   #pragma omp parallel for schedule(static) default(none) shared(theMPO, Psi0)
   for (int site=0; site<length; site++){
      const int dimLsmall = Psi0->gDimAtBound(site);
      const int dimRsmall = Psi0->gDimAtBound(site+1);
      const int dimLmpo = theMPO->dimL(site);
      const int dimRmpo = theMPO->dimR(site);
      
      for (int MPOleft=0; MPOleft<dimLmpo; MPOleft++){
         for (int MPOright=0; MPOright<dimRmpo; MPOright++){
         
            Operator * theOp = theMPO->gOperator(site,MPOleft,MPOright);
            if (!(theOp->AmIOp0())){
               const double factor = theMPO->gPrefactor(site,MPOleft,MPOright);
               for (int phys_up=0; phys_up<phys_d; phys_up++){
                  double * target = theTensors[site]->gStorage(phys_up);
                  if (theOp->AmIOpI()){
                     double * source = Psi0->gMPStensor(site)->gStorage(phys_up);
                     for (int row=0; row<dimLsmall; row++){
                        for (int col=0; col<dimRsmall; col++){
                           target[ (MPOleft * dimLsmall + row) + (dimLsmall * dimLmpo) * (MPOright * dimRsmall + col) ] = factor * source[ row + dimLsmall * col ];
                        }
                     }
                  } else {
                     for (int phys_down=0; phys_down<phys_d; phys_down++){
                        const double factorbis = factor * (*theOp)(phys_up,phys_down);
                        if (factorbis != 0.0){
                           double * source = Psi0->gMPStensor(site)->gStorage(phys_down);
                           for (int row=0; row<dimLsmall; row++){
                              for (int col=0; col<dimRsmall; col++){
                                 target[ (MPOleft * dimLsmall + row) + (dimLsmall * dimLmpo) * (MPOright * dimRsmall + col) ] += factorbis * source[ row + dimLsmall * col ];
                              }
                           }
                        }
                     }
                  }
               }
            }
            
         }
      }
   }

}

void MPSstate::ApplyMPOterm(MPO * theMPO, const int SelectedTerm){
   
   for (int site=0; site<length; site++){
      Operator * theOp = theMPO->gRN_operator(SelectedTerm, site);
      if (!(theOp->AmIOpI())){
         int sizeBlock = VirtualD[site] * VirtualD[site+1];
         int size = sizeBlock * phys_d;
         checkWork1(size);
         for (int cnt=0; cnt<size; cnt++){ work1[cnt] = 0.0; }
         for (int phys_up=0; phys_up<phys_d; phys_up++){
            for (int phys_down=0; phys_down<phys_d; phys_down++){
               double OperatorValue = (*theOp)(phys_up,phys_down);
               if (OperatorValue!=0.0){
                  int inc=1;
                  daxpy_(&sizeBlock, &OperatorValue, theTensors[site]->gStorage(phys_down), &inc, work1 + phys_up*sizeBlock, &inc);
               }
            }
         }
         int inc = 1;
         dcopy_(&size, work1, &inc, theTensors[site]->gStorage(), &inc);
      }
   }
   
   ScalarMultiplication( theMPO->gRN_prefactor(SelectedTerm) );

}

/**
 * Apply a specific two-site trotter term, generated by doing the svd of e^{S_iS_j J_[ij]} to *this.
 * basically, we have in the TrotterHeisenbergclass: e^{S_i.S_j J_[ij]} = \sum_k (OL)^i_k (OR)^j_k: We apply in this function: (OL)^i_leftSVDindex to site i and 
 * (OR)^j_rightSVDindex to site j
 */
void MPSstate::ApplyTwoSiteTrotterTerm(TrotterHeisenberg * theTrotter, const int firstSite, const int secondSite, const int leftSVDindex, const int rightSVDindex, const bool doHC){

   //get J_[ij]
   const double coupling = theTrotter->gCoupling(firstSite, secondSite);
   
   {

      //First Site
      int sizeBlock = VirtualD[firstSite] * VirtualD[firstSite+1];
      int size = sizeBlock * phys_d;

      checkWork1(size);

      double SingularValueSqrt = sqrt( theTrotter->gTwoSitePropSVD_Sing(coupling, leftSVDindex) );

      for (int cnt=0; cnt<size; cnt++)
         work1[cnt] = 0.0;

      for(int phys_up = 0;phys_up < phys_d;phys_up++)//loop over the upper physical index of the operator
         for(int phys_down = 0;phys_down < phys_d;phys_down++){//loop over the lower physical index of the operator

            double OperatorValue = SingularValueSqrt * theTrotter->gTwoSitePropSVD_Left(coupling, leftSVDindex, (doHC)?phys_down:phys_up, (doHC)?phys_up:phys_down);

            if(OperatorValue != 0.0){

               int inc = 1;
               daxpy_(&sizeBlock, &OperatorValue, theTensors[firstSite]->gStorage(phys_down), &inc, work1 + phys_up*sizeBlock, &inc);

            }

         }

      int inc = 1;
      dcopy_(&size, work1, &inc, theTensors[firstSite]->gStorage(), &inc);

   }

   {
      //Second Site
      int sizeBlock = VirtualD[secondSite] * VirtualD[secondSite+1];
      int size = sizeBlock * phys_d;

      checkWork1(size);

      double SingularValueSqrt = sqrt( theTrotter->gTwoSitePropSVD_Sing(coupling, rightSVDindex) );

      for(int cnt = 0;cnt < size;cnt++)
         work1[cnt] = 0.0;

      for(int phys_up = 0;phys_up < phys_d;phys_up++)
         for(int phys_down = 0;phys_down < phys_d;phys_down++){

            double OperatorValue = SingularValueSqrt * theTrotter->gTwoSitePropSVD_Right(coupling, rightSVDindex, (doHC)?phys_down:phys_up, (doHC)?phys_up:phys_down);

            if(OperatorValue!=0.0){

               int inc = 1;
               daxpy_(&sizeBlock, &OperatorValue, theTensors[secondSite]->gStorage(phys_down), &inc, work1 + phys_up*sizeBlock, &inc);

            }

         }

      int inc = 1;
      dcopy_(&size, work1, &inc, theTensors[secondSite]->gStorage(), &inc);

   }

}

void MPSstate::ApplyTwoSiteTrotterTerm(TrotterHeisenberg * theTrotter, const int firstSite, const int secondSite, GridGenerator * theGrid, const int gridPoint){

   const double coupling = theTrotter->gCoupling(firstSite, secondSite);

   {
      //First Site
      int sizeBlock = VirtualD[firstSite] * VirtualD[firstSite+1];
      int size = sizeBlock * phys_d;

      checkWork1(size);

      for (int cnt=0;cnt < size;cnt++)
         work1[cnt] = 0.0; 

      for (int SVDindex=0; SVDindex<theGrid->gDim(); SVDindex++){

         const double prefactor = sqrt( theTrotter->gTwoSitePropSVD_Sing(coupling, SVDindex) * theGrid->gDim() ) * theGrid->gCoOfPoint(gridPoint, SVDindex);

         if (prefactor!=0.0){
            for (int phys_up=0; phys_up<phys_d; phys_up++){
               for (int phys_down=0; phys_down<phys_d; phys_down++){
                  double OperatorValue = prefactor * theTrotter->gTwoSitePropSVD_Left(coupling, SVDindex, phys_up, phys_down);
                  if (OperatorValue!=0.0){
                     int inc = 1;
                     daxpy_(&sizeBlock, &OperatorValue, theTensors[firstSite]->gStorage(phys_down), &inc, work1 + phys_up*sizeBlock, &inc);
                  }
               }
            }
         }
      }
      int inc = 1;
      dcopy_(&size, work1, &inc, theTensors[firstSite]->gStorage(), &inc);
   }

   {
      //Second Site
      int sizeBlock = VirtualD[secondSite] * VirtualD[secondSite+1];
      int size = sizeBlock * phys_d;
      checkWork1(size);
      for (int cnt=0; cnt<size; cnt++){ work1[cnt] = 0.0; }
      for (int SVDindex=0; SVDindex<theGrid->gDim(); SVDindex++){
         const double prefactor = sqrt( theTrotter->gTwoSitePropSVD_Sing(coupling, SVDindex) * theGrid->gDim() ) * theGrid->gCoOfPoint(gridPoint, SVDindex);
         if (prefactor!=0.0){
            for (int phys_up=0; phys_up<phys_d; phys_up++){
               for (int phys_down=0; phys_down<phys_d; phys_down++){
                  double OperatorValue = prefactor * theTrotter->gTwoSitePropSVD_Right(coupling, SVDindex, phys_up, phys_down);
                  if (OperatorValue!=0.0){
                     int inc = 1;
                     daxpy_(&sizeBlock, &OperatorValue, theTensors[secondSite]->gStorage(phys_down), &inc, work1 + phys_up*sizeBlock, &inc);
                  }
               }
            }
         }
      }
      int inc = 1;
      dcopy_(&size, work1, &inc, theTensors[secondSite]->gStorage(), &inc);
   }

}

void MPSstate::ApplyOneSiteTrotterTermEverywhere(TrotterHeisenberg * theTrotter){

   if (theTrotter->gIsMagneticField()){
      for (int site=0; site<length; site++){

         int sizeBlock = VirtualD[site] * VirtualD[site+1];
         int size = sizeBlock * phys_d;
         checkWork1(size);
         for (int cnt=0; cnt<size; cnt++){ work1[cnt] = 0.0; }
         for (int phys_up=0; phys_up<phys_d; phys_up++){
            for (int phys_down=0; phys_down<phys_d; phys_down++){
               double OperatorValue = theTrotter->gSingleSiteProp( phys_up, phys_down);
               if (OperatorValue!=0.0){
                  int inc = 1;
                  daxpy_(&sizeBlock, &OperatorValue, theTensors[site]->gStorage(phys_down), &inc, work1 + phys_up*sizeBlock, &inc);
               }
            }
         }
         int inc = 1;
         dcopy_(&size, work1, &inc, theTensors[site]->gStorage(), &inc);

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

      double *storage = mps.gMPStensor(i)->gStorage();

      for(int j = 0;j < storsize;++j)
         output << i << "\t" << j << "\t" << storage[j] << endl;

   }

   return output;

}
