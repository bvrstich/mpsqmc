#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "DMRG.h"
#include "Lapack.h"
#include "Operator.h"

/*  Written by Sebastian Wouters <sebastianwouters@gmail.com> on August 16, 2013 */

using namespace std;

DMRG::DMRG(MPSstate * MPStrial, MPO * theMPO){

   this->MPSsolution = MPStrial;
   this->theMPO = theMPO;
   this->length = MPSsolution->gLength();
   
   MPOMPS = new MPSstate(length, MPSsolution->gDtrunc(), MPSsolution->gPhys_d(), MPSsolution->gRN());
   MPSsolution->LeftNormalize();
   MPOMPS->ApplyMPO(theMPO,MPSsolution);
   
   boundaryTerms = new double*[length-1];
   for (int cnt=0; cnt<length-1; cnt++){ boundaryTerms[cnt] = new double[MPSsolution->gDimAtBound(cnt+1) * MPOMPS->gDimAtBound(cnt+1)]; }
   the2siteObject = new TwoSiteObject(MPSsolution->gDtrunc(),MPSsolution->gDtrunc(),MPSsolution->gPhys_d());
   SolveHeffAllocated = false;

}

DMRG::~DMRG(){
   
   for (int cnt=0; cnt<length-1; cnt++){ delete [] boundaryTerms[cnt]; }
   delete [] boundaryTerms;
   
   delete MPOMPS;
   delete the2siteObject;
   
   if (SolveHeffAllocated){
      delete [] which;
      delete [] resid;
      delete [] v;
      delete [] iparam;
      delete [] ipntr;
      delete [] workd;
      delete [] workl;
   }
   
}

double DMRG::Solve(){

   //Construct boundaryTerms all OK
   for (int bound=1; bound<length-1; bound++){ ConstructBoundaryTermMovingRight(bound); }
   
   double EnergyPrevious = 1.0;
   double Energy = 0.0;
   
   while (fabs(EnergyPrevious-Energy)>1e-8){
      EnergyPrevious = Energy;
      Energy = sweepleft();
      Energy = sweepright();
      cout << "***** DMRG :: Energy at right edge = " << Energy << endl;
      cout << "***** DMRG :: Energy difference with previous sweep is = " << fabs (Energy-EnergyPrevious) << endl;
   }
   
   return Energy;

}

double DMRG::sweepleft(){

   double Energy = 0.0;
   for (int index=length-2; index>0; index--){
      the2siteObject->Compose(MPSsolution->gMPStensor(index),MPSsolution->gMPStensor(index+1));
      Energy = SolveHeff(index);
      cout << "DMRG :: Energy at site " << index << " = " << Energy << endl;
      the2siteObject->Decompose(MPSsolution->gMPStensor(index),MPSsolution->gMPStensor(index+1),MPSsolution->gDimAtBound(index+1),false,false);
      MPOMPS->ApplyMPO(theMPO,MPSsolution); //Can be done more efficiently, yes I know.
      ConstructBoundaryTermMovingLeft(index+1);
   }
   return Energy;

}

double DMRG::sweepright(){

   double Energy = 0.0;
   for (int index=0; index<length-2; index++){
      the2siteObject->Compose(MPSsolution->gMPStensor(index),MPSsolution->gMPStensor(index+1));
      Energy = SolveHeff(index);
      cout << "DMRG :: Energy at site " << index << " = " << Energy << endl;
      the2siteObject->Decompose(MPSsolution->gMPStensor(index),MPSsolution->gMPStensor(index+1),MPSsolution->gDimAtBound(index+1),true,false);
      MPOMPS->ApplyMPO(theMPO,MPSsolution); //Can be done more efficiently, yes I know.
      ConstructBoundaryTermMovingRight(index+1);
   }
   return Energy;

}

void DMRG::ConstructBoundaryTermMovingRight(const int bound){

   int dimUpR = MPSsolution->gDimAtBound(bound);
   int dimDownR = MPOMPS->gDimAtBound(bound);
   int dimUpL = MPSsolution->gDimAtBound(bound-1);
   int dimDownL = MPOMPS->gDimAtBound(bound-1);
   
   double one = 1.0;
   double zero = 0.0;
   double * memPrevious = (bound==1)? &one : boundaryTerms[bound-2];

   for (int cnt=0; cnt<dimUpR*dimDownR; cnt++){ boundaryTerms[bound-1][cnt] = 0.0; }
   
   char notrans = 'N';
   char trans = 'T';
   
   MPOMPS->checkWork1(dimUpL*dimDownR);
   double * workspace = MPOMPS->gWork1();
   
   for (int i_phys=0; i_phys<theMPO->gPhys_d(); i_phys++){
      dgemm_(&notrans,&notrans,&dimUpL,&dimDownR,&dimDownL,&one,memPrevious,&dimUpL,MPOMPS->gMPStensor(bound-1)->gStorage(i_phys),&dimDownL,&zero,workspace,&dimUpL);
      dgemm_(&trans,&notrans,&dimUpR,&dimDownR,&dimUpL,&one,MPSsolution->gMPStensor(bound-1)->gStorage(i_phys),&dimUpL,workspace,&dimUpL,&one,boundaryTerms[bound-1],&dimUpR);
   }

}

void DMRG::ConstructBoundaryTermMovingLeft(const int bound){

   int dimUpL = MPSsolution->gDimAtBound(bound);
   int dimDownL = MPOMPS->gDimAtBound(bound);
   int dimUpR = MPSsolution->gDimAtBound(bound+1);
   int dimDownR = MPOMPS->gDimAtBound(bound+1);
   
   double one = 1.0;
   double zero = 0.0;
   double * memPrevious = (bound==length-1)? &one : boundaryTerms[bound];

   for (int cnt=0; cnt<dimUpL*dimDownL; cnt++){ boundaryTerms[bound-1][cnt] = 0.0; }
   
   char notrans = 'N';
   char trans = 'T';
   
   MPOMPS->checkWork1(dimUpL*dimDownR);
   double * workspace = MPOMPS->gWork1();
   
   for (int i_phys=0; i_phys<theMPO->gPhys_d(); i_phys++){
      dgemm_(&notrans,&notrans,&dimUpL,&dimDownR,&dimUpR,&one,MPSsolution->gMPStensor(bound)->gStorage(i_phys),&dimUpL,memPrevious,&dimUpR,&zero,workspace,&dimUpL);
      dgemm_(&notrans,&trans,&dimUpL,&dimDownL,&dimDownR,&one,workspace,&dimUpL,MPOMPS->gMPStensor(bound)->gStorage(i_phys),&dimDownL,&one,boundaryTerms[bound-1],&dimUpL);
   }

}

double DMRG::SolveHeff(const int index){

   int ido = 0;  // internal communication flag
   char bmat = 'I';  // standard eigenvalue problem : A*x = lambda*x
   int n = MPSsolution->gDimAtBound(index) * MPSsolution->gPhys_d() * MPSsolution->gPhys_d() * MPSsolution->gDimAtBound(index+2);  // dimension of the eigenvalue problem
   
   int nev = 1; // only the smallest one
   double tol = 0.0; // tolerance for the calculation; 0.0 is machine precision
   int inc = 1;
   const int MAGICNUMBER = 32;
   int ncv = (n>MAGICNUMBER)?MAGICNUMBER:n; //number of Arnoldi vectors; capped to be tractable
   int ldv = n; // leading dimension of the v-mx
   int lworkl = (MAGICNUMBER + 8) * MAGICNUMBER; // size of work space nr. 2
   
   if (!SolveHeffAllocated){
      SolveHeffAllocated = true;
      which = new char[2]; // which eigenvalues to compute : the nev smallest ones
      which[0] = 'S';
      which[1] = 'A';
      SolveHeffVectorSize = n;
      resid = new double[SolveHeffVectorSize];  //residual vector; ini guess = previous A
      v = new double[SolveHeffVectorSize*MAGICNUMBER];  //space for the arnoldi vectors
      iparam = new int[11];
      ipntr = new int[11]; // pointer to memory
      workd = new double[3*SolveHeffVectorSize]; // work space for the package, not to be used by users during the optim.
      workl = new double[lworkl]; // work space nr. 2
   } else {
      if (n>SolveHeffVectorSize){
         SolveHeffVectorSize = n;
         delete [] resid;
         delete [] v;
         delete [] workd;
         resid = new double[SolveHeffVectorSize];
         v = new double[SolveHeffVectorSize*MAGICNUMBER];
         workd = new double[3*SolveHeffVectorSize];
      }
   }
   dcopy_(&n,the2siteObject->gStorage(),&inc,resid,&inc);

   iparam[0] = 1; // use exact shifts
   iparam[2] = 1000; // maximum number of iterations
   iparam[6] = 1; // standard eigenvalue problem A*x = l*x

   int info = 1; // we use a good initial guess @ resid 1=good guess

   dsaupd_(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

   while (ido != 99){ // ido 99 means convergence has been reached
   
             //vec                 //matvec
      matvec(workd + ipntr[0] - 1, workd + ipntr[1] - 1, index);
      dsaupd_(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

   }

   int rvec = 0; // compute eigenvectors
   bool * rvec2 = (bool *) &rvec; //fortran <-> C++ mismatch of integer size
   *rvec2 = true;
   char howmny = 'A'; // calculate all of the nev eigenvectors
   int select[ncv]; // workspace for reordering the eigenvalues
   bool * select2 = (bool *) select; //fortran <-> C++ mismatch of integer size
   double d; // array containing the nev (here one) eigenvalues

   int ldz = n; //leading dimension Z-array = n x nev array containing the ritz vectors --> store it in A
   double sigma; //unreferenced pm

   dseupd_(rvec2,&howmny,select2,&d,the2siteObject->gStorage(),&ldz,&sigma,&bmat,&n,which,&nev,&tol,resid,&ncv,v,&ldv,iparam,ipntr,workd,workl,&lworkl,&info);

   if (info != 0){
      cout << "dseupd has exited with an error, program aborted" << endl;
      cout << "check the dseupd manual for error code: " << info << endl;
      exit(1);
   }
   
   return d;

}

void DMRG::matvec(double * vec, double * resultvec, const int index){

   int dimLsmall = MPSsolution->gDimAtBound(index);
   int dimRsmall = MPSsolution->gDimAtBound(index+2);
   const int dimLmpo = theMPO->dimL(index);
   const int dimMmpo = theMPO->dimR(index);
   const int dimRmpo = theMPO->dimR(index+1);
   const int phys_d = theMPO->gPhys_d();
   
   int totalDimTemp = phys_d*phys_d*dimLsmall*dimRsmall;
   MPOMPS->checkWork2(totalDimTemp);
   MPOMPS->checkWork3(totalDimTemp);
   
   #pragma omp parallel for schedule(static) default(none) shared(vec, resultvec, dimLsmall, dimRsmall)
   for (int combined=0; combined<phys_d*phys_d; combined++){
      int physLup = combined % phys_d;
      int physRup = combined / phys_d;
      
      int dimBlk = dimLsmall * dimRsmall;
      double * target = resultvec + dimBlk*combined;
      for (int cnt=0; cnt<dimBlk; cnt++){ target[cnt] = 0.0; }
      double * work2 = MPOMPS->gWork2() + combined*dimBlk;
      double * work3 = MPOMPS->gWork3() + combined*dimBlk;
      
      for (int MPOleft=0; MPOleft<dimLmpo; MPOleft++){
         for (int MPOmiddle=0; MPOmiddle<dimMmpo; MPOmiddle++){
            Operator * theOpLeft = theMPO->gOperator(index, MPOleft, MPOmiddle);
            if (!(theOpLeft->AmIOp0())){
               for (int MPOright=0; MPOright<dimRmpo; MPOright++){
                  Operator * theOpRight = theMPO->gOperator(index+1, MPOmiddle, MPOright);
                  if (!(theOpRight->AmIOp0())){
                     double factor = theMPO->gPrefactor(index, MPOleft, MPOmiddle) * theMPO->gPrefactor(index+1, MPOmiddle, MPOright);
                     for (int combined2=0; combined2<phys_d*phys_d; combined2++){
                        int physLdown = combined2 % phys_d;
                        int physRdown = combined2 / phys_d;
                        double OpValue = (*theOpLeft)(physLup,physLdown) * (*theOpRight)(physRup,physRdown) * factor;
                        if (OpValue != 0.0){
                           double * source = vec + dimBlk*combined2;
                           
                           //OpValue * leftedge[block] * source[phys_downL,phys_downR] --> work2 (set)
                           if (index==0){
                              int inc = 1;
                              dcopy_(&dimBlk, source, &inc, work2, &inc);
                              dscal_(&dimBlk, &OpValue, work2, &inc);
                           } else {
                              char notrans = 'N';
                              double one = 1.0;
                              double zero = 0.0;
                              double * boundaryBlock = boundaryTerms[index-1] + dimLsmall*dimLsmall*MPOleft;
                              dgemm_(&notrans,&notrans,&dimLsmall,&dimRsmall,&dimLsmall,&OpValue,boundaryBlock,&dimLsmall,source,&dimLsmall,&zero,work2,&dimLsmall);
                           }

                           //1 * work2 * rightedge[block] --> target (add)                 
                           if (index==length-2){
                              int inc = 1;
                              double one = 1.0;
                              daxpy_(&dimBlk, &one, work2, &inc, target, &inc);
                           } else {
                              char notrans = 'N';
                              char trans = 'T';
                              double one = 1.0;
                              double * boundaryBlock = boundaryTerms[index+1] + dimRsmall*dimRsmall*MPOright;
                              dgemm_(&notrans,&trans,&dimLsmall,&dimRsmall,&dimRsmall,&one,work2,&dimLsmall,boundaryBlock,&dimRsmall,&one,target,&dimLsmall);
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





