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
   
   while (fabs(EnergyPrevious-Energy)>1e-10){
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
      /*Energy = SolveHeff(index);*/
      Energy = SolveDAVIDSON(index);
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
      /*Energy = SolveHeff(index);*/
      Energy = SolveDAVIDSON(index);
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

/*double DMRG::SolveHeff(const int index){

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

}*/

double DMRG::SolveDAVIDSON(const int index){

   //From people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter11.pdf : algorithm 11.1, with instead of line (16), equation (11.3).
   const int DAVIDSON_NUM_VEC            = 32;
   const int DAVIDSON_NUM_VEC_KEEP       = 3;
   const double DAVIDSON_PRECOND_CUTOFF  = 1e-12;
   const double DAVIDSON_RTOL_BASE       = 1e-10;

   int length_vec = MPSsolution->gDimAtBound(index) * MPSsolution->gPhys_d() * MPSsolution->gPhys_d() * MPSsolution->gDimAtBound(index+2);
   int num_vec = 0;
   double ** vecs  = new double*[DAVIDSON_NUM_VEC];
   double ** Hvecs = new double*[DAVIDSON_NUM_VEC];
   int num_allocated = 0;
   
   double * mxM = new double[DAVIDSON_NUM_VEC * DAVIDSON_NUM_VEC];
   double * mxM_eigs = new double[DAVIDSON_NUM_VEC];
   double * mxM_vecs = new double[DAVIDSON_NUM_VEC * DAVIDSON_NUM_VEC];
   int mxM_lwork = 3*DAVIDSON_NUM_VEC-1;
   double * mxM_work = new double[mxM_lwork];
   
   double rtol = DAVIDSON_RTOL_BASE * sqrt(length_vec);
   double rnorm = 10*rtol;
   
   double * t_vec = new double[length_vec];
   double * u_vec = new double[length_vec];
   double * work_vec = new double[length_vec];
   int inc1 = 1;
   dcopy_(&length_vec,the2siteObject->gStorage(),&inc1,t_vec,&inc1); //starting vector for Davidson is the current state of the Sobject in symmetric conventions.
   
   //Checking whether the S-object contains anything
   double Sobjectnorm = 0.0;
   for (int cnt=0; cnt<length_vec; cnt++){ Sobjectnorm += t_vec[cnt]*t_vec[cnt]; }
   Sobjectnorm = sqrt(Sobjectnorm);
   if (Sobjectnorm==0.0){
      for (int cnt=0; cnt<length_vec; cnt++){ t_vec[cnt] = ((double) rand())/RAND_MAX; }
   }
   //End checking whether the S-object contains anything
   
   double * HeffDiag = new double[length_vec];
   fillDiag(HeffDiag, index);
   
   double * Reortho_Lowdin = NULL;
   double * Reortho_Overlap_eigs = NULL;
   double * Reortho_Overlap = NULL;
   double * Reortho_Eigenvecs = NULL;
   bool Reortho_Allocated = false;
   
   int nIterations = 0;
   
   while (rnorm > rtol){

      //1. Orthogonalize the new t_vec w.r.t. the old basis
      for (int cnt=0; cnt<num_vec; cnt++){
         double min_overlap = - ddot_(&length_vec,t_vec,&inc1,vecs[cnt],&inc1);
         daxpy_(&length_vec,&min_overlap,vecs[cnt],&inc1,t_vec,&inc1);
      }
   
      //2. Normalize the t_vec
      char norm = 'F';
      double alpha = 1.0/dlange_(&norm,&length_vec,&inc1,t_vec,&length_vec,t_vec); //work not referenced as Frobenius norm
      dscal_(&length_vec,&alpha,t_vec,&inc1);
      
      //3. T_vec becomes part of vecs
      if (num_vec<num_allocated){
         double * temp = vecs[num_vec];
         vecs[num_vec] = t_vec;
         t_vec = temp;
      } else {
         vecs[num_allocated] = t_vec;
         Hvecs[num_allocated] = new double[length_vec];
         t_vec = new double[length_vec];
         num_allocated++;
      }
      matvec(vecs[num_vec], Hvecs[num_vec], index);
      nIterations++;
      
      //4. mxM contains the Hamiltonian in the basis "vecs"
      for (int cnt=0; cnt<num_vec; cnt++){
         mxM[cnt + DAVIDSON_NUM_VEC * num_vec] = ddot_(&length_vec,vecs[num_vec],&inc1,Hvecs[cnt],&inc1);
         mxM[num_vec + DAVIDSON_NUM_VEC * cnt] = mxM[cnt + DAVIDSON_NUM_VEC * num_vec];
      }
      mxM[num_vec + DAVIDSON_NUM_VEC * num_vec] = ddot_(&length_vec,vecs[num_vec],&inc1,Hvecs[num_vec],&inc1);
      
      //5. When t-vec was added to vecs, the number of vecs was actually increased by one. For convenience (doing 4.), only now the number is incremented.
      num_vec++;
      
      //6. Calculate the eigenvalues and vectors of mxM
      char jobz = 'V';
      char uplo = 'U';
      int info;
      for (int cnt1=0; cnt1<num_vec; cnt1++){
         for (int cnt2=0; cnt2<num_vec; cnt2++){
            mxM_vecs[cnt1 + DAVIDSON_NUM_VEC * cnt2] = mxM[cnt1 + DAVIDSON_NUM_VEC * cnt2];
         }
      }
      int lda = DAVIDSON_NUM_VEC;
      dsyev_(&jobz,&uplo,&num_vec,mxM_vecs,&lda,mxM_eigs,mxM_work,&mxM_lwork,&info); //ascending order of eigs
      
      //7. Calculate u and r. r is stored in t_vec, u in u_vec.
      for (int cnt=0; cnt<length_vec; cnt++){
         t_vec[cnt] = 0.0;
         u_vec[cnt] = 0.0;
      }
      for (int cnt=0; cnt<num_vec; cnt++){
         double alpha = mxM_vecs[cnt]; //eigenvector with lowest eigenvalue, hence mxM_vecs[cnt + DAVIDSON_NUM_VEC * 0]
         daxpy_(&length_vec,&alpha,Hvecs[cnt],&inc1,t_vec,&inc1);
         daxpy_(&length_vec,&alpha, vecs[cnt],&inc1,u_vec,&inc1);
      }
      alpha = -mxM_eigs[0];
      daxpy_(&length_vec,&alpha,u_vec,&inc1,t_vec,&inc1);
      
      //8. Calculate the norm of r
      rnorm = dlange_(&norm,&length_vec,&inc1,t_vec,&length_vec,t_vec);
      
      //9. In case convergence is not yet reached: prepare for the following iteration
      if (rnorm > rtol){
      
         //9a. Calculate the new t_vec based on the residual of the lowest eigenvalue, to add to the vecs.
         for (int cnt=0; cnt<length_vec; cnt++){
            if (fabs(HeffDiag[cnt] - mxM_eigs[0])> DAVIDSON_PRECOND_CUTOFF ){
               work_vec[cnt] = u_vec[cnt]/(HeffDiag[cnt] - mxM_eigs[0]); // work_vec = K^(-1) u_vec
            } else {
               work_vec[cnt] = u_vec[cnt]/DAVIDSON_PRECOND_CUTOFF ;
               if (true) cout << "|(HeffDiag[" << cnt << "] - mxM_eigs[0])| = " << fabs(HeffDiag[cnt] - mxM_eigs[0]) << endl;
            }
         }
         alpha = - ddot_(&length_vec,work_vec,&inc1,t_vec,&inc1)/ddot_(&length_vec,work_vec,&inc1,u_vec,&inc1); // alpha = - (u^T K^(-1) r) / (u^T K^(-1) u)
         daxpy_(&length_vec,&alpha,u_vec,&inc1,t_vec,&inc1); // t_vec = r - (u^T K^(-1) r) / (u^T K^(-1) u) u
         for (int cnt=0; cnt<length_vec; cnt++){
            if (fabs(HeffDiag[cnt] - mxM_eigs[0])> DAVIDSON_PRECOND_CUTOFF ){
               t_vec[cnt] = - t_vec[cnt]/(HeffDiag[cnt] - mxM_eigs[0]); //t_vec = - K^(-1) (r - (u^T K^(-1) r) / (u^T K^(-1) u) u)
            } else {
               t_vec[cnt] = - t_vec[cnt]/ DAVIDSON_PRECOND_CUTOFF ;
            }
         }
         
         // 9b. When the maximum number of vectors is reached: construct the one with lowest eigenvalue & restart
         if (num_vec == DAVIDSON_NUM_VEC){
         
            if (DAVIDSON_NUM_VEC_KEEP<=1){
            
               alpha = 1.0/dlange_(&norm,&length_vec,&inc1,u_vec,&length_vec,u_vec); //work not referenced as Frobenius norm
               dscal_(&length_vec,&alpha,u_vec,&inc1);
               dcopy_(&length_vec,u_vec,&inc1,vecs[0],&inc1);
               matvec(vecs[0], Hvecs[0], index);
               nIterations++;
               mxM[0] = ddot_(&length_vec,vecs[0],&inc1,Hvecs[0],&inc1);
            
               num_vec = 1;
            
            } else {
            
               //Construct the lowest DAVIDSON_NUM_VEC_KEEP eigenvectors
               if (!Reortho_Allocated) Reortho_Eigenvecs = new double[length_vec * DAVIDSON_NUM_VEC_KEEP];
               dcopy_(&length_vec,u_vec,&inc1,Reortho_Eigenvecs,&inc1);
               for (int cnt=1; cnt<DAVIDSON_NUM_VEC_KEEP; cnt++){
                  for (int irow=0; irow<length_vec; irow++){
                     Reortho_Eigenvecs[irow + length_vec * cnt] = 0.0;
                     for (int ivec=0; ivec<DAVIDSON_NUM_VEC; ivec++){
                        Reortho_Eigenvecs[irow + length_vec * cnt] += vecs[ivec][irow] * mxM_vecs[ivec + DAVIDSON_NUM_VEC * cnt];
                     }
                  }
               }
               
               //Reorthonormalize them
               //Reortho: Calculate the overlap matrix
               if (!Reortho_Allocated) Reortho_Overlap = new double[DAVIDSON_NUM_VEC_KEEP * DAVIDSON_NUM_VEC_KEEP];
               char trans = 'T';
               char notr = 'N';
               int DVDS_KEEP = DAVIDSON_NUM_VEC_KEEP;
               double one = 1.0;
               double zero = 0.0; //set
               dgemm_(&trans,&notr,&DVDS_KEEP,&DVDS_KEEP,&length_vec,&one,Reortho_Eigenvecs,&length_vec, Reortho_Eigenvecs,&length_vec,&zero,Reortho_Overlap,&DVDS_KEEP);
               
               //Reortho: Calculate the Lowdin tfo
               if (!Reortho_Allocated) Reortho_Overlap_eigs = new double[DVDS_KEEP];
               dsyev_(&jobz,&uplo,&DVDS_KEEP,Reortho_Overlap,&DVDS_KEEP,Reortho_Overlap_eigs,mxM_work,&mxM_lwork,&info); //ascending order of eigs
               for (int icnt=0; icnt<DVDS_KEEP; icnt++){
                  Reortho_Overlap_eigs[icnt] = pow(Reortho_Overlap_eigs[icnt],-0.25);
                  dscal_(&DVDS_KEEP, Reortho_Overlap_eigs+icnt, Reortho_Overlap+DVDS_KEEP*icnt, &inc1);
               }
               if (!Reortho_Allocated) Reortho_Lowdin = new double[DVDS_KEEP*DVDS_KEEP];
               dgemm_(&notr,&trans,&DVDS_KEEP,&DVDS_KEEP,&DVDS_KEEP,&one,Reortho_Overlap,&DVDS_KEEP,Reortho_Overlap,&DVDS_KEEP,&zero,Reortho_Lowdin,&DVDS_KEEP);
               
               //Reortho: Put the Lowdin tfo eigenvecs in vecs
               for (int ivec=0; ivec<DAVIDSON_NUM_VEC_KEEP; ivec++){
                  dscal_(&length_vec,&zero,vecs[ivec],&inc1);
                  for (int ivec2=0; ivec2<DAVIDSON_NUM_VEC_KEEP; ivec2++){
                     daxpy_(&length_vec,Reortho_Lowdin + ivec2 + DAVIDSON_NUM_VEC_KEEP * ivec,Reortho_Eigenvecs + length_vec * ivec2, &inc1, vecs[ivec], &inc1);
                  }
               }
               
               if (!Reortho_Allocated) Reortho_Allocated = true;
               
               //Construct the H*vecs
               for (int cnt=0; cnt<DAVIDSON_NUM_VEC_KEEP; cnt++){
                  matvec(vecs[cnt], Hvecs[cnt], index);
                  nIterations++;
               }
               
               //Build MxM
               for (int ivec=0; ivec<DAVIDSON_NUM_VEC_KEEP; ivec++){
                  for (int ivec2=ivec; ivec2<DAVIDSON_NUM_VEC_KEEP; ivec2++){
                     mxM[ivec + DAVIDSON_NUM_VEC * ivec2] = ddot_(&length_vec, vecs[ivec], &inc1, Hvecs[ivec2], &inc1);
                     mxM[ivec2 + DAVIDSON_NUM_VEC * ivec] = mxM[ivec + DAVIDSON_NUM_VEC * ivec2];
                  }
               }
               
               //Set num_vec
               num_vec = DAVIDSON_NUM_VEC_KEEP;
            
            }

         }
      }
      
   }
   
   if (true) cout << "   Stats: nIt(DAVIDSON) = " << nIterations << endl;
   
   double eigenvalue = mxM_eigs[0]; //mxM_eigs[0] and u_vec are eigenvalue and vector
   dcopy_(&length_vec,u_vec,&inc1,the2siteObject->gStorage(),&inc1);
   
   for (int cnt=0; cnt<num_allocated; cnt++){
      delete [] vecs[cnt];
      delete [] Hvecs[cnt];
   }
   delete [] vecs;
   delete [] Hvecs;
   delete [] t_vec;
   delete [] u_vec;
   delete [] work_vec;
   delete [] mxM;
   delete [] mxM_eigs;
   delete [] mxM_vecs;
   delete [] mxM_work;
   delete [] HeffDiag;
   
   if (Reortho_Allocated){
      delete [] Reortho_Eigenvecs;
      delete [] Reortho_Overlap;
      delete [] Reortho_Overlap_eigs;
      delete [] Reortho_Lowdin;
   }
   
   return eigenvalue;

}

void DMRG::fillDiag(double * diag, const int index){

   const int dimLsmall = MPSsolution->gDimAtBound(index);
   const int dimRsmall = MPSsolution->gDimAtBound(index+2);
   const int dimLmpo = theMPO->dimL(index);
   const int dimMmpo = theMPO->dimR(index);
   const int dimRmpo = theMPO->dimR(index+1);
   const int phys_d = theMPO->gPhys_d();
   
   #pragma omp parallel for schedule(static) default(none) shared(diag)
   for (int combined=0; combined<phys_d*phys_d; combined++){
      const int physL = combined % phys_d;
      const int physR = combined / phys_d;
      
      int dimBlk = dimLsmall * dimRsmall;
      double * target = diag + dimBlk*combined;
      for (int cnt=0; cnt<dimBlk; cnt++){ target[cnt] = 0.0; }
      
      for (int MPOleft=0; MPOleft<dimLmpo; MPOleft++){
         for (int MPOmiddle=0; MPOmiddle<dimMmpo; MPOmiddle++){
            Operator * theOpLeft = theMPO->gOperator(index, MPOleft, MPOmiddle);
            const double leftOperatorElement = (*theOpLeft)(physL,physL);
            if (leftOperatorElement != 0.0){
               for (int MPOright=0; MPOright<dimRmpo; MPOright++){
                  Operator * theOpRight = theMPO->gOperator(index+1, MPOmiddle, MPOright);
                  const double prefactor = theMPO->gPrefactor(index, MPOleft, MPOmiddle) * theMPO->gPrefactor(index+1, MPOmiddle, MPOright)
                                         * leftOperatorElement * (*theOpRight)(physR,physR);

                  if (prefactor != 0.0){
                     for (int alpha_left=0; alpha_left<dimLsmall; alpha_left++){
                        double leftFactor = 1.0;
                        if (index>0){ leftFactor = boundaryTerms[index-1][dimLsmall*dimLsmall*MPOleft + alpha_left*(dimLsmall+1)]; }
                        if (leftFactor != 0.0){
                           for (int alpha_right=0; alpha_right<dimRsmall; alpha_right++){
                              double rightFactor = 1.0;
                              if (index<length-2){ rightFactor = boundaryTerms[index+1][dimRsmall*dimRsmall*MPOright + alpha_right*(dimRsmall+1)]; }
                              target[alpha_left + dimLsmall*alpha_right] += prefactor * leftFactor * rightFactor;
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
