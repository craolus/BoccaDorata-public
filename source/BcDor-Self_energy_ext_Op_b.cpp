//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//
//
//  the following is to be placed in a second file and be adapted to the new
//  Wide_self-enerrgy class (...and tested!!).  It will eventually replace
// the old  Self-enerrgy class.
//
//


#include <cstdlib>
//#include <cstdio>
//#include <cmath>
 #include <iostream>
// #include <fstream>
// #include <iomanip>
using namespace std;


#include "BcDor-SlfEn_classes.hh"
#include "BcDor-Utilities.hh"



static int i,nr,nc;
static int ndim_loc;

// resets the Sig_HF to zero (and allocates it, if necessary):
void Wide_self_energy_t::Make_empty_1bd_SelfEn(void ) {
  if ( NULL == Sig_1b ) Sig_1b = new double[NTOT_ORB*NTOT_ORB];
  for(i=0; i<NTOT_ORB*NTOT_ORB; i++) Sig_1b[i]=0.;
  return;
  }


//  initializes Sig_1b and allocates it (if necessary). If no external
// poential is provided, sets it to zero:
void Wide_self_energy_t::Make_1bd_SelfEn(double Uab[], int NDIM) {
  Make_empty_1bd_SelfEn();
  Add_1bd_SelfEn(Uab, NDIM);
  return;
  }


//  Add an external poential to Sig_1b--if it is not initialized, it does so:
void Wide_self_energy_t::Add_1bd_SelfEn(double Uab[], int NDIM) {

  if ( NULL == Sig_1b ) Make_empty_1bd_SelfEn();

  if ((NULL == Uab ) || (NDIM < 1) ) return; // no 1-body potential is provided

  ndim_loc = (NTOT_ORB < NDIM) ? NTOT_ORB : NDIM;
  for(nc=0; nc<ndim_loc; ++nc)
   for(nr=0; nr<ndim_loc; ++nr)
     Sig_1b[nr*NTOT_ORB+nc] += Uab[nr*NDIM+nc];

  return;
  }


// resets the Sig_HF to zero (and allocates it, if necessary):
void Wide_self_energy_t::Make_empty_ExtendedHartreeFock(void ) {
  if ( NULL == Sig_HF ) Sig_HF = new double[NTOT_ORB*NTOT_ORB];
  for(i=0; i<NTOT_ORB*NTOT_ORB; i++) Sig_HF[i]=0.;
  return;
  }


void Wide_self_energy_t::Build_ExtendedHartreeFock(void ) {
  Make_empty_ExtendedHartreeFock();  // reset the Sig_HF to zero
  i = gin->ExtendedHartreeFock_Pot(Sig_HF, NTOT_ORB, Vpp, I_SH);
  if (i) cerr << "\nWarning: the calculation of the (E)HF potential gave an error (" << i << ")\n\n";
  if ((0 != Lambda) || (0 != parity) || (0 != charge)) {
    cout << "Error: cannot (yet) perform (B)HF on a nucleus that has not 0+ g.s.!!!!\n";
    exit(1000);
  }
  return;
  }

// 
// static int i1, i2;
// static double *dptr_1, *dptr_2;
// static double dx;
// static double *Ufw, *Ubk;
// 
// void Wide_self_energy_t::Evaluate_Dyn_SelfEn(double Uab[], int NDIM, double xe,
//                                         int ifb/*==0*/,double xelim/*==+1.e100*/) {  //, int isort/*==0*/) {
// 
//   if ((NULL == Uab ) || (NDIM < 1) ) return; // no 1-body potential is provided
// 
//   //if (i_sort) {
//     Sort_d_double2dim(N_PLS_fw, NTOT_DIM, 0, Sig_dyn_fw);
//     Sort_u_double2dim(N_PLS_bk, NTOT_DIM, 0, Sig_dyn_bk);
//   //}
// 
//   i1 = ((NTOT_ORB+1)*NTOT_ORB)/2;
//   Ufw = new double[i1];
//   Ubk = new double[i1];
// 
//   for(i2=0; i2<i1; ++i2) {Ufw[i2]=0.0; Ubk[i2]=0.0;}
// 
// forward:
//   if (ifb < 0) goto backward;
//   dptr_1 = Sig_dyn_fw;
//   for(i1=0; i1<N_PLS_fw; ++i1) {
//     //if ( (0 < ifb) && ((*dptr_1) < xelim) ) continue;
//     dx = xe - (*dptr_1);
// 
//     dptr_2 = dptr_1 + NTOT_DEN;
//     i2 = 0;
//     for(nc=0; nc<NTOT_ORB; ++nc)
//      for(nr=nc; nr<NTOT_ORB; ++nr) {
//        Ufw[i2] += dptr_2[nc] * dptr_2[nr] / dx;
//        ++i2;
//        }
// 
//     dptr_1 += (NTOT_DIM);
//   }
// 
// backward:
//   if (ifb > 0) goto collect;
//   dptr_1 = Sig_dyn_bk;
//   for(i1=0; i1<N_PLS_bk; ++i1) {
//     dx = xe - (*dptr_1);
// 
//     dptr_2 = dptr_1 + NTOT_DEN;
//     i2 = 0;
//     for(nc=0; nc<NTOT_ORB; ++nc)
//      for(nr=nc; nr<NTOT_ORB; ++nr) {
//        Ubk[i2] += dptr_2[nc] * dptr_2[nr] / dx;
//        ++i2;
//        }
// 
//     dptr_1 += (NTOT_DIM);
//   }
// 
// collect:
//   i2 = 0;
//   for(nc=0; nc<NTOT_ORB; ++nc)
//    for(nr=nc; nr<NTOT_ORB; ++nr) {
//      Uab[nr*NDIM+nc] = Ufw[i2] + Ubk[i2];
//      Uab[nc*NDIM+nr] = Ufw[i2] + Ubk[i2];
//      ++i2;
//      }
// 
//   delete [] Ufw; Ufw = NULL;
//   delete [] Ubk; Ubk = NULL;
//   return; }
// 
// 
// void Wide_self_energy_t::Evaluate_Der_SelfEn(double Uab[], int NDIM, double xe,
//                                         int ifb/*==0*/,double xelim/*==+1.e100*/) {  //, int isort/*==0*/) {
// 
//   if ((NULL == Uab ) || (NDIM < 1) ) return; // no 1-body potential is provided
// 
//   //if (i_sort) {
//     Sort_d_double2dim(N_PLS_fw, NTOT_DIM, 0, Sig_dyn_fw);
//     Sort_u_double2dim(N_PLS_bk, NTOT_DIM, 0, Sig_dyn_bk);
//   //}
// 
//   i1 = ((NTOT_ORB+1)*NTOT_ORB)/2;
//   Ufw = new double[i1];
//   Ubk = new double[i1];
// 
//   for(i2=0; i2<i1; ++i2) {Ufw[i2]=0.0; Ubk[i2]=0.0;}
// 
// forward:
//   if (ifb < 0) goto backward;
//   dptr_1 = Sig_dyn_fw;
//   for(i1=0; i1<N_PLS_fw; ++i1) {
//     //if ( (0 < ifb) && ((*dptr_1) < xelim) ) continue;
//     dx = xe - (*dptr_1);
//     dx *= dx;
// 
//     dptr_2 = dptr_1 + NTOT_DEN;
//     i2 = 0;
//     for(nc=0; nc<NTOT_ORB; ++nc)
//      for(nr=nc; nr<NTOT_ORB; ++nr) {
//        Ufw[i2] -= dptr_2[nc] * dptr_2[nr] / dx;
//        ++i2;
//        }
// 
//     dptr_1 += (NTOT_DIM);
//   }
// 
// backward:
//   if (ifb > 0) goto collect;
//   dptr_1 = Sig_dyn_bk;
//   for(i1=0; i1<N_PLS_bk; ++i1) {
//     dx = xe - (*dptr_1);
//     dx *= dx;
// 
//     dptr_2 = dptr_1 + NTOT_DEN;
//     i2 = 0;
//     for(nc=0; nc<NTOT_ORB; ++nc)
//      for(nr=nc; nr<NTOT_ORB; ++nr) {
//        Ubk[i2] -= dptr_2[nc] * dptr_2[nr] / dx;
//        ++i2;
//        }
// 
//     dptr_1 += (NTOT_DIM);
//   }
// 
// collect:
//   i2 = 0;
//   for(nc=0; nc<NTOT_ORB; ++nc)
//    for(nr=nc; nr<NTOT_ORB; ++nr) {
//      Uab[nr*NDIM+nc] = Ufw[i2] + Ubk[i2];
//      Uab[nc*NDIM+nr] = Ufw[i2] + Ubk[i2];
//      ++i2;
//      }
// 
//   delete [] Ufw; Ufw = NULL;
//   delete [] Ubk; Ubk = NULL;
//   return; }
// 
