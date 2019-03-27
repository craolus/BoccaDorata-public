/*
 *  BcDor-Self_energy_ext_Op_e.cpp
 *  
//
//   Construct and store the self energy (both HF and dynamc part) for a given
//  channel taking as input:
//      model spacce, interaction and sp propagator.
//
 *
 *  Created by Barbieri Carlo on 10/24/10.
 *  Copyright 2010 RIKEN. All rights reserved.
 *
 */


#include <iostream>
//#include <fstream>
//#include <iomanip>
using namespace std;

#include <cstdlib>
//#include <cstdio>
//#include <cmath>

//#include "BcDor-Global_variables.hh"
#include "BcDor-Run_vars.hh"
#include "BcDor-SlfEn_classes.hh"
#include "BcDor-Utilities.hh"


static int n1;

//static int i_pvt,j_next_pvt;

static int NDIMHF, NFWP, NBKP, NLCZSE;  // HFE, fw poles, bk poles dims. and size of the LancSelfEn
//static int NSOLDIM;                     // dims of sols. arrays.
//static double *LancSlE_Fw, *LancSlE_Bk; //  Lanczcos mtx. els and overlap w/
//// the EHF self-energy.
//static double *LancVect;        // Store Lanczos vectors

//static double *piv;   // pivots
//static int *i_last_fpvt;//, *i_used_fpvt;
//static int *i_last_bpvt;//, *i_used_bpvt;

//static double *ptr_i, *ptr_ms, *ptr_j, *ptr_mlnc_e, *ptr_mlnc_m;
static double *ptr_ms;
//static double *ptr1, *ptr2, *ptr3, *ptr_up, *ptr_lo;
//static double xa, xb, xg; // xe
//static double x1, x2, x3, x4;

//static double EFermi;



int Wide_self_energy_t::Reduce_SelfEn_Lncz(Wide_self_energy_t *SlfE_in,
                          int nt_pvts, double** pivots, int* n_itr_fw, int* n_itr_bk,
                          int i_copyMF/*=0*/) {
  //
  //  This makes a Lanczos reduction fo the input self-energy SlfE_in based on the
  // pivots and number of iterations passed. Results are then stored in 'this'
  // object.
  //
  //  NOTE the 'this'  and  SlfE_in  **must** be two different objects or the algorithm will crash
  //
  //  nt_pvts    number of pivots
  //
  //  pivots:    pointer to arrays for each pivot. element (v^n)_i is pivots[n][i]  where  0<= n < nt_pvts
  //
  //  n_itr_fw,  n_itr_bk:  number of fw and bk iteration of each pivot,  n_itr_fw[n][i]  with  0<= n < nt_pvts, etc...
  //
  //  i_copyMF, if != 0  then copy all the MF parts to the new (Lanczos reduced) self-energy
  //

  if (SlfE_in == this) {
    //  Check to avoid crashes crash...
    cerr << " ERROR in 'Reduce_SelfEn_Lncz'";
    cerr << ": input and output self-energies must be different...  -->STOP.";
    exit(100);
  }

  if ((NTOT_ORB != SlfE_in->NTOT_ORB) || (MdSp != SlfE_in->MdSp) ||
      (  Vpp    != SlfE_in->Vpp)      || (gin  != SlfE_in->gin)  || (I_SH != SlfE_in->I_SH)) {

    cout << " WARNING in Wide_self_energy_t (Reduce_SelfEn_Lncz): The two self-energies do dot agree"
            "in the mod. sp., partial wave, or other parameters...\n";
    cout << "\n  -- I will reset the output self-energy completely! (and also copy the MF)\n";
    reset(SlfE_in->I_SH, SlfE_in->Lambda, SlfE_in->parity, SlfE_in->charge,
          SlfE_in->MdSp, SlfE_in->Vpp, SlfE_in->gin);
    i_copyMF = 1;
  }



  NDIMHF =  SlfE_in->NTOT_ORB;
  NFWP   =  SlfE_in->N_PLS_fw;  //NDIMFW = NDIMHF + NFWP;
  NBKP   =  SlfE_in->N_PLS_bk;  //NDIMBK = NDIMHF + NBKP;
  NLCZSE =  NDIMHF + nt_pvts + 1;



  //////////////////////////////////////////////////////////////////////////////
  //
  //   preprocess the self energy...
  //
  //if ( NULL == SlfE_in->Sig_HF ) SlfE_in->Build_ExtendedHartreeFock();
  if ( NULL == SlfE_in->Sig_HF ) SlfE_in->Make_empty_ExtendedHartreeFock();
  if ( NULL == SlfE_in->Sig_1b ) SlfE_in->Make_1bd_SelfEn(NULL, -10);

  NDIMHF =  SlfE_in->NTOT_ORB;
  if (NDIMHF <= 0) cout << "Error: how can I build a self energy with no sp orbitals!?!?!?\n";

  //
  // order fw poles in slfen:
  //
  // select fw poles in slfen:
  //
  // final ordering of fw poles in slfen (according to ascending energies):
  //
  ptr_ms = SlfE_in->Sig_dyn_fw;
  Sort_u_double2dim(NFWP, SlfE_in->NTOT_DIM, 0, ptr_ms);
  
  //
  // order bk poles in slfen:
  //
  // select bk poles in slfen:
  //
  // final ordering of bk poles in slfen (according to descending energies):
  //
  ptr_ms = SlfE_in->Sig_dyn_bk;
  Sort_d_double2dim(NBKP, SlfE_in->NTOT_DIM, 0, ptr_ms);
  //
  //////////////////////////////////////////////////////////////////////////////


  // Plot and/or chop the input self energy (for TESTS):
  //
  ////  NFWP = 5;
  //  ptr_ms = SlfE_in->Sig_dyn_fw;
  //  for(i= (30 < NFWP) ? 29 : NFWP-1;  i >= 0; --i) {
  //    ptr_ms = SlfE_in->Sig_dyn_fw + i*(SlfE_in->NTOT_DIM);
  //    cout << i << "   ";
  //    for(j=0; j<SlfE_in->NTOT_DIM; ++j) {cout << (*ptr_ms) << "   "; ++ptr_ms;}
  //    cout << endl;
  //    }
  ////  NBKP = 5;
  //  ptr_ms = SlfE_in->Sig_dyn_bk;
  //  for(i=0;  (i<30)&&(i<NBKP); ++i) {
  //    cout << i << "   ";
  //    for(j=0; j<SlfE_in->NTOT_DIM; ++j) {cout << (*ptr_ms) << "   "; ++ptr_ms;}
  //    cout << endl;
  //    }


  //
  // Compute the tot. # of fw and bk vectors:
  
  static int nLancFw, nLancBk;
  static int nTotItrFw, nTotItrBk;
  

  nLancFw = 0;
  nLancBk = 0;
  for(n1=0; n1<nt_pvts; ++n1) {
	if (n_itr_fw[n1] > 0) nLancFw += n_itr_fw[n1];
	if (n_itr_bk[n1] > 0) nLancBk += n_itr_bk[n1];
  }
	
  nTotItrFw = nLancFw;
  nTotItrBk = nLancBk;
	
  if (nLancFw < 0) nLancFw=abs(SlfE_in->MdSp->MSp_nqp[SlfE_in->I_SH]);
  if (nLancBk < 0) nLancBk=abs(SlfE_in->MdSp->MSp_nqh[SlfE_in->I_SH]);
	
	
  if (nLancFw > NFWP) nLancFw = NFWP;
  if (nLancBk > NBKP) nLancBk = NBKP;

  cout << "nLancFw  = " << nLancFw << endl;
  cout << "nLancBk  = " << nLancBk << endl;
  cout << "NDIMHF   = " << NDIMHF  << endl;
  cout << "NFWP     = " << NFWP    << endl;
  cout << "NBKP     = " << NBKP    << endl;
  cout << "NLCZSE   = " << NLCZSE  << endl;

  cout << "Estimate of memory needs:\n";
  cout << "  Fw Lanczos vectors: " << (  NFWP*nLancFw/128)+1 << " kbytes\n";
  cout << "  Fw Lanczos self-en: " << (NLCZSE*nLancFw/128)+1 << " kbytes\n";
  cout << "  Bk Lanczos vectors: " << (  NBKP*nLancBk/128)+1 << " kbytes\n";
  cout << "  Bk Lanczos self-en: " << (NLCZSE*nLancBk/128)+1 << " kbytes\n";
	

  this->Allocate_Dynamic_SlfEn(nLancFw, nLancBk, NLCZSE-NDIMHF);


L_fw:
  if ((NFWP < 1)  || (nLancFw < 1)) goto L_bk;
  //if ((NFWP > 0)  && (nLancFw > 0)) {

  if (nTotItrFw > NFWP) {
    cout << "\n\n  Too many Fw iterations requested, I will rescale in such a way all pivots are iterated";
    double x1 = double(NFWP)/ double(nTotItrFw);
    for(n1=0; n1<nt_pvts; ++n1) n_itr_fw[n1] = int(x1 * n_itr_fw[n1] + 1.0);
  }

  cout << "\nRunning forward Lanczos:\n========================\n";
  N_PLS_fw = Reduce_Mtx_Lncz(NTOT_DEN, NTOT_ORB, N_PLS_FW_alloc, Sig_dyn_fw, i_piv_fw,
							 SlfE_in->NTOT_DEN, SlfE_in->NTOT_ORB, SlfE_in->N_PLS_fw,
							 SlfE_in->Sig_dyn_fw, nt_pvts, pivots, n_itr_fw);
  n_piv_fw = i_piv_fw[0];


L_bk:
  if ((NBKP < 1)  || (nLancBk < 1)) goto end;
  //if ((NBKP > 0)  && (nLancBk > 0)) {
  
  if (nTotItrBk > NBKP) {
    cout << "\n\n  Too many Bk iterations requested, I will rescale in such a way all pivots are iterated";
    double x1 = double(NBKP)/double(nTotItrBk);
    for(n1=0; n1<nt_pvts; ++n1) n_itr_bk[n1] = int(x1 * n_itr_bk[n1] + 1.0);
  }
  
  cout << "\nRunning backward Lanczos:\n=========================\n";
  N_PLS_bk = Reduce_Mtx_Lncz(NTOT_DEN, NTOT_ORB, N_PLS_BK_alloc, Sig_dyn_bk, i_piv_bk,
							 SlfE_in->NTOT_DEN, SlfE_in->NTOT_ORB, SlfE_in->N_PLS_bk,
							 SlfE_in->Sig_dyn_bk, nt_pvts, pivots, n_itr_bk);
  n_piv_bk = i_piv_bk[0];


end:
  if (i_copyMF) Copy_MF(SlfE_in);

	return 0;}


int Wide_self_energy_t::Copy_MF(Wide_self_energy_t *SlfE_in) {
  
  this->free_1b_mem();
  
  Make_1bd_SelfEn(SlfE_in->Sig_1b, SlfE_in->NTOT_ORB);
  Make_empty_ExtendedHartreeFock();
  int ndim_loc = (NTOT_ORB < SlfE_in->NTOT_ORB) ? NTOT_ORB : SlfE_in->NTOT_ORB;
  for(int nc=0; nc<ndim_loc; ++nc)
    for(int nr=0; nr<ndim_loc; ++nr)
      Sig_HF[nr*NTOT_ORB+nc] += SlfE_in->Sig_HF[nr*SlfE_in->NTOT_ORB+nc];
  
  N_GMTX = SlfE_in->N_GMTX;
  if (N_GMTX > 0) {
    //should copy BHF pot.
  }
  
  return 0;}
