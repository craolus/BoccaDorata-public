/*
 *  BcDor-CCD.C
 *  
 *
 *  Routined for calculating the CCD amplitudes to be later used
 * for generating the ADC(3) and the ADC(CCD) vertices
 *
 *
 *  Created by Carlo Barbieri on 02/10/2011.
 *  Copyright 2011 University of Surrey. All rights reserved.
 *
 */


#include <cstdlib>
#include <iostream>
//#include <fstream>
//#include <iomanip>  //TST//
using namespace std;

//#include <cstdlib>
//#include <cstdio>
#include <cmath>

//#define DPI   3.14159265358979323846
#define SQR2  1.414213562373095048801688


#include "BcDor-Ang_momenta.hh"
//#include "BcDor-HO_radial_me.hh"

#include "BcDor-Global_variables.hh"
#include "BcDor-Run_vars.hh"
#include "BcDor-CCD_classes.hh"



CoupledCluster_t::CoupledCluster_t(ModSpace_t *MdSpin, SpProp_t *gin, VppInt_t *vin) {
  MSp_CC = MdSpin;
  gsp_CC = gin;
  Vpp_CC = vin;
  init_NULL();
  return; }


CoupledCluster_t::~CoupledCluster_t() {
  free_CCD_aplitudes();
  return; }


static int i1;
static int n1, n2, n3;
const char chip[2] = { '+' , '-'};


void CoupledCluster_t::init_NULL(void ) {
  T2_Ld = NULL;
  T2_Pi = NULL;
  T2_Ld_dim = NULL;
  T2_Pi_dim = NULL;
  NT2_CHANNELS = -1;
  //npw_Ld = -1;
  //npw_Pi = -1;
  //NDIM_2bmsp_MX = -1;
  return;}

int  CoupledCluster_t::free_CCD_aplitudes(void ) {
  cout <<" Deleting CCD amplitudes...\n";
  if (NULL != T2_Ld) {
    for(i1=0; i1<NT2_CHANNELS; ++i1) if (NULL != T2_Ld[i1]) delete T2_Ld[i1];
    delete [] T2_Ld;
    T2_Ld = NULL;
  }
  if (NULL != T2_Ld_dim) delete [] T2_Ld_dim;
  //
  if (NULL != T2_Pi) {
    for(i1=0; i1<NT2_CHANNELS; ++i1) if (NULL != T2_Pi[i1]) delete T2_Pi[i1];
    delete [] T2_Pi;
    T2_Pi = NULL;
  }
  if (NULL != T2_Pi_dim) delete [] T2_Pi_dim;

  return 0;}


static gII_prop_t* gIIloc;
static Pol_prop_t* Piloc;

static int J, J2, dT, pi;
static int nf_dim, nb_dim, ndim_tot;
static int i_ld_njpt, i_pi_njpt;
static int n4,n5;

static double *dptr1, *dptr2;


void Correct_Jph_phase(Pol_prop_t* Piloc) {

int nf = Piloc->n_f_basis;
int nb = Piloc->n_b_basis;
int lda = Piloc->XY_nvects;

double *ptr_fb = Piloc->XYmtx + Piloc->NdataXY + nf*lda;
double *ptr_bf = Piloc->XYmtx + Piloc->NdataXY          + nf;
double *ptr_bb = Piloc->XYmtx + Piloc->NdataXY + nf*lda + nf;

if (abs(Piloc->Jf)%2) {
  for (int ifw=0; ifw<nf; ++ifw)
	for (int ibk=0; ibk<nb; ++ibk) {
	  ptr_fb[lda*ibk + ifw] = - ptr_fb[lda*ibk + ifw];
	  ptr_bf[lda*ifw + ibk] = - ptr_bf[lda*ifw + ibk];
  }
}
if (abs(2 * Piloc->Jf)%2) {
  for (int ibk=0;   ibk<nb; ++ibk)
	for (int ibk2=0; ibk2<nb; ++ibk2) {
	  ptr_bb[lda*ibk + ibk2] = - ptr_bb[lda*ibk + ibk2];
	}
}
return;}


int  CoupledCluster_t::create_CCD_aplitudes(int i_plot/*=1*/) {

  free_CCD_aplitudes();

  cout <<"\n Generating new CCD amplitudes...\n";
  
  MSp_CC->Calculate_Jch_bounds(i_plot);
  
  MULT_CH_PI = (MSp_CC->ch_ph_mx - MSp_CC->ch_ph_mn) + 1;
  MULT_CH_LD = (MSp_CC->ch_pp_mx - MSp_CC->ch_pp_mn) + 1;
  MULT_CH = (MULT_CH_PI > MULT_CH_LD) ? MULT_CH_PI : MULT_CH_LD;
  CH_OFFS_PI = MSp_CC->ch_ph_mn;
  CH_OFFS_LD = MSp_CC->ch_pp_mn;
  
  NT2_CHANNELS = (1 + MSp_CC->Jmax)*2 * MULT_CH;

  if (i_plot) {
	cout << "MULT_CH_PI=" << MULT_CH_PI << "     MULT_CH_LD="<< MULT_CH_LD << "    MULT_CH="<< MULT_CH <<endl;
	cout << "CH_OFFS_PI=" << CH_OFFS_PI << "     CH_OFFS_LD="<< CH_OFFS_LD <<endl;
	cout << "NT2_CHANNELS="<< NT2_CHANNELS <<endl;
  }


  if (NULL != T2_Ld) delete [] T2_Ld;  T2_Ld = new gII_prop_t* [NT2_CHANNELS];
  if (NULL != T2_Pi) delete [] T2_Pi;  T2_Pi = new Pol_prop_t* [NT2_CHANNELS];
  if (NULL != T2_Ld_dim) delete [] T2_Ld_dim;  T2_Ld_dim = new int[NT2_CHANNELS];
  if (NULL != T2_Pi_dim) delete [] T2_Pi_dim;  T2_Pi_dim = new int[NT2_CHANNELS];
  for(i1=0; i1<NT2_CHANNELS; ++i1) {
	T2_Ld[i1] = NULL;
	T2_Pi[i1] = NULL;
	T2_Ld_dim[i1] = -1;
	T2_Pi_dim[i1] = -1;
  }


  for(J=0;             J<=MSp_CC->Jmax;  ++J)
  for(pi=0;            pi< 2;            ++pi) {
	for(dT=MSp_CC->ch_pp_mn; dT<=MSp_CC->ch_pp_mx; ++dT) {
	  
	  J2 = 2*J;
	  i_ld_njpt = (J2+pi)*MULT_CH + dT - CH_OFFS_LD;
	  
	  if (NULL != T2_Ld[i_ld_njpt]) delete T2_Ld[i_ld_njpt];
	  T2_Ld[i_ld_njpt] = new gII_prop_t(MSp_CC,gsp_CC,Vpp_CC);
	  T2_Ld_dim[i_ld_njpt] = -1;
	  gIIloc = T2_Ld[i_ld_njpt];
	  gIIloc->Jf  = J;
	  gIIloc->dT  = dT;
	  gIIloc->pif = pi;
	  gIIloc->Count_LaddRPA_basis(J ,dT, pi, &n4, &n5);
	  if (i_plot) cout << "J\\pi,dT = "<<J<<chip[pi]<<" , "<<dT<<"     i_ld_njpt="<<i_ld_njpt<<"     "
					   << "     nf = "<<n4<<"     nb = "<<n5<<endl<<flush;
	  if ((n4 > 0) && (n5 > 0)) {
		  gIIloc->Build_LaddRPA_basis(J ,dT, pi);
		  nf_dim = gIIloc->n_f_basis;
		  nb_dim = gIIloc->n_b_basis;
		  ndim_tot = gIIloc->n_f_basis+gIIloc->n_b_basis;
		  // XYmtx is allocated as:
		  //   NdataXY is set to 0 here (won't be needed)
		  //   XY_nvects == (n_f_basis+n_b_basis)+NdataXY = dimenson of the vectors;
		  //   XY_nalloc == ndim_tot + 1    (one more vector is kept for the pole's energies)
		  gIIloc->NdataXY = 0;
		  gIIloc->allocate_XY(ndim_tot+1,ndim_tot);
		  dptr1 = gIIloc->XYmtx + gIIloc->NdataXY;
		  dptr2 = dptr1 + gIIloc->XY_nvects *  ndim_tot ;
		  gIIloc->Fill_LaddRPA_mtx(3, dptr1, gIIloc->XY_nvects, dptr2, false);
		  //
		  // Corrects for the weird sings used in the Ladd-DRPA calculation:
		  // (e_hh is already correct! so, sum sonly to n1<ndim_tot)
		  dptr1 = gIIloc->XYmtx + gIIloc->NdataXY + nf_dim;
		  n3 = gIIloc->XY_nvects;
		  for (n1=0; n1<ndim_tot; ++n1) {
			for (n2=0; n2<nb_dim; ++n2) dptr1[n2] = - dptr1[n2];
			dptr1 += n3;
		  }

		T2_Ld_dim[i_ld_njpt] = ndim_tot;
	  }
	}
	for(dT=MSp_CC->ch_ph_mn; dT<=MSp_CC->ch_ph_mx; ++dT) {
	  
	  J2 = 2*J;
	  i_pi_njpt = (J2+pi)*MULT_CH + dT - CH_OFFS_PI;

	  if (NULL != T2_Pi[i_pi_njpt]) delete T2_Pi[i_pi_njpt];
	  T2_Pi[i_pi_njpt] = new Pol_prop_t(MSp_CC,gsp_CC,Vpp_CC);
	  T2_Pi_dim[i_pi_njpt] = -1;
	  Piloc = T2_Pi[i_pi_njpt];
	  Piloc->Jf  = J;
	  Piloc->dT  = dT;
	  Piloc->pif = pi;
	  Piloc->Count_phRPA_basis(J ,dT, pi, &n4, &n5);
	  if (i_plot) cout << "J\\pi,dT = "<<J<<chip[pi]<<" , "<<dT<<"     i_pi_njpt="<<i_pi_njpt<<"     "
					   << "     nf = "<<n4<<"     nb = "<<n5<<endl<<flush;
	  if ((n4 > 0) && (n5 > 0)) {
		  Piloc->Build_phRPA_basis(J ,dT, pi);
		  nf_dim = Piloc->n_f_basis;
		  nb_dim = Piloc->n_b_basis;
		  ndim_tot = nf_dim + nb_dim;
		  // XYmtx is allocated as:
		  //   NdataXY is set to 0 here (won't be needed)
		  //   XY_nvects == (n_f_basis+n_b_basis)+NdataXY = dimenson of the vectors;
		  //   XY_nalloc == ndim_tot + 1    (one more vector is kept for the pole's energies)
		  Piloc->NdataXY = 0;
		  Piloc->allocate_XY(ndim_tot+1,ndim_tot);
		  dptr1 = Piloc->XYmtx + Piloc->NdataXY;
		  dptr2 = dptr1 + Piloc->XY_nvects *  ndim_tot ;
		  Piloc->Fill_phRPA_mtx(3, dptr1, Piloc->XY_nvects, dptr2, false);
		  //
		  // Corrects for the weird sings used in the ph-DRPA calculation:
		  // (e_hp has also the wrong sign! so, sum al the way to n1<(ndim_tot+1) )
		  dptr1 = Piloc->XYmtx + Piloc->NdataXY + nf_dim;
		  n3 = Piloc->XY_nvects;
		  for (n1=0; n1<(ndim_tot+1); ++n1) {
			for (n2=0; n2<nb_dim; ++n2) dptr1[n2] = - dptr1[n2];
			dptr1 += n3;
		  }
		  Correct_Jph_phase(Piloc);
		  
		T2_Pi_dim[i_pi_njpt] = ndim_tot;
	  }
	}
  }  
  
  return 0;}
  



static int nf;
static int nb, nvec_all;
static double *T2_ptr, *E2_ptr;
static double epp, ehh, eph, ehp;

static double Eladd, EPi;
static double x1, x2, x3, x4;

int CoupledCluster_t::Divide_deltaE(void ) {

  for (i1=0; i1<NT2_CHANNELS; ++i1) {
	// Ladd amplitudes:
	if (0 < T2_Ld_dim[i1]) {
	  gIIloc = T2_Ld[i1];
	  nf_dim = gIIloc->n_f_basis;
	  nb_dim = gIIloc->n_b_basis;
	  ndim_tot = nf_dim+nb_dim;
	  if (ndim_tot != T2_Ld_dim[i1]) {cerr << "Trouble with ndim_tot != T2_Ld_dim[i1] \n"; exit(-100);}
	  nvec_all = gIIloc->XY_nvects;
	  T2_ptr = gIIloc->XYmtx + gIIloc->NdataXY;
	  E2_ptr = T2_ptr + gIIloc->XY_nvects *  ndim_tot;
	  for (nf=0; nf<nf_dim; ++nf) {
		epp = E2_ptr[nf];
		for (nb=nf_dim; nb<ndim_tot; ++nb) {
		  ehh = E2_ptr[nb];
		  T2_ptr[nvec_all*nb + nf] /= (ehh - epp);
		}
	  }
	}
	
	// ph amplitudes:
	if (0 < T2_Pi_dim[i1]) {
	  Piloc = T2_Pi[i1];
	  nf_dim = Piloc->n_f_basis;
	  nb_dim = Piloc->n_b_basis;
	  ndim_tot = nf_dim+nb_dim;
	  if (ndim_tot != T2_Pi_dim[i1]) {cerr << "Trouble with ndim_tot != T2_Pi_dim[i1] \n"; exit(-100);}
	  nvec_all = Piloc->XY_nvects;
	  T2_ptr = Piloc->XYmtx + Piloc->NdataXY;
	  E2_ptr = T2_ptr + Piloc->XY_nvects *  ndim_tot;
	  for (nf=0; nf<nf_dim; ++nf) {
		eph = E2_ptr[nf];
		for (nb=nf_dim; nb<ndim_tot; ++nb) {
		  ehp = E2_ptr[nb];
		  T2_ptr[nvec_all*nb + nf] /= (ehp - eph);
		}
	  }
	}
  }
  
  return 0;}

int CoupledCluster_t::Calc_E2(double *Eladd , double *EPi/*=NULL*/) {

  //bool plot = false;
  if (NULL != Eladd) {
	//if (-10000.0 > (*Eladd)) plot = true;
	(*Eladd) = 0.0;
    for (i1=0; i1<NT2_CHANNELS; ++i1) {
	
	  // Ladd amplitudes:
	  x1 = 0.0;x3=0.0;
	  if (0 < T2_Ld_dim[i1]) {
		gIIloc = T2_Ld[i1];
		nf_dim = gIIloc->n_f_basis;
		nb_dim = gIIloc->n_b_basis;
		ndim_tot = nf_dim+nb_dim;
		if (ndim_tot != T2_Ld_dim[i1]) {cerr << "Trouble with ndim_tot != T2_Ld_dim[i1] \n"; exit(-100);}
		nvec_all = gIIloc->XY_nvects;
		T2_ptr = gIIloc->XYmtx + gIIloc->NdataXY;
		E2_ptr = T2_ptr + gIIloc->XY_nvects *  ndim_tot;
		x1 = 0.0;x3 = 0.0;
		for (nf=0; nf<nf_dim; ++nf) 
		  for (nb=nf_dim; nb<ndim_tot; ++nb) {
			x1 += T2_ptr[nvec_all*nf + nb] * T2_ptr[nvec_all*nb + nf];
			x3 += pow(T2_ptr[nvec_all*nf + nb], 2) / (E2_ptr[nb] - E2_ptr[nf]);
		  }
		(*Eladd) += x1 * double(2*gIIloc->Jf +1);
	  }
	  //cout << i1 << "   Ladd:   "<<x1<<" -- "<<x3<<endl;
	  //if (plot) cout << i1 << "   Ladd:   "<<x1 * double(2*gIIloc->Jf +1)<<endl;
	}
  }


  // ph amplitudes:
  if (NULL != EPi) {
	(*EPi) = 0.0;
	for (i1=0; i1<NT2_CHANNELS; ++i1) {
	  
	  x2 = 0.0; x4=0.0;
	  if (0 < T2_Pi_dim[i1]) {
		Piloc = T2_Pi[i1];
		nf_dim = Piloc->n_f_basis;
		nb_dim = Piloc->n_b_basis;
		ndim_tot = nf_dim+nb_dim;
		if (ndim_tot != T2_Pi_dim[i1]) {cerr << "Trouble with ndim_tot != T2_Pi_dim[i1] \n"; exit(-100);}
		nvec_all = Piloc->XY_nvects;
		T2_ptr = Piloc->XYmtx + Piloc->NdataXY;
		E2_ptr = T2_ptr + Piloc->XY_nvects *  ndim_tot;
		x2 = 0.0; x4=0.0;
		for (nf=0; nf<nf_dim; ++nf) 
		  for (nb=nf_dim; nb<ndim_tot; ++nb) {
			x2 += T2_ptr[nvec_all*nf + nb] * T2_ptr[nvec_all*nb + nf];
			x4 += pow(T2_ptr[nvec_all*nf + nb], 2) / (E2_ptr[nb] - E2_ptr[nf]);
		  }
		(*EPi) += x2 * double(2*Piloc->Jf +1);
	  }
	  //
	  //cout << i1 << "   Pi:   "<<x2<<" -- "<<x4<<endl;
	}
	(*EPi) /= 4.0;
  }

  //cout <<  "\n\n Total,   Ladd:   "<<*Eladd<< "     Pi:    "<<(*EPi)<<endl;

  return 0;}


static int  na,  ish_a, j2_a, ip_a, ch_a, nJab;
static int  nnb, ish_b, j2_b, ip_b, ch_b;
static int  kg,  ish_g, j2_g, ip_g, ch_g, nJgd;
static int  kd,  ish_d, j2_d, ip_d, ch_d;

static int pi_Pi, dT_Pi, pi_Ld, dT_Ld;
static int J2_Ld, J2_Pi, J2_min, J2_max;

static int nf2, nb2, nvec_all2, iph_f, iph_b;

static double t_pp, t_ph;
static double xerr,xerr_mx,xerr_mx2;

static int i2;

static int ierr;

int CoupledCluster_t::Pandya_Add_Pi_to_Ld(void ) {
  
  xerr_mx2 = 0.0;

  for (i_ld_njpt=0; i_ld_njpt<NT2_CHANNELS; ++i_ld_njpt) {

	xerr_mx = 0.0;
	// Ladd amplitudes:
	if (0 >= T2_Ld_dim[i_ld_njpt]) continue;

	gIIloc = T2_Ld[i_ld_njpt];
	nf_dim = gIIloc->n_f_basis;
	nb_dim = gIIloc->n_b_basis;
	ndim_tot = nf_dim+nb_dim;
	if (ndim_tot != T2_Ld_dim[i_ld_njpt]) {cerr << "Trouble with ndim_tot != T2_Ld_dim[i_ld_njpt] \n"; exit(-100);}
	nvec_all = gIIloc->XY_nvects;
	T2_ptr = gIIloc->XYmtx + gIIloc->NdataXY + nvec_all*nf_dim;

	ierr = 0;
	J2_Ld = 2*gIIloc->Jf;
	for (nf=0; nf<nf_dim; ++nf) 
	for (nb=0; nb<nb_dim; ++nb) {

	  na = gIIloc->n1_pp[nf]; ish_a = gIIloc->gsp->n_clj_p[na]; j2_a = MSp_CC->MSp_2j[ish_a]; ip_a = MSp_CC->MSp_ip[ish_a]; ch_a = MSp_CC->MSp_ch[ish_a];
	  nnb= gIIloc->n2_pp[nf]; ish_b = gIIloc->gsp->n_clj_p[nnb];j2_b = MSp_CC->MSp_2j[ish_b]; ip_b = MSp_CC->MSp_ip[ish_b]; ch_b = MSp_CC->MSp_ch[ish_b];
	  kg = gIIloc->k1_hh[nb]; ish_g = gIIloc->gsp->n_clj_h[kg]; j2_g = MSp_CC->MSp_2j[ish_g]; ip_g = MSp_CC->MSp_ip[ish_g]; ch_g = MSp_CC->MSp_ch[ish_g];
	  kd = gIIloc->k2_hh[nb]; ish_d = gIIloc->gsp->n_clj_h[kd]; j2_d = MSp_CC->MSp_2j[ish_d]; ip_d = MSp_CC->MSp_ip[ish_d]; ch_d = MSp_CC->MSp_ch[ish_d];


	  pi_Pi = (ip_a + ip_g)%2;	if ((ip_b + ip_d)%2 != pi_Pi) {cerr << "Trouble with pi_Pi...   STOP!\n"; exit(-100);}
	  dT_Pi = (ch_a - ch_g);	if ((ch_d - ch_b)   != dT_Pi) {cerr << "Trouble with dT_Pi...   STOP!\n"; exit(-100);}
	  //
	  i1 = abs(j2_a - j2_g); i2 = abs(j2_b - j2_d); J2_min = (i1 > i2) ? i1 : i2;
	  i1 =     j2_a + j2_g;  i2 =     j2_b + j2_d ; J2_max = (i1 < i2) ? i1 : i2;
	  //
	  t_pp = 0.0;
	  for(J2_Pi=J2_min; J2_Pi<=J2_max; J2_Pi+=2) {
		//
		//  direct-direct term:
		i_pi_njpt = (J2_Pi+pi_Pi)*MULT_CH + dT_Pi - CH_OFFS_PI;
		x1 = x2 = 0.0;
		if (0 < T2_Pi_dim[i_pi_njpt]) {
		  Piloc = T2_Pi[i_pi_njpt];
		  nvec_all2 = Piloc->XY_nvects;
		  dptr1 = Piloc->XYmtx + Piloc->NdataXY + nvec_all2*Piloc->n_f_basis;
		  //
		  ierr += Piloc->get_Pi_f_indx(&nf2, na,  kg);
		  ierr += Piloc->get_Pi_b_indx(&nb2, nnb, kd);  // OCCHIO all'ordine (n e poi k)!!
		  x1 = dptr1[nvec_all2*nb2 + nf2];
		}
		//
		//  exchange-exchange  term:
		i_pi_njpt = (J2_Pi+pi_Pi)*MULT_CH - dT_Pi - CH_OFFS_PI; // NB:  dT_Pi_ag = - dT_Pi_bd  (!!!)
		if (0 < T2_Pi_dim[i_pi_njpt]) {
		  Piloc = T2_Pi[i_pi_njpt];
		  nvec_all2 = Piloc->XY_nvects;
		  dptr1 = Piloc->XYmtx + Piloc->NdataXY + nvec_all2*Piloc->n_f_basis;
		  //
		  ierr += Piloc->get_Pi_b_indx(&nb2, na,  kg);
		  ierr += Piloc->get_Pi_f_indx(&nf2, nnb, kd);  // OCCHIO all'ordine (n e poi k)!!
		  x2 = dptr1[nvec_all2*nb2 + nf2];
		  if ( abs(j2_a - j2_g)%2 ) x2 = -x2;
		}
		
		x1 += x2;
		if ( abs((j2_b + j2_g + J2_Ld + J2_Pi)/2 )%2 ) x1 = -x1;
		t_pp += double(J2_Pi+1) * am.s6j(j2_a,j2_b,J2_Ld,j2_d,j2_g,J2_Pi) * x1;
	  }
	  
	  pi_Pi = (ip_a + ip_d)%2;	if ((ip_b + ip_g)%2 != pi_Pi) {cerr << "Trouble with pi_Pi...   STOP!\n"; exit(-100);}
	  dT_Pi = (ch_a - ch_d);	if ((ch_g - ch_b)   != dT_Pi) {cerr << "Trouble with dT_Pi...   STOP!\n"; exit(-100);}
	  //
	  i1 = abs(j2_a - j2_d); i2 = abs(j2_b - j2_g); J2_min = (i1 > i2) ? i1 : i2;
	  i1 =     j2_a + j2_d;  i2 =     j2_b + j2_g ; J2_max = (i1 < i2) ? i1 : i2;
	  x3 = 0.0;
	  for(J2_Pi=J2_min; J2_Pi<=J2_max; J2_Pi+=2) {
		//
		//  direct-exchange term:
		i_pi_njpt = (J2_Pi+pi_Pi)*MULT_CH + dT_Pi - CH_OFFS_PI;
		x1 = x2 = 0.0;
		if (0 < T2_Pi_dim[i_pi_njpt]) {
		  Piloc = T2_Pi[i_pi_njpt];
		  nvec_all2 = Piloc->XY_nvects;
		  dptr1 = Piloc->XYmtx + Piloc->NdataXY + nvec_all2*Piloc->n_f_basis;
		  //
		  ierr += Piloc->get_Pi_f_indx(&nf2, na,  kd);
		  ierr += Piloc->get_Pi_b_indx(&nb2, nnb, kg);  // OCCHIO all'ordine (n e poi k)!!
		  x1 = dptr1[nvec_all2*nb2 + nf2];
		}
		//
		//  exchange-direct  term:
		i_pi_njpt = (J2_Pi+pi_Pi)*MULT_CH - dT_Pi - CH_OFFS_PI; // NB:  dT_Pi_ad = - dT_Pi_bg  (!!!)
		if (0 < T2_Pi_dim[i_pi_njpt]) {
		  Piloc = T2_Pi[i_pi_njpt];
		  nvec_all2 = Piloc->XY_nvects;
		  dptr1 = Piloc->XYmtx + Piloc->NdataXY + nvec_all2*Piloc->n_f_basis;
		  //
		  ierr += Piloc->get_Pi_b_indx(&nb2, na,  kd);
		  ierr += Piloc->get_Pi_f_indx(&nf2, nnb, kg);  // OCCHIO all'ordine (n e poi k)!!
		  x2 = dptr1[nvec_all2*nb2 + nf2];
		  if ( abs(j2_a - j2_d)%2 ) x2 = -x2;
		}
		
		x1 += x2;
		if ( abs(1 + j2_d + (j2_b + j2_g + J2_Pi)/2)%2 ) x1 = -x1;
		x3 += double(J2_Pi+1) * am.s6j(j2_a,j2_b,J2_Ld,j2_g,j2_d,J2_Pi) * x1;
	  }

	  t_pp += x3;
	  if (na == nnb) t_pp /= SQR2;
	  if (kg == kd ) t_pp /= SQR2;


	  /* / TST
	  x1 = 4*T2_ptr[nvec_all*nb + nf];
	  xerr = x1-t_pp; if (1.e-10 < abs(t_pp)) xerr /= t_pp;
	  xerr = fabs(xerr);
	  xerr_mx = (xerr > xerr_mx) ? xerr : xerr_mx;
	  xerr_mx2 = (xerr > xerr_mx2) ? xerr : xerr_mx2;
	  if (xerr > 1.e-8) cout << nf << "  " << nb << "   -->  " << x1 << "  " << t_pp << endl;
	  //  */

	  T2_ptr[nvec_all*nb + nf] += t_pp;

	  }

	//cout << i_ld_njpt << "   Ladd (xerr_mx):   "<<xerr_mx<<" -- "<<xerr_mx2<<  "       ierr: "<<ierr <<endl;
	if (ierr) {cerr << "Trouble (CoupledCluster_t::Pandya_Copy_Ld_to_Pi):  ierr= " << ierr<< "   --> STOP.\n"<<flush; exit(100);}
  }

  //cout <<  "   Ladd (total Pandya error):   "<<xerr_mx2<<endl;

  return 0;}

int CoupledCluster_t::Pandya_Copy_Ld_to_Pi(void ) {
  
  xerr_mx2 = 0.0;
  
  for (i_pi_njpt=0; i_pi_njpt<NT2_CHANNELS; ++i_pi_njpt) {
	
	xerr_mx = 0.0;
	// Ladd amplitudes:
	if (0 >= T2_Pi_dim[i_pi_njpt]) continue;

	Piloc = T2_Pi[i_pi_njpt];
	nf_dim = Piloc->n_f_basis;
	nb_dim = Piloc->n_b_basis;
	ndim_tot = nf_dim+nb_dim;
	if (ndim_tot != T2_Pi_dim[i_pi_njpt]) {cerr << "Trouble with ndim_tot != T2_pi_dim[i_pi_njpt] \n"; exit(-100);}
	nvec_all = Piloc->XY_nvects;
	T2_ptr = Piloc->XYmtx + Piloc->NdataXY + nvec_all*nf_dim;
	
	ierr = 0;
	J2_Pi = 2*Piloc->Jf;
	for (nf=0; nf<nf_dim; ++nf) 
	  for (nb=0; nb<nb_dim; ++nb) {
		
		na = Piloc->n_ph[nf]; ish_a = Piloc->gsp->n_clj_p[na]; j2_a = MSp_CC->MSp_2j[ish_a]; ip_a = MSp_CC->MSp_ip[ish_a]; ch_a = MSp_CC->MSp_ch[ish_a];
		kg = Piloc->k_ph[nf]; ish_g = Piloc->gsp->n_clj_h[kg]; j2_g = MSp_CC->MSp_2j[ish_g]; ip_g = MSp_CC->MSp_ip[ish_g]; ch_g = MSp_CC->MSp_ch[ish_g];
		nnb= Piloc->n_hp[nb]; ish_b = Piloc->gsp->n_clj_p[nnb];j2_b = MSp_CC->MSp_2j[ish_b]; ip_b = MSp_CC->MSp_ip[ish_b]; ch_b = MSp_CC->MSp_ch[ish_b];
		kd = Piloc->k_hp[nb]; ish_d = Piloc->gsp->n_clj_h[kd]; j2_d = MSp_CC->MSp_2j[ish_d]; ip_d = MSp_CC->MSp_ip[ish_d]; ch_d = MSp_CC->MSp_ch[ish_d];
		
		pi_Ld = (ip_a + ip_b)%2;  if ((ip_g + ip_d)%2 != pi_Ld) {cerr << "Trouble with pi_Ld...   STOP!\n"; exit(-100);}
		dT_Ld = (ch_a + ch_b);	  if ((ch_g + ch_d)   != dT_Ld) {cerr << "Trouble with dT_Ld...   STOP!\n"; exit(-100);}
		
		i1 = abs(j2_a - j2_b); i2 = abs(j2_g - j2_d); J2_min = (i1 > i2) ? i1 : i2;
		i1 =     j2_a + j2_b;  i2 =     j2_g + j2_d ; J2_max = (i1 < i2) ? i1 : i2;
		t_ph = 0.0;
		for(J2_Ld=J2_min; J2_Ld<=J2_max; J2_Ld+=2) {
		  nJab = abs((j2_a + j2_b - J2_Ld + 2)/2)%2; if ((na == nnb) && nJab) continue;
		  nJgd = abs((j2_g + j2_d - J2_Ld + 2)/2)%2; if ((kg == kd ) && nJgd) continue;
		  i_ld_njpt = (J2_Ld+pi_Ld)*MULT_CH + dT_Ld - CH_OFFS_LD;
		  if (0 >= T2_Ld_dim[i_ld_njpt]) continue;
		  gIIloc = T2_Ld[i_ld_njpt];
		  nvec_all2 = gIIloc->XY_nvects;
		  dptr1 = gIIloc->XYmtx + gIIloc->NdataXY + nvec_all2*gIIloc->n_f_basis;
		  ierr += gIIloc->get_ld_f_indx(&nf2, na, nnb, &iph_f);
		  ierr += gIIloc->get_ld_b_indx(&nb2, kg, kd,  &iph_b);  // OCCHIO all'ordine (kd e poi kg)!!

		  x1 = dptr1[nvec_all2*nb2 + nf2];
		  if ( abs( (j2_b + j2_g + J2_Ld + J2_Pi)/2 + iph_f + iph_b)%2 ) x1 = - x1;
		  t_ph += double(J2_Ld+1) * am.s6j(j2_a,j2_b,J2_Ld,j2_d,j2_g,J2_Pi) * x1;
		}
		if (na == nnb) t_ph *= SQR2;
		if (kg == kd ) t_ph *= SQR2;


		/* / TST
		x1 = T2_ptr[nvec_all*nb + nf];
		xerr = x1-t_ph; if (1.e-10 < abs(t_ph)) xerr /= t_ph;
		xerr = fabs(xerr);
		xerr_mx = (xerr > xerr_mx) ? xerr : xerr_mx;
		xerr_mx2 = (xerr > xerr_mx2) ? xerr : xerr_mx2;
		//  */

		T2_ptr[nvec_all*nb + nf] = t_ph;

	  }

	//cout << i_pi_njpt << "   Pi (xerr_mx):   "<<xerr_mx<<" -- "<<xerr_mx2<<  "       ierr: "<<ierr <<endl;
	if (ierr) {cerr << "Trouble (CoupledCluster_t::Pandya_Copy_Ld_to_Pi):  ierr= " << ierr<< "   --> STOP.\n"<<flush; exit(100);}

  }
  
  //cout <<  "   Pi (total Pandya error):   "<<xerr_mx2<<endl;
  
  return 0;}

