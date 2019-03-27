/*
 *  BcDor-CCD.C
 *  
 *
 *  Routined for calculating the CCD amplitudes to be later udes
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
#include <iomanip>  //TST//
using namespace std;

//#include <cstdlib>
//#include <cstdio>
#include <cmath>

//#define DPI   3.14159265358979323846
#define SQR2  1.414213562373095048801688


//#include "BcDor-Ang_momenta.hh"
//#include "BcDor-HO_radial_me.hh"

#include "BcDor-Global_variables.hh"
#include "BcDor-CCD_classes.hh"


// LAPACK driver for diagonalizaton
extern "C" void dgemm_(char*, char*, int*, int*, int*, double*, 
                        double*, int*, double*, int*, double*, double*, int* );


void Mxt_multiply_wrt(double *A, int ndimA_r, int ndimA_c, int veclenA, double *B, int veclenB, int ndimBC, double *C, int veclenC) {
  register int	  imtx1,   imtx2,  imtx3;
  register double xmtx1;
  register double *dptrA, *dptrB, *dptrC;
  for (imtx2=0; imtx2<ndimA_c; ++imtx2) {
	dptrA = A + imtx2*veclenA;
	for (imtx1=0; imtx1<ndimA_r; ++imtx1) {
	  //
	  dptrC = C + imtx2*veclenC;
	  dptrB = B + imtx1;
	  xmtx1 = 0.0;
	  for (imtx3=0; imtx3<ndimBC; ++imtx3) {
		xmtx1 += (*dptrB)*(*dptrC);
		++dptrC;
		dptrB += veclenB;
	  }
	  (*dptrA) = xmtx1;
	  //
	  ++dptrA;
	}
  }
  return;}
void Mxt_multiply_add(double *A, int ndimA_r, int ndimA_c, int veclenA, double *B, int veclenB, int ndimBC, double *C, int veclenC, double xadd=1.0) {
  register int	  imtx1,   imtx2,  imtx3;
  register double xmtx1;
  register double *dptrA, *dptrB, *dptrC;
  for (imtx2=0; imtx2<ndimA_c; ++imtx2) {
	dptrA = A + imtx2*veclenA;
	for (imtx1=0; imtx1<ndimA_r; ++imtx1) {
	  //
	  dptrC = C + imtx2*veclenC;
	  dptrB = B + imtx1;
	  xmtx1 = 0.0;
	  for (imtx3=0; imtx3<ndimBC; ++imtx3) {
		xmtx1 += (*dptrB)*(*dptrC);
		++dptrC;
		dptrB += veclenB;
	  }
	  (*dptrA) += xmtx1*xadd;
	  //
	  ++dptrA;
	}
  }
  return;}

void Mxt_Copy_plus_transpose(double *A, int ndimA_r, int ndimA_c, int veclenA, double *B, int veclenB, double *C, int veclenC, double *Eff, double *Ebb, double aLM_in, double alpha_in) {
  register int	  imtx1,   imtx2;
  register double *dptrA, *dptrB, *dptrC, x1; //, x2, x3;
  //x2 = 1.0-aLM_in;
  //x3 = aLM_in*alpha_in;
  for (imtx2=0; imtx2<ndimA_c; ++imtx2) {
	dptrA = A + imtx2*veclenA;
	dptrB = B + imtx2*veclenB;
	dptrC = C + imtx2;
	for (imtx1=0; imtx1<ndimA_r; ++imtx1) {
	  //
    x1 = (1.0-aLM_in)*(*dptrA) + aLM_in*(alpha_in * (*dptrB) + (*dptrC))/(Ebb[imtx2]-Eff[imtx1]);
	  (*dptrA) = x1;
	  //
	  ++dptrA;
	  ++dptrB;
	  dptrC += veclenC;
	}
  }
  return;}
void Mxt_Copy(double *A, int ndimA_r, int ndimA_c, int veclenA, double *B, int veclenB, double *Eff, double *Ebb, double aLM_in, double alpha_in) {
  register int	  imtx1,   imtx2;
  register double *dptrA, *dptrB; //, x3;
  //x3 = aLM_in*alpha_in;
  for (imtx2=0; imtx2<ndimA_c; ++imtx2) {
	dptrA = A + imtx2*veclenA;
	dptrB = B + imtx2*veclenB;
	for (imtx1=0; imtx1<ndimA_r; ++imtx1) {
	  //
	  (*dptrA) = aLM_in*alpha_in * (*dptrB)/(Ebb[imtx2]-Eff[imtx1]);
	  //
	  ++dptrA;
	  ++dptrB;
	}
  }
  return;}


static int i1;

static int nf_ld_max, nb_ld_max, ntot_ld_max;
static int nf_pi_max, nb_pi_max, ntot_pi_max;

static int nf_dim,    nb_dim,    ntot_dim,  nveclen;

static double  *T2nnkk,  *Vnnnn,  *Vkkkk,  *Vkknn;
static double  *T2nkkn,  *Vnknk,  *Vknkn,  *Vknnk;

static gII_prop_t* gIIloc;
static Pol_prop_t* Piloc;

static int i_chnl;

static double  *Ld_wrk, *Pi_wrk, *Ld_wrk2, *Pi_wrk2;


//
//  Memory for storinf the Sn1n5 and Sk7k4 matrices... to be moved to
// the class' definition or some better place later on...
//
static double **Sn1n5_a=NULL, **Sn1n5_b=NULL, **Sk7k4_a=NULL, **Sk7k4_b=NULL;
static int    *Snn_dim=NULL,  *Skk_dim=NULL;
static int    NCHNLloc;

void CoupledCluster_t::Create_Snn_Skk(void ) {
  int n1;
  Destroy_Snn_Skk();
  NCHNLloc = MSp_CC->nsubsh;
  Snn_dim = new int[NCHNLloc];
  Skk_dim = new int[NCHNLloc];
  Sn1n5_a = new double*[NCHNLloc];
  Sn1n5_b = new double*[NCHNLloc];
  Sk7k4_a = new double*[NCHNLloc];
  Sk7k4_b = new double*[NCHNLloc];
  for (int ish=0; ish<NCHNLloc; ++ish) {
	n1 = gsp_CC->Sp_clj_np[ish];
	Snn_dim[ish] = n1;
	Sn1n5_a[ish] = new double[n1*n1+40];
	Sn1n5_b[ish] = new double[n1*n1+40];
	n1 = gsp_CC->Sp_clj_nh[ish];
	Skk_dim[ish] = n1;
	Sk7k4_a[ish] = new double[n1*n1+40];
	Sk7k4_b[ish] = new double[n1*n1+40];
  }
  return;}

void CoupledCluster_t::Destroy_Snn_Skk(void ) {
  for (int ish=0; ish<NCHNLloc; ++ish) {
	if (NULL != Sn1n5_a[ish]) delete [] Sn1n5_a[ish];
	if (NULL != Sn1n5_b[ish]) delete [] Sn1n5_b[ish];
	if (NULL != Sk7k4_a[ish]) delete [] Sk7k4_a[ish];
	if (NULL != Sk7k4_b[ish]) delete [] Sk7k4_b[ish];
  }
  if (NULL != Snn_dim) delete [] Snn_dim; Snn_dim = NULL;
  if (NULL != Skk_dim) delete [] Skk_dim; Skk_dim = NULL;
  if (NULL != Sn1n5_a) delete [] Sn1n5_a; Sn1n5_a = NULL;
  if (NULL != Sn1n5_b) delete [] Sn1n5_b; Sn1n5_b = NULL;
  if (NULL != Sk7k4_a) delete [] Sk7k4_a; Sk7k4_a = NULL;
  if (NULL != Sk7k4_b) delete [] Sk7k4_b; Sk7k4_b = NULL;
  return;}

void CoupledCluster_t::Clear_Snn_Skk(double **Snn, double **Skk ) {
  int n1, n2;
  for (int ish=0; ish<NCHNLloc; ++ish) {
	if (NULL != Snn) {
	  n1 = Snn_dim[ish];
	  for (n2=0; n2<n1*n1; ++n2) Snn[ish][n2] = 0.0;
	}
	if (NULL != Skk) {
	  n1 = Skk_dim[ish];
	  for (n2=0; n2<n1*n1; ++n2) Skk[ish][n2] = 0.0;
	}
  }
  return;}


int CoupledCluster_t::Solve_CCD(double aLM_in, int nStp_in) {


  //
  // seek for dimension of the required mtx.
  //
  nf_ld_max = 0;
  nb_ld_max = 0;
  nf_pi_max = 0;
  nb_pi_max = 0;
  for (i_chnl=0; i_chnl<NT2_CHANNELS; ++i_chnl) {
	// Ladd amplitudes:
	if (0 < T2_Ld_dim[i_chnl]) {
	  gIIloc = T2_Ld[i_chnl];
	  nf_dim = gIIloc->n_f_basis;
	  nb_dim = gIIloc->n_b_basis;
	  nf_ld_max = (nf_ld_max > nf_dim) ? nf_ld_max : nf_dim;
	  nb_ld_max = (nb_ld_max > nb_dim) ? nb_ld_max : nb_dim;
	}

	// ph amplitudes
	if (0 < T2_Pi_dim[i_chnl]) {
	  Piloc = T2_Pi[i_chnl];
	  nf_dim = Piloc->n_f_basis;
	  nb_dim = Piloc->n_b_basis;
	  nf_pi_max = (nf_pi_max > nf_dim) ? nf_pi_max : nf_dim;
	  nb_pi_max = (nb_pi_max > nb_dim) ? nb_pi_max : nb_dim;
	}
  }
  ntot_ld_max = nf_ld_max*nb_ld_max;
  ntot_pi_max = nf_pi_max*nb_pi_max;
  
  cout << " Max dimensions:\n     Ld:   "<<nf_ld_max<<"(fw)  x  "<<nb_ld_max<<"(bk)     ->  "<<ntot_ld_max<<"(total)\n";
  cout <<                   "     Pi:   "<<nf_pi_max<<"(fw)  x  "<<nb_pi_max<<"(bk)     ->  "<<ntot_pi_max<<"(total)\n";
  if (NULL == Ld_wrk  ) delete [] Ld_wrk;
  if (NULL == Pi_wrk  ) delete [] Pi_wrk;
  if (NULL == Ld_wrk2 ) delete [] Ld_wrk2;
  if (NULL == Pi_wrk2 ) delete [] Pi_wrk2;
  Ld_wrk  = new double [nf_ld_max*nb_ld_max];
  Pi_wrk  = new double [nf_pi_max*nb_pi_max];
  Ld_wrk2 = new double [nb_ld_max*nb_ld_max];
  Pi_wrk2 = new double [nb_pi_max*nb_pi_max];

  double **Snn_in, **Snn_out, **Skk_in, **Skk_out, **dummy;

  Create_Snn_Skk();
  Snn_in  = Sn1n5_a;
  Snn_out = Sn1n5_b;
  Skk_in  = Sk7k4_a;
  Skk_out = Sk7k4_b;



  double E2_ld, E2_ld_old, E2_ph, E2_ph_old;

  Calc_E2(&E2_ld, &E2_ph);
  cout << setfill(' ') << setprecision(12);
  cout <<  "\n Initial energies:      Ladd: " << setw(16) <<E2_ld<< "     Pi:    "<<E2_ph<<flush<<endl;

  Clear_Snn_Skk(Snn_in, Skk_in); //nothing added at first round

  int ierr = 0;
  int itr = 0;
  int i_alpha = 1;
  double alpha = 1.0;
  do {
	++itr;
	E2_ld_old = E2_ld;
	E2_ph_old = E2_ph;
	alpha = double(i_alpha)/double(nStp_in);
  if (i_alpha == nStp_in) alpha = 1.0;

	Pandya_Copy_Ld_to_Pi();

	Clear_Snn_Skk(Snn_out, Skk_out); // to store the new matrices

	Calc_New_T2(Ld_wrk, Pi_wrk, Ld_wrk2, Pi_wrk2, Snn_in, Snn_out, Skk_in, Skk_out, aLM_in, alpha);

	Pandya_Add_Pi_to_Ld();

	//Divide_deltaE();

	Calc_E2(&E2_ld, &E2_ph);
	cout <<  " Iteration number "<<itr<<" (alpha="<<setw(4)<<alpha<<"):    Ladd: " << setw(16) <<E2_ld<< "     Pi:    "<<E2_ph/alpha<<flush<<endl;

	dummy = Snn_out;
	Snn_out = Snn_in;
	Snn_in = dummy;
	dummy = Skk_out;
	Skk_out = Skk_in;
	Skk_in = Skk_out;

	if (( fabs(E2_ld-E2_ld_old) < 1.e-7) && (i_alpha < nStp_in)) {
    ++i_alpha;
    i_alpha = (i_alpha < nStp_in) ? i_alpha : nStp_in;
    cout << "   aLM_in="<<aLM_in<< "  alpha="<<alpha<<"  i_alpha="<<i_alpha<<"   nStp_in="<<nStp_in<<"\n";
    if (system( NULL )) system( "date" );
    //this->save_CCD_amplitudes();
  }

    if (fabs(E2_ld-E2_ld_old) > 1.e10) {ierr = 1; break;}

  } while(fabs(E2_ld-E2_ld_old) > 1.e-12);
  //
  Pandya_Copy_Ld_to_Pi();  // so that also Pi is correct...

  Calc_E2(&E2_ld, &E2_ph);

  if (ierr) {
    cerr <<  "\n\n CCD equations DID NOT CONVERGE!!!!\n\n";
    cout <<  "\n\n CCD equations DID NOT CONVERGE!!!!\n\n";
  } else {
  cout <<  " Converged result:       Ladd: " << setw(16) <<E2_ld<< "     Pi:    "<<E2_ph<<flush<<endl;
  }

  if (NULL == Ld_wrk  ) delete [] Ld_wrk;   Ld_wrk  = NULL;
  if (NULL == Pi_wrk  ) delete [] Pi_wrk;   Pi_wrk  = NULL;
  if (NULL == Ld_wrk2 ) delete [] Ld_wrk2;  Ld_wrk2 = NULL;
  if (NULL == Pi_wrk2 ) delete [] Pi_wrk2;  Pi_wrk2 = NULL;

  return ierr;}


static int k3, ish_3, i_k3, ikk34, ikk74, iph_74;
static int k4, ish_4, i_k4,        ikk37, iph_37;
static int k7, ish_7, i_k7, ndim_k7;

static int n1, ish_1, i_n1, inn12, inn16, iph_16;
static int n2, ish_2, i_n2,        inn62, iph_62;
static int n6, ish_6, i_n6, ndim_n6;

static double *Skk73, *Skk74;
static double *Snn26, *Snn16;

static double x1, x2, x3,  *dptr;

static int ikn_67, ikn_14;


void CoupledCluster_t::plot_nn(double **Snn_in, double **Snn_out) {
  for (int ish=0; ish<NCHNLloc; ++ish) {
	//if (5 != ish) continue;
	int ndim = Snn_dim[ish];
	cout << MSp_CC->MSp_name[ish] <<"   (dim "<< Snn_dim[ish] <<"):\n";
	for(int ir=0; ir<ndim; ++ir) {
	  for (int ic=0; ic<ndim; ++ic) cout <<Snn_in[ish][ir+ ic*ndim]<<"       ";
	  cout << "     ---     ";
	  for (int ic=0; ic<ndim; ++ic) cout <<Snn_out[ish][ir+ ic*ndim]<<"       ";
	  cout << endl;
	}
	cout << endl;
  }
  return;}

void CoupledCluster_t::plot_kk(double **Skk_in, double **Skk_out) {
  for (int ish=0; ish<NCHNLloc; ++ish) {
	//if (5 != ish) continue;
	int ndim = Skk_dim[ish];
	cout << MSp_CC->MSp_name[ish] <<"   (dim "<< Skk_dim[ish] <<"):\n";
	for(int ir=0; ir<ndim; ++ir) {
	  for (int ic=0; ic<ndim; ++ic) cout <<Skk_in[ish][ir+ ic*ndim]<<"       ";
	  cout << "     ---     ";
	  for (int ic=0; ic<ndim; ++ic) cout <<Skk_out[ish][ir+ ic*ndim]<<"       ";
	  cout << endl;
	}
	cout << endl;
  }
  return;}

int CoupledCluster_t::Calc_New_T2(double *Ld_wrk,  double *Pi_wrk,
								  double *Ld_wrk2, double *Pi_wrk2,
								  double **Snn_in, double **Snn_out,
								  double **Skk_in, double **Skk_out, double aLM_in, double alpha ) {

  double eph,ehp,epp,ehh,*E2_ptr;

  double dONE  = 1.0;
  //double dMONE  = -1.0; 
  double dZERO = 0.0; 
  double dHALF = 0.5; 
  //
  // Ladd amplitudes:
  //
  //TST//plot_nn(Snn_in, Snn_out);
  //TST//plot_kk(Skk_in, Skk_out);
  for (i_chnl=0; i_chnl<NT2_CHANNELS; ++i_chnl) {
	// Ladd amplitudes:
	if (0 < T2_Ld_dim[i_chnl]) {
	  gIIloc = T2_Ld[i_chnl];
	  nf_dim = gIIloc->n_f_basis;
	  nb_dim = gIIloc->n_b_basis;
	  ntot_dim = nf_dim + nb_dim;
	  if (ntot_dim != T2_Ld_dim[i_chnl]) {cerr << "Ohh, pleeeease...  not another bug!\n"; exit(100);}
	  nveclen = gIIloc->XY_nvects;
	  Vnnnn  = gIIloc->XYmtx + gIIloc->NdataXY;
	  T2nnkk = gIIloc->XYmtx + gIIloc->NdataXY  + nveclen *nf_dim;
	  Vkknn  = gIIloc->XYmtx + gIIloc->NdataXY                     + nf_dim;
	  Vkkkk  = gIIloc->XYmtx + gIIloc->NdataXY  + nveclen *nf_dim  + nf_dim;
	  E2_ptr = gIIloc->XYmtx + gIIloc->NdataXY  +                            nveclen * ntot_dim;

    
    for (i1=0; i1<nb_dim*nb_dim; ++i1)  Ld_wrk2[i1] = 0.0;  // this one not necessary...

	  for (i1=0; i1<nf_dim*nb_dim; ++i1)  Ld_wrk[i1] = 0.0; // this must be set to zero to start 
	  //
	  //  Qc(D3c):
	  for (ikk34=0; ikk34<nb_dim; ++ikk34) {
		k3 = gIIloc->k1_hh[ikk34]; ish_3 = gsp_CC->n_clj_h[k3]; i_k3 = k3 - (gsp_CC->Sp_clj_eh[ish_3] - gsp_CC->eh);
		k4 = gIIloc->k2_hh[ikk34]; ish_4 = gsp_CC->n_clj_h[k4]; i_k4 = k4 - (gsp_CC->Sp_clj_eh[ish_4] - gsp_CC->eh);

		for (k7 = 0; k7<gsp_CC->n_qh; ++k7) {
		  ish_7 = gsp_CC->n_clj_h[k7]; i_k7 = k7 - (gsp_CC->Sp_clj_eh[ish_7] - gsp_CC->eh);
		  ndim_k7 = gsp_CC->Sp_clj_nh[ish_7];
		  if (1 > ndim_k7) continue;
		  if (ndim_k7 != Skk_dim[ish_7]) {cerr << " AAAAAAAA---- \n"; exit(-100);}
		  Skk73 = Skk_in[ish_7] + i_k7 + i_k3*ndim_k7;
		  Skk74 = Skk_in[ish_7] + i_k7 + i_k4*ndim_k7;
		  if ((ish_7 == ish_3) && ( ! gIIloc->get_ld_b_indx(&ikk74, k7, k4, &iph_74) )) {
			/*/
			if ((i_k7 >= ndim_k7) || (i_k3 >= ndim_k7) || (i_k3 < 0) || (i_k7 < 0))
			  cout << "ish_7a:  "<<ish_7 <<"   "<<i_k3<<"   "<<i_k7<<"   "<<Skk_dim[ish_7]<<"  "<<ndim_k7<<endl;
			/*/
			x1 = (*Skk73);
			if (iph_74%2) x1 = -x1;
			if (k7 == k4) x1 *= SQR2;
			if (k3 == k4) x1 /= SQR2;
			for (inn12=0; inn12<nf_dim; ++inn12)
			  Ld_wrk[ikk34*nf_dim + inn12] += T2nnkk[ikk74*nveclen + inn12] * x1;
		  }
		  if ((ish_7 == ish_4) && ( ! gIIloc->get_ld_b_indx(&ikk37, k3, k7, &iph_37) )) {
			/*/
			if ((i_k7 >= ndim_k7) || (i_k4 >= ndim_k7) || (i_k4 < 0) || (i_k7 < 0))
			  cout << "ish_7b:  "<<ish_7 <<"   "<<i_k4<<"   "<<i_k7<<"   "<<Skk_dim[ish_7]<<"  "<<ndim_k7<<endl;
			/*/
			x1 = (*Skk74);
			if (iph_37%2) x1 = -x1;
			if (k3 == k7) x1 *= SQR2;
			if (k3 == k4) x1 /= SQR2;
			for (inn12=0; inn12<nf_dim; ++inn12) 
			  Ld_wrk[ikk34*nf_dim + inn12] += T2nnkk[ikk37*nveclen + inn12] * x1;
		  }
		}
	  }
	  //
	  //  Qd(D3d):
	  for (inn12=0; inn12<nf_dim; ++inn12) {
		n1 = gIIloc->n1_pp[inn12]; ish_1 = gsp_CC->n_clj_p[n1]; i_n1 = n1 - (gsp_CC->Sp_clj_ep[ish_1] - gsp_CC->ep);
		n2 = gIIloc->n2_pp[inn12]; ish_2 = gsp_CC->n_clj_p[n2]; i_n2 = n2 - (gsp_CC->Sp_clj_ep[ish_2] - gsp_CC->ep);

		for (n6=0; n6<gsp_CC->n_qp; ++n6) {
		  ish_6 = gsp_CC->n_clj_p[n6]; i_n6 = n6 - (gsp_CC->Sp_clj_ep[ish_6] - gsp_CC->ep);
		  ndim_n6 = gsp_CC->Sp_clj_np[ish_6];
		  if (1 > ndim_n6) continue; 
		  if (ndim_n6 != Snn_dim[ish_6]) {cerr << " BBBBBBBBB---- \n"; exit(-100);}
		  Snn16 = Snn_in[ish_6] + i_n1 + i_n6*ndim_n6;
		  Snn26 = Snn_in[ish_6] + i_n2 + i_n6*ndim_n6;
		  if ((ish_6 == ish_1) && (! gIIloc->get_ld_f_indx(&inn62, n6, n2, &iph_62) )) {
			/*/
			if ((i_n1 >= ndim_n6) || (i_n6 >= ndim_n6) || (i_n1 < 0) || (i_n6 < 0))
			  cout << "ish_6a:  "<<ish_6 <<"   "<<i_n1<<"   "<<i_n6<<"   "<<Snn_dim[ish_6]<<"  "<<ndim_n6<<endl;
			/*/
			x1 = (*Snn16);
			if (iph_62%2) x1 = -x1;
			if (n6 == n2) x1 *= SQR2;
			if (n1 == n2) x1 /= SQR2;
			for (ikk34=0; ikk34<nb_dim; ++ikk34) 
			  Ld_wrk[ikk34*nf_dim + inn12] += T2nnkk[ikk34*nveclen + inn62] * x1;
		  }
		  if ((ish_6 == ish_2) && (! gIIloc->get_ld_f_indx(&inn16, n1, n6, &iph_16) )) {
			/*/
			if ((i_n2 >= ndim_n6) || (i_n6 >= ndim_n6) || (i_n2 < 0) || (i_n6 < 0))
			  cout << "ish_6b:  "<<ish_6 <<"   "<<i_n2<<"   "<<i_n6<<"   "<<Snn_dim[ish_6]<<"  "<<ndim_n6<<endl;
			/*/
			x1 = (*Snn26);
			if (iph_16%2) x1 = -x1;
			if (n1 == n6) x1 *= SQR2;
			if (n1 == n2) x1 /= SQR2;
			for (ikk34=0; ikk34<nb_dim; ++ikk34)
			  Ld_wrk[ikk34*nf_dim + inn12] += T2nnkk[ikk34*nveclen + inn16] * x1;
		  }  
		}
	  }

	  //
	  // L2a(D2c):
	  //Mxt_multiply_wrt(Ld_wrk,  nf_dim, nb_dim, nf_dim, Vnnnn,  nveclen, nf_dim, T2nnkk,  nveclen);
	  //NO-TSTX//dgemm_("No transp","No transp", &nf_dim, &nb_dim, &nf_dim, &dONE, Vnnnn,  &nveclen, T2nnkk,  &nveclen, &dZERO, Ld_wrk, &nf_dim);
	  dgemm_("No transp","No transp", &nf_dim, &nb_dim, &nf_dim, &dONE, Vnnnn,  &nveclen, T2nnkk,  &nveclen, &dONE, Ld_wrk, &nf_dim);
	  //
	  // L2b(D2d):
	  //Mxt_multiply_add(Ld_wrk,  nf_dim, nb_dim, nf_dim, T2nnkk, nveclen, nb_dim, Vkkkk,   nveclen);
	  dgemm_("No transp","No transp", &nf_dim, &nb_dim, &nb_dim, &dONE, T2nnkk,  &nveclen, Vkkkk,  &nveclen, &dONE, Ld_wrk, &nf_dim);
	  //
	  // Qa(D3a):
	  //Mxt_multiply_wrt(Ld_wrk2, nb_dim, nb_dim, nb_dim, Vkknn,  nveclen, nf_dim, T2nnkk,  nveclen);
	  dgemm_("No transp","No transp", &nb_dim, &nb_dim, &nf_dim, &dONE, Vkknn,  &nveclen, T2nnkk,  &nveclen, &dZERO, Ld_wrk2, &nb_dim);
	  //Mxt_multiply_add(Ld_wrk,  nf_dim, nb_dim, nf_dim, T2nnkk, nveclen, nb_dim, Ld_wrk2, nb_dim );
	  dgemm_("No transp","No transp", &nf_dim, &nb_dim, &nb_dim, &dONE, T2nnkk,  &nveclen, Ld_wrk2,  &nb_dim, &dONE, Ld_wrk, &nf_dim);	  

	  //
	  //  Finally replace into the new amplitude
	  //
	  Mxt_Copy_plus_transpose(T2nnkk, nf_dim, nb_dim, nveclen, Ld_wrk, nf_dim, Vkknn, nveclen, E2_ptr, E2_ptr+nf_dim, aLM_in, alpha);
	  
	} // end of: if (T2_Ld_dim>0)...
	
	// ph amplitudes:
	if (0 < T2_Pi_dim[i_chnl]) {
	  Piloc = T2_Pi[i_chnl];
	  nf_dim = Piloc->n_f_basis;
	  nb_dim = Piloc->n_b_basis;
	  ntot_dim = nf_dim + nb_dim;
	  if (ntot_dim != T2_Pi_dim[i_chnl]) {cerr << "Ohh, pleeeease...  not another bug!\n"; exit(100);}
	  nveclen = Piloc->XY_nvects;
	  Vnknk  = Piloc->XYmtx + Piloc->NdataXY;
	  T2nkkn = Piloc->XYmtx + Piloc->NdataXY  + nveclen *nf_dim;
	  Vknnk  = Piloc->XYmtx + Piloc->NdataXY                     + nf_dim;
	  Vknkn  = Piloc->XYmtx + Piloc->NdataXY  + nveclen *nf_dim  + nf_dim;
    E2_ptr = Piloc->XYmtx + Piloc->NdataXY  +                            nveclen * ntot_dim;
	  for (i1=0; i1<nf_dim*nb_dim; ++i1)  Pi_wrk[i1] = 0.0;
	  for (i1=0; i1<nb_dim*nb_dim; ++i1)  Pi_wrk2[i1] = 0.0;
	  //
	  // L2c(D2b):
	  //Mxt_multiply_wrt(Pi_wrk,  nf_dim, nb_dim, nf_dim, T2nkkn, nveclen, nb_dim, Vknkn,   nveclen);
	  dgemm_("No transp","No transp", &nf_dim, &nb_dim, &nb_dim, &dONE, T2nkkn,  &nveclen, Vknkn,  &nveclen, &dZERO, Pi_wrk, &nf_dim);
	  //
	  //
	  // Qb(D3b):
	  //Mxt_multiply_wrt(Pi_wrk2, nb_dim, nb_dim, nb_dim, Vknnk,  nveclen, nf_dim, T2nkkn,  nveclen);
	  dgemm_("No transp","No transp", &nb_dim, &nb_dim, &nf_dim, &dONE, Vknnk,  &nveclen, T2nkkn,  &nveclen, &dZERO, Pi_wrk2, &nb_dim);
	  //Mxt_multiply_add(Pi_wrk,  nf_dim, nb_dim, nf_dim, T2nkkn, nveclen, nb_dim, Pi_wrk2, nb_dim, 0.5);
	  dgemm_("No transp","No transp", &nf_dim, &nb_dim, &nb_dim, &dHALF, T2nkkn,  &nveclen, Pi_wrk2,  &nb_dim, &dONE, Pi_wrk, &nf_dim);


	  //
	  //  Finally replace into the new amplitude
	  //
	  // This is good ONLY in the case one wants to summ ONLY the Ph rings
	  //Mxt_Copy_plus_transpose(T2nkkn, nf_dim, nb_dim, nveclen, Pi_wrk, nf_dim, Vknnk, nveclen, 1.0, 1.0);
	  //
	  // the correct one, for full CCD/CCSD is:
	  Mxt_Copy(T2nkkn, nf_dim, nb_dim, nveclen, Pi_wrk, nf_dim, E2_ptr, E2_ptr+nf_dim, aLM_in, alpha);
	  

	  //
	  //
	  // Use Pi_wrk2  to generate the hole and particle self-energy insetrions
	  //
	  x1 = double(2*Piloc->Jf + 1) / 2.0;
	  dptr = Pi_wrk2;
	  for (ikn_14=0; ikn_14<nb_dim; ++ikn_14) {
		n1 = Piloc->n_hp[ikn_14]; ish_1 = gsp_CC->n_clj_p[n1]; i_n1 = n1 - (gsp_CC->Sp_clj_ep[ish_1] - gsp_CC->ep);
		k4 = Piloc->k_hp[ikn_14]; ish_4 = gsp_CC->n_clj_h[k4]; i_k4 = k4 - (gsp_CC->Sp_clj_eh[ish_4] - gsp_CC->eh);
		x2 = x1/double(MSp_CC->MSp_2j[ish_4] + 1);
		x3 = x1/double(MSp_CC->MSp_2j[ish_1] + 1);
		for (ikn_67=0; ikn_67<nb_dim; ++ikn_67) {
		  n6 = Piloc->n_hp[ikn_67]; ish_6 = gsp_CC->n_clj_p[n6]; i_n6 = n6 - (gsp_CC->Sp_clj_ep[ish_6] - gsp_CC->ep);
		  k7 = Piloc->k_hp[ikn_67]; ish_7 = gsp_CC->n_clj_h[k7]; i_k7 = k7 - (gsp_CC->Sp_clj_eh[ish_7] - gsp_CC->eh);
		  //
		  if ((n1 == n6) && (ish_7 == ish_4)) {
			//if ((i_k4 >= Skk_dim[ish_4]) || (i_k7 >= Skk_dim[ish_4]) || (i_k4 < 0) || (i_k7 < 0))
			//  cout << "ish_4:  "<<ish_4 <<"   "<<i_k4<<"   "<<i_k7<<"   "<<Skk_dim[ish_4]<<endl;
			Skk_out[ish_4][i_k7 + i_k4*Skk_dim[ish_4]] -= (*dptr) * x2;
		  }
		  //
		  if ((k7 == k4) && (ish_1 == ish_6)) {
			//if ((i_n1 >= Snn_dim[ish_1]) || (i_n6 >= Snn_dim[ish_1]) || (i_n1 < 0) || (i_n6 < 0))
			//  cout << "ish_1:  "<<ish_1 <<"   "<<i_n1<<"   "<<i_n6<<"   "<<Snn_dim[ish_1]<<endl;
			Snn_out[ish_1][i_n1 + i_n6*Snn_dim[ish_1]] -= (*dptr) * x3;
		  }

		++dptr;
		}
	  }

	} // end of: if (T2_Pi_dim>0)...

  } // end i_chnl
  //TST//plot_nn(Snn_in, Snn_out);
  //TST//plot_kk(Skk_in, Skk_out);

  return 0;}  



