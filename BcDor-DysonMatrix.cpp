//
//
//   Takes the self-energy in a given partial wave,
//  construct the Matrix representation of Dyson eq.
//  and then diagonalizes it.
//
//   The result is stored directly in the output s.p.
//  propagator. The reoutine also returns the partial
//  contribiution to the Koltun sum rule due to the 
//  partial wave being calculated...
//
//
//   C. Barbieri  --  RIKEN, May 2010.
//

#include <iostream>
#include <cstring>

using namespace std;

#include <cstdlib>
//#include <cstdio>
#include <cmath>

#include "BcDor-Run_vars.hh"
#include "BcDor-Global_variables.hh"


static  ModSpace_t    *DMdSp=NULL;
//static  VppInt_t      *DVpp=NULL;
static  SpProp_t      *Dgin=NULL;
static  Wide_self_energy_t *DSlfE=NULL;

// LAPACK driver for diagonalizaton
extern "C" void dsyevd_(char*, char*, int*, double*, int*, 
                        double*, double*, int*, int*, int*,  int* );




static int nf_lost, nb_lost;
static int i,j;//,l;
static int nr,nc,ne;
static int ndim,ndim_fw,ndim_bk;
static int ndim_se_tot;
// 
static int NDIMHF, NFWP, NBKP;          // HFE, fw poles, and bk poles dims.
static int NSOLDIM;                     // dims of sols arrays.
// 
static double *ptr_ms, *ptr_sf_e, *ptr_sf_o;
static double *ptr1, *ptr2, *ptr3, *ptr4, *ptr5;
static double x1, x2, x3, x4;

static double EFermi;


static bool   HFB_pair;
static double mu_ext;






int DysonMatrix(Wide_self_energy_t *SlfE, SpProp_t *Dout,
              double* xe_t, double* xe_tVV, double* xA) {


  DSlfE = SlfE;
  DMdSp = SlfE->MdSp;
//  DVpp  = SlfE->Vpp;
  Dgin  = SlfE->gin;

  NDIMHF =  DSlfE->NTOT_ORB;
  NFWP   =  DSlfE->N_PLS_fw;  // # of forward  poles
  NBKP   =  DSlfE->N_PLS_bk;  // # of backward poles

  NSOLDIM = NDIMHF+2;  // dim. of X or Y vector + s.p. energy + spect. fact.


  // Get the Fermi energy
  EFermi = Seek_EFermi(Dgin , SlfE->I_SH);


 

  cout << "Dimension of E-HF matrix = " << NDIMHF  << endl;
  cout << "Number of FW poles       = " << NFWP    << endl;
  cout << "Number of BK poles       = " << NBKP    << endl;


  //
  // now build the Dyson matrix...
  //


  //  Se la parte dinamica della self-energia e` vuota (o non e`
  // stata allocata, NFWP e NBKP potrebbero essere < 0.
  //  Qui ndim_fw e ndim_bk sono forzate ad essere >=0 e vengono
  // usate come dimensioni della parte 2p1h/2h1p della matrice.
  // se ndim_fw==ndim_bk==0 vuol dire che non ci sono diagrammi
  // oltre al primo ordine e la matrice si trasforma nella
  // hamiltoniana di Hartree-Fock.
  ndim_fw = (0 > NFWP) ? 0 : NFWP;
  ndim_bk = (0 > NBKP) ? 0 : NBKP;

  ndim = NDIMHF + ndim_fw + ndim_bk;
  ndim = (ndim >   0 ) ? ndim :  0 ;  // just for extra "safety"...



  //
  //  Terms for adding an external pairing potential with given \mu
  // 
  HFB_pair = (3 == i_Ext_U1) || ((4 == i_Ext_U1) && EFermi_Is_Set(SlfE->I_SH)); // whether to add an ext. pairing term
  if (HFB_pair){
	ndim += NDIMHF;
	ndim = (ndim >   0 ) ? ndim :  0 ;  // just for extra "safety"...
	mu_ext = Ext_U1_Uch[DMdSp->MSp_ch[SlfE->I_SH]];
	EFermi = mu_ext;
  }

// The stuff needed by the LAPACK eigenvalue pakage:
// For 'DSYEV':
  int INFO;
  int LDA = ndim+1;
//int LWORK  = 3*LDA+2;          // For 'DSYEV'  only...
  int LWORK  = 2+(6+2*ndim)*(ndim);   // For 'DSYEVD' only...
  int LIWORK = 3+5*ndim;          // For 'DSYEVD' only...

  cout << " LDA    = " << LDA    << endl;
  cout << " LWORK  = " << LWORK  << endl;
  cout << " LIWORK = " << LIWORK << endl;

  double *MTX, *W, *WORK;
  int    *IWORK;
  W = new double[LDA];
  WORK = new double[LWORK];
  IWORK = new int[LIWORK];   

  int LWopt=0, LIWopt=0;

  MTX = new double[LDA*LDA];
  for(ptr1= MTX; ptr1<MTX+(LDA*LDA); ++ptr1) (*ptr1)=0.0;



//
//  the structure of the Dyson Matrix is as follow:
//
//         /     |                 |                \
//         | HF  |   M2p1h^dag     |    M2h1p^dag   |
//         |     |                 |                |
//         |-----+-----------------+----------------|
//         |     |                 |                |
//         | M   |                 |                |
//         |2p1h |     diag{       |                |
//         |     |     e_2p1h      |                |
//   MTX = {     |       }         |      0         |
//         |     |                 |                |
//         |     |                 |                |
//         |-----+-----------------+----------------|
//         |     |                 |                |
//         |     |                 |                |
//         |  M  |                 |     diag{      |
//         |2h1p |        0        |    e_2h1p      |
//         |     |                 |       }        |
//         |     |                 |                |
//         \     |                 |                /
//
//
//

// Construct the E-HF part of the self energy (HF, BHF, HFB, or wathever...)
  ptr1=MTX;
  for(nr=0 ; nr<NDIMHF; ++nr) {
   for(nc=0 ; nc<NDIMHF; ++nc) {
    ptr1[nc] = DSlfE->Sig_1b[nr*NDIMHF+nc] + DSlfE->Sig_HF[nr*NDIMHF+nc];
    }
   ptr1 += LDA;
   }
	  

// Construct the forward part of the matrix
  ndim_se_tot = DSlfE->NTOT_DIM;
//ndim_se_den = DSlfE->NTOT_DEN;
//ndim_se_orb = DSlfE->NTOT_ORB;
  if (NFWP > 0) {
    ptr_sf_e = DSlfE->Sig_dyn_fw;
    ptr_sf_o = ptr_sf_e + DSlfE->NTOT_DEN;
    ptr1 = MTX + NDIMHF*(LDA+1);
    ptr2 = MTX + NDIMHF*LDA;
    ptr3 = MTX + NDIMHF;

    for(i=0; i<ndim_fw;++i) {
      (*ptr1) = (*ptr_sf_e);
      for(nr=0 ; nr<NDIMHF; ++nr) {
        ptr2[nr]     = ptr_sf_o[nr];
        ptr3[nr*LDA] = ptr_sf_o[nr];
        }
     ptr_sf_e += ndim_se_tot;
     ptr_sf_o += ndim_se_tot;
     ptr1 +=(LDA+1);
     ptr2 += LDA;
     ++ptr3;
     }
    //
    // 
    if (DSlfE->n_piv_fw > 0) {
	  ptr_sf_e = DSlfE->Sig_dyn_fw + ndim_se_tot + 1;
	  ptr1 = MTX + NDIMHF*(LDA+1) + 1; // upper sub-diaginal
	  ptr2 = ptr1 + (LDA-1);           // lower sub-iagonal
	  nc = 1;  // will label the g_i branches
	  //
	  for (i=1; i<ndim_fw; ++i) {
		  (*ptr1) = (*ptr_sf_e);
		  (*ptr2) = (*ptr_sf_e);
		if (i == DSlfE->i_piv_fw[nc+1]) { // a new branch starts from here
		  ptr3 = ptr1;
		  ptr4 = ptr2;
		  ptr5 = ptr_sf_e + nc;
		  for(nr=i ; nr<ndim_fw; ++nr) {
			(*ptr3) = (*ptr5);
			(*ptr4) = (*ptr5);
			++ptr3;
			ptr4 += LDA;
			ptr5 += ndim_se_tot;
		  }
		  ++nc;
		}
		ptr_sf_e += ndim_se_tot;
		ptr1 +=(LDA+1);
		ptr2 +=(LDA+1);
	  }
	  if (nc > DSlfE->n_piv_fw) {
		cout << " Something went wrong while building the fw Dys. matrix.\n--> exit.\n";
		exit(-300);
	  }
    }
  }


// Construct the backward part of the matrix
  if (NBKP > 0) {
    ptr_sf_e = DSlfE->Sig_dyn_bk;
    ptr_sf_o = ptr_sf_e + DSlfE->NTOT_DEN;;
    ptr1 = MTX + (NDIMHF+ndim_fw)*(LDA+1);
    ptr2 = MTX + (NDIMHF+ndim_fw)*LDA ;
    ptr3 = MTX +  NDIMHF+ndim_fw;

    for(i=0; i<ndim_bk;++i) {
      (*ptr1) = (*ptr_sf_e);  
      for(nr=0 ; nr<NDIMHF; ++nr) {
        ptr2[nr]     = ptr_sf_o[nr];
        ptr3[nr*LDA] = ptr_sf_o[nr];
        }
     ptr_sf_e += ndim_se_tot;
     ptr_sf_o += ndim_se_tot;
     ptr1 +=(LDA+1);
     ptr2 += LDA;
     ++ptr3;
     }
    //
    // 
    if (DSlfE->n_piv_bk > 0) {
	  ptr_sf_e = DSlfE->Sig_dyn_bk + ndim_se_tot + 1;
	  ptr1 = MTX + (NDIMHF+ndim_fw)*(LDA+1) + 1; // upper sub-diaginal
	  ptr2 = ptr1 + (LDA-1);           // lower sub-iagonal
	  nc = 1;  // will label the g_i branches
	  //
	  for (i=1; i<ndim_bk; ++i) {
		(*ptr1) = (*ptr_sf_e);
		(*ptr2) = (*ptr_sf_e);
		if (i == DSlfE->i_piv_bk[nc+1]) { // a new branch starts from here
		  ptr3 = ptr1;
		  ptr4 = ptr2;
		  ptr5 = ptr_sf_e + nc;
		  for(nr=i ; nr<ndim_bk; ++nr) {
			(*ptr3) = (*ptr5);
			(*ptr4) = (*ptr5);
			++ptr3;
			ptr4 += LDA;
			ptr5 += ndim_se_tot;
		  }
		  ++nc;
		}
		ptr_sf_e += ndim_se_tot;
		ptr1 +=(LDA+1);
		ptr2 +=(LDA+1);
	  }
	  if (nc > DSlfE->n_piv_bk) {
		cout << " Something went wrong while building the bk Dys. matrix.\n--> exit.\n";
		exit(-300);
	  }
	}
  }


  
  if (HFB_pair) {
	Add_AnomHF(MTX, LDA, ndim-NDIMHF, NDIMHF, DMdSp, SlfE->I_SH, mu_ext);
	//ndim -= NDIMHF;
	//ndim = (ndim >   0 ) ? ndim :  0 ;  // just for extra "safety"...
  }
  

//// TEST STUFF --------------
//  if ( HFB_pair ) {
//	double x1 = Ext_U1_Uch[DMdSp->MSp_ch[SlfE->I_SH]];
//	for(nr=0 ; nr<NDIMHF; ++nr) MTX[nr*(LDA+1)] -= x1;
//	cout << "\n SHIFT by  "<<x1<<"\n\n";
//  }
////-------------------

// LAPACK library (w/ DSYEVD):
       INFO = 0;
/*408*/	
       dsyevd_("Vectors","Upper",&ndim,MTX,&LDA,W,WORK, &LWORK,IWORK,&LIWORK,&INFO);
       if (0 != INFO) cout<< "\nDyson (DSYEV), Wrong value of IFAIL: INFO= "<<INFO<<endl;
       if (0 == INFO) {
         ne = int(WORK[0] + 0.01);
         LWopt = (ne > LWopt) ? ne : LWopt;
         LIWopt = (IWORK[0] > LIWopt) ? IWORK[0] : LIWopt;
/*410*/} else {
        cout << " 'DEEGV' gave INFO = "   << INFO << endl;
        if (INFO < 0) cout << "\n The " <<-INFO <<"-th aggument in line 408 was illegal:\n";
        cout << "\n\n Program has been aborted since IERR != 0."
             <<   "\nAborted at the line 410, Dyson.f!"
             <<   "\n   --> stop.\n\n";
        exit(1);
       }



  double *solf, *solb, *solDys, *dptr;
  int    nsolf, nsolb;

  bool   qparticle;


  solDys = new double[ndim*NSOLDIM];
  solf = solDys;
  solb = solDys + ndim*NSOLDIM;
  nsolf = 0;
  nsolb = 0;
  

  for(ne=0; ne<ndim; ++ne) {
    ptr1=MTX+ne*LDA;
    ptr2=ptr1;
	//// TEST STUFF --------------
	//if (HFB_pair) W[ne] += Ext_U1_Uch[DMdSp->MSp_ch[SlfE->I_SH]];
	////-------------------	
    x1=0.0; x2=0.0; x3=0.0; x4=0.0;
    for (i=0; i<ndim; ++i) {
      x1 += (*ptr2)*(*ptr2);
      if (i < NDIMHF)              {x2 += (*ptr2)*(*ptr2);}
      else if (i < NDIMHF+ndim_fw) {x3 += (*ptr2)*(*ptr2);}
      else                         {x4 += (*ptr2)*(*ptr2);}
      ++ptr2;
      }
    qparticle = false; if (W[ne] > EFermi) qparticle = true;

    if (qparticle) {
      solf[0] = W[ne];
      solf[1] = x2/x1;
      for(i=0; i<NDIMHF; ++i) solf[i+2] = ptr1[i];
      solf += NSOLDIM;
      ++nsolf;
    } else {
      solb -= NSOLDIM;
      ++nsolb;
      solb[0] = W[ne];
      solb[1] = x2/x1;
      for(i=0; i<NDIMHF; ++i) solb[i+2] = ptr1[i];
    }
  } // en loop `ne'
  
  solf = solDys;
  //solb = solDys + (ndim-nsolb)*NSOLDIM;  // it sould be already set like this....


  delete [] W;
  delete [] WORK;
  delete [] IWORK;   
  delete [] MTX;


  Sort_d_double2dim(nsolf, NSOLDIM, 0, solf);
  Sort_d_double2dim(nsolb, NSOLDIM, 0, solb);

  
  if (NULL!=Dout) nf_lost=Dout->add_qp(SlfE->I_SH, nsolf, NSOLDIM, solf);
  if (NULL!=Dout) nb_lost=Dout->add_qh(SlfE->I_SH, nsolb, NSOLDIM, solb);
  //
  // se non vi e` abbastanza memoria allocata per il propagatore (Dout), le
  // routines 'add_qp' e 'add_qh' ritornano il numer disoluzioni perse...
  if (nf_lost) cout << "# fw poles not saved="<< nf_lost << endl;
  if (nb_lost) cout << "# bk poles not saved="<< nb_lost << endl;


  //
  // Partial contributions to the Kultun sum rule
  //
  *xe_t = *xe_tVV = *xA = 0.0;
  Calc_En(SlfE, nsolb, NSOLDIM, solb, xe_t, xe_tVV, xA);
  cout << "\n 1-body=" << *xe_t << " ,  <Vpp>=" << (*xe_tVV-*xe_t)/2.0
       << " ,   total en.=" << (*xe_tVV+*xe_t)/2.0 << endl;
  cout << "\n partial occupation:  " << *xA << endl << endl;


  x3=0.0;
  cout << endl << endl;
  cout << " Number of forward solutions found:  "  << nsolf << endl;
  for(j=0; j<nsolf; ++j) {
    if ((nsolf-j > 20) && (0.01 > solf[1+NSOLDIM*j])) continue;
    cout << j << " , xe=" << solf[  NSOLDIM*j]
    << "   Z(%)=" << 100.*solf[1+NSOLDIM*j] << "     ";
    for(i=0; i<NSOLDIM; ++i) cout << solf[i+NSOLDIM*j] << "  ";
    cout << endl;
    x3+=solf[1+NSOLDIM*j];
    }


  cout << " Number of backward solutions found:  "  << nsolb << endl;
  for(j=0; j<nsolb; ++j) {
    if ((j > 20) && (0.01 > solb[1+NSOLDIM*j])) continue;
    cout << j << " , xe=" << solb[  NSOLDIM*j]
    << "   Z(%)=" << 100.*solb[1+NSOLDIM*j] << "     ";
    for(i=0; i<NSOLDIM; ++i) cout << solb[i+NSOLDIM*j] << "  ";
    cout << endl;
    x3+=solb[1+NSOLDIM*j];
    }
  cout << " Sum=:  "  << x3 << endl << endl;

  delete [] solDys;

  return (nf_lost+nb_lost);}  //numero di soluzioni eventualmente perse


