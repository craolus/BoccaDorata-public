//
//
//   Construct and store the self energy (both HF and dynamc part) for a given
//  channel taking as input:
//      model spacce, interaction and sp propagator.
//
//

//#include <iostream>
//#include <fstream>
//#include <iomanip>
using namespace std;

#include <cstdlib>
//#include <cstdio>
//#include <cmath>

#include "BcDor-Global_variables.hh"
#include "BcDor-SlfEn_classes.hh"


void Wide_self_energy_t::set_ishMVG_NULL(void ) {
  MdSp = NULL;
  Vpp  = NULL;
  gin  = NULL;
  I_SH  = -100;
  NTOT_DEN = 2;
  NTOT_ORB = -100;
  NTOT_DIM = -100;
  N_EXT_PWV = -100;;
  sprintf(STOREDIR, BcDorWorkFolder);
  return;}  

void Wide_self_energy_t::init_1b_NULL(void ) {
  Sig_1b   = NULL;
  Sig_HF   = NULL;
  Sig_BHF  = NULL;
  Gmtx_ste = NULL;
  N_GMTX = 0;
  return;}

void Wide_self_energy_t::init_ext_orb_NULL(void ) {
  i_pwv    = NULL;
  j2_pwv   = NULL;
  ip_pwv   = NULL;
  ch_pwv   = NULL;
  nst_pwv  = NULL;
  norb_pwv = NULL;
  nsp_loc  = NULL;
  isp_msp  = NULL;
  return;}

void Wide_self_energy_t::init_dyn_NULL(void ) {
  Sig_dyn_fw = NULL;
  Sig_dyn_bk = NULL;
  N_PLS_fw = 0;
  N_PLS_bk = 0;
  N_PLS_FW_alloc = -100;
  N_PLS_BK_alloc = -100;
  n_piv_fw = -1;
  n_piv_bk = -1;
  i_piv_fw = NULL;
  i_piv_bk = NULL;
  return;}

void Wide_self_energy_t::init_ALL_NULL(void ) {
  init_1b_NULL();
  init_ext_orb_NULL();
  init_dyn_NULL();
  set_ishMVG_NULL();
  return;}

void Wide_self_energy_t::free_1b_mem(void ) {
  if (NULL != Sig_1b  ) delete [] Sig_1b;
  if (NULL != Sig_HF  ) delete [] Sig_HF;
  if (NULL != Sig_BHF ) delete [] Sig_BHF;
  if (NULL != Gmtx_ste) delete [] Gmtx_ste;
  init_1b_NULL();
  return;}  

void Wide_self_energy_t::free_ext_orb_mem(void ) {
  if (NULL != i_pwv    ) delete [] i_pwv;
  if (NULL != j2_pwv   ) delete [] j2_pwv;
  if (NULL != ip_pwv   ) delete [] ip_pwv;
  if (NULL != ch_pwv   ) delete [] ch_pwv;
  if (NULL != nst_pwv  ) delete [] nst_pwv;
  if (NULL != norb_pwv ) delete [] norb_pwv;
  if (NULL != nsp_loc  ) delete [] nsp_loc;
  if (NULL != isp_msp  ) delete [] isp_msp;
  init_ext_orb_NULL();
  return;}  

void Wide_self_energy_t::free_dyn_mem(void ) {
  if (NULL != Sig_dyn_fw) delete [] Sig_dyn_fw;
  if (NULL != Sig_dyn_bk) delete [] Sig_dyn_bk;
  if (NULL != i_piv_fw) delete [] i_piv_fw;
  if (NULL != i_piv_bk) delete [] i_piv_bk;
  init_dyn_NULL(); // so that other routines know that no slef-en is allocated
  return;}  

void Wide_self_energy_t::free_mem(void ) {
  free_1b_mem();
  free_dyn_mem();
  free_ext_orb_mem();
  return;}  


void Wide_self_energy_t::Allocate_ext_orb_mem(int n_part_waves/*=-10*/, int n_orbits/*=-10*/) {
  if (  n_orbits   > 0) NTOT_ORB  = n_orbits;
  if (n_part_waves > 0) N_EXT_PWV = n_part_waves;
  free_ext_orb_mem();
  if (N_EXT_PWV > 0) {
    i_pwv    = new int[N_EXT_PWV];
    j2_pwv   = new int[N_EXT_PWV];
    ip_pwv   = new int[N_EXT_PWV];
    ch_pwv   = new int[N_EXT_PWV];
    nst_pwv  = new int[N_EXT_PWV];
    norb_pwv = new int[N_EXT_PWV];
  }
  if (NTOT_ORB > 0) {
    nsp_loc  = new int[NTOT_ORB];
    isp_msp  = new int[NTOT_ORB];
  }
  return;}

static int n1;

void Wide_self_energy_t::Allocate_Dynamic_SlfEn_Fw(int npl_fw, int N_DEN/*=-100*/) {
                                    //
                                    //  If the self energy is competely reallocated,
                                    // NTOT_DEN is set to max(2, N_DEN)
                                    //

  if ( (NTOT_DEN < 2) || (NTOT_ORB < 1) || (NTOT_DIM != NTOT_DEN + NTOT_ORB) ) {
    Allocate_Dynamic_SlfEn(npl_fw, 0, N_DEN);
    return;
  }

  N_PLS_FW_alloc = (1 > npl_fw) ? 1 : npl_fw;
  if (NULL != Sig_dyn_fw) delete [] Sig_dyn_fw; Sig_dyn_fw = NULL;
  Sig_dyn_fw=new double[N_PLS_FW_alloc*NTOT_DIM];
  N_PLS_fw = 0;
  if (NULL != i_piv_fw) delete [] i_piv_fw; i_piv_fw = NULL;
  i_piv_fw = new int[NTOT_DIM];
  n_piv_fw = 0; for(n1=0; n1<NTOT_DIM; ++n1) i_piv_fw[n1] = -100;

  return;
  }

void Wide_self_energy_t::Allocate_Dynamic_SlfEn_Bk(int npl_bk, int N_DEN/*=-100*/) {
                                    //
                                    //  If the self energy is competely reallocated,
                                    // NTOT_DEN is set to max(2, N_DEN)
                                    //
                                    
  if ( (NTOT_DEN < 2) || (NTOT_ORB < 1) || (NTOT_DIM != NTOT_DEN + NTOT_ORB) ) {
    Allocate_Dynamic_SlfEn(0, npl_bk, N_DEN);
    return;
  }

  N_PLS_BK_alloc = (1 > npl_bk) ? 1 : npl_bk;
  if (NULL != Sig_dyn_bk) delete [] Sig_dyn_bk; Sig_dyn_bk = NULL;
  Sig_dyn_bk=new double[N_PLS_BK_alloc*NTOT_DIM];
  N_PLS_bk = 0;
  if (NULL != i_piv_bk) delete [] i_piv_bk; i_piv_bk = NULL;
  i_piv_bk = new int[NTOT_DIM];
  n_piv_bk = 0; for(n1=0; n1<NTOT_DIM; ++n1) i_piv_bk[n1] = -100;

  return;
  }

void Wide_self_energy_t::Allocate_Dynamic_SlfEn(int npl_fw, int npl_bk, int N_DEN/*=-100*/) {
                          //
                          //  The input N_DEN it sets the number of mtx. els.
                          // that can be set for each pole. It must be at
                          // least == 2, to store the pole and strength.
                          // NDEN > 2 is for storing the inverse mass operator
                          // in 'fishbone' form, when it is not diagonalised.
                          //

  free_dyn_mem();
  NTOT_DEN = ( N_DEN   > 2) ?  N_DEN   : 2;
  NTOT_ORB = (NTOT_ORB > 1) ? NTOT_ORB : 1;
  NTOT_DIM = NTOT_DEN + NTOT_ORB;

  N_PLS_FW_alloc = (1 > npl_fw) ? 1 : npl_fw;
  N_PLS_BK_alloc = (1 > npl_bk) ? 1 : npl_bk;
  //
  Sig_dyn_fw=new double[N_PLS_FW_alloc*NTOT_DIM];
  Sig_dyn_bk=new double[N_PLS_BK_alloc*NTOT_DIM];
  //
  i_piv_fw = new int[NTOT_DIM];
  i_piv_bk = new int[NTOT_DIM];

  // Other possibility:
  /*
  if (npl_fw>0) {Sig_dyn_fw=new double[N_PLS_FW_alloc*NTOT_DIM]; N_PLS_FW_alloc = npl_fw; i_piv_fw = new int[NTOT_DIM];}
  if (npl_bk>0) {Sig_dyn_bk=new double[N_PLS_BK_alloc*NTOT_DIM]; N_PLS_BK_alloc = npl_bk; i_piv_bk = new int[NTOT_DIM];}
  */

  N_PLS_fw = 0;
  N_PLS_bk = 0;
  n_piv_fw = 0; for(n1=0; n1<NTOT_DIM; ++n1) i_piv_fw[n1] = -100;
  n_piv_bk = 0; for(n1=0; n1<NTOT_DIM; ++n1) i_piv_bk[n1] = -100;

  return;
  }


void Wide_self_energy_t::set_MdSpVppGsp (ModSpace_t *msp, VppInt_t *pot, SpProp_t *gsp) {
  free_mem();
  MdSp = msp;
  Vpp  = pot;
  gin  = gsp;
  return;}


static int n2, n3;

static int /*ish_fd,*/ j2_fd, ip_fd, dt_fd;
static int   ish_a,    j2_a,  ip_a,  dt_a;

void Wide_self_energy_t::set_clj (int ish_fd, int lmbd/*=0*/, int d_ip/*=0*/, int d_t/*=0*/) {
  if (NULL == MdSp) {
    cerr << "\n\n Error in class' function 'Wide_self_energy_t::set_clj':" 
       << " it is not possible to initialise the self-energy to channel " << ish_fd
       << " if no model space has been specified.\n\n";
    }
  free_mem();

  I_SH  = ish_fd;
  j2_fd = MdSp->MSp_2j[ish_fd];
  ip_fd = MdSp->MSp_ip[ish_fd];
  dt_fd = MdSp->MSp_ch[ish_fd];

  Lambda = lmbd;
  parity = d_ip;
  charge = d_t;

  n1 = 0;
  n2 = 0;
  for(ish_a = 0; ish_a < MdSp->nsubsh; ++ish_a) {
    if (( abs(ip_fd + parity + MdSp->MSp_ip[ish_a])%2 )    ||
        ( dt_fd + charge != MdSp->MSp_ch[ish_a]  )         ||
        ( am.TriIneq(MdSp->MSp_2j[ish_a], 2*Lambda, j2_fd) )   ) continue;
    ++n1;
    n2 += MdSp->MSp_no[ish_a];
  }

  Allocate_ext_orb_mem(n1, n2); // # of partial wavwe, # of orbits

  n1 = 0;
  n2 = 0;
  for(ish_a = 0; ish_a < MdSp->nsubsh; ++ish_a) {
    j2_a = MdSp->MSp_2j[ish_a];
    ip_a = MdSp->MSp_ip[ish_a];
    dt_a = MdSp->MSp_ch[ish_a];

    if (( abs(ip_fd + parity + ip_a)%2 )    ||
        ( dt_fd + charge != dt_a )          ||
        ( am.TriIneq(j2_a, 2*Lambda, j2_fd) )   ) continue;

    i_pwv   [n1] = ish_a;
    j2_pwv  [n1] = j2_a;
    ip_pwv  [n1] = ip_a;
    ch_pwv  [n1] = dt_a;
    nst_pwv [n1] = n2;
    norb_pwv[n1] = MdSp->MSp_no[ish_a];
    for(n3=0; n3<MdSp->MSp_no[ish_a]; ++n3) {
      nsp_loc [n2] = n1;
      isp_msp [n2] = (MdSp->MSp_n[ish_a]  -  MdSp->Mor_n) + n3;
      ++n2;
    }
    ++n1;

  }  // ish_a loop

  if ((n2 != NTOT_ORB) || (n1 != N_EXT_PWV)) {
    cerr << "Trouble with NTOT_ORB or N_EXT_PWV...\n\n";
    exit(100);
  }
  return;}

void Wide_self_energy_t::reset(int i, ModSpace_t *msp, VppInt_t *pot, SpProp_t *gsp) {
  free_mem();
  init_ALL_NULL();
  set_MdSpVppGsp (msp, pot, gsp);
  set_clj(i);
  cout << "\nSelf-energy object reset to new mod. space, AND to the partial wave "
       << MdSp->MSp_name[I_SH] << endl;
  return;}

void Wide_self_energy_t::reset(int i, int lmbd, int d_ip, int d_t,
                               ModSpace_t *msp, VppInt_t *pot, SpProp_t *gsp) {
  free_mem();
  init_ALL_NULL();
  set_MdSpVppGsp (msp, pot, gsp);
  set_clj(i, lmbd, d_ip, d_t);
  cout << "\nSelf-energy object reset to new mod. space, AND to the partial wave "
       << MdSp->MSp_name[I_SH]
       << ". The q.# for the external field are also set for L="
       << Lambda << ",ip=" << parity <<",dt= " << charge << "." << endl;
  return;}

void Wide_self_energy_t::reset(ModSpace_t *msp, VppInt_t *pot, SpProp_t *gsp) {
  free_mem();
  init_ALL_NULL();
  set_MdSpVppGsp (msp, pot, gsp);
  cout << "\nSelf-energy object reset to new mod.space\n";
  return;}


Wide_self_energy_t::~Wide_self_energy_t () {
  cout << "\nSelf-energy object deallocated.\n\n";
  free_mem();
  //init_ALL_NULL();
  return;}


Wide_self_energy_t::Wide_self_energy_t () {
  init_ALL_NULL();
  cout << "\nSelf-energy object created without any mod.sp./intercation/prop specified\n\n";
  return;}

Wide_self_energy_t::Wide_self_energy_t (ModSpace_t *msp, VppInt_t *pot, SpProp_t *gsp) {
  init_ALL_NULL();
  set_MdSpVppGsp (msp, pot, gsp);
  cout << "\nSelf-energy object created.\n";
  return;}

Wide_self_energy_t::Wide_self_energy_t (int i, ModSpace_t *msp, VppInt_t *pot, SpProp_t *gsp) {
  init_ALL_NULL();
  set_MdSpVppGsp (msp, pot, gsp);
  set_clj(i);
  cout << "\nSelf-energy object created and set to partial wave"
       << MdSp->MSp_name[I_SH] <<endl;
  return;}

Wide_self_energy_t::Wide_self_energy_t (int i, int lmbd, int d_ip, int d_t,
                               ModSpace_t *msp, VppInt_t *pot, SpProp_t *gsp) {
  init_ALL_NULL();
  set_MdSpVppGsp (msp, pot, gsp);
  set_clj(i, lmbd, d_ip, d_t);
  cout << "\nSelf-energy object created and set to partial wave"
       << MdSp->MSp_name[I_SH]
       << ". The q.# for the external field are also set to L="
       << Lambda << ",ip=" << parity <<",dt= " << charge << "." <<endl;
  return;}

