//
//
//   Class to create the self energy (both HF and dynamc part) for a given
//  channel taking as input:
//      model spacce, interaction and sp propagator.
//
//

#ifndef  _SLFEN_CLASSES_H
#define  _SLFEN_CLASSES_H

#include "BcDor-Sp_classes.hh"
#include "BcDor-Int_classes.hh"


class Wide_self_energy_t {
  public:
  VppInt_t   *Vpp;
  SpProp_t   *gin;
  ModSpace_t *MdSp;
  double *Sig_1b;
  double *Sig_HF;
  double *Sig_BHF, *Gmtx_ste;  int N_GMTX;
  double *Sig_dyn_fw;
  double *Sig_dyn_bk;
  int I_SH;
  int NTOT_DEN, NTOT_ORB, NTOT_DIM;
  int Lambda, parity, charge;
  int N_EXT_PWV;
  int *i_pwv, *j2_pwv, *ip_pwv, *ch_pwv, *nst_pwv, *norb_pwv;
  int *nsp_loc, *isp_msp;
  //
  int N_PLS_fw,N_PLS_bk;             // number of stored fw- and bk-going poles
  int N_PLS_FW_alloc,N_PLS_BK_alloc; // size of allocated memory
  //
  int  n_piv_fw,  n_piv_bk;
  int *i_piv_fw, *i_piv_bk;


  char STOREDIR[100]; // working directory (where to read+write self energies)


  Wide_self_energy_t();
  Wide_self_energy_t(      ModSpace_t*,  VppInt_t*,  SpProp_t* );
  Wide_self_energy_t(int,  ModSpace_t*,  VppInt_t*,  SpProp_t* );
  Wide_self_energy_t(int,  int, int, int, ModSpace_t*,  VppInt_t*,  SpProp_t* );
 ~Wide_self_energy_t ();

  void  set_ishMVG_NULL(void );
  void  init_1b_NULL(void );
  void  init_ext_orb_NULL(void );
  void  init_dyn_NULL(void );
  void  init_ALL_NULL(void );
  //
  void  free_1b_mem(void );
  void  free_dyn_mem(void );
  void  free_ext_orb_mem(void );
  void  free_mem(void );
  //
  void  Allocate_ext_orb_mem(int n_part_waves=-100, int n_orbits=-100);
  void  Allocate_Dynamic_SlfEn(int, int, int N_DEN=-100/*den. size*/  );
  void  Allocate_Dynamic_SlfEn_Fw(int, int N_DEN=-100/*den. size*/  );
  void  Allocate_Dynamic_SlfEn_Bk(int, int N_DEN=-100/*den. size*/  );
  //
  void  set_MdSpVppGsp (ModSpace_t *, VppInt_t *, SpProp_t *);
  void  set_clj (int, int lmbd=0, int d_ip=0, int d_t=0);
  void  reset(      ModSpace_t*,  VppInt_t*,  SpProp_t* );
  void  reset(int,  ModSpace_t*,  VppInt_t*,  SpProp_t* );
  void  reset(int,  int, int, int,  ModSpace_t*,  VppInt_t*,  SpProp_t* );


  void  Make_empty_ExtendedBHF(int );
  void  Build_ExtendedBHF(int, int );
  void  Make_empty_ExtendedHartreeFock(void );
  void  Build_ExtendedHartreeFock(void );
//  void  Build_DynSelfEn(void );
  void  Make_empty_1bd_SelfEn(void );
  void  Make_1bd_SelfEn(double [], int);
  void   Add_1bd_SelfEn(double [], int);
  void  Evaluate_Dyn_SelfEn(double [], int, double, int ifb=0 , double xelim=+1.e100 );  //, int isort=0);
  void  Evaluate_Der_SelfEn(double [], int, double, int ifb=0 , double xelim=+1.e100 );  //, int isort=0);
//  void  Set_this_as_WorkSlE(void );
  void  Count_2p1h_2h1p_poles( int* , int* , int icut=0, int ntcut=0);
  void  Build_Dynamic_SlfEn_2ndOrder(int* , int* , int inew=1, int icut=0, int ntcut=0);
  void  Build_Dynamic_SlfEn_2ndOrder_slow(int* , int* , int inew=1, int icut=0, int ntcut=0);


  int Reduce_SelfEn_Lncz(Wide_self_energy_t*, int , double** , int* ,  int*, int i_copyMF=0);
  int Diagonalize_SelfEn(Wide_self_energy_t*, int , double** , int* ,  int*, int i_copyMF=0);
  int Copy_MF(Wide_self_energy_t* );

  int save_SelfEn_bin(int,  char *Label="SelfEn");
  int load_SelfEn_bin(char *Label="SelfEn", int nFw_min=-1, int nBk_min=-1 );
private:
  int  load_SelfEn_bin_a(ifstream *, int[] );
  void load_SelfEn_bin_b(ifstream *, double*, int, int, int);
  void load_SelfEn_bin_c(ifstream *, int*, int*, int);
  void get_filename(char*, char*, char*, char*);
  // To read the old file format (for compatibility):
  public:
  int save_SelfEn_bin_short(int );
  int load_SelfEn_bin_short(int nFw_min=-1, int nBk_min=-1 );
  private:
  int  load_SelfEn_bin_a_short(ifstream *, int[] );
  void load_SelfEn_bin_b_short(ifstream *, double*, int, int, int);
  
public:
  void plot_SelfEn(char*, char*, double *EFermi=NULL);
  void plot_SelfEn_stc(char*, double *EFermi=NULL);
  void plot_SelfEn_dyn(char*, double *EFermi=NULL);
  };

#endif

