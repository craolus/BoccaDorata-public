//
// Class 'ModSpace_t'   stores the information on the single particle basis
//                     and creates two sets of arrays. The one named 'Msp_***'
//                    group the s.p. orbits according to their partial wave
//                    (subshell), i.e. the q.#s j, l , ip, ch. In each group,
//                    the orbits differ by their principal q.# 'n'.
//                     The arrays named 'Mor_***' provide a (one dimensional)
//                    list of all the orbits.
//
//
// Class 'SpProp_t'    Stores the Lehmann representation of the single particle
//                    propagattor. That is, the ovlerlap functions and energies
//                    of a --discrete-- set of quasiparticles and quasiholes.
//                     NOTE that a discretized system is assumed, the initial 
//                    N-body wave function is assumed in a spherical state
//                   (J==O), and the ovelrap orbits are expanded in the s.p.
//                   orbits provided by the associated model space `*MdSp'
//                   (pointing to an object of class ModSpace_t).
//
//
// C.Barbieri, GSI, December 2006.  (C.Barbieri@gsi.de).
//
//


#ifndef  _SP_CLASSES_H
#define  _SP_CLASSES_H

#include <fstream>


class ModSpace_t {
   char  *MSp_chall;
  public:
   int    nsubsh, n_spb_tot;
   int    *MSp_l,   *MSp_2j, *MSp_ip, *MSp_ch, *MSp_no;
   int    *MSp_nqp, *MSp_nqh;
   int    **MSp_n, **MSp_F;
   double **MSp_xe; 
   char   **MSp_name;

   double *Mor_xe; 
   int    *Mor_n,   *Mor_sh, *Mor_F;
   int    *Mor_l, *Mor_2j, *Mor_ip, *Mor_ch;

   int Jmax;                    // bounds on angular momenta and charge
   int ch_pp_mn,  ch_pp_mx;
   int ch_hh_mn,  ch_hh_mx;
   int ch_ph_mn,  ch_ph_mx;
   int ch_pph_mn, ch_pph_mx;
   int ch_hhp_mn, ch_hhp_mx;

   double htom,bho,bsq,mass;  //  Physical paramerter in the case of
                              // an harmonic oscillator basis.

   ModSpace_t();
   ModSpace_t(int, int, int );
   ModSpace_t(char* );
  ~ModSpace_t();
   void Reset_phys_const(void );
   void Set_phys_const_bHO(double );
   void Set_phys_const_hwHO(double );
   void Reset_init_values(void );
   void Init_NULL_shells(void );
   void Init_NULL_orbits(void );
   void Clear_all(void);
   void Free_mem_shells(void );
   void Free_mem_orbits(void );
   void Allocate_shells(int );
   void Allocate_orbits(int );
   int write(char* );
// int read(char* );  //  This is not needed, really, one can deallocate
                      // object and read from the file using the constructor.
   int cout_msp();
   int get_orbit(int, int, int, int );
   int get_subsh(     int, int, int );
   void Calculate_Jch_bounds(int iplot=0);
   void count_2b_confs(int*, int*, int iplot=0);
  };


//#include "Int_classes.hh"
class VppInt_t;


// NOTE: The follofing needs a mod space to expand the overlap orbits
class SpProp_t {
   int NTOT_QP,NTOT_QH; // size of allocated memeory
  public:
   ModSpace_t *MdSp; //  Pointer to the associated model space
 //  int    nsubsh;    //  Gets this from the assoc. mod. sp. Need later 
 //                    // to deallocate the whole thing...

   int    n_qh,n_qp,N_SPB_MX; // tot number of qh, qp, and max # of components
   int    *n_clj_p, *n_clj_h;
   double      *ep,      *eh;
   double      *zp,      *zh;
   double     *fqp,     *fqh;
   //
   int        *n1p,     *n1h; // ...used to pass prcessing information, such as
   int        *n2p,     *n2h; // which poles to include in the Faddeev calculations
   double     *d1p,     *d1h; // or which qp/qh to osed as Krylov pivots...

   double **Sp_clj_ep,  **Sp_clj_eh;  // pointers to qp energies
   double **Sp_clj_zp,  **Sp_clj_zh;  // pointers to spec. factors
   double **Sp_clj_fqp, **Sp_clj_fqh; // pointers to overlap w.f.s 
   //
   int     *Sp_clj_np,   *Sp_clj_nh;  // # of ovelap functs. for each subshell
   //
   int    **Sp_clj_n1p, **Sp_clj_n1h; // pointers to n1 data
   int    **Sp_clj_n2p, **Sp_clj_n2h; // pointers to n2 data
   double **Sp_clj_d1p, **Sp_clj_d1h; // pointers to d1 data

   SpProp_t();
   SpProp_t(ModSpace_t *, int i_xe=0);
   SpProp_t(char*, ModSpace_t *, int *i_rd_flag=NULL);
   SpProp_t(int, int, ModSpace_t *);
  ~SpProp_t();
   void set_bounds_zero(void );
   void set_alloc_size_zero(void );
   void init_NULL(void );
   void init_ALL_values(void );
   void free_mem(void );
   void Clean_up(void );
   int Allocate_mem(int , int );
   int write(char* , int wcode=0);
   int cout_qplist(int wcode=0);
   int cout_byshell(int wcode=0);
   void count_qp_confs(int*, int*, int*, int*, int*, int iplot=1, int t_2qp=-1, int t_3qp=-1);
   int add_qp(int, int, int, double[]);
   int add_qh(int, int, int, double[]);
   int scale_qp_amplitude(int, double );
   int scale_qh_amplitude(int, double );
   int set_qp_amplitude(int, double );
   int set_qh_amplitude(int, double );
   int Fill_with_charge(int );
   void set_n1(int, int, int* vp=NULL, int* vh=NULL);
   void set_n2(int, int, int* vp=NULL, int* vh=NULL);
   void set_d1(double, double );
   void set_d1(int, int, double* , double* );
   //
   int    ExtendedBHF_Pot(double* , int, VppInt_t*, int, int, int, double *StrEnList=NULL);
   int    ExtendedHartreeFock_Pot(double* , int, VppInt_t*, int);
   double Compute_FirstOrderExpValue(VppInt_t*, int ist=0, int iend=10000);
   double quick_n_dirty_diff(SpProp_t*);

   void DysPivots(char*, int n_bk=1, int n_fw=1, double z_cut=60.0);

   //double Seek_EFermi_tj2p(int, int, int, int iplot=1);
   double Seek_EFermi_tj2p(int, int, int, double *Ef_up=NULL, double *Ef_dn=NULL, int iplot=1);
   int Reduce(SpProp_t*, int, int itype_calc=0);
   int RedDensMtx(SpProp_t * );
   int Build_inverted_moments(double[], int, double, int,  int, int, int);
   int    Build_plain_moments(double[], int, double, int,  int, int, int);
   int Tune_poles(int, int, int, int, int, double*, int);
   int Evaluate_SR_moments(double*, double*, int, int k_st=0, int k_max=1, double omega0=0.0);
   int Evaluate(double[], int, int, double, int iside=0, int n_der=0);

   void Plot_3D_SpectFnct(char* );
  };



class Save_g1_Poles {
 public:
  int       N_Pls_Saved, N_Pls_Subst, N_Dim_Ampls;
  int      *v_ish, *v_fb, *v_no;
  double   *v_ampls;
  int       status;
  SpProp_t *GSP_curr;

  Save_g1_Poles();
  Save_g1_Poles(SpProp_t*, char*, int *stout=NULL);
 ~Save_g1_Poles();
  void init_NULL(void );
  void free_mem(void );
  void alloc_mem(int, int );

  int load_info(SpProp_t*, char*);
  int separate_pls(int ); 
  int restore_pls(int );
  int add_fw_pls(SpProp_t*, int );
  int add_bk_pls(SpProp_t*, int );
};


class Extnl_msp {
  int *n,  *twj,  *ip,  *ch;
 public:
  ModSpace_t *MdSp; //  Pointer to the associated model space
  double bHO;
  int  Atot;
  int Ntot_orbs;
  int *i_ext, *i_msp;
  double *Esp;

  ~Extnl_msp();
   Extnl_msp();
   Extnl_msp(char*);
  void init_NULL(void );
  void free_mem(void );
  void alloc_mem(int );

  int load_oslo_msp(char* );
  int load_oslo_msp(FILE* );
  int Sync(ModSpace_t*);

};


#endif
