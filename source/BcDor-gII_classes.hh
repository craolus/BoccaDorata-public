//
//
//  Classes to build the two-body Propagator gII through the pp/hh-FGE
//
//  C.Barbieri, GSI, December 2006.  (C.Barbieri@gsi.de)
//

#ifndef _TWOBODY_PROP_H
#define _TWOBODY_PROP_H

#define NDATALDXY   7


class gII_prop_t {
  
  ModSpace_t *MdSp;
  VppInt_t   *vpp;
//  VphInt_t   *vph;
  SpProp_t   *gsp;
 public:
  int    Jf, dT, pif;
  double *e_pp,          *e_hh;
  int    *n1_pp, *n2_pp, *k1_hh, *k2_hh;
  int     n_pp_alloc,     n_hh_alloc;
  int     n_pp_basis,     n_hh_basis;

  //
  // NOTE: "n_f(b)_basis"  are meant to store the dimensions of the X and Y
  //    amplitudes, while "n_pp(hh)_basis" are the bases' dimensions. In 
  //    practice whenever these are both larger zero, they *must* always equal.
  //    Tere are cases in which "n_f(b)_basis" could be zero even if the
  //    2qp/2qh basis is not vanishing, for example when only a part of the 
  //    full RPA matrix is diagonalized (only fw, or bk, TDA).
  //
  //     NOTE that the "n_f(b)_basis" *must* always be >= zero.
  //     However, if the basis is not allocated, the "n_pp(hh)_basis" can
  //    be negative.
  //
  int n_f_basis,  n_b_basis;    //  same as n_pp_basis & n_hh_basis but forced
                                // to be >= 0 (by the Solve_LaddhRPA fnct.)

  int n_fw_sols,  n_bk_sols,  n_tot_sols;  // actual number of sols in memory

  int XY_nalloc,XY_nvects,NdataXY;
  double *XYmtx;

  friend class CoupledCluster_t;

  //
  // constructors/destructors:
  gII_prop_t(void );
  gII_prop_t(ModSpace_t*, SpProp_t*, VppInt_t* );
  gII_prop_t(ModSpace_t*, SpProp_t* );
  gII_prop_t(SpProp_t* );
  gII_prop_t(ModSpace_t* );
  gII_prop_t(VppInt_t* );
 ~gII_prop_t(void );

  void set_space(ModSpace_t*, SpProp_t*, VppInt_t* );
  void set_space(ModSpace_t*, SpProp_t* );
  void set_space(SpProp_t * );
  void set_MdSp(ModSpace_t* );
  void set_ppint(VppInt_t* );

  //
  // other functions:
 private:
  // init_NULL
  void Ndata_init(void );
  void MVG_init_NULL(void ); // ..._init_NULL fncts. should stay private 
  void bases_init_NULL(void );
  void basis_2qp_init_NULL(void );
  void matrices_init_NULL(void );
  void XY_init_NULL(void );
  //
  // free(clear) allocated memory:
  void free_bases(void );
  void free_2qp_basis(void );
  void free_2bmsp_basis(void );
  void free_matrices(void );
  void Clear_XYmtx(void );
  //
  // allocate memory:
  int  allocate_2qp_basis(int, int );
  int  allocate_XY(int, int );
  
 public:
  
  // 2qp/2qh basis
  int  Count_LaddRPA_basis(int, int, int, int*, int*, int icut=0, int ntcut=0);
  int  Build_LaddRPA_basis(int, int, int, int icut=0, int ntcut=0, int iplot=0);
  int  Fill_LaddRPA_mtx(int, double*, int , double *Free_poles=NULL, bool Add_free_poles=true);
  //
  int get_ld_f_indx(int*, int, int, int*);
  int get_ld_b_indx(int*, int, int, int*);
  //int init_2qp_search(void );

};

#endif
