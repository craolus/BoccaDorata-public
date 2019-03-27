//
//
//  Classes to build the Polarization Propagator through the ph-BSE
//
//  C.Barbieri, GSI, December 2006.  (C.Barbieri@gsi.de)
//

#ifndef _POL_PROP_H
#define _POL_PROP_H

#define NDATAPIXY   9

//#include "OneBd_mtx.hh"
class OBmtx_t;

class Pol_prop_t {
  
  ModSpace_t *MdSp;
  VppInt_t   *vpp;
  SpProp_t   *gsp;
 public:
  int    Jf, dT, pif;
  double *e_ph,        *e_hp;
  int    *n_ph, *k_ph, *n_hp, *k_hp;
  int     n_ph_alloc,   n_hp_alloc;
  int     n_ph_basis,   n_hp_basis;

  //
  // NOTE: "n_f(b)_basis"  are meant to store the dimensions of the X and Y
  //    amplitudes, while "n_ph(hp)_basis" are the bases' dimensions. In 
  //    practice whenever these are both larger zero, they *must* always equal.
  //    Tere are cases in which "n_f(b)_basis" could be zero even if the
  //    2qp/2qh basis is not vanishing, for example when only a part of the 
  //    full RPA matrix is diagonalized (only fw TDA).
  //
  //     NOTE that the "n_f(b)_basis" *must* always be >= zero.
  //     However, if the basis is not allocated, the "n_ph(hp)_basis" can
  //    be negative.
  //
  int n_f_basis,  n_b_basis;   //  same as n_ph_basis  &  n_hp_basis but
                               // forced to be >= 0 (by the Solve_phRPA fnct.)

  int n_fw_sols,  n_bk_sols,  n_tot_sols;  // actual number of sols in memory

  int XY_nalloc,XY_nvects,NdataXY;
  double *XYmtx;

  friend class CoupledCluster_t;

  //
  // constructors/destructors:
  Pol_prop_t(void );
  Pol_prop_t(ModSpace_t*, SpProp_t*, VppInt_t* );
  Pol_prop_t(ModSpace_t*, SpProp_t* );
  Pol_prop_t(SpProp_t * );
  Pol_prop_t(ModSpace_t* );
  Pol_prop_t(VppInt_t* );
 ~Pol_prop_t(void );

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
  void free_matrices(void );
  void Clear_XYmtx(void );

  //
  // allocate memory:
  int  allocate_2qp_basis(int, int );
  int  allocate_XY(int, int );
  
 public:
  // qp-qh basis
  int  Count_phRPA_basis(int, int, int, int*, int*, int icut=0, int ntcut=0);
  int  Build_phRPA_basis(int, int, int, int icut=0, int ntcut=0, int iplot=0);
  int  Fill_phRPA_mtx(int, double*, int , double *Free_poles=NULL, bool Add_free_poles=true);
  //
  int get_Pi_f_indx(int*, int, int);
  int get_Pi_b_indx(int*, int, int);
  //int init_2qp_search(void );

 };


#endif
