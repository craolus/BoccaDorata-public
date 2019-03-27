//
//
//  Classes to build the Polarization Propagator through the ph-BSE
//
//  C.Barbieri, GSI, December 2006.  (C.Barbieri@gsi.de)
//

#ifndef _POL_PROP_H
#define _POL_PROP_H

#define NDATAPIXY   9
#define NDATAPIZAB  7
#define NDATAPIDAB  3

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

//XDIST//  int n_RPA_unst;       // counts the unstable (==imaginary) RPA solutions

  int XY_nalloc,XY_nvects,NdataXY;
  double *XYmtx;

//XDIST//
  /*
   int Zab_nalloc,Zab_nvects,NdataZ;
  double *Zabmtx;

  int Dab_nalloc,Dab_nvects,NdataD;
  double *Dabmtx;

  int *ia_2b, *ib_2b;
  int N_2BMSP_alloc, n_2bmsp;

  double *FaddBx_f1, *FaddBx_f2;
  double *FaddBx_b1, *FaddBx_b2;
*/
  
//XDIST//  char STOREDIR[100];

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
//void set_phint(VppInt_t* );

  //
  // other functions:
 private:
  // init_NULL
  void Ndata_init(void );
  void MVG_init_NULL(void ); // ..._init_NULL fncts. should stay private 
  void bases_init_NULL(void );
  void basis_2qp_init_NULL(void );
//XDIST//  void basis_2bmsp_init_NULL(void );
  void matrices_init_NULL(void );
  void XY_init_NULL(void );
//XDIST//
  /*
  void Zab_init_NULL(void );
  void Dab_init_NULL(void );
  void FaddBx_Fw_init_NULL(void );
  void FaddBx_Bk_init_NULL(void );
   */
  //
  // free(clear) allocated memory:
  void free_bases(void );
  void free_2qp_basis(void );
//XDIST//  void free_2bmsp_basis(void );
  void free_matrices(void );
  void Clear_XYmtx(void );
//XDIST//
  /*
   void Clear_Zabmtx(void );
  void Clear_Dabmtx(void );
  void Clear_FaddBx_Fw(void );
  void Clear_FaddBx_Bk(void );
*/
  //
  // allocate memory:
  int  allocate_2qp_basis(int, int );
//XDIST//  int  allocate_2bmsp_basis(int );
  int  allocate_XY(int, int );
//XDIST//
  /*
  int  allocate_Zab(void );
  int  allocate_Dab(void );
  int  allocate_FaddBx_Fw(void );
  int  allocate_FaddBx_Bk(void );
*/
  
 public:
//XDIST//
  /*
  // 2-body basis
  int  Count_2bmsp_basis(int, int, int, int* );
  int  Build_2bmsp_basis(int, int, int, int iplot=0);
  void save_2bmsp_basis(ofstream *);
  int  load_2bmsp_basis(ifstream *);
  //
  //int get_ab_indx(int*, int, int);
*/
  // qp-qh basis
  int  Count_phRPA_basis(int, int, int, int*, int*, int icut=0, int ntcut=0);
  int  Build_phRPA_basis(int, int, int, int icut=0, int ntcut=0, int iplot=0);
//XDIST//  int  Solve_phRPA(int i_RPAcalc=3, bool check_imm_sol=false);
  int  Fill_phRPA_mtx(int, double*, int , double *Free_poles=NULL, bool Add_free_poles=true);
//XDIST//  void save_2qp_basis(ofstream *);
//XDIST//  int  load_2qp_basis(ifstream *);
  //
  int get_Pi_f_indx(int*, int, int);
  int get_Pi_b_indx(int*, int, int);
  //int init_2qp_search(void );
  //int get_Pi_f_indx_2(int*, int, int);
  //int get_Pi_b_indx_2(int*, int, int);

  //XDIST//
  /*
  // computation of the 'Zab' and 'Dab' vectors: 
  int  Conv_XY_to_Zab(void );
  int  Conv_XY_to_Dab(void );
  int  Conv_Zab_to_Dab(void );
  int  XY_to_Zab_matrix(double[]);
  int  XY_to_Dab_matrix(double[]);
  int  Zab_to_Dab_matrix(double[]);

  //
  // Matrices for the Faddeev eqs.
  int Make_Faddeev_Boxes(void );
  int load_Faddeev_Boxes(int iload=0  );
  inline double *getfw1_c( int i_c ) {return (FaddBx_f1 + i_c*n_fw_sols);}
  inline double *getfw2_c( int i_c ) {return (FaddBx_f2 + i_c*n_fw_sols);}
  inline double *getbk1_c( int i_c ) {return (FaddBx_b1 + i_c*n_bk_sols);}
  inline double *getbk2_c( int i_c ) {return (FaddBx_b2 + i_c*n_bk_sols);}
  inline double getfw1(int i_r, int i_c) {return *(FaddBx_f1 + i_c*n_fw_sols + i_r);}
  inline double getfw2(int i_r, int i_c) {return *(FaddBx_f2 + i_c*n_fw_sols + i_r);}
  inline double getbk1(int i_r, int i_c) {return *(FaddBx_b1 + i_c*n_bk_sols + i_r);}
  inline double getbk2(int i_r, int i_c) {return *(FaddBx_b2 + i_c*n_bk_sols + i_r);}
  inline void getfw(double *f1, double *f2, int i_r, int i_c)
          {i_c=i_c*n_fw_sols+i_r;
          (*f1) = *(FaddBx_f1 + i_c); (*f2) = *(FaddBx_f2 + i_c);
          return;}
  inline void getbk(double *b1, double *b2, int i_r, int i_c)
          {i_c=i_c*n_bk_sols+i_r;
          (*b1) = *(FaddBx_b1 + i_c); (*b2) = *(FaddBx_b2 + i_c);
          return;}


  void save_Zabmtx(ofstream *);
  int  load_Zabmtx(ifstream *);
  void save_Dabmtx(ofstream *);
  int  load_Dabmtx(ifstream *);
  void save_XYmtx(ofstream *);
  int  load_XYmtx(ifstream *);
  int  save_Zab(void );
  int  load_Zab(void );
  int  save_Dab(void );
  int  load_Dab(void );
  int  save_XY(void );
  int  load_XY(void );
*/
  
  //XDIST//
  /*
  int  Build_Strength_Op(double[], int, int , double charge0=1.0, double charge1=1.0);
  int  Compute_response(int, int , double charge0=1.0, double charge1=1.0 );
  int  Compute_responseOLD(int );
  int  Plot_XY(char* , int , double , int , int n_out=1);
  int  XY_to_Zab_matrixOLD(double[]);
  int  Compute_Zab_str(int );
  int  Plot_PiPoles(int nf=-10, int nb=-10, int nplot=3, int nout=1);
  int  Plot_Ampl(char* , char* , int );
  int  Compute_density_prof(double[], double[], int, int, int, 
                            double charge0=1.0, double charge1=1.0 );
 private:
  int  Ylm_op_file_rd(OBmtx_t*, double[], int, double[], int, int, double, double );
  int  Ylm_op(double[], int, double[], int, int, double, double );
  int  Gamow_Teller_op(double[]);
  int  Fermi_op(double[]);
  int  Density_op(double[], double, double[] );
   */
 };


#endif
