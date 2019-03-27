//
//
//   Classes for the solution of the coupled cluster CCD equations
//
//
//  C.Barbieri, Surrey, October 2011.  (C.Barbieri@surrey.ac.uk)
//


#ifndef _CCD_H
#define _CCD_H


#include "BcDor-Sp_classes.hh"
#include "BcDor-Int_classes.hh"
#include "BcDor-gII_classes.hh"
#include "BcDor-PolProp_classes.hh"

class CoupledCluster_t {
  
public:
  ModSpace_t *MSp_CC;
  VppInt_t   *Vpp_CC;
  //  VphInt_t   *Vph_CC;
  SpProp_t   *gsp_CC;
  
  int MULT_CH_PI,  MULT_CH_LD,  MULT_CH;
  int CH_OFFS_PI,  CH_OFFS_LD;
  int NT2_CHANNELS;
  
  gII_prop_t **T2_Ld;
  Pol_prop_t **T2_Pi;
  int *T2_Ld_dim, *T2_Pi_dim;
  //int npw_Ld, npw_Pi;
  
  CoupledCluster_t(ModSpace_t*, SpProp_t*, VppInt_t* );
  ~CoupledCluster_t(void );
  void init_NULL(void );
  
  int free_CCD_aplitudes(void );
  int create_CCD_aplitudes(int i_plot=1); // iplot is the level of verbosing
  //
  int Divide_deltaE(void );
  int Calc_E2(double * , double *EPi=NULL);

  int Pandya_Copy_Ld_to_Pi(void );
  int Pandya_Add_Pi_to_Ld(void );

  int Solve_CCD(double, int );

private:
  int Calc_New_T2(double*, double*, double*, double*, double**, double**, double**, double**, double, double );
  void Create_Snn_Skk(void );
  void Destroy_Snn_Skk(void );
  void Clear_Snn_Skk(double**, double** );

  //
  // for testing:
  void plot_nn(double **, double **);
  void plot_kk(double **, double **);
};

#endif
