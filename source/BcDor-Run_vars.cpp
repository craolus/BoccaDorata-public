//  
//
//

//#include <fstream>
//#include <iostream>
using namespace std;

#include "BcDor-Global_variables.hh"
#include "BcDor-Run_vars.hh"
#include "BcDor-PhysConsts.hh"


VppInt_t  *Rrms =  NULL;
ModSpace_t *MS =  NULL;
VppInt_t  *Vpp =  NULL;


//===================================================================
// run parameters:
int AddVc = 0;                              char Vcpp_file[BUFFER1];
int RrmsCalc = 0;                           char Rrmspp_file[BUFFER1];

int i_Ext_U1 = 0;                           char Ext_U1_file[BUFFER1];
double Ext_U1_Uo = 1.0;
double Ext_U1_Uch[2] = {0.0 , 0.0};

int IU1body = 3;  double U1body_fac = 1.0;  char Trel_1_file[BUFFER1];
double a_Tcm = 0.0;                         char Trel_2_file[BUFFER1];

                                            char Vpp_file[BUFFER1];
int iAddVpp = 0, AddVppMx=3;                char Vpp_add_file[3][BUFFER1];
double AddVppMult[3]={1.0, 1.0, 1.0};
                                            char Vpp_out_file[BUFFER1];
                                            char SpProp_file[BUFFER1];
                                            char MdSp_file[BUFFER1];
                                            char MdSp_oslo_file[BUFFER1];
                                            char DysPivots_file[BUFFER1];
                                            char gout_file[BUFFER1];
                                            char gen_fout[BUFFER1];
int i_gout_str = 0;                         char gout_str[BUFFER1];
                                            char Monopole_valence_file[BUFFER1];

                                            char temp_str[BUFFER1];


// Setting of Conventional/Unconventional baryon masses:
double mass_electron     = ELECTRONmass;
double mass_proton       = PROTONmass;
double mass_neutron      = NEUTRONmass;
double mass_particle_ave = NUCLEONmass;
int i_set_new_mass = 0; // ==0 no change, ==1 use ave nucl. mass


int A = -1, N=-1, Z=-1;
double bHO =-1.0;
double hwHO =-1.0;
int wAll = 0;
bool b_plot_gsp = false;
int ItrMax = 0, nDfw = -100, nDbk = -100;
int Lanczos = 0; // N.B. must be zero, sothat Lanczos is not set
                 // with no optoion specified.
int i_2nd_ord = 0, i_HF = 1, i_ExtSE=0;
int i_MeanField=0, i_MeanField_FwBk=0;
//int i_stren=0,  i_SelInt=-100;//, i_BHFnorm=0;
int i_subsh, i_FwBk;

int i_SelInt=-100; // Selects a particular interaction


int  i_sel_charge = -100;
bool sel_charge_flag = false;


int icut_ld=0, ntcut_ld=0;   //  icut_ld = 1;  ntcut_ld = 2;
int icut_Pi=0, ntcut_Pi=0;   //  icut_Pi = 1;  ntcut_Pi = 2;
int icut_fd=0, ntcut_fd=0;   //  icut_fd = 1;  ntcut_fd = 3;
int i_partSE=0,  i_AddSE=0;
int WelcheRec = 0;

// const int NMX_SET_Z_FWBK = 10;
// int SetNQp = 0, SetNQh = 0;
// int    SetNQp_i[NMX_SET_Z_FWBK], SetNQh_i[NMX_SET_Z_FWBK];
// double SetNQp_f[NMX_SET_Z_FWBK], SetNQh_f[NMX_SET_Z_FWBK];
// 
// int i_Pull = 0;
// const int NMX_PULL_ORB = 20;
// int n_PullOrb = 0;
// int    PullOrb_ish[NMX_PULL_ORB], PullOrb_n[NMX_PULL_ORB];
// double PullOrb_de[NMX_PULL_ORB];
// 


//
// Control of coupled cluster CCD iterations
//
double aLM_CCD = 0.5;
int    Stp_CCD = 1;


double Conv_check = SC_CHOP;

int iCheck = 0;

int i_gsp_rw = -1;
 
int i_npvts = 1;

int ntest = 0;


//-------------------
//  Test stuff:
int     i_hfbshift[2]={ 0 ,  0 };
double mu_hfbshift[2]={0.0, 0.0};
//-------------


//===================================================================
