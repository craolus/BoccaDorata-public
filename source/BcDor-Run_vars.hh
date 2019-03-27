//
//  
//
//

#include "BcDor-Global_variables.hh"

#define VERSION   "1.0"
#define SC_CHOP 1.e-7
#define  BUFFER1  100


extern VppInt_t  *Rrms;
extern ModSpace_t *MS;
extern VppInt_t  *Vpp;


//===================================================================
// run parameters:
extern int AddVc;                              extern char Vcpp_file[];
extern int RrmsCalc;                           extern char Rrmspp_file[];

extern int i_Ext_U1;                           extern char Ext_U1_file[];
extern double Ext_U1_Uo;
extern double Ext_U1_Uch[];

extern int IU1body; extern double U1body_fac;  extern char Trel_1_file[];
extern double a_Tcm;                           extern char Trel_2_file[];

                                               extern char Vpp_file[];
extern int iAddVpp, AddVppMx;                  extern char Vpp_add_file[][BUFFER1];
extern double AddVppMult[];
                                               extern char Vpp_out_file[];
                                               extern char SpProp_file[];
                                               extern char MdSp_file[];
                                               extern char MdSp_oslo_file[];
                                               extern char DysPivots_file[];
                                               extern char gout_file[];
                                               extern char gen_fout[];
extern int i_gout_str;                         extern char gout_str[];
                                               extern char Monopole_valence_file[];


                                               extern char temp_str[];

// Setting of Conventional/Unconventional baryon masses:
extern double mass_electron;
extern double mass_proton;
extern double mass_neutron;
extern double mass_particle_ave;
extern int i_set_new_mass;


extern int A, N, Z;
extern double bHO;
extern double hwHO;
extern int wAll;
extern bool b_plot_gsp;// = false;
extern int ItrMax, nDfw, nDbk;
extern int Lanczos; // N.B. must be zero, sothat Lanczos is not set
                 // with no optoion specified.
extern int i_2nd_ord, i_HF, i_ExtSE;
extern int i_MeanField, i_MeanField_FwBk;
//extern int i_stren,  i_SelInt;//, i_BHFnorm=0;
extern int i_subsh, i_FwBk;

extern int i_SelInt; // Selects a particular interaction


extern int  i_sel_charge;
extern bool sel_charge_flag;


extern int icut_ld, ntcut_ld;   //  icut_ld = 1;  ntcut_ld = 2;
extern int icut_Pi, ntcut_Pi;   //  icut_Pi = 1;  ntcut_Pi = 2;
extern int icut_fd, ntcut_fd;   //  icut_fd = 1;  ntcut_fd = 3;
extern int i_partSE,  i_AddSE;
extern int WelcheRec;

// extern const int NMX_SET_Z_FWBK;
// extern int SetNQp, SetNQh;
// extern int    SetNQp_i[], SetNQh_i[];  // declared to [NMX_SET_Z_FWBK]
// extern double SetNQp_f[], SetNQh_f[];  // declared to [NMX_SET_Z_FWBK]
// 
// extern int i_Pull;
// extern const int NMX_PULL_ORB;
// extern int n_PullOrb;
// extern int    PullOrb_ish[], PullOrb_n[];  // declared to [NMX_PULL_ORB]
// extern double PullOrb_de[];                // declared to [NMX_PULL_ORB]
//


//
// Control of coupled cluster CCD iterations
//
extern double aLM_CCD;
extern int    Stp_CCD;


extern double Conv_check; // = SC_CHOP;

extern int iCheck;

extern int i_gsp_rw;
 
extern int i_npvts;

extern int ntest;


//-------------------
//  Test stuff:
extern int     i_hfbshift[];
extern double mu_hfbshift[];
//-------------


//==============================================================================



void Second_order_diag(void );

void MGKoltun_sumrule(void );


void Load_Rrms_TwoBodyPart(ModSpace_t*, char*, VppInt_t*, char*);

double Calc_Rrms_TwoBody(SpProp_t*, VppInt_t*, bool plot=false);

double Calc_Rrms_TwoBody(SpProp_t*, char*, bool plot=false);

double Calc_Rrms_OneBody_rirj(SpProp_t*, VppInt_t*, double* COM_in=NULL, bool plot=false);

double Calc_Rrms_OneBody_rirj(SpProp_t*, char*, bool plot=false);


void Pull_Orbit(double[], int, int, ModSpace_t* );

void Retrive_SelfEn(Wide_self_energy_t*, int, ModSpace_t*);




//==============================================================================

//
//  File 'BcDor_Initialize.C'
//

void Set_defaults(void );

void get_part_number(void );

void Set_nucleus_data(void );

void Load_ModSp(void );

void Load_interaction(void );

void Load_SpProp(void);

void Unset_all(void);

//==============================================================================


//void verbose(void );

//void print_usage(void );

int parse_cmd_line(int , char** );

void Set_external_EFermi(int , double );

double Seek_EFermi(SpProp_t* , int , int iplot=1);

bool EFermi_Is_Set(int );



//==============================================================================

//
//
// Headres for the Lanczos and Arnoldi routines (Krylov reductions)...
//

void Make_LancPivots(ModSpace_t*, int*, double**, int[], int[], int, int iplot=0);
//void Make_LancPivots(SpProp_t *gin, int*, double**, int[], int[], int); // an OLD version that used gin too


int Reduce_Mtx_Lncz(int, int, int, double *, int*,
					int, int, int, double *,
					int, double**, int* );


//============================================================


//==============================================================================
