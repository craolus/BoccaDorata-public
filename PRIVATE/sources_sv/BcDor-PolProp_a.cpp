//
//
//   Build the Polarozation propagator based on the successive approximation
//  to the Bethe-Salpeter eqaution (BSE):
//     (D)RPA,  2Pi-RPA, ...
//
// PolProp_a.C: Initialization and utility functions
//
// C.Barbieri, GSI, December 2006.  (C.Barbieri@gsi.de)
//

#include <iostream>
#include <fstream>
//#include <iomanip>
using namespace std;

//#include <cstdlib>
//#include <cstdio>
//#include <cmath>

//#define DPI   3.14159265358979323846
//#define SQR2  1.414213562373095048801688


//#include "BcDor-Ang_momenta.hh"
//#include "BcDor-HO_radial_me.hh"

#include "BcDor-Global_variables.hh"
#include "BcDor-Sp_classes.hh"
#include "BcDor-Int_classes.hh"
#include "BcDor-PolProp_classes.hh"
//#include "Utilities.hh"

//
// Constructors:
//

Pol_prop_t::Pol_prop_t() {
  //cout << "\n\n Allocating a Polarization Propagator...    ";
  Ndata_init();
  MVG_init_NULL();
  bases_init_NULL();
  matrices_init_NULL();
  //cout << "  ...no additional information provided.\n\n";
  return; }

Pol_prop_t::Pol_prop_t(ModSpace_t *MdSpin, SpProp_t *gin, VppInt_t *vin) {
  //cout << "\n\n Allocating a Polarization Propagator...    ";
  Ndata_init();
  MVG_init_NULL();
  bases_init_NULL();
  matrices_init_NULL();
  set_space(MdSpin, gin, vin);
  //cout << "  ...mod.sp., sp.prop. and vpp int. were provided.\n\n";
  return; }

Pol_prop_t::Pol_prop_t(ModSpace_t *MdSpin, SpProp_t *gin) {
  //cout << "\n\n Allocating a Polarization Propagator...    ";
  Ndata_init();
  MVG_init_NULL();
  bases_init_NULL();
  matrices_init_NULL();
  set_space(MdSpin, gin);
  //cout << "  ...mod.sp. and sp.prop. were provided.\n\n";
  return; }

Pol_prop_t::Pol_prop_t(SpProp_t *gin) {
  //cout << "\n\n Allocating a Polarization Propagator...    ";
  Ndata_init();
  MVG_init_NULL();
  bases_init_NULL();
  matrices_init_NULL();
  set_space(gin);
  //cout << "  ...the sp.prop. was provided.\n\n";
  return; }

Pol_prop_t::Pol_prop_t(ModSpace_t *MdSpin) {
  //cout << "\n\n Allocating a Polarization Propagator...    ";
  Ndata_init();
  MVG_init_NULL();
  bases_init_NULL();
  matrices_init_NULL();
  set_MdSp(MdSpin);
  //cout << "  ...the mod.sp. was provided.\n\n";
  return; }

Pol_prop_t::Pol_prop_t(VppInt_t *vin) {
  //cout << "\n\n Allocating a Polarization Propagator...    ";
  Ndata_init();
  MVG_init_NULL();
  bases_init_NULL();
  matrices_init_NULL();
  set_ppint(vin);
  //cout << "  ...the vpp int. was provided.\n\n";
  return; }

Pol_prop_t::~Pol_prop_t() {
  //cout << "\n\n Deallocating a Polarization Propagator. \n\n";
  //cout << "Deallocating a Polarization Propagator.\n";
  free_bases();
  free_matrices();
  bases_init_NULL();
  matrices_init_NULL();
  return; }

void Pol_prop_t::Ndata_init(void ) {
//XDIST//  sprintf(STOREDIR, BcDorWorkFolder);
  NdataXY = NDATAPIXY;
//XDIST//  NdataZ  = NDATAPIZAB;
//XDIST//  NdataD  = NDATAPIDAB;
  return;}

void Pol_prop_t::MVG_init_NULL(void ) {
  MdSp = NULL;
  vpp  = NULL;
  //vph  = NULL;
  gsp  = NULL;
  return; }


//------------------------------------------------------
//  Bases for the 2qp-2qh space and the antisymmetrized
// two-fermion mod. space.

// NOTE: the 'xxx_init_NULL()' functions are to set the poiters to NULL values
//  when no memory is allocated. These should ONLY be called only when no
//  vectors are yet allocated (i.e. at the beginning of the constructor) and
//  right after cleaning up all the allocated memory--in 'free_xxx()'.
void Pol_prop_t::bases_init_NULL(void ) {
  basis_2qp_init_NULL();
//XDIST//  basis_2bmsp_init_NULL();
  return; }

void Pol_prop_t::basis_2qp_init_NULL(void ) {
  e_ph = NULL;
  n_ph = NULL;
  k_ph = NULL;
  n_ph_alloc = -1;
  n_ph_basis = -1;
  e_hp = NULL;
  n_hp = NULL;
  k_hp = NULL;
  n_hp_alloc = -1;
  n_hp_basis = -1;
  return; }

//XDIST//
/*
void Pol_prop_t::basis_2bmsp_init_NULL(void ) {
  ia_2b = NULL;
  ib_2b = NULL;
  n_2bmsp = -1;
  N_2BMSP_alloc = -1;
  return; }
*/


void Pol_prop_t::free_bases(void ) {
  free_2qp_basis();
//XDIST//  free_2bmsp_basis();
  return; }

void Pol_prop_t::free_2qp_basis(void ) { // (int i=0) {
  if (NULL != e_ph) delete [] e_ph;
  if (NULL != n_ph) delete [] n_ph;
  if (NULL != k_ph) delete [] k_ph;
  if (NULL != e_hp) delete [] e_hp;
  if (NULL != n_hp) delete [] n_hp;
  if (NULL != k_hp) delete [] k_hp;
  basis_2qp_init_NULL();
  // if (i) Clear_XYmtx();
  return; }

//XDIST//
/*
void Pol_prop_t::free_2bmsp_basis(void ) {// (int i=0) {
  if (NULL != ia_2b) delete [] ia_2b;
  if (NULL != ib_2b) delete [] ib_2b;
  basis_2bmsp_init_NULL();
  // if (i) Clear_Zabmtx();
  // if (i) Clear_Dabmtx();
  return; }
*/


int Pol_prop_t::allocate_2qp_basis(int iph, int ihp) {
  if ((n_ph_alloc > 0) || (n_hp_alloc > 0)) free_2qp_basis();
  if (iph > 0) {
    n_ph_alloc = iph;
    e_ph = new double[n_ph_alloc];
    n_ph = new int[n_ph_alloc];
    k_ph = new int[n_ph_alloc];
    n_ph_basis = -1;
    iph = 0;
    } else {iph = 1;}
  if (ihp > 0) {
    n_hp_alloc = ihp;
    e_hp = new double[n_hp_alloc];
    n_hp = new int[n_hp_alloc];
    k_hp = new int[n_hp_alloc];
    n_hp_basis = -1;
    ihp = 0;
    } else {ihp = 1;}
  return (iph*ihp); } // ==1 if nothing was allocated

//XDIST//
/*
int Pol_prop_t::allocate_2bmsp_basis(int iab) {
  if (N_2BMSP_alloc > 0) free_2bmsp_basis();
  if (iab > 0) {
    N_2BMSP_alloc = iab;
    ia_2b = new int[N_2BMSP_alloc];
    ib_2b = new int[N_2BMSP_alloc];
    n_2bmsp = -1;
    return 0;
  }
  return 1;} // ==1 if they were not allocated
*/
//------------------------------------------------------

//------------------------------------------------------
//  Lehmann representations of the two-body propagator
// in XY, santdard(Zab) and int. box (Dab) forms:

// NOTE: the 'xxx_init_NULL()' functions are to set the poiters to NULL values
//  when no memory is allocated. These should ONLY be called only when no
//  vectors are yet allocated (i.e. at the beginning of the constructor) and
//  right after cleaning up all the allocated memory--in 'free_xxx()'.
void Pol_prop_t::matrices_init_NULL(void ) { 
  XY_init_NULL();
//XDIST//  Zab_init_NULL();
//XDIST//  Dab_init_NULL();
//XDIST//  FaddBx_Fw_init_NULL();
//XDIST//  FaddBx_Bk_init_NULL();
  return; }

void Pol_prop_t::XY_init_NULL(void ) {
  XYmtx = NULL;
  XY_nalloc = -1;
  XY_nvects = -1;
  return; }

//XDIST//
/*
void Pol_prop_t::Zab_init_NULL(void ) {
  Zabmtx = NULL;
  Zab_nalloc = -1;
  Zab_nvects = -1;
  return; }

void Pol_prop_t::Dab_init_NULL(void ) {
  Dabmtx = NULL;
  Dab_nalloc = -1;
  Dab_nvects = -1;
  return; }

void Pol_prop_t::FaddBx_Fw_init_NULL(void ) {
  FaddBx_f1 = NULL;
  FaddBx_f2 = NULL;
  return;}

void Pol_prop_t::FaddBx_Bk_init_NULL(void ) {
  FaddBx_b1 = NULL;
  FaddBx_b2 = NULL;
  return;}
*/


void Pol_prop_t::free_matrices(void ) {
  Clear_XYmtx();
//XDIST//  Clear_Zabmtx();
//XDIST//  Clear_Dabmtx();
//XDIST//  Clear_FaddBx_Fw();
//XDIST//  Clear_FaddBx_Bk();
  return;}

void Pol_prop_t::Clear_XYmtx(void ) {
  if (NULL != XYmtx) delete [] XYmtx;
  XY_init_NULL();
  return; }

//XDIST//
/*
void Pol_prop_t::Clear_Zabmtx(void ) {
  if (NULL != Zabmtx) delete [] Zabmtx;
  Zab_init_NULL();
  return; }

void Pol_prop_t::Clear_Dabmtx(void ) {
  if (NULL != Dabmtx) delete [] Dabmtx;
  Dab_init_NULL();
  return; }

void Pol_prop_t::Clear_FaddBx_Fw(void ) {
  if (NULL != FaddBx_f1) delete [] FaddBx_f1;
  if (NULL != FaddBx_f2) delete [] FaddBx_f2;
  FaddBx_Fw_init_NULL();
  return;}
  
void Pol_prop_t::Clear_FaddBx_Bk(void ) {
  if (NULL != FaddBx_b1) delete [] FaddBx_b1;
  if (NULL != FaddBx_b2) delete [] FaddBx_b2;
  FaddBx_Bk_init_NULL();
  return;}
*/

int Pol_prop_t::allocate_XY(int n_exp_sol, int n_bas_vec ) {
  Clear_XYmtx(); // clear memory if something was previously allocated
  XY_nalloc = n_exp_sol;
  XY_nvects = n_bas_vec + NdataXY;
  XYmtx = new double[(XY_nvects*XY_nalloc)+NdataXY+1];
  //double *ptr2 = XYmtx + (XY_nvects*XY_nalloc);
  //for(double *ptr1=XYmtx; ptr1<ptr2; ++ptr1) (*ptr1)=0.0; //clean up array
  return 0;}

//XDIST//
/*
int Pol_prop_t::allocate_Zab(void ) {
  if (n_2bmsp    < 1) return (1 -n_2bmsp);    // tells why nothing has
  if (n_tot_sols < 1) return (1001-n_tot_sols); // been allocated
  Clear_Zabmtx();
  Zab_nalloc = n_tot_sols;
  Zab_nvects = n_2bmsp + NdataZ;
  Zabmtx = new double[(Zab_nvects*Zab_nalloc)+NdataZ+1];
  return 0;}

int Pol_prop_t::allocate_Dab(void ) {
  if (n_2bmsp    < 1) return (1 -n_2bmsp);    // tells why nothing has
  if (n_tot_sols < 1) return (1001-n_tot_sols); // been allocated
  Clear_Dabmtx();
  Dab_nalloc = n_tot_sols;
  Dab_nvects = n_2bmsp + NdataD;
  Dabmtx = new double[(Dab_nvects*Dab_nalloc)+NdataD+1];
  return 0;}

int Pol_prop_t::allocate_FaddBx_Fw(void ) {
  if (n_fw_sols > 0) {
    Clear_FaddBx_Fw();
    FaddBx_f1 = new double[n_fw_sols*n_fw_sols];
    FaddBx_f2 = new double[n_fw_sols*n_fw_sols];
    return 0;
  }
  return 1;}

int Pol_prop_t::allocate_FaddBx_Bk(void ) {
  if (n_bk_sols > 0) {
    Clear_FaddBx_Bk();
    FaddBx_b1 = new double[n_bk_sols*n_bk_sols];
    FaddBx_b2 = new double[n_bk_sols*n_bk_sols];
    return 0;
  }
  return 1;}
*/
//------------------------------------------------------

//------------------------------------------------------
//  Set the associated mod. space, s.p. propagator, 
// interaction, etc...

void Pol_prop_t::set_space(ModSpace_t *MdSpin, SpProp_t *gin, VppInt_t *vin) {
  set_MdSp(MdSpin);  // this before bacause it may reset everything (gsp included!!!!)
  gsp  = gin;
  if (MdSp != gsp->MdSp) {cout << "\n\n WARNING (Pol_prop_t::set_space): ambiguity in the Mod. space (of gsp).... \n\n";}// cin;}
  set_ppint(vin);
  return; }

void Pol_prop_t::set_space(ModSpace_t *MdSpin, SpProp_t *gin) {
  set_MdSp(MdSpin);
  gsp  = gin;
  if (MdSp != gsp->MdSp) {cout << "\n\n WARNING (Pol_prop_t::set_space): ambiguity in the Mod. space (of gsp).... \n\n";}// cin;}
  return; }

void Pol_prop_t::set_space(SpProp_t *gin) {
  set_MdSp(gin->MdSp);  // this before bacause it may reset everything (gsp included!!!!)
  gsp  = gin;
  return; }

void Pol_prop_t::set_MdSp(ModSpace_t *MdSpin) {
  if (MdSp != MdSpin) {
     // reset all
     free_bases();
     //free_matrices();  ??
     MVG_init_NULL();
     }
  MdSp = MdSpin;
  return; }

void Pol_prop_t::set_ppint(VppInt_t *vin) {
  vpp  = vin;
  if (MdSp != vpp->MdSp) {cout << "\n\n WARNING (Pol_prop_t::set_ppint): ambiguity in the Mod. space (of vpp)... \n\n";}// cin;}
  return; }

//void Pol_prop_t::set_phint(VphInt_t *vin) {
//  vph  = vin;
//  if (MdSp != vph->MdSp) {cout << "\n\n WARNING (Pol_prop_t::set_phint): ambiguity in the Mod. space (of vph)... \n\n"; cin;}
//  return; }
