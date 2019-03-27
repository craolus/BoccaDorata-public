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
  NdataXY = NDATAPIXY;
  return;}

void Pol_prop_t::MVG_init_NULL(void ) {
  MdSp = NULL;
  vpp  = NULL;
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



void Pol_prop_t::free_bases(void ) {
  free_2qp_basis();
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
  return; }

void Pol_prop_t::XY_init_NULL(void ) {
  XYmtx = NULL;
  XY_nalloc = -1;
  XY_nvects = -1;
  return; }



void Pol_prop_t::free_matrices(void ) {
  Clear_XYmtx();
  return;}

void Pol_prop_t::Clear_XYmtx(void ) {
  if (NULL != XYmtx) delete [] XYmtx;
  XY_init_NULL();
  return; }


int Pol_prop_t::allocate_XY(int n_exp_sol, int n_bas_vec ) {
  Clear_XYmtx(); // clear memory if something was previously allocated
  XY_nalloc = n_exp_sol;
  XY_nvects = n_bas_vec + NdataXY;
  XYmtx = new double[(XY_nvects*XY_nalloc)+NdataXY+1];
  //double *ptr2 = XYmtx + (XY_nvects*XY_nalloc);
  //for(double *ptr1=XYmtx; ptr1<ptr2; ++ptr1) (*ptr1)=0.0; //clean up array
  return 0;}

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

