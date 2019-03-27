//
//
//   Build the pp/hh propagator based on the successive approximations
//  to the Faynman-Galitski equation (FGE):
//     (D)RPA,  2Pi-RPA, ...
//
// gII_a.C: Initialization and utility functions
//
// C.Barbieri, GSI-RIKEN-Surrey (2006-2015)  (C.Barbieri@surrey.ac.uk)
//

#include <iostream>
#include <fstream>
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
#include "BcDor-gII_classes.hh"


//
// Different types of constructors:
//

gII_prop_t::gII_prop_t() {
  //cout << "\n\n Creating a Two-Particle--Two-Holes Propagator...    ";
  Ndata_init();
  MVG_init_NULL();
  bases_init_NULL();
  matrices_init_NULL();
  //cout << "  ...no additional information provided.\n\n";
  return; }

gII_prop_t::gII_prop_t(ModSpace_t *MdSpin, SpProp_t *gin, VppInt_t *vin) {
  //cout << "\n\n Creating a Two-Particle--Two-Holes Propagator...    ";
  Ndata_init();
  MVG_init_NULL();
  bases_init_NULL();
  matrices_init_NULL();
  set_space(MdSpin, gin, vin);
  //cout << "  ...mod.sp., sp.prop. and vpp int. were provided.\n\n";
  return; }

gII_prop_t::gII_prop_t(ModSpace_t *MdSpin, SpProp_t *gin) {
  //cout << "\n\n Creating a Two-Particle--Two-Holes Propagator...    ";
  Ndata_init();
  MVG_init_NULL();
  bases_init_NULL();
  matrices_init_NULL();
  set_space(MdSpin, gin);
  //cout << "  ...mod.sp. and sp.prop. were provided.\n\n";
  return; }

gII_prop_t::gII_prop_t(SpProp_t *gin) {
  //cout << "\n\n Creating a Two-Particle--Two-Holes Propagator...    ";
  Ndata_init();
  MVG_init_NULL();
  bases_init_NULL();
  matrices_init_NULL();
  set_space(gin);
  //cout << "  ...the sp.prop. was provided.\n\n";
  return; }

gII_prop_t::gII_prop_t(ModSpace_t *MdSpin) {
  //cout << "\n\n Creating a Two-Particle--Two-Holes Propagator...    ";
  Ndata_init();
  MVG_init_NULL();
  bases_init_NULL();
  matrices_init_NULL();
  set_MdSp(MdSpin);
  //cout << "  ...the mod.sp. was provided.\n\n";
  return; }

gII_prop_t::gII_prop_t(VppInt_t *vin) {
  //cout << "\n\n Creating a Two-Particle--Two-Holes Propagator...    ";
  Ndata_init();
  MVG_init_NULL();
  bases_init_NULL();
  matrices_init_NULL();
  set_ppint(vin);
  //cout << "  ...the vpp int. was provided.\n\n";
  return; }

gII_prop_t::~gII_prop_t() {
  //cout << "\n\n Deallocating a Two-Particle--Two-Holes Propagator. \n\n";
  //cout << "Deallocating a Two-Particle--Two-Holes Propagator.\n";
  //TST//cout << "Deallocating a Two-Particle--Two-Holes Propagator."<<Jf<<","<<pif<<","<<dT<<"\n";
  free_bases();
  free_matrices();
  bases_init_NULL();
  matrices_init_NULL();
  return; }

void gII_prop_t::Ndata_init(void ) {
  NdataXY  = NDATALDXY;
  return;}

void gII_prop_t::MVG_init_NULL(void ) {
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
void gII_prop_t::bases_init_NULL(void ) {
  basis_2qp_init_NULL();
  return; }

void gII_prop_t::basis_2qp_init_NULL(void ) {
  e_pp  = NULL;
  n1_pp = NULL;
  n2_pp = NULL;
  n_pp_alloc = -1;
  n_pp_basis = -1;
  e_hh  = NULL;
  k1_hh = NULL;
  k2_hh = NULL;
  n_hh_alloc = -1;
  n_hh_basis = -1;
  return; }


void gII_prop_t::free_bases(void ) {
  free_2qp_basis();
  return; }

void gII_prop_t::free_2qp_basis(void ) { // (int i=0) {
  if (NULL != e_pp ) delete [] e_pp ;
  if (NULL != n1_pp) delete [] n1_pp;
  if (NULL != n2_pp) delete [] n2_pp;
  if (NULL != e_hh ) delete [] e_hh ;
  if (NULL != k1_hh) delete [] k1_hh;
  if (NULL != k2_hh) delete [] k2_hh;
  basis_2qp_init_NULL();
  // if (i) Clear_XYmtx();
  return; }


int gII_prop_t::allocate_2qp_basis(int ipp, int ihh) {
  if ((n_pp_alloc > 0) || (n_hh_alloc > 0)) free_2qp_basis();
  if (ipp > 0) {
    n_pp_alloc = ipp;
    e_pp  = new double[n_pp_alloc];
    n1_pp = new int[n_pp_alloc];
    n2_pp = new int[n_pp_alloc];
    n_pp_basis = -1;
    ipp = 0;
    } else {ipp = 1;}
  if (ihh > 0) {
    n_hh_alloc = ihh;
    e_hh  = new double[n_hh_alloc];
    k1_hh = new int[n_hh_alloc];
    k2_hh = new int[n_hh_alloc];
    n_hh_basis = -1;
    ihh = 0;
    } else {ihh = 1;}
  return ipp*ihh; } // ==1 if nothing is allocated

//------------------------------------------------------

//------------------------------------------------------
//  Lehmann representations of the two-body propagator
// in XY, santdard(Zab) and int. box (Dab) forms:

// NOTE: the 'xxx_init_NULL()' functions are to set the poiters to NULL values
//  when no memory is allocated. These should ONLY be called only when no
//  vectors are yet allocated (i.e. at the beginning of the constructor) and
//  right after cleaning up all the allocated memory--in 'free_xxx()'.
void gII_prop_t::matrices_init_NULL(void ) { 
  XY_init_NULL();
  return; }

void gII_prop_t::XY_init_NULL(void ) {
  XYmtx = NULL;
  XY_nalloc = -1;
  XY_nvects = -1;
  return; }


void gII_prop_t::free_matrices(void ) {
  Clear_XYmtx();
  return;}

void gII_prop_t::Clear_XYmtx(void ) {
  if (NULL != XYmtx) delete [] XYmtx;
  XY_init_NULL();
  return; }


int gII_prop_t::allocate_XY(int n_exp_sol, int n_bas_vec ) {
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

void gII_prop_t::set_space(ModSpace_t *MdSpin, SpProp_t *gin, VppInt_t *vin) {
  set_MdSp(MdSpin);  // this before bacause it may reset everything (gsp included!!!!)
  gsp  = gin;
  if (MdSp != gsp->MdSp) {cerr << "\n\n WARNING (gII_prop_t::set_space): ambiguity in the Mod. space (of gsp).... \n\n";}// cin;}
  set_ppint(vin);
  return; }

void gII_prop_t::set_space(ModSpace_t *MdSpin, SpProp_t *gin) {
  set_MdSp(MdSpin);  // this before bacause it may reset everything (gsp included!!!!)
  gsp  = gin;
  if (MdSp != gsp->MdSp) {cerr << "\n\n WARNING (gII_prop_t::set_space): ambiguity in the Mod. space (of gsp).... \n\n";}// cin;}
  return; }

void gII_prop_t::set_space(SpProp_t *gin) {
  set_MdSp(gin->MdSp);  // this before bacause it may reset everything (gsp included!!!!)
  gsp  = gin;
  return; }

void gII_prop_t::set_MdSp(ModSpace_t *MdSpin) {
  if (MdSp != MdSpin) {
     // reset all:
     free_bases();
     //free_matrices();  ??
     MVG_init_NULL();
     }
  MdSp = MdSpin;
  return; }

void gII_prop_t::set_ppint(VppInt_t *vin) {
  vpp  = vin;
  if (MdSp != vpp->MdSp) {cerr << "\n\n WARNING (gII_prop_t::set_ppint): ambiguity in the Mod. space (of vpp)... \n\n";}// cin;}
  return; }
