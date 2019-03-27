//
//
//   Routines to calculate the rms radii corrected from the 
//   center of mass motion
//
//
//  C.Barbieri, Surrey, December 2011.
//


#include <iostream>
//#include <fstream>
#include <iomanip> 
using namespace std;

#include <math.h>

#include "BcDor-Global_variables.hh"
#include "BcDor-Run_vars.hh"
#include "BcDor-PhysConsts.hh"


void Load_Rrms_TwoBodyPart(ModSpace_t *MdSp_in, char *filein, VppInt_t *Rrms_2b, char* filetype) {
  //
  //  Loads the (R_ij)^2 matrix elements from file to be used for the two-body-only
  // formula and scales them accordingly...
  // 
  //
  cout << "\n\n Loading the matrix elements of "<<filetype<<" from file '"<<filein <<"'...\n"<<flush;
  Rrms_2b->Clear_all(MdSp_in);
  Rrms_2b->read(filein);
  //
  cout << "...will use  A=="<<A<<" for calculating the c.o.m. corrections.\n\n";
  double x1 = MdSp_in->bho/double(A); x1*=x1;
  Rrms_2b->scale(x1);
  return;}


double Calc_Rrms_TwoBody(SpProp_t *gin, VppInt_t *Rrms_2b, bool plot/*=false*/) {
  //
  //  This simply claculates the rms radius at 1st order based on the input
  // propagator and (R_ij)^2 matrix elements, both given as input.
  //
  //  Here Rrms_2b  is supposed to be already scaled by (bHO/A)^2
  //
  double x1 = gin->Compute_FirstOrderExpValue(Rrms_2b);
  x1 = sqrt(x1);
  if (plot) cout << "\n r_rms (1st order--2b)= " << x1 <<"\n";
  return x1;
}

double Calc_Rrms_TwoBody(SpProp_t *gin, char *filein, bool plot/*=false*/) {
  //
  //  Calculates R_rms w/ the two-body-only formula, prints the results,
  // gives it as output and delete the Rrms_2b matrix elemts from memory.
  //
  VppInt_t Rrms_2b(gin->MdSp);
  //
  Load_Rrms_TwoBodyPart(gin->MdSp, filein, &Rrms_2b, "(r_ij)^2");
  //
  double x1 = Calc_Rrms_TwoBody(gin, &Rrms_2b, plot);
  return x1;}


double Calc_Rrms_OneBody_rirj(SpProp_t *gin, VppInt_t *Rrms_1b, double* COM_in, bool plot/*=false*/) {
  //
  //  This simply claculates the rms radius at 1st order based on the input
  // propagator the (r_i)^2 exp. vale and the (r_i.r_j)  matrix elements, both given as input.
  //
  //  Here Rrms_1b  is supposed to be already scaled by (bHO/A)^2  (note that the file mtx els
  //                                                               must contsin a -2.0 factor!)
  //
  int ish_a,nsh_a,nsh_b, ka;
  double rsq_1b, rsq_2b, r_rms, x1;
  
  //
  //  Two-body contribution:
  rsq_2b = gin->Compute_FirstOrderExpValue(Rrms_1b);
 
  
  
  rsq_1b = 0.0;
  for(ka=0; ka<gin->n_qh; ++ka) {
    ish_a = gin->n_clj_h[ka];
    x1 = 0.0;
    for(nsh_a=0; nsh_a<gin->MdSp->MSp_no[ish_a]; ++nsh_a) {
      for(nsh_b=0; nsh_b<gin->MdSp->MSp_no[ish_a]; ++nsh_b) {
        x1 += gin->fqh[ka*gin->N_SPB_MX+nsh_a] * gin->fqh[ka*gin->N_SPB_MX+nsh_b] *
          fab_rsq(gin->MdSp->MSp_n[ish_a][nsh_a] , gin->MdSp->MSp_n[ish_a][nsh_b] , gin->MdSp->MSp_l[ish_a]);
      }
    }
    //cout << "\n r.m.s. radius of quasihole ka="<<ka<<" (eh="<<gin->eh[ka]<<",  tlj="<<gin->MdSp->MSp_name[ish_a]
    //             <<") is: "<< sqrt(x1)*(gin->MdSp->bho) <<"  (no c.o.m. correction!)\n";
    
    x1 *= double(gin->MdSp->MSp_2j[ish_a] +1);
    if (NULL != COM_in) x1 *= COM_in[gin->MdSp->MSp_ch[ish_a]];
    rsq_1b += x1;
  }
  
  x1 = gin->MdSp->bho/double(A);
  x1 *= x1*(double(A) - 1.0);
  rsq_1b *= x1;
  x1 = sqrt(rsq_1b + rsq_2b);
  
  if (plot) cout << "\n r_rms, point matter (1st order--rirj)= " << x1 <<"\n";
  if (plot) cout << "\n r_rms, point matter (1-body with COM, no 2-body )= " << sqrt(rsq_1b) <<"\n";
  if (plot) cout << "\n r_rms, point matter (1-body pure, no COM )= " << sqrt(rsq_1b*double(A)/(double(A) - 1.0)) <<"\n";
  return x1;}


double Calc_Rrms_OneBody_rirj(SpProp_t *gin, char* filein, bool plot/*=false*/) {
  //
  //  Calculates R_rms w/ the c.o.m corrected one-body-plus-two-boty formula,
  // prints the results, gives it as output and delete the Rrms_1b matrix elemts
  // from memory.
  //
  //  It also calculates the point proton and neutron radii and the charge radii.
  //
  double rA, rZ, rN, rToT;
  //
  VppInt_t Rrms_1b(gin->MdSp);
  //
  Load_Rrms_TwoBodyPart(gin->MdSp, filein, &Rrms_1b, "(r_i . r_j)^2");
  //
  //COM_factor[0] = 1.0;
  //COM_factor[1] = 1.0;
  //rA = Calc_Rrms_OneBody_rirj(gin, &Rrms_1b, COM_factor, plot);
  rA = Calc_Rrms_OneBody_rirj(gin, &Rrms_1b, NULL, plot); // NULL means no extra 1-body factors
  //
  //
  //  The calculation of proton and neutron radii proceeds by
  // taking the pure matrix elements and then rescaling them
  // appropriately. Also the 1-body sums are to be rescaled. After
  // this the same formula for the matter radii yields the proton
  // and neutron point radii.
  //
  //  Rather than reloading the same Rrms_1b matrix elements three times,
  // we first calculate protons by rescalin the original and then calculate
  // the neutons by inverting the transformation made for the protons
  //

  for(int ipw=0; ipw=gin->MdSp->nsubsh; ++ipw)
    if ( (0 > gin->MdSp->MSp_ch[ipw]) || (1 < gin->MdSp->MSp_ch[ipw]) ){
      cerr << "\n\n WARNING the model space seems to have states with charge different from protons and"
      << "    neutrons. Therefore I am olnly givng results for the point matter radii and will skip"
      << "    calculatinr neutron and charge radii...\n\n";
      return rA;
    }

  double COM_factor[2]= {0.0, 0.0};

  // Charge radii:
  //Load_Rrms_TwoBodyPart(gin->MdSp, filein, &Rrms_1b, "(r_i . r_j)^2");
  // Rescale matrix elements to calculate point-proton radii:
  Rrms_1b.scale_chrg(       -1.0              , 0);
  Rrms_1b.scale_chrg(double( A - Z )/double(Z), 1);
  Rrms_1b.scale_chrg(double(2*A - Z)/double(Z), 2);
  //
  COM_factor[0] =            1.0                   / double(A - 1);
  COM_factor[1] = double( A*(A-2) + Z) / double(Z) / double(A - 1);
  rZ = Calc_Rrms_OneBody_rirj(gin, &Rrms_1b, COM_factor, false);
  rToT = rZ*rZ * double(Z) / double(A); // This is for test, r_p^2 + r_n^2 must reproduce r_matter
  if (plot) cout << "\n r_rms, point proton (1st order--rirj)= " << rZ <<"\n";
  //
  // Apply corrections to obtain the charge radius:
  rZ *= rZ;
  rZ += PROTON_CH_RADIUS*PROTON_CH_RADIUS;   // proton ch raius
  rZ -= NEUTRON_ABS_CH_RADIUS*NEUTRON_ABS_CH_RADIUS*double(N)/double(Z);  // (neutron ch radius)^2 is negative!
  rZ += 0.75*pow( (MeVfm/PROTONmass), 2);        //  Darwin-Foldy term
  rZ = sqrt(rZ);
  if (plot) cout << "\n r_rms, charge (1st order--rirj)= " << rZ <<"\n";
  //
  //
  // Point-neutron radii:
  //Load_Rrms_TwoBodyPart(gin->MdSp, filein, &Rrms_1b, "(r_i . r_j)^2");
  //Rrms_1b.scale_chrg(double(2*A - N)/double(N), 0);
  //Rrms_1b.scale_chrg(double( A - N )/double(N), 1);
  //Rrms_1b.scale_chrg(       -1.0              , 2);
  Rrms_1b.scale_chrg(double(N - 2*A)/double(N)              , 0);
  Rrms_1b.scale_chrg(double( (A - N)*Z )/double( (A - Z)*N ), 1);
  Rrms_1b.scale_chrg(double(Z)/double(Z - 2*A)              , 2);
  //
  COM_factor[0] = double( A*(A-2) + N) / double(N) / double(A - 1);
  COM_factor[1] =            1.0                   / double(A - 1);
  rN = Calc_Rrms_OneBody_rirj(gin, &Rrms_1b, COM_factor, false);
  if (plot) cout << "\n r_rms , point neutron (1st order--rirj)= " << rN <<"\n";
  //
  // This is for test, r_p^2 + r_n^2 must reproduce r_matter
  rToT += rN*rN * double(N) / double(A);
  if (plot) cout << "\n r_rms, point neutron + point proton (1st order--rirj)= " << sqrt(rToT) <<"\n";
  //
  return rA;}

