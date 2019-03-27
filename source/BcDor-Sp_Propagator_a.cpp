//
//   The class "SpProp_t" stores the one-body Green's function (or
//  sp propagator) of a finite many-particle sytem.
//
//   The Lehmann representation is used in this implementation, so sp
//  energies and (the components in the sp basis of the) overlap wave
//  funtcions are stored.
//   A sp model space (of type ModSpace_t)nmust be associated with each
//  object, to provide a discrete basis to expend the overla amplitudes.
//
//   Here, the case of a spherical system is implemented (J_tot == 0).
//
//
//
//   The file 'Sp_Proppagator_a.C' contains routines for generating the
//  HF and BHF potentials
//
//
// C.Barbieri, RIKEN, April 2010.  (barbieri@riken.jp).
//
//

#include <cstdlib>
// #include <cstdio>
// #include <cmath>

#include <iostream>
using namespace std;

#include "BcDor-Sp_classes.hh"
#include "BcDor-Int_classes.hh"


#define SQR2 1.4142135623730950488016887




int SpProp_t::ExtendedHartreeFock_Pot(double *Uab_EHF, int NDIM, VppInt_t *Vppin, int I_SH) {

  double x1,x2,x3,x6,x7;
  double U_ehf;

  int    J2,J2u,J2l,nJ,indx;
  int    ia,im,ib,iu,nmsh,nvsh;
  int    iv_clj;

  int    Norb_mx;

  double *qh_ptr;


  if (MdSp != Vppin->MdSp) {
     cerr << "SpProp_t::ExtendedHartreeFock_Pot, the model space of the"
          << " interaction do not  correspond to the one of the propagator.\n\n";
     return 100;
     }
                       

  Norb_mx  = MdSp->MSp_no[I_SH];

  for(int nc=0; (nc<Norb_mx)&&(nc<NDIM); nc++) {
   for(int nr=nc; (nr<Norb_mx)&&(nr<NDIM); nr++) {

     //
     // Compute the generalized HF diagram:
     //
     U_ehf = 0.0;
     for(iv_clj=0; iv_clj<MdSp->nsubsh ; ++iv_clj) {

      J2l = abs(MdSp->MSp_2j[I_SH]-MdSp->MSp_2j[iv_clj]);
      J2u =     MdSp->MSp_2j[I_SH]+MdSp->MSp_2j[iv_clj] ;
      for(J2=J2l; J2<=J2u; J2+=2) {
      
        nJ = (MdSp->MSp_2j[I_SH]+MdSp->MSp_2j[iv_clj]+J2)/2;
        nJ = nJ%2;

        x1 = double(J2+1) / double(MdSp->MSp_2j[I_SH]+1);

        x7 = 0.0;
        for(nmsh=0; nmsh< MdSp->MSp_no[iv_clj]; nmsh++) {
        if ((iv_clj == I_SH) && (nmsh == nr) && (nJ != 1)) continue;
        
         x2 = 1.0;
         if ((iv_clj == I_SH) && (nmsh == nr)) x2 = SQR2;

         ia = (MdSp->MSp_n[I_SH]   - MdSp->Mor_n) + nr;
         im = (MdSp->MSp_n[iv_clj] - MdSp->Mor_n) + nmsh;

        for(nvsh=0; nvsh< MdSp->MSp_no[iv_clj]; nvsh++) {
        if ((iv_clj == I_SH) && (nvsh == nc) && (nJ != 1)) continue;
        
         x3 = 1.0;
         if ((iv_clj == I_SH) && (nvsh == nc)) x3 = SQR2;

         ib = (MdSp->MSp_n[I_SH]   - MdSp->Mor_n) + nc;
         iu = (MdSp->MSp_n[iv_clj] - MdSp->Mor_n) + nvsh;

         x6 = 0.0;
         qh_ptr = Sp_clj_fqh[iv_clj];
         for(int kg=0; kg<Sp_clj_nh[iv_clj]; ++kg) {
                              x6 += qh_ptr[nvsh]*qh_ptr[nmsh];
                              qh_ptr += N_SPB_MX;
                              }
                // NOTE: Here, gn_ptr[nvsh]  is equivalent to
                //   Sp_clj_fqh[ iv_clj ][ kg*N_SPB_MX + nvsh]
                //
         
         x7 = x7 + x2 * x3 * x6 * Vppin->get(ia,im,ib,iu,J2/2,&indx);

        //TST//if ((2==J2)&&(0==iv_clj)) cout << nmsh << "   "  << nvsh << "   --   " << x7<<endl;
        }
        }

        U_ehf = U_ehf + x1 * x7;

      }   // J2 , iv_clj
     }

     // finally:
     Uab_EHF[nr*NDIM+nc] = U_ehf;
     Uab_EHF[nc*NDIM+nr] = U_ehf;

        //TST//cout << "--------" <<nr << "   "  << nc << " -----\n";      
   }  // nr, nc
  }

  if (Norb_mx > NDIM) return (Norb_mx-NDIM);
  return 0;
  }




double SpProp_t::Compute_FirstOrderExpValue(VppInt_t *Vppin, int ist, int iend){
                                      // NOTE: defaults ist=0 , iend=10000 are
                                      //   set in the class initialization.

  int    NDIM, ish, test, nr, nc, kg;
  double x1,   x2,  U1b;
  double *qh_ptr;

  ist = (ist > 0) ? ist : 0;
  if (iend >= MdSp->nsubsh) iend = (MdSp->nsubsh)-1;

  NDIM=0;
  for(ish=0; ish<MdSp->nsubsh; ++ish) if (NDIM<MdSp->MSp_no[ish]) NDIM=MdSp->MSp_no[ish];

  double *Uab = new double[NDIM*NDIM];
  for(ish= 0; ish<NDIM*NDIM; ++ish) Uab[ish]=0.0;

  U1b = 0.0;
  for(ish=ist; ish<=iend; ++ish) {
    test = ExtendedHartreeFock_Pot(Uab, NDIM, Vppin, ish);
    if (test) cerr << "\nWarning: the calculation of the (E)HF potential gave an error (" << test << ")\n\n";

    qh_ptr = Sp_clj_fqh[ish] + N_SPB_MX*Sp_clj_nh[ish];
    x1=0.0;
    for(kg=0; kg<Sp_clj_nh[ish]; ++kg) {

     qh_ptr -= N_SPB_MX;
     x2 = 0.0;
     for(nr=0; nr<MdSp->MSp_no[ish]; ++nr) {
      for(nc=0; nc<MdSp->MSp_no[ish]; ++nc) {
       x2 += Uab[nr*NDIM+nc] * qh_ptr[nc] * qh_ptr[nr];
       }
      }
     x1 += x2;

     }
    x1 *= double(MdSp->MSp_2j[ish]+1)/2.0; // 1/2 is the symm. factor, then one
                       // still needs to account for the multiplicity of the 
                       // otbit contracted with Uab_EHF.

    U1b += x1;

    //cout << MdSp->MSp_name[ish] << " :  " << x1 << endl;
    }

  delete [] Uab; //Uab=NULL;

  return (U1b);
  }

