
#include <iostream>
//#include <fstream>
//#include <iomanip>  //TST//
using namespace std;

#include <math.h>


#include "BcDor-Global_variables.hh"
#include "BcDor-Run_vars.hh"


void MGKoltun_sumrule(void ) {
     //
     //  Copute Galiski-Migdal-Koltun sum rule and other basic
     // quantities for the sop. propagator.
     //

     SpProp_t propV(SpProp_file,MS);
     //propV.write("sp_prop-wt");
     if (propV.n_qp + propV.n_qh < 300) propV.cout_qplist();  // if too many poles, better not to output them...
     //propV.cout_byshell();
     //propV.count_qp_confs(&n1,&n2,&n3,&n4,&n5);

     int ish_a,nsh_a,nsh_b;
     
     double   U1body_tot[2],E_tot[2],A_tot[2],Vtot_1st_ord,r_rms_1st_ord;

     double  x1,x2;
     
     int i, ka;


     for(i=0; i<2; ++i) {U1body_tot[i]=0.0; E_tot[i]=0.0;  A_tot[i]=0.0;}
    

     for(ka=0; ka<propV.n_qh; ++ka) {
      ish_a = propV.n_clj_h[ka];
      x1 = x2 = 0.0;
      for(nsh_a=0; nsh_a<MS->MSp_no[ish_a]; ++nsh_a) {
        x1 += propV.fqh[ka*propV.N_SPB_MX+nsh_a]*
              propV.fqh[ka*propV.N_SPB_MX+nsh_a];
        for(nsh_b=0; nsh_b<MS->MSp_no[ish_a]; ++nsh_b) {
        x2 += propV.fqh[ka*propV.N_SPB_MX+nsh_a]*
              propV.fqh[ka*propV.N_SPB_MX+nsh_b]*
              fab_tkin(MS->MSp_n[ish_a][nsh_a] , MS->MSp_n[ish_a][nsh_b] , MS->MSp_l[ish_a]) ;
          }
        }
      x1 *= double(MS->MSp_2j[ish_a] +1);
      x2 *= double(MS->MSp_2j[ish_a] +1)*MS->htom*U1body_fac;
      U1body_tot[MS->MSp_ch[ish_a]] += x2;
      A_tot[MS->MSp_ch[ish_a]]      += x1;
      E_tot[MS->MSp_ch[ish_a]]      += (x1*propV.eh[ka]+x2)/2.0;
      }

      cout << "\nneutrons:";
      cout << "\nU1body_tot_n = " << U1body_tot[0];
      cout << "\nV_tot_n      = " << E_tot[0]-U1body_tot[0];
      cout << "\nE_tot_n      = " << E_tot[0];
      cout << "\nA_tot_n      = " << A_tot[0]           << endl;
      cout << "\nprotons:";
      cout << "\nU1body_tot_p = " << U1body_tot[1];
      cout << "\nV_tot_p      = " << E_tot[1]-U1body_tot[1];
      cout << "\nE_tot_p      = " << E_tot[1];
      cout << "\nA_tot_p      = " << A_tot[1]           << endl;
      cout << "\nU1body_tot = " << U1body_tot[0]+U1body_tot[1];              
      cout << "\nV_tot      = " << E_tot[0]-U1body_tot[0]+E_tot[1]-U1body_tot[1];
      cout << "\nE_tot      = " << E_tot[0]+E_tot[1];
      cout << "\nE_tot/A    = " << (E_tot[0]+E_tot[1])/double(A);
      cout << "\nA_tot      = " << A_tot[0]+A_tot[1]           << endl;

      Vtot_1st_ord=propV.Compute_FirstOrderExpValue(Vpp);
      cout << "\n <V_tot> at 1st ord. = " << Vtot_1st_ord <<"\n";

  //
  // if (RrmsCalc) is set calculates the rms radius
  //
  switch (RrmsCalc) {
	case 1:
	  x1 = Calc_Rrms_OneBody_rirj(&propV, Rrms, NULL, true);
	  break;
	case 2:
	  x1 = Calc_Rrms_TwoBody(&propV, Rrms, true);
	  break;
	default:
	  break;
  }
  

  return;}


//
//  Calculates the second order self-energy (for all
// partial waves and contracts it with a propagator
// to generate the 2nd order MBPT diagrams.
//
//  This gives the MBPT2 energy ONLY if the input
// propagator is a HF one.
//
void Second_order_diag(void ) {

/*
  Set_nucleus_data();
  Load_ModSp();
  Load_interaction();
  //Load_SpProp();
*/

  double x1,x2,x3,x4,x5,x6,x7,x8,x9;
  double y7,y8,y9;

  double  *fqp, *fqh,  *msf;
  double   eqp,  eqh,   esf;
  int ish,nsh_a,nsh_b,norb;

  int ka,nb;

  int n1,n3,n4;

  int NDSH;
  int NDIM=0;
  for(ish=0; ish<MS->nsubsh; ++ish) if (NDIM<MS->MSp_no[ish]) NDIM=MS->MSp_no[ish];

  double *Uab = new double[NDIM*NDIM];
  for(int i= 0; i<NDIM*NDIM; ++i) Uab[i]=0.0;

  SpProp_t propV(SpProp_file,MS);
  if (propV.n_qp + propV.n_qh < 300) propV.cout_qplist();  // if too many poles, better not to output them...

  Wide_self_energy_t sfV(MS,Vpp,&propV);


  double s1[2],s2[2],s3[2],s4[2],s5[2],s6[2],s7[2],s8[2],s9[2];

  for(int i=0; i<2; ++i) {s1[i]=0.0;s2[i]=0.0;s3[i]=0.0;s4[i]=0.0;
                          s5[i]=0.0;s6[i]=0.0;s7[i]=0.0;s8[i]=0.0;s9[i]=0.0;}
    

    ish = -100;
    for(ka=0; ka<propV.n_qh; ++ka) {
      n1 = propV.n_clj_h[ka];
      
      if (n1 != ish) {
        ish = n1;
        sfV.set_clj(ish);
        for(int i= 0; i<NDIM*NDIM; ++i) Uab[i]=0.0;
        Get_Tab_kin_ho(Uab, NDIM, ish, MS);
        for(int i= 0; i<NDIM*NDIM; ++i) Uab[i]*=U1body_fac;
        sfV.Make_1bd_SelfEn(Uab,NDIM);
        sfV.Build_ExtendedHartreeFock();
        sfV.Build_Dynamic_SlfEn_2ndOrder(&n3, &n4);
        NDSH = sfV.NTOT_ORB;
      }
      eqh = propV.eh[ka];
      fqh = propV.fqh + (ka*propV.N_SPB_MX);
      norb = MS->MSp_no[ish];
      
      x1 = x2 = x3 = x6 = 0.0;
      for(nsh_a=0; nsh_a<norb; ++nsh_a) {
       for(nsh_b=0; nsh_b<norb; ++nsh_b) {
        x1 += fqh[nsh_a]*fqh[nsh_b]*sfV.Sig_1b[nsh_a*NDSH + nsh_b];
        x2 += fqh[nsh_a]*fqh[nsh_b]*sfV.Sig_HF[nsh_a*NDSH + nsh_b] ;
       }
       x3 += fqh[nsh_a]*fqh[nsh_a]*eqh;
       x6 += fqh[nsh_a]*fqh[nsh_a];
      }

      x4 = 0.0;
      msf= sfV.Sig_dyn_fw + 2;
      for(int i=0; i<sfV.N_PLS_fw; ++i) {
        esf = *(msf-2);

        y7 = 0.0;
        for(nsh_a=0; nsh_a<norb; ++nsh_a) y7 += fqh[nsh_a]*msf[nsh_a];

        x4 += y7*y7/(eqh-esf);
        msf += (sfV.NTOT_ORB + 2);
      }

      x5 = 0.0;
      msf = sfV.Sig_dyn_bk + 2;
      for(int i=0; i<sfV.N_PLS_bk; ++i) {
        esf = *(msf-2);

        y7 = 0.0;
        for(nsh_a=0; nsh_a<norb; ++nsh_a) y7 += fqh[nsh_a]*msf[nsh_a];

        x5 += y7*y7/(eqh-esf);
        msf += (sfV.NTOT_ORB + 2);
      }


    x8 = x9 = 0.0;
    for(nb=0; nb<propV.n_qp; ++nb) {
      if (ish != propV.n_clj_p[nb]) continue;
      eqp = propV.ep[nb];
      fqp = propV.fqp + (nb*propV.N_SPB_MX);
      norb = MS->MSp_no[ish];
      
      y8 = y9 = 0.0; 
      for(nsh_a=0; nsh_a<norb; ++nsh_a) {
       for(nsh_b=0; nsh_b<norb; ++nsh_b) {
        y8 += fqh[nsh_a]*fqp[nsh_b]*sfV.Sig_1b[nsh_a*NDSH + nsh_b];
        y9 += fqh[nsh_a]*fqp[nsh_b]*sfV.Sig_HF[nsh_a*NDSH + nsh_b] ;
       }
      }
      x8 += y8*y8 / (eqp - eqh);
      x9 += y9*y9 / (eqp - eqh);
    
    }


    x1 *= double(MS->MSp_2j[ish]+1);
    x2 *= double(MS->MSp_2j[ish]+1)/2.0;
    x3 *= double(MS->MSp_2j[ish]+1);
    x4 *= double(MS->MSp_2j[ish]+1)/2.0;
    x5 *= double(MS->MSp_2j[ish]+1)/2.0;
    x6 *= double(MS->MSp_2j[ish]+1);
    x8 *= double(MS->MSp_2j[ish]+1);
    x9 *= double(MS->MSp_2j[ish]+1);

    cout << "\n Shell: "<<MS->MSp_name[ish]
         << "    " << x1 << "    " << x2 << "    " << x3
         << "    " << x4 << "    " << x5 << "    " << x6
         << "    " << x8 << "    " << x9;

    n1 = MS->MSp_ch[ish];


    s1[n1] += x1;
    s2[n1] += x2;
    s3[n1] += x3;
    s4[n1] += x4;
    s5[n1] += x5;
    s6[n1] += x6;
    s8[n1] += x8;
    s9[n1] += x9;

    }


    ish = -100;
    for(ka=0; ka<propV.n_qp; ++ka) {
      n1 = propV.n_clj_p[ka];
      
      if (n1 != ish) {
        ish = n1;
        
        sfV.set_clj(ish);
        //sfV.Build_ExtendedHartreeFock();
        //Get_Tab_kin_ho(Uab, NDIM, ish, MS);
        //for(int i= 0; i<NDIM*NDIM; ++i) Uab[i]*=U1body_fac;
        //sfV.Make_1bd_SelfEn(Uab,NDIM);
        sfV.Build_Dynamic_SlfEn_2ndOrder(&n3, &n4);
      }
      eqp = propV.ep[ka];
      fqp = propV.fqp + (ka*propV.N_SPB_MX);
      norb = MS->MSp_no[ish];
      
      x7 = 0.0;
      msf = sfV.Sig_dyn_bk + 2;
      for(int i=0; i<sfV.N_PLS_bk; ++i) {
        esf = *(msf-2);

        x1 = 0.0;
        for(nsh_a=0; nsh_a<norb; ++nsh_a) x1 += fqp[nsh_a]*msf[nsh_a];

        x7 -= x1*x1/(eqp-esf);
        msf += (sfV.NTOT_ORB + 2);
      }

      x7 *= double(MS->MSp_2j[ish]+1)/2.0;

      cout << "\n Shell: "<<MS->MSp_name[ish] << "    " << x7 << flush;

      n1 = MS->MSp_ch[ish];

      s7[n1] += x7;
    }



    cout << "\n\nneutrons:";
    cout << "\nT_n          = " << s1[0];
    cout << "\nV1_n         = " << s2[0];
    cout << "\nV20_n        = " << s7[0];
    cout << "\nV2f_n        = " << s4[0];
    cout << "\nV2b_n        = " << s5[0];
    cout << "\nV2cV_n       = " << s9[0];
    cout << "\nV2cT_n       = " << s8[0];
    cout << "\nVt           = " << s2[0]+s4[0];
    cout << "\nVt_Klt_n     = " << (s3[0]-s1[0])/2.0;
    cout << "\nE_Kolt_n     = " << (s3[0]+s1[0])/2.0;
    cout << "\nE1_n         = " << s1[0]+s2[0];
    cout << "\nE2_n         = " << s4[0];
    cout << "\nE_tot_n      = " << s1[0]+s2[0]+s4[0];
    cout << "\nA_tot_n      = " << s6[0];

    cout << "\n\nprotons:";
    cout << "\nT_p          = " << s1[1];
    cout << "\nV1_p         = " << s2[1];
    cout << "\nV20_p        = " << s7[1];
    cout << "\nV2f_p        = " << s4[1];
    cout << "\nV2b_p        = " << s5[1];
    cout << "\nV2cV_p       = " << s9[1];
    cout << "\nV2cT_p       = " << s8[1];
    cout << "\nVt           = " << s2[1]+s4[1];
    cout << "\nVt_Klt_p     = " << (s3[1]-s1[1])/2.0;
    cout << "\nE_Kolt_p     = " << (s3[1]+s1[1])/2.0;
    cout << "\nE1_p         = " << s1[1]+s2[1];
    cout << "\nE2_p         = " << s4[1];
    cout << "\nE_tot_p      = " << s1[1]+s2[1]+s4[1];
    cout << "\nA_tot_p      = " << s6[1];

    cout << "\n\ntotal:";
    cout << "\nT            = " << s1[0]+s1[1];
    cout << "\nV1           = " << s2[0]+s2[1];
    cout << "\nV20          = " << s7[0]+s7[1];
    cout << "\nV2f          = " << s4[0]+s4[1];
    cout << "\nV2b          = " << s5[0]+s5[1];
    cout << "\nV2cV         = " << s9[0]+s9[1];
    cout << "\nV2cT         = " << s8[0]+s8[1];
    cout << "\nVt           = " << s2[0]+s2[1]+s4[0]+s4[1];
    cout << "\nVt_Klt       = " << (s3[0]+s3[1]-s1[0]-s1[1])/2.0;
    cout << "\nE_Kolt       = " << (s3[0]+s3[1]+s1[0]+s1[1])/2.0;
    cout << "\nE1           = " << s1[0]+s1[1]+s2[0]+s2[1];
    cout << "\nE2           = " << s4[0]+s4[1];
    cout << "\nE_tot        = " << s1[0]+s1[1]+s2[0]+s2[1]+s4[0]+s4[1];
    cout << "\nA_tot        = " << s6[0]+s6[1];

  delete [] Uab; Uab = NULL;

  return;}
  
