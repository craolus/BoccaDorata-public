//
//
//   Build the pp/hh propagator based on the successive approximations
//  to the Bethe-Salpeter-like equation (BSE):
//     (D)RPA,  2Pi-RPA, ...
//
// gII_b.C: 2qp/2qh bases  &  solution of (D)RPA for X and Y vectors
//
// C.Barbieri, GSI, December 2006.  (C.Barbieri@gsi.de)
//

#include <iostream>
#include <iomanip>
using namespace std;

#include <cstdlib>
//#include <cstdio>
#include <cmath>

//#define DPI   3.14159265358979323846
#define SQR2  1.414213562373095048801688


#include "BcDor-Ang_momenta.hh"
//#include "BcDor-HO_radial_me.hh"


#include "BcDor-Global_variables.hh"
#include "BcDor-Sp_classes.hh"
#include "BcDor-Int_classes.hh"
#include "BcDor-gII_classes.hh"
//#include "Utilities.hh"


// LAPACK driver for diagonalizaton
extern "C" void dgeev_(char*,   char*, int*, double*, int*, double*, double*,
                       double*, int*,  double*, int*, double*, int*, int* );

//
// global variables:
//
static int ipp, ihh;//, iab;
static int nr, nc;
static double xy_r,xy_c;
static double x1;//, x2, x3;
//
static double *fqp_a , *fqp_b, *fqp_g, *fqp_d;
static int ish_a, ia1, ia2, na, ka, ia;
static int ish_b, ib1, ib2, nb, kb, ib;
static int ish_g, ig1, ig2, ng, ig;// , kg
static int ish_d, id1, id2, nd, id;// , kd
//
static int i1;//,i2;
static double *dptr1, *dptr2;
const char chip[2] = { '+' , '-'};




static bool iALL;

int gII_prop_t::Count_LaddRPA_basis(int Jf_in, int dT_in, int pif_in,
                         int *npp, int *nhh, int icut/*=0*/, int ntcut/*=0*/) {
          //
          // icut == 0,  include all posible poles
          //       > 0,  choose poles according to (n1+n2)*|icut|  >= ntcut 
          //       < 0,  choose poles according to (n1+n2)*|icut|  <  ntcut 

  iALL = true;
  if ( 0 != icut ) iALL=false;
  if ( 0 >  icut ) ntcut = 1- ntcut;

  *npp = 0;
  *nhh = 0;

  for(na=0; na<gsp->n_qp; ++na) {
   for(nb=na; nb<gsp->n_qp; ++nb) {
    ish_a = gsp->n_clj_p[na];
    ish_b = gsp->n_clj_p[nb];

    if ((!iALL) && ( (gsp->n1p[na]+gsp->n1p[nb])*icut < ntcut) ) continue;

    if ((pif_in + MdSp->MSp_ip[ish_a] + MdSp->MSp_ip[ish_b])%2) continue;
    if (am.TriIneq(MdSp->MSp_2j[ish_a],MdSp->MSp_2j[ish_b],2*Jf_in)) continue;
    if ((na==nb) && (abs( (MdSp->MSp_2j[ish_a]+MdSp->MSp_2j[ish_b])/2 - Jf_in +1)%2) ) continue;

    if ( MdSp->MSp_ch[ish_a] + MdSp->MSp_ch[ish_b] == dT_in ) ++(*npp);
   } }  // end p-p loop

  for(ka=0; ka<gsp->n_qh; ++ka) {
   for(kb=ka; kb<gsp->n_qh; ++kb) {
    ish_a = gsp->n_clj_h[ka];
    ish_b = gsp->n_clj_h[kb];

    if ((!iALL) && ( (gsp->n1h[ka]+gsp->n1h[kb])*icut < ntcut) ) continue;

    if ((pif_in + MdSp->MSp_ip[ish_a] + MdSp->MSp_ip[ish_b])%2) continue;
    if (am.TriIneq(MdSp->MSp_2j[ish_a],MdSp->MSp_2j[ish_b],2*Jf_in)) continue;
    if ((ka==kb) && (abs( (MdSp->MSp_2j[ish_a]+MdSp->MSp_2j[ish_b])/2 - Jf_in +1)%2) ) continue;

    if ( MdSp->MSp_ch[ish_a] + MdSp->MSp_ch[ish_b] == dT_in ) ++(*nhh);
   } }  // end h-h loop


//  cout << " Ladd-(D)RPA dimensions for Jf/pif = " << Jf_in << chip[pif_in] << " :" 
//       << "       dT=" <<  dT_in << " --> npp=" << *npp
//       << "  ,    dT=" << -dT_in << " --> nhh=" << *nhh;
  if (!iALL) {
    cout << "\n ( the selection rule  " << abs(icut) <<"*[n1(a) + n1(b)] ";
    if (icut > 0) cout << ">= " << ntcut;
             else cout << "< "  << 1- ntcut;
    cout << "  was used).\n";
  }
//  cout << endl;

  return 0;}

int gII_prop_t::Build_LaddRPA_basis(int Jf_in, int dT_in, int pif_in,
                            int icut/*=0*/, int ntcut/*=0*/, int iplot/*=0*/) {
          //
          // icut == 0,  include all posible poles
          //       > 0,  choose poles according to (n1+n2)*|icut|  >= ntcut 
          //       < 0,  choose poles according to (n1+n2)*|icut|  <  ntcut 
          //
          //  iplot == 0,  do not print anything except the dimensions
          //  iplot == 1,  print minor info but do not print the basis
          //  iplot >= 2,  send the basis' details to 'cout'
          //

  Jf  = Jf_in;
  dT  = dT_in;
  pif = pif_in;
  free_2qp_basis();

  if (iplot) cout << "\n\n Generating 2qp and 2qh bases... \n";
  Count_LaddRPA_basis(Jf, dT, pif, &ipp, &ihh, icut, ntcut);
  allocate_2qp_basis(ipp, ihh);
  if (iplot) cout << " ...memory allocated...\n";

  // this goes here because the count routine also uses iALL...
  iALL = true;
  if ( 0 != icut ) iALL=false;
  if ( 0 >  icut ) ntcut = 1- ntcut;

  // NOTE: the routines 'get_ld_f_indx' and 'get_ld_b_indx' rely on the
  //   following ordering of the of na-nb and ka-kb quantum numbers (obtained
  //   by keeping nb/kb in the internal loops). Thus, the loops below should not
  //   be changed to avoid errors...
  ipp = 0;
  for(na=0; na<gsp->n_qp; ++na) {
   for(nb=na; nb<gsp->n_qp; ++nb) {
    ish_a = gsp->n_clj_p[na];
    ish_b = gsp->n_clj_p[nb];

    if ((!iALL) && ( (gsp->n1p[na]+gsp->n1p[nb])*icut < ntcut) ) continue;

    if ((pif + MdSp->MSp_ip[ish_a] + MdSp->MSp_ip[ish_b])%2) continue;
    if (am.TriIneq(MdSp->MSp_2j[ish_a],MdSp->MSp_2j[ish_b],2*Jf)) continue;
    if ((na==nb) && (abs( (MdSp->MSp_2j[ish_a]+MdSp->MSp_2j[ish_b])/2 - Jf +1)%2) ) continue;

    if ( MdSp->MSp_ch[ish_a] + MdSp->MSp_ch[ish_b] == dT ) {
       e_pp[ipp] = gsp->ep[na]+gsp->ep[nb];
      n1_pp[ipp] = na;
      n2_pp[ipp] = nb;
      ++ipp;
      }
  } }  // end p-p loop

  ihh=0;
  for(ka=0; ka<gsp->n_qh; ++ka) {
   for(kb=ka; kb<gsp->n_qh; ++kb) {
    ish_a = gsp->n_clj_h[ka];
    ish_b = gsp->n_clj_h[kb];

    if ((!iALL) && ( (gsp->n1h[ka]+gsp->n1h[kb])*icut < ntcut) ) continue;

    if ((pif + MdSp->MSp_ip[ish_a] + MdSp->MSp_ip[ish_b])%2) continue;
    if (am.TriIneq(MdSp->MSp_2j[ish_a],MdSp->MSp_2j[ish_b],2*Jf)) continue;
    if ((ka==kb) && (abs( (MdSp->MSp_2j[ish_a]+MdSp->MSp_2j[ish_b])/2 - Jf +1)%2) ) continue;

    if ( MdSp->MSp_ch[ish_a] + MdSp->MSp_ch[ish_b] == dT ) {
       e_hh[ihh] = gsp->eh[ka]+gsp->eh[kb];
      k1_hh[ihh] = ka;
      k2_hh[ihh] = kb;
      ++ihh;
      }
  } }  // end h-h loop

  if ((ipp > 0) && (ipp > n_pp_alloc)) cerr << "\n\n p-p arrays are out of boundaries!!!!!!!!!!!!\n\n";
  if ((ihh > 0) && (ihh > n_hh_alloc)) cerr << "\n\n h-h arrays are out of boundaries!!!!!!!!!!!!\n\n";

  n_pp_basis   = ipp;
  n_hh_basis   = ihh;

  if (iplot > 1) {
    cout << " ...2qp basis:\n";
    for(ipp =0; ipp<n_pp_basis; ++ ipp) 
      cout << "na,nb =  " << n1_pp[ipp] << MdSp->MSp_name[gsp->n_clj_p[n1_pp[ipp]]]
                 << " , " << n2_pp[ipp] << MdSp->MSp_name[gsp->n_clj_p[n2_pp[ipp]]]
                 << "   e_pp=" << e_pp[ipp] << endl;
    cout << endl;
    cout << " ...2qh basis:\n";
    for(ihh =0; ihh<n_hh_basis; ++ ihh) 
      cout << "ka,kb =  " << k1_hh[ihh] << MdSp->MSp_name[gsp->n_clj_h[k1_hh[ihh]]]
           << " , " << k2_hh[ihh] << MdSp->MSp_name[gsp->n_clj_h[k2_hh[ihh]]]
           << "   e_hh=" << e_hh[ihh] << endl;
    cout << endl;
  }

  if (iplot) cout << " ...2qp and 2qh bases are done!";
  if (!iALL) {
    cout << "\n ( the selection rule  " << abs(icut) <<"*[n1(a) + n1(b)] ";
    if (icut > 0) cout << ">= " << ntcut;
             else cout << "< "  << 1- ntcut;
  cout << "  was used).\n";
  }
  if (iplot) cout<<"\n\n";

  n_f_basis = (0 > n_pp_basis) ? 0 : n_pp_basis;
  n_b_basis = (0 > n_hh_basis) ? 0 : n_hh_basis;

  return n_pp_basis+n_hh_basis;} // returns tot. dim. of the Ladd-DRPA matrix

  //
  // NOTE: the routines 'get_ld_f_indx' and 'get_ld_b_indx' rely on the
  //   ordering of the of na-nb and ka-kb quantum numbers done by keeping
  //   nb(kb) in the internal loop. These will require chages if the routine
  //   'Build_LaddRPA_basis' is modified.
  //


int gII_prop_t::get_ld_f_indx(int *imid, int na, int nb, int *iph) {
  int iseek,iup,idown,iab_mid,nmult;
  nmult = gsp->n_qp;
  if (na > nb) {
    iseek = (nb*nmult) +na;
    *iph = abs((MdSp->MSp_2j[gsp->n_clj_p[na]] +
                MdSp->MSp_2j[gsp->n_clj_p[nb]])/2 - Jf + 1);
  } else {
    iseek = (na*nmult) +nb;
    *iph = 0;
  }
  idown = 0;
  iup   = n_pp_basis-1;
  while(idown <= iup) {
    *imid  = (iup+idown)/2;
    iab_mid = (n1_pp[*imid]*nmult) +n2_pp[*imid];
    if (iseek == iab_mid) return 0;
    if (iseek < iab_mid) {iup   = (*imid)-1;}
    else {idown = (*imid)+1;}
  }
  
  *imid = -1000;
  return 1;}  // the 2qp configuration was not found

int gII_prop_t::get_ld_b_indx(int *imid, int ka, int kb, int *iph) {
  int iseek,iup,idown,iab_mid,nmult;
  nmult = gsp->n_qh;
  if (ka > kb) {
    iseek = (kb*nmult) +ka;
    *iph = abs((MdSp->MSp_2j[gsp->n_clj_h[ka]] +
                MdSp->MSp_2j[gsp->n_clj_h[kb]])/2 - Jf + 1);
  } else {
    iseek = (ka*nmult) +kb;
    *iph = 0;
  }
  idown = 0;
  iup   = n_hh_basis-1;
  while(idown <= iup) {
    *imid  = (iup+idown)/2;
    iab_mid = (k1_hh[*imid]*nmult) +k2_hh[*imid];
    if (iseek == iab_mid) return 0;
    if (iseek < iab_mid) {iup   = (*imid)-1;}
    else {idown = (*imid)+1;}
  }
  
  *imid = -1000;
  return 1;}  // the 2qh configuration was not found





static int ndimf,ndimb,ndim_tot;
static int NSPBMXloc;
static double *MTX, *MTX2;
//static int itest;
//static int  ndim_alloc,ndim_vects;

static bool  Make_ppladd_matrix, Make_hhladd_matrix, Make_DRPA_B_matrices, Free_ph_propagator;


int gII_prop_t::Fill_LaddRPA_mtx(int i_RPAcalc, double *Sab_LaddRPA, int ndim_vec_alloc, double *Free_poles/*=NULL*/, bool Add_free_poles/*=true*/) {
//
//  i_RPAcalc == 0   full (D)TDA: compute both forward and backward (i.e. pp
//                 and hh) amplitudes
//
//  i_RPAcalc == 1  simple pp(D)TDA (diagonalize only the forwardgoing part)
//
//  i_RPAcalc == 2  simple hh(D)TDA (diagonalize only the backwardgoing part)
//
//  i_RPAcalc == 3  (D)RPA (both forward and backward amplitudes must be
//                     computed in this case)
//
//  i_RPAcalc == 4  get the free propagator
//


  if ((i_RPAcalc < 0) || (i_RPAcalc > 4)) {cout << "gII_prop_t::Fill_LaddRPA_mtx: i_RPAcalc==" <<i_RPAcalc<< "is not a valid option \n"; return -100;}


  //itest = 1;

  //Clear_XYmtx(); 

  /////cout << " i_RPAcalc(0,1,2)?"; cin >>  i_RPAcalc;

  //cout << " i_RPAcalc=" << i_RPAcalc << endl;

  //
  // set up dimensions:
  //
  // NOTE: "n_f(b)_basis"  are meant to store the dimensions of the X and Y
  //    amplitudes, while "n_pp(hh)_basis" are the bases' dimensions. In 
  //    practice whenever these are both larger zero, they *must* always equal.
  //    Tere are cases in which "n_f(b)_basis" could be zero even if the
  //    2qp/2qh basis is not vanishing, for example when only a part of the 
  //    full RPA matrix is diagonalized (only fw, or bk, TDA).
  //
  //     NOTE that the "n_f(b)_basis" *must* always be >= zero.
  //     However, if the basis is not allocated, the "n_pp(hh)_basis" can
  //    be negative.
  //
  n_f_basis = (0 > n_pp_basis) ? 0 : n_pp_basis;
  n_b_basis = (0 > n_hh_basis) ? 0 : n_hh_basis;
  ndimf = n_f_basis;
  ndimb = n_b_basis;
     //
     // NOTE: 'ndimf', 'ndimb' and 'ndim' CANNOT be < 0 or the algorithm
     //   below will chrash!

  if (ndim_vec_alloc < 1) {cout << "gII_prop_t::Fill_LaddRPA_mtx: either there are no solutions for this parial wave or you forgot to initialize the 2qp/2qh bases! \n"; return -101;}

  
  Make_ppladd_matrix   = true;
  Make_hhladd_matrix   = true;
  Make_DRPA_B_matrices = true;
  if (0 == i_RPAcalc) {Make_DRPA_B_matrices = false;}
  if (1 == i_RPAcalc) {ndimb=0; Make_hhladd_matrix = false; Make_DRPA_B_matrices = false;}
  if (2 == i_RPAcalc) {ndimf=0; Make_ppladd_matrix = false; Make_DRPA_B_matrices = false;}
  if (4 == i_RPAcalc) {Free_ph_propagator = true;  Add_free_poles = true;}
  ndim_tot  = ndimf+ndimb;
  //
  if (ndim_tot > ndim_vec_alloc) {
	cerr << " ERROR (gII::Fill_LaddRPA_mtx): value of ndim_vec(="<<ndim_vec_alloc<<") too small. -- STOP\n";
	exit(-100);
  }


  //cout << "Dimensions of the pp/hh(D)RPA mtx:  fw=" << ndimf << " ,   bk="   << ndimb << " ,  tot=" << ndim_tot << endl;  


  // clean up the mtx:
  dptr2 = Sab_LaddRPA + (ndim_tot-1)*ndim_vec_alloc + ndim_tot;
  for(dptr1=Sab_LaddRPA; dptr1<dptr2; ++dptr1) (*dptr1)=0.0; //clean up array

  NSPBMXloc = gsp->N_SPB_MX;

  //
  // Buid mtx:
  // =========

//pp_loop:
  if (Free_ph_propagator) goto Add_energies;
  if ((!Make_ppladd_matrix) && (!Make_DRPA_B_matrices)) goto hh_loop;
  //if (ndimf < 1) goto hh_loop;


  for(nc=0; nc<n_pp_basis; ++nc) {
   ng = n1_pp[nc];
   ish_g = gsp->n_clj_p[ng];
   ig1 = (MdSp->MSp_n[ish_g]   - MdSp->Mor_n);
   ig2 = ig1 + MdSp->MSp_no[ish_g];
   //
   nd = n2_pp[nc];
   ish_d = gsp->n_clj_p[nd];
   id1 = (MdSp->MSp_n[ish_d]   - MdSp->Mor_n);
   id2 = id1 + MdSp->MSp_no[ish_d];

   //
   // Af matrix:
   MTX  = Sab_LaddRPA +         nc*ndim_vec_alloc;
   //
   if (4 != i_RPAcalc) // only diag. energies are needed for the free propagator
   for(nr=0; nr<n_pp_basis ; ++nr) {
    na = n1_pp[nr];                              
    ish_a = gsp->n_clj_p[na];                   
    ia1 = (MdSp->MSp_n[ish_a]   - MdSp->Mor_n);
    ia2 = ia1 + MdSp->MSp_no[ish_a];
    //
    nb = n2_pp[nr];
    ish_b = gsp->n_clj_p[nb];
    ib1 = (MdSp->MSp_n[ish_b]   - MdSp->Mor_n);
    ib2 = ib1 + MdSp->MSp_no[ish_b];

    x1 = 0.0;
    fqp_a = gsp->fqp + (na* NSPBMXloc);
    for(ia=ia1; ia<ia2; ++ia) {
     fqp_b = gsp->fqp + (nb* NSPBMXloc);
     for(ib=ib1; ib<ib2; ++ib) {
      xy_r = (*fqp_a) * (* fqp_b);
      if (ia==ib) xy_r *= SQR2;
    
      fqp_g = gsp->fqp + (ng* NSPBMXloc);
      for(ig=ig1; ig<ig2; ++ig) {
       fqp_d = gsp->fqp + (nd* NSPBMXloc);
       for(id=id1; id<id2; ++id) {
         xy_c =  (*fqp_g) * (* fqp_d);
         if (ig==id) xy_c *= SQR2;
         x1 += xy_r * vpp->get(ia,ib,ig,id,Jf,&i1) * xy_c;
         //x1 += xy_r * vpp->get(ia,ib,ig,id,Jf,&i1) * (*fqp_g) * (* fqp_d);
  
        ++fqp_d;}
       ++fqp_g;}

      ++fqp_b;}
     ++fqp_a;}
     if (na==nb) x1 /= SQR2;
     if (ng==nd) x1 /= SQR2;

     (*MTX)  =  x1;

    ++MTX;}  // nr loop

   if (NULL != Free_poles) Free_poles[nc] = e_pp[nc];

   if (!Make_DRPA_B_matrices) continue;

   //if (ndimb < 1) continue;

   //
   // B-rpa matrix:
   MTX   = Sab_LaddRPA +  ndimf+         nc*ndim_vec_alloc;
   MTX2  = Sab_LaddRPA +         nc + ndimf*ndim_vec_alloc;
   for(nr=0; nr<n_hh_basis ; ++nr) {
    na = k1_hh[nr];                              
    ish_a = gsp->n_clj_h[na];                   
    ia1 = (MdSp->MSp_n[ish_a]   - MdSp->Mor_n);
    ia2 = ia1 + MdSp->MSp_no[ish_a];
    //
    nb = k2_hh[nr];
    ish_b = gsp->n_clj_h[nb];
    ib1 = (MdSp->MSp_n[ish_b]   - MdSp->Mor_n);
    ib2 = ib1 + MdSp->MSp_no[ish_b];

    x1 = 0.0;
    fqp_a = gsp->fqh + (na* NSPBMXloc);
    for(ia=ia1; ia<ia2; ++ia) {
     fqp_b = gsp->fqh + (nb* NSPBMXloc);
     for(ib=ib1; ib<ib2; ++ib) {
      xy_r = (*fqp_a) * (* fqp_b);
      if (ia==ib) xy_r *= SQR2;
    
      fqp_g = gsp->fqp + (ng* NSPBMXloc);
      for(ig=ig1; ig<ig2; ++ig) {
       fqp_d = gsp->fqp + (nd* NSPBMXloc);
       for(id=id1; id<id2; ++id) {
         xy_c =  (*fqp_g) * (* fqp_d);
         if (ig==id) xy_c *= SQR2;
         x1 += xy_r * vpp->get(ia,ib,ig,id,Jf,&i1) * xy_c;
  
        ++fqp_d;}
       ++fqp_g;}

      ++fqp_b;}
     ++fqp_a;}
     if (na==nb) x1 /= SQR2;
     if (ng==nd) x1 /= SQR2;

     (*MTX)   =  -x1;
     (*MTX2)  =   x1;

    ++MTX;
    MTX2 += ndim_vec_alloc;
    }  // nr loop

   } // nc loop

hh_loop:
  if (!Make_hhladd_matrix) goto Add_energies;
  //if (ndimb < 1) goto hh_loop;


  for(nc=0; nc<n_hh_basis; ++nc) {
   ng = k1_hh[nc];                             
   ish_g = gsp->n_clj_h[ng];                  
   ig1 = (MdSp->MSp_n[ish_g]   - MdSp->Mor_n);
   ig2 = ig1 + MdSp->MSp_no[ish_g];
   //
   nd = k2_hh[nc];
   ish_d = gsp->n_clj_h[nd];
   id1 = (MdSp->MSp_n[ish_d]   - MdSp->Mor_n);
   id2 = id1 + MdSp->MSp_no[ish_d];

   //
   // D-rpa matrix:
   MTX   = Sab_LaddRPA +  ndimf+ (ndimf+nc)*ndim_vec_alloc;
   //
   //if (4 != i_RPAcalc) // only diag. energies are needed for the free propagator
   for(nr=0; nr<n_hh_basis ; ++nr) {
    na = k1_hh[nr];                              
    ish_a = gsp->n_clj_h[na];                   
    ia1 = (MdSp->MSp_n[ish_a]   - MdSp->Mor_n);
    ia2 = ia1 + MdSp->MSp_no[ish_a];
    //
    nb = k2_hh[nr];
    ish_b = gsp->n_clj_h[nb];
    ib1 = (MdSp->MSp_n[ish_b]   - MdSp->Mor_n);
    ib2 = ib1 + MdSp->MSp_no[ish_b];

    x1 = 0.0;
    fqp_a = gsp->fqh + (na* NSPBMXloc);
    for(ia=ia1; ia<ia2; ++ia) {
     fqp_b = gsp->fqh + (nb* NSPBMXloc);
     for(ib=ib1; ib<ib2; ++ib) {
      xy_r = (*fqp_a) * (* fqp_b);
      if (ia==ib) xy_r *= SQR2;
    
      fqp_g = gsp->fqh + (ng* NSPBMXloc);
      for(ig=ig1; ig<ig2; ++ig) {
       fqp_d = gsp->fqh + (nd* NSPBMXloc);
       for(id=id1; id<id2; ++id) {
         xy_c =  (*fqp_g) * (* fqp_d);
         if (ig==id) xy_c *= SQR2;
         x1 += xy_r * vpp->get(ia,ib,ig,id,Jf,&i1) * xy_c;
  
        ++fqp_d;}
       ++fqp_g;}

      ++fqp_b;}
     ++fqp_a;}
     if (na==nb) x1 /= SQR2;
     if (ng==nd) x1 /= SQR2;

     (*MTX)   =  -x1;

    ++MTX;
    }  // nr loop

   if (NULL != Free_poles) Free_poles[nc+ndimf] = e_hh[nc];

   } // nc loop


Add_energies:

  if (Add_free_poles) {
	if (Make_ppladd_matrix) for(nc=0; nc<n_pp_basis; ++nc)  *(Sab_LaddRPA + nc *(1+ndim_vec_alloc)) += e_pp[nc];
	if (Make_hhladd_matrix) for(nc=0; nc<n_hh_basis; ++nc)  *(Sab_LaddRPA + (nc+ndimf) *(1+ndim_vec_alloc)) += e_hh[nc];
  }

  
  
  return 0; }

