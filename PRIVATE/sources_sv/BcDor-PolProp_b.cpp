//
//
//   Build the Polarozation propagator based on the successive approximation
//  to the Bethe-Salpeter eqaution (BSE):
//     (D)RPA,  2Pi-RPA, ...
//
// PolProp_b.C: 2qp/2qh bases  &  solution of (D)RPA for X and Y vectors
//
// C.Barbieri, GSI, December 2006.  (C.Barbieri@gsi.de)
//

#include <iostream>
#include <fstream>
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
#include "BcDor-PolProp_classes.hh"
//#include "Utilities.hh"


double Pandya(int, int, int, int, int, VppInt_t*);


//
// global variables:
//
static int iph, ihp;//, iab;
static int ndim_tot;
static int nr, nc;
static double xy_r;//,xy_c;
static double x1, x2, x3;
//
//
static double *fqp_a , *fqh_b, *fqp_g, *fqh_d;
static int ish_a, ia, ia1, ia2, na;
static int ish_b, ib, ib1, ib2, kb;
static int ish_g, ig, ig1, ig2, ng;
static int ish_d, id, id1, id2, kd;
//
static int j2a;//la, l2a,  ich_a, nsh_a, 
static int j2b;//lb, l2b,  ich_b, nsh_b, 
//
static double *dptr1, *dptr2;
//static char chip[2] = { '+' , '-'};

static bool iALL;


int Pol_prop_t::Count_phRPA_basis(int Jf_in, int dT_in, int pif_in,
								  int *nph, int *nhp, int icut/*=0*/, int ntcut/*=0*/) {
  //
  // icut == 0,  include all posible poles
  //       > 0,  choose poles according to (n1+n3)*|icut|  >= ntcut 
  //       < 0,  choose poles according to (n1+n3)*|icut|  <  ntcut 
  
  iALL = true;
  if ( 0 != icut ) iALL=false;
  if ( 0 >  icut ) ntcut = 1- ntcut;
  
  *nph = 0;
  *nhp = 0;
  
  for(na=0; na<gsp->n_qp; ++na) {
	for(kb=0; kb<gsp->n_qh; ++kb) {
	  ish_a = gsp->n_clj_p[na];
	  ish_b = gsp->n_clj_h[kb];
	  
	  if ((!iALL) && ( (gsp->n1p[na]+gsp->n1h[kb])*icut < ntcut) ) continue;
	  
	  if ((pif_in + MdSp->MSp_ip[ish_a] + MdSp->MSp_ip[ish_b])%2) continue;
	  if (am.TriIneq(MdSp->MSp_2j[ish_a],MdSp->MSp_2j[ish_b],2*Jf_in)) continue;
	  
	  
	  if ( MdSp->MSp_ch[ish_a] - MdSp->MSp_ch[ish_b] == dT_in ) ++(*nph);
	  if ( MdSp->MSp_ch[ish_b] - MdSp->MSp_ch[ish_a] == dT_in ) ++(*nhp);
	  
	} }  // end p-h loop
  
  //  cout << " phRPA dimensions for Jf/pif = " << Jf_in << chip[pif_in] << " :" 
  //       << "       dT=" <<  dT_in << " --> nph=" << *nph
  //       << "  ,    dT=" << -dT_in << " --> nhp=" << *nhp;
  if (!iALL) {
    cout << "\n ( the selection rule  " << abs(icut) <<"*[n1p(a) + n1h(b)] ";
    if (icut > 0) cout << ">= " << ntcut;
	else cout << "< "  << 1- ntcut;
    cout << "  was used).\n";
  }
  //  cout << endl;
  
  return 0;}

int Pol_prop_t::Build_phRPA_basis(int Jf_in, int dT_in, int pif_in,
								  int icut/*=0*/, int ntcut/*=0*/, int iplot/*=0*/) {
  //
  // icut == 0,  include all posible poles
  //       > 0,  choose poles according to (n1+n3)*|icut|  >= ntcut 
  //       < 0,  choose poles according to (n1+n3)*|icut|  <  ntcut 
  //
  //  iplot == 0,  do not print anything except the dimensions
  //  iplot == 1,  print minor info but do not print the basis
  //  iplot >= 2,  send the basis' details to 'cout'
  //
  
  Jf  = Jf_in;
  dT  = dT_in;
  pif = pif_in;
  free_2qp_basis();
  
  if (iplot) cout << "\n\n Generating qp-qh and qh-qp bases... \n";
  Count_phRPA_basis(Jf, dT, pif, &iph, &ihp);
  allocate_2qp_basis(iph, ihp);
  if (iplot) cout << " ...memory allocated...\n";
  
  // this goes here because the count routine also uses iALL...
  iALL = true;
  if ( 0 != icut ) iALL=false;
  if ( 0 >  icut ) ntcut = 1- ntcut;
  
  // NOTE: the routines 'get_Pi_f_indx' and 'get_Pi_b_indx' rely on the
  //   following ordering of the of na and kb quantum numbers (obtained by
  //   keeping kb in the internal loop). Thus, the loops below should not
  //   be changed to avoid errors...
  iph = 0;
  ihp = 0;
  for(na=0; na<gsp->n_qp; ++na) {
	for(kb=0; kb<gsp->n_qh; ++kb) {
	  ish_a = gsp->n_clj_p[na];
	  ish_b = gsp->n_clj_h[kb];
	  
	  if ((!iALL) && ( (gsp->n1p[na]+gsp->n1h[kb])*icut < ntcut) ) continue;
	  
	  if ((pif + MdSp->MSp_ip[ish_a] + MdSp->MSp_ip[ish_b])%2) continue;
	  if (am.TriIneq(MdSp->MSp_2j[ish_a],MdSp->MSp_2j[ish_b],2*Jf)) continue;
	  
	  
	  if ( MdSp->MSp_ch[ish_a] - MdSp->MSp_ch[ish_b]  == dT) {
		e_ph[iph] = gsp->ep[na]-gsp->eh[kb];
		n_ph[iph] = na;
		k_ph[iph] = kb;
		++iph;
      }
	  if ( MdSp->MSp_ch[ish_b] - MdSp->MSp_ch[ish_a]  == dT) {
		e_hp[ihp] = gsp->ep[na]-gsp->eh[kb];
		n_hp[ihp] = na;
		k_hp[ihp] = kb;
		++ihp;
      }
	  
    } }  // end p-h loop
  
  if ((iph > 0) && (iph > n_ph_alloc)) cerr << "\n\n WARNING (Pol_prop_t::Build_phRPA_basis):"
	" p-h arrays are out of bounduaries!!!!!!!!!!!!\n\n";
  if ((ihp > 0) && (ihp > n_hp_alloc)) cerr << "\n\n WARNING (Pol_prop_t::Build_phRPA_basis):"
	" h-p arrays are out of bounduaries!!!!!!!!!!!!\n\n";
  
  n_ph_basis   = iph;
  n_hp_basis   = ihp;
  
  if (iplot > 1) {
    cout << " ...qp-qh basis:\n";
    for(iph =0; iph<n_ph_basis; ++ iph)
      cout << "na,kb =  " << n_ph[iph] << MdSp->MSp_name[gsp->n_clj_p[n_ph[iph]]]
	  << " , " << k_ph[iph] << MdSp->MSp_name[gsp->n_clj_h[k_ph[iph]]]
	  << "   e_ph=" <<  e_ph[iph] << endl;
    cout << endl;
    cout << " ...qh-qp basis:\n";
    for(ihp =0; ihp<n_hp_basis; ++ ihp)
      cout << "ka,nb =  " << k_hp[ihp] << MdSp->MSp_name[gsp->n_clj_h[k_hp[ihp]]]
	  << " , " << n_hp[ihp] << MdSp->MSp_name[gsp->n_clj_p[n_hp[ihp]]]
	  << "   e_hp=" << -e_hp[ihp] << endl;
    cout << endl;
  }
  
  if (iplot) cout << " ...qp-qh and qh-qp bases are done!";
  if (!iALL) {
    cout << "\n ( the selection rule  " << abs(icut) <<"*[n1(a) + n1(b)] ";
    if (icut > 0) cout << ">= " << ntcut;
	else cout << "< "  << 1- ntcut;
	cout << "  was used).\n";
  }
  if (iplot) cout<<"\n\n";

  n_f_basis = (0 > n_ph_basis) ? 0 : n_ph_basis;
  n_b_basis = (0 > n_hp_basis) ? 0 : n_hp_basis;

  return n_ph_basis+n_hp_basis;} // total dimension of the DRPA matrix

//
// NOTE: the routines 'get_Pi_f_indx' and 'get_Pi_b_indx' rely on the
//   ordering of the of na and kb quantum numbers done by keeping kb
//   in the internal loop. These will require chages if the routine
//   'Build_phRPA_basis' is modified.
//


int Pol_prop_t::get_Pi_f_indx(int *imid, int na, int kb) {
  int iseek,iup,idown,iab_mid,nmult;
  nmult = (gsp->n_qp > gsp->n_qh) ? gsp->n_qp : gsp->n_qh;
  // nmult=gsp->n_qh; is also OK...
  iseek = (na*nmult) +kb;
  idown = 0;
  iup   = n_ph_basis-1;
  while(idown <= iup) {
    *imid  = (iup+idown)/2;
    iab_mid = (n_ph[*imid]*nmult) +k_ph[*imid];
    if (iseek == iab_mid) return 0;
    if (iseek < iab_mid) {iup   = (*imid)-1;}
    else {idown = (*imid)+1;}
  }
  
  *imid = -1000;
  return 1;}  // the qp-qh configuration was not found

int Pol_prop_t::get_Pi_b_indx(int *imid, int na, int kb) {
  int iseek,iup,idown,iab_mid,nmult;
  nmult = (gsp->n_qp > gsp->n_qh) ? gsp->n_qp : gsp->n_qh;
  // nmult=gsp->n_qh; is also OK...
  iseek = (na*nmult) +kb;
  idown = 0;
  iup   = n_hp_basis-1;
  while(idown <= iup) {
    *imid  = (iup+idown)/2;
    iab_mid = (n_hp[*imid]*nmult) +k_hp[*imid];
    if (iseek == iab_mid) return 0;
    if (iseek < iab_mid) {iup   = (*imid)-1;}
    else {idown = (*imid)+1;}
  }
  
  *imid = -1000;
  return 1;}  // the qh-qp configuration was not found




static int NSPBMXloc;
static double *MTXf,*MTXb, *MTXf2,*MTXb2;
//static bool shift_imm_sol;//, found_imm_sol;
//static double ef_shift = 0.0;

static bool Make_direct_Ab_matrix, Make_copied_Ab_matrix, Make_DRPA_B_matrices, Free_ph_propagator, calc_T0_rpa;

int Pol_prop_t::Fill_phRPA_mtx(int i_RPAcalc, double *Sab_phRPA, int ndim_vec_alloc, double *Free_poles/*=NULL*/, bool Add_free_poles/*=true*/) {
//
//  i_RPAcalc == 0   full (D)TDA: compute both forward and backward (==time
//                 reversed) amplitudes
//
//  i_RPAcalc == 1  simple (D)TDA (diagonalize A without its time reverse)
//
//  i_RPAcalc == 3  (D)RPA (both forward and backward amplitudes must be
//                     computed in any case)
//
//  i_RPAcalc == 4  generates the free propagator
//

  if ((i_RPAcalc < 0) || (i_RPAcalc > 4) || (2 == i_RPAcalc)) {cout << "Pol_prop_t::Fill_phRPA_mtx: i_RPAcalc==" <<i_RPAcalc<< "is not a valid option \n"; return -100;}


  //
  // set up dimensions:
  //
  // NOTE: "n_f(b)_basis"  are meant to store the dimensions of the X and Y
  //    amplitudes, while "n_ph(hp)_basis" are the bases' dimensions. In 
  //    practice whenever these are both larger zero, they *must* always equal.
  //    Tere are cases in which "n_f(b)_basis" could be zero even if the
  //    2qp/2qh basis is not vanishing, for example when only a part of the 
  //    full RPA matrix is diagonalized (only fw TDA).
  //
  //     NOTE that the "n_f(b)_basis" *must* always be >= zero.
  //     However, if the basis is not allocated, the "n_ph(hp)_basis" can
  //    be negative.
  //
  n_f_basis = (0 > n_ph_basis) ? 0 : n_ph_basis;
  n_b_basis = (0 > n_hp_basis) ? 0 : n_hp_basis;  if (1 == i_RPAcalc) n_b_basis=0;
  ndim_tot = n_f_basis+n_b_basis;
  //
  if (ndim_tot > ndim_vec_alloc) {
	cerr << " ERROR (Pol_prop_t::Fill_phRPA_mtx): value of ndim_vec(="<<ndim_vec_alloc<<") too small. -- STOP\n";
	exit(-100);
  }
  
     //
     // NOTE: 'n_f_basis', 'n_b_basis' and 'ndim_tot' CANNOT be < 0 or the algorithm
     //   below will chrash!

  //cout << "Dimensions of the ph(D)RPA mtx:  fw=" << n_f_basis << " ,   bk="   << n_b_basis << " ,  tot=" << ndim_tot << endl;  

  if (ndim_vec_alloc < 1) {
	 cout << "Pol_prop_t::Solve_phRPA: either there are no solutions for this parial wave or you forgot to initialize the qp-qh bases!\n";
	return -101;
  }


  Make_DRPA_B_matrices = true;
  Make_direct_Ab_matrix = true;
  Make_copied_Ab_matrix = false;
  if   (0 == dT)      {Make_direct_Ab_matrix = false; Make_copied_Ab_matrix = true; }
  if (1 == i_RPAcalc) {Make_direct_Ab_matrix = false; Make_copied_Ab_matrix = false;}
  if (3 != i_RPAcalc) {Make_DRPA_B_matrices = false;}
  if (4 == i_RPAcalc) {Free_ph_propagator = true;  Add_free_poles = true;}
  //
  calc_T0_rpa = (Make_copied_Ab_matrix) && (Make_DRPA_B_matrices);

  NSPBMXloc = gsp->N_SPB_MX;


  // clean up the mtx:
  dptr2 = Sab_phRPA + (ndim_tot-1)*ndim_vec_alloc + ndim_tot;
  for(dptr1=Sab_phRPA; dptr1<dptr2; ++dptr1) (*dptr1)=0.0; //clean up array
  
  //
  // Buid mtx:
  // =========
  
begin_phRPA:

  if (Free_ph_propagator) goto Add_energies;


  //
  //  Af mtx (and Bmtx for dT==0, note that for dT==0 
  // the ph and hp bases are equal!):
  for(nc=0; nc<n_ph_basis; ++nc) {
   //cout << nc << endl;
   ng = n_ph[nc];                             
   ish_g = gsp->n_clj_p[ng];                  
   ig1 = (MdSp->MSp_n[ish_g]   - MdSp->Mor_n);
   ig2 = ig1 + MdSp->MSp_no[ish_g];
   //
   kd = k_ph[nc];
   ish_d = gsp->n_clj_h[kd];
   id1 = (MdSp->MSp_n[ish_d]   - MdSp->Mor_n);
   id2 = id1 + MdSp->MSp_no[ish_d];
   //
   //i_phas = (abs(MdSp->MSp_2j[ish_g]-MdSp->MSp_2j[ish_d])/2)%2;
   //x3 = double(1 - 2*i_phas);
   x3= 1.0; if ( (abs(MdSp->MSp_2j[ish_g]-MdSp->MSp_2j[ish_d])/2)%2 ) x3= -1.0;
   MTXf  = Sab_phRPA +            nc*ndim_vec_alloc;
   MTXb  = Sab_phRPA + (n_f_basis+nc)*ndim_vec_alloc;
   //
   if (4 != i_RPAcalc) // only diag. energies are needed for the free propagator
    for(nr=0; nr<n_ph_basis ; ++nr) {
    na = n_ph[nr];                              
    ish_a = gsp->n_clj_p[na];                   
    ia1 = (MdSp->MSp_n[ish_a]   - MdSp->Mor_n);
    ia2 = ia1 + MdSp->MSp_no[ish_a];
    //
    kb = k_ph[nr];
    ish_b = gsp->n_clj_h[kb];
    ib1 = (MdSp->MSp_n[ish_b]   - MdSp->Mor_n);
    ib2 = ib1 + MdSp->MSp_no[ish_b];

    x1 = 0.0;
    x2 = 0.0;
    fqp_a = gsp->fqp + (na* NSPBMXloc);
    for(ia=ia1; ia<ia2; ++ia) {
     fqh_b = gsp->fqh + (kb* NSPBMXloc);
     for(ib=ib1; ib<ib2; ++ib) {
      xy_r = (*fqp_a) * (* fqh_b);
    
      fqp_g = gsp->fqp + (ng* NSPBMXloc);
      for(ig=ig1; ig<ig2; ++ig) {
       fqh_d = gsp->fqh + (kd* NSPBMXloc);
       for(id=id1; id<id2; ++id) {
        //xy_c =  (*fqp_g) * (* fqh_d)

        // later it will be: 
        //     x1 += xy_r * vph->get(ia,ib,ig,id,Jf,&i1) * (*fqp_g) * (* fqh_d);

            x1 += xy_r * Pandya(ia,ib,ig,id,Jf,vpp) *      (*fqp_g) * (* fqh_d);
        if (calc_T0_rpa) 
            x2 += xy_r * Pandya(ia,ib,id,ig,Jf,vpp) * x3 * (*fqp_g) * (* fqh_d);
  
        ++fqh_d;}
       ++fqp_g;}

      ++fqh_b;}
     ++fqp_a;}

    if (!calc_T0_rpa) x2 = 0.0;

    if (Make_copied_Ab_matrix) {
      MTXb[n_f_basis] = -x1;
      MTXf[n_f_basis] = -x2;
      MTXb[0]         =  x2;   ++MTXb;
    }
    MTXf[0]  =  x1;   ++MTXf;
        
     
    }  // nr loop

	if (NULL != Free_poles) {
	  Free_poles[nc] = e_ph[nc];
	  if (Make_copied_Ab_matrix) Free_poles[n_f_basis+nc] = e_ph[nc];
	}
	
//               *(Sab_phRPA +        nc *(1+ndim_vec_alloc)) += e_ph[nc];
//   if ((0==dT) && (1 != i_RPAcalc))
//              {*(Sab_phRPA + (n_f_basis+nc)*(1+ndim_vec_alloc)) -= e_ph[nc];}
   } // nc loop

  //
  //  note that if dT!=0 then we need to compute the Ab and B matrices
  // separatedly: 
  //
  //cout << "\n\n FIRST DONE \n\n";  

if (Make_DRPA_B_matrices  &&  Make_direct_Ab_matrix) {

  //
  // B mtx:
  for(nc=0; nc<n_hp_basis; ++nc) {
   ng = n_hp[nc];                             
   ish_g = gsp->n_clj_p[ng];                  
   ig1 = (MdSp->MSp_n[ish_g]   - MdSp->Mor_n);
   ig2 = ig1 + MdSp->MSp_no[ish_g];
   //
   kd = k_hp[nc];
   ish_d = gsp->n_clj_h[kd];
   id1 = (MdSp->MSp_n[ish_d]   - MdSp->Mor_n);
   id2 = id1 + MdSp->MSp_no[ish_d];
   //
   //i_phas = (abs(MdSp->MSp_2j[ish_g]-MdSp->MSp_2j[ish_d])/2)%2;
   //x3 = double(1 - 2*i_phas);
   x3= 1.0; if ( (abs(MdSp->MSp_2j[ish_g]-MdSp->MSp_2j[ish_d])/2)%2 ) x3= -1.0;
   MTXb  = Sab_phRPA + (n_f_basis+nc)*ndim_vec_alloc;
   MTXb2 = Sab_phRPA +  n_f_basis+nc;
   //
   for(nr=0; nr<n_ph_basis ; ++nr) {
    na = n_ph[nr];                              
    ish_a = gsp->n_clj_p[na];                   
    ia1 = (MdSp->MSp_n[ish_a]   - MdSp->Mor_n);
    ia2 = ia1 + MdSp->MSp_no[ish_a];
    //
    kb = k_ph[nr];
    ish_b = gsp->n_clj_h[kb];
    ib1 = (MdSp->MSp_n[ish_b]   - MdSp->Mor_n);
    ib2 = ib1 + MdSp->MSp_no[ish_b];
    //
    //x4//x4= 1.0; if ( (abs(MdSp->MSp_2j[ish_a]-MdSp->MSp_2j[ish_b])/2)%2 ) x4= -1.0;

    x2 = 0.0;
    //x4//x4 = 0.0;
    fqp_a = gsp->fqp + (na* NSPBMXloc);
    for(ia=ia1; ia<ia2; ++ia) {
     fqh_b = gsp->fqh + (kb* NSPBMXloc);
     for(ib=ib1; ib<ib2; ++ib) {
      xy_r = (*fqp_a) * (* fqh_b);
    
      fqp_g = gsp->fqp + (ng* NSPBMXloc);
      for(ig=ig1; ig<ig2; ++ig) {
       fqh_d = gsp->fqh + (kd* NSPBMXloc);
       for(id=id1; id<id2; ++id) {
        //xy_c =  (*fqp_g) * (* fqh_d)

         x2 += xy_r * Pandya(ia,ib,id,ig,Jf,vpp) * x3 * (*fqp_g) * (* fqh_d);
  
         //x4//x4 += (*fqp_g) * (* fqh_d) * Pandya(ig,id,ib,ia,Jf,vpp) * x4 * xy_r;
  
        ++fqh_d;}
       ++fqp_g;}

      ++fqh_b;}
     ++fqp_a;}

    (*MTXb ) =  x2;  ++MTXb;
    (*MTXb2) = -x2;  MTXb2+=ndim_vec_alloc;
    //x4//(*MTXb2) = -x4;  MTXb2+=ndim_vec_alloc;

    } // nr loop
   } // nc loop

  //cout << "\n\n SECOND DONE \n\n";  


} // end of 'if (Make_DRPA_B_matrices  &&  Make_direct_Ab_matrix)'
if (Make_direct_Ab_matrix) {


  //
  // Ab mtx:
  for(nc=0; nc<n_hp_basis; ++nc) {
   ng = n_hp[nc];                             
   ish_g = gsp->n_clj_p[ng];                  
   ig1 = (MdSp->MSp_n[ish_g]   - MdSp->Mor_n);
   ig2 = ig1 + MdSp->MSp_no[ish_g];
   //
   kd = k_hp[nc];
   ish_d = gsp->n_clj_h[kd];
   id1 = (MdSp->MSp_n[ish_d]   - MdSp->Mor_n);
   id2 = id1 + MdSp->MSp_no[ish_d];
   //
   MTXf2 = Sab_phRPA + ((nc+n_f_basis)*ndim_vec_alloc + n_f_basis);
   //
   if (4 != i_RPAcalc) // only diag. energies are needed for the free propagator
    for(nr=0; nr<n_hp_basis ; ++nr) {
    na = n_hp[nr];                              
    ish_a = gsp->n_clj_p[na];                   
    ia1 = (MdSp->MSp_n[ish_a]   - MdSp->Mor_n);
    ia2 = ia1 + MdSp->MSp_no[ish_a];
    //
    kb = k_hp[nr];
    ish_b = gsp->n_clj_h[kb];
    ib1 = (MdSp->MSp_n[ish_b]   - MdSp->Mor_n);
    ib2 = ib1 + MdSp->MSp_no[ish_b];

    x1 = 0.0;
    fqp_a = gsp->fqp + (na* NSPBMXloc);
    for(ia=ia1; ia<ia2; ++ia) {
     fqh_b = gsp->fqh + (kb* NSPBMXloc);
     for(ib=ib1; ib<ib2; ++ib) {
      xy_r = (*fqp_a) * (* fqh_b);
    
      fqp_g = gsp->fqp + (ng* NSPBMXloc);
      for(ig=ig1; ig<ig2; ++ig) {
       fqh_d = gsp->fqh + (kd* NSPBMXloc);
       for(id=id1; id<id2; ++id) {
        //xy_c =  (*fqp_g) * (* fqh_d)
         x1 += xy_r * Pandya(ia,ib,ig,id,Jf,vpp) * (*fqp_g) * (* fqh_d);
        ++fqh_d;}
       ++fqp_g;}

      ++fqh_b;}
     ++fqp_a;}

    (*MTXf2) = -x1;  ++MTXf2;
     
    } // nr loop
    if (NULL != Free_poles)	Free_poles[nc+n_f_basis] = e_hp[nc];
   } // nc loop

  //cout << "\n\n THIRD DONE \n\n";  

} // end of 'if (Make_direct_Ab_matrix)'

Add_energies:

  for(nc=0; nc<n_ph_basis; ++nc) {
	if (Add_free_poles) {
	  *(Sab_phRPA +            nc *(1+ndim_vec_alloc)) += e_ph[nc];
	  if ((0==dT) && (1 != i_RPAcalc)) *(Sab_phRPA + (n_f_basis+nc)*(1+ndim_vec_alloc)) -= e_ph[nc];
	}
  }
  //
  if (Make_direct_Ab_matrix) {
	// Ab mtx:
	for(nc=0; nc<n_hp_basis; ++nc) {
	  if (Add_free_poles)	 *(Sab_phRPA + (nc+n_f_basis)*(ndim_vec_alloc+1)) -= e_hp[nc];
	}
  }
  
  return 0;}



//
//
//  Pandya transformation -- it will be moved to the Vph class...
//

static double vph = 0.0;
static int j2g,j2d,J2,JJ2,J2f,J2str,J2end; //j2a,j2b,
static int i1,i2;


double Pandya(int ia, int ib, int ig, int id, int Jf, VppInt_t *vpp) {
  
  ModSpace_t *MdSp = vpp->MdSp;

  if ((MdSp->MSp_ip[MdSp->Mor_sh[ia]]+MdSp->MSp_ip[MdSp->Mor_sh[ib]]
      +MdSp->MSp_ip[MdSp->Mor_sh[ig]]+MdSp->MSp_ip[MdSp->Mor_sh[id]])%2)
       { cout << "\n\n Warning: Pandya tranf. routine called with wrong parity.\n\n";
         cout << ia << "  " << ib << "  " << ig << "  " << id << "  " << Jf <<endl;
         return 0.0;
         }

  j2a = MdSp->MSp_2j[MdSp->Mor_sh[ia]];
  j2b = MdSp->MSp_2j[MdSp->Mor_sh[ib]];
  j2g = MdSp->MSp_2j[MdSp->Mor_sh[ig]];
  j2d = MdSp->MSp_2j[MdSp->Mor_sh[id]];

  vph = 0.0;
  J2f = Jf+Jf;
  i1 = abs(j2a-j2d);
  i2 = abs(j2b-j2g);
  J2str = (i1 > i2) ? i1 : i2;
  i1 = j2a+j2d;
  i2 = j2b+j2g;
  J2end = (i1 < i2) ? i1 : i2;
  for(JJ2=J2str; JJ2<=J2end ; JJ2+=2) {
   J2 = JJ2 / 2;

    //vph -= double(JJ2+1) * vpp->get(id,ia,ib,ig,J2,&i1) *
    vph -= double(JJ2+1) * vpp->get(ia,id,ig,ib,J2,&i1) *
                           am.s6j(j2a,j2b,J2f,j2g,j2d,JJ2);

   }

  if (ia==id) vph *= SQR2;
  if (ig==ib) vph *= SQR2;
  return vph;}
