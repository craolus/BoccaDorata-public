//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//
//  the following is to be placed in a second file and be adapted to the new
//  Wide_self-enerrgy class (...and tested!!).  It will eventually replace
// the old  Self-enerrgy class.
//
//


#include <cstdlib>
//#include <cstdio>
#include <cmath>
#include <iostream>
// #include <fstream>
// #include <iomanip>
using namespace std;


#include "BcDor-SlfEn_classes.hh"
#include "BcDor-Ang_momenta.hh"


#define SQR2 1.4142135623730950488016887


static bool test1, test2;
static int ia,ib,ig;
static int na,nb,ng, nst , nt;
static int Jbg,Jab,nJab,isum;

static bool iALL;

void Wide_self_energy_t::Count_2p1h_2h1p_poles(int* npl_fw, int* npl_bk, 
                                          int icut, int ntcut) {
                                          //defaults: icut==0, ntcut==0

  iALL = true;
  if ( 0 != icut ) iALL=false;
  if ( 0 >  icut ) ntcut = 1- ntcut;

  if (MdSp->Jmax < 0) MdSp->Calculate_Jch_bounds();

  *npl_fw = 0;
  *npl_bk = 0;
  for(Jab=0 ; Jab<=MdSp->Jmax; ++Jab) {
   for(ia=0  ; ia<MdSp->nsubsh; ++ia) {
    for(ib=ia ; ib<MdSp->nsubsh; ++ib) {
// */ TST-amTri //
     test1 = (MdSp->MSp_2j[ia]+MdSp->MSp_2j[ib] <    2*Jab) ||
             (MdSp->MSp_2j[ib]+    2*Jab        < MdSp->MSp_2j[ia]) ||
             (    2*Jab       +MdSp->MSp_2j[ia] < MdSp->MSp_2j[ib]);
     test2 = am.TriIneq(MdSp->MSp_2j[ia] , MdSp->MSp_2j[ib] , 2*Jab);
     if ( (test1 && (!test2) ) || ((!test1) && test2) ) {cerr << "ERROR w/ a.TriIneq SlfEn ----- stop!\n";exit(10);}
// TST-amTri *//
     if ( am.TriIneq(MdSp->MSp_2j[ia] , MdSp->MSp_2j[ib] , 2*Jab) ) continue;

     nJab = ( (MdSp->MSp_2j[ia]+MdSp->MSp_2j[ib])/2 + Jab + 1)%2;

     for(ig=0  ; ig<MdSp->nsubsh; ++ig) {

       if  (MdSp->MSp_ch[ia]+MdSp->MSp_ch[ib]-MdSp->MSp_ch[ig] != MdSp->MSp_ch[I_SH]) continue;
       if ((MdSp->MSp_ip[ia]+MdSp->MSp_ip[ib]+
            MdSp->MSp_ip[ig]+MdSp->MSp_ip[I_SH])%2 != 0)  continue;
// */ TST-amTri //
     test1 = (MdSp->MSp_2j[I_SH] +     2*Jab          < MdSp->MSp_2j[ig]  ) ||
             (      2*Jab        + MdSp->MSp_2j[ig]   < MdSp->MSp_2j[I_SH]) ||
             (MdSp->MSp_2j[ig]   + MdSp->MSp_2j[I_SH] < 2*Jab);
     test2 = am.TriIneq(MdSp->MSp_2j[ig] , MdSp->MSp_2j[I_SH] , 2*Jab);
     if ( (test1 && (!test2) ) || ((!test1) && test2) ) {cerr << "ERROR w/ a.TriIneq SlfEn ----- stop!\n";exit(10);}
// TST-amTri *//
       if ( am.TriIneq(MdSp->MSp_2j[ig] , MdSp->MSp_2j[I_SH] , 2*Jab) ) continue;

       if (iALL) {
         na=gin->Sp_clj_np[ia];
         nb=gin->Sp_clj_np[ib];
         ng=gin->Sp_clj_nh[ig];
         isum = na*nb*ng;
         if (ia == ib) {nb = na+1; if (1 == nJab) nb -=2; isum=(na*nb/2)*ng;}
         (*npl_fw) += isum;

         na=gin->Sp_clj_nh[ia];
         nb=gin->Sp_clj_nh[ib];
         ng=gin->Sp_clj_np[ig];
         isum = na*nb*ng;
         if (ia == ib) {nb = na+1; if (1 == nJab) nb -=2; isum=(na*nb/2)*ng;}
         (*npl_bk) += isum;
       } else {
         // here we are forced to loop over all qp and qh indices...
         for(na=0;   na<gin->Sp_clj_np[ia]; ++na) {
                        nst=0; if (ia == ib) nst = na + nJab;
         for(nb=nst; nb<gin->Sp_clj_np[ib]; ++nb) {
         for(ng=0;   ng<gin->Sp_clj_nh[ig]; ++ng) {
           nt = gin->Sp_clj_n1p[ia][na] + gin->Sp_clj_n1p[ib][nb] + gin->Sp_clj_n1h[ig][ng];
           if (nt*icut >= ntcut) ++(*npl_fw);
         } } }
         
         for(na=0;   na<gin->Sp_clj_nh[ia]; ++na) {
                        nst=0; if (ia == ib) nst = na + nJab;
         for(nb=nst; nb<gin->Sp_clj_nh[ib]; ++nb) {
         for(ng=0;   ng<gin->Sp_clj_np[ig]; ++ng) {
           nt = gin->Sp_clj_n1h[ia][na] + gin->Sp_clj_n1h[ib][nb] + gin->Sp_clj_n1p[ig][ng];
           if (nt*icut >= ntcut) ++(*npl_bk);
         } } }
       }

     }
    }
   }
  }
  return;
  }


static int ifw,ibk,indx;
static int     nash,nbsh,ngsh,ndsh;
static int      iia, iib, iig, iid;
static double  *epa, *epb, *ehg;
static int    *n1pa,*n1pb,*n1hg;
static double *fqpa,*fqpb,*fqhg;
static double *Sig_fw_ptr,*Sig_bk_ptr;
static double x1,x2,x3,x4,x5,xJab;
static int    nJgd;
static int    N_SPB_MX_loc;


void Wide_self_energy_t::Build_Dynamic_SlfEn_2ndOrder_slow(int* nfw_tot, int* nbk_tot,
                int inew/*=1*/, int icut/*=0*/, int ntcut/*=0*/) {
          //
          // Defaults: inew==1, icut==0, ntcut==0
          // inew == 0, do not (re)allocate the s.e.
          //      == 1, (re)allocate the largest amount of mem. possible
          //      == 2 or otherwise, (re)allocate only the memory
          //          needed (according to icut & ntcut) 
          //
          // icut == 0,  include all posible poles
          //       > 0,  choose poles according to (n1+n2+n3)*|icut|  >= ntcut 
          //       < 0,  choose poles according to (n1+n2+n3)*|icut|  <  ntcut 
            

  if (0 != inew) {
    if (1==inew) {Count_2p1h_2h1p_poles(&ifw , &ibk);}
      else {Count_2p1h_2h1p_poles(&ifw , &ibk, icut, ntcut);}
    free_dyn_mem();
    Allocate_Dynamic_SlfEn(ifw , ibk);
  }

  //  do this now and not above because `Count_2p1h_2h1p_poles'
  // could change it:
  iALL = true;
  if ( 0 != icut ) iALL=false;
  if ( 0 >  icut ) ntcut = 1- ntcut;


  N_SPB_MX_loc = gin->N_SPB_MX;

  *nfw_tot = 0;
  *nbk_tot = 0;
  ifw = N_PLS_fw;
  ibk = N_PLS_bk;
  Sig_fw_ptr = Sig_dyn_fw + ifw*NTOT_DIM;
  Sig_bk_ptr = Sig_dyn_bk + ibk*NTOT_DIM;
  for(Jab=0 ; Jab<=MdSp->Jmax; ++Jab) {
   for(ia=0  ; ia<MdSp->nsubsh; ++ia) {
    for(ib=ia ; ib<MdSp->nsubsh; ++ib) {
     if ( am.TriIneq(MdSp->MSp_2j[ia] , MdSp->MSp_2j[ib] , 2*Jab) ) continue;

     nJab = abs( (MdSp->MSp_2j[ia]+MdSp->MSp_2j[ib])/2 - Jab + 1)%2;
     xJab = double(1 - 2*nJab);//1.0; if (0 != nJab) xJab=-1.0;

     for(ig=0  ; ig<MdSp->nsubsh; ++ig) {
       nJgd = abs( (MdSp->MSp_2j[ig]+MdSp->MSp_2j[I_SH])/2 - Jab + 1)%2;

       if  (MdSp->MSp_ch[ia]+MdSp->MSp_ch[ib]-MdSp->MSp_ch[ig] != MdSp->MSp_ch[I_SH]) continue;
       if ((MdSp->MSp_ip[ia]+MdSp->MSp_ip[ib]+
            MdSp->MSp_ip[ig]+MdSp->MSp_ip[I_SH])%2 != 0)  continue;
       if (am.TriIneq(MdSp->MSp_2j[ig] , MdSp->MSp_2j[I_SH] , 2*Jab) ) continue;


//t2p1h:
       if (ifw >= N_PLS_FW_alloc) goto t2h1p;

        epa = gin->Sp_clj_ep[ia];
       n1pa = gin->Sp_clj_n1p[ia];
       fqpa = gin->Sp_clj_fqp[ia];
       for(na=0;   na<gin->Sp_clj_np[ia]; ++na) {
                        nst=0; if (ia == ib) nst = na + nJab;
                         epb = (gin->Sp_clj_ep[ib] )+nst;
                        n1pb = (gin->Sp_clj_n1p[ib])+nst;
                        fqpb = (gin->Sp_clj_fqp[ib])+(nst*N_SPB_MX_loc);
       for(nb=nst; nb<gin->Sp_clj_np[ib]; ++nb) {
                         ehg = (gin->Sp_clj_eh[ig]);
                        n1hg = (gin->Sp_clj_n1h[ig]);
                        fqhg = (gin->Sp_clj_fqh[ig]);
       for(ng=0;   ng<gin->Sp_clj_nh[ig]; ++ng) {
         if (ifw >= N_PLS_FW_alloc) goto t2h1p;
         if ((iALL) || ( ((*n1pa)+(*n1pb)+(*n1hg))*icut >= ntcut)) {
           Sig_fw_ptr[0]=(*epa)+(*epb)-(*ehg);
           Sig_fw_ptr[1]=0.0;

           x3 = double(Jab+Jab+1)/double(1 + MdSp->MSp_2j[I_SH]);
           x3 = sqrt(x3);
           if ((ia == ib)&&(na == nb)) x3 /= SQR2;

           for(ndsh=0  ; ndsh<MdSp->MSp_no[I_SH]; ++ndsh) {
             iid = (MdSp->MSp_n[I_SH]   - MdSp->Mor_n) + ndsh;

             x1 = 0.0;
             for(nash=0  ; nash<MdSp->MSp_no[ia]; ++nash) {nst=0; if (ia == ib) nst = nash + nJab;
             for(nbsh=nst; nbsh<MdSp->MSp_no[ib]; ++nbsh) {
              x4 = x3;
              x2 = fqpa[nash] * fqpb[nbsh];
              if (ia == ib) {
                 x2 += xJab * fqpa[nbsh] * fqpb[nash];
                 if (nash==nbsh) x4 /= SQR2;
                 }
             for(ngsh=0  ; ngsh<MdSp->MSp_no[ig]; ++ngsh) {
               if ((ig==I_SH)&&(ngsh==ndsh)&&(0!=nJgd)) continue;
               x5 = x4; if ((ig==I_SH)&&(ngsh==ndsh)) x5 *= SQR2;
               iia = (MdSp->MSp_n[ia]   - MdSp->Mor_n) + nash;
               iib = (MdSp->MSp_n[ib]   - MdSp->Mor_n) + nbsh;
               iig = (MdSp->MSp_n[ig]   - MdSp->Mor_n) + ngsh;
               //
               x1 +=  x5 * x2 * fqhg[ngsh] * Vpp->get(iia,iib,iid,iig,Jab,&indx);
             
               } } }
                 
             Sig_fw_ptr[NTOT_DEN+ndsh] = x1;
             Sig_fw_ptr[1] += x1*x1;
             }

           ++ifw; ++(*nfw_tot);
           Sig_fw_ptr += NTOT_DIM;
         }

       ++ehg; ++n1hg; fqhg+=N_SPB_MX_loc;}  // ng
       ++epb; ++n1pb; fqpb+=N_SPB_MX_loc;}  // nb
       ++epa; ++n1pa; fqpa+=N_SPB_MX_loc;}  // na


t2h1p:  if (ibk >= N_PLS_BK_alloc) continue;

        epa = gin->Sp_clj_eh[ia];
       n1pa = gin->Sp_clj_n1h[ia];
       fqpa = gin->Sp_clj_fqh[ia];
       for(na=0;   na<gin->Sp_clj_nh[ia]; ++na) {
                        nst=0; if (ia == ib) nst = na + nJab;
                         epb = (gin->Sp_clj_eh[ib] )+nst;
                        n1pb = (gin->Sp_clj_n1h[ib])+nst;
                        fqpb = (gin->Sp_clj_fqh[ib])+(nst*N_SPB_MX_loc);
       for(nb=nst; nb<gin->Sp_clj_nh[ib]; ++nb) {
                         ehg = (gin->Sp_clj_ep[ig]);
                        n1hg = (gin->Sp_clj_n1p[ig]);
                        fqhg = (gin->Sp_clj_fqp[ig]);
       for(ng=0;   ng<gin->Sp_clj_np[ig]; ++ng) {
         if (ibk >= N_PLS_BK_alloc) goto cont;
         if ((iALL) || ( ((*n1pa)+(*n1pb)+(*n1hg))*icut >= ntcut)) {
           Sig_bk_ptr[0]=(*epa)+(*epb)-(*ehg);
           Sig_bk_ptr[1]=0.0;

           x3 = double(Jab+Jab+1)/double(1 + MdSp->MSp_2j[I_SH]);
           x3 = sqrt(x3);
           if ((ia == ib)&&(na == nb)) x3 /= SQR2;

           for(ndsh=0  ; ndsh<MdSp->MSp_no[I_SH]; ++ndsh) {
             iid = (MdSp->MSp_n[I_SH]   - MdSp->Mor_n) + ndsh;

             x1 = 0.0;
             for(nash=0  ; nash<MdSp->MSp_no[ia]; ++nash) {nst=0; if (ia == ib) nst = nash + nJab;
             for(nbsh=nst; nbsh<MdSp->MSp_no[ib]; ++nbsh) {
              x4 = x3;
              x2 = fqpa[nash] * fqpb[nbsh];
              if (ia == ib) {
                 x2 += xJab * fqpa[nbsh] * fqpb[nash];
                 if (nash==nbsh) x4 /= SQR2;
                 }
             for(ngsh=0  ; ngsh<MdSp->MSp_no[ig]; ++ngsh) {
               if ((ig==I_SH)&&(ngsh==ndsh)&&(0!=nJgd)) continue;
               x5 = x4; if ((ig==I_SH)&&(ngsh==ndsh)) x5 *= SQR2;
               iia = (MdSp->MSp_n[ia]   - MdSp->Mor_n) + nash;
               iib = (MdSp->MSp_n[ib]   - MdSp->Mor_n) + nbsh;
               iig = (MdSp->MSp_n[ig]   - MdSp->Mor_n) + ngsh;
               //
               x1 +=  x5 * x2 * fqhg[ngsh] * Vpp->get(iia,iib,iid,iig,Jab,&indx);
             
               } } }
                 
             Sig_bk_ptr[NTOT_DEN+ndsh] = x1;
             Sig_bk_ptr[1] += x1*x1;
             }
             

           ++ibk; ++(*nbk_tot);
           Sig_bk_ptr += NTOT_DIM;
         }
       ++ehg; ++n1hg; fqhg+=N_SPB_MX_loc;}
       ++epb; ++n1pb; fqpb+=N_SPB_MX_loc;}
       ++epa; ++n1pa; fqpa+=N_SPB_MX_loc;}

cont: ;
     }
    }
   }
  }

  if ( (ifw - N_PLS_fw != (*nfw_tot) ) || (ibk - N_PLS_bk != (*nbk_tot) ) ) {
     cerr << " Some problems occured while counting the # of new self-energy poles...! \n";
     exit(20);
  }

  cout << endl << *nfw_tot << " fw poles added to the " << N_PLS_fw << " existing ones.\n";
  cout         << *nbk_tot << " bk poles added to the " << N_PLS_bk << " existing ones.\n";

  N_PLS_fw = ifw;
  N_PLS_bk = ibk;

  return;
  }



void Wide_self_energy_t::Build_Dynamic_SlfEn_2ndOrder(int* nfw_tot, int* nbk_tot,
                int inew/*=1*/, int icut/*=0*/, int ntcut/*=0*/) {
          //
          // Defaults: inew==1, icut==0, ntcut==0
          // inew == 0, do not (re)allocate the s.e.
          //      == 1, (re)allocate the largest amount of mem. possible
          //      == 2 or otherwise, (re)allocate only the memory
          //          needed (according to icut & ntcut) 
          //
          // icut == 0,  include all posible poles
          //       > 0,  choose poles according to (n1+n2+n3)*|icut|  >= ntcut 
          //       < 0,  choose poles according to (n1+n2+n3)*|icut|  <  ntcut 


  if (0 != inew) {
    if (1==inew) {Count_2p1h_2h1p_poles(&ifw , &ibk);}
      else {Count_2p1h_2h1p_poles(&ifw , &ibk, icut, ntcut);}
    free_dyn_mem();
    Allocate_Dynamic_SlfEn(ifw , ibk);
  }

  //  do this now and not above because `Count_2p1h_2h1p_poles'
  // could change it:
  iALL = true;
  if ( 0 != icut ) iALL=false;
  if ( 0 >  icut ) ntcut = 1- ntcut;


/*
  i =0;
  for(ifw=0; ifw<N_PLS_fw; ++ifw) {
    Sig_dyn_fw[i] = double(ifw); ++i;
   for(nr=1; nr<NTOT_DIM; ++nr) {
    Sig_dyn_fw[i] = cos(double(ifw*nr)/100.);
     ++i; }  }

  i =0;
  for(ibk=0; ibk<N_PLS_bk; ++ibk) {
    Sig_dyn_bk[i] = -15. - double(ibk); ++i;
   for(nr=1; nr<NTOT_DIM; ++nr) {
    Sig_dyn_bk[i] =  cos(double(ibk*nr)/100.);
     ++i; }  }
*/


  N_SPB_MX_loc = gin->N_SPB_MX;

  double *Vpp_dptr, *dptr_v;
  int iVpp_mx  = N_SPB_MX_loc*N_SPB_MX_loc*N_SPB_MX_loc*N_SPB_MX_loc;
  Vpp_dptr = new double[iVpp_mx];

  *nfw_tot = 0;
  *nbk_tot = 0;
  ifw = N_PLS_fw;
  ibk = N_PLS_bk;
  Sig_fw_ptr = Sig_dyn_fw + ifw*NTOT_DIM;
  Sig_bk_ptr = Sig_dyn_bk + ibk*NTOT_DIM;
  for(Jab=0 ; Jab<=MdSp->Jmax; ++Jab) {
   for(ia=0  ; ia<MdSp->nsubsh; ++ia) {
    for(ib=ia ; ib<MdSp->nsubsh; ++ib) {
     if ( am.TriIneq(MdSp->MSp_2j[ia] , MdSp->MSp_2j[ib] , 2*Jab) ) continue;

     nJab = abs( (MdSp->MSp_2j[ia]+MdSp->MSp_2j[ib])/2 - Jab + 1)%2;
     xJab = double(1 - 2*nJab);//1.0; if (0 != nJab) xJab=-1.0;

     for(ig=0  ; ig<MdSp->nsubsh; ++ig) {
       nJgd = abs( (MdSp->MSp_2j[ig]+MdSp->MSp_2j[I_SH])/2 - Jab + 1)%2;

       if  (MdSp->MSp_ch[ia]+MdSp->MSp_ch[ib]-MdSp->MSp_ch[ig] != MdSp->MSp_ch[I_SH]) continue;
       if ((MdSp->MSp_ip[ia]+MdSp->MSp_ip[ib]+
            MdSp->MSp_ip[ig]+MdSp->MSp_ip[I_SH])%2 != 0)  continue;
       if (am.TriIneq(MdSp->MSp_2j[ig] , MdSp->MSp_2j[I_SH] , 2*Jab) ) continue;

       //
       //
       // Precompute the Vpp mtx el part:
       x3 = double(Jab+Jab+1)/double(1 + MdSp->MSp_2j[I_SH]);
       x3 = sqrt(x3);

       dptr_v = Vpp_dptr;
       for(ndsh=0  ; ndsh<MdSp->MSp_no[I_SH]; ++ndsh) {
         iid = (MdSp->MSp_n[I_SH]   - MdSp->Mor_n) + ndsh;
         for(nash=0  ; nash<MdSp->MSp_no[ia]; ++nash) {nst=0; if (ia == ib) nst = nash + nJab;
           iia = (MdSp->MSp_n[ia]   - MdSp->Mor_n) + nash;
           for(nbsh=nst; nbsh<MdSp->MSp_no[ib]; ++nbsh) {
             iib = (MdSp->MSp_n[ib]   - MdSp->Mor_n) + nbsh;
             for(ngsh=0  ; ngsh<MdSp->MSp_no[ig]; ++ngsh) {
               (*dptr_v) = 0.0;
               iig = (MdSp->MSp_n[ig]   - MdSp->Mor_n) + ngsh;
               if ((iig==iid)&&(nJgd)) continue;
               (*dptr_v) = x3 * Vpp->get(iia,iib,iid,iig,Jab,&indx);
               if (iig==iid) (*dptr_v) *= SQR2;
               ++dptr_v;
           } } } }
         if ((dptr_v - Vpp_dptr) > iVpp_mx) {
            cerr << "\n Trouble wit pointers dptr_v,Vpp_dptr: " << dptr_v << "  " << Vpp_dptr
                 << "   "  << dptr_v - Vpp_dptr <<"  "<< iVpp_mx << endl;
            exit(100);
         }



//t2p1h:
        if (ifw >= N_PLS_FW_alloc) goto t2h1p;

        epa = gin->Sp_clj_ep[ia];
       n1pa = gin->Sp_clj_n1p[ia];
       fqpa = gin->Sp_clj_fqp[ia];
       for(na=0;   na<gin->Sp_clj_np[ia]; ++na) {
                        nst=0; if (ia == ib) nst = na + nJab;
                         epb = (gin->Sp_clj_ep[ib] )+nst;
                        n1pb = (gin->Sp_clj_n1p[ib])+nst;
                        fqpb = (gin->Sp_clj_fqp[ib])+(nst*N_SPB_MX_loc);
       for(nb=nst; nb<gin->Sp_clj_np[ib]; ++nb) {
                         ehg = (gin->Sp_clj_eh[ig]);
                        n1hg = (gin->Sp_clj_n1h[ig]);
                        fqhg = (gin->Sp_clj_fqh[ig]);
       for(ng=0;   ng<gin->Sp_clj_nh[ig]; ++ng) {
         if (ifw >= N_PLS_FW_alloc) goto t2h1p;
         if ((iALL) || ( ((*n1pa)+(*n1pb)+(*n1hg))*icut >= ntcut)) {
           Sig_fw_ptr[0]=(*epa)+(*epb)-(*ehg);
           Sig_fw_ptr[1]=0.0;

           x4 = 1.0;
           if ((ia == ib)&&(na == nb)) x4 /= SQR2;

           dptr_v = Vpp_dptr;
           for(ndsh=0  ; ndsh<MdSp->MSp_no[I_SH]; ++ndsh) {
             iid = (MdSp->MSp_n[I_SH]   - MdSp->Mor_n) + ndsh;

             x1 = 0.0;
             for(nash=0  ; nash<MdSp->MSp_no[ia]; ++nash) {nst=0; if (ia == ib) nst = nash + nJab;
             for(nbsh=nst; nbsh<MdSp->MSp_no[ib]; ++nbsh) {
              x2 = fqpa[nash] * fqpb[nbsh];
              if (ia == ib) {
                 x2 += xJab * fqpa[nbsh] * fqpb[nash];
                 if (nash==nbsh) x2 /= SQR2;
                 }
             for(ngsh=0  ; ngsh<MdSp->MSp_no[ig]; ++ngsh) {
               iig = (MdSp->MSp_n[ig]   - MdSp->Mor_n) + ngsh;
               if ((iig==iid)&&(nJgd)) continue;

               x1 +=  (*dptr_v) * x2 * fqhg[ngsh];

               ++dptr_v;
             } } }
             x1 *= x4;
                 
             Sig_fw_ptr[NTOT_DEN+ndsh] = x1;
             Sig_fw_ptr[1] += x1*x1;
             }

           ++ifw; ++(*nfw_tot);
           Sig_fw_ptr += NTOT_DIM;
         }

       ++ehg; ++n1hg; fqhg+=N_SPB_MX_loc;}  // ng
       ++epb; ++n1pb; fqpb+=N_SPB_MX_loc;}  // nb
       ++epa; ++n1pa; fqpa+=N_SPB_MX_loc;}  // na


t2h1p:  if (ibk >= N_PLS_BK_alloc) continue; // or 'goto cont';

        epa = gin->Sp_clj_eh[ia];
       n1pa = gin->Sp_clj_n1h[ia];
       fqpa = gin->Sp_clj_fqh[ia];
       for(na=0;   na<gin->Sp_clj_nh[ia]; ++na) {
                        nst=0; if (ia == ib) nst = na + nJab;
                         epb = (gin->Sp_clj_eh[ib] )+nst;
                        n1pb = (gin->Sp_clj_n1h[ib])+nst;
                        fqpb = (gin->Sp_clj_fqh[ib])+(nst*N_SPB_MX_loc);
       for(nb=nst; nb<gin->Sp_clj_nh[ib]; ++nb) {
                         ehg = (gin->Sp_clj_ep[ig]);
                        n1hg = (gin->Sp_clj_n1p[ig]);
                        fqhg = (gin->Sp_clj_fqp[ig]);
       for(ng=0;   ng<gin->Sp_clj_np[ig]; ++ng) {
         if (ibk >= N_PLS_BK_alloc) goto cont;
         if ((iALL) || ( ((*n1pa)+(*n1pb)+(*n1hg))*icut >= ntcut)) {
           Sig_bk_ptr[0]=(*epa)+(*epb)-(*ehg);
           Sig_bk_ptr[1]=0.0;

           x4 = 1.0;
           if ((ia == ib)&&(na == nb)) x4 /= SQR2;

           dptr_v = Vpp_dptr;
           for(ndsh=0  ; ndsh<MdSp->MSp_no[I_SH]; ++ndsh) {
             iid = (MdSp->MSp_n[I_SH]   - MdSp->Mor_n) + ndsh;

             x1 = 0.0;
             for(nash=0  ; nash<MdSp->MSp_no[ia]; ++nash) {nst=0; if (ia == ib) nst = nash + nJab;
             for(nbsh=nst; nbsh<MdSp->MSp_no[ib]; ++nbsh) {
              x2 = fqpa[nash] * fqpb[nbsh];
              if (ia == ib) {
                 x2 += xJab * fqpa[nbsh] * fqpb[nash];
                 if (nash==nbsh) x2 /= SQR2;
                 }
             for(ngsh=0  ; ngsh<MdSp->MSp_no[ig]; ++ngsh) {
               iig = (MdSp->MSp_n[ig]   - MdSp->Mor_n) + ngsh;
               if ((iig==iid)&&(nJgd)) continue;

               x1 +=  (*dptr_v) * x2 * fqhg[ngsh];

               ++dptr_v;
             } } }
             x1 *= x4;
                 
             Sig_bk_ptr[NTOT_DEN+ndsh] = x1;
             Sig_bk_ptr[1] += x1*x1;
             }
             

           ++ibk; ++(*nbk_tot);
           Sig_bk_ptr += NTOT_DIM;
         }
       ++ehg; ++n1hg; fqhg+=N_SPB_MX_loc;}
       ++epb; ++n1pb; fqpb+=N_SPB_MX_loc;}
       ++epa; ++n1pa; fqpa+=N_SPB_MX_loc;}

cont: ;
     }
    }
   }
  }

  if ( (ifw - N_PLS_fw != (*nfw_tot) ) || (ibk - N_PLS_bk != (*nbk_tot) ) ) {
     cerr << " Some problems occured while counting the # of new self-energy poles...! \n";
     exit(20);
  }

  cout << endl << *nfw_tot << " fw poles added to the " << N_PLS_fw << " existing ones.\n";
  cout         << *nbk_tot << " bk poles added to the " << N_PLS_bk << " existing ones.\n";

  N_PLS_fw = ifw;
  N_PLS_bk = ibk;

  delete [] Vpp_dptr;

  return;
  }
