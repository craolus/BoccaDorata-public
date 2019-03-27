//
//  Routines for the 3-j, 6-j and 9-j coefficients
//
//   These were taken from the CENS package (of Hjorth-Jensen):
//
//            !    Program block heff-modules.f90
//            !
//            !    Author:   Morten Hjorth-Jensen
//            !    ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
//            !    E-MAIL:   morten.hjorth-jensen@fys.uio.no
//            !    LANGUAGE: F90/F95
//            !
//
//  and translated into C++ by Carlo Barbieri (C.Barbieri@gsi.de)
//
//  August 20th, 2006.
//

#include <cstdlib>
//#include <cstdio>
#include <cmath>

#include <iostream>
//using namespace std;

#include "BcDor-Ang_momenta.hh"

//
// Initialize global variables:
//
// constructor
static int    l, k, i;//, j(sph. hamonics routine);
//static double a, sq_pi, fj, tfj, fk;// (mosh & sph.har. routine)
//
// 3-j symbols
static int ja, jb, jc, ma, mb, mc, la, lb, lc, lt, k0, k1;//ja2, jb2, jc2,  i, ld, k, ip,
static double x, fn;
//
// 6-j symbols
static int mt,na,nb,nc,ka,kb,kc,l0,l1;//ja,jb,jc,la,lb,lc,ma,mb,mc,i,l,
static double fs, fss;// x, 
//
// 9-j symbols
static int ld,le,lf,
       lg,me,mf,mg,ne,nf,ng,lx,mx,nx,jsi,jsf, js,
       ly, my,ny,m0,n0,m1,n1,m,n;
       // ja,jb,jc, i,la, ma,mc,na,nb,lb, lc,mb, l,l0,l1,
static double   fd, u, y, z, ud; // x, fn, fs, p, ps, 




ang_mom_functions::ang_mom_functions() {
// commons_to_angmom:
// ==================
//   factorials for 3j,6j and 9j symbols            
//  for moshinsky trans brackets and for           
//  vector brackets                                
//
    //    3j, 6j and 9j symbols
    //----------------------------
    //
    // used to check triangular inequalities...
    for(i=0; i<SIZEK; ++i) kha[i]=1;
    l = SIZEK/2;
    kh = kha + l;
    for(i=0; i<(SIZEK-l); i+=2) kh[i]=0;
    //
    //
    // NOTE:              ( ia )
    //    now q[ia][ib] = (    )
    //                    ( ib )
    // remember that it must be 0 <= ib <= ia <SIZEJ, or you'll get junk!
    for(l=0; l<SIZEJ; ++l) { 
      q[l][0]=1.0;
      q[l][l]=1.0;
      }
    for(l=1; l<(SIZEJ-1); ++l)
     for(k=1; k<=l;        ++k)
         q[l+1][k]=q[l][k-1]+q[l][k];

    //cout << "\n\nANGULAR MOMENTA INITIALIZED!!\n\n";

  return;
  }


double ang_mom_functions::ClG(int j_a, int j_b, int j_c,
                              int m_a, int m_b, int m_c) {
//
//     Calculates Clebsch-Gordan coefficients
//
  x = s3j(j_a, j_b, j_c, m_a, m_b, -m_c);
  //
  //  NOTE that the following correct a very old version of this routine that
  // had a wrong phase. The bug was found by A. Cipollone in year 2011-12 but
  // then I have found it copied in much later versions. THE LINE OF CODE
  // BELOW HERE IS THE CORRECT ONE: beware that if in the future you don't find
  // this note you may be running the bugged version [CB 8.04.15].
  return ( x * double( 1 - ( abs(j_a-j_b+m_c)%4 ) ) * sqrt(double(j_c+1)) );
  }



double ang_mom_functions::s3j(int j_a, int j_b, int j_c,
                              int m_a, int m_b, int m_c) {
//
//     Calculates 3j-symbols           
//
  //
  // check triangluar inequalities etc.
  ja=(j_a+m_a)/2;
  ma=(j_a-m_a)/2;
  jb=(j_b+m_b)/2;
  mb=(j_b-m_b)/2;
  jc=(j_c+m_c)/2;
  mc=(j_c-m_c)/2;
  la=(j_b+j_c-j_a)/2;
  lb=(j_c+j_a-j_b)/2;
  lc=(j_a+j_b-j_c)/2;
  lt=(j_a+j_b+j_c)/2;
  //ld=MIN(ja,jb,jc,ma,mb,mc,la,lb,lc)
  //IF(((m_a+m_b+m_c) > 0).AND.(ld <= 0)) RETURN
  if (0 != (m_a+m_b+m_c) ) return 0.0;   // it should be (m_a+m_b+m_b != 0) ==> s3j=0. ???
  if ((ja < 0) || (jb < 0) || (jc < 0)) return 0.0;
  if ((ma < 0) || (mb < 0) || (mc < 0)) return 0.0;
  if ((la < 0) || (lb < 0) || (lc < 0)) return 0.0;
  //ja2=j_a+m_a;   //note that ja2,jb2,jc2 are all >= 0 at this point
  //jb2=j_b+m_b;
  //jc2=j_c+m_c;
  //i=ja2+jb2+jc2-ja2/2*2-jb2/2*2-jc2/2*2
  //IF(i != 0) RETURN
  i = ((j_a+m_a)%2) + ((j_b+m_b)%2) + ((j_c+m_c)%2);
  if (i) return 0.0;

  //
  // now compute the thing:
  fn=q[ja+ma][lc]*q[jb+mb][lc]/(q[lt][jc+mc]*q[lt+1][1] 
       *q[ja+ma][ja]*q[jb+mb][jb]*q[jc+mc][jc]);
  k0 = (0  > lc-ja) ? 0  : lc-ja;
  k0 = (k0 > lc-mb) ? k0 : lc-mb;
  k1 = (lc < ma) ? lc : ma;
  k1 = (k1 < jb) ? k1 : jb;
  x=0.0;
  for(k=k0; k<=k1; ++k) x=-x-q[lc][k]*q[lb][ma-k]*q[la][jb-k];

  return (double(1-2*((k1+lb+jc+1)%2))*x*sqrt(fn));
  }

double ang_mom_functions::s6j(int ja, int jb, int jc,
                              int la, int lb, int lc) {
//
//    calculates 6j-symbols           
//
  i=kh[ja+jb-jc]+kh[jb+jc-ja]+kh[jc+ja-jb]+kh[ja+lb-lc]
       +kh[lb+lc-ja]+kh[lc+ja-lb]+kh[la+jb-lc]+kh[jb+lc-la]
       +kh[lc+la-jb]+kh[la+lb-jc]+kh[lb+jc-la]+kh[jc+la-lb];
  if (i) return 0.0;

  mt=(ja+jb+jc)/2 + 1;
  ma=(ja+lb+lc)/2 + 1;
  mb=(la+jb+lc)/2 + 1;
  mc=(la+lb+jc)/2 + 1;
  na=mt-ja;
  nb=mt-jb;
  nc=mt-jc;
  ka=ma-lc;
  kb=mb-lc;
  kc=mc-jc;
  fss=q[mt][ja+1]*q[ja][nc-1]/(q[ma][ja+1]*q[ja][ka-1]*q[mb][la+1]*
       q[la][kb-1]*q[mc][la+1]*q[la][kc-1]);
  fs=sqrt(fss)/(la + 1.0);
  l0 = (mt > ma) ? mt : ma;
  l0 = (l0 > mb) ? l0 : mb;
  l0 = (l0 > mc) ? l0 : mc;
  l1 = (ma+na < mb+nb) ? ma+na : mb+nb;
  l1 = (  l1  < mc+nc) ?  l1   : mc+nc;
  x=0.0;
  for(l=l0; l<l1; ++l) x=-x+q[l][mt]*q[na-1][l-ma]*q[nb-1][l-mb]*q[nc-1][l-mc];

  return (double(1-2*(l1%2))*fs*x);
  }


double ang_mom_functions::s9j(int ja, int jb, int je,
                              int jc, int jd, int jf,
                              int jg, int jh, int jt) {
//
//     Calculates 9j-symbols
//
//           | A B E \
//     snj = { C D F }
//           \ G H T |
//
//      with ia=2*A,  ib=2*B,  ic=2*C,  ...

  i =  kh[ja+jb-je]+kh[jb+je-ja]+kh[je+ja-jb]
      +kh[jc+jd-jf]+kh[jd+jf-jc]+kh[jf+jc-jd]
      +kh[jg+jh-jt]+kh[jh+jt-jg]+kh[jt+jg-jh]
      +kh[ja+jc-jg]+kh[jc+jg-ja]+kh[jg+ja-jc]
      +kh[jb+jd-jh]+kh[jd+jh-jb]+kh[jh+jb-jd]
      +kh[je+jf-jt]+kh[jf+jt-je]+kh[jt+je-jf];
  if (i) return 0.0;

  la=(je+jf+jt)/2+1;
  ld=(jg+jh+jt)/2+1;
  ma=(ja+jc+jg)/2+1;
  mc=(jf+jc+jd)/2+1;
  na=(jb+jd+jh)/2+1;
  nb=(jb+je+ja)/2+1;
  le=(je+jf-jt)/2;
  lf=(jf+jt-je)/2;
  lg=(jt+je-jf)/2;
  me=(ja+jc-jg)/2;
  mf=(jc+jg-ja)/2;
  mg=(jg+ja-jc)/2;
  ne=(jb+jd-jh)/2;
  nf=(jd+jh-jb)/2;
  ng=(jh+jb-jd)/2;
  lx=(jt+jg-jh)/2;
  mx=(jc+jd-jf)/2;
  nx=(jb+je-ja)/2;
  fn=q[la][jt+1]*q[jt][lg]*q[ma][jc+1]*q[jc][mf]*q[na][jb+1]*q[jb][ne];
  fd=q[ld][jt+1]*q[jt][lx]*q[mc][jc+1]*q[jc][mx]*q[nb][jb+1]*q[jb][nx];
  
  jsi = (abs(je-jh) > abs(jg-jf)) ? abs(je-jh) : abs(jg-jf);
  jsi = (  jsi      > abs(ja-jd)) ?   jsi      : abs(ja-jd);
  jsf = (je+jh < jg+jf) ? je+jh : jg+jf;
  jsf = ( jsf  < ja+jd) ?  jsf  : ja+jd;
  fs=double(1-2*(jsi%2))*sqrt(fn/fd)/double((jg+1)*(je+1));
  u=0.0;
  for(js=jsi; js<(jsf+2); js+=2) {
    lb=(je+jh+js)/2+1;
    lc=(jg+jf+js)/2+1;
    mb=(ja+jd+js)/2+1;
    ly=(je+jh-js)/2;
    my=(jg+jf-js)/2;
    ny=(ja-jd+js)/2;
    ud=q[lb][je+1]*q[je][ly]*q[lc][jg+1]*q[jg][my]*q[mb][js+1]*q[js][ny];
    //
    l0 = (la > lb) ? la : lb;
    l0 = (l0 > lc) ? l0 : lc;
    l0 = (l0 > ld) ? l0 : ld;
    l1 = (le+ld < lf+lb) ? le+ld : lf+lb;
    l1 = (  l1  < lg+lc) ?   l1  : lg+lc;
    //
    m0 = (ma > mb) ? ma : mb;
    m0 = (m0 > mc) ? m0 : mc;
    m0 = (m0 > lc) ? m0 : lc;
    m1 = (me+lc < mf+mb) ? me+lc : mf+mb;
    m1 = (  m1  < mg+mc) ?   m1  : mg+mc;
    //
    n0 = (na > nb) ? na : nb;
    n0 = (n0 > mb) ? n0 : mb;
    n0 = (n0 > lb) ? n0 : lb;
    n1 = (ne+lb < nf+nb) ? ne+lb : nf+nb;
    n1 = (  n1  < ng+mb) ?   n1  : ng+mb;
    //
    x=0.0;
    for(l=l0; l<=l1; ++l) x=-x-q[l][la]*q[le][l-ld]*q[lf][l-lb]*q[lg][l-lc];
    y=0.0;
    for(m=m0; m<=m1; ++m) y=-y-q[m][ma]*q[me][m-lc]*q[mf][m-mb]*q[mg][m-mc];
    z=0.0;
    for(n=n0; n<=n1; ++n) z=-z-q[n][na]*q[ne][n-lb]*q[nf][n-nb]*q[ng][n-mb];
    //
    u += double(1-2*((l1+m1+n1)%2))*x*y*z/ud;
  }

  return (u*fs);
  }

