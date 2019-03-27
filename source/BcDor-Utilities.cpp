//
//
//   Utility stuff
// =================
//
//



#include <cstdlib>
//#include <cstdio>
#include <cmath>

#include <iostream>
using namespace std;

#include "BcDor-Sp_classes.hh"
#include "BcDor-Utilities.hh"

//
// sorting routines:
//
static double *dblex;

void Sort_u_double2dim(int ntot, int NDIM, int isort,double dbl[]) {
  int    i,ir,j,l;
  dblex = new double[NDIM];
  int n;

  if (ntot < 2) return;
  l = ntot/2;
  ir = ntot-1;
  while(1) {
    if (l > 0) {
       l--;
       for(n=0; n<NDIM; ++n) dblex[n] = dbl[l*NDIM+n];
    } else {
       for(n=0; n<NDIM; ++n) dblex[n]         = dbl[ir*NDIM+n];
       for(n=0; n<NDIM; ++n) dbl[ir*NDIM+n]   =         dbl[n];
       if ((--ir) == 0) {for(n=0; n<NDIM; ++n) dbl[n] = dblex[n];
                         break; }
    }
    i = l;
    j = l+l+1;
    while (j <= ir) {
       if ((j < ir) && (dbl[j*NDIM+isort] < dbl[(j+1)*NDIM+isort])) j++;
       if (dblex[isort]  < dbl[j*NDIM+isort]) {
         for(n=0; n<NDIM; ++n) dbl[i*NDIM+n] = dbl[j*NDIM+n];
         i = j;
         j = j+j+1;
         } else break;
       }
    for(n=0; n<NDIM; ++n) dbl[i*NDIM+n] = dblex[n];
    }

  delete [] dblex;  dblex = NULL;
  return;
  }


void Sort_d_double2dim(int ntot, int NDIM, int isort,double dbl[]) {
  int    i,ir,j,l;
  dblex = new double[NDIM];
  int n;

  if (ntot < 2) return;
  l = ntot/2;
  ir = ntot-1;
  while(1) {
    if (l > 0) {
       l--;
       for(n=0; n<NDIM; ++n) dblex[n] = dbl[l*NDIM+n];
    } else {
       for(n=0; n<NDIM; ++n) dblex[n]         = dbl[ir*NDIM+n];
       for(n=0; n<NDIM; ++n) dbl[ir*NDIM+n]   =         dbl[n];
       if ((--ir) == 0) {for(n=0; n<NDIM; ++n) dbl[n] = dblex[n];
                         break; }
    }
    i = l;
    j = l+l+1;
    while (j <= ir) {
       if ((j < ir) && (dbl[j*NDIM+isort] > dbl[(j+1)*NDIM+isort])) j++;
       if (dblex[isort]  > dbl[j*NDIM+isort]) {
         for(n=0; n<NDIM; ++n) dbl[i*NDIM+n] = dbl[j*NDIM+n];
         i = j;
         j = j+j+1;
         } else break;
       }
    for(n=0; n<NDIM; ++n) dbl[i*NDIM+n] = dblex[n];
    }

  delete [] dblex;  dblex = NULL;
  return;
  }

//=======================================================================================================


//
// Tinetic energy mtx. el. for h.o. basis
//

static int    ndim, n_orb;
static double x1;

int Get_Tab_kin_ho(double Uab[], int NDIM, int ish, ModSpace_t *MdSp) {

  //ndim = (NDIM  < MdSp->MSp_no[ish]) ? NDIM : MdSp->MSp_no[ish];
  ndim = NDIM;
  n_orb = MdSp->MSp_no[ish];

  if ((NULL == Uab ) || (NDIM < 1) ) return -100; // no 1-body potential is provided

  int na,nb,l;
  for(int nc=0; nc<ndim; nc++)
   for(int nr=0; nr<ndim; nr++) {
     na = nr; if (nr < n_orb ) na = MdSp->MSp_n[ish][nr];
     nb = nc; if (nc < n_orb ) nb = MdSp->MSp_n[ish][nc];
     l  = MdSp->MSp_l[ish];
     Uab[nr*NDIM+nc] = (MdSp->htom)*fab_tkin(na,nb,l);
     }

  return ndim;
  }


double fab_rsq(int na, int nb, int l) {
  // NOTE: the quantum numbers na,nb >= 0 (i.e. they start from zero)
  //       fab_rsq is given in units of bho^2
  x1 = 0.0;
  if (na == nb)        x1 = double(na+nb+l) + 1.5;
  if (abs(na-nb) == 1) x1 = -sqrt(double(na+nb+1)*double(na+nb+l+l+2)/4.0);
  return x1;
  }

double fab_tkin(int na, int nb, int l) {
  // NOTE: the quantum numbers na,nb >= 0 (i.e. they start from zero)
  //       fab_tkin is given in units of hbar*omega
  x1 = 0.0;
  if (na == nb)        x1 = (double(na+nb+l) + 1.5) / 2.0;
  if (abs(na-nb) == 1) x1 = sqrt(double(na+nb+1)*double(na+nb+l+l+2)/4.0)/2.0;
  return x1;
  }

