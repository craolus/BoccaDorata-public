//
//
//  
//
//  C.Barbieri, GSI, December 2008.  (C.Barbieri@gsi.de)
//

#include <iostream>

#ifndef _UTIL_LISTS_H
#define _UTIL_LISTS_H


struct VectList {

 public:
  int NDIM, NVECT, n_vcts; // n_dim;
  int itype;  // to pass, e.g., the type of Krylov algorithm
  int *ia, *ib;
  double *data;
  double **vect;

  //
  // constructors/destructors:
  inline  VectList(void ) {init_NULL();}
  inline  VectList(int nd, int nv, double val=0.0, int nval=0, int nvaltype=0)
                                  {init_NULL(); set(nd, nv,val,nval,nvaltype);}
  inline ~VectList(void ) {clear();}


  inline void init_NULL(void ) {NDIM=-1; NVECT=-1; n_vcts=0; itype=0;
                                data=NULL; vect=NULL; ia=NULL; ib=NULL;}
  inline void clear(void ) {if (NULL != vect) delete [] vect;
                            if (NULL != data) delete [] data;
                            if (NULL != ia) delete [] ia;
                            if (NULL != ib) delete [] ib;
                            init_NULL();}
  inline void set(int nd, int nv, double val=0.0, int nval=0, int nvaltype=0) {
                            clear();
                            if ((nv<1) || (nd<1)) return;
                            NVECT = nv;
                            NDIM  = nd;
                            data = new double  [NVECT*NDIM];
                            vect = new double* [NVECT];
                            ia   = new int     [NVECT];
                            ib   = new int     [NVECT];
                            for(nd=0; nd<NVECT*NDIM; ++nd) data[nd] = val;
                            for(nv=0; nv<NVECT; ++nv) {
                               vect[nv] = data + nv*NDIM;
                               ia[nv] =nval;
                               ib[nv] =nval;}
                            itype = nvaltype;
                            }

};


class Execution_clock {

 int NTOT_CLOCK;
 clock_t*  ora;
 int n_ora;

public:
  Execution_clock(void ) {ora=NULL; n_ora=0; init(20);}
  Execution_clock(int i) {ora=NULL; n_ora=0; init(i); }
 ~Execution_clock(void ) {clear(); }

 void init(int i) {clear(); NTOT_CLOCK=i;  ora = new clock_t[NTOT_CLOCK]; resettime(); return;}
 void clear(void) {if (NULL != ora) delete [] ora; ora = NULL; NTOT_CLOCK=-100;}

 double marktime(void ) {
   if (n_ora < NTOT_CLOCK)  {ora[n_ora] = clock(); ++n_ora;}
   return double(ora[n_ora-1]-ora[n_ora-2])/CLOCKS_PER_SEC;}

  void resettime() {
    ora[0] = clock();
    n_ora=1;
    return;}
  void plottime(void ) {
    for(int i=1; i<n_ora; ++i)
      cout << "time:  " << double(ora[i] - ora[i-1])/CLOCKS_PER_SEC << endl;
      cout << "Total: " << double(ora[n_ora-1] - ora[0])/CLOCKS_PER_SEC << endl;
    cout << flush;
    //n_ora = 0;
    return;}

};

extern Execution_clock exclock;

#endif

