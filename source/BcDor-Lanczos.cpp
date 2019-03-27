//
//
//   Solve the Dyson Equation by first performing Lanczos iteration on the
//  fw- and bk- parts of the self-energy. Multiple pivots are used.
//
//

//#include <iostream>
using namespace std;

#include <cstdlib>
//#include <cstdio>
#include <cmath>

#include "BcDor-Global_variables.hh"
#include "BcDor-Run_vars.hh"
//#include "BcDor-Utilities.hh"



static int i,j,l;

static int i_pvt, n_pvt, j_next_pvt;

static int NTOT_in, NTOT_out;

static double *LancVect = NULL;    // Store Lanczos vectors
                                   // set NULL @ beginning for safety

static int *i_last_pvt  = NULL;//, *i_used_pvt;

static double *ptr_i, *ptr_ms, *ptr_j, *ptr_mlnc_e, *ptr_mlnc_m;
static double *ptr1, *ptr2; // , *ptr3;
static double *ptr_up, *ptr_lo;
static double xa, xb, xg; // xe
static double x1, x2; //, x3, x4;

static int nLanc;

static bool new_vector;



int Reduce_Mtx_Lncz(int NA_out, int NB_out, int NPLS_out, double *Mout, int* i_piv_out,
                    int NA_in,  int NB_in,  int NPLS_in,  double *Min,
                    int nt_pvts, double** pivots, int* n_itr) {
  //
  //  If everything goes fine, it will return the number of pivots actually processed.
  //  If there are errors, it will return an integer < 0.
  //
		
  if ((NPLS_in < 1)  || (NPLS_out < 1)) return -1;
  
  if (NB_out != NB_in) {
    cout << "We have no time, and quite some problem here!\n";
    exit(-100);
  }
  
  nLanc = 0;  j=0;
  for(i=0; i<nt_pvts; ++i) if (n_itr[i] > 0) {nLanc += n_itr[i]; ++j;}
  if (nLanc < 1) return -2;
  if (nLanc > NPLS_in) nLanc = NPLS_in;  // the number of its. already saturated Lanczos
  if (nLanc > NPLS_out) {
    //  The calling routine did not provide enough space in Min
    // to do all iterations...
    cout << "Who is this broken man?!\n";
    exit(-100);
  }
  //
  if ((j+1 > NA_out) && (nLanc < NPLS_in)) {
    cout << j<< "   "<< NA_out<< "   "<<nLanc << "   "<<NPLS_in << "  <<<-\n";
    cout << "Yes, I'm a cpative fan. I'm dying to be shown that you are just not any man!\n";
    exit(-100);
  }
  
  
  NTOT_in  = NA_in  + NB_in ;
  NTOT_out = NA_out + NB_out;
  
  
  if (NULL != i_last_pvt) delete [] i_last_pvt;
  i_last_pvt = new int[nt_pvts]; //, i_used_pvt[nt_pvts];
  
  if (NULL != LancVect  ) delete [] LancVect; // redundant, but why not...
  LancVect = new double[NPLS_in * nLanc];
  
  // clean up the memory to store the Lanczos self-energy:
  for(ptr_i= Mout; ptr_i<Mout+(NTOT_out * nLanc);  ++ptr_i) (*ptr_i)=0.0;
  //for(ptr_i= Mout; ptr_i<Mout+NTOT_out*NPLS_out; ++ptr_i) (*ptr_i)=0.0;
  //  NOTE, here the dimension of Mout is NTOT_out*NPLS_out but
  // only the first nLanc rows will be calculated and we have
  // required nLanc <= NPLS_out  above, for safety.
  
  
  ptr_mlnc_e = Mout;
  ptr_mlnc_m = Mout + NTOT_out - NB_out;
  for(i=0; i<nt_pvts; ++i) {i_last_pvt[i]=-1000;}// i_used_pvt[i]=-1000;}
  i_pvt = -1;
  n_pvt = -1;
  j_next_pvt = 0;
  for(j=0; j<nLanc; ++j) {
    
    ptr_up=LancVect+(j+1)*NPLS_in;
    ptr_lo=LancVect+    j*NPLS_in;
    
    //
    //make new vector
    //
    new_vector = (j == j_next_pvt);
    if (new_vector) {
      // build from a new pivot
      ++i_pvt; while(n_itr[i_pvt] < 1) ++i_pvt; //seek for the first allowed pivot
      if (i_pvt >= nt_pvts) {// this should never happen though...
        cout << "\n\n ERROR (Reduce_Mtx_Lncz): there is some trouble with the number of  iterations"
        <<      " and tne number of pivots!\n"
        << "\n i_pvt, nt_pvts  =" << i_pvt << "   " << nt_pvts
        << "\n j,     nLanc =" << j     << "   " << nLanc
        << "\n ---> EXIT(-100).\n";
        return (-100);
      }
      ++n_pvt;
      j_next_pvt += n_itr[i_pvt];
      i_last_pvt[n_pvt] = j_next_pvt -1;
      //i_used_pvt[n_pvt] = i_pvt;
      //
      //
      // NOTE: there is no need to normalize the input pivots here since
      //  the resulting Lanczos vector will be normalized afterward anyway
      x1 = 0.0;
      for(i=0; i<NB_in; ++i)  x1 += (pivots[i_pvt][i])*(pivots[i_pvt][i]);
      x1 = sqrt(x1);
      for(i=0; i<NB_in; ++i)  pivots[i_pvt][i]  /= x1;
      cout << "  pivot n. " << i_pvt << ":  (" << n_itr[i_pvt] << " itrs)  ";// endl;
      for(i=0; i<NB_in; ++i)  cout << "    " << pivots[i_pvt][i];
      cout << endl;
      
      ptr_ms = Min + NA_in;
      for(ptr_j = ptr_lo; ptr_j<ptr_up; ++ptr_j) {
        x1=0.0;
        for(l=0; l<NB_in; ++l) x1 += ptr_ms[l] * pivots[i_pvt][l];
        (*ptr_j) = x1;
        ptr_ms += NTOT_in;
      }
      
    } else {
      // Lanczos iteration
      ptr_j = LancVect + j*NPLS_in;
      ptr_i = ptr_j - NPLS_in;
      ptr_ms = Min;
      for(i=0; i<NPLS_in; ++i) {ptr_j[i] = ptr_i[i]*(*ptr_ms); ptr_ms+=NTOT_in;}
      //      ptr_ms = Min;
      //      ptr_i = ptr_lo - NPLS_in;
      //      for(ptr_j = ptr_lo; ptr_j<ptr_up; ++ptr_j) {(*ptr_j) = (*ptr_i)*(*ptr_ms);
      //                                                  ptr_ms+=NTOT_in;
      //                                                  ++ptr_i;}
      
    } // end of making the new vector
    
    
    //
    //  Orthogonalize to previous vectors using
    // the Gram-Schmidt algorithm:
    ptr_j = LancVect + j*NPLS_in;
    for(i=j-1; i>=0; --i) { // <-- NOTE going backward do not sum up the
      // numerical errors in successive iterations!!!
      // So, DON'T do "for(i=0; i<j; ++i) {  }"
      ptr_i = LancVect + i*NPLS_in;
      x1 = 0.0;  x2 = 0.0;
      for(l=0; l<NPLS_in; ++l) {x1 += ptr_i[l]*ptr_j[l]; x2 += ptr_i[l]*ptr_i[l];}
      for(l=0; l<NPLS_in; ++l) ptr_j[l] -= x1*ptr_i[l]/x2;
      //  note, there would be no need to divide by x2 (in principle)
      // since all preceding vectors with i<j should already be normalized...
    }

    //
    // Normalize and compute the componets of
    // the Lanczos self-energy (a, b, g_i, m_i)
    //-----------------------------------------------
    
    // Compute lanc. a,b mtx. els. and the contraction to EHF self. en.
    xa = 0.0;
    xb = 0.0;
    ptr_ms = Min;
    for(ptr_j=ptr_lo; ptr_j< ptr_up; ++ptr_j) {
      x1 = (*ptr_j);
      xa += (*ptr_ms) * x1 * x1;
      xb +=       x1 * x1;
      ptr2 = ptr_ms + NA_in;
      for(ptr1=ptr_mlnc_m; ptr1<(ptr_mlnc_m+NB_out); ++ptr1) {
        (*ptr1) += (*ptr2) * x1;
        ++ptr2;
      }
      ptr_ms += NTOT_in;
    }
    
    //
    //normalize (and generate the g_i):
    xa /= xb;
    xb = sqrt(xb);
    ptr_mlnc_e[0] = xa;
    ptr_mlnc_e[1] = xb;
    x1 = 0.0;
    for(l=0; l<NB_out; ++l) ptr_mlnc_m[l] /= xb;
    for(i=0; i<n_pvt; ++i) {
      ptr_i = LancVect + i_last_pvt[i]*NPLS_in;
      ptr_ms = Min;
      xg = 0.0;
      for(ptr_j=ptr_lo; ptr_j< ptr_up; ++ptr_j) {
        xg += (*ptr_ms) * (*ptr_i)*(*ptr_j);
        ptr_ms += NTOT_in;
        ++ptr_i;
      }
      ptr_mlnc_e[2+i] = xg/xb;
    }
    for(ptr_j=ptr_lo; ptr_j<ptr_up; ++ptr_j) (*ptr_j)/=xb;
    if (new_vector) ptr_mlnc_e[1] = ptr_mlnc_e[1+n_pvt];
    
    //cout << "\n---> TEST beta0: " << x1 << "  " << xb << endl;
    
    if (abs(xb) < 1.e-8) {
      cout << " ItrLncz n. " << j+1 << "  ,   a=" << xa << "      b=" << xb;
      if (n_pvt > 0) {cout << "     g_i=";
        for(i=0; i<n_pvt; ++i) cout << ptr_mlnc_e[2+i] << "   ";}
      cout << endl;
    }
    
    if (n_pvt + 2 + NB_out > NTOT_out) {cout << "\n\n ERROR n_pvt=" << n_pvt
      << " has exceeded the max value allowed by NTOT_out-NB_out (?)"; exit(-200);}
    //this should never happen...


    ptr_mlnc_m += NTOT_out;
    ptr_mlnc_e += NTOT_out;
  } // end of lanczos iterations
  
  ++n_pvt;
  
  // needed later on to reconstruct the fishbone-like matrix:
  i_piv_out[0] = n_pvt;
  i_piv_out[1] = 1;   // the beta_i start here
  for(i=0; i<n_pvt; ++i) i_piv_out[2+i] = i_last_pvt[i] + 1; // the i-th g_i branch starts right
  //after the end of the i-th pivot
  
  
  if (i_last_pvt[n_pvt-1]+1 != nLanc) {
    if ((i_last_pvt[n_pvt-1]+1 > nLanc) && (NPLS_in == nLanc)) {
      cout << " max fw dimensions have been reached.\n";
    } else {cout << "\n\nERROR with n_pvt!!!!!!\n\n"; exit(-200);}
  }
  cout << endl << nLanc << " Lanczos vectors have been computed, " << n_pvt << " and pivots were used\n\n ";
  
  /*
   //
   // Compute the prodocuts of between all vectors to check
   // the orthonormalization:
   //
   for(j=0; j<nLanc; ++j) {
   ptr_j = LancVect + j*NPLS_in;
   cout << "\n" << j << ":";
   for(i=0; i<=j; ++i) {
   ptr_i = LancVect + i*NPLS_in;
   x1 = 0.0;
   for(l=0; l<NPLS_in; ++l) x1 += ptr_i[l]*ptr_j[l];
   cout << "   " << x1;//+1;
   } }
   */
  
  delete [] LancVect;  LancVect = NULL;
  delete [] i_last_pvt;  i_last_pvt = NULL;
  //delete [] i_used_pvt;  i_used_pvt = NULL;
  
  /*
   //
   // Print out the coponents of the Lanczos self-en (for inspection):
   //
   ptr_ms = Mout;
   for(i=0;  (i<30)&&(i<nLanc); ++i) {
   for(j=0; j<NTOT_out; ++j) {cout << (*ptr_ms) << "   "; ++ptr_ms;}
   cout << endl;
   }
   */
  
  //
  // MUST return the number of vectors actually calculated to
  // the calling function!!
  return nLanc;}

