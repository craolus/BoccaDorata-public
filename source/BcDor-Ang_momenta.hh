//
//  routines for the 3-j, 6-j and 9-j coefficients
//
//   These have been taken from the CENS package (of Hjorth-Jensen):
//
//            !    Program block heff-modules.f90
//            !
//            !    Author:   Morten Hjorth-Jensen
//            !    ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
//            !    E-MAIL:   morten.hjorth-jensen@fys.uio.no
//            !    LANGUAGE: F90/F95
//            !
//
//  and translated into c++ by Carlo Barbieri (C.Barbieri@gsi.de)
//
//  July, ??  2006.
//


#ifndef _ANG_MOM_H
#define _ANG_MOM_H


//
// size of the table of factorial for 3n-j symbols
#define    SIZEJ  60    // note it must be:
#define    SIZEK  240   //   4*SIZEJ = SIZEK > 4*jmax is a good choice


     
class ang_mom_functions {
 private:
  //double f_mb[SIZEJ],g_mb[SIZEJ],w_mb[SIZEJ];  //for Moshinsky brakets
  int kha[SIZEK+1];
  int* kh;
  double q[SIZEJ][SIZEJ];  // binomial funcion
  //double cn(0:51,0:51) // for the spherical harmonincs

 public:
   ang_mom_functions(void );
//  ~ang_mom_functions(void );
  int inline TriIneq(int twoj1,int twoj2,int twoj3) { //triangular inequality
                          return ( kh[twoj1+twoj2-twoj3]
                                  +kh[twoj2+twoj3-twoj1]
                                  +kh[twoj3+twoj1-twoj2]);}
  double ClG(int, int, int, int, int, int );
  double s3j(int, int, int, int, int, int );
  double s6j(int, int, int, int, int, int );
  double s9j(int, int, int, int, int, int, int, int, int );

  };

extern ang_mom_functions am;

#endif
