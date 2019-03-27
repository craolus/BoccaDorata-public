//
//
//   Utility stuff
// =================
//
//

#ifndef _UTILITIES_H
#define _UTILITIES_H

#include "BcDor-Sp_classes.hh"

//
// sorting routines:
//

extern void Sort_u_double2dim(int , int , int , double []);

extern void Sort_d_double2dim(int , int , int , double []);

//============================================================


//
// radial integrals of spheric h.o. functions
//

extern int Get_Tab_kin_ho(double[] , int , int , ModSpace_t*);

extern double  fab_rsq(int , int , int);

extern double fab_tkin(int , int , int);


//============================================================


#endif
