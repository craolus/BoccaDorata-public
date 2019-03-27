//#include <stdio.h>

#define BUFFER_SIZE 300

#include <cstdlib>
#include <cstring>
//#include <cstdio>
#include <cmath>

#include <iostream>
#include <fstream>
using namespace std;

#include "BcDor-Sp_classes.hh"
#include "BcDor-Int_classes.hh"
//#include "BcDor-Ang_momenta.hh"


static char line[BUFFER_SIZE];
static char *lptr;
static int icount, count_outside, count_err;//, count_bad, count_discarded, count_added;

static int i,n1,n2;
static int ia,ib,ic,id;
static double xe;

int VppInt_t::read_Oslofmt(char *fname, int n_all_new_ints, int n_me_min) {
                                  // defualts: n_all_new_ints=1
                                  //           n_me_min=0

  double BetaCM=0.0;
  int norb_in = 0;
  int n_sten = -1;
  int nints_rd = -1, iVeff= -10;
  int iCOM=0;
  double* Est=NULL;
  char Int_name[20];

  int nme1, nme2, nme3, nme4;


  int BuildMdSp = 0;
  if (NULL == MdSp) BuildMdSp = 1;

  cout << " Get Vpp mtx elements from file '" << fname << "'...   \n";

  bool test_sep_files = false;
  char fname1[BUFFER_SIZE], fname2[BUFFER_SIZE];
  sprintf(fname2,"%s", fname);
  FILE* infile;
  infile = fopen(fname, "rb");
  if (NULL == infile) {
   test_sep_files= true;
   sscanf(fname,"%s.mhj", fname1);
   sprintf(fname2,"%s_msp.mhj", fname1);
   infile = fopen(fname2, "rb");
  }

  Extnl_msp OsloMSP;
  OsloMSP.load_oslo_msp(infile);

  if (!BuildMdSp) OsloMSP.Sync(MdSp);
  

  if (BuildMdSp) {
    cerr << "\n\n  No model space was specified. In this case one could build it using the"
         <<   "\n single particle basis given in `" << fname2 << "'...   "
         << "\n\n  However, this option has not been coded yet. So I will quit! :-P";
     exit(8);
     // build basis...
    //ModSpace_t  NewMdSp;
    //NewMdSp.Build_MdSp(OsloMSP);
    //...
  }


  // get # of sp states that are also in the basis:
  norb_in = 0;
  for(i=1; i<=OsloMSP.Ntot_orbs; ++i) if (OsloMSP.i_msp[i] >= 0) ++norb_in;

  if (OsloMSP.Ntot_orbs > norb_in) cerr << "\n WARNING (VppInt_t::read_Oslofmt): the input file "
                           << "contains some sp orbits that are not in the model space."
                           << " --> the mtx els outside the mod. sp. will be idscarded!\n";
  
  // check bHO
  if ((fabs(OsloMSP.bHO-(MdSp->bho)) > 1.e-5)&&(fabs(MdSp->bho)>0.001)) {cerr << "\n BAD bHO!!!!!" << OsloMSP.bHO << "  " << MdSp->bho; exit(9);}



  //  If separate files are being used, closes the model space file
  // and opens the one for the interaction:
  if (test_sep_files) {
   fclose(infile);
   sscanf(fname,"%s.mhj", fname1);
   sprintf(fname2,"%s_int.mhj", fname1);
   infile = fopen(fname2, "rb");
  }


  //
  // Starts reading the interaciotn:
  //

// "   ----> Interaction part                                                 "1
// "Nucleon-Nucleon interaction model:n3lo                                    "2
// "Type of calculation:no-core                                               "3
// "Number and value of starting energies:****                                "4
// "0.498492+160                                                              "5
// "Total number of twobody matx elements:   24198    4060    16078     4060  "6
// "Matrix elements with the following legend, NOTE that the interaction ...  "7
// "  Tz Par  2J   a   b   c   d          <ab|V|cd>                           "8

  fgets(line,sizeof(line), infile); //1
  fgets(line,sizeof(line), infile); sscanf(line,"Nucleon-Nucleon interaction model:%s",&Int_name);
  fgets(line,sizeof(line), infile); //3
  iVeff = -10;
  if (strstr(line,"no-core")!=NULL) {
    iVeff = 0;
    cout << "\nLee-Suzuki case not coded!\n"; return(-1);
  } else if (strstr(line,"g-matrix")!=NULL) {
    iVeff = 1;
  } else if (strstr(line,"vlowk")!=NULL) {
    iVeff = 2;
  } else if (strstr(line,"vsrg")!=NULL) {
    iVeff = 2;
  } else if (strstr(line,"bare")!=NULL) {
    iVeff = 2;
  } else {
    iVeff = -10;
    cout << "\nline:" << line << endl;
    cout << "cannot recognize the type of effective interaction!\n"; return(-1);
  }
  //
  //
  iCOM = 0;
  if (strstr(line,"com-t1")!=NULL) {
    iCOM = 1; /* kin en. in t1 form */
    cout << "\n\n WILL correct the 2b-me by  -(p_i*p_j)/mA  (with A="<<OsloMSP.Atot<<")...\n\n";
  }
  else if (strstr(line,"com-t2")!=NULL) { iCOM = 2; /* kin en. in t2 form */}
  else if (strstr(line,"AddHcm")!=NULL) { iCOM = 3; /* Add H_com */
    if (1 > sscanf(line,"AddHcm=%lf",&BetaCM)) BetaCM=10.0;}
  else { iCOM = 0; /* do nothing */}
  
  
  /*TST*/ cout << "\n(iVeff = "<<iVeff<<" ,       iCOM = "<<iCOM<<" (BetaCM="<<BetaCM<<") , Interaction name: -->"<<Int_name<<"<--)\n\n"<<flush;
  

  nints_rd = 1;
  if (1 == iVeff) {
    // g-matrix
    fgets(line,sizeof(line), infile); sscanf(line+40,"%i",&n_sten);
    fgets(line,sizeof(line), infile); //5
      Est = new double[n_sten];
      lptr = line;
      for(i=0; i<n_sten; ++i) {
        sscanf(lptr,"%lf%n",Est+i,&n1);
        lptr += n1;}
    nints_rd = n_sten;
    cout << "Starting energies: " << line;
  } else {
    fgets(line,sizeof(line), infile); //4
    fgets(line,sizeof(line), infile); //5
  }
  
  fgets(line,sizeof(line), infile); sscanf(line+40,"%i%i%i%i",&nme1,&nme2,&nme3,&nme4);
  fgets(line,sizeof(line), infile); //7
  fgets(line,sizeof(line), infile); //8




  if (n_all_new_ints >=0) {n_all_new_ints += nints_rd;} // add extra space
                     else {n_all_new_ints = nints_rd;}  // just the # of st. ens.

  if (n_me_min < nme1) n_me_min = nme1;

  cout << "Declaring a max # of " << n_all_new_ints << " interations and allocating "
       << nints_rd << " of them, the tot. # of me to be read is " << n_me_min << "...       ";
  //  set for up to `n_all_new_ints' forces
  // and allocate `nints_rd' of them
  Allocate_Vpp_table(n_all_new_ints, n_me_min, nints_rd);
  select(0);
  cout << " ok!\n";



  int iTz,iPar,J2,na,nb,nc,nd;//,ifct;
  double *xVs = new double[nints_rd];

  // set the st. ens:
  if (1 == iVeff) {
    for(i=0; i<n_sten; ++i) if (i < n_ints) StrEn[i] = Est[i];
    //cout << "Starting energies:\n";
    //for(i=0; i<n_sten; ++i) cout << "i_sel = " << i << " ,   Est = " << StrEn[i] << endl;
 }


  cout << " " << nme1 << " lines to be read, with " << nints_rd << " me(s) per line.\n";

  icount=0; count_outside=0; count_err=0;
  for(i=0; i<nme1; ++i) {
  
    // read the conf.
    fgets(line,sizeof(line), infile);
    sscanf(line,"%i%i%i %i%i%i%i%n",&iTz,&iPar,&J2,&na,&nb,&nc,&nd,&n1);
    lptr = line + n1;
    //  check basis
    ia = OsloMSP.i_msp[na];
    ib = OsloMSP.i_msp[nb];
    ic = OsloMSP.i_msp[nc];
    id = OsloMSP.i_msp[nd];
    if ((ia <0) || (ib <0) || (ic <0) || (id <0)) {++count_outside; continue;}

    // read the mtx els.
    for(n2=0; n2<nints_rd; ++n2) {
        sscanf(lptr,"%lf%n",xVs+n2,&n1);
        lptr += n1;}


    if (1 == iCOM) {
      sscanf(lptr,"%lf%n",&xe,&n1); lptr += n1;  // read p.p/m interaction
      for(n2=0; n2<nints_rd; ++n2) xVs[n2] -= xe/double(OsloMSP.Atot);
    }
    else if (2 == iCOM) {cout << "\n\n com-t2 case is not coded yet (will stop)!\n"; return(-1);}
    else if (3 == iCOM) {
      sscanf(lptr,"%lf%n",&xe,&n1); lptr += n1;  // read p.p/m interaction
      sscanf(lptr,"%lf%n",&xe,&n1); lptr += n1;  // read Hcom interaction
      cout << "\n\n AddHcm case is not coded yet (will stop)!\n"; return(-1);
    }


    //if (0==iTz) {
    //  ifct = 4;
    //  if ((na/2)+(na%2) != (nb/2)+(nb%2)) ifct /=2;
    //  if ((nc/2)+(nc%2) != (nd/2)+(nd%2)) ifct /=2;
    //  if ( 1 != ifct ) {
    //    for(n2=0; n2<nints_rd; ++n2) xVs[n2] /= sqrt(double(ifct));
    //  }
    //}

    //store them:
    if (add(ia,ib,ic,id,J2/2,xVs,nints_rd,0)) {++count_err;} else {++icount;}

  }
  if (icount != n_me_stored) 
                  {cerr << " \n\n SOMETHING wrong with the # of stored me !!!!!!"; exit(8);}
  if ((icount+count_err+count_outside) != nme1) 
                  {cerr << " \n\n SOMETHING wrong with the # of read me !!!!!!"; exit(8);}

  fclose(infile);


  prepare_for_reading();

  if (count_outside) cout << count_outside << " lines were out of the model space\n";
  if (count_err) cout << count_err << " lines had incompatible quantum numbers\n";
  cout << "  --> " << n_me_stored << " Vpp matrix elements have been loaded successfully!\n\n";


  if (NULL != Est) delete [] Est;

  delete []  xVs; xVs = NULL;

  return nints_rd;
  }

