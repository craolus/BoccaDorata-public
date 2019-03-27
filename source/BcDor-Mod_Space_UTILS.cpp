

#define BUFFER_SIZE 300


#include <cstdlib>
//#include <cstdio>
#include <cstring>
//#include <cmath>

#include <iostream>
using namespace std;

#include "BcDor-Sp_classes.hh"



Extnl_msp::Extnl_msp(char* fname) {
  init_NULL();
  bHO = 0.0;
  Atot = 0;
  MdSp = NULL;
  load_oslo_msp(fname);
  return;
  }

Extnl_msp::Extnl_msp() {
  init_NULL();
  bHO = 0.0;
  Atot = 0;
  MdSp = NULL;
  return;
  }

Extnl_msp::~Extnl_msp() {
  free_mem();
  return;
  }


void Extnl_msp::free_mem(void ) {
  if (NULL != i_ext) delete [] i_ext; i_ext = NULL;
  if (NULL !=  Esp ) delete [] Esp;    Esp  = NULL;
  if (NULL !=  ch  ) delete [] ch;     ch   = NULL;
  if (NULL !=  ip  ) delete [] ip;     ip   = NULL;
  if (NULL !=  twj ) delete [] twj;    twj  = NULL;
  if (NULL !=   n  ) delete [] n;       n   = NULL;
  if (NULL != i_msp) delete [] i_msp; i_msp = NULL;
  init_NULL();
  return;
  }


void Extnl_msp::init_NULL(void ) {
  i_ext = NULL;
  i_msp = NULL;
      n = NULL;
    twj = NULL;
     ip = NULL;
     ch = NULL;
    Esp = NULL;
  Ntot_orbs = -1;
  return;
  }


void Extnl_msp::alloc_mem(int norb) {
  free_mem();
  if (norb < 1) return;
  Ntot_orbs = norb;
  // // i_ext = new int[MdSp->nsubsh];
  i_msp  = new int[Ntot_orbs+1];
  n   = new int[Ntot_orbs+1];
  twj  = new int[Ntot_orbs+1];
  ip  = new int[Ntot_orbs+1];
  ch  = new int[Ntot_orbs+1];
  Esp = new double[Ntot_orbs+1];
  return;
  }

int Extnl_msp::load_oslo_msp(char* fname) {
  FILE* OlsoSpFile;
  OlsoSpFile = fopen(fname, "rb");
  if (NULL == OlsoSpFile) {
    cerr << "ERROR (Extnl_msp::load_oslo_msp): Cannot open file '"<<fname<<"'...";
    return 100;
    }
  int i_ret = load_oslo_msp(OlsoSpFile);
  fclose(OlsoSpFile);
  return i_ret;
  }


static int n1,n2,n3;
static char line[BUFFER_SIZE];
static char *lptr;
static int nbs,j2bs,ipbs,chbs;
static double xe;

int Extnl_msp::load_oslo_msp(FILE* infile) {
//  The Header of the file must be as follows:
//
//  0         1         2         3         4         5         6         7
//  01234567890123456789012345678901234567890123456789012345678901234567890
// "   ----> Oscillator parameters, Model space and single-particle data      "1
// "Mass number A of chosen nucleus (important for CoM corrections):        2 "2
// "Oscillator length and energy: 0.172110E+01  0.140000E+02                  "3
// " Min and max value of partial wave ang. mom           0           6       "4
// " Max value of relative orb mom or cm orb mom,  l or L=            7       "5
// " Max value of relative n:          50                                     "6
// " Max value of 2*n + l+ cm 2*N +L for large space:         100             "7
// " Max value of 2*n + l+ cm 2*N +L for model space:           6             "8
// " Total number of single-particle orbits          56                       "9
// "Legend:         n   l  2j  tz  energy  particle or hole                   "0
// "Number:   1     0     0     1    -1    0.210000E+02  hole                 "1
// "    :    "
// "    :    "
// "Number:  56     3     0     1     1    0.105000E+03  particle             "


  fgets(line,sizeof(line), infile); //1
  fgets(line,sizeof(line), infile); sscanf(line+65,"%i",&Atot);
  fgets(line,sizeof(line), infile); sscanf(line+29,"%lf",&bHO);
  fgets(line,sizeof(line), infile); //4
  fgets(line,sizeof(line), infile); //5
  fgets(line,sizeof(line), infile); //6
  fgets(line,sizeof(line), infile); //7
  fgets(line,sizeof(line), infile); //8
  fgets(line,sizeof(line), infile); sscanf(line+44,"%i",&n1);
  fgets(line,sizeof(line), infile); //10

  alloc_mem(n1);

  for(n2=1; n2<=Ntot_orbs; ++n2) {
    fgets(line,sizeof(line), infile);
      lptr = line;
      if (strstr(line,"Number:") != NULL) lptr += 7;
      sscanf(lptr,"%i%i%i%i%i%i%lf",&n1,&nbs,&ipbs,&j2bs,&chbs,&n3,&xe);
    ipbs %= 2;          // l  -->  ip
    chbs = (1-chbs)/2;  // tz (=-1 for protons)  -->  ch (=0 neut. , =1 prot.)
    i_msp[n2] = -1;
    n[n2]  = nbs;
    twj[n2]= j2bs;
    ip[n2] = ipbs;
    ch[n2] = chbs;
    Esp[n2]= xe;
  }

  return 0;
  }

int Extnl_msp::Sync(ModSpace_t *msp_in) {
  MdSp  = msp_in;
  i_msp[0] = -100;
  for(n2=1; n2<=Ntot_orbs; ++n2) {
    i_msp[n2] = MdSp->get_orbit(n[n2],twj[n2],ip[n2]%2,ch[n2]);
    //cout << n2 << " ---  "<<flush;
    //cout << n[n2] << "   "  << twj[n2  << "   "  << ip[n2]%2   << "   "  <<ch[n2]  << "  ---   " << flush; 
    //cout << n2 << "  " << i_msp[n2] << "    "  << MdSp->Mor_n[i_msp[n2]] 
    //                  << MdSp->MSp_name[MdSp->Mor_sh[i_msp[n2]]]
    //                  << "   Esp=" << Esp[n2] << endl << flush;
  }
  if (NULL != i_ext) delete [] i_ext; i_ext = NULL;
  i_ext = new int[MdSp->nsubsh];
  for(n2=0; n2<MdSp->nsubsh; ++n2) {
    i_ext[n2] = -100;
    for(n1=1; n1<=Ntot_orbs; ++n1) if (n2 == i_msp[n1]) {i_ext[n2] = n1; break;}
  }
  return 0;
  }

//int ModSpace_t::Build_MdSp(Extnl_msp Omsp_in) {
//  return 0;}

