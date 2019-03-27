//
//
//  Class to store the 2-body matrix elements of the (effective) nuclear
// interaction.
//  The matrix elements stored are in p-n formalism.
//
//  The function 'Vpp(ia,ib,ic,id,J)' returns:
//
//      < phi_a, phi_b ; J | V | phi_c, phi_d ; J >
//
//
//   ia, ib, ic and id are the sp orbits of the model space associated
//  with the object (== the interaction). The matrix elements are
//  antisymmetrized and properly normalized.
//
//
// C.Barbieri, RIKEN, April 2010.  (barbieri@riken.jp).
//


//#include <stdio.h>

#define CHOFF     -2
#define CHMULT     5
#define JMULT    128
#define ABMULT  1024

#define BUFFER_SIZE 300
#define CHOP   1.e-8


#include <cstdlib>
#include <cstring>
//#include <cstdio>
#include <cmath>

#include <iostream>
#include <fstream>
using namespace std;

#include "BcDor-Sp_classes.hh"
#include "BcDor-Int_classes.hh"

static void Asc_iii_iiiheapsort(           int, int[], int[], int[]);
static void Asc_iiid_iiiheapsort_singleV(  int, int[], int[], int[], double[]);
static void Asc_iiid_iiiheapsort_multipleV(int, int[], int[], int[], int, double* []);

//
// global variable for this file (made static so that they pertain to this file only)
//
  // used by write() and other routines;
static int iabJ,icdJ,J,ish,k;
static int iorb_a,iorb_b,iorb_c,iorb_d,iorb[4];
  // added for read();
static int nbs,j2bs,ipbs,chbs,n1,iph;
static double xe;
  // added for prepare_for_reading();
static int i1,i2,ish1,ish2,ch,ip;
static int itmp;
  // added for get();
static int isha,ishb,ishc,ishd;
static int iup,idown,itr;
static int iab_mid,icd_mid;
static int test;
static double xph;
  // added for add_file
static double x1;



//
// Encoding of charge, parity and tot J:
// 
inline int ind_chipJ(int ch, int ip, int J) {
                                return ( ((ch-CHOFF)*2+ip)*JMULT + J ); }

inline int J_chip_ind(int *ch, int *ip, int itmp) {
  *ch = itmp/JMULT; *ip = (*ch)%2; *ch = (*ch/2)+CHOFF;
  return (itmp%JMULT);}

inline int ind_i1i2J (int i1, int i2, int J) {return ((i1*ABMULT+i2)*JMULT+J);}

inline int J_i1i2_ind(int *i1, int *i2, int ind) {
    *i1 = ind/JMULT; *i2 = (*i1)%ABMULT;  
    *i1 /= ABMULT;   *i1 = (*i1)%ABMULT;
     return (ind%JMULT); }



// No matrix elements
VppInt_t::VppInt_t() {
  cout << "\n\nVppInt_t: you need to pass the model space before allocating Vpp matrix elements.\n\n";
  init_NULL();
  set_NULL();
  }  

// only mod. space
VppInt_t::VppInt_t(ModSpace_t *inmds) {
  init_NULL();
  set_NULL();
  MdSp = inmds;
  }

// both mod. space and int. file
VppInt_t::VppInt_t(char *fname, ModSpace_t *inmds, int n_all_new_ints, int n_me_min) {
                                  // defualts: n_all_new_ints=1
                                  //           n_me_min=0
  init_NULL();
  set_NULL();
  MdSp = inmds;
  read(fname, n_all_new_ints, n_me_min);
  }


// Destructor
VppInt_t::~VppInt_t() {
  cout << "\n\nDeallocating " << n_me_stored << " Vpp matrix elements...\n\n";
  Clear_all();
  return;}

void VppInt_t::Clear_all(void ) {
  free_mem();
  set_NULL();
  return;}

void VppInt_t::Clear_all(ModSpace_t *inmds) {
  free_mem();
  set_NULL();
  MdSp = inmds;
  return;}

void VppInt_t::set_NULL(void ) {
  MdSp      = NULL;
  n_me_stored = -1;
  n_me_alloc  = -1;
  n_ints       = -1;
  n_ints_alloc = -1;
  return;}  

void VppInt_t::init_NULL(void ) {
  Int_list  = NULL;
  StrEn     = NULL;
  Vpp_wrk   = NULL;
  iabJ_list = NULL;
  icdJ_list = NULL;
  Jipch_str = NULL;
  Jipch_end = NULL;
  return;}

void VppInt_t::free_mem(void ) {
  if (NULL != Int_list ) {
     for(int i=0; i<n_ints; ++i) {if (Vpp_wrk == Int_list[i]) Vpp_wrk = NULL;
                                  if (NULL != Int_list[i]) delete [] Int_list[i];}
     n_ints = -1;
     delete [] Int_list ;
     }
  if (NULL != StrEn    ) delete [] StrEn;
  if (NULL != Vpp_wrk  ) delete [] Vpp_wrk ;
  if (NULL != iabJ_list) delete [] iabJ_list;
  if (NULL != icdJ_list) delete [] icdJ_list;
  if (NULL != Jipch_str) delete [] Jipch_str;
  if (NULL != Jipch_end) delete [] Jipch_end;
  init_NULL();
  return;}

void  VppInt_t::Allocate_Vpp_table(int max_ints, int max_me, int n_new_ints) {
                                                   // default: n_new_ints = -1
  free_mem();
  n_ints_alloc = max_ints;
  n_me_alloc   = max_me;
  n_me_stored  = 0;
  n_ints       = 0;
  Int_list  = new double*[n_ints_alloc];
  StrEn     = new  double[n_ints_alloc];
  iabJ_list = new     int[n_me_alloc];
  icdJ_list = new     int[n_me_alloc];
  for(int i=0; i<n_ints_alloc; ++i) {Int_list[i] = NULL;  StrEn[i] = 0.0;}
  for(int i=0; i<n_me_alloc; ++i)   {iabJ_list[i] = -1;  icdJ_list[i] = -1;}
  
  if (n_new_ints > 0) for(int i=0; (i<n_new_ints)&&(i<n_ints_alloc); ++i) Allocate_single_Vpp();
  return;}

void  VppInt_t::Allocate_single_Vpp(void ) {
  if (n_ints >= n_ints_alloc) return;
  Int_list[n_ints] = new double[n_me_alloc];
  StrEn[n_ints] = 0.0;
  for(int i=0; i<n_me_alloc; ++i) Int_list[n_ints][i] = 0.0;
  ++n_ints;
  return;}

void  VppInt_t::Deallocate_single_Vpp(int i_cut){
  if ((i_cut < 0)||(i_cut >= n_ints_alloc)) {
    cerr <<"\nARNING (VppInt_t::Deallocate_single_Vpp): could not remove i="<<i_cut
         <<" (only "<<n_ints_alloc<<"are allocated)\n\n";
    return;}

  if (Vpp_wrk == Int_list[i_cut]) Vpp_wrk = NULL;
  delete [] Int_list[i_cut];
  --n_ints;
  for(int i=i_cut; i<n_ints; ++i) {Int_list[i]=Int_list[i+1]; StrEn[i]=StrEn[i+1];}
  StrEn[n_ints] = 0.0;
  Int_list[n_ints] = NULL;
  return;}


int VppInt_t::read(char *fname, int n_all_new_ints, int n_me_min) {
                                  // defualts: n_all_new_ints=1
                                  //           n_me_min=0
  //cout << fname+(sizeof(fname)-3) << endl;
  if (strstr(fname,".bin")!=NULL) return read_bin(fname);
  if (strstr(fname,".mhj")!=NULL) return read_Oslofmt(fname, n_all_new_ints-1, n_me_min);
  //if (strstr(fname+sizeof(fname)-5,".bin")!=NULL) return read_bin(fname);
  //if (strstr(fname,".ant")!=NULL) return read_antoine(fname);
  return read_ascii(fname, n_all_new_ints, n_me_min);
  }

int VppInt_t::write(char *fname) {
  if (strstr(fname,".bin")!=NULL) return write_bin(fname);
  //if (strstr(fname+sizeof(fname)-5,".bin")!=NULL) return write_bin(fname);
  return write_ascii(fname);
  }


int VppInt_t::write_ascii(char* fname) {
  cout << " Writing Vpp mtx elements to file '" << fname << "'... \n";

  FILE* outfile;
  outfile = fopen(fname, "w");

  fprintf(outfile," %i  %i  # tot number of mtx el. stored and to be allocated; format is: na 2ja ipa cha nb... chd J Vpp\n",n_me_stored, n_me_alloc);

  for(int i=0; i<n_me_stored; ++i) {
    iabJ = iabJ_list[i];
    icdJ = icdJ_list[i];
    
    J = J_i1i2_ind( iorb,  (iorb+1),  iabJ);
    if (J != J_i1i2_ind( (iorb+2),  (iorb+3),  icdJ))
                           cout << "Error with the coding of J \n\n";

    for(i2=0; i2<4; ++i2) {
      k = iorb[i2];
      ish = MdSp->Mor_sh[k];
      fprintf(outfile,"%i %i %i %i ",
         MdSp->Mor_n[k],MdSp->MSp_2j[ish],MdSp->MSp_ip[ish],MdSp->MSp_ch[ish]);
      if (1 == i2) fprintf(outfile," ");
    }
    fprintf(outfile," %i %.10lf\n",J,Vpp_wrk[i]);
  }
  fclose(outfile);

  cout << " --> Vpp matrix elements have been successfully written!\n\n";

  return 0;
  }

int VppInt_t::read_line_Vpp(char *lptr, int *ia, int *ib, int *ic, int *id, int *iJ, double *xV) {
  
  sscanf(lptr,"%i%i%i%i%n",&nbs,&j2bs,&ipbs,&chbs,&n1);
  lptr += n1;
  if ((*ia = MdSp->get_orbit(nbs,j2bs,ipbs%2,chbs)) < 0 ) return -1;
  //
  sscanf(lptr,"%i%i%i%i%n",&nbs,&j2bs,&ipbs,&chbs,&n1);
  lptr += n1;
  if ((*ib = MdSp->get_orbit(nbs,j2bs,ipbs%2,chbs)) < 0 ) return -1;
  //
  sscanf(lptr,"%i%i%i%i%n",&nbs,&j2bs,&ipbs,&chbs,&n1);
  lptr += n1;
  if ((*ic = MdSp->get_orbit(nbs,j2bs,ipbs%2,chbs)) < 0 ) return -1;
  //
  sscanf(lptr,"%i%i%i%i%n",&nbs,&j2bs,&ipbs,&chbs,&n1);
  lptr += n1;
  if ((*id = MdSp->get_orbit(nbs,j2bs,ipbs%2,chbs)) < 0 ) return -1;
  //
  sscanf(lptr,"%i%lf",iJ,xV);

  return 0;}

static char line[BUFFER_SIZE];
static int icount, count_outside, count_err, count_bad, count_discarded, count_added;

int VppInt_t::read_ascii(char *fname, int n_all_new_ints, int n_me_min) {
                                  // defualts: n_all_new_ints=1
                                  //           n_me_min=0

  if (NULL == MdSp) cerr << "\nVppInt_t: Cannot read Vpp mtx els. if no mod. space is specified.\n";

  cout << " Get Vpp mtx elements from file '" << fname << "'...   ";

  FILE* infile;
  infile = fopen(fname, "rb");

  fgets(line,sizeof(line), infile);  //" %i  # tot number of mtx els. on file"
  sscanf(line,"%i",&n1);

  cout << " " << n1 << " lines to be read.\n";


  if (n_all_new_ints < 1) n_all_new_ints = 1;
  if (n_me_min > n1) n1 = n_me_min;
  Allocate_Vpp_table(n_all_new_ints, n1, 1); // set for n_all_new_ints 2b-me
                                             // (with max. n1 mtx.els) but
                                             // allocate only 1 ints.
  select(0);
  n_me_stored = n1;
  //n_me_alloc = n_me_stored;
  //free_mem();
  //iabJ_list = new    int[n_me_alloc];
  //icdJ_list = new    int[n_me_alloc];
  //Vpp_wrk   = new double[n_me_alloc];



  icount=0; count_outside=0; count_err=0;
  for(int i=0; i<n_me_stored; i++) {
     fgets(line,sizeof(line), infile);  //" n_a  2j_a  ip_a  ch_a  nb ...
     if (read_line_Vpp(line,&iorb_a,&iorb_b,&iorb_c,&iorb_d,&J,&xe) ) 
                    {//cerr << "ERROR:  The following line:\n\n" << line
                     //     << "\ncontains an orbit not in the model space\n";
                     ++count_outside;
                     continue;}

   if (( (iorb_a==iorb_b)&&(abs(MdSp->MSp_2j[MdSp->Mor_sh[iorb_a]]-J)%2 ==0) )||
       ( (iorb_c==iorb_d)&&(abs(MdSp->MSp_2j[MdSp->Mor_sh[iorb_c]]-J)%2 ==0) )){
     cout << "ERROR:  The following line:\n\n" << line
          << "\ncontains a pauli forbidden state for fermions\n"
          << " --> it will be discaded\n\n";
     ++count_err;
     continue;
     }

     iph = 0;
     if (iorb_a <= iorb_b) {
       iabJ_list[icount]=ind_i1i2J(iorb_a, iorb_b, J); }
     else {
       ish = MdSp->Mor_sh[iorb_a];
       nbs = MdSp->MSp_2j[ish];
       ish = MdSp->Mor_sh[iorb_b];
       nbs += MdSp->MSp_2j[ish];
       iph += (nbs/2 - J +1);
       iabJ_list[icount]=ind_i1i2J(iorb_b, iorb_a, J); }

     if (iorb_c <= iorb_d) {
       icdJ_list[icount]=ind_i1i2J(iorb_c, iorb_d, J); }
     else {
       ish = MdSp->Mor_sh[iorb_c];
       nbs = MdSp->MSp_2j[ish];
       ish = MdSp->Mor_sh[iorb_d];
       nbs += MdSp->MSp_2j[ish];
       iph += (nbs/2 - J +1);
       icdJ_list[icount]=ind_i1i2J(iorb_d, iorb_c, J); }

     if (iabJ_list[icount] > icdJ_list[icount]) {        ipbs=iabJ_list[icount];
                                       iabJ_list[icount]=icdJ_list[icount];
                                       icdJ_list[icount]=ipbs;}

     if (abs(iph)%2 != 0) xe=-xe;
     Vpp_wrk[icount]=xe;
     ++icount;
  }
  n_me_stored=icount;

  fclose(infile);

  //for(int i=0; i<n_me_alloc; ++i)
  //   cout << i << "  " << iabJ_list[i] << "  " << icdJ_list[i] << "  " <<Vpp_wrk[i] << endl;

  prepare_for_reading();

  if (count_outside) cout << count_outside << " lines were out of the model space\n";
  if (count_err) cout << count_err << " lines had incompatible quantum numbers\n";
  cout << "  --> " << n_me_stored << " Vpp matrix elements have been loaded successfully!\n\n";

  return 0;
  }


int VppInt_t::write_bin(char* fname) {
  cout << " Write out Vpp mtx elements in binary format to the file '" << fname << "'...";

  int isel=-100,isel3NI=-100;
  for(int i=0; i<n_ints; ++i)  if (Vpp_wrk == Int_list[i]) isel=i;
  cout << "   (isel="<<isel<<")    ";
  // Later on must substitute wit this keep both iself...  it will imply remaking all the .bin files... so I'll avoid it for now.
  //cout << "   (isel="<<isel<<" ,  isel3NI="<<isel3NI<<")    ";
  
  ofstream fout(fname, ios::out|ios::trunc|ios::binary);
  fout.write((const char*) &n_me_stored,  sizeof(int));
  fout.write((const char*) &n_me_alloc,   sizeof(int));
  fout.write((const char*) &n_ints,       sizeof(int));
  fout.write((const char*) &n_ints_alloc, sizeof(int));
  fout.write((const char*) &isel,         sizeof(int));
  // Teh following is here only for future compatibility...
  fout.write((const char*) &isel3NI, sizeof(int));
  fout.write((const char*) StrEn,       n_ints*sizeof(double));
  fout.write((const char*) iabJ_list,   n_me_stored*sizeof(int));
  fout.write((const char*) icdJ_list,   n_me_stored*sizeof(int));
  for(int i=0; i<n_ints; ++i) fout.write((const char*) Int_list[i], n_me_stored*sizeof(double));
  fout.close();

  cout << " ...done!\n\n";

  return 0;
  }


int VppInt_t::read_bin(char *fname) {
  if (NULL == MdSp) cerr << "\nVppInt_t: Cannot read Vpp mtx els if not mod. space is specified.\n";

  cout << " Get Vpp mtx elements from file '" << fname << "' in binary format...   ";

  int dims[6];

  ifstream fin(fname, ios::in|ios::binary);
  fin.read((char*) dims,      6*sizeof(int));  // The last of these numbers is unsed for now but it still need to be read...

  Allocate_Vpp_table(dims[3], dims[1], dims[2]);
  n_me_stored = dims[0];
  select(dims[4]); // this already checks that dims[4] is within the bounds...

  fin.read((char*) StrEn,     n_ints*sizeof(double));
  fin.read((char*) iabJ_list, n_me_stored*sizeof(int));
  fin.read((char*) icdJ_list, n_me_stored*sizeof(int));
  for(int i=0; i<n_ints; ++i) fin.read((char*) Int_list[i], n_me_stored*sizeof(double));
  fin.close();

  cout << "    isel="  << get_wrk();
  if (get_wrk() < 0) {cout << " (-->reset to isel=0)"; select(0);}
  cout << "   ...done!\n\n";

  /* /TST//
  for(int i=0; i<n_ints; ++i) cout << Int_list[i] << "  ";
  cout << "  ---  " << Vpp_wrk << endl;
  //TST/ */

  return 0;
  }


// orders the mtx. els. and builds the ch,ip,J table
void VppInt_t::prepare_for_reading(void ) {

  int* i_chipJ; 
  i_chipJ = new    int[n_me_stored]; // this memory must be freed leaving the routine

  for(int i=0; i<n_me_stored; i++) {
    iabJ = iabJ_list[i];
    J = J_i1i2_ind(&i1, &i2, iabJ);
    ish1 = MdSp->Mor_sh[i1];
    ish2 = MdSp->Mor_sh[i2];
    
    ch =  MdSp->MSp_ch[ish1] + MdSp->MSp_ch[ish2];
    ip = (MdSp->MSp_ip[ish1] + MdSp->MSp_ip[ish2])%2;

      // test that the tot q.#s of bra and ket correspond (never hurts...):
      icdJ = icdJ_list[i];
      n1 = J_i1i2_ind(&i1, &i2, icdJ);
      ish1 = MdSp->Mor_sh[i1];
      ish2 = MdSp->Mor_sh[i2];
      if ((J  != n1) ||
          (ch != ( MdSp->MSp_ch[ish1] + MdSp->MSp_ch[ish2]   ) ) ||
          (ip != ((MdSp->MSp_ip[ish1] + MdSp->MSp_ip[ish2])%2) ) ) {
           cout << "\n\nERROR: J , ip and/or ch are not conserved for some of"
                << "the given Vpp matrix elements!\n"
                << "Type something to continue:";
           cin >> i1; // as a 'pause' statement...
           }
      //         

    i_chipJ[i] =  ind_chipJ(ch,ip,J) ;     
    }


  if (n_ints < 1)                            Asc_iii_iiiheapsort(           n_me_stored, i_chipJ, iabJ_list, icdJ_list);
  else if ((n_ints == 1)&&(NULL != Vpp_wrk)) Asc_iiid_iiiheapsort_singleV(  n_me_stored, i_chipJ, iabJ_list, icdJ_list, Vpp_wrk);
                                        else Asc_iiid_iiiheapsort_multipleV(n_me_stored, i_chipJ, iabJ_list, icdJ_list, n_ints, Int_list);
 

  if (NULL != Jipch_str) delete [] Jipch_str;
  if (NULL != Jipch_end) delete [] Jipch_end;
  Jipch_str = new int[CHMULT*JMULT*2];
  Jipch_end = new int[CHMULT*JMULT*2];
  for(int i=0; i<CHMULT*JMULT*2; i++) {Jipch_str[i]=-10;Jipch_end[i]=-10;}

  int i = 0;
  itmp = i_chipJ[i];
  J = J_chip_ind(&ch, &ip, itmp);
  Jipch_str[ ind_chipJ(ch,ip,J) ] = i;
  for(i=1; i<n_me_stored; i++) {
    if (itmp != i_chipJ[i]) {
      Jipch_end[ ind_chipJ(ch,ip,J) ] = i-1;
      itmp = i_chipJ[i];
      J = J_chip_ind(&ch, &ip, itmp);
      Jipch_str[ ind_chipJ(ch,ip,J) ] = i;
      }
    }
    Jipch_end[ ind_chipJ(ch,ip,J) ] = n_me_stored-1;

  delete [] i_chipJ;

  return;
  }




//int VppInt_t::add(int ia, int ib, int ic, int id, int J, double xV) {
int VppInt_t::add(int ia, int ib, int ic, int id, int J, double *xV, int nxe, int ist) {
                            // defaults: nxe=-1 int ist=0
                            //           nxe < 1  :  add to Vpp_wrk only
                            //           otherwise:  add nxe mtxels starting
                            //                      at Int_list[ist]
  if (n_me_stored >= n_me_alloc) {
    cout << "\n ERROR in VppInt_t::add, the Vpp array are already full and cannot add new mtx. els. !!!!\n";
    cout << "  --> WILL discard this entry....\n";
    return 1;
    }

  // make sure not to exceed the number of ints.
  nxe = (nxe < n_ints-ist) ? nxe : n_ints-ist;

  isha = MdSp->Mor_sh[ia];
  ishb = MdSp->Mor_sh[ib];
  ishc = MdSp->Mor_sh[ic];
  ishd = MdSp->Mor_sh[id];
  if ( (MdSp->MSp_ip[isha] + MdSp->MSp_ip[ishb] +  MdSp->MSp_ip[ishc] + MdSp->MSp_ip[ishd])%2 ) return -1;
  if (  MdSp->MSp_ch[isha] + MdSp->MSp_ch[ishb] != MdSp->MSp_ch[ishc] + MdSp->MSp_ch[ishd]   )  return -2;

  if ((ia == ib) && ( abs(MdSp->MSp_2j[isha] - J + 1)%2 ) ) return -3;
  if ((ic == id) && ( abs(MdSp->MSp_2j[ishc] - J + 1)%2 ) ) return -4;

  iph = 0;
  if (ia <= ib) {iabJ = ind_i1i2J(ia, ib, J);}
           else {iabJ = ind_i1i2J(ib, ia, J);
                 iph += (MdSp->MSp_2j[isha] + MdSp->MSp_2j[ishb])/2 - J +1;}

  if (ic <= id) {icdJ = ind_i1i2J(ic, id, J);}
           else {icdJ = ind_i1i2J(id, ic, J);
                 iph += (MdSp->MSp_2j[ishc] + MdSp->MSp_2j[ishd])/2 - J +1;}

  if (icdJ < iabJ) {itr=iabJ; iabJ=icdJ; icdJ=itr;}

  iabJ_list[n_me_stored] = iabJ;
  icdJ_list[n_me_stored] = icdJ;
  for(int l=0; l<n_ints; ++l) Int_list[l][n_me_stored] =  0.0;
  if (nxe < 1) {
    if (NULL == Vpp_wrk) {return -5;} else {
             Vpp_wrk[n_me_stored] = (*xV) * double(1 - 2*(abs(iph)%2) );}
  } else {
  for(int l=0; l<nxe; ++l)
    Int_list[ist+l][n_me_stored] =  xV[l]  * double(1 - 2*(abs(iph)%2) );
  }
  ++n_me_stored;

  return 0;}


double VppInt_t::get(int ia, int ib, int ic, int id, int J, int *imid) {
//
// NOTE if this routine is changed, the same corrections should be made
//   to 'VppInt_t::get_ph' as well.
//

  isha = MdSp->Mor_sh[ia];
  ishb = MdSp->Mor_sh[ib];
  ip = (MdSp->MSp_ip[isha] + MdSp->MSp_ip[ishb])%2;
  ch =  MdSp->MSp_ch[isha] + MdSp->MSp_ch[ishb];


  iph = 0;
  if (ia <= ib) {iabJ = ind_i1i2J(ia, ib, J);}
           else {iabJ = ind_i1i2J(ib, ia, J);
                 iph += (MdSp->MSp_2j[isha] + MdSp->MSp_2j[ishb])/2 - J +1;}

  if (ic <= id) {icdJ = ind_i1i2J(ic, id, J);}
           else {icdJ = ind_i1i2J(id, ic, J);
                 ishc = MdSp->Mor_sh[ic];
                 ishd = MdSp->Mor_sh[id];
                 iph += (MdSp->MSp_2j[ishc] + MdSp->MSp_2j[ishd])/2 - J +1;}

  xph = double(1 - 2*(abs(iph)%2) );
  if (icdJ < iabJ) {itr=iabJ; iabJ=icdJ; icdJ=itr;}


  idown = Jipch_str[ ind_chipJ(ch,ip,J) ];
  iup   = Jipch_end[ ind_chipJ(ch,ip,J) ];
  if ((idown < 0)||(iup < 0)) {*imid = -1100; return 0.0;}

  itr = 0;
  *imid = -1000;
  while(itr < 100) {// 100 iterations are enough for up to about 1024^10(=10^30) mtx elements
    *imid  = (iup+idown)/2;
    iab_mid = iabJ_list[*imid];
  //icd_mid = icdJ_list[*imid];
    test = (iabJ < iab_mid);
    if (iab_mid == iabJ) {icd_mid = icdJ_list[*imid];
                          if (icd_mid == icdJ) {return (Vpp_wrk[*imid] * xph);}
                            else               {test = (icdJ < icd_mid);}
                         }
    if (test) {iup   = *imid;} 
         else {idown = *imid; if (idown<iup) idown++;}
                                              
    itr++;
    }

  *imid = -1000;
  return 0.0;
  }


double VppInt_t::get_ph(int ia, int ib, int ic, int id, int J, int *imid) {
//
// NOTE if this routone is changed, the same corrections should be made
//   to 'VppInt_t::get' as well.
//

  isha = MdSp->Mor_sh[ia];
  ishb = MdSp->Mor_sh[ib];
  ip = (MdSp->MSp_ip[isha] + MdSp->MSp_ip[ishb])%2;
  ch =  MdSp->MSp_ch[isha] + MdSp->MSp_ch[ishb];


  iph = 0;
  if (ia <= ib) {iabJ = ind_i1i2J(ia, ib, J);}
           else {iabJ = ind_i1i2J(ib, ia, J);
                 iph += (MdSp->MSp_2j[isha] + MdSp->MSp_2j[ishb])/2 - J +1;}

  if (ic <= id) {icdJ = ind_i1i2J(ic, id, J);}
           else {icdJ = ind_i1i2J(id, ic, J);
                 ishc = MdSp->Mor_sh[ic];
                 ishd = MdSp->Mor_sh[id];
                 iph += (MdSp->MSp_2j[ishc] + MdSp->MSp_2j[ishd])/2 - J +1;}

  xph = double(1 - 2*(abs(iph)%2) );
  if (icdJ < iabJ) {itr=iabJ; iabJ=icdJ; icdJ=itr;}


  idown = Jipch_str[ ind_chipJ(ch,ip,J) ];
  iup   = Jipch_end[ ind_chipJ(ch,ip,J) ];
  if ((idown < 0)||(iup < 0)) {*imid = -1100; return xph;}

  itr = 0;
  *imid = -1000;
  while(itr < 100) {// 100 iterations are enough for up to about 1024^10(=10^30) mtx elements
    *imid  = (iup+idown)/2;
    iab_mid = iabJ_list[*imid];
  //icd_mid = icdJ_list[*imid];
    test = (iabJ < iab_mid);
    if (iab_mid == iabJ) {icd_mid = icdJ_list[*imid];
                          if (icd_mid == icdJ) {return xph;}
                            else               {test = (icdJ < icd_mid);}
                         }
    if (test) {iup   = *imid;} 
         else {idown = *imid; if (idown<iup) idown++;}
                                              
    itr++;
    }

  *imid = -1000;
  return xph;
  }


int VppInt_t::check_for_double_entries(void ) {
  int irep;
  cout << "\n Checking for doubly defined Vpp mtx. els....\n";

  prepare_for_reading();

  //cout << "n_me_stored=" << n_me_stored << " ,    n_me_alloc=" << n_me_alloc << endl;

  irep = 0;
  for(int i=0; i < n_me_stored-1; ++i) {
        //if (i%10000 == 0)
        // cout << n_me_stored << "  " << i << endl;
    if ((iabJ_list[i] == iabJ_list[i+1]) &&
        (icdJ_list[i] == icdJ_list[i+1])) {
        ++irep;
        // cout << iabJ_list[i] << "  " << icdJ_list[i] << "  " <<Vpp_wrk[i] << endl;
        if (fabs(Vpp_wrk[i]-Vpp_wrk[i+1]) > fabs(Vpp_wrk[i]/1.e7) )
          cerr << "WARNING: Vpp mtx. el. with Jab=" << iabJ_list[i]
               << " and Jcd=" << icdJ_list[i]
               << " has been given twice\n with values: " << Vpp_wrk[i]
               << " and " << Vpp_wrk[i+1]
               << "\n\n --> The second entry will be neglected.\n\n";
        for(int j=i+1; j < n_me_stored-1; ++j) {
          iabJ_list[j] = iabJ_list[j+1];
          icdJ_list[j] = icdJ_list[j+1];
          for(int l=0; l<n_ints; ++l) Int_list[l][j] =  Int_list[l][j+1]; // old: Vpp_wrk[j] =  Vpp_wrk[j+1];
          }
        --i;
        --n_me_stored;
        }
    }

  if (irep) cout << " Warning: " << irep << " 2-body matrix element(s) was(were)repeated, "
                 << "\n  --> the tot. number of mtx. els. is now set to = " << n_me_stored << endl << endl;


  prepare_for_reading();


  return irep;}

static int nme_file, ips;

int VppInt_t::comp_Vpp_file(char *fname) {
  if (NULL == MdSp) {
        cerr << "\nVppInt_t::comp_Vpp_file('" <<  fname
        << "') : the object is not even initialized (no mod. space is specified).\n";
        return -1;
        }

  prepare_for_reading();


  cout << " Comparing Vpp mtx elements with those contained in the file '" << fname << "'...   ";

  FILE* infile;
  infile = fopen(fname, "rb");



  fgets(line,sizeof(line), infile);  //" %i  # tot number of mtx el. stored"
  sscanf(line,"%i",&nme_file);

  cout << " " << nme_file << " lines to be read.\n";

  icount=0; count_outside=0; count_err=0; count_bad=0;
  for(int i=0; i<nme_file; i++) {
     fgets(line,sizeof(line), infile);  //" n_a  2j_a  ip_a  ch_a  nb ...
     if (read_line_Vpp(line,&iorb_a,&iorb_b,&iorb_c,&iorb_d,&J,&xe) ) 
                    {//cerr << "ERROR:  The following line:\n\n" << line
                     //     << "\ncontains an orbit not in the model space\n";
                     ++count_outside;
                     continue;}
     
 if (( (iorb_a == iorb_b) && (abs(MdSp->MSp_2j[MdSp->Mor_sh[iorb_a]]-J)%2 ==0) )
  ||( (iorb_c == iorb_d) && (abs(MdSp->MSp_2j[MdSp->Mor_sh[iorb_c]]-J)%2 ==0) ))
   {cerr << "ERROR:  The following line:\n\n" << line
         << "\ncontains a pauli forbidden state for fermions\n"
         << " --> it will be duscaded\n\n";
     ++count_err;
     continue;
     }

     if (fabs(xe-get(iorb_a,iorb_b,iorb_c,iorb_d,J,&ips)) > CHOP) {
     //if (xe != get(iorb_a,iorb_b,iorb_c,iorb_d,J,&ips))
       if (ips >= 0) {cout << "the following matrix elments do not correspond (ips=" << ips << "): ";}
       else {cout << "the following matrix elment is missing in memory: ";}
       cout << MdSp->Mor_n[iorb_a] << MdSp->MSp_name[MdSp->Mor_sh[iorb_a]] << "  " 
            << MdSp->Mor_n[iorb_b] << MdSp->MSp_name[MdSp->Mor_sh[iorb_b]] << "  " 
            << MdSp->Mor_n[iorb_c] << MdSp->MSp_name[MdSp->Mor_sh[iorb_c]] << "  " 
            << MdSp->Mor_n[iorb_d] << MdSp->MSp_name[MdSp->Mor_sh[iorb_d]] << "  J="
            << J << "  xe: " << xe;
       if (ips>=0) cout << " != " << get(iorb_a,iorb_b,iorb_c,iorb_d,J,&ips);
       cout << endl;
       ++count_bad;}

     icount++;
  }

  fclose(infile);


  cout << icount << " lines out of " << nme_file << " Vpp matrix elements present in file '"
       << fname  << "' and have been compared,\n";
  if (count_bad) cout << count_bad << " differeces were found\n"; 
  if (count_outside) cout << count_outside << " lines were out of the model space\n";
  if (count_err) cout << count_err << " lines had incompatible quantum numbers\n";
  cout << endl;

  return count_bad;
  }



int VppInt_t::add_file(double C, char *fname, int n_xe/*=-1*/, int ist/*=0*/) {
                            //           n_xe < 1  :  add to Vpp_wrk only
                            //           otherwise:  add n_xe mtxels starting
                            //                      at Int_list[ist]
  if (NULL == MdSp) {
        cerr << "\nVppInt_t::add_file('" <<  fname
        << "') : the object is not even initialized (no mod. space is specified).\n";
        return -1;
        }
  if (NULL == Vpp_wrk) {cerr << "\n(VppInt_t::add_file): CANNOT WORK WITH Vpp_wrk==NULL, consider updating the routine.\n\n-->STOP.\n";
                        exit(100);}

  // make sure not to exceed the number of ints.
  n_xe = (n_xe < n_ints-ist) ? n_xe : n_ints-ist;
  
  cout << " Adding Vpp mtx elements from file '" << fname << "' (C=" << C << ") ...    " << flush;

  double xVpp;
  FILE* infile;
  infile = fopen(fname, "rb");
  if (NULL == infile) {
    cout << "\n WARNING: cannot file '"<<fname<<"'...   will skip it!!\n";
    cerr << "\n WARNING: cannot file '"<<fname<<"'...   will skip it!!\n";
    return 1;
  }

  fgets(line,sizeof(line), infile);  //" %i  # tot number of mtx el. stored"
  sscanf(line,"%i",&nme_file);

  cout << " " << nme_file << " lines to be read.\n";

  icount=0; count_outside=0; count_err=0; count_discarded=0; count_added=0;
  for(int i=0; i<nme_file; i++) {
     fgets(line,sizeof(line), infile);  //" n_a  2j_a  ip_a  ch_a  nb ...
     if (read_line_Vpp(line,&iorb_a,&iorb_b,&iorb_c,&iorb_d,&J,&xe) ) 
                    {//cerr << "ERROR:  The following line:\n\n" << line
                     //     << "\ncontains an orbit not in the model space\n";
                     ++count_outside;
                     continue;}
     
     if (( (iorb_a == iorb_b)&&((MdSp->MSp_2j[MdSp->Mor_sh[iorb_a]]+J)%2 ==0) )
       ||( (iorb_c == iorb_d)&&((MdSp->MSp_2j[MdSp->Mor_sh[iorb_c]]+J)%2 ==0) ))
         {cerr << "ERROR:  The following line:\n\n" << line
               << "\ncontains a pauli forbidden state for fermions\n"
               << " --> it will be duscaded\n\n";
          ++count_err;
          continue;
          }


     xVpp = get(iorb_a,iorb_b,iorb_c,iorb_d,J,&ips);
     x1 = C * xe;
     if (ips>=0) {
       if (0.0 == xVpp) {
         // if xVpp == 0.0, set it momentarily to a non-zero
         // number to extract the phase...
         Vpp_wrk[ips] = 1.0;
         xVpp = get(iorb_a,iorb_b,iorb_c,iorb_d,J,&ips);
         if (xVpp * Vpp_wrk[ips] < 0) x1 = -x1;
         Vpp_wrk[ips] = 0.0;
       } 
       else if (xVpp * Vpp_wrk[ips] < 0) x1 = -x1;

       if (n_xe < 0) {Vpp_wrk[ips] += x1;}
       else { for(int l=0; l<n_xe; ++l) Int_list[ist+l][ips] += x1;}
       

       icount++;
     } else {
       // it's a good mtx el. contained in the file but not in VppInt_t
       if (add(iorb_a,iorb_b,iorb_c,iorb_d,J, &x1, n_xe, ist) != 0) {++count_discarded;} else {++icount; ++count_added;}
     }

  }

  fclose(infile);


  if (count_added) prepare_for_reading();

  cout << icount << " lines (out of " << nme_file << " Vpp matrix elements present in file '"
       << fname  << "') have been summed to the existing interaction ("<< n_me_stored <<" stored),\n";
  if (count_added) cout << count_added << " corresponed to new configurations\n"; 
  if (count_discarded) cout << count_discarded << " were discarded because could not be added the existing object\n"; 
  if (count_outside) cout << count_outside << " lines were out of the model space\n";
  if (count_err) cout << count_err << " lines had incompatible quantum numbers\n";
  cout << endl;

  return 0;
  }


//
// loop ofver all possible mtx els and chek if they are in memory
//
int VppInt_t::seek_missing_Vpp(int iadd_zeroes) {
                        // iadd_zeroes == 0 [default] , just check 
                        //             == 1 add zeroes for the missing configs.

  if (NULL == MdSp) {
     cerr << "\nVppInt_t::seek_missing_Vpp(): interaction is not even initialized (no mod. space is specified).\n";
     return -1;
     }

/*
  *n_ph   = 0;
  *n_pp   = 0;
  *n_hh   = 0;
  *n_2p1h = 0;
  *n_2h1p = 0;
*/
//  int n_chs = 0;

//  int nst,nJ,i2b;


  MdSp->Calculate_Jch_bounds();

  cout << " Seeking for missing two-body matrix elements...    ";


  int ia,  ib,  ic,  id,  nJab,nJcd,ipp_ab,ipp_cd,nst;
  int nash,nbsh,ncsh,ndsh;
  int J,ip,ch;
  int n_2bme;
  double xzero = 0.0;


  n_2bme=0;
  count_bad = count_added = 0;
  for(J=0;         J<=MdSp->Jmax;       ++J) {
  for(ip=0;        ip<2;         ++ip) {
  for(ch=MdSp->ch_pp_mn; ch<=MdSp->ch_pp_mx; ++ch) {

    ipp_ab=0;
    for(ia=0 ; ia < MdSp->nsubsh; ++ia) {
    for(ib=ia; ib < MdSp->nsubsh; ++ib) {

      nJab = ((MdSp->MSp_2j[ia]+MdSp->MSp_2j[ib])/2 + J + 1)%2;

      if ( MdSp->MSp_ch[ia]+MdSp->MSp_ch[ib]       != ch) continue;
      if ((MdSp->MSp_ip[ia]+MdSp->MSp_ip[ib]+ip)%2 != 0 ) continue;
//      if (AM.triang(MdSp->MSp_2j[ia],MdSp->MSp_2j[ib],2*J) <= 0) continue;
      if ((MdSp->MSp_2j[ia]+MdSp->MSp_2j[ib] < 2*J) ||
          (MdSp->MSp_2j[ib]+    2*J    < MdSp->MSp_2j[ia]) ||
          (    2*J   +MdSp->MSp_2j[ia] < MdSp->MSp_2j[ib]) ) continue;

      for(nash=0;   nash<MdSp->MSp_no[ia]; ++nash) {nst=0; if (ia==ib) nst=nash+nJab;
      for(nbsh=nst; nbsh<MdSp->MSp_no[ib]; ++nbsh) {
        ++ipp_ab;
        iorb_a = (MdSp->MSp_n[ia] - MdSp->Mor_n) + nash;
        iorb_b = (MdSp->MSp_n[ib] - MdSp->Mor_n) + nbsh;


        ipp_cd=0;
        for(ic=0 ; ic < MdSp->nsubsh; ++ic) {
        for(id=ic; id < MdSp->nsubsh; ++id) {

           nJcd = ((MdSp->MSp_2j[ic]+MdSp->MSp_2j[id])/2 + J + 1)%2;

           if ( MdSp->MSp_ch[ic]+MdSp->MSp_ch[id]       != ch) continue;
           if ((MdSp->MSp_ip[ic]+MdSp->MSp_ip[id]+ip)%2 != 0 ) continue;
//           if (AM.triang(MdSp->MSp_2j[ic],MdSp->MSp_2j[id],2*J) <= 0) continue;
           if ((MdSp->MSp_2j[ic]+MdSp->MSp_2j[id] < 2*J) ||
               (MdSp->MSp_2j[id]+    2*J    < MdSp->MSp_2j[ic]) ||
               (    2*J   +MdSp->MSp_2j[ic] < MdSp->MSp_2j[id]) ) continue;

           for(ncsh=0;   ncsh<MdSp->MSp_no[ic]; ++ncsh) {nst=0; if (ic==id) nst=ncsh+nJcd;
           for(ndsh=nst; ndsh<MdSp->MSp_no[id]; ++ndsh) {
             ++ipp_cd;
             if (ipp_cd < ipp_ab) continue;

             ++n_2bme;

             iorb_c = (MdSp->MSp_n[ic] - MdSp->Mor_n) + ncsh;
             iorb_d = (MdSp->MSp_n[id] - MdSp->Mor_n) + ndsh;

             get(iorb_a,iorb_b,iorb_c,iorb_d,J,&ips);
             if (ips < 0) {
               ++count_bad;
               if ( iadd_zeroes && (n_me_stored < n_me_alloc)) {
                 xzero = 0.0;
                 add(iorb_a,iorb_b,iorb_c,iorb_d,J,&xzero);
                 ++count_added;
               }
             }

           } } // ncsh,ndsh
         } } // ic , id


       } } // nash,nbsh
     } } // ia , ib
  
  //cout << "J=" << J << "  ip=" << ip << "  ch=(+/-)" << ch << " --> i2b= " << i2b << " / i2me=" << (i2b*(i2b+1))/2 << endl;
  } } } // J,ip,ch

  if (count_added) prepare_for_reading();

  if (count_bad)   cout << endl << endl << count_bad   << " 2b mtx els. are missing, out of a total of" << n_2bme  << endl  << endl;
  if (count_added) cout << endl << endl << count_added << " where added as zeroes" << endl;

  cout << "  ...check done.\n";

  return count_bad;}


void VppInt_t::scale_chrg(double C, int chrg_in){
  //
  //  Rescale matrix elements with a given charge and for the working interaction.
  //
  //  NOTE: ANY CHANGE made here MUST be made to the function VppInt_t::scaleall_chrg as well!!!
  //
  if (NULL == Vpp_wrk) return;
  //
  int ia,ib, chrg_loc;
  int ig,id, Jtst;
  for(int i=0; i<n_me_stored; ++i) {
    Jtst = J_i1i2_ind( &ia,  &ib,  iabJ_list[i]);
    chrg_loc = MdSp->MSp_ch[MdSp->Mor_sh[ia]]+MdSp->MSp_ch[MdSp->Mor_sh[ib]];
    
    // The following are extra safety tests... will slow down and are
    // a bit unnecessary but they do not hurt otherwise [CB 22.2.2015]
    if (Jtst != J_i1i2_ind( &ig,  &id,  icdJ_list[i])) cerr << "ERROR (VppInt_t::scale_chrg) with the coding of J \n\n";
    if (MdSp->MSp_ch[MdSp->Mor_sh[ig]]+MdSp->MSp_ch[MdSp->Mor_sh[id]] != chrg_loc)
      cerr << "ERROR (VppInt_t::scale_chrg) with the charges \n\n";
    
    if (chrg_loc == chrg_in) Vpp_wrk[i] *= C;
    
  }
  return;}

void VppInt_t::scaleall_chrg(double C, int chrg_in){
  //
  //  Rescale matrix elements with a given charge and for ALL interactions.
  //
  //  NOTE: ANY CHANGE made here MUST be made to the function VppInt_t::scale_chrg as well!!!
  //
  int ia,ib, chrg_loc;
  int ig,id, Jtst, j;
  for(int i=0; i<n_me_stored; ++i) {
    Jtst = J_i1i2_ind( &ia,  &ib,  iabJ_list[i]);
    chrg_loc = MdSp->MSp_ch[MdSp->Mor_sh[ia]]+MdSp->MSp_ch[MdSp->Mor_sh[ib]];
    
    // The following are extra safety tests... will slow down and are
    // a bit unnecessary but they do not hurt otherwise [CB 22.2.2015]
    if (Jtst != J_i1i2_ind( &ig,  &id,  icdJ_list[i])) cerr << "ERROR (VppInt_t::scaleall_chrg) with the coding of J \n\n";
    if (MdSp->MSp_ch[MdSp->Mor_sh[ig]]+MdSp->MSp_ch[MdSp->Mor_sh[id]] != chrg_loc)
      cerr << "ERROR (VppInt_t::scaleall_chrg) with the charges \n\n";
    
    if (chrg_loc == chrg_in) for(j=0; j<n_me_stored; ++j) Int_list[j][i] *= C;
    
  }
  return;}




static void Asc_iii_iiiheapsort(int ntot, int ia[],int ib[],int ic[]) {
  int    i,ir,j,l;
  int    iaex, ibex, icex;
  int test;

  if (ntot < 2) return;
  l = ntot/2;
  ir = ntot-1;
  while(1) {
    if (l > 0) {
       l--;
       iaex  = ia[l];
       ibex  = ib[l];
       icex  = ic[l];
    } else {
       iaex  = ia[ir];
       ibex  = ib[ir];
       icex  = ic[ir];
       ia[ir]  = ia[0];
       ib[ir]  = ib[0];
       ic[ir]  = ic[0];
       if ((--ir) == 0) {ia[0]  = iaex;
                         ib[0]  = ibex;
                         ic[0]  = icex;
                         break; }
    }
    i = l;
    j = l+l+1;
    while (j <= ir) {
       test=ia[j]  < ia[j+1];
       if  (ia[j] == ia[j+1]) {test=ib[j]  < ib[j+1];
                               if ( ib[j] == ib[j+1] ) test= ic[j] < ic[j+1];}
       if ((j < ir) && test) j++;
       test= iaex  < ia[j];
       if   (iaex == ia[j]) {test= ibex  < ib[j];
                             if   (ibex == ib[j]) test= icex  < ic[j];}
       if (test) {
         ia[i]  = ia[j];
         ib[i]  = ib[j];
         ic[i]  = ic[j];
         i = j;
         j = j+j+1;
         } else break;
       }
    ia[i]  = iaex;
    ib[i]  = ibex;
    ic[i]  = icex;
    }

  return;
  }


static void Asc_iiid_iiiheapsort_singleV(int ntot, int ia[],int ib[],int ic[],double dbl[]) {
  int    i,ir,j,l;
  int    iaex, ibex, icex;
  double dblex;
  int test;

  if (ntot < 2) return;
  l = ntot/2;
  ir = ntot-1;
  while(1) {
    if (l > 0) {
       l--;
       iaex  = ia[l];
       ibex  = ib[l];
       icex  = ic[l];
       dblex = dbl[l];
    } else {
       iaex  = ia[ir];
       ibex  = ib[ir];
       icex  = ic[ir];
       dblex = dbl[ir];
       ia[ir]  = ia[0];
       ib[ir]  = ib[0];
       ic[ir]  = ic[0];
       dbl[ir] = dbl[0];
       if ((--ir) == 0) {ia[0]  = iaex;
                         ib[0]  = ibex;
                         ic[0]  = icex;
                         dbl[0] = dblex;
                         break; }
    }
    i = l;
    j = l+l+1;
    while (j <= ir) {
       test=ia[j]  < ia[j+1];
       if  (ia[j] == ia[j+1]) {test=ib[j]  < ib[j+1];
                               if ( ib[j] == ib[j+1] ) test= ic[j] < ic[j+1];}
       if ((j < ir) && test) j++;
       test= iaex  < ia[j];
       if   (iaex == ia[j]) {test= ibex  < ib[j];
                             if   (ibex == ib[j]) test= icex  < ic[j];}
       if (test) {
         ia[i]  = ia[j];
         ib[i]  = ib[j];
         ic[i]  = ic[j];
         dbl[i] = dbl[j];
         i = j;
         j = j+j+1;
         } else break;
       }
    ia[i]  = iaex;
    ib[i]  = ibex;
    ic[i]  = icex;
    dbl[i] = dblex;
    }

  return;
  }


static void Asc_iiid_iiiheapsort_multipleV(int ntot, int ia[],int ib[],int ic[], int n_list, double* dbl_list[]) {
  int    i,ir,j,l, k;
  int    iaex, ibex, icex;
  double *dblex = new double[n_list];
  int test;

  if (ntot < 2) return;
  l = ntot/2;
  ir = ntot-1;
  while(1) {
    if (l > 0) {
       l--;
       iaex  = ia[l];
       ibex  = ib[l];
       icex  = ic[l];
       for(k=0; k<n_list; ++k) dblex[k] = dbl_list[k][l];
    } else {
       iaex  = ia[ir];
       ibex  = ib[ir];
       icex  = ic[ir];
       for(k=0; k<n_list; ++k) dblex[k] = dbl_list[k][ir];
       ia[ir]  = ia[0];
       ib[ir]  = ib[0];
       ic[ir]  = ic[0];
       for(k=0; k<n_list; ++k) dbl_list[k][ir] = dbl_list[k][0];
       if ((--ir) == 0) {ia[0]  = iaex;
                         ib[0]  = ibex;
                         ic[0]  = icex;
                         for(k=0; k<n_list; ++k) dbl_list[k][0] = dblex[k];
                         break; }
    }
    i = l;
    j = l+l+1;
    while (j <= ir) {
       test=ia[j]  < ia[j+1];
       if  (ia[j] == ia[j+1]) {test=ib[j]  < ib[j+1];
                               if ( ib[j] == ib[j+1] ) test= ic[j] < ic[j+1];}
       if ((j < ir) && test) j++;
       test= iaex  < ia[j];
       if   (iaex == ia[j]) {test= ibex  < ib[j];
                             if   (ibex == ib[j]) test= icex  < ic[j];}
       if (test) {
         ia[i]  = ia[j];
         ib[i]  = ib[j];
         ic[i]  = ic[j];
         for(k=0; k<n_list; ++k) dbl_list[k][i] = dbl_list[k][j];
         i = j;
         j = j+j+1;
         } else break;
       }
    ia[i]  = iaex;
    ib[i]  = ibex;
    ic[i]  = icex;
    for(k=0; k<n_list; ++k) dbl_list[k][i] = dblex[k];
    }

  delete [] dblex;
  return;
  }

