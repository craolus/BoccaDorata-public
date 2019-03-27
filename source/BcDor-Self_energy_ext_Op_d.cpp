//
//
//   Construct and store the self energy (both HF and dynamc part) for a given
//  channel taking as input:
//      model space, interaction and sp propagator.
//
//

//#include <iostream>
//#include <fstream>
#include <iomanip>

#include <cstdlib>
//#include <cstdio>
//#include <cmath>
using namespace std;

#include "BcDor-Global_variables.hh"
#include "BcDor-SlfEn_classes.hh"
//#include "gII_classes.hh"
//#include "BcDor-Utilities.hh"
//#include "BcDor-Ang_momenta.hh"


static int n1, n2, n3;

//----------------------------------------------------------------------------
//
//  Functions for storing the self-energy in binary format
//

static char filename[200];
static char filename1[200];
static char filename2[200];
static char strtmp[200];
static int  iload[40];  // need at least 20+15, but more do not hurt
static int *i_ptr;
static double *d1_ptr;
static char chip[2] = { '+' , '-'};
static int ierr, ierr1, ierr2, ierr3;
static int ifw, ibk;

void Wide_self_energy_t::get_filename(char *fname, char *Lbl, char *Ord, char *fb) {
  if ((0 == Lambda) && (0 == parity) && (0 == charge)) {
    sprintf(strtmp, "%s_%s", Lbl, Ord);
  } else {
    sprintf(strtmp, "%s_%s_T=%i%c_%i", Lbl, Ord, Lambda, chip[parity], charge);
  }
  sprintf(fname, "%s/%s_J2pi=%i%c_dt=%i%s.bin" , STOREDIR, strtmp,
          MdSp->MSp_2j[I_SH], chip[MdSp->MSp_ip[I_SH]], MdSp->MSp_ch[I_SH], fb);
  //cout << "File name set to: -->"<< fname<<"<--"<<endl<<flush;
  return; }


int Wide_self_energy_t::save_SelfEn_bin(int i_FwBk, char *Label/*="SelfEn"*/) {
  get_filename(filename, Label, "Dyn", "");
  cout << "\nOpening file '" << filename << "' for writing...\n" << flush;
  ofstream file(filename, ios::out|ios::trunc|ios::binary);
  if (file.bad() ) {
    cerr << "ERROR in opening the file '" << filename <<"' for writing... !!!" << flush;
    return 1;
  }
  
  int n1, n2, n3, n4;
  
  n1 = N_PLS_fw; if (i_FwBk < 0) n1 = 0;
  n2 = N_PLS_bk; if (i_FwBk > 0) n2 = 0;
  file.write((const char*) &n1            ,      sizeof(int));     // 0//
  file.write((const char*) &n2            ,      sizeof(int));
  file.write((const char*) &N_PLS_FW_alloc,      sizeof(int));
  file.write((const char*) &N_PLS_BK_alloc,      sizeof(int));
  file.write((const char*) &I_SH       ,         sizeof(int));     // 4//
  file.write((const char*) &NTOT_DEN   ,         sizeof(int));
  file.write((const char*) &NTOT_ORB   ,         sizeof(int));
  file.write((const char*) &NTOT_DIM   ,         sizeof(int));
  file.write((const char*) &N_EXT_PWV  ,         sizeof(int));     // 8//
  file.write((const char*) (MdSp->MSp_2j + I_SH ),  sizeof(int));  // 9//
  file.write((const char*) (MdSp->MSp_ip + I_SH ),  sizeof(int));  //10//
  file.write((const char*) (MdSp->MSp_ch + I_SH ),  sizeof(int));  //11//
  //
  //
  // Partial waves table
  file.write((const char*) i_pwv,         N_EXT_PWV*sizeof(int));
  file.write((const char*) nst_pwv,       N_EXT_PWV*sizeof(int));
  file.write((const char*) norb_pwv,      N_EXT_PWV*sizeof(int));
  file.write((const char*) isp_msp,       NTOT_ORB *sizeof(int));
  for(n1=0; n1<NTOT_ORB;  ++n1) {n2 = MdSp->Mor_n[isp_msp[n1]];
                                 file.write((const char*) &n2,    sizeof(int));}
  //
  //
  // Informaton on Lanczos blocks dimensions (all of this is for compatibility
  // with advanced verisons of BcDor):
  n1 = 0; n2 = 0; n3 = 0;  // This means a diagnal self-energy
  if (0 < n_piv_fw) {n2 = 1;  n4 = i_piv_fw[0] + 3;  n1 = (n1 > n4) ? n1 : n4;}
  if (0 < n_piv_bk) {n3 = 1;  n4 = i_piv_bk[0] + 3;  n1 = (n1 > n4) ? n1 : n4;}
  //
  file.write((const char*) &n1,     sizeof(int));  // max length of each piv array
  file.write((const char*) &n2,     sizeof(int));  // # of fw- blocks (==0 means diag. SE)
  file.write((const char*) &n3,     sizeof(int));  // # of bk- blocks (==0 means diag. SE)
  
  /* /TST//
  // Print-out for testing purposes
  cout << " xx- " << N_PLS_fw << "  " << N_PLS_bk << endl;
  cout << " xy- " << NTOT_DIM << "  " << n_piv_fw << "  " << n_piv_bk <<endl;
  cout << " xz- " << n1 << "  " << n2 << "  " << n3 <<endl;
  for (int i=0; i<NTOT_DIM; ++i) cout << i_piv_fw[i] << "  "; cout << endl;
  for (int i=0; i<NTOT_DIM; ++i) cout << i_piv_bk[i] << "  "; cout << endl;
  // */
  
  //
  // Forward self-energy:
  if (0 <= i_FwBk) {
    file.write((const char*) Sig_dyn_fw,   N_PLS_fw*NTOT_DIM*sizeof(double));
    if (0 < n_piv_fw) {
      n1 = i_piv_fw[0] + 2;
      file.write((const char*) &N_PLS_fw,                    sizeof(int));
      file.write((const char*) i_piv_fw,                  n1*sizeof(int));
    }
  }
  //
  // Backward self-energy:
  if (0 >= i_FwBk) {
    file.write((const char*) Sig_dyn_bk,   N_PLS_bk*NTOT_DIM*sizeof(double));
    if (0 < n_piv_bk) {
      n1 = i_piv_bk[0] + 2;
      file.write((const char*) &N_PLS_bk,                    sizeof(int));
      file.write((const char*) i_piv_bk,                  n1*sizeof(int));
    }
  }

  file.close();

  return 0;}


int  Wide_self_energy_t::load_SelfEn_bin_a(ifstream *filein, int ildin[]) {
  filein->read((char*) ildin, 12*sizeof(int));
  int i1 = 0;
  int n1, n2;
  if (ildin[4] != I_SH) {
    cout << "different subshell numbers...  "
         << I_SH << "(mod. sp.) <--> " << ildin[4] << "(file)\n";
    i1 = 1;
  }
  if ( (ildin[9]  != MdSp->MSp_2j[I_SH]) ||
       (ildin[10] != MdSp->MSp_ip[I_SH]) ||
       (ildin[11] != MdSp->MSp_ch[I_SH]) ) {
    cerr << "ERROR (Wide_self_energy_t::load_SelfEn_bin_a): j2fd, pifd, or chfd do not correspond to the ones on file\n";
    cout << ildin[9] << "  " << ildin[10] << "  " << ildin[11] << endl;
    cout << MdSp->MSp_2j[I_SH] << "  " << MdSp->MSp_ip[I_SH] << "  " << MdSp->MSp_ch[I_SH] << endl;
    i1 = 1;
  }
  //if (ildin[5] != NTOT_DEN ) {cout << "Mismatch for NTOT_DEN : "<< NTOT_DEN  << "(mod. sp.) <--> " << ildin[5] << "(file)\n"; i1=1;}
  if (ildin[6] != NTOT_ORB ) {cout << "Mismatch for NTOT_ORB : "<< NTOT_ORB  << "(mod. sp.) <--> " << ildin[6] << "(file)\n"; i1=1;}
  //if (ildin[7] != NTOT_DIM ) {cout << "Mismatch for NTOT_DIM : "<< NTOT_DIM  << "(mod. sp.) <--> " << ildin[7] << "(file)\n"; i1=1;}
  if (ildin[8] != N_EXT_PWV) {cout << "Mismatch for N_EXT_PWV: "<< N_EXT_PWV << "(mod. sp.) <--> " << ildin[8] << "(file)\n"; i1=1;}
  for(n1=0; n1<ildin[8]; ++n1) {
    filein->read((char*) &n2, sizeof(int));
    if (n2 != i_pwv[n1]) {
            cout << "Mismatch for the part. wave list (i_pwv): "<< i_pwv[n1]
                 << "(mod. sp.) <--> " << n2 << "(file)\n"; i1=1;}
  }
  for(n1=0; n1<ildin[8]; ++n1) {
    filein->read((char*) &n2, sizeof(int));
    if (n2 != nst_pwv[n1]) {
            cout << "Mismatch for the part. wave list (nst_pwv): "<< nst_pwv[n1]
                 << "(mod. sp.) <--> " << n2 << "(file)\n"; i1=1;}
  }
  for(n1=0; n1<ildin[8]; ++n1) {
    filein->read((char*) &n2, sizeof(int));
    if (n2 != norb_pwv[n1]) {
            cout << "Mismatch for the part. wave list (norb_pwv): "<< norb_pwv[n1]
                 << "(mod. sp.) <--> " << n2 << "(file)\n"; i1=1;}
  }
  for(n1=0; n1<ildin[6]; ++n1) {
    filein->read((char*) &n2, sizeof(int));
    if (n2 != isp_msp[n1]) {
            cout << "Mismatch for the orbit list (isp_msp): "<< isp_msp[n1]
                 << "(mod. sp.) <--> " << n2 << "(file)\n"; i1=1;}
  }
  for(n1=0; n1<ildin[6]; ++n1) {
    n3 = MdSp->Mor_n[isp_msp[n1]];
    filein->read((char*) &n2, sizeof(int));
    if (n2 != n3) {
            cout << "Mismatch for the 'n' q.# list (isp_msp): "<< n3
                 << "(mod. sp.) <--> " << n2 << "(file)\n"; i1=1;}
  }
  
  filein->read((char*) (ildin+12), 3*sizeof(int));
  //if (ildin[12] >= NTOT_PIV ) {cout << "Mismatch for NTOT_PIV : "<< NTOT_PIV  << "(mod. sp.) <--> " << ildin[12] << "(file)\n"; i1=1;}

  return i1;}

void Wide_self_energy_t::load_SelfEn_bin_b(ifstream *filein, double *dload, int nread, int nclear, int nspan) {
  //
  double *d2_ptr = dload + nclear*nspan;
  for(double *d1_ptr=dload; d1_ptr<d2_ptr; ++d1_ptr) (*d1_ptr)=0.0;
  filein->read((char*) dload,    nread*nspan*sizeof(double));
  return;}

void Wide_self_energy_t::load_SelfEn_bin_c(ifstream *filein, int *nrd1, int *ird2, int nclear) {
  //
  int nread;
  (*nrd1) = 0;
  for(int n1=0; n1<nclear; ++n1) ird2[n1]=0;
  filein->read((char*) nrd1,            sizeof(int));
  filein->read((char*) &nread,            sizeof(int));
  if ((nread+2) > nclear) {cerr << " BAD self-energy file! --> EXIT"; exit(100);}
  ird2[0] = nread;
  filein->read((char*) (ird2+1),  (nread+1)*sizeof(int));
  return;}

int  Wide_self_energy_t::load_SelfEn_bin(char *Label/*="SelfEn"*/,
                                     int nFw_min/*=-1*/, int nBk_min/*=-1*/) {
                         //  `nFw_min' and `nBk_min' are the minimum number of
                         // self-energy poles to be allocated (if the file
                         // contains less some emptyspace is left that can
                         // be filled out later...).
                         //
  int n1, n2;
  ierr = ierr1 = ierr2 = ierr3 = 0;
  get_filename(filename, Label, "Dyn", "");
  cout << "\nOpening file '" << filename << "' for reading...\n" << flush;
  ifstream file(filename, ios::in|ios::binary);
  ierr1 = file.bad();
  if ( ! ierr1 ) {
    //
    ierr = load_SelfEn_bin_a(&file,  iload);
    if (ierr) {cerr << "ERROR (Wide_self_energy_t::load_SelfEn_bin, file='"<<filename
                    <<"'): too many differencies, I'll skip reading the dynamic self-energy...\n\n";}
    if (ierr) {cout << "ERROR (Wide_self_energy_t::load_SelfEn_bin, file='"<<filename
                    <<"'): too many differencies, I'll skip reading the dynamic self-energy...\n\n";}
    //
    if ( ! ierr) {
      ifw = (iload[0] > iload[2]) ? iload[0] : iload[2];
      ibk = (iload[1] > iload[3]) ? iload[1] : iload[3];
      ifw = (ifw > nFw_min) ? ifw : nFw_min;
      ibk = (ibk > nBk_min) ? ibk : nBk_min;
      //
      Allocate_Dynamic_SlfEn(ifw, ibk, iload[5]);
      if (iload[5] != NTOT_DEN ) {cerr << "ERROR (Wide_self_energy_t::load_SelfEn_bin): mismatch for NTOT_DEN : "<< NTOT_DEN  << "(mod. sp.) <--> " << iload[5] << "(file)\n"; exit(100);}
      if (iload[7] != NTOT_DIM ) {cerr << "ERROR (Wide_self_energy_t::load_SelfEn_bin): mismatch for NTOT_DIM : "<< NTOT_DIM  << "(mod. sp.) <--> " << iload[7] << "(file)\n"; exit(100);}
      //
      if ( ( (0 != iload[13]) && (1 != iload[13]) ) ||
           ( (0 != iload[14]) && (1 != iload[14]) ) ||
           ( (0 >= iload[12]) && (0 != iload[13]) && (0 != iload[14]) ) ||
           ( (0  < iload[12]) && (1 != iload[13]) && (1 != iload[14]) ) ) {
        cout << "ERROR (Wide_self_energy_t::Load_Dynamic_SelfEn): cannot read self-energies with more than a lanczos block (input were "
             << iload[12]  << " , " << iload[13] << " , " << iload[13] << ") -- Will EXIT\n"; exit(100);
      }
      //
      N_PLS_fw = iload[0];
      N_PLS_bk = iload[1];
      load_SelfEn_bin_b(&file, Sig_dyn_fw, N_PLS_fw, N_PLS_FW_alloc, NTOT_DIM);
      if (0 < iload[13]) {load_SelfEn_bin_c(&file, &n1, i_piv_fw, NTOT_DIM); n_piv_fw = i_piv_fw[0];}
      load_SelfEn_bin_b(&file, Sig_dyn_bk, N_PLS_bk, N_PLS_BK_alloc, NTOT_DIM);
      if (0 < iload[14]) {load_SelfEn_bin_c(&file, &n2, i_piv_bk, NTOT_DIM); n_piv_bk = i_piv_bk[0];}
    }

    file.close();

  } else {
    file.close();
    //
    get_filename(filename1, Label, "Dyn", "-Fw");
    //sprintf(filename1, "%s/SelfEn_dyn_J2pi=%i%c_dt=%i-Fw.bin" , STOREDIR, MdSp->MSp_2j[I_SH], chip[MdSp->MSp_ip[I_SH]], MdSp->MSp_ch[I_SH]);
    cout << "\nOpening file '" << filename1 << "' for reading...\n" << flush;
    ifstream file1(filename1, ios::in|ios::binary);
    ierr2 = file1.bad();
    ierr  = load_SelfEn_bin_a(&file1,  iload);
    if (ierr) {cout << "ERROR (Wide_self_energy_t::load_SelfEn_bin, file='"<<filename1
                    <<"'): too many differencies, I'll skip reading the dynamic self-energy...\n\n";}
    //
    get_filename(filename2, Label, "Dyn", "-Bk");
    //sprintf(filename2, "%s/SelfEn_dyn_J2pi=%i%c_dt=%i-Fw.bin" , STOREDIR, MdSp->MSp_2j[I_SH], chip[MdSp->MSp_ip[I_SH]], MdSp->MSp_ch[I_SH]);
    cout << "\nOpening file '" << filename2 << "' for reading...\n" << flush;
    ifstream file2(filename2, ios::in|ios::binary);
    ierr3 = file2.bad();
    int i1    = load_SelfEn_bin_a(&file2,  iload+20);
    if ( i1 ) {cout << "ERROR (Wide_self_energy_t::load_SelfEn_bin, file='"<<filename2
                    <<"'): too many differencies, I'll skip reading the dynamic self-energy...\n\n";}
    
    ierr = ierr && i1;
    if ( ! ierr ) {
      ifw = (iload[ 0] > iload[ 2]) ? iload[ 0] : iload[ 2];
      ibk = (iload[21] > iload[23]) ? iload[21] : iload[23];
      ifw = (ifw > nFw_min) ? ifw : nFw_min;
      ibk = (ibk > nBk_min) ? ibk : nBk_min;
      //
      Allocate_Dynamic_SlfEn(ifw, ibk);
      N_PLS_fw = iload[ 0];
      N_PLS_bk = iload[21];
      //
      load_SelfEn_bin_b(&file1, Sig_dyn_fw, N_PLS_fw, N_PLS_FW_alloc, NTOT_DIM);
      if (0 < iload[13]) {load_SelfEn_bin_c(&file1, &n_piv_fw, i_piv_fw, NTOT_DIM); n_piv_fw = i_piv_fw[0];}
      load_SelfEn_bin_b(&file2, Sig_dyn_bk, N_PLS_bk, N_PLS_BK_alloc, NTOT_DIM);
      if (0 < iload[34]) {load_SelfEn_bin_c(&file2, &n_piv_bk, i_piv_bk, NTOT_DIM); n_piv_bk = i_piv_bk[0];}
      ierr = ierr1 && ierr2;
    }
    file1.close();
    file2.close();
  }
  //
  if (ierr) {
   if (ierr1)  cerr << "ERROR in opening file '" << filename  << "'!!\n";
   if (ierr2)  cerr << "ERROR in opening file '" << filename1 << "'!!\n";
   if (ierr3)  cerr << "ERROR in opening file '" << filename2 << "'!!\n";
   cerr << flush;
  }
  return ierr;}

//----------------------------------------------------------------------------
//
//  Functions for plotting the self-energy in ascii format
//
//----------------------------------------------------------------------------


static double x1;

void Wide_self_energy_t::plot_SelfEn(char* string1, char* string2, double *EFermi/*=NULL*/) {
  plot_SelfEn_stc(string1);
  plot_SelfEn_dyn(string2);
  return;}


void Wide_self_energy_t::plot_SelfEn_stc(char* string, double *EFermi/*=NULL*/) {
  ofstream file(string, ios::trunc|ios::out);

  file << "# Self-energy for channel Jf, pif, dT:"
       << "  " << MdSp->MSp_2j[I_SH] << "   " << MdSp->MSp_ip[I_SH] << "   " << MdSp->MSp_ch[I_SH] << endl
   //  << "#   Jf="<<Jf<<" ,  pif="<<pif<<" , dT="<<dT
       << "#\n#";
  if (NULL != EFermi) file << "   (E_Fermi =  "<< *EFermi <<" ,  from propagator)";
  file << "\n#\n#  number of orbits and of starting energies (if a Gmatrix is used):\n"
       << "#     " << NTOT_ORB<< "    " << N_GMTX<< "\n#\n";

  if (NULL != Sig_HF) {
    file << "# Static self-energy (EHF):\n";
    for(int nr=0;  nr<NTOT_ORB; ++nr) {
      file << MdSp->MSp_n[I_SH][nr] << MdSp->MSp_name[I_SH];
      for(int nc=0;  nc<NTOT_ORB; ++nc) {
        x1 = Sig_HF[nr*NTOT_ORB+nc];
        if (NULL != Sig_1b) x1 += Sig_1b[nr*NTOT_ORB+nc];
        file << "  " << setw(14) << setfill(' ') << setprecision(10) << x1;
      }
      file << endl;
    }
  }
  if (NULL != Sig_BHF) {
    file << "# Brueckner (w depenedent) EHF potential (NGMTX = "<<N_GMTX<<"):\n";
    for(int ne=0; ne<N_GMTX; ++ne) {
      file << "w= "<<Gmtx_ste[ne]<<endl;
      for(int nr=0;  nr<NTOT_ORB; ++nr) {
        file << MdSp->MSp_n[I_SH][nr] << MdSp->MSp_name[I_SH];
        for(int nc=0;  nc<NTOT_ORB; ++nc) {
          x1 = Sig_BHF[(ne*NTOT_ORB+nr)*NTOT_ORB+nc];
          if (NULL != Sig_1b) x1 += Sig_1b[nr*NTOT_ORB+nc];
          file << "  " << setw(14) << setfill(' ') << setprecision(10) << x1;
        }
        file << endl;
      }
    }
  }

  file.close();

  return;}

void Wide_self_energy_t::plot_SelfEn_dyn(char* string, double *EFermi/*=NULL*/) {

  ofstream file(string, ios::trunc|ios::out);

  double x_ph;

  file << "# Self-energy for channel Jf, pif, dT:"
       << "  " << MdSp->MSp_2j[I_SH] << "/2   " << MdSp->MSp_ip[I_SH] << "   " << MdSp->MSp_ch[I_SH] << endl;
  //
  if ((0 != Lambda) | (0 != parity) || (0 != charge)) 
    file << "#\n# Lambda, d_parity, d_charge :  " << Lambda << "   " << parity<< "   " << charge << endl;
  //
  file << "#\n#";
  if (NULL != EFermi) file << "   (E_Fermi =  "<< *EFermi <<" ,  from propagator)";
  file << "\n#\n#  number of fw- and bk- poles:\n"
       << "#     " << N_PLS_fw<< "    " << N_PLS_bk<< "   ("<< N_PLS_FW_alloc <<"  " << N_PLS_BK_alloc<< ", allocated)\n#\n";
/*
  if (N_EXT_PWV > 1) { 
    file << "#\n#  partial waves (from MdSp file):\n#";
    for(int nr=0;  nr<NTOT_ORB; ++nr)  file << "       "   << MdSp->Mor_sh[isp_msp[nr]];
//    for(int nr=0;  nr<NTOT_ORB; ++nr)  file << "       "   <<  MdSp->MSp_name[MdSp->Mor_sh[isp_msp[nr]]];
  }
  file << "#\n#  principal quantum numbers:\n#";
    for(int nr=0;  nr<NTOT_ORB; ++nr)  file << "       "   << MdSp->Mor_n[isp_msp[nr]];
     file << "\n#\n";
*/

    file << "#\n#  orbits :\n#";
    for(int nr=0;  nr<NTOT_ORB; ++nr)  file << "       "    << MdSp->Mor_n[isp_msp[nr]]
                   << MdSp->MSp_name[MdSp->Mor_sh[isp_msp[nr]]];
    file << "\n#\n";

  for(ifw=0; ifw<N_PLS_fw; ++ifw) {
    d1_ptr = Sig_dyn_fw + ifw*NTOT_DIM;
    file  << setw(14) << setfill(' ') << setprecision(10) <<  d1_ptr[0];
    file << setfill(' ') << "  " << setw(14) << setprecision(10) << d1_ptr[1] << "  ---  ";
    x_ph = 1.0; // cannot do this when comapring to TrOp!!! if (d1_ptr[NTOT_DEN] < 0) x_ph = -1.0;
    for(int i1=NTOT_DEN; i1<NTOT_DIM; ++i1)
      file << setfill(' ') << "  " << setw(16) << setprecision(12) << x_ph*d1_ptr[i1];
    file << endl;
  }

  file << "# ---\n";

  for(ibk=0; ibk<N_PLS_bk; ++ibk) {
    d1_ptr = Sig_dyn_bk + ibk*NTOT_DIM;
    file  << setw(14) << setfill(' ') << setprecision(10) <<  d1_ptr[0];
    file << setfill(' ') << "  " << setw(14) << setprecision(10) << d1_ptr[1] << "  ---  ";
    x_ph = 1.0; // cannot do this when comapring to TrOp!!!  if (d1_ptr[NTOT_DEN] < 0) x_ph = -1.0;
    for(int i1=NTOT_DEN; i1<NTOT_DIM; ++i1)
      file << setfill(' ') << "  " << setw(16) << setprecision(12) << x_ph*d1_ptr[i1];
    file << endl;
  }

  file.close();

  return;}
