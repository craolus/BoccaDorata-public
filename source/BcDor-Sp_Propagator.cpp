//
//   The class "SpProp_t" stores the one-body Green's function (or
//  sp propagator) of a finite many-particle sytem.
//
//   The Lehmann representation is used in this implementation, so sp
//  energies and (the components in the sp basis of the) overlap wave
//  funtcions are stored.
//   A sp model space (of type ModSpace_t)nmust be associated with each
//  object, to provide a discrete basis to expend the overla amplitudes.
//
//   Here, the case of a spherical system is implemented (J_tot == 0).
//
//
// C.Barbieri, GSI, December 2006.  (C.Barbieri@gsi.de).
//
//

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <iostream>
using namespace std;

#include "BcDor-Sp_classes.hh"
//#include "BcDor-Int_classes.hh"
//#include "BcDor-Utilities.hh"


#define SIZE_INP_LINE 300 // max line length when reading the sp.prop. from file 

#define DFLT_N1P     0  // default values for n1p, n2p, etc...
#define DFLT_N2P    -1
#define DFLT_D1P   0.0
#define DFLT_N1H     0
#define DFLT_N2H    -1
#define DFLT_D1H   0.0

#define SQR2 1.4142135623730950488016887




//
// Flags for reading/writing the sp. propagator:
static int irw_n1 = 0;
static int irw_n2 = 0;
static int irw_d1 = 0;

inline int irw_code(void )  {return ((irw_n1*2)+irw_n2)*2+irw_d1;} // uses the static variables
inline void irw_flags(int icode) {
  irw_d1 = icode%2; icode/=2; irw_n2 = icode%2; icode/=2; irw_n1 = icode%2;
  return;}
//
//inline int irw_code(int i1, int i2, int i3) {return ((i1*2)+i2)*2+i3;}
//inline void irw_flags(int *i1, int *i2, int *i3, int icode) {
//  *i3 = icode%2; icode/=2; *i2 = icode%2; icode/=2; *i1 = icode%2;
//  return;}


void SpProp_t::set_bounds_zero(void ) {
  //
  // NOTE: MdSp <- NULL should **NOT** be set here, otherwise it'll interfere
  //    with the situation in which the sp-prop is reinitialized but the 
  //    mod. space does not change...
  //     --> set it in 'init_ALL_values()' instead
  //
  n_qh = 0;
  n_qp = 0;
  return;}

void SpProp_t::set_alloc_size_zero(void ) {
  N_SPB_MX = -1;
  NTOT_QP = -1;
  NTOT_QH = -1;
  return;}


void SpProp_t::init_ALL_values(void ) {
  MdSp = NULL;
  init_NULL();
  set_alloc_size_zero();
  set_bounds_zero();
  return;
  }

void SpProp_t::init_NULL(void ) {
  // grouping by subshells
  Sp_clj_ep  = NULL;
  Sp_clj_eh  = NULL;
  Sp_clj_zp  = NULL;
  Sp_clj_zh  = NULL;
  Sp_clj_fqp = NULL;
  Sp_clj_fqh = NULL;
  Sp_clj_np  = NULL;
  Sp_clj_nh  = NULL;
  //
  Sp_clj_n1p = NULL;
  Sp_clj_n2p = NULL;
  Sp_clj_d1p = NULL;
  Sp_clj_n1h = NULL;
  Sp_clj_n2h = NULL;
  Sp_clj_d1h = NULL;
  //
  // quasiparticles       // quasiholes
  n_clj_p = NULL;         n_clj_h = NULL;
  ep  = NULL;             eh  = NULL;
  zp  = NULL;             zh  = NULL;
  fqp = NULL;             fqh = NULL;
  //
  n1p = NULL;             n1h = NULL;
  n2p = NULL;             n2h = NULL;
  d1p = NULL;             d1h = NULL;
  return;}

void SpProp_t::free_mem(void ) {
  //cout << "\nDeallocating single particle propagator... \n";
  //cout << "  max. number of s.p. orbits for a given partial wave: " << N_SPB_MX << endl;
  //cout << "  total deallocated size for quasiparticles:           " << n_qp << endl;
  //cout << "  total deallocated size forquasiholes:                " << n_qh << endl;
  cout << "SpProp_t:  freeing all the allocated memory" << endl;

  if (NULL != Sp_clj_ep  ) delete [] Sp_clj_ep;
  if (NULL != Sp_clj_eh  ) delete [] Sp_clj_eh;
  if (NULL != Sp_clj_zp  ) delete [] Sp_clj_zp;
  if (NULL != Sp_clj_zh  ) delete [] Sp_clj_zh;
  if (NULL != Sp_clj_fqp ) delete [] Sp_clj_fqp;
  if (NULL != Sp_clj_fqh ) delete [] Sp_clj_fqh;
  if (NULL != Sp_clj_np  ) delete [] Sp_clj_np;
  if (NULL != Sp_clj_nh  ) delete [] Sp_clj_nh;
  //
  if (NULL != Sp_clj_n1p  ) delete [] Sp_clj_n1p;
  if (NULL != Sp_clj_n2p  ) delete [] Sp_clj_n2p;
  if (NULL != Sp_clj_d1p  ) delete [] Sp_clj_d1p;
  if (NULL != Sp_clj_n1h  ) delete [] Sp_clj_n1h;
  if (NULL != Sp_clj_n2h  ) delete [] Sp_clj_n2h;
  if (NULL != Sp_clj_d1h  ) delete [] Sp_clj_d1h;

  // quasiparticles
  if (NULL != n_clj_p ) delete [] n_clj_p;
  if (NULL != ep      ) delete [] ep;
  if (NULL != zp      ) delete [] zp;
  if (NULL != fqp     ) delete [] fqp;
  //
  if (NULL != n1p     ) delete [] n1p;
  if (NULL != n2p     ) delete [] n2p;
  if (NULL != d1p     ) delete [] d1p;

  // quasiholes
  if (NULL != n_clj_h ) delete [] n_clj_h;
  if (NULL != eh      ) delete [] eh;
  if (NULL != zh      ) delete [] zh;
  if (NULL != fqh     ) delete [] fqh;
  //
  if (NULL != n1h     ) delete [] n1h;
  if (NULL != n2h     ) delete [] n2h;
  if (NULL != d1h     ) delete [] d1h;

  init_NULL();
  set_alloc_size_zero();
  set_bounds_zero();
  return;}


int SpProp_t::Allocate_mem(int nqp, int nqh) {
  free_mem();

  N_SPB_MX = 0;
  for(int i=0; i <MdSp->nsubsh; i++)
    if (N_SPB_MX < MdSp->MSp_no[i]) N_SPB_MX = MdSp->MSp_no[i];
  if (N_SPB_MX <= 0) {
    cerr << "\n An empty model space was assigned. Cannot allocate memory for the propagator!\n";
    set_bounds_zero();
    free_mem(); // it has actually been done above...  it doesn't hhurt though.
    return 1; // <-- to tell the calling routine that nothing was allocated...
    }

  //
  NTOT_QP = (1 > nqp) ? 1 : nqp;
  NTOT_QH = (1 > nqh) ? 1 : nqh;
  //
  cout << " Allocating memory..." << endl;
  //
  // grouping by subshells
  Sp_clj_ep  = new double*[MdSp->nsubsh]; // qp energies by subshell
  Sp_clj_eh  = new double*[MdSp->nsubsh];
  Sp_clj_zp  = new double*[MdSp->nsubsh]; // spec. factors by subshell
  Sp_clj_zh  = new double*[MdSp->nsubsh];
  Sp_clj_fqp = new double*[MdSp->nsubsh]; // overlap w.f.s by subshell
  Sp_clj_fqh = new double*[MdSp->nsubsh];
  Sp_clj_np  = new int[MdSp->nsubsh];     // # of overlap fnct.s by subshell
  Sp_clj_nh  = new int[MdSp->nsubsh];
  //
  Sp_clj_n1p  = new int*[MdSp->nsubsh];
  Sp_clj_n1h  = new int*[MdSp->nsubsh];
  Sp_clj_n2p  = new int*[MdSp->nsubsh];
  Sp_clj_n2h  = new int*[MdSp->nsubsh];
  Sp_clj_d1p  = new double*[MdSp->nsubsh];
  Sp_clj_d1h  = new double*[MdSp->nsubsh];
  //
  // quasiparticles
  n_clj_p = new int[NTOT_QP];       // partial wave the qp belongs to
  ep      = new double[NTOT_QP];
  zp      = new double[NTOT_QP];
  fqp     = new double[NTOT_QP*N_SPB_MX];
  //
  n1p     = new int[NTOT_QP];
  n2p     = new int[NTOT_QP];
  d1p     = new double[NTOT_QP];
  //
  // quasiholes
  n_clj_h = new int[NTOT_QH];       // partial wave the qh belongs to
  eh      = new double[NTOT_QH];
  zh      = new double[NTOT_QH];
  fqh     = new double[NTOT_QH*N_SPB_MX];
  //
  n1h     = new int[NTOT_QH];
  n2h     = new int[NTOT_QH];
  d1h     = new double[NTOT_QH];
  
  cout << "  max. number of s.p. orbits for a given partial wave: " << N_SPB_MX << endl;
  cout << "  total number of quasiparticles ALLOCATED:            " << NTOT_QP << endl;
  cout << "  total number of quasiholes ALLOCATED:                " << NTOT_QH << endl;
  cout << " --> arrays allocated." << endl << endl;

  Clean_up();
  return 0;}
  
void SpProp_t::Clean_up(void) {
  for(int i=0; i<MdSp->nsubsh; ++i) {
    Sp_clj_ep[i]  = NULL;    Sp_clj_eh[i]  = NULL;
    Sp_clj_zp[i]  = NULL;    Sp_clj_zh[i]  = NULL;
    Sp_clj_fqp[i] = NULL;    Sp_clj_fqh[i] = NULL;
    Sp_clj_np[i]  = 0;      Sp_clj_nh[i]  = 0;
    //
    Sp_clj_n1p[i] = NULL;    Sp_clj_n1h[i] = NULL;
    Sp_clj_n2p[i] = NULL;    Sp_clj_n2h[i] = NULL;
    Sp_clj_d1p[i] = NULL;    Sp_clj_d1h[i] = NULL;
    }
  for(int i=0; i<NTOT_QP; ++i) n_clj_p[i]=-1;
  for(int i=0; i<NTOT_QH; ++i) n_clj_h[i]=-1;
  n_qp = 0;
  n_qh = 0;
  return;}


// Error mesage if no mod. sp. information is given
SpProp_t::SpProp_t() {
  cerr << "\n ERROR: No information was given to build the s.p. propagator \n";
  init_ALL_values();
  return;
  }

// Deallocate the s.p. propgator to free memory space
SpProp_t::~SpProp_t() {
  cout << "\nDeallocating single particle propagator... \n";
  cout << "  max. number of s.p. orbits for a given partial wave: " << N_SPB_MX << endl;
  cout << "  total number of quasiparticles:                      " << n_qp << endl;
  cout << "  total number of quasiholes:                          " << n_qh << endl;

  free_mem();
  init_ALL_values(); // not really needed here...

  return;
  }


//
//  Create and empty propagator but allocate memory as indicated. If both
// iqp and iqh are negative, does not allocate anything.
SpProp_t::SpProp_t(int iqp, int iqh, ModSpace_t * inmds) {

  init_ALL_values();

  MdSp =inmds;

  if ((iqp >0) || (iqh > 0)) {
    iqp = (iqp>0) ? iqp : 0;
    iqh = (iqh>0) ? iqh : 0;
    cout << "  Allocating empty memory for a sp propagator..." << endl;
    if (Allocate_mem(iqp, iqh)) return;
    }

  return;
  }

static double xe;
static int ip,ih;
static int i, j, k;

// Builds a bare propagator, based on the mod. sp. given
SpProp_t::SpProp_t(ModSpace_t * inmds, int i_xe/*=0*/) {
                                   // defaults i_xe=0 (use 2n+l)hw for spe
                                   //     if   i_xe=1, get spe from Md. sp.

  init_ALL_values();

  MdSp = inmds;

  ih = 0;
  ip = 0;
//NSPBMXmsp = 0;
  for(i=0; i <MdSp->nsubsh; i++) {
//  if (NSPBMXmsp < MdSp->MSp_no[i]) NSPBMXmsp = MdSp->MSp_no[i];
    for (j=0; j < MdSp->MSp_no[i]; j++) 
      if (1 == MdSp->MSp_F[i][j]) ih++; else ip++;
    }
  
  cout << "\nGenerating the bare single particle propagtor (from the given model space)...\n";

  //
  // Allocate the relevant arrays:
  // -----------------------------
  //
  if (Allocate_mem(ip, ih)) return; // nothing was allocated

  
  ip=0, ih=0;
  for(i=0; i <MdSp->nsubsh; i++) {
   Sp_clj_n1p[i] = n1p + ip;
   Sp_clj_n2p[i] = n2p + ip;
   Sp_clj_d1p[i] = d1p + ip;
   Sp_clj_ep[i]  = ep  + ip;
   Sp_clj_zp[i]  = zp  + ip;
   Sp_clj_fqp[i] = fqp + ip*N_SPB_MX;
   Sp_clj_np[i]  = 0;
   //
   Sp_clj_n1h[i] = n1h + ih;
   Sp_clj_n2h[i] = n2h + ih;
   Sp_clj_d1h[i] = d1h + ih;
   Sp_clj_eh[i]  = eh  + ih;
   Sp_clj_zh[i]  = zh  + ih;
   Sp_clj_fqh[i] = fqh + ih*N_SPB_MX;
   Sp_clj_nh[i]  = 0;
   for (j=0; j<MdSp->MSp_no[i]; j++) {

     xe = (MdSp->MSp_n[i][j] * 2  + MdSp->MSp_l[i]) * MdSp->htom;
     if (i_xe) xe = MdSp->MSp_xe[i][j];
     
     if (1 == MdSp->MSp_F[i][j]) {
       n1h[ih] = DFLT_N1H; n1h[ih] = 1; // if (Sp_clj_nh[i] < 1) n1h[ih] = 1;
       n2h[ih] = DFLT_N2H;
       d1h[ih] = DFLT_D1H;
       eh[ih] = xe - 10.0;
       zh[ih] = 1.0;
       n_clj_h[ih] = i;
       for (int k=0; k<N_SPB_MX; ++k) {
         fqh[ih*N_SPB_MX+k] = 0.0;
         if ((k == MdSp->MSp_n[i][j]) && (k < MdSp->MSp_no[i] )) fqh[ih*N_SPB_MX+k] = 1.0;
       }
       ++ih;
       ++(Sp_clj_nh[i]);
     }
     else  {
       n1p[ip] = DFLT_N1P;  if (Sp_clj_np[i] < 2) n1p[ip] = 1;
       n2p[ip] = DFLT_N2P;
       d1p[ip] = DFLT_D1P;
       ep[ip] = xe;
       zp[ip] = 1.0;
       n_clj_p[ip] = i;
       for (int k=0; k<N_SPB_MX; ++k) {
         fqp[ip*N_SPB_MX+k] = 0.0;
         if ((k == MdSp->MSp_n[i][j]) && (k < MdSp->MSp_no[i] )) fqp[ip*N_SPB_MX+k] = 1.0;
       }
       ++ip;
       ++(Sp_clj_np[i]);
     }

    }
   }

  n_qp = ip;   n_qh = ih;

  cout << "  total number of quasiparticles:                     " << n_qp << endl;
  cout << "  total number of quasiholes:                         " << n_qh << endl;

  return;
  }

static int n1,n2,n3;
static double xz;
static char line[SIZE_INP_LINE];// ,orbit_name[45];
static char *line_ptr;
static int NSPBMXmsp,NSPBMXfile;

// Read the sp propagator from a file, the mod. sp. must also be given as input
SpProp_t::SpProp_t(char* fname, ModSpace_t * inmds, int *in_rw_flag/*=NULL*/) {
  cout << "\nReading the single particle propagtor form file '"
       << fname << "'..." <<endl;

  init_ALL_values();

  MdSp = inmds;

  n_qh = 0;
  n_qp = 0;
  
  NSPBMXmsp = 0;
  for(i=0; i <MdSp->nsubsh; i++)
    if (NSPBMXmsp < MdSp->MSp_no[i]) NSPBMXmsp = MdSp->MSp_no[i];


  FILE *infile;
  infile = fopen(fname, "rb");
  if (infile == NULL) {
    cerr << " Error in opening file '" << fname << "'\n";
    set_bounds_zero();
    return;
    }
     

  fgets(line,sizeof(line), infile);  //"#  Quasi- particle and hole...
  fgets(line,sizeof(line), infile);  //"#--------------------------...
   

  // read the total number of j and \pi subshells:
  fgets(line,sizeof(line), infile);
  fgets(line,sizeof(line), infile);  //"\n# number of (ilj\\pi) sub...
  fgets(line,sizeof(line), infile);
  if (sscanf(line,  "# %i %i %i\n", &n1,&NSPBMXfile,&n3) < 3) n3=0;
  if (n1 != MdSp->nsubsh) {
    cerr << "\n The number of shells given in file '" << fname 
         << "' does not seem to be the same as that of the model space:\n";
    cerr << " --> I'll skip reading the propagator!!\n\n";
    set_bounds_zero();
    return;
    }
  if (NSPBMXfile != NSPBMXmsp) {
    cout << "\n\n WARNING: It seems that--at least for some shells--there is a difference\n"
         <<         "  between the number of oribits in the propagator file and in the model space.\n"
         <<         "  If so, the s.p. propagator will be truncated to fit the model space.\n";
  }
  if (NULL != in_rw_flag) (*in_rw_flag) = n3;
  irw_flags(n3); // set input flags

  // read the total number of quasiparticles and quasiholes:
  fgets(line,sizeof(line), infile);
  fgets(line,sizeof(line), infile);  // "\n# Total numbers of qp and...
  fgets(line,sizeof(line), infile);
  sscanf(line,  "#%i%i", &n_qp,&n_qh);


  if ((n_qp <= 0) && (n_qh <= 0)) {
    cerr << "\n WARNING: it seems that an empty propagator has been assigned(\?\?)...\n";
    cerr << " --> I'll skip reading it!!\n\n";
    set_bounds_zero();
    return;
    }


    {// Head down count of the # of quasi-particles and holes:
    
    fpos_t pos; fgetpos(infile,&pos);// cout << " Position: " << pos << endl;
    ip=0, ih=0;
    for(i=0; i <MdSp->nsubsh; i++) {
      fgets(line,sizeof(line), infile);
      fgets(line,sizeof(line), infile);  // "\n# Subshell:\n"
      fgets(line,sizeof(line), infile);
      //  sscanf(line,  "# %s",orbitname); ...if one wants to check this too...
    
      fgets(line,sizeof(line), infile);
      sscanf(line,  "#%i",&n1);
      ip += n1;
      for (j=0; j<n1; j++) fgets(line,sizeof(line), infile);

      fgets(line,sizeof(line), infile);
      sscanf(line,  "#%i",&n1);
      ih += n1;
      for (j=0; j<n1; j++) fgets(line,sizeof(line), infile);
      }

    if ((ip != n_qp) || (ih != n_qh)) {
      cout << "\n\n WARNING: the number of quasiparticles and/or quasiholes detected";
      cout << "\n   is not the same as that given by the file '" << fname << "'.";
      cout << "\n   will reset n_qp and n_qh to these new values:\n";
      n_qp = ip;
      n_qh = ih;
      cout << "  max. number of s.p. orbits for a given partial wave: " << NSPBMXmsp << endl;
      cout << "  total number of quasiparticles:                      " << n_qp << endl;
      cout << "  total number of quasiholes:                          " << n_qh << endl;
      }

    // resume the position where the count has started:
    fsetpos(infile, &pos);
    //fgetpos(infile, &pos); cout << " Position re set to: " << pos << endl;
    }


  //
  // Allocate the relevant arrays:
  // -----------------------------
  //
  if (Allocate_mem(n_qp, n_qh)) return; // nothing was allocated


  ip=0, ih=0;
  for(i=0; i <MdSp->nsubsh; i++) {
    fgets(line,sizeof(line), infile);
    fgets(line,sizeof(line), infile);  // "\n# Subshell:\n"
    fgets(line,sizeof(line), infile);
    //  sscanf(line,  "# %s",orbitname);  ...if one wants to check this too...
    
    fgets(line,sizeof(line), infile);
    sscanf(line,  "#%i",&n1);
    Sp_clj_np[i] = n1;
    Sp_clj_ep[i] = ep+ip;
    Sp_clj_zp[i] = zp+ip;
    Sp_clj_fqp[i] = fqp + ip*N_SPB_MX;
    Sp_clj_n1p[i] = n1p+ip;
    Sp_clj_n2p[i] = n2p+ip;
    Sp_clj_d1p[i] = d1p+ip;
    for (j=0; j<Sp_clj_np[i]; j++) {
      fgets(line,sizeof(line), infile);
      line_ptr = line;
      //sscanf(line_ptr,  "%lf%lf%n",&xe,&xz,&n1);
      //line_ptr += n1;
      //ep[ip] = xe;
      //zp[ip] = xz / 100.;
      sscanf(line_ptr,  "%lf%n",&xe,&n1);
      ep[ip] = xe;
      line_ptr += n1;
      sscanf(line_ptr,  "%lf%n",&xz,&n1);
      zp[ip] = xz / 100.0;
      line_ptr += n1;
      if (irw_n1) {sscanf(line_ptr,  "%i%n",&n2,&n1);
                   n1p[ip] = n2;
                   line_ptr += n1;}
      if (irw_n2) {sscanf(line_ptr,  "%i%n",&n2,&n1);
                   n2p[ip] = n2;
                   line_ptr += n1;}
      if (irw_d1) {sscanf(line_ptr,  "%lf%n",&xz,&n1);
                   d1p[ip] = xz;
                   line_ptr += n1;}
      n_clj_p[ip] = i;
      for (int k=0; k<MdSp->MSp_no[i]; k++) {n2 = sscanf(line_ptr,"%le%n",&xz,&n1);
                                             if (n2 < 1) {xz=0.0; n1=0;}
                                             fqp[ip*N_SPB_MX+k] = xz;
                                             line_ptr += n1;
                                             }
      ip++;
      }

    fgets(line,sizeof(line), infile);
    sscanf(line,  "#%i",&n1);
    Sp_clj_nh[i] = n1;
    Sp_clj_eh[i] = eh+ih;
    Sp_clj_zh[i] = zh+ih;
    Sp_clj_fqh[i] = fqh + ih*N_SPB_MX;
    Sp_clj_n1h[i] = n1h+ih;
    Sp_clj_n2h[i] = n2h+ih;
    Sp_clj_d1h[i] = d1h+ih;
    for (j=0; j<Sp_clj_nh[i]; j++) {
      fgets(line,sizeof(line), infile);
      line_ptr = line;
      //sscanf(line_ptr,  "%lf%lf%n",&xe,&xz,&n1);
      //line_ptr += n1;
      //eh[ih] = xe;
      //zh[ih] = xz / 100.;
      sscanf(line_ptr,  "%lf%n",&xe,&n1);
      eh[ih] = xe;
      line_ptr += n1;
      sscanf(line_ptr,  "%lf%n",&xz,&n1);
      zh[ih] = xz / 100.0;
      line_ptr += n1;
      if (irw_n1) {sscanf(line_ptr,  "%i%n",&n2,&n1);
                   n1h[ih] = n2;
                   line_ptr += n1;}
      if (irw_n2) {sscanf(line_ptr,  "%i%n",&n2,&n1);
                   n2h[ih] = n2;
                   line_ptr += n1;}
      if (irw_d1) {sscanf(line_ptr,  "%lf%n",&xz,&n1);
                   d1h[ih] = xz;
                   line_ptr += n1;}
      n_clj_h[ih] = i;
      for (int k=0; k<MdSp->MSp_no[i]; k++) {n2 = sscanf(line_ptr,"%le%n",&xz,&n1);
                                             if (n2 < 1) {xz=0.0; n1=0;}
                                             fqh[ih*N_SPB_MX+k] = xz;
                                             line_ptr += n1;
                                             }
      ih++;
      }

    }


  n_qp = ip;  n_qh = ih;

  cout << "  total number of quasiparticles:                      " << n_qp << endl;
  cout << "  total number of quasiholes:                          " << n_qh << endl;


  if ((n_qp > NTOT_QP) || (n_qh > NTOT_QH)) {
    cout << "\n\n WARNING: the number of quasiparticles and/or quasiholes actually read";
    cout << "\n   is larger that the space allocated on memory. One of the following things";
    cout << "\n   may have happnend: given by the file '" << fname << "'.";
    cout << "\n\n  1) arrays have been allocated wtih the WRONG size! (maybe the first line in file `"
         << fname << "' is meesed up.)";    
    cout <<   "\n  2) you have decleard a wrong number od qp of qh for some subshell (still in  file `"
         << fname << "'.)";
    cout <<   "\n  3) soemthing else..." << endl;
    cout << "   total number of quasiparticles (read):    " << ip << endl;
    cout << "   total number of quasiholes (read):        " << ih << endl << endl;
    }

  if (fclose(infile) != 0) {
    cout << "Error in closing file '" << fname << "'...   " << endl;
    }

  return;
  }

// Write the sp propagator into a file
int SpProp_t::write(char* fname , int wcode) {  // wcode==0 (default)
  cout << "\nWriting the single particle propagtor to file '"
       << fname << "'..." <<endl; 
  cout << "  max. number of s.p. orbits for a given partial wave: " << N_SPB_MX << endl;
  cout << "  total number of quasiparticles:                     " << n_qp << endl;
  cout << "  total number of quasiholes:                         " << n_qh << endl;

  irw_flags(wcode); // set input flags
  wcode = irw_code(); // set input flags

  FILE *outfile;
  outfile = fopen(fname, "w");

  fprintf(outfile,"#  Quasi- particle and hole fragments of the sp propagator\n");
  fprintf(outfile,"#----------------------------------------------------------\n");
   

  // write the total number of j and \pi subshells:
  fprintf(outfile,"\n# number of (ilj\\pi) subshells, max n. of radial orbitals:\n");
  fprintf(outfile,  "# %7i %7i %7i\n", MdSp->nsubsh,N_SPB_MX,wcode);

  // write the total number of quasiparticles and quasiholes:
  fprintf(outfile,"\n# Total numbers of qp and qh stored here:\n");
  fprintf(outfile,  "# %7i %7i\n", n_qp,n_qh);

  for(i=0; i <MdSp->nsubsh; i++) {
    fprintf(outfile,"\n# Subshell:\n");
    fprintf(outfile,  "# %s \n",MdSp->MSp_name[i]);
    
    fprintf(outfile,  "# %7i  # -> tot n. of quasiparticles\n",Sp_clj_np[i]);
    for (j=0; j<Sp_clj_np[i]; j++) {
      fprintf(outfile,  " %15.6lf ",Sp_clj_ep[i][j]);
      fprintf(outfile,  " %12.3lf   ",100.*Sp_clj_zp[i][j]);
      if (irw_n1) fprintf(outfile,  " %4i ",Sp_clj_n1p[i][j]);
      if (irw_n2) fprintf(outfile,  " %4i ",Sp_clj_n2p[i][j]);
      if (irw_d1) fprintf(outfile,  " %7.2lf  ",Sp_clj_d1p[i][j]);
      for (int k=0; k<MdSp->MSp_no[i]; k++)
                          fprintf(outfile,  " %18.6le ",Sp_clj_fqp[i][j*N_SPB_MX+k]);
      fprintf(outfile,  "\n");
      }

    fprintf(outfile,  "# %7i  # -> tot n. of quasiholes\n",    Sp_clj_nh[i]);
    for (j=0; j<Sp_clj_nh[i]; j++) {
      fprintf(outfile,  " %15.6lf ",Sp_clj_eh[i][j]);
      fprintf(outfile,  " %12.3lf   ",100.*Sp_clj_zh[i][j]);
      if (irw_n1) fprintf(outfile,  " %4i ",Sp_clj_n1h[i][j]);
      if (irw_n2) fprintf(outfile,  " %4i ",Sp_clj_n2h[i][j]);
      if (irw_d1) fprintf(outfile,  " %7.2lf  ",Sp_clj_d1h[i][j]);
      for (int k=0; k<MdSp->MSp_no[i]; k++)
                          fprintf(outfile,  " %18.6le ",Sp_clj_fqh[i][j*N_SPB_MX+k]);
      fprintf(outfile,  "\n");
      }

    }


  if (fclose(outfile) != 0) {
    cout << "Error in closing file '" << fname << "'...   " << endl;
    }

  return 0;
  }


static double *ptr1,*ptr2;
static int nst, i_orb;
static double x1;

//
// Add nin quasiparticles with partial wave ish to the propagator
int SpProp_t::add_qp(int ish, int nin, int NDIM, double spectf[]) {
   // NOTE: for now this can be done ONLY if the subshell
   // is STILL EMPTY or if it is the last one that has been filled.
   // Also, NO allocation of further memory is done and the poles
   // in excess are discarded (it will return the # of lost poles).

  nin = (nin > 0) ? nin : 0;

  if (Sp_clj_np[ish] < 1) {
    j = n_qp;  // == starting point of the new subshell
    j = (j > 0) ? j : 0;
    Sp_clj_np[ish]  = 0;
    Sp_clj_ep[ish]  = ep+j;
    Sp_clj_zp[ish]  = zp+j;
    Sp_clj_n1p[ish] = n1p+j;
    Sp_clj_n2p[ish] = n2p+j;
    Sp_clj_d1p[ish] = d1p+j;
    Sp_clj_fqp[ish] = fqp+(j*N_SPB_MX);
  } else {
    if (n_qp < 1) {cerr << "\nAborted from SpProp_t::add_qp : nqp and Sp_clj_np are incinsistent\n"; exit(100);}
    if (ish != n_clj_p[n_qp-1]) return nin;
  }

  i = NTOT_QP - n_qp;
  nst = (nin < i) ? nin : i;// limit to the max allocated mem.available
  nst = (nst > 0) ? nst : 0;// make sure it is not negative since it is used
                            // below (should never happen though...)

  i_orb = MdSp->MSp_no[ish];
  i_orb = (i_orb < (NDIM-2)) ? i_orb : (NDIM-2);
  if ( (Sp_clj_ep[ish] - ep) + Sp_clj_np[ish] != n_qp) {
    cerr << "\nAborted from SpProp_t::add_qp : nqp, Sp_clj_np, Sp_clj_ep for the last subshell are incinsistent\n";
    exit(100);
  }
  j=n_qp; // starting point where to insert the new poles
  ptr1 = spectf;
  for(i=0;  i<nst; ++i) {
    n_clj_p[j] = ish;
    n1p[j] = DFLT_N1P;
    n2p[j] = DFLT_N2P;
    d1p[j] = DFLT_D1P;
    ep[j] = (*ptr1);
  //zp[j] = (*(ptr_1+1));  get it from imput data...
    ptr2 = ptr1+2;
    x1 = 0.0;
    for(k=0; k<i_orb; ++k) {
      fqp[j*N_SPB_MX+k] = *ptr2;
      x1 += (*ptr2)*(*ptr2);
      ++ptr2;
      }
    for(k=i_orb; k<N_SPB_MX; ++k) fqp[j*N_SPB_MX+k] = 0.0;
    zp[j] = x1;  // recomputes it...
    
    ptr1 += NDIM;
    ++j;
    }

  Sp_clj_np[ish] += nst; // increment by the # of new qparticles
  n_qp += nst;

  return (nin-nst); // ==0 if all sols added, >0 if (nin-nst) were discarded
  }

//
// Add nin quasiholes with partial wave ish to the propagator
int SpProp_t::add_qh(int ish, int nin, int NDIM, double spectf[]) {
   // NOTE: for now this can be done ONLY if the subshell
   // is STILL EMPTY or if it is the last one that has been filled.
   // Also, NO allocation of further memory is done and the poles
   // in excess are discarded (it will return the # of lost poles).

  nin = (nin > 0) ? nin : 0;

  if (Sp_clj_nh[ish] < 1) {
    j = n_qh;  // == starting point of the new subshell
    j = (j > 0) ? j : 0;
    Sp_clj_nh[ish]  = 0;
    Sp_clj_eh[ish]  =  eh+j;
    Sp_clj_zh[ish]  =  zh+j;
    Sp_clj_n1h[ish] = n1h+j;
    Sp_clj_n2h[ish] = n2h+j;
    Sp_clj_d1h[ish] = d1h+j;
    Sp_clj_fqh[ish] = fqh+(j*N_SPB_MX);
  } else {
    if (n_qh < 1) {cerr << "\nAbortet from SpProp_t::add_qh : nqh and Sp_clj_nh are incinsistent\n"; exit(100);}
    if (ish != n_clj_h[n_qh-1]) return nin;
  }

  i = NTOT_QH - n_qh;
  nst = (nin < i) ? nin : i;// limit to the max allocated mem.available
  nst = (nst > 0) ? nst : 0;// make sure it is not negative since it is used
                            // below (should never happen though...)

  i_orb = MdSp->MSp_no[ish];
  i_orb = (i_orb < (NDIM-2)) ? i_orb : (NDIM-2);
  if ( (Sp_clj_eh[ish] - eh) + Sp_clj_nh[ish] != n_qh) {
    cerr << "\nAborted from SpProp_t::add_qh : nqh, Sp_clj_nh, Sp_clj_eh for the last subshell are incinsistent\n";
    exit(100);
  }
  j=n_qh; // starting point where to insert the new poles
  ptr1 = spectf;
  for(i=0;  i<nst; ++i) {
    n_clj_h[j] = ish;
    n1h[j] = DFLT_N1H;
    n2h[j] = DFLT_N2H;
    d1h[j] = DFLT_D1H;
    eh[j] = (*ptr1);
  //zh[j] = (*(ptr_1+1));  get it from imput data...
    ptr2 = ptr1+2;
    x1 = 0.0;
    for(k=0; k<i_orb; ++k) {
      fqh[j*N_SPB_MX+k] = *ptr2;
      x1 += (*ptr2)*(*ptr2);
      ++ptr2;
      }
    for(k=i_orb; k<N_SPB_MX; ++k) fqh[j*N_SPB_MX+k] = 0.0;
    zh[j] = x1;  // recomputes it...
    
    ptr1 += NDIM;
    ++j;
    }

  Sp_clj_nh[ish] += nst; // increment by the # of new qholes
  n_qh += nst;

  return (nin-nst); // ==0 if all sols added, >0 if (nin-nst) were discarded
  }

int SpProp_t::Fill_with_charge(int dT_in) {
  //
  //  This is somewaht quick and dirty but it works, needs to be
  // coded better later on...
  //
  int ierr = 0;
  for(int i=0; i <MdSp->nsubsh; i++) {
	if (dT_in == MdSp->MSp_ch[i]) continue;

	int j2 = MdSp->MSp_2j[i];
	int ip = MdSp->MSp_ip[i];
	int ish_orig = MdSp->get_subsh(j2, ip, dT_in);
	if (0 > ish_orig) {
	  cout << " there's not channel with j,pi,ch="<<j2<<"/2,"<<ip<<","<<dT_in;
	  continue;
	}
	if (Sp_clj_nh[i] < 1) {
	  int j = Sp_clj_eh[ish_orig] - eh;
	  double *dptr = new double[(2+N_SPB_MX)*Sp_clj_nh[ish_orig]];
	  double *dptr2 = dptr;
	  for(int n=j; n<j+Sp_clj_nh[ish_orig]; ++n) {
		dptr2[0]=eh[n];
		dptr2[1]=zh[n];
		for (int l=0; l<N_SPB_MX; ++l) dptr2[l+2] = fqh[l + n*N_SPB_MX];
		dptr2 += (2+N_SPB_MX);
	  }
	  ierr += add_qh(i, Sp_clj_nh[ish_orig], 2+N_SPB_MX, dptr);
	  delete [] dptr;
	} else {
	  cout << "qh: ish_orig="<<ish_orig<<"    "<<Sp_clj_nh[ish_orig]<<endl;
	}

	if (Sp_clj_np[i] < 1) {
	  int j = Sp_clj_ep[ish_orig] - ep;
	  double *dptr = new double[(2+N_SPB_MX)*Sp_clj_np[ish_orig]];
	  double *dptr2 = dptr;
	  for(int n=j; n<j+Sp_clj_np[ish_orig]; ++n) {
		dptr2[0]=ep[n];
		dptr2[1]=zp[n];
		for (int l=0; l<N_SPB_MX; ++l) dptr2[l+2] = fqp[l + n*N_SPB_MX];
		dptr2 += (2+N_SPB_MX);
	  }
	  ierr += add_qp(i, Sp_clj_np[ish_orig], 2+N_SPB_MX, dptr);
	  delete [] dptr;
	} else {
	  cout << "qp: ish_orig="<<ish_orig<<"    "<<Sp_clj_np[ish_orig]<<endl;
	}

  }
  return 0;}

//-----------------------------------------------------------------
//
//  these routined to set the values of the arrays n1, n2 and d1
// are called as follows:
//   1) set_xx(v_p, v_h) : sets ALL qp to v_p and ALL qhs to v_h
//
//   2) set_xx(ip, ih, ptr_p, ptr_h) : sets the first ip qps to the values
//                                    conteined in the array ptr_p[0,..,ip-1].
//                                    Same for ih, ptr_h and the qhs.
//
//   3) set_xx(0, 0, NULL+1, NULL+1) : sets all to defaults
//
void SpProp_t::set_n1(int ip, int ih, int* vp, int* vh) {
  if (NULL == vp) {
    for(i=0; i<n_qp; ++i) n1p[i]=ip;
  } else {
    for(i=0; i<n_qp; ++i)
      if (i<ip) {n1p[i]=vp[i];} else {n1p[i]=DFLT_N1P;}
  }
  if (NULL == vh) {
    for(i=0; i<n_qh; ++i) n1h[i]=ih;
  } else {
    for(i=0; i<n_qh; ++i)
      if (i<ih) {n1h[i]=vh[i];} else {n1h[i]=DFLT_N1H;}
  }
  return;}

void SpProp_t::set_n2(int ip, int ih, int* vp, int* vh) {
  if (NULL == vp) {
    for(i=0; i<n_qp; ++i) n2p[i]=ip;
  } else {
    for(i=0; i<n_qp; ++i)
      if (i<ip) {n2p[i]=vp[i];} else {n2p[i]=DFLT_N2P;}
  }
  if (NULL == vh) {
    for(i=0; i<n_qh; ++i) n2h[i]=ih;
  } else {
    for(i=0; i<n_qh; ++i)
      if (i<ih) {n2h[i]=vh[i];} else {n2h[i]=DFLT_N2H;}
  }
  return;}

void SpProp_t::set_d1(double vp, double vh) {
    for(i=0; i<n_qp; ++i) d1p[i]=vp;
    for(i=0; i<n_qh; ++i) d1h[i]=vh;
  return;}

void SpProp_t::set_d1(int ip, int ih, double* vp, double* vh) {
    for(i=0; i<n_qp; ++i) if (i<ip) {d1p[i]=vp[i];} else {d1p[i]=DFLT_D1P;}
    for(i=0; i<n_qh; ++i) if (i<ih) {d1h[i]=vh[i];} else {d1h[i]=DFLT_D1H;}
  return;}
//
//----------------------- set of n1, n2, d1 ------------------------


int SpProp_t::cout_qplist(int wcode) {  // wcode==0 (default)

  irw_flags(wcode); // set input flags
  wcode = (irw_n1) || (irw_n2) || (irw_d1);

  cout << "\nQuasiparticles: \n";
  for(int i=0; i<n_qp; i++) {
    cout << i << "  " << MdSp->MSp_name[n_clj_p[i]] << "  " << ep[i] << "  " << zp[i];
    if (wcode) cout << "  (";
    if (irw_n1) cout << n1p[i] << " ";
    if (irw_n2) cout << n2p[i] << " ";
    if (irw_d1) cout << d1p[i];
    if (wcode) cout << ")";
    cout << "  -- ";
    for (int j=0; j<MdSp->MSp_no[n_clj_p[i]]; j++) cout << fqp[i*N_SPB_MX+j] << "  ";
    cout << endl;
    }
  cout << "\nQuasiholes: \n";
  for(int i=0; i<n_qh; i++) {
    cout << i << "  " << MdSp->MSp_name[n_clj_h[i]] << "  " << eh[i] << "  " << zh[i];
    if (wcode) cout << "  (";
    if (irw_n1) cout << n1h[i] << " ";
    if (irw_n2) cout << n2h[i] << " ";
    if (irw_d1) cout << d1h[i];
    if (wcode) cout << ")";
    cout << "  -- ";
    for (int j=0; j<MdSp->MSp_no[n_clj_h[i]]; j++) cout << fqh[i*N_SPB_MX+j] << "  ";
    cout     << endl;
    }

  return 0;
  }


int SpProp_t::cout_byshell(int wcode) {  // wcode==0 (default)

  irw_flags(wcode); // set input flags
  wcode = (irw_n1) || (irw_n2) || (irw_d1);

  for(int i=0; i <MdSp->nsubsh; i++) {
    cout << endl << "#   " << MdSp->MSp_name[i] << "   # subshell" << endl;
    cout         << "#   " << Sp_clj_np[i] << "   # tot n. of quasiparticles" << endl;
    for (int j=0; j<Sp_clj_np[i]; j++) {
      cout << j << "  " << Sp_clj_ep[i][j] << "  " << Sp_clj_zp[i][j];
      if (wcode) cout << "  (";
      if (irw_n1) cout << Sp_clj_n1p[i][j] << " ";
      if (irw_n2) cout << Sp_clj_n2p[i][j] << " ";
      if (irw_d1) cout << Sp_clj_d1p[i][j];
      if (wcode) cout << ")";
      cout << "  -- ";
      for (int k=0; k<MdSp->MSp_no[i]; k++) cout << Sp_clj_fqp[i][j*N_SPB_MX+k] << "  ";
      cout     << endl;
      }
    cout         << "#   " << Sp_clj_nh[i] << "   # tot n. of quasiholes"     << endl;
    for (int j=0; j<Sp_clj_nh[i]; j++) {
      cout << j << "  " << Sp_clj_eh[i][j] << "  " << Sp_clj_zh[i][j];
      if (wcode) cout << "  (";
      if (irw_n1) cout << Sp_clj_n1h[i][j] << " ";
      if (irw_n2) cout << Sp_clj_n2h[i][j] << " ";
      if (irw_d1) cout << Sp_clj_d1h[i][j];
      if (wcode) cout << ")";
      cout << "  -- ";
      for (int k=0; k<MdSp->MSp_no[i]; k++) cout << Sp_clj_fqh[i][j*N_SPB_MX+k] << "  ";
      cout     << endl;
      }  
    }

  return 0;
  }



static int n_pp, n_hh, n_ph, n_2p1h, n_2h1p;

void SpProp_t::count_qp_confs(int *n_pp_mx, int *n_hh_mx, int *n_ph_mx,
                              int *n_2p1h_mx, int *n_2h1p_mx,
                              int iplot/*=1*/,
                              int t_2qp/*=-1*/, int t_3qp/*=-1*/) {

  *n_ph_mx   = 0;
  *n_pp_mx   = 0;
  *n_hh_mx   = 0;
  *n_2p1h_mx = 0;
  *n_2h1p_mx = 0;
  int n_chs = 0;

  int n1,n2,n3,nst,nJ,nt;

  double x1,x2;
  double ph_Mb, pp_Mb, hh_Mb;   //, 2p1h_Mb, 2h1p_Mb;
  ph_Mb = pp_Mb = hh_Mb = 0.0;  //  2p1h_Mb = 2h1p_Mb = 0.0;

  MdSp->Calculate_Jch_bounds();


  //
  // coupling of pp and hh fragments:
  //
  pp_Mb = 0.0; hh_Mb = 0.0;
  n_chs=0;
  for(int J=0;                J<=MdSp->Jmax;      J++) {
  for(int ip=0;               ip<2;              ip++) {
  for(int ch=MdSp->ch_pp_mn; ch<=MdSp->ch_pp_mx; ch++) {

    n_pp = 0; n_hh = 0;
    for(int ish1=0 ; ish1 <MdSp->nsubsh; ish1++) {
    for(int ish2=ish1; ish2 <MdSp->nsubsh; ish2++) {

       nJ = ((MdSp->MSp_2j[ish1]+MdSp->MSp_2j[ish2])/2 + J + 1)%2;

       if ( MdSp->MSp_ch[ish1]+MdSp->MSp_ch[ish2]       != ch) continue;
       if ((MdSp->MSp_ip[ish1]+MdSp->MSp_ip[ish2]+ip)%2 != 0 ) continue;
//       if (AM.triang(MdSp->MSp_2j[ish1],MdSp->MSp_2j[ish2],2*J) <= 0) continue;
       if ((MdSp->MSp_2j[ish1]+MdSp->MSp_2j[ish2] < 2*J) ||
           (MdSp->MSp_2j[ish2]+    2*J          < MdSp->MSp_2j[ish1]) ||
           (          2*J   +MdSp->MSp_2j[ish1] < MdSp->MSp_2j[ish2]) ) continue;


       // pp frags.
       for(n1=0;   n1<Sp_clj_np[ish1]; n1++) {nst=0; if (ish1==ish2) nst=n1+nJ;
       for(n2=nst; n2<Sp_clj_np[ish2]; n2++) {
         nt = Sp_clj_n1p[ish1][n1]+Sp_clj_n1p[ish2][n2];
         //if ((t_2qp < 0) || (nt==t_2qp)) ++n_pp;
         if ((t_2qp < 0) || (nt>=t_2qp)) ++n_pp;
       } }

       // hh frags.
       for(n1=0;   n1<Sp_clj_nh[ish1]; n1++) {nst=0; if (ish1==ish2) nst=n1+nJ;
       for(n2=nst; n2<Sp_clj_nh[ish2]; n2++) {
         nt = Sp_clj_n1h[ish1][n1]+Sp_clj_n1h[ish2][n2];
         //if ((t_2qp < 0) || (nt==t_2qp)) ++n_hh;
         if ((t_2qp < 0) || (nt>=t_2qp)) ++n_hh;
       } }

     } } // ish1, ish2
     if ( n_pp + n_hh > 0) n_chs++;
     *n_pp_mx = (*n_pp_mx > n_pp) ? *n_pp_mx : n_pp;
     *n_hh_mx = (*n_hh_mx > n_hh) ? *n_hh_mx : n_hh;

     x1 = pow((n_pp/1.e3),2)*8.0;
     x2 = pow((n_hh/1.e3),2)*8.0;
     pp_Mb += x1;
     hh_Mb += x2;
  
  if (iplot) cout << "J=" << J << "  ip=" << ip << "  ch=(+/-)" << ch 
                  << " --> n_pp= " << n_pp << " / n_hh=" << n_hh
                  << "    tot Mb:   " << x1 << " and "<<  x2 << endl;
  } } } // J,ip,ch
  if (iplot) cout << "\n number of pp/hh channels: " << n_chs << endl;
  if (iplot) cout << "\n required memory (in Mb): " << 2.0*pp_Mb << "(pp) ,   " << 2.0*hh_Mb << "(hh)" << endl << endl;

  //
  // coupling of ph fragments:
  //
  ph_Mb = 0.0;
  n_chs=0;
  for(int J=0;                J<=MdSp->Jmax;      J++) {
  for(int ip=0;               ip<2;              ip++) {
  for(int ch=MdSp->ch_ph_mn; ch<=MdSp->ch_ph_mx; ch++) {

    n_ph = 0;
    for(int ish1=0; ish1 <MdSp->nsubsh; ish1++) {
    for(int ish3=0; ish3 <MdSp->nsubsh; ish3++) {

       if ( MdSp->MSp_ch[ish1]-MdSp->MSp_ch[ish3]       != ch) continue;
       if ((MdSp->MSp_ip[ish1]+MdSp->MSp_ip[ish3]+ip)%2 != 0 ) continue;
//       if (AM.triang(MdSp->MSp_2j[ish1],MdSp->MSp_2j[ish2],2*J) <= 0) continue;

       if ((MdSp->MSp_2j[ish1]+MdSp->MSp_2j[ish3] < 2*J) ||
           (MdSp->MSp_2j[ish3]+    2*J          < MdSp->MSp_2j[ish1]) ||
           (          2*J   +MdSp->MSp_2j[ish1] < MdSp->MSp_2j[ish3]) ) continue;


       // ph frags.
       for(n1=0; n1<Sp_clj_np[ish1]; n1++) {
       for(n3=0; n3<Sp_clj_nh[ish3]; n3++) {
         nt = Sp_clj_n1p[ish1][n1]+Sp_clj_n1h[ish3][n3];
         //if ((t_2qp < 0) || (nt==t_2qp)) ++n_ph;
         if ((t_2qp < 0) || (nt>=t_2qp)) ++n_ph;
       } }

     } } // ish1, ish2
     if ( n_ph > 0) n_chs++;
     *n_ph_mx = (*n_ph_mx > n_ph) ? *n_ph_mx : n_ph;

     x1 = pow((n_ph/1.e3),2)*8.0;
     ph_Mb += x1;
  
  if (iplot) cout << "J=" << J << "  ip=" << ip << "  ch=" << ch << " --> n_ph=" << n_ph
                  << "    tot Mb:   " << x1 << endl;
  } } } // J,ip,ch
  if (iplot) cout << "\n number of ph/hp channels: " << n_chs << endl;
  if (iplot) cout << "\n required memory (in Mb): " << 2.0*ph_Mb << "(ph)" << endl << endl;




  n_chs = 0;
  for(int ish=0; ish<MdSp->nsubsh; ish++) {
    int J2 = MdSp->MSp_2j[ish];
    int ch = MdSp->MSp_ch[ish];
    int ip = MdSp->MSp_ip[ish];

    x1 = 0.0;    x2 = 0.0;
    n_2p1h = 0;  n_2h1p = 0;

    for(int Jab=0; Jab<=MdSp->Jmax; Jab++) {
        
      for(int ish1=0 ;   ish1 <MdSp->nsubsh; ish1++) {
      for(int ish2=ish1; ish2 <MdSp->nsubsh; ish2++) {

       if ((MdSp->MSp_2j[ish1]+MdSp->MSp_2j[ish2] < 2*Jab) ||
           (MdSp->MSp_2j[ish2]+    2*Jab          < MdSp->MSp_2j[ish1]) ||
           (          2*Jab   +MdSp->MSp_2j[ish1] < MdSp->MSp_2j[ish2]) ) continue;

       nJ = ((MdSp->MSp_2j[ish1]+MdSp->MSp_2j[ish2])/2 + Jab + 1)%2;

       for(int ish3=0; ish3 <MdSp->nsubsh; ish3++) {

        if  (MdSp->MSp_ch[ish1]+MdSp->MSp_ch[ish2]-MdSp->MSp_ch[ish3]       != ch) continue;
        if ((MdSp->MSp_ip[ish1]+MdSp->MSp_ip[ish2]+MdSp->MSp_ip[ish3]+ip)%2 != 0)  continue;
        if ((     2*Jab         + MdSp->MSp_2j[ish3] < J2) ||
            (MdSp->MSp_2j[ish3] +    J2              < 2*Jab) ||
            (     J2            +    2*Jab       < MdSp->MSp_2j[ish3]) ) continue;


//            nij++;

          // pph frags.
         for(n1=0;   n1<Sp_clj_np[ish1]; n1++) {nst=0;if (ish1==ish2) nst=n1+nJ;
         for(n2=nst; n2<Sp_clj_np[ish2]; n2++) {
         for(n3=0;   n3<Sp_clj_nh[ish3]; n3++) {
           nt = Sp_clj_n1p[ish1][n1]+Sp_clj_n1p[ish2][n2]+Sp_clj_n1h[ish3][n3];
           //if ((t_3qp < 0) || (nt==t_3qp)) ++n_2p1h;
           if ((t_3qp < 0) || (nt>=t_3qp)) ++n_2p1h;
         } } }

          // hhp frags.
         for(n1=0;   n1<Sp_clj_nh[ish1]; n1++) {nst=0;if (ish1==ish2) nst=n1+nJ;
         for(n2=nst; n2<Sp_clj_nh[ish2]; n2++) {
         for(n3=0;   n3<Sp_clj_np[ish3]; n3++) {
           nt = Sp_clj_n1h[ish1][n1]+Sp_clj_n1h[ish2][n2]+Sp_clj_n1p[ish3][n3];
           //if ((t_3qp < 0) || (nt==t_3qp)) ++n_2h1p;
           if ((t_3qp < 0) || (nt>=t_3qp)) ++n_2h1p;
         } } }

        }   // ish3
       } }  // ish1, ish2
    } // Jab
    if (n_2p1h + n_2h1p > 0) n_chs++;
    *n_2p1h_mx = (*n_2p1h_mx > n_2p1h) ? *n_2p1h_mx : n_2p1h;
    *n_2h1p_mx = (*n_2h1p_mx > n_2h1p) ? *n_2h1p_mx : n_2h1p;

     x1 = pow((n_2p1h/1.e3),2)*8.0;
     x2 = pow((n_2h1p/1.e3),2)*8.0;
     //Mb_2p1h += x1;
     //Mb_2h1p += x2;
  
  if (iplot) cout << "J=" << J2 << "/2   ip=" << ip << "  ch=" << ch
                  << " --> n_pph=" << n_2p1h << "  ,  n_hhp=" << n_2h1p
                  << "    tot Mb:   " << x1 << " and "<<  x2<< endl;

  } // ish

  if (iplot) cout << "\n number of pph/hhp channels: " << n_chs<< endl << endl;


  return;
  }



double SpProp_t::quick_n_dirty_diff(SpProp_t* gin) {
  double diff = 0.0;
  double x1,x2;
  int imax;

  // qps
  imax= (n_qp > gin->n_qp) ? n_qp : gin->n_qp;
  for(i=0; i<imax; ++i) {
   if (i<     n_qp) {x1=     ep[i];} else {x1=0.0;}
   if (i<gin->n_qp) {x2=gin->ep[i];} else {x2=0.0;}
   diff += fabs(x1-x2);
   //if (fabs(x1-x2) > .01) cout << "qp: " << i << "  " << x1 << "  " << x2 << endl;
   }

  // qhs
  imax= (n_qh > gin->n_qh) ? n_qh : gin->n_qh;
  for(i=0; i<imax; ++i) {
   if (i<     n_qh) {x1=     eh[i];} else {x1=0.0;}
   if (i<gin->n_qh) {x2=gin->eh[i];} else {x2=0.0;}
   diff += fabs(x1-x2);
   //if (fabs(x1-x2) > .01) cout << "qh: " << i << "  " << x1 << "  " << x2 << endl;
   }

  return diff;
  }



//
// Utilities sections....  better to move it to a new file soon or later..
//
//

// 
//  Make a pivots file to be used to be Used with the Dyson-BAGEL algorythm
//
void SpProp_t::DysPivots(char* fname, int n_bk, int n_fw, double z_cut) {
  cout << "\nGenerating a file named '"
       << fname << "' that included BAGEL pivots absed in this propagator..." << endl; 

  int itest;
  double x1,x2;
  FILE *outfile;
  outfile = fopen(fname, "w");

  fprintf(outfile,"#   Template pivots file\n");
  fprintf(outfile,"#\n");
  fprintf(outfile,"# number of (ilj\\pi) subshells, max n. of radial orbitals:\n");
  fprintf(outfile,"# %7i %7i\n\n", MdSp->nsubsh,N_SPB_MX);


  for(i=0; i <MdSp->nsubsh; i++) {
    fprintf(outfile,"\n\nshell,npvts:   %i   %i  # %s \n", i ,MdSp->MSp_no[i] , MdSp->MSp_name[i]);

    x2 = -1.0e+10;
    itest = 1; if (Sp_clj_nh[i] < 1) itest = -100;
    while(itest >= 0) {
      itest = -100;
      x1 = +1.0e+10;
      for (j=0; j<Sp_clj_nh[i]; j++) 
        if ((Sp_clj_eh[i][j] < x1) && (Sp_clj_eh[i][j] > x2) && (100.*Sp_clj_zh[i][j] >= z_cut)) {
          x1 = Sp_clj_eh[i][j];
          itest = j;}
      if (itest >= 0) {
          x2 = x1;
          fprintf(outfile,  " %3i  %3i  ", n_fw, n_bk );
          for (int k=0; k<MdSp->MSp_no[i]; k++)
                          fprintf(outfile,  " %18.6le ",Sp_clj_fqh[i][itest*N_SPB_MX+k]);
      fprintf(outfile,  "   #    e- = %13.6lf\n",Sp_clj_eh[i][itest]);
      }
    }

    x2 = -1.0e+10;
    itest = 1; if (Sp_clj_np[i] < 1) itest = -100;
    while(itest >= 0) {
      itest = -100;
      x1 = +1.0e+10;
      for (j=0; j<Sp_clj_np[i]; j++) 
        if ((Sp_clj_ep[i][j] < x1) && (Sp_clj_ep[i][j] > x2) && (100.*Sp_clj_zp[i][j] >= z_cut)) {
          x1 = Sp_clj_ep[i][j];
          itest = j;}
      if (itest >= 0) {
          x2 = x1;
          fprintf(outfile,  " %3i  %3i  ", n_fw, n_bk );
      for (int k=0; k<MdSp->MSp_no[i]; k++)
                          fprintf(outfile,  " %18.6le ",Sp_clj_fqp[i][itest*N_SPB_MX+k]);
      fprintf(outfile,  "   #    e+ = %13.6lf\n",Sp_clj_ep[i][itest]);
      }
    }
    
  } // end loop over subshells


  if (fclose(outfile) != 0) {
    cout << "Error in closing file '" << fname << "'...   " << endl;
  } else {
    cout << "  ...pivots done." << endl;
  }

  return;
  }

