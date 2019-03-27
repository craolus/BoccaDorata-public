// 
//
//  The class "ModSpace_t" contains the single particle orbits that form
// the basis for a spherical nucleus in p-n formalism
//
//
//
//   The model space can be read from or stored into a file.
//
//   The input format is as follow:
//   ------------------------------
//
//	/*   #
//	     #
//	     #
//	     #
//	     n1
//
//	     #  			\
//	     l    2j	ip		 |
//	     nmx  npmx  nhmx		 |
//	     #  			 |
//	     0    shF(0)    xe(0)	 }   1st shell 
//	     1    shF(1)    xe(1)	 |
//	     :  			 |
//	     :  			 |
//	     n    shF(n)    xe(n)	 |
//	     :  			 |
//	     nmx  shF(nmx)  xe(nmx)	 |
//					/
//	     #  			\
//	     l    2j	ip               |
//	     nmx  npmx  nhmx		 |
//	     #  			 |
//	     0    shF(0)    xe(0)	 }   2nd shell
//	     1    shF(1)    xe(1)	 |
//	     :  			 |
//	     :  			 |
//	     nmx  shF(nmx)  xe(nmx)	 |
//					/
//		     :
//		     :
//		     :
//	     #  			\
//	     l    2j	ip      	 |
//	     nmx  npmx  nhmx		 |
//	     #  			 |
//	     0    shF(0)    xe(0)	 }   n1-th shell
//	     1    shF(1)    xe(1)	 |
//	     :  			 |
//	     :  			 |
//	     nmx  shF(nmx)  xe(nmx)	 |
//					/              */
//
//
//   The meanings of the parameters are as follow:
//   ---------------------------------------------
//
//             n1 :  total number of s.p. shells with different l, j, ip
//                  and ch (charge) quantum numbers.
//
//             l  , 2j  , ip , ch  : orbital angular momentum, total angular
//                               momentum (MULTIPLIED by 2!), parity and
//                               charge of 
//                               of the subshell (ip = 0 or 1 , 0==even and
//                               1==odd, ch=-1,0,1,..).  These
//                               are all integer numbers. 
//                              
//             nmx  :  total number of s.p. shells with given l, j and ip. These
//                    shells differ only by the value of the principal q.# 'n'
//
//             npmx :  total number of PARTICLE fragments to be kept for this
//                    shell when iterating.
//
//             npmx :  total number of HOLE fragments to be kept for this
//                    shell when iterating.
//
//	       n    shF(n)    xe(n) :
//                       'n' is the principal q.# of the s.p. state.  The
//                          complete set of q.#s is n,l,j,{m_j,m_t}
//
//                       shF  tells if the state is to be considered above
//                          (shF=1) or below (shF=0) the fermi level. (Actually
//                          this info is not used in practice).
//
//                       xe   Hartee-Fock (or MF) energy of the single particle
//                          state. It could be used as to pass the eigenvalue of
//                          a MF potential U (as was done with the 5-lvs code).
//                          In the O16 case, this info is not used in practice.
//
//
//  Here is a (wrong) example of the input file:
//  ------------------------------------
//
//
//	 /*  # input model space quantum numbers
//	     #-----------------------------------
//
//	     # number of j and \pi  subshells:
//	     8
//
//	     # subshell:  s1/2
//	     0 1 0
//	     2 1 2
//	      Fermi level and sp energies
//	     0 0 -3.7	    #  0s1/2
//	     1 1 -4.8	    #  1s1/2
//
//	     # subshell:  p3/2
//	     1 3 1
//	     2 1 3
//	      Fermi level and sp energies
// 	     0 0 -3.8	    #  0p3/2
// 	     1 1 -7.78      #  1p3/2
//		     :
//		     :
//		     :                             */
//
//
//

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <iostream>
using namespace std;

#define STRL 9

#include "BcDor-PhysConsts.hh"

#include "BcDor-Sp_classes.hh"

//
//  (Re)set the intial values of some important variables. To be called 
// by the constructors.
void ModSpace_t::Reset_init_values(void ) {
  Init_NULL_shells();
  Init_NULL_orbits();
  nsubsh = -1;
  n_spb_tot = -1;
  //
  Jmax=-100;
  ch_pp_mn = 100 , ch_pp_mx = -100 ;
  ch_hh_mn = 100 , ch_hh_mx = -100 ;
  ch_ph_mn = 100 , ch_ph_mx = -100 ;
  ch_pph_mn= 100 , ch_pph_mx= -100 ;
  ch_hhp_mn= 100 , ch_hhp_mx= -100 ;
  //
  Reset_phys_const();
  return;
  }


void ModSpace_t::Reset_phys_const(void ) {
  htom = 0.0;
  bho  = 0.0;
  bsq  = 0.0;
  mass = NUCLEONmass;
  return;
  }

void ModSpace_t::Set_phys_const_bHO(double xb) {
  bho  = xb;
  bsq  = bho*bho;
  //mass = NUCLEONmass;
  htom = MeVfm/(bsq * mass / MeVfm);
  cout << "\n Oscillator frequency set to  b = " << bho
       << " fm ,     htom = " << htom << " MeV\n\n";
  return;
  }


void ModSpace_t::Set_phys_const_hwHO(double hw) {
  htom = hw;
  //mass = NUCLEONmass;
  bsq  = MeVfm/(htom * mass / MeVfm);
  bho  = sqrt(bsq);
  cout << "\n Oscillator frequency set to  b = " << bho
       << " fm ,     htom = " << htom << " MeV\n\n";
  return;
  }


// Initialize and empty model space
ModSpace_t::ModSpace_t() {
  Reset_init_values();
  cout << "\nAn empty model space was allocated... what for?\n" << endl;
  return;
  }


//  Destructor: free the allocated memory
ModSpace_t::~ModSpace_t() {
  cout << "\nDeallocating single particle basis..." << endl;
  Clear_all();
  return;}


void ModSpace_t::Init_NULL_shells(void ) {
  MSp_l      = NULL;
  MSp_2j     = NULL;
  MSp_ip     = NULL;
  MSp_ch     = NULL;
  MSp_chall  = NULL;
  MSp_name   = NULL;
  MSp_no     = NULL;
  //
  MSp_n      = NULL;
  MSp_F      = NULL;
  MSp_xe     = NULL;
  //
  MSp_nqp    = NULL;
  MSp_nqh    = NULL;
  return;}

void ModSpace_t::Init_NULL_orbits(void ) {
  Mor_sh     = NULL;
  Mor_n      = NULL;
  Mor_F      = NULL;
  Mor_xe     = NULL;
  return;}

void ModSpace_t::Allocate_shells(int n_of_subshells ) {
  nsubsh = n_of_subshells;
  //
  // List of all the partial waves (and charge states)
  MSp_l  = new int[nsubsh];
  MSp_2j = new int[nsubsh];
  MSp_ip = new int[nsubsh];
  MSp_ch = new int[nsubsh];
  MSp_chall = new char[nsubsh*STRL];
  MSp_name  = new char*[nsubsh];
  MSp_no = new int[nsubsh];
  for(int i=0; i<nsubsh; ++i) MSp_name[i] = MSp_chall+(STRL*i);
  //
  // pointers to pointers for the loops over other q.#s
  MSp_n  = new int*[nsubsh];
  MSp_F  = new int*[nsubsh];
  MSp_xe = new double*[nsubsh];
  //
  MSp_nqp  = new int[nsubsh];  // tot # of quasipart. (to be kept iterating?)
  MSp_nqh  = new int[nsubsh];  // tot # of quasiholes (to be kept iterating?)
  //
  cout << " --> total number of allocated subshells:     " << nsubsh << endl;
  return;}

void ModSpace_t::Allocate_orbits(int n_of_orbits ) {
  n_spb_tot = n_of_orbits;
  //
  // Complete list of the orbitals present int the mod space: 
  Mor_sh = new int[n_spb_tot];       // partial wave it belongs to
  Mor_n  = new int[n_spb_tot];       // principal q.# 'n'
  Mor_F  = new int[n_spb_tot];       // can be used to store an integer q.#
  Mor_xe = new double[n_spb_tot];    // can be used to store a double q.#
  //
  //  These are not strictly necessary but help later
  // when constructing pp and ph babses:
  Mor_l   = new int[n_spb_tot];       // orbital ang. momentum q.# 'l'
  Mor_2j  = new int[n_spb_tot];       // twice the tot. angular momentum q.# '2j'
  Mor_ip  = new int[n_spb_tot];       // parity: 0<-->+ , 1<-->-
  Mor_ch  = new int[n_spb_tot];       // charge  (usually, -1<-->e,  0<-->n, 1<-->p)
  //
  cout << " --> total number of allocated orbits:        " << n_spb_tot << endl;
  return;}


void ModSpace_t::Clear_all() {
  Free_mem_shells();
  Free_mem_orbits();
  return;}

void ModSpace_t::Free_mem_shells(void ) {
  cout << "  number of deallocated subshells:    " << nsubsh << endl;
  if (NULL != MSp_l     ) delete [] MSp_l;
  if (NULL != MSp_2j    ) delete [] MSp_2j;
  if (NULL != MSp_ip    ) delete [] MSp_ip;
  if (NULL != MSp_ch    ) delete [] MSp_ch;
  if (NULL != MSp_chall ) delete [] MSp_chall;
  if (NULL != MSp_name  ) delete [] MSp_name;
  if (NULL != MSp_no    ) delete [] MSp_no;

  if (NULL != MSp_n     ) delete [] MSp_n;
  if (NULL != MSp_F     ) delete [] MSp_F;
  if (NULL != MSp_xe    ) delete [] MSp_xe;

  if (NULL != MSp_nqp   ) delete [] MSp_nqp;
  if (NULL != MSp_nqh   ) delete [] MSp_nqh;
  Init_NULL_shells();
  nsubsh = -1;
  return;}

void ModSpace_t::Free_mem_orbits(void ) {
  cout << "  number of deallocated s.p. orbits:  " << n_spb_tot << endl;
  if (NULL != Mor_sh ) delete [] Mor_sh;
  if (NULL != Mor_n  ) delete [] Mor_n;
  if (NULL != Mor_F  ) delete [] Mor_F;
  if (NULL != Mor_xe ) delete [] Mor_xe;
  //
  if (NULL != Mor_l  ) delete [] Mor_l ;
  if (NULL != Mor_2j ) delete [] Mor_2j;
  if (NULL != Mor_ip ) delete [] Mor_ip;
  if (NULL != Mor_ch ) delete [] Mor_ch;
  Init_NULL_orbits();
  n_spb_tot = -1;
  return;}

ModSpace_t::ModSpace_t(int norb, int nemx, int nspc) {
  //
  //  norb = total number of partial waves to be generated (starts with the
  //                                                    lowest values of lj).
  //
  //  nemx = max value of  rho==2*n+l
  //
  //  nspc = number of fermionic species:
  //                             if nspc==2,    protons and neutrons in equal #
  //                             if nspc==1 (or anything else),       electrons
  //


  Reset_init_values();

  int nlj,l,nmx,n1,n2,nch;

  //
  // dimension of the arrays to be allocated by 'head down' counting...
  //
  nsubsh    = 0;  //total number of subshells
  n_spb_tot = 0;  // total number of s.p. orbitals (in all subshells)
  for (int i=0; i<norb; i++) {
     nlj = i;
     if ((2 == nspc) && (i >= norb/2)) nlj=i-(norb/2);
     l = (nlj+1)/2;
     if (nemx-l < 0) continue;  // we don't want orbitals beyond nemx: if
                               // we don't skip here the following satements
                               // will give troubles...!!!
          nmx = 1 + (nemx-l)/2; 
          n_spb_tot += nmx;
          nsubsh++;   //the first half is neutron the second is proton
     }
  cout << "\nGenerating single particle basis (from truncation constraints)..." << endl;

  Allocate_shells(nsubsh);
  Allocate_orbits(n_spb_tot);


  n1 = 0;
  n2 = 0;
  for (int i=0; i<norb; i++) {
     nlj = i; nch = -1;  // note that the default for nspc is set here
     switch(nspc) {
        case 1:  nch = -1; nlj=i; break;   // electrons
        case 2:  nch = 0;  nlj=i;                            // neutrons
                 if (i >= norb/2) {nlj=i-(norb/2); nch=1;}   // protons
                 break;
        }

     l = (nlj+1)/2;
     if (nemx-l < 0) continue;

        MSp_l[n1]  = l;
        MSp_2j[n1] = 2*(nlj-l)+1;
        MSp_ip[n1] = l%2;
        MSp_ch[n1] = nch;
        MSp_no[n1] = 1 + (nemx-l)/2;
        MSp_nqh[n1] = MSp_no[n1]/2;
        MSp_nqp[n1] = MSp_no[n1]-MSp_nqh[n1];
         {MSp_name[n1] = MSp_chall+(STRL*n1);
          int n=0;
          if (nch == -1) MSp_chall[STRL*n1+n]='e';
          if (nch ==  0) MSp_chall[STRL*n1+n]='v';
          if (nch ==  1) MSp_chall[STRL*n1+n]='p';
          n++;
          MSp_chall[STRL*n1+n]='_';
          n++;
          if (l == 0) MSp_chall[STRL*n1+n]='s';
          if (l == 1) MSp_chall[STRL*n1+n]='p';
          if (l == 2) MSp_chall[STRL*n1+n]='d';
          if (l >= 3) MSp_chall[STRL*n1+n]='f'+(l - 3);
          n++;
          int jj = MSp_2j[n1];
          if (jj > 9) {MSp_chall[STRL*n1+n]='0'+(jj/10);
                       jj %= 10;
                       n++; }
          MSp_chall[STRL*n1+n]='0'+jj;
          MSp_chall[STRL*n1+n+1]='/';
          MSp_chall[STRL*n1+n+2]='2';
          MSp_chall[STRL*n1+n+3]='\0';
          }
          
          
          MSp_n[n1]  = Mor_n+n2;
          MSp_F[n1]  = Mor_F+n2;
          MSp_xe[n1] = Mor_xe+n2;
          for (int ii=0; ii<MSp_no[n1]; ii++) {
            Mor_sh[n2] = n1;
            Mor_n[n2]  = ii;
            Mor_F[n2]  = 0;
            Mor_xe[n2] = 0.e0;
			
			Mor_l[n2]   = MSp_l[n1];       
			Mor_2j[n2]  = MSp_2j[n1];   
			Mor_ip[n2]  = MSp_ip[n1];      
			Mor_ch[n2]  = MSp_ch[n1];       
			
            n2++;
            }

        n1++;
     }
  }


ModSpace_t::ModSpace_t(char* fname) {

  Reset_init_values();

  cout << "\nReading single particle basis from file '"
       << fname << "'..." <<endl;
  nsubsh = 0;
  n_spb_tot = 0;

  char line[132],orbit_name[45];
  FILE *infile;
  int n1,n2,n3,n4,k;
  double x1;


  infile = fopen(fname, "rb");
  if (infile == NULL) {
    cerr << " Error in opening file '" << fname << "'";
    return;
    }
     

  fgets(line,sizeof(line), infile);  //"# Definition of the model space by ...
  fgets(line,sizeof(line), infile);  //"#----------------------------------...
   
  // read the total number of j and \pi subshells:
  fgets(line,sizeof(line), infile); 
  fgets(line,sizeof(line), infile);  //"# number of (ilj\\pi) subshells and...

  fgets(line,sizeof(line), infile);
  sscanf(line," %i%i\n", &nsubsh,&n_spb_tot);

  Allocate_shells(nsubsh);
  Allocate_orbits(n_spb_tot);
  

  k = 0;  // counts the  of orbitals, at the end it must be == n_spb_tot

  for(int i=0; i<nsubsh; i++) {

    fgets(line,sizeof(line), infile);
    fgets(line,sizeof(line), infile);
    sscanf(line,"# subshell:  %8s ",orbit_name);
    MSp_name[i] = MSp_chall+(STRL*i);
    for(n1=0; n1<STRL; n1++) MSp_chall[STRL*i+n1]=orbit_name[n1];

    // read quantum numbers of the i-th subshell
    fgets(line,sizeof(line), infile);
    sscanf(line,"%i%i%i%i",  &n1,&n2,&n3,&n4);
    MSp_l[i]  = n1;
    MSp_2j[i] = n2;
    MSp_ip[i] = n3;
    MSp_ch[i] = n4;
    // read the number of model space shells, p fragments and h fragments
    fgets(line,sizeof(line), infile);
    sscanf(line,"%i%i%i", &n1,&n2,&n3); 
    MSp_no[i]  = n1;
    MSp_nqp[i] = n2;
    MSp_nqh[i] = n3;


    // for each sp state, give position w.r.t. Fermi level and sp energy
    fgets(line,sizeof(line), infile);  // "# principal q.#,  p/h character and...

    MSp_n[i] = Mor_n + k;
    MSp_F[i] = Mor_F + k;
    MSp_xe[i] = Mor_xe + k;
    for(int j=0; j<MSp_no[i]; j++) {
      fgets(line,sizeof(line), infile);
      sscanf(line,"%i%i%lf ",&n1,&n2,&x1);
      Mor_sh[k] = i;
      Mor_n[k]  = n1;
      Mor_F[k]  = n2;
      Mor_xe[k] = x1;
	  
	  Mor_l[k]   = MSp_l[i];       
	  Mor_2j[k]  = MSp_2j[i];   
	  Mor_ip[k]  = MSp_ip[i];      
	  Mor_ch[k]  = MSp_ch[i];       
	  
      k++;
    }

  }
  fclose(infile);

  cout << " --> the model space has been read successfully!" << endl << endl;
  
  return;
  }


int ModSpace_t::write(char* fname) {
  cout << "\nWriting single particle basis to file '"
       << fname << "'..." <<endl;
  cout << " total number of subshells:     " << nsubsh << endl;
  cout << " total number of s.p. orbitals: " << n_spb_tot << endl;

  FILE *outfile;
  outfile = fopen(fname, "w");


  fprintf(outfile,"# Definition of the model space by listing its s.p. orbitals\n");
  fprintf(outfile,"#------------------------------------------------------------\n\n");
   

  // total number of j and \pi subshells:
  fprintf(outfile,"# number of (ilj\\pi) subshells and total number of orbitals:\n");
  fprintf(outfile," %7i %7i\n", nsubsh,n_spb_tot);

  for(int i=0; i<nsubsh; i++) {

    fprintf(outfile,"\n# subshell:  %s \n",MSp_name[i]);
    // quantum numbers of the i-th subshell
    fprintf(outfile," %7i %7i %7i %7i\n",  MSp_l[i],MSp_2j[i],MSp_ip[i],MSp_ch[i]);
    // number of model space shells, p fragments and h fragments
    fprintf(outfile," %7i %7i %7i\n",  MSp_no[i],MSp_nqp[i],MSp_nqh[i]);

    // for each sp state, give position w.r.t. Fermi level and sp energy
    fprintf(outfile,"# principal q.#,  p/h character and sp energies:\n");
    for(int n3=0; n3<MSp_no[i]; n3++) {
      fprintf(outfile," %7i %7i %8.3f     # %i%s\n",
            MSp_n[i][n3],MSp_F[i][n3],MSp_xe[i][n3],MSp_n[i][n3],MSp_name[i]);
    }

  }
  
  if (fclose(outfile) != 0) {
    cout << "Error in closing file '" << fname << "'...   " << endl;
    }

  return 0;
  }


int ModSpace_t::cout_msp(void ) {
  cout << "\n\nSingle particle basis:\n";
  for(int i=0; i<nsubsh; i++) {
    cout << " " << MSp_name[i] 
         << " :  l="    << MSp_l[i]  << " ,  j="  << MSp_2j[i]
         << "/2 ,  ip=" << MSp_ip[i] << " ,  ch=" << MSp_ch[i];
    
    cout << "  (no,nqp,nqh=" << MSp_no[i]  << ", "
                             << MSp_nqp[i] << ", "
                             << MSp_nqh[i] << ")\n";
    }
  cout << endl;

  return 0;
  }


int ModSpace_t::get_orbit(int nor, int j2or, int ipor, int chor) {
  int ish;
  for(int i=0; i<n_spb_tot; i++) {
    ish = Mor_sh[i];
    if (  (nor  == Mor_n[i])    && (ipor == MSp_ip[ish])
       && (j2or == MSp_2j[ish]) && (chor == MSp_ch[ish]) ) return i;
    }

  return -100;
  }


int ModSpace_t::get_subsh(int j2sh, int ipsh, int chsh) {
  for(int ish=0; ish<nsubsh; ish++) {
    if ( (ipsh == MSp_ip[ish])
       && (j2sh == MSp_2j[ish]) && (chsh == MSp_ch[ish]) ) return ish;
    }

  return -100;
  }


void ModSpace_t::Calculate_Jch_bounds(int iplot) { // default, iplot==0
  int n1,n2,n3;

  // get bounds on tot. ang. mom.:
  Jmax=0;
  for(int i=0; i <nsubsh; i++) if (Jmax < MSp_2j[i]) Jmax=MSp_2j[i];
  
  // bounds on tot. charge:
  ch_pp_mn = 2*(MSp_ch[0]), ch_pp_mx = 2*(MSp_ch[0]);
  ch_hh_mn =-2*(MSp_ch[0]), ch_hh_mx =-2*(MSp_ch[0]);
  ch_ph_mn = 0                  , ch_ph_mx = 0;
  ch_pph_mn=    MSp_ch[0] , ch_pph_mx=    MSp_ch[0] ;
  ch_hhp_mn=   -MSp_ch[0] , ch_hhp_mx=   -MSp_ch[0] ;
  //
  for(int ish1=0 ; ish1 <nsubsh; ish1++) {n1=MSp_ch[ish1];
  for(int ish2=ish1; ish2 <nsubsh; ish2++) {n2=MSp_ch[ish2];
   if (ch_pp_mn >  n1+n2) ch_pp_mn= n1+n2;
   if (ch_pp_mx <  n1+n2) ch_pp_mx= n1+n2;
   if (ch_hh_mn > -n1-n2) ch_hh_mn=-n1-n2;
   if (ch_hh_mx < -n1-n2) ch_hh_mx=-n1-n2;
   for(int ish3=0; ish3 <nsubsh; ish3++) {n3=MSp_ch[ish3];
     if (ch_ph_mn >  n1-n3) ch_ph_mn= n1-n3;
     if (ch_ph_mx <  n1-n3) ch_ph_mx= n1-n3;
     if (ch_pph_mn > n1+n2-n3) ch_pph_mn = n1+n2-n3;
     if (ch_pph_mx < n1+n2-n3) ch_pph_mx = n1+n2-n3;
     if (ch_hhp_mn > n3-n2-n1) ch_hhp_mn = n3-n2-n1;
     if (ch_hhp_mx < n3-n2-n1) ch_hhp_mx = n3-n2-n1;
     }
   }
   }
  
  if (iplot) {
    cout << "\n\nCalculated bounds on total ang. mom. and charge:\n";
    cout << " Jmax= " << Jmax << endl;
    cout << " ch_pp_mn =" << ch_pp_mn  << "  ,  ch_pp_mx = " << ch_pp_mx  << endl;
    cout << " ch_hh_mn =" << ch_hh_mn  << "  ,  ch_hh_mx = " << ch_hh_mx  << endl;
    cout << " ch_ph_mn =" << ch_ph_mn  << "  ,  ch_ph_mx = " << ch_ph_mx  << endl;
    cout << " ch_pph_mn=" << ch_pph_mn << "  ,  ch_pph_mx= " << ch_pph_mx << endl;
    cout << " ch_hhp_mn=" << ch_hhp_mn << "  ,  ch_hhp_mx= " << ch_hhp_mx
         << endl << endl;
  }

  return;}


//void ModSpace_t::count_2b_confs(int *n_pp, int *n_hh, int *n_ph, int *n_2p1h, int *n_2h1p) {
void ModSpace_t::count_2b_confs(int *n_2b, int *n_2bme, int iplot/*=0*/) {

/*
  *n_ph   = 0;
  *n_pp   = 0;
  *n_hh   = 0;
  *n_2p1h = 0;
  *n_2h1p = 0;
*/
  int n_chs = 0;

  int nst,nJ,i2b;


  Calculate_Jch_bounds();

  //
  // coupling of two 1-body orbitals:
  //
  *n_2b   = 0;
  *n_2bme = 0;
  n_chs=0;
  for(int J=0;         J<=Jmax;       ++J) {
  for(int ip=0;        ip<2;         ++ip) {
  for(int ch=ch_pp_mn; ch<=ch_pp_mx; ++ch) {

    i2b=0;
    for(int ia=0 ; ia < nsubsh; ++ia) {
    for(int ib=ia; ib < nsubsh; ++ib) {

       nJ = abs((MSp_2j[ia]+MSp_2j[ib])/2 - J + 1)%2;

       if ( MSp_ch[ia]+MSp_ch[ib]       != ch) continue;
       if ((MSp_ip[ia]+MSp_ip[ib]+ip)%2 != 0 ) continue;
//       if (AM.triang(MSp_2j[ia],MSp_2j[ib],2*J) <= 0) continue;
       if ((MSp_2j[ia]+MSp_2j[ib] < 2*J) ||
           (MSp_2j[ib]+    2*J    < MSp_2j[ia]) ||
           (    2*J   +MSp_2j[ia] < MSp_2j[ib]) ) continue;

       // pp frags.
       for(int nash=0;   nash<MSp_no[ia]; ++nash) {nst=0; if (ia==ib) nst=nash+nJ;
       for(int nbsh=nst; nbsh<MSp_no[ib]; ++nbsh) {
         ++i2b;
       } }

     } } // ia , ib
     if ( i2b > 0) n_chs++;

  *n_2b   += i2b;
  *n_2bme += (i2b*(i2b+1))/2;
  
  if (iplot > 1) cout << "J=" << J << "  ip=" << ip << "  ch=(+/-)" << ch << " --> i2b= " << i2b << " / i2me=" << (i2b*(i2b+1))/2 << endl;
  } } } // J,ip,ch

  if (iplot) {
    cout << "\n total number of 2-body configuations: " << *n_2b   << endl;
    cout <<   " total number of 2-body mateix els.  : " << *n_2bme << endl << endl;
    cout <<   " tot. number of non empty 2-body channels: " << n_chs << endl << endl;
    }

  return;}
