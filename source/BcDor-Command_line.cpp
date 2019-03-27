#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
using namespace std;

#include "BcDor-Run_vars.hh"


static int n1,n2;
static double x1,x2;

void verbose(void ) {
  cout << "\n BoccaDorata package for many-body Green's functions (version "
       << VERSION << ", compiled on " << __DATE__ << " at " << __TIME__ << ").\n\n";
  return;}

void print_usage(void ) {
  cout << "\n\n USAGE:\n";
  cout << "\n  -h                          This help \n";
  cout << "\n  -v                          Version \n\n";
  cout << "\n  WkFold[=bcdwk]              Set the working folder\n";
  cout << "\n  A=...                       Mass # (if not specified, the program will ask for it)";
  cout << "\n  N=...                       Neutron # (if not specified, the program will ask for it)";
  cout << "\n  Z=...                       Proton # (if not specified, the program will ask for it)\n";
  cout << "\n  bHO=...      hwHO=...       H.o. length or energy (if none is specified, the program will ask for bHO)\n";
  cout << "\n  SpProp[=sp_prop]            Input s.p. propagator file\n";
//  cout << "\n  OutProp=...                 Passes a string that will be included in the name of each output s.p. propagator files\n";
  cout << "\n  MdSp[=input_msp]            Model space file\n";
  cout << "\n  Vpp[=Vpp.bcd]               Interaction file\n";
  cout << "\n  SelInt=i                    When more than one interaction is loaded, it selects the interaction number 'i' (must be i>=0)\n";
//  cout << "\n  ck                         set 'iCheck==1' (check for doubly def. mtx els of Vpp --normally this is not needed)\n";
  cout << "\n  Vc[=Vpp_Coul.bcd]           Add Coulomb from file (it will divide it by the h.o. length)  \n";
  cout << "\n  AddVpp=C,Vpp_file.bcd       Add the interaction from the given file to all alredy stored 2-body force. The matrix"
       << "\n                             elements will be multiplied by the factor 'C'.  In case of a G-matrix (more than a starting"
       << "\n                             energy) or when several sets of interaction are present, 'Vpp_file.bcd' will be added"
       << "\n                             to all of them.\n";
  cout << "\n  Trel[=3]                    Options for adding the kinetic energy in the **harmonic oscillator** basis:"
       << "\n                               0  no correction for c.m. motion"
       << "\n                               1  ti*(A-1)/A + 2b(pi*pj)  kinetic energy  (the 2b part must be already included in the interaction file)"
       << "\n                               2  2b((pi-pj)^2) form of   the kin. energy (the 2b part must be already included in the interaction file)"
       << "\n                               3  ti*(A-1)/A + 2b(pi*pj)  kinetic energy     (read the 2b part from the file passesd using 'Trel_pipj')"
       << "\n                               4  2b((pi-pj)^2)  form of  the kinetic energy (read the 2b part from the file passesd using 'Trel_2body')"
//       << "\n                               5  ti*(A-1)/A + 2b(pi*pj) + a_Tcm*T_cm  kinetic energy (read the 2b part from a file, as for Trel==3)"
//       << "\n                               6  2b((pi-pj)^2) + a_Tcm*T_cm  form of the kinetic energy (read the 2b part from a file, as for Trel==4')"
       << "\n                               7  kinetic term is SUPPRESSED\n";
  cout << "\n  Trel_pipj [=Trelpp_1.bcd]   Matrix elements of p_i*p_j\n";
  cout << "\n  Trel_2body[=Trelpp_2.bcd]   Matrix elements ot the kin. energy in 2b form\n";
  cout << "\n  ExtU1[=i[,U[,file]]]        External 1-body potential:"
	   << "\n  ExtU1[=3,U1,U2]              i=0  no 1-b potential [default]"
	   << "\n                                 1  diagonal form, taken from spe in the model space file"
	   << "\n                                 2  the potential from the specified file and multiply it by U [Defaults are U=1.0 file='Usp_ext.bcd']"
	   << "\n                                 3  shitfs charge ==0 by U0 and charge ==1 by U1";
//	   << "\n                                 4  other options are not available yet...\n";
  //  cout << "\n  a_Tcm[=0.0]                 Fraction of T_cm to be left in (with Trel = 5, 6 only)\n";
  cout << "\n  Radii                       Calculates the matter radii bsed on given sp propagator and Rrmspp_file (defaults: 1+2-body"
  	   << "\n                             and 'Rrmspp_1.bcd').\n";
  cout << "\n  Rrms[=i,[Rrmspp_file]]      Options for calculating the center-of-mas corrected r.m.s. matter raidus and input file with"
	     << "\n                             required two-body matrix elements:"
	     << "\n                              i = 1    1- and 2-body  form with 2-body   mtx els.  (r1*r2)   , default file is 'Rrmspp_1.bcd'"
	     << "\n                                = 2     2-body-only   form with reative  mtx els.  (r1-r2)^2 , default file is 'Rrmspp_2.bcd'"
	     << "\n                             If only the Rrms option is given r.m.s. radii are calculated at each iteration."
  	   << "\n                             Defaults:  i==1   and  Rrmspp_file == 'Rrmspp_1.bcd'.\n";

  cout << "\n  mass_p[=938.918...]          Set the particle mass to be used for calculations:\n";
  cout <<   "  mass_n[=938.918...]            - mass_p,  mass_n,  and mass_e  are the proton, neutron and electron masses \n";
  cout <<   "  mass_e[=0.511...]                     For all of these, if the input value is missing then it is set to\n";
  cout <<   "                                        the BcDor defaults values found in 'BcDot-PhysConsts.hh'.\n";
  cout <<   "  mass_ave[=938.918,0.511]       - mass_ave is used for all other charges then p, n, and e  and  also for setting\n";
  cout <<   "                                        the h.o. parameters in the model space.  NOTE that THIS IS THE ONE USED FOR\n";
  cout <<   "                                        SETTING THE MODEL SPACE AND THE KINETIC ENERGY (up to April 2015).\n";
  cout <<   "                                 - Setting either mass_p and mass_n, it also changes the mass_ave ot the agerage.\n";
  cout <<   "                                        of the two.\n";
  cout <<   "                                 - Setting either mass_e, it also changes  mass_ave the same as mass_e.\n";

  cout << "\n  ItrMax[=0]                  Max # of iterations\n";
  cout << "\n  nDfw[=-100]                 Approximate # of particle solutions to look for when solving the Dyson eq. (<0 means all!)\n";
  cout << "\n  nDbk[=-100]                 Approximate # of   hole   solutions to look for when solving the Dyson eq. (<0 means all!)\n";
  cout << "\n  HF=[i]                      Hartree-Fock calculation, it sets ItrMax=i [default==0]\n";
  cout << "\n  2nd[=i]                     Second order self-energy, it sets ItrMax=i [default==0]\n";
  cout << "\n  ExtSE[=i]                   Reads the self-energy form a file in forlder 'bcdwk', it sets ItrMax=i [default==0]\n";
  cout << "\n  CCD[=aLM[,nStp]]            Solves the CCD equations ASSUMING AN HF PROPAGATOR AS INPUT. It useslinear mixing with mixing\n"
       <<   "                             parameter aLM [default value is aLM == 0.5] and it adiabatically increases the contributions beyond\n"
       <<   "                             1st order to the CC amplitudes in nStp steps (default nStp == 1). This last feature sometimes helps\n"
       <<   "                             reaching convergence but it slows down iterations; setting nStp==1 means to switch this off.\n";
  cout << "\n  SelCharge[=i_ch]            When calculating or manipulating a sefl-energy ('MakeSelfEn', 'PlotSelfEn' and 'PlotSelfEn'), only the\n"
       <<   "                             the channels with charge == to 'i_ch' are done.\n"
       <<   "                             If solving the Dyson equation, only channles with charge == 'i_ch' are calculated and then copied into\n"
       <<   "                             the other channels (assuming cherge independence).\n";
  cout << "\n  MakeSelfEn[Fw,Bk][=i]       Builds the dinamic self-energy for the subshell i (and for all of them if i<0) [default: i=-1]\n"
       << "                               Adding 'Fw' computes only the forward part (~2p1h), 'Bk' only the backward part (~2h1p), otherwise\n"
       << "                               both are done.  The outputs are stored in folder 'bcdwk'. \n";
  cout << "\n  PlotSelfEn[=i]              Write the full self-energy for the subshell i (and for all of them if i<0)to text files. The quantum"
       << "\n                             numbers Jf, pi, dt and EHF or Fw/Bw description will be added to the file name to distiguish between"
       << "\n                             different partial waves.\n";
//  cout << "\n  PlotSelfEn[=i]             plot the full self-energy for the subshell i (and for all of them if i<0). Use 'OutProp=...' to indicate"
//       << "\n                             the file name and/or directory (Jf,pi,dt and EHF or Fw/Bw description will be added).\n";
  cout << "\n  Lanczos                     Uses a Krylov subspace algorhitm (BAGEL-like) to solve the Dyson Eq.\n";
  cout << "\n  LanDysPiv[=DysPivots.bcd]   File with the pivots to be used by the dyson-BAGEL algorhytm (if the exists, reads the pivots from there)\n";
  cout << "\n  SetEf=i,Ef                  The Fermi level for the subshell n. 'i' will be forced to 'Ef'\n";
  cout << "\n  Koltun                      Extract information from the sp propagator (# of particles, GMK sum rule etc..)\n";
  cout << "\n  MBPT2                       Second order many-body PT\n";
  cout << "\n  MakeMdSp                    Builds a model space file. The output will be in 'input_msp-wt'\n";
  cout << "\n  MakeSpProp[=i]              Builds a template sp_prop file based in the given model space (output is 'sp_prop-wt')\n"
       <<   "                             i = 4*irw_n1 + 2*irw_n2 + irw_d1 [defaults i=0] controls the flags to write the output.\n";
  cout << "\n  ConvSpProp[=i]              Rewrite the input sp_prop file with the flags specified by i (output is 'sp_prop-wt')\n"
       <<   "                             i = 4*irw_n1 + 2*irw_n2 + irw_d1 [defaults i=0].\n";
  cout << "\n  SpPropStat                  Prints out information from the input propagator (such has the dimensions of the ladd- and ring-DRPA\n";
  cout <<   "                             eigenvalue problems or the number of 2qp1qh/2qh1qp poles for each channel).\n";
  cout << "\n  MakePivots=n,file           Creates a pivots file based on the input propagator, for use with the BAGEL algorythm.\n";
  cout << "\n  ConvVpp[=outfile]           Loads the interaction (using option 'Vpp') and converts it into a new format. This could be either\n"
       <<   "                             in native .bcd format or a binary file .bin for fast reading. The type of output file depends on thr\n"
       <<   "                             extension of the name given for'outfile'.\n";
  cout << "\n\n" << flush;
  return; }

int parse_cmd_line(int argc, char **argv) {

  if(argc < 2) {
    print_usage();
    return(1);}

  cout << "\nInput line: ";
  for(int iarg=1; iarg<argc; ++iarg) cout << "  " << argv[iarg];
  cout << endl;

  for(int iarg=1; iarg<argc; ++iarg) {

    if(strstr(argv[iarg],"-h")!=NULL) {
         verbose();
         print_usage();
         return(1);}

    if(strstr(argv[iarg],"-v")!=NULL) {
         verbose();
         return(1);}

    if(strstr(argv[iarg],"WkFold")!=NULL) {
	  if (1 == sscanf(argv[iarg],"WkFold=%s", temp_str)) 
		sprintf(BcDorWorkFolder,temp_str);
	  continue;}
	
    if(strstr(argv[iarg],"MakeMdSp")!=NULL) {
         WelcheRec = 31;
         continue;}

    if(strstr(argv[iarg],"MakeSpProp")!=NULL) {
         WelcheRec = 32;
         i_gsp_rw =  0;
         if (1 != sscanf(argv[iarg],"MakeSpProp=%i",&i_gsp_rw )) i_gsp_rw =  0;
         continue;}

    if(strstr(argv[iarg],"ConvSpProp")!=NULL) {
         WelcheRec = 33;
         i_gsp_rw =  -1;
         if (1 != sscanf(argv[iarg],"ConvSpProp=%i",&i_gsp_rw )) i_gsp_rw =  -1;
         continue;}

    if(strstr(argv[iarg],"SpPropStat")!=NULL) {
         WelcheRec = 35;
         continue;}

    if(strstr(argv[iarg],"MakePivots")!=NULL) {
         WelcheRec = 40;
         if (2 != sscanf(argv[iarg],"MakePivots=%i,%s", &i_npvts,gen_fout)) {
           i_npvts = 1;
           sprintf(gen_fout,"DysPivots.bcd-wt");
           }
         continue;}

    if(strstr(argv[iarg],"ConvVpp")!=NULL) {
         WelcheRec = 45;
         if (1 != sscanf(argv[iarg],"ConvVpp=%s", Vpp_out_file)) 
                  sprintf(Vpp_out_file,"Vpp.bin-wt");
         continue;}

    if(strstr(argv[iarg],"Vc")!=NULL) {
         AddVc = 1;
         if (1 != sscanf(argv[iarg],"Vc=%s", Vcpp_file)) 
                  sprintf(Vcpp_file,"Vpp_Coul.bcd");
         continue;}

    if(strstr(argv[iarg],"Radii")!=NULL) {
	  //
	  WelcheRec = 44;
	  // Defaults:
	  if (!RrmsCalc) { 
		RrmsCalc = 1;
		sprintf(Rrmspp_file,"Rrmspp_1.bcd");
	  }
	  continue;}
	
    if(strstr(argv[iarg],"Rrms")!=NULL) {
	  //
	  // Defaults:
	  RrmsCalc = 1;
	  sprintf(Rrmspp_file,"Rrmspp_1.bcd");
	  //
	  if (2 == sscanf(argv[iarg],"Rrms=%i,%s", &RrmsCalc, Rrmspp_file)) continue;
	  if (1 != sscanf(argv[iarg],"Rrms=%i", &RrmsCalc) ) RrmsCalc = 1;  // Default
	  switch (RrmsCalc) {
		case 2:
		  sprintf(Rrmspp_file,"Rrmspp_2.bcd");  // Default
		  break;
		default:
		case 1:
		  sprintf(Rrmspp_file,"Rrmspp_1.bcd");  // Default
		  break;
	  }
	  continue;}
	
	
    if(strstr(argv[iarg],"SpProp")!=NULL) {
         if (1 != sscanf(argv[iarg],"SpProp=%s", SpProp_file)) 
                  sprintf(SpProp_file,"sp_prop");
         continue;}

//     if(strstr(argv[iarg],"OutProp=")!=NULL) {
//          if (1 == sscanf(argv[iarg],"OutProp=%s", gout_str)) 
//                   i_gout_str=1;
//          continue;}

    if(strstr(argv[iarg],"MdSp")!=NULL) {
         if (1 != sscanf(argv[iarg],"MdSp=%s", MdSp_file)) 
                  sprintf(MdSp_file,"input_msp");
         continue;}

    if(strstr(argv[iarg],"Trel_pipj")!=NULL) {
         if (1 != sscanf(argv[iarg],"Trel_pipj=%s", Trel_1_file)) 
                  sprintf(Trel_1_file,"Trelpp_1.bcd");
         continue;}

    if(strstr(argv[iarg],"Trel_2body")!=NULL) {
         if (1 != sscanf(argv[iarg],"Trel_2body=%s", Trel_2_file)) 
                  sprintf(Trel_2_file,"Trelpp_2.bcd");
         continue;}

    if(strstr(argv[iarg],"Trel")!=NULL) {
         IU1body = 3;
         if (1 != sscanf(argv[iarg],"Trel=%i", &IU1body)) 
                  IU1body=3;
         continue;}

//     if(strstr(argv[iarg],"a_Tcm")!=NULL) {
//          a_Tcm = 0.0;
//          if (1 != sscanf(argv[iarg],"a_Tcm=%lf", &a_Tcm)) a_Tcm=0.0;
//          continue;}

	if(strstr(argv[iarg],"ExtU1")!=NULL) {
	  i_Ext_U1 = 0;
	  if (1 != sscanf(argv[iarg],"ExtU1=%i", &i_Ext_U1)) {i_Ext_U1=0;  continue;}
	  else {
		Ext_U1_Uch[0] = 0.0;
		Ext_U1_Uch[1] = 0.0;
		switch (i_Ext_U1) {
		  case 3:
		  case 4:
			Ext_U1_Uo = 0.0;
			if (3 == sscanf(argv[iarg],"ExtU1=%i,%lf,%lf", &n1, &x1, &x2)) {Ext_U1_Uch[0] = x1; Ext_U1_Uch[1] = x2;}
			break;
		  default:
			Ext_U1_Uo = 1.0;
			if (1 < sscanf(argv[iarg],"ExtU1=%i,%lf,%s", &n1, &x1, Ext_U1_file)) Ext_U1_Uo = x1;
			break;
		}
	  }
	  continue;}
	
	
    if(strstr(argv[iarg],"AddVpp")!=NULL) {
         if (iAddVpp >= AddVppMx) continue;
         if (2 == sscanf(argv[iarg],"AddVpp=%lf,%s", AddVppMult+iAddVpp, Vpp_add_file[iAddVpp])) ++iAddVpp;
         continue;}

    if(strstr(argv[iarg],"Vpp")!=NULL) {
         if (1 != sscanf(argv[iarg],"Vpp=%s", Vpp_file)) 
                  sprintf(Vpp_file,"Vpp.bcd");
         continue;}

    if(strstr(argv[iarg],"ItrMax")!=NULL) {
         ItrMax = 0;
         if (1 != sscanf(argv[iarg],"ItrMax=%i", &ItrMax)) 
                 ItrMax = 0;
         continue;}

    if(strstr(argv[iarg],"Lanczos")!=NULL) {
         Lanczos = 2;
         if (1 != sscanf(argv[iarg],"Lanczos=%i", &Lanczos)) 
                  Lanczos = 2;
         continue;}

    if(strstr(argv[iarg],"SelCharge")!=NULL) {
      sel_charge_flag = true;
      if (1 != sscanf(argv[iarg],"SelCharge=%i", &i_sel_charge)) sel_charge_flag=false;
      continue;}
	
    if(strstr(argv[iarg],"SelInt")!=NULL) {
      i_SelInt = -1;
      if (1 == sscanf(argv[iarg],"SelInt=%i", &n1)) i_SelInt = n1;
      continue;}
    

	
    if(strstr(argv[iarg],"LanDysPiv")!=NULL) {
         if (1 != sscanf(argv[iarg],"LanDysPiv=%s", DysPivots_file)) 
                  sprintf(DysPivots_file,"DysPivots.bcd");
         continue;}

    if(strstr(argv[iarg],"2nd")!=NULL) {
         i_HF = 0;
         i_2nd_ord = 1;
         i_ExtSE = 0;
         ItrMax = 0;
         if (1 == sscanf(argv[iarg],"2nd=%i", &n1)) ItrMax = n1;
         continue;}

    if(strstr(argv[iarg],"HF")!=NULL) {
         i_HF = 1;
         i_2nd_ord = 0;
         i_ExtSE = 0;
         ItrMax = 0;
         if (1 == sscanf(argv[iarg],"HF=%i", &n1)) ItrMax = n1;
         continue;}

	
	if(strstr(argv[iarg],"nDfw")!=NULL) {
	  nDfw = -100;
	  if (1 != sscanf(argv[iarg],"nDfw=%i", &nDfw)) 
		nDfw = -100;
	  continue;}
	
    if(strstr(argv[iarg],"nDbk")!=NULL) {
	  nDbk = -100;
	  if (1 != sscanf(argv[iarg],"nDbk=%i", &nDbk)) 
		nDbk = -100;
	  continue;}
	
	
    if(strstr(argv[iarg],"ExtSE")!=NULL) {
         i_HF = 0;
         i_2nd_ord = 0;
         i_ExtSE = 1;
         ItrMax = 0;
         if (1 == sscanf(argv[iarg],"ExtSE=%i", &n1)) ItrMax = n1;
         continue;}

     if(strstr(argv[iarg],"SetEf=")!=NULL) {
          if (2 == sscanf(argv[iarg],"SetEf=%i,%lf", &n1,&x1)) 
            Set_external_EFermi(n1,x1);
          continue;}

    if(strstr(argv[iarg],"MakeSelfEn")!=NULL) {
         WelcheRec = 6;
         i_subsh = -1;
         i_FwBk = 0;
         if (strstr(argv[iarg],"MakeSelfEnFw")!=NULL) i_FwBk = +1;
         if (strstr(argv[iarg],"MakeSelfEnBk")!=NULL) i_FwBk = -1;
         if (1 == sscanf(argv[iarg],"MakeSelfEnFw=%i",&n1 )) {i_subsh = n1; i_FwBk = +1; continue;}
         if (1 == sscanf(argv[iarg],"MakeSelfEnBk=%i",&n1 )) {i_subsh = n1; i_FwBk = -1; continue;}
         if (1 == sscanf(argv[iarg],"MakeSelfEn=%i"  ,&n1 )) {i_subsh = n1; continue;}
         continue;}

    if(strstr(argv[iarg],"PlotSelfEn")!=NULL) {
         WelcheRec = 7;
         i_subsh = -1;
         if (1 == sscanf(argv[iarg],"PlotSelfEn=%i"  ,&n1 )) {i_subsh = n1; continue;}
         continue;}
         
    if(strstr(argv[iarg],"Koltun")!=NULL) {
         WelcheRec = 10;
         continue;}

    if(strstr(argv[iarg],"MBPT2")!=NULL) {
         WelcheRec = 12;
         continue;}

    if(strstr(argv[iarg],"CCD")!=NULL) {
      WelcheRec = 20;
      aLM_CCD = 0.5;
      Stp_CCD = 1;
      if (2 == sscanf(argv[iarg],"CCD=%lf,%i",&x1, &n1 )) {aLM_CCD=x1; Stp_CCD = n1; continue;}
      if (1 == sscanf(argv[iarg],"CCD=%lf"   ,&x1      )) {aLM_CCD=x1; continue;}
    continue;}
    

    
    if(strstr(argv[iarg],"bHO=")!=NULL) {
         //bHO = -1.0;
         hwHO = -1.0;
         if (1 == sscanf(argv[iarg],"bHO=%lf", &bHO)) continue;}

    if(strstr(argv[iarg],"hwHO=")!=NULL) {
         bHO = -1.0;
         //hWHO = -1.0;
         if (1 == sscanf(argv[iarg],"hwHO=%lf", &hwHO)) continue;}

    if(strstr(argv[iarg],"mass_p")!=NULL) {
      if (1 == sscanf(argv[iarg],"mass_p=%lf", &mass_proton)) {
        mass_particle_ave = 2.0*mass_proton*mass_neutron/(mass_proton+mass_neutron);
        i_set_new_mass = 1;
        continue;}
    }
    
    if(strstr(argv[iarg],"mass_n")!=NULL) {
      if (1 == sscanf(argv[iarg],"mass_n=%lf", &mass_neutron)) {
        mass_particle_ave = 2.0*mass_proton*mass_neutron/(mass_proton+mass_neutron);
        i_set_new_mass = 1;
        continue;}
    }
    
    if(strstr(argv[iarg],"mass_e")!=NULL) {
      if (1 == sscanf(argv[iarg],"mass_e=%lf", &mass_electron)) {
        mass_particle_ave = mass_electron;
        i_set_new_mass = 1;
        continue;}
    }
    
    if(strstr(argv[iarg],"mass_ave")!=NULL) {
      if (1 == sscanf(argv[iarg],"mass_ave=%lf", &mass_particle_ave)) {
        i_set_new_mass = 1;
        continue;}
    }

    
    
    
    if(strstr(argv[iarg],"A=")!=NULL) {
         //A = -1;
         if (1 == sscanf(argv[iarg],"A=%i",&A )) continue;}

    if(strstr(argv[iarg],"Z=")!=NULL) {
         //Z = -1;
         if (1 == sscanf(argv[iarg],"Z=%i",&Z )) continue;}

    if(strstr(argv[iarg],"N=")!=NULL) {
         //N = -1;
         if (1 == sscanf(argv[iarg],"N=%i",&N )) continue;}

//     if(strstr(argv[iarg],"ck")!=NULL) {iCheck = 1; continue;}


//-------------------
//  Test stuff:
	if(strstr(argv[iarg],"hfbshift")!=NULL) {
	  if (2 == sscanf(argv[iarg],"hfbshift=%lf", &x1)) {
		i_hfbshift[0]=1;   mu_hfbshift[0]=x1;
		i_hfbshift[1]=1;   mu_hfbshift[1]=x1;
	  }
	  if (2 == sscanf(argv[iarg],"hfbshift_n=%lf", &x1)) {
		i_hfbshift[0]=1;   mu_hfbshift[0]=x1;
	  }
	  if (2 == sscanf(argv[iarg],"hfbshift_p=%lf", &x1)) {
		i_hfbshift[1]=1;   mu_hfbshift[1]=x1;
	  }
	  continue;}
//-------------


    cout << "\n\n BcDor; WARNING: the option '" << argv[iarg] << "' is not recognized.\n\n";
    cerr << "\n\n BcDor; WARNING: the option '" << argv[iarg] << "' is not recognized.\n\n";

  }

  return 0;
  }
