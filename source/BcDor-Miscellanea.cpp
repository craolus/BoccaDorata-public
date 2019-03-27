//
//
//  Side programs for generatin template input files, converting 
// interaction files, saving/printing slef-energies and so on...
//
//  C. Barbieri -- RIKEN, May 2010
//


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
//#include <iomanip>  //TST//
using namespace std;

#include <math.h>

#include "BcDor-Global_variables.hh"
#include "BcDor-Run_vars.hh"



static double   x1,x2,x3;
static int n1,n2,n3,n4,n5;

static char chip[2] = { '+' , '-'};


//==============================================================================



int Miscellanea(int ICalcType) {

   if (31==ICalcType) {

    // Read in the model space:
    int norb,eMax,n_fluid=-1;
    cout << "\n norb? "; cin >> norb;
    cout << "\n eMax? "; cin >> eMax;
    while(n_fluid<1 || n_fluid>2) {
       cout << "\n # of fermions (1 or 2)? "; cin >> n_fluid;}
    ModSpace_t MSn(norb,eMax,n_fluid);
    MSn.write("input_msp-wt");
    MSn.cout_msp();
    MSn.count_2b_confs(&n1, &n2);

    }
   else if (32==ICalcType) {
    ModSpace_t MS(MdSp_file);
    cout << "\n\n ...will use the following mod. sp. to create the propagator:\n\n";
    MS.cout_msp();
    MS.Set_phys_const_hwHO(1.0);
    SpProp_t gab_n(&MS,1);
    gab_n.write("sp_prop-wt",i_gsp_rw);
    //gab_n.cout_qplist(i_gsp_rw);
    gab_n.cout_byshell(i_gsp_rw);

    }
   else if (33==ICalcType) {
    //
    // This will need some clean up at some point...
    //  if the reduce opiotn is set only the input
    //  prop is not written on file. When the write flag
    //  is not specified in the input line, it should not
    //  be changed.
    //
    SpProp_t *g_ptr1, *g_ptr2, *g_ptr3;
    ModSpace_t MS(MdSp_file);
    SpProp_t *propV = new SpProp_t(SpProp_file,&MS,&n1);
    if (0 > i_gsp_rw) i_gsp_rw=n1;
    g_ptr1 = propV;  // g_ptr1 will be the output prop.
    //
    if (g_ptr1->n_qp + g_ptr1->n_qh < 300) g_ptr1->cout_qplist();  // if too many poles, better not to output them...
    //g_ptr1->cout_byshell();
    //g_ptr1->count_qp_confs(&n1,&n2,&n3,&n4,&n5);
    g_ptr1->write("sp_prop-wt",i_gsp_rw);
    if (NULL != g_ptr1) delete g_ptr1; g_ptr1 = NULL;
    if (NULL != g_ptr2) delete g_ptr2; g_ptr2 = NULL;
    if (NULL != g_ptr3) delete g_ptr3; g_ptr3 = NULL;
    if (NULL != propV ) delete propV;  propV  = NULL;
   }
   else if (35==ICalcType) {
    ModSpace_t MS(MdSp_file);
    SpProp_t propV(SpProp_file,&MS);
    //propV.write("sp_prop-wt");
    //propV.cout_qplist(i_gsp_rw);
    //propV.cout_byshell();
    propV.count_qp_confs(&n1,&n2,&n3,&n4,&n5,1,ntcut_ld,ntcut_fd);
    cout << "\nMaximum dimensions for the 2qp and 3qp bases:\n"
      << " nmx_pp  = " << n1 << endl
      << " nmx_hh  = " << n2 << endl
      << " nmx_ph  = " << n3 << endl
      << " nmx_pph = " << n4 << endl
      << " nmx_hhp = " << n5 << endl;
   }
   //
   else if (40==ICalcType) {
    ModSpace_t MS(MdSp_file);
    SpProp_t propV(SpProp_file,&MS);
    //propV.write("sp_prop-wt");
    //propV.cout_qplist(i_gsp_rw);
    //propV.cout_byshell();
    //propV.count_qp_confs(&n1,&n2,&n3,&n4,&n5);
    propV.DysPivots(gen_fout,i_npvts,i_npvts);

   }
   else if (44==WelcheRec) {
     Load_ModSp();
     SpProp_t propV(SpProp_file,MS);
	 switch (RrmsCalc) {
	   case 1:
		 Calc_Rrms_OneBody_rirj(&propV, Rrmspp_file, true);
		 break;
	   case 2:
		 Calc_Rrms_TwoBody(&propV, Rrmspp_file, true);
		 break;
	   case 3:
		 //
		 break;
	   default:
		 break;
	 }
	 
   }
   else if (45==ICalcType) {
    //
    if ((AddVc) || ((2 < IU1body)&&(IU1body < 7)) || (RrmsCalc) || (iAddVpp)) Set_nucleus_data();
    Load_ModSp();
    Load_interaction();

     //cout << "Write into file " << Vpp_out_file << "...";
     Vpp->write(Vpp_out_file);
     //cout << "   done!\n";
   }

  return 0;}
  // end of Miscellanea


int Miscellanea2(int ICalcType) {

  if (7==ICalcType) {

    int ish;
    double EFermi;

    SpProp_t gsp(SpProp_file,MS);
    //gsp.cout_qplist();

    Wide_self_energy_t sfV(MS,Vpp,&gsp);

    if (!i_gout_str) sprintf(gout_str,"SelfEn");

    for(ish=0; ish<MS->nsubsh; ++ish) {
       if ((i_subsh >= 0) && (ish != i_subsh)) continue;
      cout << " ---------------------------------------\n Subshell: "
           << ish << "   " << MS->MSp_name[ish] << endl << endl;
      sfV.set_clj(ish);
      Retrive_SelfEn(&sfV, ish, MS);

      EFermi = Seek_EFermi(&gsp , ish);

      sprintf(temp_str,"%s_J2pi=%i%c_dt=%i-static.dat"  , gout_str, MS->MSp_2j[ish], chip[MS->MSp_ip[ish]], MS->MSp_ch[ish]);
      sfV.plot_SelfEn_stc(temp_str,&EFermi);


      Sort_d_double2dim(sfV.N_PLS_fw, sfV.NTOT_ORB+2, 0, sfV.Sig_dyn_fw);
      Sort_d_double2dim(sfV.N_PLS_bk, sfV.NTOT_ORB+2, 0, sfV.Sig_dyn_bk);
      sprintf(temp_str,"%s_J2pi=%i%c_dt=%i-dynamic.dat" , gout_str, MS->MSp_2j[ish], chip[MS->MSp_ip[ish]], MS->MSp_ch[ish]);
      sfV.plot_SelfEn_dyn(temp_str,&EFermi);
      }

     }
  else if (6==ICalcType) {
    //
    //  Generates the self-energy (2nd ord.)
    //
    SpProp_t gin(SpProp_file,MS);
    if (gin.n_qp + gin.n_qh < 300) gin.cout_qplist();  // if too many poles, better not to output them...
    //gin.cout_byshell();
    //gin.count_qp_confs(&n1,&n2,&n3,&n4,&n5);
    
    Wide_self_energy_t sfV(MS,Vpp,&gin);
    //TST//Wide_self_energy_t sfT(MS,Vpp,&gin);
    
    
    int i1,i2,i3;
    int j2fd, ipfd, chfd;
    int ndimse;
    double *dptr1;
    
    
    VectList  *Pivots = NULL;
    Wide_self_energy_t *sfV_lncz = NULL;
    if (0 != Lanczos) {
      int NDIM = 0;
      for(int ish=0; ish<MS->nsubsh; ++ish) if (NDIM < MS->MSp_no[ish]) NDIM = MS->MSp_no[ish];
      int n_pvt_alloc = (10 > NDIM) ? 10 : NDIM;
      Pivots = new VectList(NDIM,n_pvt_alloc);
      Pivots->itype = Lanczos;
      sfV_lncz = new Wide_self_energy_t(MS,Vpp,&gin);
    }
    
    
    
    // loop over the subshells:
    for(int ish=0; ish<MS->nsubsh; ++ish) {
      if ((i_subsh >= 0) && (ish != i_subsh)) continue;
      
      sfV.reset(ish,MS,Vpp,&gin);
      
      if (i_2nd_ord) {
        //
        //  2nd order self energy:
        //
        sfV.Build_Dynamic_SlfEn_2ndOrder(&i1, &i2, 1, icut_fd, ntcut_fd);
        //cout << " number of 2p1h/2h1p poles(V): " <<i1 << " / " <<i2 << endl;
        
      } else {
        cout << "\n\n\n Don't know which approximation to use for constraxting the self energy.\n";
        cout<<        " Did you forget to add '2nd' in the command line??? \n\n Will STOP for now.\n\n";
        exit(100);
      }
      
      
      if (0 != Lanczos) {
        
        cout << "\n\n Will now perform a Lanczos reduction based on Lanczos="<< Lanczos<<" input fleag:";
        cout <<   "\n========================================================================\n";
        
        // legge i pivots dal file:
        if (NULL == Pivots) {cerr << "Supertrouble #29:  EXIT!\n";  exit(100);}
        Pivots->n_vcts = 0;
        //Make_LancPivots(&gin, &(Pivots->n_vcts), Pivots->vect, Pivots->ia, Pivots->ib, ish); //same line as in BcDor_FRPA_v1 but the first argument changes
        Make_LancPivots( MS,  &(Pivots->n_vcts), Pivots->vect, Pivots->ia, Pivots->ib, ish, 1);
        
        // sfV_lncz  ...was already declared above
        sfV_lncz->reset(ish,MS,Vpp,&gin);
        
        sfV_lncz->Reduce_SelfEn_Lncz(&sfV,Pivots->n_vcts, Pivots->vect, Pivots->ia, Pivots->ib, 0); // ..., i_copyMF=0); 'cause is not needed here.
        
        sfV_lncz->save_SelfEn_bin(i_FwBk);
        //sfV_lncz->plot_SelfEn_dyn("SE_cmpd_sts");
        
      } else {
        //
        // if no Lanczos, just save the full self-energy
        //
        
        
        cout << "Total # of poles in the self-energy: n_fw="<<sfV.N_PLS_fw <<" ,   n_bk="<<sfV.N_PLS_bk<<endl;
        
        //if (sfV.N_PLS_fw+sfV.N_PLS_bk > 0) sfV.save_SelfEn_bin(i_FwBk);
        sfV.save_SelfEn_bin(i_FwBk);
        //sfV.plot_SelfEn_dyn("SE_cmpd_sts");
        
      }
      
      
    }  // end 'ish' loop
    
    if (NULL != Pivots  ) delete Pivots;    Pivots   = NULL;
    if (NULL != sfV_lncz) delete sfV_lncz;  sfV_lncz = NULL;
    
  }
  else if (12==ICalcType) {
     //
     // Calculate the MBPT2 contribution to the correlation energy!
     //
     Second_order_diag();
     }
  else if (10==ICalcType) {
     //
     //  Copute Galiski-Migdal-Koltun sum rule and other basic
     // quantities for the input sp. propagator.
     //
     MGKoltun_sumrule();
     }
  else if (20==ICalcType) {
    //
    //  Solve the CCD or equations starting form the given
    // input sp. propagator. NOTE that this algoritm expect
    // the propagator (taken as reference state) is the
    // converged HF one.
    //
    
    cout << "\n\n  Solution of the coupled cluster CCD equations. NOTE that this algoritm\n"
          <<    " expects that the reference state (given by the input propagator) is a\n"
          <<    " CONVERGED Hartree-Fock state. Using a different reference wuold neclect\n"
          <<    " the f_ab contributions in the CCD equartions.\n";

    cout << "\n Will read the reference state fro file '"<< SpProp_file << "'...\n\n" << flush;

    SpProp_t* gin = new SpProp_t (SpProp_file,MS);
    if (gin->n_qp + gin->n_qh < 300) gin->cout_qplist();  // if too many poles, better not to output them...
    
    CoupledCluster_t CCD(MS, gin, Vpp);
    
    CCD.create_CCD_aplitudes();
    CCD.Divide_deltaE();
    
    cout << "\n\n Will solve the CCD equations using aLM_CCD="<<aLM_CCD<<"   Stp_CCD="<< Stp_CCD<<"\n\n";
    
    int ierr = CCD.Solve_CCD(aLM_CCD, Stp_CCD);
    
    
  }

  return 0;}

void Retrive_SelfEn(Wide_self_energy_t *sfV, int ish, ModSpace_t *MSin) {

      sfV->Build_ExtendedHartreeFock();

      if (i_2nd_ord) {
        sfV->Build_Dynamic_SlfEn_2ndOrder(&n3, &n4, 1, icut_fd, ntcut_fd);
        cout << " number of 2p1h/2h1p poles(V): " << n3 << " / " << n4 << endl;
      //
      } else if (i_ExtSE) {
        n3 = n4 = -1;
        sfV->load_SelfEn_bin();//n3 , n4);
        
        cout << " number of Fd-2p1h/Fd-2h1p poles(V): " << sfV->N_PLS_fw << " / " << sfV->N_PLS_bk << endl;
        //
      }

  return; }


static int nfb;

void Make_LancPivots(ModSpace_t *MdSp_in, int* n_pvt, double** pvt_list, int n_Litr_fw[], int n_Litr_bk[], int ish, int iplot/*=0*/) {

  if (1 < Lanczos) {
	int n1,n2,n3;
	double x1;
	char line[300], *line_ptr;
	FILE *infile;
	infile = fopen(DysPivots_file, "rb");
	if (infile == NULL) {
	  cerr << "  Cannot open the pivots file ('" << DysPivots_file << "')\n";
	  return;
	}
	(*n_pvt) = 0;
	
	while(1) {
      n3 = -1; while ((n3<2)&& (!feof(infile)) ) {
		fgets(line,sizeof(line), infile);
		n3=sscanf(line,  "shell,npvts:%i%i\n", &n1,&n2);}
      
	  if ( feof(infile) ) {cerr <<"\nEND OF FILE encountered while reading the pivots!!!\n"; break;}
      
	  if (n1 != ish) continue;
	  
	  (*n_pvt) = 0;
	  n4 = n2;
	  for(int i=0; i<n4; ++i)  {
		fgets(line,sizeof(line), infile);
		line_ptr = line;
		if (2 > sscanf(line_ptr,  "%i %i%n",&n2,&n3,&n1) ) continue;
		
		n_Litr_fw[(*n_pvt)] = n2;
		n_Litr_bk[(*n_pvt)] = n3;
		line_ptr += n1;
		for(nfb=0; nfb<MdSp_in->MSp_no[ish]; ++nfb) {
		  sscanf(line_ptr,  "%lf%n",&x1,&n1);
		  pvt_list[(*n_pvt)][nfb] = x1;
		  line_ptr += n1;
		}
		++(*n_pvt);
		if (iplot) {
		  cout  << "  " << i << "   -->" << n2 << "   " << n3 << "  --";
		  for(nfb=0; nfb<MdSp_in->MSp_no[ish]; ++nfb) cout  << "   " << pvt_list[i][nfb];
		  cout << endl;
		}
		//cout << "        " << line;
		
	  }
	  fclose(infile);
	  return;
      
	}
	
	fclose(infile);
	return;
  }
  
  return;}

