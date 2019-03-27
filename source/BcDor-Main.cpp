//
//   This is the main Green's function program at 2nd order. It uses routines
//  from the BcDor library to do:
//     - standard manipulation of propagator, mod-sp, and interaction files
//     - an Hartree-Fock calculation
//     - compute the Midgal-Galitski-Koltun sum rule for a given propagator
//     - a calculation with 2nd order self-energy 
//  Still to be included:
//     - solve the Dyson Eq. in the BAGEL approximation, hence allowing
//         iterations
//     - solve the dressed ph-DRPA equations (with and without charge ex.)
//
//
//   Only harmonic oscillator wave function are supported here (however
//  the library can in principle handle any basis...).
//
//   This program runs getting all the input information from the command
//  line. Type
//    "BcDor.exe -h"
//  for help.
//
//
//                                  C.Barbieri, RIKEN, May 2010.
//


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <cmath>
using namespace std;

#include "BcDor-Global_variables.hh"
#include "BcDor-Run_vars.hh"


int Run_itrs(ModSpace_t*, VppInt_t*, SpProp_t*);
int Miscellanea(int );
int Miscellanea2(int );
int DysonMatrix(Wide_self_energy_t*, SpProp_t*, double*, double*, double* );

//
// void plot_str_dist(char*, double[], double[], int);
// // 

static double   x1,x2,x3;
static int n1,n2,n3,n4,n5;

static char chip[2] = { '+' , '-'};

//==============================================================================



int main(int argc, char **argv) {

  //
  // Initializaions
  BcDor_Init();

  //
  // some defaults values:
  Set_defaults();

  if (n1=parse_cmd_line(argc, argv)) return (n1-1);


  if (30<WelcheRec && WelcheRec<60) {Miscellanea(WelcheRec);  Unset_all(); exit(0);}

  Set_nucleus_data();
  Load_ModSp();
  Load_interaction();
  //Load_SpProp();

  if (0<WelcheRec && WelcheRec<=30) {Miscellanea2(WelcheRec); Unset_all(); exit(0);}


  SpProp_t propV(SpProp_file,MS);
  if (propV.n_qp + propV.n_qh < 300) propV.cout_qplist();  // if too many poles, better not to output them...

  sprintf(gout_file,"sp_prop-SC-A=%i_itr%%i-wt", A);
  if (i_gout_str) sprintf(gout_file,"sp_prop-%s_finalitr%%i-wt", gout_str);

  cout << Run_itrs(MS, Vpp, &propV) << " iterations have been executed.\n\n";

  
  Unset_all(); // deallocates Rrms, ...  

  return 0;
  }



static SpProp_t *g_ptr1= NULL, *g_ptr2=NULL, *g_ptr3=NULL;

static double   U1body_tot[2],E_tot[2],A_tot[2],Vtot_1st_ord;


//
// This is the main loop routine that iterates the Dyson equation
//  If I_ExtSE or i_2nd_ord, then only a Hartree-Foch is done.
//
int Run_itrs(ModSpace_t *MSin, VppInt_t *Vppin, SpProp_t *gin) {

  double prop_diff;
  int itr, ish;
  int NDIM=0;
  for(ish=0; ish<MSin->nsubsh; ++ish) if (NDIM<MSin->MSp_no[ish]) NDIM=MSin->MSp_no[ish];
  //  NDIM is the max dimension of the HF matrix among all possible
  // partial waves

  

  //  The class VectList can be used to contain some set of vectors...
  //  Here, an object is created to sotre the set of Lanczos pivots (in (B)HF(B) space)
  // that are read in from file for ever new partial wave (ish).
  int n_pvt_alloc = (10 > NDIM) ? 10 : NDIM; //usually, one needs NDIM pvts
  VectList Pivots(NDIM,n_pvt_alloc); // decleare memory for the Dyson's pivots
  Pivots.itype = Lanczos;
  
  

  double *Uab = new double[NDIM*NDIM];
  for(int i= 0; i<NDIM*NDIM; ++i) Uab[i]=0.0;


  //
  // Bounds of the tot number of poles to allocate for the propagators
  if (i_HF) {nDfw=nDbk=NDIM;}
  cout << "\n\n\nNDIM,nDfw,nDbk=" << NDIM << " ,  " << nDfw << " ,  " << nDbk << endl << endl;
  

  //
  // NOTE The following is just to estimate the number of qps and qhs
  //  that need to be allocated. It was a quick try and will need to be
  //  coded better at some point...
  if (Lanczos) {
    //---
	// questa non e` una stima sicura al 100%...
	//int nLczItrs = 200; // questo e` il numero di iterazioni x ogni pivot... 
    //                    // ... di solito 200 bastano.
    //for(ish=0; ish<MSin->nsubsh; ++ish) n3 += nLczItrs * MSin->MSp_no[ish];
    //n4 = n3;
    //---

	// meglio calcolare le dimensioni esatte:
	n3 = n4 = 0; n5 = NDIM;
	for(ish=0; ish<MSin->nsubsh; ++ish) {
      // reads # of qps / qhs for the pivots file: 
	  Pivots.n_vcts = 0;
	  Make_LancPivots(MSin, &Pivots.n_vcts, Pivots.vect, Pivots.ia, Pivots.ib, ish);
 	  for (n1=0; n1<Pivots.n_vcts; ++n1 ) {if (0 < Pivots.ia[n1]) n3 +=Pivots.ia[n1];
	  									   if (0 < Pivots.ib[n1]) n4 +=Pivots.ib[n1];}
	  
	  // then add space for the HF part of the dyson matrix:
	  n3 += MSin->MSp_no[ish];
	  n4 += MSin->MSp_no[ish];
   }

 	n3 += 30;  //  add extra 30 spaces, they will stay 
    n4 += 30;  //  empty bu JUST for safety...
  } else {
    n3 = n4 = NDIM; n5 =0;
    //
    //  dichiara un self energia solo per usare la routine che 
    // calcola il numero di poli..
    Wide_self_energy_t tmp(MSin,Vppin,gin);
	bool i_calc = true;
	if (2 * MSin->n_spb_tot < (gin->n_qh + gin->n_qp)) i_calc = false;
    if (i_calc) for(ish=0; ish<MSin->nsubsh; ++ish) {
	  tmp.set_clj(ish);
	  tmp.Count_2p1h_2h1p_poles(&n1 , &n2);
	  n3 += n1 + MSin->MSp_no[ish];
	  n4 += n2 + MSin->MSp_no[ish];
	}
	if (! i_calc) {n3 = 200000; n4 = 30000;} 
  } // 'tmp' viene deallocata qui.
  cout << " Max number of expected qps =" << n3 << endl;
  cout << " Max number of expected qhs =" << n4 << endl;

  //
  if (nDfw < 0) {n1=n3;} else {n1 = (nDfw+n5)* MSin->nsubsh;} //   NOTE: n5 ==0 if no BAGEL is
  if (nDbk < 0) {n2=n4;} else {n2 = (nDbk+n5)* MSin->nsubsh;} //  done, otherwise n5==NDIM.
  if (n1 < gin->n_qp) n1 = 10 + gin->n_qp;
  if (n2 < gin->n_qh) n2 = 10 + gin->n_qh;
  if (n1 > 100000000) n1 = 10 + gin->n_qp;
  if (n2 > 100000000) n2 = 10 + gin->n_qh;
  //
  //
  //  gin e` il propagatore di input viene usato solo alla
  // prime iterazione e poi non viene piu` toccato. 'gin' e`
  // ancora presente quando si ritorna alla main routine
  //
  //  g_ptr2 e g_ptr3 vengono dichiarati in modo a contenere
  // abbastanza spazio per il n. di poli generati durante le
  // iterazioni
  g_ptr1 = gin;
  g_ptr2 = new SpProp_t(n1,n2,MSin);
  g_ptr3 = new SpProp_t(n1,n2,MSin);


  //Wide_self_energy_t sfV(MSin,Vppin,g_ptr1); //it has to be reset anyway...
  Wide_self_energy_t sfV, sfV_lncz; // sfV       is the self-energy
	                                // sfV_lncz  will be used for the reduced Krylov one
						            //           and will remain empty if no Lanczo in made

  if ( (ItrMax>0) && (!i_HF) && (!Lanczos)&&(!i_ExtSE) ) {
    cout << "\n\n WARNING: one shold not make itrations beyond HF without at least Lanzcos!\n\n";
    ItrMax=0;
    }


  // loop delle itarazioni
  for(itr=0; (itr<=ItrMax)||(ItrMax<0); ++itr) {
    //cout << "\n START itr. n. "<<itr<<" of " <<ItrMax<<"\n==================================\n";

    //
    // Inizializzaztioni
    //

    sfV.reset(MSin,Vppin,g_ptr1); // **must** reset here because gptr1 changes!
	sfV_lncz.reset(MSin,Vppin,g_ptr1);

    g_ptr2->Clean_up(); // this will be new output

    // Koltun sum rule
    for(int i=0; i<2; ++i) {
      U1body_tot[i]=0.0;
      E_tot[i]=0.0;
      A_tot[i]=0.0;
      }

   
    //
    // Nuovo propagatore: self-energia e Dyson sono calcolati onda per onda)
    //

    // loop su tutte le onde parziali (
    for(ish=0; ish<MSin->nsubsh; ish++) {
      cout << " ---------------------------------------\n Subshell: "
           << ish << "   " << MSin->MSp_name[ish] << endl << endl;

	  if ( sel_charge_flag && (i_sel_charge != MS->MSp_ch[ish]) ) {
		cout << " -- skip this partial wave (it will be copied from the one with charge="<<i_sel_charge<<")...\n\n";
		continue;
	  }

      sfV.set_clj(ish);
	  sfV_lncz.set_clj(ish);

      //
      // Construct the self energy:
      //

      //  HF(B) diagram
      sfV.Build_ExtendedHartreeFock();

      //
      // Add the kinetic energy (the self-energy is usually defined as non
      //  containing the kin. en., but T_kin must be included in the Dyson
      //  eq.  and we do it here by including it as an external 1-body
      //  potential).
      //
      if (7 != IU1body) {
        Get_Tab_kin_ho(Uab, NDIM, ish, MSin);
        for(int i= 0; i<NDIM*NDIM; ++i) Uab[i]*=U1body_fac;
        sfV.Add_1bd_SelfEn(Uab,NDIM);
      }

	  //
      // Add the external potential:
      //
      if (i_Ext_U1) {  // all this stuff to be moved to a 'retrive self-en' file...
        for(int i= 0; i<NDIM*NDIM; ++i) Uab[i] = 0.0;
        //
        switch(i_Ext_U1) {
          case(1):
		  {int ndim = (NDIM  < MSin->MSp_no[ish]) ? NDIM : MSin->MSp_no[ish];
            for(int i= 0; i<ndim; ++i) {
              Uab[i*(NDIM+1)] = MSin->MSp_xe[ish][i];
              cout << i << "    - -   " << Uab[i*(NDIM+1)] << endl;
			} }
            break;
          case(3):
			break;
		  {int ndim = (NDIM  < MSin->MSp_no[ish]) ? NDIM : MSin->MSp_no[ish];
			double x1 =0;
			x1 = Ext_U1_Uch[MSin->MSp_ch[ish]];
            for(int i= 0; i<ndim; ++i) {
              Uab[i*(NDIM+1)] += x1;
			  cout << i << "    - -   " << Uab[i*(NDIM+1)] << "  ( +  "<< x1<<")\n";
			} }
            break;
          case(2):
            char line[300], *line_ptr;
            int ext_ndim, ext_2j, ext_ip, ext_ch;
            int n1, n2, n3;
            double x1;
            FILE *infile;
            infile = fopen(Ext_U1_file, "rb");
            if (infile == NULL) {
              cerr << "ERROR: Cannot open the the external potential file ('" << Ext_U1_file << "')\n";
              break;
			}
            ext_ndim = -100;
            while(!feof(infile)) {
              fgets(line,sizeof(line), infile);
			  //if (4 == sscanf(line_ptr,  "dim=%i%i%i%i%n",&ext_ndim,&ext_2j,&ext_ip,&ext_ch,&n1) ) {...};
              if (4 == sscanf(line,  "dim=%i%i%i%i",&n1,&ext_2j,&ext_ip,&ext_ch) ) 
				if ((ext_2j==MSin->MSp_2j[ish]) && (ext_ip==MSin->MSp_ip[ish])
					&& (ext_ch==MSin->MSp_ch[ish])) {ext_ndim=n1; break;}
            }
            if (ext_ndim < 1) {
              cout << "\nWARNING: Cannot find the external potential for channel 2j,ip,ch="
			  << MSin->MSp_2j[ish] <<"/2"<< chip[MSin->MSp_ip[ish]]  <<","
			  << MSin->MSp_ch[ish] << "  (from file '" << Ext_U1_file << "')\n";
              break;
			}
            ext_ndim = (ext_ndim < NDIM) ? ext_ndim : NDIM;
            cout << "\n Adding the following external potential (J\\pi="<<ext_2j <<"/2"<< chip[ext_ip]
			<<", ch=" << ext_ch << "  and  with Uo="<<Ext_U1_Uo<<"):\n";
            for(n1=0; n1<ext_ndim; ++n1) {
              fgets(line,sizeof(line), infile);
              line_ptr = line;
              cout << n1;
              for(n2=0; n2<ext_ndim; ++n2)
                if (0 < sscanf(line_ptr,  "%lf%n",&x1,&n3) ) {
                  line_ptr += n3;
                  x1 *= Ext_U1_Uo;
                  Uab[n1*NDIM + n2] = x1;
                  cout << "     " << x1;
				}
              cout << endl;
            }
            fclose(infile);
            break;
        } // end of switch(i_Ext_U1)
        //
        sfV.Add_1bd_SelfEn(Uab,NDIM);
      }
	  
      if (i_2nd_ord) {
        //
        // Slef-enegy at 2nd order
        sfV.Build_Dynamic_SlfEn_2ndOrder(&n3, &n4, 1, icut_fd, ntcut_fd);
        cout << " number of 2p1h/2h1p poles in the self-en.: " << n3 << " / " << n4 << endl;

      } else if (i_ExtSE) {
        //
        // Reads the self-energy frmom file (instead of calculating one)...
        sfV.load_SelfEn_bin();//n3 , n4);
        cout << " number of Ext-2p1h/2h1p poles in the self-en.: " << sfV.N_PLS_fw << " / " << sfV.N_PLS_bk << endl;
      //
      }


      //
      // Solve Dyson:
      //
      if (Lanczos) {
        //
        // Usa prima una proiezione nello spazio di Krylov
	    // per ridurre i poli della self-energia:
	
	    cout << "\n Lanczos "<< Lanczos<<":\n============\n";
	
  	    // legge i pivots dal file:
		Pivots.n_vcts = 0;
		Make_LancPivots(g_ptr1->MdSp, &Pivots.n_vcts, Pivots.vect, Pivots.ia, Pivots.ib, ish, 1);

		// 	sfV_lncz  ...was already declared above
	    sfV_lncz.Reduce_SelfEn_Lncz(&sfV,Pivots.n_vcts, Pivots.vect, Pivots.ia, Pivots.ib, 1);

		n3 = DysonMatrix(&sfV_lncz,g_ptr2,&x1,&x2,&x3);
		
      } else  {

        //
		// conto di HFB oppure di Gorkov al 2nd-ord senza lanczos:
		//
		
		n3 = DysonMatrix(&sfV,g_ptr2,&x1,&x2,&x3);

      }
      U1body_tot[MSin->MSp_ch[ish]] += x1;
      E_tot[MSin->MSp_ch[ish]]   += (x1+x2)/2.0;
      A_tot[MSin->MSp_ch[ish]]    += x3;

    } // end of loop on 'ish' --- onde parziali

	if (sel_charge_flag) {
      U1body_tot[1-i_sel_charge] = U1body_tot[i_sel_charge];
      E_tot[1-i_sel_charge]	   = E_tot[i_sel_charge];
      A_tot[1-i_sel_charge]      = A_tot[i_sel_charge];
	  g_ptr2->Fill_with_charge(i_sel_charge);
	}

  //
  //  Risultati parziali e totali della regola di somma di Koltun
  //
  cout << "\nneutrons:";
  cout << "\nU1body_tot_n     = " << setw(16) << U1body_tot[0];
  cout << "\nV_tot_n          = " << setw(16) << E_tot[0]-U1body_tot[0];
  cout << "\nE_tot_n          = " << setw(16) << E_tot[0];
  cout << "\nA_tot_n          = " << setw(16) << A_tot[0]           << endl;
  cout << "\nprotons:";
  cout << "\nU1body_tot_p     = " << setw(16) << U1body_tot[1];
  cout << "\nV_tot_p          = " << setw(16) << E_tot[1]-U1body_tot[1];
  cout << "\nE_tot_p          = " << setw(16) << E_tot[1];
  cout << "\nA_tot_p          = " << setw(16) << A_tot[1]           << endl;
  cout << "\nU1body_tot     = " << setw(16) << U1body_tot[0]+U1body_tot[1];              
  cout << "\nV_tot          = " << setw(16) << E_tot[0]-U1body_tot[0]+E_tot[1]-U1body_tot[1];
  cout << "\nE_tot          = " << setw(16) << E_tot[0]+E_tot[1];
  cout << "\nE_tot/A        = " << setw(16) << (E_tot[0]+E_tot[1])/double(A);
  cout << "\nA_tot          = "  << setw(16) << A_tot[0]+A_tot[1]           << endl;



  Vtot_1st_ord=g_ptr2->Compute_FirstOrderExpValue(Vppin);
  cout << "\n <V_tot> at 1st ord. = " << Vtot_1st_ord <<"\n";

	//
	//  Calculate and write out the r.m.s. radius if RrmsCalc i set:
	switch (RrmsCalc) {
	  case 1:
		Calc_Rrms_OneBody_rirj(g_ptr2, Rrms, NULL, true);
		break;
	  case 2:
		Calc_Rrms_TwoBody(g_ptr2, Rrms, true);
		break;
	  default:
		break;
	}
	
	
  cout << "\n\n Itr."<< itr << ": a rough estimate of the summed difference of all sp energies of the 2 propagators is: "
       << (prop_diff = g_ptr1->quick_n_dirty_diff(g_ptr2)) << " MeV.\n\n"; 
  

  // ogni volta che una nuova iterazione inizia, g_ptr1 sara` in propagatore
  // di input mentre g_ptr2 contiene quello di 2 iterazioni precedenti (che
  // viene quindi svuotato e usato per il nuovo output).
  if (itr!=0) g_ptr3 = g_ptr1;
  g_ptr1 = g_ptr2;
  g_ptr2 = g_ptr3;
  g_ptr3 = NULL;  // won't try to cancel it later.


  if (prop_diff < Conv_check) break;
 } // end of loop on 'itr'


  //
  //  Write the last two propagators to disk -- this is the output!
  //
  char fmn[BUFFER1];
  sprintf(fmn, gout_file , itr-1);
  g_ptr2->write(fmn);
  sprintf(fmn, gout_file , itr);
  if (itr > 1) sprintf(fmn, gout_file , -999);
  g_ptr1->write(fmn);


  //  You can send the final propagator to stdout (but beyond Hartree-Fock
  // you'll have zilions of output lines printed...).
  if ((i_HF) || (g_ptr1->n_qp + g_ptr1->n_qh < 300)) g_ptr1->cout_qplist();


  if ((g_ptr1==g_ptr2) || (g_ptr1==g_ptr3)) g_ptr1=NULL;
  if  (g_ptr2==g_ptr3)                      g_ptr2=NULL;
  if (NULL!=g_ptr1) delete g_ptr1; g_ptr1=NULL;
  if (NULL!=g_ptr2) delete g_ptr2; g_ptr2=NULL;
  if (NULL!=g_ptr3) delete g_ptr3; g_ptr3=NULL;

  delete [] Uab;  Uab = NULL;

  return itr;}
