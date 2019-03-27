


#include <iostream>
#include <cstdlib>
using namespace std;

#include "BcDor-Run_vars.hh"


static int n1,n2;
static double x1;

void Set_defaults(void ){

  //
  // Set the default names of the I/O files:
  //
  sprintf(Ext_U1_file,"Usp_ext.bcd");
  sprintf(Vpp_file,"Vpp.bcd");
  sprintf(Vpp_out_file,"Vpp.bin-wt");
  sprintf(SpProp_file,"sp_prop");
  sprintf(MdSp_file,"input_msp");
  sprintf(gout_file,"sp_prop-SC-itr%000i-wt");
  sprintf(Trel_1_file,"Trelpp_1.bcd");
  sprintf(Trel_2_file,"Trelpp_2.bcd");
  sprintf(DysPivots_file,"DysPivots.bcd");

  return;}


void get_part_number(void) {

  while( (0 > N) || (0 > Z) || (N+Z != A) ) {
   A=N=Z=-1;
   while(A<=0) {cout << "\n Mass number?    "; cin >> A;}
   while(N<0)  {cout << "\n Neutron number? "; cin >> N;}
   while(Z<0)  {cout << "\n Proton number?  "; cin >> Z;}
   if (A == N+Z) {
   //cout << "\n Mass number = " << A << endl;
   break;
   }
  }

  return;}

void Set_nucleus_data(void) {

  while(1){
   while(A<=0) {cout << "\n Mass number?    "; cin >> A;}
   while(N<0)  {cout << "\n Neutron number? "; cin >> N;}
   while(Z<0)  {cout << "\n Proton number?  "; cin >> Z;}
   if (A == N+Z) break;
   cout << "\n\n A="<<A<<",  N="<<N<<", and Z="<<Z <<" are not consistent!\n";
   A=N=Z=-1;
  }

  while(bHO<=0.0 && hwHO<=0.0) {cout << "\n H.o. parameter [fm]? "; cin >> bHO;}

  while(IU1body<0 || 7<IU1body) {
    cout << "\n Type of kin. energy (0==no c.m. corr., 1(3)==Trel w/ 1b+2b(ext. pipj), 2(4)==Trel w/ 2b(ext. 2b) only)? ";
    cin >> IU1body;
  }

  switch(IU1body) {
    case(0): U1body_fac=1.0;                   break;
    case(1): U1body_fac=double(A-1)/double(A); break;
    case(2): U1body_fac=0.0;                   break;
    case(3): U1body_fac=double(A-1)/double(A); break;
    case(4): U1body_fac=0.0;                   break;
    case(5): U1body_fac=(double(A-1)+a_Tcm)/double(A); break;
    case(6): U1body_fac=a_Tcm;                 break;
    case(7): U1body_fac=0.0;                   break;
    }


//
//  print_input parameters:
//

  cout << "\n Mass number = " << A << " ,  Itkin = " << IU1body
                      << " ,  U1body_fac = " << U1body_fac << endl;

  if (bHO  > 0.0) cout << "\n bHO = "  << bHO  << endl;
  if (hwHO > 0.0) cout << "\n hwHO = " << hwHO << endl;

  cout << "\n Number of expected fw poles = " << nDfw;
  cout << "\n Number of expected bk poles = " << nDbk << endl << endl;

  return;}


void Load_ModSp(void ) {
  //
  //  Get Mod.Sp.
  // ============

  // Read in the model space:
  if ( NULL != MS) delete MS;
  MS = new ModSpace_t(MdSp_file);
  if (0 != i_set_new_mass) MS->mass = mass_particle_ave;
//MS->write("input_msp-wt");
  MS->cout_msp();
  MS->count_2b_confs(&n1, &n2, 1);
  if (bHO > 0.0) { MS->Set_phys_const_bHO(bHO);    // set b=...fm
          } else { MS->Set_phys_const_hwHO(hwHO);  // set hw=...MeV
          }
  return;}

void Load_interaction(void ) {
  // "Build" the  Interaction:
  // ================

  if ( NULL != Vpp) delete Vpp;
  Vpp = new VppInt_t(Vpp_file,MS);
    if (n1=Vpp->check_for_double_entries()) 
      cerr << "\n WARNING:" << n1 << " Vpp mtx els. have been given twice!\n\n";
  if (iCheck) {
    if (n1=Vpp->seek_missing_Vpp())
      cerr << "\n WARNING:" << n1 << " Vpp mtx els. are missing! (maybe they're zero?)\n\n";
  }

  if (-1 < i_SelInt) {
    cout << "\n\n Selcet interaction # " << i_SelInt << "...    ";
    Vpp->select(i_SelInt);
    cout << "  StrEn[ " << Vpp->get_wrk() << " ] = " << Vpp->StrEn[Vpp->get_wrk()] << endl<< endl;
  }

  if (AddVc) Vpp->add_file((1.0/MS->bho),Vcpp_file, Vpp->n_ints, 0);
  //
  for(int iv=0; iv<iAddVpp; ++iv) Vpp->add_file(AddVppMult[iv], Vpp_add_file[iv], Vpp->n_ints, 0);
  //
  if (3 == IU1body) Vpp->add_file(            (2.0*MS->htom/double(A)),Trel_1_file, Vpp->n_ints, 0);
  if (4 == IU1body) Vpp->add_file(            (2.0*MS->htom/double(A)),Trel_2_file, Vpp->n_ints, 0);
  if (5 == IU1body) Vpp->add_file(((1.0-a_Tcm)*2.0*MS->htom/double(A)),Trel_1_file, Vpp->n_ints, 0);
  if (6 == IU1body) Vpp->add_file(((1.0-a_Tcm)*2.0*MS->htom/double(A)),Trel_2_file, Vpp->n_ints, 0);


  //
  // To calculate the rms intrinsic radius:
  if (RrmsCalc) {
	if ( NULL != Rrms) delete Rrms;
	Rrms = new VppInt_t(MS);
	switch (RrmsCalc) {
	  case 1:
		Load_Rrms_TwoBodyPart(MS, Rrmspp_file, Rrms, "(r_i . r_j)^2");
		break;
	  case 2:
		Load_Rrms_TwoBodyPart(MS, Rrmspp_file, Rrms, "(r_ij)^2");
		break;
	  default:
		cerr << "  Bad RrmsCalc="<<RrmsCalc<<endl; exit(100);
		break;
	}
  }

  return;}

void Load_SpProp(void) {
  return;}


void Unset_all(void) {
  //
  //  Deallocate things that might have been initialized along the way...
  //
  //  THIS FUNCTION IS STILL TO BE COMPLETED....  (C.B. 28/12/2011)

  if (NULL!=Rrms) delete Rrms;
  
  return;}


