
#include <fstream>
#include <iostream>
using namespace std;

#include "BcDor-Global_variables.hh"


// j-coupling routines (public to all functions)
ang_mom_functions am;


// code wide clock functions
Execution_clock exclock;


//
// Folder were self-energies and other
// intermediate propagators are stored.
//
// This is set to 'bcdwk' in file by
// the routtine just below
//
char BcDorWorkFolder[FNAMESTACK];


void BcDor_Init(void) {

  //
  //  SET DEFAULTS
  //

  // Default working directory:
  sprintf(BcDorWorkFolder, "bcdwk");

  return;
  }

