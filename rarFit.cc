/************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarFit.cc,v 1.20 2014/09/14 17:33:51 fwilson Exp $
 *          Main function
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 ************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This is the main program for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This is the main program for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "TROOT.h"
#include "TApplication.h"
#include "TStopwatch.h"

#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "rarDatasets.hh"
#include "rarMLFitter.hh"

extern Int_t doBanner();  // reference to RooFit's banner

/// \brief Usage page
///
/// Print a short help page.
void rarFitUsage(TString myCommand)
{
  cout <<"Usage for RooRarFit ("<< RARFIT_VERSION <<") :" << endl << endl
       <<myCommand<<" [-options] <RooRarFit_Config_file>"<<endl
       <<"\t-h this help page"<<endl
       <<"\t-D <data input section>"
       <<" (default \"Dataset Input\")"<<endl
       <<"\t-C <mlFitter config section>"
       <<" (default \"mlFitter Config\")"<<endl
       <<"\t-A <fitter action section>"
       <<" (default \"Fitter Action\")"<<endl
       <<"\t-t <toy job id> (default 0)"<<endl
       <<"\t-n <toyNexp> (default 0, use config)"<<endl
       <<"\t-d <toy dir> (default .toyData)"<<endl
       <<"e.g."<<endl
       <<"\tTo run "<<myCommand<<" from config file demo.config"<<endl
       <<myCommand<<" demo.config"<<endl
       <<"\tand with mlfitter config from section [myMLFitter]"<<endl
       <<myCommand<<" -C myMLFitter demo.config"<<endl
       <<"\tand with mlfitter action from section"
       <<" [my Fit Action]"<<endl
       <<myCommand<<" -C myMLFitter"
       <<" -A \"my Fit Action\" demo.config"<<endl
       <<endl;
}

/// \brief Main program of the mlFitter
///
/// After it parses all the command line options,
/// with normal operation,
/// it first creates a rarDatasets object through which all the datasets
/// are read in, then it instantiates mlFitter class, rarMLFitter,
/// so that all the PDFs are created,
/// and finally it sets toyID (random seed) and
/// calls rarMLFitter::run() of the fitter to finish the job.
/// If successfully compiled and linked,
/// type rarFit in workdir to see short help page.
int main(int argc, char **argv)
{
  Int_t optFlag;
  TString dataInputSec="Dataset Input";
  TString fitterConfigSec="mlFitter Config";
  TString fitterActionSec="Fitter Action";
  Int_t randomBase(0);
  Int_t toyID(0);
  Int_t toyNexp(0);
  TString toyDir(".toyData");
  while (EOF!=(optFlag=getopt(argc, argv, "hD:C:A:t:n:d:"))) {
    switch (optFlag) {
    case 'h' :
      rarFitUsage(argv[0]);
      return 0;
      break;
    case 'D' :
      dataInputSec=optarg;
      break;
    case 'C' :
      fitterConfigSec=optarg;
      break;
    case 'A' :
      fitterActionSec=optarg;
      break;
    case 't' :
      toyID=atoi(optarg);
      break;
    case 'n' :
      toyNexp=atoi(optarg);
      break;
    case 'd' :
      toyDir=optarg;
      break;
    }
  }

  if (optind >= argc) {
    rarFitUsage(argv[0]);
    return 0;
  }
  
  // always run with at least one argument, the config file
  TString ConfigFile=argv[optind];
  ifstream ifs(ConfigFile);
  if(ifs.fail()) {
    cout<<argv[0]<<": can not open config file "<<ConfigFile<<endl;
    exit(-1);
  }
  ifs.close();
  
  // started
  TStopwatch timer;
  timer.Start();

  cout<<endl<<"  "<<argv[0]<<"    S T A R T E D  ("<<RARFIT_VERSION<<")"
      <<endl<<endl
      <<"Config File: "<<ConfigFile<<endl
      <<"Dataset  Input  Section: \""<<dataInputSec<<"\""<<endl
      <<"mlFitter Config Section: \""<<fitterConfigSec<<"\""<<endl
      <<"mlFitter Action Section: \""<<fitterActionSec<<"\""<<endl
      <<endl<<endl;

  static bool dummy = false;
  if ( dummy ) {
    doBanner();  // force linker to leave it in; will print the RooFit banner
  }
  
  // due to inconsistency between RooRealVar::format() and parser,
  // use printScientific method
  RooRealVar::printScientific(kTRUE);
  
  // set random seed
  if (getenv("RANDOMSEEDBASE")) randomBase=atoi(getenv("RANDOMSEEDBASE"));
  Int_t randomSeed=randomBase+toyID;
  if (randomSeed) {
    cout<<" Set random seed to "<<randomSeed<<endl;
    RooRandom::randomGenerator()->SetSeed(randomSeed);
  }
  
  // read in the datasets (for sig, bkg, MC, Onpeak data, etc.)
  rarDatasets theDatasets(ConfigFile, dataInputSec, fitterActionSec);
  // first set the mlFitter/action section name
  rarConfig::setMasterSec(fitterConfigSec);
  rarConfig::setRunSec(fitterActionSec);
  // instantiate the ML fitter
  rarMLFitter theFitter(ConfigFile, fitterConfigSec,
			"mlFitter MLFitter \"ML Function\"",
			&theDatasets, 0, "mlFitter", "ML Function");
  // set other initial values from user for the fitter
  TString paramDir=".params";
  if (getenv("PARAMDIR")) paramDir=getenv("PARAMDIR");
  theFitter.setParamDir(paramDir);
  TString resultDir="results";
  if (getenv("RESULTDIR")) resultDir=getenv("RESULTDIR");
  theFitter.setResultDir(resultDir);
  theFitter.setToyID(toyID);
  theFitter.setToyNexp(toyNexp);
  theFitter.setToyDir(toyDir);
  // then run it with configs
  theFitter.run();
  
  timer.Stop();
  cout<<endl<<endl
      <<"  "<<argv[0]<<"    F I N I S H E D"
      << " RealTime = " << timer.RealTime() 
      << " secs. CpuTime = " << timer.CpuTime() << " secs." << endl;

  return 0;
  // isn't it simple?
}
