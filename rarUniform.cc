/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarUniform.cc,v 1.2 2014/09/14 17:33:56 fwilson Exp $
 * Authors: F. wilson
 * History:
 * 
 * Copyright (C) 2005-2012, RAL/STFC
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides Polynomial/Chebychev Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides Polynomial/Chebychev Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooFormulaVar.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

// not available in analysis-52 version of babar root
#ifdef USENEWROOT
#include "RooUniform.h"
#endif

#include "rarUniform.hh"

ClassImp(rarUniform)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarUniform::rarUniform()
  : rarBasePdf()
{
  init();
}

/// \brief Default ctor
///
/// \param configFile The config file
/// \param configSec The config section
/// \param configStr The config string
/// \param theDatasets Available datasets
/// \param theData Default dataset for this PDF
/// \param name The name
/// \param title The title
///
/// The default ctor first initializes data members,
/// and then calls #init.
rarUniform::rarUniform(const char *configFile, const char *configSec,
		 const char *configStr, rarDatasets *theDatasets,
		 RooDataSet *theData, const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title)
{
  init();
}

rarUniform::~rarUniform()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates the parameters by calling #createAbsReal,
/// and finally it builds RooUniformnomial/RooChebychev PDF
/// with #_pdfType being Uniformnomial/Chebychev, respectively.
void rarUniform::init()
{
  cout<<"init of rarUniform for "<<GetName()<<":"<<endl;
  
   // read the obs string
  rarStrParser obsStrParser=readConfStr("obs", "", getVarSec());

  // get obs list
  getFormulaArgs(obsStrParser);

  // print out the observables
  cout<<" Obs in pdf: "<<endl;
  _obsSet.Print("v");
 
  // create pdf
#ifdef USENEWROOT
  _thePdf=new RooUniform(Form("the_%s",GetName()),_pdfType+" "+GetTitle(), _obsSet);
#else
  cout << "Uniform not implemented in old root versions." << endl;
  exit(-1);
#endif
}
