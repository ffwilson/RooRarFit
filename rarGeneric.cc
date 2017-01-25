/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarGeneric.cc,v 1.7 2014/09/14 17:33:51 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides Generic Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides Generic Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "rarGeneric.hh"

ClassImp(rarGeneric)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarGeneric::rarGeneric()
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
rarGeneric::rarGeneric(const char *configFile, const char *configSec,
		       const char *configStr,
		       rarDatasets *theDatasets, RooDataSet *theData,
		       const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title)
{
  init();
}

rarGeneric::~rarGeneric()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first reads in the formula info from config item \p formula,
/// gets the first token as formula string,
/// and use #getFormulaArgs to get ArgList of the PDF,
/// and finally it builds RooGenericPdf.
void rarGeneric::init()
{
  cout<<"init of rarGeneric for "<<GetName()<<":"<<endl;
  
  // read the formula string
  TString formulaStr=readConfStr("formula", "1", getVarSec());
  rarStrParser formulaStrParser=formulaStr;
  // get formula
  TString formula=formulaStrParser[0];
  formulaStrParser.Remove();
  // get formula dep list
  RooArgList *depVarList=getFormulaArgs(formulaStrParser);
  _params.Print("v");
  cout<<"formula string:\t"<<formulaStr<<endl;
  
  // create the generic pdf
  _thePdf=new RooGenericPdf(Form("the_%s", GetName()), _pdfType+" "+GetTitle(),
			    formula, *depVarList);
}
