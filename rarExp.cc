/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarExp.cc,v 1.6 2014/09/14 17:33:51 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides Exponential Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides Exponential Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "RooExponential.h"

#include "rarExp.hh"

ClassImp(rarExp)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarExp::rarExp()
  : rarBasePdf(),
    _x(0), _c(0)
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
rarExp::rarExp(const char *configFile, const char *configSec,
	       const char *configStr,
	       rarDatasets *theDatasets, RooDataSet *theData,
	       const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _c(0)
{
  init();
}

rarExp::~rarExp()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates the parameters by calling #createAbsReal,
/// and finally it builds the RooExponential PDF.
void rarExp::init()
{
  cout<<"init of rarExp for "<<GetName()<<":"<<endl;
  
  // first get its obs
  _x=createAbsReal("x", "observable"); assert(_x);
  // Config pdf params
  _c=createAbsReal("c", "c", 0, -10, 10);
  _params.Print("v");
  
  // create pdf
  _thePdf=new RooExponential(Form("the_%s", GetName()),_pdfType+" "+GetTitle(),
			     *_x, *_c);
}
