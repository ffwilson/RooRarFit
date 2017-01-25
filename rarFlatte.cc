/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarFlatte.cc,v 1.4 2014/09/14 17:33:51 fwilson Exp $
 * Authors: Fergus Wilson
 * History:
 * 
 * Copyright (C) 2005-2012, RAL
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides Flatte Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides Flatte Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "RooFlatte.hh"
#include "rarFlatte.hh"

ClassImp(rarFlatte);

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarFlatte::rarFlatte()
  : rarBasePdf(),
    _x(0), _mean(0), _g0(0), _m0a(0), _m0b(0), _g1(0), _m1a(0), _m1b(0)
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
rarFlatte::rarFlatte(const char *configFile, 
				     const char *configSec,
				     const char *configStr,
				     rarDatasets *theDatasets, 
				     RooDataSet *theData,
				     const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _mean(0), _g0(0), _m0a(0), _m0b(0),_g1(0), _m1a(0), _m1b(0)
{
  init();
}

rarFlatte::~rarFlatte()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates the parameters by calling #createAbsReal,
/// and finally it builds the RooFlatte PDF.
void rarFlatte::init()
{

  cout <<"init of rarFlatte for " << GetName() << ":"<<endl;
  
  // first get its observable
  _x = createAbsReal("x", "observable"); assert(_x);

  // Config pdf params
  _mean = createAbsReal("mean", "mean", 0.980, -100, 100);
  _g0  = createAbsReal("g0", "g0", 0.1108, 0, 100);
  _m0a = createAbsReal("m0a", "m0a",  0.13957, 0, 100);
  _m0b = createAbsReal("m0b", "m0b", 0.13957, 0, 100);
  _g1  = createAbsReal("g1", "g1", 0.42, 0, 100);
  _m1a = createAbsReal("m1a", "m1a", 0.49368, 0, 100);
  _m1b = createAbsReal("m1b", "m1b", 0.49368, 0, 100);

  // masses will be made constant even if floated in configuration
  _m0a->setAttribute("Constant", kTRUE);
  _m0b->setAttribute("Constant", kTRUE);
  _m1a->setAttribute("Constant", kTRUE);
  _m1b->setAttribute("Constant", kTRUE);

  _params.Print("v");
  
  // create pdf
  _thePdf=new RooFlatte(Form("the_%s", GetName()),_pdfType+" "+GetTitle(),
			*_x, *_mean, *_g0, *_m0a, *_m0b, *_g1, *_m1a, *_m1b);
}
