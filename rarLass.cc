/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarLass.cc,v 1.5 2014/09/14 17:33:52 fwilson Exp $
 * Authors: F. Wilson
 * History:
 * 
 * Copyright (C) 2005-2012, RAL
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides parameterisation of LASS lineshape Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides parameterisation of LASS lineshape Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

// Include the RooFit Pdf
#include "RooLass.hh"

// Include this classes header
#include "rarLass.hh"

ClassImp(rarLass)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarLass::rarLass()
  : rarBasePdf(),
    _x(0), _mean(0), _width(0), _effRange(0), _scatlen(0), _turnOffVal(0)
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
rarLass::rarLass(const char *configFile, const char *configSec,
		     const char *configStr,
		     rarDatasets *theDatasets, RooDataSet *theData,
		     const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _mean(0), _width(0), _effRange(0), _scatlen(0), _turnOffVal(0)
{
  init();
}

rarLass::~rarLass()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates the parameters by calling #createAbsReal,
/// and finally it builds the LASS  PDF.
void rarLass::init()
{
  cout<<"init of rarLass for "<<GetName()<<":"<<endl;
  
  _x=createAbsReal("x", "observable"); assert(_x);
  RooRealVar *x=(RooRealVar *)RooArgList(_obsSet).at(0); assert(x);
 
  _mean=createAbsReal("mean", "#mu", (x->getMin()+x->getMax())/2,
                      x->getMin(), x->getMax(), _x->getUnit());
  // initialise with values from Bill's 34pt fit
  // http://www.slac.stanford.edu/~wmd/kpi_swave/kpi_swave_fit.note
  // They should replace those in Nucl Phys B296, 493 (1988)
  _width      = createAbsReal("width", "#Lambda", 0.279, 0, 1000);
  _effRange   = createAbsReal("effRange", "R", 1.76, 0, 1000);
  _scatlen    = createAbsReal("scatlen", "s", 1.95, 0, 1000);
  // LASS not tested beyond about 1.65 GeV
  _turnOffVal = createAbsReal("turnOffVal", "t", 1.825, 0, 1000);
  _params.Print("v");

  _thePdf = new RooLass(Form("the_%s", GetName()),_pdfType+" "+GetTitle(),
			       *_x, *_mean, *_width, *_effRange, *_scatlen, *_turnOffVal);

}
