/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarThreshold.cc,v 1.2 2014/09/14 17:33:55 fwilson Exp $
 * Authors: Fergus Wilson
 * History:
 * 
 * Copyright (C) 2005-2012, RAL
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides Threshold Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides Threshold Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"
//#include "RooNumber.h"

#include "RooThreshold.hh"
#include "rarThreshold.hh"

using namespace RooFit;

ClassImp(rarThreshold);

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarThreshold::rarThreshold()
  : rarBasePdf(),
    _x(0), _m0(0), _power(0), _nOrder(1)
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
rarThreshold::rarThreshold(const char *configFile, 
				     const char *configSec,
				     const char *configStr,
				     rarDatasets *theDatasets, 
				     RooDataSet *theData,
				     const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _m0(0), _power(0), _nOrder(1)
{
  init();
}

rarThreshold::~rarThreshold()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates the parameters by calling #createAbsReal,
/// and finally it builds the RooThreshold PDF.
void rarThreshold::init()
{

  cout <<"init of rarThreshold for " << GetName() << ":"<<endl;
  
  // first get its observable
  _x = createAbsReal("x", "observable"); assert(_x);

  // Config pdf params
  _m0 = createAbsReal("m0", "threshold", 0, -100000, 100000);
  _power = createAbsReal("power", "power", 1, -1000, 1000);

  // read in pdf params
  _nOrder=atoi(readConfStr("nOrder", "1", getVarSec()));
  if (_nOrder<=0) {_nOrder=1;}
  _nOrder++;
  for (Int_t i=1; i<_nOrder; i++) {
    RooAbsReal *P=createAbsReal(Form("P%02d",i),Form("p_{%d}", i),
				0, -10000, +10000);
    _coeffs.add(*P);
  }

  _params.Print("v");
  
  // create pdf
  _thePdf=new RooThreshold(Form("the_%s", GetName()),_pdfType+" "+GetTitle(),
			   *_x, *_m0, *_power, _coeffs);
}
