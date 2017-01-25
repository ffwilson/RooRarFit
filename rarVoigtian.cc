/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarVoigtian.cc,v 1.3 2014/09/14 17:33:56 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides Voigtian Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides Voigtian Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "RooVoigtian.h"

#include "rarVoigtian.hh"

ClassImp(rarVoigtian)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarVoigtian::rarVoigtian()
  : rarBasePdf(),
    _x(0), _mean(0), _width(0), _sigma(0)
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
rarVoigtian::rarVoigtian(const char *configFile, const char *configSec,
			 const char *configStr,
			 rarDatasets *theDatasets, RooDataSet *theData,
			 const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _mean(0), _width(0), _sigma(0)
{
  init();
}

rarVoigtian::~rarVoigtian()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates the parameters by calling #createAbsReal,
/// and finally it builds RooVoigtian PDF
/// with #_pdfType being Voigtian, respectively.
void rarVoigtian::init()
{
  cout<<"init of rarVoigtian for "<<GetName()<<":"<<endl;
  
  // first get its dependent/observable
  _x=createAbsReal("x", "observable"); assert(_x);
  RooRealVar *x=(RooRealVar *)RooArgList(_obsSet).at(0); assert(x);

  // Config pdf params
  _mean=createAbsReal("mean", "#mu", (x->getMin()+x->getMax())/2,
		      x->getMin(), x->getMax(), _x->getUnit());

  _width=createAbsReal("width", "#mu", (x->getMin()+x->getMax())/8,
		      x->getMin(), x->getMax(), _x->getUnit());

  _sigma=createAbsReal("sigma", "#sigma", 0.0, 0., 1., _x->getUnit());

  _params.Print("v");
  
  // create pdf
  Bool_t doFast = kTRUE;
  _thePdf=new RooVoigtian(Form("the_%s", GetName()), _pdfType+" "+GetTitle(),
			  *_x, *_mean, *_width, *_sigma, doFast);
}
