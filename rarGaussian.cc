/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarGaussian.cc,v 1.8 2014/09/14 17:33:51 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides Gaussian/BreitWigner Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides Gaussian/BreitWigner Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "RooBreitWigner.h"
#include "RooGaussian.h"
#include "RooLandau.h"

#include "rarGaussian.hh"

ClassImp(rarGaussian)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarGaussian::rarGaussian()
  : rarBasePdf(),
    _x(0), _mean(0), _sigma(0),
    _scale(0), _shift(0)
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
rarGaussian::rarGaussian(const char *configFile, const char *configSec,
			 const char *configStr,
			 rarDatasets *theDatasets, RooDataSet *theData,
			 const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _mean(0), _sigma(0),
    _scale(0), _shift(0)
{
  init();
}

rarGaussian::~rarGaussian()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates the parameters by calling #createAbsReal,
/// and finally it builds RooGaussian/RooBreitWigner PDF
/// with #_pdfType being Gaussian/BreitWigner, respectively.
void rarGaussian::init()
{
  cout<<"init of rarGaussian for "<<GetName()<<":"<<endl;
  
  // first get its dependent/observable
  _x=createAbsReal("x", "observable"); assert(_x);
  RooRealVar *x=(RooRealVar *)RooArgList(_obsSet).at(0); assert(x);
  // Config pdf params
  _mean=createAbsReal("mean", "#mu", (x->getMin()+x->getMax())/2,
		      x->getMin(), x->getMax(), _x->getUnit());
  if ("notSet"!=readConfStr("shift", "notSet", getVarSec())) { // use shift
    _shift=createAbsReal("shift", "d#mu", 0.);
    _mean=(RooAbsReal*)
      createAbsVar(Form("%s %s", "meanS RooFormulaVar","@0+@1 mean shift"));
  }
  _sigma=createAbsReal("sigma", "#sigma", .1, 0., 1., _x->getUnit());
  if ("notSet"!=readConfStr("scale", "notSet", getVarSec())) { // use scale
    _scale=createAbsReal("scale", "S#sigma_{C}", 1.);
    _sigma=(RooAbsReal*)
      createAbsVar(Form("%s %s","sigmaS RooFormulaVar","@0*@1 sigma scale"));
  }
  _params.Print("v");
  
  // create pdf
  if("BreitWigner"==_pdfType) {
    _thePdf=new RooBreitWigner(Form("the_%s", GetName()),
			       _pdfType+" "+GetTitle(), *_x, *_mean, *_sigma);
  } else if ("Landau"==_pdfType) {
    _thePdf=new RooLandau(Form("the_%s", GetName()),
			  _pdfType+" "+GetTitle(), *_x, *_mean, *_sigma);
  } else { // default
    _thePdf=new RooGaussian(Form("the_%s", GetName()), _pdfType+" "+GetTitle(),
			    *_x, *_mean, *_sigma);
  }
}
