/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarGaussModel.cc,v 1.7 2014/09/14 17:33:51 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides Gaussian Resolution Model class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides Gaussian Resolution Model class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "RooGaussModel.h"

#include "rarGaussModel.hh"

ClassImp(rarGaussModel)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarGaussModel::rarGaussModel()
  : rarBasePdf(),
    _x(0), _mean(0), _sigma(0), _msSF(0), _meanSF(0), _sigmaSF(0)
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
rarGaussModel::rarGaussModel(const char *configFile, const char *configSec,
			     const char *configStr,
			     rarDatasets *theDatasets, RooDataSet *theData,
			     const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _mean(0), _sigma(0), _msSF(0), _meanSF(0), _sigmaSF(0)
{
  init();
}

rarGaussModel::~rarGaussModel()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates the parameters by calling #createAbsReal,
/// and finally it builds RooGaussModel PDF.
void rarGaussModel::init()
{
  cout<<"init of rarGaussModel for "<<GetName()<<":"<<endl;
  
  // first get its dependent/observable
  _x=(RooRealVar*)createAbsReal("x", "observable"); assert(_x);
  RooRealVar *x=(RooRealVar *)RooArgList(_obsSet).at(0); assert(x);
  // make sure theDep is not derived
  if (_x!=x) {
    cout <<"No derived dependent allowed in rarGaussModel"<<endl;
    exit(-1);
  }
  // Config pdf params
  _mean=createAbsReal("mean", "#mu", (x->getMin()+x->getMax())/2,
		      x->getMin(), x->getMax(), _x->getUnit());
  _sigma=createAbsReal("sigma", "#sigma", .1, 0., 1., _x->getUnit());
  // param Level
  Int_t paramLevel(0);
  TString msSFStr=readConfStr("msSF", "notSet", getVarSec());
  TString meanSFStr=readConfStr("meanSF", "notSet", getVarSec());
  TString sigmaSFStr=readConfStr("sigmaSF", "notSet", getVarSec());
  if("notSet"!=msSFStr) { // level 1
    paramLevel=1;
    _msSF=createAbsReal("msSF", "SF_{#mu#sigma}", 1.);
  } else if (("notSet"!=meanSFStr) && ("notSet"!=sigmaSFStr)) { // level 2
    paramLevel=2;
    _meanSF=createAbsReal("meanSF", "SF_{#mu}", 1.);
    _sigmaSF=createAbsReal("sigmaSF", "SF_{#sigma}", 1.);
  }
  _params.Print("v");
  
  // create pdf
  if (1==paramLevel) {
    _thePdf=new RooGaussModel(Form("the_%s",GetName()),_pdfType+" "+GetTitle(),
			      *_x, *_mean, *_sigma, *_msSF);
  } else if (2==paramLevel) {
    _thePdf=new RooGaussModel(Form("the_%s",GetName()),_pdfType+" "+GetTitle(),
			      *_x, *_mean, *_sigma, *_meanSF, *_sigmaSF);
  } else { //default level 0
    _thePdf=new RooGaussModel(Form("the_%s",GetName()),_pdfType+" "+GetTitle(),
			      *_x, *_mean, *_sigma);
  }
  // FlatScaleFactorIntegral?
  Bool_t FSFIntegral(kTRUE);
  if ("no"==readConfStr("FlatScaleFactorIntegral", "yes", getVarSec()))
    FSFIntegral=kFALSE;
  ((RooGaussModel*)_thePdf)->advertiseFlatScaleFactorIntegral(FSFIntegral);
  
}
