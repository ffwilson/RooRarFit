/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarTriGauss.cc,v 1.8 2014/09/14 17:33:55 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides TriGauss related Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides TriGauss related Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "RooAddModel.h"
#include "RooDecay.h"
#include "RooGaussian.h"
#include "RooGaussModel.h"

#include "rarTriGauss.hh"

ClassImp(rarTriGauss)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarTriGauss::rarTriGauss()
  : rarBasePdf(),
    _x(0), _tau(0), _msSF(0),
    _meanC(0), _sigmaC(0), _meanT(0), _sigmaT(0), _meanO(0), _sigmaO(0),
    _fracC(0), _fracO(0)
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
rarTriGauss::rarTriGauss(const char *configFile, const char *configSec,
			 const char *configStr, rarDatasets *theDatasets,
			 RooDataSet *theData, const char*name,const char*title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _tau(0), _msSF(0),
    _meanC(0), _sigmaC(0), _meanT(0), _sigmaT(0), _meanO(0), _sigmaO(0),
    _fracC(0), _fracO(0)
{
  init();
}

rarTriGauss::~rarTriGauss()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates the parameters by calling #createAbsReal,
/// and finally it builds RooTripleGaussian PDF/Model
/// based on #_pdfType.
void rarTriGauss::init()
{
  cout<<"init of rarTriGauss for "<<GetName()<<":"<<endl;
  
  // first get its dependent/observable
  _x=createAbsReal("x", "observable"); assert(_x);
  RooRealVar *x=(RooRealVar *)RooArgList(_obsSet).at(0); assert(x);
  TString msSFStr=readConfStr("msSF", "notSet", getVarSec());
  if (_x!=x) {
    // make sure theDep is not derived for most of pdfTypes
    if ("TriGauss"!=_pdfType) {
      cout <<"No derived dependent allowed for "<<_pdfType<<endl;
      exit(-1);
    } else if ("notSet"!=msSFStr) {
      cout <<"No derived dependent allowed for TriGauss with scaling"<<endl;
      exit(-1);
    }
  }
  // Config pdf params
  if (("notSet"!=msSFStr)||("TriGauss"!=_pdfType))
    _msSF=createAbsReal("msSF", "SF_{#mu#sigma}", 1.);
  if ("GexpShape"==_pdfType) _tau=createAbsReal("tau","#tau", 1.5, 0,10, "ps");
  _meanC=createAbsReal("meanC", "#mu_{C}",
		       (x->getMin()+x->getMax())/2,
		       x->getMin(), x->getMax(), _x->getUnit());
  _sigmaC=createAbsReal("sigmaC", "#sigma_{C}", .05, 0.,1., _x->getUnit());
  _meanT=createAbsReal("meanT", "#mu_{T}",
		       (x->getMin()+x->getMax())/2,
		       x->getMin(), x->getMax(), _x->getUnit());
  _sigmaT=createAbsReal("sigmaT", "#sigma_{T}", .09, 0.,1., _x->getUnit());
  _meanO=createAbsReal("meanO", "#mu_{O}",
		       (x->getMin()+x->getMax())/2,
		       x->getMin(), x->getMax(), _x->getUnit());
  _sigmaO=createAbsReal("sigmaO", "#sigma_{O}", .20, 0.,1., _x->getUnit());
  _fracC=createAbsReal("fracC", "f_{C}", .80, 0., 1.);
  _fracO=createAbsReal("fracO", "f_{O}", .05, 0., 1.);
  
  _params.Print("v");
  
  // create pdf
  RooAbsPdf *corePdf(0), *tailPdf(0), *outlPdf(0);
  if (("TriGauss"==_pdfType)&&("notSet"==msSFStr)) {
    // core, tail, and outliner
    corePdf=new RooGaussian(Form("core_%s",GetName()),
			    Form("Core Gaussian %s", GetTitle()),
			    *_x, *_meanC, *_sigmaC);
    tailPdf=new RooGaussian(Form("tail_%s",GetName()),
			    Form("Tail Gaussian %s", GetTitle()),
			    *_x, *_meanT, *_sigmaT);
    outlPdf=new RooGaussian(Form("outl_%s",GetName()),
			    Form("Outl Gaussian %s", GetTitle()),
			    *_x, *_meanO, *_sigmaO);
    _thePdf=new
      RooAddPdf(Form("the_%s", GetName()), _pdfType+" "+GetTitle(),
		RooArgList(*corePdf, *outlPdf, *tailPdf),
		RooArgList(*_fracC, *_fracO));
  } else {
    // FlatScaleFactorIntegral?
    Bool_t FSFIntegral(kTRUE);
    if ("no"==readConfStr("FlatScaleFactorIntegral", "yes", getVarSec()))
      FSFIntegral=kFALSE;
    // core, tail, and outliner
    corePdf=new RooGaussModel(Form("core_%s",GetName()),
			      Form("Core Gaussian %s", GetTitle()),
			      *x, *_meanC, *_sigmaC, *_msSF);
    ((RooGaussModel*)corePdf)->advertiseFlatScaleFactorIntegral(FSFIntegral);
    tailPdf=new RooGaussModel(Form("tail_%s",GetName()),
			      Form("Tail Gaussian %s", GetTitle()),
			      *x, *_meanT, *_sigmaT, *_msSF);
    ((RooGaussModel*)tailPdf)->advertiseFlatScaleFactorIntegral(FSFIntegral);
    outlPdf=new RooGaussModel(Form("outl_%s",GetName()),
			      Form("Outl Gaussian %s", GetTitle()),
			      *x, *_meanO, *_sigmaO);
    ((RooGaussModel*)outlPdf)->advertiseFlatScaleFactorIntegral(FSFIntegral);
    
    if ("TriGauss"==_pdfType) {
      _thePdf=new
	RooAddPdf(Form("the_%s", GetName()), _pdfType+" "+GetTitle(),
		  RooArgList(*corePdf, *outlPdf, *tailPdf),
		  RooArgList(*_fracC, *_fracO));
    } else if ("TriGaussModel"==_pdfType) {
      _thePdf=new
	RooAddModel(Form("the_%s", GetName()), _pdfType+" "+GetTitle(),
		    RooArgList(*corePdf, *outlPdf, *tailPdf),
		    RooArgList(*_fracC, *_fracO));
    } else if ("GexpShape"==_pdfType) {
      // parse the decay type
      RooDecay::DecayType decayType=RooDecay::DoubleSided;
      TString decayTypeStr=readConfStr("decayType","DoubleSided",getVarSec());
      if ("SingleSided"==decayTypeStr) {
	decayType=RooDecay::SingleSided;
      } else if ("DoubleSided"==decayTypeStr) {
	decayType=RooDecay::DoubleSided;
      } else if ("Flipped"==decayTypeStr) {
	decayType=RooDecay::Flipped;
      }
      RooAbsPdf *expDecay=new RooDecay(Form("decay_%s",GetName()),
				       Form("decay %s", GetTitle()),
				       *x, *_tau,
				       *((RooGaussModel*)corePdf), decayType);
      _subPdfs.add(*expDecay);
      _thePdf=new
	RooAddPdf(Form("the_%s", GetName()), _pdfType+" "+GetTitle(),
		  RooArgList(*expDecay, *outlPdf, *tailPdf),
		  RooArgList(*_fracC, *_fracO));
    }
  }
  if ("GexpShape"!=_pdfType) _subPdfs.add(*corePdf);
  _subPdfs.add(*tailPdf);
  _subPdfs.add(*outlPdf);
  // by default do comp plot
  setControlBit("CompsOnPlot", "compsOnPlot");
}
