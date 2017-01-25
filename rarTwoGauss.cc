/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarTwoGauss.cc,v 1.8 2014/09/14 17:33:55 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides TwoGauss Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides TwoGauss Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "RooGaussian.h"

#include "rarTwoGauss.hh"

ClassImp(rarTwoGauss)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarTwoGauss::rarTwoGauss()
  : rarBasePdf(),
    _x(0), _meanC(0), _sigmaC(0), _meanT(0), _sigmaT(0), _fracC(0),
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
rarTwoGauss::rarTwoGauss(const char *configFile, const char *configSec,
			 const char *configStr, rarDatasets *theDatasets,
			 RooDataSet *theData, const char*name,const char*title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _meanC(0), _sigmaC(0), _meanT(0), _sigmaT(0), _fracC(0),
    _scale(0), _shift(0)
{
  init();
}

rarTwoGauss::~rarTwoGauss()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates the parameters by calling #createAbsReal,
/// and finally it builds Double-Gaussian PDF by
/// creating two RooGaussian PDFs and combining them
/// using RooAddPdf.
void rarTwoGauss::init()
{
  cout<<"init of rarTwoGauss for "<<GetName()<<":"<<endl;
  
  // first get its dependent/observable
  _x=createAbsReal("x", "observable"); assert(_x);
  RooRealVar *x=(RooRealVar *)RooArgList(_obsSet).at(0); assert(x);
  // Config pdf params
  _meanC=createAbsReal("meanC", "#mu_{C}",
		       (x->getMin()+x->getMax())/2,
		       x->getMin(), x->getMax(), _x->getUnit());
  _meanT=createAbsReal("meanT", "#mu_{T}",
		       (x->getMin()+x->getMax())/2,
		       x->getMin(), x->getMax(), _x->getUnit());
  if ("notSet"!=readConfStr("shift", "notSet", getVarSec())) { // use shift
    _shift=createAbsReal("shift", "d#mu", 0.);
    _meanC=(RooAbsReal*)
      createAbsVar(Form("%s %s", "meanCS RooFormulaVar","@0+@1 meanC shift"));
    _meanT=(RooAbsReal*)
      createAbsVar(Form("%s %s", "meanTS RooFormulaVar","@0+@1 meanT shift"));
  }
  _sigmaC=createAbsReal("sigmaC", "#sigma_{C}", .1, 0., 1., _x->getUnit());
  if ("notSet"!=readConfStr("scale", "notSet", getVarSec())) { // use scale
    _scale=createAbsReal("scale", "S#sigma_{C}", 1.);
    _sigmaC=(RooAbsReal*)
      createAbsVar(Form("%s %s","sigmaCS RooFormulaVar","@0*@1 sigmaC scale"));
  }
  _sigmaT=createAbsReal("sigmaT", "#sigma_{T}", .1, 0., 1., _x->getUnit());
  _fracC=createAbsReal("fracC", "f_{C}", .6, 0., 1.);
  _params.Print("v");
  
  // create pdf
  RooAbsPdf *corePdf=
    new RooGaussian(Form("core_%s",GetName()),
		    Form("Core Gaussian %s", GetTitle()),
		    *_x, *_meanC, *_sigmaC);
  _subPdfs.add(*corePdf);
  RooAbsPdf *tailPdf=
    new RooGaussian(Form("tail_%s",GetName()),
		    Form("Tail Gaussian %s", GetTitle()),
		    *_x, *_meanT, *_sigmaT);
  _subPdfs.add(*tailPdf);
  _thePdf=new RooAddPdf(Form("the_%s", GetName()), _pdfType+" "+GetTitle(),
			*corePdf, *tailPdf, *_fracC);
  // by default do comp plot
  setControlBit("CompsOnPlot", "compsOnPlot");
}
