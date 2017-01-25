/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarBifurGauss.cc,v 1.7 2014/09/14 17:33:49 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides BifurGauss Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides BifurGauss Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooFormulaVar.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "RooBifurGauss.h"

#include "rarBifurGauss.hh"

ClassImp(rarBifurGauss)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarBifurGauss::rarBifurGauss()
  : rarBasePdf(),
    _x(0), _parSymLevel(0), _peak(0), _sigL(0), _sigR(0)
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
rarBifurGauss::rarBifurGauss(const char *configFile,
			     const char *configSec,
			     const char *configStr,
			     rarDatasets *theDatasets,
			     RooDataSet *theData,
			     const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _parSymLevel(0), _peak(0), _sigL(0), _sigR(0)
{
  init();
}

rarBifurGauss::~rarBifurGauss()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates the parameters by calling #createAbsReal,
/// and finally it builds the RooBifurGauss PDF.
void rarBifurGauss::init()
{
  cout<<"init of rarBifurGauss for "<<GetName()<<":"<<endl;
  
  // first get its dependent/observable
  _x=createAbsReal("x", "observable"); assert(_x);
  RooRealVar *x=(RooRealVar *)RooArgList(_obsSet).at(0); assert(x);
  // default _parSymLevel for BGGauss pdfType
  if ("BGGauss"==_pdfType) _parSymLevel=1;
  // read in _parSymLevel from config section
  _parSymLevel=atoi(readConfStr("parSymLevel", Form("%d", _parSymLevel),
				getVarSec()));
  if (_parSymLevel<0) _parSymLevel=0;
  if (_parSymLevel>3) _parSymLevel=3;
  // Config pdf params
  if (0==_parSymLevel) {
    _peak=createAbsReal("peak", "#mu", (x->getMin()+x->getMax())/2,
			x->getMin(), x->getMax(), _x->getUnit());
    _sigL=createAbsReal("sigL", "#sigma_{L}", .1, 0.,
			(x->getMax()-x->getMin())/2,_x->getUnit());
    _sigR=createAbsReal("sigR", "#sigma_{R}", .1, 0.,
			(x->getMax()-x->getMin())/2,_x->getUnit());
  } else {
    createAbsReal("mean", "#mu", (x->getMin()+x->getMax())/2,
		  x->getMin(), x->getMax(), _x->getUnit());
    createAbsReal("rms", "#sigma", .1, 0.,(x->getMax()-x->getMin())/2,
		  _x->getUnit());
    createAbsReal("asym", "A", 0., -1., 1.);
  }
  if (1==_parSymLevel) {
    createAbsReal("Cmean", "sqrt(8/pi)", sqrt(8./3.14159265), "RooConstVar");
    _peak=(RooAbsReal*)
      createAbsVar("peak RooFormulaVar @0-@3*@1*@2 mean rms asym Cmean");
    _sigL=(RooAbsReal*)createAbsVar("sigL RooFormulaVar @0*(1-@1) rms asym");
    _sigR=(RooAbsReal*)createAbsVar("sigR RooFormulaVar @0*(1+@1) rms asym");
  } else if (2==_parSymLevel) {
    createAbsReal("Cmean", "sqrt(8/pi)", sqrt(8./3.14159265), "RooConstVar");
    createAbsReal("CA2", "3/2 - 4/pi", .5*(3. - 8./3.14159265), "RooConstVar");
    _peak=(RooAbsReal*)
      createAbsVar(Form("%s %s", "peak RooFormulaVar",
			"@0-@3*@1*@2*(1-@4*@2*@2) mean rms asym Cmean CA2"));
    _sigL=(RooAbsReal*)
      createAbsVar("sigL RooFormulaVar @0*(1-@2*@1*@1)*(1-@1) rms asym CA2");
    _sigR=(RooAbsReal*)
      createAbsVar("sigR RooFormulaVar @0*(1-@2*@1*@1)*(1+@1) rms asym CA2");
  } else if (3==_parSymLevel) {
    createAbsReal("CA2x3", "(3*pi-8)/8", (3.*3.14159265-8.)/8., "RooConstVar");
    createAbsReal("CAx3", "sqrt(pi/8)", sqrt(3.14159265/8.), "RooConstVar");
    _peak=(RooAbsReal*)
      createAbsVar(Form("%s %s %s", "peak RooFormulaVar",
			"@0-@2/@1/@1/sqrt(1+@3*(@2/@1/@1/@1)*(@2/@1/@1/@1))",
			"mean rms asym CA2x3"));
    _sigL=(RooAbsReal*)
      createAbsVar(Form("%s %s%s %s", "sigL RooFormulaVar",
			"@0*(1-@2*@1/@0/@0/@0)/sqrt(1+@3*(@1/@0/@0/@0)",
			"*(@1/@0/@0/@0))", "rms asym CAx3 CA2x3"));
    _sigR=(RooAbsReal*)
      createAbsVar(Form("%s %s%s %s", "sigR RooFormulaVar",
			"@0*(1+@2*@1/@0/@0/@0)/sqrt(1+@3*(@1/@0/@0/@0)",
			"*(@1/@0/@0/@0))", "rms asym CAx3 CA2x3"));
  }
  _params.Print("v");
  
  // create pdf
  _thePdf=new RooBifurGauss(Form("the_%s", GetName()), _pdfType+" "+GetTitle(),
			    *_x, *_peak, *_sigL, *_sigR);
}
