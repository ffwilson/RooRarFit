/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarBinned.cc,v 1.5 2014/09/14 17:33:49 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides Binned Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides Binned Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "RooBinnedPdf.hh"

#include "rarBinned.hh"

ClassImp(rarBinned)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarBinned::rarBinned()
  : rarBasePdf(),
    _x(0), _nBins(0), _limits(0)
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
rarBinned::rarBinned(const char *configFile, const char *configSec,
		 const char *configStr, rarDatasets *theDatasets,
		 RooDataSet *theData, const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _nBins(0), _limits(0)
{
  init();
}

rarBinned::~rarBinned()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first reads in #_nBins and #_limits info from its param config section,
/// and creates \p H00 ... \p H\<nBins-1\> parameters
/// by calling #createAbsReal,
/// and finally it builds \p Binned PDF.
void rarBinned::init()
{
  cout<<"init of rarBinned for "<<GetName()<<":"<<endl;
  
  // first get its dependent/observable
  _x=createAbsReal("x", "observable"); assert(_x);
  RooRealVar *x=(RooRealVar *)RooArgList(_obsSet).at(0); assert(x);
  // Config pdf params
  _nBins=atoi(readConfStr("nBins", "0", getVarSec()));
  if (_nBins<=0) _nBins= // from plotBins
    atoi(readConfStr(Form("plotBins_%s", x->GetName()),"0", getVarSec()));
  if (_nBins<=0) _nBins=x->getBins(); // get bins from its obs
  _limits= new TArrayD(_nBins+1);
  // read in limit info
  TString limitStr=readConfStr("limits", "", getVarSec());
  rarStrParser limitStrParser=limitStr;
  Int_t nLimits=limitStrParser.nArgs();
  if (_nBins+1==nLimits) {
    for (Int_t i=0; i<_nBins+1; i++) {
      (*_limits)[i]=atof(limitStrParser[i]);
    }
  } else { // set equal interval
    Double_t binned=(x->getMax() - x->getMin())/_nBins;
    for (Int_t i=0; i<_nBins+1; i++) {
      (*_limits)[i]=x->getMin()+i*binned;
    }
  }
  // set lo and hi limits
  //(*_limits)[0]=x->getMin();
  //(*_limits)[_nBins]=x->getMax();
  // read in coeffs
  for (Int_t i=0; i<_nBins-1; i++) {
    RooAbsReal *H=createAbsReal(Form("H%02d",i),Form("h_{%d}", i),.1,0,1000);
    _coeffs.add(*H);
  }
  _params.Print("v");
  
  // create pdf
  _thePdf=
    new RooBinnedPdf(Form("the_%s", GetName()),
		     _pdfType+" "+GetTitle(),
		     *_x, _coeffs, *_limits);
}
