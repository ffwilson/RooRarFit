/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarSimPdf.cc,v 1.18 2014/09/14 17:33:54 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides Simultaneous Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides Simultaneous Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"
//#include <sstream>

#include "RooAbsPdf.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooSimPdfBuilder.h"
#include "RooStringVar.h"

#include "rarMLFitter.hh"
#include "rarSimPdf.hh"

ClassImp(rarSimPdf)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarSimPdf::rarSimPdf()
  : rarCompBase(),
    _simBuilder(0), _simConfig(0)
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
/// #_configFile, #_configSec, #_configStr, name, title,
/// and then calls #init.
rarSimPdf::rarSimPdf(const char *configFile, const char *configSec,
		     const char *configStr,
		     rarDatasets *theDatasets, RooDataSet *theData,
		     const char *name, const char *title)
  : rarCompBase(configFile, configSec, configStr,
		theDatasets, theData, name, title, kFALSE),
    _simBuilder(0), _simConfig(0)
{
  init();
}

rarSimPdf::~rarSimPdf()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// The initialization work to read in component config info
/// is done with rarCompBase::init.
/// \p init here then reads in category info for SimPdf,
/// then it builds datasets for each category type
/// so that each SimPdf component built has its own default fit dataset,
/// and finally it builds RooSimultaneous PDF with configured components.
void rarSimPdf::init()
{
  // it should always have fitData
  if (!_theData) {
    cout<<" No dataset for Sim Pdf"<<endl;
    exit(-1);
  }
  // create SimPdfBuilder
  _simBuilder=new RooSimPdfBuilder(_subPdfs);
  _simConfig=_simBuilder->createProtoBuildConfig();
  _simConfig->readFromFile(_configFile, 0, getVarSec());
  cout<<"simPdfBuilder configs in the config file:"<<endl;
  _simConfig->Print("v");
  // do we need to use _simBuilder->addSpecializations?
  // Yes, let's do it
  if (getFitter()->getSimBuilder()) {
    RooSimPdfBuilder *theFitterBuilder=getFitter()->getSimBuilder();
    _simBuilder->addSpecializations(theFitterBuilder->splitLeafList());
  }
  // build SimPdf
  _thePdf=(RooAbsPdf*)_simBuilder->buildPdf(*_simConfig, _theData, 0, kTRUE);
  _thePdf->SetName(Form("the_%s", GetName()));
  // fit the total pdf for pdfFit as well
  setControlBits("UseAlsoBasePdfFit");
  
  //_thePdf->Print();
  //_thePdf->Print("v");
  //_thePdf->getComponents()->Print("v");
  cout<<"done init of rarSimPdf for "<<GetName()<<endl<<endl;
}

/// \brief Get a PDF not depending on specified observables for sPlot
/// \param ignoredObs The observables to check
/// \return The Pdf required
///
/// It checks of the pdf depends on the observables.
/// if no, it returns the default pdf;
/// if yes, it builds the non-dependent pdf.
RooAbsPdf *rarSimPdf::getPdfWOvar(RooArgList ignoredObs)
{
  RooAbsPdf *thePdf=_thePdf;
  cout<<"This pdf can not be used in model"<<endl;
  exit(-1);
  ignoredObs.getHashTableSize() ; // avoid "unused" warning
  return thePdf;
}

/// \brief Return prototype var generator for toy study
/// \return The generator
///
/// It constrcuts prototype var generator for toy study
RooAbsPdf *rarSimPdf::getProtGen()
{
  cout<<"This pdf can not be used in model"<<endl;
  exit(-1);
  return _theProtGen;
}
