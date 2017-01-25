/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarCompBase.cc,v 1.20 2014/09/14 17:33:49 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides component base class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides component base class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooStringVar.h"

#include "rarCompBase.hh"

ClassImp(rarCompBase)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarCompBase::rarCompBase()
  : rarBasePdf(),
    _nComp(0)
{
  setControlBits("UseBasePdfFit");
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
/// \param useBasePdfFit If set, use #doPdfFit of base class,
///                      otherwise, use its own #doPdfFit.
///
/// The default ctor first initializes data members,
/// #_configFile, #_configSec, #_configStr, name, title,
/// #_nComp, UseBasePdfFit,
/// and then calls #init.
rarCompBase::rarCompBase(const char *configFile, const char *configSec,
			 const char *configStr,
			 rarDatasets *theDatasets, RooDataSet *theData,
			 const char *name, const char *title,
			 Bool_t useBasePdfFit)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _nComp(0)
{
  if (useBasePdfFit) setControlBits("UseBasePdfFit");
  else setControlBits("noUseBasePdfFit");
  init();
}

rarCompBase::~rarCompBase()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It calls #createPdfs to create all its components.
void rarCompBase::init()
{
  cout<<"init of rarCompBase for "<<GetName()<<":"<<endl;
  
  createPdfs("Comps", &_pdfList, &_subPdfs);
  _nComp=_pdfList.GetSize();
  if (_nComp<=0) {
    cout<<"no components defined in RooStringVar \"Comps\" in section \""
	<<_configSec<<"\""<<endl;
    exit(-1);
  }
  
  // add protDataVars
  addProtVars();
  
}

/// \brief Set SimPdf for composite PDF
///
/// \param simPdf The target SimPdf to set
///
/// It sets SimPdf to \p simPdf for itself
/// and also for each its component.
void rarCompBase::setSimPdf(RooSimultaneous *simPdf)
{
  rarBasePdf::setSimPdf(simPdf);
  for (Int_t i=0; i<_nComp; i++) {
    rarBasePdf *thePdf=(rarBasePdf*)_pdfList.At(i);
    thePdf->setSimPdf(simPdf);
  }
}

/// \brief Return RooRealVar ArgSet associated with composite PDF
///
/// It adds #_params from its components onto its own #_params
/// to form ArgSet to return.
RooArgSet rarCompBase::getParams()
{
  RooArgSet paramSet(rarBasePdf::getParams());
  // then add params in each of the comps
  for (Int_t i=0; i<_nComp; i++) {
    rarBasePdf *thePdf=(rarBasePdf*)_pdfList.At(i);
    paramSet.add(thePdf->getParams());
  }
  return paramSet;  
}

/// \brief Actions right after every RooRarFit Pdf is created
/// and before any other action is taken.
///
/// It first calls #preAction of its base class,
/// then it calls #preAction of every component.
void rarCompBase::preAction()
{
  cout<<" In rarCompBase preAction for "<<GetName()<<endl;
  // first run the function in base class
  rarBasePdf::preAction();
  // then run the function in each of the comps
  for (Int_t i=0; i<_nComp; i++) {
    rarBasePdf *thePdf=(rarBasePdf*)_pdfList.At(i);
    thePdf->preAction();
  }
}

/// \brief Constructs ArgSet of params with given names
///
/// \param paramNames The param names or config item name
/// \param useRead To use #readConfStr or not
/// \param fullSet Argset to get additional args
/// \return The ArgSet created
///
/// It first calls rarBasePdf::getArgSet to get ArgSet of its own,
/// then it calls rarBasePdf::getArgSet of each component to get
/// ArgSet from them, and return the total ArgSet.
RooArgSet rarCompBase::getArgSet(TString paramNames, Bool_t useRead,
				 RooArgSet *fullSet)
{
  // first for its own (base class)
  RooArgSet retVal(rarBasePdf::getArgSet(paramNames, useRead, fullSet));
  // then run the function in each of the comps if #useRead is set
  if (useRead) {
    for (Int_t i=0; i<_nComp; i++) {
      rarBasePdf *thePdf=(rarBasePdf*)_pdfList.At(i);
      retVal.add(thePdf->getArgSet(paramNames, useRead, fullSet));
    }
  }
  return retVal;
}

/// \brief RooRealVar ArgSet of correlation coefficient of params
///   associated with composite PDF
/// \return Correlation coefficients
///
/// This function returns RooRealVar ArgSet of correlation coefficient
/// of params associated with this RooRarFitPdf and from its components
RooArgSet rarCompBase::getCorrCoeffs()
{
  RooArgSet theSet(rarBasePdf::getCorrCoeffs());
  // then add correlation coeffs in each of the comps
  for (Int_t i=0; i<_nComp; i++) {
    rarBasePdf *thePdf=(rarBasePdf*)_pdfList.At(i);
    theSet.add(thePdf->getCorrCoeffs());
  }
  return theSet;
}

/// \brief Do pdfFit for given PDFs
///
/// \param pdfList Pdfs need to do pdfFit
///
/// If UseBasePdfFit is set,
/// it calls its base class' rarBasePdf::doPdfFit,
/// if not, it first calls #doXPdfFit for extra pdfs,
/// then it calls rarBasePdf::doPdfFit for each of its components.
void rarCompBase::doPdfFit(TString pdfList)
{
  // first for xtraPfs
  doXPdfFit(pdfList);
  
  cout<<endl<<" In rarCompBase doPdfFit for "<<GetName()<<endl;
  if (!getControlBit("PdfFit")) return;
  if (getControlBit("UseBasePdfFit")) return rarBasePdf::doPdfFit(pdfList);
  // go through all its components
  for (Int_t i=0; i<_nComp; i++) {
    rarBasePdf *thePdf=(rarBasePdf*)_pdfList.At(i);
    thePdf->doPdfFit(pdfList);
  }
  if (getControlBit("UseAlsoBasePdfFit")) return rarBasePdf::doPdfFit(pdfList);
}

/// \brief Do pdfPlot for given PDFs
///
/// \param plotList List of plots
/// \param pdfList List of PDFs to be plotted
/// \return The last RooPlot created
///
/// It first calls #doXPdfPlot to get PDF plots of extra PDFs,
/// then it calls rarBasePdf::doPdfPlot for each of its components
/// to get their PDFs plotted.
RooPlot *rarCompBase::doPdfPlot(TList &plotList, TString pdfList)
{
  RooPlot *frame(0);
  // first for xtraPdfs
  frame=doXPdfPlot(plotList, pdfList);
  
  cout<<endl<<" In rarCompBase doPdfPlot for "<<GetName()<<endl;
  
  //if (!getControlBit("PdfFit")) return frame;
  if (!getControlBit("PdfPlot")) return frame;
  if (getControlBit("PdfPlotDone")) return frame;
  
  // save obs and param
  string obsSaveStr, coeffParamSSaver;
  RooArgSet *theParams(0);
  if (getControlBit("UseBasePdfFit")) {
    writeToStr(_fObsSet, obsSaveStr);
    // restore the fit value for plotting
    theParams=_thePdf->getParameters(_theData);
    writeToStr(*theParams, coeffParamSSaver);
    readFromStr(*theParams, _afterFitSaverStr);
  }
  
  for (Int_t i=0; i<_nComp; i++) {
    rarBasePdf *thePdf=(rarBasePdf*)_pdfList.At(i);
    frame=thePdf->doPdfPlot(plotList, pdfList);
  }
  
  // restore params and obs
  if (getControlBit("UseBasePdfFit")) {// restore params and obs
    readFromStr(*theParams, coeffParamSSaver);
    readFromStr(_fObsSet, obsSaveStr);
  }
  setControlBits("PdfPlotDone");
  
  return frame;
}

/// \brief Get ArgList of PDFs not depending on specified observables for sPlot
/// \param ignoredObs The observables to check
/// \return ArgList of PDFs required
///
/// It loops over its components to get the required PDF list
/// by calling #rarBasePdf::getPdfWOvar.
RooArgList rarCompBase::getPdfsWOvar(RooArgList ignoredObs)
{
  RooArgList thePdfs;
  for (Int_t i=0; i<_nComp; i++) {
    rarBasePdf *thePdf=(rarBasePdf*)_pdfList.At(i);
    thePdfs.add(*thePdf->getPdfWOvar(ignoredObs));
  }
  //thePdfs.Print("v");
  //thePdfs.Print();
  return thePdfs;
}

/// \brief Return prototype var generator for toy study
/// \return The generator
///
/// It constrcuts prototype var generator for toy study
RooAbsPdf *rarCompBase::getProtGen()
{
  if (_theProtGen) return _theProtGen;
  if (_subProtGenPdfs.getSize()<1) {
    for (Int_t i=0; i<_nComp; i++) {
      rarBasePdf *thePdf=(rarBasePdf*)_pdfList.At(i);
      _subProtGenPdfs.add(*thePdf->getProtGen());
    }
  }
  for (Int_t i=0; i<_nComp; i++)
    if (!((rarBasePdf*)_pdfList.At(i))->protGenIsDummy())return _theProtGen;
  // if we go this far, we should have dummyPdf for all subGenPdfs
  setControlBits("dummySubGenPdfs");
  
  return _theProtGen;
}

/// \brief Attach dataset to the pdf
/// \param data The dataset to be attached to the Pdf
///
/// It attaches dataset to this PDF
void rarCompBase::attachDataSet(const RooAbsData &data)
{
  // first call base class'
  rarBasePdf::attachDataSet(data);
  // then call its sub pdfs
  for (Int_t i=0; i<_nComp; i++) {
    rarBasePdf *thePdf=(rarBasePdf*)_pdfList.At(i);
    thePdf->attachDataSet(data);
  }
}

/// \brief Check if the current value of PDF is negative
/// \return true if negative
Bool_t rarCompBase::isNegativeValue()
{

  Bool_t isNegative(kFALSE);
  // first call its base class'
  isNegative=rarBasePdf::isNegativeValue();
  if (isNegative) return isNegative;
  // then call its sub pdfs
  for (Int_t i=0; i<_nComp; i++) {
    rarBasePdf *thePdf=(rarBasePdf*)_pdfList.At(i);
    isNegative=thePdf->isNegativeValue();
    if (isNegative) return isNegative;
  }
  return isNegative;
}
