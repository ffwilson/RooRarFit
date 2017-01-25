/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarAdd.cc,v 1.12 2014/09/14 17:33:48 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides AddPdf/AddModel class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides AddPdf/AddModel class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"
#include "RooGlobalFunc.h"
using namespace RooFit;

#include "RooAddPdf.h"
#include "RooAddModel.h"

#include "rarAdd.hh"

ClassImp(rarAdd)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarAdd::rarAdd()
  : rarCompBase(),
    _nCoeff(0)
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
/// \param useBasePdfFit If set, use #doPdfFit of base class,
///                      otherwise, use its own #doPdfFit.
/// \param buildAddPdf If set, the RooFit pdf will be built by rarAdd,
///                      otherwise, the pdf is supposed to be created
///                      by sub-class.
///
/// The default ctor first initializes data members,
/// and then calls #init.
rarAdd::rarAdd(const char *configFile, const char *configSec,
	       const char *configStr,
	       rarDatasets *theDatasets, RooDataSet *theData,
	       const char *name, const char *title,
	       Bool_t useBasePdfFit, Bool_t buildAddPdf)
  : rarCompBase(configFile, configSec, configStr,
		theDatasets, theData, name, title, useBasePdfFit),
    _nCoeff(0), _buildAddPdf(buildAddPdf)
{
  init();
}

rarAdd::~rarAdd()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first reads in the info of coeffs,
/// and it prints out the coeffs configured in the config section,
/// then it creates those coeffs by calling #createAbsVar and adds them
/// into the coeff list, #_coeffs,
/// and finally it builds the RooAddPdf/RooAddModel if needed.
void rarAdd::init()
{
  // create coeffs
  createAbsVars("Coeffs", 0, &_coeffs);
  _nCoeff=_coeffs.getSize();
  if (_nCoeff>_nComp) {
    cout<<"You should define at most "<<_nComp
	<<" coeffs in config \"Coeffs\" in section \""
	<<getVarSec()<<"\""<<endl;
    exit(-1);
  }
  if (_nCoeff<_nComp-1) {
    cout<<"You should define at least "<<_nComp-1
	<<" coeffs in config \"Coeffs\" in section \""
	<<getVarSec()<<"\""<<endl;
    exit(-1);
  }
  // make sure AddModel is non-extended
  if (("AddModel"==_pdfType) && (_nCoeff==_nComp)) {
    cout<<"You should define "<<_nComp-1
	<<" coeffs in config \"Coeffs\" in section \""
	<<getVarSec()<<"\""<<endl;
    exit(-1);
  }
  _params.Print("v");
  
  if (_buildAddPdf) {
    // build pdf function
    if ("AddModel"==_pdfType) {
      _thePdf=new RooAddModel(Form("the_%s", GetName()),
			      _pdfType+" "+GetTitle(), _subPdfs, _coeffs);
    } else {
      _thePdf=new RooAddPdf(Form("the_%s", GetName()), _pdfType+" "+GetTitle(),
			    _subPdfs, _coeffs);
    }
  }
  // by default do comp plot
  setControlBit("CompsOnPlot", "compsOnPlot");
  
  cout<<"done init of rarAdd for "<<GetName()<<endl<<endl;
}

/// \brief Do pdfPlot for AddPdf (Overlay component plots)
///
/// \param plotList List of plots
/// \param pdfList List of PDFs to be plotted
/// \return The last RooPlot created
///
/// Use rarBasePdf::doPdfPlot
RooPlot *rarAdd::doPdfPlot(TList &plotList, TString pdfList)
{
  return rarBasePdf::doPdfPlot(plotList, pdfList);
}

/// \brief Get a PDF not depending on specified observables for sPlot
/// \param ignoredObs The observables to check
/// \return The Pdf required
///
/// It checks of the pdf depends on the observables.
/// if no, it returns the default pdf;
/// if yes, it builds the non-dependent pdf.
RooAbsPdf *rarAdd::getPdfWOvar(RooArgList ignoredObs)
{
  RooAbsPdf *thePdf=_thePdf;
  // check if thePdf depends on theVar
  if (!(thePdf->dependsOn(ignoredObs))) return thePdf;
  // create individual Pdf
  TString theName=Form("the_%s_%s", ignoredObs[0].GetName(), GetName());
  TString theTitle=Form("%s w/o %s", GetTitle(), ignoredObs[0].GetName());
  RooArgList myPdfList=getPdfsWOvar(ignoredObs);
  if ("AddModel"==_pdfType) {
    thePdf=new RooAddModel(theName, theTitle, myPdfList, _coeffs);
  } else {
    thePdf=new RooAddPdf(theName, theTitle, myPdfList, _coeffs);
  }
  //thePdf->Print("v");
  //thePdf->Print();
  if (!(thePdf->dependsOn(ignoredObs))) return thePdf;
  delete thePdf;
  // create a dummy for it
  thePdf=new RooGenericPdf(theName, theTitle, "1",
			   *_thePdf->getParameters(_fullObs));
  //thePdf->Print("v");
  //thePdf->Print();
  return thePdf;
}

/// \brief Return prototype var generator for toy study
/// \return The generator
///
/// It constrcuts prototype var generator for toy study
RooAbsPdf *rarAdd::getProtGen()
{
  if (_theProtGen) return _theProtGen;
  // call super class' getProtGen
  rarCompBase::getProtGen();
  _theProtGen=_myDummyPdf;
  if (!getControlBit("dummySubGenPdfs")) {
    if ("AddModel"==_pdfType) {
      _theProtGen=new RooAddModel
	(Form("genAdd_%s", GetName()), Form("GenAdd for %s", GetName()),
	 _subProtGenPdfs, _coeffs);
    } else {
      _theProtGen=new RooAddPdf
	(Form("genAdd_%s", GetName()), Form("GenAdd for %s", GetName()),
	 _subProtGenPdfs, _coeffs);
    }
  }
  if (_protGenPdfs.getSize()>0) {
    RooArgList pdfList(_protGenPdfs);
    pdfList.add(*_theProtGen);
    _theProtGen = new RooProdPdf
      (Form("protGen_%s",GetName()),Form("Prot Gen for %s",GetName()),pdfList);
  }
  
  return _theProtGen;
}
