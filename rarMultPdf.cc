/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarMultPdf.cc,v 1.6 2014/09/14 17:33:53 fwilson Exp $
 * Authors: Lei Zhang, Wolfgang Gradl
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides MultPdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides MultPdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"
#include <string>

#include "RooArgList.h"
#include "RooDataSet.h"

#include "RooRealVar.h"
#include "RooStringVar.h"
#include "RooGlobalFunc.h"
using namespace RooFit;

#include "rarMultPdf.hh"

ClassImp(rarMultPdf)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarMultPdf::rarMultPdf()
  : rarCompBase()
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
rarMultPdf::rarMultPdf(const char *configFile, const char *configSec,
		       const char *configStr,
		       rarDatasets *theDatasets, RooDataSet *theData,
		       const char *name, const char *title)
  : rarCompBase(configFile, configSec, configStr,
		theDatasets, theData, name, title, kTRUE)
{
  init();
}

rarMultPdf::~rarMultPdf()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first reads in the info of comps,
/// and it prints out the components configured in the config section,
/// then it creates those component PDFs and multiplies them
/// to build the RooMultPdf.
void rarMultPdf::init()
{
  
  TString format = "@0";
  for(Int_t i = 1; i < _subPdfs.getSize(); ++i) {
    format += Form (" * @%d", i);
  }
  
  _thePdf=new RooGenericPdf(Form("the_%s", GetName()), _pdfType+" "+GetTitle(),
			    format, _subPdfs);
  
  
  cout<<"done init of rarMultPdf for "<<GetName()<<endl<<endl;
}

/// \brief Do pdfPlot for MultPdf (Overlay component plots)
///
/// \param plotList List of plots
/// \param pdfList List of PDFs to be plotted
/// \return The last RooPlot created
///
/// Use rarBasePdf::doPdfPlot
RooPlot *rarMultPdf::doPdfPlot(TList &plotList, TString pdfList)
{
  return rarBasePdf::doPdfPlot(plotList, pdfList);
}
