/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarMultPdf.rdl,v 1.6 2014/09/14 17:33:53 fwilson Exp $
 * Authors: Lei Zhang, Wolfgang Gradl
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_MULTPDF
#define RAR_MULTPDF

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarCompBase.hh"

/// \brief MultPdfPdf builder.
///
/// Build composite pdfs through
/// <a href="http://roofit.sourceforge.net/docs/classref/RooGenericPdf.html"
/// target=_blank>RooMultGeneric</a>
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_MultPdf">See doc for MultPdf configs.</a>
class rarMultPdf : public rarCompBase {
  
public:
  rarMultPdf();
  rarMultPdf(const char *configFile, const char *configSec,
	     const char *configStr,
	     rarDatasets *theDatasets, RooDataSet *theData,
	     const char *name, const char *title);
  virtual ~rarMultPdf();
  
  virtual RooPlot *doPdfPlot(TList &plotList, TString pdfList="");
  //virtual RooAbsPdf *getPdfWOvar(RooArgList ignoredObs);
  //virtual RooAbsPdf *getDPdfWvar(RooRealVar *theVar);
  /// \todo should have its own getProtGen(), similar to that for prodpdf
  //virtual RooAbsPdf *getProtGen();
  
protected:
  void init();
  
 private:
  rarMultPdf(const rarMultPdf&);
  ClassDef(rarMultPdf, 0) // RooRarFit Multiple class
    ;
};

#endif
