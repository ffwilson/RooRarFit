/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarSimPdf.rdl,v 1.10 2014/09/14 17:33:55 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_SIMPDF
#define RAR_SIMPDF

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarCompBase.hh"

class RooSimPdfBuilder;

/// \brief RooSimultaneous PDF builder
///
/// Build
/// <a href="http://roofit.sourceforge.net/docs/classref/RooSimultaneous.html"
/// target=_blank>RooSimultaneous</a>
/// Pdf.
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_SimPdf">See doc for SimPdf configs.</a>
class rarSimPdf : public rarCompBase {
  
public:
  rarSimPdf();
  rarSimPdf(const char*configFile, const char*configSec, const char*configStr,
	    rarDatasets *theDatasets, RooDataSet *theData,
	    const char *name, const char *title);
  virtual ~rarSimPdf();
  
  virtual RooAbsPdf *getPdfWOvar(RooArgList ignoredObs);
  virtual RooAbsPdf *getProtGen();
  
protected:
  void init();
  
  RooSimPdfBuilder *_simBuilder; ///< SimPdf builder
  RooArgSet *_simConfig; ///< SimPdf config ArgSet
  
private:
  rarSimPdf(const rarSimPdf&);
  ClassDef(rarSimPdf, 0) // RooRarFit Simultaneous Pdf class
    ;
};

#endif
