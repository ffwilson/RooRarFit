/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarHistPdf.rdl,v 1.5 2014/09/14 17:33:52 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012,, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_HISTPDF
#define RAR_HISTPDF

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

class RooDataHist;

/// \brief HistPdf PDF builder
///
/// Build
/// <a href="http://roofit.sourceforge.net/docs/classref/RooHistPdf.html"
/// target=_blank>RooHistPdf</a> Pdf.
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_HistPdf">See doc for HistPdf PDF configs.</a>
class rarHistPdf : public rarBasePdf {
  
public:
  rarHistPdf();
  rarHistPdf(const char*configFile, const char*configSec, const char*configStr,
	     rarDatasets *theDatasets, RooDataSet *theData,
	     const char *name, const char *title);
  virtual ~rarHistPdf();
  
protected:
  void init();
  
  RooDataHist *_theHist; ///< RooDataHist for the PDF
  
private:
  rarHistPdf(const rarHistPdf&);
  ClassDef(rarHistPdf, 0) // RooRarFit HistPdf class
    ;
};

#endif
