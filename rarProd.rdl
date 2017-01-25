/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarProd.rdl,v 1.9 2014/09/14 17:33:54 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_PROD
#define RAR_PROD

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarCompBase.hh"

/// \brief ProdPdf builder.
///
/// Build composite pdfs through
/// <a href="http://roofit.sourceforge.net/docs/classref/RooProdPdf.html"
/// target=_blank>RooProdPdf</a>.
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_Prod">See doc for ProdPdf configs.</a>
class rarProd : public rarCompBase {
  
public:
  rarProd();
  rarProd(const char *configFile, const char *configSec, const char *configStr,
	  rarDatasets *theDatasets, RooDataSet *theData,
	  const char *name, const char *title);
  virtual ~rarProd();
  
  virtual RooAbsPdf *getPdfWOvar(RooArgList ignoredObs);
  virtual RooAbsPdf *getDPdfWvar(RooRealVar *theVar);
  virtual RooAbsPdf *getProtGen();
  
protected:
  void init();
  
  TList _condPdfList; ///< RooRarFit Pdf list
  RooArgList _condPdfs; ///< condPdf ArgList

  
private:
  rarProd(const rarProd&);
  ClassDef(rarProd, 0) // RooRarFit ProdPdf class
    ;
};

#endif
