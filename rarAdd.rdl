/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarAdd.rdl,v 1.8 2014/09/14 17:33:48 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_ADD
#define RAR_ADD

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarCompBase.hh"

/// \brief AddPdf/AddModel builder.
///
/// Build composite pdfs through
/// <a href="http://roofit.sourceforge.net/docs/classref/RooAddPdf.html"
/// target=_blank>RooAddPdf</a> /
/// <a href="http://roofit.sourceforge.net/docs/classref/RooAddModel.html"
/// target=_blank>RooAddModel</a>.
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_Add">See doc for AddPdf/AddModel configs.</a>
class rarAdd : public rarCompBase {
  
public:
  rarAdd();
  rarAdd(const char *configFile, const char *configSec,
	 const char *configStr,
	 rarDatasets *theDatasets, RooDataSet *theData,
	 const char *name, const char *title,
	 Bool_t useBasePdfFit=kTRUE, Bool_t buildAddPdf=kTRUE);
  virtual ~rarAdd();
  
  /// return coeff List
  virtual RooArgList getCoeffList() {return _coeffs;}
  virtual RooPlot *doPdfPlot(TList &plotList, TString pdfList="");
  virtual RooAbsPdf *getPdfWOvar(RooArgList ignoredObs);
  virtual RooAbsPdf *getProtGen();
  
protected:
  void init();
  
  Int_t _nCoeff; ///< Number of coeffs
  
private:
  Bool_t _buildAddPdf; ///< If build addpdf within this rarAdd
  
private:
  rarAdd(const rarAdd&);
  ClassDef(rarAdd, 0) // RooRarFit AddPdf/AddModel class
    ;
};

#endif
