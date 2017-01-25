/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarCompBase.rdl,v 1.10 2014/09/14 17:33:49 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_COMPBASE
#define RAR_COMPBASE

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief Composite base class
///
/// It is the base class for composite PDFs, like AddPdf and ProdPdf.
/// It defines common data and functions for composite PDF
/// upon #rarBasePdf.
/// \par Config Directives:
/// \verbatim
/// Comps = <name1> <name2> ... <nameN>\endverbatim
/// \p Comps defines names of all its components.
/// Each component has its own config section, namely,
/// comp \p \<name1\> is configured in config section named
/// \p "<name1> Config".
///
/// \p rarCompBase serves as base class for all
/// composite classes and it is not supposed to be referred
/// directly in the config file.
class rarCompBase : public rarBasePdf {
  
public:
  rarCompBase();
  rarCompBase(const char*configFile,const char*configSec,const char*configStr,
	      rarDatasets *theDatasets, RooDataSet *theData,
	      const char *name, const char *title,
	      Bool_t useBasePdfFit=kTRUE);
  virtual ~rarCompBase();
  
  virtual void setSimPdf(RooSimultaneous *simPdf);
  virtual RooArgSet getParams();
  
  virtual void preAction();
  virtual RooArgSet getArgSet(TString paramNames, Bool_t useRead=kFALSE,
			      RooArgSet *fullSet=0);
  virtual void doPdfFit(TString pdfList="");
  virtual void attachDataSet(const RooAbsData &data);
  virtual Bool_t isNegativeValue();
  virtual RooPlot *doPdfPlot(TList &plotList, TString pdfList="");
  virtual RooArgSet getCorrCoeffs();
  /// Return the pdfList
  virtual TList *getPdfList() {return &_pdfList;}
  virtual RooArgList getPdfsWOvar(RooArgList ignoredObs);
  virtual RooAbsPdf *getProtGen();
  
protected:
  void init();
  
  Int_t _nComp; ///< Number of RooRarFit Pdf components
  TList _pdfList; ///< RooRarFit Pdf list
  RooArgList _subProtGenPdfs; ///< subProtGenPdf ArgList
  
private:
  rarCompBase(const rarCompBase&);
  ClassDef(rarCompBase, 0) // RooRarFit Component base class
    ;
};

#endif
