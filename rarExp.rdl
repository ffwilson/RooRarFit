/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarExp.rdl,v 1.5 2014/09/14 17:33:51 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_EXP
#define RAR_EXP

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief Exponential PDF builder
///
/// Build <a href="http://roofit.sourceforge.net/docs/classref/RooExponential.html"
/// target=_blank>RooExponential</a> Pdf.
/// \par Config Directives:
/// \verbatim
/// configStr = Exp ["<Optional Title>"]
/// x = AbsReal Def
/// c = AbsReal Def\endverbatim
/// \p x is the default observable.
/// \p c is the exponent of the pdf.
/// All the variables can be \p RooRealVar or \p RooFormulaVar.
class rarExp : public rarBasePdf {
  
public:
  rarExp();
  rarExp(const char *configFile, const char *configSec, const char *configStr,
	 rarDatasets *theDatasets, RooDataSet *theData,
	 const char *name, const char *title);
  virtual ~rarExp();
  
protected:
  void init();
  
  RooAbsReal *_x; ///< Default obs
  RooAbsReal *_c; ///< Exponent of the PDF
  
private:
  rarExp(const rarExp&);
  ClassDef(rarExp, 0) // RooRarFit Exponential Pdf class
    ;
};

#endif
