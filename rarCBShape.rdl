/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarCBShape.rdl,v 1.4 2014/09/14 17:33:49 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_CBSHAPE
#define RAR_CBSHAPE

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief CBShape PDF builder
///
/// Build <a href="http://roofit.sourceforge.net/docs/classref/RooCBShape.html"
/// target=_blank>RooCBShape</a> Pdf.
/// \par Config Directives:
/// \verbatim
/// configStr = CBShape ["<Optional Title>"]
/// x = AbsReal Def
/// mean = AbsReal Def
/// sigma = AbsReal Def
/// alpha = AbsReal Def
/// n = AbsReal Def\endverbatim
/// \p x is the default observable.
/// \p mean is the mean of the pdf,
/// \p sigma is the sigma of the pdf,
/// \p alpha is the alpha of the pdf,
/// \p n is the n of the pdf.
/// All the variables can be \p RooRealVar or \p RooFormulaVar.
class rarCBShape : public rarBasePdf {
  
public:
  rarCBShape();
  rarCBShape(const char *configFile, const char *configSec,
	     const char *configStr,
	     rarDatasets *theDatasets, RooDataSet *theData,
	     const char *name, const char *title);
  virtual ~rarCBShape();
  
protected:
  void init();
  
  RooAbsReal *_x; ///< Default obs
  RooAbsReal *_mean; ///< Pdf mean
  RooAbsReal *_sigma; ///< Pdf sigma
  RooAbsReal *_alpha; ///< Pdf alpha
  RooAbsReal *_n; ///< Pdf n
  
private:
  rarCBShape(const rarCBShape&);
  ClassDef(rarCBShape, 0) // RooRarFit CBShape Pdf class
    ;
};

#endif
