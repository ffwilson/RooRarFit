/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarArgusBG.rdl,v 1.4 2014/09/14 17:33:48 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_ARGUSBG
#define RAR_ARGUSBG

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief ArgusBG PDF builder
///
/// Build
/// <a href="http://roofit.sourceforge.net/docs/classref/RooArgusBG.html"
/// target=_blank>RooArgusBG</a> Pdf.
/// \par Config Directives:
/// \verbatim
/// configStr = ArgusBG ["<Optional Title>"]
/// x = AbsReal Def
/// max = AbsReal Def
/// c = AbsReal Def
/// pow = AbsReal Def\endverbatim
/// \p x is the default observable.
/// \p max is the end point of the pdf, usually a constant.
/// \p c is the slope parameter of the pdf.
/// \p pow is the power of the multiplying term (normally 1/2).
/// All the variables can be \p RooRealVar or \p RooFormulaVar.
class rarArgusBG : public rarBasePdf {
  
public:
  rarArgusBG();
  rarArgusBG(const char *configFile, const char *configSec,
	     const char *configStr,
	     rarDatasets *theDatasets, RooDataSet *theData,
	     const char *name, const char *title);
  virtual ~rarArgusBG();
  
protected:
  void init();
  
  RooAbsReal *_x; ///< Default obs
  RooAbsReal *_max; ///< End point
  RooAbsReal *_c; ///< Slope parameter
  RooAbsReal *_pow; ///< Power (normally 1/2)
  
private:
  rarArgusBG(const rarArgusBG&);
  ClassDef(rarArgusBG, 0) // RooRarFit ArgusBG Pdf class
    ;
};

#endif
