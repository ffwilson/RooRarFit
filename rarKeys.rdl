/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarKeys.rdl,v 1.4 2014/09/14 17:33:52 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_KEYS
#define RAR_KEYS

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief 1/2D Keys PDF builder
///
/// Build
/// <a href="http://roofit.sourceforge.net/docs/classref/RooKeysPdf.html"
/// target=_blank>RooKeysPdf</a> /
/// <a href="http://roofit.sourceforge.net/docs/classref/Roo2DKeysPdf.html"
/// target=_blank>Roo2DKeysPdf</a>
/// Pdf.
/// \par Config Directives:
/// \verbatim
/// configStr = Keys ["<Optional Title>"]
/// configStr = 2DKeys ["<Optional Title>"]
/// x = AbsReal Def
/// y = AbsReal Def
/// rho = Double_t
/// keysOption = Options\endverbatim
/// \p x and \p y are the default observables (\p y for 2D only).
/// \p rho is width scale factor.
/// \p keysOption is 1/2D Keys opitons.
/// For 1D, the options can be any enum value of
/// <a href="http://roofit.sourceforge.net/docs/classref/RooKeysPdf.html#RooKeysPdf:Mirror"
/// target=_blank>RooKeysPdf::Mirror</a>,
/// for 2D, the options can be any valid characters of
/// <a href="http://roofit.sourceforge.net/docs/classref/Roo2DKeysPdf.html#Roo2DKeysPdf:setOptions"
/// target=_blank>Roo2DKeysPdf::setOptions</a>.
/// All the \p AbsReal parameters can be \p RooRealVar or \p RooFormulaVar.
class rarKeys : public rarBasePdf {
  
public:
  rarKeys();
  rarKeys(const char *configFile, const char *configSec, const char *configStr,
	  rarDatasets *theDatasets, RooDataSet *theData,
	  const char *name, const char *title);
  virtual ~rarKeys();
  
  void setFitData(RooDataSet *theData=0);
  
protected:
  void init();
  
  RooAbsReal *_x; ///< Default obs
  RooAbsReal *_y; ///< Default obs
  Double_t _rho; ///< Width scale factor
  TString _keysOption; ///< Options
  
private:
  rarKeys(const rarKeys&);
  ClassDef(rarKeys, 0) // RooRarFit 1D/2D Keys Pdf class
    ;
};

#endif
