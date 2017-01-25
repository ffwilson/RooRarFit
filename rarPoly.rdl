/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarPoly.rdl,v 1.4 2014/09/14 17:33:54 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_POLY
#define RAR_POLY

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief RooPolynomial/RooChebychev PDF builder
///
/// Build
/// <a href="http://roofit.sourceforge.net/docs/classref/RooPolynomial.html"
/// target=_blank>RooPolynomial</a> /
/// <a href="http://roofit.sourceforge.net/docs/classref/RooChebychev.html"
/// target=_blank>RooChebychev</a>
/// Pdf.
/// \par Config Directives:
/// \verbatim
/// configStr = Polynomial ["<Optional Title>"]
/// configStr = Chebychev ["<Optional Title>"]
/// x = AbsReal Def
/// nOrder = <orderVal>
/// P01 = AbsReal Def
/// ...
/// P<orderVal> = AbsReal Def\endverbatim
/// \p x is the default observable.
/// \p nOrder is an integer representing the order of the polynomial.
/// \p P01 ... \c P\<orderVal\> define all the coeffs.
/// \p P00 is not needed.
/// All the \p AbsReal variables can be \p RooRealVar or \p RooFormulaVar.
class rarPoly : public rarBasePdf {
  
public:
  rarPoly();
  rarPoly(const char *configFile, const char *configSec, const char *configStr,
	  rarDatasets *theDatasets, RooDataSet *theData,
	  const char *name, const char *title);
  virtual ~rarPoly();
  
protected:
  void init();
  
  RooAbsReal *_x; ///< Default obs
  Int_t _nOrder; ///< Order of polynomial
  
private:
  rarPoly(const rarPoly&);
  ClassDef(rarPoly, 0) // RooRarFit Polynomial/Chebychev Pdf class
    ;
};

#endif
