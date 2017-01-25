/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: RooRarFit
 *    File: $Id: rarThreshold.rdl,v 1.2 2014/09/14 17:33:55 fwilson Exp $
 * Authors: Fergus Wilson
 * History:
 *
 * Copyright (C) 2005-2012, RAL
 *****************************************************************************/
#ifndef RAR_THRESHOLD
#define RAR_THRESHOLD

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief Threshold PDF builder
///
/// Build <a href="http://roofit.sourceforge.net/docs/classref/RooThreshold.html"
/// target=_blank>RooThreshold</a> Pdf.
/// \par Config Directives:
/// \verbatim
/// configStr = Threshold ["<Optional Title>"]
/// x = AbsReal Def
/// m0  = AbsReal Def
/// power = AbsReal Def
/// nOrder = <orderval>
/// P01 = AbsReal Def
/// P02 = ...
/// \endverbatim
/// \p x is the default observable.
/// \p m0 is the threshold. PDF= for x<m0.
/// \p power is the power
/// \p nOrder is an integer representing the order of the polynomial.
/// \p P01 ... \c P\<orderVal\> define all the coeffs.
/// \p P00 is not needed.
/// All the floating variables can be \p RooRealVar or \p RooFormulaVar.
///
class rarThreshold : public rarBasePdf {

public:
  rarThreshold();
  rarThreshold(const char *configFile, const char *configSec, 
		    const char *configStr,
		    rarDatasets *theDatasets, RooDataSet *theData,
		    const char *name, const char *title);
  virtual ~rarThreshold();
  
protected:
  void init();
  
  RooAbsReal *_x ;     ///< Default obs (mass)
  RooAbsReal *_m0 ;  //< square of the coupling constant to channel 1
  RooAbsReal *_power ; //< resonance mass
  Int_t _nOrder ; //< order of polynomial
 
private:
  rarThreshold(const rarThreshold&);
 
  ClassDef(rarThreshold, 0) // RooRarFit Threshold PDF class
    ;
};

#endif
