/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarDecay.rdl,v 1.8 2014/09/14 17:33:50 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_DECAY
#define RAR_DECAY

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief BDecay/Decay Model builder
///
/// Build
/// <a href="http://roofit.sourceforge.net/docs/classref/RooBCPGenDecay.html"
/// target=_blank>RooBCPGenDecay</a> /
/// <a href="http://roofit.sourceforge.net/docs/classref/RooBDecay.html"
/// target=_blank>RooBDecay</a> /
/// <a href="http://roofit.sourceforge.net/docs/classref/RooDecay.html"
/// target=_blank>RooDecay</a>
/// model.
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_Decay">See doc for Decay configs.</a>
class rarDecay : public rarBasePdf {
  
public:
  rarDecay();
  rarDecay(const char*configFile, const char*configSec, const char*configStr,
	   rarDatasets *theDatasets, RooDataSet *theData,
	   const char *name, const char *title);
  virtual ~rarDecay();
  
  virtual RooAbsPdf *getProtGen();
  
protected:
  void init();
  
  RooRealVar *_x; ///< Default obs (deltaT)
  RooAbsReal *_tau; ///< \e B0 lifetime
  rarBasePdf *_model; ///< resolution model
  TString _decayType; ///< Decay type
  RooAbsReal *_dm; ///< mixing frequency
  RooAbsReal *_dgamma; ///< DeltaGamma
  RooAbsCategory *_tag; ///< tagFlav
  RooAbsReal *_w; ///< avg mistag rate
  RooAbsReal *_dw; ///< diff mistag rate
  RooAbsReal *_mu; ///< diff tagging efficiency
  TString _blindStatus; ///< Blind status
  TString _blindString; ///< Blind string
  RooAbsReal *_Cb; ///< C blinded
  RooAbsReal *_Sb; ///< S blinded
  RooAbsReal *_f0; ///< f0
  RooAbsReal *_f1; ///< f1
  RooAbsReal *_f2; ///< f2
  RooAbsReal *_f3; ///< f3
  
private:
  rarDecay(const rarDecay&);
  ClassDef(rarDecay, 0) // RooRarFit BDecay/Decay model class
    ;
};

#endif
