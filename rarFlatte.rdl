/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: RooRarFit
 *    File: $Id: rarFlatte.rdl,v 1.4 2014/09/14 17:33:51 fwilson Exp $
 * Authors: Fergus Wilson
 * History:
 *
 * Copyright (C) 2005-2012, RAL
 *****************************************************************************/
#ifndef RAR_FLATTE
#define RAR_FLATTE

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief Flatte PDF builder
///
/// Build <a href="http://roofit.sourceforge.net/docs/classref/RooFlatte.html"
/// target=_blank>RooFlatte</a> Pdf.
/// \par Config Directives:
/// \verbatim
/// configStr = Flatte ["<Optional Title>"]
/// x = AbsReal Def
/// mean  = AbsReal Def
/// g0 = AbsReal Def
/// m0a = AbsReal Def
/// m0b = AbsReal Def
/// g1 = AbsReal Def
/// m1a = AbsReal Def
/// m1b = AbsReal Def
/// \endverbatim
/// \p x is the default observable.
/// \p mean is the peak position of the pdf.
/// \p g0 is the square of the coupling constant to the first channel (default = 0.1108 GeV)
/// \p m0a is the mass of one of the final state particles in the first channel (default = 0.1108 GeV)
/// \p m0b is the mass of other final state particle in the first channel (default = 0.13957 GeV)
/// \p g1 is the square of the coupling constant to the first channel (default = 0.4229 GeV)
/// \p m1a is the mass of one of the final state particles in the second channel (default = 0.49368 GeV)
/// \p m1b is the mass of other final state particle in the second channel (default = 0.49368 GeV)
/// After being defined, the four masses are held constant in the fit.
/// All the floating variables can be \p RooRealVar or \p RooFormulaVar.
///
class rarFlatte : public rarBasePdf {

public:
  rarFlatte();
  rarFlatte(const char *configFile, const char *configSec, 
		    const char *configStr,
		    rarDatasets *theDatasets, RooDataSet *theData,
		    const char *name, const char *title);
  virtual ~rarFlatte();
  
protected:
  void init();
  
  RooAbsReal *_x;     ///< Default obs (mass)
  RooAbsReal *_mean ; // resonance mass
  RooAbsReal *_g0 ;  // square of the coupling constant to channel 1
  RooAbsReal *_m0a ; // mass of first final state particle in channel 1
  RooAbsReal *_m0b ; // mass of second final state particle in channel 1
  RooAbsReal *_g1 ;  // square of the coupling constant to channel 2
  RooAbsReal *_m1a ; // mass of first final state particle in channel 2
  RooAbsReal *_m1b ; // mass of second final state particle in channel 2
 
private:
  rarFlatte(const rarFlatte&);
 
  ClassDef(rarFlatte, 0) // RooRarFit Flatte PDF class
    ;
};

#endif
