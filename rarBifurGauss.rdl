/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarBifurGauss.rdl,v 1.5 2014/09/14 17:33:49 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_BIFURGAUSS
#define RAR_BIFURGAUSS

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief BifurGauss PDF builder
///
/// Build <a href="http://roofit.sourceforge.net/docs/classref/RooBifurGauss.html"
/// target=_blank>RooBifurGauss</a> Pdf
/// \par Config Directives:
/// \verbatim
/// configStr = BifurGauss ["<Optional Title>"]
/// configStr = BGGauss ["<Optional Title>"]
/// parSymLevel = <0|1|2|3>
/// x = AbsReal Def
/// peak = AbsReal Def
/// sigL = AbsReal Def
/// sigR = AbsReal Def
/// mean = AbsReal Def
/// rms = AbsReal Def
/// asym = AbsReal Def\endverbatim
/// If \p parSymLevel = 0, use the default BifurGauss;
/// \p mean, \p rms, and \p asym are not required.
/// If \p parSymLevel = 1, use lowest-order mapping;
/// if = 2, keep O(A^2, A^3) terms also;
/// if = 3, take \p asym to mean 3rd moment.
/// With \p parSymLevel > 0, \p mean, \p rms, and \p asym are used;
/// \p peak, \p sigL, and \p sigR are hard-coded to
/// the corresponding \p RooFormulaVar.
/// \p x is the default observable.
/// Pdf type \p BifurGauss and \p BGGauss are interchangeable,
/// with the only difference being that
/// if \p parSymLevel not specified,
/// the default value of it for \p BifurGauss is 0,
/// while that for \p BGGauss is 1.
/// All the \p AbsReal variables can be \p RooRealVar or \p RooFormulaVar.
class rarBifurGauss : public rarBasePdf {
  
public:
  rarBifurGauss();
  rarBifurGauss(const char *configFile, const char *configSec,
		const char *configStr,
		rarDatasets *theDatasets, RooDataSet *theData,
		const char *name, const char *title);
  virtual ~rarBifurGauss();
  
protected:
  void init();
  
  RooAbsReal *_x; ///< Default obs
  Int_t _parSymLevel; ///< parSymLevel
  RooAbsReal *_peak; ///< Pdf mean
  RooAbsReal *_sigL; ///< Left side sigma
  RooAbsReal *_sigR; ///< Right side sigma
  
private:
  rarBifurGauss(const rarBifurGauss&);
  ClassDef(rarBifurGauss, 0) // RooRarFit BifurGauss Pdf class
    ;
};

#endif
