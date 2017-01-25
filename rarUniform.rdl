/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarUniform.rdl,v 1.2 2014/09/14 17:33:56 fwilson Exp $
 * Authors: F. Wilson
 * History:
 *
 * Copyright (C) 2005-2012, RAL/STFC
 *****************************************************************************/
#ifndef RAR_UNIFORM
#define RAR_UNIFORM

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief RooUniform PDF builder
///
/// Build
/// <a href="http://roofit.sourceforge.net/docs/classref/RooUniform.html"
/// target=_blank>RooUniform</a> /
/// Pdf.
/// \par Config Directives:
/// \verbatim
/// configStr = Uniform ["<Optional Title>"]
/// obs = AbsReal Def
// \p obs is the list of observables.
/// All the \p AbsReal variables can be \p RooRealVar or \p RooFormulaVar.
class rarUniform : public rarBasePdf {
  
public:
  rarUniform();
  rarUniform(const char *configFile, const char *configSec, const char *configStr,
	  rarDatasets *theDatasets, RooDataSet *theData,
	  const char *name, const char *title);
  virtual ~rarUniform();
  
protected:
  void init();
  
  RooAbsReal *_x; ///< Default obs
  
private:
  rarUniform(const rarUniform&);
  ClassDef(rarUniform, 0) // RooRarFit Uniform Pdf class
    ;
};

#endif
