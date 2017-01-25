/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarCruijff.rdl,v 1.7 2014/09/14 17:33:50 fwilson Exp $
 * Authors: Karsten Koeneke
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

#ifndef RAR_CRUIJFF
#define RAR_CRUIJFF

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief User-defined PDF builder
///
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_UsrPdf">See doc for UsrPdf configs.</a>
class rarCruijff : public rarBasePdf {
  
public:
  rarCruijff();
  rarCruijff(const char *configFile,
	     const char *configSec, const char *configStr,
	     rarDatasets *theDatasets, RooDataSet *theData,
	     const char *name, const char *title);
  virtual ~rarCruijff();
  
protected:
  void init();
  
  RooAbsReal *_x; ///< Default obs
  RooAbsReal *_m0; ///< mean of the distribution
  RooAbsReal *_sigmaL; ///< Left handed width
  RooAbsReal *_sigmaR; ///< Right handed width
  RooAbsReal *_alphaL; ///< Left handed alpha
  RooAbsReal *_alphaR; ///< Right haded alpha
  
private:
  rarCruijff(const rarCruijff&);
  ClassDef(rarCruijff, 0) // RooRarFit Cruijff Pdf class
    ;
};

#endif
