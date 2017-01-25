/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarUsrPdf.rdl,v 1.5 2014/09/14 17:33:56 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

#ifndef RAR_USRPDF
#define RAR_USRPDF

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief User-defined PDF builder
///
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_UsrPdf">See doc for UsrPdf configs.</a>
class rarUsrPdf : public rarBasePdf {
  
public:
  rarUsrPdf();
  rarUsrPdf(const char *configFile,
	    const char *configSec, const char *configStr,
	    rarDatasets *theDatasets, RooDataSet *theData,
	    const char *name, const char *title);
  virtual ~rarUsrPdf();
  
protected:
  void init();
  
  RooAbsReal *_x; ///< Default obs
  RooAbsReal *_a; ///< Param a
  RooAbsReal *_b; ///< Param b
  RooAbsReal *_c; ///< Param c
  RooAbsReal *_d; ///< Param d
  RooAbsReal *_e; ///< Param e
  
private:
  rarUsrPdf(const rarUsrPdf&);
  ClassDef(rarUsrPdf, 0) // RooRarFit User-defined Pdf class
    ;
};

#endif
