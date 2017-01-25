/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarGeneric.rdl,v 1.6 2014/09/14 17:33:52 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_GENERIC
#define RAR_GENERIC

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief Generic PDF builder
///
/// Build
/// <a href="http://roofit.sourceforge.net/docs/classref/RooGenericPdf.html"
/// target=_blank>RooGenericPdf</a> Pdf.
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_Generic">See doc for Generic PDF configs.</a>
class rarGeneric : public rarBasePdf {
  
public:
  rarGeneric();
  rarGeneric(const char*configFile, const char*configSec, const char*configStr,
	     rarDatasets *theDatasets, RooDataSet *theData,
	     const char *name, const char *title);
  virtual ~rarGeneric();
  
protected:
  void init();
  
private:
  rarGeneric(const rarGeneric&);
  ClassDef(rarGeneric, 0) // RooRarFit Generic Pdf class
    ;
};

#endif
