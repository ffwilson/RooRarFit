/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarNovosibirsk.rdl,v 1.5 2014/09/14 17:33:53 fwilson Exp $
 * Authors: Karsten Koeneke, Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, Massachusetts Institute of Technology, Cambridge
 *  and          2005 University of California, Riverside
 *****************************************************************************/
#ifndef RAR_NOVOSIBIRSK
#define RAR_NOVOSIBIRSK

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief Exponential PDF builder
///
/// Build <a href="http://roofit.sourceforge.net/docs/classref/RooNovosibirsk.html"
/// target=_blank>RooNovosibirsk</a> Pdf.
/// \par Config Directives:
/// \verbatim
/// configStr = Exp ["<Optional Title>"]
/// x = AbsReal Def
/// peak = AbsReal Def
/// width = AbsReal Def
/// tail  = AbsReal Def\endverbatim
/// \p x is the default observable.
/// \p peak is the peak position of the pdf.
/// \p width is the width of the pdf.
/// \p tail is the tail parameter of the pdf.
/// All the variables can be \p RooRealVar or \p RooFormulaVar.
class rarNovosibirsk : public rarBasePdf {
  
public:
  rarNovosibirsk();
  rarNovosibirsk(const char *configFile, const char *configSec, const char *configStr,
	         rarDatasets *theDatasets, RooDataSet *theData,
 	         const char *name, const char *title);
  virtual ~rarNovosibirsk();
  
protected:
  void init();
  
  RooAbsReal *_x;     ///< Default obs
  RooAbsReal *_peak;  ///< Peak position of the PDF
  RooAbsReal *_width; ///< Width of the PDF
  RooAbsReal *_tail;  ///< Tail parameter of the PDF
  
private:
  rarNovosibirsk(const rarNovosibirsk&);
  ClassDef(rarNovosibirsk, 0) // RooRarFit Novosibirsk Pdf class
    ;
};

#endif
