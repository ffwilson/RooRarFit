/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarBallack.rdl,v 1.5 2014/09/14 17:33:49 fwilson Exp $
 * Authors: Karsten Koeneke, Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, Massachusetts Institute of Technology, Cambridge
 *  and          2005 University of California, Riverside
 *****************************************************************************/
#ifndef RAR_BALLACK
#define RAR_BALLACK

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief Exponential PDF builder
///
/// Build <a href="http://roofit.sourceforge.net/docs/classref/RooBallack.html"
/// target=_blank>RooBallack</a> Pdf.
/// \par Config Directives:
/// \verbatim
/// configStr = Exp ["<Optional Title>"]
/// x = AbsReal Def
/// mean  = AbsReal Def
/// width = AbsReal Def
/// tail  = AbsReal Def
/// alpha = AbsReal Def
/// n     = AbsReal Def\endverbatim
/// \p x is the default observable.
/// \p mean is the peak position of the pdf.
/// \p width is the width of the pdf.
/// \p tail is the tail parameter of the pdf.
/// \p alpha is the atachment point of the polynomial tail.
/// \p n is the power of the polynomial of the pdf.
/// All the variables can be \p RooRealVar or \p RooFormulaVar.
class rarBallack : public rarBasePdf {
  
public:
  rarBallack();
  rarBallack(const char *configFile, const char *configSec, const char *configStr,
	     rarDatasets *theDatasets, RooDataSet *theData,
 	     const char *name, const char *title);
  virtual ~rarBallack();
  
protected:
  void init();
  
  RooAbsReal *_x;     ///< Default obs
  RooAbsReal *_mean;  ///< Peak position of the PDF
  RooAbsReal *_width; ///< Width of the PDF
  RooAbsReal *_tail;  ///< Tail parameter of the PDF
  RooAbsReal *_alpha; ///< Attaching point
  RooAbsReal *_n;     ///< Attached polynomial power
  
private:
  rarBallack(const rarBallack&);
  ClassDef(rarBallack, 0) // RooRarFit Ballack Pdf class
    ;
};

#endif
