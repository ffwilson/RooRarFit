/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarGaussian.rdl,v 1.5 2014/09/14 17:33:51 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_GAUSSIAN
#define RAR_GAUSSIAN

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief RooGaussian/RooBreitWigner PDF builder
///
/// Build
/// <a href="http://roofit.sourceforge.net/docs/classref/RooGaussian.html"
/// target=_blank>RooGaussian</a> /
/// <a href="http://roofit.sourceforge.net/docs/classref/RooBreitWigner.html"
/// target=_blank>RooBreitWigner</a>
/// Pdf.
/// \par Config Directives:
/// \verbatim
/// configStr = Gaussian ["<Optional Title>"]
/// configStr = BreitWigner ["<Optional Title>"]
/// x = AbsReal Def
/// mean = AbsReal Def
/// sigma = AbsReal Def\endverbatim
/// \p x is the default observable.
/// \p mean is the mean of the PDF.
/// \p sigma is the sigma of the PDF.
/// All the variables can be \p RooRealVar or \p RooFormulaVar.
class rarGaussian : public rarBasePdf {
  
public:
  rarGaussian();
  rarGaussian(const char*configFile,const char*configSec,const char*configStr,
	      rarDatasets *theDatasets, RooDataSet *theData,
	      const char *name, const char *title);
  virtual ~rarGaussian();
  
protected:
  void init();
  
  RooAbsReal *_x; ///< Default obs
  RooAbsReal *_mean; ///< Mean of PDF
  RooAbsReal *_sigma; ///< Sigma of PDF
  RooAbsReal *_scale; ///< Scale of the sigma
  RooAbsReal *_shift; ///< Shift of the mean
  
private:
  rarGaussian(const rarGaussian&);
  ClassDef(rarGaussian, 0) // RooRarFit Gaussian/BreitWigner Pdf class
    ;
};

#endif
