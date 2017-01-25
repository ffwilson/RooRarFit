/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarVoigtian.rdl,v 1.3 2014/09/14 17:33:56 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_VOIGTIAN
#define RAR_VOIGTIAN

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief RooVoigtian PDF builder
///
/// Build
/// <a href="http://roofit.sourceforge.net/docs/classref/RooVoigtian.html"
/// target=_blank>RooVoigtian</a> Pdf.
/// \par Config Directives:
/// \verbatim
/// configStr = Voigtian ["<Optional Title>"]
/// x = AbsReal Def
/// mean = AbsReal Def
/// sigma = AbsReal Def\endverbatim
/// \p x is the default observable.
/// \p mean is the mean of the Breit-Wigner PDF.
/// \p width is the width of the Breit-Wigner PDF.
/// \p sigma is the width of gaussian that is convoluted with the Breit-Wigner PDF.
/// All the variables can be \p RooRealVar or \p RooFormulaVar.
class rarVoigtian : public rarBasePdf {
  
public:
  rarVoigtian();
  rarVoigtian(const char*configFile,const char*configSec,const char*configStr,
	      rarDatasets *theDatasets, RooDataSet *theData,
	      const char *name, const char *title);
  virtual ~rarVoigtian();
  
protected:
  void init();
  
  RooAbsReal *_x; ///< Default obs
  RooAbsReal *_mean; ///< Mean of Breit-Wigner PDF
  RooAbsReal *_width; ///< Sigma of Breit-Wigner PDF
  RooAbsReal *_sigma; ///< Sigma of Gaussian convolution function PDF
  
private:
  rarVoigtian(const rarVoigtian&);
  ClassDef(rarVoigtian, 0) // RooRarFit Voigtian Pdf class
    ;
};

#endif
