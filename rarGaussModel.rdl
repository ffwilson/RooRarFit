/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarGaussModel.rdl,v 1.4 2014/09/14 17:33:51 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_GAUSSMODEL
#define RAR_GAUSSMODEL

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief RooGaussModel PDF builder
///
/// Build
/// <a href="http://roofit.sourceforge.net/docs/classref/RooGaussModel.html"
/// target=_blank>RooGaussModel</a>
/// Pdf.
/// \par Config Directives:
/// \verbatim
/// configStr = GaussModel ["<Optional Title>"]
/// x = RealVar Def
/// mean = AbsReal Def
/// sigma = AbsReal Def
/// msSF = AbsReal Def
/// meanSF = AbsReal Def
/// sigmaSF = AbsReal Def
/// FlatScaleFactorIntegral = <yes|no>\endverbatim
/// \p x is the default observable (deltaT) and should be \p RooRealVar.
/// \p mean is the mean of the PDF.
/// \p sigma is the sigma of the PDF.
/// \p msSF is the scale factor for both mean and sigma.
/// \p meanSF is the mean scale factor of the PDF.
/// \p sigmaSF is the sigma scale factor of the PDF.
/// One can choose to specify \p msSF for both mean and sigma,
/// or choose \p meanSF for mean, and/or \p sigmaSF for sigma,
/// or if none is chosen, the default values are 1.
/// All the \p AbsReal parameters can be \p RooRealVar or \p RooFormulaVar.
/// If optional \p FlatScaleFactorIntegral is set to \p no (default \p yes),
/// \p RooGaussModel::advertiseFlatScaleFactorIntegral will be called
/// with \p kFALSE, otherwise with \p kTRUE.
class rarGaussModel : public rarBasePdf {
  
public:
  rarGaussModel();
  rarGaussModel(const char *configFile, const char *configSec,
		const char *configStr,
		rarDatasets *theDatasets, RooDataSet *theData,
		const char *name, const char *title);
  virtual ~rarGaussModel();
  
protected:
  void init();
  
  RooRealVar *_x; ///< Default obs (deltaT)
  RooAbsReal *_mean; ///< Mean of PDF
  RooAbsReal *_sigma; ///< Sigma of PDF
  RooAbsReal *_msSF; ///< Scale factor for mean and sigma
  RooAbsReal *_meanSF; ///< Scale factor of the mean
  RooAbsReal *_sigmaSF; ///< Scale factor of the sigma
  
private:
  rarGaussModel(const rarGaussModel&);
  ClassDef(rarGaussModel, 0) // RooRarFit Gaussian Resolution Model class
    ;
};

#endif
