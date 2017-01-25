/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarTriGauss.rdl,v 1.6 2014/09/14 17:33:55 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_TRIGAUSS
#define RAR_TRIGAUSS

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief TripleGaussian PDF/Model/GexpShape builder
///
/// Build TripleGaussian PDF/Model/GexpShape.
/// This class is for convenience of three commonly used
/// composite pdfs related to triple-Gaussian.
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_TriGauss">See doc for TriGauss configs.</a>
class rarTriGauss : public rarBasePdf {
  
public:
  rarTriGauss();
  rarTriGauss(const char*configFile,const char*configSec,const char*configStr,
	      rarDatasets *theDatasets, RooDataSet *theData,
	      const char *name, const char *title);
  virtual ~rarTriGauss();
  
protected:
  void init();
  
  RooAbsReal *_x; ///< Default obs
  RooAbsReal *_tau; ///< \e B0 lifetime
  RooAbsReal *_msSF; ///< Scale factor for mean and sigma
  RooAbsReal *_meanC; ///< Mean of the core Gaussian
  RooAbsReal *_sigmaC; ///< Sigma of the core Gaussian
  RooAbsReal *_meanT; ///< Mean of the tail Gaussian
  RooAbsReal *_sigmaT; ///< Sigma of the tail Gaussian
  RooAbsReal *_meanO; ///< Mean of the outliner Gaussian
  RooAbsReal *_sigmaO; ///< Sigma of the outliner Gaussian
  RooAbsReal *_fracC; ///< Fraction of the core Gaussian
  RooAbsReal *_fracO; ///< Fraction of the outliner Gaussian
  
private:
  rarTriGauss(const rarTriGauss&);
  ClassDef(rarTriGauss, 0) // RooRarFit TriGauss related Pdf class
    ;
};

#endif
