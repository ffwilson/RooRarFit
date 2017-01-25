/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarTwoGauss.rdl,v 1.6 2014/09/14 17:33:56 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_TWOGAUSS
#define RAR_TWOGAUSS

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief Double-Gaussian PDF builder
///
/// Build Double-Gaussian PDF.
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_TwoGauss">See doc for TwoGauss configs.</a>
class rarTwoGauss : public rarBasePdf {
  
public:
  rarTwoGauss();
  rarTwoGauss(const char*configFile,const char*configSec,const char*configStr,
	      rarDatasets *theDatasets, RooDataSet *theData,
	      const char *name, const char *title);
  virtual ~rarTwoGauss();
  
protected:
  void init();
  
  RooAbsReal *_x; ///< Default obs
  RooAbsReal *_meanC; ///< Mean of core Gaussian
  RooAbsReal *_sigmaC; ///< Sigma of core Gaussian
  RooAbsReal *_meanT; ///< Mean of tail Gaussian
  RooAbsReal *_sigmaT; ///< Sigma of tail Gaussian
  RooAbsReal *_fracC; ///< Fraction of core Gaussian
  RooAbsReal *_scale; ///< Scale of the core sigma
  RooAbsReal *_shift; ///< Shift of the two means
  
private:
  rarTwoGauss(const rarTwoGauss&);
  ClassDef(rarTwoGauss, 0) // RooRarFit TwoGauss Pdf class
    ;
};

#endif
