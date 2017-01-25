/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarNLL.rdl,v 1.9 2014/09/14 17:33:53 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_NLL
#define RAR_NLL

#include "TList.h"
#include "TMatrixD.h"
#include "TString.h"
#include "TArrayD.h"
#include "TObject.h"

#include "RooCurve.h"

/// \brief RooRarFit class for NLL manipulation
///
/// It calculate NLL values, likelihood integrals, etc.
/// based on NLL curve.
class rarNLL : public TNamed {
  
public:
  rarNLL();
  rarNLL(RooCurve *curve,
	 const char *name="theRARNLL", const char *title="the RARNLL",
         const Bool_t verbose=kFALSE);
  virtual ~rarNLL();
  
  /// \brief Set verbose mode
  /// \param verbose Verbose mode
  virtual void setVerbose(Bool_t verbose=kTRUE) {_verbose=verbose;}
  
  virtual void init(RooCurve *curve=0);
  Double_t getNLL(Double_t x);
  Double_t getY(Double_t x);
  TArrayD getX(Double_t y);
  void getMin(Double_t &x, Double_t &y);
  TArrayD getMin();
  Double_t getLIntegral(Double_t x);
  Double_t getLIntegral(Double_t xl, Double_t xh);
  Double_t getLIntegral();
  Double_t getLIntegralInverse(Double_t xl, Double_t iVal);
  Double_t getLIntegralInverse(Double_t iVal);
  Double_t lIntegralFunc(Double_t x, TMatrixD &A, TMatrixD &lA);
  Double_t lIntegralFuncInverse(Double_t x0, Double_t x2, Int_t iter,
                                Double_t &thisXI,
                                Double_t &la, Double_t &lb, Double_t &lc);

protected:
  void getMin(TArrayD &xy, Double_t x, Double_t a, Double_t b, Double_t c);
  TArrayD getMin(Double_t x0, Double_t x1, Double_t x2,
		 Double_t a, Double_t b, Double_t c);
  
  Int_t _nPoints; ///< Number of points in NLL curve
  TArrayD _xs; ///< NLL curve x values
  TArrayD _ys; ///< NLL curve y values
  Int_t _mIdx; ///< NLL point index for minY
  Int_t _nSteps; ///< Number of calculation steps
  TList _coeffMList; ///< List of coeff matrices
  TList _lCoefMList; ///< List of coeff matrices for likelihood curve fit
  TArrayD _x0s; ///< Starting point for each step
  TArrayD _x1s; ///< Middle point for each step
  TArrayD _x2s; ///< Ending point for each step
  TArrayD _iLIntegrals; ///< Lower likelihood integrals within each step
  TArrayD _tLIntegrals; ///< Total likelihood integrals before each step
  
  Bool_t _verbose; ///< Boolean to control debug output
  
private:
  rarNLL(const rarNLL&);
  ClassDef(rarNLL, 0) // RooRarFit NLL class
    ;
};

#endif
