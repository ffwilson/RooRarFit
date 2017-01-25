/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: RooBinnedPdf.rdl,v 1.5 2014/09/14 17:33:46 fwilson Exp $
 * Authors:                                                                  *
 *    Aaron Roodman, Stanford Linear Accelerator Center, Stanford University *
 *    Modified by Wouter                                                     *
 *                                                                           *
 * Copyright (C) 2005-2012, Stanford University. All rights reserved.        *
 *           
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_BINNEDPDF_HH
#define ROO_BINNEDPDF_HH

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "TArrayD.h"
//#include <iostream>
using namespace std;


class RooRealVar;
class RooArgList ;

class RooBinnedPdf : public RooAbsPdf {

private:
  Double_t localEval(const Double_t) const;

public:

  RooBinnedPdf(const char *name, const char *title,
		RooAbsReal& x, const RooArgList& coefList, const TArrayD& limits) ;

  RooBinnedPdf(const RooBinnedPdf& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooBinnedPdf(*this, newname); }
  virtual ~RooBinnedPdf() ;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
  Int_t getnBins();
  Double_t* getLimits();

protected:

  RooRealProxy _x;
  RooListProxy _coefList ;
  TArrayD _limits;
  Int_t _nBins ;
  Double_t evaluate() const;

  ClassDef(RooBinnedPdf,1) // Parametric Step Function Pdf
};

#endif
