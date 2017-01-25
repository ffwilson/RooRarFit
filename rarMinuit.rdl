/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarMinuit.rdl,v 1.9 2014/09/14 17:33:53 fwilson Exp $
 * Author:                                                                   *
 *   WTF, Bill Ford, U. of Colorado, wtford@pizero.colorado.edu              *
 *                                                                           *
 * Copyright (C) 2005-2012, U. of Colorado
 *****************************************************************************/

// -- CLASS DESCRIPTION [PDF] --
// This class derived from RooMinuit overloads the contour() method to
// produce a RooPlot.

#ifndef RAR_MINUIT
#define RAR_MINUIT

#include "TObject.h"
class TGraph;
class TH2F ;  // Needed because missing from RooMinuit.rdl
#include "RooMinuit.h"

class RooAbsReal ;
class RooRealVar ;
class RooPlot ;

class rarMinuit : public RooMinuit {
public:

  rarMinuit(RooAbsReal& function) ;
  RooPlot* contour(RooRealVar& var1, RooRealVar& var2,
		     Double_t n1=1, Double_t n2=2, Double_t n3=0,
		     Double_t n4=0, Double_t n5=0, Double_t n6=0);
  void fixGraph(TGraph *graph, Int_t lineStyle=1);
  
protected:

private:

  RooArgList* _floatParamList ;
  RooAbsReal* _func ;

protected:

  ClassDef(rarMinuit,0)   // RooMinuit derivative with contour RooPlot
} ;

#endif

