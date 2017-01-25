/*****************************************************************************
* Package: RooRarFit
 *    File: $Id: RooBallack.rdl,v 1.4 2014/09/14 17:33:46 fwilson Exp $     *
 * Authors:                                                                  *
 *    Karsten Koeneke, Massachusetts Institute of Technology, Cambridge, USA *
 *                                                                           *
 * Copyright (C) 2005-2012, Massachsetts Institute of Technology, Cambridge, USA  *
 *****************************************************************************/

#ifndef ROO_BALLACK
#define ROO_BALLACK

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;

class RooBallack : public RooAbsPdf {
public:
  RooBallack(const char *name, const char *title, 
	     RooAbsReal& _x,
	     RooAbsReal& _mean,
	     RooAbsReal& _width,
	     RooAbsReal& _tail,
	     RooAbsReal& _alpha,
	     RooAbsReal& _n);
  
  RooBallack(const RooBallack& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const { 
    return new RooBallack(*this,newname); }

  inline virtual ~RooBallack() { }

protected:
  RooRealProxy x;
  RooRealProxy mean;
  RooRealProxy width;
  RooRealProxy tail;
  RooRealProxy alpha;
  RooRealProxy n;

  Double_t evaluate() const;

private:
  ClassDef(RooBallack,0)
};

#endif
