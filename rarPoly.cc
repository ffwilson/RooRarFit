/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarPoly.cc,v 1.7 2014/09/14 17:33:54 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides Polynomial/Chebychev Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides Polynomial/Chebychev Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooFormulaVar.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "RooChebychev.h"
#include "RooPolynomial.h"

#include "rarPoly.hh"

ClassImp(rarPoly)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarPoly::rarPoly()
  : rarBasePdf(),
    _x(0), _nOrder(1)
{
  init();
}

/// \brief Default ctor
///
/// \param configFile The config file
/// \param configSec The config section
/// \param configStr The config string
/// \param theDatasets Available datasets
/// \param theData Default dataset for this PDF
/// \param name The name
/// \param title The title
///
/// The default ctor first initializes data members,
/// and then calls #init.
rarPoly::rarPoly(const char *configFile, const char *configSec,
		 const char *configStr, rarDatasets *theDatasets,
		 RooDataSet *theData, const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _nOrder(1)
{
  init();
}

rarPoly::~rarPoly()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates the parameters by calling #createAbsReal,
/// and finally it builds RooPolynomial/RooChebychev PDF
/// with #_pdfType being Polynomial/Chebychev, respectively.
void rarPoly::init()
{
  cout<<"init of rarPoly for "<<GetName()<<":"<<endl;
  
  // first get its obs
  _x=createAbsReal("x", "observable"); assert(_x);
  // read in pdf params
  _nOrder=atoi(readConfStr("nOrder", "1", getVarSec()));
  if (_nOrder<=0) _nOrder=1;
  _nOrder++;
  for (Int_t i=1; i<_nOrder; i++) {
    RooAbsReal *P=createAbsReal(Form("P%02d",i),Form("p_{%d}", i),
				0, -10000, +10000);
    _coeffs.add(*P);
  }
  _params.Print("v");
  
  // create pdf
  if ("Chebychev"==_pdfType) {
    _thePdf=new RooChebychev(Form("the_%s", GetName()),_pdfType+" "+GetTitle(),
			     *_x, _coeffs);
  } else { // default
    _thePdf=new RooPolynomial(Form("the_%s",GetName()),_pdfType+" "+GetTitle(),
			      *_x, _coeffs);
  }
}
