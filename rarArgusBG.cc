/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarArgusBG.cc,v 1.6 2014/09/14 17:33:48 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides ArgusBG Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides ArgusBG Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "RooArgusBG.h"

#include "rarArgusBG.hh"

ClassImp(rarArgusBG)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarArgusBG::rarArgusBG()
  : rarBasePdf(),
    _x(0), _max(0), _c(0), _pow(0)
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
rarArgusBG::rarArgusBG(const char *configFile, const char *configSec,
		       const char *configStr,
		       rarDatasets *theDatasets, RooDataSet *theData,
		       const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _max(0), _c(0), _pow(0)
{
  init();
}

rarArgusBG::~rarArgusBG()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates variables by calling #createAbsReal,
/// and finally it builds RooArgusBG PDF.
void rarArgusBG::init()
{
  cout<<"init of rarArgusBG for "<<GetName()<<":"<<endl;
  
  // first get its dependent/observable
  _x=createAbsReal("x", "observable"); assert(_x);
  RooRealVar *x=(RooRealVar *)RooArgList(_obsSet).at(0); assert(x);
  // Config pdf params
  _max=createAbsReal("max", "E_{end}", x->getMax(), _x->getUnit());
  _c=createAbsReal("c", "#xi", -23, -80, -1);
  _pow=createAbsReal("pow", "n", 0.5);
  _params.Print("v");
  
  // create pdf
  _thePdf=new RooArgusBG(Form("the_%s",GetName()), _pdfType+" "+GetTitle(),
			 *_x, *_max, *_c, *_pow);
}
