/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarBallack.cc,v 1.5 2014/09/14 17:33:48 fwilson Exp $
 * Authors: Karsten Koeneke, Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, Massachusetts Institute of Technology, Cambridge
 *  and          2005 University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides Ballack Pdf class for RooRarFit
// It is a Novosibirsk core with a polynomial attached
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides Ballack Pdf class for RooRarFit
// It is a Novosibirsk core with a polynomial attached
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "RooBallack.hh"

#include "rarBallack.hh"

ClassImp(rarBallack)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarBallack::rarBallack()
  : rarBasePdf(),
    _x(0), _mean(0), _width(0), _tail(0), _alpha(0), _n(0)
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
rarBallack::rarBallack(const char *configFile, const char *configSec,
		       const char *configStr,
		       rarDatasets *theDatasets, RooDataSet *theData,
		       const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _mean(0), _width(0), _tail(0), _alpha(0), _n(0)
{
  init();
}

rarBallack::~rarBallack()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates the parameters by calling #createAbsReal,
/// and finally it builds the RooBallack PDF.
void rarBallack::init()
{
  cout<<"init of rarBallack for "<<GetName()<<":"<<endl;
  
  // first get its obs
  _x=createAbsReal("x", "observable"); assert(_x);
  // Config pdf params
  _mean=createAbsReal("mean", "mean", 0, -10, 10);
  _width=createAbsReal("width", "width", 0, -10, 10);
  _tail=createAbsReal("tail", "tail", 0, -10, 10);
  _alpha=createAbsReal("alpha", "alpha", 0, -10, 10);
  _n=createAbsReal("n", "n", 0, -10, 10);
  _params.Print("v");
  
  // create pdf
  _thePdf=new RooBallack(Form("the_%s", GetName()),_pdfType+" "+GetTitle(),
			 *_x, *_mean, *_width, *_tail, *_alpha, *_n);
}
