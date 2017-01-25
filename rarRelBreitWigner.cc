/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarRelBreitWigner.cc,v 1.6 2014/09/24 17:09:29 fwilson Exp $
 * Authors: Fergus Wilson
 * History:
 * 
 * Copyright (C) 2005-2012, RAL
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides RelBreitWigner Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides RelBreitWigner Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "RooRelBreitWigner.hh"
#include "rarRelBreitWigner.hh"

ClassImp(rarRelBreitWigner);

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarRelBreitWigner::rarRelBreitWigner()
  : rarBasePdf(),
    _x(0), _mean(0), _width(0), _radius(0), _mass_a(0), _mass_b(0), _spin(0)
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
rarRelBreitWigner::rarRelBreitWigner(const char *configFile, 
				     const char *configSec,
				     const char *configStr,
				     rarDatasets *theDatasets, 
				     RooDataSet *theData,
				     const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _mean(0), _width(0),  _radius(0), _mass_a(0), _mass_b(0), _spin(0)
{
  init();
}

rarRelBreitWigner::~rarRelBreitWigner()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates the parameters by calling #createAbsReal,
/// and finally it builds the RooRelBreitWigner PDF.
void rarRelBreitWigner::init()
{
  // read in _parSpin from config section
  //_parSpin=atoi(readConfStr("spin", Form("%d", _parSpin), getVarSec()));

  cout <<"init of rarRelBreitWigner for "<<GetName()<<endl;
  
  // first get its obs
  _x=createAbsReal("x", "observable"); assert(_x);

  // Config pdf params
  _mean=createAbsReal("mean", "mean", 0, -10, 10);
  _width=createAbsReal("width", "width", 0, 0, 100);
  //  _radius=createAbsReal("radius", "radius", 3.1, "RooConstVar");
  _radius=createAbsReal("radius", "radius", 3.1, 0, 10.0);
  _mass_a=createAbsReal("mass_a", "mass_a", 0.493677, "RooConstVar");
  _mass_b=createAbsReal("mass_b", "mass_b", 0.139570, "RooConstVar");
  _spin=createAbsReal("spin", "spin", 1, "RooConstVar");

  _params.Print("v");
  
  // create pdf
  _thePdf=new RooRelBreitWigner(Form("the_%s", GetName()),_pdfType+" "+GetTitle(),
				*_x, *_mean, *_width, *_radius, *_mass_a, *_mass_b, *_spin);
}
