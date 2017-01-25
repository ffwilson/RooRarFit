/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarGounarisSakurai.cc,v 1.3 2014/09/14 17:33:52 fwilson Exp $
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides GounarisSakurai Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides GounarisSakurai Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include <iostream>
#include <fstream>
using namespace std;

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "RooGounarisSakurai.hh"
#include "rarGounarisSakurai.hh"

ClassImp(rarGounarisSakurai);

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarGounarisSakurai::rarGounarisSakurai()
  : rarBasePdf(),
    _x(0), _mean(0), _width(0), _spin(0), 
    _radius(0), _mass_a(0), _mass_b(0)
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
rarGounarisSakurai::rarGounarisSakurai(const char *configFile, 
				     const char *configSec,
				     const char *configStr,
				     rarDatasets *theDatasets, 
				     RooDataSet *theData,
				     const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _mean(0), _width(0), _spin(0),
    _radius(0), _mass_a(0), _mass_b(0)
{
  init();
}

rarGounarisSakurai::~rarGounarisSakurai()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates the parameters by calling #createAbsReal,
/// and finally it builds the RooGounarisSakurai PDF.
void rarGounarisSakurai::init()
{
  
  // first get its obs
  _x=createAbsReal("x", "observable"); assert(_x);

  // Config pdf params
  _mean   = createAbsReal("mean", "mean", 0, -10, 10);
  _width  = createAbsReal("width", "width", 0, 0, 10);
  _spin   = createAbsReal("spin", "spin", 1, 0, 3);
  _radius = createAbsReal("radius", "radius", 3.1, 0, 10); // GeV-1
  _mass_a = createAbsReal("mass_a", "mass_a", 0.1359, 0, 10); // pi mass
  _mass_b = createAbsReal("mass_b", "mass_b", 0.1359, 0, 10);

  // set spin and masses to fixed ?
  _spin->setAttribute("Constant", kTRUE);
  _mass_a->setAttribute("Constant", kTRUE);
  _mass_b->setAttribute("Constant", kTRUE);

  _params.Print("v");
  
  // create pdf
  _thePdf=new RooGounarisSakurai(Form("the_%s", GetName()),_pdfType+" "+GetTitle(),
				 *_x, *_mean, *_width, *_spin, 
				 *_radius, *_mass_a, *_mass_b);
}
