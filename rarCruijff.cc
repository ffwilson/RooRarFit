/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarCruijff.cc,v 1.5 2014/09/14 17:33:50 fwilson Exp $
 * Authors: Karsten Koeneke, Lei Zhang
 * History:
 * 
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides the Cruijff Pdf class for RooRarFit
//
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides the Cruijff Pdf class for RooRarFit
//
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

// Include the RooFit Pdf
#include "RooCruijff.hh"

// Include this classes header
#include "rarCruijff.hh"

ClassImp(rarCruijff)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarCruijff::rarCruijff()
  : rarBasePdf(),
    _x(0), _m0(0), _sigmaL(0), _sigmaR(0), _alphaL(0), _alphaR(0)
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
rarCruijff::rarCruijff(const char *configFile, const char *configSec,
		       const char *configStr,
		       rarDatasets *theDatasets, RooDataSet *theData,
		       const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _m0(0), _sigmaL(0), _sigmaR(0), _alphaL(0), _alphaR(0)
{
  init();
}

rarCruijff::~rarCruijff()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates the parameters by calling #createAbsReal,
/// and finally it builds the User-defined PDF.
void rarCruijff::init()
{
  cout<<"init of rarCruijff for "<<GetName()<<":"<<endl;
  
  // first get its obs
  _x=createAbsReal("x", "observable"); assert(_x);
  // Config pdf params
  // for example, a is now mean, b sigma, etc.
  // [myPdf Config]
  // configStr = UsrPdf
  // x = AbsReal Def
  // mean = AbsReal Def
  // sigmaL = AbsReal Def
  // sigmaR = AbsReal Def
  // alphaL = AbsReal Def
  // alphaR = AbsReal Def
  
  // Default param creation
  _m0     = createAbsReal("mean", "#mu", 0, -10, 10);
  _sigmaL = createAbsReal("sigmaL", "#sigma_{L}", 0, -10, 10);
  _sigmaR = createAbsReal("sigmaR", "#sigma_{R}", 0, -10, 10);
  _alphaL = createAbsReal("alphaL", "#alpha_{L}", 0, -10, 10);
  _alphaR = createAbsReal("alphaR", "#alpha_{R}", 0, -10, 10);
  _params.Print("v");

  // create pdf
  _thePdf=new RooCruijff(Form("the_%s", GetName()),_pdfType+" "+GetTitle(),
			 *_x, *_m0, *_sigmaL, *_sigmaR, *_alphaL, *_alphaR);

}
