/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarUsrPdf.cc,v 1.6 2014/09/14 17:33:56 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides User-defined Pdf class for RooRarFit
//
// Please change this cc file to make use of your PDF.
// In many cases, you just need to change the code inside marks:
//====================================================================>
//                                                                    v
//                                                                    v
//  and
//                                                                    ^
//                                                                    ^
//====================================================================>
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides User-defined Pdf class for RooRarFit
//
// Please change this cc file to make use of your PDF.
// In many cases, you just need to change the code inside marks:
//====================================================================>
//                                                                    v
//                                                                    v
//  and
//                                                                    ^
//                                                                    ^
//====================================================================>
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"
//using namespace std;

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

//==================================================================>
// Please include your RooFit Pdf header here                       v
// #include "mydir/myPdf.hh"                                        v

//                                                                  ^
//                                                                  ^
//==================================================================>


#include "rarUsrPdf.hh"

ClassImp(rarUsrPdf)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarUsrPdf::rarUsrPdf()
  : rarBasePdf(),
    _x(0), _a(0), _b(0), _c(0), _d(0), _e(0)
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
rarUsrPdf::rarUsrPdf(const char *configFile, const char *configSec,
		     const char *configStr,
		     rarDatasets *theDatasets, RooDataSet *theData,
		     const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _a(0), _b(0), _c(0), _d(0), _e(0)
{
  init();
}

rarUsrPdf::~rarUsrPdf()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates the parameters by calling #createAbsReal,
/// and finally it builds the User-defined PDF.
void rarUsrPdf::init()
{
  cout<<"init of rarUsrPdf for "<<GetName()<<":"<<endl;
  
  //==================================================================>
  // change the lines in between the marks as you want                v
  //                                                                  v
  // first get its obs
  _x=createAbsReal("x", "observable"); assert(_x);
  // Config pdf params
  // instead of a b c etc, you can give them more meaningful names and titles
  // for example
  // _a=createAbsReal("mean", "#mu", 0, -10, 10);
  // _b=createAbsReal("sigma", "#sigma", 0, -10, 10);
  // if you give them different names, please use those names
  // in the PDF config sections,
  // for example, a is now mean, b sigma, etc.
  // [myPdf Config]
  // configStr = UsrPdf
  // x = AbsReal Def
  // mean = AbsReal Def
  // sigma = AbsReal Def
  
  // Default param creation
  _a=createAbsReal("a", "a", 0, -10, 10);
  _b=createAbsReal("b", "b", 0, -10, 10);
  _c=createAbsReal("c", "c", 0, -10, 10);
  _d=createAbsReal("d", "d", 0, -10, 10);
  _e=createAbsReal("e", "e", 0, -10, 10);
  _params.Print("v");
  
  // YOU MUST CREATE YOUR PDF AND SET IT TO _thePdf
  // create pdf
  //_thePdf=new myPdf(Form("the_%s", GetName()),_pdfType+" "+GetTitle(),
  // *_x, *_a, *_b, *_c, *_d, *_e);
  //                                                                  ^
  // change the lines in between the marks as you want                ^
  //==================================================================>

}
