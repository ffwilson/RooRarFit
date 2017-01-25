/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarKeys.cc,v 1.8 2014/09/14 17:33:52 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides 1D/2D Keys Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides 1D/2D Keys Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "Roo2DKeysPdf.h"
#include "RooKeysPdf.h"

#include "rarKeys.hh"

ClassImp(rarKeys)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarKeys::rarKeys()
  : rarBasePdf(),
    _x(0), _y(0), _rho(1.), _keysOption("")
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
rarKeys::rarKeys(const char *configFile, const char *configSec,
		 const char *configStr,
		 rarDatasets *theDatasets, RooDataSet *theData,
		 const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr,
	       theDatasets, theData, name, title),
    _x(0), _y(0), _rho(1.), _keysOption("")
{
  init();
}

rarKeys::~rarKeys()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
///
/// There is no need to do pdfFit,
/// so the first step is to set #_pdfFit to false.
/// It reads in #_rho and #_keysOption,
/// and finally it builds RooKeysPdf/Roo2DKeysPdf PDF
/// with #_pdfType being Keys/2DKeys, respectively.
void rarKeys::init()
{
  cout<<"init of rarKeys for "<<GetName()<<":"<<endl;
  
  // no need for pdf fit
  setControlBits("noPdfFit");
  
  // get its obs
  _x=createAbsReal("x", "observable"); assert(_x);
  if ("2DKeys"==_pdfType) {
    _y=createAbsReal("y", "observable"); assert(_y);
  }
  // it should always have fitData
  if (!_theData) {
    cout<<" No dataset for Keys Pdf"<<endl;
    exit(-1);
  }
  // read in options
  _rho=atof(readConfStr("rho", "1.", getVarSec()));
  _keysOption=readConfStr("keysOption", "", getVarSec());
  
  // create pdf
  if ("2DKeys"==_pdfType) {
    _thePdf=new Roo2DKeysPdf(Form("the_%s", GetName()),_pdfType+" "+GetTitle(),
			     *_x, *_y, *_theData, _keysOption, _rho);
  } else { // default
    // parse the option
    RooKeysPdf::Mirror mirror=RooKeysPdf::NoMirror;
    if ("NoMirror"==_keysOption) mirror=RooKeysPdf::NoMirror;
    else if ("MirrorLeft"==_keysOption) mirror=RooKeysPdf::MirrorLeft;
    else if ("MirrorRight"==_keysOption) mirror=RooKeysPdf::MirrorRight;
    else if ("MirrorBoth"==_keysOption) mirror=RooKeysPdf::MirrorBoth;
    else if ("MirrorAsymLeft"==_keysOption) mirror=RooKeysPdf::MirrorAsymLeft;
    else if ("MirrorAsymLeftRight"==_keysOption)
      mirror=RooKeysPdf::MirrorAsymLeftRight;
    else if ("MirrorAsymRight"==_keysOption)mirror=RooKeysPdf::MirrorAsymRight;
    else if ("MirrorLeftAsymRight"==_keysOption)
      mirror=RooKeysPdf::MirrorLeftAsymRight;
    else if ("MirrorAsymBoth"==_keysOption) mirror=RooKeysPdf::MirrorAsymBoth;
    else mirror=RooKeysPdf::NoMirror;
    _thePdf=new RooKeysPdf(Form("the_%s", GetName()), _pdfType+" "+GetTitle(),
			   *_x, *_theData, mirror, _rho);
  }
}

/// \brief Load Keys dataset
///
/// \param theData The dataset to be loaded
///
/// Load \p theData as Keys dataset for pdf modeling.
/// There is no need to do pdfFit,
/// instead, 1/2D Keys PDFs call their LoadDataSet/loadDataSet functions
/// to load initial values (datasets).
void rarKeys::setFitData(RooDataSet *theData)
{
  RooDataSet *oldData=_theData;
  // set dataset
  rarBasePdf::setFitData(theData);
  if (_theData==oldData) return;
  // recalculate keys pdf
  if ("Keys"==_pdfType) ((RooKeysPdf*)_thePdf)->LoadDataSet(*_theData);
  if ("2DKeys"==_pdfType)
    ((Roo2DKeysPdf*)_thePdf)->loadDataSet(*_theData, _keysOption);
}
