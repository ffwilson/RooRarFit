/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarDecay.cc,v 1.13 2014/09/14 17:33:50 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides BDecay/Decay model for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides BDecay/Decay model for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "RooBCPGenDecay.h"
#include "RooBDecay.h"
#include "RooDecay.h"

#include "rarDecay.hh"

ClassImp(rarDecay)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarDecay::rarDecay()
  : rarBasePdf(),
    _x(0), _tau(0), _model(0), _decayType("DoubleSided"),
    _dm(0), _dgamma(0), _tag(0), _w(0), _dw(0), _mu(0),
    _Cb(0), _Sb(0),
    _f0(0), _f1(0), _f2(0), _f3(0)
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
rarDecay::rarDecay(const char *configFile, const char *configSec,
		   const char *configStr,
		   rarDatasets *theDatasets, RooDataSet *theData,
		   const char *name, const char *title)
  : rarBasePdf(configFile, configSec, configStr, theDatasets, theData,
	       name, title),
    _x(0), _tau(0), _model(0), _decayType("DoubleSided"),
    _dm(0), _dgamma(0), _tag(0), _w(0), _dw(0), _mu(0),
    _Cb(0), _Sb(0),
    _f0(0), _f1(0), _f2(0), _f3(0)
{
  init();
}

rarDecay::~rarDecay()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first creates the parameters by calling #createAbsReal,
/// and finally it builds RooArgusBG PDF.
void rarDecay::init()
{
  cout<<"init of rarDecay for "<<GetName()<<":"<<endl;
  
  // first get its obs
  _x=(RooRealVar*)createAbsReal("x", "observable"); assert(_x);
  RooRealVar *x=(RooRealVar *)RooArgList(_obsSet).at(0); assert(x);
  // make sure theDep is not derived
  if (_x!=x) {
    cout <<"No derived dependent allowed in rarDecay"<<endl;
    exit(-1);
  }
  // Config pdf params
  _tau=createAbsReal("tau", "#tau", 1.5, 0, 10, "ps");
  if ("BCPGenDecay"==_pdfType||"BDecay"==_pdfType) {
    _dm=createAbsReal("dm", "#Delta m", 0.502, 0, 1, "ps^{-1}");
    _dgamma=createAbsReal("dgamma", "#Delta#Gamma", 0, "RooConstVar");
    //_dgamma=createAbsReal("dgamma", "#Delta#Gamma", 0);
    _tag=(RooAbsCategory*)createAbsReal("tag", "tag");
    _w=createAbsReal("w", "w", 0.2, 0., 1.);
    _dw=createAbsReal("dw", "dw", 0., -2., 2.);
    _mu=createAbsReal("mu", "#mu", 0., -2., 2.);
    createAbsReal("S", "S_{CP}", .7, -2., 2.);
    createAbsReal("C", "C_{CP}", 0., -2., 2.);
    // get blind status and blind string
    _blindStatus=readConfStr("blindStatus", "blind", getVarSec());
    _blindString=readConfStr("blindString", "notSet", getVarSec());
    if ("notSet"==_blindString) {
      cout<<"You must specify blinding string in section \""
	  <<getVarSec()<<"\""<<endl;
      exit(-1);
    }
    // let's remove any "
    _blindString.ReplaceAll("\"", " ");
    // get blinding values
    TString bVStr=readConfStr("blindValues", ".2 .2 .2 .2", getVarSec());
    rarStrParser bVStrParser=bVStr;
    if (bVStrParser.nArgs()<4) {
      cout<<"Please specify 4 blind values"<<endl;
      exit(-1);
    }
    Double_t Cv=atof(bVStrParser[0]);
    Double_t Cs=atof(bVStrParser[1]);
    Double_t Sv=atof(bVStrParser[2]);
    Double_t Ss=atof(bVStrParser[3]);
    // create blind state cat
    RooCategory *blindCat=(RooCategory*)
      createAbsVar(Form("%s %s %s", "blindCat RooCategory",
			"\"CP blinding category\"",
			"useIdx    Unblind 0  Blind 1"));
    blindCat->setLabel("Blind");
    if ("unblind"==_blindStatus) blindCat->setLabel("Unblind");
    // put blindCat into ignore list
    addToConfStr("Ignored", blindCat->GetName(), getVarSec());
    // create blind variables for C and S
    _Cb=(RooAbsReal*)createAbsVar(Form("%s \"%s\" %f %f C blindCat kTRUE",
                                       "Cb RooUnblindPrecision \"C blind\"",
                                       _blindString.Data(), Cv, Cs));
    _Sb=(RooAbsReal*)createAbsVar(Form("%s \"%s\" %f %f S blindCat kTRUE",
                                       "Sb RooUnblindPrecision \"S blind\"",
                                       _blindString.Data(), Sv, Ss));
    if ("BDecay"==_pdfType) {
      cout<<" Due to bugs in PDF RooBDecay,"
	  <<" please use RooBCPGenDecay with"<<endl
	  <<" configStr = BCPGenDecay"<<endl
	  <<" in config section ["<<_configSec<<"]"<<endl;
      exit(-1);
      
      // creat f0, f1, f2, f3
      _f0=(RooAbsReal*)
        createAbsVar(Form("%s %s %s", "f0 RooFormulaVar",
                          "\"(1-@0*@2+@3*@0*(1.-2.*@1))\"", "tag w dw mu"));
      _f1=(RooAbsReal*)createAbsVar("f1 RooConstVar \"f1\" 0.");
      //_f1=createAbsReal("f1", "f1", 0);
      _f2=(RooAbsReal*)
        createAbsVar(Form("%s %s %s", "f2 RooFormulaVar",
                          "\"-1.*(@0*(1-2*@1)+@3*(1.-@0*@2))*@4\"",
                          "tag w dw mu Cb"));
      _f3=(RooAbsReal*)
        createAbsVar(Form("%s %s %s", "f3 RooFormulaVar",
                          "\"(@0*(1-2*@1)+@3*(1.-@0*@2))*@4\"",
                          "tag w dw mu Sb"));
    }
  }
  _params.Print("v");
  
  // resolution model
  _model=createPdfs("model", &_xPdfList);
  // let's make sure the default for model is noPdfFit noPdfPlot
  _model->setControlBit("noPdfFit", "pdfFit");
  _model->setControlBit("noPdfPlot", "pdfPlot");
    
  // parse the decay type
  RooBCPGenDecay::DecayType bcpGenDecayType=RooBCPGenDecay::DoubleSided;
  RooBDecay::DecayType bDecayType=RooBDecay::DoubleSided;
  RooDecay::DecayType decayType=RooDecay::DoubleSided;
  _decayType=readConfStr("decayType", "DoubleSided", getVarSec());
  if ("SingleSided"==_decayType) {
    decayType=RooDecay::SingleSided;
    bDecayType=RooBDecay::SingleSided;
    bcpGenDecayType=RooBCPGenDecay::SingleSided;
  } else if ("DoubleSided"==_decayType) {
    decayType=RooDecay::DoubleSided;
    bDecayType=RooBDecay::DoubleSided;
    bcpGenDecayType=RooBCPGenDecay::DoubleSided;
  } else if ("Flipped"==_decayType) {
    decayType=RooDecay::Flipped;
    bDecayType=RooBDecay::Flipped;
    bcpGenDecayType=RooBCPGenDecay::Flipped;
  }
  
  // create pdf
  if ("BCPGenDecay"==_pdfType) {
    _thePdf=new
      RooBCPGenDecay(Form("the_%s",GetName()), _pdfType+" "+GetTitle(),
                     *_x, *_tag, *_tau, *_dm, *_w, *_Cb, *_Sb, *_dw, *_mu,
                     *((RooResolutionModel*)_model->getPdf()),bcpGenDecayType);
  } else if ("BDecay"==_pdfType) {
    _thePdf=new RooBDecay(Form("the_%s",GetName()), _pdfType+" "+GetTitle(),
			  *_x, *_tau, *_dgamma, *_f0, *_f1, *_f2, *_f3, *_dm,
			  *((RooResolutionModel*)_model->getPdf()),bDecayType);
  } else {
    _thePdf=new RooDecay(Form("the_%s",GetName()), _pdfType+" "+GetTitle(),
			 *_x, *_tau,
			 *((RooResolutionModel*)_model->getPdf()), decayType);
  }
}

/// \brief Return prototype var generator for toy study
/// \return The generator
///
/// It constrcuts prototype var generator for toy study
RooAbsPdf *rarDecay::getProtGen()
{
  if (_theProtGen) return _theProtGen;
  // do we have model generator?
  if (_model) {
    RooAbsPdf *modelGen=_model->getProtGen();
    if (!_model->protGenIsDummy()) _protGenPdfs.add(*modelGen);
  }
  // call super class' getProtGen
  return rarBasePdf::getProtGen();
}
