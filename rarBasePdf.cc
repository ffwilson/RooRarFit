/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarBasePdf.cc,v 1.65 2014/09/14 17:33:49 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides a base Pdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides a base Pdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"
#include <fstream>
#include <sstream>
using namespace std;

#include "TArrayI.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TH1F.h"

#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooProdPdf.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooStringVar.h"
#include "RooSuperCategory.h"
#include "RooGlobalFunc.h"
using namespace RooFit;

#include "rarBasePdf.hh"
#include "rarMLFitter.hh"

ClassImp(rarBasePdf)
  ;

RooGenericPdf rarBasePdf::_dummyPdf("dummyPdf","dummy Pdf","1",RooArgList());
RooConstVar rarBasePdf::_dummyExpEvt("dummyExpEvt", "dummyExpEvt", 0);
RooExtendPdf rarBasePdf::_dummyExtPdf("dummyExtPdf","dummy extend Pdf",
			     rarBasePdf::_dummyPdf, rarBasePdf::_dummyExpEvt);
RooCategory rarBasePdf:: _compCat("compCat", "compCat");
rarMLFitter *rarBasePdf::_theFitter(0);
// plot colors
Int_t rarBasePdf::_rarColors[NCOLORS]={
  8, // light green
  6, // magenta
  9, // light blue
  49, // dark brown
  42, // brown
  //11, // light brown
  30, // dark green
  13, // dark gray
  7, // cyan
  2, // red
  3, // green
  5, // yellow
};
TString rarBasePdf::_fracNames;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarBasePdf::rarBasePdf()
  : rarConfig(),
    _pdfType("notSet"), _datasets(0), _theData(0),
    _thePdf(0), _theProtGen(0),
    _thisSimPdf(0), _thisSimPdfWOP(0), _theSimPdf(0),
    _nxPdf(0), _myDummyPdf(0), _controlStr("")
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
/// The default ctor initializes several common data members,
/// and then calls #init.
rarBasePdf::rarBasePdf(const char *configFile, const char *configSec,
		       const char *configStr,
		       rarDatasets *theDatasets, RooDataSet *theData,
		       const char *name, const char *title)
  : rarConfig(configFile, configSec, configStr, name, title),
    _pdfType("notSet"), _datasets(theDatasets), _theData(theData),
    _thePdf(0), _theProtGen(0),
    _thisSimPdf(0), _thisSimPdfWOP(0), _theSimPdf(0),
    _nxPdf(0), _myDummyPdf(0), _controlStr("")
{
  init();
}

rarBasePdf::~rarBasePdf()
{
  // _dataSets.Delete();
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first sets the pdf type using pdf info in #_configStr,
/// then it sets the param config section, #_paramSec,
/// to value in config item named <tt>paramSec_\<#_pdfType\></tt>,
/// if not found, sets to value of config item named \p paramSec,
/// if still not found,  #_paramSec is set to #_configSec.
/// All other config items should be in section #_paramSec.
/// #_fullObs is set to that found in #_datasets.
/// #_configStr is a string specifying the name, type and (optional)
/// title of the pdf. For example, for config section
/// \verbatim
/// [myPdf Config]
/// configStr=Gaussian "my Gaussian"
/// ...\endverbatim
/// The #_configStr will be
/// \verbatim
/// myPdf Gaussian "my Gaussian"\endverbatim
///
/// It also sets the default dataset by calling #setFitData,
/// and finally it builds extra pdfs,
/// if config \p xtraPdfs is specified.
void rarBasePdf::init()
{
  // set pdf type
  rarStrParser configStrParser=_configStr;
  configStrParser.Remove(); // name has been set
  _pdfType=configStrParser[0];
  configStrParser.Remove();
  // if this is the main fitter? if none set yet, it should be
  if (("MLFitter"==_pdfType)&&(!getFitter())) setFitter((rarMLFitter*)this);
  
  // set paramSec
  TString paramSec=readConfStr("paramSec_"+_pdfType, "notSet");
  if ("notSet"==paramSec) {
    if ("notSet"==(paramSec=readConfStr("paramSec", "notSet"))) {
      paramSec=_configSec;
    }
  }
  setVarSec(paramSec);
  cout<<"init of rarBasePdf:"<<endl
      <<" pdf Type: "<<_pdfType<<endl
      <<" pdfTitle: "<<GetTitle()<<endl
      <<" paramSec: "<<getVarSec()<<endl;
  
  // set control bits
  setControlBit("PdfFit", "pdfFit");
  setControlBit("PdfPlot", "pdfPlot");
  setControlBit("ParamsOnPlot", "paramsOnPlot");
  setControlBit("Chi2OnPlot", "chi2OnPlot");
  setControlBit("FirstFitOnly", "firstFitOnly");
  // default no
  setControlBit("noCompsOnPlot", "compsOnPlot");
  setControlBit("noCompsDataOnPlot", "compsDataOnPlot");
  // for SimPdf
  setControlBit("noSimFit", "simultaneousFit", getMasterSec());
  
  //setControlBit("noUseGenOnly", "useGeneratorOnly", getMasterSec());
  
  // get full obs
  _fullObs=_datasets->getFullFObs();
  // get the default dataset
  setFitData(_theData);
  
  // creates extra params if needed
  // The extra params created may or my not be directly used by this Pdf
  createAbsVars("xtraParams", &_xParams);
  
  // creates extra pdfs if needed
  // It creates extra PDFs associated with this RooRarFitPdf.
  // The extra pdfs themselves might not be direct components
  // of the (final) PDFs, but their parameters, components, etc.,
  // can be part of the (final) PDFs.
  createPdfs("xtraPdfs", &_xPdfList);
  // creates extra generators if needed
  createPdfs("xtraGenerators", &_xPdfList, &_protGenPdfs);
  
  // add protDataVars
  addProtVars();

}

/// \brief Add protDataVars (conditionalObs)
/// to #_protDataVars (#_conditionalObs)
///
/// It adds protDataVars (conditionalObs)
/// to #_protDataVars (#_conditionalObs)
void rarBasePdf::addProtVars()
{
  addProtVars("conditionalObs", _conditionalObs);
  addProtVars("protDataVars", _protDataVars);
  _protDataVars.add(_conditionalObs, kTRUE);
  // add _conditionalObs for pdfFit
  _condObss.add(_conditionalObs);
}

/// \brief Add obs in a config to argSet
/// \param configName Config name having obs
/// \param protVars Reference to argSet
///
/// It adds obs in given config to give argSet
void rarBasePdf::addProtVars(TString configName, RooArgSet &protVars)
{
  // get protDataVars defined in config file
  RooArgSet allProtVars(getArgSet(configName, kTRUE));
  allProtVars.add(getArgSet(readConfStr(configName,"", _runSec),kFALSE));
  if (allProtVars.getSize()>0) {
    //cout<<"protDataVars (conditionalObs) defined:"<<endl;
    //allProtVars.Print("v");
    RooArgList allProtVarList(allProtVars);
    for (Int_t i=0; i<allProtVarList.getSize(); i++) {
      TString varName=allProtVarList[i].GetName();
      RooAbsArg *theVar=_fullObs->find(varName);
      if (!theVar) {
        cout<<" Can not find obs named "<<varName<<endl;
      }
      protVars.add(*theVar);
    }
  }
}

/// \brief Set CondObss for conditional Pdf
///
/// \param condObsSet The CondObss for conditional Pdf
///
/// It sets the CondObss for conditional Pdf
void rarBasePdf::setCondObss(RooArgSet condObsSet)
{
  if (condObsSet.getSize()<=0) return;
  if (!_thePdf) return;
  if (!_theData) return;
  // first get all observables
  _condObss.add(*(_thePdf->getObservables(_theData)));
  // remove CondObss
  _condObss.remove(condObsSet);
  cout<<" conditionalObs for pdffit"<<endl;
  _condObss.Print("v");
}

/// \brief Set the default fit dataset.
///
/// \param theData The dataset
///
/// It sets the default dataset, #_theData, to its parameter,
/// and then it looks into param section to see if
/// \p fitData config item is defined, if so,
/// it sets the default dataset to the dataset defined in \p fitData.
void rarBasePdf::setFitData(RooDataSet *theData)
{
  _theData=theData;
  // check if fitData is defined in param section
  TString dsName=readConfStr("fitData", "notSet", getVarSec());
  if ("notSet"!=dsName) {
    RooDataSet *myData=_datasets->getData(dsName);
    if (myData) _theData=myData;
    else {
      cout<<" W A R N I N G !"<<endl
          <<" Dataset named "<<dsName<<" can not be found"<<endl
	  <<" fitData = "<<dsName<<endl
	  <<" defined in section "<<getVarSec()<<endl;
      _theData=0;
    }
  }
  
  if (_theData) {
    cout<<" default Dataset for "<<GetName()<<": "
	<<_theData->GetName()<<endl;
  } else {
    cout<<" no default dataset for "<<GetName()<<endl;
  }
  return;
}

/// \brief Return RooRealVar ArgSet associated with this RooRarFitPdf
/// \return params associated with this RooRarFitPdf
///
/// This function returns RooRealVar ArgSet associated with this RooRarFitPdf,
/// ie, #_params, #_params of its extra pdfs, and #_xParams.
RooArgSet rarBasePdf::getParams()
{
  RooArgSet theSet(_params);
  for (Int_t i=0; i<_nxPdf; i++) {
    rarBasePdf *thePdf=(rarBasePdf*)_xPdfList.At(i);
    theSet.add(thePdf->getParams());
  }
  theSet.add(_xParams);
  
  return theSet;
}

/// \brief Add param to #_params
///
/// \param theVar The var needs to be added
///
/// It adds \p theVar into #_params
void rarBasePdf::addToParams(RooRealVar *theVar)
{
  _params.add(*theVar);
}

/// \brief Add obs to #_obsSet
///
/// \param theVar The var needs to be added
///
/// It adds \p theVar into #_obsSet
void rarBasePdf::addToObs(RooRealVar *theVar)
{
  _obsSet.add(*theVar);
}

/// \brief Get RooFormulaVar ArgList
///
/// \param fStrParser The parsed tokens of formula Args
/// \return The formula ArgList
///
/// It creates a RooArgList from \p fStrParser.
/// It first checks if the token in \p fStrParser is in #_obsSet,
/// if not, checks if in #_params, then if in #_fullObs,
/// if not, it creats one by calling #createAbsReal,
/// and finally, it adds the var into the ArgList.
/// It repeats until all the tokens are scanned.
RooArgList *rarBasePdf::getFormulaArgs(rarStrParser fStrParser)
{
  RooArgList *depList=new RooArgList;
  Int_t nArgs=fStrParser.nArgs();
  if (nArgs<=0) return depList;
  for (Int_t i=0; i<nArgs; i++) {
    RooAbsArg *theDep(0);
    if (theDep=_obsSet.find(fStrParser[i]));
    else if (theDep=_params.find(fStrParser[i]));
    else if (theDep=_fullObs->find(fStrParser[i])) {
      _obsSet.add(*theDep);
    } else { // create it
      if (isNumber(fStrParser[i])) break;
      theDep=createAbsReal(fStrParser[i], fStrParser[i]);
    }
    depList->add(*theDep);
  }
  
  return depList;
}

/// \brief Return a Double_t value of a string as formula var
///
/// \param varStr The string to be evaluated
/// \return The evaluation
///
/// It calculates the value of \p varStr
/// which is parsed as formula and returns that number.
/// The formula string looks like
/// \li ``<tt>\<formula\> \<var1\> \<var2\> ...</tt>''
Double_t rarBasePdf::getFormulaVal(TString varStr)
{
  Double_t retVal(0);
  rarStrParser varStrParser=varStr;
  TString formula=varStrParser[0]; // get the formula
  varStrParser.Remove();
  // we need to prserve _params RooArgSet
  RooArgSet params(_params);
  RooFormulaVar myVar("myVar", formula, *getFormulaArgs(varStrParser));
  retVal=myVar.getVal();
  // restore original _params
  _params.removeAll();
  _params.add(params);
  return retVal;
}

/// \brief Creates pdfs according to component string
/// \param Comps A String of pdf names
/// \param pdfList A TList to store the created RooRarFitPdfs
/// \param PDFs A RooArgList to store the created RooFit Pdfs
/// \param secName Config section name
/// \return The last RooRarFitPdf created
///
/// This function creates all pdfs specified by a string
/// containing names of pdfs.
rarBasePdf *rarBasePdf::createPdfs(TString Comps, TList *pdfList,
				   RooAbsCollection *PDFs, TString secName)
{
  // secName
  if (""==secName) secName=getVarSec();
  // get number of components
  rarStrParser compStrParser=readConfStr(Comps, "", secName);
  Int_t nComp=compStrParser.nArgs();
  // read in comp info
  RooArgSet pdfCompSet("PDF Comp Set");
  for (Int_t i=0; i<nComp; i++) {
    TString configStr=
      readConfStr("configStr", "notSet", compStrParser[i]+" Config");
    RooStringVar *pdf=
      new RooStringVar(compStrParser[i], compStrParser[i], configStr);
    pdfCompSet.addOwned(*pdf);
  }
  //pdfCompSet.Print("v");
  // list configed pdfs
  if (nComp>0) {
    cout <<endl<<"Pdfs defined with config \""<<Comps<<"\""
	 <<" in config file for "<<GetName()<<":"<<endl;
  }
  for (Int_t i=0; i<nComp; i++) {
    TString pdfStr=(RooStringVar&)pdfCompSet[compStrParser[i]];
    cout<<Form(" #%02d ",i)<<compStrParser[i]<<" "<<pdfStr<<endl;
  }
  cout<<endl;
  // now build individual pdf
  rarBasePdf *thePdf(0);
  for (Int_t i=0; i<nComp; i++) {
    TString pdfStr=(RooStringVar&)pdfCompSet[compStrParser[i]];
    thePdf=createPdf(compStrParser[i]+" "+pdfStr);
    if (PDFs) PDFs->add(*thePdf->getPdf());
    if (pdfList) pdfList->Add(thePdf);
  }
  
  return thePdf;
}

/// \brief Actions right after every RooRarFitPdf is created
/// and before any other action is taken
///
/// It first calls \p preAction of all the RooRarFitPdfs in #_xPdfList.
/// It makes a list of full observables, #_fObsSet,
/// reads in common config options for every RooRarFitPdf,
/// \p pdfFit, \p pdfPlot,
/// \p paramsOnPlot, \p chi2OnPlot, etc.
/// Unlike config options in rarMLFitter::run function,
/// those options affect this RooRarFitPdf only.
/// It also builds a SimPdf of its own for pdfFit purpose
/// if the final pdf is a SimPdf (disabled for now).
void rarBasePdf::preAction()
{
  // first for xtraPdf
  _nxPdf=_xPdfList.GetSize();
  for (Int_t i=0; i<_nxPdf; i++) {
    rarBasePdf *thePdf=(rarBasePdf*)_xPdfList.At(i);
    thePdf->preAction();
  }
  cout<<" In rarBasePdf preAction for "<<GetName()<<endl;
  // first get full obs including those from sub pdfs
  _fObsSet.add(*_thePdf->getDependents(*_fullObs));
  // create _myDummyPdf
  if (!_myDummyPdf) {
    RooArgSet myFullVars(*_thePdf->getParameters(_theData));
    _myDummyPdf=new RooGenericPdf
      (Form("myDummyPdf_%s", GetName()), "1", myFullVars);
  }
  // build simPdf for this Pdf
  _thisSimPdf=(RooSimultaneous*)getSimPdf(_theSimPdf);
  _thisSimPdfWOP=_thisSimPdf;
}

/// \brief Return the created SimPdf for the PDF specified
/// \param simPdf The reference SimPdf
/// \param srcPdf The pdf to get SimPdf
/// \return The simPdf created
///
/// This function creates simPdf based on reference SimPdf for srcPdf
RooAbsPdf *rarBasePdf::getSimPdf(RooSimultaneous *simPdf, RooAbsPdf *srcPdf)
{
  if (!simPdf) return _thisSimPdf;
  RooSimultaneous *theSimPdf(0);
  if (TString(simPdf->ClassName())!="RooSimultaneous") return theSimPdf;
  if (!srcPdf) srcPdf=_thePdf;
  // simPdf for master pdf is itself
  if ((simPdf==_theSimPdf)&&(srcPdf==_theSimPdf)) return _theSimPdf;
  // build simPdf for the pdf based on full SimPdf
  RooAbsCategoryLValue *simCats=(RooAbsCategoryLValue *)(&simPdf->indexCat());
  //simCats->Print("v");
  // get pdf name
  theSimPdf=new RooSimultaneous(Form("simPdf_%s", srcPdf->GetName()),
				Form("simPdf for %s", srcPdf->GetTitle()),
				*simCats);
  Bool_t isFound=kFALSE;
  Int_t nCats=simCats->numTypes();
  for(Int_t i=0; i<nCats; i++) {
    RooCatType *catType=(RooCatType*)simCats->lookupType(i);
    TString catName=catType->GetName();
    TString catName2=catName;
    // do we have physCat?
    TString physCat=getFitter()->getPhysCat();
    if ((""!=physCat)&&(catName2.Contains(";"))) {
      TString physCatType=catName2(1,catName2.First(';')-1);
      catName2.ReplaceAll(physCatType+";", "");
      catName2.ReplaceAll("}", ";"+physCatType+"}");
    }
    // check if we can find cloned pdf
    RooAbsPdf *thisCatFullPdf=simPdf->getPdf(catName.Data());
    if (!thisCatFullPdf) {
      theSimPdf->addPdf(_dummyPdf, catName);
    } else {
      //thisCatFullPdf->Print("v");
      RooArgSet *thisCatCompSet=thisCatFullPdf->getComponents();
      RooAbsPdf *thisCatPdf=(RooAbsPdf*)thisCatCompSet->
	find(Form("%s_%s", srcPdf->GetName(), catName.Data()));
      if (!thisCatPdf) thisCatPdf=(RooAbsPdf*)thisCatCompSet->
	find(Form("%s_%s", srcPdf->GetName(), catName2.Data()));
      if (thisCatPdf) {
	theSimPdf->addPdf(*thisCatPdf, catName);
        isFound=kTRUE;
      } else {
        theSimPdf->addPdf(_dummyPdf, catName);
      }
    }
  }
  // do we really find sim pdf for this pdf
  if (!isFound) {
    // this Pdf does not have SimPdf counter-part
    //cout<<" No need to build simPdf for "<<GetName()<<endl;
    delete theSimPdf;
    theSimPdf=0;
  }
  if (theSimPdf) {
    cout<<" simPdf created for "<<srcPdf->GetName()<<endl;
    theSimPdf->Print();
    //theSimPdf->Print("v");
  }
  return theSimPdf;
}

/// \brief Constructs ArgSet of params with given names
///
/// \param paramNames The param names or config item name
/// \param useRead To use #readConfStr or not
/// \param fullSet Argset to get additional args
/// \return The ArgSet created
///
/// It constructs ArgSet according to \p paramNames from its #_params.
/// If \p useRead is set, \p paramNames is actually config item name
/// and #readConfStr is used to read in the actaully names.
/// If a name is actually the name of this RooRarFitPdf itself,
/// it adds all params in #_params,
/// if not, it checks if the name is in #_params,
/// if yes, it adds the param into the ArgSet.
/// It loops through all the names and return the created ArgSet.
RooArgSet rarBasePdf::getArgSet(TString paramNames, Bool_t useRead,
				RooArgSet *fullSet)
{
  RooArgSet retVal;
  rarStrParser paramNameStrParser="";
  if (!useRead) {
    paramNameStrParser=paramNames;
  } else {
    paramNameStrParser=readConfStr(paramNames, "", getVarSec());
  }
  Int_t nParams=paramNameStrParser.nArgs();
  RooAbsArg *theParam(0);
  rarBasePdf *thePdf(0);
  for (Int_t i=0; i<nParams; i++) {
    TString paramName=paramNameStrParser[i];
    // if fullSet have the name
    if (fullSet) {
      RooArgSet *inFS=(RooArgSet*)fullSet->selectByName(paramName);
      retVal.add(*inFS);
      delete inFS;
    }
    // check if the paramName is the pdf Name
    if(paramName==GetName()) {
      retVal.add(_params);
      // have any split param?
      if (fullSet) {
	RooArgList paramList(_params);
	for (Int_t j=0; j<paramList.getSize(); j++) {
          RooArgSet *inFS=(RooArgSet*)
            fullSet->selectByName(Form("%s_*",paramList[j].GetName()));
          retVal.add(*inFS);
          delete inFS;
        }
      }
    } else if (theParam=getAbsVar(paramName)) {
      //RooAbsArg *theParam=_params.find(paramName);
      retVal.add(*theParam);
      if (fullSet) {
        RooArgSet *inFS=(RooArgSet*)
          fullSet->selectByName(Form("%s_*",theParam->GetName()));
        retVal.add(*inFS);
        delete inFS;
      }
    } else if (thePdf=(rarBasePdf*)(_rarPdfs.FindObject(paramName))) {
      retVal.add(thePdf->getArgSet(thePdf->GetName(), kFALSE, fullSet));
    }
  }
  // now for xtraPdf
  for (Int_t i=0; i<_nxPdf; i++) {
    rarBasePdf *thePdf=(rarBasePdf*)_xPdfList.At(i);
    retVal.add(thePdf->getArgSet(paramNames, useRead, fullSet));
  }
  return retVal;
}

/// \brief Check if a new catType is part of old catType
/// \param catTypeN New catType
/// \param catTypeO Old catType
/// \return If catTypeN matches part of catTypeO, true
Bool_t rarBasePdf::matchCatType(RooCatType *catTypeN, RooCatType *catTypeO)
{
  Bool_t retVal(kTRUE);
  TString typeNameO=catTypeO->GetName();
  typeNameO.ReplaceAll("{", "");
  typeNameO.ReplaceAll("}", "");
  typeNameO.ReplaceAll(";", " ");
  rarStrParser typeParserO=typeNameO;
  
  TString typeNameN=catTypeN->GetName();
  typeNameN.ReplaceAll("{", "");
  typeNameN.ReplaceAll("}", "");
  typeNameN.ReplaceAll(";", " ");
  rarStrParser typeParserN=typeNameN;
  for (Int_t i=0; i<typeParserN.nArgs(); i++) {
    if (!typeParserO.Have(typeParserN[i])) {
      retVal=kFALSE;
      break;
    }
  }
  
  return retVal;
}

/// \brief Save frac name for splitting
/// \param fracName The frac name to save
///
/// It saves frac name for splitting
void rarBasePdf::saveFracName(TString fracName)
{
  if (isFracName(fracName)) return;
  _fracNames+=" "+fracName;
  //cout<<_fracNames<<endl;
}

/// \brief Check if a name is frac name or not
/// \param fracName The name to check
/// \return kTRUE if the name is frac Name
///
/// It checks if a name is frac name or not (for splitting)
Bool_t rarBasePdf::isFracName(TString fracName)
{
  Bool_t retVal(kFALSE);
  rarStrParser fracNamesParser=_fracNames;
  if (fracNamesParser.Have(fracName)) retVal=kTRUE;
  
  return retVal;
}

/// \brief Attach dataset to the pdf
/// \param data The dataset to be attached to the Pdf
///
/// It attaches dataset to this PDF
void rarBasePdf::attachDataSet(const RooAbsData &data)
{
  RooAbsPdf *thePdf(0);
  if (_thePdf) thePdf=_thePdf;
  if (_thisSimPdf) thePdf=_thisSimPdf;
  if (!thePdf) return;
  thePdf->attachDataSet(data);
}

/// \brief Check if the current value of PDF is negative
/// \return true if negative
Bool_t rarBasePdf::isNegativeValue()
{
  Bool_t isNegative(kFALSE);
  RooAbsPdf *thePdf(0);
  if (_thePdf) thePdf=_thePdf;
  if (_thisSimPdf) thePdf=_thisSimPdf;
  if (!thePdf) return isNegative;
  if (thePdf->getVal()<0) {
    isNegative=kTRUE;
    cout<<" Pdf "<<thePdf->GetName()<<" has negative value"<<endl;
    thePdf->Print();
    thePdf->Print("v");
    return isNegative;
  }
  return isNegative;
}


/// \brief Get var ranges from config file
/// \param theVar The var to get ranges from
/// \param rPrefix The config prefix for ranges
/// \param min The double var to store min
/// \param max The double var to store max
/// \param sec Config section name
/// \param nBins Number of bins
/// \return The RooBinning var if created
///
/// It returns min and max for a var from config file.
/// It nBins is set, it will return RooBinning.
RooBinning *rarBasePdf::getRange(RooRealVar *theVar, TString rPrefix,
                                 Double_t &min, Double_t &max,
                                 const Char_t *sec, Int_t *nBins)
{
  RooBinning *abins(0);
  if (!sec) sec=getVarSec();
  min=theVar->getMin();
  max=theVar->getMax();
  rarStrParser rangeParser=
    readConfStr(rPrefix+theVar->GetName(),"", sec);
  if (rangeParser.nArgs()<=0) { // try runSec then
    rangeParser=readConfStr(rPrefix+theVar->GetName(),"", _runSec);
  }
  if (rangeParser.nArgs()>0) {
    min=atof(rangeParser[0]);
    rangeParser.Remove();
    if (min<theVar->getMin()) min=theVar->getMin();
  }
  if (rangeParser.nArgs()>0) {
    max=atof(rangeParser[rangeParser.nArgs()-1]);
    rangeParser.Remove(rangeParser.nArgs()-1);
    if (max>theVar->getMax()) max=theVar->getMax();
  }
  if (max<min) {
    Double_t v=max;
    max=min;
    min=v;
  }
  if (nBins) {
    abins=new RooBinning(*nBins, min, max);
    if (rangeParser.nArgs()>0) {
      delete abins;
      abins=new RooBinning(min, max);
      while (rangeParser.nArgs()>0) {
        Double_t boundary=atof(rangeParser[0]);
        rangeParser.Remove();
        if ((boundary<min) || (boundary>max)) continue;
        abins->addBoundary(boundary);
      }
    }
  }
  return abins;
}

/// \brief RooRealVar ArgSet of correlation coefficient of params
///   associated with this RooRarFitPdf.
/// \return Correlation coefficients
///
/// This function returns RooRealVar ArgSet of correlation coefficient
/// of params associated with this RooRarFitPdf.
RooArgSet rarBasePdf::getCorrCoeffs()
{
  RooArgSet theSet(_corrCoeffs);
  for (Int_t i=0; i<_nxPdf; i++) {
    rarBasePdf *thePdf=(rarBasePdf*)_xPdfList.At(i);
    theSet.add(thePdf->getCorrCoeffs());
  }
  
  return theSet;
}

/// \brief Return correlation coefficient of two params
/// \param pn1 Name of the first param
/// \param pn2 Name of the second param
/// \return Correlation coefficient
///
/// This function returns correlation coefficient
/// of two params specified in args
Double_t rarBasePdf::getCorrCoeff(const TString pn1, const TString pn2)
{
  // trivial case, but need to be considered specially
  if (pn1==pn2) return 1;
  // do we have it in argset?
  RooArgSet theCorrSet(getCorrCoeffs());
  return theCorrSet.getRealValue(getCorrCoefName(pn1, pn2), 0);
}

/// \brief Save the correlation coefficients for error calculation
/// \param fr Fit results
/// 
/// Save the correlation coefficients for (systematic) error calculation
void rarBasePdf::saveCorrCoeffs(RooFitResult *fr)
{
  // add correlation matrix
  Int_t nFloat=fr->floatParsFinal().getSize();
  for (Int_t iF=0; iF<nFloat; iF++) {
    RooAbsArg &iArg=fr->floatParsFinal().operator[](iF);
    for(Int_t jF=0; jF<iF; jF++) {
      RooAbsArg &jArg=fr->floatParsFinal().operator[](jF);
      TString corrName=getCorrCoefName(iArg.GetName(), jArg.GetName());
      saveCorrCoeff(corrName, fr->correlation(iArg, jArg));
    }
  }
}

/// \brief Save the correlation coefficient for error calculation
/// \param corrName Correlation coefficient name
/// \param corrCoef Correlation coefficient value
/// \param saveTrivial Save trivial values
/// \return true if successful
/// 
/// Save the correlation coefficient for (systematic) error calculation
Bool_t rarBasePdf::saveCorrCoeff(TString corrName, Double_t corrCoef,
				 Bool_t saveTrivial)
{
  Bool_t retVal(kFALSE);
  if (fabs(corrCoef)<1e-4) {
    if (!saveTrivial)
      return retVal;
  }
  if (fabs(corrCoef)>1) return retVal;
  if (_corrCoeffs.setRealValue(corrName, corrCoef)) {
    // create corr coeff RRV
    RooRealVar *theCorrCoefVar=new RooRealVar(corrName, corrName, corrCoef);
    _corrCoeffs.add(*theCorrCoefVar);
  }
  retVal=kTRUE;
  
  return retVal;
}

/// \brief Get correlation coefficient name of two params
/// \param pn1 Name of the first param
/// \param pn2 Name of the second param
/// \return Correlation coefficient Name
TString rarBasePdf::getCorrCoefName(const TString pn1, const TString pn2) const
{
  // construct corrcoef var name
  TString corrName=Form("corrCoef_%s_%s", pn1.Data(), pn2.Data());
  if (pn1.CompareTo(pn2)>0)
    corrName=Form("corrCoef_%s_%s", pn2.Data(), pn1.Data());
  return corrName;
}

/// \brief Pdf fit for extra Pdfs
/// \param pdfList Pdfs need to do pdfFit
///
/// Run pdfFit for extra pdfs.
void rarBasePdf::doXPdfFit(TString pdfList)
{
  for (Int_t i=0; i<_nxPdf; i++) {
    rarBasePdf *thePdf=(rarBasePdf*)_xPdfList.At(i);
    thePdf->doPdfFit(pdfList);
  }
}

/// \brief Do pdfFit for given PDFs
///
/// \param pdfList Pdfs need to do pdfFit
///
/// It first calls #doXPdfFit for extra pdfs.
/// Then it checks various booleans, etc., to see if it needs to do pdfFit
/// for itself. If \p pdfList is not null it will run pdfFit only
/// when its name is in \p pdfList.
/// It fixes any params in #_prePdfFixParamSet and calls
/// \p fitTo of RooFit to do the fit.
/// If there is SimPdf for this pdf, it also fits #_thisSimPdf.
/// After fitting, it will store the fit params in a string,
/// #_afterFitSaverStr, for later use, like plotting, etc.
void rarBasePdf::doPdfFit(TString pdfList)
{
  // first for xtraPfs
  doXPdfFit(pdfList);
  
  cout<<endl<<" In rarBasePdf doPdfFit for "<<GetName()<<endl;
  
  if (!getControlBit("PdfFit")) return;
  TString pdfFitStr=readConfStr("pdfFit","yes",getVarSec());
  // is fit been done?
  if (getControlBit("FirstFitOnly")) {
    if (getControlBit("FirstFitDone")) return;
  }
  // are we in the list for plotting
  rarStrParser pdfListParser=pdfList;
  if (pdfListParser.nArgs()>0)
    if (!pdfListParser.Have(GetName())) return;
  
  // first check if it has its pdf and data
  assert(_thePdf);
  assert(_theData);
  // display dataset to fit
  cout<<" pdfFit to dataset: "<<_theData->GetName()<<endl;
  // fitoption
  TString fitOption="hrq";
  //if (_thePdf->isExtended()) fitOption="emhr";
  // convert to newer fitTo format for steering options
  Bool_t pdfExtended = fitOption.Contains("e");
  Bool_t pdfMinos    = fitOption.Contains("m");
  Bool_t pdfHesse    = fitOption.Contains("h");
  Bool_t pdfVerbose  = !fitOption.Contains("q");
  Bool_t pdfSave     = fitOption.Contains("r"); // return results
 
  // get number of cpus option
  Int_t pdfFitNumCPU=atoi(readConfStr("useNumCPU", "1", getMasterSec()));

  // save obs
  string obsSaveStr;
  writeToStr(_fObsSet, obsSaveStr);
  // set fit ranges
  RooArgList fObsList(_fObsSet);
  for (Int_t i=0; i<fObsList.getSize(); i++) {
    // get obs
    RooRealVar *theObs=dynamic_cast<RooRealVar*>(fObsList.at(i));
    if (!theObs) continue;
    Double_t fitMin=theObs->getMin();
    Double_t fitMax=theObs->getMax();
    getRange(theObs, "fitRange_", fitMin, fitMax);
    theObs->setRange(fitMin, fitMax);
  }
  
  // get all params for this pdf and simPdf
  RooArgSet thisFParams(*_thePdf->getParameters(_theData));
  if (_thisSimPdf) thisFParams.add(*_thisSimPdf->getParameters(_theData));
  // first fix params specified with prePdfFix
  RooArgSet prePdfFixParamSet=
    getArgSet(readConfStr("prePdfFix","", getVarSec()), kFALSE, &thisFParams);
  // set those params to values specified
  rarStrParser prePdfFixParser=readConfStr("prePdfFix", "", getVarSec());
  while (prePdfFixParser.nArgs()>0) {
    TString varName=prePdfFixParser[0];
    prePdfFixParser.Remove();
    RooRealVar *theVar=(RooRealVar*)prePdfFixParamSet.find(varName);
    if (!theVar) continue;
    // check if the next is number
    TString varVal=prePdfFixParser[0];
    if (isNumber(varVal)) {
      prePdfFixParser.Remove();
      if (theVar->ClassName()==TString("RooRealVar"))
	theVar->setVal(atof(varVal));
    }
  }
  if (prePdfFixParamSet.getSize()>0) {
    cout<<" prePdfFix:"<<endl;
    prePdfFixParamSet.Print("v");
  }
  // then float params specified w/ prePdfFloat
  RooArgSet prePdfFloatParamSet=
    getArgSet(readConfStr("prePdfFloat","", getVarSec()),kFALSE,&thisFParams);
  // merge those params to prePdfFix for saving
  prePdfFixParamSet.add(prePdfFloatParamSet);
  string prePdfFixParamSetSSaver;
  writeToStr(prePdfFixParamSet, prePdfFixParamSetSSaver);
  // fix or float
  prePdfFixParamSet.setAttribAll("Constant"); 
  prePdfFloatParamSet.setAttribAll("Constant", kFALSE);
  
  // fit to default dataset
  if ("simFitOnly"!=pdfFitStr) {
    cout<<endl<<" In rarBasePdf doPdfFit for "<<GetName()<< " Options: " << fitOption 
      << " using " << pdfFitNumCPU << " CPUs (if available)." << endl;
    RooFitResult *fitResult=
      _thePdf->fitTo(*_theData,ConditionalObservables(_condObss),
		     Save(pdfSave), Extended(pdfExtended), 
		     Verbose(pdfVerbose), Hesse(pdfHesse),
		     Minos(pdfMinos), NumCPU(pdfFitNumCPU));

    saveCorrCoeffs(fitResult);
  }
  // for simPdf
  if (_thisSimPdf&&pdfFitStr.Contains("simFit")) {
    cout<<"Should also have simPdf fit"<<endl;
    //// first fix all other params
    //string coeffSSaver;
    //writeToStr(_params, coeffSSaver);
    //_params.setAttribAll("Constant");
    /// \todo Create_thisSimPdfWOP through another function
    // let's find if there is any protCat?
    RooArgSet protDataEVars(getFitter()->getProtDataEVars());
    RooAbsCategoryLValue *simCats=(RooAbsCategoryLValue *)
      (&_thisSimPdf->indexCat());
    RooSuperCategory *sCatO=dynamic_cast<RooSuperCategory*>(simCats);
    if (sCatO) {
      RooArgSet inputCats(sCatO->inputCatList());
      Int_t nCatsO=inputCats.getSize();
      inputCats.remove(protDataEVars, kFALSE, kTRUE);
      if ((inputCats.getSize()!=nCatsO)&&(inputCats.getSize()>0)) {
        RooSuperCategory *sCatN=new
          RooSuperCategory("sCatN", "sCatN", inputCats);
        _thisSimPdfWOP=new RooSimultaneous
          (Form("%sWOP", _thisSimPdfWOP->GetName()), "simPDFWOP", *sCatN);
        Int_t nTypesN=sCatN->numTypes();
        Int_t nTypesO=sCatO->numTypes();
        for(Int_t i=0; i<nTypesN; i++) {
          RooCatType *catTypeN=(RooCatType*)sCatN->lookupType(i);
          for (Int_t j=0; j<nTypesO; j++) {
            RooCatType *catTypeO=(RooCatType*)sCatO->lookupType(j);
            if (matchCatType(catTypeN, catTypeO)) {
              RooAbsPdf *typePdf=_thisSimPdf->getPdf(catTypeO->GetName());
              if (TString(typePdf->GetName())==_dummyPdf.GetName()) continue;
              _thisSimPdfWOP->addPdf(*typePdf, catTypeN->GetName());
              break;
            }
          }
        }
      }
    }
    //RooFitResult *fitResult=
    //  _thisSimPdfWOP->fitTo(*_theData,_condObss,fitOption);
    RooFitResult *fitResult= _thisSimPdfWOP->fitTo(*_theData,ConditionalObservables(_condObss),
						   Save(pdfSave), Extended(pdfExtended), 
						   Verbose(pdfVerbose), Hesse(pdfHesse),
						   Minos(pdfMinos), NumCPU(pdfFitNumCPU));
    saveCorrCoeffs(fitResult);
    //// restore params
    //readFromStr(_params, coeffSSaver);
  }
  // after fit, save the params for later use, like plotting, etc
  RooArgSet *theParams=_thePdf->getParameters(_theData);
  writeToStr(*theParams, _afterFitSaverStr);
  // restore params in prePdfFixParamSet
  readFromStr(prePdfFixParamSet, prePdfFixParamSetSSaver);
  // restore obs saved
  readFromStr(_fObsSet, obsSaveStr);
  setControlBits("FirstFitDone");
}

/// \brief Pdf plot for extra Pdfs
///
/// \param plotList List of plots
/// \param pdfList List of PDFs to be plotted
/// \return The last RooPlot created
///
/// It loops through all extra PDFs to get pdf plots by
/// calling their #doPdfPlot.
RooPlot *rarBasePdf::doXPdfPlot(TList &plotList, TString pdfList)
{
  RooPlot *frame(0);
  for (Int_t i=0; i<_nxPdf; i++) {
    rarBasePdf *thePdf=(rarBasePdf*)_xPdfList.At(i);
    frame=thePdf->doPdfPlot(plotList, pdfList);
  }
  return frame;
}

/// \brief Do pdfPlot for given PDFs
///
/// \param plotList List of plots
/// \param pdfList List of PDFs to be plotted
/// \return The last RooPlot created
///
/// It first calls #doXPdfPlot to get PDF plots of extra PDFs.
/// Before it gets plot for each observable, it restore parameters
/// to the fit values right after pdfFit.
/// It then loops through all observables to get plot for each of them.
/// For each plot, it puts data points, pdf contour,
/// and parameters and chi square (if their booleans set) on.
/// Before it returns the last plot, it restores parameters back
/// to values before plotting.
RooPlot *rarBasePdf::doPdfPlot(TList &plotList, TString pdfList)
{
  RooPlot *frame(0);
  // first for xtraPdfs
  frame=doXPdfPlot(plotList, pdfList);
  
  cout<<endl<<" In rarBasePdf doPdfPlot for "<<GetName()<<endl;
  
  //if (!getControlBit("PdfFit")) return frame;
  TString pdfFitStr=readConfStr("pdfFit","yes",getVarSec());
  if (!getControlBit("PdfPlot")) return frame;
  if (getControlBit("PdfPlotDone")) return frame;
  // are we in the list for plotting
  rarStrParser pdfListParser=pdfList;
  if (pdfListParser.nArgs()>0)
    if (!pdfListParser.Have(GetName())) return frame;
  
  // check if it has its pdf and data
  if ((!_thePdf)||(!_theData)) return frame;
  
  // save obs
  string obsSaveStr;
  writeToStr(_fObsSet, obsSaveStr);
  // restore the fit value for plotting
  string coeffParamSSaver;
  RooArgSet *theParams=_thePdf->getParameters(_theData);
  writeToStr(*theParams, coeffParamSSaver);
  readFromStr(*theParams, _afterFitSaverStr);
  
  // plot for each obs
  RooArgList fObsList(_fObsSet);
  for (Int_t i=0; i<fObsList.getSize(); i++) {
    // get obs
    RooRealVar *theObs=dynamic_cast<RooRealVar*>(fObsList.at(i));
    if (!theObs) continue;
    TString nBinStr="0";
    // if obs is in _protDataVars or _condObss,
    // do not plot unless explicitly requested
    if (_protDataVars.find(theObs->GetName())) nBinStr="-1";
    if (_condObss.find(theObs->GetName())) nBinStr="-1";
    // get nBins
    Int_t nBins=atoi(readConfStr(Form("plotBins_%s", theObs->GetName()),
				 nBinStr, getVarSec()));
    if (nBins<0) continue;
    if (0==nBins) nBins=theObs->getBins();
    // plot range
    Double_t plotMin=theObs->getMin();
    Double_t plotMax=theObs->getMax();
    getRange(theObs, "plotRange_", plotMin, plotMax);
    // normalization CmdArg
    // see post on HN from Wouter for another way
    // to handle fit and plot for sub-ranges
    // using named ranges, like:
    // mes.setRange("signal", 5.27, 5.29);
    // de.setRange("signal",-0.1,0.1) ;
    // data->plotOn(frame2,CutRange("signal")) ;
    // model.plotOn(frame2,ProjectionRange("signal")) ;
    RooCmdArg plotCutNormCmd;
    TString rangeCuts="1";
    if (plotMin>theObs->getMin())
      rangeCuts+=Form("&&(%s>%f)",theObs->GetName(), plotMin);
    if (plotMax<theObs->getMax())
      rangeCuts+=Form("&&(%s<%f)",theObs->GetName(), plotMax);
    if ("1"!=rangeCuts) {
      Double_t d=_theData->sumEntries();
      Double_t n=_theData->sumEntries(rangeCuts);
      Double_t plotRangeNorm=n/d;
      cout<<" Based on total evts "<<d<<" and fitted evts "<<n<<endl
          <<" norm. scale factor for plotting: "<<plotRangeNorm<<endl
          <<endl;
      if (plotRangeNorm!=1) {
        cout<<" Using relative norm. factor "<<plotRangeNorm<<endl;
        plotCutNormCmd=Normalization(plotRangeNorm, RooAbsReal::Relative);
      }
    }
    theObs->setRange(plotMin, plotMax);
    // do we have projWData for plotting
    RooDataSet *projWData(0);
    TString projWDataStr=
      readConfStr(Form("projWData_%s", theObs->GetName()), "no", getVarSec());
    if ("yes"==projWDataStr) projWData=_theData;
    else if ("no"!=projWDataStr)
      projWData=getDatasets()->getData(projWDataStr);
    if (projWData) {
      cout<<" Get plot of "<<theObs->GetName()<<" for Pdf "<<_thePdf->GetName()
	  <<" using projWData "<<projWData->GetName()<<endl;
    }
    // any plots for simPdf?
    if (_thisSimPdfWOP&&pdfFitStr.Contains("simFit")) {
      // all sim comps
      RooArgSet *simComps=(RooArgSet*)_thisSimPdfWOP->getComponents()->
        selectByName(Form("%s_*",_thePdf->GetName()));
      RooAbsCategoryLValue *simCats=(RooAbsCategoryLValue *)
        (&((RooSimultaneous*)_thisSimPdfWOP)->indexCat());
      RooArgSet inputCats;
      RooSuperCategory *simSuperCats=dynamic_cast<RooSuperCategory*>(simCats);
      if (simSuperCats) inputCats.add(simSuperCats->inputCatList());
      else inputCats.add(*simCats);
      TIterator *cIter=inputCats.createIterator();
      RooCategory *theCat(0);
      while(theCat=(RooCategory*)cIter->Next()) {
	TIterator *tIter=theCat->typeIterator();
	RooCatType *theType(0);
	while (theType=(RooCatType*)tIter->Next()) {
	  // do we have simPdf comp for this cat?
	  TString compStr="";
	  RooArgSet compParams;
	  TIterator *compIter=simComps->createIterator();
	  RooAbsPdf *thePdf(0);
	  while(thePdf=(RooAbsPdf*)compIter->Next()) {
	    TString thePdfName=thePdf->GetName();
	    thePdfName.Remove(0, strlen(_thePdf->GetName()));
	    if (!thePdfName.Contains(theType->GetName())) continue;
	    compStr+=thePdf->GetName();
	    compStr+=",";
	    compParams.add(*thePdf->getParameters(_theData));
	  }
	  delete compIter;
	  if (""==compStr) continue;
	  // do we have dataset for this cat?
	  TString catCut=Form("%s==%s::%s",theCat->GetName(),
			      theCat->GetName(),theType->GetName());
	  RooDataSet *subData=(RooDataSet*)_theData->reduce(catCut);
	  if (subData->numEntries()<=0) {
	    delete subData;
	    continue;
	  }
	  // now plot the subData and pdf
	  RooPlot *subFrame=theObs->frame(plotMin, plotMax, nBins);
	  if (getControlBit("ParamsOnPlot"))
	    doParamsOnPlot(subFrame, &compParams);
	  RooLinkedList plotOpts;
	  // X Error size relative to bin width, default 1.0
	  Double_t xerrorscale = atof(readConfStr("XErrorSize", "1.0", 
						  getVarSec()));
	  RooCmdArg xerrArg = XErrorSize(xerrorscale);
	  plotOpts.Add((TObject*)&xerrArg);
	  RooCmdArg datErrArg = DataError(RooAbsData::SumW2);
	  if (subData->isWeighted()) // use sumw2
	    plotOpts.Add((TObject*)&datErrArg);
	  subData->plotOn(subFrame, plotOpts);
	  _thisSimPdfWOP->plotOn(subFrame, plotCutNormCmd, Components(compStr),
				 ProjWData(*subData));
	  if (getControlBit("Chi2OnPlot")) doChi2OnPlot(subFrame);
	  subFrame->SetNameTitle(Form("%s_%s_%s",theObs->GetName(),
				      _thePdf->GetName(), catCut.Data()),
				 Form("%s %s %s",theObs->GetTitle(),
				      _thePdf->GetTitle(), catCut.Data()));
	  //#ifndef RAR_USE_ROOT5
	  subFrame->SetTitleSize(0.05);
	  //#endif
	  plotList.Add(subFrame);
	  
	  delete subData;
	}
	delete tIter;
      }
      delete cIter;
      delete simComps;
      // now plot for full simPdf
      {
	RooPlot *subFrame=theObs->frame(plotMin, plotMax, nBins);
	if (getControlBit("ParamsOnPlot")) doParamsOnPlot(subFrame);
	RooLinkedList plotOpts;
	Double_t xerrorscale = atof(readConfStr("XErrorSize", "1.0", 
						getVarSec()));
	RooCmdArg xerrArg = XErrorSize(xerrorscale);
	plotOpts.Add((TObject*)&xerrArg);
	RooCmdArg datErrArg = DataError(RooAbsData::SumW2);
	if (_theData->isWeighted()) // use sumw2
	  plotOpts.Add((TObject*)&datErrArg);
	_theData->plotOn(subFrame, plotOpts);
	_thisSimPdfWOP->plotOn(subFrame, ProjWData(*_theData),plotCutNormCmd);
	if (getControlBit("Chi2OnPlot")) doChi2OnPlot(subFrame);
	subFrame->SetNameTitle(Form("%s_%s",theObs->GetName(),
				    _thisSimPdfWOP->GetName()),
			       Form("%s %s",theObs->GetTitle(),
				    _thisSimPdfWOP->GetTitle()));
	//#ifndef RAR_USE_ROOT5
	subFrame->SetTitleSize(0.05);
	//#endif
	plotList.Add(subFrame);
      }
    }
    if ("simFitOnly"==pdfFitStr) continue;
    // now categorized plot for the total pdf
    TString plotWCatStr=
      readConfStr(Form("plotWCat_%s",theObs->GetName()),"notSet",getVarSec());
    if ("notSet"==plotWCatStr)
      plotWCatStr=readConfStr("plotWCat","no",getVarSec());
    if ("no"!=plotWCatStr) {
      rarStrParser plotCatParser=plotWCatStr;
      while (plotCatParser.nArgs()>0) {
        TString catName=plotCatParser[0];
        plotCatParser.Remove();
        RooAbsCategory *plotCat=(RooAbsCategory *)_fullObs->find(catName);
        if (!plotCat) {
          cout<<" W A R N I N G !"<<endl
              <<" category "<<catName<<" does not exist!"<<endl;
          continue;
        }
        TIterator *cIter=plotCat->typeIterator();
        RooCatType *cType(0);
        while(cType=(RooCatType*)cIter->Next()) {
          TString catCut=catName+"=="+catName+"::"+cType->GetName();
          RooDataSet *catSliceData=(RooDataSet*)_theData->reduce(catCut);
          RooPlot *subFrame=theObs->frame(plotMin, plotMax, nBins);
          if (catSliceData->isWeighted()) // use sumw2
            catSliceData->plotOn(subFrame, DataError(RooAbsData::SumW2));
          else
            catSliceData->plotOn(subFrame);
          _thePdf->plotOn(subFrame, ProjWData(*catSliceData), plotCutNormCmd);
          subFrame->SetNameTitle(Form("%s_%s_%s",theObs->GetName(),
                                      _thePdf->GetName(), catCut.Data()),
                                 Form("%s %s %s",theObs->GetTitle(),
                                      _thePdf->GetTitle(), catCut.Data()));
          plotList.Add(subFrame);
          delete catSliceData;
        }
        delete cIter;
      }
    }
    // pdfPlot for total pdf
    // get plot frame
    frame=theObs->frame(plotMin, plotMax, nBins);
    // get params
    if (getControlBit("ParamsOnPlot")) doParamsOnPlot(frame);
    // data points
    RooLinkedList plotOpts;
    Double_t xerrorscale = atof(readConfStr("XErrorSize", "1.0", getVarSec()));
    RooCmdArg xerrArg = XErrorSize(xerrorscale);
    plotOpts.Add((TObject*)&xerrArg);
    RooCmdArg datErrArg = DataError(RooAbsData::SumW2);
    if (_theData->isWeighted()) // use sumw2
      plotOpts.Add((TObject*)&datErrArg);
    _theData->plotOn(frame, plotOpts);
    // pdf
    _thePdf->plotOn(frame, ProjWData(*projWData), plotCutNormCmd);
    // put chi2 on
    if (getControlBit("Chi2OnPlot")) doChi2OnPlot(frame);
    // any sub plots?
    if (getControlBit("CompsOnPlot")) {
      Int_t lineStyle=4;
      Int_t nComp=_subPdfs.getSize();
      for (Int_t j=0; j<nComp; j++) {
	if (lineStyle>2) lineStyle=2;
	else lineStyle=4;
	RooAbsPdf *theSubPdf=(RooAbsPdf*)_subPdfs.at(j);
	_thePdf->plotOn(frame, ProjWData(*projWData),LineWidth(2),
			plotCutNormCmd, Components(theSubPdf->GetName()),
			Name(theSubPdf->GetName()),
			LineStyle(lineStyle),LineColor(getColor(j)));
      }
    }
    // any sub data points
    if (getControlBit("CompsOnPlot")&&getControlBit("CompsDataOnPlot")) {
      RooAddPdf *refPdf(0);
      // do we have ref pdf in the config?
      TString refName=readConfStr("compsDataOnPlot", "no", getVarSec());
      if ("yes"==refName) refPdf=(RooAddPdf*)_thePdf;
      else {
	// ref to RooRarFit PDF?
	rarBasePdf *rarRef=(rarBasePdf*)(_rarPdfs.FindObject(refName));
	if (rarRef) refPdf=(RooAddPdf*)rarRef->getPdf();
	else {
	  cout<<"Can not find ref pdf named "<<refName<<endl
	      <<"Use the default one: "<<GetName()<<endl;
	  refPdf=(RooAddPdf*)_thePdf;
	}
      }
      // make sure the refPdf is RooAddPdf
      if (refPdf->ClassName()!=TString("RooAddPdf")) {
        cout <<" The reference pdf "<<refPdf->GetName()
             <<" for data point plotting is not RooAddPdf"<<endl;
        exit(-1);
      }
      Int_t nComp=refPdf->pdfList().getSize();
      // now create nComp histograms for data points
      TList hList;
      for (Int_t j=0; j<nComp; j++) {
        TH1F *h=new
          TH1F(Form("subDP_%s_%s", theObs->GetName(),
                    refPdf->pdfList()[j].GetName()), "sub data points",
               frame->GetXaxis()->GetNbins(),
               frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax());
        h->Sumw2();
        h->SetMarkerStyle(8+j);
        hList.Add(h);
      }
      // 
      refPdf->attachDataSet(*_theData);
      // now get data point plots
      Int_t nEvts=_theData->numEntries();
      RooArgSet *normSet=(RooArgSet *)_theData->get(0);
      for (Int_t evtIdx=0; evtIdx<nEvts; evtIdx++) {
        _theData->get(evtIdx);
        Double_t val=((RooAbsReal*)normSet->find(theObs->GetName()))->getVal();
        Double_t totProb=refPdf->getVal(normSet);
        Double_t lastW=1.;
        for (Int_t j=0; j<nComp; j++) {
          TH1F *h=(TH1F *)hList.At(j);
          if (j<nComp-1) {
            RooAbsPdf *subPdf=(RooAbsPdf*)refPdf->pdfList().at(j);
            RooAbsReal *subCoef=(RooAbsReal*)refPdf->coefList().at(j);
            Double_t w=subPdf->getVal(normSet)*subCoef->getVal()/totProb;
            lastW-=w;
            h->Fill(val, w);
          } else {
            h->Fill(val, lastW);
          }
        }
      }
      // add the histograms into sub frames
      for (Int_t j=0; j<nComp; j++) {
        RooPlot *subframe=theObs->frame(nBins);
        subframe->SetName(Form("sub_%s_%s", theObs->GetName(),
                               refPdf->pdfList()[j].GetName()));
        RooAbsPdf *theSubPdf=(RooAbsPdf*)_subPdfs.at(j);
        TH1F *h=(TH1F *)hList.At(j);
        // convert TH1 to RooHist
        RooHist *rHist=new RooHist(*h);
        subframe->addPlotable(rHist, "P");
	theSubPdf->plotOn(subframe, ProjWData(*projWData));
        plotList.Add(subframe);
      }
    }
    // set frame's attrs
    frame->SetNameTitle(Form("%s_%s", theObs->GetName(), _thePdf->GetName()),
			Form("%s %s",theObs->GetTitle(),_thePdf->GetTitle()));
    //#ifndef RAR_USE_ROOT5
    frame->SetTitleSize(0.05);
    //#endif
    plotList.Add(frame);
  }
  // restore params
  readFromStr(*theParams, coeffParamSSaver);
  // restore obs saved
  readFromStr(_fObsSet, obsSaveStr);
  setControlBits("PdfPlotDone");
  return frame;
}

/// \brief Put parameters on the plot
///
/// \param frame The plot frame
/// \param params Params to draw
/// \param sigDigits Number of significant digits
/// \param options Plot options
/// \param xmin xmin of the param frame
/// \param xmax xmax of the param frame
/// \param ymax ymax of the param frame
/// \return The plot frame
///
/// It is modified from
/// <a href="http://roofit.sourceforge.net/docs/classref/RooAbsPdf.html#RooAbsPdf:paramOn"
/// target=_blank>RooAbsPdf::paramOn</a> to fix some bug and
/// also allow chi square of the plot on the param frame.
/// The chi square here has not be calculated yet,
/// so the text line for it is just a place holder.
/// #doChi2OnPlot is the function to fill that field.
RooPlot *rarBasePdf::doParamsOnPlot(RooPlot* frame, RooArgSet *params,
				    Int_t sigDigits, Option_t *options,
				    Double_t xmin,Double_t xmax, Double_t ymax)
{
  // parse the options
  TString opts = options;
  opts.ToLower();
  Bool_t showConstants= opts.Contains("c");
  
  // calculate the box's size, adjusting for constant parameters
  if (!params) params = _thePdf->getParameters(_theData);
  TIterator* pIter = params->createIterator();
  
  Int_t nLine(0);
  Int_t nFreeParam(0);
  Double_t ymin(ymax), dy(0.06);
  RooRealVar *var = 0;
  while(var=(RooRealVar*)pIter->Next()) {
    if (var->ClassName()!=TString("RooRealVar")) continue;
    if (!var->isConstant()) nFreeParam++;
    if(showConstants || !var->isConstant()) {
      ymin-= dy;
      nLine++;
    }
  }
  if(getControlBit("Chi2OnPlot")) {
    ymin-= dy;
    nLine++;
  }
  if (nLine<=0) nLine=1;
  Double_t step=1./nLine;
  
  // create the box and set its options
  TPaveText *box= new TPaveText(xmin,ymin,xmax,ymax,"BRNDC");
  if(!box) return frame;
  box->SetFillColor(0);
  box->SetBorderSize(1);
  box->SetTextAlign(12);
  box->SetTextSize(0.04F);
  box->SetFillStyle(1001);
  box->SetFillColor(0);
  TText *text(0); 
  pIter->Reset() ;
  Double_t y=1-.4*step;
  // add chi2 if specified
  if(getControlBit("Chi2OnPlot")) {
    text= box->AddText(0, y, Form("%d", nFreeParam));
    y-=step;
  }
  while(var=(RooRealVar*)pIter->Next()) {
    if (var->ClassName()!=TString("RooRealVar")) continue;
    if(var->isConstant() && !showConstants) continue;
    if(var->GetName()==TString(var->getPlotLabel())) {
      TString label=var->GetTitle();
      label.ReplaceAll(" ({", " _{");
      label.ReplaceAll(" (", " _{");
      label.ReplaceAll("})", "}");
      label.ReplaceAll(")", "}");
      var->setPlotLabel(label);
    }
    TString *formatted= var->format(sigDigits, opts.Data());
    text= box->AddText(0, y, formatted->Data());
    y-=step;
    delete formatted;
  }
  
  // Add box to frame 
  frame->addObject(box) ;
  
  delete pIter ;
  return frame ;
}

/// \brief Put chi square on the plot
///
/// \param frame The plot frame
/// \return The plot frame
///
/// It replaces the place holder with
/// the actual chi square value on the plot.
/// The chisquare number has to be put on in this two step manner,
/// because param frame has to be the first on in the frame
/// so it will not block any other plots,
/// but chisquare has to be calculated after all the plots are drawn,
/// which means the value before the plots is just place holder.
RooPlot *rarBasePdf::doChi2OnPlot(RooPlot* frame)
{
  if (!getControlBit("ParamsOnPlot")) return frame;
  TPaveText* pbox = (TPaveText*) frame->findObject("TPave");
  if (!pbox) return frame;
  if (!pbox->GetSize()) return frame;
  TLatex *text=(TLatex*)pbox->GetLine(0);
  if ("nbin"==readConfStr("chi2OnPlot", "dof", getVarSec())) {
    text->SetTitle(Form("%s%.3f", "#chi^{2}/n = ", frame->chiSquare()));
  } else { // DOF
    Int_t nFreeParam=atoi(text->GetTitle());
    //text->SetLimitIndiceSize(1);
    //text->SetTitle(Form("%s%.3f", "#chi^{2}/_{^{DOF}} = ",
    text->SetTitle(Form("%s%.3f", "#chi^{2} / ndf = ",
			frame->chiSquare(nFreeParam)));
  }
  
  return frame;
}

/// \brief Return prototype var generator for toy study
/// \return The generator
///
/// It constrcuts prototype var generator for toy study
RooAbsPdf *rarBasePdf::getProtGen()
{
  if (_theProtGen) return _theProtGen;
  _theProtGen=_myDummyPdf;
  if (_protGenPdfs.getSize()>0) {
    RooArgList pdfList(_protGenPdfs);
    pdfList.add(*_theProtGen);
    _theProtGen=new RooProdPdf
      (Form("protGen_%s",GetName()),Form("Gen for %s",GetName()), pdfList);
  }
  
  return _theProtGen;
}

/// \brief Get a PDF not depending on specified observable for sPlot
/// \param ignoredObs The observables to check
/// \return The Pdf required
///
/// It checks if the pdf depends on the observables.
/// if yes, it returns dummy constant pdf;
/// if no, it returns the default pdf
RooAbsPdf *rarBasePdf::getPdfWOvar(RooArgList ignoredObs)
{
  RooAbsPdf *thePdf=_thePdf;
  // check if thePdf depends on theVar
  if (!(thePdf->dependsOn(ignoredObs))) return thePdf;
  // create a dummy for it
  TString theName=Form("the_%s_%s", ignoredObs[0].GetName(), GetName());
  TString theTitle=Form("%s w/o %s", GetTitle(), ignoredObs[0].GetName());
  thePdf=new RooGenericPdf(theName, theTitle, "1",
			   *_thePdf->getParameters(_fullObs));
  //thePdf->Print("v");
  //thePdf->Print();
  return thePdf;
}

/// \brief Return the pdf if it depends on the var
/// \param theVar The var to check
/// \return This (sim)Pdf if it depends on the var
///
/// It checks if the pdf depends on the obs;
/// if yes, it returns this (sim)Pdf
RooAbsPdf *rarBasePdf::getDPdfWvar(RooRealVar *theVar)
{
  RooAbsPdf *retVal(0);
  if (!_thePdf->dependsOn(*theVar)) return retVal;
  retVal=_thePdf;
  if (_thisSimPdf) retVal=_thisSimPdf;
  return retVal;
}

/// \brief Set PDF control bit in the control string by reading config
///
/// \param controlBitStr String of control bit
/// \param bitConfigStr Config name for the control bit
/// \param configSec Config section
///
/// It saves the control bit specified by \p controlBitStr
/// to #_controlStr.
void rarBasePdf::setControlBit(TString controlBitStr, TString bitConfigStr,
			       TString configSec)
{
  // bit
  TString bit="yes";
  // first get bit name
  if (controlBitStr.BeginsWith("no")) {
    controlBitStr.Replace(0, 2, "");
    bit="no";
  }
  // then construct config name
  if (""==bitConfigStr) bitConfigStr=controlBitStr;
  if (""==configSec) configSec=getVarSec();
  bit=readConfStr(bitConfigStr, bit, configSec);
  if ("no"==bit) setControlBits("no"+controlBitStr);
  else setControlBits(controlBitStr);
}

/// \brief Set PDF control bit in the control string
///
/// \param controlBitsStr String of control bits
///
/// It saves all the control bits in \p controlBitsStr
/// to #_controlStr.
/// If the control name is \p PdfPlot,
/// valid bits in the control string can be
/// \p PdfPlot, which means the bit is set to \p true,
/// or \p noPdfPlot, which means the bit is set to \p false.
void rarBasePdf::setControlBits(TString controlBitsStr)
{
  rarStrParser controlBitsStrParser=controlBitsStr;
  for (Int_t i=0; i<controlBitsStrParser.nArgs(); i++) {
    TString myBit=controlBitsStrParser[i];
    // remove any such bit
    rarStrParser masterStrParser=_controlStr;
    _controlStr="";
    for (Int_t j=0; j<masterStrParser.nArgs(); j++) {
      if (myBit==masterStrParser[j]) continue;
      if ("no"+myBit==masterStrParser[j]) continue;
      if (myBit=="no"+masterStrParser[j]) continue;
      _controlStr+=" "+masterStrParser[j];
    }
    // add the bit
    _controlStr+=" "+myBit;
  }
}

/// \brief Return a specified PDF control bit in the control string
/// \param controlBitStr String of control bit
/// \return The value of the control bit
///
/// It gets the value of the specified control bit in #_controlStr.
/// It returns \p true if controlBitStr matches a bit inside #_controlStr,
/// otherwise \p false.
/// If the bit does not exist, it reutrns \p false.
Bool_t rarBasePdf::getControlBit(TString controlBitStr)
{
  rarStrParser controlStrParser=_controlStr;
  for (Int_t i=0; i<controlStrParser.nArgs(); i++) {
    if (controlBitStr==controlStrParser[i]) return kTRUE;
    if ("no"+controlBitStr==controlStrParser[i]) return kFALSE;
    if (controlBitStr=="no"+controlStrParser[i]) return kFALSE;
  }
  // not exist
  if (controlBitStr.BeginsWith("no")) return kTRUE;
  return kFALSE;
}

/// \brief Return color from repository
/// \param i Index for color
/// \return Color returned
Int_t rarBasePdf::getColor(Int_t i)
{
  return _rarColors[i%NCOLORS];
}
