/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarMLPdf.cc,v 1.28 2014/09/14 17:33:53 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides ML model class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides ML model class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "Roo1DTable.h"
#include "RooAbsPdf.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "RooAddPdf.h"

#include "rarMLFitter.hh"
#include "rarMLPdf.hh"

ClassImp(rarMLPdf)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarMLPdf::rarMLPdf()
  : rarAdd(),
    _specialStr("")
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
rarMLPdf::rarMLPdf(const char *configFile, const char *configSec,
		   const char *configStr,
		   rarDatasets *theDatasets, RooDataSet *theData,
		   const char *name, const char *title)
  : rarAdd(configFile, configSec, configStr,
	      theDatasets, theData, name, title, kFALSE, kFALSE),
    _specialStr("")
{
  init();
}

rarMLPdf::~rarMLPdf()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// The coeff info is read in with rarAdd::init(),
/// then it re-define the coeffs using the splitting info
/// of the final ML model,
/// and finally it builds the RooAddPdf.
///
/// Essentially, the splitting can be done manually
/// by fitting on to separate sub-datasets for individual categories.
/// So for each category, the yields, nSig, nBkg, etc.
/// need to be split, which means for n number of category types,
/// there should be n sub-yields.
/// To get the final yield, just sum up all these values.
/// It is then convenient to have the total yield,
/// and n-1 sub-yield fractions as free parameters:
/// \verbatim
/// subY1=Ntot*(1-frac2...-fracn),
/// subY2=Ntot*frac2,
/// ...
/// subYn=Ntot*fracn.\endverbatim
/// That is the basic idea for #_coeffs splitting of the extended AddPdf.
/// For each category, you will need similar splitting setttings.
/// Depending on the auto splitting method used,
/// the unity constraint can be imposed by \p RooFit or \p RooRarFit itself.
void rarMLPdf::init()
{
  // make sure it is extended
  if (_nCoeff!=_nComp) {
    cout<<"Pdf model should be extended pdf"<<endl;
    exit(-1);
  }
  // add compCat type, the last type for embedded events
  for (Int_t i=0; i<=_nComp; i++) _compCat.defineType(Form("%02d",i));
  // save the original coeff List anyway
  _sCoeffs.add(_coeffs);
  // to deal with simu fit and splitting of coeffs
  TString masterSec=getMasterSec();
  Bool_t simFit=getControlBit("SimFit")&&
    ("simple"!=readConfStr("simultaneousFit", "yes",masterSec));
  // for split specials
  if (simFit) {
    createAbsVars("splitSpecials", &_specialSet);
    // run master sec's
    TString myVarSec=getVarSec();
    setVarSec(masterSec);
    createAbsVars("splitSpecials", &_specialSet);
    setVarSec(myVarSec);
  }
  TString yieldSplitMethod=readConfStr("yieldSplitMethod", "auto", masterSec);
  if (simFit && ("manual"!=yieldSplitMethod)) {
    // use master sec config
    TString myVarSec=getVarSec();
    setVarSec(masterSec);
    // get coeff list, coeff param set, param list
    RooArgSet coeffParamSet("Coeff Parameter Set");
    for (Int_t i=0; i<_nCoeff; i++) {
      coeffParamSet.add(*(_coeffs[i].getParameters(*_fullObs)));
    }
    RooArgList coeffParamList(coeffParamSet);
    // let's make sure coeffs are not split in the mlFitter
    // coz we are going to split them here
    TString pdfSplitStr=readConfStr(Form("the_%s", GetName()), "", masterSec);
    if ("semi"==yieldSplitMethod) {
      for (Int_t i=0; i<_coeffs.getSize(); i++) {
	if (pdfSplitStr.Contains(_coeffs[i].GetName())) {
	  cout<<"Please do not include "<<_coeffs[i].GetName()
	      <<" or its dependents in your simPdfBuilder configs"<<endl
	      <<"in section \""<<masterSec<<"\":"<<endl
	      <<pdfSplitStr<<endl
	      <<"We will split it for you."<<endl;
	  exit(-1);
	}
      }
      for (Int_t i=0; i<coeffParamList.getSize(); i++) {
	if (pdfSplitStr.Contains(coeffParamList[i].GetName())) {
	  cout<<"Please do not include "<<coeffParamList[i].GetName()
	      <<" or its dependents in your simPdfBuilder configs"<<endl
	      <<"in section \""<<masterSec<<"\":"<<endl
	      <<pdfSplitStr<<endl
	      <<"We will split it for you."<<endl;
	  exit(-1);
	}
      }
    }
    // get splitting cats
    TString physCat=getFitter()->getPhysCat();
    TString splitCats=getFitter()->getSplitCats();
    // build split coeffs
    _coeffs.removeAll();
    for (Int_t i=0; i<_nCoeff; i++) {
      // formula for the final coeff
      TString coeffArgsStr=_sCoeffs[i].GetName();
      TString formulaStr="@0";
      // get frac split rules
      TString fracRuleStr=readConfStr("fracRule", splitCats, masterSec);
      fracRuleStr=readConfStr("fracRule_"+coeffArgsStr,fracRuleStr, masterSec);
      rarStrParser catsStrParser=fracRuleStr;
      Int_t nSCat=catsStrParser.nArgs();
      // create multiple cat coeffs
      RooArgSet splitCatSet(getFitter()->getSplitCatSet());
      for (Int_t j=0; j<nSCat; j++) {
	TString catName=catsStrParser[j];
	RooAbsCategory *theCat=getFitter()->getSplitCat(splitCatSet, catName);
	if (!theCat) {
	  cout<<"Can not find cat "<<catName<<endl
	      <<"Please make sure the cat you are using exists"<<endl;
	  exit(-1);
	}
	// rename the catName
	catName=theCat->GetName();
	// cat ceoff name
	TString catCoeffName=Form("Frac_%s_%s",
				  _sCoeffs[i].GetName(), catName.Data());
	TString catCoeffStr=readConfStr(catCoeffName, "notSet", masterSec);
	if ("notSet"==catCoeffStr)
	  catCoeffStr=Form("%s %s RooRealVar \"cat coeff\" C 1. 0. 1.",
			   catCoeffName.Data(), catCoeffName.Data());
	else
	  catCoeffStr=catCoeffName+" "+catCoeffStr;
	RooRealVar *catCoeff=(RooRealVar*)createAbsVar(catCoeffStr);
	catCoeffName=catCoeff->GetName();
	// save the var name
	saveFracName(catCoeffName);
	// addToConfStr("Ignored", catCoeff->GetName(), myVarSec);
	coeffArgsStr+=Form(" %s ", catCoeff->GetName());
	formulaStr+=Form("*@%d",j+1);
	// let's make sure catCoeff is created as RooRealVar
	if ("RooRealVar"!=TString(catCoeff->ClassName())) {
	  cout<<catCoeffName<<" is created in config file as "
	      <<catCoeff->ClassName()<<","<<endl
	      <<"but it should be RooRealVar"<<endl;
	  exit(-1);
	}
	// now, we create split catCoeffs for this catCoeff
	// first we find out how many types for this cat
	TString catTypesStr="";
	TIterator* catTypeIter = theCat->typeIterator();
	RooCatType *theType(0);
	while(theType=(RooCatType*)catTypeIter->Next()) {
	  catTypesStr+=Form(" %s ", theType->GetName());
	}
	delete catTypeIter;
	//cout<<catTypesStr<<endl;
	rarStrParser catTypesStrParser=catTypesStr;
	if (catTypesStrParser.nArgs()<1) continue; // nothing in the cat
	TString cat1FormulaStr="1.";
	TString cat1ArgSet="";
	Roo1DTable* theTable(0);
	if (_theData) theTable = _theData->table(*theCat);
        // check if we have fracSrc config
        Bool_t useFracSrc=kFALSE;
        TString fracSrcStr=
          readConfStrCnA(Form("fracSrc_%s_%s", _sCoeffs[i].GetName(),
                              catName.Data()),"notSet");
        if ("notSet"==fracSrcStr)
          fracSrcStr=
            readConfStrCnA(Form("fracSrc_%s",_sCoeffs[i].GetName()), "notSet");
        if ("notSet"==fracSrcStr)
          fracSrcStr=readConfStrCnA("fracSrc", "notSet");
        if ("notSet"!=fracSrcStr) {
          RooDataSet *theData=getDatasets()->getData(fracSrcStr);
          if (theData) {
            theTable = theData->table(*theCat);
            useFracSrc=kTRUE;
            cout<<" Using table from "<<theData->GetName()
                <<" for "<<_sCoeffs[i].GetName()<<" "<<catName<<endl;
            theTable->Print();
          }
        }
	// now we need to find out which type to be the non-free
	Int_t nonIdx=0;
	// first for auto mode
	if ("auto"==yieldSplitMethod) {
	  Bool_t notFound=(catTypesStrParser.nArgs()>1); // if #types>1
	  for (Int_t k=0; k<catTypesStrParser.nArgs(); k++) {
	    TString myNonFreeName=catCoeffName+"["+catTypesStrParser[k]+"]";
	    Int_t strIdx=pdfSplitStr.Index(myNonFreeName);
	    if (strIdx<0) continue;
	    // check cat seq
	    TString subpdfSplitStr=pdfSplitStr(0,strIdx);
	    Int_t colonIdx=subpdfSplitStr.Last(':');
	    if (colonIdx<=0) break;
	    rarStrParser pdfSplitStrParser=TString(subpdfSplitStr(0,colonIdx));
	    //cout<<pdfSplitStr<<endl<<myNonFreeName<<endl<<strIdx<<endl;
	    // get cats
	    TString splitCatName=
	      pdfSplitStrParser[pdfSplitStrParser.nArgs()-1];
	    splitCatName.ReplaceAll(",", "_");
	    //cout<<splitCatName<<endl;
	    if (splitCatName==catName) {
	      nonIdx=k;
	      notFound=kFALSE;
	    }
	    break;
	  }
	  if (notFound) {
	    cout<<" Can not find splitting specif. for"<<endl
		<<" "<<catCoeffName<<" wrt types:"<<endl
		<<" "<<catTypesStr<<endl
		<<" in the_"<<GetName()<<endl<<pdfSplitStr<<endl
		<<" from section ["<<getVarSec()<<"]"<<endl;
	    exit(-1);
	  }
	} else { // semi mode
	  TString nonFreeStr=readConfStr("nonFreeCatTypes", "", getVarSec());
	  for (Int_t k=0; k<catTypesStrParser.nArgs(); k++) {
	    if (nonFreeStr.
		Contains(catCoeffName+"["+catTypesStrParser[k]+"]")) {
	      nonIdx=k;
	      break;
	    }
	  }
	}
	for (Int_t k=0; k<catTypesStrParser.nArgs(); k++) {
	  // make sure the name is a cat type
	  const RooCatType *theType=theCat->lookupType(catTypesStrParser[k]);
	  if (!theType) {
	    cout <<"The cat type "<<catTypesStrParser[k]
		 <<" does not exist in cat "<<theCat->GetName()<<endl;
	    exit(-1);
	  }
	  if (k==nonIdx) continue;
	  TString catCoeffTypeName=catCoeffName+"_"+catTypesStrParser[k]+"";
	  Double_t typeFrac=1./catTypesStrParser.nArgs();
	  if (theTable) typeFrac=theTable->getFrac(catTypesStrParser[k]);
	  RooAbsReal *catCoeffType=(RooAbsReal*)
	    createAbsVar(Form("%s %s RooRealVar \"cat type coeff\" C %f 0 1",
			      catCoeffTypeName.Data(), catCoeffTypeName.Data(),
			      typeFrac));
	  // save the var name
	  saveFracName(catCoeffTypeName);
	  if (useFracSrc)
            addToConfStr("Ignored", catCoeffType->GetName(), myVarSec);
	  _specialSet.add(*catCoeffType);
	  if (2==catTypesStrParser.nArgs()) _asymSet.add(*catCoeffType);
	  Int_t fvIdx = k<nonIdx ? k : k-1;
	  cat1FormulaStr+=Form("-@%d", fvIdx);
	  cat1ArgSet+=Form(" %s ", catCoeffType->GetName());
	}
	// for nonfree cat type
	TString catCoeffTypeName=catCoeffName+"_"+catTypesStrParser[nonIdx]+"";
	// old way, ie RoorarFit way, RooRealVar
	RooAbsReal *catCoeffType=(RooAbsReal*)
	  createAbsVar(Form("%s %s RooFormulaVar %s %s",
			    catCoeffTypeName.Data(), catCoeffTypeName.Data(),
			    cat1FormulaStr.Data(), cat1ArgSet.Data()));
	// make sure it is RooFormulaVar
	if ("RooFormulaVar"!=TString(catCoeffType->ClassName())) {
	  cout<<catCoeffTypeName<<" is created in config file as "
	      <<catCoeffType->ClassName()<<","<<endl
	      <<"but it should be RooFormulaVar"<<endl;
	  exit(-1);
	}
	_specialSet.add(*catCoeffType);
	if (catTypesStrParser.nArgs()>1) _cat1Set.add(*catCoeffType);
	if ("semi"==yieldSplitMethod) {
	  _specialStr+=Form(" %s : %s ", catName.Data(), catCoeffName.Data());
	}
      }
      // do we still have any left over?
      if (splitCatSet.getSize()>0) {
	cout<<" No splitting rules for "<<_sCoeffs[i].GetName()<<" wrt:"<<endl;
	splitCatSet.Print();
	cout<<" Please check config fracRule or fracRule_"
	    <<_sCoeffs[i].GetName()<<" in section \""<<masterSec<<"\":"<<endl;
	exit(-1);
      }
      TString coeffName=Form("theSim_%s", _sCoeffs[i].GetName());
      //cout<<coeffName<<endl
      //<<formulaStr<<endl
      //<<coeffArgsStr<<endl;
      RooAbsReal *coeff=(RooAbsReal*)
	createAbsVar(Form("%s %s RooFormulaVar %s %s",
			  coeffName.Data(), coeffName.Data(),
			  formulaStr.Data(), coeffArgsStr.Data()));
      //coeff->Print("v");
      _coeffs.add(*coeff);
    }
    //_specialSet.Print("v");
    // restore myVarSec
    setVarSec(myVarSec);
  }
  // build pdf function
  _thePdf=new RooAddPdf(Form("the_%s", GetName()), _pdfType+" "+GetTitle(),
			_subPdfs, _coeffs);
  //_thePdf->Print("v");
  cout<<"done init of rarMLPdf for "<<GetName()<<endl<<endl;
}

/// \brief Return the special ArgSet for splitting
///
/// \param setName The name of the ArgSet
/// \return The special ArgSet for splitting
///
/// It returns the named special ArgSet for splitting
RooArgSet rarMLPdf::getSpecialSet(TString setName)
{
  if ("cat1Set"==setName) return _cat1Set;
  if ("asymSet"==setName) return _asymSet;
  return _specialSet; // default
}

/// \brief Return prototype var generator for toy study
/// \return The generator
///
/// It constrcuts prototype var generator for toy study
RooAbsPdf *rarMLPdf::getProtGen()
{
  if (_theProtGen) return _theProtGen;
  // call super class' getProtGen
  rarCompBase::getProtGen();
  // construct prototype var generator from components
  RooArgList compList;
  for (Int_t i=0; i<_nComp; i++) {
    if (_protGenPdfs.getSize()>0) {
      RooArgList pdfList(_protGenPdfs);
      pdfList.add(_subProtGenPdfs[i]);
      TString compName=
	Form("prodProtGen_%s_%s",GetName(),_subPdfs[i].GetName());
      TString compTitle=
	Form("Prod Prot Gen %s %s",GetName(),_subPdfs[i].GetName());
      RooAbsPdf *thePdf=new RooProdPdf (compName, compTitle, pdfList);
      compList.add(*thePdf);
    } else compList.add(_subProtGenPdfs[i]);
  }
  _theProtGen=new RooAddPdf
    (Form("gen_%s", GetName()), Form("Gen Add for %s", GetName()),
     compList, _coeffs);
  
  cout<<"getProtGen: Expected number of events "
      <<_theProtGen->expectedEvents(0)<<endl;
  
  cout<<"protGen created for "<<GetName()<<endl;
  _theProtGen->Print();
  //_theProtGen->Print("v");
  return _theProtGen;
}

/// \brief Do pdfPlot for given PDFs
///
/// \param plotList List of plots
/// \param pdfList List of PDFs to be plotted
/// \return The last RooPlot created
///
/// The function of #rarCompBase should be used
/// to get pdfPlot instead of that of #rarAdd, its direct super-class,
/// because it is more useful to get individual plots of its components
/// rather than the total plot of its own.
RooPlot *rarMLPdf::doPdfPlot(TList &plotList, TString pdfList)
{
  return rarCompBase::doPdfPlot(plotList, pdfList);
}
