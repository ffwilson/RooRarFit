/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarMLFitter.cc,v 1.148 2014/09/14 17:33:52 fwilson Exp $
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
#include <sstream>
#include <map>
#include <set>
#include <vector>
using namespace std;

#include <libgen.h>
#include "TFile.h"
#include "TTree.h"
#include "TObjString.h"
#include "TStopwatch.h"

#include "Roo1DTable.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooChi2Var.h"
#include "RooHistError.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include "RooNLLVar.h"
#include "RooPlot.h"
#include "RooProdPdf.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooSimPdfBuilder.h"
#include "RooStringVar.h"
#include "RooSuperCategory.h"
#include "RooGlobalFunc.h"
#include "RooPullVar.h"

using namespace RooFit;

#include "rarMinuit.hh"
#include "rarMLPdf.hh"
#include "rarNLL.hh"
#include "rarSPlot.hh"

#include "rarMLFitter.hh"
#include "rarToyList.hh"

ClassImp(rarMLFitter)
  ;

TString rarMLFitter::_physCatStr="";
TString rarMLFitter::_splitCatStr="";
RooArgSet rarMLFitter::_splitCatSet;
RooArgSet rarMLFitter::_splitCatSet2;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarMLFitter::rarMLFitter()
  : rarCompBase(),
    _simBuilder(0), _simConfig(0), _theGen(0), _protGenLevel(0),
    _protDataset(0), _theToyParamGen(0), _theSPdf(0), _theBPdf(0),
    _toyID(0), _toyNexp(0)
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
rarMLFitter::rarMLFitter(const char *configFile, const char *configSec,
			 const char *configStr, rarDatasets *theDatasets,
			 RooDataSet *theData, const char*name,const char*title)
  : rarCompBase(configFile, configSec, configStr,
		theDatasets, theData, name, title, kFALSE),
    _simBuilder(0), _simConfig(0), _theGen(0), _protGenLevel(0),
    _protDataset(0), _theToyParamGen(0), _theSPdf(0), _theBPdf(0),
    _toyID(0), _toyNexp(0)
{
  init();
}

rarMLFitter::~rarMLFitter()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It cheks to see if it is needed to build SimPdf for the final ml model.
/// If yes, it will build the SimPdf using RooSimPdfBuilder.
void rarMLFitter::init()
{
  //if (!getFitter()) setFitter(this);
  if (!getControlBit("SimFit")) // use the first in the prototype pdf list
    _thePdf=(RooAbsPdf*)_subPdfs.at(0);
  else if (""!=readConfStr("category", "", getVarSec())) {
    TString catStr=readConfStr("category", "", getVarSec());
    _category=(RooAbsCategoryLValue *)_fullObs->find(catStr);
    if (!_category) {
      cout<<"Can not find category named "<<catStr<<endl;
      exit(-1);
    }
    _thePdf=new RooSimultaneous(Form("the_%s",GetName()),
				_pdfType+" "+GetTitle(),
				_subPdfs, *_category);
  } else {
    // create SimPdfBuilder
    _simBuilder=new RooSimPdfBuilder(_subPdfs);
    _simConfig=_simBuilder->createProtoBuildConfig();
    _simConfig->readFromFile(_configFile, 0, getVarSec());
    cout<<"simPdfBuilder configs in the config file:"<<endl;
    _simConfig->Print("v");
    // read in special config string
    setSpecialStr(_simConfig);
    cout<<"simPdfBuilder configs after reading in special strings:"<<endl;
    _simConfig->Print("v");
    // do we need to use _simBuilder->addSpecializations?
    // Yes, let's do it
    _simBuilder->addSpecializations(getSpecialSet());
    cout<<"special Set"<<endl;
    getSpecialSet().Print("v");
    if (!_theData) {
      cout<<" No dataset defined to build simPdf."
          <<" Please add/edit"<<endl
          <<"fitData=<datasetName>"<<endl
          <<" in master pdf section ["<<getMasterSec()<<"]"<<endl;
      exit(-1);
    }
    _thePdf=(RooAbsPdf*)_simBuilder->buildPdf(*_simConfig, _theData, 0, kTRUE);
    if (!_thePdf) {
      cout<<" simPdf not built!!"<<endl
          <<" Please check your simPdf config in section ["
          <<getVarSec()<<"]"<<endl;
      exit(-1);
    }
    //_thePdf->Print();
    //_thePdf->Print("v");
    setSimPdf((RooSimultaneous*)_thePdf);
  }
  { // create extra pdfs as toy param randomizer
    Int_t xPdfs=_xPdfList.GetSize();
    createPdfs("preToyRandGenerators",&_xPdfList,&_preToyRandGenerators,
	       _runSec);
    Int_t xPdfsAfter=_xPdfList.GetSize();
    if (xPdfsAfter>xPdfs) {//set the default pdfFit behavior for those xtraPdfs
      for (Int_t i=xPdfs; i<xPdfsAfter; i++) {
	rarBasePdf *thePdf=(rarBasePdf*)_xPdfList.At(i);
	thePdf->setControlBit("noPdfFit", "pdfFit");
	thePdf->setControlBit("noPdfPlot", "pdfPlot");
      }
      // create the randomizer
      _theToyParamGen=new RooProdPdf
	(Form("paraRand_%s", GetName()),"paraRand",_preToyRandGenerators);
    }
  }
  // create extra pdfs as embed obs randomizer
  rarStrParser postEmbdRandObsParser=readConfStr("postEmbdRandObs","",_runSec);
  TString embdObsGensStr="";
  while (postEmbdRandObsParser.nArgs()>0) {
    // the first is data src name
    TString dataSrcName=postEmbdRandObsParser[0];
    postEmbdRandObsParser.Remove();
    if (_embdObsRandSet.find(dataSrcName)) {
      cout<<dataSrcName<<" has been defined in postEmbdRandObs"<<endl;
      exit(-1);
    }
    RooStringVar *srcRandStr=new RooStringVar(dataSrcName,dataSrcName,"",8192);
    _embdObsRandSet.addOwned(*srcRandStr);
    while (postEmbdRandObsParser.nArgs()>0) {
      TString obsName=postEmbdRandObsParser[0];
      postEmbdRandObsParser.Remove();
      srcRandStr->setVal(Form("%s %s",srcRandStr->getVal(), obsName.Data()));
      // check if obsName is obs or pdf
      if (_thePdf->getDependents(*_fullObs)->find(obsName)) continue;
      embdObsGensStr+=" "+obsName;
      break;
    }
  }
  if (""!=embdObsGensStr) {
    setConfStr("embdObsGensStr", embdObsGensStr);
    Int_t xPdfs=_xPdfList.GetSize();
    createPdfs("embdObsGensStr",&_xPdfList,&_embdObsGens);
    Int_t xPdfsAfter=_xPdfList.GetSize();
    if (xPdfsAfter>xPdfs) {//set the default pdfFit behavior for those xtraPdfs
      for (Int_t i=xPdfs; i<xPdfsAfter; i++) {
	rarBasePdf *thePdf=(rarBasePdf*)_xPdfList.At(i);
	thePdf->setControlBit("noPdfFit", "pdfFit");
	thePdf->setControlBit("noPdfPlot", "pdfPlot");
      }
    }
  }
  
  if (_protDataVars.getSize()>0) { // display protDataVars if any
    cout<<"protDataVars (including conditionalObs) defined:"<<endl;
    _protDataVars.Print("v");
  }
  if (_conditionalObs.getSize()>0) { // display conditionalObs if any
    cout<<"conditionalObs defined:"<<endl;
    _conditionalObs.Print("v");
  }
  
  { // get protDataEVars defined in config file
    RooArgSet allProtEVars(getArgSet("protDataEVars", kTRUE));
    allProtEVars.add(getArgSet(readConfStr("protDataEVars","",_runSec),kFALSE));
    if (allProtEVars.getSize() > 0) {
      cout<<"protDataEVars defined:"<<endl;
      allProtEVars.Print("v");
    } else {
      cout << "No protDataEVars defined." << endl;
    }
    RooArgList allProtEVarList(allProtEVars);
    for (Int_t i=0; i<allProtEVarList.getSize(); i++) {
      RooAbsArg *theEVar=_fullObs->find(allProtEVarList[i].GetName());
      if (theEVar) _protDataEVars.add(*theEVar);
    }
  }
  
  cout<<"done init of rarMLFitter for "<<GetName()<<endl<<endl;
}

/// \brief Set special config string for yield splitting
/// \param simConfigSet Config string ArgSet
///
/// It calls rarMLPdf::getSpecialStr to get special config string
/// for all its components and append it to appropriate config string.
void rarMLFitter::setSpecialStr(RooArgSet *simConfigSet)
{
  for (Int_t i=0; i<_nComp; i++) {
    rarMLPdf *model=(rarMLPdf*)_pdfList.At(i);
    RooStringVar* modelSplitStr=(RooStringVar*)
      simConfigSet->find(Form("the_%s",model->GetName()));
    modelSplitStr->setVal(Form("%s %s", modelSplitStr->getVal(),
			       (model->getSpecialStr()).Data()));
    // because of bugs in RooFit, let's remove [] stuff
    TString splitVar=modelSplitStr->getVal();
    Int_t lIdx, rIdx;
    while (((lIdx=splitVar.First("["))>-1)&&
	   ((rIdx=splitVar.First("]"))>-1)) {
      splitVar.Replace(lIdx, rIdx-lIdx+1, "");
    }
    modelSplitStr->setVal(splitVar.Data());
  }
}

/// \brief Return the special ArgSet for splitting
/// \param setName The name of the ArgSet
/// \return The special ArgSet for splitting
///
/// It calls rarMLPdf::getSpecialSet to get special ArgSet
/// for all its components.
RooArgSet rarMLFitter::getSpecialSet(TString setName)
{
  RooArgSet specialSet("special ArgSet for splitting");
  for (Int_t i=0; i<_nComp; i++) {
    specialSet.add(((rarMLPdf*)_pdfList.At(i))->getSpecialSet(setName));
  }
  return specialSet;
}

/// \brief Calculate values of splitting coeff fraction
///
/// \param coeffList ArgList of splitting coeff fraction
/// \param valType value type for calculation (Frac, Asym)
/// \param o The output stream
///
/// It calculates values for these splitting coeff fraction.
/// - Fraction. One coeff fraction is \p RooFormulaVar
///   so it is desirable to show its value for each cat.
/// - Asymmetry. If the number of types for a cat is two,
///   asym is also calculated.
void rarMLFitter::getSplitCoeffValues(RooArgList coeffList, TString valType,
				      ostream &o)
{
  for (Int_t i=0; i<coeffList.getSize(); i++) {
    RooRealVar *theCoeffFrac=(RooRealVar*)(&coeffList[i]);
    Double_t val=theCoeffFrac->getVal();
    if ("Frac"==valType) {
      o<<theCoeffFrac->GetName()<<" = "<<val<<endl;
    } else if ("Asym"==valType) {
      Double_t asym=2*val-1;
      Double_t err=2*theCoeffFrac->getError();
      Double_t errHi=2*theCoeffFrac->getAsymErrorHi();
      Double_t errLo=2*theCoeffFrac->getAsymErrorLo();
      o<<"Asym wrt "<<theCoeffFrac->GetName()<<endl
       <<" 1-2*("<<theCoeffFrac->GetName()<<") = "
       <<-asym<<" +/- "<<err<<" (+"<<errLo<<", -"<<errHi<<")"<<endl;
      //<<" 2*("<<theCoeffFrac->GetName()<<")-1 = "
      //<< asym<<" +/- "<<err<<" (+"<<errHi<<", -"<<errLo<<")"<<endl;
    }
  }
}

/// \brief Return the phys cat name
/// \return The cat name for phys model splitting
TString rarMLFitter::getPhysCat()
{
  if (""!=_physCatStr) return _physCatStr;
  rarStrParser physModelsParser=readConfStr("physModels", "", getMasterSec());
  _physCatStr=physModelsParser[0];
  _physCatStr.ReplaceAll(":", "");
  RooAbsCategory *theCat=(RooAbsCategory*)(getCats()->find(_physCatStr));
  if (!theCat) _physCatStr="";
  return _physCatStr;
}

/// \brief Return splitCat string (w/o physCat)
/// \return The splitCat string
TString rarMLFitter::getSplitCats()
{
  if (""!=_splitCatStr) return _splitCatStr;
  // first get splitCatSet
  getSplitCatSet();
  // loop splitCat set to fill the string
  Int_t nCat=_splitCatSet.getSize();
  if (nCat<=0) {
    _splitCatStr=" ";
    return _splitCatStr;
  }
  RooArgList splitCatList(_splitCatSet);
  for(Int_t i=0; i<nCat; i++) {
    _splitCatStr+=splitCatList[i].GetName();
    _splitCatStr+=" ";
  }
  _splitCatStr.Remove(_splitCatStr.Length()-1);
  
  return _splitCatStr;
}

/// \brief Return splitCats (w/o physCat)
/// \return The splitCat set
RooArgSet rarMLFitter::getSplitCatSet()
{
  if (_splitCatSet.getSize()>0) return _splitCatSet;
  // first remove physCat from splitCat
  TString splitCats=readConfStr("splitCats", "", getMasterSec());
  TString physCat=getPhysCat();
  if (""!=physCat) splitCats.ReplaceAll(physCat, "");
  rarStrParser catsStrParser=splitCats;
  Int_t nCat=catsStrParser.nArgs();
  for (Int_t i=0; i<nCat; i++) {
    TString catName=catsStrParser[i];
    // we do not want anything in between "(" and ")"
    Int_t lIdx, rIdx;
    while (((lIdx=catName.First("("))>-1)&&
	   ((rIdx=catName.First(")"))>-1)) {
      catName.Replace(lIdx, rIdx-lIdx+1, "");
    }
    if (""==catName) continue;
    RooAbsCategory *theCat=(RooAbsCategory*)(getCats()->find(catName));
    if (!theCat) {
      cout<<"Can not find cat "<<catName<<endl
	  <<"Please make sure the cat you are using exists"<<endl;
      exit(-1);
    }
    // do we only use part of the cat?
    lIdx=catsStrParser[i].First("(");
    rIdx=catsStrParser[i].First(")");
    if (lIdx<0 && rIdx<0) { // no
      _splitCatSet.add(*theCat);
      continue;
    }
    if (lIdx>-1 && rIdx<0) {
      cout<<" Invalid specif. w/ config splitCats in sec. "<<getMasterSec()
	  <<endl;
      exit(-1);
    }
    if (lIdx<0 && rIdx>-1) {
      cout<<" Invalid specif. w/ config splitCats in sec. "<<getMasterSec()
	  <<endl;
      exit(-1);
    }
    // now create a new cat based on it and add it to _splitCatSet
    TString catTypesStr=catsStrParser[i](lIdx+1, rIdx-lIdx-1);
    catTypesStr.ReplaceAll(",", " ");
    rarStrParser catTypesStrParser=catTypesStr;
    if (catTypesStrParser.nArgs()<1) continue; // nothing in the cat
    // create a new cat
    RooCategory *thePCat=new RooCategory(theCat->GetName(),theCat->GetTitle());
    // add types into the cat
    for (Int_t k=0; k<catTypesStrParser.nArgs(); k++) {
      // make sure the name is a cat type
      const RooCatType *theType=theCat->lookupType(catTypesStrParser[k]);
      if (!theType) {
	cout <<"The cat type "<<catTypesStrParser[k]
	     <<" does not exist in cat "<<theCat->GetName()<<endl;
	exit(-1);
      }
      thePCat->defineType(theType->GetName(), theType->getVal());
    }
    _splitCatSet.add(*thePCat);
  }
  
  //RooArgList l(_splitCatSet);
  //for (Int_t i=0; i<l.getSize(); i++)
  //l[i].Print("v");
  return _splitCatSet;
}

/// \brief Return splitCat (maybe derived)
/// \param splitCatSet ArgSet to extract the cats
/// \param catName String for cats
/// \return The splitCat
RooAbsCategory *rarMLFitter::getSplitCat(RooArgSet &splitCatSet,
					 TString catName)
{
  TString retName="";
  RooAbsCategory *retVal(0);
  RooArgSet thisSet;
  rarStrParser catNameParser=catName;
  while (catNameParser.nArgs()>0) {
    TString theName=catNameParser[0];
    catNameParser.Remove();
    if (!splitCatSet.find(theName)) {
      cout<<" Can not find "<<theName<<" in split cats "<<endl;
      return retVal;
    }
    thisSet.add(*splitCatSet.find(theName));
    retName+=theName+" ";
  }
  retName.Remove(retName.Length()-1);
  retName.ReplaceAll(" ", "_");
  // remove found cats
  splitCatSet.remove(thisSet);
  
  // do we need to create super cat?
  if (1==thisSet.getSize()) {
    retVal=(RooAbsCategory*)thisSet.find(retName);
    return retVal;
  }
  // do we have the super cat?
  retVal=(RooAbsCategory*)_splitCatSet2.find(retName);
  if (retVal) return retVal;
  // create super cat
  retVal=new RooSuperCategory(retName, retName, thisSet);
  // add it to the argset
  _splitCatSet2.add(*retVal);

  return retVal;
}

/// \brief Return a root file name
///
/// \param aType Action type
/// \param configToken config token
/// \return Root file name constructed
TString rarMLFitter::getRootFileName(TString aType, TString configToken)
{
  TString rootFile;
  if (("yes"==configToken)||("default"==configToken)) // default name
    rootFile=getFullFileName(_resultDir,aType,_runSec,"none", getMasterSec());
  else {
    char *bName=basename((char*)configToken.Data());
    rootFile=_resultDir+"/"+bName;
  }
  if (!rootFile.EndsWith(".root")) rootFile+=".root";
  if (_toyID) {
    rootFile.Replace(rootFile.Length()-5, 5, ".%03d.root");
    rootFile=Form(rootFile.Data(), _toyID);
  }
  
  return rootFile;
}

/// \brief Return a param file name
///
/// \param aType Action type
/// \param configToken config token
/// \param dirName Common dir
/// \return Param file name constructed
TString rarMLFitter::getParamFileName(TString aType, TString configToken,
				      TString dirName)
{
  TString paramFile;
  TString name="none";
  TString dsName="";
  TString msName="";
  TString cfName=_configFile;
  rarStrParser paramIDParser=configToken;
  while (paramIDParser.nArgs()>1) {
    char flag=paramIDParser[0][0];
    paramIDParser.Remove();
    switch (flag) {
    case 'F' :
      cfName=paramIDParser[0];
      paramIDParser.Remove();
      break;
    case 'D' :
      dsName=paramIDParser[0];
      paramIDParser.Remove();
      break;
    case 'C' :
      msName=paramIDParser[0];
      paramIDParser.Remove();
      break;
    case 'A' :
      aType=paramIDParser[0];
      paramIDParser.Remove();
      break;
    case 'N' :
      name=paramIDParser[0];
      paramIDParser.Remove();
      break;
    }
  }
  
  if (""==dirName) dirName=_paramDir;
  return getFullFileName(dirName, aType, name, dsName, msName, cfName);
}

/// \brief Read in/write out params
/// \param params Argset to IO
/// \param paramFile param file
/// \param In Input/ouput
///
/// This is param IO function
void rarMLFitter::paramFileIO(RooArgSet params, TString paramFile, Bool_t In)
{
  // add correlation coeffs
  RooArgSet corrCoeffs(getCorrCoeffs());
  TString paramOrder=readConfStr("outParamOrder", "ascend", getMasterSec());
  if (("ascend"==paramOrder) || ("descend"==paramOrder)) {
    Bool_t reverse=("descend"==paramOrder);
    RooArgList paramList(corrCoeffs);
    corrCoeffs.removeAll();
    paramList.sort(reverse);
    corrCoeffs.add(paramList);
  }
  if (!In) { // remove trivial corr coeffs for output
    RooArgSet trivSet;
    TIterator* iter = corrCoeffs.createIterator();
    RooRealVar *theCorrCoef(0);
    while(theCorrCoef=(RooRealVar*)iter->Next()) {
      Double_t theVal=theCorrCoef->getVal();
      if ((fabs(theVal)>1)||(fabs(theVal)<1e-4))
	trivSet.add(*theCorrCoef);
    }
    delete iter;
    corrCoeffs.remove(trivSet);
  }
  params.add(corrCoeffs);
  // construct param info
  RooStringVar configFile("configFile", "configFile", _configFile);
  params.add(configFile);
  RooStringVar DataInputSec("DataInputSec", "DataInputSec",
			    getDatasets()->getVarSec());
  params.add(DataInputSec);
  RooStringVar MasterPdfSec("MasterPdfSec", "MasterPdfSec", getMasterSec());
  params.add(MasterPdfSec);
  RooStringVar ActionSec("ActionSec", "ActionSec", _runSec);
  params.add(ActionSec);
  
  if (!In) params.writeToFile(paramFile);
  else {
    Bool_t fail=params.readFromFile(paramFile);
    if (fail) {
      cout<<"Fail to read in from "<<paramFile<<endl;
      exit(-1);
    }
    cout<<"Info of params you just read in"<<endl
	<<"Config  file  : "<<configFile.getVal()<<endl
	<<"Data input Sec: "<<DataInputSec.getVal()<<endl
	<<"Master pdf Sec: "<<MasterPdfSec.getVal()<<endl
	<<"Action     Sec: "<<ActionSec.getVal()<<endl;
  }
}

/// \brief Read in params
/// \param params Argset to input
/// \param paramFile param file
///
/// This is param input function
void rarMLFitter::paramFileI(RooArgSet params, TString paramFile)
{
  paramFileIO(params, paramFile);
}

/// \brief Write out params
/// \param params Argset to output
/// \param paramFile param file
///
/// This is param output function
void rarMLFitter::paramFileO(RooArgSet params, TString paramFile)
{
  paramFileIO(params, paramFile, kFALSE);
}

/// \brief Check if the dataset is blind or not
/// \param dsName The name of the dataset
///
/// It check if the dataset is blind.
void rarMLFitter::chkBlind(TString dsName)
{
  if (_datasets->isBlind(dsName)) {
    cout<<endl<<" Dataset \""<<dsName<<"\" is still blind."<<endl
	<<" (any action other than pdfFit or toyStudy on the dataset"
	<<" requires its unblind)"<<endl
	<<" (unblinding dataset does not necessarily"
	<<" unblind your results,)"<<endl
	<<" (as long as you use RooUnblindPrecision or RooUnblindOffset"
	<<" and the state is set to blind,)"<<endl
	<<" (those variables they blind still remain blind.)"<<endl
	<<" If you really want to UNBLIND "<<dsName<<","<<endl
	<<" please specify in dsi section ["
	<<_datasets->getVarSec()<<"]:"<<endl<<endl
	<<"  ub_"<<dsName<<" = [<OldIDsIfAny>] "<<_datasets->ubStr(dsName)
	<<endl<<endl<<" and re-run the job"<<endl;
    //    exit(-1);
  }
  return;
}

/// \brief ML fit driver
/// \param mlFitData Dataset to fit for mlFit action
/// \param opt Fitting options (default: "ehr").
///            Changeable by config item \p mlFitOption in run section
/// \param o The output stream
/// \param ncpus Number of CPUs (cores) to use in fit (default=1)
/// \return The fit result object
///
/// It is a wrapper for final ml fit.
RooFitResult *rarMLFitter::doMLFit(RooDataSet *mlFitData,
				   TString opt, ostream &o, Int_t ncpus)
{
  cout<<endl<<" In rarMLFitter doMLFit for "<<GetName()<< " Options: " << opt 
      << " using " << ncpus << " CPUs (if available)." << endl;
  RooFitResult *fitResult(0);
  // first check if it has its pdf
  assert(_thePdf);
  
  // convert to newer fitTo format for steering options
  Bool_t mlFitExtended = opt.Contains("e");
  Bool_t mlFitMinos    = opt.Contains("m");
  Bool_t mlFitHesse    = opt.Contains("h");
  Bool_t mlFitVerbose  = !opt.Contains("q");
  Bool_t mlFitSave     = opt.Contains("r"); // return results

  //_thePdf->Print();
  //fit to its dataset
  fitResult=_thePdf->fitTo(*mlFitData, ConditionalObservables(_conditionalObs),
  			   Save(mlFitSave), Extended(mlFitExtended), 
			   Verbose(mlFitVerbose), Hesse(mlFitHesse),
			   Minos(mlFitMinos), NumCPU(ncpus));

  // needed to fill "GblCorr." column of printout (now redundant?) FFW
  fitResult->globalCorr();
  // save correlation coeffs from final fit
  saveCorrCoeffs(fitResult);
  // output the result
  o<<endl<<"The mlFit results (wrt "<<mlFitData->GetName()
   <<" with " << mlFitData->numEntries() <<" events):"<<endl;
#ifndef USENEWROOT 
  fitResult->printToStream(o, RooPrintable::Verbose);
#else
  Int_t contents(0);
  fitResult->printMultiline(o, contents, kTRUE);
#endif
  // calculate frac in cat1Set and asym in asymSet
  getSplitCoeffValues(getSpecialSet("cat1Set"), "Frac", o);
  getSplitCoeffValues(getSpecialSet("asymSet"), "Asym", o);
  
  return fitResult;
}

/// \brief Wrapper for ML fit calls
/// \param pdf Dataset to fit for mlFit action
/// \param fitData Dataset to fit for mlFit action
/// \param fitOptions Fitting options (default: "ehr").
/// \param ncpus Number of CPUs (cores) to use in fit (default=1)
/// \return The fit result object
///
/// It is a wrapper for calls to fiTo() method.
RooFitResult *rarMLFitter::doTheFit(RooAbsPdf *pdf, RooDataSet *fitData, TString fitOptions, Int_t ncpus)
{
  cout<<endl<<"In rarMLFitter doTheFit for "<<GetName()<< " Options: " << fitOptions 
      << " using " << ncpus << " CPUs (if available)." << endl;

  // first check if it has its pdf
  assert(pdf);
  
  // convert to newer fitTo format for steering options
  Bool_t fitExtended = fitOptions.Contains("e");
  Bool_t fitMinos    = fitOptions.Contains("m");
  Bool_t fitHesse    = fitOptions.Contains("h");
  Bool_t fitVerbose  = !fitOptions.Contains("q");
  Bool_t fitSave     = fitOptions.Contains("r"); // return results

  TStopwatch timer;
  timer.Start();
  // fit to its dataset
  RooFitResult *fitResult=_thePdf->fitTo(*fitData, ConditionalObservables(_conditionalObs),
					 Save(fitSave), Extended(fitExtended), 
					 Verbose(fitVerbose), Hesse(fitHesse),
					 Minos(fitMinos), NumCPU(ncpus));
  timer.Stop();
  // output the time
  cout<<endl<<"The doTheFit RealTime= " << timer.RealTime() << " CpuTime= "<< timer.CpuTime() << endl;
  
  return fitResult;
}

/// \brief Significance calculation
/// \param mlFitData Dataset to fit for mlFit action
/// \param signfStr String of params to calculate significance
/// \param fitResult The fit result from nominal fit
/// \param fullParams ArgSet of all params
/// \param o The output stream
///
/// It performs significance calculation
void rarMLFitter::doSignf(RooDataSet *mlFitData, TString signfStr,
			  RooFitResult *fitResult, RooArgSet fullParams,
			  ostream &o)
{
  cout<<endl<<" In rarMLFitter doSignf for "<<GetName()<<endl;
  
  // save params
  string fParamSStr0;
  writeToStr(fullParams, fParamSStr0);
  
  TString signfFitOpt="qemhr";
  // convert to newer fitTo format for steering options
  Bool_t signfFitExtended = signfFitOpt.Contains("e");
  Bool_t signfFitMinos    = signfFitOpt.Contains("m");
  Bool_t signfFitHesse    = signfFitOpt.Contains("h");
  Bool_t signfFitVerbose  = !signfFitOpt.Contains("q");
  Bool_t signfFitSave     = signfFitOpt.Contains("r"); // return results

  // get nll
  if (!fitResult) {
    cout<<"No fit results from mlFit"<<endl;
    exit(-1);
    return;
  }
  Double_t nll=2*fitResult->minNll();
  o<<endl;
  rarStrParser signfStrParser=signfStr;
  while (signfStrParser.nArgs()>0) {
    TString paramName=signfStrParser[0];
    signfStrParser.Remove();
    RooRealVar *theParam=(RooRealVar*)fullParams.find(paramName);
    if (!theParam) continue;
    // restore params
    readFromStr(fullParams, fParamSStr0);
    Double_t zSignfVal(0);
    if ((signfStrParser.nArgs()>0)&&(isNumber(signfStrParser[0]))) {
      zSignfVal = atof(signfStrParser[0]);
      signfStrParser.Remove();
    }
    Double_t sVal=theParam->getVal();
    theParam->setVal(zSignfVal);
    theParam->setConstant();

    // fit again
    //Int_t ncpus(1);

    RooFitResult *theResult=
      _thePdf->fitTo(*mlFitData, ConditionalObservables(_conditionalObs), 
    		     Save(signfFitSave), Extended(signfFitExtended), 
    		     Verbose(signfFitVerbose), Hesse(signfFitHesse),
		     Minos(signfFitMinos));
    //RooFitResult *theResult = doTheFit(_thePdf, mlFitData, signfFitOpt, ncpus);

    o<<" Signf. of "<<theParam->GetName()<<" being "<<sVal
     <<" wrt "<<zSignfVal<<" is "<<sqrt(2*theResult->minNll()-nll)
     <<" (sigma)"<<endl;
  }
  
  // restore saved params
  readFromStr(fullParams, fParamSStr0);
  return;
}

/// \brief Systematic error study driver
/// \param mlFitData Dataset to fit for mlFit action
/// \param paramsStr String of params to vary
/// \param varsStr String of params to study
/// \param fullParams ArgSet of all params
/// \param o The output stream
///
/// It performs systematic error study.
void rarMLFitter::doSysStudy(RooDataSet *mlFitData, TString paramsStr,
			     TString varsStr, RooArgSet fullParams,
			     ostream &o)
{
  cout<<endl<<" In rarMLFitter doSysStudy for "<<GetName()<<endl;
  
  // save params
  string fParamSStr0;
  writeToStr(fullParams, fParamSStr0);
  // first find out max length for param name
  Int_t paramNameLen=12;
  {
    rarStrParser paramsStrParser=paramsStr;
    for (Int_t i=0; i<paramsStrParser.nArgs(); i++) {
      if (paramNameLen<paramsStrParser[i].Length())
	paramNameLen=paramsStrParser[i].Length();
    }
  }
  // ArgSet of study vars
  RooArgSet studyVars;
  rarStrParser varsStrParser=varsStr;
  while (varsStrParser.nArgs()>0) {
    TString varName=varsStrParser[0];
    varsStrParser.Remove();
    RooRealVar *theVar=(RooRealVar*)fullParams.find(varName);
    if (theVar) studyVars.add(*theVar);
  }
  if (studyVars.getSize()<1) {
    cout<<" No vars in postMLSysVars found!"<<endl;
    return;
  }
  cout<<" postMLSysVars:"<<endl;
  studyVars.Print("v");
  // fit again with option emhr
  TString sysFitOpt="qemhr";

  //Int_t ncpus(1);

  // convert to newer fitTo format for steering options
  Bool_t sysFitExtended = sysFitOpt.Contains("e");
  Bool_t sysFitMinos    = sysFitOpt.Contains("m");
  Bool_t sysFitHesse    = sysFitOpt.Contains("h");
  Bool_t sysFitVerbose  = !sysFitOpt.Contains("q");
  Bool_t sysFitSave     = sysFitOpt.Contains("r"); // return results

  _thePdf->fitTo(*mlFitData, ConditionalObservables(_conditionalObs), 
		 Save(sysFitSave), Extended(sysFitExtended), 
		 Verbose(sysFitVerbose), Hesse(sysFitHesse),
		 Minos(sysFitMinos));
  //  doTheFit(_thePdf, mlFitData, sysFitOpt, ncpus);

  string fParamSStr;
  writeToStr(fullParams, fParamSStr);
  // clone it
  RooArgSet *cStudyVars=(RooArgSet*)studyVars.snapshot(kTRUE);
  // vary params
  rarStrParser paramsStrParser=paramsStr;
  // default variant unit
  Double_t defVU(1);
  if (paramsStrParser.nArgs()>0) {
    if (isNumber(paramsStrParser[0])) {
      defVU=fabs(atof(paramsStrParser[0]));
      paramsStrParser.Remove();
    }
  }
  Int_t nSysParams(0);
  // string for all params varied
  TString vParams="";
  TArrayD pV(1), mV(1);
  TArrayD pArray(studyVars.getSize()), mArray(studyVars.getSize());
  TArrayD aArray(studyVars.getSize()); // avg error
  // loop over all params
  while (paramsStrParser.nArgs()>0) {
    TString theParamName=paramsStrParser[0];
    paramsStrParser.Remove();
    RooArgSet theParamSet;
    TIterator* iter = fullParams.createIterator();
    RooRealVar *theParam(0);
    Bool_t foundvParam=kFALSE;
    while(theParam=(RooRealVar*)iter->Next()) {
      TString theName=theParam->GetName();
      if ((theName==theParamName)||(theName.BeginsWith(theParamName+"_"))) {
	theParamSet.add(*theParam);
        if (!foundvParam) {
          foundvParam=kTRUE;
          vParams=vParams+" \""+theParamName+"\"";
        }
      }
    }
    delete iter;
    if (theParamSet.getSize()<1) continue;
    nSysParams++;
    // reset arrays
    pV.Set(nSysParams);
    mV.Set(nSysParams);
    pArray.Set(nSysParams*studyVars.getSize());
    mArray.Set(nSysParams*studyVars.getSize());
    aArray.Set(nSysParams*studyVars.getSize());
    // for this param
    cout<<" SysStudy for "<<theParamName<<endl;
    theParamSet.Print("v");
    // for plus variation
    Bool_t useErr=kTRUE;
    pV[nSysParams-1]=defVU;
    // do we have plusV specified?
    if (paramsStrParser.nArgs()>0) {
      if (isNumber(paramsStrParser[0])) {
	TString myVStr=paramsStrParser[0];
	paramsStrParser.Remove();
	pV[nSysParams-1]=fabs(atof(myVStr));
	if (myVStr.EndsWith("V")||(myVStr.EndsWith("v"))) useErr=kFALSE;
      }
    }
    // restore params
    readFromStr(fullParams, fParamSStr);
    // set variation
    setVariation(theParamSet, pV[nSysParams-1], useErr, kTRUE);
    // fit

    _thePdf->fitTo(*mlFitData, ConditionalObservables(_conditionalObs), 
   		   Save(sysFitSave), Extended(sysFitExtended), 
		   Verbose(sysFitVerbose), Hesse(sysFitHesse),
		   Minos(sysFitMinos));
   // doTheFit(_thePdf, mlFitData, sysFitOpt, ncpus);

    // calculation errors
    calSysErrors(nSysParams-1, *cStudyVars, studyVars, pArray);
    // for minus variation
    // do we have minusV specified?
    mV[nSysParams-1]=pV[nSysParams-1];
    if (paramsStrParser.nArgs()>0) {
      if (isNumber(paramsStrParser[0])) {
	TString myVStr=paramsStrParser[0];
	paramsStrParser.Remove();
	mV[nSysParams-1]=fabs(atof(myVStr));
	if (myVStr.EndsWith("V")||(myVStr.EndsWith("v"))) useErr=kFALSE;
      }
    }
    // restore params
    readFromStr(fullParams, fParamSStr);
    // set variation
    setVariation(theParamSet, mV[nSysParams-1], useErr, kFALSE);
    // fit
    _thePdf->fitTo(*mlFitData, ConditionalObservables(_conditionalObs), 
    		   Save(sysFitSave), Extended(sysFitExtended), 
		   Verbose(sysFitVerbose), Hesse(sysFitHesse),
		   Minos(sysFitMinos));
    //doTheFit(_thePdf, mlFitData, sysFitOpt, ncpus);
    // calculation errors
    calSysErrors(nSysParams-1, *cStudyVars, studyVars, mArray);
  }
  // get correlation matrix
  TMatrixD corrM1=getCorrMatrix(nSysParams, vParams, kTRUE);
  TMatrixD corrM=getCorrMatrix(nSysParams, vParams);
  // output error table
  o<<endl<<" Systematic Error Table:"<<endl;
  o<<setw(paramNameLen+11)<<" ";
  // output obs names
  RooArgList studyVarList(studyVars);
  for (Int_t i=0; i<studyVarList.getSize(); i++) {
    o<<setw(40)<<studyVarList[i].GetName();
  }
  // output correlation matrix index
  o<<setw(15)<<"corr matrix:";
  rarStrParser vParamsParser=vParams;
  for (Int_t i=0; i<nSysParams; i++) {
    o<<setw(paramNameLen+1)<<vParamsParser[i];
  }
  o<<endl;
  // get avg errros, and output param variations
  avgSysErrors(o, vParams, paramNameLen, pV,pArray, mV,mArray, aArray, corrM);
  // final error output
  outSysErrors(o, "(w/o corr):", paramNameLen, corrM1, aArray);
  outSysErrors(o, "(w/ corr):", paramNameLen, corrM, aArray);
  // restore saved params
  readFromStr(fullParams, fParamSStr0);
  return;
}

/// \brief Set variation for all params specified
/// \param theParams ArgSet of params to vary
/// \param myV Variant
/// \param useErr If the variant is unit of error or not
/// \param isPlus Plus variation
void rarMLFitter::setVariation(RooArgSet theParams, Double_t myV,
			       Bool_t useErr, Bool_t isPlus)
{
  TIterator* iter = theParams.createIterator();
  RooRealVar *theParam(0);
  while(theParam=(RooRealVar*)iter->Next()) {
    Double_t myVariant=myV;
    if (useErr) {
      Double_t myErr(0);
      if (theParam->hasAsymError()) {
	if (isPlus) myErr=theParam->getAsymErrorHi();
	else myErr=theParam->getAsymErrorLo();
      }
      if ((0==myErr)&&(theParam->hasError()))
	myErr=theParam->getError();
      if (0==myErr) {
	cout<<" No err specified for "<<theParam->GetName()<<endl;
	continue;
      }
      myVariant=myV*myErr;
    }
    myVariant=fabs(myVariant);
    if (!isPlus) myVariant=-myVariant;
    cout<<" Variation for "<<theParam->GetName()<<": "<<myVariant<<endl;
    theParam->setVal(theParam->getVal()+myVariant);
  }
  delete iter;
  return;
}

/// \brief Calculate systematic errors
/// \param iParam Index of param currently being studied
/// \param cstudyVars Vars to study (original values)
/// \param studyVars Vars to study
/// \param eArray Array to store errors
void rarMLFitter::calSysErrors(Int_t iParam, RooArgSet &cstudyVars,
                               RooArgSet &studyVars, TArrayD &eArray)
{
  TIterator* iter = studyVars.createIterator();
  RooRealVar *theVar(0), *theIVar(0);
  Int_t vIdx(iParam*studyVars.getSize());
  while(theVar=(RooRealVar*)iter->Next()) {
    theIVar=(RooRealVar*)cstudyVars.find(theVar->GetName());
    if (!theIVar) {
      cout<<" Can not find "<<theVar->GetName()<<endl;
      exit(-1);
    }
    eArray[vIdx++]=theVar->getVal()-theIVar->getVal();
  }
  
  delete iter;
  return;
}

/// \brief Average the plus and minus variations and output systematical errors
/// \param o Output stream
/// \param vParams Variant param names
/// \param pnLen Max length for param name
/// \param pV Positve variations
/// \param pArray Positve errors
/// \param mV Negative variations
/// \param mArray Negative errors
/// \param aArray Average errors
/// \param corrM Correlation matrix
void rarMLFitter::avgSysErrors(ostream &o, TString vParams, Int_t pnLen,
                               TArrayD &pV, TArrayD &pArray,
                               TArrayD &mV, TArrayD &mArray,
                               TArrayD &aArray, TMatrixD corrM)
{
  Int_t nParam=pV.GetSize();
  rarStrParser vParamsParser=vParams;
  if (nParam!=vParamsParser.nArgs()) {
    cout<<vParams<<endl
        <<"does not match # of params: "<<nParam<<endl;
    o<<vParams<<endl
        <<"does not match # of params: "<<nParam<<endl;
  }
  Int_t nVar=pArray.GetSize()/nParam;
  for (Int_t j=0; j<nParam; j++) {
    o.unsetf(ios_base::scientific);
    o.setf(ios_base::showpos);
    o.precision(4);
    o<<setw(pnLen)<<vParamsParser[j]<<setw(5)<<pV[j]<<setw(5)<<-mV[j]<<":";
    o.setf(ios_base::scientific);
    for (Int_t i=0; i<nVar; i++) {
      // output +/- errors
      o<<setw(16)<<pArray[j*nVar+i]
       <<setw(12)<<mArray[j*nVar+i];
      // calculate the avg error
      Double_t aError(0);
      if (pV[j]>0) aError+=pArray[j*nVar+i];
      if (mV[j]>0) aError-=mArray[j*nVar+i];
      if ((pV[j]>0)&&(mV[j]>0)) {
        aError/=2.;
        if (pArray[j*nVar+i]*mArray[j*nVar+i]>0) { // same sign
          cout<<" W A R N I N G !"<<endl
              <<" + and - variation have the same sign!"<<endl
              <<" Using the larger one only!"<<endl;
          aError=pArray[j*nVar+i];
	  if (fabs(aError)<fabs(mArray[j*nVar+i]))
	    aError=mArray[j*nVar+i];
	}
      }
      aArray[j*nVar+i]=aError;
      o.unsetf(ios_base::showpos);
      o<<setw(12)<<aError;
    }
    // now correlation matrix elements
    o<<setw(15)<<" ";
    for (Int_t i=0; i<nParam; i++) {
      o<<setw(pnLen+1)<<corrM(j, i);
    }
    o<<endl;
  }
  
  return;
}

/// \brief Output the final errors
/// \param o Output stream
/// \param rowName Row name
/// \param pnLen Max length for param name
/// \param corrM Correlation matrix for error calculation
/// \param aArray Average errors
void rarMLFitter::outSysErrors(ostream &o, TString rowName, Int_t pnLen,
			       TMatrixD corrM, TArrayD &aArray)
{
  // final error output for the row
  o<<setw(pnLen)<<rowName<<setw(11)<<" ";
  o.setf(ios_base::scientific);
  // # of params
  Int_t nSysParams=corrM.GetNcols();
  Int_t nSysVars=aArray.GetSize()/nSysParams;
  for (Int_t i=0; i<nSysVars; i++) {
    // get error matrics
    TMatrixD errM(nSysParams,1);
    for (Int_t j=0; j<nSysParams; j++) {
      errM(j,0)=aArray[j*nSysVars+i];
    }
    // get transpose of errM
    TMatrixD errMT(1,nSysParams);
    errMT.Transpose(errM);
    // get final error matrix
    TMatrixD errorSq(errMT*corrM*errM);
    o<<setw(40)<<sqrt(errorSq(0,0));
    //errM.Print();
  }
  o<<endl;
  //corrM.Print();
}

/// \brief Return the correlation matrix
/// \param nSysParams Number of varying parameters
/// \param vParams Names of the varying parameters
/// \param diagOnly Diagonal terms only
TMatrixD rarMLFitter::getCorrMatrix(Int_t nSysParams, TString vParams,
				    Bool_t diagOnly)
{
  TMatrixD retM(nSysParams, nSysParams);
  retM.Zero();
  rarStrParser vParamsParser=vParams;
  if (vParamsParser.nArgs()!=nSysParams) {
    cout<<" Names and number of names ("<<nSysParams<<") do not match"<<endl
	<<vParams<<endl;
    return retM;
  }
  for (Int_t i=0; i<nSysParams; i++) {
    for (Int_t j=0; j<nSysParams; j++) {
      if (i==j) {
	retM(i,j)=1;
	continue;
      }
      if (diagOnly) continue;
      retM(i,j)=getCorrCoeff(vParamsParser[i], vParamsParser[j]);
    }
  }
  return retM;
}

/// \brief Chisq GOF study
/// \param mlFitData Dataset to fit for mlFit action
/// \param o The output stream
/// \param plotList Plot list
/// \return GOF chisq
///
/// It does chisq GOF study for mlFit.
Double_t rarMLFitter::doGOFChisq(RooDataSet *mlFitData, ostream &o,
                             TList *plotList)
{
  // first get obs
  RooArgSet GOFObsSet;
  addProtVars("postMLGOFObs", GOFObsSet);
  GOFObsSet.add(*_thePdf->getDependents(mlFitData), kTRUE);
  if (GOFObsSet.getSize()<=0) return 0;
  // save obs
  string obsStr;
  writeToStr(GOFObsSet, obsStr);
  // set # of bins for obs
  RooArgList GOFObsList(GOFObsSet);
  for (Int_t i=0; i<GOFObsList.getSize(); i++) {
    RooRealVar *theVar=dynamic_cast<RooRealVar*>(GOFObsList.at(i));
    if (!theVar) continue;
    Double_t min, max;
    getRange(theVar, "plotRange_", min, max, _runSec);
    Int_t nBins=atoi
      (readConfStr(Form("plotBins_%s",theVar->GetName()),"-1",_runSec));
    if (nBins>0) theVar->setBins(nBins);
  }
  if (plotList) {
    cout<<"GOFObsSet:"<<endl;
    GOFObsSet.Print("v");
  }
  //_thePdf->Print("v");
  //abort();
  // fill data histogram
  RooDataHist dHist("dHist", "dHist", GOFObsSet, *mlFitData);
  if (plotList) {
    //dHist.dump();
    dHist.dump2();
  }
  // calculate chisq using RooChi2Var
  //RooChi2Var chisqVar("chisqVar", "chisqVar", *_thePdf, dHist,
  //_conditionalObs, kTRUE);
  //o<<"GOF chisq: "<<chisqVar.getVal()<<endl;
  // if simpdf, get super cat
  RooSuperCategory *theSuperCat(0);
  RooArgSet inputCats;
  if (getControlBit("SimFit")) {
    theSuperCat=(RooSuperCategory*) &((RooSimultaneous*)_thePdf)->indexCat();
    inputCats.add(theSuperCat->inputCatList());
  }
  // fill pdf histogram
  RooDataHist pHist("pHist", "pHist", GOFObsSet);
  Int_t nBins=pHist.numEntries();
  Double_t sum(0);
  for (Int_t i=0; i<nBins; i++) {
    const RooArgSet *theBinSet=pHist.get(i);
    // set bin for each obs
    RooArgList theBinList(*theBinSet);
    for (Int_t j=0; j<theBinList.getSize(); j++) {
      RooRealVar *theBin=dynamic_cast<RooRealVar*>(theBinList.at(j));
      if (!theBin) {
        RooAbsCategory *theCat=dynamic_cast<RooAbsCategory*>
          (theBinList.at(j));
        if (!theCat) continue;
        RooAbsCategoryLValue *theInCat=(RooAbsCategoryLValue*)
          inputCats.find(theCat->GetName());
        if (!theInCat) continue;
        theInCat->setLabel(theCat->getLabel());
        continue;
      }
      RooRealVar *theVar=dynamic_cast<RooRealVar*>(GOFObsList.at(j));
      assert(theVar);
      Double_t hbSize=(theBin->getMax()-theBin->getMin())/theBin->getBins()/2.;
      theVar->setRange("subrange",
		       theBin->getVal()-hbSize, theBin->getVal()+hbSize);
      theVar->setVal(theBin->getVal());
      //theBin->Print();
      //theVar->Print();
    }
    // the pdf to integrate
    RooAbsPdf *theIntPdf(_thePdf);
    if (theSuperCat) {
      theIntPdf=((RooSimultaneous*)_thePdf)->getPdf(theSuperCat->getLabel());
      //theIntPdf->Print();
      //theIntPdf->Print("v");
    }
    // get integral
    RooAbsReal *pdfIntegral(0);
    if (theIntPdf) pdfIntegral=
      theIntPdf->createIntegral(GOFObsSet, GOFObsSet, "subrange");
    //pdfIntegral->Print("v");
    Double_t thisIntegral(0);
    RooArgSet nullDS;
    if (pdfIntegral) thisIntegral=
      theIntPdf->expectedEvents(&nullDS)*pdfIntegral->getVal();
    if (theSuperCat) thisIntegral/=theSuperCat->numTypes();
    sum+=thisIntegral;
    //cout<<"pdfIntegral="<<thisIntegral<<endl;
    // fill the bin
    pHist.set(*theBinSet, thisIntegral);
    //theBin->Print("v");
  }
  if (plotList) {
    cout<<"sum="<<sum<<endl;
    pHist.dump2();
  }
  // create histograms
  TH1F *dataHist(0);
  TH1F *pdfHist(0);
  if (plotList) {
    dataHist=new TH1F("dataHist_GOF", "dataHist_GOF", nBins, 0, 1);
    pdfHist=new TH1F("pdfHist_GOF", "pdfHist_GOF", nBins, 0, 1);
    plotList->Add(dataHist);
    plotList->Add(pdfHist);
  }
  // calculate chisq by myself
  Double_t chisq(0);
  // fill histograms
  for (Int_t i=0; i<nBins; i++) {
    const RooArgSet *theBinSet=dHist.get(i);
    //theBinSet->Print("v");
    Double_t dw=dHist.weight(*theBinSet,0);
    Double_t ym,yp;
    RooHistError::instance().getPoissonInterval(Int_t(dw+0.5),ym,yp,1);
    Double_t del=dw-ym;
    Double_t deh=yp-dw;
    Double_t pw=pHist.weight(*theBinSet,0);
    if (plotList) {
      dataHist->SetBinContent(i, dw);
      pdfHist->SetBinContent(i, pw);
    }
    Double_t pwmdw=pw-dw;
    Double_t err=(pwmdw>0)?deh:del;
    if (0==err) continue;
    chisq+=pwmdw*pwmdw/(err*err);
    
    if (plotList)
      cout<<"err="<<err<<" dw="<<dw<<" pw="<<pw<<endl;
    
  }
  o<<"GOF chisq: "<<chisq<<endl;
  
  // restore saved obs
  readFromStr(GOFObsSet, obsStr);
  return chisq;
}

/// \brief Toy study driver
/// \param fullParams Full params
/// \return RooDataSet storing fitPar results
///
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_toyStudy">See doc for toy study action.</a>
///
/// The function reads in config options from action section.
/// It finds out prototype dataset,
/// number of experiments, number of events, and if fluctuation is needed.
/// Then it finds out how to generate for each component.
/// And finally it generates, fits, and returns the fitting results
///
/// \todo Make protGen, etc, independent of toy study
/// so they can be used by other routines.
RooDataSet *rarMLFitter::doToyStudy(RooArgSet fullParams)
{
  rarToyList toylist; // keeps tarck of requested events 

  cout<<endl<<" In rarMLFitter doToyStudy for "<<GetName()<<endl;
  // verbose/quiet
  TString temp=readConfStr("toyVerbose", "yes", _runSec);
  Bool_t toyVerbose = kFALSE;
  if (temp == "yes") {toyVerbose = kTRUE;}

  RooDataSet *toyResults(0);
  // get protDatasets from config
  TString protDataStr=readConfStr("protToyData", "", _runSec);
  if (""==protDataStr) protDataStr=readConfStr("protDatasets", "", _runSec);
  rarStrParser protDataStrParser=protDataStr;
  // master protDataset
  if (protDataStrParser.nArgs()>0) {
    TString mpdsName=protDataStrParser[0];
    protDataStrParser.Remove();
    _protDataset=_datasets->getData(mpdsName);
    if (!_protDataset) {
      cout<<" Toy Error: Can not find dataset named "<<mpdsName<<endl
	  <<" for your prototype dataset."<<endl
	  <<" Please make sure you have the right Dataset name"<<endl;
      exit(-1);
    }
  } else if (_theData) _protDataset=_theData;
  else {
    cout<<" Toy Error: No prototype datasets defined!!"<<endl
	<<" You could have problem with your toy study"<<endl
	<<" See config `protDatasets' in online doc for more info"<<endl;
    exit(-1);
  }
  cout<<" Toy Using "<<_protDataset->GetName()<<" as master protDataset"<<endl;

  // string for param randomization
  TString preToyRandParams="";
  // find param randomizer obs
  RooArgSet paramRandSet;
  if (_theToyParamGen) {
    RooArgSet paramNobs(*_thePdf->getParameters(_protDataset));
    paramNobs.add(*_thePdf->getDependents(_protDataset));
    RooArgList allOtherParams(*_theToyParamGen->getParameters(paramNobs));
    for (Int_t i=0; i<allOtherParams.getSize(); i++) {
      RooRealVar *theParam=(RooRealVar*)allOtherParams.at(i);
      if (theParam->hasMin()&&theParam->hasMax()) paramRandSet.add(*theParam);
    }
    // params to be randomized
    preToyRandParams=readConfStr("preToyRandParams", "", _runSec);
  }
  rarStrParser preToyRandParser=preToyRandParams;
  if ((preToyRandParser.nArgs()>0)&&(preToyRandParser[0]=="no"))
    preToyRandParser="";
  
  // get number of toy experiment
  Int_t toyNexp=atoi(readConfStr("toyNexp", "1", _runSec));
  if (_toyNexp>0) toyNexp=_toyNexp;
  // get #evt option
  rarStrParser toyNevtParser=readConfStr("toyNevt", "0 fixed", _runSec);
  Double_t toyNevt=atoi(toyNevtParser[0]);
#ifndef USENEWROOT 
  if (toyNevt<=0) toyNevt=_protDataset->numEntries(kTRUE); // use weight
#else
  if (toyNevt<=0) toyNevt=_protDataset->numEntries(); // use weight
#endif
  toylist.setDataset(_protDataset->GetName(), toyNevt);

  Bool_t extendedGen=
    ("floated"==toyNevtParser[1])||
    ("extended"==toyNevtParser[1])||
    ("notfixed"==toyNevtParser[1]);
  cout<<"Toy Total Number of events: "<<toyNevt << " to be generated ";
  if (extendedGen) cout<<"with";
  else cout<<"without";
  cout<<" Poisson fluctuation"<<endl;
  // find all parameters of nEvt
  RooArgSet compCoeffSet("Comp Coeff Set");
  RooArgSet coeffParamSet("Comp Coeff Parameter Set");
  Int_t nComp=1;
  if (getControlBit("SimFit")) nComp=_nComp;
  for (Int_t i=0; i<nComp; i++) {
    rarAdd *compPdf=(rarAdd*)_pdfList.At(i);
    RooArgList thisCompCoeffList(compPdf->getCoeffList());

    compCoeffSet.add(thisCompCoeffList);
    cout << "Toy Size thisCompCoeffList " << thisCompCoeffList.getSize() << endl;
    thisCompCoeffList.Print("v");
    for (Int_t j=0; j<thisCompCoeffList.getSize(); j++) {

      // this line seems to work differently  in standalone - FFW
      RooArgList thisCompCoeffParams(*(thisCompCoeffList[j].getParameters(*_fullObs)));
      //cout << "Toy Size " << j << " " << thisCompCoeffParams.getSize() << endl;
 
      if (thisCompCoeffParams.getSize() == 0) {
	coeffParamSet.add(thisCompCoeffList[j]); // needed in standalone
      } else {
	for (Int_t k=0; k<thisCompCoeffParams.getSize(); k++) {
	  TString thisParamName=thisCompCoeffParams[k].GetName();
	  //cout << "Toy thisParamName " << thisParamName << " " << thisCompCoeffParams.getSize() << endl;
	  if (isFracName(thisParamName)) continue;
	  coeffParamSet.add(thisCompCoeffParams[k]);
	}
      }
    }
  }
  Int_t nCoeff=compCoeffSet.getSize();
  Int_t nCoeffParam=coeffParamSet.getSize();
  RooArgList compCoeffList(compCoeffSet);
  RooArgList coeffParamList(coeffParamSet);

  if (toyVerbose) {
    cout << "Toy: List of coefficients: " << nCoeff << endl;
    compCoeffList.Print("v");
    cout << "Toy: List of parameters: " << nCoeffParam << endl;
    coeffParamList.Print("v");
  }

  // find number of events for each type of data source
  RooArgSet evtTypeVarSet("Event Type Var Set");
  RooArgSet evtTypeStrSet("Event Type Str Set"); // for dynamic embed gen
  for (Int_t i=0; i<nCoeffParam; i++) {
    RooAbsReal *theCoeff=dynamic_cast<RooAbsReal*>(coeffParamList.at(i));
    if (!theCoeff) {
      cout<<endl<<"Toy W A R N I N G ! ! !"<<endl
	  <<" "<<coeffParamList[i].GetName()
	  <<" is not a (sub-)class of RooAbsReal."<<endl
	  <<" If it is a blind category var,"
	  <<" please set it to unblind for toy study."<<endl
	  <<" ie, by adding in section ["<<_runSec<<"]"<<endl
	  <<" "<<coeffParamList[i].GetName()<<" = 0"<<endl
	  <<endl;
      continue;
    }
    //
    TString coeffName = theCoeff->GetName();
    toylist.setInitial(coeffName, theCoeff->getVal());

    rarStrParser toySrcParser=
      readConfStr(Form("toySrc_%s", theCoeff->GetName()),
		  Form("pdf %f",theCoeff->getVal()), _runSec);
    Double_t coeffParamVal(0);
    Int_t nSrcType=toySrcParser.nArgs()/2; // type of sources for this Nevt
    for (Int_t j=0; j<nSrcType; j++) {
      TString evtSrc=toySrcParser[0];
      toySrcParser.Remove();
      TString evtSrcStr=toySrcParser[0];
      toySrcParser.Remove();
      Double_t evtSrcVal=getFormulaVal(evtSrcStr);
      coeffParamVal+=evtSrcVal;
      TString evtTypeName=
	Form("%s \"%s\"", coeffParamList[i].GetName(), evtSrc.Data());
      RooRealVar *evtTypeVar=(RooRealVar*)evtTypeVarSet.find(evtTypeName);
      toylist.setMethod(coeffName, evtSrc.Data());

      if (!evtTypeVar) {
	evtTypeVar=new RooRealVar(evtTypeName, evtTypeName, 0);
	evtTypeVarSet.addOwned(*evtTypeVar);
      }
      evtTypeVar->setVal(evtTypeVar->getVal()+evtSrcVal);
      // save the string for dynamic toy generating
      RooStringVar *evtTypeStr=(RooStringVar*)evtTypeStrSet.find(evtTypeName);
      if (!evtTypeStr) {
	evtTypeStr=new RooStringVar(evtTypeName, evtTypeName, "");
	evtTypeStrSet.addOwned(*evtTypeStr);
      }
      evtTypeStr->
	setVal(Form("%s \"%s\"",evtTypeStr->getVal(),evtSrcStr.Data()));
    }
    RooRealVar *coeffParamVar=(RooRealVar*)coeffParamList.at(i);
    // seg faults with RooFormula
    if (!coeffParamVar->inRange(coeffParamVal, "range")) {
      cout<<coeffParamVal<<" not within the range of "
	  <<coeffParamVar->GetName()<<": ("
	  <<coeffParamVar->getMin()<<","
	  <<coeffParamVar->getMax()<<")"<<endl;
      exit(-1);
    }
    coeffParamVar->setVal(coeffParamVal);
    //cout<<" "<<coeffParamVar->GetName()<<"="<<coeffParamVal<<endl;
    toylist.setRequested(coeffName, coeffParamVal);
  } // for (Int_t i=0; i<nCoeffParam; i++)

  if (toyVerbose) {
    cout << "Toy: List of generating methods: " << nCoeff << endl;
    evtTypeVarSet.Print("v");
    evtTypeStrSet.Print("v");
  }

  RooArgList evtTypeVarList(evtTypeVarSet);
  RooArgList evtTypeStrList(evtTypeStrSet);
  Int_t nEvtType=evtTypeVarList.getSize();
  { // now, let's figure out what need to be adjusted

    // methods for adjusting events: total/largest/prorata/smallest
    TString toyAdjustMethod=readConfStr("toyAdjustMethod", "Largest", _runSec);
    toyAdjustMethod.ToLower();
    if (toyAdjustMethod.Contains("largest") ||
	toyAdjustMethod.Contains("total") ||
	toyAdjustMethod.Contains("smallest") ||
	toyAdjustMethod.Contains("prorata")) {
      //cout << "toyAdjustMethod: using method " << toyAdjustMethod 
      //	   << " to adjust events in toy." << endl;
    } else {
      cout << "Toy toyAdjustMethod: expected token [Largest|Total|Smallest|Prorata]" 
	   << " found : " << toyAdjustMethod << "." << endl;
      //exit(-1);
    }

    // identify which dataset to adjutst if necessary. Use the dataset
    // that has been set to 0 events e.g. nSig = pdf 0 first; if not set
    // use largest.
    TString zEvtName("notSet"), lEvtName("notSet");
    Int_t zEvtIdx(-1),lEvtIdx(-1); //0 evt comp index, largest # evt comp index
    Double_t lEvt(0); // largest #s evts
    for (Int_t i=0; i<nEvtType; i++) {
      RooRealVar *evtTypeVar=(RooRealVar *)evtTypeVarList.at(i);
      TString evtTypeVarName=evtTypeVar->GetName();
      Double_t evtTypeVarVal=evtTypeVar->getVal();
      if(0==evtTypeVarVal) {
	if (zEvtIdx!=-1) {
	  cout<<"Toy Can not have more than one coeff=0"<<endl;
	  exit(-1);
	} else {
	  //zEvtIdx=i; // commented out to disable 0 event adjustment
	  zEvtName=evtTypeVarName;
	}
      }
      if(evtTypeVarVal>lEvt) {
	lEvt=evtTypeVarVal;
	lEvtIdx=i;
	lEvtName=evtTypeVarName;
      }
    }
    if (toyVerbose) {
      if (zEvtName != "notSet") {
	cout << "Toy Dataset with #evts = 0 specified:     " << zEvtName << endl;
      }
      cout << "Toy Dataset with largest #evts specified: " << lEvtName << endl;
    }
    // see if they can be adjusted
    Bool_t zParamAdj(kFALSE), lParamAdj(kFALSE);
    while (1) {
      if (-1==zEvtIdx) break;
      if (!compCoeffSet.find(rarStrParser(zEvtName)[0])) {
	if (!compCoeffSet.find("theSim_"+rarStrParser(zEvtName)[0])) {//simuPdf
	  // do not allow this
	  cout<<"Toy " <<zEvtName<<" is set to zero but is not adjustable"<<endl;
	  exit(-1);
	  break;
	}
      }
      zParamAdj=kTRUE;
      break;
    }
    while (1) {
      if (-1==lEvtIdx) break;
      if (!compCoeffSet.find(rarStrParser(lEvtName)[0])) {
	if (!compCoeffSet.find("theSim_"+rarStrParser(lEvtName)[0])) {
	  // warning about this
	  cout<<"Toy "<<lEvtName<<" is the largest but might not be adjustable"<<endl;
	}
      }
      lParamAdj=kTRUE;
      break;
    }
    // check total # evt
    Double_t totEvt(0);
    for (Int_t i=0; i<nCoeff; i++) { totEvt+=((RooAbsReal&)compCoeffList[i]).getVal(); }

    if (totEvt!=toyNevt) {
      Double_t nAdjust=toyNevt-totEvt;
      Int_t adjIdx(-1);
      // we prefer to adjust zEvtIdx
      if (zParamAdj) adjIdx=zEvtIdx;
      else if (lParamAdj) adjIdx=lEvtIdx;
      else {
	cout << "Toy #totEvts = " << totEvt << " #toyNevt = " << toyNevt << endl;
	cout<<"Need to adjust number of evenst but none of the toySrc_* can be adjusted" 
	    <<endl;
	exit(-1);
      }
      // adjust that evt type for the 'var', not for the 'str',
      // hopefully, the 'str' is a formula so the adj is done automatically
      RooRealVar *evtTypeVar=(RooRealVar *)evtTypeVarList.at(adjIdx);
      evtTypeVar->setVal(evtTypeVar->getVal()+nAdjust);
      // adjust that coeff
      RooRealVar *coeffAdjust=(RooRealVar*)
	compCoeffSet.find(rarStrParser(evtTypeVar->GetName())[0]);
      if (!coeffAdjust) coeffAdjust=(RooRealVar*)
	coeffParamSet.find(rarStrParser(evtTypeVar->GetName())[0]);
      if (!coeffAdjust) {
	cout<<"Toy Can not find var "<<rarStrParser(evtTypeVar->GetName())[0]
	    <<endl;
	exit(-1);
      }
      nAdjust+=coeffAdjust->getVal();
      cout<<"Toy "<<coeffAdjust->GetName()<<" needs to be adjusted from "
	  <<coeffAdjust->getVal()<<" to "<<nAdjust<<endl;
      if (!coeffAdjust->inRange(nAdjust, "range")) {
	cout<<"Toy "<<nAdjust<<" not within the range of "<<coeffAdjust->GetName()
	    <<": ("<<coeffAdjust->getMin()<<","
	    <<coeffAdjust->getMax()<<")"<<endl;
	exit(-1);
      }
      coeffAdjust->setVal(nAdjust);
      { // double-check if we have adjusted
	Double_t totEvt(0);
	for (Int_t i=0; i<nCoeff; i++) {totEvt+=((RooAbsReal&)compCoeffList[i]).getVal();}

	if (fabs(totEvt-toyNevt)>.01) {
	  cout<<endl<<"Toy Adjusting "<<evtTypeVar->GetName()
	      <<" will not set totEvtGen "<<totEvt<<" to "<<toyNevt<<endl
	      <<" If you are using blind var for yields,"<<endl
	      <<" please make sure you set it to unblind"<<endl
	      <<" (or temporarily do not use blind var for toy study)"<<endl
	      <<" (search for `W A R N I N G ! ! !\' in your log file)"<<endl;
	  exit(-1);
	}
      }
    }

    // update table with number to be used
    for (Int_t i=0; i<nCoeff; i++) {
      Double_t nevts = ((RooAbsReal&)compCoeffList[i]).getVal();
      TString coeffName = ((RooAbsReal&)compCoeffList[i]).getTitle();
      //cout << "coeffname " << coeffName << endl;
      toylist.setUsed(coeffName, nevts);
    }

    //if (toyVerbose) { evtTypeVarSet.Print("v"); }
  }
  toylist.print();

  // get data output file prefix
  TString toyFilePrefix=readConfStrCnA("toyDataFilePrefix", "default");
  if (("default"==toyFilePrefix)||
      ("no"==toyFilePrefix))
    toyFilePrefix="no";
  else {
    if ("yes"==toyFilePrefix) toyFilePrefix=Form("N \"%s\"", _runSec.Data());
    toyFilePrefix=getParamFileName("toySample", toyFilePrefix, _toyDir);
    toyFilePrefix+=Form(".%03d", _toyID);
    toyFilePrefix+=".%03d";
  }
  // generate and fit options
  TString genOpt="r"; // always randomize
  if (extendedGen) genOpt+="e";
  TString fitOpt=readConfStr("toyFitOption", "mhq", _runSec);
  // force option to be extended and return results
  fitOpt+="er";
  // convert to newer fitTo format for steering options
  //Bool_t toyFitExtended = fitOpt.Contains("e"); // must be true
  Bool_t toyFitMinos    = fitOpt.Contains("m");
  Bool_t toyFitHesse    = fitOpt.Contains("h");
  Bool_t toyFitVerbose  = !fitOpt.Contains("q");
  //Bool_t toyFitSave     = fitOpt.Contains("r"); // must be true

  // get number of cpus option
  Int_t toyFitNumCPU=atoi(readConfStr("useNumCPU", "1", getMasterSec()));

  // now see if we have toyFitMinos (do Minos only for some parameters)
  RooArgSet toyFitMinosAS;
  TString toyFitMinosStr=readConfStr("toyFitMinos", "notSet", _runSec);
  if ("notSet"!=toyFitMinosStr) {
    toyFitMinosAS.add(getArgSet(toyFitMinosStr,kFALSE,&fullParams));
  }
  
  // print full obs
  if (toyVerbose) {
    cout<<"Toy full obs for toy study"<<endl;
    _fObsSet.Print("v");
  }

  // get protGen
  if (!getProtGen()) {
    cout<<"Toy Can not create protGen!"<<endl;
    exit(-1);
  }
  
  // get comp-cat-ed datasets
  getCompCatDS(&_protDatasetsM, _protDataset);
  // for individual protDataset
  const char *dummy("");
  for (Int_t i=0; i<_nComp; i++) {
    rarMLPdf *model=(rarMLPdf*)_pdfList.At(i);
    TList *modelPdfList=model->getPdfList();
    Int_t nModelComp=modelPdfList->GetSize();
    for (Int_t j=0; j<nModelComp; j++) {
      rarBasePdf *comp=(rarBasePdf*)modelPdfList->At(j);
      RooDataSet *theData(0);
      if (protDataStrParser.nArgs()>0) {
	theData=_datasets->getData(protDataStrParser[0]);
	protDataStrParser.Remove();
      } else theData=comp->getData(dummy);
      _compCat.setIndex(j);
      getCompCatDS(&_protDatasets, theData, &_compCat);
    }
  }
  
  // determine protGenLevel
  _protGenLevel=atoi(readConfStr("protDataGenLevel", "1", _runSec));
  if (_protGenLevel>3) _protGenLevel=3;
  if (_protGenLevel<=0) _protGenLevel=0;
  // let's find fullProtVars
  RooArgSet fullProtVars(_protDataVars);
  { // let's find fullProtVars
    RooArgSet protGenVars(*_theProtGen->getDependents(_protDataset));
    fullProtVars.add(protGenVars);
    if (fullProtVars.getSize()<=0) _protGenLevel=0; // no need for protData
    else {
      if (0==_protGenLevel) _protGenLevel=2;
    }
    //_protDataVars.add(fullProtVars); // get full protVars//disable it for now
    if (toyVerbose) {
      if (fullProtVars.getSize() > 0) {
	cout<<"Toy Full protDataVars defined:"<<endl;
	fullProtVars.Print("v");
      } else {
	cout<<"Toy No Full protDataVars defined."<<endl;
      }
      cout<<"Toy protGenLevel: "<<_protGenLevel<<endl;
    }
  }
  
  if (!getGenerator()) {
    cout<<"Toy Can not create toy Gen!"<<endl;
    exit(-1);
  }
  // clone it to keep the initial values for fitting
  RooArgSet pdfNparamFSet(*_theGen, *_theProtGen);
  pdfNparamFSet.add(_subPdfs); // add original params, pdfs
  RooArgSet *fCloneSet=(RooArgSet*)pdfNparamFSet.snapshot(kTRUE);
  RooAbsPdf *theGen=(RooAbsPdf*)fCloneSet->find(_theGen->GetName());
  RooAbsPdf *theProtGen=(RooAbsPdf*)fCloneSet->find(_theProtGen->GetName());
  
  // define protData
  RooDataSet *protData(0);
  if (1==_protGenLevel) { // get default protData
    protData=(RooDataSet*)_protDataset->reduce(fullProtVars);
    if (toyVerbose) {
      cout<<"Toy The prototype dataset for toy study:"<<endl;
      protData->Print("v");
      //protData->write("/tmp/pd.text");
      cout<<endl;
    }
  }
  
  // save params just before toy study begins
  string fParamSStr;
  writeToStr(fullParams, fParamSStr);
  // try to figure out how many loops
  Int_t nLoops=1; // one loop
  Int_t nExpPerLoop=toyNexp;
  Bool_t tooManyEvts=kFALSE;
  if (nExpPerLoop*toyNevt>1000000) tooManyEvts=kTRUE;
  if (tooManyEvts||(_protGenLevel>1)||(preToyRandParser.nArgs()>0)) {
    nLoops=toyNexp;
    nExpPerLoop=1;
  }
  // loop to do toy study
  Bool_t firstToy(kTRUE);
  //for (Int_t expIdx=0; expIdx<toyNexp; expIdx++) {
  for (Int_t expIdx=0; expIdx<nLoops; expIdx++) {
    //if (expIdx>0) {RooMsgService::instance().setSilentMode(kTRUE); }
    // get toy data file
    TString toyFileName=Form(toyFilePrefix.Data(), expIdx);
    if (nExpPerLoop>1) toyFileName+=".%03d";
    toyFileName+=".text";
    // restore param for each toy loop
    readFromStr(fullParams, fParamSStr);
    readFromStr(*fCloneSet, fParamSStr);
    // do we need to randomize params?
    if ((preToyRandParser.nArgs()>0)&&_theToyParamGen) { // yes
      // generate an set of randomizers
      RooDataSet *theRandDataSet=_theToyParamGen->generate(paramRandSet, 1);
      RooArgSet *theRandSet=(RooArgSet*)theRandDataSet->get(0);
      stringstream theRandSetStr;
      theRandSet->writeToStream((ostream&)theRandSetStr, kFALSE);
      delete theRandDataSet;
      paramRandSet.readFromStream((istream&)theRandSetStr, kFALSE);
      cout<<"Toy Randomizer ArgSet:"<<endl;
      paramRandSet.Print("v");
      Int_t nParams=preToyRandParser.nArgs()/2; // # of params
      RooArgSet rParams, cParams;
      for (Int_t i=0; i<nParams; i++) {
	RooRealVar *theParam=dynamic_cast<RooRealVar*>
	  (getAbsVar(preToyRandParser[2*i]));
	if (!theParam) {
	  cout<<"Toy Can not find param "<<preToyRandParser[2*i]<<endl;
	  continue;
	}
	rParams.add(*theParam);
	// find cloned param
	theParam=(RooRealVar*)fCloneSet->find(theParam->GetName());
	if (!theParam) {
	  cout<<"Toy Can not find cloned param "<<theParam->GetName()<<endl;
	  exit(-1);
	}
	cParams.add(*theParam);
	cout<<"Toy "<<theParam->GetName()<<" = "<<theParam->getVal();
	theParam->setVal(getFormulaVal(preToyRandParser[2*i+1]));
	cout<<" --> "<<theParam->getVal()<<endl;
      }
      // now save cParams to rParams
      string cTorStr;
      writeToStr(cParams, cTorStr);
      readFromStr(rParams, cTorStr);
    }
    // save params after possible randomization
    string randParamSStr;
    writeToStr(fullParams, randParamSStr);
    // now check if there are any data need to be generated from dataset
    TString genStr="";
    Double_t toyNevtGen=toyNevt;
    for (Int_t i=0; i<nEvtType; i++) {
      RooStringVar *evtTypeStr=(RooStringVar *)evtTypeStrList.at(i);
      rarStrParser toySrcParser=evtTypeStr->GetName();
      rarStrParser toySrcStrParser=evtTypeStr->getVal();
      if ("pdf"==toySrcParser[1]) continue; // otherwise, from other src
      // get evtTypeStr's value
      Double_t evtTypeStrVal(0);
      while (toySrcStrParser.nArgs()>0) {
	evtTypeStrVal+=getFormulaVal(toySrcStrParser[0]);
	toySrcStrParser.Remove();
      }
      // get #evt for this src
      Double_t nEvtWSrc(0);
      for (Int_t j=0; j<nCoeff; j++) {
	TString coefName=compCoeffList[j].GetName();
	nEvtWSrc+=((RooAbsReal*)fCloneSet->find(coefName))->getVal();
      }
      RooRealVar *paramVar=(RooRealVar*)fCloneSet->find(toySrcParser[0]);
      paramVar->setVal(paramVar->getVal()-evtTypeStrVal);
      Double_t nEvtWoSrc(0);
      for (Int_t j=0; j<nCoeff; j++) {
	TString coefName=compCoeffList[j].GetName();
	nEvtWoSrc+=((RooAbsReal*)fCloneSet->find(coefName))->getVal();
      }
      Double_t nEvtGen=nEvtWSrc-nEvtWoSrc;
      if (nEvtGen<0) {
	cout<<"Toy W A R N I N G !"<<endl
            <<" Can not generate negative event from src "<<toySrcParser[1]
	    <<" for parameter "<<toySrcParser[0]<<endl;
	continue;
      }
      genStr+=Form(" \"%s\" %f ", toySrcParser[1].Data(), nEvtGen);
      toyNevtGen-=nEvtGen;
    }
    if (""!=genStr) {
      cout<<"Toy Using generation string: \""<<genStr<<"\" for embedding"<<endl;
    }
    // do we need to generate protData?
    if (_protGenLevel>=2) {
      RooArgSet protDeps(fullProtVars);
      // add compCat as well
      protDeps.add(_compCat);
      // get expected number of event
      Int_t nEvt=_protDataset->numEntries();
      if (protData) delete protData; // delete old one
      if (_protGenLevel>=3) protData=theProtGen->generate(protDeps, nEvt);
      else { // from individual protDatasets
	// first generate cats
	RooArgSet fCatSet(*((RooSimultaneous*)theProtGen)
			 ->indexCat().getDependents(_protDataset));
	fCatSet.add(_compCat);
	RooDataSet *fCatData=theProtGen->generate(fCatSet, nEvt);
	// cat w/o protEVars
	RooArgSet catSet(fCatSet);
	// remove protEVars
	catSet.remove(_protDataEVars, kFALSE, kTRUE);
	RooSuperCategory sCat("sCat", "sCat", catSet);
	sCat.attachDataSet(*fCatData);
	protData=new RooDataSet("protData", "protData", protDeps);
	for (Int_t i=0; i<nEvt; i++) {
	  const RooArgSet *theCat=fCatData->get(i);
	  TString dsName=sCat.getLabel();
	  dsName.ReplaceAll("{", "");
	  dsName.ReplaceAll("}", "");
	  // find the dataset
	  RooDataSet *theData=(RooDataSet*)_protDatasets.FindObject(dsName);
	  if (!theData) {
	    if (dsName.Last(';')<0) theData=(RooDataSet*)_protDatasetsM.At(0);
	    else {
	      dsName.Remove(dsName.Last(';'), dsName.Length());
	      theData=(RooDataSet*)_protDatasetsM.FindObject(dsName);
	    }
	  }
	  if (!theData) theData=(RooDataSet*)_protDatasetsM.At(0);
	  Int_t I=RooRandom::randomGenerator()->Integer(theData->numEntries());
	  RooArgSet protSet(*theCat);
	  protSet.add(*theData->get(I), kTRUE);
	  protData->add(protSet);
	}
	delete fCatData;
      }
      cout<<"Toy protData from protGen:"<<endl;
      protData->Print("v");
      //protData->write("/tmp/pd.text");
      //TFile f("/tmp/pd.root", "recreate");
      //protData->tree().Write();
      //f.Close();
      RooArgList protDepList(protDeps);
      for (Int_t i=0; i<protDepList.getSize(); i++) {
	RooCategory *theCat=(RooCategory *)&(protDepList[i]);
	if (theCat->ClassName()!=TString("RooCategory")) continue;
	Roo1DTable *theTable(0);
	if (theCat&&(theTable=protData->table(*theCat))) {
	  theTable->Print();
	  // get frac table
	  TIterator* catTypeIter = theCat->typeIterator();
	  RooCatType *theType(0);
	  while(theType=(RooCatType*)catTypeIter->Next()) {
	    TString typeName=theType->GetName();
	    cout<<"    "<<typeName<<"\t"<<theTable->getFrac(typeName)<<endl;
	  }
	  cout<<endl;
	  delete catTypeIter;
	}
      }
    }
    // toy dependents
    RooArgSet toyDeps(*theGen->getDependents(_protDataset));
    if (protData) toyDeps.remove(*theGen->getDependents(protData));
    //protData->Print("v");
    toyDeps.Print();
    RooMCStudy *theToy(0);
    if (toyFitMinosAS.getSize()<1) {
      theToy=new RooMCStudy
	(*theGen, toyDeps, FitModel(*_thePdf),
	 Extended(extendedGen), ProtoData(*protData, kTRUE),
	 ProjectedObservables(_conditionalObs),
	 FitOptions(Save(kTRUE), Extended(kTRUE),
		    Verbose(toyFitVerbose), Hesse(toyFitHesse),
		    Minos(toyFitMinos), NumCPU(toyFitNumCPU)));
    } else {
      theToy=new RooMCStudy
	(*theGen, toyDeps, FitModel(*_thePdf),
	 Extended(extendedGen), ProtoData(*protData, kTRUE),
	 ProjectedObservables(_conditionalObs),
	 FitOptions(Save(kTRUE), Extended(kTRUE),
		    Verbose(toyFitVerbose), Hesse(toyFitHesse),
		    Minos(toyFitMinos), Minos(toyFitMinosAS)), FitOptions(NumCPU(toyFitNumCPU)));
    }
    if (firstToy) {
      cout<<"Toy The toyDeps"<<endl;
      toyDeps.Print("v");
      cout<<"Toy The protVars"<<endl;
      _protDataVars.Print("v");
    }
    // do we want to check negative pdf value
    Bool_t toyChkNegativePdf(kFALSE);
    if ("yes"==readConfStr("toyChkNegativePdf", "no", _runSec)) {
      toyChkNegativePdf=kTRUE;
    }
    RooDataSet *chkNPdfDS(0);
    // generate
    if ("no"!=readConfStr("toyGenerate", "yes", _runSec)) {
      if (firstToy) {
	cout<<endl<<"Toy RooMCStudy params for generating:"<<endl;
	theGen->getParameters(_protDataset)->Print("v");
      }
      // first generate dataset to check negative pdf value
      if (toyChkNegativePdf) {
        chkNPdfDS=_dummyPdf.generate(toyDeps, *protData, (Int_t)toyNevt);
      }
      // generate through standard RooMCStudy
      theToy->generate(nExpPerLoop, randInt(toyNevtGen), kTRUE);
      RooArgSet toyObs(toyDeps);
      if (protData) toyObs.add(*protData->get(),kTRUE);
      if (firstToy) {
	cout<<"Toy Full obs generated:"<<endl;
	toyObs.Print();
      }
      if (kTRUE) { // now we need to find if we need two-step generation
        //generate(theToy, theGen, genOpt, nExpPerLoop);
      }
      if (""!=genStr) {
	// find obs to get distribution from prototype dataset
	RooArgList etoyDeps(toyObs);
	etoyDeps.remove(_protDataEVars, kFALSE, kTRUE);
	if (firstToy) {
	  cout<<"Toy etoyDeps"<<endl;
	  etoyDeps.Print("v");
	  cout<<"Toy protDataEVars"<<endl;
	  _protDataEVars.Print("v");
	  cout<<"Toy genStr "<<genStr<<endl;
	}
	generate(theToy, etoyDeps, genStr, genOpt, nExpPerLoop);
      }
      if (!toyFileName.BeginsWith("no")) { // output toy data per request
        Int_t nSamples=nExpPerLoop;
        while (nSamples--) {
	  TString thisToyFileName=Form(toyFileName.Data(),nSamples);
          ((RooDataSet*)theToy->genData(nSamples))->write(thisToyFileName);
          //cout<<"toy sample structure"<<endl;
          //((RooDataSet*)theToy->genData(nSamples))->Print();
	  thisToyFileName.ReplaceAll(".text", ".root");
#ifndef USENEWROOT
	  TFile f(thisToyFileName, "recreate");
	  ((RooDataSet*)theToy->genData(nSamples))->tree().Write();
	  f.Close();
#else
	  RooDataSet *theSet = (RooDataSet*) theToy->genData(nSamples);
	  saveAsRootFile(theSet, thisToyFileName, kTRUE);
#endif
	}
      }
      if (firstToy) { // save sample 
	RooDataSet *toySample=(RooDataSet*)theToy->genData(0);
	toySample=(RooDataSet*)toySample->Clone();
	toySample->SetName(_datasets->getDSName("toySample"));
	_datasets->getDatasetList()->Add(toySample);
	_datasets->ubStr(_datasets->getDSName("toySample"), "Unblinded");
      }
    }
    // fit
    if ("no"!=readConfStr("toyFit", "yes", _runSec)) {
      if (firstToy) {
	cout<<endl<<"Toy fitting coeffParamSet:"<<endl;
	coeffParamSet.Print("v");
      }
      // check negative pdf values
      if (toyChkNegativePdf&&chkNPdfDS) {
        // first attach the dataset
        attachDataSet(*chkNPdfDS);
        Int_t nChkEvt=chkNPdfDS->numEntries();
        for (Int_t iChkEvt=0; iChkEvt<nChkEvt; iChkEvt++) {
          RooArgSet *theEvt=(RooArgSet*)chkNPdfDS->get(iChkEvt);
          if (isNegativeValue()) {
            cout<<"Toy The PDF have negative value for current event"<<endl;
            theEvt->Print();
            theEvt->Print("v");
            exit(-1);
          }
        }
      }
      if ("no"!=readConfStr("toyGenerate", "yes", _runSec)) {
        // construct dataset list
        TList genSamples;
        Int_t nSamples=nExpPerLoop;
        while (nSamples--) genSamples.Add(theToy->genData(nSamples)->Clone());
        theToy->fit(nExpPerLoop, genSamples);
      } else {
        theToy->fit(nExpPerLoop, toyFileName);
      }
      // pulls for embedded toy are not right, re-calculate
      if (""!=genStr) {
	cout<<"Toy Re-calculate pulls for embedded toys"<<endl;
	// restore the randomized params
	readFromStr(fullParams, randParamSStr);
	RooDataSet &fitParData=(RooDataSet &)theToy->fitParDataSet();
	const RooArgSet *fitParDataSet=fitParData.get();
	TIterator* iter = coeffParamSet.createIterator();
	RooRealVar *par(0);
	while(par=(RooRealVar*)iter->Next()) {
	  TString pullName=Form("%spull", par->GetName());
	  RooPullVar *origPull=(RooPullVar*)fitParDataSet->find(pullName);
	  if (!origPull) continue;
	  RooAbsReal* tPar=(RooAbsReal*)par->Clone("truth");
	  RooPullVar pull(pullName+"_embd", origPull->GetTitle(),*par,*tPar);
	  fitParData.addColumn(pull);
	  delete tPar;
	}
	delete iter;
	// add vars for # of embedded events
	rarStrParser genStrParser=genStr;
	Int_t nSubset=genStrParser.nArgs()/2;
	for(Int_t i=0; i<nSubset; i++) {
	  TString genSrcName=genStrParser[0];
	  genStrParser.Remove();
	  Double_t nEvtGen=atof(genStrParser[0]);
	  genStrParser.Remove();
	  RooRealVar embdEvt("embdEvt_"+genSrcName, "embd Evt", nEvtGen);
	  fitParData.addColumn(embdEvt);
	}
      }
      { // add toy ID's
	RooDataSet &fitParData=(RooDataSet &)theToy->fitParDataSet();
	RooRealVar ID_COL0("ID_COL0", "Toy ID", _toyID);
	fitParData.addColumn(ID_COL0);
	RooRealVar ID_COL1("ID_COL1", "Toy Loop ID", expIdx);
	fitParData.addColumn(ID_COL1);
	RooRealVar ID_COL2("ID_COL2", "Toy Loop #", 0);
	fitParData.addColumn(ID_COL2);
	// add correlation matrix for floating params
	const RooFitResult *fr=theToy->fitResult(0);
	Int_t nFloat=fr->floatParsFinal().getSize();
	for (Int_t iF=1; iF<=nFloat; iF++) {
	  RooRealVar CorrG(Form("CorrG_%i",iF), "Global Corr",0);
	  fitParData.addColumn(CorrG);
	  for(Int_t jF=1; jF<iF; jF++) {
	    RooRealVar Corr(Form("Corr_%i_%i", iF, jF), "Corr",0);
	    fitParData.addColumn(Corr);
	  }
	}
	// add more result info
	RooRealVar covQual("covQual", "covQual", 0);
	fitParData.addColumn(covQual);
	RooRealVar numInvalidNLL("numInvalidNLL", "numInvalidNLL", 0);
	fitParData.addColumn(numInvalidNLL);
	RooRealVar edm("edm", "edm", 0);
	fitParData.addColumn(edm);
	
        // add gof chisq
        RooRealVar GOFChisq("GOFChisq", "GOFChisq", 0);
        fitParData.addColumn(GOFChisq);
      }
      if (!toyResults) {
	toyResults=new
          RooDataSet("toyResults","toyResults",*theToy->fitParDataSet().get());
      }
      // merge w/ toy IDs and check if all toy fits converged
      Int_t ii=0;
      for (Int_t i=0; i<nExpPerLoop; i++) {
	RooFitResult *fr=(RooFitResult*)theToy->fitResult(i);
	if (fr->status()) {
	  cout<<"Toy Fit status for experiment #"<<expIdx<<"-"<<nExpPerLoop-i
	      <<": "<<fr->status()<<endl;
	  continue;
	}
	RooArgSet *fitResultSet=(RooArgSet *)theToy->fitParams(ii);
	((RooRealVar*)fitResultSet->find("ID_COL2"))->setVal(nExpPerLoop-i);
	// add correlation matrix
	Int_t nFloat=fr->floatParsFinal().getSize();
	for (Int_t iF=0; iF<nFloat; iF++) {
	  RooAbsArg &iArg=fr->floatParsFinal().operator[](iF);
	  ((RooRealVar*)fitResultSet->find(Form("CorrG_%i",iF+1)))->
	    setVal(fr->globalCorr(iArg));
	   for(Int_t jF=0; jF<iF; jF++) {
	     RooAbsArg &jArg=fr->floatParsFinal().operator[](jF);
	     ((RooRealVar*)fitResultSet->find(Form("Corr_%i_%i", iF+1,jF+1)))->
	       setVal(fr->correlation(iArg, jArg));
	   }
	}
	// add more result info
	((RooRealVar*)fitResultSet->find("covQual"))->setVal(fr->covQual());
	((RooRealVar*)fitResultSet->find("numInvalidNLL"))
	  ->setVal(fr->numInvalidNLL());
	((RooRealVar*)fitResultSet->find("edm"))->setVal(fr->edm());
	
        // gof chisq
        Double_t gofChisq(0);
	TString postMLGOFChisq=readConfStrCnA("postMLGOFChisq", "no");
        if (!postMLGOFChisq.BeginsWith("no"))
          gofChisq=doGOFChisq((RooDataSet*)theToy->genData(i), cout);
        ((RooRealVar*)fitResultSet->find("GOFChisq"))->setVal(gofChisq);

	toyResults->add(*fitResultSet);
	ii++;
      }
    }
    delete theToy;
    firstToy=kFALSE;
  }
  
  //return theToy;
  return toyResults;
}

/// \brief Get comp-cat-ed datasets
/// \param ds List of comp-cat-ed datasets
/// \param iData Dataset to be comp-cat-ed
/// \param compCat Optional component cat
///
/// It splits the input dataset into comp-cat-ed datasets
void rarMLFitter::getCompCatDS(TList*ds,RooDataSet *iData,RooCategory *compCat)
{
  RooArgSet catSet(*((RooSimultaneous*)_theProtGen)
		   ->indexCat().getDependents(iData));
  catSet.remove(_compCat, kFALSE, kTRUE);
  catSet.remove(_protDataEVars, kFALSE, kTRUE);
  if (catSet.getSize()<=0) {
    RooDataSet *theData=(RooDataSet*)iData->reduce("1");
    theData->SetName("noCompCatDS");
    if (compCat) {
      theData->addColumn(*compCat);
      theData->SetName(compCat->getLabel());
    }
    ds->Add(theData);
    return;
  }
  // we do have split cats, let's create datasets for each cat-type
  RooSuperCategory sCat("sCat", "sCat", catSet);
  sCat.attachDataSet(*iData);
  Int_t nEvt=iData->numEntries();
  for (Int_t i=0; i<nEvt; i++) {
    RooArgSet *theEvt=(RooArgSet*)iData->get(i);
    TString dsName=sCat.getLabel();
    if (compCat) dsName=Form("%s;%s", dsName.Data(), compCat->getLabel());
    dsName.ReplaceAll("{", "");
    dsName.ReplaceAll("}", "");
    RooDataSet *theData=(RooDataSet*)ds->FindObject(dsName);
    if (!theData) {
      theData=new RooDataSet(dsName, dsName, *theEvt);
      if (compCat) theData->addColumn(*compCat);
      ds->Add(theData);
    }
    theData->add(*theEvt);
  }
}

/// \brief Generate a random integer around a real number
/// \param iNumber Input real number;
/// \return Generated integer close
///
/// It generates a integer around input real number (within +/- .5)
/// so that the average of the generated integers will be the input number
Int_t rarMLFitter::randInt(Double_t iNumber)
{
  Int_t retVal=(Int_t) iNumber;
  if (RooRandom::randomGenerator()->Uniform()+retVal<iNumber) retVal++;
  return retVal;
}

/// \brief Two-step generator
void rarMLFitter::generate(RooMCStudy *theToy, RooAbsPdf *theGen,
			   const TString genOpt, const Int_t toyNexp)
{
  cout<<endl<<" In rarMLFitter generate (two-step) for "<<GetName()<<endl;
  
  // generate toys
  Int_t nSamples=toyNexp;
  while(nSamples--) {
    cout<<endl;
    RooDataSet *genSample=(RooDataSet *)theToy->genData(nSamples);
    //genSample->Print();
    //genSample->Print("v");
    // remove d-cosheli
    RooArgSet subSet(*genSample->get());
    subSet.remove(*subSet.selectByName("hzCat,hcCat"),
                  kTRUE, kTRUE);
    // reduce the dataset
    RooDataSet *subData=(RooDataSet*)genSample->reduce(subSet);
    // now add dcoshel0 and dcoshelC
    RooArgSet *addOns=getAddOnCols();
    addOns=(RooArgSet*)(addOns->selectByName("hzCat,hcCat"));
    subData->addColumns(*addOns);
    //subData->Print();
    //subData->Print("v");
    genSample->reset();
    genSample->append(*subData);
    //genSample->Print("v");
    delete subData;
  }
  if (genOpt.Contains("safasse")) {theGen->Print();} // avoid "unused" warning
}

/// \brief Generate toy sample from datasets
///
/// \param theToy RooMCStudy object for toy study
/// \param etoyDeps Observables to be embedded from dataset
/// \param genStr Generation string
/// \param genOpt Generation option
/// \param toyNexp Number of toy experiments
///
/// It generates all the events need to be generated from data sources
/// and appends those events to the data generated by PDFs.
/// #doToyStudy forms the generation string and it passes
/// the string to #generate to get this part of the toy
/// sample generated.
void rarMLFitter::generate(RooMCStudy *theToy,
			   const RooArgList& etoyDeps,
			   const TString genStr, const TString genOpt,
			   const Int_t toyNexp)
{
  cout<<endl<<" In rarMLFitter generate (embed) for "<<GetName()<<endl;
  
  // generate toys
  Int_t nSamples=toyNexp;
  while(nSamples--) {
    cout<<endl;
    RooDataSet *genSample=(RooDataSet *)theToy->genData(nSamples);
    // generate each sub sample from genStr
    rarStrParser genStrParser=genStr;
    Int_t nSubset=genStrParser.nArgs()/2;
    for(Int_t i=0; i<nSubset; i++) {
      TString genSrcName=genStrParser[0];
      genStrParser.Remove();
      Double_t nEvtGen=atof(genStrParser[0]);
      genStrParser.Remove();
      RooDataSet *embedSample=generate(etoyDeps, genSrcName, nEvtGen, genOpt);
      if (embedSample) {
	genSample->append(*embedSample);
	delete embedSample;
      }
    }
  }
}

/// \brief Generate toy dataset from source datasets
///
/// \param dependents Observables to be generated from dataset
/// \param genSrcName Data src name
/// \param nEvtGen Evt to be embedded
/// \param genOpt Generation option
/// \return Generated toy dataset
///
/// \p genStr has form like
/// \li <tt>\<srcDataName1\> \<nEvt1\> [\<srcDataName2\> \<nEvt2\> ...]</tt>
///
/// It loops through all such data sources and retrieve specified
/// number of events for each of them.
/// If \p genOpt contains option \p e, which means extended generation,
/// the number of events will be Poisson fluctuated.
RooDataSet *rarMLFitter::generate(const RooArgList& dependents,
				  const TString genSrcName, Double_t nEvtGen,
				  const TString genOpt)
{
  //cout<<" In rarMLFitter generate (embed) 2 for "<<GetName()<<endl;
  RooDataSet *theSample(0);
  RooArgSet fullDeps(dependents);
  RooArgSet randObs;
  RooAbsPdf *theObsGen(0);
  // first find out if need to randomize obs
  RooStringVar *srcRandStr=(RooStringVar *)_embdObsRandSet.find(genSrcName);
  if (srcRandStr) {
    rarStrParser srcRandStrParser=srcRandStr->getVal();
    //cout<<srcRandStr->getVal()<<endl;
    while (srcRandStrParser.nArgs()>0) {
      TString obsName=srcRandStrParser[0];
      srcRandStrParser.Remove();
      if (_fObsSet.find(obsName)) randObs.add(*_fObsSet.find(obsName));
      else theObsGen=(RooAbsPdf*)_embdObsGens.find("the_"+obsName);
    }
    if (theObsGen) {
      cout<<" Using "<<theObsGen->GetName()<<" as embd randomizer"<<endl;
      theObsGen->Print();
      //RooArgSet(*theObsGen).snapshot(kTRUE)->Print("v");
    } else {
      cout<<" Can not find randomizer defined in \""<<srcRandStr->getVal()
	  <<"\""<<endl;
    }
    if (randObs.getSize()>0) {
      cout<<" The following params will be randomized:"<<endl;
      randObs.Print("v");
    } else {
      cout<<" Can not find params to randomize in \""<<srcRandStr->getVal()
	  <<"\""<<endl;
    }
  }
  if (theObsGen) fullDeps.add(*theObsGen->getDependents(_protDataset), kTRUE);
  // remove compCat from fullDeps
  fullDeps.remove(_compCat, kFALSE, kTRUE);
  // get source dataset
  RooDataSet *genSrcData=_datasets->getData(genSrcName);
  if (!genSrcData) {
    cout<<"toy src data \""<<genSrcName<<"\" does not exist"<<endl;
    exit(-1);
  }

  if (genSrcData->numEntries()<=0) {
    cout<<"toy src data \""<<genSrcName<<"\" does not contain any entries"
	<<endl;
    exit(-1);
  }
  
  if (genOpt.Contains("e")) { // extended
    cout<<"Extended generating nEvt "<<nEvtGen;
    nEvtGen = RooRandom::randomGenerator()->Poisson(nEvtGen);
    cout<<" -> "<<nEvtGen<<endl;
  }
  nEvtGen=randInt(nEvtGen);
  cout<<"Generate "<<nEvtGen<<" events from dataset "
      <<genSrcData->GetName()<<endl;
  if (nEvtGen<.1) {
    return theSample;
  }
  // generate each uncorrelated sub sample
  rarStrParser embdUnCorrParser=readConfStr("toyEmbdUnCorrelate","",_runSec);
  for (Int_t i=0; i<=embdUnCorrParser.nArgs(); i++) {
    RooDataSet *subSample(0);
    if ((i==embdUnCorrParser.nArgs())&&(fullDeps.getSize()>0)) {
      //cout<<" Remaining obs:"<<endl;
      //fullDeps.Print();
      subSample=new RooDataSet("subSample", "subSample", fullDeps);
    } else {
      RooArgSet subSampleSet;
      rarStrParser groupObsParser=embdUnCorrParser[i];
      for (Int_t j=0; j<groupObsParser.nArgs(); j++) {
	RooAbsArg *theObs=fullDeps.find(groupObsParser[j]);
	if (theObs) {
	  subSampleSet.add(*theObs);
	  fullDeps.remove(*theObs);
	}
      }
      if (subSampleSet.getSize()<1) continue;
      cout<<" Sub obs:"<<endl;
      subSampleSet.Print();
      subSample=new RooDataSet("subSample", "subSample", subSampleSet);
    }
    if(!subSample) continue;
    for (Int_t j=0; j<nEvtGen; j++) {
      Int_t I=RooRandom::randomGenerator()->Integer(genSrcData->numEntries());
      const RooArgSet* row=genSrcData->get(I);
      if(!row) {
	cout<<"skipping null row entry I="<<I<<endl;
	continue;
      }
      subSample->add(*row);
    }
    if (!theSample)
      theSample=subSample;
    else {
      theSample->merge(subSample);
      delete subSample;
    }
  }
  
  // merge those columns from prototype dataset
  if (_protDataEVars.getSize()>0) {
    RooDataSet *protSample=
      new RooDataSet("protEDataSet", "prot embedded Dataset", _protDataEVars);
    Int_t nEvt=theSample->numEntries();
    for (Int_t i=0; i<nEvt; i++) {
      Int_t Idx=RooRandom::randomGenerator()->
	Integer(_protDataset->numEntries());
      const RooArgSet* row=_protDataset->get(Idx);
      protSample->add(*row);
    }
    cout<<"Merge "<<nEvt<<" events from prototype dataset "
	<<_protDataset->GetName()<<endl;
    theSample->merge(protSample);
    delete protSample;
  }
  { // merge with compCat set to max
    _compCat.setIndex(_compCat.numTypes()-1);
    theSample->addColumn(_compCat);
  }
  // do we need to randomize obs
  if (theObsGen&&randObs.getSize()>0) {
    RooArgSet embProtSet;
    embProtSet.add(*theSample->get(0));
    embProtSet.remove(randObs, kFALSE, kTRUE);
    RooDataSet *protDS=(RooDataSet*)theSample->reduce(embProtSet);
    RooDataSet *oldembedSample=theSample;
    theSample=theObsGen->generate(randObs,*protDS, (Int_t)nEvtGen);
    //for (Int_t iEvt=0; iEvt<nEvtGen; iEvt++) {
    //  oldembedSample->get(iEvt)->Print("v");
    //  theSample->get(iEvt)->Print("v");
    //}
    //theSample->write("/tmp/new.txt");
    //oldembedSample->write("/tmp/old.txt");
    delete oldembedSample;
    delete protDS;
  }
  
  return theSample;
}

/// \brief Return prototype var generator for toy study
/// \return The generator
///
/// It constrcuts prototype var generator for toy study
RooAbsPdf *rarMLFitter::getProtGen()
{
  if (_theProtGen) return _theProtGen;
  // call super class' getProtGen
  rarCompBase::getProtGen();
  // now tricky part, different treatment for simFit or not
  if (getControlBit("noSimFit")) { // no simFit
    _theProtGen=(RooAbsPdf*)_subProtGenPdfs.at(0);
  } else { // now, create simProtGen
    RooSimPdfBuilder *simBuilder=new RooSimPdfBuilder(_subProtGenPdfs);
    RooArgSet *simConfig=simBuilder->createProtoBuildConfig();
    // get configs from #_simConfig
    ((RooStringVar*)simConfig->find("splitCats"))->
      setVal(((RooStringVar*)_simConfig->find("splitCats"))->getVal());
    TString modelStr=((RooStringVar*)_simConfig->find("physModels"))->getVal();
    for (Int_t i=0; i<_subProtGenPdfs.getSize(); i++) {
      TString genName = _subProtGenPdfs[i].GetName();
      TString pdfName =_subPdfs[i].GetName();
      // construct splitting configStr
      TString configStr=((RooStringVar*)_simConfig->find(pdfName))->getVal();
      ((RooStringVar*)simConfig->find(genName))->setVal(configStr);
      modelStr.ReplaceAll(pdfName, genName);
    }
    ((RooStringVar*)simConfig->find("physModels"))->setVal(modelStr);
    cout<<"simPdfBuilder configs for protGen:"<<endl;
    simConfig->Print("v");
    simBuilder->addSpecializations(_simBuilder->splitLeafList());
    _theProtGen=(RooAbsPdf*)
      simBuilder->buildPdf(*simConfig, _protDataset, 0, kTRUE);
    // now for leftovers
    RooAbsCategoryLValue *simCats=(RooAbsCategoryLValue *)
      (&((RooSimultaneous*)_theProtGen)->indexCat());
    Int_t nCats=simCats->numTypes();
    for(Int_t i=0; i<nCats; i++) {
      TString catName=((RooCatType*)simCats->lookupType(i))->GetName();
      if (!((RooSimultaneous*)_theProtGen)->getPdf(catName.Data())) {
	((RooSimultaneous*)_theProtGen)->addPdf(_dummyExtPdf, catName);
      }
    }
    
    _theProtGen->SetName(Form("simProtGen_%s", GetName()));
    _theProtGen->Print();
    //_theProtGen->Print("v");
  }
  // now component-categorized
  _theProtGen=compGen(_theProtGen, _subProtGenPdfs, _compCat);
  
  return _theProtGen;
}

/// \brief Return generator for toy study
/// \return The generator
///
/// It constrcuts toy event generator
RooAbsPdf *rarMLFitter::getGenerator()
{
  if (_theGen) return _theGen;
  _theGen=_thePdf;
  if (getControlBit("SimFit")) { // get fitter cats
    RooAbsCategoryLValue *simCats=(RooAbsCategoryLValue *)
      (&((RooSimultaneous*)_thePdf)->indexCat());
    Int_t nCats=simCats->numTypes();
    // construct toy generator from fitter
    _theGen=new RooSimultaneous
      (Form("simGen_%s",GetName()),"simGen",*simCats);
    for(Int_t i=0; i<nCats; i++) {
      TString catName=((RooCatType*)simCats->lookupType(i))->GetName();
      RooAbsPdf *catPdf=((RooSimultaneous*)_thePdf)->getPdf(catName.Data());
      if (!catPdf) catPdf=&_dummyExtPdf;
      ((RooSimultaneous*)_theGen)->addPdf(*catPdf, catName);
    }
  }
  // now component-categorized
  if (_protGenLevel>1) _theGen=compGen(_theGen, _subPdfs, _compCat);
  
  cout<<"Generator created from fitter:"<<endl;
  _theGen->Print();
  return _theGen;
}

/// \brief Return the simPdf component with the name and cat label
/// \param simSet ArgSet of simPdf components
/// \param argName Component name
/// \param catName Cat label name
/// \param srcPdf Source Pdf
/// \return The pointer of the component found
///
/// It returns the component wrt the name and cat label
RooAbsArg *rarMLFitter::findSimed(RooArgSet &simSet,
				  TString argName, TString catName,
				  RooAbsPdf *srcPdf)
{
  RooAbsArg *theComp(0);
  // the default comb
  theComp=simSet.find(argName+"_"+catName);
  if (theComp) return theComp;
  // maybe we can put the 1st cat to the end
  if (catName.Contains(";")) {
    // find the first ;
    TString phyCat=catName(1, catName.First(";")-1);
    catName.Replace(1, phyCat.Length()+1,"");
    catName.Replace(catName.Length()-1, 0, ";"+phyCat);
    theComp=simSet.find(argName+"_"+catName);
    if (theComp) return theComp;
  }
  // now for the name itself
  theComp=simSet.find(argName);
  if (""==getPhysCat()||" "==getPhysCat()) return theComp; // no physCat
  if (!srcPdf) return theComp;
  RooSimultaneous *sPdf=dynamic_cast<RooSimultaneous*>(srcPdf);
  if (!sPdf) return theComp;
  // now let's find the arg in sPdf, set to default first
  theComp=0;
  RooAbsPdf *mlPdf=sPdf->getPdf(catName);
  if (!mlPdf) return theComp;
  RooArgSet *mlSet=mlPdf->getComponents();
  theComp=mlSet->find(argName);
  return theComp; // return it anyway
}

/// \brief Construct component-specific generator
/// \param gen The orignal generator
/// \param subGens Sub-mode generator
/// \param compCat Component category
///
/// It constructs component-categorized generator.
///
/// \todo Make sure it works for multi physics cats
RooSimultaneous *rarMLFitter::compGen(RooAbsPdf *gen, RooArgList subGens,
				      RooCategory &compCat)
{
  RooSimultaneous *theCompGen(0);
  if (getControlBit("noSimFit")) { // the first mode only
    Int_t nCompTypes=compCat.numTypes();
    theCompGen=new RooSimultaneous
      (Form("simComp_%s", gen->GetName()), "comp-cated gen", compCat);
    RooAddPdf &theMode=(RooAddPdf&)subGens[0];
    RooArgList pdfList(theMode.pdfList());
    RooArgList coefList(theMode.coefList());
    Int_t nComps=pdfList.getSize();
    for (Int_t i=0; i<nComps; i++) {
      RooAddPdf *thePdf=new RooAddPdf
	(Form("comp_%s_%s", theMode.GetName(), pdfList[i].GetName()),
	 Form("comp %s %s", theMode.GetName(), pdfList[i].GetName()),
	 pdfList[i], coefList[i]);
      theCompGen->addPdf(*thePdf, Form("%02d", i));
    }
    // any leftover?
    for (Int_t i=nComps; i<nCompTypes; i++)
      theCompGen->addPdf(_dummyExtPdf, Form("%02d", i));
    
    return theCompGen;
  }
  // now for sim generator
  RooAbsCategoryLValue *simCats=(RooAbsCategoryLValue *)
    (&((RooSimultaneous*)gen)->indexCat());
  RooArgSet *genSimSet=gen->getComponents();
  // get superCat inputList
  RooArgSet myInputCats(compCat);
  // get genCat inputCat
  RooSuperCategory *simSuperCats=dynamic_cast<RooSuperCategory*>(simCats);
  if (simSuperCats) myInputCats.add(simSuperCats->inputCatList());
  else myInputCats.add(*simCats);
  // create superCat
  RooSuperCategory *compGenCat=new RooSuperCategory
    ("compGenCat", "compGenCat", myInputCats);
  // create compGen
  theCompGen=new RooSimultaneous
    (Form("simComp_%s", gen->GetName()), "comp-cated gen", *compGenCat);
  Int_t nCats=compGenCat->numTypes();
  for (Int_t i=0; i<nCats; i++) { // each component
    TString catName=((RooCatType*)compGenCat->lookupType(i))->GetName();
    TString compLabel=catName(1,catName.First(";")-1);
    Int_t compIdx=atoi(compLabel);
    // find genCat name
    TString gCatName=catName
      (catName.First(";")+1, catName.Last('}')-catName.First(";")-1);
    if (simSuperCats) gCatName="{"+gCatName+"}";
    // construct extended add
    RooAbsPdf *thisCatGen(0);
    for (Int_t j=0; j<_nComp; j++) { // sub-mode
      RooAddPdf &theMode=(RooAddPdf&)subGens[j];
      if (compIdx>=theMode.pdfList().getSize()) continue;
      // find comp pdf and coef
      TString thePdfName=theMode.pdfList()[compIdx].GetName();
      TString theCoefName=theMode.coefList()[compIdx].GetName();
      RooAbsPdf *thisPdf=(RooAbsPdf*)findSimed(*genSimSet,thePdfName,gCatName,
					       gen);
      RooAbsReal *thisCoef=(RooAbsReal*)
	findSimed(*genSimSet, theCoefName, gCatName, gen);
      if (thisPdf&&thisCoef) {
	thisCatGen=new RooAddPdf
	  ("comp_"+catName+"_"+thePdfName+"_"+theCoefName,
	   "comp "+catName+" "+thePdfName+" "+theCoefName, *thisPdf,*thisCoef);
	break;
      }
    }
    if (thisCatGen) theCompGen->addPdf(*thisCatGen, catName);
    else theCompGen->addPdf(_dummyExtPdf, catName);
  }
  
  theCompGen->Print();
  //theCompGen->Print("v");
  RooArgSet nullDS;
  cout << "Expected number of events generated : " 
       << theCompGen->expectedEvents(&nullDS)<<endl;
  for (Int_t i=0; i<nCats; i++) {
    TString catName=((RooCatType*)compGenCat->lookupType(i))->GetName();
    compGenCat->setLabel(catName);
    cout<<"\t"<<catName<<" "
	<<theCompGen->expectedEvents(&nullDS)<<endl;
  }

  return theCompGen;
}

/// \brief Get extended component pdf
/// \param thePdf Component pdf
/// \param theCoef Component coef (yield)
/// \return The extended (Sim)Pdf.
RooAbsPdf *rarMLFitter::getExtCompPdf(RooAbsPdf *thePdf, RooAbsReal *theCoef)
{
  RooAbsPdf *theComp(0);
  TString thePdfName=thePdf->GetName();
  TString theCoefName=theCoef->GetName();
  TString theCompName="ExtComp_"+thePdfName+"_"+theCoefName;
  if (theComp=(RooAbsPdf*)_extCompPdfs.find(theCompName)) return theComp;
  // if simFit?
  if (!getControlBit("SimFit")) {
    theComp=new RooAddPdf
      (theCompName, theCompName, RooArgList(*thePdf), RooArgList(*theCoef));
    _extCompPdfs.add(*theComp);
    return theComp;
  }
  // create extended sim comp
  RooAbsCategoryLValue *theCat=(RooAbsCategoryLValue *)
    (&((RooSimultaneous*)_thePdf)->indexCat());
  RooArgSet *theSimSet=_thePdf->getComponents();
  theComp=new RooSimultaneous(theCompName, theCompName, *theCat);
  _extCompPdfs.add(*theComp);
  TIterator* catTypeIter = theCat->typeIterator();
  RooCatType *theType(0);
  while(theType=(RooCatType*)catTypeIter->Next()) {
    TString theTypeName=theType->GetName();
    RooAbsPdf *thisPdf=(RooAbsPdf*)
      findSimed(*theSimSet, thePdfName, theTypeName);
    RooAbsReal *thisCoef=(RooAbsReal*)
      findSimed(*theSimSet, theCoefName, theTypeName);
    RooAbsPdf *thisCatPdf=&_dummyExtPdf;
    if (thisPdf&&thisCoef) {
      thisCatPdf=new RooAddPdf
	("extComp_"+theTypeName+"_"+thePdfName+"_"+theCoefName,
	 "extComp "+theTypeName+" "+thePdfName+" "+theCoefName,
	 *thisPdf, *thisCoef);
    }
    //thisCatPdf->Print();
    //thisCatPdf->Print("v");
    ((RooSimultaneous *)theComp)->addPdf(*thisCatPdf, theTypeName);
  }
  
  delete catTypeIter;
  return theComp;
}

/// \brief Get S and B LL
///
/// It constructs S and B likelihood for LLR plot and proj plot.
void rarMLFitter::getSnB()
{
  if (_theSPdf) return;
  // get proj components
  rarStrParser compStrParser=readConfStr("projComps", "", _runSec);
  if (compStrParser.nArgs()<=0) {
    cout<<"no projection component defined with config \"projComps\""
        <<" in section \""<<_runSec<<"\""<<endl;
    exit(-1);
  }
  // S and B argList
  RooArgList sList, bList;
  Int_t nComp=1;
  if (getControlBit("SimFit")) nComp=_nComp;
  for (Int_t i=0; i<nComp; i++) {
    rarMLPdf *mlPdf=(rarMLPdf*)_pdfList.At(i);
    TList *pdfList=mlPdf->getPdfList();
    Int_t pdfSize=pdfList->GetSize();
    RooArgList coefList(mlPdf->getCoeffList());
    for (Int_t j=0; j<pdfSize; j++) {
      rarBasePdf *thePdf=(rarBasePdf*)pdfList->At(j);
      RooAbsReal *theCoef=(RooAbsReal*)coefList.at(j);
      _fCompList.Add(thePdf);
      _fCoefList.add(*theCoef);
      RooAbsPdf *theExtCompPdf=getExtCompPdf(thePdf->getPdf(), theCoef);
      // S or B?
      if (compStrParser.Have(thePdf->GetName())) {
	sList.add(*theExtCompPdf);
	_sCompList.Add(thePdf);
	_sCoefList.add(*theCoef);
      } else bList.add(*theExtCompPdf);
    }
  }
  if ((sList.getSize()<1)||(bList.getSize()<1)) {
    cout<<"Can not find any comp for projPlot"<<endl;
    exit(-1);
  }
  // create S and B
  _theSPdf=new RooAddPdf("S_Pdf", "S Pdf", sList);
  _theBPdf=new RooAddPdf("B_Pdf", "B Pdf", bList);
  
  return;
}

/// \brief Get dataset with LLR added
/// \param theData The dataset to add LLR
/// \param LLR LLR function to add
/// \return The dataset created
///
/// It creates dataset with LLR added
RooDataSet *rarMLFitter::getLLRDataset(RooDataSet *theData, RooFormulaVar *LLR)
{
  RooDataSet *retDS=new
    RooDataSet(*theData,Form("%s_wLLR", theData->GetName()));
  
  // add limits to the newly created dataset,
  // because the derived column, LLR, needs those limits to get normalization
  setColLimits(retDS);
  // now add the derived column
  //LLR->Print();
  //LLR->Print("v");
  //retDS->Print();
  //retDS->Print("v");
  retDS->addColumn(*LLR);
  // now remove the limits
  setColLimits(retDS, kFALSE);
  // if removing limits of derived columns is keeping causing problems,
  // we should not remove the limits any longer,
  // because if the limits are causing the problems,
  // it usually should be fixed by changing the deriving formula,
  // or other method. It is reasonable to have limits.
  
  return retDS;
}

/// \brief Get LLR plot
/// \param projData Dataset to fit for mlFit action
/// \param plotList Plot list
///
/// It is a simplified version of #doProjPlot to get LLR plots
/// based on full models.
void rarMLFitter::doLLRPlot(RooDataSet *projData, TList &plotList)
{
  cout<<endl<<" In rarMLFitter doLLRPlot for "<<GetName()<<endl;
  TString llrStr=readConfStr("projLLRPlots", "yes", _runSec);
  if ("no"==llrStr) return;
  // dummy projection settings
  RooArgSet projVarSet;
  RooArgSet depVarSet(_fObsSet);
  depVarSet.remove(_conditionalObs,kFALSE,kTRUE);
  // construct lRatio
  const RooAbsReal *theSProj=_theSPdf->createPlotProjection(depVarSet, projVarSet);
  const RooAbsReal *theBProj=_theBPdf->createPlotProjection(depVarSet, projVarSet);
  RooArgList projPdfList(*theSProj, *theBProj, "Projected Pdf List");
  TString formulaStr="@0/(@0+@1)";
  TString funcName="LLR";
  RooFormulaVar *lRatioFunc=new
    RooFormulaVar(funcName, "LRatio Func", formulaStr, projPdfList);
  // get bins
  Int_t nBins=atoi(readConfStr("plotBins_LLR", "100", _runSec));
  // now for projection dataset
  TH1F *LLRP0=new TH1F(Form("LLR_%s", projData->GetName()),"LLR ds",nBins,0,1);
  projData=getLLRDataset(projData, lRatioFunc);
  for (Int_t idxEvt=0; idxEvt<projData->numEntries(); idxEvt++) {
    RooArgSet *dEntry=(RooArgSet*)projData->get(idxEvt);
    LLRP0->Fill(((RooRealVar*)dEntry->find(funcName))->getVal(),
                projData->weight());
  }
  LLRP0->Print("v");
  plotList.Add(LLRP0);
  RooDataSet *protData(0);
  /// \todo this _protDataVars should be replaced with fullProtVars
  if (_protDataVars.getSize()>0)
    protData=(RooDataSet*)projData->reduce(_protDataVars);
  RooArgSet theDeps(*_thePdf->getDependents(projData));
  if (protData) theDeps.remove(*_thePdf->getDependents(protData));
  
  // now for total model
  TString projLLRScale=
    readConfStr(Form("projLLRScale_%s", GetName()), "-1", _runSec);
  if ("-1"==projLLRScale)
    projLLRScale=readConfStr("projLLRScale", "10", _runSec);
  Double_t scale=atof(projLLRScale);
  Int_t nEvt=(Int_t)(scale*projData->numEntries());
  if (nEvt>0) {
    getLLRPlot(plotList, GetName(), _thePdf, nEvt,
	       theDeps, protData, nBins, lRatioFunc);
  } else {
    cout<<" Do not generate sample for LLR plot of "<<GetName()<<endl;
  }
  // for comp in projComps and projLLRPlots
  TString llrCompStr="";
  if ("yes"!=llrStr) llrCompStr=llrStr;
  else for (Int_t i=0; i<_fCompList.GetSize(); i++)
    llrCompStr+=Form(" \"%s\"", _fCompList.At(i)->GetName());
  rarStrParser llrStrParser=llrCompStr;
  while (llrStrParser.nArgs()>0) {
    TString compName=llrStrParser[0];
    llrStrParser.Remove();
    rarBasePdf *compPdf=(rarBasePdf*)_fCompList.FindObject(compName);
    if (!compPdf) continue;
    Int_t compIdx=_fCompList.IndexOf(compPdf);
    RooAbsPdf *thePdf=compPdf->getSimPdf();
    if (!thePdf) thePdf=compPdf->getPdf();
    // scale
    TString projLLRScale=
      readConfStr(Form("projLLRScale_%s", compPdf->GetName()), "-1", _runSec);
    if ("-1"==projLLRScale)
      projLLRScale=readConfStr("projLLRScale", "10", _runSec);
    Double_t scale=atof(projLLRScale);
    // find #evt
    Int_t nEvt=(Int_t)(scale*((RooAbsReal&)_fCoefList[compIdx]).getVal());
    if (nEvt<=0) {
      cout<<" Negative yield for component "<<compPdf->GetName()
	  <<", Ignored!!!"<<endl;
      continue;
    }
    getLLRPlot(plotList, compPdf->GetName(), thePdf, nEvt,
	       theDeps, protData, nBins, lRatioFunc);
  }
  // delete the created projData
  delete projData;
  
  return;
}

/// \brief Get LLR histogram
/// \param plotList Plot list
/// \param plotName Plot Name
/// \param thePdf Pdf to get LLR plot
/// \param nEvt Number of event to generate for the plot
/// \param theDeps ArgSet to generate
/// \param protData Prototype dataset
/// \param nBins Number of bins
/// \param LLRFunc LLR function
///
/// It creates LLR histogram from pdf specified.
void rarMLFitter::getLLRPlot(TList &plotList, TString plotName,
			     RooAbsPdf *thePdf, Int_t nEvt,
			     RooArgSet theDeps, RooDataSet *protData,
			     Int_t nBins, RooAbsReal *LLRFunc)
{
  TString histName="LLR_"+plotName;
  if (plotList.FindObject(histName)) return;
  TH1F *theHist(0);
  // generate dataset with the non-fancy generator
  RooDataSet *theData(0);
  if (protData)
    theData=thePdf->generate(theDeps,*protData,nEvt,kFALSE,kTRUE);
  else theData=thePdf->generate(theDeps, nEvt);
  theData->addColumn(*LLRFunc);
  theHist=new TH1F(histName, "LLR hist", nBins,0,1);
  for (Int_t idxEvt=0; idxEvt<theData->numEntries(); idxEvt++) {
    RooArgSet *dEntry=(RooArgSet*)theData->get(idxEvt);
    theHist->Fill(((RooRealVar*)dEntry->find(LLRFunc->GetName()))->getVal(),
                  theData->weight());
  }
  theHist->Print("v");
  plotList.Add(theHist);
  delete theData;
  return;
}


/// \brief Get projection plot
/// \param projData Dataset to fit for mlFit action
/// \param theVar obs to plot
/// \param plotList Plot list
/// \return The last plot created
///
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_projPlot">See doc for projection plot action.</a>
///
/// It creates projection plots of \p projVar
/// for default dataset, total pdf, and each component, on one frame.
RooPlot *rarMLFitter::doProjPlot(RooDataSet *projData, RooRealVar *theVar,
				 TList &plotList)
{
  cout<<endl<<" In rarMLFitter doProjPlot of "<<theVar->GetName()
      <<" for "<<GetName()<<endl;
  RooPlot *frame(0);
  // theVar name
  TString varName=theVar->GetName();
  // get proj plot data
  TString ppdStr=readConfStrCnA("projPlotData_"+varName,"");
  RooDataSet *projPlotData(0);
  if (""!=ppdStr) projPlotData=_datasets->getData(ppdStr);
  if(projPlotData) chkBlind(projPlotData->GetName());
  if(!projPlotData) projPlotData=projData;
  cout<<" Using dataset "<<projPlotData->GetName()<<endl;
  RooArgSet projVarSet;
  //projVarSet.add(*theVar); // JA
  RooArgSet depVarSet(_fObsSet);
  //depVarSet.remove(projVarSet, kFALSE, kTRUE); // JA
  depVarSet.remove(*theVar, kFALSE, kTRUE); // JA
  depVarSet.remove(_conditionalObs, kFALSE, kTRUE);
  // construct lRatio
  //_theSPdf->Print();
  //_theSPdf->Print("v");
  //_theBPdf->Print();
  //_theBPdf->Print("v");

  TString LRatioCut=readConfStr("projLRatioCut_"+varName, "-1", _runSec);
  if ("-1"==LRatioCut) LRatioCut=readConfStr("projLRatioCut", "-1", _runSec);
  // see if we need to find the optimal cut
  TString findOptCut=readConfStr("projFindOptimCut_"+varName,"notSet",_runSec);
  if ("notSet"==findOptCut)
    findOptCut=readConfStr("projFindOptimCut", "no", _runSec);

  RooDataSet* sliceData = 0;
  bool delDS = false;

  if ( LRatioCut!="-1" || findOptCut!="no" ) {
    const RooAbsReal *theSProj=_theSPdf->createPlotProjection(depVarSet, projVarSet);
    const RooAbsReal *theBProj=_theBPdf->createPlotProjection(depVarSet, projVarSet);
    RooArgList projPdfList(*theSProj, *theBProj, "Projected Pdf List");
    TString formulaStr="@0/(@0+@1)";
    TString funcName="lRatioFunc_"+varName;
    RooFormulaVar *lRatioFunc=new
      RooFormulaVar(funcName, "LRatio Func", formulaStr, projPdfList);
  
    projPlotData=getLLRDataset(projPlotData, lRatioFunc);
    delDS = true;
    // get the cut
    if ("no"!=findOptCut) {
      cout<<"Finding optimal cut"<<endl;
      // set steps
      TString projOptStep=readConfStr("projOptimStep_"+varName, "-1", _runSec);
      if ("-1"==projOptStep)
	projOptStep=readConfStr("projOptimStep", ".005", _runSec);
      Int_t sStep=(Int_t)(.5+1./atof(projOptStep));
      if (sStep<=0) sStep=10;
      cout<<"Searching step: "<<1./sStep<<", ie "<<sStep<<" steps"<<endl;
      // do we have any cut from projPlotData
      TString ppdCut="1";
      rarStrParser cutParser=ppdStr;
      if (cutParser.nArgs()>=2) ppdCut=cutParser[1];
      // do we have any range cuts
      TString optRangeCut=readConfStr("projOptimRange_"+varName, "-1", _runSec);
      if ("-1"==optRangeCut)
	optRangeCut=readConfStr("projOptimRange", "1", _runSec);
      if (("1"!=optRangeCut)||("1"!=ppdCut)) {
	// remove any " in cut string
	optRangeCut.ReplaceAll("\"", "");
	optRangeCut="("+optRangeCut+")&&("+ppdCut+")";
	cout<<"With cuts "<<optRangeCut<<endl;
      }
      // reduced total sample
      cout<<" Reducing "<<projPlotData->GetName()<<endl;
      RooDataSet *theData=(RooDataSet*)projPlotData->reduce(optRangeCut);
      // create histogram for ratio searching
      TH1F *nTotHist=new TH1F("nTotHist_"+varName,"nTotHist_"+varName,sStep,0,1);
      cout<<" Creating nTotHist_"+varName<<endl;
      for (Int_t idxEvt=0; idxEvt<theData->numEntries(); idxEvt++) {
	RooArgSet *dEntry=(RooArgSet*)theData->get(idxEvt);
	nTotHist->Fill(((RooRealVar*)dEntry->find(funcName))->getVal(),
		       theData->weight());
      }
      delete theData;
      nTotHist->Print("v");
      plotList.Add(nTotHist);
      // get datasets for efficiency calculation
      rarStrParser optDataStrParser=readConfStr("projOptimData", "", _runSec);
      // get max size of the dataset
      Int_t projDataSize=atoi(readConfStr("projOptimDataLimit","40000",_runSec));
      if (projDataSize<=1000) projDataSize=1000;
      // get histograms for components
      TList hList;
      RooArgSet fullObsWOEVars(*_fullObs);
      fullObsWOEVars.remove(_protDataEVars, kFALSE, kTRUE);
      const char *dummy("");
      for (Int_t i=0; i<_sCompList.GetSize(); i++) {
	rarBasePdf *rarPdf=(rarBasePdf*)_sCompList.At(i);
	RooDataSet *theData=rarPdf->getData(dummy);
	if (optDataStrParser.nArgs()>0) {
	  TString theDataName=optDataStrParser[0];
	  optDataStrParser.Remove();
	  theData=_datasets->getData(theDataName);
	  if (!theData) {
	    cout<<"Can not find dataset named "<<theDataName<<endl
		<<"Use the default dataset"<<endl;
	    theData=rarPdf->getData(dummy);
	  }
	}
	Bool_t reducedData(kFALSE);
	// now get a smaller subset of the dataset
	if ((theData->numEntries()>projDataSize)||(_protDataEVars.getSize()>0)) {
	  cout<<" Reducing "<<theData->GetName()<<endl;
	  Int_t nEvt=theData->numEntries();
	  if (nEvt>projDataSize) nEvt=projDataSize;
	  theData=(RooDataSet*)
	    theData->reduce(SelectVars(fullObsWOEVars),EventRange(0, nEvt-1));
	  reducedData=kTRUE;
	}
	// make sure protDataEVars in theData are from prototype data
	if (_protDataEVars.getSize()>0) {
	  RooDataSet protEData("protEData", "protEData", _protDataEVars);
	  Int_t nEvt=theData->numEntries();
	  Int_t nProtEvt=projPlotData->numEntries();
	  cout<<" Getting protDataEVars from "<<projPlotData->GetName()<<endl;
	  for (Int_t j=0; j<nEvt; j++) {
	    Int_t I=RooRandom::randomGenerator()->Integer(nProtEvt);
	    protEData.add(*projPlotData->get(I));
	  }
	  cout<<" Merging protDataEVars to "<<theData->GetName()<<endl;
	  theData->merge(&protEData);
	}
	cout<<" Adding LRatio to "<<theData->GetName()<<endl;
	RooDataSet *tmpData=theData;
	RooDataSet *theDataNew=getLLRDataset(theData, lRatioFunc);
	theData=theDataNew;
	Double_t sigScale=theData->numEntries();
	theData=(RooDataSet*)theData->reduce(optRangeCut);
	if (reducedData) delete tmpData; // delete if created here
	sigScale=theData->numEntries()/sigScale;
	TString hName=Form("hist_%s_%s", rarPdf->GetName(), varName.Data());
	TH1F *h=new TH1F(hName, hName, sStep, 0, 1);
	cout<<" Creating histogram "<<h->GetName()<<endl;
	for (Int_t idxEvt=0; idxEvt<theData->numEntries(); idxEvt++) {
	  RooArgSet *dEntry=(RooArgSet*)theData->get(idxEvt);
	  h->Fill(((RooRealVar*)dEntry->find(funcName))->getVal(),
		  theData->weight());
	}
	plotList.Add(h);
	h->Scale(sigScale/h->GetEntries());
	h->Print("v");
	hList.Add(h);
	delete theData;
	delete theDataNew;
      }
      // set optimization formula, default (n1+n2+...)^2/ntotal
      TString optimFormulaStr=
	readConfStr("projOptimFormula_"+varName, "notSet", _runSec);
      if ("notSet"==optimFormulaStr)
	optimFormulaStr=readConfStr("projOptimFormula", "notSet", _runSec);
      if ("notSet"==optimFormulaStr) optimFormulaStr="@0*@0/@1";
      optimFormulaStr.ReplaceAll("\"", "");
      RooRealVar NprojVar("NprojVar", "NprojVar", 0);
      RooRealVar NtotalVar("NtotalVar", "NtotalVar", 1);
      RooFormulaVar optimFormula("optimFormula", optimFormulaStr,
				 RooArgList(NprojVar, NtotalVar));
      cout<<" Optimization formula: "<<optimFormulaStr<<endl;
      Double_t bestRatio(0.), bestCut(0.);
      for (Int_t i=0; i<sStep+1; i++) {
	//Double_t nTot=nTotHist->Integral(i, sStep);
	//if (nTot<=0) continue;
	NtotalVar.setVal(nTotHist->Integral(i, sStep));
	if (NtotalVar.getVal()<=0) continue;
	Double_t nSig(0.);
	for (Int_t j=0; j<_sCoefList.getSize(); j++) {
	  nSig+=((RooAbsReal &)_sCoefList[j]).getVal()*
	    ((TH1F*)hList.At(j))->Integral(i, sStep);
	}
	NprojVar.setVal(nSig);
	//Double_t ratio=nSig*nSig*nSig/nTot;
	//cout<<i<<" nSig:"<<nSig<<" nTot:"<<nTot<<" ratio:"<<ratio<<endl;
	Double_t ratio=optimFormula.getVal();
	cout<<i<<" nSig:"<<NprojVar.getVal()<<" nTot:"<<NtotalVar.getVal()
	    <<" ratio:"<<ratio<<endl;
	if (ratio>bestRatio) {
	  bestRatio=ratio;
	  bestCut=i;
	}
      }
      bestCut=bestCut/((Double_t)sStep);
      LRatioCut=Form("%f", bestCut);
      cout<<"optimal cut @ "<<bestCut<<" with n^2/ntot="<<bestRatio<<endl;
    }
    LRatioCut=funcName+">"+LRatioCut;
    // reduce the dataset
    sliceData=(RooDataSet*)projPlotData->reduce(LRatioCut);
  } else {
    sliceData = (RooDataSet*)projPlotData->Clone(); // to avoid SEGV on delete
  }
  // get plot bins
  Int_t nBins=atoi(readConfStr(Form("plotBins_%s", theVar->GetName()),
			       "0", _runSec));
  if (nBins<=0) nBins=theVar->getBins();
  // get X error scaling
  Double_t xerrorscale = atof(readConfStr("XErrorSize", "1.0", _runSec));
  // plot range
  Double_t plotMin=theVar->getMin();
  Double_t plotMax=theVar->getMax();
  RooBinning *abins=
    getRange(theVar, "plotRange_", plotMin, plotMax, _runSec, &nBins);
  abins->Print();
  // plot name and title
  TString frameName=Form("proj_%s", theVar->GetName());
  TString frameTitle;
  if ( LRatioCut != "" ) {
    frameTitle =
      Form("projection of %s with %s",theVar->GetTitle(), LRatioCut.Data());
  } else {
    frameTitle =  Form("projection of %s",theVar->GetTitle());
  }
  // do we need to have asym plot?
  TString projAsymPlot=readConfStr("projAsymPlot_"+varName,"notSet", _runSec);
  if ("notSet"==projAsymPlot)
    projAsymPlot=readConfStr("projAsymPlot","no", _runSec);
  // Lists for proj datasets
  TList projDS;
  TList projNames;
  projDS.Add(sliceData);
  projNames.Add(new TObjString(frameName));
  // proj cat
  TString projPlotCat=readConfStr("projPlotCat_"+varName,"notSet", _runSec);
  if ("notSet"==projPlotCat)
    projPlotCat=readConfStr("projPlotCat","no", _runSec);
  if ("no"!=projPlotCat) {
    rarStrParser projPlotCatParser=projPlotCat;
    while (projPlotCatParser.nArgs()>0) {
      TString catName=projPlotCatParser[0];
      projPlotCatParser.Remove();
      RooAbsCategory *projCat=(RooAbsCategory *)_fullObs->find(catName);
      if (!projCat) {
	cout<<" W A R N I N G !"<<endl
	    <<" category "<<catName<<" does not exist!"<<endl;
	continue;
      }
      if (catName==projAsymPlot) {
	cout<<" "<<catName<<" is used for asymPlot"<<endl
	    <<" Ignored"<<endl;
	continue;
      }
      TIterator *cIter=projCat->typeIterator();
      RooCatType *cType(0);
      while(cType=(RooCatType*)cIter->Next()) {
	TString catCut=catName+"=="+catName+"::"+cType->GetName();
	RooDataSet *catSliceData=(RooDataSet*)sliceData->reduce(catCut);
	projDS.Add(catSliceData);
	projNames.Add(new TObjString(frameName+"_"+catCut));
      }
      delete cIter;
    }
  }
  // projection plot for each dataset
  for (Int_t pIdx=0; pIdx<projDS.GetSize(); pIdx++) {
    sliceData=(RooDataSet*)projDS.At(pIdx);
    frameName=projNames.At(pIdx)->GetName();
    if ("no"==projAsymPlot) {
      frame=getProjPlot(theVar, plotMin, plotMax, nBins, sliceData,
		     frameName, frameTitle, plotList, RooCmdArg(), RooCmdArg(),
		      XErrorSize(xerrorscale));
    } else { // asym plot
      // check cat
      RooCategory *cat=(RooCategory *)_fullObs->find(projAsymPlot);
      if (!cat) {
	cout<<"Can not find category: "<<projAsymPlot<<endl;
	exit(-1);
      }
      // make sure it is two-type cat
      if (2!=cat->numTypes()) {
	cout<<"Asym plot is for two-type category only"<<endl;
	exit(-1);
      }
      TIterator *catI=cat->typeIterator();
      RooCatType *theCatType=(RooCatType*)catI->Next();
      TString catCutA=
        projAsymPlot+"=="+projAsymPlot+"::"+theCatType->GetName();
      Int_t catIdxA=theCatType->getVal();
      theCatType=(RooCatType*)catI->Next();
      TString catCutB=
        projAsymPlot+"=="+projAsymPlot+"::"+theCatType->GetName();
      Int_t catIdxB=theCatType->getVal();
      TString catCutMinus=catIdxA<catIdxB?catCutA:catCutB;
      TString catCutPlus =catIdxA>catIdxB?catCutA:catCutB;
      RooDataSet *catSliceData=(RooDataSet*)sliceData->reduce(catCutMinus);
      cout<<"Dataset (-) with cut: "<<catCutMinus<<endl;
      catSliceData->Print();
      RooPlot *frameM=
	getProjPlot(theVar, plotMin, plotMax, nBins, catSliceData,
		    frameName+"_Minus", frameTitle+" Minus", plotList,
		    RooCmdArg(), RooCmdArg(), XErrorSize(xerrorscale));
      delete catSliceData;
      catSliceData=(RooDataSet*)sliceData->reduce(catCutPlus);
      cout<<"Dataset (+) with cut: "<<catCutPlus<<endl;
      catSliceData->Print();
      RooPlot *frameP=
	getProjPlot(theVar, plotMin, plotMax, nBins, catSliceData,
		    frameName+"_Plus", frameTitle+" Plus", plotList,
		     RooCmdArg(), RooCmdArg(), XErrorSize(xerrorscale));
      delete catSliceData;
      delete catI;
      // now asym plot
      frame=getProjPlot(theVar, plotMin, plotMax, nBins, sliceData,
			frameName+"_Asym", frameTitle+" Asym", plotList,
			Asymmetry(*cat), Binning(*abins), 
			XErrorSize(xerrorscale), frameM, frameP);
    }
    sliceData->Print("v");
  }
  projDS.Delete();
  projNames.Delete();
  
  // do we need to save the dataset as well
  TString saveDS=readConfStrCnA("projPlotSaveLLR", "no");
  if ("no"==saveDS) {
    if ( delDS ) delete projPlotData;
  } else {
    TString dsName=projPlotData->GetName();
    dsName+="_";
    dsName+=theVar->GetName();
    projPlotData->SetName(dsName);

    TTree *theDSTree = createTreeFromDataset(projPlotData, kFALSE);
    plotList.Add(theDSTree);
  }
  
  cout<<endl<<"lRatio Cut: "<<LRatioCut<<endl;
  return frame;
}

/// \brief Get projection plot frame
/// \param theVar obs to plot
/// \param plotMin Plot min
/// \param plotMax Plot max
/// \param nBins Number of bins
/// \param sliceData Reduced dataset
/// \param frameName Frame name
/// \param frameTitle Frame title
/// \param plotList Plot list
/// \param asymCat RooCmdArg for asym plot
/// \param nuBins Non-uniform binning
/// \param xerrorscale Scale of X errors with respect to bin width
/// \param frameM Frame of Minus type
/// \param frameP Frame of Plus type
/// \return The frame created
///
/// This the actual function to draw the projection plots.
/// It is separated frome #doProjPlot for more flexibility
/// in plotting.
RooPlot *rarMLFitter::getProjPlot(RooRealVar *theVar, Double_t plotMin,
				  Double_t plotMax, Int_t nBins,
				  RooDataSet *sliceData, TString frameName,
				  TString frameTitle, TList &plotList,
				  const RooCmdArg &asymCat,
				  const RooCmdArg &nuBins,
				  const RooCmdArg &xerrorscale,
				  RooPlot *frameM, RooPlot *frameP)
{
  // first add limits to addOns for sliceData
  setColLimits(sliceData);
  RooPlot *frame(0);
  frame=theVar->frame(plotMin, plotMax, nBins);
  RooLinkedList plotOpts;
  plotOpts.Add((TObject*)&asymCat);
  plotOpts.Add((TObject*)&nuBins);
  plotOpts.Add((TObject*)&xerrorscale);
  RooCmdArg datErrArg = DataError(RooAbsData::SumW2);
  if (sliceData->isWeighted())
    plotOpts.Add((TObject*)&datErrArg);
  sliceData->plotOn(frame, plotOpts);
  Int_t lineStyle=1;
  Int_t lineWidth=2;
  Int_t lineColor=kBlue;
  if (!(frameM&&frameP)) { // plot pdf and components
    //_thePdf->plotOn(frame, ProjWData(_conditionalObs, *sliceData),asymCat,
    _thePdf->plotOn(frame, ProjWData(*sliceData),asymCat,
                    LineStyle(lineStyle), LineColor(lineColor));
    for (Int_t i=0; i<_fCompList.GetSize(); i++) {
      if (lineStyle>2) lineStyle=2;
      else lineStyle=4;
      rarBasePdf *compPdf=(rarBasePdf*)_fCompList.At(i);
      TString pdfName=compPdf->getPdf()->GetName();
      TString compStr=pdfName+","+pdfName+"_*";
      cout<<"Comps to proj: "<<compStr<<endl;
      _thePdf->plotOn(frame, Components(compStr), ProjWData(*sliceData),
		      asymCat, LineWidth(lineWidth), LineStyle(lineStyle),
		      LineColor(getColor(i)), Name(pdfName));
    }
  } else { // plot asym between two sets of pdf components
    Int_t nCurves=(Int_t)frameP->numItems();
    for (Int_t i=1; i<nCurves; i++) {
      RooCurve *asymCurve(0);
      RooCurve *B0Curve(0), *B0barCurve(0);
      B0Curve= (RooCurve*) frameP->getObject(i);
      B0barCurve = (RooCurve*) frameM->getObject(i);
      if (B0Curve&&B0barCurve) {
	TString curveName=Form("asymCurve_%s", B0Curve->GetName());
	Int_t lineWidth=B0Curve->GetLineWidth();
	Int_t lineStyle=B0Curve->GetLineStyle();
	Int_t lineColor=B0Curve->GetLineColor();
	asymCurve=combCurves(B0Curve, B0barCurve, "(@0-@1)/(@0+@1)", "@0+@1>0",
			     curveName, lineWidth, lineStyle, lineColor);
      }
      if (asymCurve)
	frame->addPlotable(asymCurve);
    }
  }
  
  frame->SetNameTitle(frameName, frameTitle);
  //#ifndef RAR_USE_ROOT5
  frame->SetTitleSize(0.05);
  //#endif
  plotList.Add(frame);
  // remove limits
  setColLimits(sliceData, kFALSE);
  return frame;
}

/// \brief Combine two RooCurves by a formula
/// \param crv1 The first RooCurve
/// \param crv2 The second RooCurve
/// \param formula Combination formula
/// \param valid Test formula
/// \param crvName Name of the new curve
/// \param lineWidth Line width of the new curve
/// \param lineStyle Line style of the new curve
/// \param lineColor Line color of the new curve
/// \return Combined curve
///
/// This funciton combines two RooCurves with formula.
/// It is taken from Q2BPlots::combCurves()
RooCurve *rarMLFitter::combCurves(RooCurve *crv1, RooCurve *crv2,
				  const char *formula, const char *valid,
				  const char *crvName, Int_t lineWidth,
				  Int_t lineStyle, Int_t lineColor)
{
  RooCurve* newCurve = new RooCurve();
  if (0 != crvName) newCurve->SetName(crvName);
  newCurve->SetLineWidth(lineWidth);
  newCurve->SetLineStyle(lineStyle);
  newCurve->SetLineColor(lineColor);
  Double_t x1 = -1.E30;
  Double_t x2 = -1.E30;
  Double_t y1, y2;
  Int_t n1 = crv1->GetN();
  Int_t n2 = crv2->GetN();
  RooRealVar y1RV("y1RV", "", 0.), y2RV("y2RV", "", 0.);
  RooFormulaVar yComb("yComb", "", formula, RooArgSet(y1RV, y2RV));
  RooFormulaVar testOK("testOK", "", valid, RooArgSet(y1RV, y2RV));
  
  bool needPt = true;
  int i(0), j(0);
  while (i < n1-1 && j < n2-1) {
    while (i < n1-1 && (needPt || x1 < x2)) {
      needPt = false;
      crv1->GetPoint(++i, x1, y1);
      y1RV.setVal(y1);
    }
    while (j < n2-1 && (needPt || x2 < x1)) {
      needPt = false;
      crv2->GetPoint(++j, x2, y2);
      y2RV.setVal(y2);
    }
    if (x1 == x2 && 0 != testOK.getVal()) {
      //cout << i << " " << j << " " << x1 << " " << x2
      //<< "  " << y1 << "  " << y2 << endl;
      newCurve->addPoint(x1, yComb.getVal());
    }
    needPt = true;
  }
  //  newCurve->Print("v");
  
  return newCurve;
}

/// \brief Shift the fixed values back to the nominal mlFit values
/// \param scanVars Vars to scan
/// \param scanVarDiff Array of diffs
///
/// It shifts the fixed values back to the nominal mlFit values
void rarMLFitter::scanVarShiftToNorm(RooArgList scanVars, TArrayD &scanVarDiff)
{
  if ("no"==readConfStr("scanVarShiftToNorm", "no", _runSec)) return;
  for (Int_t i=0; i<scanVars.getSize(); i++) {
    RooRealVar *theVar=dynamic_cast<RooRealVar*>(scanVars.at(i));
    if (!theVar) continue;
    theVar->setVal(theVar->getVal()+scanVarDiff[i]);
    //theVar->Print("v");
  }
}

/// \brief Get scan plot
/// \param plotList Plot list
/// \return The plot, if any, created
///
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_scanPlot">See doc for scan plot action.</a>
///
/// The function scans the allowed ranges for specified vars randomly
/// to get NLL points for scan plot.
RooPlot *rarMLFitter::doScanPlot(TList &plotList)
{
  cout<<endl<<" In rarMLFitter doScanPlot for "<<GetName()<<endl;
  RooPlot *frame(0);
  RooCurve *curve(0);
  // get fit option
  TString fitOption=readConfStr("scanPlotFitOption", "qemhr", _runSec);
  cout<<"scanPlot fit option: \""<<fitOption<<"\""<<endl;

  // convert to newer fitTo format for steering options
  Bool_t scanPlotExtended = fitOption.Contains("e");
  Bool_t scanPlotMinos    = fitOption.Contains("m");
  Bool_t scanPlotHesse    = fitOption.Contains("h");
  Bool_t scanPlotVerbose  = !fitOption.Contains("q");
  Bool_t scanPlotSave     = fitOption.Contains("r"); // return results

  // get scan plot data
  RooDataSet *scanPlotData=
    _datasets->getData(readConfStrCnA("scanPlotData", "notSet"));
  if(!scanPlotData) {
    cout<<endl<<" Please specify in section ["<<_runSec<<"]"<<endl
	<<" the dataset you want to scan with config "<<endl<<endl
	<<"   scanPlotData = <datasetName>"<<endl
	<<endl<<" and ru-run your job"<<endl;
    exit(-1);
  }
  chkBlind(scanPlotData->GetName());
  cout<<" Using dataset "<<scanPlotData->GetName()<<endl;
  // full paramters
  RooArgSet fullParams(*_thePdf->getParameters(scanPlotData));
  // save them first
  string paramSStr0;
  writeToStr(fullParams, paramSStr0);
  // get scan plot vars
  RooArgSet scanVars;
  rarStrParser scanVarsParser=readConfStr("scanVars", "notSet", _runSec);
  while (scanVarsParser.nArgs()>0) {
    TString varName=scanVarsParser[0];
    scanVarsParser.Remove();
    RooRealVar *theVar=(RooRealVar*)fullParams.find(varName);
    if (!theVar) continue;
    scanVars.add(*theVar);
    // new min?
    Double_t min=theVar->getMin();
    Double_t max=theVar->getMax();
    // check if the next two are limits for it
    if (isNumber(scanVarsParser[0])) {
      min=atof(scanVarsParser[0]);
      scanVarsParser.Remove();
    }
    if (isNumber(scanVarsParser[0])) {
      max=atof(scanVarsParser[0]);
      scanVarsParser.Remove();
    }
    if (min>max) {
      Double_t v=min;
      min=max;
      max=v;
    }
    theVar->setRange(min,max);
    theVar->setConstant();
  }
  if (scanVars.getSize()<1) {
    cout<<" Please specify var to scan with config scanVars"<<endl;
    exit(-1);
  }
  // if only one var to scan, create RooCurve for it
  if (scanVars.getSize()==1) {
    RooRealVar *theVar=(RooRealVar*)RooArgList(scanVars).at(0);
    frame=theVar->frame(theVar->getMin(), theVar->getMax(), 100);
    frame->SetNameTitle(Form("NLLScanPlot_%s", theVar->GetName()),
			Form("NLLScanPlot_%s", theVar->GetName()));
    frame->SetYTitle("-2 ln (L/L_{0})");
    curve=new RooCurve();
    curve->SetName("NLL_curve");
    frame->addPlotable(curve);
    frame->SetMinimum(0);
    plotList.Add(frame);
  }
  // save params
  string paramSStr;
  writeToStr(fullParams, paramSStr);
  // number of points
  Int_t nPoints=atoi(readConfStr("nScanPoints", "100", _runSec));
  if (nPoints<1) nPoints=1;
  RooRealVar NLL("NLL", "NLL", 0);
  RooArgSet scanSet(scanVars);
  scanSet.add(NLL);
  // Dataset for scan points
  RooDataSet *theDS=new RooDataSet("scanDS", "scanDS", scanSet);
  // save the old values for scanVars
  TArrayD scanVarDiff(scanVars.getSize());
  RooArgList scanList(scanVars);
  for (Int_t i=0; i<scanList.getSize(); i++) {
    scanVarDiff[i]=((RooAbsReal&)scanList[i]).getVal();
  }
  // do an initial fit to find min nll
  Double_t mNLL(0), maxNLL(0);
  scanVars.setAttribAll("Constant", kFALSE);
  cout<<" Refit to find mins for scanPlot"<<endl;
  RooFitResult *fr=_thePdf->fitTo(*scanPlotData, ConditionalObservables(_conditionalObs),
  				  Save(scanPlotSave),Extended(scanPlotExtended), 
				  Verbose(scanPlotVerbose), Hesse(scanPlotHesse),
				  Minos(scanPlotMinos));

  //Int_t ncpus(1);
  //RooFitResult *fr=doTheFit(_thePdf, scanPlotData, fitOption, ncpus);

  if (fr->status()) {
    cout<<" Fit status for inital fit: "<<fr->status()<<endl;
    exit(-1);
  }
  mNLL=2*fr->minNll();
  // set the difference
  for (Int_t i=0; i<scanList.getSize(); i++) {
    scanVarDiff[i]-=((RooAbsReal&)scanList[i]).getVal();
    cout<<" scanVarDiff["<<scanList[i].GetName()<<"]="<<scanVarDiff[i]<<endl;
  }
  // save the min point
  Double_t theVarNormVal=((RooRealVar*)RooArgList(scanVars).at(0))->getVal();
  Double_t theVarNormNLL=2*fr->minNll();
  // fix the obs and refit again for scan points
  {
    scanVars.setAttribAll("Constant");
    RooFitResult *fr=_thePdf->fitTo(*scanPlotData, ConditionalObservables(_conditionalObs),
    				    Save(scanPlotSave),Extended(scanPlotExtended), 
				    Verbose(scanPlotVerbose), Hesse(scanPlotHesse),
				    Minos(scanPlotMinos));
    //Int_t ncpus(1);
    //RooFitResult *fr=doTheFit(_thePdf, scanPlotData, fitOption, ncpus);

    if (fr->status()) {
      cout<<" Fit status for fit: "<<fr->status()<<endl;
      exit(-1);
    }
    // shift fixed values
    scanVarShiftToNorm(scanVars, scanVarDiff);
    // reset the mins
    theVarNormNLL=2*fr->minNll();
    theVarNormVal=((RooRealVar*)RooArgList(scanVars).at(0))->getVal();
    // save the point
    NLL.setVal(theVarNormNLL-mNLL);
    theDS->add(scanSet);
  }
  // first find nScanSegments for 1D
  Int_t nSegs=atoi(readConfStr("nScanSegments", "1", _runSec));
  if (nSegs<1) {
    cout<<" nScanSegments = "<<nSegs<<endl
	<<" reset to 1"<<endl;
    nSegs=1;
  }
  Int_t segIdx=_toyID%nSegs;
  for(Int_t i=0; i<nPoints; i++) {
    // restore params
    readFromStr(fullParams, paramSStr);
    // loop over scanVars to set initial values
    RooArgList scanList(scanVars);
    for(Int_t j=0; j<scanList.getSize(); j++) {
      RooRealVar *theVar=(RooRealVar*)scanList.at(j);
      Double_t min=theVar->getMin();
      Double_t max=theVar->getMax();
      if (scanList.getSize()>1) {
	theVar->setVal(min+(max-min)*RooRandom::randomGenerator()->Uniform());
      } else {
	// first find seg size
	Double_t segSize=(max-min)/nSegs;
	// reset the max and min
	min=min+segSize*segIdx;
	max=min+segSize;
	theVar->setVal(min+(max-min)*i/nPoints);
      }
      cout<<" Set scan var "<<theVar->GetName()<<" to "
	  <<theVar->getVal()<<endl;
    }
    RooFitResult *fr=_thePdf->fitTo(*scanPlotData, ConditionalObservables(_conditionalObs),
				    Save(scanPlotSave),Extended(scanPlotExtended), 
				    Verbose(scanPlotVerbose), Hesse(scanPlotHesse),
				    Minos(scanPlotMinos));
    //Int_t ncpus(1);
    //RooFitResult *fr=doTheFit(_thePdf, scanPlotData, fitOption, ncpus);

    if (fr->status()) {
      cout<<" Fit status for point #"<<i<<": "<<fr->status()<<endl;
      scanVars.Print("v");
      continue;
    }
    // shift fixed values
    scanVarShiftToNorm(scanVars, scanVarDiff);
    // save NLL
    NLL.setVal(2*fr->minNll()-mNLL);
    theDS->add(scanSet);
    if (curve) {
      RooRealVar *theVar=(RooRealVar*)RooArgList(scanVars).at(0);
      Double_t x=theVar->getVal();
      Double_t min=theVar->getMin();
      Double_t max=theVar->getMax();
      // first find seg size
      Double_t segSize=(max-min)/nSegs;
      // reset the max and min
      min=min+segSize*segIdx;
      max=min+segSize;
      curve->addPoint(x, NLL.getVal());
      if ((x<theVarNormVal)&&(theVarNormVal<x+(max-min)/nPoints)) {
	curve->addPoint(theVarNormVal, theVarNormNLL-mNLL);
      }
      if (NLL.getVal()>maxNLL) maxNLL=NLL.getVal();
    }
  }
  // save the TTree;
  TTree *theDSTree = createTreeFromDataset(theDS, kFALSE);
  plotList.Add(theDSTree);
  // set max for frame
  if (frame) frame->SetMaximum(maxNLL);
  
  if (scanVars.getSize()==1) {
    // do we have uncorrelated error?
    TString unCorrStr=readConfStr("scanUnCorrErr", "no", _runSec);
    if (!unCorrStr.BeginsWith("no")) {
      cout<<"scanUnCorrErr obsoleted, please use combine.C"<<endl
	  <<"remove the config to re-run the job"<<endl;
      exit(-1);
      Double_t ucErr=atof(unCorrStr);
      if (ucErr<=0) {
	cout<<" Uncorrelated error should be > 0 ("<<unCorrStr<<")"<<endl;
	exit(-1);
      }
      // create a new curve for NLL w/ uncorrelated syst error
      RooCurve *uccurve=frame->getCurve("NLL_curve");
      if (!uccurve) {
	cout<<" Can not find curve: NLL_curve"<<endl;
	exit(-1);
      }
      uccurve=(RooCurve*)uccurve->Clone();
      uccurve->SetName("NLL_curve_unCorr");
      frame->addPlotable(uccurve);
      addErrToCurve(uccurve, ucErr, &maxNLL);
      frame->SetMaximum(maxNLL);
    }
    
    // do we have correlated error?
    TString corrStr=readConfStr("scanCorrErr", "no", _runSec);
    if (!corrStr.BeginsWith("no")) {
      cout<<"scanCorrErr obsoleted, please use combine.C"<<endl
	  <<"remove the config to re-run the job"<<endl;
      exit(-1);
      Double_t cErr=atof(corrStr);
      if (cErr<=0) {
	cout<<" Correlated error should be > 0 ("<<corrStr<<")"<<endl;
	exit(-1);
      }
      // create a new curve for NLL w/ correlated syst error
      RooCurve *ccurve=frame->getCurve("NLL_curve_unCorr");
      if (!ccurve) ccurve=frame->getCurve("NLL_curve");
      if (!ccurve) {
	cout<<" Can not find curve named NLL_curve_unCorr or NLL_curve"<<endl;
	exit(-1);
      }
      ccurve=(RooCurve*)ccurve->Clone();
      ccurve->SetName("NLL_curve_Corr");
      frame->addPlotable(ccurve);
      addErrToCurve(ccurve, cErr, &maxNLL);
      frame->SetMaximum(maxNLL);
    }
  }
  
  // restore original params
  readFromStr(fullParams, paramSStr0);
  return frame;
}

/// \brief Add error to scan plot
/// \param curve Curve to change
/// \param errLo Low error to add
/// \param errHi High error to add
/// \param maxNLL The max NLL of the plot
///
/// It widens the NLL curve by low and high errors
/// to take systematic error into account
void rarMLFitter::addErrToCurve(RooCurve *curve, Double_t errLo,Double_t errHi,
				Double_t *maxNLL)
{
  if ((0==errLo)&&(0==errHi)) return;
  // first get the min point
  Double_t mean(0), yMin(0), err(errLo);
  rarNLL nll(curve);
  nll.getMin(mean, yMin);
  // shift yMin to zero
  shiftNLLCurve(curve, 0, -yMin);
  // add err onto curve
  Int_t nPoints=curve->GetN();
  for(Int_t i=0; i<nPoints; i++) {
    Double_t x, chi2Stat;
    curve->GetPoint(i, x, chi2Stat);
    err=(x<=mean) ? errLo : errHi;
    if (0==err) continue;
    Double_t chi2Syst=(x-mean)*(x-mean)/(err*err);
    if (chi2Stat+chi2Syst==0.) continue;
    Double_t chi2Tot = chi2Stat*chi2Syst/(chi2Stat+chi2Syst);
    curve->SetPoint(i, x, chi2Tot);
    if (maxNLL) if (*maxNLL<chi2Tot) *maxNLL=chi2Tot;
  }
}

/// \brief Add error to scan plot
/// \param curve Curve to change
/// \param err Error to add
/// \param maxNLL The max NLL of the plot
///
/// It widens the NLL curve by the error
/// to take systematic error into account
void rarMLFitter::addErrToCurve(RooCurve *curve, Double_t err,Double_t *maxNLL)
{
  addErrToCurve(curve, err, err, maxNLL);
}

/// \brief Shift the curve w/ specified shifts
/// \param curve The curve to shift
/// \param dx X shift
/// \param dy Y shift
/// \return The Curve shifted
///
/// It shifts the NLL curve with specified amount
RooCurve *rarMLFitter::shiftNLLCurve(RooCurve *curve, Double_t dx, Double_t dy)
{
  int nPoints=curve->GetN();
  for(Int_t i=0; i<nPoints; i++) {
    Double_t x(0), y(0);
    curve->GetPoint(i, x, y);
    curve->SetPoint(i, x+dx, y+dy);
  }
  return curve;
}

/// \brief Combine NLL curves together
/// \param curves List of NLL curve to combine
/// \param shiftToZero To shift new curve to zero or not
/// \param maxNLL The max NLL of the plot
/// \return Curve combined
///
/// It combines input NLL curves
/// \todo Combine 2D plots (and more-D if possible)
RooCurve *rarMLFitter::combineNLLCurves(TList &curves, Bool_t shiftToZero,
					Double_t *maxNLL)
{
  // first find all x values
  set<Double_t> xs;
  TList nlls;
  for (Int_t i=0; i<curves.GetSize(); i++) {
    RooCurve *curve=(RooCurve*)curves.At(i);
    rarNLL *nll=new rarNLL(curve);
    nlls.Add(nll);
    Int_t nPoints=curve->GetN();
    for (Int_t j=0; j<nPoints; j++) {
      Double_t x, y;
      curve->GetPoint(j, x, y);
      xs.insert(x);
    }
  }
  //nlls.Print();
  // Max y
  Double_t yMax(0);
  // build the combined curve
  RooCurve *theCurve=new RooCurve();
  set<Double_t>::iterator the_iterator=xs.begin();
  while (the_iterator!=xs.end()) {
    Double_t y=0;
    for (Int_t i=0; i<curves.GetSize(); i++) {
      rarNLL *nll=(rarNLL*)nlls.At(i);
      y+=nll->getNLL(*the_iterator);
    }
    if (yMax<y) yMax=y;
    theCurve->addPoint(*the_iterator, y);
    the_iterator++;
  }
  // shift min back to zero
  if (shiftToZero) {
    // find the new min
    rarNLL nll(theCurve);
    Double_t xMin(0), yMin(0);
    nll.getMin(xMin, yMin);
    if (yMin<-1e-05) { // should be >= 0
      cout<<" Min of combined NLL is "<<yMin<<endl
          <<" Should be equal to or greater than 0"<<endl;
    }
    // shift it back to zero
    shiftNLLCurve(theCurve, 0, -yMin);
    yMax-=yMin;
  }
  // return the max y
  if (maxNLL) if (*maxNLL<yMax) *maxNLL=yMax;
  
  return theCurve;
}

/// \brief Combine NLL curves with correlated errors
/// \param curves List of NLL curve to combine
/// \param errs Array of correlated errors
/// \param maxNLL The max NLL of the plot
/// \return Curve combined
///
/// It combines input NLL curves and also `add' correlated errors properly
RooCurve *rarMLFitter::combineNLLCurves(TList &curves, Double_t errs[],
                                        Double_t *maxNLL)
{
  // First build the combined curve
  RooCurve *theCurve=combineNLLCurves(curves, kTRUE, maxNLL);
  // Then add correlated errors
  TList lCurves, rCurves; // left and right shifted curves by errors
  for(Int_t i=0; i<curves.GetSize(); i++) {
    RooCurve *curve=(RooCurve*)curves.At(i);
    RooCurve *lCurve=(RooCurve*)curve->Clone();
    RooCurve *rCurve=(RooCurve*)curve->Clone();
    shiftNLLCurve(lCurve, -fabs(errs[i]), 0);
    shiftNLLCurve(rCurve, +fabs(errs[i]), 0);
    lCurves.Add(lCurve);
    rCurves.Add(rCurve);
  }
  RooCurve *lCurve=combineNLLCurves(lCurves);
  RooCurve *rCurve=combineNLLCurves(rCurves);
  Double_t lx(0), rx(0), yMin(0);
  rarNLL lnll(lCurve), rnll(rCurve);
  lnll.getMin(lx, yMin);
  rnll.getMin(rx, yMin);
  // the overall correlated error
  Double_t err=fabs(rx-lx)/2.;
  // add the overall correlated error to plot
  addErrToCurve(theCurve, err, maxNLL);
  
  return theCurve;
}

/// \brief Get mean and error associated with NLL curve
/// \param curve NLL curve
/// \param errLo Low error
/// \param errHi High error
/// \return Mean with the NLL curve
Double_t rarMLFitter::getMeanErrs(RooCurve *curve,
				  Double_t *errLo, Double_t *errHi)
{
  Double_t theX(0), theY(0);
  rarNLL theNLL(curve);
  theNLL.getMin(theX, theY);
  if (fabs(theY)>1e-5) {
    cout<<" The NLL curve is minimized @ y="<<theY<<endl
	<<" (should be y=0)"<<endl;
  }
  TArrayD errs=theNLL.getX(1+theY);
  if (errs.GetSize()!=2) {
    cout<<" W A R N I N G ! !"<<endl
	<<" The number of points for error from NLL curve is not TWO"<<endl
	<<" The first and last points (if any)"
	<<" will be used for error calculation"<<endl
	<<" Your results might be problematic"<<endl
	<<" The points found:"<<endl;
    for (Int_t i=0; i<errs.GetSize(); i++) {
      cout<<" "<<errs[i];
    }
    cout<<endl;
  }
  if (errs.GetSize()>0) {
    errs[0]-=theX;
    errs[errs.GetSize()-1]-=theX;
    if (errLo) *errLo=errs[0];
    if (errHi) *errHi=errs[errs.GetSize()-1];
  }
  return theX;
}

/// \brief Get significance from NLL curve
/// \param curve NLL curve
/// \param refVal 0 significance reference point
/// \return significance
Double_t rarMLFitter::getSignf(RooCurve *curve, Double_t refVal)
{
  rarNLL nll(curve);
  Double_t signf=nll.getNLL(refVal);
  if (signf<0) {
    cout<<" NLL("<<refVal<<")="<<signf<<"<0"<<endl;
    signf=0;
  }
  signf=sqrt(signf);
  //cout<<" Signf = "<<signf<<" ("<<curve->GetName()<<" @ "<<refVal<<")"<<endl;
  return signf;
}

/// \brief Get UL from NLL curve
/// \param curve NLL curve
/// \param CL Confidence limit (default .90)
/// \return UL
Double_t rarMLFitter::getUL(RooCurve *curve, Double_t CL)
{
  // first instantiate rarNLL object
  rarNLL theNLL(curve);
  Double_t theIntegral=theNLL.getLIntegral();
  Double_t UL=theNLL.getLIntegralInverse(theIntegral*CL);
  return UL;
}

/// \brief Get contour plot
/// \param plotList Plot list
/// \return The last plot created
///
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_contourPlot">See doc for contour plot action.</a>
///
/// The function reads in the config options from action section.
/// It then draws a 2D contour plot for specified floating parameters
/// with number of contour given by config \p nContours in action section,
/// and finally it returns the frame containing the plot.
RooPlot *rarMLFitter::doContourPlot(TList &plotList)
{
  cout<<endl<<" In rarMLFitter doContourPlot for "<<GetName()<<endl;
  RooPlot *frame(0);
  
  // get contour plot data
  RooDataSet *contourPlotData=
    _datasets->getData(readConfStr("contourPlotData", "notSet", _runSec));
  if(!contourPlotData) {
    cout<<endl<<" Please specify in section ["<<_runSec<<"]"<<endl
	<<" the dataset you want to do contourPlot with config "<<endl<<endl
	<<"   contourPlotData = <datasetName>"<<endl
	<<endl<<" and ru-run your job"<<endl;
    exit(-1);
  }
  chkBlind(contourPlotData->GetName());
  cout<<" Using dataset "<<contourPlotData->GetName()<<endl;
  
  // full paramters
  RooArgSet fullParams(*_thePdf->getParameters(contourPlotData));
  // save them first
  string paramSSaver;
  writeToStr(fullParams, paramSSaver);
  
  // get contour plot vars
  rarStrParser contourVarsParser=readConfStr("contourVars", "notSet", _runSec);
  RooRealVar *contX(0), *contY(0);
  Double_t xMin(0), xMax(0), yMin(0), yMax(0);
  // for x
  contX=(RooRealVar*)fullParams.find(contourVarsParser[0]);
  contourVarsParser.Remove();
  if (!contX) {
    cout<<"Can not find parameter "<<contourVarsParser[0]<<endl;
    exit(-1);
  }
  xMin=contX->getMin();
  xMax=contX->getMax();
  // check if the next two are limits for it
  if (isNumber(contourVarsParser[0])) {
    xMin=atof(contourVarsParser[0]);
    contourVarsParser.Remove();
  }
  if (isNumber(contourVarsParser[0])) {
    xMax=atof(contourVarsParser[0]);
    contourVarsParser.Remove();
  }
  if (xMin>xMax) {
    Double_t v=xMin;
    xMin=xMax;
    xMax=v;
  }
  contX->setRange(xMin,xMax);
  // for y
  contY=(RooRealVar*)fullParams.find(contourVarsParser[0]);
  contourVarsParser.Remove();
  if (!contY) {
    cout<<"Can not find parameter "<<contourVarsParser[0]<<endl;
    exit(-1);
  }
  yMin=contY->getMin();
  yMax=contY->getMax();
  // check if the next two are limits for it
  if (isNumber(contourVarsParser[0])) {
    yMin=atof(contourVarsParser[0]);
    contourVarsParser.Remove();
  }
  if (isNumber(contourVarsParser[0])) {
    yMax=atof(contourVarsParser[0]);
    contourVarsParser.Remove();
  }
  if (yMin>yMax) {
    Double_t v=yMin;
    yMin=yMax;
    yMax=v;
  }
  contY->setRange(yMin,yMax);
  // change floating ranges for free parameters
  rarStrParser restParser=
    readConfStr("contourRestrictFloatParams", "no", _runSec);
  if ("no"!=restParser[0]) {
    Double_t dScale=2;
    if (isNumber(restParser[0])) {
      dScale=atof(restParser[0]);
      restParser.Remove();
    }
    // go through all parameters
    RooArgList fullParamList(fullParams);
    for (Int_t i=0; i<fullParamList.getSize(); i++) {
      RooRealVar *theVar=(RooRealVar*)fullParamList.at(i);
      if ((contX==theVar)||(contY==theVar)) continue;
      if (TString("RooRealVar")!=theVar->ClassName()) continue;
      if (theVar->isConstant()) continue;
      if ((!theVar->hasError()) && (!theVar->hasAsymError())) continue;
      Double_t min=theVar->getMin();
      Double_t max=theVar->getMax();
      Double_t myScale=dScale;
      Int_t myIdx=restParser.Index(theVar->GetName());
      if ((myIdx>=0)&&isNumber(restParser[myIdx+1]))
	myScale=atof(restParser[myIdx+1]);
      Double_t nMin(min), nMax(max);
      if (theVar->hasError()&&(theVar->getError()>0)) {
	nMin=theVar->getVal()-theVar->getError()*myScale;
	nMax=theVar->getVal()+theVar->getError()*myScale;
      }
      if (theVar->hasAsymError()) {
	if (theVar->getAsymErrorLo()<0)
	  nMin=theVar->getVal()+theVar->getAsymErrorLo()*myScale;
	if (theVar->getAsymErrorHi()>0)
	  nMax=theVar->getVal()+theVar->getAsymErrorHi()*myScale;
      }
      if (nMin>nMax) {
	Double_t v=nMin;
	nMin=nMax;
	nMax=v;
      }
      if (nMin>min) min=nMin;
      if (nMax<max) max=nMax;
      theVar->setRange(min, max);
      cout<<" Reset limits: "
	  <<theVar->GetName()<<" ("<<min<<", "<<max<<")"<<endl;
    }
  }

  // number of contour
  Int_t nContours=atoi(readConfStr("nContours", "2", _runSec));
  if (nContours<1) nContours=1;
  if (nContours>6) nContours=6;
  cout<<" Drawing "<<nContours<<" contour(s) of "
      <<contX->GetName()<<" vs "<<contY->GetName()<<endl<<endl;
  RooNLLVar nll("nll","nll",*_thePdf, *contourPlotData, kTRUE);
  rarMinuit min(nll);
  if      (nContours<=1) frame = min.contour(*contX, *contY, 1, 0);
  else if (2==nContours) frame = min.contour(*contX, *contY, 1, 2);
  else if (3==nContours) frame = min.contour(*contX, *contY, 1, 2, 3);
  else if (4==nContours) frame = min.contour(*contX, *contY, 1, 2, 3, 4);
  else if (5==nContours) frame = min.contour(*contX, *contY, 1, 2, 3, 4, 5);
  else if (6<=nContours) frame = min.contour(*contX, *contY, 1, 2, 3, 4, 5, 6);
  TString myName=Form("contour_%s_%s", contX->GetName(), contY->GetName());
  if (frame) {
    frame->SetNameTitle(myName, myName);
    plotList.Add(frame);
  }
  
  // restore params
  readFromStr(fullParams, paramSSaver);
  return frame;
}

/// \brief Get combine plot
/// \param plotList Plot list
/// \return The plot, if any, created
///
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_combinePlot">See doc for combine plot action.</a>
///
/// The function combines NLL curves creatted by scan action and convolutes
/// them with systematic errors to get Central Values, Upper limits and/or significance
RooPlot *rarMLFitter::doCombinePlot(TList &plotList)
{
  cout<<endl<<" In rarMLFitter doCombinePlot for "<<GetName()<<endl;
  RooPlot *frame(0);

  // number of curves to combine
  Int_t ncurves = atoi(readConfStr("combineNcurves", "0", _runSec));

  if (ncurves <= 0) {
    cout << "combineNcurves must be set to 1 or more curves." << endl;
    return (0);
  }

  //cout << "Ncurves : " << ncurves << endl;

  vector<Double_t> addSystErrLo, addSystErrHi, fitBias;
  vector<Double_t> uncorrSystErrLo, uncorrSystErrHi, corrSystErr;
  vector<TString> filenames, plotnames;

  // set default for systematic errors
  Bool_t doUL = kTRUE;
  Double_t ConfidenceLevel(0.9);
  for (Int_t i=0; i < ncurves ; i++) {
    addSystErrLo.push_back(0.0);  addSystErrHi.push_back(0.0);
    uncorrSystErrLo.push_back(0.0); uncorrSystErrHi.push_back(0.0);
    corrSystErr.push_back(0.0);
    fitBias.push_back(0.0);
  }

  TString tempstr("");
  // read in file names
  tempstr = readConfStr("combineFilenames", "notSet", _runSec);
  if ("notSet" == tempstr) {
    cout << "combineFilenames must be specified." << endl;
    return(0);
  } else {
    rarStrParser rarFilenames = tempstr;
    if (rarFilenames.nArgs() < ncurves) {
      cout << "combineFilenames: Expected " << ncurves << " filenames, found " 
	   << rarFilenames.nArgs() << "." << endl;
      return(0);
    }
    for (Int_t i=0; i < rarFilenames.nArgs() ; i++) {
      filenames.push_back(rarFilenames[i]);
      //cout << filenames[i] << endl;
    }
  }

  // read in RooPlot names from Scan action
  tempstr = readConfStr("combinePlotnames", "notSet", _runSec);
  if ("notSet" == tempstr) {
    cout << "combinePlotnames must be specified." << endl;
    return(0);
  } else {
    rarStrParser rarPlotnames = tempstr;
    if (rarPlotnames.nArgs() < ncurves) {
      cout << "combinePlotnames: Expected " << ncurves << " RooPlot scan plot names, found " 
	   << rarPlotnames.nArgs() << "." << endl;
      return(0);
    }
    for (Int_t i=0; i < rarPlotnames.nArgs() ; i++) {
      plotnames.push_back(rarPlotnames[i]);
      //  cout << plotnames[i] << endl;
    }
  }

  // read in Fit bias
  tempstr = readConfStr("combineFitBias", "notSet", _runSec);
  if ("notSet" != tempstr) {
    rarStrParser rarFitBias = tempstr;
    if (rarFitBias.nArgs() < ncurves) {
      cout << "combineFitBias: Expected " << ncurves << " fit biases, found " 
	   << rarFitBias.nArgs() << "." << endl;
      return(0);
    }
    for (Int_t i=0; i < rarFitBias.nArgs() ; i++) {
      fitBias[i]    = atof(rarFitBias[i]);
    }
  }
  // read in Additive systematic errors. If nArgs = ncurves, assume symmetric errors
  // if nArgs = 2*ncurves assume asymmetric (-ve/+ve).
  tempstr = readConfStr("combineAdditive", "notSet", _runSec);
  if ("notSet" != tempstr) {
    rarStrParser rarAddErrs = tempstr;
    Bool_t sym   = (rarAddErrs.nArgs() == ncurves);
    Bool_t assym = (rarAddErrs.nArgs() == (2*ncurves));

    if (!sym && !assym) {
      cout << "combineAdditive: Expected either " << ncurves << " (symmetric) errors or " 
	   << ncurves*2 << " (asymmetric) errors. Found " << rarAddErrs.nArgs() << "." << endl;
      return (0);
    }
    for (Int_t i=0; i < ncurves ; i++) {
      if (assym) {
	addSystErrLo[i] = atof(rarAddErrs[2*i]);
	addSystErrHi[i] = atof(rarAddErrs[2*i+1]);
      } else {
	addSystErrHi[i] = atof(rarAddErrs[i]);
	addSystErrLo[i] = -1 * addSystErrHi[i];
      }
    }
    
  }

  // read in Uncorrelated systematic errors. 
  tempstr = readConfStr("combineMultiplicativeUncorrelated", "notSet", _runSec);
  if ("notSet" != tempstr) {
    rarStrParser rarSystErrs = tempstr;
    Bool_t sym   = (rarSystErrs.nArgs() == ncurves);
    Bool_t assym = (rarSystErrs.nArgs() == (2*ncurves));

    if (!sym && !assym) {
      cout << "combineMultiplicativeUncorrelated: Expected either " 
	   << ncurves << " (symmetric) errors or " 
	   << ncurves*2 << " (asymmetric) errors. Found "
	   << rarSystErrs.nArgs() << "." << endl;
      return (0);
    }
    for (Int_t i=0; i < ncurves ; i++) {
      //     cout << "ASSYM : " << i << " " << rarSystErrs[2*i] << endl;
      if (assym) {
	//cout << rarSystErrs[2*i] << " " << rarSystErrs[2*i+1] << endl;
	uncorrSystErrLo[i] = atof(rarSystErrs[2*i]);
	uncorrSystErrHi[i] = atof(rarSystErrs[2*i+1]);
      } else {
	uncorrSystErrHi[i] = atof(rarSystErrs[i]);
	uncorrSystErrLo[i] = -1 * uncorrSystErrHi[i];
      }
    }
  }

  // read in correlated systematic errors. 
  tempstr = readConfStr("combineMultiplicativeCorrelated", "notSet", _runSec);
  if ("notSet" != tempstr) {
    rarStrParser rarSystMultErrs = tempstr;
    Bool_t sym   = (rarSystMultErrs.nArgs() == ncurves);

    if (!sym) {
      cout << "combineMultiplicativeCorrelated: Expected either " << ncurves 
	   << " errors. Found " << rarSystMultErrs.nArgs() << "." << endl;
      return (0);
    }
    for (Int_t i=0; i < rarSystMultErrs.nArgs() ; i++) {
      corrSystErr[i] = atof(rarSystMultErrs[i]);
    }
  }

  // read in upper limit test
  tempstr = readConfStr("combineUpperLimit", "notSet", _runSec);
  if ("notSet" != tempstr) {
    rarStrParser rarULtest = tempstr;
    if (rarULtest.nArgs() != 2) {
      cout << "combineUpperLimit: Expected <yes/no> <CL>." << endl;
      return (0);
    }
    if (rarULtest[0] != "yes") {doUL = kFALSE;}
    ConfidenceLevel = atof(rarULtest[1]);
    if (ConfidenceLevel <= 0.0 || ConfidenceLevel >= 100) {
      cout << "combineUpperLimit: Expected CL between 0 and 100%. Found "
	   <<  ConfidenceLevel << "." << endl;
      return (0);
    }
    ConfidenceLevel *=  0.01; // fraction rather than %
  }

  // read in upper limit test
  TString xAxisTitle = "";
  tempstr = readConfStr("combineXaxisTitle", "notSet", _runSec);
  rarStrParser xtitle = tempstr;
  if ("notSet" != tempstr) {xAxisTitle = xtitle[0];}

  Double_t uncorr_var(0);

  for (Int_t j=0; j < ncurves ; j++) {
    cout << endl << "Curve/mode: " << j << endl;
    cout << "ScanPlot File: " << filenames[j] << endl;
    cout << "ScanPlot Name: " << plotnames[j] << endl;
    cout << "Additive uncorrelated Error         : " << addSystErrLo[j] << "/" << addSystErrHi[j] << endl;
    cout << "Multiplicative correlated Error     : +/- " << corrSystErr[j] << endl;
    cout << "Multiplicative uncorrelated Error   : " << uncorrSystErrLo[j] 
	 << "/" << uncorrSystErrHi[j] << endl;

    // Add uncorrelated errors together
    uncorr_var = pow(uncorrSystErrHi[j],2) + pow(addSystErrHi[j],2);
    uncorrSystErrHi[j] =  sqrt(uncorr_var);

    uncorr_var = pow(uncorrSystErrLo[j],2) + pow(addSystErrLo[j],2);
    uncorrSystErrLo[j] =  -1 * sqrt(uncorr_var);

    cout << "Total Uncorrelated (add+Mult) Error : " << uncorrSystErrLo[j] 
	 << "/" << uncorrSystErrHi[j] << endl;

    Double_t tothi = pow(uncorrSystErrHi[j],2) + pow(corrSystErr[j],2);
    Double_t totlo = pow(uncorrSystErrLo[j],2) + pow(corrSystErr[j],2);
    cout << "Total Error                         : " << -1 * sqrt(totlo) << "/" << sqrt(tothi) << endl;
  }
  cout << endl;

  // now create the plots for each curved separately showing the total and stat. only
  // NLL distributions
  Bool_t doSignf = kTRUE;
  vector<Double_t> addLo, addHi, fb, unSystLo, unSystHi, corrSyst;
  vector<TString> fn, pn;
  RooPlot *thePlot(0);
  for (Int_t j=0; j < ncurves ; j++) {
    // reset the vectors
    fn.clear(); pn.clear(); addLo.clear(); addHi.clear();
    unSystLo.clear();unSystHi.clear(); corrSyst.clear(); fb.clear();
    // copy to a one element vector
    fn.push_back(filenames[j]); pn.push_back(plotnames[j]);
    addLo.push_back(addSystErrLo[j]); addHi.push_back(addSystErrHi[j]); 
    unSystLo.push_back(uncorrSystErrLo[j]); unSystHi.push_back(uncorrSystErrHi[j]);
    corrSyst.push_back(corrSystErr[j]); fb.push_back(fitBias[j]);

    // create the plot
    cout << endl << "Curve/mode: " << j << endl;
    thePlot = combine(1, fn, pn, fb, addLo, addHi, unSystLo, unSystHi,
		      corrSyst, xAxisTitle, doSignf, doUL, ConfidenceLevel);
    // update title
    thePlot->SetNameTitle(Form("combinePlot_Mode%d", j),"NLL with syst+stat and stat. errs only");
    plotList.Add(thePlot);
  }

  // Now combine the curves
  cout << endl << "Combined curves: " << endl;
  thePlot = combine(ncurves, filenames, plotnames, fitBias,
			     addSystErrLo, addSystErrHi,
			     uncorrSystErrLo, uncorrSystErrHi, corrSystErr, 
			     xAxisTitle, doSignf, doUL, ConfidenceLevel);

  thePlot->SetNameTitle("combinePlot","Combined curves");
  cout << endl;
  plotList.Add(thePlot);
  
  return (frame);
}

// combine NLL curves (w sys errs), draw the plots,
// calculate signf. and UL based on the final curve
RooPlot *rarMLFitter::combine(Int_t nModes, 
			      const vector<TString> fileNames,
			      const vector<TString> plotNames,
			      const vector<Double_t> fitBias,
			      const vector<Double_t> addSystErrLo,
			      const vector<Double_t> addSystErrHi,
			      const vector<Double_t> uncorrSystErrLo,
			      const vector<Double_t> uncorrSystErrHi,
			      const vector<Double_t> corrSystErr,
			      const TString xAxisTitle,
			      Bool_t doSignf, Bool_t doUL, Double_t CL)
{
  RooPlot *thePlot(0);
  // List of RooCurves
  TList curves; // original individual NLL curves
  TList curvesWadds; // individual NLL curves with additive errors
  TList curvesWerrs; // individual NLL curves with uncorrelated errors

  // RooPlot properties
  TString yAxisTitle="-2 ln (L/L_{0})";
  Double_t xMin(0), xMax(0), yMax(0);
  // read in all the NLL Curves

  for (Int_t i=0; i<nModes; i++) {
    TFile *f=new TFile(fileNames[i]);
    if (!f) {
      cout<<" Can not read from root file: "<<fileNames[i]<<endl;
      return thePlot;
    }
    RooPlot *scanPlot=(RooPlot*)f->Get(plotNames[i]);
    if (!scanPlot) {
      cout<<" Can not read in RooPlot "<<plotNames[i]
	  <<" from "<<fileNames[i]<<endl;
      return thePlot;
    }
    Double_t theYmax=scanPlot->GetMaximum();
    if (yMax<theYmax) yMax=theYmax;
    TAxis *a=scanPlot->GetXaxis();
    Double_t theXmin=a->GetXmin();
    if (xMin>theXmin) xMin=theXmin;
    Double_t theXmax=a->GetXmax();
    if (xMax<theXmax) xMax=theXmax;
    //if (!xAxisTitle) xAxisTitle=(Char_t*)a->GetTitle();

    // Get curve
    RooCurve *nllCurve=scanPlot->getCurve("NLL_curve");
    if (!nllCurve) {
      cout<<" Can not read NLL curve NLL_curve from RooPlot "<<plotNames[i]
	  <<" in root file "<<fileNames[i]<<endl;
      return thePlot;
    }
    // rename it so one can tell which is which
    nllCurve->SetNameTitle(Form("NLL_curve_Mode%d", i),
                           Form("curve for sub-mode %d", i));
    { // now shift the curve to min=0
      Double_t mean(0), yMin(0);
      rarNLL nll(nllCurve);
      nll.getMin(mean, yMin);
      shiftNLLCurve(nllCurve, -fitBias[i], -yMin);
    }
    // Clone it
    RooCurve *nllCurveWadd=(RooCurve*)nllCurve->Clone();
    RooCurve *nllCurveWerr=(RooCurve*)nllCurve->Clone();
    // `add' errors onto the curve
    addErrToCurve(nllCurveWadd, addSystErrLo[i], addSystErrHi[i]);
    addErrToCurve(nllCurveWerr, uncorrSystErrLo[i], uncorrSystErrHi[i]);
    // add the curve into lists
    curves.Add(nllCurve);
    curvesWadds.Add(nllCurveWadd);
    curvesWerrs.Add(nllCurveWerr);
  }
  //curves.Print();
  //curvesWadds.Print();
  //curvesWerrs.Print();
  
  // combine the curve together
  RooCurve *tCurve=combineNLLCurves(curves, kTRUE, &yMax);
  tCurve->SetNameTitle("NLL_curve_total", "total curve w/o syst errors");
  // combine the curve (w additive errors) together
  RooCurve *aCurve=combineNLLCurves(curvesWadds, kTRUE, &yMax);
  aCurve->SetNameTitle("NLL_curve_additive", "total curve w/ additve errors");
  // combine the curve (w uncorrelated errors) together
  RooCurve *uCurve=combineNLLCurves(curvesWerrs, kTRUE, &yMax);
  uCurve->SetNameTitle("NLL_curve_unCorr", "total curve w/ unCorr. errors");
  // combine the curves and add correlated errors
  // silly hack
  Double_t corrSystErrArray[20];
  for (Int_t i=0; i< (Int_t) corrSystErr.size(); i++) {corrSystErrArray[i] = corrSystErr[i];}
  RooCurve *cCurve=combineNLLCurves(curvesWerrs, corrSystErrArray, &yMax);
  cCurve->SetNameTitle("NLL_curve_Corr", "total curve w/ ALL syst errors");
  // Draw the RooPlot
  thePlot=new RooPlot(xMin, xMax, 0, yMax);
  //  thePlot->SetNameTitle(plotNames[0], plotNames[0]);
  thePlot->SetNameTitle("combinePlot","Combined curves");
  thePlot->SetYTitle(yAxisTitle);
  thePlot->SetXTitle(xAxisTitle);
  // total curve without syst. error
  thePlot->addPlotable(tCurve);
  tCurve->SetLineWidth(2);
  tCurve->SetLineStyle(kDashed);
  thePlot->setInvisible(tCurve->GetName()); // invisible
  // total curve w/ additive errors
  thePlot->addPlotable(aCurve);
  aCurve->SetLineWidth(2);
  aCurve->SetLineStyle(kDashDotted);
  aCurve->SetLineColor(kYellow);
  thePlot->setInvisible(aCurve->GetName()); // invisible
  // total curve w/ uncorrelated errors
  thePlot->addPlotable(uCurve);
  uCurve->SetLineWidth(2);
  uCurve->SetLineStyle(kDashDotted);
  uCurve->SetLineColor(kGreen);
  thePlot->setInvisible(uCurve->GetName()); // invisible
  // total curve w/ all errors
  thePlot->addPlotable(cCurve);
  cCurve->SetLineWidth(2);
  // add individual NLL curves
  for (Int_t i=0; i<nModes; i++) {
    RooCurve *theCCurve=(RooCurve*)curves.At(i)->Clone();
    thePlot->addPlotable(theCCurve);
    theCCurve->SetLineWidth(1);
    theCCurve->SetLineColor(rarBasePdf::getColor(i));
    theCCurve->SetLineStyle(0==(i%2)?2:4); // kDashed or KDashDotted
  }
  // draw it
  thePlot->Draw();
  thePlot->Print("V");
  cout<<endl
      <<" To show hidden curves (those with Option \":I\"),"
      <<" type, for example"<<endl
      <<"   NLL_curve_total->SetDrawOption(\"\")"
      <<endl
      <<" To hide a curve, type, for example"<<endl
      <<"   NLL_curve_Mode0->SetDrawOption(\"I\")"
      <<endl
      <<" Then redraw the plot: gPad->Modified()"
      <<endl;
  // get new mean, stat. error, syst. errors
  Double_t theMean(0), systErrLo(0), systErrHi(0), statErrLo(0), statErrHi(0);
  // find stat. errors
  getMeanErrs(tCurve, &statErrLo, &statErrHi);
  // find final mean and syst. errors
  theMean=getMeanErrs(cCurve, &systErrLo, &systErrHi);
  // get pure syst. error
  systErrLo=(systErrLo*systErrLo-statErrLo*statErrLo);
  if (systErrLo>=0) systErrLo=-sqrt(systErrLo);
  else {
    cout<<" systErrLo^2="<<systErrLo<<endl
	<<" Can not be true!!!! Please check!!!"<<endl;
  }
  systErrHi=(systErrHi*systErrHi-statErrHi*statErrHi);
  if (systErrHi>=0) systErrHi=sqrt(systErrHi);
  else {
    cout<<" systErrHi^2="<<systErrHi<<endl
	<<" Can not be true!!!! Please check!!!"<<endl;
  }
  // signf and UL
  Double_t signfStat(0), ULStat(0), signf(0), UL(0);
  if (doSignf) {
    signfStat=getSignf(tCurve);
    signf=getSignf(aCurve);
  }
  if (doUL) {
    ULStat=getUL(tCurve, CL);
    UL=getUL(cCurve, CL);
  }
  
  // final output for combined results
  cout<<endl;
  cout<<"Combined results based on "<<cCurve->GetName()
      <<" (and "<<tCurve->GetName()<<" for stat.):"<<endl
      <<" mean = "<<theMean
      <<" ("<<statErrLo<<", +"<<statErrHi<<")(stat)"
      <<" ("<<systErrLo<<", +"<<systErrHi<<")(syst)"
      <<endl;
  if (doSignf) cout<<" Signf = "<<signf<<" (wrt 0)"
		   <<" (stat only = "<<signfStat<<")"<<endl;
  if (doUL) cout<<" UL = "<<UL<<" (CL="<<CL<<")"
		<<" (stat only = "<<ULStat<<")"<<endl;
  
  return thePlot;
}

/// \brief Get sPlot
/// \param theVar obs to plot
/// \param plotList Plot list
/// \return The last plot created
///
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_sPlot">See doc for sPlot action.</a>
///
/// The function reads in the config options from action section.
/// It then draws sPlots for given components by calling #rarSPlot.
RooPlot *rarMLFitter::doSPlot(RooRealVar *theVar, TList &plotList)
{
  cout<<endl<<" In rarMLFitter doSPlot of "<<theVar->GetName()
      <<" for "<<GetName()<<endl;
  RooPlot *frame(0);
  // theVar name
  TString varName=theVar->GetName();
  // get sPlot data
  TString spdStr=readConfStrCnA("sPlotData_"+varName,"notSet");
  if ("notSet"==spdStr) spdStr=readConfStrCnA("sPlotData", "");
  RooDataSet *sPlotData(0);
  if (""!=spdStr) sPlotData=_datasets->getData(spdStr);
  if(!sPlotData) {
    cout<<endl<<" Please specify in section ["<<_runSec<<"]"<<endl
	<<" the dataset you want to project with config "<<endl<<endl
	<<"   sPlotData = <datasetName>  // and/or (see online doc)"<<endl
	<<"   sPlotData_<obsName> = <datasetName>"<<endl
	<<endl<<" and ru-run your job"<<endl;
    exit(-1);
  }
  chkBlind(sPlotData->GetName());
  cout<<" Using dataset "<<sPlotData->GetName()<<endl;
  
  // get plot bins
  Int_t nBins=atoi(readConfStr(Form("plotBins_%s", theVar->GetName()),
			       "0", _runSec));
  if (nBins<=0) nBins=theVar->getBins();
  // plot range
  Double_t plotMin=theVar->getMin();
  Double_t plotMax=theVar->getMax();
  getRange(theVar, "plotRange_", plotMin, plotMax,  _runSec);
  // construct ignored obs
  RooArgSet ignoredObsSet(*theVar);
  // any more?
  TString ignStr=readConfStr("sPlotIgnoredVars_"+varName, "notSet", _runSec);
  if ("notSet"==ignStr) ignStr=readConfStr("sPlotIgnoredVars", "", _runSec);
  rarStrParser ignoredParser=ignStr;
  for (Int_t i=0; i<ignoredParser.nArgs(); i++) {
    RooRealVar *ignored=(RooRealVar*)_fullObs->find(ignoredParser[i]);
    if (ignored) ignoredObsSet.add(*ignored);
  }
  cout<<" Observables not in sPlot PDF:"<<endl;
  ignoredObsSet.Print();
  // first get pdfs w/o ignored
  RooArgList thePdfsWOvar=getPdfsWOvar(ignoredObsSet);
  //RooArgList thePdfsWOvar=getPdfsWOvar(theVar);
  if (thePdfsWOvar.getSize()<1) {
    cout<<" Can not build pdf for sPlot"<<endl;
    exit(-1);
  }
  // do we use simultaneousFit?
  Bool_t simFit=getControlBit("SimFit");
  // no simPdf for multiple physCats for now
  if ((thePdfsWOvar.getSize()>1) && simFit) {
    cout<<" Currently no sPlot support for multi physCat in RooRarFit"<<endl;
    exit(-1);
  }
  // build Pdf for sPlot
  RooAbsPdf *theFullPdfWOvar=(RooAbsPdf*)thePdfsWOvar.at(0);
  // do we have simultaneousFit
  if (simFit) { // build SimPdf for sPlot
    RooSimPdfBuilder *simBuilder=new RooSimPdfBuilder(thePdfsWOvar);
    RooArgSet *simConfig=simBuilder->createProtoBuildConfig();
    // get configs from #_simConfig
    ((RooStringVar*)simConfig->find("splitCats"))->
      setVal(((RooStringVar*)_simConfig->find("splitCats"))->getVal());
    TString modelStr=((RooStringVar*)_simConfig->find("physModels"))->getVal();
    for (Int_t i=0; i<thePdfsWOvar.getSize(); i++) {
      TString pdfNameWOvar= thePdfsWOvar[i].GetName();
      TString pdfNameWvar=_subPdfs[i].GetName();
      ((RooStringVar*)simConfig->find(pdfNameWOvar))->
	setVal(((RooStringVar*)_simConfig->find(pdfNameWvar))->getVal());
      modelStr.ReplaceAll(pdfNameWvar, pdfNameWOvar);
    }
    ((RooStringVar*)simConfig->find("physModels"))->setVal(modelStr);
    cout<<"simPdfBuilder configs for sPlot:"<<endl;
    simConfig->Print("v");
    simBuilder->addSpecializations(_simBuilder->splitLeafList());
    theFullPdfWOvar=(RooAbsPdf*)
      simBuilder->buildPdf(*simConfig, sPlotData, 0, kTRUE);
    //simBuilder->splitLeafList().Print();
  }
  theFullPdfWOvar->Print("v");
  // get full paramters
  RooArgSet fullParams(*theFullPdfWOvar->getParameters(sPlotData));
  // save them first
  string paramSSaver;
  writeToStr(fullParams, paramSSaver);
  // get free yield, pdfsWvar, and pdfWOvar for each component
  RooArgList yields, pdfsWvar, pdfsWOvar, yield0s, pdf0sWOvar;
  TList rarPdfsWvar;
  Int_t nModel=1;
  if (simFit) nModel=_nComp;
  for (Int_t i=0; i<nModel; i++) {
    rarMLPdf *thisModel=(rarMLPdf*)_pdfList.At(i);
    TList *modelPdfList=thisModel->getPdfList();
    RooArgList modelCoeffList(thisModel->getSCoeffList());
    RooAddPdf *thisModelWOvar=(RooAddPdf*)thePdfsWOvar.at(i);
    RooArgList thisModelWOvarList(thisModelWOvar->pdfList());
    for (Int_t j=0; j<modelCoeffList.getSize(); j++) {
      // get yield
      if (TString("RooRealVar")!=modelCoeffList[j].ClassName()) {
	cout<<" Yield "<<modelCoeffList[j].GetName()
	    <<" must be RooRealVar"<<endl;
	exit(-1);
      }
      RooRealVar &theYield=(RooRealVar&)modelCoeffList[j];
      if (theYield.isConstant()) {
	cout<<" "<<theYield.GetName()<<" is constant, ignored in sPlot"<<endl;
        yield0s.add(theYield);
      } else {
        // force fit limits large enough
        if (theYield.getMin()>-1e6) {
          theYield.setMin(-1e6);
          cout<<" The lower limit of "<<theYield.GetName()
              <<" changed to -1e6"<<endl;
        }
        if (theYield.getMax()<1e6) {
          theYield.setMax(1e6);
          cout<<" The upper limit of "<<theYield.GetName()
              <<" changed to 1e6"<<endl;
        }
        yields.add(theYield);
        // get pdfsWvar
        rarBasePdf *modelCompPdf=(rarBasePdf*)modelPdfList->At(j);
        RooAbsPdf *thePdfWvar=modelCompPdf->getSimPdf();
        if (!thePdfWvar) thePdfWvar=modelCompPdf->getPdf();
        pdfsWvar.add(*thePdfWvar);
        rarPdfsWvar.Add(modelCompPdf);
      }
      // get pdfsWOvar
      RooAbsPdf *thePdfWOvar=(RooAbsPdf*)thisModelWOvarList.at(j);
      RooAbsPdf *theSimPdfWOvar=
	getSimPdf((RooSimultaneous*)theFullPdfWOvar, thePdfWOvar);
      if (!theSimPdfWOvar) theSimPdfWOvar=thePdfWOvar;
      if (theYield.isConstant()) pdf0sWOvar.add(*theSimPdfWOvar);
      else pdfsWOvar.add(*theSimPdfWOvar);
    }
  }
  // fix all params
  fullParams.setAttribAll("Constant");
  // float those in yields
  yields.setAttribAll("Constant", kFALSE);
  // construct PdfOverlay array
  TString sPlotPdfOverlay=
    readConfStr("sPlotPdfOverlay_"+varName,"notSet", _runSec);
  if ("notSet"==sPlotPdfOverlay)
    sPlotPdfOverlay=readConfStr("sPlotPdfOverlay", "direct", _runSec);
  rarStrParser pdfOverlayParser=sPlotPdfOverlay;
  TObjArray pdfsOverlay(yields.getSize());
  for (Int_t i=0; i<yields.getSize(); i++) {
    RooAbsPdf *OpdfI=(RooAbsPdf*)pdfsWvar.at(i);
    RooAbsPdf *OpdfD=((rarBasePdf*)rarPdfsWvar.At(i))->getDPdfWvar(theVar);
    if (!OpdfD) OpdfD=OpdfI;
    RooAbsPdf *theOverlayPdf=OpdfI;
    if ("no"==pdfOverlayParser[0]) theOverlayPdf=0;
    else if ("yes"!=pdfOverlayParser[0]) theOverlayPdf=OpdfD;
    // now check if any specific config pairs
    Int_t idx=pdfOverlayParser.Index(yields[i].GetName())+1;
    if (idx>0) {
      if ("no"==pdfOverlayParser[idx]) theOverlayPdf=0;
      else if("yes"==pdfOverlayParser[idx]) theOverlayPdf=OpdfI;
      else if("direct"==pdfOverlayParser[idx]) theOverlayPdf=OpdfD;
      else { // get its pdf
	rarBasePdf *thePdf=(rarBasePdf*)
	  (_rarPdfs.FindObject(pdfOverlayParser[idx]));
	if (!thePdf)
	  cout<<" Can not find pdf named "<<pdfOverlayParser[idx]
	      <<" for yield "<<yields[i].GetName()
	      <<" in config sPlotPdfOverlay"<<endl;
	else {
	  RooAbsPdf *thisPdf=thePdf->getSimPdf();
	  if (!thisPdf) thisPdf=thePdf->getPdf();
	  theOverlayPdf=thisPdf;
	}
      }
    }
    pdfsOverlay[i]=theOverlayPdf;
  }
  
  yields.Print();
  pdfsWvar.Print();
  pdfsWOvar.Print();
  /// \todo add _protDataVars into projDeps for sPlot
  /// \todo standardize protDataVars, protDataEVars
  // find sPlotNormIgnoredObs
  RooArgSet projDeps(getArgSet("sPlotNormIgnoredObs", kTRUE));
  projDeps.add
    (getArgSet(readConfStr("sPlotNormIgnoredObs","",_runSec),kFALSE));
  // add _conditionalObs
  projDeps.add(_conditionalObs);
  // fit the dataset again
  //RooFitResult *fitStat=
  //  theFullPdfWOvar->fitTo(*sPlotData, _conditionalObs,"ehmr");

  TString fitOption("qemhr");

  // convert to newer fitTo format for steering options
  Bool_t sPlotExtended = fitOption.Contains("e");
  Bool_t sPlotMinos    = fitOption.Contains("m");
  Bool_t sPlotHesse    = fitOption.Contains("h");
  Bool_t sPlotVerbose  = !fitOption.Contains("q");
  Bool_t sPlotSave     = fitOption.Contains("r"); // return results
  RooFitResult *fitStat=
    theFullPdfWOvar->fitTo(*sPlotData, ConditionalObservables(_conditionalObs),
  			   Save(sPlotSave), Extended(sPlotExtended), 
			   Verbose(sPlotVerbose), Hesse(sPlotHesse),
			   Minos(sPlotMinos));
  //Int_t ncpus(1);
  //RooFitResult *fitStat=doTheFit(theFullPdfWOvar, sPlotData, fitOption, ncpus);

  // create sPlot obj
  TString sPlotName=Form("sPlot_%s", theVar->GetName());
  rarSPlot mySPlot(sPlotName, sPlotName, theVar, *sPlotData, fitStat,
                   pdfsWOvar, yields, pdf0sWOvar, yield0s, projDeps);
  // plot comps
  rarStrParser sPlotComps=readConfStr("sPlotComps", "all", _runSec);
  // do we need to have asym plot?
  TString sPlotAsym=readConfStr("sPlotAsym_"+varName,"notSet", _runSec);
  if ("notSet"==sPlotAsym)
    sPlotAsym=readConfStr("sPlotAsym","no", _runSec);
  RooCategory *sPlotAsymCat(0);
  if ("no"!=sPlotAsym) {
    sPlotAsymCat=(RooCategory *)_fullObs->find(sPlotAsym);
    if (!sPlotAsymCat) {
      cout<<"Can not find category: "<<sPlotAsym<<endl;
      exit(-1);
    }
    // make sure it is two-type cat
    if (2!=sPlotAsymCat->numTypes()) {
      cout<<"Asym plot is for two-type category only"<<endl;
      exit(-1);
    }
  }
  // sPlot cat
  TList sPlotCatList;
  TString sPlotCatName=readConfStr("sPlotCat_"+varName,"notSet", _runSec);
  if ("notSet"==sPlotCatName)
    sPlotCatName=readConfStr("sPlotCat","no", _runSec);
  if ("no"!=sPlotCatName) {
    rarStrParser sPlotCatParser=sPlotCatName;
    while (sPlotCatParser.nArgs()>0) {
      TString catName=sPlotCatParser[0];
      sPlotCatParser.Remove();
      RooAbsCategory *sPlotCat=(RooAbsCategory *)_fullObs->find(catName);
      if (!sPlotCat) {
	cout<<" W A R N I N G !"<<endl
	    <<" category "<<catName<<" does not exist!"<<endl;
	continue;
      }
      if (catName==sPlotAsym) {
	cout<<" "<<catName<<" is used for sPlotAsym"<<endl
	    <<" Ignored"<<endl;
	continue;
      }
      sPlotCatList.Add(sPlotCat);
    }
  }
  RooLinkedList plotOpts;
  Double_t xerrorscale = atof(readConfStr("XErrorSize", "1.0", _runSec));
  RooCmdArg xerrArg = XErrorSize(xerrorscale);
  plotOpts.Add((TObject*)&xerrArg);
  RooCmdArg datErrArg = DataError(RooAbsData::SumW2);
  RooCmdArg asymArg;
  if(sPlotAsymCat) {
    asymArg = Asymmetry(*sPlotAsymCat);
    plotOpts.Add((TObject*)&asymArg);
  }
  else
    plotOpts.Add((TObject*)&datErrArg);
  for (Int_t i=0; i<yields.getSize(); i++) {
    RooRealVar &yield=(RooRealVar&)yields[i];
    TString yieldName=yield.GetName();
    TString pdfName=pdfsWvar[i].GetName();
    pdfName.Remove(0,4);
    if (("all"!=sPlotComps[0])&&(!sPlotComps.Have(yieldName))
	&&(!sPlotComps.Have(pdfName))) continue;
    cout<<" For "<<yieldName<<" "<<varName<<endl;
    TString mySPlotName=sPlotName+"_"+yieldName;
    RooDataSet *fillData=mySPlot.fill(yield, nBins, plotMin, plotMax);
    // Set up splitting
    TList sPlotDS;
    TList sPlotSlice;
    TList sPlotNames;
    sPlotDS.Add(new RooDataSet(*fillData));
    sPlotSlice.Add(new RooDataSet(*sPlotData));
    sPlotNames.Add(new TObjString(mySPlotName));
    for( int sPlotCatIdx=0; sPlotCatIdx < sPlotCatList.GetSize(); 
	 sPlotCatIdx++) {
      RooAbsCategory *sPlotCat=(RooAbsCategory*)sPlotCatList.At(sPlotCatIdx);
      TIterator *cIter=sPlotCat->typeIterator();
      RooCatType *cType(0);
      while(cType=(RooCatType*)cIter->Next()) {
	TString catName=sPlotCat->GetName();
	TString typeName=cType->GetName();
	TString catCut=catName+"=="+catName+"::"+typeName;
	RooDataSet *catFillData=(RooDataSet*)fillData->reduce(catCut);
	RooDataSet *catSliceData=(RooDataSet*)sPlotData->reduce(catCut);
	sPlotDS.Add(catFillData);
	sPlotSlice.Add(catSliceData);
	sPlotNames.Add(new TObjString(mySPlotName+"_"+typeName));
      }
      delete cIter;
    }
    for (Int_t sPlotIdx=0; sPlotIdx<sPlotDS.GetSize(); sPlotIdx++) {
      frame=theVar->frame(plotMin, plotMax, nBins);
      TString frameName=sPlotNames.At(sPlotIdx)->GetName();
      RooDataSet* localData = (RooDataSet*)sPlotDS.At(sPlotIdx);
      RooDataSet* localSlice = (RooDataSet*)sPlotSlice.At(sPlotIdx);
      frame->SetNameTitle(frameName,frameName);
      // Check for Hist
      if(sPlotIdx==0 && 
	 !(localData&&"yes"!=readConfStr("sPlotHist", "no", _runSec)))
	frame->addTH1(mySPlot.getSPlotHist(), "E");
      // Plot the data
      if (localData&&"yes"!=readConfStr("sPlotHist", "no", _runSec)) {	
	if (!sPlotAsymCat)
	  localData->plotOn(frame, plotOpts);
	else
	  // We just have to live with bin rounding here:
	  localData->plotOn(frame, plotOpts);
      }
      RooAbsPdf *theOverlayPdf=(RooAbsPdf*)pdfsOverlay[i];
      if (theOverlayPdf) {
	if (theOverlayPdf->dependsOn(*theVar)) {
	  cout<<" Overlay "<<theOverlayPdf->GetName()<<" on sPlot"<<endl;
	  if (!sPlotAsymCat)
	    theOverlayPdf->plotOn(frame, ProjWData(projDeps, *localSlice));
	  else {      	    
	    TIterator *catI=sPlotAsymCat->typeIterator();
	    RooCatType *theCatType=(RooCatType*)catI->Next();
	    TString catCutA=
	      sPlotAsym+"=="+sPlotAsym+"::"+theCatType->GetName();
	    Int_t catIdxA=theCatType->getVal();
	    theCatType=(RooCatType*)catI->Next();
	    TString catCutB=
	      sPlotAsym+"=="+sPlotAsym+"::"+theCatType->GetName();
	    Int_t catIdxB=theCatType->getVal();
	    TString catCutMinus=catIdxA<catIdxB?catCutA:catCutB;
	    TString catCutPlus =catIdxA>catIdxB?catCutA:catCutB;
	    RooDataSet *catLocalData=(RooDataSet*)localData->reduce(catCutMinus);
	    RooDataSet *catLocalSlice=(RooDataSet*)localSlice->reduce(catCutMinus);
	    cout<<"Dataset (-) with cut: "<<catCutMinus<<endl;
	    catLocalData->Print();
	    RooPlot *frameM=theVar->frame(plotMin, plotMax, nBins);
	    catLocalData->plotOn(frameM, DataError(RooAbsData::SumW2));
	    theOverlayPdf->plotOn(frameM, ProjWData(projDeps,*catLocalSlice));
	    RooCurve *B0barCurve = (RooCurve*) frameM->getObject(1);
	    delete catLocalData; delete catLocalSlice;
	    catLocalData=(RooDataSet*)localData->reduce(catCutPlus);
	    catLocalSlice=(RooDataSet*)localSlice->reduce(catCutPlus);
	    cout<<"Dataset (+) with cut: "<<catCutPlus<<endl;
	    catLocalData->Print();
	    RooPlot *frameP=theVar->frame(plotMin, plotMax, nBins);
	    catLocalData->plotOn(frameP, DataError(RooAbsData::SumW2));
	    theOverlayPdf->plotOn(frameP, ProjWData(projDeps,*catLocalSlice));
	    RooCurve *B0Curve = (RooCurve*) frameP->getObject(1);
	    delete catLocalData;
	    delete catI;
	    TString curveName=Form("asymCurve_%s", B0Curve->GetName());
	    RooCurve* asymCurve=combCurves(B0Curve,B0barCurve,
					   "(@0-@1)/(@0+@1)", "@0+@1>0", 
					   curveName, 3, kSolid, kBlue);
	    if(asymCurve)
	      frame->addPlotable(asymCurve);
	  }
	  cout<<endl;
	} else {
	  cout<<" Overlaying pdf "<<theOverlayPdf->GetName()
	      <<" does not depend on "<<theVar->GetName()<<endl
	      <<" Overlaying ignored!!"<<endl;
	}
      }
      // add the sPlot for output
      plotList.Add(frame);
    }
    sPlotDS.Delete();
    sPlotSlice.Delete();
    sPlotNames.Delete();
    if (("yes"==readConfStr("sPlotSaveSWeight", "no", _runSec))&&fillData) {
      TString dName=frame->GetName();
      dName.ReplaceAll("sPlot_", "sWeight_");

      TTree *dataTree = createTreeFromDataset(fillData, kFALSE);
      plotList.Add(dataTree);
    } else {
      delete fillData;
    }    
  }
  
  // restore params
  readFromStr(fullParams, paramSSaver);
  return frame;
}

/// \brief To initialize param values for each action
/// \param act The action name
/// \param fullParams full param to display
/// \param fullParamsWOI full param set to initialize
/// \param readParams Read or not from param file
/// \param paramFileID Default fileID for reading params
/// \param readSecParams Read or not from section
/// \param useFloatFix Float/fix param after initialization
/// \param postPdfFloatSet ArgSet for postPdfFloat
/// \param preMLFixSet ArgSet for preMLFix
/// \param preMLFloatSet ArgSet for preMLFloat
/// \return kTRUE if successful, kFALSE if not.
///
/// It initializes param values for each action.
Bool_t rarMLFitter::initParams(TString act, RooArgSet fullParams,
			       RooArgSet fullParamsWOI,
			       TString readParams, TString paramFileID,
			       TString readSecParams, Bool_t useFloatFix,
			       RooArgSet postPdfFloatSet,
			       RooArgSet preMLFixSet, RooArgSet preMLFloatSet)
{
  Bool_t retVal(kTRUE);
  // read from param file
  TString readParamsConfig=readConfStrCnA("pre"+act+"ReadParams", readParams);
  if ("no"!=readParamsConfig) {
    TString paramFile=getParamFileName(paramFileID, readParamsConfig);
    cout<<endl<<"Reading in params before "<<act<<" from"<<endl
	<<paramFile<<endl;
    paramFileI(fullParamsWOI, paramFile);
  }
  // read from section
  TString readSecParamsConfig=
    readConfStr("pre"+act+"ReadSecParams", readSecParams, _runSec);
  if ("no"!=readSecParamsConfig) {
    if ("yes"!=readSecParamsConfig) {
      cout<<endl<<"Reading in params for "<<act<<" from section \""
	  <<readSecParamsConfig<<"\""<<endl;
      fullParams.readFromFile(_configFile, 0, readSecParamsConfig);
    }
    cout<<endl<<"Reading in params for "<<act<<" from section \""
	<<_runSec<<"\""<<endl;
    fullParams.readFromFile(_configFile, 0, _runSec);
  }
  if (useFloatFix) {
    // make sure postPdfFloat vars floated
    postPdfFloatSet.setAttribAll("Constant", kFALSE);
    // make sure preMLFix vars fixed
    preMLFixSet.setAttribAll("Constant");
    // make sure preMLFloat vars floated
    preMLFloatSet.setAttribAll("Constant", kFALSE);
  }
  cout<<endl<<" The init values of variables for "<<act<<":"<<endl;
  fullParams.Print("v");
  
  return retVal;
}

/// \brief release-independent code for saving dataset as an ntuple
/// \param ds The dataset
/// \param rootfilename Name of file to contain the ntuple
/// \return void
//----------------------------------------
void rarMLFitter::saveAsRootFile(const RooDataSet *ds, 
				 const TString rootfilename,
				 Bool_t withErrors = kTRUE) {

  cout << "Saving tree for "  << ds->GetName() << endl;
  TFile *f = new TFile(rootfilename,"RECREATE");
  if (f != 0) {
    TTree *tree = createTreeFromDataset(ds, withErrors);
    if (tree != 0) {
      tree->Write();
    }  else {
      cout << "Warning: empty tree for dataset " << ds->GetName() << endl;
    }
    f->Close();
    //    delete tree;
  }
  return;
}

/// \brief release-independent code for creating ntuple from a dataset
/// \param ds The dataset
/// \return The pointer to the newly created tree
//------------------------------------------------
TTree *rarMLFitter::createTreeFromDataset(const RooDataSet *ds,
					  Bool_t withErrors = kTRUE) {

  const Int_t NVARS = 5000; // max variables to save
  Double_t array[NVARS]; 

  const Bool_t debug = kFALSE;

  // loop over rows in dataset
  const Int_t nrows = ds->numEntries();

  if (debug) {
    cout << "Number of nrows " << nrows << " for ds " << ds->GetName() << endl;
  }

  if (nrows < 1) {
    cout << "No entries in dataset " << ds->GetName() << endl;
    return(0);
  }

  // get vector of variable names
  const RooArgSet *names = ds->get(0);
  const Int_t nvars = names->getSize();

  if (nvars > NVARS) {
    cout << "Too many observables )"<<nvars<<") to save. "
	 << "Maximum is " << NVARS << endl;
  }

  // create the ntuple
  TTree *tree = new TTree(ds->GetName(), ds->GetTitle());

  // iterate over observable names and add names for the errors
  Int_t icol(0);
  TIterator *iter = names->createIterator();

  // map variable name to an element in the ntuple array
  typedef map<TString, Int_t> MapType;
  MapType my_map;

  RooRealVar *rvar(0);
  while ((rvar = (RooRealVar*) iter->Next())) {
    if (rvar->ClassName() == TString("RooRealVar")) {
      TString varname = rvar->GetName();
      my_map[varname] = icol; icol++;
      if (withErrors) {

	// create names for error variables
	// older versions (<=v5.26) of hasError/hasAsymError do not have a parameter
	if (rvar->hasError()) {
	  if (rvar->getError() != 0) {my_map[varname+"_err"] = icol; icol++;}
	}
	if (rvar->hasAsymError()) {
	  if (rvar->getAsymErrorLo() != 0 || rvar->getAsymErrorHi() != 0) {
	    my_map[varname+"_aerr_lo"] = icol; icol++;
	    my_map[varname+"_aerr_hi"] = icol; icol++;
	  }
	}
      }
      if (debug) {cout << "VARR  " << icol << " " << varname << endl;}
    }
  }

  // look for any remaining non-RooRealVar types
  RooAbsArg *var(0);
  iter->Reset();
  while ((var = (RooAbsArg*) iter->Next())) {
    if (var->ClassName() != TString("RooRealVar")) {
      TString varname = var->GetName();
      if (debug) {cout << "VAR C " << icol << " " << varname << endl;}
      my_map[varname] = icol;
      icol++;
    }
  }

  // now create branches of tree
  MapType::const_iterator it;
  for (it = my_map.begin(); it != my_map.end(); ++it) {
    TString varname = it->first;
    Int_t icol   = it->second;
    tree->Branch(varname.Data(), &array[icol]);
  }

  // loop over each row in dataset and fill array
  for (Int_t irow=0; irow < nrows; irow++) {
    // extract each field value in row and put in tree
    const RooArgSet *row = ds->get(irow);
    icol = -1; // reset column index - not necessary now
    iter->Reset() ; // reset to start
    const Int_t    idefval = -992746521;
    const Double_t ddefval = -992746521;

    while ((var = (RooAbsArg*) iter->Next())) {
      icol++;
      TString varname = var->GetName();
      MapType::const_iterator iter1 = my_map.find(varname);
      if (iter1 != my_map.end()) {icol = iter1->second;}

      // look for RooCategory type events
      Int_t itemp = row->getCatIndex(varname, idefval);
      if (itemp != idefval) {
	array[icol] = (Double_t) itemp;
	continue; // found a RooCategory
      }

      Double_t temp = row->getRealValue(varname, ddefval);
      if (temp != ddefval) {

	array[icol] = temp;

	// look for errors again
	RooRealVar *rvar = dynamic_cast<RooRealVar*>(var);
	Double_t err(0);
	iter1 = my_map.find(varname+"_err");
	if (iter1 != my_map.end()) {
	  icol = iter1->second;
	  array[icol] = rvar->getError();
	}

	// lower asymmetric error
	iter1 = my_map.find(varname+"_aerr_lo");
	if (iter1 != my_map.end()) {
	  icol = iter1->second;
	  if (rvar->hasAsymError()) {
	    err = rvar->getAsymErrorLo();
	  } else {
	    err = 1;
	  }
	  array[icol] = err;
	}

	// upper asymmetric error
	iter1 = my_map.find(varname+"_aerr_hi");
	if (iter1 != my_map.end()) {
	  icol = iter1->second;
	  array[icol] = rvar->getAsymErrorHi();
	  if (rvar->hasAsymError()) {
	    err = rvar->getAsymErrorHi();
	  } else {
	    err = -1;
	  }
	  array[icol] = err;
	}
	continue; // found a RooRealVar
      }
      cout << "Unexpected type " << var->ClassName() << endl;
      
    }
    tree->Fill();
  }

  return (tree);
}

/// \brief The main driver to finish fitting actions
/// \par Config Directives:
/// - <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_pdfFit">pdfFit Action</a>
/// - <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_mlFit">mlFit Action</a>
/// - <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_toyStudy">toyStudy Action</a>
/// - <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_projPlot">projPlot Action</a>
/// - <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_contourPlot">contourPlot Action</a>
/// - <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_sPlot">sPlot Action</a>
///
/// This function can be divided into several different sections,
/// and each section finishes particular jobs, like pdfFit, mlFit,
/// toy studies, projection plot, etc.
/// It is advisable to divide the action config section accordingly
/// and run the program with different config section (option -A)
/// for each type of job.
///
/// \todo It is possible and more elegant to have this function
///       divided into several different functions according to
///       fitting job type.
void rarMLFitter::run()
{
  // first have pre-actions done for each RooRarFitPdf
  cout<<endl<<" Pre-Actions for each RooRarFitPdf"<<endl<<endl;
  preAction();
  // then see if we need to compute correlation matrix
  rarStrParser compCorrParser=
    readConfStr("computeCorrelations", "yes", getDatasets()->getVarSec());
  // compute correlation matrix
  for (Int_t i=0; i<getDatasets()->getDatasetList()->GetSize(); i++) {
    RooDataSet *theData=(RooDataSet*)getDatasets()->getDatasetList()->At(i);
    if (("yes"==compCorrParser[0])||compCorrParser.Have(theData->GetName()))
      computeCorrelations(_fObsSet, theData);
  }
  // then do actions in action section
  cout<<endl<<endl<<" Run rarMLFitter with configs from section \""
      <<_runSec<<"\""<<endl<<endl;
  
  // get the full params
  RooArgSet fullParams(getParams());
  // add parameters in _thePdf
  fullParams.add(*_thePdf->getParameters(_protDataset));
  // let's sorted it as desired
  TString paramOrder=readConfStr("outParamOrder","ascend",getMasterSec());
  if (("ascend"==paramOrder) || ("descend"==paramOrder)) {
    Bool_t reverse=("descend"==paramOrder);
    RooArgList paramList(fullParams);
    fullParams.removeAll();
    paramList.sort(reverse);
    fullParams.add(paramList);
  }
  //fullParams.Print();
  RooArgSet fullParamsWOI(fullParams); // fullParamsWithoutIgnored
  fullParamsWOI.remove(getArgSet("Ignored", kTRUE, &fullParams)); // ignored
  // postPdfFloat
  RooArgSet postPdfFloatSet(getArgSet("postPdfFloat", kTRUE, &fullParams));
  postPdfFloatSet.add(getArgSet(readConfStr("postPdfFloat","", _runSec),kFALSE,
				&fullParams));
  // preMLFix
  RooArgSet preMLFixSet(getArgSet("preMLFix", kTRUE, &fullParams));
  preMLFixSet.add(getArgSet(readConfStr("preMLFix","", _runSec),kFALSE,
			    &fullParams));
  // preMLFloat
  RooArgSet preMLFloatSet(getArgSet("preMLFloat", kTRUE, &fullParams));
  preMLFloatSet.add(getArgSet(readConfStr("preMLFloat","", _runSec),kFALSE,
			      &fullParams));
  
  // do pdf fit
  if ("no"!=readConfStr("pdfFit", "no", _runSec)) {
    cout<<endl<<" Do pdfFit Action"<<endl<<endl;
    // read in params before pdf fit
    initParams("Pdf", fullParams, fullParamsWOI, "no", "pdfFit", "no");
    // do we have a pdf list to fit?
    TString pdfToFit=readConfStr("pdfToFit", "", _runSec);
    if (""!=pdfToFit) cout<<" Pdf(s) to fit: "<<pdfToFit<<endl;
    // do pdf fit
    doPdfFit(pdfToFit);
    // plot pdf
    TString postPdfMakePlot=readConfStr("postPdfMakePlot", "no", _runSec);
    if ("no"!=postPdfMakePlot) {
      // create plot list
      TList *plotList=new TList();
      doPdfPlot(*plotList, pdfToFit);
      TString postPdfPlotFile=getRootFileName("pdfPlot", postPdfMakePlot);
      TFile plotFile(postPdfPlotFile, "recreate");
      plotList->Print();
      plotList->Write();
      plotFile.Close();
      cout<<endl<<"Writing out Pdf plots after pdf fit to "
	  <<postPdfPlotFile<<endl;
    }
    // read in any values in the action section
    TString postPdfReadSecParams="no";
    postPdfReadSecParams= // the old way
      readConfStr("postPdfReadParams", postPdfReadSecParams, _runSec);
    postPdfReadSecParams= // the current way
      readConfStr("postPdfReadSecParams", postPdfReadSecParams, _runSec);
    if ("no"!=postPdfReadSecParams) {
      if ("yes"!=postPdfReadSecParams) {
	cout<<endl<<"Reading in params after pdfFit from section \""
	    <<postPdfReadSecParams<<"\""<<endl;
	fullParams.readFromFile(_configFile, 0, postPdfReadSecParams);
      }
      cout<<endl<<"Reading in params after pdfFit from section \""
	  <<_runSec<<"\""<<endl;
      fullParams.readFromFile(_configFile, 0, _runSec);
    }
    // fix params before output
    //fullParams.Print("v");
    cout<<endl<<"Fix params after pdfFit"<<endl;
    fullParams.setAttribAll("Constant");
    // write params after pdf fit
    TString postPdfWriteParams=readConfStrCnA("postPdfWriteParams", "no");
    if ("no"!=postPdfWriteParams) {
      // do we need to warn?
      if (""!=pdfToFit) {
	fullParams.writeToStream(cout, kFALSE);
	cout<<endl<<"W A R N I N G !"<<endl
	    <<"You do not fit all PDFs!"<<endl
	    <<"No output for the params!"<<endl;
      } else {
	// then write them out
	TString postPdfParamFile=
	  getParamFileName("pdfFit", postPdfWriteParams);
	cout<<endl<<"Writing out params after pdfFit to "
	    <<postPdfParamFile<<endl;
	paramFileO(fullParams, postPdfParamFile);
      }
    }
  } //done pdfFit
  
  // do toys
  if ("no"!=readConfStr("toyStudy", "no", _runSec)) {
    cout<<endl<<" Do toyStudy Action"<<endl<<endl;
    // read in params before toy study
    initParams("Toy", fullParams, fullParamsWOI, "yes", "pdfFit", "yes",
	       kTRUE, postPdfFloatSet, preMLFixSet, preMLFloatSet);
    // call the function to do toys
    RooDataSet *theFitParDataSet=doToyStudy(fullParams);
    TString postToyWriteParams=readConfStr("postToyWriteParams","yes",_runSec);
    if (("no"!=postToyWriteParams) && theFitParDataSet) {
      TString toyRootFile=getRootFileName("toyPlot", postToyWriteParams);
      cout<<endl<<"Writing out param pulls, etc after toy study to "
	  <<toyRootFile<<endl;
      saveAsRootFile(theFitParDataSet, toyRootFile, kTRUE);
    }
  } // done toy
  
  // do ml fit
  if ("no"!=readConfStr("mlFit", "no", _runSec)) {
    cout<<endl<<" Do mlFit Action"<<endl<<endl;
    // to get correlation coeffs?
    TString postMLSysParams=readConfStrCnA("postMLSysParams", "no");
    TString postMLSysVars=readConfStrCnA("postMLSysVars", "no");
    if((!postMLSysParams.BeginsWith("no"))&&(!postMLSysVars.BeginsWith("no")))
    {
      TString vParams="";
      rarStrParser paramsStrParser=postMLSysParams;
      while (paramsStrParser.nArgs()>0) {
	TString theParamName=paramsStrParser[0];
	paramsStrParser.Remove();
	TIterator* iter = fullParams.createIterator();
	RooRealVar *theParam(0);
	while(theParam=(RooRealVar*)iter->Next()) {
	  TString theName=theParam->GetName();
	  if ((theName==theParamName)||(theName.BeginsWith(theParamName+"_"))){
	    vParams=vParams+" \""+theParamName+"\"";
	    break;
	  }
	}
	delete iter;
      }
      rarStrParser vParamsParser=vParams;
      for (Int_t i=0; i<vParamsParser.nArgs(); i++) {
	for (Int_t j=0; j<i; j++) {
	  saveCorrCoeff(getCorrCoefName(vParamsParser[i],vParamsParser[j]),0,
			kTRUE);
	}
      }
      //getCorrCoeffs().Print("v");
    }
    // read in params before pdf fit
    initParams("ML", fullParams, fullParamsWOI, "yes", "pdfFit", "yes",
	        kTRUE, postPdfFloatSet, preMLFixSet, preMLFloatSet);
    // get fit option
    TString mlFitOption=readConfStr("mlFitOption", "hq", _runSec);
    // force option to be extended and return results
    mlFitOption+="er";
    cout<<"mlFit option: \""<<mlFitOption<<"\""<<endl;
    // string to save important results
    stringstream o;
    // mlFitData
    RooDataSet *mlFitData=
      _datasets->getData(readConfStrCnA("mlFitData", "notSet"));
    if(!mlFitData) {
      cout<<endl<<" Please specify in section ["<<_runSec<<"]"<<endl
	  <<" the dataset you want to fit on with config "<<endl<<endl
	  <<"    mlFitData = <datasetName>"<<endl<<endl
	  <<" and ru-run your job"<<endl;
      exit(-1);
    }
    //cout << "Skipping blind test" << endl;
    chkBlind(mlFitData->GetName());
    cout<<" Using dataset "<<mlFitData->GetName()<<" in mlAction"<<endl;
    // get number of cpus option
    Int_t mlFitNumCPU=atoi(readConfStr("useNumCPU", "1", getMasterSec()));

    // create plot list
    TList *plotList=new TList();
    // do mlFit
    RooFitResult *fitResult=doMLFit(mlFitData, mlFitOption, o, mlFitNumCPU);

    // do signf calculation
    TString postMLSignf=readConfStrCnA("postMLSignf", "no");
    if (!postMLSignf.BeginsWith("no"))
      doSignf(mlFitData, postMLSignf, fitResult, fullParams, o);
    // systematic study
    if ((!postMLSysParams.BeginsWith("no"))&&(!postMLSysVars.BeginsWith("no")))
      doSysStudy(mlFitData, postMLSysParams, postMLSysVars, fullParams, o);
    // gof study
    TString postMLGOFChisq=readConfStrCnA("postMLGOFChisq", "no");
    if (!postMLGOFChisq.BeginsWith("no"))
      doGOFChisq(mlFitData, o, plotList);
    // save plot if any to root file
    TString mlPlotFile=readConfStr("mlPlotFile", "default", _runSec);
    mlPlotFile=getRootFileName("mlPlot", mlPlotFile);
    TFile plotFile(mlPlotFile,"recreate");
    plotList->Print();
    plotList->Write();
    TString postMLWriteFitResult=readConfStrCnA("postMLWriteFitResult", "no");
    if("no"!=postMLWriteFitResult) {
      // Save the RooFitResult
      fitResult->Write();      
    }
    plotFile.Close();
    cout<<"Writing out ml plots to "<<mlPlotFile<<endl;
    // total ml output
    cout<<o.str()<<endl;
    // write params after mlfit
    TString postMLWriteParams=readConfStrCnA("postMLWriteParams", "yes");
    if ("no"!=postMLWriteParams) {
      TString postMLParamFile=getParamFileName("mlFit", postMLWriteParams);
      cout<<endl<<"Writing out params after mlFit to "<<postMLParamFile<<endl;
      // make sure all vars fixed before output
      fullParams.setAttribAll("Constant");
      paramFileO(fullParams, postMLParamFile);
    }
  } // done ml fit
  
  // do projection plot
  if ("no"!=readConfStr("projPlot", "no", _runSec)) {
    cout<<endl<<" Do projPlot Action"<<endl<<endl;
    // read in params before projection plot
    initParams("ProjPlot", fullParams, fullParamsWOI, "yes", "mlFit", "no",
	       kTRUE, postPdfFloatSet, preMLFixSet, preMLFloatSet);
    // need to float some/all parameters for plotting to work in standalone mode - FFW
    cout << "Floating all parameters so that projection plots display correctly." << endl;
    fullParams.setAttribAll("Constant",kFALSE);
    // get projData
    RooDataSet *projData=
      _datasets->getData(readConfStrCnA("projPlotData", "notSet"));
    if(!projData) {
      cout<<endl<<" Please specify in section ["<<_runSec<<"]"<<endl
	  <<" the dataset you want to project with config "<<endl<<endl
	  <<"   projPlotData = <datasetName>  // and (see online doc)"<<endl
	  <<"   projPlotData_<obsName> = <datasetName>"<<endl
	  <<endl<<" and ru-run your job"<<endl;
      exit(-1);
    }
    chkBlind(projData->GetName());
    // create plot list
    TList *plotList=new TList();
    // get S and B LL
    getSnB();
    // loop over all projVars
    TString varNames=readConfStr("projVars", "", _runSec);
    if (""==varNames) varNames=readConfStr("projVar","", _runSec);
    rarStrParser varNamesParser=varNames;
    for (Int_t i=0; i<varNamesParser.nArgs(); i++) {
      // get proj variable
      RooRealVar *theVar=(RooRealVar*)_fullObs->find(varNamesParser[i]);
      if(!theVar) {
	cout<<"Can not find obs named "<<varNamesParser[i]<<endl;
	continue;
      }
      doProjPlot(projData, theVar, *plotList);
    }
    // if any LLR plots?
    doLLRPlot(projData, *plotList);
    // output to root file
    TString projPlotFile=readConfStr("projPlotFile", "default", _runSec);
    projPlotFile=getRootFileName("projPlot", projPlotFile);
    TFile plotFile(projPlotFile,"recreate");
    plotList->Print();
    plotList->Write();
    plotFile.Close();
    cout<<"Writing out proj plots to "<<projPlotFile<<endl;
  } // done projection plot

  // do scan plot
  if ("no"!=readConfStr("scanPlot", "no", _runSec)) {
    cout<<endl<<" Do scanPlot Action"<<endl<<endl;
    // read in params before contour plot
    initParams("ScanPlot", fullParams, fullParamsWOI, "yes", "mlFit", "no",
	       kTRUE, postPdfFloatSet, preMLFixSet, preMLFloatSet);
    // create plot list
    TList *plotList=new TList();
    doScanPlot(*plotList);
    TString scanPlotFile=readConfStr("scanPlotFile", "default", _runSec);
    scanPlotFile=getRootFileName("scanPlot", scanPlotFile);
    TFile plotFile(scanPlotFile,"recreate");
    plotList->Print();
    plotList->Write();
    plotFile.Close();
    cout<<endl<<"Writing out scan plots (dataset) to "<<scanPlotFile<<endl;
  } // done scan plot
  
  // do contour plot
  if ("no"!=readConfStr("contourPlot", "no", _runSec)) {
    cout<<endl<<" Do contourPlot Action"<<endl<<endl;
    // read in params before contour plot
    initParams("ContourPlot", fullParams, fullParamsWOI, "yes", "mlFit", "no",
	       kTRUE, postPdfFloatSet, preMLFixSet, preMLFloatSet);
    // create plot list
    TList *plotList=new TList();
    doContourPlot(*plotList);
    TString contourPlotFile=readConfStr("contourPlotFile", "default", _runSec);
    contourPlotFile=getRootFileName("contPlot", contourPlotFile);
    TFile plotFile(contourPlotFile,"recreate");
    plotList->Print();
    plotList->Write();
    plotFile.Close();
    cout<<"Writing out contour plots to "<<contourPlotFile<<endl;
  } // done contour plot
  
  // do sPlot
  if ("no"!=readConfStr("sPlot", "no", _runSec)) {
    cout<<endl<<" Do sPlot Action"<<endl<<endl;
    // read in params before sPlot
    initParams("SPlot", fullParams, fullParamsWOI, "yes", "mlFit", "no",
	       kTRUE, postPdfFloatSet, preMLFixSet, preMLFloatSet);
    // create plot list
    TList *plotList=new TList();
    // loop over all sPlotVars
    TString varNames=readConfStr("sPlotVars", "", _runSec);
    if (""==varNames) varNames=readConfStr("sPlotVar","", _runSec);
    rarStrParser varNamesParser=varNames;
    for (Int_t i=0; i<varNamesParser.nArgs(); i++) {
      // get proj variable
      RooRealVar *theVar=(RooRealVar*)_fullObs->find(varNamesParser[i]);
      if(!theVar) {
	cout<<"Can not find obs named "<<varNamesParser[i]<<endl;
	continue;
      }
      doSPlot(theVar, *plotList);
    }
    TString sPlotFile=readConfStr("sPlotFile", "default", _runSec);
    sPlotFile=getRootFileName("sPlot", sPlotFile);
    TFile plotFile(sPlotFile, "recreate");
    plotList->Print();
    plotList->Write();
    plotFile.Close();
    cout<<"Writing out sPlots to "<<sPlotFile<<endl;
  }// done sPlot

  // do combine plot
  if ("no"!=readConfStr("combinePlot", "no", _runSec)) {
    cout<<endl<<" Do combinePlot Action"<<endl<<endl;
    // read in params before contour plot
    //initParams("CombinePlot", fullParams, fullParamsWOI, "yes", "mlFit", "no",
    //	       kTRUE, postPdfFloatSet, preMLFixSet, preMLFloatSet);
    // create plot list
    TList *plotList=new TList();
    doCombinePlot(*plotList);
    TString combinePlotFile=readConfStr("combinePlotFile", "default", _runSec);
    combinePlotFile=getRootFileName("combinePlot", combinePlotFile);
    TFile plotFile(combinePlotFile,"recreate");
    plotList->Print();
    plotList->Write();
    plotFile.Close();
    cout<<"Writing out combine plots to "<<combinePlotFile<<endl;
  } // done combine plot
}
