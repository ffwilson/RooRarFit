/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarDatasetDef.cc,v 1.14 2014/09/14 17:33:50 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides dataset definition class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides dataset definition class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "rarDatasetDef.hh"

ClassImp(rarDatasetDef)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarDatasetDef::rarDatasetDef()
  : rarConfig(),
    _primaryObs(0), _addonCols(0)
{
  init();
}

/// \brief Default ctor
///
/// \param configFile The config file
/// \param configSec The config section
///
/// The default ctor initializes common data members,
/// and because there is only one instantiation of this class,
/// #_configStr, name and title are set to constant values.
/// It then calls #init.
rarDatasetDef::rarDatasetDef(const char *configFile, const char *configSec)
  : rarConfig(configFile, configSec, "null", "DatasetDef", "Dataset Def"),
    _primaryObs(0), _addonCols(0)
{
  init();
}

rarDatasetDef::~rarDatasetDef()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// It first reads in the info of fields,
/// and it prints out the fields configured int the config section,
/// then it creates those field variables by calling #createAbsVar
/// and adds them to the field list, #_fullObs,
/// as final observables.
void rarDatasetDef::init()
{
  std::cout << "init of rarDatasetDef:" << std::endl;
  // first set #_createFundamental and #_fullNameSchema
  _createFundamental=kTRUE;
  _fullNameSchema="self";
  
  // create the full, primary, and addon obs sets.
  _fullObs=new RooArgSet("theFullObs");
  _primaryObs=new RooArgSet("thePrimaryObs");
  _addonCols=new RooArgSet("theAddOnCols");
  
  // read in field info
  createAbsVars("Fields", _fullObs, _primaryObs);
  // check number of fields
  Int_t nField=_primaryObs->getSize();
  if (nField<=0) {
    std::cout<<"no fields defined in RooStringVar \"Fields\" in section \""
	     <<_configSec<<"\""<< std::endl;
    exit(-1);
  }
  // now for addon cols
  createAbsVars("AddOns", _fullObs, _addonCols);
  // print cout full obs
  _fullObs->Print("v");
  std::cout << std::endl;
}

/// \brief Get RooFormulaVar ArgList
///
/// \param fStrParser The parsed tokens of formula Args
/// \return The formula ArgList
///
/// It creates a RooArgList from \p fStrParser.
/// It checks if the token in \p fStrParser is in #_fullObs,
/// if not, it creats one by calling #createAbsReal,
/// and finally, it adds the var into the ArgList.
/// It repeats until all the tokens are scanned.
RooArgList *rarDatasetDef::getFormulaArgs(rarStrParser fStrParser)
{
  RooArgList *depList=new RooArgList;
  Int_t nArgs=fStrParser.nArgs();
  if (nArgs<=0) return depList;
  for (Int_t i=0; i<nArgs; i++) {
    if (isNumber(fStrParser[i])) break;
    RooAbsArg *theDep=createAbsReal(fStrParser[i], fStrParser[i]); // create it
    depList->add(*theDep);
  }
  
  return depList;
}

/// \brief Set value for var
/// \param var Name of var
/// \param val Value to set
///
/// This function sets value of var defined in #_primaryObs
void rarDatasetDef::setVal(TString var, Double_t val)
{
  RooAbsArg *theArg=_primaryObs->find(var);
  { // RooRealVar
    RooRealVar *theVar=dynamic_cast<RooRealVar*>(theArg);
    if (theVar) {
      theVar->setVal(val);
      return;
    }
  }
  { // RooCategory
    RooCategory *theVar=dynamic_cast<RooCategory*>(theArg);
    if (theVar) {
      theVar->setIndex((Int_t)val);
      return;
    }
  }
  { // RooStringVar
    RooStringVar *theVar=dynamic_cast<RooStringVar*>(theArg);
    if (theVar) {
      theVar->setVal(Form("%8x", (UInt_t) val));
      return;
    }
  }
  
  static Int_t counter(0);
  counter++;
  if (counter<100)
    std::cout <<"Can not find "<<var<< std::endl;
}

/// \brief Set value for var
/// \param var Name of var
/// \param val Value to set
///
/// This function sets value of var defined in #_primaryObs
void rarDatasetDef::setVal(TString var, TString val)
{
  RooAbsArg *theArg=_primaryObs->find(var);
  { // RooStringVar
    RooStringVar *theVar=dynamic_cast<RooStringVar*>(theArg);
    if (theVar) {
      theVar->setVal(val);
      return;
    }
  }
  { // RooCategory
    RooCategory *theVar=dynamic_cast<RooCategory*>(theArg);
    if (theVar) {
      theVar->setLabel(val);
      return;
    }
  }
  
  static Int_t counter(0);
  counter++;
  if (counter<100)
    std::cout <<"Can not find "<<var<< std::endl;
}
