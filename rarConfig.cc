/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarConfig.cc,v 1.62 2014/09/14 17:33:50 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides config class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides config class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"
#include <sstream>
#include <vector>

#include <ctype.h>

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TMatrixDSym.h"

#include "Roo1DTable.h"
#include "RooAbsCategory.h"
#include "RooArgList.h"
#include "RooCategory.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooMappedCategory.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooStringVar.h"
#include "RooThresholdCategory.h"
#include "RooGlobalFunc.h"

#include "RooSuperCategory.h"
#include "RooUnblindOffset.h"
#include "RooUnblindPrecision.h"

#include "rarAdd.hh"
#include "rarArgusBG.hh"
#include "rarDecay.hh"
#include "rarBasePdf.hh"
#include "rarBifurGauss.hh"
#include "rarCBShape.hh"
#include "rarDatasets.hh"
#include "rarExp.hh"
#include "rarGaussModel.hh"
#include "rarGaussian.hh"
#include "rarOsipDisc.hh"
#include "rarGeneric.hh"
#include "rarHistPdf.hh"
#include "rarKeys.hh"
#include "rarMLPdf.hh"
#include "rarNovosibirsk.hh"
#include "rarBallack.hh"
#include "rarPoly.hh"
#include "rarProd.hh"
#include "rarSimPdf.hh"
#include "rarStep.hh"
#include "rarBinned.hh"
#include "rarTriGauss.hh"
#include "rarTwoGauss.hh"
#include "rarUsrPdf.hh"
#include "rarCruijff.hh"
#include "rarMultPdf.hh"
#include "rarLass.hh"
#include "rarRelBreitWigner.hh"
#include "rarVoigtian.hh"
#include "rarGounarisSakurai.hh"
#include "rarFlatte.hh"
#include "rarUniform.hh"
#include "rarThreshold.hh"

#include "rarConfig.hh"

ClassImp(rarConfig)
  ;

TList rarConfig::_rarPdfs;
TList rarConfig::_rarVars;
TList rarConfig::_rarOVars;
RooArgSet rarConfig::_rarCats;
TString rarConfig::_masterSec;
TString rarConfig::_runSec;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarConfig::rarConfig()
  : TNamed ("",""),
    _configFile(""),_configSec(""), _configStr(""),
    _createFundamental(kFALSE), _fullObs(0)
{
  init();
}

/// \brief Default ctor
///
/// \param configFile The config file
/// \param configSec The config section
/// \param configStr The config string
/// \param name The name
/// \param title The title
///
/// Default ctor to set several common data members,
/// and then call #init.
rarConfig::rarConfig(const char *configFile, const char *configSec,
			 const char *configStr,
			 const char *name, const char *title)
  : TNamed(name, title),
    _configFile(configFile), _configSec(configSec), _configStr(configStr),
    _createFundamental(kFALSE), _fullObs(0)
{
  init();
}

/// \brief Trivial dtor
rarConfig::~rarConfig()
{
}

/// \brief Initial function called by ctor
///
/// \p init will be called by every ctor,
/// and every sub-class has its own init function.
/// init in base class rarConfig will just print out config settings
/// and check if the config file, #_configFile, is readable.
void rarConfig::init()
{
  cout<<endl
      <<GetName()<<":"<<endl
      <<" init of rarConfig"<<endl
      <<" configSec: "<<_configSec<<endl
      <<" configStr: "<<_configStr<<endl;
  
  // default name schema
  _fullNameSchema="prefix";
  
  // check if the config file is readable
  ifstream ifs(_configFile);
  if(ifs.fail()) {
    cout<<"rarConfig::init() can not open config file "<<_configFile<<endl;
    exit(-1);
  }
  ifs.close();
}

/// \brief Create a i/o file name based on various config info
///
/// \param dir Directory name of the file
/// \param aType Action type name
/// \param name Conceptual name
/// \param dsName Dataset input section name
/// \param msName Master section name
/// \param cfName Config file name
/// \return The filename constructed
///
/// This function returns a filename based on given arguments,
/// #_configFile, #_masterSec, Dataset input section.
TString rarConfig::getFullFileName(const TString dir, const TString aType, const TString name,
                                   const TString dsName, const TString msName,
                                   const TString cfName)
{
  TString fileName=dir;
  // check dir
  if (gSystem->AccessPathName(dir)&&gSystem->mkdir(dir)) {
    cout << "getFullFileName: Can not access dir "<<dir<<endl;
    exit(-1);
  }
  // add configFile
  TString configFile=cfName;
  if (""==configFile) configFile=_configFile;
  configFile=gSystem->BaseName(configFile);
  // remove .config
  if (configFile.EndsWith(".config"))
    configFile.Remove(configFile.Length()-7, 7);
  fileName+="/"+getFileName(configFile);

  // add dataset sec
  TString datasetName = dsName;
  if (""==datasetName) {
    if (getDatasets()) datasetName=getDatasets()->getVarSec();
    else datasetName=getVarSec();
  }
  fileName+=getFileName(datasetName);

  // add master sec
  TString masterName = msName;
  if (""==masterName) masterName=getMasterSec();
  fileName+=getFileName(masterName);

  // add action type
  fileName+=getFileName(aType);

  // add action name
  fileName+=getFileName(name);
  // remove the last dot
  if (fileName.EndsWith(".")) fileName.Remove(fileName.Last('.'), 1);
  
  return fileName;
}

/// \brief Return a string w/o special chars for file name
///
/// \param in_name Input name for file
/// \return The filename created
///
/// This function returns a filename based on input string
/// with special chars removed.
TString rarConfig::getFileName(const TString in_name) const
{
  TString name = in_name;
  if (("none"==name)||("no"==name)||("yes"==name)||("default"==name))
    return "";
  
  for (Int_t i=0; i<name.Length(); i++)
    if (!isalnum(name[i])) name[i]='_';
  name+=".";
  
  return name;
}

/// \brief Read in a config string
///
/// \param name The name of the config item to read in
/// \param val The default value of the string
/// \param secName The config section from which to read in
/// \return The string it reads in
///
/// A utility function to read in a config item as a string.
/// It actually uses RooFit's
/// <a href="http://roofit.sourceforge.net/docs/classref/RooArgSet.html#RooArgSet:readFromFile" target=_blank>RooArgSet::readFromFile</a>
/// function.
TString rarConfig::readConfStr(const char *name, const char *val,
			       const char *secName)
{
  // construct section name
  TString configSec=_configSec;
  if (secName) configSec=secName;
  TString secVarName=Form("%s_%s", name, configSec.Data());
  // first check if it has been read in
  RooStringVar *theStr=(RooStringVar*)(_configStrSet.find(secVarName));
  if (theStr) return theStr->getVal();
  // not read in yet, try to get it from config file
  RooArgSet strList("Read Config String List");
  RooStringVar strVar(name, "config string", val, 40960);
  strList.add(strVar);
  strList.readFromFile(_configFile, 0, configSec);
  
  Int_t strLen=strlen(strVar.getVal())+1;
  if (strLen<1024) strLen=1024;
  if (strVar.getVal()!=TString(val)) _configStrSet.addOwned
    (*(new RooStringVar(secVarName, secVarName, strVar.getVal(), strLen)));
  
  return strVar;
}

/// \brief Set a config string
/// \param name The config name
/// \param val The value to set
/// \param secName The config section from which to read in
///
/// A utility function to set config string for a config
void rarConfig::setConfStr(const char*name, const char*val, const char*secName)
{
  // construct section name
  TString configSec=_configSec;
  if (secName) configSec=secName;
  TString secVarName=Form("%s_%s", name, configSec.Data());
  // first check if it has been read in
  RooStringVar *theStr=(RooStringVar*)(_configStrSet.find(secVarName));
  if (!theStr) { // create it
    Int_t strLen=strlen(val)+1;
    if (strLen<1024) strLen=1024;
    theStr=new RooStringVar(secVarName, secVarName, val, strLen);
    _configStrSet.addOwned(*theStr);
  }
  theStr->setVal(val);
  if (!val) { // we should remove this config
    _configStrSet.remove(*theStr);
    delete theStr;
  }
}

/// \brief Add a string to a config string
/// \param name The config name
/// \param val The value to add
/// \param secName The config section from which to read in
///
/// A utility function to add config string for a config
void rarConfig::addToConfStr(const char*name,const char*val,const char*secName)
{
  TString configStr=readConfStr(name, "", secName);
  setConfStr(name, configStr+" "+val, secName);
}

/// \brief Read config string from master section and action section
/// \param configStr Config string
/// \param defVal Default value
/// \return The string read in
///
/// It read in config string from master section and action section.
/// Master section has higher priority.
TString rarConfig::readConfStrCnA(TString configStr, TString defVal)
{
  TString retVal=readConfStr(configStr, "notSet", getVarSec());
  if ("notSet"==retVal) retVal=readConfStr(configStr, defVal, _runSec);
  return retVal;
}

/// \brief Check if a string can be interpreted as a real number
///
/// \param numStr String to check
/// \retval true If the first character is a digit, `.', `+', or `-',
///              or if it can be converted to non-zero number by atof
/// \retval false Otherwise
///
/// It checks if a string is actually a number.
Bool_t rarConfig::isNumber(TString numStr)
{
  Bool_t retVal(kTRUE);
  char c=numStr[0];
  if (isdigit(c)) return retVal;
  if (0!=atof(numStr)) return retVal;
  if ('.'==c) return retVal;
  if (('+'==c)||('-'==c)) return retVal;
  
  retVal=kFALSE;
  return retVal;
}

/// \brief Check if a string is a varType
///
/// \param typeStr String to check
/// \retval true If string is defined varType
/// \retval false Otherwise
///
/// It checks if a string is a varType
Bool_t rarConfig::isVarType(TString typeStr)
{
  if (
      "RooRealVar"==typeStr ||
      "RooConstVar"==typeStr ||
      "RooUnblindOffset"==typeStr ||
      "RooUnblindPrecision"==typeStr ||
      "RooCategory"==typeStr ||
      "RooMappedCategory"==typeStr ||
      "RooThresholdCategory"==typeStr ||
      "RooSuperCategory"==typeStr ||
      "RooStringVar"==typeStr ||
      "RooFormulaVar"==typeStr
      ) return kTRUE;
  return kFALSE;
}

/// \brief Save RooArgSet to a string
/// \param aSet RooArgSet to save
/// \param aStr string to save the RooArgSet
///
/// It saves RooArgSet to a string
void rarConfig::writeToStr(RooArgSet &aSet, string &aStr)
{
  stringstream aSaver;
  aSet.writeToStream((ostream&)aSaver, kFALSE);
  aStr=aSaver.str();  
}

/// \brief Restore values of RooArgSet from a string
/// \param aSet RooArgSet to restore values
/// \param aStr string to restore the RooArgSet
///
/// It restores RooArgSet from a string
void rarConfig::readFromStr(RooArgSet &aSet, string &aStr)
{
  stringstream aSaver;
  aSaver.str(aStr);
  aSet.readFromStream((istream&)aSaver, kFALSE);
}

/// \brief Get the `header' part of varStr
/// \param varStr The varStr to parse
/// \param option Optional string to parse
/// \param fullName Full name constructed
/// \param fullTitle Full Title constructed
/// \param varType Type of the var
/// \param Name Short name in varStr
/// \param Title Short title in varStr/option
/// \param Unit Unit in varStr/option
/// \return The leftover of the string parser
///
/// This function parses the header part of a varStr
/// to determine the var's full name, title, type, etc.,
/// and return all the tokens left.
rarStrParser rarConfig::getVarTNTU(TString varStr, TString option,
				   TString *fullName, TString *fullTitle,
				   TString *varType, TString *Name,
				   TString *Title, TString *Unit)
{
  // get the parsers
  rarStrParser varStrParser=varStr;
  rarStrParser optionParser=option;
  TString myFullName(""), myFullTitle(""), myVarType(""),
    myName(""), myTitle(""), myUnit("");
  Bool_t fNameSpecified(kFALSE);
  // first for short name
  if (varStrParser.nArgs()>0) {
    myName=varStrParser[0];
    varStrParser.Remove();
  }
  if (varStrParser.nArgs()>0) {
    // check if the next token is name or not
    if (!isVarType(varStrParser[0])) { // it is full name, or new schema
      if ((varStrParser.nArgs()>1)&&isVarType(varStrParser[1])){// the old way
	myFullName=varStrParser[0];
	varStrParser.Remove();
	fNameSpecified=kTRUE; // full name specified
      } else if (varStrParser.nArgs()==1) { // full name only
	myFullName=varStrParser[0];
	varStrParser.Remove();
	fNameSpecified=kTRUE; // full name specified
      } else { // the new schema
	myVarType="nsRooRealVar"; // new schema RRV
	myFullName=myName;
	if (optionParser.Have("N")) {
	  myFullName=optionParser[optionParser.Index("N")+1];
	  fNameSpecified=kTRUE; // full name specified
	}
	if (optionParser.Have("T"))
	  myTitle=optionParser[optionParser.Index("T")+1];
	if (optionParser.Have("U"))
	  myUnit=optionParser[optionParser.Index("U")+1];
	// remember current index 0 for N, 1 for T and 2 for U
	Int_t cvIdx(0);
	while((varStrParser.nArgs()>0)&&(!isNumber(varStrParser[0]))) {
	  TString flag=varStrParser[0];
	  varStrParser.Remove();
	  if (varStrParser.nArgs()>0) {
	    if ("N"==flag) {
	      myFullName=varStrParser[0];
	      varStrParser.Remove();
	      fNameSpecified=kTRUE; // full name specified
	    } else if ("T"==flag) {
	      myTitle=varStrParser[0];
	      varStrParser.Remove();
	    } else if ("U"==flag) {
	      myUnit=varStrParser[0];
	      varStrParser.Remove();
	    } else { // OK flag is the one
	      if (0==cvIdx) {
		myFullName=flag;
		fNameSpecified=kTRUE; // full name specified
	      } else if(1==cvIdx) myTitle=flag;
	      else if(2==cvIdx) myUnit=flag;
	    }
	  } else {
	    if (0==cvIdx) {
	      myFullName=flag;
	      fNameSpecified=kTRUE; // full name specified
	    } else if(1==cvIdx) myTitle=flag;
	    else if(2==cvIdx) myUnit=flag;
	  }
	  cvIdx++;
	}
      }
    }
    // let's make sure the name is not notSet
    if ("notSet"==myName) {
      cout<<"var name should not be notSet, please give the var a name"<<endl
	  <<"varStr: "<<varStr<<endl;
      exit(-1);
    }
    // now, find varType
    if ((varStrParser.nArgs()>0)&&(""==myVarType)) {
      if (!isVarType(varStrParser[0])) { // wrong
	cout<<"This is supposed to be a varType "<<varStrParser[0]<<endl;
	exit(-1);
      }
      myVarType=varStrParser[0];
      varStrParser.Remove();
      // then the next token is title
      myTitle=varStrParser[0];
      varStrParser.Remove();
    }
  }
  if (!fNameSpecified) {
    myFullName=getFullVarName(myName);
    if (myFullName==myName) myFullTitle=myTitle;
    else myFullTitle=myTitle; //+" "+GetTitle();
  } else {
    myFullTitle=myTitle;
  }
  
  if (fullName) *fullName=myFullName;
  if (fullTitle) *fullTitle=myFullTitle;
  if (varType) *varType=myVarType;
  if (Name) *Name=myName;
  if (Title) *Title=myTitle;
  if (Unit) *Unit=myUnit;
  return varStrParser;
}

/// \brief Get a full var name based on the input and naming schema
///
/// \param nameStr The initial name
/// \return The generated full name
///
/// It returns a full name based on the input string and naming schema.
/// The \p self method should be disabled in base definition.
TString rarConfig::getFullVarName(TString nameStr)
{
  TString myName=GetName();
  TString myFullName;
  if ("prefix"==_fullNameSchema) myFullName=myName+"_"+nameStr;
  else if ("suffix"==_fullNameSchema) myFullName=nameStr+"_"+myName;
  else myFullName=nameStr;
  return myFullName;
}

/// \brief Get a full var name based on the input and naming schema
///
/// \param nameStr The initial name
/// \param varStr Var def string
/// \param option option def string
/// \return The generated full name
///
/// It returns a full name based on the input string and naming schema.
/// The \p self method should be disabled in base definition.
TString rarConfig::getFullVarName(TString nameStr,
				  TString varStr, TString option)
{
  TString myFullName;
  if ((""==varStr)&&(""==option)) return getFullVarName(nameStr);
  getVarTNTU(varStr, option, &myFullName);
  return myFullName;
}

/// \brief To create a RooAbsReal var
///
/// \param name The name of the var to create
/// \param title Its title
/// \param val Initial value of the var
/// \param min Lower limit
/// \param max Higher limit
/// \param unit Optional unit name for the var,
///             and/or "RooConstVar" for constant var
/// \return The RooAbsReal var created
///
/// It is mainly called to create pdf parameters.
/// It constructs a string, \p varStr, based on its input parameters.
/// It also check the config section to see
/// if the var's config item exists there.
/// If so, \p varStr is set to that value.
/// Then, it calls #createAbsVar with \p varStr passed to it
/// to get the var created.
/// If \p varStr is a number (checked by #isNumber),
/// #createAbsVar is called with parameter \p rescan set to true.
///
/// When parameters \p min and \p max are both 0, \p varStr is constructed as
/// - ``<tt>\<name\> RooRealVar "<title>" \<val\></tt>''
/// .
/// otherwise, it is created as
/// - ``<tt>\<name\> RooRealVar "<title>" \<val\> \<min\> \<max\></tt>''
/// .
/// If the char pointer \p unit is not null, it will be appended to
/// \p varStr,
/// and if \p unit contains "RooConstVar",
/// it will construct a stirng to create \p RooConstVar.
///
/// After \p varStr is constructed, the function will look into the config
/// section to see if there is config item named \p name in that section
/// using #readConfStr,
/// if it is true, it will replace \p varStr with that string,
/// and if that string is actually a number (checked by #isNumber),
/// it will set \p rescan parameter to true.
/// Finally, it calls #createAbsVar to create the var.
RooAbsReal *rarConfig::createAbsReal(const char *name, const char *title,
				     const Double_t val,
				     const Double_t min, const Double_t max,
				     const char *unit)
{
  RooAbsReal *theVar(0);
  if (max<min) {
    cout<<"max="<<max<<" must >= "<<"min="<<min<<endl;
    exit(-1);
  }
  // check option to see if it is for RooConstVar
  TString option=unit;
  Bool_t forRooConstVar(kFALSE);
  if (option.Contains("RooConstVar")) {
    forRooConstVar=kTRUE;
    option.ReplaceAll("RooConstVar", "");
  }
  
  // create the construct string
  TString varStr=Form("%s RooRealVar \"%s\"", name, title);
  if (forRooConstVar) varStr=Form("%s RooConstVar \"%s\"", name, title);
  if ((0==min) && (0==max)) {
    varStr=Form("%s %f", varStr.Data(), val);
  } else {
    varStr=Form("%s %f %f %f", varStr.Data(), val, min, max);
  }
  if (""!=option) varStr=Form("%s \"%s\"", varStr.Data(), option.Data());
  TString varConfigString=readConfStr(name, "notSet", getVarSec());
  if ("notSet"!=varConfigString) {
    if (isNumber(varConfigString)&&forRooConstVar) {
      varStr=
	Form("%s RooConstVar \"%s\" %s", name, title, varConfigString.Data());
    } else {
      // construct new varStr
      varStr=Form("%s %s", name, varConfigString.Data());
    }
  }
  // whatever in option right now is actually unit
  if (""!=option) option ="U "+option;
  // now add title into option
  option.Append(Form(" T %s ", title));
  // the actual function to create variable
  theVar=(RooAbsReal*)createAbsVar(varStr, option);
  return theVar;
}

/// \brief To create a RooAbsReal var
///
/// \param name The name of the var to create
/// \param title Its title
/// \param val Initial value of the var
/// \param unit Optional unit name for the var,
///             and/or "RooConstVar" for constant var
/// \return The RooAbsReal var created
///
/// It calls #createAbsReal will \p min and \p max set to 0.
/// It is mainly used to create constant parameters.
/// If the option given is "RooConstVar",
/// the created variable will be \p RooConstVar;
/// otherwise it is the unit of the variable if specified.
RooAbsReal *rarConfig::createAbsReal(const char *name, const char *title,
				       const Double_t val, const char *unit)
{
  return createAbsReal(name, title, val, 0, 0, unit);
}

/// \brief To create a RooAbsArg var
///
/// \param varStr The config string
/// \param option Option to create the var.
///               You can put optional units, etc., here.
/// \return The \p RooAbsArg var created
///
/// It creates a \p RooAbsArg var according to its config string, \p varStr.
/// The config string is of form like:
/// - ``<tt>\<varName\> [\<finalName\>] \<varType\>
///     [\<token for this var\>...]</tt>''
/// .
/// where \p varName is the initial name of the var to be created,
/// (and it is the name of config item in the config section),
/// an optional \p finalName, if set, will be the created var's final
/// name, instead of generated from \p varName and RooRarFitPdf's name,
/// \p varType is the actual RooAbsArg type to be created, including
/// \p \b RooRealVar,
/// \p \b RooConstVar,
/// \p \b RooUnblindOffset,
/// \p \b RooUnblindPrecision,
/// \p \b RooCategory,
/// \p \b RooMappedCategory,
/// \p \b RooThresholdCategory,
/// \p \b RooSuperCategory,
/// \p \b RooStringVar,
/// and \p \b RooFormulaVar.
/// The tokens of the rest of the string depend on the type of the var
/// and how it is created.
/// Usually there is a title token right after \p varType for \p RooRealVar.
/// The config string \p varStr is constructed by the caller in various ways.
/// It can be formed by #createAbsReal or directly by RooRarFit Pdf ctor,
/// #init,
/// where, in the latter cases, it is read in
/// with config item named \p varName from config file.
/// For example, if the config item is
/// - myVarName=blah1 blah2 blah3 ...
/// .
/// then the config string formed will be
/// - ``myVarName blah1 blah2 blah3 ...''
/// .
/// The function will first check if the var with the same final name
/// has been created,
/// if so, it returns that var instead of creating a new one.
/// The final name, if not specified in \p varStr, will be generated
/// based on \p varName and the name of RooRarPdf object who creates
/// this var. For example, if \p varName is ``myVarName'', and the RooRarPdf
/// object's name is ``myPdfName'', the final name for the var will be
/// ``myVarName_myPdfName''.
///
/// \todo Print Unblind Status after job
///
/// \par Config String Structure:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#item_createVar">See doc for general rules to create variable.</a>
RooAbsArg *rarConfig::createAbsVar(const char *varStr, const char *option)
{
  RooAbsArg *theVar(0);
  RooAbsArg *theFVar(0); // fundamental counterpart
  TString fullName, fullTitle, varType, myName, myTitle, Unit;
  rarStrParser varStrParser=getVarTNTU
    (varStr, option, &fullName,&fullTitle, &varType, &myName,&myTitle, &Unit);
  
  // check if it has been created
  theVar=getAbsVar(fullName);
  // now check if the short name is actually created
  if (!theVar) theVar=getAbsVar(myName);
  if (theVar) {
    if ("RooRealVar"==TString(theVar->ClassName())) {
      // add it to the params list if it is not in obs list
      if (_fullObs->find(theVar->GetName())) {
	addToObs((RooRealVar*)theVar);
      } else {
	addToParams((RooRealVar*)theVar);
      }
    }
    return theVar;
  }
  
  if ("nsRooRealVar"==varType) { // new schema
    theVar=new RooRealVar(myName, fullTitle, 0, Unit);
    // set label
    ((RooRealVar*)theVar)->setPlotLabel(myTitle);
    // construct the RRV output stream
    TString rrvStr=myName+"=";
    while (varStrParser.nArgs()>0) {
      rrvStr+=" "+varStrParser[0];
      varStrParser.Remove();
    }
    stringstream rrvStream(rrvStr.Data());
    RooArgSet varSet(*theVar);
    varSet.readFromStream((istream&)rrvStream, kFALSE);
    addToParams((RooRealVar*)theVar);
  } else if ("RooRealVar"==varType) {  // dealing with RooRealVar
    Double_t min(0), max(0), val(0);
    Int_t nBins(100);
    TString unit="";
    Bool_t setToConstant(kFALSE);
    if (("C"==varStrParser[0])||("c"==varStrParser[0])) {
      varStrParser.Remove();
      setToConstant=kTRUE;
    }
    if (varStrParser.nArgs()>=5) { // #1
      val=atof(varStrParser[0]);
      min=atof(varStrParser[1]);
      max=atof(varStrParser[2]);
      TString binStr=varStrParser[3](2, varStrParser[3].Length());
      nBins=atoi(binStr);
      unit=varStrParser[4];
    } else if (4==varStrParser.nArgs()) {
      if (isNumber(varStrParser[2])) {
	val=atof(varStrParser[0]);
	min=atof(varStrParser[1]);
	max=atof(varStrParser[2]);
	if (varStrParser[3].Contains("B(")) { // #2
	  TString binStr=varStrParser[3](2, varStrParser[3].Length());
	  nBins=atoi(binStr);
	} else { // #3
	  unit=varStrParser[3];
	}
      } else { // #4
	min=atof(varStrParser[0]);
	max=atof(varStrParser[1]);
	TString binStr=varStrParser[2](2, varStrParser[2].Length());
	nBins=atoi(binStr);
	val=(max+min)/2.;
	unit=varStrParser[3];
      }
    } else if (3==varStrParser.nArgs()) {
      if (varStrParser[2].Contains("B(")) { // #5
	min=atof(varStrParser[0]);
	max=atof(varStrParser[1]);
	TString binStr=varStrParser[2](2, varStrParser[2].Length());
	nBins=atoi(binStr);
	val=(max+min)/2.;
      } else if (!isNumber(varStrParser[2])) { // #6
	min=atof(varStrParser[0]);
	max=atof(varStrParser[1]);
	val=(max+min)/2.;
	unit=varStrParser[2];
      } else { // #7
	val=atof(varStrParser[0]);
	min=atof(varStrParser[1]);
	max=atof(varStrParser[2]);
      }
    } else if (2==varStrParser.nArgs()) {
      // check if varStrParser[1] is unit
      if (isNumber(varStrParser[1])) { // #8
	min=atof(varStrParser[0]);
	max=atof(varStrParser[1]);
	val=(max+min)/2.;
      } else { // #9
	val=atof(varStrParser[0]);
	unit=varStrParser[1];
      }
    } else if (1==varStrParser.nArgs()) { // #10
      val=atof(varStrParser[0]);
    }
    if (min==max) {
      theVar=new RooRealVar(myName, fullTitle, val, unit);
    } else {
      theVar=new RooRealVar(myName, fullTitle, val, min, max, unit);
    }
    // set #bins
    ((RooRealVar*)theVar)->setBins(nBins);
    // set label
    ((RooRealVar*)theVar)->setPlotLabel(myTitle);
    // set error
    Double_t theError(0);
    if (theError<=0) theError=(max-min)/2.;
    if (theError>100) theError=100;
    if (min!=max) ((RooRealVar*)theVar)->setError(theError);
    if (setToConstant) ((RooRealVar*)theVar)->setConstant(kTRUE);
    // add it to the params list
    addToParams((RooRealVar*)theVar);
  } else if ("RooConstVar"==varType) {// RooConstVar
    Double_t val=0;
    if (varStrParser.nArgs()>0) val=atof(varStrParser[0]);
    theVar=new RooConstVar(myName, fullTitle, val);
  } else if ("RooUnblindOffset"==varType) {// RooUnblindOffset
    if (varStrParser.nArgs()<=2) {
      cout<<"\""<<varStr<<"\""<<endl
	  <<"At least 3 more parameters needed"<<endl;
      exit(-1);
    }
    TString blindString=varStrParser[0];
    varStrParser.Remove();
    Double_t scale=atof(varStrParser[0]);
    varStrParser.Remove();
    TString blindValueName=varStrParser[0];
    varStrParser.Remove();
    RooAbsReal *blindValue=(RooAbsReal *)
      createAbsReal(blindValueName, blindValueName);
    RooAbsCategory *blindState(0);
    if (varStrParser.nArgs()>1) { // have blindState cat
      TString blindStateName=varStrParser[0];
      varStrParser.Remove();
      blindState=(RooAbsCategory *)
	createAbsReal(blindStateName, blindStateName);
    }
    if (blindState) {
      theVar=new RooUnblindOffset
	(myName, fullTitle, blindString, scale, *blindValue, *blindState);
    } else {
      theVar=new RooUnblindOffset
	(myName, fullTitle, blindString, scale, *blindValue);
    }
  } else if ("RooUnblindPrecision"==varType) {// RooUnblindPrecision
    if (varStrParser.nArgs()<=4) {
      cout<<"\""<<varStr<<"\""<<endl
	  <<"At least 5 more parameters needed"<<endl;
      exit(-1);
    }
    TString blindString=varStrParser[0];
    varStrParser.Remove();
    Double_t centralValue=atof(varStrParser[0]);
    varStrParser.Remove();
    Double_t scale=atof(varStrParser[0]);
    varStrParser.Remove();
    TString blindValueName=varStrParser[0];
    varStrParser.Remove();
    RooAbsReal *blindValue=(RooAbsReal *)
      createAbsReal(blindValueName, blindValueName);
    RooAbsCategory *blindState(0);
    if (varStrParser.nArgs()>1) { // have blindState cat
      TString blindStateName=varStrParser[0];
      varStrParser.Remove();
      blindState=(RooAbsCategory *)
	createAbsReal(blindStateName, blindStateName);
    }
    Bool_t sin2betaMode(kFALSE);
    TString sin2betaModeStr=varStrParser[0];
    varStrParser.Remove();
    if (("kTRUE"==sin2betaModeStr) ||
	("yes"==sin2betaModeStr) ||
	(0!=atoi(sin2betaModeStr))) {
      sin2betaMode=kTRUE;
    }
    if (blindState) {
      theVar=new
	RooUnblindPrecision(myName, fullTitle, blindString, centralValue,
			    scale, *blindValue, *blindState, sin2betaMode);
    } else {
      theVar=new
	RooUnblindPrecision(myName, fullTitle, blindString, centralValue,
			    scale, *blindValue, sin2betaMode);
    }
  } else if ("RooCategory"==varType) { // dealing with RooCategory
    if (varStrParser.nArgs()<=1) {
      cout<<"\""<<varStr<<"\""<<endl
	  <<"At least 2 more parameters needed"<<endl;
      exit(-1);
    }
    TString idxOption=varStrParser[0];
    varStrParser.Remove();
    Bool_t useIdx(kFALSE);
    if (isNumber(idxOption)) {
      Int_t nCat=atoi(idxOption);
      Int_t nArgs=varStrParser.nArgs();
      if (nArgs==nCat*2) useIdx=kTRUE;
    } else {
      if ("noIdx"!=idxOption) useIdx=kTRUE;
    }
    if (useIdx && (varStrParser.nArgs()<=1)) {
      cout<<"Please give also index of cat"<<endl;
      exit(-1);
    }
    theVar=new RooCategory(myName, fullTitle);
    Int_t nTokenPerCat=1;
    if (useIdx) nTokenPerCat=2;
    Int_t nCat=varStrParser.nArgs()/nTokenPerCat;
    for (Int_t j=0; j<nCat; j++) {
      if (useIdx) { // #1
	((RooCategory*)theVar)->defineType(varStrParser[j*nTokenPerCat],
					 atoi(varStrParser[1+j*nTokenPerCat]));
      } else { // #2
	((RooCategory*)theVar)->defineType(varStrParser[j*nTokenPerCat]);
      }
    }
    _rarCats.add(*theVar); // remember this category
  } else if ("RooMappedCategory"==varType) {
    if (varStrParser.nArgs()<=2) {
      cout<<"\""<<varStr<<"\""<<endl
	  <<"At least 3 more parameters needed"<<endl;
      exit(-1);
    }
    TString idxOption=varStrParser[0];
    varStrParser.Remove();
    Bool_t useIdx(kFALSE);
    if ("noIdx"!=idxOption) useIdx=kTRUE;
    if (useIdx && (varStrParser.nArgs()<=2)) {
      cout<<"Please give default index"<<endl;
      exit(-1);
    }
    TString inputCatName=varStrParser[0];
    varStrParser.Remove();
    RooAbsCategory *inputCat=dynamic_cast<RooAbsCategory*>
      (createAbsVar(inputCatName+" "+
                    readConfStr(inputCatName, "", getVarSec())));
    if (!inputCat) {
      cout<<" Can not find cat named "<<inputCatName<<endl;
      exit(-1);
    }
    TString defCatName=varStrParser[0];
    varStrParser.Remove();
    Int_t defCatIdx=RooMappedCategory::NoCatIdx;
    if (useIdx) {
      defCatIdx=atoi(varStrParser[0]);
      varStrParser.Remove();
    }
    theVar=new RooMappedCategory(myName, fullTitle, *inputCat,
				 defCatName, defCatIdx);
    Int_t nTokenPerMap=2;
    if (useIdx) nTokenPerMap=3;
    Int_t nMap=varStrParser.nArgs()/nTokenPerMap;
    for (Int_t j=0; j<nMap; j++) {
      Int_t catIdx=RooMappedCategory::NoCatIdx;
      if (useIdx) catIdx=atoi(varStrParser[2+j*nTokenPerMap]);
      ((RooMappedCategory*)theVar)->map(varStrParser[j*nTokenPerMap],
					varStrParser[1+j*nTokenPerMap],catIdx);
    }
    if (_createFundamental) {
      theFVar=(RooCategory*) theVar->createFundamental();
      _rarCats.add(*theFVar); // remember this category
    } else _rarCats.add(*theVar); // remember this category
  } else if ("RooThresholdCategory"==varType) {
    if (varStrParser.nArgs()<=2) {
      cout<<"\""<<varStr<<"\""<<endl
	  <<"At least 3 more parameters needed"<<endl;
      exit(-1);
    }
    TString idxOption=varStrParser[0];
    varStrParser.Remove();
    Bool_t useIdx(kFALSE);
    if ("noIdx"!=idxOption) useIdx=kTRUE;
    if (useIdx && (varStrParser.nArgs()<=2)) {
      cout<<"Please give default index"<<endl;
      exit(-1);
    }
    TString inputVarName=varStrParser[0];
    varStrParser.Remove();
    RooAbsReal *inputVar=(RooAbsReal *)
      createAbsVar(inputVarName+" "+
		   readConfStr(inputVarName, "", getVarSec()));
    TString defCatName=varStrParser[0];
    varStrParser.Remove();
    Int_t defCatIdx=0;
    if (useIdx) {
      defCatIdx=atoi(varStrParser[0]);
      varStrParser.Remove();
    }
    theVar=new RooThresholdCategory(myName, fullTitle, *inputVar,
				    defCatName, defCatIdx);
    Int_t nTokenPerThres=2;
    if (useIdx) nTokenPerThres=3;
    Int_t nThres=varStrParser.nArgs()/nTokenPerThres;
    for (Int_t j=0; j<nThres; j++) {
      Int_t catIdx=-99999;
      if (useIdx) catIdx=atoi(varStrParser[2+j*nTokenPerThres]);
      ((RooThresholdCategory*)theVar)->
	addThreshold(atof(varStrParser[j*nTokenPerThres]),
		     varStrParser[1+j*nTokenPerThres],catIdx);
    }
    if (_createFundamental) {
      theFVar=(RooCategory*) theVar->createFundamental();
      _rarCats.add(*theFVar); // remember this category
    } else _rarCats.add(*theVar); // remember this category
  } else if ("RooSuperCategory"==varType) {
    if (varStrParser.nArgs()<=1) {
      cout<<"\""<<varStr<<"\""<<endl
	  <<"At least 2 more parameters needed"<<endl;
      exit(-1);
    }
    RooArgSet inputCats;
    while (varStrParser.nArgs()>0) {
      TString inputCatName=varStrParser[0];
      varStrParser.Remove();
      RooAbsCategory *inputCat=dynamic_cast<RooAbsCategory*>
        (createAbsVar(inputCatName+" "+
                      readConfStr(inputCatName, "", getVarSec())));
      if (!inputCat) {
        cout<<" Can not find cat named "<<inputCatName<<endl;
        exit(-1);
      }
      inputCats.add(*inputCat);
    }
    theVar=new RooSuperCategory(myName, fullTitle, inputCats);
    if (_createFundamental) {
      theFVar=(RooCategory*) theVar->createFundamental();
      _rarCats.add(*theFVar); // remember this category
    } else _rarCats.add(*theVar); // remember this category
  } else if ("RooStringVar"==varType) { // dealing with RooStringVar
    TString val="";
    if (varStrParser.nArgs()) val=varStrParser[0];
    theVar=new RooStringVar(myName, fullTitle, val);
  } else if ("RooFormulaVar"==varType) { // dealing with RooFormulaVar
    theVar=new RooFormulaVar(myName, myTitle, *getFormulaArgs(varStrParser));
    if (_createFundamental) {
      theFVar=(RooRealVar*) theVar->createFundamental();
      // check if we specify range
      Bool_t gotMin(kFALSE);
      Double_t min(0), max(0);
      while (varStrParser.nArgs()>0) {
	if (!isNumber(varStrParser[0])) {
	  varStrParser.Remove();
	  continue;
	}
	if (!gotMin) {
	  min=atof(varStrParser[0]);
	  gotMin=kTRUE;
	} else {
	  max=atof(varStrParser[0]);
	}
	varStrParser.Remove();
      }
      if (gotMin) {
	if (min>max) {
	  Double_t v=min;
	  min=max;
	  max=v;
	}
	((RooRealVar*)theFVar)->setRange(min, max);
	((RooRealVar*)theFVar)->setVal((min+max)/2.);
      }
    }
  } else { // not implemented for this type
    cout<<"\""<<varStr<<"\""<<endl
	<<"Var type \""<<varType<<"\" not implemented"<<endl;
    exit(-1);
  }
  
  theVar->SetName(fullName);
  if (theFVar) {
    theFVar->SetName(fullName);
    _rarVars.Add(theFVar);
    _rarOVars.Add(theVar);
  } else {
    _rarVars.Add(theVar);
  }
  
  cout<<varType<<" "<<fullName<<" created:\t"<<varStr<<endl;
  return theVar;
}

/// \brief Find a variable
///
/// \param varName The name of the variable
/// \return The AbsVar created, or null if not found
///
/// \p varName can be the var's full name,short name, or config name.
/// For the later two to work, the var must be defined in the same section.
RooAbsArg *rarConfig::getAbsVar(TString varName)
{
  RooAbsArg *theVar(0);
  if (_createFundamental) {
    // for _rarOVars;
    // first check the full name
    theVar=(RooAbsArg*)_rarOVars.FindObject(getFullVarName(varName));
    if (theVar) return theVar;
    // then check for varName as full name
    theVar=(RooAbsArg*)_rarOVars.FindObject(varName);
    if (theVar) return theVar;
  }
  // for _rarVars
  // first check the full name
  theVar=(RooAbsArg*)_rarVars.FindObject(getFullVarName(varName));
  if (theVar) return theVar;
  // then check for varName as full name
  theVar=(RooAbsArg*)_rarVars.FindObject(varName);
  if (theVar) return theVar;
  // now check if varName is actually pre-defined config
  // for example, mean, sigma, etc.
  // in this case, the AbsVar returned should be what varName refer to
  TString varConfigStr=readConfStr(varName, "notSet", getVarSec());
  if ("notSet"==varConfigStr) return theVar;
  TString refName=getFullVarName(varName, varName+" "+varConfigStr);
  if (varName==refName) return theVar; // safeguard for self refer
  if (refName==getFullVarName(varName)) return theVar;
  return getAbsVar(refName);
}

/// \brief Creates AbsVars according to config string
///
/// \param configName A String of AbsVar names
/// \param argCollA A RooAbsCollection to store the created RooAbsArg
/// \param argCollB A RooAbsCollection to store the created RooAbsArg
/// \return The last RooAbsArg created
///
/// This function creates AbsVars specified by a string
/// containing param config names.
RooAbsArg *rarConfig::createAbsVars(TString configName,
				      RooAbsCollection *argCollA,
				      RooAbsCollection *argCollB)
{
  RooAbsArg *theVar(0);
  // get number of vars
  TString varsStr=readConfStr(configName, "notSet", getVarSec());
  if ("notSet"==varsStr) return theVar;
  rarStrParser varsStrParser=varsStr;
  // is the first arg fullNamed
  Bool_t fullNamed(kFALSE);
  if ("fullNamed"==varsStrParser[0]) {
    varsStrParser.Remove();
    fullNamed=kTRUE;
  }
  Int_t nVar=varsStrParser.nArgs();
  if (nVar>0) {
    for (Int_t i=0; i<nVar; i++)
      readConfStr(varsStrParser[i], "notSet", getVarSec());
    cout<<"Variables defined with config \""<<configName<<"\""
	<<" in config file for "<<GetName();
    if (fullNamed) cout<<" (fullNamed)";
    cout<<":"<<endl;
    for (Int_t i=0; i<nVar; i++) {
      cout<<Form(" #%02d ", i)<<varsStrParser[i]<<" "
	  <<readConfStr(varsStrParser[i], "notSet", getVarSec())<<endl;
    }
  }
  for (Int_t i=0; i<nVar; i++) {
    if (!fullNamed) theVar=createAbsReal(varsStrParser[i], varsStrParser[i]);
    else { // created as if the fullName is specified
      TString varConfigStr=
	readConfStr(varsStrParser[i], varsStrParser[i], getVarSec());
      rarStrParser varConfigStrParser=varConfigStr;
      TString fullName=varsStrParser[i];
      TString option ="N "+fullName;
      if (isVarType(varConfigStrParser[0])) fullName+=" "+fullName;
      theVar=createAbsVar(fullName+" "+varConfigStr, option);
    }
    if (argCollA) argCollA->add(*theVar);
    if (argCollB) argCollB->add(*theVar);
  }
  
  return theVar;
}

/// \brief To set derived columns' limits
/// \param data The dataset to set derived column limits
/// \param setLimits To set addOns' limits
///
/// it sets addOns limits (for plotting)
void rarConfig::setColLimits(RooDataSet *data, Bool_t setLimits)
{
  if (!data) return;
  RooArgSet *addOns=getAddOnCols();
  if (!addOns) return;
  RooArgList addOnsList(*addOns);
  RooArgSet theRow(*data->get());
  Int_t nAddOns=addOnsList.getSize();
  for (Int_t i=0; i<nAddOns; i++) {
    TString theName=addOnsList[i].GetName();
    RooRealVar *theVar=dynamic_cast<RooRealVar*>(theRow.find(theName));
    if (!theVar) continue;
    RooRealVar *fVar=dynamic_cast<RooRealVar*>(_rarVars.FindObject(theName));
    if (!fVar) {
      cout<<" No fundamental var created for "<<theName<<endl;
      exit(-1);
    }
    if (setLimits) {
      Double_t theMin=fVar->getMin();
      Double_t theMax=fVar->getMax();
      Double_t dMin, dMax;
      data->getRange(*theVar, dMin, dMax);
      if (dMin<theMin) {
	cout<<" W A R N I N G ! ! !"<<endl
	    <<" The allowable min value ("<<theMin<<") for "<<theVar->GetName()
	    <<" is larger than the min value ("<<dMin<<")"
	    <<" in dataset "<<data->GetName()<<endl
	    <<" The min of "<<theVar->GetName()<<" is set to "<<dMin
	    <<" for "<<data->GetName()<<endl
	    <<" Please make sure the range for "<<theVar->GetName()
	    <<" is valid"<<endl
	    <<" Your results based on "<<data->GetName()
	    <<" may NOT be accurate!"<<endl
	    <<endl;
	theMin=dMin;
      }
      if (dMax>theMax) {
	cout<<" W A R N I N G ! ! !"<<endl
	    <<" The allowable max value ("<<theMax<<") for "<<theVar->GetName()
	    <<" is smaller than the max value ("<<dMax<<")"
	    <<" in dataset "<<data->GetName()<<endl
	    <<" The max of "<<theVar->GetName()<<" is set to "<<dMax
	    <<" for "<<data->GetName()<<endl
	    <<" Please make sure the range for "<<theVar->GetName()
	    <<" is valid"<<endl
	    <<" Your results based on "<<data->GetName()
	    <<" may NOT be accurate!"<<endl
	    <<endl;
	theMax=dMax;
      }
      theVar->setRange(theMin, theMax);
    } else
      theVar->removeRange();
    //cout<<" Derived obs "<<theName <<" set to ("<<theVar->getMin()<<", "
    //<<theVar->getMax()<<")"<<endl;
  }
}

/// \brief To add derived columns to dataset
/// \param data The dataset to add derived columns
/// \param addColmns To add addOns
/// \param setLimits To set addOns' limits
///
/// It adds derived columns to the dataset and
/// try to set the correct obs ranges for those derived obs
void rarConfig::addColumns(RooDataSet *data, Bool_t addColmns,
			   Bool_t setLimits)
{
  if (!data) return;
  RooArgSet *addOns=getAddOnCols();
  if (!addOns) return;
  // add the columns
  if (addColmns) {
    // data->addColumns(*addOns); // gives errors "can only add to owned list"
    TIterator* iter = addOns->createIterator();
    RooAbsArg *temp(0);
    while ( temp = (RooAbsArg*)(iter->Next()) ) {
      data->addColumn(*temp);
    }
  }
  // then set the limits for those obs if it is RRV
  setColLimits(data, setLimits);
}

/// \brief To create a dataset
///
/// \param dsStr The dataset config string
/// \param isUB Returned boolean for the ub status
/// \return The dataset created
///
/// It creates a \p RooDataSet object according to its config string,
/// \p dsStr. The config string has the form like
/// - ``<tt>\<dsName\> \<dsType\> "<dsTitle>" ...</tt>''
/// .
/// where \p dsName is the name of the dataset to be created,
/// and it is the name of config item in the config section.
/// The next token, \p dsType, is the dataset type, i.e.,
/// how the dataset will be created.
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_dsi">See doc for Dataset Input Section</a> for more details.
RooDataSet *rarConfig::createDataSet(const char *dsStr, Bool_t &isUB, TString wgtVarName)
{
  isUB=kFALSE;
  rarStrParser datasetStrParser=dsStr;
  TString myName=datasetStrParser[0];
  datasetStrParser.Remove();
  TString dsType=datasetStrParser[0];
  datasetStrParser.Remove();
  TString myTitle=datasetStrParser[0];
  datasetStrParser.Remove();
  
  // cout<<endl<<" For "<<myName<<endl;
  RooDataSet *data(0);
  if ("ascii"==dsType) { // dealing with ascii text file method
    RooDataSet *theData(0); // temporary data set
    if (1==datasetStrParser.nArgs()) { // #1
      theData=RooDataSet::read(datasetStrParser[0], *getPrimaryObs());
    } else if(2==datasetStrParser.nArgs()) { // #2, with options
      theData=RooDataSet::read(datasetStrParser[0], *getPrimaryObs(),
			    datasetStrParser[1]);
    } else if(3==datasetStrParser.nArgs()) { // #3, with multi text files
      theData=RooDataSet::read(datasetStrParser[0], *getPrimaryObs(),
			    datasetStrParser[1], datasetStrParser[2]);
    } else if(4==datasetStrParser.nArgs()) { // #4, with indexCatName
      theData=RooDataSet::read(datasetStrParser[0], *getPrimaryObs(),
			    datasetStrParser[1], datasetStrParser[2],
			    datasetStrParser[3]);
    } else {
      cout<<"Not implemented for "<<datasetStrParser.nArgs()<<" args"<<endl;
      exit(-1);
    }
    // if weight
    // RooArgSet primaryObs = *getPrimaryObs();
    // primaryObs.add(evtWgt); // add the weight to the list
    // data = new RooDataSet("theData", "theData", tree, primaryObs, cut, wgtVarName);
    // TFile f("/tmp/temporay.root");
    // data->Write();
    // f.ls();
    // f.Close();
    // RooDataSet *d = (RooDataSet*) f.FindObject("d");
    addColumns(theData);
    // copy and add weight variable
    data=new RooDataSet("asciiData", "asciiData", theData, *_fullObs, "1", wgtVarName);
  } else if ("root"==dsType) { // dealing with root file
    if (datasetStrParser.nArgs()<2) {
      cout<<"At least two more args needed"<<endl;
      exit(-1);
    }
    TString cut="";
    if (datasetStrParser.nArgs()>2) cut=datasetStrParser[2];
 
    // open the file
    TFile f(datasetStrParser[0]);

    TTree *tree = (TTree*) f.Get(datasetStrParser[1]);
    if (!tree) {
      cout<< "Can not find tree "<<datasetStrParser[1] <<" in file "<<datasetStrParser[0]<<endl;
      exit(-1);
    }
    
    // now create the dataset and add weight variable at same time
    RooArgSet primaryObs = *getPrimaryObs();
    if (wgtVarName != "") {
      RooRealVar evtWgt(wgtVarName, "weight", 1.0);
      primaryObs.add(evtWgt); // add the weight to the list
      data = new RooDataSet("theData", "theData", tree, primaryObs, cut, wgtVarName);
    } else {
      data = new RooDataSet("theData", "theData", tree, primaryObs, cut);
    }

    addColumns(data);
    //data->Print("v");
  } else if ("add"==dsType) { // dealing with add method
    // create dataset by adding existing datasets
    cout<<"RooDataSet::add: adding";
    for (Int_t i=0; i<datasetStrParser.nArgs(); i++) {
      cout<<" "<<datasetStrParser[i];
    }
    cout<<endl;
    data=new RooDataSet("addData", "add dataset", *_fullObs);
    Int_t nDS=datasetStrParser.nArgs()/2;
    for (Int_t i=0; i<nDS; i++) {
      // get dataset #i
      RooDataSet *theData=getData(datasetStrParser[0]);
      datasetStrParser.Remove();
      if (!theData) continue;
      Int_t srcNevt=(Int_t)theData->numEntries();
      Double_t nEvtScale=atof(datasetStrParser[0]);
      datasetStrParser.Remove();
      Int_t nEvt=(Int_t)nEvtScale;
      if (nEvtScale<1) { // if less than 1, it is scaling factor
	nEvt=(Int_t)(.5+nEvtScale*srcNevt);
      }
      if ((nEvt<=0) || (nEvt>srcNevt)) nEvt=srcNevt;
      if ((nEvt>0)&&(nEvt<srcNevt)) isUB=kTRUE;
      cout<<" Adding "<<nEvt<<" events from "<<theData->GetName()<<endl;
      // create index vector
      vector<Int_t> indexVector;
      for (Int_t j=0; j<srcNevt; j++) {
	indexVector.push_back(j);
      }
      for (Int_t j=0; j<nEvt; j++) {
	// random access to index to make sure the sample added
	// is not sequence dependent to the src data
	Int_t randomIdx=
	  RooRandom::randomGenerator()->Integer(indexVector.size());
	RooArgSet *theRow=(RooArgSet *)theData->get(indexVector[randomIdx]);
	data->add(*theRow);
	// remove the index so it won't be selected again
	vector<Int_t>::iterator the_iterator;
	the_iterator=indexVector.begin();
	the_iterator+=randomIdx;
	indexVector.erase(the_iterator);
      }
    }
    addColumns(data);
    cout<<"RooDataSet::add: added "<<data->numEntries()<<" events"<<endl;
  } else if ("reduce"==dsType) { // dealing with reduce method
    if (datasetStrParser.nArgs()<=0) {
      cout<<"need more args"<<endl;
      exit(-1);
    }
    TString srcDataName=datasetStrParser[0];
    datasetStrParser.Remove();
    RooDataSet *theData=getData(srcDataName);
    if (!theData) {
      cout<<"Can not find dataset named \""<<srcDataName<<"\""<<endl;
      exit(-1);
    }
    TString cuts="1";
    if (datasetStrParser.nArgs()>=1) {
      cuts=datasetStrParser[0];
      datasetStrParser.Remove();
    }
    data=(RooDataSet*)theData->reduce(cuts);

    cout<<"RooDataSet::reduce: reduced "<<data->numEntries()<<" events"<<endl;
 
  } else if ("hist"==dsType) { // dealing with RooDataHist
    if (datasetStrParser.nArgs()<=0) {
      cout<<"need more args"<<endl;
      exit(-1);
    }
    TString srcDataName=datasetStrParser[0];
    datasetStrParser.Remove();
    RooDataSet *theData=getData(srcDataName);
    if (!theData) {
      cout<<"Can not find dataset named \""<<srcDataName<<"\""<<endl;
      exit(-1);
    }
    data=(RooDataSet*)
      new RooDataHist("theData", "theData", *getPrimaryObs(), *theData);
  } else {
    cout<<"Not implemented for "<<datasetStrParser[1]<<" dataset type"<<endl;
    exit(-1);
  }
  
  // change its name and title
  data->SetName(myName);
  data->SetTitle(myTitle);
  
  return data;
}

/// \brief To compute correlation matrix for dataset
/// \param varList The vars to compute
/// \param data The dataset to compute
///
/// It computes the correlation matrix for given dataset.
void rarConfig::computeCorrelations(RooArgList varList,
				    const RooDataSet *data) {
  RooArgSet theFVars(*data->get());
  // remove any categories
  RooArgSet theVars;
  for (Int_t i=0; i<varList.getSize(); i++) {
    RooAbsArg *theVar=varList.at(i);
    if (_rarOVars.FindObject(theVar->GetName())) continue;
    if (getCats()->find(theVar->GetName())) continue;
    theVars.add(*(theFVars.find(theVar->GetName())));
  }
  // make sure there are vars for correlation matrix
  if (theVars.getSize()<=1) return;
  
  // the following codes are mainly from Q2BTimeMLFit::computeCorrelations
  unsigned i,j,k;
  unsigned nEvt = data->numEntries();
  unsigned nPar =  theVars.getSize();
  double* average = new double[nPar];
  RooAbsReal **arg = new RooAbsReal*[nPar];
  TIterator* iter = theVars.createIterator();
  i=0;
  RooAbsArg *y(0);
  while ( (y=(RooAbsArg*)(iter->Next())) ) {
    arg[i] =  dynamic_cast<RooAbsReal*>(y);
    average[i]=0;
    i++;
  }
  for (i=0;i<nEvt;++i) {
    data->get(i);
    for (unsigned j=0;j<nPar;++j)
      if (arg[j]!=0) average[j]+=(arg[j]->getVal())/nEvt;
  }
  double** cor = new double*[nPar];
  for (i=0;i<nPar;++i) {
    cor[i] = new double[nPar];
    for (unsigned j=0;j<nPar;++j) cor[i][j]=0;
  }
  for (i=0;i<nEvt;++i) {
    data->get(i);
    for (j=0;j<nPar;++j) for (k=0;k<nPar;++k) {
      if (arg[j]!=0)
	cor[j][k]+=(arg[j]->getVal()-average[j])*(arg[k]->getVal()-average[k]);
    }
  }
  for (j=0;j<nPar;++j) for (k=0;k<nPar;++k)  {
    if (j!=k && arg[j]!=0 && arg[k]!=0)
      cor[j][k]/=sqrt(cor[j][j]*cor[k][k]);
  }
  
  cout << endl << " Correlation matrix for " << data->GetName()
       << ":" << endl;
  cout << "        ";
  char chrbuf[9], label[9];
  for (i=0; i< nPar-1;i++) {
    strncpy(label, arg[i]->GetName(), 8);  label[8] = '\0';
    sprintf(chrbuf, "%8s", label);
    cout << "  " << chrbuf;
  }
  cout << endl;
  for (i=0;i<nPar;++i) {
    if (i>0) {
      strncpy(label, arg[i]->GetName(), 8);  label[8] = '\0';
      sprintf(chrbuf, "%8s", label);
      cout << chrbuf;
      for (j=0;j<i;++j) {
	sprintf(chrbuf, "%8.4f", cor[i][j]);
	cout << "  " << chrbuf;
      }
      cout << endl;
    }
    delete[] cor[i];
  }
  cout << endl;
  delete[] cor;
  delete[] average;
  delete[] arg;

  // FFW but only outside babar
  /*
  const TMatrixDSym* corn = data->correlationMatrix(theVars) ;
  //const TMatrixDSym* covn = data->covarianceMatrix(theVars) ;

  // Print correlation, covariance matrix
  cout << "correlation matrix for " << data->GetName() << endl ;
  corn.Print() ;
  //cout << "covariance matrix" << endl ;
  //covn.Print() ;
  */
}

/// \brief To create a RooRarFit Pdf object
///
/// \param configStr The config string
/// \return The Pdf created
///
/// It creates Pdf object according to its config string, \p configStr.
/// The config string has form like
/// \verbatim
/// <pdfName> <pdfType> ["<Optional pdf Title>"]\endverbatim
/// where \p \<pdfName\> is the name of the #rarBasePdf object to be created,
/// and its config section is [\<pdfName\> Config].
/// The second token is the type of the Pdf,
/// and at the end of the config string, it can have an optional
/// title to replace the default one. 
/// It first check if the pdf named \p pdfName has been created,
/// if yes, it just returns that pdf,
/// if not, it creates the pdf and returns it.
///
/// \par Add New Type of PDF
/// One has to do two things to make new type of PDF available.
/// First he needs to create a new class inherited from #rarBasePdf,
/// second he adds an entry here so the new class can be instantiated
/// through the standard creation mechanism in RooRarFit.
rarBasePdf *rarConfig::createPdf(const char *configStr)
{
  rarStrParser configStrParser=configStr;
  if(configStrParser.nArgs()<1) {
    cout<<"\""<<configStr<<"\""<<endl
	<<" At least 1 arg needed for rarConfig::createPdf"<<endl;
    exit(-1);
  }
  rarBasePdf *thePdf(0);
  // get the name of the pdf
  TString theName=configStrParser[0];
  configStrParser.Remove();
  // check if the pdf has already been created
  if (thePdf=(rarBasePdf*)_rarPdfs.FindObject(theName)) return thePdf;
  TString theConfigSec=theName+" Config";
  // get the pdf type
  if(configStrParser.nArgs()<1) {
    cout<<" You have to specify the pdf type"<<endl;
    exit(-1);
  }
  TString thePdfType=configStrParser[0];
  configStrParser.Remove();
  // get the pdf title
  TString theTitle(thePdfType+" "+GetTitle()); // the default title
  if(configStrParser.nArgs()>0) {
    theTitle=configStrParser[0]; // this is the title
    configStrParser.Remove();
  }
  const char *dummy("");
  RooDataSet *_theData=getData(dummy);
  rarDatasets *_datasets=getDatasets();
  // create different types of PDFs
  if ("MLPdf"==thePdfType) {
    if (getPdfType()!="MLFitter") {
      cout<<"Only the final mlFitter can create "<<thePdfType<<endl;
      exit(-1);
    }
    thePdf=new rarMLPdf(_configFile, theConfigSec, configStr,
			_datasets, _theData, theName, theTitle);
  } else if ("ProdPdf"==thePdfType) {
    thePdf=new rarProd(_configFile, theConfigSec, configStr,
		       _datasets, _theData, theName, theTitle);
  } else if (("AddPdf"==thePdfType)||("AddModel"==thePdfType)) {
    thePdf=new rarAdd(_configFile, theConfigSec, configStr,
		      _datasets, _theData, theName, theTitle);
  } else if ("Simultaneous"==thePdfType) {
    thePdf=new rarSimPdf(_configFile, theConfigSec, configStr,
			 _datasets, _theData, theName, theTitle);
  } else if ("Exp"==thePdfType) {
    thePdf=new rarExp(_configFile, theConfigSec, configStr,
		      _datasets, _theData, theName, theTitle);
  } else if (("Gaussian"==thePdfType)||("BreitWigner"==thePdfType) ||
	     ("Landau"==thePdfType)) {
    thePdf=new rarGaussian(_configFile, theConfigSec, configStr,
			   _datasets, _theData, theName, theTitle);
  } else if ("TwoGauss"==thePdfType) {
    thePdf=new rarTwoGauss(_configFile, theConfigSec, configStr,
			   _datasets, _theData, theName, theTitle);
  } else if ("OsipDisc"==thePdfType) {
    thePdf=new rarOsipDisc(_configFile, theConfigSec, configStr,
			   _datasets, _theData, theName, theTitle);
  } else if (("TriGauss"==thePdfType)||
	     ("TriGaussModel"==thePdfType)||("GexpShape"==thePdfType)) {
    thePdf=new rarTriGauss(_configFile, theConfigSec, configStr,
			   _datasets, _theData, theName, theTitle);
  } else if (("BifurGauss"==thePdfType)||("BGGauss"==thePdfType)) {
    thePdf=new rarBifurGauss(_configFile, theConfigSec, configStr,
			     _datasets, _theData, theName, theTitle);
  } else if ("Novosibirsk"==thePdfType) {
    thePdf=new rarNovosibirsk(_configFile, theConfigSec, configStr,
			      _datasets, _theData, theName, theTitle);
  } else if ("Ballack"==thePdfType) {
    thePdf=new rarBallack(_configFile, theConfigSec, configStr,
			  _datasets, _theData, theName, theTitle);
  } else if ("CBShape"==thePdfType) {
    thePdf=new rarCBShape(_configFile, theConfigSec, configStr,
			  _datasets, _theData, theName, theTitle);
  } else if (("Polynomial"==thePdfType)||("Chebychev"==thePdfType)) {
    thePdf=new rarPoly(_configFile, theConfigSec, configStr,
		       _datasets, _theData, theName, theTitle);
  } else if ("ArgusBG"==thePdfType) {
    thePdf=new rarArgusBG(_configFile, theConfigSec, configStr,
			  _datasets, _theData, theName, theTitle);
  } else if ("Step"==thePdfType) {
    thePdf=new rarStep(_configFile, theConfigSec, configStr,
		       _datasets, _theData, theName, theTitle);
  } else if ("Binned"==thePdfType) {
    thePdf=new rarBinned(_configFile, theConfigSec, configStr,
			 _datasets, _theData, theName, theTitle);
  } else if (("Keys"==thePdfType)||("2DKeys"==thePdfType)) {
    thePdf=new rarKeys(_configFile, theConfigSec, configStr,
		       _datasets, _theData, theName, theTitle);
  } else if ("Generic"==thePdfType) {
    thePdf=new rarGeneric(_configFile, theConfigSec, configStr,
			  _datasets, _theData, theName, theTitle);
  } else if ("HistPdf"==thePdfType) {
    thePdf=new rarHistPdf(_configFile, theConfigSec, configStr,
			  _datasets, _theData, theName, theTitle);
  } else if ("GaussModel"==thePdfType) {
    thePdf=new rarGaussModel(_configFile, theConfigSec, configStr,
			     _datasets, _theData, theName, theTitle);
  } else if (("BCPGenDecay"==thePdfType)||("BDecay"==thePdfType)||
             ("Decay"==thePdfType)) {
    thePdf=new rarDecay(_configFile, theConfigSec, configStr,
			_datasets, _theData, theName, theTitle);
  } else if ("Cruijff"==thePdfType) {
    thePdf=new rarCruijff(_configFile, theConfigSec, configStr,
			  _datasets, _theData, theName, theTitle);
  } else if ("MultPdf"==thePdfType) {
    thePdf=new rarMultPdf(_configFile, theConfigSec, configStr,
			 _datasets, _theData, theName, theTitle);
  } else if ("UsrPdf"==thePdfType) {
    thePdf=new rarUsrPdf(_configFile, theConfigSec, configStr,
			 _datasets, _theData, theName, theTitle);
  } else if ("RelBreitWigner"==thePdfType) {
    thePdf=new rarRelBreitWigner(_configFile, theConfigSec, configStr,
                         _datasets, _theData, theName, theTitle);
  } else if ("Lass"==thePdfType) {
    thePdf=new rarLass(_configFile, theConfigSec, configStr,
    		       _datasets, _theData, theName, theTitle);
  } else if ("Voigtian"==thePdfType) {
    thePdf=new rarVoigtian(_configFile, theConfigSec, configStr,
    		       _datasets, _theData, theName, theTitle);
  } else if ("GounarisSakurai"==thePdfType) {
    thePdf=new rarGounarisSakurai(_configFile, theConfigSec, configStr,
				  _datasets, _theData, theName, theTitle);
  } else if ("Flatte"==thePdfType) {
    thePdf=new rarFlatte(_configFile, theConfigSec, configStr,
			 _datasets, _theData, theName, theTitle);
  } else if ("Uniform"==thePdfType) {
    thePdf=new rarUniform(_configFile, theConfigSec, configStr,
			 _datasets, _theData, theName, theTitle);
  } else if ("Threshold"==thePdfType) {
    thePdf=new rarThreshold(_configFile, theConfigSec, configStr,
			 _datasets, _theData, theName, theTitle);
  } else {
    cout<<"\""<<configStr<<"\""<<endl
	<<"Pdf \""<<thePdfType<<"\" not implemented"<<endl;
    exit(-1);
  }
  _rarPdfs.Add(thePdf);
  
  cout<<thePdfType<<" created"<<endl<<endl;
  return thePdf;
}
