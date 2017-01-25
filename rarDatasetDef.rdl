/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarDatasetDef.rdl,v 1.6 2014/09/14 17:33:50 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_DATASETDEF
#define RAR_DATASETDEF

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarConfig.hh"

/// \brief Dataset definition class
///
/// This class defines the dataset structure for the fitter.
/// \par Config Directives:
/// \verbatim
/// Fields = <field01> <field02> ... <fieldN>
/// <field01> = <field01> <varType> ...
/// <field02> = <field02> <varType> ...
/// ...
/// <fieldN> = <fieldN> <varType> ...
/// 
/// AddOns = <addon01> ... <addonM>
/// <addon01> = <addon01> <derivedVarType> ...
/// ...
/// <addonM> = <addonM> <derivedVarType> ...\endverbatim
/// \p Fields is a string containing all the names of
/// the variables in data files.
/// The order of the variables should be the same as in the ascii data file.
/// \p \<varType\> can be \p RooRealVar, \p RooCategory,
/// or \p RooStringVar.
/// Optional \p AddOns specifies derived columns for a dataset after it is read
/// in from data file.
/// \p \<addon01\>, \p ..., \p \<addonM\> must
/// be derived from other variables declared in \p Fields.
/// Because all these variables declared here are globally accessible,
/// they should have unique and meaningful names,
/// and it is advisable to have their full names explicitly specified
/// in the config items.
class rarDatasetDef : public rarConfig {
  
public:
  rarDatasetDef();
  rarDatasetDef(const char *configFile, const char *configSec);
  virtual ~rarDatasetDef();
  
  /// \brief Return primary observables in dataset files
  /// \return The primary observables defined
  virtual RooArgSet *getPrimaryObs() {return _primaryObs;}
  
  /// \brief Return addon columns for datasets
  /// \return The addon columns defined
  virtual RooArgSet *getAddOnCols() {return _addonCols;}
  
  virtual RooArgList *getFormulaArgs(rarStrParser fStrParser);
  
  virtual void setVal(TString var, Double_t val);
  virtual void setVal(TString var, Int_t val) { setVal(var, (Double_t)val);}
  virtual void setVal(TString var, TString val);
  
protected:
  void init();
  
  RooArgSet *_primaryObs; ///< Primary obs in dataset file
  RooArgSet *_addonCols; ///< Addon columns for dataset
  
private:
  rarDatasetDef(const rarDatasetDef&);
  ClassDef(rarDatasetDef, 0) // RooRarFit dataset definition class
    ;
};

#endif
