/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarDatasets.rdl,v 1.16 2014/09/14 17:33:50 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_DATASETS
#define RAR_DATASETS

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarConfig.hh"
#include "rarDatasetDef.hh"

/// \brief Dataset holder class
///
/// This class instantiates a #rarDatasetDef class for dataset definition,
/// reads in and holds all the datasets from ascii/root files.
/// It also holds datasets derived from those primary datasets.
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_dsi">See doc for dataset input section.</a>
class rarDatasets : public rarConfig {
  
public:
  rarDatasets();
  rarDatasets(const char*configFile,const char*configSec,const char*actionSec);
  virtual ~rarDatasets();
  
  /// \brief Return primary observables in dataset files
  /// \return The primary observables defined
  virtual RooArgSet *getPrimaryObs() {return _dsd->getPrimaryObs();}
  
  /// \brief Return addon columns for datasets
  /// \return The addon columns defined
  virtual RooArgSet *getAddOnCols() {return _dsd->getAddOnCols();}
  
  /// \brief Return full observables
  /// \return The full (fundamental) observables defined
  virtual RooArgSet *getFullFObs() {return _fullFObs;}
  
  virtual TString getDSName(TString name);
  virtual RooDataSet *getData(const char *name=0);
  
  /// \brief Get dataset list
  /// \return the dataset list
  virtual TList *getDatasetList() {return &_dataSets;}
  
  virtual TString ubStr(TString dsName, const char *ubStrVal=0);
  virtual Bool_t isBlind(TString dsName);
  
protected:
  void init();
  //virtual void setWeightVar(); // not needed in new versions of root
  virtual void tabulateDatasets(const char *dsName=0);
  
  TString getWeightVarName(TString datasetName);

  TString _actionSec; ///< Action config section name
  rarDatasetDef *_dsd; ///< Dataset definition object
  TList _dataSets; ///< Defined datasets
  RooArgSet *_fullFObs; ///< Full set of fundamental observables
  RooArgSet _UBs; ///< Unblind strings for datasets
  
private:
  rarDatasets(const rarDatasets&);
  ClassDef(rarDatasets, 0) // RooRarFit dataset class
    ;
};

#endif
