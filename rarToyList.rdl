/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarToyList.rdl,v 1.2 2011/08/26 17:54:18 fwilson Exp $
 * Authors: Fergus Wilson
 * History:
 *
 * Copyright (C) 2005-2012, RAL
 *****************************************************************************/
//
//
// A helper class to keep track of events used in Toys
//
//
#ifndef RARTOYLIST_HH
#define RARTOYLIST_HH

#include "Riostream.h"
#include "TString.h"

#include <map>

using namespace std;

class rarToyList {

public:

  rarToyList() {datasetName= "unknown"; datasetEvts = -1;}

  virtual ~rarToyList() {}

  inline void setDataset(TString name, Double_t value) {
    datasetName = name;
    datasetEvts = value;
  }

  inline void setInitial(TString parameter, Double_t value) {
    mapInitial[parameter] = value;
  }

  inline void setRequested(TString parameter, Double_t value) {
    mapRequested[parameter] = value;
  }

  inline void setFound(TString parameter, Double_t value) {
    mapFound[parameter] = value;
  }

  inline void setUsed(TString parameter, Double_t value) {
    mapUsed[parameter] = value;
  }

  inline void setMethod(TString parameter, TString value) {
    mapMethod[parameter] = value;
  }

  void print() const;

  inline void reset() {
    mapInitial.clear(); 
    mapRequested.clear();
    mapFound.clear(); 
    mapUsed.clear(); 
    mapMethod.clear();
  }

private:

  // " Ind     Observable  Initial Request   Found      Used Adj  Method"
  // mapping with observable as key
  map<TString, Double_t> mapInitial;   // #events in PdfAct constructor
  map<TString, Double_t> mapRequested; // #events requested from ToyAct
  map<TString, Double_t> mapFound;     // #events available from pdf/data
  map<TString, Double_t> mapUsed;      // #events finally used in Toy
  map<TString, TString>  mapMethod;    // method type pdf/data ...

  TString  datasetName; // events in the protodataset
  Double_t datasetEvts; // events in the protodataset

  ClassDef(rarToyList,0);

};

#endif
