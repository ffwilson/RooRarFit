/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarToyList.cc,v 1.3 2014/09/14 17:33:55 fwilson Exp $
 * Authors: Fergus Wilson
 * History:
 * 
 * Copyright (C) 2005-2012, RAL
 *****************************************************************************/
//
// BEGIN_HTML
// This is a helper class for displaying Toy events numbers
// END_HTML
//

#include "Riostream.h"
#include <map>
#include "rarVersion.hh"
#include "rarToyList.hh"

using namespace std;

ClassImp(rarToyList);

//--------------------------------------------
void rarToyList::print() const {
  
  if (mapInitial.empty()) {
    cout << "No observables defined in Toy so nothing to print" <<endl;
    return;
  }

  typedef map<TString, Double_t> MapType;
  typedef map<TString, TString>  MapType2;
  
  Double_t nInit_tot(0), nReq_tot(0), nFound_tot(0), nUsed_tot(0);
  Int_t index(0);

  cout << endl;
  cout << "Summary of Toy. Dataset : " << datasetName << " with " << datasetEvts << " events " << endl;
  cout << " Ind     Observable  Initial  Request      Used Adj      Method" << endl;
  cout << " ---     ----------  -------  -------      ---- ---      ------" << endl;

  MapType::const_iterator it;
  for (it = mapInitial.begin(); it != mapInitial.end(); ++it) {
    TString obs = it->first;
    Double_t nInit = it->second;
    nInit_tot += nInit;

    // now look for other entries
    MapType::const_iterator iter1 = mapRequested.find(obs);
    Double_t nReq(0);
    if (iter1 != mapRequested.end()) {nReq = iter1->second;}
    nReq_tot += nReq;

    //
    MapType::const_iterator iter2 = mapFound.find(obs);
    Double_t nFound(0);
    if (iter2 != mapFound.end()) {nFound = iter2->second;}
    nFound_tot += nFound;

    //
    MapType::const_iterator iter3 = mapUsed.find(obs);
    Double_t nUsed(0);
    if (iter3 != mapUsed.end()) {nUsed = iter3->second;}
    nUsed_tot += nUsed;

    //
    MapType2::const_iterator iter4 = mapMethod.find(obs);
    TString method("unk");
    if (iter4 != mapMethod.end()) {method = iter4->second;}

    //
    TString adj("");
    if (nUsed != nReq) {adj = "*";}

    index++;
    TString line(Form("%4d %14s %8d %8d  %8d   %1s %11s", 
		      index, obs.Data(), 
		      (Int_t) nInit, (Int_t) nReq, (Int_t) nUsed, 
		      adj.Data(), method.Data()));

    cout << line << endl;
  }

  TString totals(Form("Totals: %20d %8d  %8d", (Int_t) nInit_tot, (Int_t) nReq_tot, (Int_t) nUsed_tot));
  cout << " ---     ----------  -------  -------      ---- ---      ------" << endl;
  cout << totals << endl;
  cout << endl;
  return;
}
