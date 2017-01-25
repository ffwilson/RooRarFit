/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarSPlot.cc,v 1.13 2014/09/14 17:33:54 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides sPlot class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides sPlot class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "rarSPlot.hh"
#include "RooAbsPdf.h"

#include "TMatrixD.h"
using std::cout;
using std::endl;

ClassImp(rarSPlot) ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarSPlot::rarSPlot() :
  TNamed(),
  _obs(0), _fitRes(0),
  _nComps(0), _nComp0s(0), _idxMap(0), _covM(0,0), _iV(0,0), _V(0,0),
  _sHist(0), _sData(0),
  _verbose(kFALSE)
{
  // Default constructor
}

/// \brief Default ctor
/// \param name The name
/// \param title The title
/// \param obs The obs to get sPlot
/// \param data The dataset to get sPlot from
/// \param fitRes The fit result
/// \param pdfs The component pdfs
/// \param yields The component yields
/// \param pdf0s The component pdfs w/ fixed yields
/// \param yield0s The fixed component yields
/// \param projDeps Obs to be ignored for normalization
/// \param verbose Boolean to control debug output
///
/// Default ctor to set several common data members,
/// and then call #init.
rarSPlot::rarSPlot(const char* name, const char *title, RooAbsReal *obs,
                   RooDataSet &data, RooFitResult *fitRes,
                   const RooArgList &pdfs, const RooArgList &yields,
                   const RooArgList &pdf0s, const RooArgList &yield0s,
                   const RooArgSet &projDeps, const Bool_t verbose) :
  TNamed(name, title),
  _obs(obs), _data(data), _nEvts(_data.numEntries()),
  _fitRes(fitRes), _fitPars(fitRes->floatParsFinal()),
  _pdfs(pdfs), _yields(yields), _pdf0s(pdf0s), _yield0s(yield0s),
  _projDeps(projDeps), _normVars(*pdfs[0].getDependents(&data)),
  _nComps(0), _nComp0s(0), _idxMap(0), _covM(0,0), _iV(0,0), _V(0,0),
  _sHist(0), _sData(0),
  _verbose(verbose)
{
  init();
}

/// \brief Trivial dtor
rarSPlot::~rarSPlot()
{
}

/// \brief Initial function called by ctor
///
/// It will initialize sPlot Matrics, etc. based on input
void rarSPlot::init()
{
  
  if (_verbose) {
    cout<<" rarSPlot:  yields:"<<endl;
    _yields.Print();
    cout<<" pdfs:"<<endl;
    _pdfs.Print();
    cout<<" data:"<<endl;
    _data.Print();
    //_data.Print("v");
    cout<<endl;
  }
  
  _nComps=_yields.getSize();
  if (_nComps<=0) {
    cout<<" No comp found"<<endl;
    _nComps=0;
    return;
  }
  if (_pdfs.getSize()!=_nComps) {
    cout<<" The size of yields and pdfs do not match"<<endl;
    _nComps=0;
    return;
  }
  if (_fitPars.getSize()!=_nComps) {
    cout<<" The size of yields and floating params do not match"<<endl;
    _yields.Print();
    _fitPars.Print();
    _nComps=0;
    return;
  }
  _nComp0s=_yield0s.getSize();
  if (_pdf0s.getSize()!=_nComp0s) {
    cout<<" The size of yield0s and pdf0s do not match"<<endl;
    _nComps=0;
    return;
  }
  if ((_nComp0s!=0)&&_verbose) {
    cout<<" yield0s:"<<endl;
    _yield0s.Print();
    cout<<" pdf0s:"<<endl;
    _pdf0s.Print();
    cout<<endl;
  }

  // Set the right size
  _idxMap.Set(_nComps);
  _covM.ResizeTo(_nComps, _nComps);
  _iV.ResizeTo(_nComps, _nComps);
  _V.ResizeTo(_nComps, _nComps);
  
  // Remember, the order of finalPars is not necessarily the same
  // as in #_yields, and we'd better map between them
  for (Int_t i=0; i<_nComps; i++) {
    RooAbsArg *theYield(0);
    if(theYield=_yields.find(_fitPars[i].GetName()))
      _idxMap[i]=_yields.index(theYield);
    else {
      cout<<" Can not find "<<_fitPars[i].GetName()
          <<" in yields:"<<endl;
      _yields.Print();
      _nComps=0;
      return;
    }
    if (_verbose)
      cout<<" fitPar #"<<i<<" --> yield #"<<_idxMap[i]
          <<" "<<_fitPars[i].GetName()<<endl;
  }
  if (_verbose) cout<<endl;
  
  // normalization ignored obs
  // The list of variables to normalize over when calculating PDF values.
  _normVars.remove(_projDeps, kTRUE, kTRUE);
  if (_verbose) {
    cout<<" Normalization variables" << endl;
    _normVars.Print();
    cout<<" Norm ignored variables"<<endl;
    _projDeps.Print();
    cout<<endl;
  }

  // create full snapshot
  RooArgList *pdfs=(RooArgList*)_pdfs.snapshot();
  RooArgList *yields=(RooArgList*)_yields.snapshot();
  RooArgList *pdf0s=(RooArgList*)_pdf0s.snapshot();
  RooArgList *yield0s=(RooArgList*)_yield0s.snapshot();
  _pdfs.removeAll();
  _pdfs.add(*pdfs);
  _yields.removeAll();
  _yields.add(*yields);
  _pdf0s.removeAll();
  _pdf0s.add(*pdf0s);
  _yield0s.removeAll();
  _yield0s.add(*yield0s);
  // attach dataset
  for (Int_t i = 0; i<_nComps; i++) {
    RooAbsPdf *pdf = (RooAbsPdf*)_pdfs.at(_idxMap[i]);
    pdf->attachDataSet(_data);
  }
  for (Int_t i = 0; i<_nComp0s; i++) {
    RooAbsPdf *pdf = (RooAbsPdf*)_pdf0s.at(i);
    pdf->attachDataSet(_data);
  }
  
  // Loop over all the parameters to make the covariance matrix.
  for (Int_t i=0; i<_nComps; i++) {
    for (Int_t j=0; j<_nComps; j++) {
      const RooRealVar *rowVar= (const RooRealVar*)_fitPars.at(i);
      const RooRealVar *colVar= (const RooRealVar*)_fitPars.at(j);
      assert(0 != rowVar && 0 != colVar);
      _covM(i,j) = rowVar->getError()*colVar->getError()
	*_fitRes->correlation(rowVar->GetName(),colVar->GetName());
      _V(i,j)=_covM(i,j);
    }
  }
  
  // Get the inverse of the covariance matrix
  // First check if it's singular
  if (_covM.Determinant() == 0) {
    cout << "rarSPlot Error: covariance matrix is singular;"
	 << " I can't invert it!" << endl;
    _covM.Print();
    _nComps=0;
    return;
  }
  _iV=TMatrixD(TMatrixD::kInverted, _V);
  
  if (_verbose) {
    cout << "Covariance matrix:" << endl;
    _V.Print();
    cout << "Inverse of covariance matrix:" << endl;
    _iV.Print();
  }
  
  // calculate all denominators
  _dens.Set(_nEvts);
  cout<<" calculate denominators ";
  for (Int_t ievt=0; ievt<_nEvts; ievt++) {
    if (_verbose&&ievt%1000==0) {
      cout << ".";
      cout.flush();
    }
    _data.get(ievt);
    _dens[ievt]=0;
    for (Int_t k=0; k<_nComps; k++) {
      _dens[ievt]+=((RooAbsReal*)_yields.at(_idxMap[k]))->getVal()*
	((RooAbsPdf*)_pdfs.at(_idxMap[k]))->getVal(&_normVars);
      if (0==((RooAbsPdf*)_pdfs.at(_idxMap[k]))->getVal(&_normVars)) {
        _data.get(ievt);
        _dens[ievt]+=((RooAbsReal*)_yields.at(_idxMap[k]))->getVal()*
          ((RooAbsPdf*)_pdfs.at(_idxMap[k]))->getVal(&_normVars);
      }
    }
    for (Int_t k=0; k<_nComp0s; k++) {
      _dens[ievt]+=((RooAbsReal*)_yield0s.at(k))->getVal()*
	((RooAbsPdf*)_pdf0s.at(k))->getVal(&_normVars);
    }
  }
  cout<<endl;
  
  // init _sPns
  _sPns.Set((_nComps+1)*_nEvts);
  _sPnb.Set(_nComps+1);
  _sPnb.Reset();
  
  if (_nComp0s>0) {
    // now using eqn 15 in BAD 509 to calculate #_iV
    for (Int_t i=0; i<_nComps; i++) {
      RooAbsPdf *fi=(RooAbsPdf*)_pdfs.at(_idxMap[i]);
      for (Int_t j=i; j<_nComps; j++) {
	RooAbsPdf *fj=(RooAbsPdf*)_pdfs.at(_idxMap[j]);
	Double_t iVij(0);
	for(Int_t ievt=0; ievt<_nEvts; ievt++) {
	  _data.get(ievt);
	  iVij+=fi->getVal(&_normVars)*fj->getVal(&_normVars)
	    /_dens[ievt]/_dens[ievt];
	}
	_iV(i,j)=iVij;
      }
    }
    for (Int_t i=0; i<_nComps; i++) {
      for (Int_t j=0; j<i; j++) {
	_iV(i,j)=_iV(j,i);
      }
    }
    
    _V=TMatrixD(TMatrixD::kInverted, _iV);
    
    if (_verbose) {
      cout << "Covariance matrix:" << endl;
      _V.Print();
      cout << "Inverse of covariance matrix:" << endl;
      _iV.Print();
    }
    
    // get sP0
    fillsPn(_nComps);
  }
}

/// \brief Fill sPn array
/// \param compIdx Component index
void rarSPlot::fillsPn(Int_t compIdx)
{
  if ((compIdx>_nComps)||(compIdx<0)) {
    cout<<" Invalid comp index "<<compIdx<<endl;
    return;
  }
  // check if done already
  if (_sPnb[compIdx]) return;
  if (compIdx==_nComps) {
    // fill sPn's
    for (Int_t i=0; i<_nComps; i++) {
      fillsPn(i);
    }
    // fill sP0
    cout<<" fill sP0 ";
    for (Int_t i=0; i<_nEvts; i++) {
      if (_verbose&&i%1000==0) {
        cout << ".";
        cout.flush();
      }
      _sPns[compIdx*_nEvts+i]=1;
      for (Int_t j=0; j<_nComps; j++) {
	_sPns[compIdx*_nEvts+i]-=_sPns[j*_nEvts+i];
      }
    }
  } else {
    // fill sPn
    cout<<" fill sPn for fitPar #"<<compIdx
        <<" "<<_fitPars[compIdx].GetName()<<" ";
    for (Int_t i=0; i<_nEvts; i++) {
      if (_verbose&&i%1000==0) {
        cout << ".";
        cout.flush();
      }
      _data.get(i);
      Double_t numerator=0.;
      for (Int_t j=0; j<_nComps; j++) {
        numerator+=_V(compIdx,j)*
          ((RooAbsPdf*)_pdfs.at(_idxMap[j]))->getVal(&_normVars);
      }
      _sPns[compIdx*_nEvts+i]=numerator/_dens[i];
    }
  }
  cout<<endl;
  
  _sPnb[compIdx]=1;
  return;
}

/// \brief Fill the sPlot histogram and dataset for a yield
/// \param yield The yield to get sPlot for
/// \param nbins Number of bins for histogram
/// \param min Histogram min
/// \param max Histogram max
/// \param doErrors Use Sumw2() for histogram
/// \return The sPlot dataset created
///
/// It fills the sPlot histogram and dataset
RooDataSet* rarSPlot::fill(RooAbsReal &yield, Int_t nbins,
                           Double_t min, Double_t max, Bool_t doErrors)
{
  if (_nComps<=0) {
    _sHist=0;
    _sData=0;
    return 0;
  }
  
  // create sPlot histogram
  TString histName=Form("hist_%s_%s",GetName(), yield.GetName());
  _sHist=new TH1F(histName, histName, nbins, min, max);
  // Size of bins...
  Double_t binSize = (max-min)/nbins;
  if (_verbose) {
    cout <<" SPlot histogram bins: " << min << " to " << max
	 <<" with " << nbins << " bins."
	 <<" binSize = " << binSize << endl;
  }
  
  Int_t istar=_fitPars.index(_fitPars.find(yield.GetName()));
  // get sPns
  fillsPn(istar);
  if (_verbose) {
    Double_t sum = 0;
    for (Int_t j=0; j<_nComps; j++) {
      sum += _V(j,istar);
    }
    cout <<" Sum of star column in V: " << sum << endl
         <<" fitted yield "<<yield.GetName()<<" = "<<yield.getVal()<<endl;
  }
  // add a weight column to the dataset
  RooArgSet obs = *_data.get();
  RooRealVar evtWgt("evtWgt", "Weight for sPlot entry", 1.);
  obs.add(evtWgt); // add weight to list of observables
  _sData = new RooDataSet("sPlotData", "sPlotData", obs, "evtWgt");
  
  Double_t sum = 0;
  Double_t sumtmp = 0;
  // This forces the error in a bin to be calculated
  // as the sqrt of the sum of the squares of weights, as they should be.
  if (doErrors) _sHist->Sumw2();
  
  Double_t Nstar= yield.getVal();
  Double_t Vstar=0;
  for (Int_t j=0; j<_nComps; j++) {
    Vstar+=_V(istar, j);
  }
  Double_t Vsum=0;
  for(Int_t i=0; i<_nComps; i++) {
    for (Int_t j=0; j<_nComps; j++) {
      Vsum+=_V(i,j);
    }
  }
  Double_t Nratio=(Nstar-Vstar)/(_nEvts-Vsum);
  Double_t Nratio2=(Nstar-Vstar)*Nratio;
  if (_nComp0s>0) cout<<" Nratio = "<<Nratio<<" \tNratio2 = "<<Nratio2<<endl;
  cout<<" loop over dataset to fill sPlot ";
  //_data.write("/tmp/sd.txt");
  for (Int_t ievt=0; ievt<_nEvts; ievt++) {
    if (_verbose&&ievt%1000==0) {
      cout << ".";
      cout.flush();
    }
    // Read this event and find the value of x for this event.
    const RooArgSet *row=_data.get(ievt);
    Double_t obsVal = ((RooAbsReal*)row->find(_obs->GetName()))->getVal();
    Double_t p=1/Nstar*_sPns[istar*_nEvts+ievt];
    if (_nComp0s>0) {
      p=1/Nstar*(_sPns[istar*_nEvts+ievt]+Nratio*_sPns[_nComps*_nEvts+ievt]);
    }
    sumtmp += ((RooAbsPdf*)_pdfs.at(istar))->getVal(&_normVars)/_dens[ievt];
    
    //if (obsVal > min && obsVal < max)
    sum += p*Nstar;
    
    // Add weighted event to the sPlot histogram
    _sHist->Fill(obsVal, p*Nstar);
    // ... and to the dataSet
    _sData->add(*row, p*Nstar);
    
  }// end event loop
  cout<<endl;
  
  _sHist->SetEntries(sum*Double_t(binSize));
  
  if (_verbose)
    cout << " Entries should be: " << sum*Double_t(binSize)
	 << " (" << sum << " events)" << endl
	 << " Sum of likelihood ratios for nstar: " << sumtmp << endl;
  
  return _sData;
}
