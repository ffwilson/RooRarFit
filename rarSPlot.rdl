/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarSPlot.rdl,v 1.9 2014/09/14 17:33:54 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

#ifndef RAR_SPLOT
#define RAR_SPLOT

#define MAXNSPEC 100

class RooAbsReal;

#include "TH1.h"
#include "TMatrixD.h"

#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooHist.h"
#include "RooPlot.h"

/// \brief Create an sPlot
///
/// This creates an sPlot.  See BAD 509 for more information on sPlots.
/// An sPlot gives us the distribution of some variable,
/// x in our data sample for a given species,
/// for instance the distribution of signal mES.
/// The result is similar to a likelihood projection plot,
/// but no cuts are made, so every event contributes to the distribution.
///
/// To use this class, you first must perform your fit twice:
/// The first time perform your nominal fit.
/// For the second fit, fix your parameters at the minimum,
/// float only your yields, and remove any PDFs correlated with
/// the variable of interest.
/// (In particular, if you want to make a signal mES sPlot,
/// fix everything but yields, and perform the fit without the mES PDFs.)
/// Be sure to save the RooFitResult.   Now you can use rarSPlot.

class rarSPlot : public TNamed {
public:
  rarSPlot();
  rarSPlot(const char* name, const char *title, RooAbsReal *obs,
           RooDataSet &data, RooFitResult *fitRes,
           const RooArgList &pdfs, const RooArgList &yields,
           const RooArgList &pdf0s=RooArgList(),
           const RooArgList &yield0s=RooArgList(),
           const RooArgSet &projDeps=RooArgSet(), const Bool_t verbose=kTRUE);
  virtual ~rarSPlot();
  
  /// \brief Set verbose mode
  /// \param verbose Verbose mode
  virtual void setVerbose(Bool_t verbose=kTRUE) {_verbose=verbose;}
  
  RooDataSet* fill(RooAbsReal &yield, Int_t nbins, Double_t min, Double_t max,
                   Bool_t doErrors=kTRUE);
  
  /// \brief Get sPlot histogram
  /// \return The sPlot histogram
  virtual TH1F *getSPlotHist() {return _sHist;}
  /// \brief Get sPlot dataset
  /// \return The sPlot histogram
  virtual RooDataSet *getSPlotData() {return _sData;}
protected:
  virtual void init();
  virtual void fillsPn(Int_t compIdx);
  
  RooAbsReal *_obs; ///< Obs to do sPlot
  RooDataSet _data; ///< Input dataset
  Int_t _nEvts; ///< Number of events in dataset
  RooFitResult *_fitRes; ///< Fit results
  RooArgList _fitPars; ///< Fit params
  RooArgList _pdfs; ///< Component pdfs
  RooArgList _yields; ///< Component yields
  RooArgList _pdf0s; ///< Component pdfs w/ fixed yields
  RooArgList _yield0s; ///< Fixed component yields
  RooArgSet _projDeps; /// Obs to be ignored for normalization
  RooArgSet _normVars; /// Obs to be used for normalization
  
  Int_t _nComps; ///< Number of components
  Int_t _nComp0s; ///< Number of fixed components
  TArrayI _idxMap; ///< Array for index mapping between #_yields and #_fitPars
  TMatrixD _covM; ///< Covariance matrix from Minuit
  TMatrixD _iV; ///< Inverse of covM or from direct calculation
  TMatrixD _V; ///< covM or from inverse of iV
  TArrayD _dens; ///< Denominators
  TArrayD _sPns; ///< Array for sPn
  TArrayI _sPnb; ///< Array for sPn fill bit
  
  TH1F *_sHist; ///< SPlot histogram
  RooDataSet *_sData; ///< SPlot dataset
  
  Bool_t _verbose; ///< Boolean to control debug output
  
private:
  ClassDef(rarSPlot,0)   // RooRarFit sPlot class
    ;
};

#endif
