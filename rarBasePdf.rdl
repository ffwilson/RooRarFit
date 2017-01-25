/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarBasePdf.rdl,v 1.34 2014/09/14 17:33:49 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_BASEPDF
#define RAR_BASEPDF

#define NCOLORS 11

#include <string>
using namespace std;

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "RooBinning.h"
#include "RooCategory.h"
#include "RooConstVar.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
#include "RooSimultaneous.h"

#include "rarConfig.hh"
#include "rarDatasets.hh"

class rarMLFitter;

/// \brief Base class for all RooRarFit PDF classes
///
/// It is the base class for all RooRarFit PDF classes.
/// It defines common data and functions for PDF
/// upon #RooConfig.
class rarBasePdf : public rarConfig {
  
public:
  rarBasePdf();
  rarBasePdf(const char *configFile,const char*configSec,const char*configStr,
	     rarDatasets *theDatasets, RooDataSet *theData,
	     const char *name, const char *title);
  virtual ~rarBasePdf();
  
  /// Return the name of param config section
  virtual TString getVarSec() {return _paramSec;}
  
  /// \brief Return primary observables in dataset files
  /// \return The primary observables defined
  virtual RooArgSet *getPrimaryObs() {return _datasets->getPrimaryObs();}
  
  /// \brief Return addon columns for datasets
  /// \return The addon columns defined
  virtual RooArgSet *getAddOnCols() {return _datasets->getAddOnCols();}
  
  /// \brief Return the default dataset
  /// \param name The name of the dataset, dummy parameter here
  /// \return The dataset returned
  ///
  /// Return the default dataset
  virtual RooDataSet *getData(const char * /*name=0 */) {return _theData;}
  
  /// \brief Return the default rarDatasets object
  /// \return The dataset holder
  virtual rarDatasets *getDatasets() const {return _datasets;}
  
  /// \brief Return the pdf type string
  /// \return The pdf type string
  virtual TString getPdfType() {return _pdfType;}
  
  /// \brief Return the created RooFit PDF
  /// \return  The default RooFit PDF
  virtual RooAbsPdf *getPdf() {return _thePdf;}
  virtual RooAbsPdf *getPdfWOvar(RooArgList ignoredObs);
  virtual RooAbsPdf *getDPdfWvar(RooRealVar *theVar);
  virtual RooAbsPdf *getSimPdf(RooSimultaneous *simPdf=0, RooAbsPdf *srcPdf=0);

  /// \brief Set final rarMLFitter
  /// \param theFitter The final rarMLFitter
  virtual void setFitter(rarMLFitter *theFitter) {_theFitter=theFitter;}
  /// \brief Return the fianl rarMLFitter
  /// \return The final rarMLFitter
  virtual rarMLFitter *getFitter() {return _theFitter;}
  
  /// \brief Set (total) SimPdf for this RooRarFit Pdf
  /// \param simPdf SimPdf to be set
  virtual void setSimPdf(RooSimultaneous *simPdf) {_theSimPdf=simPdf;}
  virtual void setCondObss(RooArgSet condObsSet);
  virtual void setFitData(RooDataSet *theData=0);
  virtual RooArgSet getParams();
  
  virtual void preAction();
  virtual RooArgSet getArgSet(TString paramNames, Bool_t useRead=kFALSE,
			      RooArgSet *fullSet=0);
  virtual void doPdfFit(TString pdfList="");
  virtual void attachDataSet(const RooAbsData &data);
  virtual Bool_t isNegativeValue();
  virtual RooPlot *doPdfPlot(TList &plotList, TString pdfList="");
  virtual RooArgSet getCorrCoeffs();
  virtual Double_t getCorrCoeff(const TString pn1, const TString pn2);
  virtual RooAbsPdf *getProtGen();
  virtual Bool_t protGenIsDummy() {return _theProtGen==_myDummyPdf;}
  
  virtual void setControlBit(TString controlBitStr, TString bitConfigStr="",
			     TString configSec="");
  virtual void setControlBits(TString controlBitsStr);
  virtual Bool_t getControlBit(TString controlBitStr);
  
  static Int_t getColor(Int_t i);
  
protected:
  virtual void init();
  virtual void addProtVars();
  virtual void addProtVars(TString configName, RooArgSet &protVars);
  /// Set the name of param config section
  virtual void setVarSec(TString paramSec) {_paramSec=paramSec;}
  virtual void addToParams(RooRealVar *theVar);
  virtual void addToObs(RooRealVar *theVar);
  virtual RooArgList *getFormulaArgs(rarStrParser fStrParser);
  virtual Double_t getFormulaVal(TString varStr);
  virtual rarBasePdf *createPdfs(TString Comps="Comps", TList *pdfList=0,
				 RooAbsCollection *PDFs=0, TString secName="");
  virtual Bool_t matchCatType(RooCatType *catN, RooCatType *catO);
  virtual void saveFracName(TString fracName);
  virtual Bool_t isFracName(TString fracName);
  
  virtual RooBinning *getRange(RooRealVar *theVar, TString rPrefix,
                               Double_t &min, Double_t &max,
                               const Char_t *sec=0, Int_t *nBins=0);
  virtual void saveCorrCoeffs(RooFitResult *fr);
  virtual Bool_t saveCorrCoeff(TString corrCoefName, Double_t corrCoef,
			       Bool_t saveTrivial=kFALSE);
  virtual TString getCorrCoefName(const TString pn1, const TString pn2) const;
  virtual void doXPdfFit(TString pdfList="");
  virtual RooPlot *doXPdfPlot(TList &plotList, TString pdfList="");
  virtual RooPlot *doParamsOnPlot(RooPlot* frame, RooArgSet *params=0,
				  Int_t sigDigits=2, Option_t *options="NELU",
				  Double_t xmin=0.65, Double_t xmax=0.99,
				  Double_t ymax=0.95);
  virtual RooPlot *doChi2OnPlot(RooPlot *frame);
  
  TString _pdfType; ///< Pdf type string
  
  RooArgSet _obsSet; ///< Observables directly for this pdf (no sub-pdfs)
  RooArgSet _fObsSet; ///< Full obs of this pdf (including sub-pdfs')
  
  rarDatasets *_datasets; ///< Datasets holder
  RooDataSet *_theData; ///< Default dataset
  RooArgSet _protDataVars; ///< protDataVars
  RooArgSet _conditionalObs; ///< Conditional observables
  RooArgSet _condObss; ///< Conditional observables for production

  RooAbsPdf *_thePdf; ///< Default RooFit pdf for this class
  RooArgList _subPdfs; ///< subPdf ArgList
  RooAbsPdf *_theProtGen; ///< Constructed prototype var generator
  RooArgList _protGenPdfs; ///< Extra pdfs as prototype var generator
  RooSimultaneous *_thisSimPdf; ///< SimPdf associated with #_thePdf
  RooSimultaneous *_thisSimPdfWOP; ///< SimPdf w/o protCats
  RooSimultaneous *_theSimPdf; ///< Final pdf model if it is a SimPdf
  Int_t _nxPdf; ///< Number of extra pdfs built (directly) within this class
  TList _xPdfList; ///< List of extra pdfs (directly) within this class
  RooAbsPdf *_myDummyPdf; ///< Dummy constant pdf with params
  RooArgSet _corrCoeffs; ///< Correlation coefficients
  static RooGenericPdf _dummyPdf; ///< Dummy constant pdf
  static RooConstVar _dummyExpEvt; ///< Dummy expected events 
  static RooExtendPdf _dummyExtPdf; ///< Dummy extended pdf (0 evt)
  static RooCategory _compCat; ///< Component category for generation
  static Int_t _rarColors[NCOLORS]; ///< Colors for plots
  static rarMLFitter *_theFitter; ///< The final fitter
  static TString _fracNames; ///< String of frac names

  TString _paramSec; ///< Param config section name
  RooArgSet _params; ///< Param List (all created RooRealVar params)
  RooArgList _coeffs; ///< Coeff List (directly created params)
  RooArgSet _xParams; ///< ArgSet of extra params within this class
  string _afterFitSaverStr; ///< String to save params after pdf fit
  
  TString _controlStr; ///< String for pdf control booleans
  
private:
  rarBasePdf(const rarBasePdf&);
  ClassDef(rarBasePdf, 0) // RooRarFit base Pdf class
    ;
};

#endif
