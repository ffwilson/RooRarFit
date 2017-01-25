/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarMLPdf.rdl,v 1.7 2014/09/14 17:33:53 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_MLPDF
#define RAR_MLPDF

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarAdd.hh"

/// \brief Simple mlFit model builder.
///
/// Build a simple mlFit model which is basically an extended AddPdf.
/// This simple model can be used to build SimPdf for the final mlFit model,
/// or if no SimPdf is used at all, the first such simple model can be used
/// as the final model (you can build many such simple model in
/// #rarMLFitter anyway).
///
/// \par Config Directives:
/// \verbatim
/// configStr = MLPdf ["<Optional Title>"]
/// Comps = <name1> <name2>  ... <nameN>
/// Coeffs = <coeff1> <coeff2> ... <coeffN>
/// <coeff1> = AbsReal Def
/// ...
/// <coeffN> = AbsReal Def\endverbatim
/// \p Comps defines names of all its components.
/// Each component has its own config section, namely,
/// comp \p \<name1\> is configured in config section named
/// \p "<name1> Config".
/// \p Coeffs defines names of all its coefficients.
/// All the coeffs can be defined as \p RooRealVar or \p RooFormulaVar.
/// The number of coeffs should match that of components
/// since it is an extended AddPdf.
class rarMLPdf : public rarAdd {
  
public:
  rarMLPdf();
  rarMLPdf(const char *configFile, const char *configSec, const char*configStr,
	   rarDatasets *theDatasets, RooDataSet *theData,
	   const char *name, const char *title);
  virtual ~rarMLPdf();
  
  /// return sCoeff List
  virtual RooArgList getSCoeffList() {return _sCoeffs;}
  virtual RooArgSet getSpecialSet(TString setName="specialSet");
  
  /// \brief Return the special config string for splitting
  /// \return The special config string for splitting
  virtual TString getSpecialStr(){return _specialStr;}
  virtual RooAbsPdf *getProtGen();
  
  virtual RooPlot *doPdfPlot(TList &plotList, TString pdfList="");
  
protected:
  void init();
  
  RooArgList _sCoeffs; ///< Coeff before splitting
  RooArgSet _specialSet; ///< Special ArgSet for splitting
  RooArgSet _cat1Set; ///< ArgSet for the 1st type of each cat
  RooArgSet _asymSet; ///< ArgSet for two-state cat
  TString _specialStr; ///< Special config string for splitting
  
private:
  rarMLPdf(const rarMLPdf&);
  ClassDef(rarMLPdf, 0) // RooRarFit ML model class
    ;
};

#endif
