/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarBinned.rdl,v 1.5 2014/09/14 17:33:49 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_BINNED
#define RAR_BINNED

#include "TList.h"
#include "TString.h"
#include "TObject.h"
#include "TArrayD.h"

#include "rarBasePdf.hh"

/// \brief RhhBinned PDF builder
///
/// Build
/// Pdf.
/// \par Config Directives:
/// \verbatim
/// configStr = Binned ["<Optional Title>"]
/// x = AbsReal Def
/// nBins = <binVal>
/// limits = <<binVal>+1 points>
/// H00 = AbsReal Def
/// ...
/// H<binVal-1> = AbsReal Def\endverbatim
/// \p x is the default observable.
/// \p nBins is the number of bins of the binned function.
/// \p limits is a set of \p \<binVal\>+1 Double_t's setting
/// the bounds of those bins.
/// \p H00 \p ... \p H\<binVal-1\> are \p \<binVal\>-1 free
/// parameters of the binned function.
/// All the \p AbsReal variables can be \p RooRealVar or \p RooFormulaVar.
class rarBinned : public rarBasePdf {
  
public:
  rarBinned();
  rarBinned(const char *configFile, const char *configSec, const char *configStr,
	  rarDatasets *theDatasets, RooDataSet *theData,
	  const char *name, const char *title);
  virtual ~rarBinned();
  
protected:
  void init();
  
  RooAbsReal *_x; ///< Default obs
  Int_t _nBins; ///< Number of bins
  TArrayD *_limits; ///< Limits of those bins
  
private:
  rarBinned(const rarBinned&);
  ClassDef(rarBinned, 0) // RooRarFit Binned Pdf class
    ;
};

#endif
