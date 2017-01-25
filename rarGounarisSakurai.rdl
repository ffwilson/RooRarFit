/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: RooRarFit
 *    File: $Id: rarGounarisSakurai.rdl,v 1.5 2014/09/14 17:33:52 fwilson Exp $
 * Authors: Fergus Wilson
 * History:
 *
 * Copyright (C) 2005-2012, RAL
 *****************************************************************************/
#ifndef RAR_GOUNARISSAKURAI
#define RAR_GOUNARISSAKURAI

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief GounarisSakurai rho line shape PDF builder
///
/// Build <a href="http://roofit.sourceforge.net/docs/classref/RooGounarisSakurai.html"
/// target=_blank>RooGounarisSakurai</a> Pdf.
/// \par Config Directives:
/// \verbatim
/// configStr = GounarisSakurai ["<Optional Title>"]
/// spin = <0|1|2>
/// x = AbsReal Def
/// mean  = AbsReal Def
/// width = AbsReal Def
/// spin  = AbsReal Def
/// radius = AbsReal Def
/// mass_a = AbsReal Def
/// mass_b = AbsReal Def
/// \endverbatim
/// \p x is the default observable.
/// \p mean is the peak position of the pdf.
/// \p width is the width of the pdf.
/// \p spin is the spin (= 0, 1, 2) (default 1).
/// \p radius is the form factor radius (default 3.1/GeV)
/// \p mass_a is the mass of the first daughter (default pi+ mass)
/// \p mass_b is the mass of the second daughter (default pi- mass)
/// After being defined, the two daughter masses and the spin are held constant
/// in the fit.
/// All the variables can be \p RooRealVar or \p RooFormulaVar.
class rarGounarisSakurai : public rarBasePdf {

public:
  rarGounarisSakurai();
  rarGounarisSakurai(const char *configFile, const char *configSec, 
		    const char *configStr,
		    rarDatasets *theDatasets, RooDataSet *theData,
		    const char *name, const char *title);
  virtual ~rarGounarisSakurai();
  
protected:
  void init();
  
  RooAbsReal *_x;     ///< Default obs (mass)
  RooAbsReal *_mean;  ///< Peak position of the PDF
  RooAbsReal *_width; ///< Width of the PDF
  RooAbsReal *_spin;  ///< Spin (0,1, or 2)
  RooAbsReal *_radius; ///< Form factor radius
  RooAbsReal *_mass_a; ///< Mass of daughter A
  RooAbsReal *_mass_b; ///< Mass of daughter B
  
private:
  rarGounarisSakurai(const rarGounarisSakurai&);
 
  ClassDef(rarGounarisSakurai, 0) // RooRarFit GounarisSakurai PDF class
    ;
};

#endif
