/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarRelBreitWigner.rdl,v 1.5 2014/09/14 17:33:54 fwilson Exp $
 * Authors: Fergus Wilson
 * History:
 * 
 * Copyright (C) 2005-2012, RAL
 *****************************************************************************/
#ifndef RAR_RELBREITWIGNER
#define RAR_RELBREITWIGNER

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief Relativistic Breit Wigner PDF builder
///
/// Build <a href="http://roofit.sourceforge.net/docs/classref/RooRelBreitWigner.html"
/// target=_blank>RooRelBreitWigner</a> Pdf.
/// \par Config Directives:
/// \verbatim
/// configStr = RelBreitWigner ["<Optional Title>"]
/// spin = <0|1|2>
/// x = AbsReal Def
/// mean  = AbsReal Def
/// width = AbsReal Def
/// spin  = AbsReal Def
/// \endverbatim
/// \p x is the default observable.
/// \p mean is the peak position of the pdf.
/// \p width is the width of the pdf.
/// \p spin is the spin (= 0, 1, 2).
/// All the variables can be \p RooRealVar or \p RooFormulaVar.
class rarRelBreitWigner : public rarBasePdf {

public:
  rarRelBreitWigner();
  rarRelBreitWigner(const char *configFile, const char *configSec, 
		    const char *configStr,
		    rarDatasets *theDatasets, RooDataSet *theData,
		    const char *name, const char *title);
  virtual ~rarRelBreitWigner();
  
protected:
  void init();
  
  RooAbsReal *_x;     ///< Default obs (mass)
  RooAbsReal *_mean;  ///< Peak position of the PDF
  RooAbsReal *_width; ///< Width of the PDF
  RooAbsReal *_radius; ///< Width of the PDF
  RooAbsReal *_mass_a; ///< 
  RooAbsReal *_mass_b; ///< 
  RooAbsReal *_spin; ///< 

  //Int_t _parSpin;     ///< Spin (0,1, or 2)
  
private:
  rarRelBreitWigner(const rarRelBreitWigner&);
 
  ClassDef(rarRelBreitWigner, 0) // RooRarFit RelBreitWigner PDF class
    ;
};

#endif
