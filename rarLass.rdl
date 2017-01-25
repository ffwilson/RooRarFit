/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarLass.rdl,v 1.4 2014/09/14 17:33:52 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

#ifndef RAR_KMATRIXLASS
#define RAR_KMATRIXLASS

#include "TList.h"
#include "TString.h"
#include "TObject.h"

#include "rarBasePdf.hh"

/// \brief LASS PDF builder
///
/// \par Config Directives:
/// <a href="http://www.slac.stanford.edu/~zhanglei/RooRarFit/RooRarFit.html#sec_LassPdf">See doc for LassPdf configs.</a>
class rarLass : public rarBasePdf {
  
public:
  rarLass();
  rarLass(const char *configFile,
	    const char *configSec, const char *configStr,
	    rarDatasets *theDatasets, RooDataSet *theData,
	    const char *name, const char *title);
  virtual ~rarLass();
  
protected:
  void init();
  
  RooAbsReal *_x; ///< Default obs
  RooAbsReal *_mean;       ///< LASS Resonance K*0(1430) mass
  RooAbsReal *_width;      ///< LASS Resonance mass width
  RooAbsReal *_effRange;   ///< LASS S-wave Effective Range
  RooAbsReal *_scatlen;    ///< LASS S-wave Scattering Length
  RooAbsReal *_turnOffVal; ///< LASS S-wave Turn off value
  
private:
  rarLass(const rarLass&);
  ClassDef(rarLass, 0) // RooRarFit User-defined Pdf class
    ;
};

#endif
