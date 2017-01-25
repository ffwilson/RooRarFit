/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarStrParser.rdl,v 1.7 2014/09/14 17:33:55 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_STRPARSER
#define RAR_STRPARSER

#include "TList.h"
#include "TString.h"
// #include "TObject.h"
#include "RooStringVar.h"

class TObjString;

/// \brief String parser for RooRarFit
///
/// It is the string parser used widely by other RooRarFit classes
/// to direct their actions through config items in config file.
/// It breaks a string into tokens seperated by spaces.
/// Characters inside quote(") are considered one token.
/// Currently, quote (") can not be inside the parsed tokens.
class rarStrParser : public TObject {
  
public:
  rarStrParser();
  rarStrParser(const char *str);
  rarStrParser(const TString str);
  rarStrParser(const rarStrParser&);
  virtual ~rarStrParser();
  void operator=(const char *str);
  void operator=(const TString str);
  void operator=(const RooStringVar str);
  TString &operator[](const Int_t idx);
  void Remove(const Int_t idx=0);
  Int_t Index(const TString token);
  Bool_t Have(const TString token);
  
  /// \brief Return number of tokens
  /// \return Number of tokens stored
  Int_t nArgs() {return _strs.GetSize();}
  
protected:
  void init();
  TObjString *nextToken();
  
  TString _str; ///< String to parse
  TList _strs; ///< List of parsed tokens
  
private:
  Int_t _idx; ///< Current char index used for string parsing
  ClassDef(rarStrParser, 0) // RooRarFit String Parser class
    ;
};

#endif
