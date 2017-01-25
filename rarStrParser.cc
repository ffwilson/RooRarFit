/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarStrParser.cc,v 1.10 2014/09/14 17:33:55 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides String Parser class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides String Parser class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "TObjString.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooStringVar.h"

#include "rarStrParser.hh"

ClassImp(rarStrParser)
  ;

/// \brief Trivial ctor
rarStrParser::rarStrParser()
  : _str("")
{
  init();
}

/// \brief Ctor using char pointer
///
/// \param str String to parse
rarStrParser::rarStrParser(const char *str)
  : _str(str)
{
  init();
}


/// \brief Ctor using TString object
///
/// \param str String to parse
rarStrParser::rarStrParser(const TString str)
  : _str(str)
{
  init();
}

/// \brief Ctor for assignment using its own object
///
/// \param aParser Source object
rarStrParser::rarStrParser(const rarStrParser& aParser)
  : TObject(aParser)
{
  _str=aParser._str;
  // go through tokens and copy them
  for (Int_t i=0; i<((rarStrParser&)aParser).nArgs(); i++)
    _strs.Add(new TObjString(((rarStrParser&)aParser)[i]));
}

rarStrParser::~rarStrParser()
{
  _strs.Delete();
}

/// \brief Operator = with char pointer
///
/// \param str String to parse
void rarStrParser::operator=(const char *str)
{
  _str=str;
  init();
}

/// \brief Operator = with TString object
///
/// \param str String to parse
void rarStrParser::operator=(const TString str)
{
  rarStrParser::operator=(str.Data());
}

/// \brief Operator = with RooStringVar object
///
/// \param str String to parse
void rarStrParser::operator=(const RooStringVar str)
{
  rarStrParser::operator=(str.getVal());
}

/// \brief Get the indexed token
///
/// \param idx Token index in #_strs
/// \return The indexed token
TString &rarStrParser::operator[](const Int_t idx)
{
  TObjString *retVal=(TObjString*)_strs.At(idx);
  if (retVal) return retVal->String();
  return _str;
}

/// \brief Remove the indexed token
///
/// \param idx The index of token to be removed
void rarStrParser::Remove(const Int_t idx)
{
  TObject *val=_strs.At(idx);
  if (val) {
    _strs.Remove(val);
    delete val;
  }
}

/// \brief Check if have the token
/// \param token Token to check
/// \return true if found, false if not
///
/// It checks if the parser has the token
Bool_t rarStrParser::Have(const TString token) {
  return (Index(token)<0) ? kFALSE : kTRUE;
}

/// \brief Get the index for a token
/// \param token Token to check
/// \return the index of the give token, -1 if not found
///
/// It returns index for a given token, -1 if not found.
Int_t rarStrParser::Index(const TString token) {
  Int_t index(-1);
  for (Int_t i=0; i<nArgs(); i++)
    if (token==operator[](i)) return i;
  return index;
}

/// \brief Initial function to parse a string
///
/// It initializes #_idx and #_strs and
/// calls #nextToken to parse the string
void rarStrParser::init()
{
  // reset index and parsed string objs
  _idx=0;
  _strs.Delete();
  TObjString *objStr(0);
  while (objStr=nextToken()) {
    _strs.Add(objStr);
  }
}

/// \brief Get next token in the string from current position
///
/// \todo Make quote (") regular character with `\"'.
///
/// This is the actual parser.
/// It ignores any blank characters at current position
/// and return the sub-string until next blank character as the token.
/// If current non-blank character is quote ("),
/// it will take all characters after that until the next quote (")
/// as token and return it.
/// It return 0 (null) if current position, #_idx, is out of range.
TObjString *rarStrParser::nextToken()
{
  Int_t nIdx(0);
  if (_idx>=_str.Length()) return 0;
  while (' '==_str[_idx] || '\t'==_str[_idx] || '\n'==_str[_idx]) {
    _idx++;
    if (_idx>=_str.Length()) return 0;
  }
  if ('"'==_str[_idx]) {
    _idx++;
    nIdx=_str.Index("\"", _idx);
    if(nIdx<0) return 0;
  } else {
    nIdx=_str.Index(" ", _idx);
    if(nIdx<0) nIdx=_str.Length();
  }
  //cout<<"in Parser "<<nIdx<<" "<<_idx<<" "<<_str(_idx, nIdx-_idx)<<endl;
  TString substr=_str(_idx, nIdx-_idx);
  TObjString *subStr=new TObjString(substr);
  _idx=nIdx;
  if (_idx<_str.Length())
    if ('"'==_str[_idx]) _idx++;
  
  return subStr;
}
