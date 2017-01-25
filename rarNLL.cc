/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarNLL.cc,v 1.18 2014/09/14 17:33:53 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides NLL class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides NLL class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"
#include <sstream>
#include <vector>
using namespace std;

#include "TMatrixD.h"
#include "TArrayI.h"
#include "TAxis.h"
#include "TMath.h"

#include "rarNLL.hh"

ClassImp(rarNLL)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarNLL::rarNLL()
  : TNamed("",""),
    _xs(0), _ys(0), _mIdx(-1), _x0s(0), _x1s(0), _x2s(0),
    _iLIntegrals(0), _tLIntegrals(0),
    _verbose(kFALSE)
{
  init();
}

/// \brief Default ctor
///
/// \param curve NLL curve
/// \param name The name
/// \param title The title
/// \param verbose Boolean to control debug output
///
/// Default ctor to set several common data members,
/// and then call #init.
rarNLL::rarNLL(RooCurve *curve, const char *name, const char *title,
               const Bool_t verbose)
  : TNamed(name, title),
    _xs(0), _ys(0), _mIdx(-1), _x0s(0), _x1s(0), _x2s(0),
    _iLIntegrals(0), _tLIntegrals(0),
    _verbose(verbose)
{
  init(curve);
}

/// \brief Trivial dtor
rarNLL::~rarNLL()
{
}

/// \brief Initial function called by ctor
/// \param curve NLL curve to use
///
/// It will initialize NLL data members based on NLL curve
void rarNLL::init(RooCurve *curve)
{
  _nPoints=0;
  _nSteps=0;
  _coeffMList.Delete();
  _lCoefMList.Delete();
  
  if (!curve) {
    cout<<" Can not find any NLL curve"<<endl;
    return;
  }
  // now clone the curve
  curve=(RooCurve*)curve->Clone();
  // sort it
  curve->Sort();
  {
    // make sure two adjacent points are not too close
    Double_t dLimit=curve->GetXaxis()->GetXmax()-curve->GetXaxis()->GetXmin();
    dLimit=TMath::Abs(dLimit/2000.);
    if (0==dLimit) dLimit=1./2000.;
    // check points
    for (Int_t i=0; i<curve->GetN()-1; i++) {
      Double_t xlo, ylo, xhi, yhi;
      curve->GetPoint(i, xlo, ylo);
      curve->GetPoint(i+1, xhi, yhi);
      if (TMath::Abs(xhi-xlo)<dLimit) {
	Int_t rIdx=ylo<yhi?i+1:i;
	if (0==i) rIdx=i+1; // do not remove end point
	if (curve->GetN()-2==i) rIdx=i; // do not remove end point
	curve->RemovePoint(rIdx);
	if (_verbose)
	  cout<<" Point #"<<i<<" ("<<xlo<<", "<<ylo<<") and "<<endl
	      <<" point #"<<i+1<<" ("<<xhi<<", "<<yhi<<")"
	      <<" are too close,"<<endl
	      <<" remove point #"<<rIdx<<" for NLL curve calculation"<<endl;
	i--; // check that spot again
      }
    }
  }
  // first make sure the curve has at least 3 points
  _nPoints=curve->GetN();
  if (_nPoints<3) {
    cout<<" At least 3 points in the NLL curve ("<<curve->GetName()
	<<")are needed"<<endl;
    // do not need the curve any longer
    delete curve;
    return;
  }
  _nSteps=(_nPoints-1)/2;
  // set right size for arrays
  _xs.Set(_nPoints);
  _ys.Set(_nPoints);
  _x0s.Set(_nSteps);
  _x1s.Set(_nSteps);
  _x2s.Set(_nSteps);
  _iLIntegrals.Set(_nSteps);
  _tLIntegrals.Set(_nSteps+1); // The last one is the total integral
  
  // save points, and minY index
  _mIdx=0;
  for (Int_t i=0; i<_nPoints; i++) {
    curve->GetPoint(i, _xs[i], _ys[i]);
    if (_ys[i]<_ys[_mIdx]) _mIdx=i;
  }
  
  // calculate coeff matrics, integrals, etc.
  for (Int_t i=0; i<_nSteps; i++) {
    TMatrixD X(3,3), Y(3, 1, &_ys[2*i]), lY(3,1);
    // lY is y's in Likelihood space
    lY(0,0)=exp(-.5*Y(0,0));
    lY(1,0)=exp(-.5*Y(1,0));
    lY(2,0)=exp(-.5*Y(2,0));
    
    _x0s[i]=_xs[2*i+0];
    _x1s[i]=_xs[2*i+1];
    _x2s[i]=_xs[2*i+2];
    
    X(0,0)=_x0s[i]*_x0s[i];
    X(0,1)=_x0s[i];
    X(0,2)=1;
    X(1,0)=_x1s[i]*_x1s[i];
    X(1,1)=_x1s[i];
    X(1,2)=1;
    X(2,0)=_x2s[i]*_x2s[i];
    X(2,1)=_x2s[i];
    X(2,2)=1;
    
    // invert X
    TMatrixD invX=TMatrixD(TMatrixD::kInverted, X);
    // get coeff matrix
    TMatrixD *A=new TMatrixD(invX*Y);
    _coeffMList.Add(A);
    // coeff matrix in likelihood space
    TMatrixD *lA=new TMatrixD(invX*lY);
    _lCoefMList.Add(lA);
    
    // calculate integral
    Double_t a=A->operator()(0,0);
    Double_t b=A->operator()(1,0);
    Double_t c=A->operator()(2,0);
    // check if a <= 0
    if (a<=0) {
      cout<<" W A R N I N G !"<<endl
	  <<" The 3-point fit on (X, Y)"<<endl
	  <<" ("<<X(0,1)<<", "<<Y(0,0)<<")"<<endl
	  <<" ("<<X(1,1)<<", "<<Y(1,0)<<")"<<endl
	  <<" ("<<X(2,1)<<", "<<Y(2,0)<<")"<<endl
	  <<" returns a="<<a<<" (<=0)"
	  <<" b="<<b<<" c="<<c<<endl
	  <<" Please make sure the scanning fits around there are normal"
	  <<endl;
      if (i<=0) { // do not allow
        curve->RemovePoint(0);
        init(curve);
        delete curve;
        return;
      }
      if (i>=_nSteps-1) { // do not allow
        curve->RemovePoint(_nPoints-1);
        init (curve);
        delete curve;
        return;
      }
    }
    
    Double_t x=_x0s[i];
    const Double_t bignumber = 1e10;
    if (i<=0) x = -bignumber;
    // get lower integral for this step
    _iLIntegrals[i]=lIntegralFunc(x, *A, *lA);
    x=_x2s[i];
    if (i>=_nSteps-1) x = bignumber;
    // get integral for this step
    Double_t thisIntegral=lIntegralFunc(x, *A, *lA)-_iLIntegrals[i];
    if (thisIntegral<0) {
      cout<<" Negative integral="<<thisIntegral
          <<" for i="<<i<<" x0="<<_x0s[i]
          <<" x1="<<_x1s[i]<<" x2="<<_x2s[i]<<endl;
      thisIntegral=0;
    }
    // get total integral for next step
    if (i<=0) _tLIntegrals[0]=0;
    _tLIntegrals[i+1]=_tLIntegrals[i]+thisIntegral;
  }
  // printout
  if (_verbose) {
    cout<<" rarNLL based on:"<<endl;
    curve->Print("v");
    _coeffMList.Print();
    cout<<" "<<_nSteps<<" steps for integral calculation:"<<endl;
    for (Int_t i=0; i<_nSteps; i++) {
      cout<<" i="<<i<<"\t x0="<<_x0s[i]
          <<"\t I="<<_iLIntegrals[i]
          <<"\t T="<<_tLIntegrals[i]
          <<endl;
    }
  }
  
  // do not need the curve any longer
  delete curve;
  return;
}

/// \brief Return NLL for a given value
/// \param x The point for NLL
/// \return The NLL value for a given point
///
/// It returns an NLL value for a given point
Double_t rarNLL::getNLL(Double_t x)
{
  return getY(x);
}

/// \brief Return NLL (Y) for a given value
/// \param x The point for NLL
/// \return The NLL value for a given point
///
/// It returns an NLL value for a given point
Double_t rarNLL::getY(Double_t x)
{
  Double_t retVal(0);
  if (_nSteps<=0) return retVal;
  // if exact match, return cached value
  for (Int_t i=0; i<_nPoints; i++) {
    if (x==_xs[i]) {
      if (_verbose) cout<<" curve point #"<<i<<" = "<<_ys[i]<<endl;
      return _ys[i];
    }
  }
  
  // first find the bin
  Int_t i;
  for (i=0; i<_nSteps; i++) {
    if (x<=_x2s[i]) break;
  }
  if (i<0) i=0;
  if (i>=_nSteps) i=_nSteps-1;
  TMatrixD *A=(TMatrixD*)_coeffMList.At(i);
  Double_t a=A->operator()(0,0);
  Double_t b=A->operator()(1,0);
  Double_t c=A->operator()(2,0);
  retVal=a*x*x+b*x+c;
  
  if (_verbose) {
    cout<<" i="<<i
        <<" a="<<a
        <<" b="<<b
        <<" c="<<c
        <<endl;
    cout<<" NLL("<<x<<")="<<retVal<<endl;
  }
  return retVal;
}

/// \brief Return x values for a given y
/// \param y NLL value
/// \return Array having x's
///
/// It returns array of x's for a given y value
TArrayD rarNLL::getX(Double_t y)
{
  TArrayD X(0); // return array
  if (_nSteps<=0) return X;
  const Double_t bignumber = 1e10;
  // go through bins to find x's
  for (Int_t i=0; i<_nSteps; i++) {
    TArrayD s(0); // solutions for the bin
    Double_t x0=_x0s[i];
    if (0==i) x0 = -bignumber;
    Double_t x2=_x2s[i];
    if (_nSteps-1==i) x2 = bignumber;
    TMatrixD *A=(TMatrixD*)_coeffMList.At(i);
    Double_t a=A->operator()(0,0);
    Double_t b=A->operator()(1,0);
    Double_t c=A->operator()(2,0);
    if ((0==a)&&(0==b)) continue;
    Double_t b24ac=b*b-4*a*(c-y);
    if (0==a) { // for a==0
      s.Set(1);
      s[0]=(y-c)/b;
    } else if (b24ac>=0) {
      s.Set(2);
      s[0]=(-b-sqrt(b24ac))/2/a;
      s[1]=(-b+sqrt(b24ac))/2/a;
    } else {
      continue;
    }
    // save solutions within the bin
    for (Int_t j=0; j<s.GetSize(); j++) {
      if ((x0<=s[j])&&(s[j]<x2)) {
	X.Set(X.GetSize()+1);
	X[X.GetSize()-1]=s[j];
      }
    }
  }
  // sort x's
  TArrayI I(X.GetSize());
  TMath::Sort(X.GetSize(), X.GetArray(), I.GetArray(), kFALSE);
  // combine points too close
  const Double_t limit = 1e-6;
  Double_t dLimit=(_x2s[_nSteps-1]-_x0s[0]) * limit;
  if (0==dLimit) dLimit = limit;
  if (dLimit>limit) dLimit = limit;
  TArrayD XX(0);
  if (I.GetSize()>0) {
    XX.Set(1);
    XX[0]=X[I[0]];
  }
  for (Int_t i=1; i<I.GetSize(); i++) {
    Double_t x=X[I[i]];
    Int_t iXX=XX.GetSize()-1;
    if (TMath::Abs(XX[iXX]-x)<dLimit) { // combine them
      if (_verbose) {
	cout<<" Merging "<<X[iXX]<<" and "<<x
	    <<" w/ diff: "<<x-X[iXX]<<endl;
      }
      XX[iXX]=(X[iXX]+x)/2; /// \todo should use weights
    } else { // add x to XX
      XX.Set(iXX+2);
      XX[iXX+1]=x;
    }
  }
  
  if (_verbose) {
    cout<<" WRT y="<<y<<" x's are:"<<endl;
    for (Int_t i=0; i<XX.GetSize(); i++) {
      cout<<"  #"<<i<<"="<<XX[i];
    }
    cout<<endl;
  }
  
  return XX;
}

/// \brief Replace array contents if y calculated is smaller
/// \param xy Array for min (x,y)
/// \param x x to be evaluated
/// \param a 2nd order coeff from fit
/// \param b 1st order coeff from fit
/// \param c 0th order coeff from fit
///
/// It replaces the array contents if y calculated is smaller
void rarNLL::getMin(TArrayD &xy, Double_t x, Double_t a,Double_t b,Double_t c)
{
  const Double_t bignumber = 1e10;
  if (xy.GetSize()<2) {
    xy.Set(2);
    xy.Reset(bignumber);
  }
  Double_t y=a*x*x+b*x+c;
  if (y<xy[1]) {
    xy[0]=x;
    xy[1]=y;
  }
  return;
}

/// \brief Return the local min point for a given bin
/// \param x0 low x boundary
/// \param x1 middle point
/// \param x2 high x boundary
/// \param a 2nd order coeff from fit
/// \param b 1st order coeff from fit
/// \param c 0th order coeff from fit
/// \return Array of min point (x, y)
///
/// It calculates the min point on NLL curve for a given bin
TArrayD rarNLL::getMin(Double_t x0, Double_t x1, Double_t x2,
		       Double_t a, Double_t b, Double_t c)
{
  const Double_t bignumber = 1e10;
  TArrayD xy(2);
  xy.Reset(bignumber);
  // get min from boundaries
  getMin(xy, x0, a, b, c);
  getMin(xy, x2, a, b, c);
  // for normal fit
  if (a>0) {
    const Double_t limit = 1e-5;
    Double_t xmin=-b/2/a;
    // check if xmin is actually x1
    Double_t dLimit=TMath::Abs(x2-x0) * limit;
    if (dLimit>limit) dLimit = limit;
    if (TMath::Abs(xmin-x1)<dLimit) {
      if (_verbose) cout<<" Set xmin from "<<xmin<<" to "<<x1<<endl;
      xmin=x1;
    }
    if ((x0<xmin)&&(xmin<x2)) getMin(xy, xmin, a, b, c);
    else if (_verbose) {
      Double_t ymin=-b*b/4/a+c;
      cout<<"  parabolic min ("<<xmin<<", "<<ymin<<")"
	  <<" not within x ("<<x0<<", "<<x2<<")"<<endl;
    }
  }
  
  if (_verbose) cout<<"  local min: ("<<xy[0]<<", "<<xy[1]<<")"<<endl;
  return xy;
}

/// \brief Return the min point
/// \param x x var
/// \param y y var
///
/// It calculates the min point on NLL curve
void rarNLL::getMin(Double_t &x, Double_t &y)
{
  if (_nSteps<=0) return;
  const Double_t bignumber = 1e10;
  x=y=bignumber;
  // go through bins to find min bin and min
  for (Int_t i=0; i<_nSteps; i++) {
    Double_t x0=_x0s[i];
    if (0==i) x0 = -bignumber;
    Double_t x1=_x1s[i];
    Double_t x2=_x2s[i];
    if (_nSteps-1==i) x2 = bignumber;
    TMatrixD *A=(TMatrixD*)_coeffMList.At(i);
    Double_t a=A->operator()(0,0);
    Double_t b=A->operator()(1,0);
    Double_t c=A->operator()(2,0);
    if (_verbose)
      cout<<" i="<<i
          <<" x0="<<x0
          <<" x1="<<x1
          <<" x2="<<x2
          <<" a="<<a
          <<" b="<<b
          <<" c="<<c
          <<endl;
    TArrayD xy(getMin(x0, x1, x2, a, b, c));
    if (2==xy.GetSize()) { // local min
      if (xy[1]<y) {
	x=xy[0];
	y=xy[1];
      }
    }
  }
  
  if (_verbose) cout<<" global min: ("<<x<<", "<<y<<")"<<endl;
  return;
}

/// \brief Return the min point (x,y)
/// \return The array of min point (x,y)
///
/// It return the min point in array
TArrayD rarNLL::getMin()
{
  TArrayD xy(2);
  getMin(xy[0], xy[1]);
  return xy;
}

/// \brief Return the likelihood curve integral
/// \param x The limit for integration
/// \return The likelihood curve integral
///
/// It returns likelihood curve integral from -inf to the point
Double_t rarNLL::getLIntegral(Double_t x)
{
  Double_t retVal(0);
  const Double_t bignumber = 1e10;
  if (_nSteps<=0) return retVal;
  
  // first find the bin
  Int_t i(0);
  for (i=0; i<_nSteps; i++) {
    Double_t x0=_x0s[i];
    if (0==i) x0 = -bignumber;
    if (x==x0) {
      if (_verbose)
	cout<<" i="<<i<<" x="<<x<<" integral="<<_tLIntegrals[i]<<endl;
      return _tLIntegrals[i];
    }
    if (x<x0) break;
  }
  // the point falls into the previous bin
  i--;
  if (i<0) return 0; // below the lowest point
  // get x0, a, b, c
  Double_t x0=_x0s[i];
  TMatrixD *A=(TMatrixD*)_coeffMList.At(i);
  TMatrixD *lA=(TMatrixD*)_lCoefMList.At(i);
  if (_verbose) cout<<" i="<<i<<" x0="<<x0<<endl;
  // get integral for this x, a b c within the bin
  Double_t thisIntegral=lIntegralFunc(x, *A, *lA)-_iLIntegrals[i];
  retVal=_tLIntegrals[i]+thisIntegral;
  
  if (_verbose)
    cout<<" integral="<<retVal<<"\t(base+"<<thisIntegral<<")"<<endl;
  return retVal;
}

/// \brief Return the likelihood curve integral
/// \param xl The low limit for integration
/// \param xh The high limit for integration
/// \return The likelihood curve integral
///
/// It returns likelihood curve integral between the limits
Double_t rarNLL::getLIntegral(Double_t xl, Double_t xh)
{
  return getLIntegral(xh)-getLIntegral(xl);
}

/// \brief Return the likelihood curve integral
/// \return The likelihood curve integral
///
/// It returns likelihood curve integral
Double_t rarNLL::getLIntegral()
{
  const Double_t bignumber = 1e10;
  return getLIntegral(0., bignumber);
}

/// \brief Calculate x for given integral
/// \param xl The intergral starting point
/// \param iVal The integral
/// \return The x for the given integral
///
/// It calculates x for a given integral value
Double_t rarNLL::getLIntegralInverse(Double_t xl, Double_t iVal)
{
  const Double_t bignumber = 1e10;
  if (_nSteps<=0) return 0;
  // increase iVal by the integral at xl
  iVal+=getLIntegral(xl);
  
  Double_t x=_x0s[0];
  // Make sure iVal is valid
  if (iVal<=0) return (-bignumber);
  if (iVal>=_tLIntegrals[_nSteps]) return (bignumber);
  
  // first identify which step iVal belongs to
  Int_t i(0);
  for (i=0; i<_nSteps; i++) {
    if (iVal<_tLIntegrals[i+1]) break;
  }
  
  // get x0, a, b, c
  Double_t x0=_x0s[i];
  TMatrixD *A=(TMatrixD*)_coeffMList.At(i);
  Double_t a=A->operator()(0,0);
  Double_t b=A->operator()(1,0);
  Double_t c=A->operator()(2,0);
  // get iVal's residual within the bin
  Double_t residual=iVal-_tLIntegrals[i];
  if (_verbose) 
    cout<<" i="<<i<<" x0="<<x0
	<<" a="<<a
	<<" b="<<b
	<<" c="<<c
	<<" residual="<<residual
	<<endl;
  // get integral from this bin for this x
  Double_t thisXI=_iLIntegrals[i]+residual;
  if (_verbose) cout<<" thisXI="<<thisXI<<endl;
  // identify integral method as in rarNLL::lIntegralFunc().
  // for a=b=0
  if ((0==a)&&(0==b)) {
    cout<<" Using constant fit"<<endl;
    x=thisXI*exp(c/2.);
  } else if (a>0) { // normal fit
    // get erf value
    Double_t thisXerf=thisXI*sqrt(2*a/TMath::Pi())*exp(-(b*b-4*a*c)/8/a);
    if ((thisXerf<-1)||(1<thisXerf)) {
      cout<<" "<<thisXerf<<" should be between -1 and 1"<<endl;
    } else {
      Double_t y=TMath::ErfInverse(thisXerf);
      x=(y*sqrt(8*a)-b)/2/a;
    }
  } else if ((0!=b)&&(a>-1e-5)&&(TMath::Abs(a/b)<100)) { // a ~ 0
    cout<<" Using linear fit"<<endl;
    Double_t thisXExp=-thisXI*b/2.;
    if (thisXExp<=0) {
      cout<<" "<<thisXExp<<" should be >0"<<endl;
    } else {
      x=-(2*TMath::Log(thisXExp)+c)/b;
    }
  } else { // numerical integral
    TMatrixD *lA=(TMatrixD*)_lCoefMList.At(i);
    Double_t la=lA->operator()(0,0);
    Double_t lb=lA->operator()(1,0);
    Double_t lc=lA->operator()(2,0);
    x=lIntegralFuncInverse(_x0s[i], _x2s[i], 10, thisXI, la, lb, lc);
    cout<<" Using numerical integral: la="<<la<<" lb="<<lb<<" lc="<<lc<<endl
	<<" x="<<x<<endl;
  }
  
  if (_verbose) cout<<" x="<<x<<endl;
  return x;
}

/// \brief Calculate x for given integral
/// \param iVal The integral
/// \return The x for the given integral
///
/// It calculates x for a given integral value
Double_t rarNLL::getLIntegralInverse(Double_t iVal)
{
  return getLIntegralInverse(0., iVal);
}

/// \brief The likelihood curve integral function
/// \param x x value
/// \param A Coeff matrix
/// \param lA Coeff matrix in likelihood space
///
/// It is the likelihood integral function based on 2NLL curve parameters
/// The formula is based on exp(-y/2) where y is from 2NLL curve
/// \verbatim
/// Integrate[Exp[(-a x x - b x - c)/2],x]
/// 
///   2
///  b /(8 a) - c/2      Pi          b + 2 a x
/// E               Sqrt[--] Erf[-----------------]
///                      2       2 Sqrt[2] Sqrt[a]
/// -----------------------------------------------
///                     Sqrt[a]
/// \endverbatim
/// If a<0, the analytical integral is a function of imaginary erf
/// \verbatim
/// Integrate[Exp[(a x x - b x - c)/2],x]
///
///    2
///  -b /(8 a) - c/2      Pi          -b + 2 a x
/// E                Sqrt[--] Erfi[-----------------]
///                       2        2 Sqrt[2] Sqrt[a]
/// -------------------------------------------------
///                      Sqrt[a]
/// \endverbatim
/// We prefer to use numerical integral in likelihood space
/// for a<0 because of the lack of Erfi and InverseErfi in ROOT/RooFit
/// and the complicated nature of this function.
/// If a==0, we use analytical integral
/// \verbatim
/// Integrate[Exp[(- b x - c)/2], x]
///
///     -c/2 - (b x)/2
/// -2 E
/// ------------------
///         b
/// \endverbatim
/// If a==0 and b==0, we use
/// \verbatim
/// Integrate[Exp[(- c)/2], x]
///
///  x
/// ----
///  c/2
/// E
/// \endverbatim
///
/// http://mathworld.wolfram.com/Erf.html \n
/// http://mathworld.wolfram.com/Erfi.html
Double_t rarNLL::lIntegralFunc(Double_t x, TMatrixD &A, TMatrixD &lA)
{
  Double_t retVal(0);
  Double_t a=A(0,0);
  Double_t b=A(1,0);
  Double_t c=A(2,0);
  // for a=b=0
  if ((0==a)&&(0==b)) {
    cout<<" Using constant fit"<<endl;
    retVal=x/exp(c/2.);
  } else if (a>0) { // normal fit
    retVal=sqrt(TMath::Pi()/2/a)*exp((b*b-4*a*c)/8/a)*
      TMath::Erf((b+2*a*x)/sqrt(8*a));
  } else if ((0!=b)&&(a>-1e-5)&&(TMath::Abs(a/b)<100)) { // a ~ 0
    cout<<" Using linear fit"<<endl;
    retVal=-2/b*exp((-b*x-c)/2.);
  } else { // numerical integral
    Double_t la=lA(0,0);
    Double_t lb=lA(1,0);
    Double_t lc=lA(2,0);
    retVal=la/3*x*x*x+lb/2*x*x+lc*x;
    cout<<" Using numerical integral: la="<<la<<" lb="<<lb<<" lc="<<lc<<endl
	<<" x="<<x<<" retV:"<<retVal<<endl;
  }
  const Double_t bignumber = 1e200;
  if (retVal>bignumber)  { retVal = bignumber; }
  if (retVal<-bignumber) { retVal = -bignumber; }
  if (_verbose)
    cout<<"x:"<<x<<" a:"<<a<<" b:"<<b<<" c:"<<c<<" retV LIFunc:"<<retVal<<endl;
  return retVal;
}

/// \brief Calculate x for given integral using numerical evaluation
/// \param x0 Lower limit
/// \param x2 Upper limit
/// \param iter Number of iterations
/// \param thisXI The integral
/// \param la Coeff la
/// \param lb Coeff lb
/// \param lc Coeff lc
/// \return The x for the given integral
///
/// It calculates x for a given integral value using numerical evaluation
/// (divide and conquer)
Double_t rarNLL::lIntegralFuncInverse(Double_t x0, Double_t x2, Int_t iter,
                                      Double_t &thisXI,
                                      Double_t &la, Double_t &lb, Double_t &lc)
{
  Double_t x=(x0+x2)/2.;
  if (iter--<=0) return x;
  Double_t theIntegral=la/3*x*x*x+lb/2*x*x+lc*x;
  if (theIntegral==thisXI) return x;
  if (theIntegral<thisXI)
    return lIntegralFuncInverse(x, x2, iter, thisXI, la, lb, lc);
  return lIntegralFuncInverse(x0, x, iter, thisXI, la, lb, lc);
}
