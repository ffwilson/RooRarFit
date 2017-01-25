/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarMinuit.cc,v 1.12 2014/09/14 17:33:53 fwilson Exp $
 * Author:                                                                   *
 *   WTF, Bill Ford, U. of Colorado, wtford@pizero.colorado.edu              *
 *                                                                           *
 *****************************************************************************/

// -- CLASS DESCRIPTION [PDF] --
// This class is derived from RooMinuit and overloads the contour() method to
// produce a RooPlot.
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class is derived from RooMinuit and overloads the contour() method to
// produce a RooPlot.
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"
#include "TH1.h"
#include "TH2.h"
#include "TMarker.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include "TFitter.h"
#include "TMinuit.h"
#include "TDirectory.h"
#include "RooMinuit.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "rarMinuit.hh"

rarMinuit::rarMinuit(RooAbsReal& function) : RooMinuit(function)
{
  // The following are private in the base class so we have to set them up here
  _func = &function ;
  // Examine parameter list
  RooArgSet* paramSet = function.getParameters(RooArgSet()) ;
  RooArgList paramList(*paramSet) ;
  delete paramSet ;

  _floatParamList = (RooArgList*) paramList.selectByAttrib("Constant",kFALSE) ;
  if (_floatParamList->getSize()>1) {
    _floatParamList->sort() ;
  }
  _floatParamList->setName("floatParamList") ;
  
}

RooPlot* rarMinuit::contour(RooRealVar& var1, RooRealVar& var2,
			    Double_t n1, Double_t n2, Double_t n3,
			    Double_t n4, Double_t n5, Double_t n6) 
{
  // Verify that both variables are floating parameters of PDF
  Int_t index1= _floatParamList->index(&var1);
  if(index1 < 0) {
    std::cout << "rarMinuit::contour(" << GetName() 
	      << ") ERROR: " << var1.GetName()
	      << " is not a floating parameter of " 
	      << _func->GetName() << std::endl ;
    return 0;
  }
  
  Int_t index2= _floatParamList->index(&var2);
  if(index2 < 0) {
    std::cout << "rarMinuit::contour(" << GetName() 
	      << ") ERROR: " << var2.GetName()
	      << " is not a floating parameter of PDF "
	      << _func->GetName() << std::endl ;
    return 0;
  }
  
  // create and draw a frame
  RooPlot *frame = new RooPlot(var1, var2) ;
  frame->SetStats(kFALSE);
  
  // draw a point at the current parameter values
  TMarker *point= new TMarker(var1.getVal(), var2.getVal(), 8);
  
  // remember our original value of ERRDEF
  Double_t errdef= gMinuit->fUp;
  
  TGraph* graph1 = 0;
  if(n1 > 0) {
    // set the value corresponding to an n1-sigma contour
    gMinuit->SetErrorDef(n1*n1*errdef);
    // calculate and draw the contour
    graph1= (TGraph*)gMinuit->Contour(25, index1, index2);
    fixGraph(graph1);
  }
  
  TGraph* graph2 = 0;
  if(n2 > 0) {
    // set the value corresponding to an n2-sigma contour
    gMinuit->SetErrorDef(n2*n2*errdef);
    // calculate and draw the contour
    graph2= (TGraph*)gMinuit->Contour(25, index1, index2);
    fixGraph(graph2,2);
  }
  
  TGraph* graph3 = 0;
  if(n3 > 0) {
    // set the value corresponding to an n3-sigma contour
    gMinuit->SetErrorDef(n3*n3*errdef);
    // calculate and draw the contour
    graph3= (TGraph*)gMinuit->Contour(25, index1, index2);
    fixGraph(graph3,3);
  }
  
  TGraph* graph4 = 0;
  if(n4 > 0) {
    // set the value corresponding to an n4-sigma contour
    gMinuit->SetErrorDef(n4*n4*errdef);
    // calculate and draw the contour
    graph4= (TGraph*)gMinuit->Contour(25, index1, index2);
    fixGraph(graph4,4);
  }
  
  TGraph* graph5 = 0;
  if(n5 > 0) {
    // set the value corresponding to an n5-sigma contour
    gMinuit->SetErrorDef(n5*n5*errdef);
    // calculate and draw the contour
    graph5= (TGraph*)gMinuit->Contour(25, index1, index2);
    fixGraph(graph5,5);
  }
  
  TGraph* graph6 = 0;
  if(n6 > 0) {
    // set the value corresponding to an n6-sigma contour
    gMinuit->SetErrorDef(n6*n6*errdef);
    // calculate and draw the contour
    graph6= (TGraph*)gMinuit->Contour(25, index1, index2);
    fixGraph(graph6,6);
  }
  // restore the original ERRDEF
  gMinuit->SetErrorDef(errdef);
  
  
  // Add all objects to the frame
  frame->addObject(point);
  if (graph1) frame->addObject(graph1, "C");
  if (graph2) frame->addObject(graph2, "C");
  if (graph3) frame->addObject(graph3, "C");
  if (graph4) frame->addObject(graph4, "C");
  if (graph5) frame->addObject(graph5, "C");
  if (graph6) frame->addObject(graph6, "C");
  
  return frame;
}

/// \brief make it close graph (fix the graph, looks a bug in Root)
void rarMinuit::fixGraph(TGraph *graph, Int_t lineStyle) {
  if (!graph) return;
  // set line style
  graph->SetLineStyle(lineStyle);
  // add a point to form closure
  Int_t nPoint=graph->GetN();
  graph->Set(nPoint+1);
  Double_t x, y;
  graph->GetPoint(0, x, y);
  graph->SetPoint(nPoint, x, y);
  return;
}

