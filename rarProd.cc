/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarProd.cc,v 1.14 2014/09/14 17:33:54 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 * 
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/

// -- CLASS DESCRIPTION [RooRarFit] --
// This class provides ProdPdf class for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This class provides ProdPdf class for RooRarFit
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooStringVar.h"
#include "RooGlobalFunc.h"
using namespace RooFit;

#include "rarProd.hh"

ClassImp(rarProd)
  ;

/// \brief Trivial ctor
///
/// Usually the objects should be created using other ctors.
rarProd::rarProd()
  : rarCompBase()
{
  init();
}

/// \brief Default ctor
///
/// \param configFile The config file
/// \param configSec The config section
/// \param configStr The config string
/// \param theDatasets Available datasets
/// \param theData Default dataset for this PDF
/// \param name The name
/// \param title The title
///
/// The default ctor first initializes data members,
/// and then calls #init.
rarProd::rarProd(const char *configFile, const char *configSec,
		 const char *configStr,
		 rarDatasets *theDatasets, RooDataSet *theData,
		 const char *name, const char *title)
  : rarCompBase(configFile, configSec, configStr,
		theDatasets, theData, name, title, kFALSE)
{
  init();
}

rarProd::~rarProd()
{
}

/// \brief Initial function called by ctor
///
/// \p init is called by the ctor.
/// The initialization work to read in config info
/// is done with rarCompBase::init(),
/// and \p init here just builds the RooProdPdf with configured components.
void rarProd::init()
{
  // Config pdfs done by parent Pdf
  // list configed pdfs done by parent pdf
  
  // set control bits
  setControlBit("noUseBasePdfFit", "ndFit");
  
  // list of cmd
  RooLinkedList cmdList;

  // do we have conditional obseervables?
  RooArgSet condObsSet(getArgSet(readConfStr("CondObss",""),kFALSE));
  // do we have conditional Pdf ?
  rarStrParser condPdfsStr=readConfStr("CondPdfs","");
  while (condPdfsStr.nArgs()>0) {
    rarBasePdf *condPdf=(rarBasePdf*)_pdfList.FindObject(condPdfsStr[0]);
    condPdfsStr.Remove();
    if (condPdf) {
      _condPdfList.Add(condPdf);
      _condPdfs.add(*condPdf->getPdf());
      // now set CondObss to this pdf
      //RooArgSet *thisPdfSet=new RooArgSet;
      //thisPdfSet->add(*condPdf->getPdf());
      RooArgSet *thisPdfSet=new RooArgSet(*condPdf->getPdf());
      RooArgSet *thisCondObsSet=new RooArgSet(condObsSet);
      thisCondObsSet->
        add(getArgSet(readConfStr
                      (Form("CondObss_%s",condPdf->GetName()),""),kFALSE));
      condPdf->setCondObss(*thisCondObsSet);
      cmdList.Add(Conditional(*thisPdfSet,*thisCondObsSet).Clone());
      cout<<" "<<condPdf->GetName()<<" is conditional PDF"<<endl;
      cout<<" The normalization of this PDF is on "<<endl;
      thisCondObsSet->Print();
    }
  }
  // Pdfs without conditional
  RooArgList subPdfs(_subPdfs);
  subPdfs.remove(_condPdfs);
  if (_condPdfs.getSize()>0) {
    // do not use ndFit by default
    //setControlBit("UseBasePdfFit", "ndFit");
    // cout...
    cout<<"The following PDFs are conditional PDFs:"<<endl;
    _condPdfs.Print();
    cout<<endl;
    //cout<<"The normalization of those PDFs are on observables:"<<endl;
    //condObsSet.Print();
    //cout<<endl;
    cout<<"The remaining PDFs in the production are:"<<endl;
    subPdfs.Print();
    cout<<endl;
  }
  
  // build ProdPdf
  //cmdList.Print("v");
  _thePdf=new RooProdPdf(Form("the_%s", GetName()), _pdfType+" "+GetTitle(),
			 subPdfs, cmdList);
  
  cout<<"done init of rarProd for "<<GetName()<<endl<<endl;
}

/// \brief Get a PDF not depending on specified observables for sPlot
/// \param ignoredObs The observables to check
/// \return The Pdf required
///
/// It checks of the pdf depends on the observables.
/// if no, it returns the default pdf;
/// if yes, it builds the non-dependent pdf.
RooAbsPdf *rarProd::getPdfWOvar(RooArgList ignoredObs)
{
  RooAbsPdf *thePdf=_thePdf;
  // check if thePdf depends on theVar
  if (!(thePdf->dependsOn(ignoredObs))) return thePdf;
  // create individual Pdf
  TString theName=Form("the_%s_%s", ignoredObs[0].GetName(), GetName());
  TString theTitle=Form("%s w/o %s", GetTitle(), ignoredObs[0].GetName());
  RooArgList myPdfList=getPdfsWOvar(ignoredObs);
  thePdf=new RooProdPdf(theName, theTitle, myPdfList);
  //thePdf->Print("v");
  //thePdf->Print();
  if (!(thePdf->dependsOn(ignoredObs))) return thePdf;
  delete thePdf;
  // create a dummy for it
  thePdf=new RooGenericPdf(theName, theTitle, "1",
			   *_thePdf->getParameters(_fullObs));
  //thePdf->Print("v");
  //thePdf->Print();
  return thePdf;
}

/// \brief Return the pdf if it depends on the var
/// \param theVar The var to check
/// \return This (sim)Pdf if it depends on the var
///
/// It checks if the pdf depends on the obs;
/// if yes, it returns this (sim)Pdf
RooAbsPdf *rarProd::getDPdfWvar(RooRealVar *theVar)
{
  RooAbsPdf *retVal(0);
  for (Int_t i=0; i<_nComp; i++) {
    rarBasePdf *thePdf=(rarBasePdf*)_pdfList.At(i);
    retVal=thePdf->getDPdfWvar(theVar);
    if (retVal) return retVal;
  }
  return retVal;
}

/// \brief Return prototype var generator for toy study
/// \return The generator
///
/// It constrcuts prototype var generator for toy study
RooAbsPdf *rarProd::getProtGen()
{
  if (_theProtGen) return _theProtGen;
  // call super class' getProtGen
  rarCompBase::getProtGen();
  _theProtGen=_myDummyPdf;
  RooArgList pdfList(_protGenPdfs);
  if (!getControlBit("dummySubGenPdfs")) pdfList.add(_subProtGenPdfs);
  else pdfList.add(*_theProtGen);
  if ((pdfList.getSize()>1)||(!getControlBit("dummySubGenPdfs")))
    _theProtGen = new RooProdPdf (Form("protGen_%s",GetName()),
				  Form("Prot Gen for %s", GetName()), pdfList);
  
  return _theProtGen;
}
