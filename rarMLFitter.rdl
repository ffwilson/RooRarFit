/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: rarMLFitter.rdl,v 1.55 2014/09/14 17:33:53 fwilson Exp $
 * Authors: Lei Zhang
 * History:
 *
 * Copyright (C) 2005-2012, University of California, Riverside
 *****************************************************************************/
#ifndef RAR_MLFITTER
#define RAR_MLFITTER

#include "TList.h"
#include "TString.h"
#include "TObject.h"
#include "TMatrixD.h"
#include "TArrayD.h"

#include "rarCompBase.hh"

class RooFormulaVar;
class RooMCStudy;
class RooSimPdfBuilder;

/// \brief ML fitter builder
///
/// Build the final mlFit model.
/// The config section is specified with
/// command line option
/// \verbatim
/// -C "<fitter config section>"\endverbatim
/// If not specified, the default section is \p "mlFitter Config".
/// The final configStr,
/// \verbatim
/// mlFitter MLFitter "ML Function"\endverbatim
/// where \p mlFitter is its name, \p MLFitter is its pdfType,
/// and \p "ML Function" is its title,
/// is given by main() when it instantiates the ml fitter.
/// \par Config Directives:
/// <a href="http://rarfit.sourceforge.net/RooRarFit.html#sec_rarMLFitterConfig">See doc for this RooRarFit PDF configs.</a>
///
/// After all PDF components of the whole fitter are created,
/// the #main function calls rarMLFitter::run to finish fitting jobs.
/// The actions in this step is controlled by action config section
/// specified with command line option
/// \verbatim
/// -A "<fitter action section>"\endverbatim
/// If not specified, the default section is \p "Action Config".
class rarMLFitter : public rarCompBase {
  
public:
  rarMLFitter();
  rarMLFitter(const char *configFile, const char *configSec,
	      const char *configStr, rarDatasets *theDatasets,
	      RooDataSet *theData, const char *name, const char *title);
  virtual ~rarMLFitter();
  
  /// \brief Set #_paramDir
  /// \param paramDir param dir
  ///
  /// It sets the param file dir
  void setParamDir(TString paramDir) {_paramDir=paramDir;}

  /// \brief Set #_resultDir
  /// \param resultDir result dir
  ///
  /// It sets the result root file dir
  void setResultDir(TString resultDir) {_resultDir=resultDir;}
  
  /// \brief Set #_toyID (random seed)
  /// \param toyID toyID to set
  ///
  /// It sets random seed (through #_toyID)
  /// so that different toy study jobs can have different random sequences.
  void setToyID(Int_t toyID=0) {_toyID=toyID;}
  
  /// \brief Set #_toyNexp
  /// \param toyNexp Number of experiments to set
  ///
  /// It sets number of experiment from command line option
  void setToyNexp(Int_t toyNexp=0) {_toyNexp=toyNexp;}
  
  /// \brief Set #_toyDir
  /// \param toyDir toy sample dir
  ///
  /// It sets the toy sample dir
  void setToyDir(TString toyDir) {_toyDir=toyDir;}

  /// \brief Get #_protDataEVars
  virtual RooArgSet getProtDataEVars() {return _protDataEVars;}
  virtual RooAbsPdf *getProtGen();
  virtual RooAbsPdf *getGenerator();
  virtual RooSimultaneous *compGen(RooAbsPdf *gen, RooArgList subGens,
				   RooCategory &compCat);
  /// \brief Get SimPdfBuilder
  /// \return The simPdfBuilder
  virtual RooSimPdfBuilder* getSimBuilder(){return _simBuilder;}
  TString getPhysCat();
  TString getSplitCats();
  RooArgSet getSplitCatSet();
  RooAbsCategory *getSplitCat(RooArgSet &splitCatSet, TString catName);
  
  void run();

  static void addErrToCurve(RooCurve *curve, Double_t errLo, Double_t errHi,
			    Double_t *maxNLL=0);
  static void addErrToCurve(RooCurve *curve, Double_t err, Double_t *maxNLL=0);
  static RooCurve *shiftNLLCurve(RooCurve *curve, Double_t dx, Double_t dy);
  static RooCurve *combineNLLCurves(TList &curves, Bool_t shiftToZero=kTRUE,
				    Double_t *maxNLL=0);
  static RooCurve *combineNLLCurves(TList &curves, Double_t errs[],
                                    Double_t *maxNLL=0);
  static Double_t getMeanErrs(RooCurve *curve,
			      Double_t *errLo=0, Double_t *errHi=0);
  static Double_t getSignf(RooCurve *curve, Double_t refVal=0.);
  static Double_t getUL(RooCurve *curve, Double_t CL=0.90);
  static RooCurve *combCurves(RooCurve* crv1, RooCurve* crv2,
			      const char* formula, const char* valid="1",
			      const char* crvName=0, Int_t lineWidth=2,
			      Int_t lineStyle=2, Int_t lineColor=kBlue);
  
protected:
  void init();
  virtual Bool_t initParams(TString act, RooArgSet fullParams,
			    RooArgSet fullParamsWOI,
			    TString readParams="yes",
			    TString paramFileID="pdfFit",
			    TString readSecParams="yes",
			    Bool_t useFloatFix=kFALSE,
			    RooArgSet postPdfFloatSet=RooArgSet(),
			    RooArgSet preMLFixSet=RooArgSet(),
			    RooArgSet preMLFloatSet=RooArgSet());
  virtual RooArgSet getSpecialSet(TString setName="specialSet");
  virtual void setSpecialStr(RooArgSet *simConfigSet);
  virtual void getSplitCoeffValues(RooArgList coeffList, TString valType,
				   ostream &o);
  virtual void chkBlind(TString dsName);
  virtual RooFitResult *doMLFit(RooDataSet *mlFitData,TString opt,ostream &o, Int_t ncpus=1);
  virtual RooFitResult *doTheFit(RooAbsPdf *pdf, RooDataSet *mlFitData, TString opt, Int_t ncpus=1);
  virtual void doSignf(RooDataSet *mlFitData,TString signfStr,
                       RooFitResult *fitResult, RooArgSet fullParams,
                       ostream &o);
  virtual void doSysStudy(RooDataSet *mlFitData, TString paramsStr,
			  TString varsStr,RooArgSet fullParams, ostream &o);
  virtual void setVariation(RooArgSet theParams, Double_t myV,
			    Bool_t useErr, Bool_t isPlus);
  virtual void calSysErrors(Int_t iParam, RooArgSet &cstudyVars,
                            RooArgSet &studyVars, TArrayD &eArray);
  virtual void avgSysErrors(ostream &o, TString vParams, Int_t pnLen,
                            TArrayD &pV, TArrayD &pArray, TArrayD &mV,
			    TArrayD &mArray, TArrayD &aArray, TMatrixD corrM);
  virtual void outSysErrors(ostream &o, TString rowName, Int_t pnLen,
			    TMatrixD corrM, TArrayD &aArray);
  virtual TMatrixD getCorrMatrix(Int_t nSysParams, TString vParams,
				 Bool_t diagOnly=kFALSE);
  virtual Double_t doGOFChisq(RooDataSet *mlFitData, ostream &o,
                              TList *plotList=0);
  
  virtual RooDataSet *doToyStudy(RooArgSet fullParams);
  virtual void getCompCatDS(TList*ds,RooDataSet *iData,RooCategory *compCat=0);
  virtual RooAbsArg *findSimed(RooArgSet &simSet,
			       TString argName, TString catName,
			       RooAbsPdf *srcPdf=0);
  virtual Int_t randInt(Double_t iNumber);
  virtual void generate(RooMCStudy *theToy, RooAbsPdf *theGen,
                        const TString genOpt, const Int_t toyNexp);
  virtual void generate(RooMCStudy *theToy,
			const RooArgList& etoyDeps,
			const TString genStr, const TString genOpt,
			const Int_t toyNexp);
  virtual RooDataSet *generate(const RooArgList& dependents,
			       const TString genSrcName,
			       Double_t nEvtGen, const TString genOpt);
  virtual RooAbsPdf *getExtCompPdf(RooAbsPdf *thePdf, RooAbsReal *theCoef);
  virtual void getSnB();
  virtual RooDataSet *getLLRDataset(RooDataSet *theData, RooFormulaVar *LLR);
  virtual void doLLRPlot(RooDataSet *projData, TList &plotList);
  virtual void getLLRPlot(TList &plotList, TString plotName,
			  RooAbsPdf *thePdf, Int_t nEvt,
			  RooArgSet theDeps, RooDataSet *protData,
			  Int_t nBins, RooAbsReal *LLRFunc);
  virtual RooPlot *doProjPlot(RooDataSet *projData, RooRealVar *theVar,
			      TList &plotList);
  virtual RooPlot *getProjPlot(RooRealVar *theVar, Double_t plotMin,
			       Double_t plotMax, Int_t nBins,
			       RooDataSet *sliceData, TString frameName,
			       TString frameTitle, TList &plotList,
			       const RooCmdArg &asymCat=RooCmdArg(),
			       const RooCmdArg &nuBins=RooCmdArg(),
			       const RooCmdArg &xerrorscale=RooCmdArg(),
			       RooPlot *frameM=0, RooPlot *frameP=0);
  virtual void scanVarShiftToNorm(RooArgList scanVars, TArrayD &scanVarDiff);
  virtual RooPlot *doScanPlot(TList &plotList);
  virtual RooPlot *doContourPlot(TList &plotList);
  virtual RooPlot *doSPlot(RooRealVar *theVar, TList &plotList);
  
  virtual TString getRootFileName(TString aType, TString configToken="yes");
  virtual TString getParamFileName(TString aType, TString configToken="yes",
				   TString dirName="");
  virtual void paramFileIO(RooArgSet params, TString paramFile,
			   Bool_t In=kTRUE);
  virtual void paramFileI(RooArgSet params, TString paramFile);
  virtual void paramFileO(RooArgSet params, TString paramFile);
  virtual RooPlot *doCombinePlot(TList &plotList);
  virtual RooPlot *combine(Int_t nModes, 
			   const vector<TString> fileNames,
			   const vector<TString> plotNames,
			   const vector<Double_t> fitBias,
			   const vector<Double_t> addSystErrLo,
			   const vector<Double_t> addSystErrHi,
			   const vector<Double_t> uncorrSystErrLo,
			   const vector<Double_t> uncorrSystErrHi,
			   const vector<Double_t> corrSystErr,
			   const TString xAxisTitle,
			   Bool_t doSignf=kTRUE, Bool_t doUL=kFALSE, Double_t CL=0.90);
  
  virtual void saveAsRootFile(const RooDataSet *ds, 
			      const TString rootfilename, 
			      Bool_t withErrors);
  virtual TTree *createTreeFromDataset(const RooDataSet *ds, Bool_t withErrors);

  static TString _physCatStr; ///< physCat
  static TString _splitCatStr; ///< All other splitting Cat string
  static RooArgSet _splitCatSet; ///< All other splitting Cats
  static RooArgSet _splitCatSet2; ///< All other splitting Cats (derived)
  RooSimPdfBuilder *_simBuilder; ///< SimPdf builder
  RooArgSet *_simConfig; ///< SimPdf config ArgSet
  RooAbsCategoryLValue *_category; ///< Category to build SimPdf directly with
  RooAbsPdf *_theGen; ///< Constructed toy generator
  Int_t _protGenLevel; ///< Generating level for protData
  RooDataSet *_protDataset; ///< Master protDataset
  TList _protDatasetsM; ///< List of dataset for each cat in master protD
  TList _protDatasets; ///< List of dataset for each cat
  RooArgSet _protDataEVars; ///< protDataEVars
  RooArgList _preToyRandGenerators; ///< Extra pdfs as param randomizer for toy
  RooAbsPdf *_theToyParamGen; ///< Pdf to randomize params for toy
  RooArgSet _embdObsRandSet; ///< ArgSet of data src to randomize the obs
  RooArgSet _embdObsGens; ///< Extra pdfs as embed obs randomizer
  RooArgSet _extCompPdfs; ///< Extended component pdfs
  RooAbsPdf *_theSPdf; ///< S Pdf
  RooAbsPdf *_theBPdf; ///< B Pdf
  TList _fCompList; ///< Full comp list
  RooArgList _fCoefList; ///< Full coef list
  TList _sCompList; ///< S comp list
  RooArgList _sCoefList; ///< S coef list
  
  TString _paramDir; ///< Dir for param files
  TString _resultDir; ///< Dir for output root files
  Int_t _toyID; ///< Toy ID used as random seed
  Int_t _toyNexp; ///< Number of experiments from command line
  TString _toyDir; ///< Dir for toy samples
  
private:
  rarMLFitter(const rarMLFitter&);
  ClassDef(rarMLFitter, 0) // Final RooRarFit ML fitter class
    ;
};

#endif
