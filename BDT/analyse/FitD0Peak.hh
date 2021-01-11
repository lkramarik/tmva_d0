#ifndef FitD0Peak_hh
#define FitD0Peak_hh

#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TList.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TCut.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include <iostream>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include "TObject.h"

using namespace std;
using namespace TMath;

class FitD0Peak : public TObject {
public:
    FitD0Peak(TH1F*, TH1F*, Float_t, Float_t, TString);
    ~FitD0Peak();

    float getMeanError();
    float getMean();
    float getSigmaError();
    float getSigma();
    float getHeight();
    float getSignificance();
    float getRawYield();
    float getRawYieldError();
    float getRawYieldFit();
    float getRawYieldFitError();

    void makeTuple(TString, TCut, bool);
    void setMean(float);
    void setSigma(float);
    void setHeight(float);
    void addMixedEventBckg(TH1F*);
    void setFitRange(float, float);
    void plotPub();
    void plotPubWithResidual();
    void plotOrig();
    void doStuff();
    void makeTupleCalculations();

    void setTupleSignal(TNtuple*);
    void setTupleBackground(TNtuple*);

    TH1F* sigOrig;
    TH1F* bckgOrig;
    TH1F* mxdOrig;
    TH1F* bckgToWork;

private:
    TF1* fun0;
    TLatex tx2;

    Float_t mMean;
    Float_t mSigma;
    Float_t mSigmaE;
    Float_t mMeanE;
    Float_t mHeight;
    Float_t mHeightE;

    Float_t nsigma;
    Float_t binSize;

    Float_t mFitRMin;
    Float_t mFitRMax;
    Float_t ptMin;
    Float_t ptMax;

    TNtuple* mTupleSignal;
    TNtuple* mTupleBackground;

    bool scale;

    TH1F* sigSubtracted;
    TH1F* sigSubtractedResidualBckg;
    TH1F* bckgAddedResidualBckg;
    TFile* fOut;

    Float_t significanceBins;
    Float_t SoverB;
    Double_t rawYieldError;
    Float_t rawYield;

    Double_t  rawYieldFit;
    Double_t  rawYieldFitError;

    bool lines;
    TLine* rightLine;
    TLine* leftLine;

    void setHistoStyle(TH1F*, Int_t, Int_t);
    void fitFunction();
    void fitComeOn();

    const float peakMin;
    const float peakMax;

    bool isMxdEv;

    TPaveText *text1;



ClassDef(FitD0Peak,1)

};

inline float FitD0Peak::getMean() {return mMean;}
inline float FitD0Peak::getMeanError() {return mMeanE;}
inline float FitD0Peak::getHeight() {return mHeight;}
inline float FitD0Peak::getSigma() {return mSigma;}
inline float FitD0Peak::getSigmaError() {return mSigmaE;}
inline float FitD0Peak::getSignificance() {return significanceBins;}
inline float FitD0Peak::getRawYield() {return rawYield;}
inline float FitD0Peak::getRawYieldError() {return rawYieldError;}
inline float FitD0Peak::getRawYieldFit() {return rawYieldFit;}
inline float FitD0Peak::getRawYieldFitError() {return rawYieldFitError;}

inline void FitD0Peak::setMean(float a) {mMean = a;}
inline void FitD0Peak::setSigma(float a) {mSigma = a;}
inline void FitD0Peak::setHeight(float a) {mHeight = a;}
inline void FitD0Peak::setFitRange(float min, float max) {mFitRMin=min; mFitRMax=max;}
inline void FitD0Peak::setTupleSignal(TNtuple* nt) {mTupleSignal=nt;}
inline void FitD0Peak::setTupleBackground(TNtuple* nt) {mTupleBackground=nt;}


#endif
