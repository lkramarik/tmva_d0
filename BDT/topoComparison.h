#ifndef BDT_TOPOCOMPARISON_H
#define BDT_TOPOCOMPARISON_H
#include "TCut.h"
#include "TNtuple.h"
#include "TH1F.h"
#include <vector>

class topoComparison {
public:
    topoComparison();
    topoComparison(TString);
    void addCut(TCut);
    void addNtpToCompare(TNtuple*);
    void addNtpToCompare(TNtuple*, TString);
    void addNtpToCompare(TNtuple*, TString, TCut);
    void addNtpToCompareBackground(TNtuple*, TString, TNtuple*);
    void addNtpToCompareBackground(TNtuple*, TString, TCut, TNtuple*);
    void addNtpToCompare(TNtuple*, TString, TCut, TString);
    void addVarToCompare(TString, Double_t, Double_t);
    void addVarToCompare(TString, Double_t, Double_t, TString);
    void addVarToCompare(TString, Double_t, Double_t, TString, Int_t);
    void addVarToCompare(TString, Double_t, Double_t, TString, Int_t, int);
    void project();
    void setHisto(TH1D*, int);
    void setOutFileName(TString);
    void setFolderName(TString);
    void setText(TString);
    ~topoComparison();

private:
    std::vector<TCut> mCuts;
    std::vector<TCut> mCutsForTuples;
    std::vector<TNtuple*> mTuples;
    std::vector<TNtuple*> mTuplesToSubtract;
    std::vector<TString> mVars;
    std::vector<Double_t> mLimsMin;
    std::vector<Double_t> mLimsMax;
    std::vector<Int_t> mnBins;
    std::vector<TString> mTitleX;
    TString mFolder;
    TString mOutputFileName;
    float findMaximumHistos(TH1F*, int);
    TString mText;
    bool mDrawText;
    std::vector<TString> mWeightExpr;
    std::vector<int> mLogy;

    TCut connectCuts();
};

inline void topoComparison::addCut(TCut setCut) {mCuts.push_back(setCut);}
inline void topoComparison::addNtpToCompare(TNtuple * tuple) {
    mTuples.push_back(tuple);
    mCutsForTuples.push_back("");
    mWeightExpr.push_back("1");
    TNtuple* nullTP=NULL;
    mTuplesToSubtract.push_back(nullTP);
}
inline void topoComparison::addNtpToCompare(TNtuple * tuple, TString name) {
    tuple->SetName(name);
    mTuples.push_back(tuple);
    mCutsForTuples.push_back("");
    mWeightExpr.push_back("1");
    TNtuple* nullTP=NULL;
    mTuplesToSubtract.push_back(nullTP);
}
inline void topoComparison::addNtpToCompare(TNtuple * tuple, TString name, TCut cut) {
    tuple->SetName(name);
    mTuples.push_back(tuple);
    mCutsForTuples.push_back(cut);
    mWeightExpr.push_back("1");
    TNtuple* nullTP=NULL;
    mTuplesToSubtract.push_back(nullTP);
}
inline void topoComparison::addNtpToCompareBackground(TNtuple * tuple, TString name, TNtuple * bckg) {
    tuple->SetName(name);
    mTuples.push_back(tuple);
    mCutsForTuples.push_back("");
    mWeightExpr.push_back("1");
    mTuplesToSubtract.push_back(bckg);
}
inline void topoComparison::addNtpToCompareBackground(TNtuple * tuple, TString name, TCut cut, TNtuple * bckg) {
    tuple->SetName(name);
    mTuples.push_back(tuple);
    mCutsForTuples.push_back(cut);
    mWeightExpr.push_back("1");
    mTuplesToSubtract.push_back(bckg);
}
inline void topoComparison::addNtpToCompare(TNtuple * tuple, TString name, TCut cut, TString weight) {
    tuple->SetName(name);
    mTuples.push_back(tuple);
    mCutsForTuples.push_back(cut);
    mWeightExpr.push_back(weight);
    TNtuple* nullTP=NULL;
    mTuplesToSubtract.push_back(nullTP);
}
inline void topoComparison::addVarToCompare(TString var, Double_t min, Double_t max) {
    mnBins.push_back(50);
    mTitleX.push_back(var);
    mLimsMin.push_back(min);
    mLimsMax.push_back(max);
    mVars.push_back(var);
    mLogy.push_back(0);
}

inline void topoComparison::addVarToCompare(TString var, Double_t min, Double_t max, TString titleX) {
    mnBins.push_back(50);
    mTitleX.push_back(titleX);
    mLimsMin.push_back(min);
    mLimsMax.push_back(max);
    mVars.push_back(var);
    mLogy.push_back(0);
}

inline void topoComparison::addVarToCompare(TString var, Double_t min, Double_t max, TString titleX, Int_t nBins) {
    mnBins.push_back(nBins);
    mTitleX.push_back(titleX);
    mLimsMin.push_back(min);
    mLimsMax.push_back(max);
    mVars.push_back(var);
    mLogy.push_back(0);
}

inline void topoComparison::addVarToCompare(TString var, Double_t min, Double_t max, TString titleX, Int_t nBins, int logy) {
    mnBins.push_back(nBins);
    mTitleX.push_back(titleX);
    mLimsMin.push_back(min);
    mLimsMax.push_back(max);
    mVars.push_back(var);
    mLogy.push_back(logy);
}
inline void topoComparison::setOutFileName(TString name) {mOutputFileName = name;}
inline void topoComparison::setText(TString name) {
    mText = name;
    mDrawText = true;
}

#endif //BDT_TOPOCOMPARISON_H
