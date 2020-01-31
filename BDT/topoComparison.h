//
// Created by lukas on 31.01.20.
//

#ifndef BDT_TOPOCOMPARISON_H
#define BDT_TOPOCOMPARISON_H
#include "TCut.h"
#include "TNtuple.h"
#include "TH1F.h"
#include <vector>

class topoComparison {
public:
    topoComparison();
    void addCut(TCut);
    void addNtpToCompare(TNtuple*);
    void addNtpToCompare(TNtuple*, TString);
    void addVarToCompare(TString, float, float);
    void addVarToCompare(TString, float, float, TString);
    void project();
    void setHisto(TH1F*, int);
    void setOutFileName(TString);
private:
    std::vector<TCut> mCuts;
    std::vector<TNtuple*> mTuples;
    std::vector<TString> mVars;
    std::vector<float> mLimsMin;
    std::vector<float> mLimsMax;
    std::vector<int> mnBins;
    std::vector<TString> mTitleX;
    TString mFolder;
    TString mOutputFileName;
    float findMaximumHistos(TH1F*, int);

    TCut connectCuts();
};

inline void topoComparison::addCut(TCut setCut) {mCuts.push_back(setCut);}
inline void topoComparison::addNtpToCompare(TNtuple * tuple) {mTuples.push_back(tuple);}
inline void topoComparison::addNtpToCompare(TNtuple * tuple, TString name) {
    tuple->SetName(name);
    mTuples.push_back(tuple);
}

inline void topoComparison::addVarToCompare(TString var, float min, float max) {
    mnBins.push_back(100);
    mTitleX.push_back(var);
    mLimsMin.push_back(min);
    mLimsMax.push_back(max);
    mVars.push_back(var);
}

inline void topoComparison::addVarToCompare(TString var, float min, float max, TString titleX) {
    mnBins.push_back(100);
    mTitleX.push_back(titleX);
    mLimsMin.push_back(min);
    mLimsMax.push_back(max);
    mVars.push_back(var);
}
inline void topoComparison::setOutFileName(TString name) {mOutputFileName = name;}

#endif //BDT_TOPOCOMPARISON_H
