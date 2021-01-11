#ifndef BDT_EFFICIENCYCALCULATION_H
#define BDT_EFFICIENCYCALCULATION_H
#include "TH1.h"
#include "TF1.h"
#include "StAnaCuts.h"


class efficiencyCalculation {
public:
    efficiencyCalculation();
    void setPidGraphs();
    void setTpcGraphs();

    void setFolderTpc(TString);
    void setFolderTpcPid(TString);
    void setFolderTofPid(TString);
    void setFolderTofMatch(TString);

    bool isTpcReconstructed(TString, float, int, float);
    void setHistoSuffix(TString);
    void setTofMatch(int);
    void setTofPid();
    void setTpcPid();
    bool isTofmatched(TString, float, float);
    bool isTpcPid(TString, float);
    bool isTofPid(TString, float);

private:
    bool mPidSet;
    bool mTpcSet;
    bool mTofMatchSet;
    bool mTpcPidSet;
    bool mTofPidSet;
    TString mFolderTpc;
    TString mFolderTpcPid;
    TString mFolderTofPid;
    TString mFolderTofMatch;
    const int nmultEdge;
//    float const multEdge[8];
//    extern TH1D* hTpcPiPlus[]; //embedding
//    extern TH1D* hTpcPiMinus[]; //embedding
//    extern TH1D* hTpcKPlus[]; //embedding
//    extern TH1D* hTpcKMinus[]; //embedding
//
    TH1D* hTpcPiPlus[vars::nmultEdge]; //embedding
    TH1D* hTpcPiMinus[vars::nmultEdge]; //embedding
    TH1D* hTpcKPlus[vars::nmultEdge]; //embedding
    TH1D* hTpcKMinus[vars::nmultEdge]; //embedding
    TString mSuffix;

    TF1 *mfTofPidPion;
    TF1 *mfTofPidKaon;

    TF1 *mfTpcPidPion;
    TF1 *mfTpcPidKaon;

    TF1* mfTofMatch[4];

//    ClassDef(efficiencyCalculation, 1) //set to 1

};

inline void efficiencyCalculation::setFolderTpc(TString m) {mFolderTpc=m; }
inline void efficiencyCalculation::setHistoSuffix(TString m) {mSuffix=m; }
inline void efficiencyCalculation::setFolderTofMatch(TString m) {mFolderTofMatch=m; }
inline void efficiencyCalculation::setFolderTofPid(TString m) {mFolderTofPid=m; }
inline void efficiencyCalculation::setFolderTpcPid(TString m) {mFolderTpcPid=m; }

#endif //BDT_EFFICIENCYCALCULATION_H
