#include "efficiencyCalculation.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TDatime.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TPad.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include <iostream>

ClassImp(efficiencyCalculation)

efficiencyCalculation::efficiencyCalculation() :
        mPidSet(false),  mTpcSet(false),
        mFolderTpc("/home/lukas/work/embedding_dAu/analyse/tpc_eff/"),
        mFolderTofMatch("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/tofMatch/"),
        mFolderTofPid("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/tofMatch/"),
        mFolderTpcPid("/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/analyse/tofMatch/"),
        mTofMatchSet(false),
        mTofPidSet(false),
        mTpcPidSet(false),
        nmultEdge(7),
        mSuffix("")
        {
}

//________________________________________________
void efficiencyCalculation::setPidGraphs(){
    mPidSet=true;
    return;
}

//________________________________________________
void efficiencyCalculation::setTpcGraphs() {
    mTpcSet=true;

    TFile fTpcPiPlus(mFolderTpc+"piplus_tpc_eff_embedding.root");
    TFile fTpcPiMinus(mFolderTpc+"piminus_tpc_eff_embedding.root");
    TFile fTpcKPlus(mFolderTpc+"kplus_tpc_eff_embedding.root");
    TFile fTpcKMinus(mFolderTpc+"kminus_tpc_eff_embedding.root");

//    if (!fTpcPiPlus) {
//        cout<<"some TPC files are not ok. End."
//        return;
//    }
    cout<<"Getting TPC from folder "<<mFolderTpc<<endl;
    for (int iCent = 0; iCent < vars::nmultEdge; ++iCent)
    {
        TString hName = Form("TrackEffMult%i", iCent);
        hName+=mSuffix;
        cout<<"Loading "<<hName<<endl;
        hTpcPiPlus[iCent] = (TH1D*)fTpcPiPlus.Get(hName);
        hTpcPiPlus[iCent]->SetDirectory(0);
        hTpcPiMinus[iCent] = (TH1D*)fTpcPiMinus.Get(hName);
        hTpcPiMinus[iCent] ->SetDirectory(0);
        hTpcKPlus[iCent] = (TH1D*)fTpcKPlus.Get(hName);
        hTpcKPlus[iCent]->SetDirectory(0);
        hTpcKMinus[iCent] = (TH1D*)fTpcKMinus.Get(hName);
        hTpcKMinus[iCent]->SetDirectory(0);
    }

    fTpcPiPlus.Close();
    fTpcPiMinus.Close();
    fTpcKPlus.Close();
    fTpcKMinus.Close();

    return;
}

//____________________________________________________________________________________
//void efficiencyCalculation::getTpcEff(){
//    return;
//}


//____________________________________________________________________________________
bool efficiencyCalculation::isTpcReconstructed(TString particle, float charge, int cent, float pt) {
    TH1D* h = NULL;
    if (particle == "pion") {
        if (charge > 0) h = hTpcPiPlus[cent];
        else h = hTpcPiMinus[cent];
    }
    else {
        if (charge > 0) h = hTpcKPlus[cent];
        else h = hTpcKMinus[cent];
    }

    int const bin = h->FindBin(pt);
    return gRandom->Rndm() < h->GetBinContent(bin);
}


//____________________________________________________________________________________
void efficiencyCalculation::setTofMatch(int nsigma=1) {
    TFile* fTofMatch=new TFile(mFolderTofMatch, "READ");
    if (!fTofMatch) return;
    cout<<"TOF matching file loaded."<<endl;

    const int m_nParticlesCharged = 4;
    const TString m_ParticleChargedName[4] = {"PionPlus", "PionMinus", "KaonPlus", "KaonMinus"};

    TString hisName;
    for (int i = 0; i < m_nParticlesCharged; ++i) {
        hisName=Form("f_tofMatch_%s_nsigma%i", m_ParticleChargedName[i].Data(), nsigma);
        mfTofMatch[i]=(TF1*)fTofMatch->Get(hisName);
        if (!mfTofMatch[i]) cout<<hisName<<" not loaded!!!!"<<endl;
    }

    mTofMatchSet= true;
    return;
}

//____________________________________________________________________________________
bool efficiencyCalculation::isTofmatched(TString particle, float charge, float pt) {
    if (!mTofMatchSet) return false;

    TF1* h = NULL;
    if (particle == "pion") {
        if (charge > 0) h = mfTofMatch[0];
        else h = mfTofMatch[1];
    }
    else {
        if (charge > 0) h = mfTofMatch[2];
        else h = mfTofMatch[3];
    }

    return gRandom->Rndm() < h->Eval(pt);
}

//____________________________________________________________________________________
void efficiencyCalculation::setTofPid() {
    TFile* fileTofPidKaon=new TFile(mFolderTofPid+"/f_tofPidEff_K.root", "READ");
    TFile* fileTofPidPion=new TFile(mFolderTofPid+"/f_tofPidEff_pi.root", "READ");

    if (!fileTofPidKaon || !fileTofPidPion) {
        cout<<"one of tof pid files not loaded."<<endl;
        return;
    }


    mfTofPidPion = (TF1*) fileTofPidPion->Get("f_f_tofPidEff_pi");
    mfTofPidKaon = (TF1*) fileTofPidKaon->Get("f_f_tofPidEff_K");

    mTofPidSet = true;

    return;
}

//____________________________________________________________________________________
void efficiencyCalculation::setTpcPid() {
    TFile* fileKaon=new TFile(mFolderTpcPid+"/f_tpc_K.root", "READ");
    TFile* filePion=new TFile(mFolderTpcPid+"/f_tpc_pi.root", "READ");
    if (!fileKaon || !filePion) {
        cout<<"one of tof pid files not loaded."<<endl;
        return;
    }

    mfTpcPidPion = (TF1*) filePion->Get("f_f_tpc_pi");
//    mfTpcPidPion->SetDirectory(0);
    mfTpcPidKaon = (TF1*) fileKaon->Get("f_f_tpc_K");
//    mfTpcPidKaon->SetDirectory(0);

//    fileKaon->Close();
//    filePion->Close();
    mTpcPidSet = true;
    return;
}

//____________________________________________________________________________________
bool efficiencyCalculation::isTofPid(TString particle, float pt) {
    if (!mTofPidSet) return false;

    TF1 *h = NULL;
    if (particle == "pion") h = mfTofPidPion;
    if (particle == "kaon") h = mfTofPidKaon;

    return gRandom->Rndm() < h->Eval(pt);
}

//____________________________________________________________________________________
bool efficiencyCalculation::isTpcPid(TString particle, float pt){
    if (!mTpcPidSet) return false;

    TF1* h = NULL;
    if (particle == "pion") h = mfTpcPidPion;
    if (particle == "kaon") h = mfTpcPidKaon;

    return gRandom->Rndm() < h->Eval(pt);
}