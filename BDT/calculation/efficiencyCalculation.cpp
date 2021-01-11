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
        hTpcPiPlus[iCent]->Draw();
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
    TFile* fTofMatch=new TFile(mFolderTofMatch+"tofMatchEfficiency.root", "READ");
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
    TFile* fileTofPidKaon=new TFile(mFolderTofPid+"tofPidEff_K.root", "READ");
    TFile* fileTofPidPion=new TFile(mFolderTofPid+"tofPidEff_pi.root", "READ");

//    TFile* fileTofPidKaon=new TFile(mFolderTofPid+"results_KK.root", "READ");
//    TFile* fileTofPidPion=new TFile(mFolderTofPid+"results_pipi.root", "READ");

    if (!fileTofPidKaon || !fileTofPidPion) {
        cout<<"one of tof pid files not loaded."<<endl;
        return;
    }


    mfTofPidPion = (TF1*) fileTofPidPion->Get("f_tofPidEff_pi");
    mfTofPidKaon = (TF1*) fileTofPidKaon->Get("f_tofPidEff_K");

    //    TGraphErrors *gPion = (TGraphErrors*) fileTofPidPion->Get("eff");
//    TGraphErrors *gKaon = (TGraphErrors*) fileTofPidKaon->Get("eff");

    /*
    mfTofPidPion = new TF1("tof_pif_fit_pion", "[0]+[1]/x+[2]*x+[3]/x/x", 0.15, 5);
    mfTofPidPion->SetParameters(1, -0.06, -0.1, 0.02, 0.006);
    mfTofPidPion->SetParLimits(0, 0, 1);
    mfTofPidPion->SetParLimits(1, 0, 1);
    mfTofPidPion->SetParLimits(2, 0, 1);
    mfTofPidPion->SetParLimits(3, -1, 0);
    TCanvas* canPionFit=new TCanvas("cPionTof", "cPionTof", 1000, 1000);
    gPion->Fit(mfTofPidPion, "REX0W");
    gPion->Draw("ap");

    mfTofPidKaon = new TF1("tof_pif_fit_kaon", "[0]+[1]/x+[2]*x+[3]/x/x", 0.15, 5);
    mfTofPidKaon->SetParameters(1, -0.06, -0.1, 0.02, 0.006);

    TCanvas* canKaonFit=new TCanvas("cKaonFit", "cKaonFit", 1000, 1000);
    gKaon->Fit(mfTofPidKaon, "REX0W");
    gKaon->Draw("ap");

    fileTofPidKaon->Close();
    fileTofPidPion->Close();
    */

    mTofPidSet = true;

    return;
}

//____________________________________________________________________________________
void efficiencyCalculation::setTpcPid() {
    TFile* fileKaon=new TFile(mFolderTpcPid+"tpc_K.root", "READ");
    TFile* filePion=new TFile(mFolderTpcPid+"tpc_pi.root", "READ");
    if (!fileKaon || !filePion) {
        cout<<"one of tof pid files not loaded."<<endl;
        return;
    }

    mfTpcPidPion = (TF1*) filePion->Get("f_tpc_pi");
    mfTpcPidKaon = (TF1*) fileKaon->Get("f_tpc_K");


    /*
//    TGraphErrors *gPion = (TGraphErrors*) filePion->Get("eff");
//    TGraphErrors *gKaon = (TGraphErrors*) fileKaon->Get("eff");

    mfTpcPidPion = new TF1("tpc_pif_fit_pion", "[0]+[1]/x+[2]*x+[3]/x/x", 0.15, 5);
    mfTpcPidPion->SetParameters(1, -0.06, -0.1, 0.02, 0.006);
    mfTpcPidPion->SetParLimits(0, 0, 1);
    mfTpcPidPion->SetParLimits(1, 0, 1);
    mfTpcPidPion->SetParLimits(2, 0, 1);
    mfTpcPidPion->SetParLimits(3, -1, 0);
    TCanvas* canPionFit=new TCanvas("cPionTof", "cPionTof", 1000, 1000);
    gPion->Fit(mfTpcPidPion, "REX0W");
    gPion->Draw("ap");

    mOutFileTpcPidPion->cd();
    mfTpcPidPion->Write();

    mfTpcPidKaon = new TF1("tpc_pif_fit_kaon", "[0]+[1]/x+[2]*x+[3]/x/x", 0.15, 5);
    mfTpcPidKaon->SetParameters(1, -0.06, -0.1, 0.02, 0.006);

    TCanvas* canKaonFit=new TCanvas("cKaonFit", "cKaonFit", 1000, 1000);
    gKaon->Fit(mfTpcPidKaon, "REX0W");
    gKaon->Draw("ap");
    mOutFileTpcPidKaon->cd();
    mfTpcPidKaon->Write();


    fileKaon->Close();
    filePion->Close();
    */

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