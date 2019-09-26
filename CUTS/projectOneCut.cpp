/*
 * Application of trained CUTS, however, in more "physics" way - e.g. not using k_dca_max,...
*/
#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLatex.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include <iostream>
#include <ctime>
#include <cstdio>
#include <fstream>
#include "FitD0Peak.hh"
#include "TSystem.h"


using namespace std;
using namespace TMath;

TH1F* hS;
TH1F* hB;

//void project(TString input, float ptmin, float ptmax, float decayL, float dcaDaughters, float dcaD0, float cos, float kdca, float pidca) {
void project(TString input, TString inputFile, float ptmin, float ptmax, TString setCuts) {
//    TString ptCut = Form("D_pt>%f && D_pt<=%f", ptmin, ptmax);
//    TString decayLCut = Form("D_decayL>%f",decayL);
//    TString dcaDaughtersCut = Form("dcaDaughters<%f", dcaDaughters);
//    TString dcaD0Cut = Form("dcaD0ToPv<%f", dcaD0);
//    TString cosCut = Form("cosTheta>%f", cos);
//    TString kdcaCut = Form("k_dca>%f",kdca);
//    TString pidcaCut = Form("pi1_dca>%f", pidca);
//    TString setCuts = ptCut+" && "+decayLCut+" && "+dcaDaughtersCut+" && "+dcaD0Cut+" && "+cosCut+" && "+kdcaCut+" && "+pidcaCut+" && k_pt>0.15 && pi1_pt>0.15";
//    cout<<setCuts<<endl;
    gSystem->Exec(Form("mkdir %s", inputFile.Data()));

    TString nameSig = Form("hS_%.1f_%.1f", ptmin, ptmax);
    TString nameBack = Form("hB_%.1f_%.1f", ptmin, ptmax);
    cout<<nameBack<<endl;
    TFile* data = new TFile(input ,"r");
    TNtuple* ntpS = (TNtuple*)data -> Get("ntp_signal");
    TNtuple* ntpB = (TNtuple*)data -> Get("ntp_background");

    TString outFileName="res_oneCut"+inputFile+".root";

    TFile *fOut = new TFile(outFileName,"recreate");
    hS = new TH1F(nameSig, nameSig, 2000, 0.4, 2.4);
    hB = new TH1F(nameBack, nameBack, 2000, 0.4, 2.4);

    ntpS -> Project(nameSig, "D_mass", setCuts);
    ntpB -> Project(nameBack, "D_mass", setCuts);

//    hS->SetStats(0);
//    hB->SetStats(0);

    gROOT->ProcessLine(".L FitD0Peak.cpp++");

    FitD0Peak *fitmass = new FitD0Peak(hS, hB, ptmin, ptmax, outFileName);
    cout<<fitmass->getSignificance()<<endl;
    cout<<fitmass->getRawYield()<<" "<<fitmass->getRawYieldError()<<endl;

//    TCanvas* publ = (TCanvas*)fOut -> Get("2.0_3.0_publ");
//    publ->SaveAs("aaaa.png");

    delete fitmass;
    delete hS;
    delete hB;
//    hS -> Write(nameSig);
//    hB -> Write(nameBack);


    gSystem->Exec(Form("mv *.png %s", inputFile.Data()));

    data->Close();
    fOut->Close();
}



void projectOneCut(){
    TString folder="/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/";

//    TString inputEff = "/home/lukas/work/tmva_d0/CUTS/pt_3/pass0/TMVA_cuts_pt_3.0_5.0.root";
    TString inputEff = "/home/lukas/work/tmva_d0/CUTS/pt_2/pass0/TMVA_cuts_pt_2.0_3.0.root";



//    TString inputEff = "/home/lukas/work/tmva_d0/CUTS/pt_1/pass0/TMVA_cuts_pt_1.0_2.0.root";
//    TString inputPairs = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.KF.D0.1306.root";
//    TString inputPairs = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.D0.KF.dcadaughters007.costheta09.USLS.1706.root";
//    TString inputPairs = "ntp.D0.KF.dcaD007.cos09.dca009.2507.root";
//    TString inputPairs = "ntp.D0.KF.rename.root";
//    TString inputPairs = "ntp.KF.D0.dcaD007.costheta09.dca012.2307.root";
//    TString inputPairs = "ntp.KF.D0.1306.root";
//    TString inputPairs = "ntp.D0.KF.dcadaughters07.cos09.US.root";
//    TString inputPairs = "ntp.D0.KFnew.1607.root";
//    TString inputPairs = "ntp.D0.dcaDaughers007.costheta.US.1006.root";
//    TString inputPairs = "ntp.D0.dca012.prim.all.2406.root";
//    TString inputPairs = "ntp.D0.2705.root";
//    TString inputPairs = "ntp.D0.2805.dca009.root";
//    TString inputPairs = "ntp.D0.KF.dcad007.cos09.dca009.US.mass1.81-1.92.root";
//    TString inputPairs = "ntp.D0.KF.dcaD007.cos09.pairDCA008.dca012.US.mass1.81-1.92.root";
    TString inputPairs = "ntp.D0.KF.dcaD007.cos09.pairDCA008.dca010.US.mass1.81-1.92.root";
//    TString inputPairs = "ntp.D0.KF.US.dcaD007.cos09.dca012.mass1.81-1.92.root";
//    TString inputPairs = "ntp.D0.KF.dcadaughters007.costheta09.USLS.1706.root";
//    TString inputPairs = "ntp.D0.2506.root";
//    TString inputPairs = "ntp.D0.KF.dcadaughter007.costheta09.USLS.refitOnlyRem.1906.root";
//    TString inputPairs = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.D0.KF.dcadaughter007.costheta09.USLS.refitOnlyRem.1906.root";
//    TString inputPairs = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.D0.KF.dcadaughters07.cos09.US.root";
//    TString inputPairs = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.picoD0AnaMaker.0802.0415.root";
//    TString inputPairs = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.D0.dcaDaughers007.costheta.US.1006.root";
//    TString inputPairs = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.D0.KFvertex.1106.root";
//    TString inputPairs = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.D0.dca009.globalRem.root";
//    TString inputPairs = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.D0.0406.dca015.root";
//    TString inputPairs = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.D0.2805.dca009.root";

    Float_t selectedEff=0.2550002;
//    Float_t selectedEff=0.2050002;
    Float_t ptMin=2;
    Float_t ptMax=3;

    TFile* dataEff = new TFile(inputEff ,"r");

    bool cutFound=false;

//    TTree ntp[2] = {(TTree*)input->Get("ntp_signal"), (TTree*)input->Get("ntp_background")};
    TTree* ntpCuts = (TTree*)dataEff->Get("ntpCuts");

    Float_t effS, k_dca_min, k_dca_max, pi1_dca_min, pi1_dca_max, dcaDaughters_min, dcaDaughters_max, dcaD0ToPv_min, dcaD0ToPv_max, decayL_min, decayL_max, cosTheta_min;

    decayL_min=0.0005;
    decayL_max=2000;
    cosTheta_min=0.9;

    ntpCuts->SetBranchAddress("effS", &effS);
//    ntpCuts->SetBranchAddress("k_dca_min", &k_dca_min);
//    ntpCuts->SetBranchAddress("k_dca_max", &k_dca_max);
//    ntpCuts->SetBranchAddress("pi1_dca_min", &pi1_dca_min);
//    ntpCuts->SetBranchAddress("pi1_dca_max", &pi1_dca_max);
//    ntpCuts->SetBranchAddress("dcaDaughters_min", &dcaDaughters_min);
//    ntpCuts->SetBranchAddress("dcaDaughters_max", &dcaDaughters_max);
//    ntpCuts->SetBranchAddress("dcaD0ToPv_min", &dcaD0ToPv_min);
//    ntpCuts->SetBranchAddress("dcaD0ToPv_max", &dcaD0ToPv_max);

    ntpCuts->SetBranchAddress("k_dca_min", &dcaDaughters_min);
    ntpCuts->SetBranchAddress("k_dca_max", &dcaDaughters_max);
    ntpCuts->SetBranchAddress("pi1_dca_min", &dcaD0ToPv_min);
    ntpCuts->SetBranchAddress("pi1_dca_max", &dcaD0ToPv_max);
    ntpCuts->SetBranchAddress("dcaDaughters_min", &k_dca_min);
    ntpCuts->SetBranchAddress("dcaDaughters_max", &k_dca_max);
    ntpCuts->SetBranchAddress("dcaD0ToPv_min", &pi1_dca_min);
    ntpCuts->SetBranchAddress("dcaD0ToPv_max", &pi1_dca_max);


    Int_t nEffCuts=ntpCuts->GetEntries();
    for (int i = 0; i < nEffCuts; ++i) {
        ntpCuts->GetEntry(i);
        if (abs(effS-selectedEff)<0.001) {
            cout<<"effS cut found "<<effS<<endl;
            cutFound=true;

            cout<<"dcaD0 "<<dcaD0ToPv_min<<" "<<dcaD0ToPv_max<<endl;
            cout<<"piDca "<<pi1_dca_min<<" "<<pi1_dca_max<<endl;
            cout<<"kDca "<<k_dca_min<<" "<<k_dca_max<<endl;
            cout<<"dcaDaughters "<<dcaDaughters_min<<" "<<dcaDaughters_max<<endl;
            TString ptCut = Form("D_pt>%f && D_pt<%f", ptMin, ptMax);

//            TMVA:
//            TString decayLCut = Form("D_decayL>=%f && D_decayL<=%f",decayL_min, decayL_max);
//            TString dcaDaughtersCut = Form("dcaDaughters>=%f && dcaDaughters<=%f", dcaDaughters_min, dcaDaughters_max);
//            TString dcaD0Cut = Form("dcaD0ToPv>=%f && dcaD0ToPv<=%f", dcaD0ToPv_min, dcaD0ToPv_max);
//            TString cosCut = Form("cosTheta>=%f", cosTheta_min);
//            TString kdcaCut = Form("k_dca>=%f && k_dca<=%f",k_dca_min, k_dca_max);
//            TString pidcaCut = Form("pi1_dca>=%f && pi1_dca<=%f", pi1_dca_min, pi1_dca_max);

//            TString decayLCut = Form("D_decayL>%f && D_decayL<%f",decayL_min, decayL_max);
//            TString dcaDaughtersCut = Form("dcaDaughters>%f && dcaDaughters<%f", dcaDaughters_min, dcaDaughters_max);
//            TString dcaD0Cut = Form("dcaD0ToPv>%f && dcaD0ToPv<%f", dcaD0ToPv_min, dcaD0ToPv_max);
//            TString cosCut = Form("cosTheta>%f", cosTheta_min);
//            TString kdcaCut = Form("k_dca>%f && k_dca<%f",k_dca_min, k_dca_max);
//            TString pidcaCut = Form("pi1_dca>%f && pi1_dca<%f", pi1_dca_min, pi1_dca_max);

            //Physics:
            TString decayLCut = Form("D_decayL>%f",decayL_min);
            TString dcaDaughtersCut = Form("dcaDaughters<%f", dcaDaughters_max);
            TString dcaD0Cut = Form("dcaD0ToPv<%f", dcaD0ToPv_max);
            TString cosCut = Form("cosTheta>%f", cosTheta_min);
            TString kdcaCut = Form("k_dca>%f",k_dca_min);
            TString pidcaCut = Form("pi1_dca>%f", pi1_dca_min);

            TString setCuts = ptCut+" && "+decayLCut+" && "+dcaDaughtersCut+" && "+dcaD0Cut+" && "+cosCut+" && "+kdcaCut+" && "+pidcaCut+" && k_pt>0.15 && pi1_pt>0.15";
//
//            TString setCuts = kdcaCut+" && "+pidcaCut;
//            TString setCuts = dcaDaughtersCut+" && "+cosCut;

//            TString setCuts = ptCut+" && "+decayLCut+" && "+dcaDaughtersCut+" && "+dcaD0Cut+" && "+cosCut+" && "+kdcaCut+" && "+pidcaCut;
//            setCuts=pidcaCut;
//            setCuts=pidcaCut;
            cout<<setCuts<<endl;

//            project(inputPairs, 2, 3, 0.01, dcaDaughters_max, dcaD0ToPv_max, 0.7, k_dca_min, pi1_dca_min);
            project(folder+inputPairs, inputPairs, ptMin, ptMax, setCuts);
        }

    }




}


