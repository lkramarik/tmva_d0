#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TROOT.h>
#include <TH1.h>
#include <TFile.h>
#include <iostream>
#include <ctime>
#include <cstdio>
#include <fstream>
#include "FitD0Peak.hh"

void project_bdt_oneCut(Double_t ptmin = 2, Double_t ptmax = 3, Double_t nTrees = 350, Double_t maxDepth = 4, Double_t bdtRange = 0.3, TString input = "../out_local.root") {
    bool mixed=false;
    TFile *f = new TFile(Form("D0_bdt_%.3f_pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f.root", bdtRange, ptmin, ptmax, nTrees ,maxDepth),"RECREATE");

    TH1F *his[2];
    TString name[2] = {"signal", "background"};

    for (int jj = 0; jj < 2; ++jj) {
        his[jj] = new TH1F(Form("D_mass_%.4f_", bdtRange) + name[jj], "D_mass;m [GeV]", 2000, 0.4, 2.4);
        his[jj] -> Sumw2();
    }

    TFile* data = new TFile(input ,"r");
    TNtuple* ntp[2] = {(TNtuple*)data -> Get("ntp_"+name[0]), (TNtuple*)data -> Get("ntp_"+name[1])};

    float D_mass, D_pt, BDTresponse;

    for (int k = 0; k < 2; ++k) {
        ntp[k]->SetBranchAddress("D_mass", &D_mass);
        ntp[k]->SetBranchAddress("D_pt", &D_pt);
        ntp[k]->SetBranchAddress("BDTresponse", &BDTresponse);
        Long64_t nentries = ntp[k]->GetEntriesFast();

        for (Long64_t jentry=0; jentry<nentries; jentry++) {
            ntp[k] -> GetEntry(jentry);
            if (D_pt>ptmin && D_pt<ptmax) {
                if (BDTresponse >= bdtRange) his[k]->Fill(D_mass);
            }
        }

        f -> cd();
        TList *listOut = new TList();
        listOut->Add(his[k]);

        listOut->Write("hists_"+name[k], 1, 0);
        delete listOut;
    }

    if (mixed){
        TH1F *hisMxd = new TH1F(Form("D_mass_%.4f_ME", bdtRange), "D_mass;m [GeV]", 2000, 0.4, 2.4);
        hisMxd -> Sumw2();

        TFile* dataMxd = new TFile("../out_local_mix.root" ,"r");
        TNtuple* ntpMxd = (TNtuple*)dataMxd -> Get("ntp_ME");

        ntpMxd->SetBranchAddress("D_mass", &D_mass);
        ntpMxd->SetBranchAddress("D_pt", &D_pt);
        ntpMxd->SetBranchAddress("BDTresponse", &BDTresponse);
        Long64_t nentries = ntpMxd->GetEntriesFast();

        for (Long64_t jentry=0; jentry<nentries; jentry++) {
            ntpMxd -> GetEntry(jentry);
                if (D_pt>ptmin && D_pt<ptmax) {
                    if (BDTresponse >= bdtRange) hisMxd->Fill(D_mass);
                }
        }
        f -> cd();
        TList *listOut = new TList();
        listOut->Add(hisMxd);
        listOut->Write("hists_ME", 1, 0);
        delete listOut;
    }

    TH1F *hisSubtr;
    TList *listOut = new TList();
    TString outFile;

    hisSubtr = new TH1F(Form("D_mass_%.4f", bdtRange), "D_mass;m [GeV]", 2000, 0.4, 2.4);
    hisSubtr = (TH1F*)his[0]->Clone();
    Double_t stat = his[0] -> Integral(his[0] -> FindBin(1.7), his[0] -> FindBin(2));
    hisSubtr -> Add(his[1], -1);
    hisSubtr -> Rebin(10);
    listOut -> Add(hisSubtr);

    outFile=Form("signals_%.4fbdt_pt_%.1f_%.1f.root", bdtRange, ptmin, ptmax);
    FitD0Peak *fitmass = new FitD0Peak(his[0], his[1], ptmin, ptmax, outFile);
    fitmass->doStuff();
//        y[l] = fitmass->getSignificance();
//        rawYields[l] = fitmass->getRawYield();
//        rawYieldsE[l] = fitmass->getRawYieldError();
    delete fitmass;

    f -> cd();
    listOut -> Write("hists_D_mass", 1, 0);
    data -> Close();
    f -> Close();
}