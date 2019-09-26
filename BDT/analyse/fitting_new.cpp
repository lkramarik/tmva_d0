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

using namespace std;
using namespace TMath;

void fitting_new(TString input = "D0_bdt_cuts_pt_2.0_3.0_nTrees_350.0_maxDepth_4.0.root", bool scale = false, bool subtract = true,  Double_t ptminText = 1, Double_t ptmaxText = 5, float bdtcut = 0.3536) {
//void fitting_new(TString input = "/home/lukas/work/tmva_d0/BDT/pt_2_3_scan/n300_d4/analyse/D0_bdt_cuts_pt_2.0_3.0_nTrees_300.0_maxDepth_4.0.root", bool scale = false, bool subtract = true,  Double_t ptminText = 1, Double_t ptmaxText = 5, float bdtcut = 0.3536) {
    std::cout << input << endl;
    TFile *data = new TFile(input, "r");
    TList *listS = (TList *) data->Get("hists_signal");
    TList *listB = (TList *) data->Get("hists_background");
    TList *listME = (TList *) data->Get("hists_ME");

    //just number of pairs:
    std::cout << Form("hB_%.1f_%.1f", ptminText, ptmaxText) << endl;

    TH1F *hInvMassBack = (TH1F*) listB->FindObject(Form("D_mass_%.4f_background", bdtcut));
    TH1F *hInvMassSign = (TH1F*) listS->FindObject(Form("D_mass_%.4f_signal", bdtcut));
    TH1F *hInvMassMxd = (TH1F*) listME->FindObject(Form("D_mass_%.4f_ME", bdtcut));

//    fit(hInvMassSign, hInvMassBack, ptminText,  ptmaxText, scale, subtract, input, bdtcut);

    gROOT->ProcessLine(".L FitD0Peak.cpp++");

    FitD0Peak *fitmass = new FitD0Peak(hInvMassSign, hInvMassBack, 2, 3, "test.root");
    fitmass->addMixedEventBckg(hInvMassMxd);
    fitmass->doStuff();
    cout<<fitmass->getSignificance()<<endl;
    cout<<fitmass->getRawYield()<<" "<<fitmass->getRawYieldError()<<endl;
    delete fitmass;

    data->Close();

    delete data;
    delete listB;
    delete listS;
    delete hInvMassBack;
    delete  hInvMassSign;
}