//
// Created by lukas on 2. 10. 2019.
//
#include "tmvaCuts.h"
#include "TMath.h"
#include "TFile.h"
#include "TAxis.h"
#include "TF1.h"
#include "TLine.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TNtuple.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TPaveText.h"
#include "TMultiGraph.h"
#include "TCollection.h"
#include <iostream>
#include <iostream>
#include <fstream>

void project_vars() {
    Double_t ptmin=0.5;
    Double_t ptmax=6;
    TCut mycuts = Form("D_mass > 1.815 && D_mass < 1.905 && D_pt>%1.2f && D_pt<%1.2f && k_pt>%1.2f && pi1_pt>%1.2f && "
                       "D_decayL>%f && D_decayL<0.2 && "
                       "dcaDaughters<%f && "
                       "k_dca>%f && k_dca<0.2 && "
                       "pi1_dca>%f && pi1_dca<0.2 && "
                       "dcaD0ToPv<%f && "
                       "cosTheta>%f",
                       ptmin, ptmax, tmvaCuts::minPt, tmvaCuts::minPt,
                       tmvaCuts::decayLength, tmvaCuts::dcaDaughters,
                       tmvaCuts::kDca, tmvaCuts::pDca,
                       tmvaCuts::dcaV0ToPv,
                       tmvaCuts::cosTheta);

    TString input = "/home/lukas/work/dmesons/Dmaker_ndAu/Dmaker_dAu/ntp/ntp.D0.0210.root";
    TFile *inputBckg = new TFile(input);

    TNtuple *bckgTree = (TNtuple*) inputBckg->Get("ntp_background");
    TNtuple *ntpSignal = (TNtuple*) inputBckg->Get("ntp_signal");

    TCanvas *out11 = new TCanvas("out11", "out11", 900, 1100);

    TH1F *hPtPi = new TH1F("hPtPi", "hPtPi", 30, -0.005, 30);
//    TH1F *hPtPi = new TH1F("hPtPi", "hPtPi", 1000000, 17000000, 18000000);
    hPtPi->Sumw2();
    hPtPi->SetMarkerStyle(20);
    hPtPi->SetMarkerColor(9);
    hPtPi->SetLineColor(9);
    hPtPi->SetStats(0);

//    TH1F *hPtK = new TH1F("hPtK", "hPtK", 120, 0, 1200);
    TH1F *hPtK = new TH1F("hPtK", "hPtK", 30, -0.005, 30);
//    TH1F *hPtK = new TH1F("hPtK", "hPtK", 200, -0.005, 0.195);
//    TH1F *hPtK = new TH1F("hPtK", "hPtK", 100000, 170000, 180000);
    hPtK->Sumw2();
    hPtK->SetMarkerStyle(21);
    hPtK->SetMarkerColor(46);
    hPtK->SetLineColor(46);
    hPtK->SetStats(0);
    TCut outHotSpot = "hotSpot<1";
//    ntpSignal->Project("hPtK", "k_dca", mycuts);
//    ntpSignal->Project("hPtPi", "k_dca", mycuts+outHotSpot);

    ntpSignal->Project("hPtK", "refMult", mycuts);
    ntpSignal->Project("hPtPi", "refMult", mycuts+outHotSpot);

    hPtK->Scale(1/hPtK->GetEntries());
    hPtPi->Scale(1/hPtPi->GetEntries());
    cout<<ntpSignal->GetEntries(mycuts+outHotSpot)<<endl;
    cout<<ntpSignal->GetEntries(mycuts)<<endl;

    cout<<bckgTree->GetEntries(mycuts+outHotSpot)<<endl;
    cout<<bckgTree->GetEntries(mycuts)<<endl;

    ntpSignal->Draw("k_pt:pi1_pt");
    hPtK->Draw();
    hPtPi->Draw("same");
//    out11->SaveAs("ooutt.png");

    TH2F *h2d = new TH2F("h2d", "h2d", 1000, 0, 10, 1000, 0, 10);
    ntpSignal->Project("h2d","k_pt:pi1_pt", mycuts);
    gPad->SetRightMargin(0.15);
    h2d->GetYaxis()->SetTitle("kaon p_{T} (GeV/c)");
    h2d->GetXaxis()->SetTitle("pion p_{T} (GeV/c)");
    h2d->GetXaxis()->SetTitleOffset(1.1);
    h2d->SetStats(0);
    h2d->SetTitle("");
//    h2d->Draw("colz");

    TPaveText *text4 = new TPaveText(0.35, 0.848, 0.229, 0.873, "brNDC");
    text4->SetTextSize(0.03);
    text4->SetLineColor(0);
    text4->SetShadowColor(0);
    text4->SetFillColor(0);
    text4->AddText("D^{0} candidates pre-tuning cuts");
    text4->AddText(Form("%.1f < pair p_{T} < %.1f", ptmin, ptmax));
    text4->Draw("same");
}
