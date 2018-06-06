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
#include <fitting.C>

void project_bdt(Double_t ptmin = 2, Double_t ptmax = 3) {
    gROOT->LoadMacro("fitting.C++");

    TFile *f = new TFile(Form("D0_bdt_cuts_pt_%.1f_%.1f.root", ptmin, ptmax),"RECREATE");
    TString input = Form("../out_local.root", ptmin, ptmax);

    const int n_bin = 21;
    float bdtRange[] = {0., 0.025, 0.03, 0.05, 0.06, 0.075, 0.08, 0.09, 0.1, 0.115, 0.125, 0.135, 0.15, 0.175, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};

    TH1F *his[n_bin][2];
    TString name[2] = {"signal", "background"};

    for (int jj = 0; jj < 2; ++jj) {
        for (int j = 0; j < n_bin; j++) {
            his[j][jj] = new TH1F(Form("D_mass_%.4f_", bdtRange[j]) + name[jj], "D_mass;m [GeV]", 2000, 0.4, 2.4);
            his[j][jj] -> Sumw2();
        }
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
            for(int bin = 0; bin < n_bin; bin++) {
                if (BDTresponse >= bdtRange[bin])  his[bin][k]->Fill(D_mass);
            }
        }

        f -> cd();
        TList *listOut = new TList();
        for (int i = 0; i < n_bin; i++) {
            listOut->Add(his[i][k]);
        }
        listOut->Write("hists_"+name[k], 1, 0);
        delete listOut;
    }


    Double_t y[n_bin], x[n_bin], stat[n_bin];

    TH1F *hisSubtr[n_bin];
    TList *listOut = new TList();

//    for (int l = 0; l < 1; ++l) {
    for (int l = 0; l < n_bin; ++l) {
        hisSubtr[l] = new TH1F(Form("D_mass_%.4f", bdtRange[l]), "D_mass;m [GeV]", 2000, 0.4, 2.4);
        hisSubtr[l] = (TH1F*)his[l][0]->Clone();
        stat[l] = his[l][0] -> Integral(his[l][0] -> FindBin(1.7), his[l][0] -> FindBin(2));
        hisSubtr[l] -> Add(his[l][1], -1);
        hisSubtr[l] -> Rebin(10);
        listOut -> Add(hisSubtr[l]);
        //        cout<<fit(his[l][0], his[l][1], ptmin, ptmax, false, true, "bd")<<endl;
        y[l] = fit(his[l][0], his[l][1], ptmin, ptmax, false, true, "bdt.root", bdtRange[l]);
        x[l] = bdtRange[l];
    }

    f -> cd();
    listOut -> Write("hists_D_mass", 1, 0);
    data -> Close();
    f -> Close();


    gStyle->SetMarkerStyle(2);
    gStyle->SetOptFit(1);
    gStyle->SetStatY(0.899);
    gStyle->SetStatX(0.9);

    TCanvas *cSign = new TCanvas("cSign","cSign",1200,900);
    TGraph* gr = new TGraph(n_bin, x, y);
    gr -> SetMarkerColor(9);
    gr->SetTitle("");
    gr->GetXaxis()->SetLabelSize(0.04);
    gr->GetYaxis()->SetLabelSize(0.04);
    gr->GetXaxis()->SetTitleSize(0.045);
    gr->GetYaxis()->SetTitleSize(0.045);
    gr->GetYaxis()->SetTitle("Significance");
    gr->GetYaxis()->SetTitleOffset(1.1);
    gr->GetXaxis()->SetTitle("BDT cut");
    gr->SetMarkerSize(2);
//    TGraphErrors* gr = new TGraph(n_bin, x, y);
    gr -> Draw("ap");
    TLatex txR;
    txR.SetNDC();
    txR.SetTextSize(0.04);
    txR.DrawLatex(0.1,0.93,Form("p_{T}: %3.1f-%3.1f GeV/c", ptmin, ptmax));

    TPaveText *text1 = new TPaveText(0.66,0.91,0.9391,0.967,"brNDC");
    text1 -> SetTextSize(0.04);
    text1 -> SetShadowColor(0);
    text1 -> SetLineColor(0);
    text1 -> SetFillColor(0);
    text1 -> AddText("d+Au 200 GeV");
    text1 -> Draw("same");

    TCanvas *cStat = new TCanvas("cStat","cStat",1200,900);
    TGraph* grStat = new TGraph(n_bin, x, stat);
    grStat -> SetMarkerColor(9);
    grStat->SetTitle("");
    grStat->GetXaxis()->SetLabelSize(0.04);
    grStat->GetYaxis()->SetLabelSize(0.04);
    grStat->GetXaxis()->SetTitleSize(0.045);
    grStat->GetYaxis()->SetTitleSize(0.045);
    grStat->GetYaxis()->SetTitle("Statistics (1.7 - 2.0 GeV/c^{2})");
    grStat->GetYaxis()->SetTitleOffset(1.1);
    grStat->GetXaxis()->SetTitle("BDT cut");
    grStat->SetMarkerSize(2);
//    TGraphErrors* gr = new TGraph(n_bin, x, y);
    grStat -> Draw("ap");
    text1 -> Draw("same");
    txR.DrawLatex(0.1,0.93,Form("p_{T}: %3.1f-%3.1f GeV/c", ptmin, ptmax));



}
