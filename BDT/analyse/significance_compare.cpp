//script to compare selected significance vs BDT plots
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TCanvas.h"


void significance_compare(Double_t ptmin = 1, Double_t ptmax = 2) {

    Double_t nTreesArr[] =    {200, 200, 200, 300, 300, 300, 400, 400, 400};
    Double_t treeDepthArr[] = {2,     3,   4,   2,   3,   4,   2,    3,   4};
    Int_t colors[] = {46, 8, 1, 38, 7, 9, 29, 42, 41};
    TCanvas *cSign = new TCanvas("cSign","cSign",1200,900);

    TLegend *legend = new TLegend(0.1269, 0.59, 0.27, 0.9);
    legend -> SetFillStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextSize(0.03);

    TGraph* gr[9];
    TString name;
    gStyle->SetPalette(kSolar);

    for (int i = 0; i < 6; ++i) {
        cout<<"Ntrees: "<<nTreesArr[i]<<" Depth: "<<treeDepthArr[i]<<endl;
        TFile *fSign = new TFile(Form("significance_pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f.root", ptmin, ptmax, nTreesArr[i], treeDepthArr[i]));
        fSign -> GetObject(Form("gr_sign_pt_%.1f_%.1f_nTrees_%.1f_maxDepth_%.1f", ptmin, ptmax, nTreesArr[i], treeDepthArr[i]), gr[i]);
        gr[i] -> SetLineColor(colors[i]);
        gr[i] -> SetMarkerColor(colors[i]);
        gr[i] -> SetMarkerStyle(20);
        gr[i] -> SetMarkerSize(0.03);
        name = Form("n%.0f_d%.0f", nTreesArr[i], treeDepthArr[i]);
        legend -> AddEntry(gr[i], name, "pl");

//        gr[i]->Draw("C");

        if (i == 0) gr[i]->Draw();
//        if (i == 0) gr[i]->Draw("PFC PMC PLC");
        else gr[i]->Draw("same");
        fSign -> Close();

    }

    legend -> Draw("same");








}