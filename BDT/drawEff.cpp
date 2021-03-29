//
// Created by lukas on 29.1.2019.
//
#include<iostream>
#include<fstream>
#include "TGraphErrors.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TPaveText.h"
#include <stdio.h>

using namespace std;

void drawEff(){

    std::vector<TString> inputFiles;
    std::vector<TFile*> inputFilesF;
    std::vector<TGraphErrors*> graphs;
//    std::vector<TString*> legendStrings;
    std::vector<const char *> legendStrings;

    TString names[] = {"pt_1_2/n150_d3/pt_1.0_2.0_nTrees_150.0_maxDepth_3.0",
                       "pt_2_3/n150_d3/pt_2.0_3.0_nTrees_150.0_maxDepth_3.0",
                       "pt_3_5/n400_d3/pt_3.0_5.0_nTrees_400.0_maxDepth_3.0"};


    TString leg[] = {"1<p_{T}<2 GeV/c",
                       "2<p_{T}<3 GeV/c",
                       "3<p_{T}<5 GeV/c"
    };

    for (int k = 0; k < 3; ++k) {
        inputFiles.push_back(Form("significance_%s.root", names[k].Data()));
        legendStrings.push_back(leg[k]);

    }

    Int_t colors[] = {1, 46, 9, 9, 40, 41, 42, 28, 2};
    Int_t markers[] = {25, 20, 34, 9, 40, 41, 42, 28, 2};

    Double_t x,y;
    for (unsigned short i = 0; i < inputFiles.size(); ++i) {
        TFile *in = TFile::Open(inputFiles[i]);
        inputFilesF.push_back(in);
        TGraphErrors *gr = (TGraphErrors*) inputFilesF[i]->Get(Form("gr_sign_%s",names[i].Data()));
        for (int j = 0; j < gr->GetN(); ++j) {
            gr->GetPoint(j,x,y);
            if (abs(y)>20) gr->RemovePoint(j);
            if (x>1) gr->RemovePoint(j);
            if (y!=y) gr->RemovePoint(j);
        }
        graphs.push_back(gr);

    }
    TCanvas *out = new TCanvas("out", "out", 1000, 800);

    TMultiGraph *mg = new TMultiGraph();
//    mg->GetYaxis()->SetTitle("n#sigma sigma [n#sigma]");
//    mg->GetYaxis()->SetTitle("n#sigma mean [n#sigma]");
    mg->GetYaxis()->SetTitle("Significance");
//    mg->GetYaxis()->SetTitle("TOF match efficiency");
    mg->GetXaxis()->SetTitle("Cut value applied on BDT output");
    mg->SetTitle("");

//    mg->SetTitle("Pions");
    TLegend *legend = new TLegend(0.126, 0.71, 0.277, 0.85);
    legend -> SetFillStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextSize(0.035);
    cout<<graphs.size()<<endl;
    for (unsigned short j = 0; j < graphs.size(); ++j) {
        graphs[j]->SetMarkerColor(colors[j]);
        graphs[j]->SetMarkerStyle(markers[j]);
        graphs[j]->SetMarkerSize(1.2);
        graphs[j]->SetLineColor(colors[j]);
        graphs[j]->GetYaxis()->SetRangeUser(0,10);
        graphs[j]->GetYaxis()->SetTitleOffset(0.9);
        graphs[j]->SetName(Form("%i",j));
        legend -> AddEntry(graphs[j], legendStrings[j], "p");
        mg->Add(graphs[j]);
    }

//    mg->SetMinimum(0);
    graphs[0]->GetYaxis()->SetLabelSize(0.045);
    graphs[0]->GetXaxis()->SetLabelSize(0.045);
    graphs[0]->GetXaxis()->SetTitleSize(0.055);
    graphs[0]->GetYaxis()->SetTitleSize(0.055);
    graphs[0]->GetXaxis()->CenterTitle();
    graphs[0]->GetYaxis()->CenterTitle();
    graphs[0]->Draw("ap");

    for (int l = 1; l < graphs.size(); ++l) {
        graphs[l]->Draw("p same");

    }



    gPad->Modified();

    TPaveText *text5 = new TPaveText(0.232,0.86,0.30,0.88,"brNDC");
    text5->SetTextSize(0.035);
    text5->SetLineColor(0);
    text5->SetShadowColor(0);
    text5->SetFillColor(0);
    text5->SetTextFont(42);
    text5->AddText("d+Au #sqrt{s_{NN}} = 200 GeV");
    text5->Draw("same");
    legend->Draw("same");
}