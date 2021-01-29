//
// Created by lukas on 04.03.20.
//
#include "TH1.h"
#include "TRatioPlot.h"
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
#include "TGaxis.h"
#include "TDatime.h"
#include <iostream>
#include <iostream>
#include <fstream>

void projectFilesBig(TString* inputPlot, int nInputFiles);
void plot(TString* input, Int_t nInputFiles, TString* legendNames, TString* grNames, Int_t nGrNames);
TH1F* graphToHisto(TGraphErrors *gr);
void efficiencyCompareEfficiency();
void efficiencyCompareResolution();
TString folderDate;

//__________________________________________________________________________________________________________________
TH1F* graphToHisto(TGraphErrors *gr){
    Int_t const nPoints = gr->GetN();
    Double_t pT[nPoints],y[nPoints],yE[nPoints];
    Float_t xAxis[nPoints+1];

    for (int i = 0; i < nPoints; ++i) {
        gr->GetPoint(i, pT[i], y[i]);
        yE[i] = gr->GetErrorY(i);
        xAxis[i]=pT[i]-gr->GetErrorX(i);
    }
    xAxis[nPoints]=pT[nPoints-1]+gr->GetErrorX(nPoints-1);

    TString name = gr->GetName();
    name+="_histo";
    TH1F* h = new TH1F(name, gr->GetTitle(), nPoints, xAxis);
//    h ->SetDirectory(0);

    for (int j = 0; j < nPoints; ++j) {
        h->SetBinContent(j+1, y[j]);
        h->SetBinError(j+1, yE[j]);
    }

    int color = gr->GetMarkerColor();
    Int_t style = gr->GetMarkerStyle();
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->SetLineWidth(2);
    h->SetStats(0);
    h->SetTitle("");
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.9);
//    h->SetMarkerSize(gr->GetMarkerSize());
    h->Sumw2();


    return h;
};

//__________________________________________________________________________________________________________________
void plot(TString* input, const int nInputFiles, TString* legendNames, TString* grNames, const int nGrNames) {
    cout << nInputFiles << " " << nGrNames << endl;
    Int_t colors[] = {1, 46, 9, 8, 6, 40, 42, 28, 2};
    Int_t markers[] = {25, 20, 34, 9, 40, 41, 42, 28, 2};

    gSystem->Exec(Form("rm finalAnalysis/%s/final_efficiencies_SIM.root", folderDate.Data()));
    TFile *fSimEfficiency = new TFile(Form("finalAnalysis/%s/final_efficiencies_SIM.root", folderDate.Data()), "RECREATE");

    TGraphErrors *gBackground = new TGraphErrors();
    TGraphErrors *gr = new TGraphErrors();
    std::vector < TGraphErrors * > graphsSIM[nGrNames];

    for (int m = 0; m < nInputFiles; ++m) {
        TFile *fInEff = new TFile(Form("finalAnalysis/%s/%s/final_result_SIM.root", folderDate.Data(), input[m].Data()), "READ");
        if (!fInEff) continue;
        for (int k = 0; k < nGrNames; ++k) {
            gr = (TGraphErrors *) fInEff->Get(grNames[k]);
            graphsSIM[k].push_back(gr);
        }
        fInEff->Close();
    }

    TLegend *legend = new TLegend(0.136, 0.66, 0.34, 0.88);
    legend->SetFillStyle(0);
    legend->SetLineColor(0);
    legend->SetTextSize(0.04);

    for (int l = 0; l < nGrNames; ++l) {
        TCanvas *out1 = new TCanvas("out1", "out1", 1200, 1000);
        out1->SetLogy();
        out1->SetGrid();
        TMultiGraph *mg = new TMultiGraph();
        mg->SetMaximum(0.7);
        mg->GetYaxis()->SetTitle("#varepsilon");
        mg->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
        mg->GetXaxis()->CenterTitle(kTRUE);
        mg->GetXaxis()->SetTitleSize(0.96);
        //    mg->GetYaxis()->SetRangeUser(0.000001,2);
        mg->SetTitle("");
        for (unsigned short j = 0; j < graphsSIM[l].size(); ++j) {// nfiles
            graphsSIM[l][j]->SetMarkerColor(colors[j]);
            graphsSIM[l][j]->SetMarkerStyle(2);
            graphsSIM[l][j]->SetMarkerSize(1.7);
            graphsSIM[l][j]->SetLineColor(colors[j]);
            graphsSIM[l][j]->SetLineWidth(2);
            graphsSIM[l][j]->SetName(Form("%i", j));
            if (l == 0) legend->AddEntry(graphsSIM[l][j], legendNames[j], "pl");
            mg->Add(graphsSIM[l][j]);
            fSimEfficiency->cd();
            graphsSIM[l][j]->Write(grNames[l]);
        }
        mg->Draw("ap");
        TPaveText *text5 = new TPaveText(0.724, 0.925, 0.793, 0.945, "brNDC");
        text5->SetTextSize(0.035);
        text5->SetLineColor(0);
        text5->SetShadowColor(0);
        text5->SetFillColor(0);
        text5->SetTextFont(42);
        text5->AddText("d+Au #sqrt{s_{NN}} = 200 GeV");
//    mg->GetYaxis()->SetRangeUser(0.,20.);
        text5->Draw("same");
        legend->Draw("same");
        TString imgName = "finalAnalysis/efficiencyRatio/" ;
        imgName += grNames[l];
        imgName += ".png";
        out1->SaveAs(imgName);

        imgName = "finalAnalysis/efficiencyRatio/" ;
        imgName += grNames[l];
        imgName += ".eps";
        out1->SaveAs(imgName);

        //----------------------------------Average value error------------------------------------
        TLegend *legend2 = new TLegend(0.188377, 0.779, 0.34, 0.9186);
        legend2->SetFillStyle(0);
        legend2->SetLineColor(0);
        legend2->SetTextSize(0.03);
        const int nPoints = graphsSIM[0][0]->GetN();
        Double_t x[nPoints], y[nPoints], xE[nPoints], yE[nPoints], xTemp, yTemp;

        for (int i = 0; i < graphsSIM[0][0]->GetN(); ++i) { //all of the points
            x[i]=0;
            y[i]=0;
            xE[i]=0;
            yE[i]=0;
            for (int j = 0; j < nInputFiles; ++j) { //all input files
                graphsSIM[l][j]->GetPoint(i,xTemp,yTemp);
                x[i]=xTemp;
                xE[i]=graphsSIM[l][j]->GetErrorX(i);
                y[i]+=yTemp;
            }
            y[i]/=nInputFiles; //average

            for (int j = 0; j < nInputFiles; ++j) { //all input files
                graphsSIM[l][j]->GetPoint(i,xTemp,yTemp);
                yE[i]+=pow(yTemp-y[i], 2);
            }
            yE[i]/=(nInputFiles*(nInputFiles-1));
            yE[i]=sqrt(yE[i]);
            y[i]=yE[i]/y[i];
            yE[i]=0;
        }

        gBackground = new TGraphErrors(nPoints,x,y,xE,yE);
        gBackground->SetMarkerStyle(20);
        gBackground->SetMarkerSize(0.9);
        gBackground->SetMarkerColor(kBlack);
        gBackground->SetLineColor(kBlack);
        gBackground->GetYaxis()->SetTitle("#varepsilon average value error");
        gBackground->GetYaxis()->SetTitleOffset(1.55);
        gBackground->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
        gBackground->SetTitle("");
        gBackground->GetYaxis()->SetRangeUser(0,0.2);

        TCanvas *out2 = new TCanvas("out2", "out2", 1000, 1000);
        gPad->SetLeftMargin(0.15);
        gBackground->Draw("ap");
        legend2->AddEntry(gBackground, grNames[l], "pl");
        imgName= "finalAnalysis/efficiencyRatio/";
        imgName+=grNames[l];
        imgName+=".avError.png";

        legend2->Draw("same");

        out2->SaveAs(imgName);
    }


    //----------------------------------ratios------------------------------------
    int baseNumberFile = 0;
    TH1F* histo[nInputFiles][nGrNames];
    TCanvas* out[nInputFiles];
    const int nRatioPlots = nInputFiles-1;
    TRatioPlot* ratios[nInputFiles][nGrNames];

    for (int n = 0; n < nGrNames; ++n) { //gr type
        TLegend *legend1 = new TLegend(0.13, 0.85, 0.334, 0.9); //0.1,0.7,0.48,0.9
        legend1->SetFillStyle(0);
        legend1->SetLineColor(0);
        legend1->SetTextSize(0.03);

        Float_t max=0;
        Float_t min=9999;

        for (int j = 0; j < nInputFiles; ++j) { //file type
            histo[j][n] = graphToHisto(graphsSIM[n][j]);
            histo[j][n]->SetName(Form("%s_%s", grNames[n].Data(), input[j].Data()));
            histo[j][n]->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
            if (max<histo[j][n]->GetMaximum()) max=histo[j][n]->GetMaximum();
            if (min>histo[j][n]->GetMinimum()) min=histo[j][n]->GetMinimum();
            legend1->AddEntry(histo[j][n], legendNames[j], "pl");
        }

        for (int i = 0; i < nInputFiles; ++i) {
//            if (i==baseNumberFile) continue;
            out[i] = new TCanvas(Form("%i",i), Form("%i",i), 900, 1000);
            out[i]->cd();
            ratios[i][n] = new TRatioPlot(histo[i][n], histo[baseNumberFile][n],"divsym");
//            ratios[i][n] = new TRatioPlot(histo[i][n], histo[baseNumberFile][n]);
            ratios[i][n]->SetH1DrawOpt("E");
            ratios[i][n]->SetH2DrawOpt("E");
            ratios[i][n]->Draw();
            std::vector<double> lines = {0.5, 1};
            ratios[i][n]->SetGridlines(lines);
            ratios[i][n]->SetSeparationMargin(0.01);
            ratios[i][n]->GetUpperPad()->SetLogy();

            ratios[i][n] -> GetLowerRefYaxis() -> CenterTitle(kTRUE);
            ratios[i][n] -> GetLowerRefXaxis() -> CenterTitle(kTRUE);
            ratios[i][n] -> GetUpperRefYaxis() -> CenterTitle(kTRUE);
            ratios[i][n] -> GetUpperRefYaxis() -> SetRangeUser(0.7*min, 1.4*max);
//            ratios[i][n] -> GetUpperRefYaxis() -> SetRangeUser(0.0001, 1.2*max);

            ratios[i][n] -> GetLowerRefYaxis() -> SetTitle(Form("%s/%s",legendNames[i].Data(),legendNames[baseNumberFile].Data()));
            ratios[i][n] -> GetLowerRefGraph() -> SetMinimum(0.001); //if this is 0.0, problems with X axis
            ratios[i][n] -> GetLowerRefGraph() -> SetMaximum(1.4);

//                ratio -> GetLowerRefXaxis() -> SetTitleSize(0.04);
            ratios[i][n] -> GetLowerRefXaxis() -> SetTitleOffset(1.1);
            ratios[i][n] -> GetUpperRefYaxis() -> SetTitle("#varepsilon");

            out[i]->Update();

//            ratios[i]->GetLowerPad()->SetLeftMargin(0.15);
//            ratios[i]->GetLowerPad()->SetRightMargin(0.05);
//            ratios[i][n]->GetUpperPad()->cd();
//            ratios[i]->GetUpperPad()->SetLeftMargin(0.15);
//            ratios[i]->GetUpperPad()->SetRightMargin(0.05);

            legend1->Draw("same");
            if(i!=baseNumberFile) {
                TGraph *g = ratios[i][n]->GetLowerRefGraph();
                fSimEfficiency->cd();
                g->Write(Form("%s_ratio_%s_to_%s", grNames[n].Data(), input[i].Data(), input[baseNumberFile].Data()));
            }
        }

        for (int k = 0; k < nInputFiles; ++k) {
            if (k!=baseNumberFile) {
                TGraph *g = ratios[k][n]->GetLowerRefGraph();
                TH1F *hUp = (TH1F *) ratios[k][n]->GetUpperRefObject();

//                ratios[baseNumberFile][n]->GetLowerPad()->cd();
//                g->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
//                g->Draw("psame");
////            ratios[0]->GetLowerRefXaxis() -> SetTitle("D^{0} p_{T} (GeV/c)");
//                ratios[baseNumberFile][n]->GetUpperPad()->cd();
//                hUp->DrawCopy("same");
//            out[0]->Update();
                out[k]->Close();
            }
        }

        TString nameRatio = "finalAnalysis/efficiencyRatio/ratio."+grNames[n]+".png";
        out[baseNumberFile]->SaveAs(nameRatio);
        nameRatio = "finalAnalysis/efficiencyRatio/ratio."+grNames[n]+".eps";
        out[baseNumberFile]->SaveAs(nameRatio);
    }

//    fSimEfficiency->Close();
}

//________________________________________________________________
void efficiencyCompareResolution() {
    folderDate="test";

    Int_t colors[] = {1, 46, 9, 8, 6, 40, 42, 28, 2};
    Int_t markers[] = {25, 20, 34, 9, 40, 41, 42, 28, 2};
    gROOT->ProcessLine(".L analyse/FitD0Peak.cpp+");

//    TString input[] = {"ntpTMVA_full_D0.toyMc.1003.105perc.nHitsFit15.root"};
//    TString input[] = {"ntpTMVA_full_D0.toyMc.0303.nHitsFit15.root", "ntpTMVA_full_D0.toyMc.0303.95perc.nHitsFit15.root", "ntpTMVA_full_D0.toyMc.1003.105perc.root"};
    TString input[] = {"ntpTMVA_full_D0.toyMc.0303.nHitsFit15.root", "ntpTMVA_full_D0.toyMc.0303.95perc.nHitsFit15.root", "ntpTMVA_full_D0.toyMc.1003.105perc.nHitsFit15.root"};
    const int nInputFiles = sizeof(input) / sizeof(TString);
    TString legendNamesPlot[nInputFiles] = {"resolution #times 1", "resolution #times 0.95", "resolution #times 1.05"};

    cout<<"number of input files: "<<nInputFiles<<endl;

    TString grNames[] = {"grTpcAccHftPidPreCutsBDT", "grTpcAcc", "grTpcAccHftPid", "grTpcAccHft", "grTpcAccHftPidPreCuts"};
    const int nGrNames = sizeof(grNames) / sizeof(TString);

    TGraphErrors *gr = new TGraphErrors();
    std::vector < TGraphErrors * > graphsSIM[nGrNames];
//    projectFilesBig(input, nInputFiles);
//
    TString ffImg = "resolution";
    plot(input, nInputFiles, legendNamesPlot, grNames, nGrNames);
    gSystem->Exec(Form("mkdir -p calculation/%s/.", ffImg.Data()));
    gSystem->Exec(Form("mv calculation/*.png calculation/%s/.", ffImg.Data()));


}