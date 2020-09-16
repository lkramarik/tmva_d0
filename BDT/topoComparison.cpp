#include "topoComparison.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TDatime.h"
#include "TSystem.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TRatioPlot.h"
#include "TPad.h"
#include "TFile.h"
#include <iostream>

ClassImp(topoComparison)

topoComparison::topoComparison() :
                   mCuts(), mCutsForTuples(), mVars(), mTuples(), mFolder(), mText(""), mDrawText(false), mWeightExpr(){
    cout<<"ahoj"<<endl;
    mOutputFileName="out.root";

    auto* time=new TDatime();
    Int_t day=time->GetDay();
    Int_t month=time->GetMonth();
    Int_t hour=time->GetHour();
    Int_t minute=time->GetMinute();
    TString folderDate=Form("%i%i_%i%i", month, day, hour, minute);
    cout<<folderDate<<endl;
    mFolder=Form("topoProjections/%s/", folderDate.Data());
    gSystem->Exec(Form("mkdir topoProjections/%s", folderDate.Data()));
}

//____________________________________________________________________________________________
topoComparison::topoComparison(TString folderName) :
                   mCuts(), mCutsForTuples(), mVars(), mTuples(), mFolder(), mText(""), mDrawText(false), mWeightExpr(){
    mOutputFileName="out.root";

    mFolder=Form("topoProjections/%s/", folderName.Data());
    gSystem->Exec(Form("mkdir topoProjections/%s", folderName.Data()));
}

//____________________________________________________________________________________________
void topoComparison::project() {
    TCut cuts = connectCuts();
    cout<<"Cuts you set are: "<<cuts<<endl;

    TFile* outF = new TFile(mFolder+mOutputFileName ,"RECREATE");

    const int nTuples = mTuples.size();
    TString nameHisto;
    TH1D* mHists[nTuples];

    if (mVars.size()==0 || mTuples.size()==0){
        cout<<"no vars or tuples to project"<<endl;
        return;
    }

    if(mTuples.size()!=mTuplesToSubtract.size()) {
        cout<<"bckg and signal tuples not loaded correctly"<<endl;
        return;
    }


    for(unsigned int iVar = 0; iVar < mVars.size(); ++iVar) {
        auto *c = new TCanvas("c","c",1100,1100);
        auto *legend = new TLegend(0.5, 0.70, 0.75, 0.89);
        gPad->SetLeftMargin(0.15);

        for(unsigned int iNt = 0; iNt < mTuples.size(); ++iNt) {
            TCut cutToProject = cuts+mCutsForTuples[iNt];
            cout<<"you are projecting tuple "<<mTuples[iNt]->GetName()<<endl;
            cout<<"        with cuts "<<cutToProject<<endl;
            nameHisto = mVars[iVar] + "_" + mTuples[iNt]->GetName();
            cout << nameHisto << endl;
            mHists[iNt] = new TH1D(nameHisto, nameHisto, mnBins[iVar], mLimsMin[iVar], mLimsMax[iVar]);
            mHists[iNt]->SetNameTitle(nameHisto, "");
            mHists[iNt]->GetXaxis()->SetTitle(mTitleX[iVar]);
            mHists[iNt]->GetYaxis()->SetTitle("Scaled counts");
            mTuples[iNt]->Project(nameHisto, mVars[iVar], cutToProject*mWeightExpr[iNt]);
            if (mTuplesToSubtract[iNt]){
                nameHisto+="_bckg";
                TH1F* backgHisto = new TH1F(nameHisto, nameHisto, mnBins[iVar], mLimsMin[iVar], mLimsMax[iVar]);
                mTuplesToSubtract[iNt]->Project(nameHisto, mVars[iVar], cutToProject*mWeightExpr[iNt]);
                mHists[iNt]->Add(backgHisto,-1);
            }

            setHisto(mHists[iNt], iNt);
            mHists[iNt]->Write();
            cout<<"---------------- ENTRIES: "<<mHists[iNt]->GetEntries()<<endl;
            cout<<"---------------- INT: "<<mHists[iNt]->Integral()<<endl;
            mHists[iNt]->Scale(10. / mHists[iNt]->Integral());
//            mHists[iNt]->Scale(1 / mHists[iNt]->GetEntries());
            legend->AddEntry(mHists[iNt], mTuples[iNt]->GetName(), "pl");

        }
//        float maximumY = findMaximumHistos(*mHists, mTuples.size());

        legend->SetFillStyle(0);
        legend->SetLineColor(0);
        legend->SetTextSize(0.03);

        float max = 0, min = 99999999;
        for (unsigned int i1 = 0; i1 < mTuples.size(); ++i1) {
            if (mHists[i1]->GetMaximum()>max) max = mHists[i1]->GetMaximum();
            if (mHists[i1]->GetMinimum()<min) min = mHists[i1]->GetMinimum();
        }

        if (mTuples.size()==2) {
            auto *cRatio = new TCanvas("cRatio","cRatio",900,1100);
            gPad->SetLeftMargin(0.2);

            auto ratioRefMult = new TRatioPlot(mHists[1], mHists[0], "divsym");
            ratioRefMult->SetH1DrawOpt("E");
            ratioRefMult->SetH2DrawOpt("E");
            ratioRefMult->SetSeparationMargin(0.01);
            ratioRefMult -> Draw();

            Double_t lowerPadMin=0.002;
            Double_t lowerPadMax=2.5;
//            ratioRefMult->SetGridlines(lines);

            TString titleUp =mTuples[1]->GetName();
            TString titleDown =mTuples[0]->GetName();
            TString titleRatio = Form("#frac{%s}{%s}",titleUp.Data(),titleDown.Data());

            ratioRefMult -> GetLowerRefYaxis() -> SetTitle(titleRatio);
            ratioRefMult -> GetLowerRefGraph() -> SetMinimum(lowerPadMin);
            ratioRefMult -> GetLowerRefGraph() -> SetMaximum(lowerPadMax);
            ratioRefMult -> GetLowerRefYaxis() -> CenterTitle(kTRUE);
            ratioRefMult -> GetLowerRefGraph() -> SetLineColor(1);
            ratioRefMult -> GetLowerRefGraph() -> SetMarkerColor(1);
            ratioRefMult -> GetLowerRefGraph() -> SetMarkerStyle(1);

            ratioRefMult -> GetUpperRefYaxis() -> SetRangeUser(0.00001, max*1.1);
            if (mLogy[iVar]==1) ratioRefMult -> GetUpperRefYaxis() -> SetRangeUser(0.9*min, max*5.0);


            ratioRefMult->GetLowerPad()->SetLeftMargin(0.15);
            ratioRefMult->GetLowerPad()->SetRightMargin(0.05);

            ratioRefMult -> GetLowerRefXaxis() -> SetTitle(mTitleX[iVar]);
            ratioRefMult -> GetLowerRefXaxis() -> SetTitleSize(0.04);
            ratioRefMult -> GetLowerRefXaxis() -> SetTitleOffset(1.);
            ratioRefMult->GetUpperPad()->SetLogy(mLogy[iVar]);

//                        cRatio->Update();

            ratioRefMult->GetUpperPad()->cd();

            if(mDrawText){
                TPaveText *textRatioPlot = new TPaveText(0.76, 0.925, 0.9, 0.945, "brNDC");
                textRatioPlot->SetTextSize(0.035);
                textRatioPlot->SetLineColor(0);
                textRatioPlot->SetShadowColor(0);
                textRatioPlot->SetFillColor(0);
                textRatioPlot->SetTextFont(42);
                textRatioPlot->AddText(mText);
                textRatioPlot->Draw("same");
            }


            legend->Draw("same");
            cRatio->SaveAs(mFolder+mVars[iVar]+".ratio.png");
            cRatio->Close();
        }

        for (unsigned int i = 0; i < mTuples.size(); ++i) {
            mHists[i]->GetYaxis()->SetRangeUser(0.9*min, max*1.1);
            if (i==0) mHists[i]->Draw(); else mHists[i]->Draw("same");
        }

        gPad->SetLogy(mLogy[iVar]);

        if(mDrawText){
            TPaveText *text5 = new TPaveText(0.724, 0.925, 0.793, 0.945, "brNDC");
            text5->SetTextSize(0.03);
            text5->SetLineColor(0);
            text5->SetShadowColor(0);
            text5->SetFillColor(0);
            text5->SetTextFont(42);
            text5->AddText(mText);
//    mg->GetYaxis()->SetRangeUser(0.,20.);
            text5->Draw("same");
        }

        legend->Draw("same");
        c->SaveAs(mFolder+mVars[iVar]+".png");
        c->Clear();
        c->Close();
        delete c;
        delete legend;
    }
    outF->Close();
}

//____________________________________________________________________________________________
void topoComparison::setHisto(TH1D* histo, int iterator){
    Int_t colors[] = {1, 46, 8, 9, 7, 42};
    histo->SetMarkerStyle(20);
    histo->SetMarkerSize(2);
    histo->Sumw2();
    histo->SetStats(0);
    histo->SetMarkerColor(colors[iterator]);
    histo->SetLineColor(colors[iterator]);
}

//___________________________________________________________________________________________
float topoComparison::findMaximumHistos(TH1F* hists, int nhists){
//    int nHists = sizeof(hists) / sizeof(TH1F);
    float max = 0;
    for (int i = 0; i < nhists; ++i) {
        if (hists[i].GetMaximum()>max) max = hists[i].GetMaximum();
    }
    return max;
}

//____________________________________________________________________________________________
TCut topoComparison::connectCuts() {
    TCut setCuts = "";
    for(unsigned int k = 0; k < mCuts.size(); ++k) {
        setCuts += mCuts[k];
    }
    return setCuts;
}

//____________________________________________________________________________________________
void topoComparison::setFolderName(TString name) {
    mFolder+=name;
    gSystem->Exec(Form("mkdir %s", mFolder.Data()));
    return;
}

//____________________________________________________________________________________________
topoComparison::~topoComparison(){
    mCuts.clear();
    mCuts.shrink_to_fit();
    mCutsForTuples.clear();
    mCutsForTuples.shrink_to_fit();
    mTuples.clear();
    mTuples.shrink_to_fit();
    mTuplesToSubtract.clear();
    mTuplesToSubtract.shrink_to_fit();
    mVars.clear();
    mVars.shrink_to_fit();
    mLimsMin.clear();
    mLimsMin.shrink_to_fit();
    mLimsMax.clear();
    mLimsMax.shrink_to_fit();
    mnBins.clear();
    mnBins.shrink_to_fit();
    mTitleX.clear();
    mTitleX.shrink_to_fit();
    mWeightExpr.clear();
    mWeightExpr.shrink_to_fit();
}
